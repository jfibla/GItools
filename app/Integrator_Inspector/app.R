# ==============================
# app.R (Integrator_Inspector)
# ==============================

options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)
library(readr)
library(dplyr)
library(DT)
library(tibble)
library(stringr)
library(tidyr)
library(purrr)
library(shinycssloaders)
library(ComplexUpset)
library(ggplot2)
library(UpSetR)


# Gene model (hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(IRanges)

# Optional labels
library(AnnotationDbi)
library(org.Hs.eg.db)

library(igraph)
library(plotly)
library(clusterProfiler)
# Per PPI:
#library(STRINGdb)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tags <- shiny::tags
HTML <- shiny::HTML
div  <- shiny::div

# -----------------------------
# Portable config (GItools)
# -----------------------------
cfg_file <- normalizePath(file.path("..", "..", "config.R"), winslash = "/", mustWork = FALSE)
if (!file.exists(cfg_file)) {
  cfg_file <- normalizePath("config.R", winslash = "/", mustWork = FALSE)
}
stopifnot(file.exists(cfg_file))
source(cfg_file, local = TRUE)

source("R/integrator_candidates_helpers.R")

source("R/mod_gwas_hit_priority.R")

source("R/module_audit_table.R")

cfg <- gi_cfg()
gi_root        <- cfg$root
gi_shared_root <- cfg$shared
gi_res_root    <- cfg$resources
gi_pop_dir     <- cfg$pop_dir

Sys.setenv(
  GITOOLS_ROOT      = gi_root,
  GITOOLS_SHARED    = gi_shared_root,
  GITOOLS_RESOURCES = gi_res_root,
  GITOOLS_POP_DIR   = gi_pop_dir
)

cat("[Integrator_Inspector] gi_pop_dir =", gi_pop_dir, "\n")
cat("[Integrator_Inspector] dir.exists(gi_pop_dir) =", dir.exists(gi_pop_dir), "\n")

# -----------------------------
# Shared modules
# -----------------------------
dl_file <- file.path(gi_shared_root, "GItools_local_deeplinks_ALL_IN_ONE.R")
if (file.exists(dl_file)) source(dl_file, local = TRUE)

stopifnot(nzchar(gi_pop_dir))
stopifnot(dir.exists(gi_pop_dir))
stopifnot(length(list.files(gi_pop_dir, pattern = "\\.txt$", full.names = TRUE)) > 0)

ld_file <- file.path(gi_shared_root, "mod_ld_integrator.R")
stopifnot(file.exists(ld_file))
source(ld_file, local = globalenv())

stopifnot(exists("ld_integrator_module_ui"), exists("ld_integrator_module_server"))

# -----------------------------
# Robust file reader
# -----------------------------
guess_delim_by_ext <- function(filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext == "csv") return(",")
  if (ext == "tsv") return("\t")
  NULL
}

detect_delim_from_line <- function(line) {
  cands <- c("\t", ",", ";", "|")
  counts <- vapply(cands, function(d) {
    if (length(line) == 0 || is.na(line) || nchar(line) == 0) return(0L)
    length(strsplit(line, d, fixed = TRUE)[[1]]) - 1L
  }, integer(1))
  cands[which.max(counts)]
}

read_any_delim <- function(path, name, header = TRUE, user_delim = NULL) {
  delim <- user_delim
  if (is.null(delim) || !nzchar(delim)) delim <- guess_delim_by_ext(name)
  if (is.null(delim) || !nzchar(delim)) {
    first_line <- readLines(path, n = 1, warn = FALSE)
    delim <- detect_delim_from_line(first_line)
  }
  
  readr::read_delim(
    file = path,
    delim = delim,
    col_names = isTRUE(header),
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE
  )
}

read_multi_files <- function(files, header = TRUE, user_delim = NULL, label = "file") {
  req(files)
  req(nrow(files) >= 1)
  
  lst <- lapply(seq_len(nrow(files)), function(i) {
    df <- read_any_delim(
      path       = files$datapath[i],
      name       = files$name[i],
      header     = isTRUE(header),
      user_delim = user_delim
    )
    validate(need(is.data.frame(df) && nrow(df) > 0, paste0("Empty ", label, ": ", files$name[i])))
    df
  })
  
  dplyr::bind_rows(lst)
}

# -----------------------------
# Helpers
# -----------------------------
pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

chr_map_plink19 <- function(x) {
  x <- toupper(as.character(x))
  x <- trimws(x)
  x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT", "M", "MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

normalize_trait_term <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- tolower(x)
  
  # eliminar tot el que queda a partir d'un bloc tipus [123456] o [MIM:...]
  x <- gsub("\\[[^\\]]*$", "", x, perl = TRUE)
  x <- gsub("\\[[^\\]]*\\]", "", x, perl = TRUE)
  
  # neteja general
  x <- gsub("[[:space:]]+", " ", x)
  x <- gsub("[^[:alnum:] /_\\-]", "", x)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)
  
  # harmonitzacions útils
  x <- gsub("\\btype 2 diabetes mellitus\\b", "type 2 diabetes", x)
  x <- gsub("\\bt2d\\b", "type 2 diabetes", x)
  x <- gsub("\\bcoronary artery disease\\b", "coronary disease", x)
  x <- gsub("\\bb-cell acute lymphoblastic leukaemia\\b", "acute lymphoblastic leukemia", x)
  x <- gsub("\\bacute lymphoblastic leukemia in childhood b cell precursor\\b", "acute lymphoblastic leukemia", x)
  x <- gsub("\\bacute lymphoblastic leukemia childhood\\b", "acute lymphoblastic leukemia", x)
  
  # placeholders / soroll
  bad_vals <- c(
    "", "na", "n/a", "none", "unknown", "unspecified",
    "not specified", "not_specified", "not provided", "not_provided",
    "disease", "diseases", "trait", "traits", "phenotype", "phenotypes",
    "clinical significance", "association", "risk factor"
  )
  x[x %in% bad_vals] <- NA_character_
  
  # eliminar termes que encara acaben amb PMID enganxat
  x[grepl(".*[0-9]{6,}$", x)] <- NA_character_
  
  # eliminar termes numèrics purs
  x[grepl("^[0-9]+$", x)] <- NA_character_
  
  # massa curt = soroll
  x[nchar(x) < 2] <- NA_character_
  
  x
}


sanitize_dt_types <- function(df) {
  if (!is.data.frame(df)) return(df)
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  df
}

safe_top1 <- function(x) {
  x <- unique(stats::na.omit(trimws(as.character(x))))
  if (!length(x)) NA_character_ else x[1]
}

integrator_root <- file.path(gi_shared_root, "integrator_exports")
dir.create(integrator_root, recursive = TRUE, showWarnings = FALSE)

list_integrator_sessions <- function(integrator_root) {
  if (!dir.exists(integrator_root)) return(character(0))
  
  dd <- list.dirs(integrator_root, full.names = TRUE, recursive = FALSE)
  if (!length(dd)) return(character(0))
  
  dd <- dd[file.exists(file.path(dd, "manifest.rds"))]
  if (!length(dd)) return(character(0))
  
  info <- file.info(dd)
  dd <- dd[order(info$mtime, decreasing = TRUE)]
  names(dd) <- basename(dd)
  dd
}

read_manifest_safe <- function(session_dir) {
  if (is.null(session_dir) || !nzchar(session_dir)) return(NULL)
  
  mf <- file.path(session_dir, "manifest.rds")
  if (!file.exists(mf)) return(NULL)
  
  x <- tryCatch(readRDS(mf), error = function(e) NULL)
  if (!is.list(x)) return(NULL)
  
  if (is.null(x$files)) x$files <- list()
  if (is.null(x$apps_present)) x$apps_present <- character(0)
  x
}

manifest_to_display_df <- function(manifest) {
  if (is.null(manifest) || !is.list(manifest)) {
    return(tibble::tibble(
      field = c("session_id", "gwas_session_file", "cluster_method", "threshold_used", "apps_present", "created_at", "last_updated"),
      value = c(NA, NA, NA, NA, NA, NA, NA)
    ))
  }
  
  tibble::tibble(
    field = c(
      "session_id",
      "gwas_session_file",
      "cluster_method",
      "threshold_used",
      "apps_present",
      "created_at",
      "last_updated"
    ),
    value = c(
      paste(manifest$session_id %||% "", collapse = " | "),
      paste(as.character(manifest$gwas_session_file %||% ""), collapse = " | "),
      paste(as.character(manifest$cluster_method %||% ""), collapse = " | "),
      paste(as.character(manifest$threshold_used %||% ""), collapse = " | "),
      paste(as.character(manifest$apps_present %||% ""), collapse = "; "),
      paste(as.character(manifest$created_at %||% ""), collapse = " | "),
      paste(as.character(manifest$last_updated %||% ""), collapse = " | ")
    )
  )
}

validate_manifest_basic <- function(manifest) {
  if (is.null(manifest) || !is.list(manifest)) return("manifest missing or unreadable")
  
  req_fields <- c("session_id", "cluster_method", "threshold_used", "gwas_session_file")
  miss <- req_fields[!req_fields %in% names(manifest)]
  if (length(miss)) {
    return(paste("manifest missing fields:", paste(miss, collapse = ", ")))
  }
  
  if (is.null(manifest$files) || !is.list(manifest$files)) {
    return("manifest$files missing or invalid")
  }
  
  "ok"
}



empty_gene_bridge <- function() {
  tibble::tibble(
    gene = character(),
    source_app = character(),
    evidence_type = character(),
    cluster_id = character(),
    chr = integer(),
    start = integer(),
    end = integer()
  )
}

empty_term_bridge <- function() {
  tibble::tibble(
    term = character(),
    term_type = character(),
    source_app = character(),
    evidence_type = character(),
    cluster_id = character(),
    chr = integer(),
    start = integer(),
    end = integer()
  )
}

# clean terms
clean_term_bridge_records <- function(df) {
  
  if (!is.data.frame(df) || !nrow(df)) {
    return(empty_term_bridge())
  }
  
  out <- df %>%
    dplyr::mutate(
      term = clean_term_label_strict(term),
      term_type = trimws(as.character(term_type)),
      source_app = trimws(tolower(as.character(source_app))),
      evidence_type = trimws(tolower(as.character(evidence_type))),
      cluster_id = trimws(as.character(cluster_id)),
      chr = as.character(chr),
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end))
    ) %>%
    dplyr::filter(
      !is.na(term),
      nzchar(term)
    ) %>%
    dplyr::distinct()
  
  out
}

clean_term_label_strict <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # elimina fragments tipus [27798627 o [27798627]
  x <- gsub("\\[[0-9]+\\]?$", "", x)
  
  # elimina termes que siguin només números o números + ]
  x <- gsub("^[0-9]+\\]?$", "", x)
  
  # neteja espais sobrants
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  x
}

normalize_term_key <- function(x) {
  x <- clean_term_label_strict(x)
  x <- tolower(x)
  x <- trimws(x)
  x
}

################

collapse_unique_terms <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- unique(x)
  x <- sort(x)
  
  if (!length(x)) return("")
  
  paste(x, collapse = "; ")
}

make_terms_collapsed_links <- function(x, base_url) {
  vapply(x, function(cell) {
    if (is.null(cell) || is.na(cell) || !nzchar(trimws(cell))) return("")
    
    terms <- unlist(strsplit(as.character(cell), "\\s*;\\s*"))
    terms <- trimws(terms)
    terms <- terms[!is.na(terms) & nzchar(terms)]
    terms <- unique(terms)
    
    if (!length(terms)) return("")
    
    linked <- vapply(terms, function(term) {
      safe_term <- utils::URLencode(term, reserved = TRUE)
      paste0('<a href="', base_url, safe_term, '" target="_blank">', term, '</a>')
    }, character(1))
    
    linked <- linked[!is.na(linked) & nzchar(linked)]
    
    if (!length(linked)) return("")
    
    paste(linked, collapse = "; ")
  }, character(1))
}

make_omim_collapsed_links <- function(x) {
  make_terms_collapsed_links(
    x = x,
    base_url = "https://www.omim.org/search/?search="
  )
}

read_bridge_rds <- function(path, empty_df) {
  if (!file.exists(path)) return(empty_df)
  
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.data.frame(x)) return(empty_df)
  
  x
}



read_existing_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.data.frame(x) || !nrow(x)) return(NULL)
  x
}

integrated_clusters_from_rds <- function(selected_session_dir) {
  paths <- c(
    file.path(selected_session_dir, "catalog_clusters_master.rds"),
    file.path(selected_session_dir, "gtex_clusters_master.rds"),
    file.path(selected_session_dir, "nonsyn_clusters_master.rds"),
    file.path(selected_session_dir, "ewastum_clusters_master.rds"),
    file.path(selected_session_dir, "ewasdis_clusters_master.rds")
  )
  
  lst <- Filter(Negate(is.null), lapply(paths, read_existing_rds))
  if (!length(lst)) {
    return(data.frame(
      cluster_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  df <- dplyr::bind_rows(lst)
  
  chr_col <- pick_col(df, c("chr", "CHR", "cluster_chr_n", "chrom", "CHROM", "chromosome"))
  st_col  <- pick_col(df, c("cluster_start", "start", "START", "start_bp", "FROM", "from", "bp1"))
  en_col  <- pick_col(df, c("cluster_end", "end", "END", "end_bp", "TO", "to", "bp2"))
  id_col  <- pick_col(df, c("cluster_id", "CLUSTER_ID", "cluster", "id", "ID", "label"))
  
  if (is.null(chr_col) || is.null(st_col) || is.null(en_col)) {
    stop("Integrated cluster RDSs must contain chr + start + end.")
  }
  
  df %>%
    dplyr::transmute(
      cluster_id = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_,
      chr   = trimws(as.character(.data[[chr_col]])),
      start = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
      end   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]]))))
    ) %>%
    dplyr::filter(!is.na(chr), nzchar(chr), !is.na(start), !is.na(end)) %>%
    dplyr::mutate(
      cluster_id = trimws(cluster_id),
      cluster_id = ifelse(
        is.na(cluster_id) | !nzchar(cluster_id),
        paste0("cluster_", dplyr::row_number()),
        cluster_id
      ),
      start2 = pmin(start, end),
      end2   = pmax(start, end),
      start = start2,
      end   = end2
    ) %>%
    dplyr::select(cluster_id, chr, start, end) %>%
    dplyr::distinct()
}

integrated_candidates_from_rds <- function(selected_session_dir) {
  files <- c(
    catalog  = file.path(selected_session_dir, "catalog_candidates.rds"),
    gtex     = file.path(selected_session_dir, "gtex_candidates.rds"),
    nonsyn   = file.path(selected_session_dir, "nonsyn_candidates.rds"),
    ewastum  = file.path(selected_session_dir, "ewastum_candidates.rds"),
    ewasdis  = file.path(selected_session_dir, "ewasdis_candidates.rds")
  )
  
  lst <- lapply(names(files), function(app_nm) {
    x <- read_existing_rds(files[[app_nm]])
    if (is.null(x) || !is.data.frame(x) || !nrow(x)) return(NULL)
    x$source_app <- app_nm
    x
  })
  
  lst <- Filter(Negate(is.null), lst)
  
  if (!length(lst)) {
    return(data.frame(
      cluster_id = character(),
      chr = character(),
      position = integer(),
      rsid = character(),
      id_hit = character(),
      classe = character(),
      source_app = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  df <- dplyr::bind_rows(lst)
  
  id_cluster_col <- pick_col(df, c("cluster_id", "CLUSTER_ID", "cluster", "id"))
  chr_col <- pick_col(df, c("chr", "CHR", "chrom", "CHROM", "chromosome"))
  pos_col <- pick_col(df, c("position", "pos_ini", "POS", "pos", "BP", "bp", "start", "pos_start"))
  rsid_col <- pick_col(df, c("rsid", "RSID", "SNP", "snp", "marker", "MARKER"))
  id_col  <- pick_col(df, c("id_hit", "id", "ID"))
  cl_col  <- pick_col(df, c("classe", "class", "CLASS", "type", "TYPE"))
  src_col <- pick_col(df, c("source_app"))
  
  if (is.null(chr_col) || is.null(pos_col)) {
    stop("Integrated candidate RDSs must contain chr + position.")
  }
  
  out <- df %>%
    dplyr::transmute(
      cluster_id = if (!is.null(id_cluster_col)) trimws(as.character(.data[[id_cluster_col]])) else NA_character_,
      chr = trimws(as.character(.data[[chr_col]])),
      position = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
      rsid = if (!is.null(rsid_col)) trimws(as.character(.data[[rsid_col]])) else NA_character_,
      id_hit = if (!is.null(id_col)) trimws(as.character(.data[[id_col]])) else NA_character_,
      classe = if (!is.null(cl_col)) normalize_classe(trimws(as.character(.data[[cl_col]]))) else NA_character_,
      source_app = if (!is.null(src_col)) trimws(as.character(.data[[src_col]])) else NA_character_
    ) %>%
    dplyr::mutate(
      rsid = dplyr::coalesce(dplyr::na_if(rsid, ""), id_hit),
      id_hit = dplyr::coalesce(dplyr::na_if(id_hit, ""), rsid),
      chr = trimws(as.character(chr))
    ) %>%
    dplyr::filter(
      !is.na(chr), nzchar(chr),
      !is.na(position),
      !is.na(rsid), nzchar(rsid)
    ) %>%
    dplyr::distinct(
      cluster_id, chr, position, rsid, id_hit, classe, source_app,
      .keep_all = TRUE
    )
  
  out
}

read_keep_file_global <- function(pop, ld_pops_dir) {
  d <- tryCatch(normalizePath(ld_pops_dir, winslash = "/", mustWork = FALSE),
                error = function(e) ld_pops_dir)
  if (is.null(d) || !nzchar(d) || !dir.exists(d)) {
    stop("POP keep-files directory not found: ", ld_pops_dir)
  }
  
  pop <- as.character(pop %||% "")
  if (!nzchar(pop)) {
    av <- available_pops_dir(d)
    stop("Population not selected. Available: ", paste(av, collapse = ", "))
  }
  
  cand <- c(
    file.path(d, paste0(pop, ".txt")),
    file.path(d, paste0(pop, ".TXT"))
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f) || !nzchar(f)) {
    av <- available_pops_dir(d)
    stop("Missing keep file for pop='", pop, "' in ", d,
         ". Available: ", paste(av, collapse = ", "))
  }
  
  x <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 1) stop("Keep file must have at least 1 column (IID)")
  if (ncol(x) == 1) x <- cbind(x[,1], x[,1]) else x <- x[,1:2]
  colnames(x) <- c("FID","IID")
  
  out <- tempfile(pattern = paste0("keep_", pop, "_"), fileext = ".txt")
  write.table(x, out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  out
}

compute_ld_for_cluster_global <- function(
    cluster_row,
    candidates_df,
    bfile_ref,
    keep_path,
    plink_bin,
    workdir,
    pop,
    ld_metric = "R2",
    r2_min = 0.6,
    compute_blocks = TRUE,
    max_snps_interval = 400
) {
  bundle <- compute_ld_bundle_common(
    cluster_row = cluster_row,
    candidates_df = candidates_df,
    bfile_ref = bfile_ref,
    keep_path = keep_path,
    plink_bin = plink_bin,
    workdir = workdir,
    pop = pop,
    ld_metric = ld_metric,
    r2_min = r2_min,
    max_snps_interval = max_snps_interval,
    compute_blocks = compute_blocks,
    append_log = NULL
  )
  
  cid <- as.character(bundle$cluster_id)
  
  if (!is.data.frame(bundle$block_summary)) {
    bundle$block_summary <- tibble::tibble()
  }
  if (!is.data.frame(bundle$block_hits)) {
    bundle$block_hits <- tibble::tibble()
  }
  if (!is.data.frame(bundle$block_genes)) {
    bundle$block_genes <- tibble::tibble()
  }
  if (!is.data.frame(bundle$block_ranges)) {
    bundle$block_ranges <- tibble::tibble()
  }
  if (!is.data.frame(bundle$proxies)) {
    bundle$proxies <- tibble::tibble()
  }
  if (!is.data.frame(bundle$seeds)) {
    bundle$seeds <- tibble::tibble()
  }
  if (!is.data.frame(bundle$gwas_hits)) {
    bundle$gwas_hits <- tibble::tibble()
  }
  if (!is.data.frame(bundle$candidates)) {
    bundle$candidates <- tibble::tibble()
  }
  if (!is.data.frame(bundle$fl)) {
    bundle$fl <- tibble::tibble()
  }
  
  summary_row <- data.frame(
    cluster_id = cid,
    chr = as.character(cluster_row$chr[1]),
    start = as.integer(cluster_row$start[1]),
    end = as.integer(cluster_row$end[1]),
    population = as.character(pop),
    ld_metric = as.character(ld_metric),
    n_interval_snps = nrow(bundle$fl),
    n_candidate_snps = nrow(bundle$candidates),
    n_gwas_hits = nrow(bundle$gwas_hits),
    gwas_hits = if ("rsid" %in% names(bundle$gwas_hits)) {
      paste(unique(stats::na.omit(bundle$gwas_hits$rsid)), collapse = "; ")
    } else {
      ""
    },
    n_seeds_exact = if ("seed_type" %in% names(bundle$seeds)) sum(bundle$seeds$seed_type == "exact", na.rm = TRUE) else 0L,
    n_seeds_nearest = if ("seed_type" %in% names(bundle$seeds)) sum(bundle$seeds$seed_type == "nearest", na.rm = TRUE) else 0L,
    seed_snps = if ("seed_snp" %in% names(bundle$seeds)) {
      paste(unique(stats::na.omit(bundle$seeds$seed_snp)), collapse = "; ")
    } else {
      ""
    },
    n_proxy_snps = if ("proxy_snp" %in% names(bundle$proxies)) dplyr::n_distinct(bundle$proxies$proxy_snp) else 0L,
    proxy_snps = if ("proxy_snp" %in% names(bundle$proxies)) {
      paste(unique(stats::na.omit(bundle$proxies$proxy_snp)), collapse = "; ")
    } else {
      ""
    },
    max_ld_value = if (nrow(bundle$proxies) && "ld_value" %in% names(bundle$proxies)) max(bundle$proxies$ld_value, na.rm = TRUE) else NA_real_,
    mean_ld_value = if (nrow(bundle$proxies) && "ld_value" %in% names(bundle$proxies)) mean(bundle$proxies$ld_value, na.rm = TRUE) else NA_real_,
    n_blocks = if (is.data.frame(bundle$block_ranges)) nrow(bundle$block_ranges) else 0L,
    n_block_hits = if (nrow(bundle$block_summary) && "n_block_hits" %in% names(bundle$block_summary)) {
      sum(bundle$block_summary$n_block_hits, na.rm = TRUE)
    } else {
      0L
    },
    n_block_genes = if (nrow(bundle$block_summary) && "n_genes_in_block" %in% names(bundle$block_summary)) {
      sum(bundle$block_summary$n_genes_in_block, na.rm = TRUE)
    } else {
      0L
    },
    block_apps = if (nrow(bundle$block_summary) && "external_apps" %in% names(bundle$block_summary)) {
      paste(unique(stats::na.omit(bundle$block_summary$external_apps)), collapse = "; ")
    } else {
      ""
    },
    block_genes = if (nrow(bundle$block_summary) && "genes_in_block" %in% names(bundle$block_summary)) {
      paste(unique(stats::na.omit(bundle$block_summary$genes_in_block)), collapse = "; ")
    } else {
      ""
    },
    ld_has_signal = as.integer(nrow(bundle$proxies) > 0),
    ld_computed = TRUE,
    status = "ok",
    stringsAsFactors = FALSE
  )
  
  details <- list(
    cluster_id = cid,
    chr = as.character(cluster_row$chr[1]),
    start = as.integer(cluster_row$start[1]),
    end = as.integer(cluster_row$end[1]),
    population = as.character(pop),
    ld_metric = as.character(ld_metric),
    fl = bundle$fl,
    candidates = bundle$candidates,
    gwas_hits = bundle$gwas_hits,
    seeds = bundle$seeds,
    proxies = bundle$proxies,
    blocks = bundle$blocks,
    blocks_ij = bundle$blocks_ij,
    block_ranges = bundle$block_ranges,
    block_hits = bundle$block_hits,
    block_genes = bundle$block_genes,
    block_summary = bundle$block_summary,
    summary = summary_row
  )
  
  list(summary = summary_row, details = details)
}


run_ld_global_for_clusters <- function(
    clusters_df,
    candidates_df,
    out_summary_rds,
    out_details_rds,
    bfile_ref,
    keep_path,
    plink_bin,
    workdir,
    pop,
    ld_metric = "R2",
    r2_min = 0.6,
    compute_blocks = FALSE,
    max_snps_interval = 400,
    progress_fun = NULL
) {
  if (!is.data.frame(clusters_df) || !nrow(clusters_df)) {
    stop("clusters_df has 0 rows")
  }
  if (!is.data.frame(candidates_df) || !nrow(candidates_df)) {
    stop("candidates_df has 0 rows")
  }
  
  req_cols_cl <- c("cluster_id", "chr", "start", "end")
  req_cols_ca <- c("chr", "position", "rsid", "classe")
  
  miss_cl <- setdiff(req_cols_cl, names(clusters_df))
  miss_ca <- setdiff(req_cols_ca, names(candidates_df))
  
  if (length(miss_cl)) stop("clusters_df missing: ", paste(miss_cl, collapse = ", "))
  if (length(miss_ca)) stop("candidates_df missing: ", paste(miss_ca, collapse = ", "))
  
  clusters_df$chr <- as.character(clusters_df$chr)
  candidates_df$chr <- as.character(candidates_df$chr)
  
  dir.create(dirname(out_summary_rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_details_rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  
  res <- lapply(seq_len(nrow(clusters_df)), function(i) {
    cl <- clusters_df[i, , drop = FALSE]
    cid <- as.character(cl$cluster_id[1])
    
    if (is.function(progress_fun)) {
      progress_fun(i = i, n = nrow(clusters_df), cluster_id = cid, stage = "start")
    }
    
    ans <- tryCatch(
      {
        out_one <- compute_ld_for_cluster_global(
          cluster_row = cl,
          candidates_df = candidates_df,
          bfile_ref = bfile_ref,
          keep_path = keep_path,
          plink_bin = plink_bin,
          workdir = workdir,
          pop = pop,
          ld_metric = ld_metric,
          r2_min = r2_min,
          compute_blocks = compute_blocks,
          max_snps_interval = max_snps_interval
        )
        
        cat("[GLOBAL LD][OK] ", cid,
            " | details names = ",
            paste(names(out_one$details), collapse = ", "),
            "\n", sep = "")
        
        out_one
      },
      error = function(e) {
        msg <- conditionMessage(e)
        cat("[GLOBAL LD][ERROR] ", cid, " -> ", msg, "\n", sep = "")
        
        summary_row <- data.frame(
          cluster_id = cid,
          chr = as.character(cl$chr[1]),
          start = as.integer(cl$start[1]),
          end = as.integer(cl$end[1]),
          population = as.character(pop),
          ld_metric = as.character(ld_metric),
          n_interval_snps = NA_integer_,
          n_candidate_snps = NA_integer_,
          n_gwas_hits = NA_integer_,
          gwas_hits = "",
          n_seeds_exact = 0L,
          n_seeds_nearest = 0L,
          seed_snps = "",
          n_proxy_snps = 0L,
          proxy_snps = "",
          max_ld_value = NA_real_,
          mean_ld_value = NA_real_,
          n_blocks = 0L,
          n_block_hits = 0L,
          n_block_genes = 0L,
          block_apps = "",
          block_genes = "",
          ld_has_signal = 0L,
          ld_computed = FALSE,
          status = paste("error:", msg),
          stringsAsFactors = FALSE
        )
        
        list(summary = summary_row, details = NULL)
      }
    )
    
    if (is.function(progress_fun)) {
      progress_fun(i = i, n = nrow(clusters_df), cluster_id = cid, stage = "done")
    }
    
    ans
  })
  
  summary_df <- dplyr::bind_rows(lapply(res, `[[`, "summary"))
  details_list <- lapply(res, `[[`, "details")
  names(details_list) <- summary_df$cluster_id
  
  saveRDS(summary_df, out_summary_rds)
  saveRDS(details_list, out_details_rds)
  
  cat("[GLOBAL LD] summary status table:\n")
  print(table(summary_df$status, useNA = "ifany"))
  
  invisible(list(
    summary = summary_df,
    details = details_list
  ))
}

# ============================================================
# Block-level annotation helpers
# ============================================================

blocks_ij_to_ranges <- function(blocks_ij, fl, cluster_id = NA_character_, chr = NA_integer_) {
  if (!is.data.frame(blocks_ij) || !nrow(blocks_ij)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      block_id = character(),
      i = integer(),
      j = integer(),
      block_start = integer(),
      block_end = integer(),
      block_size_bp = integer()
    ))
  }
  
  stopifnot(is.data.frame(fl), all(c("BP") %in% names(fl)))
  
  out <- blocks_ij %>%
    dplyr::mutate(
      i = suppressWarnings(as.integer(i)),
      j = suppressWarnings(as.integer(j))
    ) %>%
    dplyr::filter(
      is.finite(i), is.finite(j),
      i >= 1, j >= 1,
      i <= nrow(fl), j <= nrow(fl),
      j > i
    ) %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      block_start = as.integer(fl$BP[i]),
      block_end   = as.integer(fl$BP[j]),
      block_size_bp = as.integer(block_end - block_start),
      block_id = paste0(cluster_id, "_block_", dplyr::row_number())
    ) %>%
    dplyr::select(cluster_id, chr, block_id, i, j, block_start, block_end, block_size_bp)
  
  tibble::as_tibble(out)
}

######## helpers commons

#########


truncate_for_dt <- function(x, n = 80) {
  
  x <- as.character(x)
  
  full <- htmltools::htmlEscape(x)
  
  short <- ifelse(
    nchar(x) > n,
    paste0(substr(x, 1, n), "..."),
    x
  )
  
  short <- htmltools::htmlEscape(short)
  
  ifelse(
    nchar(x) > n,
    paste0(
      "<details class='dt-details'>",
      "<summary><span class='dt-trunc'>", short, "</span></summary>",
      "<div class='dt-fulltext'>", full, "</div>",
      "</details>"
    ),
    full
  )
}


safe_int0 <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x[is.na(x)] <- 0L
  x
}

safe_num0 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

collapse_text_list_html <- function(x, max_visible = 4L, summary_label = NULL) {
  vals <- split_semicolon(x)
  vals <- unique(vals)
  vals <- vals[nzchar(vals)]
  
  if (!length(vals)) return("")
  
  if (length(vals) <= max_visible) {
    return(html_escape_basic(paste(vals, collapse = "; ")))
  }
  
  first_part <- html_escape_basic(paste(vals[seq_len(max_visible)], collapse = "; "))
  rest_part  <- html_escape_basic(paste(vals[(max_visible + 1):length(vals)], collapse = "; "))
  
  summ <- summary_label %||% paste0(length(vals), " items")
  
  paste0(
    "<details>",
    "<summary>", html_escape_basic(summ), "</summary>",
    "<div style='margin-top:6px;'>", first_part,
    if (nzchar(rest_part)) paste0("; ", rest_part) else "",
    "</div>",
    "</details>"
  )
}


read_gwas_significance_bridge <- function(session_dir) {
  f <- file.path(session_dir, "gwas_significance_bridge.rds")
  if (!file.exists(f)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      position = integer(),
      rsid = character(),
      p_value = numeric(),
      logp = numeric()
    ))
  }
  
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.data.frame(x) || !nrow(x)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      position = integer(),
      rsid = character(),
      p_value = numeric(),
      logp = numeric()
    ))
  }
  
  x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      p_value = suppressWarnings(as.numeric(p_value)),
      logp = suppressWarnings(as.numeric(logp))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(position)
    ) %>%
    dplyr::distinct()
}

summarize_gwas_significance_by_block <- function(block_df, gwas_bridge_df) {
  if (!is.data.frame(block_df) || !nrow(block_df) ||
      !is.data.frame(gwas_bridge_df) || !nrow(gwas_bridge_df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  if (!all(c("cluster_id", "chr", "block_id", "block_start", "block_end") %in% names(block_df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  bb <- block_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      block_id = as.character(block_id),
      block_start = suppressWarnings(as.integer(block_start)),
      block_end = suppressWarnings(as.integer(block_end))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(block_start), is.finite(block_end)
    )
  
  gg <- gwas_bridge_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      logp = suppressWarnings(as.numeric(logp))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(position)
    )
  
  ov <- gg %>%
    dplyr::inner_join(
      bb %>%
        dplyr::select(cluster_id, chr, block_id, block_start, block_end),
      by = c("cluster_id", "chr"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      is.finite(position),
      is.finite(block_start),
      is.finite(block_end),
      position >= block_start,
      position <= block_end
    ) %>%
    dplyr::distinct(cluster_id, block_id, rsid, position, .keep_all = TRUE)
  
  if (!nrow(ov)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  collapse_hits <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- x[nzchar(x)]
    paste(sort(x), collapse = "; ")
  }
  
  ov %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      n_gwas_sig_hits = dplyr::n_distinct(rsid[!is.na(rsid) & nzchar(rsid)]),
      max_gwas_logp = if (any(is.finite(logp))) max(logp, na.rm = TRUE) else NA_real_,
      mean_gwas_logp = if (any(is.finite(logp))) mean(logp[is.finite(logp)], na.rm = TRUE) else NA_real_,
      gwas_sig_hits = collapse_hits(rsid),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_gwas_sig_hits = safe_int0(n_gwas_sig_hits),
      max_gwas_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      mean_gwas_logp = suppressWarnings(as.numeric(mean_gwas_logp)),
      gwas_sig_hits = dplyr::coalesce(gwas_sig_hits, "")
    )
}

extract_block_order <- function(block_id, block_label = NULL) {
  x <- dplyr::coalesce(as.character(block_label), as.character(block_id))
  x <- trimws(x)
  
  out <- suppressWarnings(as.integer(sub("^.*B([0-9]+)$", "\\1", x)))
  out[is.na(out)] <- 9999L
  out
}
## Helper shared terms
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
}

split_semicolon_values <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character())
  
  out <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  out <- trimws(out)
  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

build_shared_traits_app_matrix <- function(shared_traits_tbl) {
  if (!is.data.frame(shared_traits_tbl) || !nrow(shared_traits_tbl)) {
    return(tibble::tibble())
  }
  
  need_cols <- c("term_label", "term_norm", "term_type", "apps", "n_apps", "n_clusters")
  miss <- setdiff(need_cols, names(shared_traits_tbl))
  if (length(miss)) {
    stop("Missing columns in shared_traits_tbl: ", paste(miss, collapse = ", "))
  }
  
  app_long <- shared_traits_tbl %>%
    dplyr::mutate(
      term_label = as.character(term_label),
      term_norm = as.character(term_norm),
      term_type = as.character(term_type),
      apps = as.character(apps),
      n_apps = suppressWarnings(as.integer(n_apps)),
      n_clusters = suppressWarnings(as.integer(n_clusters)),
      term_display = paste0(
        dplyr::if_else(!is.na(term_label) & nzchar(term_label), term_label, term_norm),
        " [", dplyr::coalesce(term_type, "unknown"), "]"
      )
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(app_list = list(split_semicolon_values(apps))) %>%
    tidyr::unnest(app_list) %>%
    dplyr::rename(app = app_list) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(app), nzchar(app))
  
  if (!nrow(app_long)) return(tibble::tibble())
  
  app_levels <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
  app_levels <- unique(c(app_levels, sort(unique(app_long$app))))
  
  app_long %>%
    dplyr::mutate(
      app = factor(app, levels = app_levels),
      term_display = as.character(term_display),
      term_type = dplyr::coalesce(term_type, "unknown"),
      n_apps = dplyr::coalesce(n_apps, 0L),
      n_clusters = dplyr::coalesce(n_clusters, 0L)
    ) %>%
    dplyr::distinct(term_norm, term_display, term_type, n_apps, n_clusters, app) %>%
    dplyr::arrange(dplyr::desc(n_apps), dplyr::desc(n_clusters), term_display, app)
}

build_shared_traits_heatmap_df <- function(shared_traits_tbl) {
  app_long <- build_shared_traits_app_matrix(shared_traits_tbl)
  if (!is.data.frame(app_long) || !nrow(app_long)) return(tibble::tibble())
  
  term_order <- app_long %>%
    dplyr::distinct(term_norm, term_display, term_type, n_apps, n_clusters) %>%
    dplyr::arrange(
      dplyr::desc(n_apps),
      dplyr::desc(n_clusters),
      term_type,
      term_display
    ) %>%
    dplyr::pull(term_display)
  
  app_levels <- levels(app_long$app) %||% sort(unique(as.character(app_long$app)))
  
  tidyr::expand_grid(
    term_display = unique(term_order),
    app = app_levels
  ) %>%
    dplyr::left_join(
      app_long %>%
        dplyr::transmute(
          term_display,
          app = as.character(app),
          present = 1L
        ) %>%
        dplyr::distinct(),
      by = c("term_display", "app")
    ) %>%
    dplyr::mutate(
      present = dplyr::coalesce(present, 0L),
      term_display = factor(term_display, levels = rev(unique(term_order))),
      app = factor(app, levels = app_levels)
    )
}

build_shared_traits_upset_df <- function(shared_traits_tbl) {
  app_long <- build_shared_traits_app_matrix(shared_traits_tbl)
  if (!is.data.frame(app_long) || !nrow(app_long)) return(tibble::tibble())
  
  app_levels <- levels(app_long$app) %||% sort(unique(as.character(app_long$app)))
  
  out <- app_long %>%
    dplyr::mutate(app = as.character(app)) %>%
    dplyr::distinct(term_norm, term_display, term_type, app) %>%
    dplyr::mutate(value = 1L) %>%
    tidyr::pivot_wider(
      names_from = app,
      values_from = value,
      values_fill = 0L
    )
  
  for (nm in app_levels) {
    if (!nm %in% names(out)) out[[nm]] <- 0L
  }
  
  out %>%
    dplyr::select(term_norm, term_display, term_type, dplyr::all_of(app_levels))
}

short_term_label <- function(x, max_chars = 45) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x)] <- ""
  
  out <- ifelse(
    nchar(x) > max_chars,
    paste0(substr(x, 1, max_chars - 3), "..."),
    x
  )
  
  out
}

################################################################################
# funcio per ordenar cluster
cluster_levels_order <- function(df) {
  if (!is.data.frame(df) || !nrow(df) || !"cluster_id" %in% names(df)) return(character())
  
  df %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=chr)\\d+"))),
      cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
    ) %>%
    dplyr::distinct(cluster_id, chr_num, cluster_num) %>%
    dplyr::mutate(
      chr_num = dplyr::coalesce(chr_num, 999L),
      cluster_num = dplyr::coalesce(cluster_num, 999L)
    ) %>%
    dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
    dplyr::pull(cluster_id) %>%
    rev()
}

# helper plot
consensus_candidates_apps_legend_ui <- function() {
  
  items <- list(
    list(label = "catalog", color = "#e6ab02"),
    list(label = "gtex", color =    "#f28e2b"),
    list(label = "nonsyn", color =  "#e15759"),
    list(label = "ewasdis", color = "#b24745"),
    list(label = "ewastum", color = "#7b3f98")
  )
  
  tags$div(
    style = "display:flex; flex-wrap:wrap; align-items:center; gap:14px; margin-top:8px; margin-bottom:4px;",
    lapply(items, function(it) {
      tags$div(
        style = "display:flex; align-items:center; gap:6px;",
        tags$span(
          style = paste0(
            "display:inline-block; width:14px; height:14px; ",
            "background:", it$color, "; border:1px solid white;"
          )
        ),
        tags$span(it$label, style = "font-size:12px;")
      )
    })
  )
}
# helper HTML
escape_html_basic <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

make_expandable_genecards_cell <- function(x, max_show = 5L) {
  x <- as.character(x %||% "")
  if (!nzchar(trimws(x))) return("")
  
  genes <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes <- unique(genes)
  
  if (!length(genes)) return("")
  
  links_all <- make_genecards_links(paste(genes, collapse = ";"))
  n_genes <- length(genes)
  
  if (n_genes <= max_show) {
    return(links_all)
  }
  
  genes_short <- genes[seq_len(max_show)]
  links_short <- make_genecards_links(paste(genes_short, collapse = ";"))
  
  htmltools::HTML(paste0(
    "<details>",
    "<summary style='cursor:pointer;'>", n_genes, " genes · View</summary>",
    "<div style='margin-top:6px;'>", links_short,
    "<div style='margin-top:6px;'>", links_all, "</div>",
    "</div>",
    "</details>"
  ))
}
# priority classification by tertiles 
classify_priority_tertiles <- function(score) {
  score <- suppressWarnings(as.numeric(score))
  
  q33 <- stats::quantile(score, 0.33, na.rm = TRUE, names = FALSE, type = 7)
  q66 <- stats::quantile(score, 0.66, na.rm = TRUE, names = FALSE, type = 7)
  
  dplyr::case_when(
    is.na(score) ~ NA_character_,
    q33 == q66 & score >= q66 ~ "High",
    q33 == q66 ~ "Low",
    score >= q66 ~ "High",
    score >= q33 ~ "Medium",
    TRUE ~ "Low"
  )
}

# warnnig big tables
full_table_render_note <- function() {
  tags$div(
    class = "smallNote",
    style = "margin-bottom:8px; color:#666;",
    HTML("Some tables may take a few seconds to appear. Full export of all rows is enabled.")
  )
}

# -----------------------------
# UI
# -----------------------------
ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("🧩 Integrator Inspector")
  ),
  id = "topnav",
  
  tags$script(HTML("
  Shiny.addCustomMessageHandler('filter_candidates_consensus_table_by_cluster', function(message) {
  var cid = (message && message.cluster_id) ? message.cluster_id : '';
  
  if (!window.candidatesConsensusTable || window.candidatesConsensusClusterCol === undefined) {
    return;
  }
  
  var tbl = window.candidatesConsensusTable;
  var col = window.candidatesConsensusClusterCol;
  
  tbl.search('');
  tbl.columns().search('');
  
  if (cid !== '') {
    var esc = $.fn.dataTable.util.escapeRegex(cid);
    tbl.column(col).search('^' + esc + '$', true, false).draw();
  } else {
    tbl.draw();
  }
});
# filtratge taula priority block esquema antic
    Shiny.addCustomMessageHandler('filter_block_table_by_cluster', function(message) {
      var cid = (message && message.cluster_id) ? message.cluster_id : '';
      
      if (!window.blockSummaryTable || window.blockSummaryClusterCol === undefined) {
        console.log('[filter_block_table_by_cluster] table not ready');
        return;
      }
      
      var tbl = window.blockSummaryTable;
      var col = window.blockSummaryClusterCol;
      
      tbl.search('');
      tbl.columns().search('');
      
      if (cid !== '') {
        var esc = $.fn.dataTable.util.escapeRegex(cid);
        tbl.column(col).search('^' + esc + '$', true, false).draw();
      } else {
        tbl.draw();
      }
      
      console.log('[filter_block_table_by_cluster] filtered to', cid);
    });
    
    Shiny.addCustomMessageHandler('filter_gene_table_by_cluster', function(message) {
  var cid = (message && message.cluster_id) ? message.cluster_id : '';
  
  if (!window.prioritizedGenesTable || window.prioritizedGenesClusterCol === undefined) {
    console.log('[filter_gene_table_by_cluster] table not ready');
    return;
  }
  
  var tbl = window.prioritizedGenesTable;
  var col = window.prioritizedGenesClusterCol;
  
  tbl.search('');
  tbl.columns().search('');
  
  if (cid !== '') {
    var esc = $.fn.dataTable.util.escapeRegex(cid);
    var rx = '(^|;\\\\s*)' + esc + '(;|$)';
    tbl.column(col).search(rx, true, false).draw();
  } else {
    tbl.draw();
  }
  
  console.log('[filter_gene_table_by_cluster] filtered to', cid);
});
    Shiny.addCustomMessageHandler('filter_gwasxapp_table_by_cluster', function(message) {
  var cid = (message && message.cluster_id) ? message.cluster_id : '';
  
  if (!window.gwasxappSummaryTable || window.gwasxappSummaryClusterCol === undefined) {
    console.log('[filter_gwasxapp_table_by_cluster] table not ready');
    return;
  }
  
  var tbl = window.gwasxappSummaryTable;
  var col = window.gwasxappSummaryClusterCol;
  
  tbl.search('');
  tbl.columns().search('');
  
  if (cid !== '') {
    var esc = $.fn.dataTable.util.escapeRegex(cid);
    tbl.column(col).search('^' + esc + '$', true, false).draw();
  } else {
    tbl.draw();
  }
  
  console.log('[filter_gwasxapp_table_by_cluster] filtered to', cid);
});
  ")),
  
  tags$style(HTML("
    .panelGrey{ background:#f2f2f2; padding:14px; border-radius:10px; }
    .panelGrey.padDT .dataTables_wrapper{ background:#f2f2f2; padding:10px; border-radius:8px; }
    .panelGrey.padDT table.dataTable{ background:transparent !important; }
    .panelGrey.padDT table.dataTable thead th{ background:transparent !important; }
    .panelGrey.padDT{ min-height:420px; overflow:auto; }
    .smallNote{ color:#555; font-size:13px; line-height:1.35; }

    .panelGrey.padDT .dataTables_wrapper {
      width: 100% !important;
    }

    .panelGrey.padDT table.dataTable {
      width: 100% !important;
    }

    .panelGrey.padDT .dataTables_scrollHeadInner,
    .panelGrey.padDT .dataTables_scrollHeadInner table,
    .panelGrey.padDT .dataTables_scrollBody table {
      width: 100% !important;
    }

    div.dataTables_scrollHead table.dataTable,
    div.dataTables_scrollBody table.dataTable {
      margin-top: 0 !important;
      margin-bottom: 0 !important;
    }
  ")),
  
  tags$head(
    tags$style(HTML("
      .panelGrey.padDT {
        min-height: 420px;
      }

      .dt-outer-scroll {
        overflow-x: auto;
        overflow-y: visible;
      }

      .dt-outer-scroll table.dataTable thead th {
        white-space: nowrap !important;
        vertical-align: middle !important;
      }

      .dt-trunc {
        display: inline-block;
        max-width: 260px;
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        vertical-align: top;
        cursor: pointer;
      }

      .dt-details summary {
        list-style: none;
        cursor: pointer;
      }

      .dt-details summary::-webkit-details-marker {
        display: none;
      }

      .dt-fulltext {
        margin-top: 6px;
        white-space: normal !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
        max-width: 600px;
      }

      table.dataTable tbody tr.selected,
      table.dataTable tbody td.selected {
        background: transparent !important;
      }
    "))
  ),
  
  # ============================================================
  # Nivell 1 · Evidence inputs
  # ============================================================
  tabPanel(
    title = div(
      style = "font-weight:700; font-size:22px; color:#1A4E8A;",
      HTML("🧾 Evidence inputs")
    ),
    fluidPage(
      tabsetPanel(
        
        tabPanel(
          title = "Session & consensus",
          
          # ============================================================
          # Row 1
          # Loaded integrator session | Consensus candidates by cluster
          # ============================================================
          fluidRow(
            column(
              6,
              h3("Loaded integrator session"),
              div(
                class = "panelGrey",
                uiOutput("integrator_session_ui")
              ),
              tags$br(),
              div(
                class = "panelGrey padDT",
                full_table_render_note(),
                withSpinner(DTOutput("manifest_dt")),
                verbatimTextOutput("integrator_session_status")
              )
            ),
            
            column(
              6,
              h3("Consensus candidates by cluster"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                withSpinner(plotly::plotlyOutput("consensus_candidates_apps_plot", height = "420px")),
                uiOutput("consensus_candidates_apps_legend"),
                tags$details(
                  tags$summary("Show numeric table by cluster"),
                  withSpinner(DT::DTOutput("consensus_candidates_apps_table_dt"))
                )
              )
            )
          ),
          
          tags$br(),
          
          # ============================================================
          # Row 2
          # Consensus clusters | Consensus candidates
          # ============================================================
          fluidRow(
            column(
              6,
              h3("Consensus clusters"),
              div(
                class = "panelGrey padDT",
                withSpinner(DTOutput("clusters_consensus_dt"))
              )
            ),
            
            column(
              6,
              h3("Consensus candidates"),
              div(
                class = "panelGrey padDT",
                withSpinner(DTOutput("candidates_consensus_dt"))
              )
            )
          ),
          
          tags$br(),
          
          # ============================================================
          # Row 3
          # Cluster master consistency check | -
          # ============================================================
          fluidRow(
            h3("Cluster master consistency check"),
            uiOutput("clusters_master_consistency_alert"),
            column(
              6,
              
              div(
                class = "panelGrey padDT",
                
                withSpinner(DTOutput("clusters_master_consistency_summary_dt")))
            ),
            
            column(
              6,
             
              div(
                class = "panelGrey padDT",
              withSpinner(DTOutput("clusters_master_consistency_mismatch_dt")))
            )
          )
        ),
        
        tabPanel(
          title = "Bridges",
          
          tabsetPanel(
            id = "bridges_tabs",
            
            tabPanel(
              title = "Gene bridges",
              br(),
              fluidRow(
                column(
                  12,
                  div(class = "panelGrey padDT", full_table_render_note(), DTOutput("all_gene_bridges_dt"))
                )
              )
            ),
            
            tabPanel(
              title = "Term bridges",
              br(),
              fluidRow(
                column(
                  12,
                  div(class = "panelGrey padDT", full_table_render_note(), DTOutput("all_term_bridges_dt"))
                )
              )
            )
          )
        )
      )
    )
  ),
  
  # ============================================================
  # Nivell 1 · LD evidence
  # ============================================================
  tabPanel(
    title = div(
      style = "font-weight:700; font-size:22px; color:#1A4E8A;",
      HTML("🔗 LD evidence")
    ),
    tabsetPanel(
      
      tabPanel(
        title = "Global LD",
        fluidPage(
          fluidRow(
            column(
              6,
              h3("LD/Block settings"),
              div(
                class = "panelGrey",
                fluidRow(
                  column(4, uiOutput("ld_global_pop_ui")),
                  column(
                    4,
                    radioButtons(
                      "ld_global_metric", "LD metric",
                      choices = c("r²" = "R2", "D′" = "Dprime"),
                      selected = "R2",
                      inline = TRUE
                    )
                  ),
                  column(
                    4,
                    numericInput(
                      "ld_global_r2_min",
                      "Min LD threshold",
                      value = 0.6,
                      min = 0,
                      max = 1,
                      step = 0.05
                    )
                  )
                ),
                fluidRow(
                  column(3, actionButton("run_ld_global", "Compute global LD / blocks")),
                  column(9, verbatimTextOutput("ld_global_status"))
                ),
                tags$br(),
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "Global LD / blocks may take several minutes depending on the number of clusters and population subset. ",
                      "This run also annotates hits and genes within the same haploblock as each GWAS lead SNP.",
                      "<br><b>Precomputed LD results, when available for the selected session, are loaded automatically.</b>"
                    )
                  )
                )
              )
            ),
            
            column(
              6,
              h3("Global block overview"),
              div(
                class = "panelGrey",
                tags$div(
                  style = "font-size:12px; color:#666; margin-bottom:6px;",
                  HTML("<i>Tip:</i> Click a block to inspect/filter its corresponding cluster elsewhere.")
                ),
                withSpinner(plotly::plotlyOutput("ld_global_blocks_plot", height = "300px")),
                uiOutput("ld_global_blocks_legend_ui")
              )
            )
          ),
          
          tags$br(),
          
          fluidRow(
            column(
              12,
              h3("LD cluster summary"),
              div(class = "panelGrey padDT", full_table_render_note(), withSpinner(DTOutput("ld_summary_dt")))
            )
          )
        )
      ),
      
      tabPanel(
        title = "Cluster LD",
        ld_integrator_module_ui("ld")
      )
    )
  ),
  
  # ============================================================
  # Nivell 2 · Shared evidence
  # ============================================================
  tabPanel(
    title = div(
      style = "font-weight:700; font-size:22px; color:#1A4E8A;",
      HTML("🧬 Shared evidence")
    ),
    tabsetPanel(
      
      tabPanel(
        title = "Shared genes",
        fluidPage(
          fluidRow(
            column(
              12,
              h3("Shared genes"),
              
              div(
                class = "panelGrey",
                style = "margin-bottom:10px; padding:10px;",
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "<b>Shared genes:</b> ",
                      "Visual summary of genes supported across apps. ",
                      "Use the heatmap to inspect which genes appear in each app, ",
                      "and the bubble plot to identify genes with broader cross-app support and more supporting records."
                    )
                  )
                ),
                fluidRow(
                  column(
                    4,
                    actionButton(
                      "show_shared_genes_heatmap",
                      "Show gene × app heatmap",
                      icon = icon("table-cells-large")
                    )
                  )
                )
              ),
              
              tabsetPanel(
                id = "shared_genes_view_tab",
                
                tabPanel(
                  title = "Table",
                  value = "table",
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("shared_genes_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Bubble plot",
                  value = "bubble",
                  div(
                    class = "panelGrey",
                    style = "margin-top:10px; padding:10px;",
                    
                    tags$div(
                      class = "smallNote",
                      style = "margin-bottom:10px;",
                      HTML(
                        paste0(
                          "Each bubble represents one gene from the shared evidence genes table. ",
                          "The X axis shows the number of supporting apps, ",
                          "bubble size reflects the total number of supporting records, ",
                          "color indicates the supporting app combination, ",
                          "and bubble shape reflects the gene link mode."
                        )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        4,
                        checkboxInput(
                          "shared_genes_bubble_shared_only",
                          "Show only shared genes (n_apps > 1)",
                          value = TRUE
                        )
                      ),
                      column(
                        4,
                        numericInput(
                          "shared_genes_bubble_top_n",
                          "Top N genes",
                          value = 75,
                          min = 10,
                          max = 300,
                          step = 5
                        )
                      )
                    ),
                    
                    uiOutput("shared_genes_bubble_shape_legend"),
                    
                    withSpinner(
                      plotly::plotlyOutput("shared_genes_bubble_plot", height = "800px")
                    )
                  )
                )
              )
            )
          )
        )
      ),
      
      tabPanel(
        title = "Shared terms",
        fluidPage(
          fluidRow(
            column(
              12,
              h3("Shared terms"),
              
              div(
                class = "panelGrey",
                style = "margin-bottom:10px; padding:10px;",
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "<b>Shared terms:</b> ",
                      "Visual summary of normalized terms across apps. ",
                      "Use the heatmap to inspect which terms appear in each app, ",
                      "and the bubble plot to identify the most recurrent and cross-app supported terms."
                    )
                  )
                ),
                fluidRow(
                  column(
                    4,
                    actionButton(
                      "show_traits_heatmap",
                      "Show term × app heatmap",
                      icon = icon("table-cells")
                    )
                  )
                )
              ),
              
              tabsetPanel(
                tabPanel(
                  title = "Table",
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("shared_traits_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Bubble plot",
                  div(
                    class = "panelGrey",
                    style = "margin-top:10px; padding:10px;",
                    
                    tags$div(
                      class = "smallNote",
                      style = "margin-bottom:10px;",
                      HTML(
                        paste0(
                          "Each bubble represents one term from the shared evidence terms table. ",
                          "The X axis shows the number of clusters where the term appears, ",
                          "bubble size reflects the total number of supporting records, ",
                          "and color indicates which app combination supports the term ",
                          "(for example, catalog + nonsyn or catalog + EWAS)."
                        )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        3,
                        selectInput(
                          "shared_terms_bubble_term_type",
                          "Term type",
                          choices = c("All", "disease", "trait"),
                          selected = "All"
                        )
                      ),
                      column(
                        3,
                        checkboxInput(
                          "shared_terms_bubble_shared_only",
                          "Show only shared terms (n_apps > 1)",
                          value = TRUE
                        )
                      ),
                      column(
                        3,
                        numericInput(
                          "shared_terms_bubble_top_n",
                          "Top N terms",
                          value = 75,
                          min = 10,
                          max = 300,
                          step = 5
                        )
                      ),
                      column(
                        3,
                        selectInput(
                          "shared_terms_bubble_label_col",
                          "Label field",
                          choices = c("term_label", "term_norm"),
                          selected = "term_label"
                        )
                      )
                    ),
                    
                    withSpinner(
                      plotly::plotlyOutput("shared_terms_bubble_plot", height = "800px")
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  ),
  
  # ============================================================
  # Nivell 2 · Prioritized evidence
  # ============================================================
  tabPanel(
    title = div(
      style = "font-weight:700; font-size:22px; color:#1A4E8A;",
      HTML("⭐ Prioritized evidence")
    ),
    tabsetPanel(
      
      tabPanel(
        title = "Prioritized clusters",
        fluidPage(
          fluidRow(
            column(
              6,
              uiOutput("ld_notice_clusters"),
              h3("Prioritized clusters"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "<b>Cluster priority score:</b> ",
                      "Each cluster is prioritized using the same canonical logic applied to prioritized genes and LD blocks.",
                      "<br><br>",
                      
                      "The current cluster score combines:<br>",
                      "1) <b>top block support</b>,<br>",
                      "2) <b>other block support (×0.05)</b>,<br>",
                      "3) <b>gene-content bonus</b>, and<br>",
                      "4) <b>fragmentation bonus</b>.",
                      "<br><br>",
                      
                      "Thus, the strongest prioritized block within the cluster has the greatest weight, ",
                      "while additional block-level support is downweighted and structural context contributes through gene content and fragmentation.",
                      "<br><br>",
                      
                      "All clusters remain visible in the plots and tables, including those without prioritized support.",
                      "<br><br>",
                      
                      "<b>Priority classes:</b> ",
                      "Two complementary classifications are shown.<br>",
                      "<b>Absolute priority class</b> uses fixed score thresholds and supports comparison across sessions.<br>",
                      "<b>Relative priority class</b> uses tertiles computed within the current result table and shows whether an item falls in the lower, middle, or upper part of the current score distribution."
                    )
                  )
                )
              ),
              tags$br(),
              actionButton(
                "show_cluster_priority_plot",
                "📊 Cluster priority components",
                class = "btn-primary"
              )
            ),
            
            column(
              6,
              h3("Clusters by support score"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                withSpinner(plotly::plotlyOutput("prioritized_clusters_plot", height = "300px")),
                uiOutput("prioritized_clusters_legend"),
                tags$br()
              )
            )
          ),
          
          tags$br(),
          
          fluidRow(
            column(
              12,
              tags$p(
                style = "font-size:13px; color:#555; margin-bottom:8px; margin-top:10px;",
                "One row per cluster. Cluster scores summarize the canonical prioritized block support together with gene-content and fragmentation bonuses."
              ),
              div(class = "panelGrey padDT", withSpinner(DTOutput("prioritized_cluster_dt")))
            )
          )
        )
      ),
      
      tabPanel(
        title = "Prioritized LD blocks",
        fluidPage(
          fluidRow(
            column(
              6,
              h3("Prioritized LD blocks"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "<b>Block priority score:</b> ",
                      "Each LD block is prioritized from the prioritized GWAS hits assigned to that block.",
                      "<br><br>",
                      
                      "The current block score follows the same canonical logic used for prioritized genes. ",
                      "It is interpreted as the combination of three components:<br>",
                      "1) <b>top GWAS-hit support</b>, corresponding to the highest prioritized GWAS-hit score mapped to the block,<br>",
                      "2) <b>other GWAS-hit support (×0.05)</b>, corresponding to the remaining prioritized GWAS-hit scores in the block with reduced weight, and<br>",
                      "3) <b>gene-content bonus</b>, reflecting the additional biological context provided by the presence and number of genes within the block.",
                      "<br><br>",
                      
                      "This means that the strongest support unit dominates the block score, ",
                      "while additional support contributes more softly and does not inflate the score excessively.",
                      "<br><br>",
                      
                      "All blocks remain visible in the plots and tables, including those without prioritized support.",
                      "<br><br>",
                      
                      "<b>Priority classes:</b> ",
                      "Two complementary classifications are shown.<br>",
                      "<b>Absolute priority class</b> uses fixed score thresholds and supports comparison across sessions.<br>",
                      "<b>Relative priority class</b> uses tertiles computed within the current result table and shows whether an item falls in the lower, middle, or upper part of the current score distribution."
                    )
                  )
                )
              ),
              tags$br(),
              actionButton(
                "show_block_priority_plot",
                "📊 Block priority components",
                class = "btn-primary"
              )
            ),
            
            column(
              6,
              h3("Global LD blocks"),
              h5("Click on block squares to see cluster blocks' score"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                withSpinner(plotly::plotlyOutput("global_ld_blocks_plot", height = "300px")),
                uiOutput("global_ld_blocks_score_legend")
              )
            )
          ),
          
          tags$br(),
          
          fluidRow(
            column(
              12,
              h3("Block priority scores"),
              tags$p(
                style = "font-size:13px; color:#555; margin-bottom:8px; margin-top:10px;",
                "One row per LD block. Block scores summarize the canonical prioritized GWAS-hit support linked to each block together with the block gene-content bonus."
              ),
              div(class = "panelGrey padDT", withSpinner(DTOutput("prioritized_block_dt"))),
              tags$br()
            )
          )
        )
      ),
      
      tabPanel(
        title = "Prioritized genes",
        
        tabsetPanel(
          
          # ============================================================
          # SUBTAB 1 — PRIORITIZATION (EL TEU CODI ORIGINAL)
          # ============================================================
          tabPanel(
            "Prioritization",
            
            fluidPage(
              fluidRow(
                column(
                  6,
                  uiOutput("ld_notice_genes"),
                  h3("Prioritized genes"),
                  div(
                    class = "panelGrey",
                    style = "margin-top:12px;",
                    div(
                      class = "smallNote",
                      HTML(
                        paste0(
                          "<b>Gene priority score:</b> ",
                          "A gene enters the prioritization table if it is linked to at least one prioritized cluster under a canonical gene model.",
                          "<br><br>",
                          "This canonical model combines two complementary gene-linking modes:<br>",
                          "i) <b>position-linked genes</b> (physical), corresponding to genes physically overlapping the cluster interval or structurally linked to a prioritized GWAS hit within the same LD block, and<br>",
                          "ii) <b>effect-linked genes</b> (effect), corresponding to genes supported by real regulatory or epigenetic evidence (for example GTEx or EWAS) connected to a prioritized GWAS hit, even when the gene itself lies outside the cluster interval.",
                          "<br><br>",
                          "Genes included by either route are retained in the table, including genes with score 0 when they belong to the canonical gene set for a prioritized cluster but do not receive any <b>MATCH</b> or <b>MARKER</b> support from prioritized GWAS hits.",
                          "<br><br>",
                          "The <b>gene score</b> is computed from canonical support units linking prioritized GWAS hits to the gene. ",
                          "Each support unit combines the evidence type and the priority of the linked GWAS hit. ",
                          "For scoring purposes, <b>ewasdis</b> and <b>ewastum</b> are treated as a single <b>EWAS</b> evidence family rather than as independent additive supports.",
                          "<br><br>",
                          "To reduce inflation from dense loci, <b>GTEx</b> and <b>EWAS</b> support are collapsed by LD block for both <b>MATCH</b> and <b>MARKER</b> evidence, keeping only the support linked to the highest-priority GWAS hit within each block. ",
                          "When a linked hit falls outside defined LD blocks, it is retained as an independent support unit.",
                          "<br><br>",
                          "For each retained support unit, a unit score is defined from the evidence weight and the linked prioritized GWAS-hit score. ",
                          "The final <b>gene score</b> is then computed as the <b>top support-unit score</b> plus <b>0.05 × the sum of all other support-unit scores</b>, so that the strongest signal drives prioritization while additional independent support contributes more moderately.",
                          "<br><br>",
                          "<b>Priority classes:</b> ",
                          "Two complementary classifications are shown.<br>",
                          "<b>Absolute priority class</b> uses fixed score thresholds and supports comparison across sessions.<br>",
                          "<b>Relative priority class</b> uses tertiles computed within the current result table and shows whether an item falls in the lower, middle, or upper part of the current score distribution."
                        )
                      )
                    )
                  ),
                  tags$br(),
                  tags$br()
                ),
                
                column(
                  6,
                  h3("Genes by cluster"),
                  h5("Click to see genes on cluster"),
                  div(
                    class = "panelGrey",
                    style = "margin-top:12px;",
                    withSpinner(plotly::plotlyOutput("prioritized_genes_cluster_plot", height = "300px")),
                    uiOutput("prioritized_genes_cluster_legend")
                  ),
                  tags$br(),
                  tags$br(),
                  tags$div(
                    style = "display:flex; flex-wrap:wrap; gap:8px; align-items:center;",
                    
                    actionButton(
                      "show_gene_priority_plot",
                      "📊 Gene evidence components",
                      class = "btn-primary"
                    ),
                    
                    actionButton(
                      "show_prioritized_genes_evidence_types_plot",
                      "📊 Gene evidence types",
                      class = "btn-primary"
                    ),
                    
                    actionButton(
                      "show_gene_networks",
                      "🕸️ Gene networks",
                      icon = icon("project-diagram")
                    )
                  )
                )
              ),
              
              tags$br(),
              
              fluidRow(
                h3("Gene priority score"),
                column(
                  12,
                  tags$p(
                    style = "font-size:13px; color:#555; margin-bottom:8px; margin-top:10px;",
                    "One row per gene. Gene scores summarize the prioritized GWAS-hit evidence linked to each gene."
                  ),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("prioritized_gene_dt"))
                  )
                )
              ),
              
              tags$br(),
              
              fluidRow(
                column(
                  12,
                  h4("Gene × GWAS-hit score audit"),
                  tags$p(
                    style = "font-size:13px; color:#555; margin-bottom:8px; margin-top:10px;",
                    "Audit table showing which prioritized GWAS hits contribute to each gene score."
                  ),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("gene_gwas_hit_score_audit_dt"))
                  )
                )
              )
            )
          ),
          
          # ============================================================
          # SUBTAB 2 — ENRICHMENT (NOU)
          # ============================================================
          tabPanel(
            "Enrichment",
            
            fluidPage(
              fluidRow(
                
                column(
                  4,
                  
                  h3("Functional enrichment"),
                  tags$br(),
                  
                  div(
                    class = "panelGrey",
                    style = "padding:12px;",
                    
                    uiOutput("enrich_bg_note"),
                    tags$br(),
                    
                    fluidRow(
                      column(
                        6,
                        radioButtons(
                          "func_scope",
                          "Scope",
                          choices = c("Global" = "global", "Cluster" = "cluster"),
                          selected = "global"
                        )
                      ),
                      column(
                        6,
                        radioButtons(
                          "enrich_background",
                          "Background",
                          choices = c(
                            "Dataset genes" = "dataset",
                            "Default background" = "orgdb"
                          ),
                          selected = "dataset"
                        )
                      )
                    ),
                    
                    conditionalPanel(
                      condition = "input.func_scope == 'cluster'",
                      uiOutput("func_cluster_ui")
                    ),
                    
                    fluidRow(
                      column(
                        6,
                        checkboxGroupInput(
                          "go_ontos",
                          "GO ontologies",
                          choices = c("BP", "CC", "MF"),
                          selected = c("BP", "CC", "MF")
                        )
                      ),
                      column(
                        6,
                        tags$div(
                          style = "margin-top:25px;",
                          checkboxInput("go_simplify", "Simplify GO terms", value = TRUE)
                        )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        6,
                        numericInput(
                          "enrich_pcut",
                          "FDR cutoff",
                          value = 0.05,
                          min = 0,
                          max = 1,
                          step = 0.001
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "enrich_min_gs",
                          "minGSSize",
                          value = 10,
                          min = 1,
                          step = 1
                        )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        6,
                        numericInput(
                          "enrich_max_gs",
                          "maxGSSize",
                          value = 500,
                          min = 1,
                          step = 1
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "go_topn",
                          "Top GO",
                          value = 10,
                          min = 1,
                          step = 1
                        )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        6,
                        numericInput(
                          "enrich_kegg_top",
                          "Top KEGG",
                          value = 15,
                          min = 1,
                          step = 1
                        )
                      ),
                      column(
                        6,
                        tags$div(
                          style = "margin-top:25px;",
                          actionButton(
                            "run_enrich",
                            "Run enrichment",
                            class = "btn-primary",
                            width = "100%"
                          )
                        )
                      )
                    )
                  )
                ),
                
                column(
                  8,
                  
                  tabsetPanel(
                    id = "enrich_tabs",
                    
                    tabPanel(
                      title = "GO",
                      value = "tab_enrich_go",
                      br(),
                      div(
                        class = "panelGrey padDT",
                        withSpinner(plotlyOutput("go_bar", height = "700px"))
                      ),
                      br(),
                      div(
                        class = "panelGrey padDT",
                        withSpinner(DTOutput("go_table"))
                      )
                    ),
                    
                    tabPanel(
                      title = "KEGG",
                      value = "tab_enrich_kegg",
                      br(),
                      div(
                        class = "panelGrey padDT",
                        withSpinner(plotlyOutput("kegg_bar", height = "700px"))
                      ),
                      br(),
                      div(
                        class = "panelGrey padDT",
                        withSpinner(DTOutput("kegg_table"))
                      )
                    )
                  )
                )
              )
            )
          )
          
        )
      ),
      
      tabPanel(
        title = "Prioritized GWAS hits",
        fluidPage(
          gwas_hit_priority_module_ui("gwas_hit_priority")
        )
      ),
      
      tabPanel(
        title = "Priority audit",
        fluidPage(
          
          fluidRow(
            column(
              12,
              h3("Priority audit"),
              div(
                class = "panelGrey",
                style = "margin-top:12px;",
                div(
                  class = "smallNote",
                  HTML(
                    paste0(
                      "<b>Priority audit:</b> ",
                      "This section provides a hierarchical audit of the current canonical prioritization framework across ",
                      "<b>clusters, LD blocks, genes, and GWAS hits</b>.",
                      "<br><br>",
                      
                      "The audit is designed to verify whether prioritization remains coherent across levels: ",
                      "whether strong clusters are supported by strong blocks, whether strong blocks contain strong prioritized genes, ",
                      "and whether prioritized genes are supported by strong GWAS hits. ",
                      "It also helps identify cases where an entity appears unusually strong or weak relative to its linked context.",
                      "<br><br>",
                      
                      "The audit includes three complementary components:<br>",
                      "1) <b>Audit summary</b>, providing a compact cross-level view with consistency flags across cluster, block, gene, and hit levels,<br>",
                      "2) <b>Hierarchy flow</b>, a Sankey visualization summarizing the canonical structure ",
                      "<b>cluster → block → gene → hit</b>, allowing inspection of how prioritized evidence propagates across levels, and<br>",
                      "3) <b>Pairwise coherence plots</b>, showing consistency between adjacent levels:<br>",
                      "&nbsp;&nbsp;&nbsp;• <b>Cluster ↔ Block</b><br>",
                      "&nbsp;&nbsp;&nbsp;• <b>Block ↔ Gene</b><br>",
                      "&nbsp;&nbsp;&nbsp;• <b>Gene ↔ Hit</b>",
                      "<br><br>",
                      
                      "Priority classes (absolute or relative) are used in the plots to facilitate interpretation of coherence patterns across levels.",
                      "<br><br>",
                      
                      "This audit relies on the canonical prioritized tables used throughout the app and is intentionally conservative, ",
                      "focusing on structural consistency rather than redefining scores."
                    )
                  )
                )
              )
            )
          ),
          
          tags$br(),
          
          fluidRow(
            column(
              12,
              tabsetPanel(
                id = "priority_audit_tabs",
                
                tabPanel(
                  title = "Audit summary",
                  tags$br(),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("priority_audit_summary_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Cluster ↔ Block",
                  tags$br(),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("cluster_block_audit_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Block ↔ Gene",
                  tags$br(),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("block_gene_audit_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Gene ↔ Hit",
                  tags$br(),
                  div(
                    class = "panelGrey padDT",
                    withSpinner(DTOutput("gene_hit_audit_dt"))
                  )
                ),
                
                tabPanel(
                  title = "Hierarchy flow",
                  tags$br(),
                  
                  fluidRow(
                    column(
                      8,
                      fluidRow(
                        column(
                          4,
                          selectInput(
                            "priority_audit_sankey_top_n_clusters",
                            "Top clusters to include",
                            choices = c(5, 10, 15, 20, 30, 40, 50, "All"),
                            selected = 15
                          )
                        ),
                        column(
                          4,
                          selectInput(
                            "priority_audit_sankey_top_n_genes_per_block",
                            "Top genes per block",
                            choices = c(2, 3, 4, 5, 8, 10, "All"),
                            selected = 4
                          )
                        ),
                        column(
                          4,
                          selectInput(
                            "priority_audit_sankey_top_n_hits_per_gene",
                            "Top hits per gene",
                            choices = c(2, 3, 4, 5, 8, 10, "All"),
                            selected = 4
                          )
                        )
                      ),
                      
                      tags$div(
                        class = "smallNote",
                        style = "margin-bottom:10px;",
                        HTML(
                          paste0(
                            "<b>Hierarchy flow:</b> ",
                            "This Sankey plot summarizes the canonical hierarchy ",
                            "<b>cluster → block → gene → hit</b>. ",
                            "The selectors control how many of the strongest entities are included in the flow, ",
                            "allowing either a compact overview or a broader exploration of the prioritized hierarchy."
                          )
                        )
                      ),
                      
                      div(
                        class = "panelGrey",
                        withSpinner(plotly::plotlyOutput("priority_audit_sankey_full", height = "1200px"))
                      )
                    ),
                    
                    column(
                      4,
                      radioButtons(
                        "priority_audit_class_mode",
                        "Priority class mode",
                        choices = c("ABS", "REL"),
                        selected = "ABS",
                        inline = TRUE
                      ),
                      
                      div(
                        class = "panelGrey",
                        style = "margin-bottom:12px;",
                        tags$div(
                          style = "font-weight:600; margin-bottom:6px;",
                          "Cluster ↔ Block"
                        ),
                        withSpinner(plotly::plotlyOutput("cluster_block_audit_class_plot", height = "360px"))
                      ),
                      
                      div(
                        class = "panelGrey",
                        style = "margin-bottom:12px;",
                        tags$div(
                          style = "font-weight:600; margin-bottom:6px;",
                          "Block ↔ Gene"
                        ),
                        withSpinner(plotly::plotlyOutput("block_gene_audit_class_plot", height = "360px"))
                      ),
                      
                      div(
                        class = "panelGrey",
                        tags$div(
                          style = "font-weight:600; margin-bottom:6px;",
                          "Gene ↔ Hit"
                        ),
                        withSpinner(plotly::plotlyOutput("gene_hit_audit_class_plot", height = "360px"))
                      )
                    )
                  )
                ),
                
                tabPanel(
                  title = "Audit table (temp)",
                  tags$br(),
                  fluidPage(
                    priority_audit_hierarchy_module_ui("audit_table_temp")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
  
)
# -----------------------------
# 
# -----------------------------
server <- function(input, output, session) {
  
  # exportacio taules model
  dt_buttons_full_allcols <- list(
    list(
      extend = "copyHtml5",
      text = "Copy",
      exportOptions = list(
        columns = ":all",
        modifier = list(
          page = "all",
          search = "none",
          order = "index"
        )
      )
    ),
    list(
      extend = "csvHtml5",
      text = "CSV",
      exportOptions = list(
        columns = ":all",
        modifier = list(
          page = "all",
          search = "none",
          order = "index"
        )
      )
    ),
    list(
      extend = "excelHtml5",
      text = "Excel",
      exportOptions = list(
        columns = ":all",
        modifier = list(
          page = "all",
          search = "none",
          order = "index"
        )
      )
    )
  )
  ##############
  ld_summary_version <- reactiveVal(0)
  ld_details_version <- reactiveVal(0)
  
  #---------------- Reactius base --------------------
  session_dir_r <- reactive({
    ss <- available_sessions()
    req(length(ss) > 0)
    
    sel <- input$integrator_session
    if (is.null(sel) || !nzchar(sel) || !sel %in% names(ss)) {
      sel <- names(ss)[1]
    }
    
    unname(ss[[sel]])
  })
  
  selected_session_dir <- session_dir_r
  
  ld_summary_r <- reactive({
    ld_summary_version()
    
    path <- file.path(session_dir_r(), "ld_cluster_summary.rds")
    if (!file.exists(path)) return(tibble::tibble())
    
    x <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.data.frame(x)) return(tibble::tibble())
    
    tibble::as_tibble(x)
  })
  
  ld_details_r <- reactive({
    ld_details_version()
    
    path <- file.path(session_dir_r(), "ld_cluster_details.rds")
    if (!file.exists(path)) return(list())
    
    x <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.list(x)) return(list())
    
    x
  })
  
  gwas_bridge_r <- reactive({
    read_gwas_significance_bridge(session_dir_r())
  })
  #---------------------------------------------------
  
  available_sessions <- reactive({
    list_integrator_sessions(integrator_root)
  })
  
  output$integrator_session_ui <- renderUI({
    ss <- available_sessions()
    
    if (!length(ss)) {
      return(
        div(
          style = "color:#b00020; font-weight:600;",
          paste0("No session folders with manifest.rds found in: ", integrator_root)
        )
      )
    }
    
    sel <- input$integrator_session
    if (is.null(sel) || !nzchar(sel) || !sel %in% names(ss)) {
      sel <- names(ss)[1]
    }
    
    selectInput(
      "integrator_session",
      "Integrator session folder",
      choices = names(ss),
      selected = sel
    )
  })
  
  output$integrator_session_status <- renderText({
    ss_dir <- selected_session_dir()
    mf <- selected_manifest()
    chk <- validate_manifest_basic(mf)
    
    paste0(
      "session_dir: ", ss_dir, "\n",
      "manifest_check: ", chk
    )
  })
  
  # ============================================================
  # =========     EVIDENCE INPUTS TAB       ====================
  # ============================================================

  # ============================================================
  # TABLE CACHE PER SESSION
  # Carrega i prepara una sola vegada quan canvia session_dir_r()
  # ============================================================
  
  tables_cache <- reactiveValues(
    manifest_obj = NULL,
    manifest_df = NULL,
    clusters_raw = NULL,
    clusters_dt = NULL,
    candidates_raw = NULL,
    candidates_dt = NULL,
    candidates_cluster_col_idx = NA_integer_
  )
  
  observeEvent(session_dir_r(), {
    req(session_dir_r())
    
    t0 <- Sys.time()
    message("[TIMING] table preload start")
    
    # ---------------- manifest
    t_manifest <- Sys.time()
    mf <- read_manifest_safe(session_dir_r())
    tables_cache$manifest_obj <- mf
    tables_cache$manifest_df <- manifest_to_display_df(mf)
    message("[TIMING] manifest: ", round(as.numeric(Sys.time() - t_manifest, units = "secs"), 2), " s")
    
    # ---------------- clusters read
    t_cl_read <- Sys.time()
    cl <- integrated_clusters_from_rds(session_dir_r())
    message("[TIMING] clusters read: ", round(as.numeric(Sys.time() - t_cl_read, units = "secs"), 2), " s")
    
    # ---------------- clusters sanitize
    t_cl_san <- Sys.time()
    if (is.data.frame(cl) && nrow(cl) > 0) {
      tables_cache$clusters_raw <- cl
      tables_cache$clusters_dt <- sanitize_dt_types(cl)
    } else {
      tables_cache$clusters_raw <- tibble::tibble()
      tables_cache$clusters_dt <- NULL
    }
    message("[TIMING] clusters sanitize: ", round(as.numeric(Sys.time() - t_cl_san, units = "secs"), 2), " s")
    
    # ---------------- candidates read
    t_ca_read <- Sys.time()
    ca <- integrated_candidates_from_rds(session_dir_r())
    message("[TIMING] candidates read: ", round(as.numeric(Sys.time() - t_ca_read, units = "secs"), 2), " s")
    
    # ---------------- candidates sanitize
    t_ca_san <- Sys.time()
    if (is.data.frame(ca) && nrow(ca) > 0) {
      ca2 <- sanitize_dt_types(ca)
      
      tables_cache$candidates_raw <- ca
      tables_cache$candidates_dt <- ca2
      
      idx <- which(names(ca2) == "cluster_id")
      tables_cache$candidates_cluster_col_idx <- if (length(idx)) idx else NA_integer_
    } else {
      tables_cache$candidates_raw <- tibble::tibble()
      tables_cache$candidates_dt <- NULL
      tables_cache$candidates_cluster_col_idx <- NA_integer_
    }
    message("[TIMING] candidates sanitize: ", round(as.numeric(Sys.time() - t_ca_san, units = "secs"), 2), " s")
    
    message("[TIMING] table preload total: ", round(as.numeric(Sys.time() - t0, units = "secs"), 2), " s")
  }, ignoreInit = FALSE)
  
  
  # manifest #########################
  
  manifest_r <- reactive({
    tables_cache$manifest_obj
  })
  
  selected_manifest <- manifest_r
  
  output$manifest_dt <- DT::renderDT({
    dt <- tables_cache$manifest_df
    req(is.data.frame(dt))
    
    DT::datatable(
      dt,
      rownames = FALSE,
      selection = "none",
      options = list(
        dom = "t",
        pageLength = nrow(dt),
        autoWidth = TRUE
      )
    )
  }, server = FALSE)
  
  # clusters #########################
  
  integrated_clusters_r <- reactive({
    tables_cache$clusters_raw
  })
  
  clusters_consensus <- integrated_clusters_r
  
  output$clusters_consensus_dt <- DT::renderDT({
    dt <- tables_cache$clusters_dt
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No consensus clusters available."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    DT::datatable(
      dt,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copyHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          ),
          list(
            extend = "csvHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          ),
          list(
            extend = "excelHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        ),
        pageLength = 15,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  ### validacio Cluster consensus
  
  clusters_master_consistency_df <- reactive({
    sess_dir <- selected_session_dir()
    
    validate(
      need(is.character(sess_dir) && length(sess_dir) == 1 && nzchar(sess_dir) && dir.exists(sess_dir),
           "No valid Integrator session directory available.")
    )
    
    apps <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
    
    read_cluster_master_safe <- function(app) {
      path <- file.path(sess_dir, paste0(app, "_clusters_master.rds"))
      
      if (!file.exists(path)) {
        return(list(
          app = app,
          path = path,
          exists = FALSE,
          data = tibble::tibble()
        ))
      }
      
      obj <- tryCatch(readRDS(path), error = function(e) NULL)
      
      if (!is.data.frame(obj) || !nrow(obj)) {
        return(list(
          app = app,
          path = path,
          exists = TRUE,
          data = tibble::tibble()
        ))
      }
      
      needed <- c("cluster_id", "chr", "start", "end")
      if (!all(needed %in% names(obj))) {
        return(list(
          app = app,
          path = path,
          exists = TRUE,
          data = tibble::tibble()
        ))
      }
      
      dat <- obj %>%
        dplyr::transmute(
          app = app,
          cluster_id = trimws(as.character(cluster_id)),
          chr = as.character(chr),
          start = suppressWarnings(as.numeric(start)),
          end = suppressWarnings(as.numeric(end))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(chr), nzchar(chr),
          is.finite(start), is.finite(end)
        ) %>%
        dplyr::mutate(
          start = pmin(start, end),
          end = pmax(start, end)
        )
      
      list(
        app = app,
        path = path,
        exists = TRUE,
        data = dat
      )
    }
    
    lst <- lapply(apps, read_cluster_master_safe)
    
    summary_df <- dplyr::bind_rows(lapply(lst, function(x) {
      dat <- x$data
      tibble::tibble(
        app = x$app,
        file_exists = x$exists,
        n_rows = if (is.data.frame(dat)) nrow(dat) else 0L,
        n_cluster_ids = if (is.data.frame(dat) && nrow(dat)) dplyr::n_distinct(dat$cluster_id) else 0L,
        n_duplicated_cluster_ids = if (is.data.frame(dat) && nrow(dat)) {
          sum(duplicated(dat$cluster_id))
        } else {
          0L
        },
        file_path = x$path
      )
    }))
    
    all_df <- dplyr::bind_rows(lapply(lst, `[[`, "data"))
    
    mismatch_df <- if (nrow(all_df)) {
      all_df %>%
        dplyr::mutate(
          coord = paste0("chr", chr, ":", start, "-", end)
        ) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          n_apps = dplyr::n_distinct(app),
          apps = paste(sort(unique(app)), collapse = "; "),
          n_coords = dplyr::n_distinct(coord),
          coords = paste(sort(unique(coord)), collapse = " | "),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          issue_type = dplyr::case_when(
            n_coords > 1 ~ "coordinate_mismatch",
            TRUE ~ "ok"
          )
        ) %>%
        dplyr::filter(issue_type != "ok")
    } else {
      tibble::tibble()
    }
    
    list(
      summary = summary_df,
      mismatches = mismatch_df,
      all = all_df
    )
  })
  
  output$clusters_master_consistency_alert <- renderUI({
    chk <- clusters_master_consistency_df()
    sm <- chk$summary
    mm <- chk$mismatches
    
    if (!is.data.frame(sm) || !nrow(sm)) return(NULL)
    
    n_missing_files <- sum(!sm$file_exists, na.rm = TRUE)
    n_dup <- sum(sm$n_duplicated_cluster_ids, na.rm = TRUE)
    n_mismatch <- if (is.data.frame(mm)) nrow(mm) else 0L
    
    if (n_missing_files == 0 && n_dup == 0 && n_mismatch == 0) {
      return(
        tags$div(
          style = "padding:10px 14px; border-radius:8px; background:#d9ead3; color:#1b5e20; margin-bottom:12px;",
          tags$b("Cluster master consistency check: OK. "),
          "All available _clusters_master.rds files are present and structurally consistent."
        )
      )
    }
    
    tags$div(
      style = "padding:10px 14px; border-radius:8px; background:#fff3cd; color:#7a5a00; margin-bottom:12px;",
      tags$b("Cluster master consistency check: WARNING. "),
      paste0(
        "Missing files: ", n_missing_files,
        " | duplicated cluster IDs: ", n_dup,
        " | cluster mismatches: ", n_mismatch
      )
    )
  })
  
  output$clusters_master_consistency_summary_dt <- DT::renderDT({
    chk <- clusters_master_consistency_df()
    df <- chk$summary
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster master summary available.")
    )
    
    show_df <- df %>%
      dplyr::select(-dplyr::any_of("n_rows")) %>%
      dplyr::rename(
        n_clusters_with_hits = n_cluster_ids
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5"),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  output$clusters_master_consistency_mismatch_dt <- DT::renderDT({
    chk <- clusters_master_consistency_df()
    df <- chk$mismatches
    
    if (!is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(
        data.frame(Message = "No cluster master mismatches detected."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    DT::datatable(
      df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5"),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  # candidates #########################
  
  integrated_candidates_r <- reactive({
    tables_cache$candidates_raw
  })
  
  candidates_consensus <- integrated_candidates_r
  
  output$candidates_consensus_dt <- DT::renderDT({
    dt <- tables_cache$candidates_dt
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No consensus candidates available."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    DT::datatable(
      dt,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copyHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          ),
          list(
            extend = "csvHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          ),
          list(
            extend = "excelHtml5",
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        ),
        pageLength = 15,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  #------------------------------------------------------------------  
  # ============================================================================
  # ============================ SHARED EVIDENCE TAB  ==========================
  # ============================================================================

  # ==============================
  # SHARED GENES
  # ==============================
  shared_genes_df <- reactive({
    
    ge <- gene_bridges_combined()
    if (!is.data.frame(ge) || !nrow(ge)) return(tibble::tibble())
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(sort(unique(x)), collapse = "; ")
    }
    
    ge %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        source_app = as.character(source_app),
        evidence_type = as.character(evidence_type),
        cluster_id = as.character(cluster_id),
        chr = as.character(chr)
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        n_apps = dplyr::n_distinct(source_app[!is.na(source_app) & nzchar(source_app)]),
        apps = collapse_unique_semicolon(source_app),
        n_evidence_types = dplyr::n_distinct(evidence_type[!is.na(evidence_type) & nzchar(evidence_type)]),
        evidence_types = collapse_unique_semicolon(evidence_type),
        n_clusters = dplyr::n_distinct(cluster_id[!is.na(cluster_id) & nzchar(cluster_id)]),
        clusters = collapse_unique_semicolon(cluster_id),
        n_records = dplyr::n(),
        chr_set = collapse_unique_semicolon(chr),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        has_position_link = stringr::str_detect(evidence_types, "catalog_gene|nonsyn_gene"),
        has_effect_link = stringr::str_detect(evidence_types, "eqtl_gene|ewas_nearby_gene"),
        gene_link_mode = dplyr::case_when(
          has_position_link & has_effect_link ~ "position + effect",
          has_position_link ~ "position",
          has_effect_link ~ "effect",
          TRUE ~ "other"
        ),
        gene_label = dplyr::case_when(
          gene_link_mode == "position + effect" ~ paste0("◐ ", gene),
          gene_link_mode == "position" ~ paste0("● ", gene),
          gene_link_mode == "effect" ~ paste0("◆ ", gene),
          TRUE ~ gene
        ),
        gene_link_rank = dplyr::case_when(
          gene_link_mode == "position + effect" ~ 1L,
          gene_link_mode == "position" ~ 2L,
          gene_link_mode == "effect" ~ 3L,
          TRUE ~ 4L
        )
      ) %>%
      dplyr::arrange(
        gene_link_rank,
        dplyr::desc(n_apps),
        dplyr::desc(n_evidence_types),
        dplyr::desc(n_clusters),
        gene
      ) %>%
      sanitize_dt_types()
  })
  
  output$genes_dt <- DT::renderDT({
    dt <- genes_integrated()
    if (!nrow(dt)) {
      return(DT::datatable(data.frame(Message = "No optional gene evidence loaded yet."), rownames = FALSE))
    }
    DT::datatable(
      sanitize_dt_types(dt),
      rownames = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "450px",
        autoWidth = FALSE
      ),
      class = "compact stripe hover order-column",
      fillContainer = TRUE
    )
  })
  
  shared_genes_dt_data <- reactive({
    req(input$shared_genes_view_tab == "table")
    
    dt <- shared_genes_df()
    
    validate(
      need(is.data.frame(dt), "Load optional gene evidence to populate this panel."),
      need(nrow(dt) > 0, "Load optional gene evidence to populate this panel.")
    )
    
    dt %>%
      dplyr::mutate(
        gene = make_genecards_links(gene)
      )
  }) %>%
    bindCache(shared_genes_df())
  
  output$shared_genes_dt <- DT::renderDT({
    show_dt <- shared_genes_dt_data()
    
    DT::datatable(
      show_dt,
      rownames = FALSE,
      escape = FALSE,
      selection = "none",
      extensions = "Buttons",
      width = "100%",
      class = "compact stripe hover order-column",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          ),
          list(
            extend = "csv",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          ),
          list(
            extend = "excel",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          )
        ),
        pageLength = 15,
        autoWidth = FALSE,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  outputOptions(output, "shared_genes_dt", suspendWhenHidden = TRUE)
  
  shared_genes_bubble_df <- reactive({
    df <- shared_genes_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No shared genes available.")
    )
    
    needed_cols <- c("gene", "n_apps", "n_records", "apps")
    
    validate(
      need(all(needed_cols %in% names(df)), "Shared genes table missing required columns.")
    )
    
    normalize_apps_signature <- function(x) {
      x <- as.character(x)
      
      vapply(x, function(s) {
        if (is.na(s) || !nzchar(trimws(s))) return("")
        
        parts <- unlist(strsplit(s, "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
        parts <- trimws(tolower(parts))
        parts <- parts[nzchar(parts)]
        parts <- unique(parts)
        parts <- sort(parts)
        
        paste(parts, collapse = " + ")
      }, character(1))
    }
    
    normalize_gene_link_mode <- function(x) {
      x <- as.character(x)
      x <- trimws(tolower(x))
      
      dplyr::case_when(
        x %in% c("position", "pos", "position-linked", "position linked") ~ "position",
        x %in% c("effect", "eff", "effect-linked", "effect linked") ~ "effect",
        x %in% c("position+effect", "position + effect", "pos+effect", "position/effect",
                 "position + effect-linked", "position+effect-linked",
                 "position + effect linked", "position effect") ~ "position+effect",
        TRUE ~ "other"
      )
    }
    
    out <- df %>%
      dplyr::mutate(
        gene = as.character(gene),
        apps = dplyr::coalesce(as.character(apps), ""),
        evidence_types = if ("evidence_types" %in% names(.)) dplyr::coalesce(as.character(evidence_types), "") else "",
        gene_link_mode = if ("gene_link_mode" %in% names(.)) dplyr::coalesce(as.character(gene_link_mode), "") else "",
        n_apps = suppressWarnings(as.numeric(n_apps)),
        n_records = suppressWarnings(as.numeric(n_records))
      ) %>%
      dplyr::filter(
        !is.na(n_apps),
        !is.na(n_records),
        nzchar(trimws(gene))
      )
    
    if (isTRUE(input$shared_genes_bubble_shared_only)) {
      out <- out %>%
        dplyr::filter(n_apps > 1)
    }
    
    out <- out %>%
      dplyr::mutate(
        plot_label = gene,
        apps_signature = normalize_apps_signature(apps),
        gene_link_mode_norm = normalize_gene_link_mode(gene_link_mode)
      ) %>%
      dplyr::arrange(
        dplyr::desc(n_apps),
        dplyr::desc(n_records),
        plot_label
      )
    
    top_n <- suppressWarnings(as.numeric(input$shared_genes_bubble_top_n))
    if (is.na(top_n) || top_n < 1) top_n <- 75
    
    n_levels <- length(unique(out$apps_signature))
    n_per_level <- max(1, ceiling(top_n / max(1, n_levels)))
    
    out <- out %>%
      dplyr::group_by(apps_signature) %>%
      dplyr::slice_head(n = n_per_level) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        dplyr::desc(n_apps),
        dplyr::desc(n_records),
        plot_label
      ) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::mutate(
        hover_txt = paste0(
          "<b>", plot_label, "</b>",
          "<br>Supporting apps: ", apps_signature,
          "<br>n_apps: ", n_apps,
          ifelse(nzchar(gene_link_mode), paste0("<br>Gene link mode: ", gene_link_mode), ""),
          ifelse(nzchar(evidence_types), paste0("<br>Evidence types: ", evidence_types), ""),
          "<br>Total records: ", n_records
        )
      )
    
    out$plot_label <- factor(out$plot_label, levels = rev(unique(out$plot_label)))
    
    out
  })
  
  output$shared_genes_bubble_shape_legend <- renderUI({
    tags$div(
      class = "smallNote",
      style = "margin:8px 0 12px 0; color:#444;",
      
      tags$span(
        style = "font-weight:600; margin-right:10px;",
        "Gene link mode"
      ),
      
      tags$span(
        style = "display:inline-block; margin-right:14px;",
        HTML("&#9679; position")
      ),
      
      tags$span(
        style = "display:inline-block; margin-right:14px;",
        HTML("&#9670; effect")
      ),
      
      tags$span(
        style = "display:inline-block; margin-right:14px;",
        HTML("&#9632; position+effect")
      ),
      
      tags$span(
        style = "display:inline-block; margin-right:14px;",
        HTML("&#9675; other/unknown")
      )
    )
  })
  
  output$shared_genes_bubble_plot <- plotly::renderPlotly({
    req(input$shared_genes_view_tab == "bubble")
    req(input$shared_genes_bubble_top_n)
    req(!is.null(input$shared_genes_bubble_shared_only))
    
    df <- shared_genes_bubble_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No shared genes available for plotting.")
    )
    
    rescale_size <- function(x, to = c(10, 45)) {
      x <- as.numeric(x)
      rng <- range(x, na.rm = TRUE)
      
      if (!all(is.finite(rng)) || diff(rng) == 0) {
        return(rep(mean(to), length(x)))
      }
      
      to[1] + (x - rng[1]) * (to[2] - to[1]) / diff(rng)
    }
    
    symbol_map <- c(
      "position" = "circle",
      "effect" = "diamond",
      "position+effect" = "square",
      "other" = "circle-open"
    )
    
    app_sig_levels <- unique(df$apps_signature)
    mode_levels <- c("position", "effect", "position+effect", "other")
    mode_levels <- mode_levels[mode_levels %in% unique(df$gene_link_mode_norm)]
    
    pal <- viridisLite::viridis(length(app_sig_levels))
    names(pal) <- app_sig_levels
    
    p <- plotly::plot_ly(
      source = "shared_genes_bubble_plot"
    )
    
    for (sig in app_sig_levels) {
      first_for_sig <- TRUE
      
      for (mode_i in mode_levels) {
        subdf <- df %>%
          dplyr::filter(apps_signature == sig, gene_link_mode_norm == mode_i)
        
        if (!nrow(subdf)) next
        
        p <- p %>%
          plotly::add_markers(
            data = subdf,
            x = ~n_apps,
            y = ~plot_label,
            text = ~hover_txt,
            hoverinfo = "text",
            name = sig,
            legendgroup = sig,
            showlegend = first_for_sig,
            marker = list(
              size = rescale_size(subdf$n_records),
              sizemode = "diameter",
              color = pal[[sig]],
              symbol = symbol_map[[mode_i]],
              opacity = 0.85,
              line = list(
                color = "rgba(80,80,80,0.4)",
                width = 0.7
              )
            )
          )
        
        first_for_sig <- FALSE
      }
    }
    
    p %>%
      plotly::layout(
        xaxis = list(
          title = "Number of supporting apps",
          dtick = 1,
          rangemode = "tozero",
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = levels(df$plot_label)
        ),
        legend = list(
          title = list(text = "Supporting app combination")
        ),
        margin = list(l = 220, r = 30, t = 30, b = 60),
        hoverlabel = list(align = "left")
      )
  })
  
  # ==============================
  # SHARED TRAITS
  # ==============================
  shared_traits_df <- reactive({
    tr <- term_bridges_combined()
    
    if (!is.data.frame(tr) || !nrow(tr)) {
      return(tibble::tibble())
    }
    
    tr2 <- tr %>%
      dplyr::filter(!is.na(term), nzchar(term)) %>%
      dplyr::mutate(
        term = trimws(as.character(term)),
        term_type = dplyr::coalesce(as.character(term_type), "unknown"),
        source_app = trimws(as.character(source_app)),
        evidence_type = trimws(as.character(evidence_type)),
        cluster_id = trimws(as.character(cluster_id)),
        chr = suppressWarnings(as.integer(chr)),
        term_norm = normalize_trait_term(term)
      ) %>%
      dplyr::filter(!is.na(term_norm), nzchar(term_norm))
    
    if (!nrow(tr2)) {
      return(tibble::tibble())
    }
    
    tr2 %>%
      dplyr::group_by(term_norm, term_type) %>%
      dplyr::summarise(
        term_label = safe_top1(term[order(nchar(term), term)]),
        n_apps = dplyr::n_distinct(source_app[!is.na(source_app) & nzchar(source_app)]),
        apps = paste(sort(unique(stats::na.omit(source_app[source_app != ""]))), collapse = "; "),
        n_evidence_types = dplyr::n_distinct(evidence_type[!is.na(evidence_type) & nzchar(evidence_type)]),
        evidence_types = paste(sort(unique(stats::na.omit(evidence_type[evidence_type != ""]))), collapse = "; "),
        n_clusters = dplyr::n_distinct(cluster_id[!is.na(cluster_id) & nzchar(cluster_id)]),
        clusters = paste(sort(unique(stats::na.omit(cluster_id[cluster_id != ""]))), collapse = "; "),
        n_records = dplyr::n(),
        chr_set = paste(sort(unique(stats::na.omit(chr))), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        term_label = dplyr::coalesce(term_label, term_norm)
      ) %>%
      dplyr::arrange(
        dplyr::desc(n_apps),
        dplyr::desc(n_evidence_types),
        dplyr::desc(n_clusters),
        term_type,
        term_label
      ) %>%
      dplyr::select(
        term_label, term_norm, term_type, n_apps, apps,
        n_evidence_types, evidence_types, n_clusters, clusters,
        chr_set, n_records
      ) %>%
      sanitize_dt_types()
  })
  
   output$traits_dt <- DT::renderDT({
    dt <- traits_integrated()
    if (!nrow(dt)) {
      return(DT::datatable(data.frame(Message = "No optional trait / disease evidence loaded yet."), rownames = FALSE))
    }
    DT::datatable(
      sanitize_dt_types(dt),
      rownames = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "450px",
        autoWidth = FALSE
      ),
      class = "compact stripe hover order-column",
      fillContainer = TRUE
    )
  }) 
  
   output$shared_traits_dt <- DT::renderDT({
    dt <- shared_traits_df()
    
    if (!nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "Load optional trait/disease evidence to populate this panel."),
        rownames = FALSE
      ))
    }
    
    dt_export <- dt
    dt_show <- dt
    
    long_text_cols <- intersect(
      c("term_norm", "clusters"),
      names(dt_show)
    )
    
    if ("term_label" %in% names(dt_show)) {
      dt_show$term_label_export <- dt_export$term_label
    }
    
    if ("term_norm" %in% names(dt_show)) {
      dt_show$term_norm_export <- dt_export$term_norm
    }
    
    if ("clusters" %in% names(dt_show)) {
      dt_show$clusters_export <- dt_export$clusters
    }
    
    if ("block_genes" %in% names(dt_show)) {
      dt_show$block_genes_export <- dt_export$block_genes
      
      dt_show$block_genes <- vapply(
        dt_export$block_genes,
        function(x) as.character(make_expandable_genecards_cell(x, max_show = 5L)),
        character(1)
      )
    }
    
    for (cc in long_text_cols) {
      dt_show[[cc]] <- truncate_for_dt(dt_show[[cc]], n = 90)
    }
    
    if ("term_label" %in% names(dt_show)) {
      dt_show <- dt_show %>%
        dplyr::mutate(
          term_label = make_gwascatalog_term_links(term_label)
        )
    }
    
    hidden_export_cols <- intersect(
      c("term_label_export", "term_norm_export", "clusters_export", "block_genes_export"),
      names(dt_show)
    )
    
    js_export_body <- DT::JS(
      "function(data, row, column, node) {
       var tbl = window.sharedTraitsTable;
       if (!tbl) return $('<div>').html(data).text();
       
       var rowData = tbl.row(row).data();
       var colName = tbl.column(column).header().textContent.trim();
       
       if (!rowData) return $('<div>').html(data).text();
       
       if (colName === 'term_label' && rowData.term_label_export !== undefined) {
         return String(rowData.term_label_export);
       }
       if (colName === 'term_norm' && rowData.term_norm_export !== undefined) {
         return String(rowData.term_norm_export);
       }
       if (colName === 'clusters' && rowData.clusters_export !== undefined) {
         return String(rowData.clusters_export);
       }
       if (colName === 'block_genes' && rowData.block_genes_export !== undefined) {
         return String(rowData.block_genes_export);
       }
       
       return $('<div>').html(data).text();
     }"
    )
    
    DT::datatable(
      dt_show,
      rownames = FALSE,
      escape = FALSE,
      selection = "none",
      extensions = "Buttons",
      width = "100%",
      class = "compact stripe hover order-column",
      callback = DT::JS(
        "var tbl = table;
       window.sharedTraitsTable = tbl;
       tbl.columns.adjust();
       $(window).on('resize.sharedTraitsTable', function() {
         tbl.columns.adjust();
       });"
      ),
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(
              columns = which(!names(dt_show) %in% hidden_export_cols) - 1L,
              format = list(body = js_export_body)
            )
          ),
          list(
            extend = "csv",
            exportOptions = list(
              columns = which(!names(dt_show) %in% hidden_export_cols) - 1L,
              format = list(body = js_export_body)
            )
          ),
          list(
            extend = "excel",
            exportOptions = list(
              columns = which(!names(dt_show) %in% hidden_export_cols) - 1L,
              format = list(body = js_export_body)
            )
          )
        ),
        pageLength = 15,
        autoWidth = TRUE,
        columnDefs = if (length(hidden_export_cols)) {
          lapply(hidden_export_cols, function(cc) {
            list(
              targets = which(names(dt_show) == cc) - 1L,
              visible = FALSE,
              searchable = FALSE
            )
          })
        } else {
          list()
        }
      )
    )
  }, server = FALSE)
  
  # ==============================
  # Shared evidence terms bubble plot
  # ==============================
  
  shared_terms_bubble_df <- reactive({
    df <- shared_traits_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No shared evidence terms available.")
    )
    
    needed_cols <- c(
      "term_label", "term_norm", "term_type", "n_apps",
      "n_clusters", "n_records", "apps", "evidence_types", "clusters", "chr_set"
    )
    
    validate(
      need(all(needed_cols %in% names(df)), "Shared terms table does not contain the required columns.")
    )
    
    normalize_apps_signature <- function(x) {
      x <- as.character(x)
      
      vapply(x, function(s) {
        if (is.na(s) || !nzchar(trimws(s))) return("")
        
        parts <- unlist(strsplit(s, "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
        parts <- trimws(tolower(parts))
        parts <- parts[nzchar(parts)]
        parts <- unique(parts)
        parts <- sort(parts)
        
        paste(parts, collapse = " + ")
      }, character(1))
    }
    
    out <- df %>%
      dplyr::mutate(
        term_label = as.character(term_label),
        term_norm = as.character(term_norm),
        term_type = as.character(term_type),
        apps = dplyr::coalesce(as.character(apps), ""),
        evidence_types = dplyr::coalesce(as.character(evidence_types), ""),
        clusters = dplyr::coalesce(as.character(clusters), ""),
        chr_set = dplyr::coalesce(as.character(chr_set), ""),
        n_apps = suppressWarnings(as.numeric(n_apps)),
        n_clusters = suppressWarnings(as.numeric(n_clusters)),
        n_records = suppressWarnings(as.numeric(n_records))
      ) %>%
      dplyr::filter(
        !is.na(n_apps),
        !is.na(n_clusters),
        !is.na(n_records)
      )
    
    if (isTRUE(input$shared_terms_bubble_shared_only)) {
      out <- out %>%
        dplyr::filter(n_apps > 1)
    }
    
    if (!is.null(input$shared_terms_bubble_term_type) &&
        input$shared_terms_bubble_term_type != "All") {
      out <- out %>%
        dplyr::filter(term_type == input$shared_terms_bubble_term_type)
    }
    
    label_col <- input$shared_terms_bubble_label_col %||% "term_label"
    if (!label_col %in% c("term_label", "term_norm")) {
      label_col <- "term_label"
    }
    
    if (label_col == "term_norm") {
      out <- out %>%
        dplyr::mutate(plot_label = term_norm)
    } else {
      out <- out %>%
        dplyr::mutate(plot_label = term_label)
    }
    
    out <- out %>%
      dplyr::mutate(
        apps_signature = normalize_apps_signature(apps)
      ) %>%
      dplyr::arrange(
        dplyr::desc(n_clusters),
        dplyr::desc(n_records),
        apps_signature,
        plot_label
      )
    
    top_n <- suppressWarnings(as.numeric(input$shared_terms_bubble_top_n))
    if (is.na(top_n) || top_n < 1) top_n <- 75
    
    n_levels <- length(unique(out$apps_signature))
    n_per_level <- max(1, ceiling(top_n / max(1, n_levels)))
    
    out <- out %>%
      dplyr::group_by(apps_signature) %>%
      dplyr::slice_head(n = n_per_level) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        dplyr::desc(n_clusters),
        dplyr::desc(n_records),
        apps_signature,
        plot_label
      ) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::mutate(
        hover_txt = paste0(
          "<b>", plot_label, "</b>",
          "<br>Type: ", term_type,
          "<br>Supporting apps: ", apps_signature,
          "<br>n_apps: ", n_apps,
          "<br>Evidence types: ", evidence_types,
          "<br>Clusters: ", n_clusters,
          "<br>Cluster IDs: ", clusters,
          "<br>Chromosomes: ", chr_set,
          "<br>Total records: ", n_records
        )
      )
    
    out$plot_label <- factor(out$plot_label, levels = rev(unique(out$plot_label)))
    
    out
  })
  
  output$shared_terms_bubble_plot <- plotly::renderPlotly({
    df <- shared_terms_bubble_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No shared evidence terms available for plotting.")
    )
    
    rescale_size <- function(x, to = c(10, 45)) {
      x <- as.numeric(x)
      rng <- range(x, na.rm = TRUE)
      
      if (!all(is.finite(rng)) || diff(rng) == 0) {
        return(rep(mean(to), length(x)))
      }
      
      to[1] + (x - rng[1]) * (to[2] - to[1]) / diff(rng)
    }
    
    app_sig_levels <- unique(df$apps_signature)
    
    pal <- viridisLite::viridis(length(app_sig_levels))
    names(pal) <- app_sig_levels
    
    p <- plotly::plot_ly(
      source = "shared_terms_bubble_plot"
    )
    
    for (sig in app_sig_levels) {
      subdf <- df %>% dplyr::filter(apps_signature == sig)
      
      p <- p %>%
        plotly::add_markers(
          data = subdf,
          x = ~n_clusters,
          y = ~plot_label,
          text = ~hover_txt,
          hoverinfo = "text",
          name = sig,
          marker = list(
            size = rescale_size(subdf$n_records, to = c(10, 45)),
            sizemode = "diameter",
            color = pal[[sig]],
            opacity = 0.82,
            line = list(
              color = "rgba(80,80,80,0.45)",
              width = 0.7
            )
          )
        )
    }
    
    p %>%
      plotly::layout(
        xaxis = list(
          title = "Number of clusters",
          rangemode = "tozero",
          dtick = 1,
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = levels(df$plot_label)
        ),
        legend = list(
          title = list(text = "Supporting app combination")
        ),
        margin = list(l = 220, r = 30, t = 30, b = 60),
        hoverlabel = list(align = "left")
      )
  })
  
  #--------------------------------------------
  
  gene_bridges_collapsed <- reactive({
    dt <- gene_bridges_combined()
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(empty_gene_bridge())
    }
    
    dt %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::mutate(
        dplyr::across(c(gene, source_app, evidence_type, cluster_id), as.character)
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        source_app    = paste(sort(unique(stats::na.omit(source_app[source_app != ""]))), collapse = "; "),
        evidence_type = paste(sort(unique(stats::na.omit(evidence_type[evidence_type != ""]))), collapse = "; "),
        cluster_id    = paste(sort(unique(stats::na.omit(cluster_id[cluster_id != ""]))), collapse = "; "),
        chr           = paste(sort(unique(stats::na.omit(chr))), collapse = "; "),
        start         = suppressWarnings(min(as.numeric(start), na.rm = TRUE)),
        end           = suppressWarnings(max(as.numeric(end), na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        start = ifelse(is.infinite(start), NA, start),
        end   = ifelse(is.infinite(end), NA, end)
      ) %>%
      sanitize_dt_types()
  })
  
  make_collapsed_links <- function(x, base_url, split_regex = "\\s*[,;]\\s*") {
    vals <- as.character(x)
    
    vapply(vals, function(v) {
      v <- trimws(v)
      if (is.na(v) || !nzchar(v)) return("")
      
      parts <- unlist(strsplit(v, split_regex, perl = TRUE))
      parts <- trimws(parts)
      parts <- parts[!is.na(parts) & nzchar(parts)]
      parts <- unique(parts)
      
      if (!length(parts)) return("")
      
      links <- vapply(parts, function(p) {
        lab <- htmltools::htmlEscape(p)
        href <- paste0(base_url, utils::URLencode(p, reserved = TRUE))
        paste0("<a href='", href, "' target='_blank'>", lab, "</a>")
      }, character(1))
      
      paste(links, collapse = "; ")
    }, character(1))
  }
  
  make_genecards_collapsed_links <- function(x) {
    make_collapsed_links(
      x,
      base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
      split_regex = "\\s*[,;]\\s*"
    )
  }
  
  split_pipe_terms <- function(x) {
    x <- as.character(x %||% "")
    x <- trimws(x)
    if (!nzchar(x)) return(character(0))
    
    parts <- unlist(strsplit(x, "\\s*\\|\\s*", perl = TRUE))
    parts <- trimws(parts)
    parts <- parts[!is.na(parts) & nzchar(parts)]
    unique(parts)
  }
  
  is_omim_term <- function(x) {
    grepl("^\\s*\\[MIM:\\d+\\]", x, perl = TRUE)
  }
  
  extract_omim_id <- function(x) {
    m <- regmatches(x, regexpr("\\[MIM:(\\d+)\\]", x, perl = TRUE))
    sub("^\\[MIM:(\\d+)\\]$", "\\1", m, perl = TRUE)
  }
  
  make_catalog_only_links <- function(x) {
    vals <- as.character(x)
    
    vapply(vals, function(v) {
      terms <- split_pipe_terms(v)
      terms <- terms[!is_omim_term(terms)]
      
      if (!length(terms)) return("")
      
      links <- vapply(terms, function(term) {
        href <- paste0(
          "https://www.ebi.ac.uk/gwas/search?query=",
          utils::URLencode(term, reserved = TRUE)
        )
        paste0(
          "<a href='", href, "' target='_blank'>",
          htmltools::htmlEscape(term),
          "</a>"
        )
      }, character(1))
      
      paste(links, collapse = "; ")
    }, character(1))
  }
  
  make_omim_only_links <- function(x) {
    vals <- as.character(x)
    
    vapply(vals, function(v) {
      terms <- split_pipe_terms(v)
      terms <- terms[is_omim_term(terms)]
      
      if (!length(terms)) return("")
      
      links <- vapply(terms, function(term) {
        omim_id <- extract_omim_id(term)
        href <- paste0("https://www.omim.org/entry/", omim_id)
        
        paste0(
          "<a href='", href, "' target='_blank'>",
          htmltools::htmlEscape(term),
          "</a>"
        )
      }, character(1))
      
      paste(links, collapse = "; ")
    }, character(1))
  }
  
  #########################################################################
  ####################### BRIDGES #########################################
  output$all_gene_bridges_dt <- DT::renderDT({
    dt <- gene_bridges_combined()
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(
        DT::datatable(
          data.frame(Message = "No shared gene bridges found yet."),
          rownames = FALSE
        )
      )
    }
    
    # Assegurem existència de columnes GTEx per evitar errors en apps/sessions antigues
    if (!"gtex_gene_name" %in% names(dt)) dt$gtex_gene_name <- NA_character_
    if (!"gtex_gene_id" %in% names(dt)) dt$gtex_gene_id <- NA_character_
    
    collapse_unique_comma <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = ", "), "\\s*,\\s*|\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = ", ")
    }
    
    dt_sum <- dt %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        gtex_gene_name = trimws(as.character(gtex_gene_name)),
        gtex_gene_id = trimws(as.character(gtex_gene_id)),
        source_app = trimws(as.character(source_app)),
        evidence_type = trimws(as.character(evidence_type)),
        cluster_id = trimws(as.character(cluster_id)),
        chr = as.character(chr)
      ) %>%
      dplyr::group_by(cluster_id, chr, start, end) %>%
      dplyr::summarise(
        n_bridges = dplyr::n(),
        n_genes = dplyr::n_distinct(gene[!is.na(gene) & nzchar(gene)]),
        source_apps = paste(sort(unique(source_app[!is.na(source_app) & nzchar(source_app)])), collapse = ", "),
        evidence_types = paste(sort(unique(evidence_type[!is.na(evidence_type) & nzchar(evidence_type)])), collapse = ", "),
        gene_bridges = paste(sort(unique(gene[!is.na(gene) & nzchar(gene)])), collapse = ", "),
        gtex_gene_names = collapse_unique_comma(gtex_gene_name),
        gtex_gene_ids = collapse_unique_comma(gtex_gene_id),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        chr_order = dplyr::case_when(
          is.na(suppressWarnings(as.integer(chr))) & toupper(chr) %in% c("X", "CHRX") ~ 23L,
          is.na(suppressWarnings(as.integer(chr))) & toupper(chr) %in% c("Y", "CHRY") ~ 24L,
          is.na(suppressWarnings(as.integer(chr))) & toupper(chr) %in% c("MT", "M", "CHRMT", "CHRM") ~ 26L,
          TRUE ~ suppressWarnings(as.integer(chr))
        )
      ) %>%
      dplyr::arrange(chr_order, start, end) %>%
      dplyr::select(-chr_order) %>%
      dplyr::mutate(
        gene_bridges = make_genecards_collapsed_links(gene_bridges),
        gtex_gene_names = make_collapsible_text(gtex_gene_names, max_chars = 120),
        gtex_gene_ids = make_collapsible_text(gtex_gene_ids, max_chars = 120)
      )
    
    gene_col <- which(names(dt_sum) == "gene_bridges") - 1L
    gtex_name_col <- which(names(dt_sum) == "gtex_gene_names") - 1L
    gtex_id_col <- which(names(dt_sum) == "gtex_gene_ids") - 1L
    
    escape_cols <- setdiff(
      seq_along(dt_sum) - 1L,
      c(gene_col, gtex_name_col, gtex_id_col)
    )
    
    DT::datatable(
      sanitize_dt_types(dt_sum),
      rownames = FALSE,
      escape = escape_cols,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = FALSE,
        columnDefs = list(
          list(width = "320px", targets = gene_col),
          list(width = "220px", targets = gtex_name_col),
          list(width = "220px", targets = gtex_id_col)
        )
      )
    )
  }, server = TRUE)
  
  output$all_term_bridges_dt <- DT::renderDT({
    
    dt <- term_bridges_combined()
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(
        DT::datatable(
          data.frame(Message = "No shared term bridges found yet."),
          rownames = FALSE,
          options = list(dom = "t")
        )
      )
    }
    
    # ------------------------------------------------------------
    # 1) Base neta
    # ------------------------------------------------------------
    dt_clean <- dt %>%
      dplyr::mutate(
        term = clean_term_label_strict(term),
        term_key = normalize_term_key(term),
        term_type = trimws(as.character(term_type)),
        source_app = trimws(as.character(source_app)),
        evidence_type = trimws(as.character(evidence_type)),
        cluster_id = trimws(as.character(cluster_id)),
        chr = as.character(chr)
      ) %>%
      dplyr::filter(
        !is.na(term), nzchar(term),
        !is.na(term_key), nzchar(term_key)
      )
    
    # ------------------------------------------------------------
    # 2) Un sol label representatiu per terme normalitzat dins cluster
    #    Això elimina duplicacions tipus:
    #    "Smoking initiation" vs "smoking initiation"
    # ------------------------------------------------------------
    dt_display <- dt_clean %>%
      dplyr::group_by(cluster_id, chr, start, end, term_key) %>%
      dplyr::summarise(
        term_display = sort(unique(term))[1],
        has_omim = any(
          evidence_type %in% c("mim_disease", "omim_disease") |
            source_app %in% c("omim"),
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    # ------------------------------------------------------------
    # 3) Resum per cluster
    # ------------------------------------------------------------
    dt_sum <- dt_clean %>%
      dplyr::group_by(cluster_id, chr, start, end) %>%
      dplyr::summarise(
        n_terms = dplyr::n(),
        n_unique_terms = dplyr::n_distinct(term),
        n_terms_clean = dplyr::n_distinct(term_key),
        term_types = paste(sort(unique(term_type[!is.na(term_type) & nzchar(term_type)])), collapse = "; "),
        source_apps = paste(sort(unique(source_app[!is.na(source_app) & nzchar(source_app)])), collapse = "; "),
        evidence_types = paste(sort(unique(evidence_type[!is.na(evidence_type) & nzchar(evidence_type)])), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        dt_display %>%
          dplyr::group_by(cluster_id, chr, start, end) %>%
          dplyr::summarise(
            shared_terms_raw = collapse_unique_terms(term_display),
            shared_terms_omim_raw = collapse_unique_terms(term_display[has_omim]),
            .groups = "drop"
          ),
        by = c("cluster_id", "chr", "start", "end")
      ) %>%
      dplyr::arrange(
        suppressWarnings(as.integer(chr)),
        start,
        end
      ) %>%
      dplyr::mutate(
        shared_terms = make_terms_collapsed_links(
          shared_terms_raw,
          base_url = "https://www.ebi.ac.uk/gwas/search?query="
        ),
        shared_terms_omim = make_omim_collapsed_links(shared_terms_omim_raw)
      ) %>%
      dplyr::select(
        cluster_id, chr, start, end,
        n_terms, n_unique_terms, n_terms_clean,
        term_types, source_apps, evidence_types,
        shared_terms, shared_terms_omim
      )
    
    shared_terms_col <- which(names(dt_sum) == "shared_terms") - 1L
    omim_col <- which(names(dt_sum) == "shared_terms_omim") - 1L
    
    DT::datatable(
      sanitize_dt_types(dt_sum),
      rownames = FALSE,
      escape = setdiff(seq_along(dt_sum) - 1L, c(shared_terms_col, omim_col)),
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          ),
          list(
            extend = "csv",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          ),
          list(
            extend = "excel",
            exportOptions = list(
              modifier = list(
                page = "all",
                search = "none",
                order = "index"
              )
            )
          )
        ),
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = FALSE,
        columnDefs = list(
          list(width = "420px", targets = shared_terms_col),
          list(width = "420px", targets = omim_col)
        )
      )
    )
  }, server = TRUE)
  
  all_gene_bridges <- reactive({
    idir <- selected_session_dir()
    
    dplyr::bind_rows(
      read_bridge_rds(file.path(idir, "catalog_gene_bridge.rds"), empty_gene_bridge()),
      read_bridge_rds(file.path(idir, "gtex_gene_bridge.rds"), empty_gene_bridge()),
      read_bridge_rds(file.path(idir, "nonsyn_gene_bridge.rds"), empty_gene_bridge()),
      read_bridge_rds(file.path(idir, "ewasdis_gene_bridge.rds"), empty_gene_bridge()),
      read_bridge_rds(file.path(idir, "ewastum_gene_bridge.rds"), empty_gene_bridge()),
      read_bridge_rds(file.path(idir, "ld_gene_bridge.rds"), empty_gene_bridge())
    ) %>%
      sanitize_bridge() %>%
      dplyr::distinct()
  })
  
  gene_bridges_combined <- reactive({
    manual_ge <- tryCatch(genes_integrated(), error = function(e) tibble::tibble())
    
    shared_ge <- tryCatch(all_gene_bridges(), error = function(e) empty_gene_bridge())
    
    if (!is.data.frame(manual_ge) || !nrow(manual_ge)) {
      return(shared_ge %>% sanitize_bridge() %>% dplyr::distinct())
    }
    
    # Si genes_integrated() no té exactament el schema bridge, el normalitzem
    if (!all(c("gene", "source_app", "evidence_type", "cluster_id", "chr", "start", "end") %in% names(manual_ge))) {
      manual_ge <- tibble::tibble()
    }
    
    dplyr::bind_rows(shared_ge, manual_ge) %>%
      sanitize_bridge() %>%
      dplyr::distinct()
  })
    
  all_term_bridges <- reactive({
    idir <- selected_session_dir()
    
    dplyr::bind_rows(
      read_bridge_rds(file.path(idir, "catalog_term_bridge.rds"), empty_term_bridge()),
      read_bridge_rds(file.path(idir, "gtex_term_bridge.rds"), empty_term_bridge()),
      read_bridge_rds(file.path(idir, "nonsyn_term_bridge.rds"), empty_term_bridge()),
      read_bridge_rds(file.path(idir, "ewasdis_term_bridge.rds"), empty_term_bridge()),
      read_bridge_rds(file.path(idir, "ewastum_term_bridge.rds"), empty_term_bridge())
    ) %>%
      sanitize_bridge() %>%
      dplyr::distinct()
  })
  
  term_bridges_combined <- reactive({
    
    manual_tr <- tryCatch(traits_integrated(), error = function(e) tibble::tibble())
    shared_tr <- tryCatch(all_term_bridges(), error = function(e) empty_term_bridge())
    
    if (is.data.frame(manual_tr) && nrow(manual_tr)) {
      
      if (all(c("trait", "source_app", "evidence_type", "cluster_id", "chr", "start", "end") %in% names(manual_tr))) {
        manual_tr <- manual_tr %>%
          dplyr::transmute(
            term = as.character(trait),
            term_type = "trait",
            source_app = as.character(source_app),
            evidence_type = as.character(evidence_type),
            cluster_id = as.character(cluster_id),
            chr = chr,
            start = start,
            end = end
          )
      } else if (!all(c("term", "term_type", "source_app", "evidence_type", "cluster_id", "chr", "start", "end") %in% names(manual_tr))) {
        manual_tr <- tibble::tibble()
      }
    }
    
    dplyr::bind_rows(shared_tr, manual_tr) %>%
      sanitize_bridge() %>%
      clean_term_bridge_records()
  })
  
  #########################################################################
  # nova adaptacio: 
  block_ranges_global <- reactive({
    dd <- ld_details_rds()
    
    if (!is.list(dd) || !length(dd)) {
      return(tibble::tibble(
        cluster_id = character(),
        chr = integer(),
        block_id = character(),
        block_label = character(),
        i = integer(),
        j = integer(),
        block_start = integer(),
        block_end = integer(),
        block_size_bp = integer()
      ))
    }
    
    xx <- dplyr::bind_rows(lapply(dd, function(x) {
      if (is.null(x) || !is.list(x) || !is.data.frame(x$block_ranges) || !nrow(x$block_ranges)) {
        return(NULL)
      }
      
      br <- x$block_ranges
      
      if (!"cluster_id" %in% names(br)) br$cluster_id <- x$cluster_id %||% NA_character_
      if (!"chr" %in% names(br)) br$chr <- suppressWarnings(as.integer(x$chr %||% NA))
      if (!"block_label" %in% names(br)) br$block_label <- safe_block_label_from_id(
        block_id = br$block_id,
        cluster_id = br$cluster_id
      )
      
      br %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          chr = suppressWarnings(as.integer(chr)),
          block_id = as.character(block_id),
          block_label = as.character(block_label),
          i = if ("i" %in% names(.)) suppressWarnings(as.integer(i)) else NA_integer_,
          j = if ("j" %in% names(.)) suppressWarnings(as.integer(j)) else NA_integer_,
          block_start = suppressWarnings(as.integer(block_start)),
          block_end = suppressWarnings(as.integer(block_end)),
          block_size_bp = suppressWarnings(as.integer(block_size_bp))
        )
    }))
    
    if (!is.data.frame(xx) || !nrow(xx)) {
      return(tibble::tibble(
        cluster_id = character(),
        chr = integer(),
        block_id = character(),
        block_label = character(),
        i = integer(),
        j = integer(),
        block_start = integer(),
        block_end = integer(),
        block_size_bp = integer()
      ))
    }
    
    xx %>%
      dplyr::filter(!is.na(block_id), nzchar(block_id)) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE) %>%
      dplyr::arrange(cluster_id, block_start, block_end)
  })
  
  #########################################################################
  observeEvent(input$run_ld_global, {
    
    ld_global_log_buf(character(0))
    output$ld_global_status <- renderText("Preparing global LD batch...")
    
    tryCatch({
      cl <- clusters_consensus()
      ca <- candidates_consensus()
      
      validate(
        need(is.data.frame(cl) && nrow(cl) > 0, "No integrated clusters available."),
        need(is.data.frame(ca) && nrow(ca) > 0, "No integrated candidates available.")
      )
      
      plink_bin <- gi_ld_defaults()$plink %||% ""
      bfile_ref <- gi_ld_defaults()$bfile %||% ""
      ld_pops_dir <- gi_ld_defaults()$popdir %||% ""
      
      validate(
        need(file.exists(plink_bin), paste0("PLINK binary not found: ", plink_bin)),
        need(nzchar(bfile_ref), "Reference bfile prefix is empty."),
        need(
          file.exists(paste0(bfile_ref, ".bed")) &&
            file.exists(paste0(bfile_ref, ".bim")) &&
            file.exists(paste0(bfile_ref, ".fam")),
          "Reference bfile prefix invalid (missing .bed/.bim/.fam)."
        )
      )
      
      keep_path <- read_keep_file_global(input$ld_global_pop, ld_pops_dir)
      
      withProgress(message = "Running global LD / blocks", value = 0, {
        
        append_ld_global_log("[GLOBAL LD] population = ", input$ld_global_pop)
        append_ld_global_log("[GLOBAL LD] metric = ", input$ld_global_metric %||% "R2")
        append_ld_global_log("[GLOBAL LD] min threshold = ", input$ld_global_r2_min %||% 0.6)
        append_ld_global_log("[GLOBAL LD] compute blocks = TRUE")
        append_ld_global_log("[GLOBAL LD] clusters = ", nrow(cl))
        
        
        sum_path <- file.path(selected_session_dir(), "ld_cluster_summary.rds")
        det_path <- file.path(selected_session_dir(), "ld_cluster_details.rds")
        
        idir <- selected_session_dir()
        
        
        if (file.exists(sum_path)) unlink(sum_path)
        if (file.exists(det_path)) unlink(det_path)
        
        out <- run_ld_global_for_clusters(
          clusters_df = cl,
          candidates_df = ca,
          out_summary_rds = sum_path,
          out_details_rds = det_path,
          bfile_ref = bfile_ref,
          keep_path = keep_path,
          plink_bin = plink_bin,
          workdir = tempdir(),
          pop = input$ld_global_pop,
          ld_metric = input$ld_global_metric %||% "R2",
          r2_min = as.numeric(input$ld_global_r2_min %||% 0.6),
          compute_blocks = TRUE,
          max_snps_interval = 400,
          progress_fun = function(i, n, cluster_id, stage) {
            if (stage == "start") {
              detail_txt <- paste0("Cluster ", i, "/", n, ": ", cluster_id)
              incProgress(0, detail = detail_txt)
              output$ld_global_status <- renderText(
                paste0("Running global LD... cluster ", i, " of ", n, " (", cluster_id, ")")
              )
              append_ld_global_log("[START] ", cluster_id, " (", i, "/", n, ")")
            }
            
            if (stage == "done") {
              incProgress(1 / max(1, n), detail = paste0("Completed ", cluster_id))
              append_ld_global_log("[DONE] ", cluster_id, " (", i, "/", n, ")")
            }
          }
        )
        
        ld_summary_version(ld_summary_version() + 1)
        ld_details_version(ld_details_version() + 1)
        
        append_ld_global_log("[GLOBAL LD] finished")
        
        output$ld_global_status <- renderText(
          paste0("Global LD summary generated for ", nrow(out$summary), " clusters.")
        )
      })
      
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (!nzchar(msg)) msg <- "<empty error>"
      append_ld_global_log("[ERROR] ", msg)
      output$ld_global_status <- renderText(paste("LD batch failed:", msg))
    })
  })
  
  ########### CONSENSUS LOCI TABLE ######################
  block_support_by_cluster <- reactive({
    blk <- block_overlap_summary_global_df()
    
    if (!is.data.frame(blk) || !nrow(blk)) {
      return(tibble::tibble(
        cluster_id = character(),
        block_support_any = integer(),
        block_support_blocks = integer(),
        block_support_apps = integer(),
        block_support_hits = integer(),
        block_catalog_hits = integer(),
        block_gtex_hits = integer(),
        block_nonsyn_hits = integer(),
        block_ewasdis_hits = integer(),
        block_ewastum_hits = integer()
      ))
    }
    
    for (nm in c(
      "cluster_id",
      "block_id",
      "block_support_any",
      "block_support_apps",
      "block_support_hits",
      "block_catalog_hits",
      "block_gtex_hits",
      "block_nonsyn_hits",
      "block_ewasdis_hits",
      "block_ewastum_hits"
    )) {
      if (!nm %in% names(blk)) blk[[nm]] <- 0
    }
    
    blk %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_support_any = dplyr::coalesce(as.integer(block_support_any), 0L),
        block_support_apps = dplyr::coalesce(as.integer(block_support_apps), 0L),
        block_support_hits = dplyr::coalesce(as.integer(block_support_hits), 0L),
        block_catalog_hits = dplyr::coalesce(as.integer(block_catalog_hits), 0L),
        block_gtex_hits = dplyr::coalesce(as.integer(block_gtex_hits), 0L),
        block_nonsyn_hits = dplyr::coalesce(as.integer(block_nonsyn_hits), 0L),
        block_ewasdis_hits = dplyr::coalesce(as.integer(block_ewasdis_hits), 0L),
        block_ewastum_hits = dplyr::coalesce(as.integer(block_ewastum_hits), 0L)
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        block_support_any = as.integer(any(block_support_any > 0, na.rm = TRUE)),
        block_support_blocks = dplyr::n_distinct(block_id[block_support_any > 0]),
        block_support_apps = suppressWarnings(max(block_support_apps, na.rm = TRUE)),
        block_support_hits = sum(block_support_hits, na.rm = TRUE),
        block_catalog_hits = sum(block_catalog_hits, na.rm = TRUE),
        block_gtex_hits = sum(block_gtex_hits, na.rm = TRUE),
        block_nonsyn_hits = sum(block_nonsyn_hits, na.rm = TRUE),
        block_ewasdis_hits = sum(block_ewasdis_hits, na.rm = TRUE),
        block_ewastum_hits = sum(block_ewastum_hits, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        block_support_any = dplyr::coalesce(as.integer(block_support_any), 0L),
        block_support_blocks = dplyr::coalesce(as.integer(block_support_blocks), 0L),
        block_support_apps = ifelse(
          is.infinite(block_support_apps),
          0L,
          dplyr::coalesce(as.integer(block_support_apps), 0L)
        ),
        block_support_hits = dplyr::coalesce(as.integer(block_support_hits), 0L),
        block_catalog_hits = dplyr::coalesce(as.integer(block_catalog_hits), 0L),
        block_gtex_hits = dplyr::coalesce(as.integer(block_gtex_hits), 0L),
        block_nonsyn_hits = dplyr::coalesce(as.integer(block_nonsyn_hits), 0L),
        block_ewasdis_hits = dplyr::coalesce(as.integer(block_ewasdis_hits), 0L),
        block_ewastum_hits = dplyr::coalesce(as.integer(block_ewastum_hits), 0L)
      )
  })
  
  # ============================================================================= 
  # ============================================================
  # PRIORITIZED CLUSTERS
  # ============================================================
  
  ################################################################################
  # ---------------------------
  # Prioritized clusters v2
  # ---------------------------
  prioritized_cluster_df_v2 <- reactive({
    cl  <- priority_source_clusters()
    gh  <- gwas_hit_priority_df_v2()
    blk <- prioritized_block_df_v2()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available.")
    )
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    # ============================================================
    # 1) BASE CANÒNICA: 1 fila per cluster_id
    # ============================================================
    cl_base <- cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        chr = as.character(chr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(chr), nzchar(chr),
        is.finite(start),
        is.finite(end)
      ) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        chr = dplyr::first(chr[!is.na(chr) & nzchar(chr)]),
        start = min(start, na.rm = TRUE),
        end = max(end, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        !is.na(chr), nzchar(chr),
        is.finite(start),
        is.finite(end)
      ) %>%
      dplyr::mutate(
        cluster_size_bp = end - start + 1,
        cluster_size_kb = round(cluster_size_bp / 1000, 2)
      )
    
    # ============================================================
    # 2) DESCRIPTORS DE GWAS HITS
    #    (informatius / auditoria, no entren directament al score)
    # ============================================================
    hit_sum <- if (is.data.frame(gh) && nrow(gh) > 0) {
      gh %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          gwas_hit = as.character(gwas_hit),
          gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
          gwas_hit_score = dplyr::coalesce(as.numeric(gwas_hit_score), 0),
          apps_supported = dplyr::coalesce(as.character(gwas_hit_apps), "")
        ) %>%
        dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          n_gwas_hits = dplyr::n_distinct(paste(cluster_id, gwas_hit, gwas_pos, sep = "||")),
          top_gwas_hit_score = round(max(gwas_hit_score, na.rm = TRUE), 2),
          apps_supported = collapse_unique_semicolon(apps_supported),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          n_apps_supported = dplyr::if_else(
            nzchar(apps_supported),
            lengths(strsplit(apps_supported, "\\s*;\\s*")),
            0L
          )
        )
    } else {
      tibble::tibble(
        cluster_id = character(),
        n_gwas_hits = integer(),
        top_gwas_hit_score = numeric(),
        apps_supported = character(),
        n_apps_supported = integer()
      )
    }
    
    # ============================================================
    # 3) COMPONENT CANÒNIC DES DE BLOCKS
    #    cluster_score_from_blocks =
    #      top_block_score + 0.05 * other_block_score
    # ============================================================
    block_sum <- if (is.data.frame(blk) && nrow(blk) > 0) {
      blk %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          block_score = dplyr::coalesce(as.numeric(block_score), 0),
          block_score_from_hits = dplyr::coalesce(as.numeric(block_score_from_hits), 0),
          block_gene_content_bonus = dplyr::coalesce(
            as.numeric(block_gene_content_bonus),
            as.numeric(block_bio_bonus),
            0
          ),
          n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
          genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
        ) %>%
        dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          n_blocks = dplyr::n_distinct(block_id[!is.na(block_id) & nzchar(block_id)]),
          n_blocks_with_genes = sum(n_genes_in_block > 0, na.rm = TRUE),
          block_ids = paste(sort(unique(block_id[!is.na(block_id) & nzchar(block_id)])), collapse = "; "),
          genes_in_blocks = collapse_unique_semicolon(genes_in_block),
          
          top_block_score = round(max(block_score, na.rm = TRUE), 2),
          other_block_score = round(sum(block_score, na.rm = TRUE) - max(block_score, na.rm = TRUE), 2),
          
          top_block_score_from_hits = round(max(block_score_from_hits, na.rm = TRUE), 2),
          other_block_score_from_hits = round(sum(block_score_from_hits, na.rm = TRUE) - max(block_score_from_hits, na.rm = TRUE), 2),
          
          sum_block_scores = round(sum(block_score, na.rm = TRUE), 2),
          sum_block_scores_from_hits = round(sum(block_score_from_hits, na.rm = TRUE), 2),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          other_block_score = dplyr::if_else(is.finite(other_block_score) & other_block_score > 0, other_block_score, 0),
          other_block_score_from_hits = dplyr::if_else(is.finite(other_block_score_from_hits) & other_block_score_from_hits > 0, other_block_score_from_hits, 0),
          
          n_unique_genes_in_blocks = vapply(genes_in_blocks, count_semicolon_items, integer(1)),
          
          cluster_score_from_blocks = round(
            top_block_score + 0.05 * other_block_score,
            2
          ),
          
          cluster_gene_bonus = round(
            dplyr::if_else(n_unique_genes_in_blocks > 0, 2.5, 0) +
              0.75 * log1p(n_unique_genes_in_blocks),
            2
          ),
          
          cluster_fragmentation_bonus = round(
            0.5 * log1p(n_blocks),
            2
          )
        )
    } else {
      tibble::tibble(
        cluster_id = character(),
        n_blocks = integer(),
        n_blocks_with_genes = integer(),
        block_ids = character(),
        genes_in_blocks = character(),
        top_block_score = numeric(),
        other_block_score = numeric(),
        top_block_score_from_hits = numeric(),
        other_block_score_from_hits = numeric(),
        sum_block_scores = numeric(),
        sum_block_scores_from_hits = numeric(),
        n_unique_genes_in_blocks = integer(),
        cluster_score_from_blocks = numeric(),
        cluster_gene_bonus = numeric(),
        cluster_fragmentation_bonus = numeric()
      )
    }
    
    # ============================================================
    # 4) RESULTAT FINAL
    # ============================================================
    out <- cl_base %>%
      dplyr::left_join(hit_sum, by = "cluster_id") %>%
      dplyr::left_join(block_sum, by = "cluster_id") %>%
      dplyr::mutate(
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_blocks_with_genes = dplyr::coalesce(as.integer(n_blocks_with_genes), 0L),
        block_ids = dplyr::coalesce(as.character(block_ids), ""),
        genes_in_blocks = dplyr::coalesce(as.character(genes_in_blocks), ""),
        
        top_block_score = dplyr::coalesce(as.numeric(top_block_score), 0),
        other_block_score = dplyr::coalesce(as.numeric(other_block_score), 0),
        top_block_score_from_hits = dplyr::coalesce(as.numeric(top_block_score_from_hits), 0),
        other_block_score_from_hits = dplyr::coalesce(as.numeric(other_block_score_from_hits), 0),
        sum_block_scores = dplyr::coalesce(as.numeric(sum_block_scores), 0),
        sum_block_scores_from_hits = dplyr::coalesce(as.numeric(sum_block_scores_from_hits), 0),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        
        cluster_score_from_blocks = dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0),
        cluster_gene_bonus = dplyr::coalesce(as.numeric(cluster_gene_bonus), 0),
        cluster_fragmentation_bonus = dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0)
      ) %>%
      dplyr::mutate(
        cluster_score = round(
          cluster_score_from_blocks +
            cluster_gene_bonus +
            cluster_fragmentation_bonus,
          2
        ),
        priority_class = dplyr::case_when(
          cluster_score >= 30 ~ "High",
          cluster_score >= 10 ~ "Medium",
          TRUE ~ "Low"
        ),
        priority_class_relative = classify_priority_tertiles(cluster_score)
      ) %>%
      dplyr::arrange(
        dplyr::desc(cluster_score),
        dplyr::desc(top_block_score),
        cluster_id
      )
    
    out
  })
  
  ################################################################################
  
  ################################################################################
  
  
  #---------------- bar plot---------------------------------------
  
  output$prioritized_cluster_dt <- DT::renderDT({
    df <- prioritized_cluster_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized clusters available.")
    )
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    # ------------------------------------------------------------
    # Garantir columnes mínimes esperades
    # ------------------------------------------------------------
    needed_chr <- c(
      "cluster_id",
      "priority_class",
      "priority_class_relative",
      "genes_in_blocks",
      "apps_supported",
      "support_signature"
    )
    
    needed_num <- c(
      "cluster_score",
      "cluster_score_from_blocks",
      "top_block_score",
      "other_block_score",
      "cluster_gene_bonus",
      "cluster_fragmentation_bonus",
      "top_gwas_hit_score",
      "cluster_size_kb",
      "n_gwas_hits",
      "n_blocks",
      "n_blocks_with_genes",
      "n_unique_genes_in_blocks",
      "n_apps_supported"
    )
    
    for (nm in needed_chr) {
      if (!nm %in% names(df)) df[[nm]] <- ""
    }
    
    for (nm in needed_num) {
      if (!nm %in% names(df)) df[[nm]] <- 0
    }
    
    show_df <- df %>%
      dplyr::mutate(
        cluster_id = dplyr::coalesce(as.character(cluster_id), ""),
        genes_in_blocks = dplyr::coalesce(as.character(genes_in_blocks), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        support_signature = dplyr::coalesce(as.character(support_signature), ""),
        
        cluster_score = round(dplyr::coalesce(as.numeric(cluster_score), 0), 2),
        cluster_score_from_blocks = round(dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0), 2),
        top_block_score = round(dplyr::coalesce(as.numeric(top_block_score), 0), 2),
        other_block_score = round(dplyr::coalesce(as.numeric(other_block_score), 0), 2),
        cluster_gene_bonus = round(dplyr::coalesce(as.numeric(cluster_gene_bonus), 0), 2),
        cluster_fragmentation_bonus = round(dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0), 2),
        top_gwas_hit_score = round(dplyr::coalesce(as.numeric(top_gwas_hit_score), 0), 2),
        cluster_size_kb = round(dplyr::coalesce(as.numeric(cluster_size_kb), 0), 2),
        
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_blocks_with_genes = dplyr::coalesce(as.integer(n_blocks_with_genes), 0L),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        
        n_genes_from_label = vapply(genes_in_blocks, count_semicolon_items, integer(1)),
        n_unique_genes_in_blocks = dplyr::if_else(
          n_unique_genes_in_blocks > 0L,
          n_unique_genes_in_blocks,
          as.integer(n_genes_from_label)
        ),
        
        evidence_summary = paste0(
          n_gwas_hits, " GWAS hits",
          " | ",
          n_blocks, " blocks",
          " | ",
          n_unique_genes_in_blocks, " genes",
          " | ",
          n_apps_supported, " apps"
        ),
        
        support_signature = dplyr::case_when(
          nzchar(support_signature) ~ support_signature,
          nzchar(apps_supported) ~ apps_supported,
          TRUE ~ ""
        ),
        
        score_breakdown = paste0(
          "total=", format(cluster_score, nsmall = 2),
          " | top_block=", format(top_block_score, nsmall = 2),
          " | other_blocks=", format(0.05 * other_block_score, nsmall = 2),
          " | gene_bonus=", format(cluster_gene_bonus, nsmall = 2),
          " | frag_bonus=", format(cluster_fragmentation_bonus, nsmall = 2)
        ),
        
        genes_in_blocks = make_genecards_links(genes_in_blocks),
        priority_class = make_priority_badge(priority_class),
        priority_class_relative = make_priority_badge(priority_class_relative)
      ) %>%
      dplyr::select(
        cluster_id,
        score = cluster_score,
        `priority_class (ABS)` = priority_class,
        `priority_class (REL)` = priority_class_relative,
        evidence_summary,
        support_signature,
        score_breakdown,
        genes_in_blocks,
        cluster_score_from_blocks,
        top_block_score,
        other_block_score,
        cluster_gene_bonus,
        cluster_fragmentation_bonus,
        top_gwas_hit_score,
        cluster_size_kb,
        n_gwas_hits,
        n_blocks,
        n_blocks_with_genes,
        n_genes = n_unique_genes_in_blocks,
        n_apps_supported,
        apps_supported
      )
    
    hidden_cols <- which(names(show_df) %in% c(
      "cluster_score_from_blocks",
      "top_block_score",
      "other_block_score",
      "cluster_gene_bonus",
      "cluster_fragmentation_bonus",
      "top_gwas_hit_score",
      "cluster_size_kb",
      "n_gwas_hits",
      "n_blocks",
      "n_blocks_with_genes",
      "n_genes",
      "n_apps_supported",
      "apps_supported"
    )) - 1L
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        columnDefs = list(
          list(
            targets = hidden_cols,
            visible = FALSE
          )
        ),
        buttons = list(
          list(
            extend = "copyHtml5",
            text = "Copy",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "csvHtml5",
            text = "CSV",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "excelHtml5",
            text = "Excel",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          )
        )
      )
    )
  }, server = FALSE)
  
  # ============================================================
  # PRIORITIZED CLUSTERS PLOT DF (v2)
  # ============================================================
  ##
  prioritized_clusters_plot_df <- reactive({
    df <- prioritized_cluster_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized clusters available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    ord_df <- df %>%
      dplyr::mutate(
        cluster_id_chr = as.character(cluster_id),
        chr_num = extract_chr_num(cluster_id_chr),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id_chr, "(?<=_)\\d+$"))),
        
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0),
        cluster_score_from_blocks = dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0),
        top_block_score = dplyr::coalesce(as.numeric(top_block_score), 0),
        other_block_score = dplyr::coalesce(as.numeric(other_block_score), 0),
        cluster_gene_bonus = dplyr::coalesce(as.numeric(cluster_gene_bonus), 0),
        cluster_fragmentation_bonus = dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0),
        
        priority_class = dplyr::coalesce(as.character(priority_class), ""),
        priority_class_relative = dplyr::coalesce(as.character(priority_class_relative), ""),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_blocks_with_genes = dplyr::coalesce(as.integer(n_blocks_with_genes), 0L),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        genes_in_blocks = dplyr::coalesce(as.character(genes_in_blocks), ""),
        cluster_size_kb = dplyr::coalesce(as.numeric(cluster_size_kb), 0)
      ) %>%
      dplyr::mutate(
        n_genes_from_label = vapply(genes_in_blocks, count_semicolon_items, integer(1)),
        n_unique_genes_in_blocks = dplyr::if_else(
          n_unique_genes_in_blocks > 0L,
          n_unique_genes_in_blocks,
          n_genes_from_label
        ),
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999)
      ) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id_chr)
    
    cluster_levels <- rev(unique(ord_df$cluster_id_chr))
    
    score_vals <- ord_df$cluster_score
    rng <- range(score_vals, na.rm = TRUE)
    
    if (!all(is.finite(rng)) || diff(rng) == 0) {
      fill_fun <- function(x) rep("#fe9929", length(x))
    } else {
      pal_vals <- grDevices::colorRampPalette(
        c("#d9d9d9", "#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404")
      )(100)
      
      fill_fun <- function(x) {
        idx <- floor((x - rng[1]) / (rng[2] - rng[1]) * 99) + 1
        idx <- pmax(1, pmin(100, idx))
        pal_vals[idx]
      }
    }
    
    ord_df %>%
      dplyr::mutate(
        cluster_id = factor(cluster_id_chr, levels = cluster_levels),
        fill_col = fill_fun(cluster_score),
        hover_txt = paste0(
          "Cluster: ", cluster_id_chr,
          "<br>Cluster score: ", round(cluster_score, 2),
          "<br>ABS class: ", priority_class,
          "<br>REL class: ", priority_class_relative,
          "<br><br><b>Canonical score components</b>",
          "<br>Top block support: ", round(top_block_score, 2),
          "<br>Other block support (x0.05): ", round(0.05 * other_block_score, 2),
          "<br>Gene-content bonus: ", round(cluster_gene_bonus, 2),
          "<br>Fragmentation bonus: ", round(cluster_fragmentation_bonus, 2),
          "<br><br><b>Context</b>",
          "<br>GWAS hits: ", n_gwas_hits,
          "<br>Top GWAS-hit score: ", round(top_gwas_hit_score, 2),
          "<br>Apps supported: ", n_apps_supported,
          "<br>Blocks: ", n_blocks,
          "<br>Blocks with genes: ", n_blocks_with_genes,
          "<br>Unique genes in blocks: ", n_unique_genes_in_blocks,
          ifelse(
            nzchar(apps_supported),
            paste0("<br>Apps: ", apps_supported),
            ""
          ),
          ifelse(
            nzchar(genes_in_blocks),
            paste0("<br>Block genes: ", genes_in_blocks),
            ""
          )
        )
      )
  })

  # ============================================================
  # PRIORITIZED CLUSTERS PLOT (size on X, color by score)
  # ============================================================
  output$prioritized_clusters_plot <- plotly::renderPlotly({
    df <- prioritized_clusters_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized clusters plot data available."),
      need("cluster_size_kb" %in% names(df), "cluster_size_kb column not available.")
    )
    
    df <- df %>%
      dplyr::mutate(
        cluster_size_kb = dplyr::coalesce(as.numeric(cluster_size_kb), 0)
      )
    
    score_vals <- df$cluster_score
    rng <- range(score_vals, na.rm = TRUE)
    
    if (!all(is.finite(rng)) || diff(rng) == 0) {
      df$fill_col <- "#fe9929"
    } else {
      pal_fun <- scales::col_numeric(
        palette = c("#d9d9d9", "#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
        domain = rng,
        na.color = "#d9d9d9"
      )
      df$fill_col <- pal_fun(score_vals)
    }
    
    plotly::plot_ly(
      data = df,
      source = "prioritized_clusters_plot",
      y = ~cluster_id,
      x = ~cluster_size_kb,
      type = "bar",
      orientation = "h",
      customdata = ~cluster_id_chr,
      hovertext = ~hover_txt,
      hoverinfo = "text",
      marker = list(
        color = df$fill_col,
        line = list(color = "black", width = 0.5)
      ),
      showlegend = FALSE
    ) %>%
      plotly::layout(
        clickmode = "event+select",
        xaxis = list(
          title = "Cluster size (kb)",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Cluster",
          categoryorder = "array",
          categoryarray = levels(df$cluster_id),
          automargin = TRUE
        ),
        margin = list(l = 90, r = 20, t = 20, b = 60)
      )
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "prioritized_clusters_plot"), {
    ed <- plotly::event_data("plotly_click", source = "prioritized_clusters_plot")
    req(ed)
    
    cid <- ed$customdata[[1]] %||% ""
    req(nzchar(cid))
    
    session$sendCustomMessage(
      "filter_prioritized_cluster_table_by_cluster",
      list(cluster_id = cid)
    )
  })
  
  # ============================================================
  # CLUSTER SCORE LEGEND (v2)
  # ============================================================
  prioritized_clusters_legend_ui <- function(df) {
    if (!is.data.frame(df) || !nrow(df) || !"cluster_score" %in% names(df)) {
      return(NULL)
    }
    
    vals <- suppressWarnings(as.numeric(df$cluster_score))
    vals <- vals[is.finite(vals)]
    
    if (!length(vals)) {
      return(NULL)
    }
    
    vmin <- min(vals, na.rm = TRUE)
    vmax <- max(vals, na.rm = TRUE)
    
    tags$div(
      style = "margin-top:8px; font-size:12px; color:#444;",
      tags$span(
        style = "font-weight:600; margin-right:10px;",
        "Canonical cluster priority score"
      ),
      tags$span(
        style = paste0(
          "display:inline-block; width:220px; height:12px; vertical-align:middle; ",
          "border:1px solid #999; border-radius:6px; margin-right:8px; ",
          "background: linear-gradient(to right, ",
          "#d9d9d9 0%, ",
          "#fff7bc 20%, ",
          "#fec44f 40%, ",
          "#fe9929 60%, ",
          "#d95f0e 80%, ",
          "#993404 100%);"
        )
      ),
      tags$span(formatC(vmin, format = "f", digits = 2)),
      tags$span(" – "),
      tags$span(formatC(vmax, format = "f", digits = 2))
    )
  }
  
  output$prioritized_clusters_legend <- renderUI({
    df <- prioritized_clusters_plot_df()
    prioritized_clusters_legend_ui(df)
  })
  
  ################################################################### 
  # ============================================================
  # CLUSTER PRIORITY COMPONENTS · v2
  # ============================================================
  observeEvent(input$show_cluster_priority_plot, {
    df <- prioritized_cluster_df_v2()
    req(is.data.frame(df), nrow(df) > 0)
    
    showModal(
      modalDialog(
        title = "Cluster priority plot",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$div(
          style = "font-size:13px; color:#444; margin-bottom:10px;",
          HTML(
            paste0(
              "Each horizontal bar represents one cluster. ",
              "Bar segments show the contribution of the current canonical cluster-level priority components: ",
              "<span style='display:inline-block; padding:2px 8px; border-radius:999px; ",
              "background:#1f78b4; color:white; font-weight:600; margin:0 2px;'>Top block support</span>, ",
              "<span style='display:inline-block; padding:2px 8px; border-radius:999px; ",
              "background:#6baed6; color:white; font-weight:600; margin:0 2px;'>Other block support (×0.05)</span>, ",
              "<span style='display:inline-block; padding:2px 8px; border-radius:999px; ",
              "background:#33a02c; color:white; font-weight:600; margin:0 2px;'>Gene-content bonus</span>, ",
              "and ",
              "<span style='display:inline-block; padding:2px 8px; border-radius:999px; ",
              "background:#6a3d9a; color:white; font-weight:600; margin:0 2px;'>Fragmentation bonus</span>. ",
              "The final cluster score is computed as the sum of these four components. ",
              "This keeps cluster prioritization aligned with the canonical gene and block scoring logic: ",
              "the strongest support unit dominates, while additional support contributes with reduced weight."
            )
          )
        ),
        
        plotly::plotlyOutput("prioritized_clusters_global_plot", height = "650px")
      )
    )
  })
  
  prioritized_clusters_global_plot_df <- reactive({
    df <- prioritized_cluster_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized clusters available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    needed_chr <- c(
      "cluster_id",
      "genes_in_blocks",
      "priority_class",
      "priority_class_relative"
    )
    
    needed_num <- c(
      "cluster_score",
      "cluster_score_from_blocks",
      "top_block_score",
      "other_block_score",
      "cluster_gene_bonus",
      "cluster_fragmentation_bonus",
      "top_gwas_hit_score",
      "n_gwas_hits",
      "n_blocks",
      "n_blocks_with_genes",
      "n_unique_genes_in_blocks",
      "n_apps_supported"
    )
    
    for (nm in needed_chr) {
      if (!nm %in% names(df)) df[[nm]] <- ""
    }
    
    for (nm in needed_num) {
      if (!nm %in% names(df)) df[[nm]] <- 0
    }
    
    df_ord <- df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        genes_in_blocks = dplyr::coalesce(as.character(genes_in_blocks), ""),
        chr_num = extract_chr_num(cluster_id),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$"))),
        
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0),
        cluster_score_from_blocks = dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0),
        top_block_score = dplyr::coalesce(as.numeric(top_block_score), 0),
        other_block_score = dplyr::coalesce(as.numeric(other_block_score), 0),
        cluster_gene_bonus = dplyr::coalesce(as.numeric(cluster_gene_bonus), 0),
        cluster_fragmentation_bonus = dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_blocks_with_genes = dplyr::coalesce(as.integer(n_blocks_with_genes), 0L),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        
        priority_class = dplyr::coalesce(as.character(priority_class), ""),
        priority_class_relative = dplyr::coalesce(as.character(priority_class_relative), "")
      ) %>%
      dplyr::mutate(
        n_genes_from_label = vapply(genes_in_blocks, count_semicolon_items, integer(1)),
        n_unique_genes_in_blocks = dplyr::if_else(
          n_unique_genes_in_blocks > 0L,
          n_unique_genes_in_blocks,
          n_genes_from_label
        ),
        weighted_other_block_score = 0.05 * other_block_score,
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999)
      ) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id)
    
    cluster_levels <- rev(unique(df_ord$cluster_id))
    
    df_ord %>%
      dplyr::transmute(
        cluster_id,
        cluster_score,
        top_block_score,
        weighted_other_block_score,
        cluster_gene_bonus,
        cluster_fragmentation_bonus,
        top_gwas_hit_score,
        n_gwas_hits,
        n_blocks,
        n_blocks_with_genes,
        n_unique_genes_in_blocks,
        n_apps_supported,
        priority_class,
        priority_class_relative
      ) %>%
      tidyr::pivot_longer(
        cols = c(
          top_block_score,
          weighted_other_block_score,
          cluster_gene_bonus,
          cluster_fragmentation_bonus
        ),
        names_to = "component",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        component = dplyr::recode(
          component,
          top_block_score = "Top block support",
          weighted_other_block_score = "Other block support (x0.05)",
          cluster_gene_bonus = "Gene-content bonus",
          cluster_fragmentation_bonus = "Fragmentation bonus"
        ),
        cluster_id = factor(cluster_id, levels = cluster_levels),
        cluster_id_chr = as.character(cluster_id),
        value = dplyr::coalesce(as.numeric(value), 0),
        hover_txt = paste0(
          "Cluster: ", cluster_id_chr,
          "<br>Component: ", component,
          "<br>Value: ", round(value, 2),
          "<br>Total cluster score: ", round(cluster_score, 2),
          "<br>Top GWAS hit score: ", round(top_gwas_hit_score, 2),
          "<br>ABS class: ", priority_class,
          "<br>REL class: ", priority_class_relative,
          "<br>GWAS hits: ", n_gwas_hits,
          "<br>Blocks: ", n_blocks,
          "<br>Blocks with genes: ", n_blocks_with_genes,
          "<br>Unique genes in blocks: ", n_unique_genes_in_blocks,
          "<br>Apps supported: ", n_apps_supported
        )
      ) %>%
      dplyr::select(cluster_id, cluster_id_chr, component, value, hover_txt)
  })
  
  output$prioritized_clusters_global_plot <- plotly::renderPlotly({
    df <- prioritized_clusters_global_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster priority component data available.")
    )
    
    plotly::plot_ly(
      data = df,
      y = ~cluster_id,
      x = ~value,
      color = ~component,
      colors = c(
        "Top block support" = "#1f78b4",
        "Other block support (x0.05)" = "#6baed6",
        "Gene-content bonus" = "#33a02c",
        "Fragmentation bonus" = "#6a3d9a"
      ),
      type = "bar",
      orientation = "h",
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis = list(
          title = "Cluster score components",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Cluster",
          categoryorder = "array",
          categoryarray = levels(df$cluster_id),
          automargin = TRUE
        ),
        legend = list(
          title = list(text = "Score component"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 90, r = 20, t = 20, b = 60)
      )
  })
  
  
  
  
  
  
  ########### troç modificat
  
  ld_global_log_buf <- reactiveVal(character(0))
  
  append_ld_global_log <- function(...) {
    txt <- paste(..., collapse = "")
    cur <- ld_global_log_buf()
    ld_global_log_buf(c(cur, txt))
  }
  
  output$ld_global_log <- renderText({
    paste(ld_global_log_buf(), collapse = "\n")
  })
  
  ld_summary_rds <- ld_summary_r
  
  ld_details_rds <- ld_details_r
  
  ld_block_hits_global <- reactive({
    dd <- ld_details_rds()
    
    if (!is.list(dd) || !length(dd)) {
      return(tibble::tibble())
    }
    
    pop_sel <- as.character(input$ld_global_pop %||% "")
    metric_sel <- as.character(input$ld_global_metric %||% "R2")
    
    out <- dplyr::bind_rows(lapply(dd, function(x) {
      if (is.null(x) || !is.list(x) || !is.data.frame(x$block_hits) || !nrow(x$block_hits)) {
        return(NULL)
      }
      
      bh <- tibble::as_tibble(x$block_hits)
      
      if (!"cluster_id" %in% names(bh)) bh$cluster_id <- as.character(x$cluster_id %||% NA_character_)
      if (!"chr" %in% names(bh)) bh$chr <- suppressWarnings(as.integer(x$chr %||% NA))
      if (!"population" %in% names(bh)) bh$population <- as.character(x$population %||% NA_character_)
      if (!"ld_metric" %in% names(bh)) bh$ld_metric <- as.character(x$ld_metric %||% NA_character_)
      
      bh
    }))
    
    if (!is.data.frame(out) || !nrow(out)) {
      return(tibble::tibble())
    }
    
    if (nzchar(pop_sel) && "population" %in% names(out)) {
      out <- out %>% dplyr::filter(as.character(population) == pop_sel)
    }
    
    if (nzchar(metric_sel) && "ld_metric" %in% names(out)) {
      out <- out %>% dplyr::filter(as.character(ld_metric) == metric_sel)
    }
    
    out
  })
  
  ### equival al block_overlap_summary_df del cluster
  block_overlap_summary_global_df <- reactive({
    block_overlap_summary_df(
      block_ranges = block_ranges_global(),
      block_hits = ld_block_hits_global(),
      block_genes = block_gene_overlap_summary_global(),
      gwas_bridge_df = gwas_bridge_r(),
      proxy_tbl = ld_proxy_table_global()
    )
  })
  
  ###################################### Taula amb score de blocks
  # pendent de confirmar si s'elimina
  output$block_overlap_summary_XXXXX <- DT::renderDT({
    
    df <- block_overlap_summary_global_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No global block summary available.")
    )
    
    show_df <- make_block_canonical_table_display(df)
    
    cat("[APP GLOBAL TABLE] show_df cols = ", paste(names(show_df), collapse = ", "), "\n", sep = "")
    
    show_df <- show_df %>%
      dplyr::mutate(
        gwas_min_p = dplyr::if_else(
          is.finite(gwas_min_p),
          formatC(gwas_min_p, format = "e", digits = 2),
          ""
        )
      )
    
    if ("genes_name" %in% names(show_df)) {
      genes_raw <- as.character(show_df$genes_name %||% rep("", nrow(show_df)))
      if (length(genes_raw) != nrow(show_df)) {
        genes_raw <- rep("", nrow(show_df))
      }
      
      show_df$genes_name_export <- genes_raw
      
      show_df$genes_name <- vapply(
        genes_raw,
        function(x) {
          x <- as.character(x %||% "")
          if (!nzchar(trimws(x))) return("")
          as.character(make_genecards_links(x))
        },
        character(1)
      )
    }
    
    if ("ld_proxy_snps" %in% names(show_df)) {
      ld_raw <- as.character(show_df$ld_proxy_snps %||% rep("", nrow(show_df)))
      if (length(ld_raw) != nrow(show_df)) {
        ld_raw <- rep("", nrow(show_df))
      }
      
      show_df$ld_proxy_snps_export <- ld_raw
      
      show_df$ld_proxy_snps <- vapply(
        ld_raw,
        function(x) {
          x <- as.character(x %||% "")
          if (!nzchar(trimws(x))) return("")
          build_snp_link_column_html(x)
        },
        character(1)
      )
    }
    
    if ("gwas_hits" %in% names(show_df)) {
      gwas_raw <- as.character(show_df$gwas_hits %||% rep("", nrow(show_df)))
      if (length(gwas_raw) != nrow(show_df)) {
        gwas_raw <- rep("", nrow(show_df))
      }
      
      show_df$gwas_hits_export <- gwas_raw
      
      show_df$gwas_hits <- vapply(
        gwas_raw,
        function(x) {
          x <- as.character(x %||% "")
          if (!nzchar(trimws(x))) return("")
          collapse_plain_list_html(
            x,
            max_visible = 5L,
            summary_label = paste0(length(split_semicolon(x)), " GWAS hits")
          )
        },
        character(1)
      )
    }
    
    cluster_col_idx <- which(names(show_df) == "cluster_id") - 1L
    
    js_export_block_table <- DT::JS("
function(data, row, column, node){
  var tbl = window.blockSummaryTable;
  if (!tbl) return $('<div>').html(data).text();
  
  var colName = tbl.column(column).header().textContent.trim();
  var rowData = tbl.row(row).data();
  if (!rowData) return $('<div>').html(data).text();
  
  if (colName === 'genes_name' && rowData.genes_name_export !== undefined && rowData.genes_name_export !== null) {
    return String(rowData.genes_name_export);
  }
  
  if (colName === 'ld_proxy_snps' && rowData.ld_proxy_snps_export !== undefined && rowData.ld_proxy_snps_export !== null) {
    return String(rowData.ld_proxy_snps_export);
  }
  
  if (colName === 'gwas_hits' && rowData.gwas_hits_export !== undefined && rowData.gwas_hits_export !== null) {
    return String(rowData.gwas_hits_export);
  }
  
  return $('<div>').html(data).text();
}
")
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      escape = FALSE,
      extensions = "Buttons",
      callback = DT::JS(sprintf("
      var tbl = table;
      window.blockSummaryTable = tbl;
      window.blockSummaryClusterCol = %d;
      
      tbl.columns.adjust();
      
      $(window).on('resize.blockSummaryTable', function() {
        tbl.columns.adjust();
      });
      
      $(document).on('shown.bs.tab.blockSummaryTable', 'button[data-bs-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function() {
        $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
      });
    ", cluster_col_idx)),
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = FALSE,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(
              modifier = list(page = "all"),
              format = list(body = js_export_block_table)
            )
          ),
          list(
            extend = "csv",
            exportOptions = list(
              modifier = list(page = "all"),
              format = list(body = js_export_block_table)
            )
          ),
          list(
            extend = "excel",
            exportOptions = list(
              modifier = list(page = "all"),
              format = list(body = js_export_block_table)
            )
          )
        ),
        columnDefs = c(
          if ("genes_name_export" %in% names(show_df)) {
            list(
              list(
                targets = which(names(show_df) == "genes_name_export") - 1L,
                visible = FALSE,
                searchable = FALSE
              )
            )
          } else {
            list()
          },
          if ("ld_proxy_snps_export" %in% names(show_df)) {
            list(
              list(
                targets = which(names(show_df) == "ld_proxy_snps_export") - 1L,
                visible = FALSE,
                searchable = FALSE
              )
            )
          } else {
            list()
          },
          if ("gwas_hits_export" %in% names(show_df)) {
            list(
              list(
                targets = which(names(show_df) == "gwas_hits_export") - 1L,
                visible = FALSE,
                searchable = FALSE
              )
            )
          } else {
            list()
          }
        ),
        order = list(list(which(names(show_df) == "block_mean_ld") - 1L, "desc"))
      )
    ) %>%
      DT::formatRound(
        c(
          "start", "end",
          "support_apps",
          "catalog_hits", "gtex_hits", "nonsyn_hits", "ewasdis_hits", "ewastum_hits",
          "n_ld_proxy_hits"
        ),
        0
      ) %>%
      DT::formatRound(
        c(
          "size_kb",
          "gwas_top_logp",
          "block_max_ld",
          "block_mean_ld"
        ),
        2
      )
  }, server = TRUE)
  
  ######################################
  
  output$ld_global_pop_ui <- renderUI({
    dd <- gi_pop_dir %||% ""
    d  <- tryCatch(normalizePath(dd, winslash = "/", mustWork = FALSE), error = function(e) dd)
    
    ff <- if (nzchar(d) && dir.exists(d)) list.files(d, pattern = "\\.[Tt][Xx][Tt]$", full.names = FALSE) else character(0)
    pops <- sub("\\.[Tt][Xx][Tt]$", "", ff)
    
    if (length(pops) == 0) {
      return(div(style = "color:#b00020; font-weight:600;", paste0("No .txt keep files found in: ", dd)))
    }
    
    sel <- if (!is.null(input$ld_global_pop) && input$ld_global_pop %in% pops) input$ld_global_pop else {
      if ("EUR" %in% pops) "EUR" else pops[1]
    }
    
    selectInput("ld_global_pop", "Population (global LD)", choices = sort(pops), selected = sel)
  })
  
  ld_global_ready <- reactive({
    dt <- ld_summary_rds()
    is.data.frame(dt) && nrow(dt) > 0
  })
  
  output$ld_notice_clusters <- renderUI({
    if (ld_global_ready()) return(NULL)
    
    div(
      class = "smallNote",
      style = "background:#fff3cd; border:1px solid #ffe69c; padding:10px; border-radius:8px; margin-bottom:12px;",
      HTML("LD global summary not computed yet. Go to <b>LD summary</b> and run <b>Compute global LD / blocks</b> to incorporate LD into the cluster prioritization score.")
    )
  })
  
  output$ld_notice_genes <- renderUI({
    if (ld_global_ready()) return(NULL)
    
    div(
      class = "smallNote",
      style = "background:#fff3cd; border:1px solid #ffe69c; padding:10px; border-radius:8px; margin-bottom:12px;",
      HTML("LD global summary not computed yet. Go to <b>LD summary</b> and run <b>Compute global LD / blocks</b> to incorporate LD into the gene prioritization score.")
    )
  })
  
  output$ld_global_summary_dt <- DT::renderDT({
    dt <- ld_summary_rds()
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No LD summary available yet. Run 'Compute global LD / blocks'."),
        rownames = FALSE
      ))
    }
    
    DT::datatable(
      sanitize_dt_types(dt),
      rownames = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        scrollX = TRUE,
        pageLength = 15,
        autoWidth = FALSE
      )
    )
  }, server = TRUE)
  
  output$ld_summary_dt <- DT::renderDT({
    dt <- ld_summary_rds()
    
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No LD summary available yet. Run 'Compute global LD / blocks'."),
        rownames = FALSE
      ))
    }
    
    # ------------------------------------------------------------
    # Helper: converteix cel·les amb llistes llargues en desplegables
    # i opcionalment aplica funcions de links als ítems
    # ------------------------------------------------------------
    make_collapsible_items <- function(x, max_items = 5, link_fun = NULL) {
      
      split_items <- function(txt) {
        if (is.null(txt) || is.na(txt) || !nzchar(trimws(as.character(txt)))) {
          return(character(0))
        }
        
        txt <- as.character(txt)
        parts <- unlist(strsplit(txt, "\\s*;\\s*|\\s*\\|\\s*|\\n+"))
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]
        unique(parts)
      }
      
      vapply(x, function(cell) {
        items <- split_items(cell)
        
        if (length(items) == 0) {
          return("")
        }
        
        items_out <- items
        if (is.function(link_fun)) {
          items_out <- link_fun(items)
        }
        
        if (length(items_out) <= max_items) {
          return(paste(items_out, collapse = "<br>"))
        }
        
        paste0(
          "<div class='ld-collapsible-cell'>",
          "<a href='#' class='ld-toggle'>", length(items_out), " items ▼</a>",
          "<div class='ld-collapsible-content' style='display:none; margin-top:4px;'>",
          paste(items_out, collapse = "<br>"),
          "</div>",
          "</div>"
        )
      }, character(1))
    }
    
    clean_gene_items_text <- function(x) {
      if (is.null(x) || is.na(x) || !nzchar(trimws(as.character(x)))) {
        return("")
      }
      
      txt <- as.character(x)
      parts <- unlist(strsplit(txt, "\\s*;\\s*|\\s*\\|\\s*|\\n+"))
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      
      # elimina identificadors purament numèrics
      parts <- parts[!grepl("^[0-9]+$", parts)]
      
      parts <- unique(parts)
      paste(parts, collapse = "; ")
    }
    
    dt2 <- sanitize_dt_types(dt)
    
    
    # ------------------------------------------------------------
    # Columnes amb links
    # ------------------------------------------------------------
    gene_cols <- intersect(
      c(
        "block_genes",
        "genes_in_block",
        "genes",
        "gwas_genes",
        "catalog_genes",
        "gtex_genes",
        "nonsyn_genes",
        "ewasdis_genes",
        "ewastum_genes"
      ),
      names(dt2)
    )
    
    snp_cols <- intersect(
      c(
        "gwas_hits",
        "seed_snps",
        "proxy_snps",
        "lead_snps",
        "supporting_snps"
      ),
      names(dt2)
    )
    
    trait_cols <- intersect(
      c("catalog_traits"),
      names(dt2)
    )
    
    plain_collapsible_cols <- intersect(
      c(),
      names(dt2)
    )
    
    # Gens -> GeneCards
    if (length(gene_cols)) {
      for (cc in gene_cols) {
        dt2[[cc]] <- vapply(dt2[[cc]], clean_gene_items_text, character(1))
        
        dt2[[cc]] <- make_collapsible_items(
          dt2[[cc]],
          max_items = 5,
          link_fun = make_genecards_links
        )
      }
    }
    
    # SNPs / rsid -> dbSNP
    if (length(snp_cols)) {
      for (cc in snp_cols) {
        dt2[[cc]] <- make_collapsible_items(
          dt2[[cc]],
          max_items = 5,
          link_fun = make_dbsnp_links
        )
      }
    }
    
    # Traits -> GWAS Catalog
    if (length(trait_cols)) {
      for (cc in trait_cols) {
        dt2[[cc]] <- make_collapsible_items(
          dt2[[cc]],
          max_items = 5,
          link_fun = make_gwascatalog_term_links
        )
      }
    }
    
    # Altres columnes llargues sense links
    if (length(plain_collapsible_cols)) {
      for (cc in plain_collapsible_cols) {
        dt2[[cc]] <- make_collapsible_items(
          dt2[[cc]],
          max_items = 5,
          link_fun = NULL
        )
      }
    }
    
    DT::datatable(
      dt2,
      rownames = FALSE,
      escape = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(modifier = list(page = "all"))
          ),
          list(
            extend = "csv",
            exportOptions = list(modifier = list(page = "all"))
          ),
          list(
            extend = "excel",
            exportOptions = list(modifier = list(page = "all"))
          )
        ),
        scrollX = TRUE,
        pageLength = 15,
        autoWidth = FALSE
      ),
      callback = DT::JS("
      table.on('click', 'a.ld-toggle', function(e) {
        e.preventDefault();
        var $link = $(this);
        var $content = $link.closest('.ld-collapsible-cell').find('.ld-collapsible-content');
        
        if ($content.is(':visible')) {
          $content.hide();
          $link.text($link.text().replace('▲', '▼'));
        } else {
          $content.show();
          $link.text($link.text().replace('▼', '▲'));
        }
      });
    ")
    ) %>%
      DT::formatRound("mean_ld_value", 2)
  })
  
  block_gene_overlap_summary_global <- reactive({
    dd <- ld_details_rds()
    
    empty_out <- tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_genes_in_block = integer(),
      genes_in_block = character()
    )
    
    if (!is.list(dd) || !length(dd)) {
      return(empty_out)
    }
    
    clean_gene_vector <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      
      # elimina identificadors purament numèrics tipus ENTREZID
      x <- x[!grepl("^[0-9]+$", x)]
      
      sort(unique(x))
    }
    
    clean_gene_string <- function(x) {
      x <- as.character(x %||% "")
      x <- trimws(x)
      if (!nzchar(x)) {
        return(list(
          genes_in_block = "",
          n_genes_in_block = 0L
        ))
      }
      
      parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
      parts <- clean_gene_vector(parts)
      
      list(
        genes_in_block = paste(parts, collapse = "; "),
        n_genes_in_block = length(parts)
      )
    }
    
    out_list <- lapply(dd, function(x) {
      if (is.null(x) || !is.list(x) || !is.data.frame(x$block_genes) || !nrow(x$block_genes)) {
        return(NULL)
      }
      
      bg <- tibble::as_tibble(x$block_genes)
      
      if (!"cluster_id" %in% names(bg)) {
        bg$cluster_id <- as.character(x$cluster_id %||% NA_character_)
      }
      
      if (!"block_id" %in% names(bg)) {
        return(NULL)
      }
      
      # --------------------------------------------------
      # Cas 1: ja ve resumit
      # --------------------------------------------------
      if (all(c("n_genes_in_block", "genes_in_block") %in% names(bg))) {
        out1 <- bg %>%
          dplyr::transmute(
            cluster_id = as.character(cluster_id),
            block_id = as.character(block_id),
            genes_in_block_raw = dplyr::coalesce(as.character(genes_in_block), "")
          ) %>%
          dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
        
        cleaned <- lapply(out1$genes_in_block_raw, clean_gene_string)
        
        out1$genes_in_block <- vapply(cleaned, `[[`, character(1), "genes_in_block")
        out1$n_genes_in_block <- vapply(cleaned, `[[`, integer(1), "n_genes_in_block")
        
        return(
          out1 %>%
            dplyr::select(cluster_id, block_id, n_genes_in_block, genes_in_block)
        )
      }
      
      # --------------------------------------------------
      # Cas 2: ve en format llarg amb una columna gene
      # --------------------------------------------------
      gene_col <- intersect(c("gene", "gene_symbol", "symbol", "label"), names(bg))
      gene_col <- if (length(gene_col)) gene_col[1] else NA_character_
      
      if (is.na(gene_col)) return(NULL)
      
      bg %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          gene = trimws(as.character(.data[[gene_col]]))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(block_id), nzchar(block_id),
          !is.na(gene), nzchar(gene)
        ) %>%
        dplyr::group_by(cluster_id, block_id) %>%
        dplyr::summarise(
          genes_in_block = paste(clean_gene_vector(gene), collapse = "; "),
          n_genes_in_block = length(clean_gene_vector(gene)),
          .groups = "drop"
        )
    })
    
    out <- dplyr::bind_rows(out_list)
    
    if (!is.data.frame(out) || !nrow(out)) {
      return(empty_out)
    }
    
    out %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        n_genes_in_block = safe_int0(n_genes_in_block),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
      ) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
  })
  
  output$ld_block_genes_global_dt <- DT::renderDT({
    xx <- block_gene_overlap_summary_global()
    
    if (!is.data.frame(xx) || !nrow(xx)) {
      return(DT::datatable(
        data.frame(Message = "No block-level genes available."),
        rownames = FALSE
      ))
    }
    
    DT::datatable(
      sanitize_dt_types(xx),
      rownames = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        scrollX = TRUE,
        pageLength = 15,
        autoWidth = FALSE
      )
    )
  }, server = TRUE)
  
  ########################################
  
  ld_proxy_table_global <- reactive({
    x <- ld_details_rds()
    
    empty <- tibble::tibble(
      cluster_id = character(),
      query_hit = character(),
      query_pos = integer(),
      seed_snp = character(),
      seed_pos = integer(),
      seed_type = character(),
      seed_dist_bp = integer(),
      proxy_snp = character(),
      proxy_pos = integer(),
      ld_value = numeric()
    )
    
    if (!is.list(x) || !length(x)) {
      return(empty)
    }
    
    out <- lapply(names(x), function(cid) {
      obj <- x[[cid]]
      if (!is.list(obj)) return(NULL)
      
      px <- obj$proxies
      if (!is.data.frame(px) || !nrow(px)) return(NULL)
      
      px %>%
        dplyr::mutate(
          cluster_id   = as.character(cluster_id %||% cid),
          query_hit    = as.character(query_hit),
          query_pos    = suppressWarnings(as.integer(query_pos)),
          seed_snp     = as.character(seed_snp),
          seed_pos     = suppressWarnings(as.integer(seed_pos)),
          seed_type    = as.character(seed_type),
          seed_dist_bp = suppressWarnings(as.integer(seed_dist_bp)),
          proxy_snp    = as.character(proxy_snp),
          proxy_pos    = suppressWarnings(as.integer(proxy_pos)),
          ld_value     = suppressWarnings(as.numeric(ld_value))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(proxy_snp), nzchar(proxy_snp),
          is.finite(ld_value)
        ) %>%
        dplyr::distinct()
    })
    
    dplyr::bind_rows(out)
  })
  
  output$ld_proxy_global_dt <- DT::renderDT({
    px <- ld_proxy_table_global()
    
    if (!is.data.frame(px) || !nrow(px)) {
      return(DT::datatable(
        data.frame(Message = "No proxy table available in LD details."),
        rownames = FALSE
      ))
    }
    
    DT::datatable(
      sanitize_dt_types(px),
      rownames = FALSE,
      extensions = "Buttons",
      width = "100%",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        scrollX = TRUE,
        pageLength = 15,
        autoWidth = FALSE
      )
    )
  }, server = TRUE)
  
  #########################################################################
  observeEvent(plotly::event_data("plotly_click", source = "ld_global_blocks"), {
    ed <- plotly::event_data("plotly_click", source = "ld_global_blocks")
    
    if (is.null(ed) || !nrow(ed)) return()
    
    cid <- as.character(ed$customdata[[1]] %||% "")
    if (!nzchar(cid)) return()
    
    session$sendCustomMessage(
      "filter_block_table_by_cluster",
      list(cluster_id = cid)
    )
  }, ignoreInit = TRUE)
  
  # ============================================================
  # ALL BLOCKS + PRIORITIZED SCORE OVERLAY
  # ============================================================
  global_ld_blocks_plot_df <- reactive({
    # base: TOTS els blocs
    df_all <- block_overlap_summary_global_df()
    
    validate(
      need(is.data.frame(df_all) && nrow(df_all) > 0, "No global block summary available.")
    )
    
    # scores: només blocs prioritzats
    df_score <- prioritized_block_df_v2()
    
    if (!is.data.frame(df_score)) {
      df_score <- tibble::tibble()
    }
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    size_kb_col <- NULL
    if ("block_size_kb" %in% names(df_all)) {
      size_kb_col <- "block_size_kb"
    } else if ("size_kb" %in% names(df_all)) {
      size_kb_col <- "size_kb"
    }
    
    block_label_col <- NULL
    if ("block_label" %in% names(df_all)) {
      block_label_col <- "block_label"
    }
    
    genes_col <- NULL
    if ("genes_in_block" %in% names(df_all)) {
      genes_col <- "genes_in_block"
    } else if ("genes_name" %in% names(df_all)) {
      genes_col <- "genes_name"
    }
    
    ngenes_col <- NULL
    if ("n_genes_in_block" %in% names(df_all)) {
      ngenes_col <- "n_genes_in_block"
    } else if ("n_genes" %in% names(df_all)) {
      ngenes_col <- "n_genes"
    }
    
    base_df <- df_all %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        chr_num = extract_chr_num(cluster_id),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$"))),
        block_label = if (!is.null(block_label_col)) as.character(.data[[block_label_col]]) else NA_character_,
        block_ord = extract_block_order(block_id, block_label),
        block_label2 = paste0("B", block_ord),
        block_size_kb = if (!is.null(size_kb_col)) suppressWarnings(as.numeric(.data[[size_kb_col]])) else 0,
        genes_in_block = if (!is.null(genes_col)) dplyr::coalesce(as.character(.data[[genes_col]]), "") else "",
        n_genes_in_block = if (!is.null(ngenes_col)) dplyr::coalesce(as.integer(.data[[ngenes_col]]), 0L) else 0L
      ) %>%
      dplyr::mutate(
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999),
        block_ord = dplyr::coalesce(block_ord, 999L),
        block_size_kb = dplyr::coalesce(block_size_kb, 0)
      ) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id, block_ord, block_id)
    
    score_df <- df_score %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score_v2 = suppressWarnings(as.numeric(block_score)),
        n_gwas_hits_v2 = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        top_gwas_hit_score_v2 = dplyr::coalesce(suppressWarnings(as.numeric(top_gwas_hit_score)), 0),
        n_apps_supported_v2 = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        priority_class_v2 = dplyr::coalesce(as.character(priority_class), "")
      ) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
    
    base_df %>%
      dplyr::left_join(score_df, by = c("cluster_id", "block_id")) %>%
      dplyr::mutate(
        block_score = block_score_v2,
        n_gwas_hits = dplyr::coalesce(n_gwas_hits_v2, 0L),
        top_gwas_hit_score = dplyr::coalesce(top_gwas_hit_score_v2, 0),
        n_apps_supported = dplyr::coalesce(n_apps_supported_v2, 0L),
        priority_class = dplyr::coalesce(priority_class_v2, ""),
        block_score = dplyr::coalesce(block_score, NA_real_)
      ) %>%
      dplyr::select(-dplyr::ends_with("_v2"))
  })
  
  # ============================================================
  # SCORE LEGEND
  # ============================================================
  global_ld_blocks_score_legend_ui <- function(df) {
    if (!is.data.frame(df) || !nrow(df) || !"block_score" %in% names(df)) {
      return(NULL)
    }
    
    vals <- suppressWarnings(as.numeric(df$block_score))
    vals <- vals[is.finite(vals)]
    
    if (!length(vals)) {
      return(
        tags$div(
          style = "margin-top:8px; font-size:12px; color:#444;",
          tags$span(style = "font-weight:600; margin-right:10px;", "Block priority score"),
          tags$span(style = "color:#777;", "No score values available")
        )
      )
    }
    
    vmin <- min(vals, na.rm = TRUE)
    vmax <- max(vals, na.rm = TRUE)
    
    tags$div(
      style = "margin-top:8px; font-size:12px; color:#444;",
      tags$span(style = "font-weight:600; margin-right:10px;", "Block priority score"),
      tags$span(
        style = paste0(
          "display:inline-block; width:220px; height:12px; vertical-align:middle; ",
          "border:1px solid #999; border-radius:6px; margin-right:8px; ",
          "background: linear-gradient(to right, ",
          "#d9d9d9 0%, ",
          "#fff7bc 20%, ",
          "#fec44f 40%, ",
          "#fe9929 60%, ",
          "#d95f0e 80%, ",
          "#993404 100%);"
        )
      ),
      tags$span(formatC(vmin, format = "f", digits = 2)),
      tags$span(" – "),
      tags$span(formatC(vmax, format = "f", digits = 2))
    )
  }
  
  output$global_ld_blocks_score_legend <- renderUI({
    df <- global_ld_blocks_plot_df()
    global_ld_blocks_score_legend_ui(df)
  })
  
  # ============================================================
  # BARPLOT: ALL BLOCKS, COLOR BY PRIORITIZED SCORE
  # ============================================================
  output$global_ld_blocks_plot <- plotly::renderPlotly({
    df <- global_ld_blocks_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No global block summary available.")
    )
    
    cluster_levels <- df %>%
      dplyr::distinct(cluster_id, chr_num, cluster_num) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
      dplyr::pull(cluster_id)
    
    cluster_levels <- rev(cluster_levels)
    
    vals <- suppressWarnings(as.numeric(df$block_score))
    has_score <- any(is.finite(vals))
    
    if (has_score) {
      rng <- range(vals, na.rm = TRUE)
      
      if (!all(is.finite(rng))) {
        df$fill_col <- ifelse(is.finite(vals), "#fff7bc", "#d9d9d9")
      } else if (rng[1] == rng[2]) {
        df$fill_col <- ifelse(is.finite(vals), "#fe9929", "#d9d9d9")
      } else {
        pal_fun <- scales::col_numeric(
          palette = c("#d9d9d9","#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
          domain = rng,
          na.color = "#d9d9d9"
        )
        df$fill_col <- pal_fun(vals)
      }
    } else {
      df$fill_col <- "#d9d9d9"
    }
    
    df <- df %>%
      dplyr::mutate(
        cluster_id_chr = as.character(cluster_id),
        cluster_id = factor(cluster_id, levels = cluster_levels),
        score_lab = dplyr::if_else(
          is.finite(block_score),
          formatC(block_score, format = "f", digits = 2),
          "NA"
        ),
        hover_txt = paste0(
          "Cluster: ", cluster_id_chr,
          "<br>Block: ", block_id,
          "<br>Order: ", block_label2,
          "<br>Size (kb): ", round(block_size_kb, 2),
          "<br>Priority score: ", score_lab,
          ifelse(
            nzchar(priority_class),
            paste0("<br>Priority class: ", priority_class),
            ""
          ),
          "<br>GWAS hits: ", n_gwas_hits,
          "<br>Top GWAS-hit score: ", round(top_gwas_hit_score, 2),
          "<br>Apps supported: ", n_apps_supported,
          ifelse(
            nzchar(genes_in_block),
            paste0("<br>Genes: ", genes_in_block),
            ""
          )
        )
      )
    
    ords <- sort(unique(df$block_ord))
    
    p <- plotly::plot_ly(source = "ld_global_blocks")
    
    for (oo in ords) {
      dsub <- df %>% dplyr::filter(block_ord == oo)
      if (!nrow(dsub)) next
      
      p <- p %>%
        plotly::add_trace(
          data = dsub,
          y = ~cluster_id,
          x = ~block_size_kb,
          type = "bar",
          orientation = "h",
          customdata = ~cluster_id_chr,
          text = NULL,
          textposition = "none",
          hovertext = ~hover_txt,
          hoverinfo = "text",
          showlegend = FALSE,
          marker = list(
            color = dsub$fill_col,
            line = list(color = "black", width = 0.5)
          )
        )
    }
    
    p %>%
      plotly::layout(
        barmode = "stack",
        clickmode = "event+select",
        xaxis = list(
          title = "Stacked block size (kb)",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Cluster",
          categoryorder = "array",
          categoryarray = cluster_levels,
          automargin = TRUE
        ),
        margin = list(l = 90, r = 20, t = 20, b = 60)
      )
  })
  
# observeEvent(input$clear_block_table_filter, {
#   session$sendCustomMessage(
#     "filter_block_table_by_cluster",
#     list(cluster_id = "")
#   )
# })
  
  ############################################################## 
  # ============================================================
  # Canonical prioritization sources
  # ============================================================
  
  gwas_hit_priority_mod <- gwas_hit_priority_module_server(
    id = "gwas_hit_priority",
    gene_bridge_lookup_df = gene_bridge_lookup_df,
    nonsyn_gene_hit_bridge_df = nonsyn_gene_hit_bridge_df,
    integrated_candidates_r = integrated_candidates_r,
    block_hits_global_r = reactive({
      ld_block_hits_global()
    }),
    clusters_consensus_r = integrated_clusters_r,
    canonical_prioritized_cluster_genes_r = canonical_prioritized_cluster_genes_df,
    gene_gwas_hit_score_audit_r = gene_gwas_hit_score_audit_df,
    gwas_hit_priority_df_v2_r = gwas_hit_priority_df_v2
  )
  
  # ---------------------------
  # Font canónica genes
  # ---------------------------
  
  priority_source_gwas_hits <- reactive({
    req(gwas_hit_priority_mod)
    
    df <- tryCatch(
      gwas_hit_priority_mod$gwas_hit_match_audit_summary_with_gene_df(),
      error = function(e) tibble::tibble()
    )
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0,
           "No GWAS hit × app summary with gene available.")
    )
    
    df
  })
  
  # ###
  prioritized_gwas_hits_df <- reactive({
    df <- gwas_hit_priority_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0,
           "No prioritized GWAS-hit table available.")
    )
    
    df
  })
  
  # ---------------------------
  # Font canónica blocks
  # ---------------------------
  
  
  ld_mod <- ld_integrator_module_server(
    id = "ld",
    activate_r = reactive(TRUE),
    integrator_dir_r = selected_session_dir,
    integrated_clusters_r = integrated_clusters_r,
    integrated_candidates_r = integrated_candidates_r,
    gwas_bridge_r = gwas_bridge_r
  )
  
  # ---------------------------
  # Font canónica clusters
  # ---------------------------
  
  priority_source_clusters <- reactive({
    df <- tryCatch(
      integrated_clusters_from_rds(selected_session_dir()),
      error = function(e) tibble::tibble()
    )
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0,
           "No integrated clusters available.")
    )
    
    df
  })
  
  ###############################################################################
  # Canonical GWAS-hit priority table
  ###############################################################################
  
  physical_prioritized_cluster_genes_df <- reactive({
    cl <- priority_source_clusters()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available.")
    )
    
    cl_std <- cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        chr = as.character(chr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(chr), nzchar(chr),
        is.finite(start),
        is.finite(end)
      ) %>%
      dplyr::mutate(
        chr = sub("^chr", "", chr, ignore.case = TRUE),
        chr = dplyr::case_when(
          chr == "23" ~ "X",
          chr == "24" ~ "Y",
          chr %in% c("25", "26") ~ "MT",
          TRUE ~ chr
        ),
        seqnames = paste0("chr", chr)
      )
    
    gene_ranges <- genes(txdb)
    
    gene_df <- tibble::tibble(
      entrezid = names(gene_ranges),
      seqnames = as.character(GenomicRanges::seqnames(gene_ranges)),
      gene_start = BiocGenerics::start(gene_ranges),
      gene_end = BiocGenerics::end(gene_ranges)
    ) %>%
      dplyr::mutate(
        gene = AnnotationDbi::mapIds(
          org.Hs.eg.db,
          keys = entrezid,
          keytype = "ENTREZID",
          column = "SYMBOL",
          multiVals = "first"
        ),
        gene = dplyr::coalesce(as.character(gene), as.character(entrezid)),
        gene = trimws(gene)
      ) %>%
      dplyr::filter(
        !is.na(gene),
        nzchar(gene),
        tolower(gene) != "numeric",
        !grepl("^[0-9]+$", gene)
      )
    
    cl_std %>%
      dplyr::inner_join(gene_df, by = "seqnames", relationship = "many-to-many") %>%
      dplyr::filter(gene_end >= start, gene_start <= end) %>%
      dplyr::transmute(
        cluster_id,
        chr,
        gene,
        gene_start,
        gene_end
      ) %>%
      dplyr::distinct(cluster_id, gene, .keep_all = TRUE)
  })
  
  canonical_prioritized_cluster_genes_df <- reactive({
    cl <- priority_source_clusters()
    gb <- tryCatch(gene_bridges_combined(), error = function(e) tibble::tibble())
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available.")
    )
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(sort(unique(x)), collapse = "; ")
    }
    
    # ------------------------------------------------------------
    # 1) physical genes from TxDb overlap
    # ------------------------------------------------------------
    cl_std <- cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        chr = as.character(chr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(chr), nzchar(chr),
        is.finite(start),
        is.finite(end)
      ) %>%
      dplyr::mutate(
        chr = sub("^chr", "", chr, ignore.case = TRUE),
        chr = dplyr::case_when(
          chr == "23" ~ "X",
          chr == "24" ~ "Y",
          chr %in% c("25", "26") ~ "MT",
          TRUE ~ chr
        ),
        seqnames = paste0("chr", chr)
      )
    
    gene_ranges <- genes(txdb)
    
    gene_df <- tibble::tibble(
      entrezid = names(gene_ranges),
      seqnames = as.character(GenomicRanges::seqnames(gene_ranges)),
      gene_start = BiocGenerics::start(gene_ranges),
      gene_end = BiocGenerics::end(gene_ranges)
    ) %>%
      dplyr::mutate(
        gene = AnnotationDbi::mapIds(
          org.Hs.eg.db,
          keys = entrezid,
          keytype = "ENTREZID",
          column = "SYMBOL",
          multiVals = "first"
        ),
        gene = dplyr::coalesce(as.character(gene), as.character(entrezid)),
        gene = trimws(gene)
      ) %>%
      dplyr::filter(
        !is.na(gene),
        nzchar(gene),
        tolower(gene) != "numeric",
        !grepl("^[0-9]+$", gene)
      )
    
    physical_df <- cl_std %>%
      dplyr::inner_join(gene_df, by = "seqnames", relationship = "many-to-many") %>%
      dplyr::filter(gene_end >= start, gene_start <= end) %>%
      dplyr::transmute(
        cluster_id,
        chr,
        gene,
        gene_link_mode = "physical",
        gene_evidence_types = "physical_overlap",
        gene_source_apps = ""
      ) %>%
      dplyr::distinct(cluster_id, gene, .keep_all = TRUE)
    
    # ------------------------------------------------------------
    # 2) effect / evidence-linked genes from bridges
    # ------------------------------------------------------------
    effect_df <- if (is.data.frame(gb) && nrow(gb)) {
      gb %>%
        dplyr::mutate(
          gene = trimws(as.character(gene)),
          cluster_id = trimws(as.character(cluster_id)),
          chr = as.character(chr),
          source_app = tolower(trimws(as.character(source_app))),
          evidence_type = trimws(as.character(evidence_type))
        ) %>%
        dplyr::filter(
          !is.na(gene), nzchar(gene),
          tolower(gene) != "numeric",
          !grepl("^[0-9]+$", gene),
          !is.na(cluster_id), nzchar(cluster_id),
          evidence_type %in% c("catalog_gene", "nonsyn_gene", "eqtl_gene", "ewas_nearby_gene")
        ) %>%
        dplyr::semi_join(
          cl %>% dplyr::transmute(cluster_id = as.character(cluster_id)) %>% dplyr::distinct(),
          by = "cluster_id"
        ) %>%
        dplyr::group_by(cluster_id, gene) %>%
        dplyr::summarise(
          chr = dplyr::first(as.character(chr)),
          gene_link_mode = "effect",
          gene_evidence_types = collapse_unique_semicolon(evidence_type),
          gene_source_apps = collapse_unique_semicolon(source_app),
          .groups = "drop"
        )
    } else {
      tibble::tibble(
        cluster_id = character(),
        chr = character(),
        gene = character(),
        gene_link_mode = character(),
        gene_evidence_types = character(),
        gene_source_apps = character()
      )
    }
    
    # ------------------------------------------------------------
    # 3) canonical union
    # ------------------------------------------------------------
    dplyr::bind_rows(physical_df, effect_df) %>%
      dplyr::group_by(cluster_id, gene) %>%
      dplyr::summarise(
        chr = dplyr::first(chr[!is.na(chr) & nzchar(chr)]),
        gene_link_mode = collapse_unique_semicolon(gene_link_mode),
        gene_evidence_types = collapse_unique_semicolon(gene_evidence_types),
        gene_source_apps = collapse_unique_semicolon(gene_source_apps),
        .groups = "drop"
      ) %>%
      dplyr::arrange(cluster_id, gene)
  })
  
  ###############################################################################
  # ---------------------------
  # Prioritized genes v2
  # ---------------------------
  
  prioritized_gene_df_v2 <- reactive({
    cl_base <- canonical_prioritized_cluster_genes_df() %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        gene = trimws(as.character(gene))
      ) %>%
      tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s*,\\s*") %>%
      dplyr::mutate(
        gene = trimws(gene)
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::distinct()
    
    audit_units <- gene_gwas_hit_score_audit_df()
    
    validate(
      need(is.data.frame(cl_base) && nrow(cl_base) > 0, "No canonical genes available for prioritized clusters.")
    )
    
    split_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return(character(0))
      
      parts <- unlist(strsplit(x, "\\s*;\\s*"))
      parts <- trimws(parts)
      parts <- parts[!is.na(parts) & nzchar(parts)]
      unique(parts)
    }
    
    collapse_unique_semicolon <- function(x) {
      vals <- split_unique_semicolon(x)
      if (!length(vals)) return("")
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    count_unique_semicolon <- function(x) {
      vals <- split_unique_semicolon(x)
      length(vals)
    }
    
    count_blocks_by_app_string <- function(block_keys, app_strings) {
      block_keys <- as.character(block_keys)
      app_strings <- as.character(app_strings)
      
      n <- min(length(block_keys), length(app_strings))
      if (!n) return("")
      
      df_tmp <- tibble::tibble(
        pair_id = seq_len(n),
        block_key = block_keys[seq_len(n)],
        app_str = app_strings[seq_len(n)]
      ) %>%
        dplyr::mutate(
          block_key = dplyr::coalesce(.data$block_key, ""),
          app_str = dplyr::coalesce(.data$app_str, "")
        ) %>%
        dplyr::filter(nzchar(.data$block_key), nzchar(.data$app_str)) %>%
        tidyr::separate_rows(app_str, sep = "\\s*;\\s*") %>%
        dplyr::mutate(
          app_str = trimws(.data$app_str)
        ) %>%
        dplyr::filter(nzchar(.data$app_str)) %>%
        dplyr::distinct(.data$block_key, .data$app_str)
      
      if (!nrow(df_tmp)) return("")
      
      df_tmp %>%
        dplyr::count(.data$app_str, name = "n_blocks") %>%
        dplyr::arrange(.data$app_str) %>%
        dplyr::transmute(txt = paste0(.data$app_str, " (", .data$n_blocks, ")")) %>%
        dplyr::pull(.data$txt) %>%
        paste(collapse = "; ")
    }
    
    gene_support <- if (is.data.frame(audit_units) && nrow(audit_units)) {
      audit_units %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          gene = trimws(as.character(gene)),
          gwas_hit = trimws(as.character(gwas_hit)),
          gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
          block_id = dplyr::coalesce(as.character(block_id), ""),
          source_app = dplyr::coalesce(as.character(source_app), ""),
          score_app = dplyr::coalesce(as.character(score_app), ""),
          apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
          match_apps = dplyr::coalesce(as.character(match_apps), ""),
          marker_apps = dplyr::coalesce(as.character(marker_apps), ""),
          gene_score_component = dplyr::coalesce(as.numeric(gene_score_component), 0),
          gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0),
          hit_support_score = gene_score_component * gwas_hit_priority_score,
          supporting_block_key = dplyr::case_when(
            nzchar(block_id) ~ paste(cluster_id, block_id, sep = "||"),
            TRUE ~ ""
          )
        ) %>%
        dplyr::group_by(cluster_id, gene) %>%
        dplyr::arrange(
          dplyr::desc(hit_support_score),
          dplyr::desc(gwas_hit_priority_score),
          .by_group = TRUE
        ) %>%
        dplyr::summarise(
          top_hit_support_score = {
            hs <- hit_support_score
            hs <- hs[!is.na(hs)]
            if (!length(hs)) 0 else hs[1]
          },
          other_hit_support_score = {
            hs <- hit_support_score
            hs <- hs[!is.na(hs)]
            if (length(hs) <= 1) 0 else sum(hs[-1], na.rm = TRUE)
          },
          other_hits_weight = 0.05,
          gene_score = top_hit_support_score + other_hits_weight * other_hit_support_score,
          n_supporting_blocks = dplyr::n_distinct(supporting_block_key[nzchar(supporting_block_key)]),
          top_gwas_hit_score = max(gwas_hit_priority_score, na.rm = TRUE),
          apps_supported_audit = collapse_unique_semicolon(apps_supported),
          match_support = count_blocks_by_app_string(
            supporting_block_key[nzchar(supporting_block_key)],
            match_apps[nzchar(supporting_block_key)]
          ),
          marker_support = count_blocks_by_app_string(
            supporting_block_key[nzchar(supporting_block_key)],
            marker_apps[nzchar(supporting_block_key)]
          ),
          .groups = "drop"
        )
    } else {
      tibble::tibble(
        cluster_id = character(),
        gene = character(),
        top_hit_support_score = numeric(),
        other_hit_support_score = numeric(),
        other_hits_weight = numeric(),
        gene_score = numeric(),
        n_supporting_blocks = integer(),
        top_gwas_hit_score = numeric(),
        apps_supported_audit = character(),
        match_support = character(),
        marker_support = character()
      )
    }
    
    cl_base %>%
      dplyr::left_join(
        gene_support,
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::mutate(
        top_hit_support_score = dplyr::coalesce(as.numeric(top_hit_support_score), 0),
        other_hit_support_score = dplyr::coalesce(as.numeric(other_hit_support_score), 0),
        other_hits_weight = dplyr::coalesce(as.numeric(other_hits_weight), 0.05),
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0),
        n_supporting_blocks = dplyr::coalesce(as.integer(n_supporting_blocks), 0L),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        apps_supported_audit = dplyr::coalesce(as.character(apps_supported_audit), ""),
        match_support = dplyr::coalesce(as.character(match_support), ""),
        marker_support = dplyr::coalesce(as.character(marker_support), "")
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        cluster_ids = paste(sort(unique(cluster_id)), collapse = "; "),
        top_hit_support_score = round(sum(top_hit_support_score, na.rm = TRUE), 4),
        other_hit_support_score = round(sum(other_hit_support_score, na.rm = TRUE), 4),
        other_hits_weight = max(other_hits_weight, na.rm = TRUE),
        gene_score = round(sum(gene_score, na.rm = TRUE), 2),
        n_supporting_blocks = sum(n_supporting_blocks, na.rm = TRUE),
        n_clusters = dplyr::n_distinct(cluster_id),
        top_gwas_hit_score = round(max(top_gwas_hit_score, na.rm = TRUE), 2),
        gene_link_mode = collapse_unique_semicolon(gene_link_mode),
        gene_evidence_types = collapse_unique_semicolon(gene_evidence_types),
        gene_source_apps = collapse_unique_semicolon(gene_source_apps),
        apps_supported_audit = collapse_unique_semicolon(apps_supported_audit),
        match_support = collapse_unique_semicolon(match_support),
        marker_support = collapse_unique_semicolon(marker_support),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        apps_supported = dplyr::coalesce(as.character(apps_supported_audit), ""),
        n_apps_supported = vapply(apps_supported, count_unique_semicolon, integer(1)),
        priority_class = dplyr::case_when(
          gene_score >= 40 ~ "High",
          gene_score >= 10 ~ "Medium",
          TRUE ~ "Low"
        ),
        priority_class_relative = classify_priority_tertiles(gene_score)
      ) %>%
      dplyr::filter(gene_score > 0) %>%
      dplyr::arrange(
        dplyr::desc(gene_score),
        dplyr::desc(n_supporting_blocks),
        gene
      )
  })
  
  gtex_pos_summary_r <- reactive({
    normalize_chr_label <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x
    }
    
    dt <- tryCatch(all_gene_bridges(), error = function(e) tibble::tibble())
    
    validate(
      need(is.data.frame(dt) && nrow(dt) > 0, "No gene bridges available.")
    )
    
    # Assegurem que les columnes opcionals existeixin
    if (!"gtex_gene_name" %in% names(dt)) dt$gtex_gene_name <- NA_character_
    if (!"gtex_gene_id" %in% names(dt)) dt$gtex_gene_id <- NA_character_
    if (!"gtex_tissues" %in% names(dt)) dt$gtex_tissues <- NA_character_
    if (!"gtex_variant_id" %in% names(dt)) dt$gtex_variant_id <- NA_character_
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*|\\s*,\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    dt %>%
      dplyr::mutate(
        source_app = tolower(trimws(as.character(source_app))),
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr),
        gene = trimws(as.character(gene)),
        gtex_gene_name = dplyr::coalesce(as.character(gtex_gene_name), ""),
        gtex_gene_id = dplyr::coalesce(as.character(gtex_gene_id), ""),
        gtex_tissues = dplyr::coalesce(as.character(gtex_tissues), ""),
        gtex_variant_id = dplyr::coalesce(as.character(gtex_variant_id), "")
      ) %>%
      dplyr::filter(
        source_app == "gtex",
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(chr), nzchar(chr),
        !is.na(gene), nzchar(gene)
      ) %>%
      dplyr::group_by(cluster_id, chr, gene) %>%
      dplyr::summarise(
        gtex_gene_name = collapse_unique_semicolon(gtex_gene_name),
        gtex_gene_id = collapse_unique_semicolon(gtex_gene_id),
        gtex_tissues = collapse_unique_semicolon(gtex_tissues),
        gtex_variant_id = collapse_unique_semicolon(gtex_variant_id),
        .groups = "drop"
      )
  })
  
  global_block_ranges_r <- reactive({
    normalize_chr_label <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x
    }
    
    ld_global_blocks_plot_df() %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr_num),
        block_id = trimws(as.character(block_id)),
        block_start = suppressWarnings(as.numeric(block_start)),
        block_end = suppressWarnings(as.numeric(block_end))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(chr), nzchar(chr),
        !is.na(block_id), nzchar(block_id),
        !is.na(block_start),
        !is.na(block_end)
      ) %>%
      dplyr::distinct()
  })
  
  gene_gwas_hit_score_audit_df <- reactive({
    ga <- gwas_hit_priority_mod$gwas_hit_match_audit_summary_with_gene_df()
    gh <- gwas_hit_priority_df_v2()
    canon_genes <- canonical_prioritized_cluster_genes_df()
    block_ranges_df <- global_block_ranges_r()
    gtex_src <- gtex_pos_summary_r()
    
    validate(
      need(is.data.frame(ga) && nrow(ga) > 0, "No GWAS hit × app × gene audit data available."),
      need(is.data.frame(gh) && nrow(gh) > 0, "No GWAS hit priority data available."),
      need(is.data.frame(canon_genes) && nrow(canon_genes) > 0, "No canonical prioritized genes available."),
      need(is.data.frame(block_ranges_df) && nrow(block_ranges_df) > 0, "No global LD block ranges available.")
    )
    
    normalize_chr_label <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x
    }
    
    app_weights <- c(
      catalog = 3.0,
      gtex    = 2.5,
      nonsyn  = 4.0,
      ewas    = 2.0
    )
    
    state_weights <- c(
      MATCH  = 1.0,
      MARKER = 0.5
    )
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(sort(unique(x)), collapse = "; ")
    }
    
    # ------------------------------------------------------------
    # 1) Base de gens canònics
    # ------------------------------------------------------------
    canon_gene_keys <- canon_genes %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        gene = trimws(as.character(gene))
      ) %>%
      tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s*,\\s*") %>%
      dplyr::mutate(gene = trimws(gene)) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::distinct()
    
    # ------------------------------------------------------------
    # 2) Base audit gene × GWAS hit × app
    # ------------------------------------------------------------
    ga_base <- ga %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr),
        gwas_hit = trimws(as.character(gwas_hit)),
        gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
        source_app = tolower(trimws(as.character(source_app))),
        score_app = dplyr::case_when(
          source_app %in% c("ewasdis", "ewastum") ~ "ewas",
          TRUE ~ source_app
        ),
        link_state = trimws(as.character(link_state)),
        gene = trimws(as.character(nearest_gene)),
        matched_ids = dplyr::coalesce(as.character(matched_ids), ""),
        matched_pos = dplyr::coalesce(as.character(matched_pos), ""),
        hit_label = dplyr::coalesce(as.character(hit_label), ""),
        nearest_gene_dist_bp = suppressWarnings(as.numeric(nearest_gene_dist_bp))
      ) %>%
      dplyr::filter(
        link_state %in% c("MATCH", "MARKER"),
        !is.na(source_app), nzchar(source_app),
        !is.na(score_app), nzchar(score_app),
        score_app %in% names(app_weights),
        !is.na(gene), nzchar(gene),
        tolower(gene) != "numeric",
        !grepl("^[0-9]+$", gene),
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(gwas_hit), nzchar(gwas_hit),
        !is.na(gwas_pos)
      ) %>%
      tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s*,\\s*") %>%
      dplyr::mutate(gene = trimws(gene)) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::inner_join(canon_gene_keys, by = c("cluster_id", "gene")) %>%
      dplyr::distinct(
        gene, cluster_id, chr, gwas_hit, gwas_pos,
        source_app, score_app, link_state, matched_ids, matched_pos, hit_label,
        .keep_all = TRUE
      )
    
    # ------------------------------------------------------------
    # 3) Afegir score global del GWAS hit
    # ------------------------------------------------------------
    gh_base <- gh %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr),
        gwas_hit = trimws(as.character(gwas_hit)),
        gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
        gwas_hit_priority_score = suppressWarnings(as.numeric(gwas_hit_priority_score))
      ) %>%
      dplyr::select(
        cluster_id, chr, gwas_hit, gwas_pos,
        gwas_hit_priority_score,
        priority_class,
        priority_class_relative
      ) %>%
      dplyr::distinct()
    
    ga_base <- ga_base %>%
      dplyr::left_join(
        gh_base,
        by = c("cluster_id", "chr", "gwas_hit", "gwas_pos")
      ) %>%
      dplyr::mutate(
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0)
      )
    
    # ------------------------------------------------------------
    # 4) Assignar block_id per solapament amb blocs globals
    # ------------------------------------------------------------
    ga_dt <- ga_base %>%
      dplyr::transmute(
        gene,
        cluster_id,
        chr,
        gwas_hit,
        gwas_pos,
        hit_start = gwas_pos,
        hit_end = gwas_pos,
        source_app,
        score_app,
        link_state,
        matched_ids,
        matched_pos,
        hit_label,
        nearest_gene_dist_bp,
        gwas_hit_priority_score,
        priority_class,
        priority_class_relative
      ) %>%
      data.table::as.data.table()
    
    block_dt <- block_ranges_df %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr),
        block_id = trimws(as.character(block_id)),
        block_start = suppressWarnings(as.numeric(block_start)),
        block_end = suppressWarnings(as.numeric(block_end))
      ) %>%
      data.table::as.data.table()
    
    data.table::setkey(ga_dt, cluster_id, chr, hit_start, hit_end)
    data.table::setkey(block_dt, cluster_id, chr, block_start, block_end)
    
    joined_dt <- data.table::foverlaps(
      x = ga_dt,
      y = block_dt,
      by.x = c("cluster_id", "chr", "hit_start", "hit_end"),
      by.y = c("cluster_id", "chr", "block_start", "block_end"),
      type = "within",
      nomatch = NA
    )
    
    ga_with_blocks <- joined_dt %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        block_id = dplyr::coalesce(as.character(block_id), ""),
        has_block = !is.na(block_id) & nzchar(block_id),
        score_unit_id = dplyr::case_when(
          score_app %in% c("gtex", "ewas") & has_block ~ paste(gene, cluster_id, chr, score_app, block_id, sep = "||"),
          score_app %in% c("gtex", "ewas") & !has_block ~ paste(gene, cluster_id, chr, score_app, gwas_hit, gwas_pos, sep = "||"),
          link_state == "MARKER" & has_block ~ paste(gene, cluster_id, chr, score_app, block_id, sep = "||"),
          TRUE ~ paste(gene, cluster_id, chr, score_app, gwas_hit, gwas_pos, sep = "||")
        )
      )
    
    # ------------------------------------------------------------
    # 5) MATCH funcionals: GTEx/EWAS col·lapsats per block
    #    (o per hit si no tenen block)
    #    NONSYN es manté per hit
    # ------------------------------------------------------------
    match_functional <- ga_with_blocks %>%
      dplyr::filter(
        link_state == "MATCH",
        score_app %in% c("gtex", "ewas", "nonsyn")
      ) %>%
      dplyr::group_by(score_unit_id) %>%
      dplyr::arrange(dplyr::desc(gwas_hit_priority_score), gwas_pos, gwas_hit, .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        unit_type = "MATCH",
        gene_score_component = unname(app_weights[score_app]) * state_weights["MATCH"]
      )
    
    # ------------------------------------------------------------
    # 6) MATCH catalog: només si no hi ha MATCH funcional del mateix hit
    # ------------------------------------------------------------
    match_functional_keys <- match_functional %>%
      dplyr::transmute(
        gene, cluster_id, chr, gwas_hit, gwas_pos
      ) %>%
      dplyr::distinct()
    
    match_catalog <- ga_with_blocks %>%
      dplyr::filter(
        link_state == "MATCH",
        score_app == "catalog"
      ) %>%
      dplyr::anti_join(
        match_functional_keys,
        by = c("gene", "cluster_id", "chr", "gwas_hit", "gwas_pos")
      ) %>%
      dplyr::group_by(score_unit_id) %>%
      dplyr::arrange(dplyr::desc(gwas_hit_priority_score), gwas_pos, gwas_hit, .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        unit_type = "MATCH",
        gene_score_component = unname(app_weights[score_app]) * state_weights["MATCH"]
      )
    
    # ------------------------------------------------------------
    # 7) MARKER: GTEx/EWAS col·lapsats per block o per hit si no block
    #    ALTRES apps: mantenim criteri general per score_unit_id
    # ------------------------------------------------------------
    marker_units <- ga_with_blocks %>%
      dplyr::filter(
        link_state == "MARKER"
      ) %>%
      dplyr::group_by(score_unit_id) %>%
      dplyr::arrange(dplyr::desc(gwas_hit_priority_score), gwas_pos, gwas_hit, .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        unit_type = "MARKER",
        gene_score_component = unname(app_weights[score_app]) * state_weights["MARKER"]
      )
    
    # ------------------------------------------------------------
    # 8) Unir unitats
    # ------------------------------------------------------------
    units_df <- dplyr::bind_rows(
      match_functional,
      match_catalog,
      marker_units
    ) %>%
      dplyr::mutate(
        apps_supported = source_app,
        score_family = score_app,
        n_apps_supported = 1L,
        match_apps = dplyr::if_else(unit_type == "MATCH", source_app, ""),
        marker_apps = dplyr::if_else(unit_type == "MARKER", source_app, ""),
        n_match_apps = dplyr::if_else(unit_type == "MATCH", 1L, 0L),
        n_marker_apps = dplyr::if_else(unit_type == "MARKER", 1L, 0L)
      )
    
    # ------------------------------------------------------------
    # 9) Enriquir amb metadata GTEx des del gene bridge
    # ------------------------------------------------------------
    if (!is.data.frame(gtex_src) || !nrow(gtex_src)) {
      gtex_src <- tibble::tibble(
        cluster_id = character(),
        chr = character(),
        gene = character(),
        gtex_gene_name = character(),
        gtex_gene_id = character(),
        gtex_tissues = character(),
        gtex_variant_id = character()
      )
    }
    
    units_df <- units_df %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        chr = normalize_chr_label(chr),
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::left_join(
        gtex_src %>%
          dplyr::mutate(
            cluster_id = trimws(as.character(cluster_id)),
            chr = normalize_chr_label(chr),
            gene = trimws(as.character(gene))
          ) %>%
          dplyr::group_by(cluster_id, chr, gene) %>%
          dplyr::summarise(
            gtex_gene_name = collapse_unique_semicolon(gtex_gene_name),
            gtex_gene_id = collapse_unique_semicolon(gtex_gene_id),
            gtex_tissues = collapse_unique_semicolon(gtex_tissues),
            gtex_variant_id = collapse_unique_semicolon(gtex_variant_id),
            .groups = "drop"
          ),
        by = c("cluster_id", "chr", "gene")
      ) %>%
      dplyr::mutate(
        gtex_gene_name = dplyr::if_else(score_app == "gtex", dplyr::coalesce(as.character(gtex_gene_name), ""), ""),
        gtex_gene_id = dplyr::if_else(score_app == "gtex", dplyr::coalesce(as.character(gtex_gene_id), ""), ""),
        gtex_tissues = dplyr::if_else(score_app == "gtex", dplyr::coalesce(as.character(gtex_tissues), ""), ""),
        gtex_variant_id = dplyr::if_else(score_app == "gtex", dplyr::coalesce(as.character(gtex_variant_id), ""), "")
      )
    
    # ------------------------------------------------------------
    # 10) Selecció final
    # ------------------------------------------------------------
    units_df %>%
      dplyr::select(
        gene,
        cluster_id,
        chr,
        block_id,
        gwas_hit,
        gwas_pos,
        hit_label,
        nearest_gene_dist_bp,
        matched_ids,
        matched_pos,
        gtex_gene_name,
        gtex_gene_id,
        gtex_tissues,
        gtex_variant_id,
        source_app,
        score_app = score_family,
        link_state = unit_type,
        score_unit_id,
        gene_score_component,
        gwas_hit_priority_score,
        apps_supported,
        n_apps_supported,
        match_apps,
        marker_apps,
        n_match_apps,
        n_marker_apps,
        priority_class,
        priority_class_relative
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(
        top_gwas_hit_score = if (all(is.na(gwas_hit_priority_score))) 0 else max(gwas_hit_priority_score, na.rm = TRUE),
        is_top_gwas_hit_for_gene = dplyr::if_else(
          !is.na(gwas_hit_priority_score) & gwas_hit_priority_score == top_gwas_hit_score,
          TRUE,
          FALSE
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        gene,
        dplyr::desc(gwas_hit_priority_score),
        dplyr::desc(gene_score_component),
        cluster_id,
        chr,
        block_id,
        score_app,
        source_app,
        gwas_pos,
        gwas_hit
      )
  })
  

  gene_priority_comparison_df <- reactive({
    audit_sum <- gene_priority_summary_from_audit_df()
    final_df <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(final_df) && nrow(final_df) > 0, "No final prioritized gene table available.")
    )
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    fin <- final_df %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        cluster_ids = dplyr::coalesce(as.character(cluster_ids), ""),
        gene_score = suppressWarnings(as.numeric(gene_score)),
        n_gwas_hits = suppressWarnings(as.integer(n_gwas_hits)),
        n_clusters = suppressWarnings(as.integer(n_clusters)),
        top_gwas_hit_score = suppressWarnings(as.numeric(top_gwas_hit_score)),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        n_apps_supported = suppressWarnings(as.integer(n_apps_supported))
      ) %>%
      dplyr::select(
        gene,
        cluster_ids,
        gene_score,
        n_gwas_hits,
        n_clusters,
        top_gwas_hit_score,
        apps_supported,
        n_apps_supported
      )
    
    aud <- if (is.data.frame(audit_sum) && nrow(audit_sum)) {
      audit_sum %>%
        dplyr::mutate(
          gene = trimws(as.character(gene)),
          cluster_ids_from_audit = dplyr::coalesce(as.character(cluster_ids), ""),
          gene_score_from_audit = suppressWarnings(as.numeric(gene_score_from_audit)),
          n_gwas_hits_from_audit = suppressWarnings(as.integer(n_gwas_hits_from_audit)),
          n_clusters_from_audit = suppressWarnings(as.integer(n_clusters_from_audit)),
          top_gwas_hit_score_from_audit = suppressWarnings(as.numeric(top_gwas_hit_score_from_audit)),
          apps_supported_from_audit = dplyr::coalesce(as.character(apps_supported_from_audit), ""),
          n_apps_supported_from_audit = suppressWarnings(as.integer(n_apps_supported_from_audit))
        ) %>%
        dplyr::select(
          gene,
          cluster_ids_from_audit,
          gene_score_from_audit,
          n_gwas_hits_from_audit,
          n_clusters_from_audit,
          top_gwas_hit_score_from_audit,
          apps_supported_from_audit,
          n_apps_supported_from_audit
        )
    } else {
      tibble::tibble(
        gene = character(),
        cluster_ids_from_audit = character(),
        gene_score_from_audit = numeric(),
        n_gwas_hits_from_audit = integer(),
        n_clusters_from_audit = integer(),
        top_gwas_hit_score_from_audit = numeric(),
        apps_supported_from_audit = character(),
        n_apps_supported_from_audit = integer()
      )
    }
    
    fin %>%
      dplyr::left_join(aud, by = "gene") %>%
      dplyr::mutate(
        has_audit_support = !is.na(gene_score_from_audit),
        
        cluster_ids_from_audit = dplyr::coalesce(cluster_ids_from_audit, ""),
        gene_score_from_audit = dplyr::coalesce(gene_score_from_audit, 0),
        n_gwas_hits_from_audit = dplyr::coalesce(n_gwas_hits_from_audit, 0L),
        n_clusters_from_audit = dplyr::coalesce(n_clusters_from_audit, 0L),
        top_gwas_hit_score_from_audit = dplyr::coalesce(top_gwas_hit_score_from_audit, 0),
        apps_supported_from_audit = dplyr::coalesce(apps_supported_from_audit, ""),
        n_apps_supported_from_audit = dplyr::coalesce(n_apps_supported_from_audit, 0L),
        
        diff_gene_score = round(gene_score - gene_score_from_audit, 10),
        diff_n_gwas_hits = n_gwas_hits - n_gwas_hits_from_audit,
        diff_n_clusters = n_clusters - n_clusters_from_audit,
        diff_top_gwas_hit_score = round(top_gwas_hit_score - top_gwas_hit_score_from_audit, 10),
        diff_n_apps_supported = n_apps_supported - n_apps_supported_from_audit,
        
        # ----------------------------------------------------------
        # Core consistency of score reconstruction
        # For score-0 genes, absence from audit is acceptable
        # ----------------------------------------------------------
        zero_gene_ok = dplyr::coalesce(abs(gene_score) < 1e-09, FALSE) &
          dplyr::coalesce(n_gwas_hits == 0L, FALSE),
        
        gene_score_match = abs(diff_gene_score) < 1e-9,
        n_gwas_hits_match = diff_n_gwas_hits == 0L,
        top_gwas_hit_score_match = abs(diff_top_gwas_hit_score) < 1e-9,
        
        core_score_match = zero_gene_ok | (
          gene_score_match &
            n_gwas_hits_match &
            top_gwas_hit_score_match
        ),
        
        # ----------------------------------------------------------
        # Cluster comparison
        # For score-0 genes with no audit support, cluster mismatch is expected
        # and should not fail the comparison
        # ----------------------------------------------------------
        cluster_ids_match = zero_gene_ok | (
          collapse_unique_semicolon(cluster_ids) ==
            collapse_unique_semicolon(cluster_ids_from_audit)
        ),
        
        # ----------------------------------------------------------
        # App support comparison
        # Only meaningful when the gene has audit-supported hits
        # ----------------------------------------------------------
        apps_supported_match = dplyr::case_when(
          zero_gene_ok ~ TRUE,
          !has_audit_support & gene_score > 0 ~ FALSE,
          TRUE ~ collapse_unique_semicolon(apps_supported) ==
            collapse_unique_semicolon(apps_supported_from_audit)
        ),
        
        n_apps_supported_match = dplyr::case_when(
          zero_gene_ok ~ TRUE,
          !has_audit_support & gene_score > 0 ~ FALSE,
          TRUE ~ diff_n_apps_supported == 0L
        ),
        
        # ----------------------------------------------------------
        # Cluster count comparison
        # Same rule as above
        # ----------------------------------------------------------
        n_clusters_match = dplyr::case_when(
          zero_gene_ok ~ TRUE,
          !has_audit_support & gene_score > 0 ~ FALSE,
          TRUE ~ diff_n_clusters == 0L
        ),
        
        # ----------------------------------------------------------
        # Final flags
        # ----------------------------------------------------------
        all_match_core = core_score_match & cluster_ids_match & n_clusters_match,
        all_match = all_match_core & apps_supported_match & n_apps_supported_match
      ) %>%
      dplyr::arrange(
        dplyr::desc(!all_match),
        dplyr::desc(abs(diff_gene_score)),
        gene
      )
  })
  
  output$gene_gwas_hit_score_audit_dt <- DT::renderDT({
    df <- gene_gwas_hit_score_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No gene × GWAS hit score audit data available.")
    )
    
    DT::datatable(
      df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 15,
        lengthMenu = list(c(15, 25, 50, 100), c("15", "25", "50", "100")),
        scrollX = TRUE,
        autoWidth = FALSE
      ),
      class = "compact stripe hover order-column"
    )
  }, server = FALSE)
  
  # optional to include as internal check 
  output$gene_priority_comparison_dt <- DT::renderDT({
    df <- gene_priority_comparison_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No gene priority comparison data available.")
    )
    
    show_df <- df %>%
      dplyr::select(
        gene,
        cluster_ids,
        cluster_ids_from_audit,
        gene_score,
        gene_score_from_audit,
        diff_gene_score,
        n_gwas_hits,
        n_gwas_hits_from_audit,
        diff_n_gwas_hits,
        top_gwas_hit_score,
        top_gwas_hit_score_from_audit,
        diff_top_gwas_hit_score,
        apps_supported,
        apps_supported_from_audit,
        n_apps_supported,
        n_apps_supported_from_audit,
        diff_n_apps_supported
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 15,
        lengthMenu = list(c(15, 25, 50, 100), c("15", "25", "50", "100")),
        scrollX = TRUE,
        autoWidth = FALSE
      ),
      class = "compact stripe hover order-column"
    )
  }, server = FALSE)
  
  ################################################################################
  
  # ---------------------------
  # Font canónica blocks
  # ---------------------------
  priority_source_blocks <- reactive({
    df <- tryCatch(
      block_overlap_summary_global_df(),
      error = function(e) tibble::tibble()
    )
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0,
           "No global block summary available.")
    )
    
    df
  })
  
  # ---------------------------
  # Prioritized blocks v2
  # ---------------------------
  
  prioritized_block_df_v2 <- reactive({
    
    blk <- priority_source_blocks()
    gh  <- gwas_hit_priority_df_v2()
    
    req(is.data.frame(blk), nrow(blk) > 0)
    req(is.data.frame(gh),  nrow(gh)  > 0)
    
    # ---------------------------
    # Helpers
    # ---------------------------
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      length(unique(vals))
    }
    
    # ---------------------------
    # 1) Detectar columnes reals del summary global
    # ---------------------------
    start_col  <- intersect(c("start", "block_start", "pos_ini"), names(blk))
    end_col    <- intersect(c("end", "block_end", "pos_end"), names(blk))
    score_col  <- intersect(c("block_priority_score", "priority_score", "block_score"), names(blk))
    genes_col  <- intersect(c("genes_name", "genes_in_block", "genes"), names(blk))
    ngenes_col <- intersect(c("n_genes", "n_genes_in_block"), names(blk))
    
    start_col  <- if (length(start_col))  start_col[1]  else NULL
    end_col    <- if (length(end_col))    end_col[1]    else NULL
    score_col  <- if (length(score_col))  score_col[1]  else NULL
    genes_col  <- if (length(genes_col))  genes_col[1]  else NULL
    ngenes_col <- if (length(ngenes_col)) ngenes_col[1] else NULL
    
    # ---------------------------
    # 2) Estandaritzar blocks
    # ---------------------------
    blocks_std <- blk %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        block_id   = trimws(as.character(block_id)),
        block_score_input = if (!is.null(score_col)) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_,
        block_start = if (!is.null(start_col)) suppressWarnings(as.numeric(.data[[start_col]])) else NA_real_,
        block_end   = if (!is.null(end_col))   suppressWarnings(as.numeric(.data[[end_col]]))   else NA_real_,
        genes_in_block = if (!is.null(genes_col)) as.character(.data[[genes_col]]) else "",
        n_genes_in_block = if (!is.null(ngenes_col)) suppressWarnings(as.integer(.data[[ngenes_col]])) else NA_integer_
      ) %>%
      dplyr::mutate(
        genes_in_block   = dplyr::coalesce(genes_in_block, ""),
        n_genes_in_block = dplyr::coalesce(n_genes_in_block, 0L),
        block_size_bp = dplyr::if_else(
          is.finite(block_start) & is.finite(block_end),
          block_end - block_start + 1,
          NA_real_
        ),
        block_size_kb = block_size_bp / 1000
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(block_id), nzchar(block_id)
      ) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
    
    # ---------------------------
    # 3) Assignar GWAS hits v2 als blocs per posició
    # ---------------------------
    block_hit_map <- gh %>%
      dplyr::inner_join(
        blocks_std %>%
          dplyr::select(cluster_id, block_id, block_start, block_end) %>%
          dplyr::distinct(),
        by = "cluster_id",
        relationship = "many-to-many"
      ) %>%
      dplyr::filter(
        is.finite(gwas_pos),
        is.finite(block_start),
        is.finite(block_end),
        gwas_pos >= block_start,
        gwas_pos <= block_end
      ) %>%
      dplyr::distinct(cluster_id, block_id, gwas_hit, gwas_pos, .keep_all = TRUE)
    
    validate(
      need(nrow(block_hit_map) > 0, "No prioritized GWAS hits mapped to global LD blocks.")
    )
    
    # ---------------------------
    # 4) Score canònic del bloc:
    #    top_gwas_hit_score + 0.05 * other_gwas_hit_score_sum
    # ---------------------------
    block_support_summary <- block_hit_map %>%
      dplyr::group_by(cluster_id, block_id) %>%
      dplyr::summarise(
        n_gwas_hits = dplyr::n_distinct(paste(gwas_hit, gwas_pos, sep = "||")),
        gwas_hits = paste(sort(unique(gwas_hit[!is.na(gwas_hit) & nzchar(gwas_hit)])), collapse = "; "),
        apps_supported = collapse_unique_semicolon(gwas_hit_apps),
        
        hit_scores_desc = list(sort(as.numeric(gwas_hit_score), decreasing = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        top_gwas_hit_score = {
          s <- unlist(hit_scores_desc)
          s <- s[is.finite(s)]
          if (!length(s)) 0 else s[1]
        },
        other_gwas_hit_score = {
          s <- unlist(hit_scores_desc)
          s <- s[is.finite(s)]
          if (length(s) <= 1) 0 else sum(s[-1], na.rm = TRUE)
        },
        block_score_from_hits = top_gwas_hit_score + 0.05 * other_gwas_hit_score
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        n_apps_supported = vapply(apps_supported, count_semicolon_items, integer(1))
      ) %>%
      dplyr::select(
        -hit_scores_desc
      )
    
    # ---------------------------
    # 5) Unir info estructural del bloc i afegir bonus per gene content
    # ---------------------------
    out <- block_support_summary %>%
      dplyr::left_join(
        blocks_std,
        by = c("cluster_id", "block_id")
      ) %>%
      dplyr::mutate(
        block_gene_content_bonus = round(
          1.50 * as.numeric(dplyr::coalesce(n_genes_in_block, 0L) > 0) +
            0.75 * log1p(dplyr::coalesce(n_genes_in_block, 0L)),
          2
        ),
        
        # Àlies per compatibilitat amb codi antic
        block_bio_bonus = block_gene_content_bonus,
        
        block_score = block_score_from_hits + block_gene_content_bonus,
        
        priority_class = dplyr::case_when(
          block_score >= 25 ~ "High",
          block_score >= 10 ~ "Medium",
          TRUE              ~ "Low"
        ),
        
        priority_class_relative = classify_priority_tertiles(block_score)
      ) %>%
      dplyr::select(
        cluster_id,
        block_id,
        block_score_from_hits,
        top_gwas_hit_score,
        other_gwas_hit_score,
        block_gene_content_bonus,
        block_bio_bonus,
        block_score,
        n_gwas_hits,
        apps_supported,
        n_apps_supported,
        block_start,
        block_end,
        block_size_bp,
        block_size_kb,
        genes_in_block,
        n_genes_in_block,
        gwas_hits,
        block_score_input,
        priority_class,
        priority_class_relative
      ) %>%
      dplyr::arrange(
        dplyr::desc(block_score),
        dplyr::desc(top_gwas_hit_score),
        dplyr::desc(n_gwas_hits),
        cluster_id,
        block_id
      )
    
    out
  })
  
  ################################################################################
  # ============================================================
  # PRIORITIZED GWAS HITS
  # ============================================================
  
  ns_gwas <- NS("gwas_hit_priority")
  
  gwas_hit_priority_df_v2 <- reactive({
    df <- gwas_hit_priority_mod$gwas_hit_gene_support_long_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No GWAS hit gene-support data available."),
      need(all(c(
        "cluster_id", "chr", "gwas_hit", "gwas_pos",
        "source_app", "link_state", "gene"
      ) %in% names(df)), "Required columns missing in GWAS hit gene-support data.")
    )
    
    app_weights <- c(
      catalog = 3.0,
      gtex    = 2.5,
      nonsyn  = 4.0,
      ewasdis = 2.0,
      ewastum = 2.0
    )
    
    state_weights <- c(
      MATCH  = 1.0,
      MARKER = 0.5
    )
    
    app_levels   <- c("catalog", "gtex", "nonsyn", "ewastum", "ewasdis")
    state_levels <- c("MATCH", "MARKER")
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      paste(sort(unique(vals)), collapse = "; ")
    }
    
    collapse_unique_pipe <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(unique(x), collapse = " | ")
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      as.integer(length(unique(vals)))
    }
    
    # ------------------------------------------------------------
    # GWAS significance bridge
    # ------------------------------------------------------------
    sig_path <- file.path(selected_session_dir(), "gwas_significance_bridge.rds")
    sig_bridge <- tryCatch(readRDS(sig_path), error = function(e) NULL)
    
    if (!is.data.frame(sig_bridge) || !nrow(sig_bridge)) {
      sig_bridge <- tibble::tibble(
        cluster_id = character(),
        chr = character(),
        gwas_pos = numeric(),
        gwas_hit = character(),
        gwas_p = numeric(),
        gwas_logp = numeric()
      )
    } else {
      sig_bridge <- sig_bridge %>%
        dplyr::transmute(
          cluster_id = trimws(as.character(cluster_id)),
          chr = as.character(chr),
          gwas_pos = suppressWarnings(as.numeric(position)),
          gwas_hit = trimws(as.character(rsid)),
          gwas_p = suppressWarnings(as.numeric(p_value)),
          gwas_logp = suppressWarnings(as.numeric(logp))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(gwas_pos),
          !is.na(gwas_p)
        ) %>%
        dplyr::distinct(cluster_id, chr, gwas_pos, gwas_hit, .keep_all = TRUE)
    }
    
    hit_app_support <- df %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        chr = as.character(chr),
        gwas_hit = trimws(as.character(gwas_hit)),
        gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
        source_app = tolower(trimws(as.character(source_app))),
        link_state = trimws(as.character(link_state)),
        gene = trimws(as.character(gene)),
        state_rank = dplyr::case_when(
          link_state == "MATCH"  ~ 2L,
          link_state == "MARKER" ~ 1L,
          TRUE ~ 0L
        ),
        app_weight = dplyr::coalesce(unname(app_weights[source_app]), 0),
        state_weight = dplyr::coalesce(unname(state_weights[link_state]), 0),
        support_score = app_weight * state_weight
      ) %>%
      dplyr::filter(
        source_app %in% app_levels,
        link_state %in% state_levels,
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(gwas_hit), nzchar(gwas_hit),
        !is.na(source_app), nzchar(source_app)
      ) %>%
      dplyr::group_by(cluster_id, chr, gwas_hit, gwas_pos, source_app) %>%
      dplyr::summarise(
        best_link_state = dplyr::case_when(
          max(state_rank, na.rm = TRUE) >= 2 ~ "MATCH",
          max(state_rank, na.rm = TRUE) == 1 ~ "MARKER",
          TRUE ~ ""
        ),
        best_state_rank = max(state_rank, na.rm = TRUE),
        app_weight = dplyr::first(app_weight),
        state_weight = max(state_weight, na.rm = TRUE),
        support_score = max(support_score, na.rm = TRUE),
        genes = collapse_unique_semicolon(gene),
        n_genes = dplyr::n_distinct(gene[!is.na(gene) & nzchar(gene)]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        genes = dplyr::coalesce(as.character(genes), ""),
        n_genes = dplyr::coalesce(as.integer(n_genes), 0L),
        app_state_label = paste0(source_app, ":", best_link_state),
        app_score_label = paste0(source_app, "=", format(round(support_score, 2), nsmall = 2)),
        app_gene_label = dplyr::if_else(
          nzchar(genes),
          paste0(source_app, ": ", genes),
          paste0(source_app, ": -")
        )
      )
    
    hit_level_base <- hit_app_support %>%
      dplyr::group_by(cluster_id, chr, gwas_hit, gwas_pos) %>%
      dplyr::summarise(
        raw_evidence_score = sum(support_score, na.rm = TRUE),
        
        gwas_hit_n_apps = dplyr::n_distinct(source_app),
        gwas_hit_n_match_apps = sum(best_link_state == "MATCH", na.rm = TRUE),
        gwas_hit_n_marker_apps = sum(best_link_state == "MARKER", na.rm = TRUE),
        
        gwas_hit_has_nonsyn_match = any(source_app == "nonsyn" & best_link_state == "MATCH", na.rm = TRUE),
        gwas_hit_has_catalog_match = any(source_app == "catalog" & best_link_state == "MATCH", na.rm = TRUE),
        gwas_hit_has_gtex_match = any(source_app == "gtex" & best_link_state == "MATCH", na.rm = TRUE),
        gwas_hit_has_ewasdis_match = any(source_app == "ewasdis" & best_link_state == "MATCH", na.rm = TRUE),
        gwas_hit_has_ewastum_match = any(source_app == "ewastum" & best_link_state == "MATCH", na.rm = TRUE),
        
        gwas_hit_apps = paste(sort(unique(source_app[!is.na(source_app) & nzchar(source_app)])), collapse = "; "),
        gwas_hit_match_apps = paste(sort(unique(source_app[best_link_state == "MATCH"])), collapse = "; "),
        gwas_hit_marker_apps = paste(sort(unique(source_app[best_link_state == "MARKER"])), collapse = "; "),
        
        linked_genes = collapse_unique_semicolon(genes),
        linked_gene_n = sum(vapply(unique(genes), count_semicolon_items, integer(1))),
        
        support_signature = collapse_unique_pipe(app_state_label[order(match(source_app, app_levels, nomatch = 999))]),
        score_breakdown_core = collapse_unique_pipe(app_score_label[order(match(source_app, app_levels, nomatch = 999))]),
        genes_by_app = collapse_unique_pipe(app_gene_label[order(match(source_app, app_levels, nomatch = 999))]),
        
        .groups = "drop"
      )
    
    out <- hit_level_base %>%
      dplyr::left_join(
        sig_bridge,
        by = c("cluster_id", "chr", "gwas_hit", "gwas_pos")
      ) %>%
      dplyr::mutate(
        linked_gene_n = dplyr::coalesce(as.integer(linked_gene_n), 0L),
        
        gwas_significance_bonus = dplyr::case_when(
          is.na(gwas_p) | gwas_p <= 0 ~ 0,
          TRUE ~ pmax(0, pmin(-log10(gwas_p), 20) - 6) / 3.5
        ),
        
        gwas_hit_score_old = raw_evidence_score,
        gwas_hit_score = raw_evidence_score + gwas_significance_bonus,
        gwas_hit_priority_score = gwas_hit_score,
        
        score_breakdown = dplyr::case_when(
          nzchar(score_breakdown_core) & gwas_significance_bonus > 0 ~ paste0(
            score_breakdown_core,
            " | gwas_significance=",
            format(round(gwas_significance_bonus, 2), nsmall = 2)
          ),
          nzchar(score_breakdown_core) ~ score_breakdown_core,
          gwas_significance_bonus > 0 ~ paste0(
            "gwas_significance=",
            format(round(gwas_significance_bonus, 2), nsmall = 2)
          ),
          TRUE ~ ""
        ),
        
        evidence_summary = dplyr::case_when(
          gwas_hit_n_match_apps > 0 & gwas_hit_n_marker_apps > 0 ~ paste0(
            gwas_hit_n_match_apps, " MATCH + ", gwas_hit_n_marker_apps, " MARKER"
          ),
          gwas_hit_n_match_apps > 0 ~ paste0(gwas_hit_n_match_apps, " MATCH"),
          gwas_hit_n_marker_apps > 0 ~ paste0(gwas_hit_n_marker_apps, " MARKER"),
          TRUE ~ "No support"
        ),
        
        priority_class = dplyr::case_when(
          gwas_hit_has_nonsyn_match & gwas_hit_n_match_apps >= 2 ~ "High",
          gwas_hit_n_match_apps >= 2 ~ "High",
          gwas_hit_n_match_apps >= 1 & gwas_hit_n_apps >= 2 ~ "Medium",
          gwas_hit_n_match_apps >= 1 ~ "Medium",
          gwas_hit_n_marker_apps >= 1 ~ "Low",
          TRUE ~ "Low"
        ),
        
        priority_class_relative = classify_priority_tertiles(gwas_hit_priority_score),
        genes = linked_genes,
        
        gwas_significance_label = dplyr::if_else(
          !is.na(gwas_p),
          paste0(
            format(round(gwas_significance_bonus, 2), nsmall = 2),
            " (p=",
            format(gwas_p, scientific = TRUE, digits = 3),
            ", -log10P=",
            format(round(gwas_logp, 2), nsmall = 2),
            ")"
          ),
          format(round(gwas_significance_bonus, 2), nsmall = 2)
        )
      ) %>%
      dplyr::arrange(
        dplyr::desc(gwas_hit_priority_score),
        dplyr::desc(gwas_hit_n_match_apps),
        dplyr::desc(gwas_hit_n_apps),
        cluster_id,
        gwas_pos,
        gwas_hit
      ) %>%
      dplyr::select(
        cluster_id,
        chr,
        gwas_hit,
        gwas_pos,
        gwas_p,
        gwas_logp,
        gwas_significance_bonus,
        gwas_significance_label,
        gwas_hit_priority_score,
        priority_class,
        priority_class_relative,
        evidence_summary,
        support_signature,
        score_breakdown,
        genes_by_app,
        linked_genes,
        linked_gene_n,
        genes,
        gwas_hit_n_apps,
        gwas_hit_n_match_apps,
        gwas_hit_n_marker_apps,
        gwas_hit_has_nonsyn_match,
        gwas_hit_has_catalog_match,
        gwas_hit_has_gtex_match,
        gwas_hit_has_ewasdis_match,
        gwas_hit_has_ewastum_match,
        gwas_hit_apps,
        gwas_hit_match_apps,
        gwas_hit_marker_apps,
        gwas_hit_score,
        gwas_hit_score_old,
        raw_evidence_score
      )
    
    out
  })
  
# observe({
#   req(gwas_hit_priority_df_v2())
#   cat("\n[gwas_hit_priority_df_v2] column names:\n")
#   print(names(gwas_hit_priority_df_v2()))
#   cat("\n")
# })
  
  output[[ns_gwas("gwas_hit_priority_dt")]] <- DT::renderDT({
    df <- gwas_hit_priority_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized GWAS hits available.")
    )
    
    make_yes_no_badge <- function(x) {
      vapply(x, function(val) {
        if (isTRUE(val)) {
          "<span style='display:inline-block;padding:2px 8px;border-radius:10px;background:#d9ead3;color:#1b5e20;font-weight:600;'>YES</span>"
        } else {
          "<span style='display:inline-block;padding:2px 8px;border-radius:10px;background:#f2f2f2;color:#666;font-weight:600;'>NO</span>"
        }
      }, character(1))
    }
    
    show_df <- df %>%
      dplyr::mutate(
        gwas_hit = make_dbsnp_links(gwas_hit),
        gwas_pos = as.integer(gwas_pos),
        gwas_hit_priority_score = round(gwas_hit_priority_score, 2),
        gwas_significance_bonus = round(gwas_significance_bonus, 2),
        gwas_p = dplyr::case_when(
          is.na(gwas_p) ~ NA_character_,
          TRUE ~ format(gwas_p, scientific = TRUE, digits = 3)
        ),
        linked_genes = make_genecards_links(linked_genes),
        genes = make_genecards_links(genes),
        gwas_hit_has_nonsyn_match = make_yes_no_badge(gwas_hit_has_nonsyn_match),
        gwas_hit_has_catalog_match = make_yes_no_badge(gwas_hit_has_catalog_match),
        gwas_hit_has_gtex_match = make_yes_no_badge(gwas_hit_has_gtex_match),
        gwas_hit_has_ewasdis_match = make_yes_no_badge(gwas_hit_has_ewasdis_match),
        gwas_hit_has_ewastum_match = make_yes_no_badge(gwas_hit_has_ewastum_match),
        priority_class = make_priority_badge(priority_class),
        priority_class_relative = make_priority_badge(priority_class_relative)
      ) %>%
      dplyr::select(
        cluster_id,
        chr,
        gwas_hit,
        gwas_pos,
        score = gwas_hit_priority_score,
        `priority_class (ABS)` = priority_class,
        `priority_class (REL)` = priority_class_relative,
        evidence_summary,
        support_signature,
        score_breakdown,
        genes_by_app,
        linked_genes,
        linked_gene_n,
        gwas_hit_n_apps,
        gwas_hit_n_match_apps,
        gwas_hit_n_marker_apps,
        gwas_hit_has_nonsyn_match,
        gwas_hit_has_catalog_match,
        gwas_hit_has_gtex_match,
        gwas_hit_has_ewasdis_match,
        gwas_hit_has_ewastum_match,
        gwas_p,
        gwas_significance_bonus
      )
    
    hidden_cols <- which(names(show_df) %in% c(
      "linked_gene_n",
      "gwas_hit_n_apps",
      "gwas_hit_n_match_apps",
      "gwas_hit_n_marker_apps",
      "gwas_hit_has_nonsyn_match",
      "gwas_hit_has_catalog_match",
      "gwas_hit_has_gtex_match",
      "gwas_hit_has_ewasdis_match",
      "gwas_hit_has_ewastum_match",
      "gwas_p",
      "gwas_significance_bonus"
    )) - 1L
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        columnDefs = list(
          list(
            targets = hidden_cols,
            visible = FALSE
          )
        ),
        buttons = list(
          list(
            extend = "copyHtml5",
            text = "Copy",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "csvHtml5",
            text = "CSV",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "excelHtml5",
            text = "Excel",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          )
        )
      )
    )
  }, server = FALSE)

  # ============================================================
  # PRIORITIZED GENES
  # ============================================================
  
  output$prioritized_gene_dt <- DT::renderDT({
    df <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized genes available.")
    )
    
    show_df <- df %>%
      dplyr::mutate(
        gene = make_genecards_links(gene),
        gene_score = round(gene_score, 2),
        top_gwas_hit_score = round(top_gwas_hit_score, 2),
        cluster_ids = dplyr::coalesce(as.character(cluster_ids), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        match_support = dplyr::coalesce(as.character(match_support), ""),
        marker_support = dplyr::coalesce(as.character(marker_support), ""),
        n_supporting_blocks = dplyr::coalesce(as.integer(n_supporting_blocks), 0L),
        n_clusters = dplyr::coalesce(as.integer(n_clusters), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        
        evidence_summary = paste0(
          n_supporting_blocks, " LD blocks",
          " | ",
          n_clusters, " clusters",
          " | ",
          n_apps_supported, " apps"
        ),
        
        support_signature = dplyr::case_when(
          nzchar(gene_link_mode) & nzchar(apps_supported) ~ paste0(gene_link_mode, " | ", apps_supported),
          nzchar(gene_link_mode) ~ gene_link_mode,
          nzchar(apps_supported) ~ apps_supported,
          TRUE ~ ""
        ),
        
        score_breakdown = paste0(
          "total=", format(gene_score, nsmall = 2),
          " | top_block=", format(top_gwas_hit_score, nsmall = 2)
        ),
        
        priority_class = make_priority_badge(priority_class),
        priority_class_relative = make_priority_badge(priority_class_relative)
      ) %>%
      dplyr::select(
        gene,
        cluster_ids,
        score = gene_score,
        `priority_class (ABS)` = priority_class,
        `priority_class (REL)` = priority_class_relative,
        evidence_summary,
        support_signature,
        score_breakdown,
        gene_evidence_types,
        apps_supported,
        top_gwas_hit_score,
        n_supporting_blocks,
        n_clusters,
        n_apps_supported,
        gene_link_mode,
        gene_source_apps,
        match_support,
        marker_support
      )
    
    hidden_cols <- which(names(show_df) %in% c(
      "top_gwas_hit_score",
      "n_supporting_blocks",
      "n_clusters",
      "n_apps_supported",
      "gene_link_mode",
      "gene_source_apps",
      "match_support",
      "marker_support"
    )) - 1L
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        columnDefs = list(
          list(
            targets = hidden_cols,
            visible = FALSE
          )
        ),
        buttons = list(
          list(
            extend = "copyHtml5",
            text = "Copy",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "csvHtml5",
            text = "CSV",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "excelHtml5",
            text = "Excel",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          )
        )
      )
    )
  }, server = FALSE)
 
  # ============================================================
  # PRIORITIZED BLOCKS
  # ============================================================

  output$prioritized_block_dt <- DT::renderDT({
    df <- prioritized_block_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized blocks available.")
    )
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    split_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(character(0))
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      unique(vals)
    }
    
    html_escape_basic <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x <- gsub("\"", "&quot;", x, fixed = TRUE)
      x
    }
    
    make_one_dbsnp_link <- function(x) {
      x2 <- trimws(as.character(x))
      if (!nzchar(x2)) return("")
      
      url <- paste0(
        "https://www.ncbi.nlm.nih.gov/snp/",
        utils::URLencode(x2, reserved = TRUE)
      )
      
      paste0("<a href='", url, "' target='_blank'>", html_escape_basic(x2), "</a>")
    }
    
    make_dbsnp_links_collapsible <- function(x, max_visible = 5L) {
      vals <- split_semicolon_items(x)
      
      if (!length(vals)) return("")
      
      links <- vapply(vals, make_one_dbsnp_link, character(1))
      
      if (length(links) <= max_visible) {
        return(paste(links, collapse = "; "))
      }
      
      first_part <- paste(links[seq_len(max_visible)], collapse = "; ")
      rest_part <- paste(links[(max_visible + 1):length(links)], collapse = "; ")
      
      paste0(
        "<details>",
        "<summary>", length(links), " GWAS hits</summary>",
        "<div style='margin-top:6px;'>",
        first_part,
        if (nzchar(rest_part)) paste0("; ", rest_part) else "",
        "</div>",
        "</details>"
      )
    }
    
    needed_chr <- c(
      "cluster_id",
      "block_id",
      "genes_in_block",
      "gwas_hits",
      "priority_class",
      "priority_class_relative",
      "apps_supported",
      "source_apps",
      "support_signature",
      "gene_link_mode"
    )
    
    needed_num <- c(
      "block_score",
      "top_gwas_hit_score",
      "block_size_kb",
      "n_gwas_hits",
      "n_genes",
      "n_apps_supported",
      "n_clusters"
    )
    
    for (nm in needed_chr) {
      if (!nm %in% names(df)) df[[nm]] <- ""
    }
    
    for (nm in needed_num) {
      if (!nm %in% names(df)) df[[nm]] <- 0
    }
    
    show_df <- df %>%
      dplyr::mutate(
        cluster_id = dplyr::coalesce(as.character(cluster_id), ""),
        block_id = dplyr::coalesce(as.character(block_id), ""),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), ""),
        gwas_hits = dplyr::coalesce(as.character(gwas_hits), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        source_apps = dplyr::coalesce(as.character(source_apps), ""),
        support_signature = dplyr::coalesce(as.character(support_signature), ""),
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        
        block_score = round(as.numeric(block_score), 2),
        top_gwas_hit_score = round(as.numeric(top_gwas_hit_score), 2),
        block_size_kb = round(as.numeric(block_size_kb), 2),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        n_clusters = dplyr::coalesce(as.integer(n_clusters), 0L),
        
        n_genes_from_label = vapply(genes_in_block, count_semicolon_items, integer(1)),
        n_genes = dplyr::coalesce(as.integer(n_genes), 0L),
        n_genes = dplyr::if_else(
          n_genes > 0,
          n_genes,
          as.integer(n_genes_from_label)
        ),
        
        evidence_summary = paste0(
          n_gwas_hits, " GWAS hits",
          " | ",
          n_genes, " genes",
          " | ",
          n_apps_supported, " apps"
        ),
        
        support_signature = dplyr::case_when(
          nzchar(support_signature) ~ support_signature,
          nzchar(apps_supported) ~ apps_supported,
          nzchar(source_apps) ~ source_apps,
          TRUE ~ ""
        ),
        
        score_breakdown = paste0(
          "total=", format(block_score, nsmall = 2),
          " | top_hit=", format(top_gwas_hit_score, nsmall = 2),
          " | size_kb=", format(block_size_kb, nsmall = 2)
        ),
        
        genes_in_block = make_genecards_links(genes_in_block),
        gwas_hits = vapply(gwas_hits, make_dbsnp_links_collapsible, character(1)),
        priority_class = make_priority_badge(priority_class),
        priority_class_relative = make_priority_badge(priority_class_relative)
      ) %>%
      dplyr::select(
        cluster_id,
        block_id,
        score = block_score,
        `priority_class (ABS)` = priority_class,
        `priority_class (REL)` = priority_class_relative,
        evidence_summary,
        support_signature,
        score_breakdown,
        genes_in_block,
        gwas_hits,
        top_gwas_hit_score,
        block_size_kb,
        n_gwas_hits,
        n_genes,
        n_apps_supported,
        n_clusters,
        apps_supported,
        source_apps,
        gene_link_mode
      )
    
    hidden_cols <- which(names(show_df) %in% c(
      "top_gwas_hit_score",
      "block_size_kb",
      "n_gwas_hits",
      "n_genes",
      "n_apps_supported",
      "n_clusters",
      "apps_supported",
      "source_apps",
      "gene_link_mode"
    )) - 1L
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        columnDefs = list(
          list(
            targets = hidden_cols,
            visible = FALSE
          )
        ),
        buttons = list(
          list(
            extend = "copyHtml5",
            text = "Copy",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "csvHtml5",
            text = "CSV",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          ),
          list(
            extend = "excelHtml5",
            text = "Excel",
            exportOptions = list(
              columns = ":visible",
              modifier = list(
                page = "all",
                search = "applied",
                order = "applied"
              )
            )
          )
        )
      )
    )
  }, server = FALSE)
  
  #############################################################################
  
  nonsyn_hit_gene_lookup_df <- reactive({
    req(session_dir_r())
    
    fp <- file.path(session_dir_r(), "nonsyn_hit_gene_lookup.rds")
    
    validate(
      need(file.exists(fp), "nonsyn_hit_gene_lookup.rds not found in selected session.")
    )
    
    df <- tryCatch(readRDS(fp), error = function(e) NULL)
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "NonSyn hit→gene lookup is empty.")
    )
    
    df %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        source_app = tolower(trimws(as.character(source_app))),
        hit_key = trimws(as.character(hit_key)),
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(source_app), nzchar(source_app),
        !is.na(hit_key), nzchar(hit_key),
        !is.na(gene), nzchar(gene)
      ) %>%
      dplyr::distinct()
  })
  
  # ============================================================
  # GENE PRIORITY COMPONENTS · v2
  # ============================================================
  
  prioritized_genes_global_plot_df <- reactive({
    df <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized genes available."),
      need(all(c("top_hit_support_score", "other_hit_support_score", "other_hits_weight") %in% names(df)),
           "Canonical gene-score columns are not available in prioritized_gene_df_v2().")
    )
    
    df_base <- df %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0),
        top_hit_support_score = dplyr::coalesce(as.numeric(top_hit_support_score), 0),
        other_hit_support_score = dplyr::coalesce(as.numeric(other_hit_support_score), 0),
        other_hits_weight = dplyr::coalesce(as.numeric(other_hits_weight), 0.05),
        weighted_other_support = other_hits_weight * other_hit_support_score
      ) %>%
      dplyr::filter(
        !is.na(gene),
        nzchar(gene),
        gene_score > 0
      )
    
    dplyr::bind_rows(
      df_base %>%
        dplyr::transmute(
          gene,
          gene_score,
          component = "Top support unit",
          value = top_hit_support_score,
          hover_txt = paste0(
            "<b>", gene, "</b><br>",
            "Component: Top support unit<br>",
            "Value: ", formatC(top_hit_support_score, format = "f", digits = 2), "<br>",
            "Total gene score: ", formatC(gene_score, format = "f", digits = 2)
          )
        ),
      df_base %>%
        dplyr::transmute(
          gene,
          gene_score,
          component = "Other support units (0.05x)",
          value = weighted_other_support,
          hover_txt = paste0(
            "<b>", gene, "</b><br>",
            "Component: Other support units (weighted)<br>",
            "Raw other support: ", formatC(other_hit_support_score, format = "f", digits = 2), "<br>",
            "Weight: ", formatC(other_hits_weight, format = "f", digits = 2), "<br>",
            "Weighted value: ", formatC(weighted_other_support, format = "f", digits = 2), "<br>",
            "Total gene score: ", formatC(gene_score, format = "f", digits = 2)
          )
        )
    ) %>%
      dplyr::mutate(
        component = factor(
          component,
          levels = c(
            "Top support unit",
            "Other support units (0.05x)"
          )
        )
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(
        gene_total_score = max(gene_score, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        dplyr::desc(gene_total_score),
        gene,
        component
      )
  })
  
  output$prioritized_genes_global_plot <- plotly::renderPlotly({
    df <- prioritized_genes_global_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized gene plot data available.")
    )
    
    gene_levels <- df %>%
      dplyr::distinct(gene, gene_total_score) %>%
      dplyr::arrange(dplyr::desc(gene_total_score), gene) %>%
      dplyr::pull(gene)
    
    df_plot <- df %>%
      dplyr::mutate(
        gene = factor(gene, levels = rev(gene_levels))
      )
    
    plotly::plot_ly(
      data = df_plot,
      y = ~gene,
      x = ~value,
      color = ~component,
      type = "bar",
      orientation = "h",
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis = list(
          title = "Canonical gene score",
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE
        ),
        legend = list(
          orientation = "h",
          x = 0,
          y = -0.15
        ),
        margin = list(l = 140, r = 30, t = 30, b = 90)
      )
  })
  
  observeEvent(input$show_gene_priority_plot, {
    df <- prioritized_gene_df_v2()
    req(is.data.frame(df), nrow(df) > 0)
    
    showModal(
      modalDialog(
        title = "Canonical gene score components",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$div(
          style = "font-size:13px; color:#444; margin-bottom:10px;",
          HTML("Bars show the canonical gene score split into the top support unit and the weighted contribution of all remaining support units.")
        ),
        
        plotly::plotlyOutput("prioritized_genes_global_plot", height = "650px")
      )
    )
  })
  
  #############################################################################
  
  
  prioritized_genes_evidence_types_plot_df <- reactive({
    genes_df <- prioritized_gene_df_v2()
    audit_df <- gene_gwas_hit_score_audit_df()
    
    validate(
      need(is.data.frame(genes_df) && nrow(genes_df) > 0, "No prioritized genes available."),
      need(is.data.frame(audit_df) && nrow(audit_df) > 0, "No gene audit data available.")
    )
    
    split_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return(character(0))
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      unique(vals)
    }
    
    # ------------------------------------------------------------
    # 1) Base gens prioritzats (només score > 0)
    # ------------------------------------------------------------
    gene_base <- genes_df %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        n_supporting_blocks = dplyr::coalesce(as.integer(n_supporting_blocks), 0L),
        n_clusters = dplyr::coalesce(as.integer(n_clusters), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), "")
      ) %>%
      dplyr::filter(
        !is.na(gene),
        nzchar(gene),
        gene_score > 0
      )
    
    # ------------------------------------------------------------
    # 2) Contribució real per tipus d'evidència des de l'audit canònic
    #    - usa score_app (EWAS unificat)
    #    - usa score_unit_id (unitats ja col·lapsades correctament)
    #    - hits sense block_id es mantenen com a unitats independents
    #    - reparteix el gene score real:
    #         top unit -> 100%
    #         altres -> 0.05
    # ------------------------------------------------------------
    audit_units <- audit_df %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        block_id = dplyr::coalesce(as.character(block_id), ""),
        score_app = dplyr::coalesce(as.character(score_app), ""),
        score_unit_id = dplyr::coalesce(as.character(score_unit_id), ""),
        gene_score_component = dplyr::coalesce(as.numeric(gene_score_component), 0),
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0),
        hit_support_score = gene_score_component * gwas_hit_priority_score,
        evidence_type = dplyr::case_when(
          score_app == "catalog" ~ "catalog_gene",
          score_app == "gtex" ~ "eqtl_gene",
          score_app == "nonsyn" ~ "nonsyn_gene",
          score_app == "ewas" ~ "ewas_nearby_gene",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(
        !is.na(gene), nzchar(gene),
        !is.na(score_unit_id), nzchar(score_unit_id),
        !is.na(evidence_type), nzchar(evidence_type)
      ) %>%
      dplyr::group_by(gene) %>%
      dplyr::arrange(
        dplyr::desc(hit_support_score),
        dplyr::desc(gwas_hit_priority_score),
        score_unit_id,
        .by_group = TRUE
      ) %>%
      dplyr::mutate(
        unit_rank_within_gene = dplyr::row_number(),
        unit_weight = dplyr::if_else(unit_rank_within_gene == 1L, 1, 0.05),
        weighted_unit_score = hit_support_score * unit_weight,
        support_block_or_unit = dplyr::if_else(
          nzchar(block_id),
          block_id,
          score_unit_id
        )
      ) %>%
      dplyr::ungroup()
    
    audit_long <- audit_units %>%
      dplyr::group_by(gene, evidence_type) %>%
      dplyr::summarise(
        n_supporting_blocks_type = dplyr::n_distinct(support_block_or_unit),
        score_from_blocks_type = sum(weighted_unit_score, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ------------------------------------------------------------
    # 3) Physical overlap com a component visual binari
    #    (es manté com a senyal visual, no com a score canònic)
    # ------------------------------------------------------------
    physical_df <- gene_base %>%
      dplyr::mutate(
        has_physical = vapply(
          gene_evidence_types,
          function(x) "physical_overlap" %in% split_unique_semicolon(x),
          logical(1)
        )
      ) %>%
      dplyr::filter(has_physical) %>%
      dplyr::transmute(
        gene,
        evidence_type = "physical_overlap",
        n_supporting_blocks_type = 1L,
        score_from_blocks_type = 1
      )
    
    # ------------------------------------------------------------
    # 4) Unió final
    # ------------------------------------------------------------
    binded <- dplyr::bind_rows(audit_long, physical_df) %>%
      dplyr::group_by(gene, evidence_type) %>%
      dplyr::summarise(
        n_supporting_blocks_type = max(n_supporting_blocks_type, na.rm = TRUE),
        score_from_blocks_type = sum(score_from_blocks_type, na.rm = TRUE),
        .groups = "drop"
      )
    
    evidence_label_map <- c(
      physical_overlap = "Physical overlap",
      catalog_gene     = "Catalog",
      eqtl_gene        = "GTEx eQTL",
      nonsyn_gene      = "Non-synonymous",
      ewas_nearby_gene = "EWAS nearby"
    )
    
    binded %>%
      dplyr::left_join(
        gene_base %>%
          dplyr::select(
            gene,
            gene_score,
            top_gwas_hit_score,
            n_supporting_blocks,
            n_clusters,
            n_apps_supported
          ),
        by = "gene"
      ) %>%
      dplyr::mutate(
        evidence_label = dplyr::recode(evidence_type, !!!evidence_label_map, .default = evidence_type),
        value = score_from_blocks_type,
        evidence_label = factor(
          evidence_label,
          levels = c(
            "Physical overlap",
            "Catalog",
            "GTEx eQTL",
            "Non-synonymous",
            "EWAS nearby"
          )
        ),
        hover_txt = paste0(
          "<b>", gene, "</b><br>",
          "Evidence type: ", evidence_label, "<br>",
          "Score contribution: ", formatC(score_from_blocks_type, format = "f", digits = 2), "<br>",
          "Supporting units/blocks for this type: ", n_supporting_blocks_type, "<br>",
          "Total gene score: ", formatC(gene_score, format = "f", digits = 2), "<br>",
          "Top GWAS hit score: ", formatC(top_gwas_hit_score, format = "f", digits = 2), "<br>",
          "Total supporting LD blocks: ", n_supporting_blocks, "<br>",
          "Clusters: ", n_clusters, "<br>",
          "Supporting apps: ", n_apps_supported
        )
      ) %>%
      dplyr::arrange(
        dplyr::desc(gene_score),
        gene,
        evidence_label
      )
  })
  
  output$prioritized_genes_evidence_types_plot <- plotly::renderPlotly({
    df <- prioritized_genes_evidence_types_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized gene evidence-type data available.")
    )
    
    top_n_raw <- input$prioritized_genes_evidence_types_top_n
    
    gene_order <- df %>%
      dplyr::distinct(gene, gene_score) %>%
      dplyr::arrange(dplyr::desc(gene_score), gene)
    
    if (!is.null(top_n_raw) && !identical(top_n_raw, "All")) {
      top_n <- suppressWarnings(as.integer(top_n_raw))
      if (is.finite(top_n) && top_n > 0) {
        gene_order <- utils::head(gene_order, top_n)
      }
    }
    
    df_plot <- df %>%
      dplyr::filter(gene %in% gene_order$gene) %>%
      dplyr::mutate(
        gene = factor(gene, levels = rev(gene_order$gene))
      )
    
    plotly::plot_ly(
      data = df_plot,
      x = ~value,
      y = ~gene,
      color = ~evidence_label,
      type = "bar",
      orientation = "h",
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis = list(
          title = "Canonical support contribution by evidence type",
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE
        ),
        legend = list(
          orientation = "h",
          x = 0,
          y = 1.15
        ),
        margin = list(l = 140, r = 30, t = 80, b = 40)
      )
  })
  
  # modal
  observeEvent(input$show_prioritized_genes_evidence_types_plot, {
    ns <- session$ns
    
    showModal(
      modalDialog(
        title = "Canonical gene evidence-type contributions",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        fluidRow(
          column(
            width = 3,
            
            selectInput(
              inputId = ns("prioritized_genes_evidence_types_top_n"),
              label = "Top genes to display",
              choices = c(10, 20, 30, 50, 75, 100, "All"),
              selected = 30
            ),
            
            tags$div(
              style = "font-size: 12px; color: #666; margin-top: 6px;",
              HTML("GTEx and EWAS contributions are collapsed by LD block. Physical overlap is shown as a visual component.")
            )
          ),
          
          column(
            width = 9,
            plotly::plotlyOutput(
              outputId = ns("prioritized_genes_evidence_types_plot"),
              height = "700px"
            )
          )
        )
      )
    )
  })
  
  
  #############################################################################
  # Plot shared evidence
  #############################################################################
  # ============================================================
  # Shared evidence genes heatmap helpers
  # ============================================================
  
  shared_genes_heatmap_df <- function(df) {
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    app_levels <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(sort(unique(x)), collapse = "; ")
    }
    
    plot_df <- df %>%
      dplyr::mutate(
        gene = as.character(gene),
        gene_label = dplyr::coalesce(as.character(gene_label), as.character(gene)),
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        apps = dplyr::coalesce(as.character(apps), ""),
        evidence_types = dplyr::coalesce(as.character(evidence_types), ""),
        clusters = dplyr::coalesce(as.character(clusters), ""),
        chr_set = dplyr::coalesce(as.character(chr_set), "")
      ) %>%
      dplyr::transmute(
        gene,
        gene_label,
        gene_link_mode,
        apps,
        evidence_types,
        clusters,
        chr_set,
        app_list = strsplit(apps, "\\s*;\\s*")
      ) %>%
      tidyr::unnest_longer(app_list, values_to = "app") %>%
      dplyr::mutate(
        app = trimws(as.character(app))
      ) %>%
      dplyr::filter(!is.na(app), nzchar(app)) %>%
      dplyr::mutate(value = 1L) %>%
      dplyr::distinct(gene, gene_label, gene_link_mode, apps, evidence_types, clusters, chr_set, app, value)
    
    full_grid <- tidyr::expand_grid(
      gene = unique(df$gene),
      app = app_levels
    ) %>%
      dplyr::left_join(
        df %>%
          dplyr::transmute(
            gene = as.character(gene),
            gene_label = dplyr::coalesce(as.character(gene_label), as.character(gene)),
            gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
            apps = dplyr::coalesce(as.character(apps), ""),
            evidence_types = dplyr::coalesce(as.character(evidence_types), ""),
            clusters = dplyr::coalesce(as.character(clusters), ""),
            chr_set = dplyr::coalesce(as.character(chr_set), "")
          ),
        by = "gene"
      ) %>%
      dplyr::left_join(
        plot_df %>% dplyr::select(gene, app, value),
        by = c("gene", "app")
      ) %>%
      dplyr::mutate(
        value = dplyr::coalesce(value, 0L),
        gene_label = dplyr::coalesce(gene_label, gene),
        hover_txt = paste0(
          "Gene: ", gene,
          "<br>Link mode: ", gene_link_mode,
          "<br>App: ", app,
          "<br>Present: ", ifelse(value == 1L, "yes", "no"),
          ifelse(nzchar(evidence_types), paste0("<br>Evidence types: ", evidence_types), ""),
          ifelse(nzchar(apps), paste0("<br>Apps: ", apps), ""),
          ifelse(nzchar(clusters), paste0("<br>Clusters: ", clusters), ""),
          ifelse(nzchar(chr_set), paste0("<br>Chromosomes: ", chr_set), "")
        )
      )
    
    gene_levels <- df %>%
      dplyr::mutate(
        gene_label = dplyr::coalesce(as.character(gene_label), as.character(gene))
      ) %>%
      dplyr::pull(gene_label)
    
    full_grid %>%
      dplyr::mutate(
        gene = factor(gene_label, levels = rev(unique(gene_levels))),
        app = factor(app, levels = app_levels)
      ) %>%
      dplyr::select(gene, app, value, hover_txt)
  }
  
  #############################################################################
  # shared genes plots
  #############################################################################
  
  observeEvent(input$show_shared_genes_heatmap, {
    showModal(
      modalDialog(
        title = "Shared genes · heatmap",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$p(
          "Rows are genes and columns are apps. Green = present, Grey = absent."
        ),
        tags$div(
          class = "smallNote",
          style = "margin-top:6px;",
          HTML("Gene-link mode labels: <b>●</b> position-linked, <b>◆</b> effect-linked, <b>◐</b> position + effect.")
        ),
        uiOutput("shared_genes_heatmap_ui")
      )
    )
  })
  
  output$shared_genes_heatmap_ui <- renderUI({
    df <- shared_genes_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "")
    )
    
    plot_df <- shared_genes_heatmap_df(df)
    n_genes <- dplyr::n_distinct(plot_df$gene)
    
    h <- max(500, min(1600, 18 * n_genes))
    
    plotly::plotlyOutput("shared_genes_heatmap_plot", height = paste0(h, "px"))
  })
  
  output$shared_genes_heatmap_plot <- plotly::renderPlotly({
    df <- shared_genes_df()
    
    shiny::validate(
      shiny::need(is.data.frame(df) && nrow(df) > 0, "No shared evidence genes available.")
    )
    
    plot_df <- shared_genes_heatmap_df(df)
    
    shiny::validate(
      shiny::need(nrow(plot_df) > 0, "No heatmap data available.")
    )
    
    plotly::plot_ly(
      data = plot_df,
      x = ~app,
      y = ~gene,
      z = ~value,
      type = "heatmap",
      text = ~hover_txt,
      hoverinfo = "text",
      colorscale = list(
        c(0.0, "#F2F2F2"),
        c(1.0, "green")
      ),
      zmin = 0,
      zmax = 1,
      showscale = FALSE,
      xgap = 1,
      ygap = 1
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "",
          side = "top"
        ),
        yaxis = list(
          title = "",
          automargin = TRUE
        ),
        margin = list(l = 180, r = 30, t = 55, b = 60)
      )
  })
  
  #############################################################################
  # shared terms plots
  #############################################################################
  
  observeEvent(input$show_traits_heatmap, {
    showModal(
      modalDialog(
        title = "Shared traits · heatmap",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$p(
          "Rows are normalized traits and columns are apps. Orange = present, Grey = absent."
        ),
        plotly::plotlyOutput("shared_traits_heatmap_plot", height = "700px")
      )
    )
  })
  
  output$shared_traits_heatmap_plot <- plotly::renderPlotly({
    df <- shared_traits_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No shared traits available.")
    )
    
    app_levels <- c("catalog", "nonsyn", "ewasdis", "ewastum")
    
    heat_df <- df %>%
      dplyr::mutate(
        term_show = short_term_label(term_label, max_chars = 50)
      ) %>%
      dplyr::mutate(app_list = strsplit(as.character(apps), "\\s*;\\s*")) %>%
      tidyr::unnest_longer(app_list, values_to = "source_app") %>%
      dplyr::filter(!is.na(source_app), nzchar(source_app)) %>%
      dplyr::mutate(source_app = trimws(source_app)) %>%
      dplyr::filter(source_app %in% app_levels) %>%
      dplyr::distinct(term_label, term_show, source_app)
    
    all_terms <- df %>%
      dplyr::mutate(term_show = short_term_label(term_label, max_chars = 50)) %>%
      dplyr::arrange(dplyr::desc(n_apps), dplyr::desc(n_records), term_label) %>%
      dplyr::distinct(term_label, term_show, n_apps, n_records) %>%
      dplyr::slice_head(n = 40)
    
    plot_df <- tidyr::expand_grid(
      term_label = all_terms$term_label,
      source_app = app_levels
    ) %>%
      dplyr::left_join(
        all_terms %>% dplyr::select(term_label, term_show, n_apps, n_records),
        by = "term_label"
      ) %>%
      dplyr::left_join(
        heat_df %>% dplyr::mutate(present = 1L),
        by = c("term_label", "term_show", "source_app")
      ) %>%
      dplyr::mutate(
        present = ifelse(is.na(present), 0L, 1L)
      )
    
    term_order <- all_terms %>%
      dplyr::arrange(dplyr::desc(n_apps), dplyr::desc(n_records), term_label) %>%
      dplyr::pull(term_show)
    
    plot_df <- plot_df %>%
      dplyr::mutate(
        term_show = factor(term_show, levels = rev(unique(term_order))),
        source_app = factor(source_app, levels = app_levels),
        present_label = ifelse(present == 1, "Present", "Absent"),
        hover_txt = paste0(
          "<b>Term:</b> ", htmltools::htmlEscape(term_label),
          "<br><b>App:</b> ", as.character(source_app),
          "<br><b>Status:</b> ", present_label,
          "<br><b>Apps with term:</b> ", n_apps,
          "<br><b>Total records:</b> ", n_records
        )
      )
    
    plotly::plot_ly(
      data = plot_df,
      x = ~source_app,
      y = ~term_show,
      type = "heatmap",
      z = ~present,
      text = ~hover_txt,
      hoverinfo = "text",
      colorscale = list(
        c(0.0, "#F2F2F2"),
        c(0.4999, "#F2F2F2"),
        c(0.5, "orange"),
        c(1.0, "orange")
      ),
      zmin = 0,
      zmax = 1,
      showscale = FALSE,
      xgap = 1,
      ygap = 1
    ) %>%
      plotly::layout(
        # title = list(text = "Shared terms across apps"),
        xaxis = list(
          # title = "",
          side = "top",
          tickangle = 0,
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = app_levels
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = rev(unique(term_order))
        ),
        margin = list(l = 160, r = 20, t = 60, b = 40),
        annotations = list(
          list(
            x = 0,
            y = 1.12,
            xref = "paper",
            yref = "paper",
            #   text = "<b>Legend:</b> grey = Absent, orange = Present",
            showarrow = FALSE,
            xanchor = "left",
            align = "left",
            font = list(size = 12)
          )
        )
      )
  })
  
  #############################################################################
  
  # ============================================================
  # GLOBAL LD BLOCKS PLOT DATA
  # ============================================================
  ld_global_blocks_plot_df <- reactive({
    df <- block_overlap_summary_global_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No global block summary available."),
      need("cluster_id" %in% names(df), "Column cluster_id not found."),
      need("block_id" %in% names(df), "Column block_id not found.")
    )
    
    # ----------------------------
    # helpers
    # ----------------------------
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    extract_block_index <- function(x) {
      x <- as.character(x)
      out <- stringr::str_extract(x, "B[0-9]+$")
      out <- sub("^B", "", out)
      suppressWarnings(as.numeric(out))
    }
    
    get_num_col <- function(df, candidates, default = NA_real_) {
      nm <- candidates[candidates %in% names(df)][1]
      if (length(nm) == 0 || is.na(nm)) {
        return(rep(default, nrow(df)))
      }
      suppressWarnings(as.numeric(df[[nm]]))
    }
    
    get_chr_col <- function(df, candidates, default = "") {
      nm <- candidates[candidates %in% names(df)][1]
      if (length(nm) == 0 || is.na(nm)) {
        return(rep(default, nrow(df)))
      }
      as.character(df[[nm]])
    }
    
    # ----------------------------
    # normalize columns safely
    # ----------------------------
    start_vec <- get_num_col(df, c("start", "pos_ini", "block_start", "start_bp"))
    end_vec   <- get_num_col(df, c("end", "pos_end", "block_end", "end_bp"))
    
    size_kb_vec <- if ("size_kb" %in% names(df)) {
      suppressWarnings(as.numeric(df[["size_kb"]]))
    } else {
      get_num_col(df, c("block_size_kb"))
    }
    
    if (all(!is.finite(size_kb_vec)) && any(is.finite(start_vec)) && any(is.finite(end_vec))) {
      size_kb_vec <- (end_vec - start_vec + 1) / 1000
    }
    
    block_mean_ld_vec <- get_num_col(df, c("block_mean_ld", "block_mean_ld_value", "mean_ld"))
    block_max_ld_vec  <- get_num_col(df, c("block_max_ld", "block_max_ld_value", "max_ld"))
    n_ld_proxy_hits_vec <- get_num_col(df, c("n_ld_proxy_hits"))
    n_genes_vec <- get_num_col(df, c("n_genes", "n_genes_in_block"))
    block_priority_score_vec <- get_num_col(df, c("block_priority_score"))
    support_apps_vec <- get_num_col(df, c("support_apps", "block_support_apps", "n_external_apps"))
    genes_name_vec <- get_chr_col(df, c("genes_name", "genes_in_block"), default = "")
    
    out <- df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id   = as.character(block_id),
        
        chr_num = extract_chr_num(cluster_id),
        block_index = extract_block_index(block_id),
        
        block_label = dplyr::if_else(
          !is.na(block_id) & nzchar(block_id),
          block_id,
          paste0(cluster_id, "_block")
        ),
        
        start = start_vec,
        end = end_vec,
        size_kb = size_kb_vec,
        
        block_mean_ld = block_mean_ld_vec,
        block_max_ld = block_max_ld_vec,
        n_ld_proxy_hits = n_ld_proxy_hits_vec,
        n_genes = n_genes_vec,
        genes_name = genes_name_vec,
        block_priority_score = block_priority_score_vec,
        support_apps = support_apps_vec
      ) %>%
      dplyr::arrange(chr_num, cluster_id, block_index, start, end) %>%
      dplyr::mutate(
        cluster_id = factor(cluster_id, levels = unique(cluster_id)),
        block_label = factor(block_label, levels = unique(block_label))
      )
    
    out
  })
  
  # ============================================================
  # GLOBAL LD BLOCKS STACKED BAR PLOT
  # one horizontal bar per cluster, stacked by blocks
  # ============================================================
  output$ld_global_blocks_plot <- plotly::renderPlotly({
    df <- ld_global_blocks_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No global block summary available.")
    )
    
    plot_df <- df %>%
      dplyr::mutate(
        size_kb = suppressWarnings(as.numeric(size_kb)),
        block_mean_ld = suppressWarnings(as.numeric(block_mean_ld)),
        block_max_ld = suppressWarnings(as.numeric(block_max_ld)),
        n_ld_proxy_hits = suppressWarnings(as.numeric(n_ld_proxy_hits)),
        n_genes = suppressWarnings(as.numeric(n_genes)),
        block_priority_score = suppressWarnings(as.numeric(block_priority_score)),
        support_apps = suppressWarnings(as.numeric(support_apps)),
        
        size_kb = dplyr::if_else(is.finite(size_kb), size_kb, 0),
        
        size_kb_txt = dplyr::if_else(
          is.finite(size_kb),
          formatC(size_kb, format = "f", digits = 1),
          "NA"
        ),
        mean_ld_txt = dplyr::if_else(
          is.finite(block_mean_ld),
          formatC(block_mean_ld, format = "f", digits = 3),
          "NA"
        ),
        max_ld_txt = dplyr::if_else(
          is.finite(block_max_ld),
          formatC(block_max_ld, format = "f", digits = 3),
          "NA"
        ),
        proxy_txt = dplyr::if_else(
          is.finite(n_ld_proxy_hits),
          as.character(as.integer(n_ld_proxy_hits)),
          "NA"
        ),
        genes_n_txt = dplyr::if_else(
          is.finite(n_genes),
          as.character(as.integer(n_genes)),
          "0"
        ),
        score_txt = dplyr::if_else(
          is.finite(block_priority_score),
          formatC(block_priority_score, format = "f", digits = 2),
          "NA"
        ),
        support_apps_txt = dplyr::if_else(
          is.finite(support_apps),
          as.character(as.integer(round(support_apps))),
          "NA"
        ),
        genes_txt = dplyr::if_else(
          !is.na(genes_name) & nzchar(genes_name),
          genes_name,
          "-"
        ),
        hover_txt = paste0(
          "<b>", block_id, "</b><br>",
          "Cluster: ", cluster_id, "<br>",
          "Block size (kb): ", size_kb_txt, "<br>",
          "Mean LD: ", mean_ld_txt, "<br>",
          "Max LD: ", max_ld_txt, "<br>",
          "LD proxy hits: ", proxy_txt, "<br>",
          "Support apps: ", support_apps_txt, "<br>",
          "Genes: ", genes_n_txt, "<br>",
          "Gene list: ", genes_txt, "<br>",
          "Priority score: ", score_txt
        )
      )
    
    validate(
      need(any(plot_df$size_kb > 0, na.rm = TRUE), "No valid block sizes available for plotting.")
    )
    
    vals <- plot_df$block_mean_ld
    vals_ok <- vals[is.finite(vals)]
    
    if (!length(vals_ok)) {
      vals <- rep(0, nrow(plot_df))
      vmin <- 0
      vmax <- 1
    } else {
      vmin <- min(vals_ok, na.rm = TRUE)
      vmax <- max(vals_ok, na.rm = TRUE)
      
      if (!is.finite(vmin)) vmin <- 0
      if (!is.finite(vmax)) vmax <- 1
      if (identical(vmin, vmax)) vmax <- vmin + 1e-6
      
      vals[!is.finite(vals)] <- vmin
    }
    
    norm_vals <- (vals - vmin) / (vmax - vmin)
    norm_vals[!is.finite(norm_vals)] <- 0
    
    plot_df$fill_col <- grDevices::rgb(
      grDevices::colorRamp(c(
        "#d9d9d9",  # gris clar
        "#a6bddb",  # blau grisós
        "#3690c0",  # blau
        "#756bb1",  # violeta
        "#54278f"   # morat fosc
      ))(norm_vals),
      maxColorValue = 255
    )
    
    cluster_levels <- levels(plot_df$cluster_id)
    
    p <- plotly::plot_ly()
    
    for (i in seq_len(nrow(plot_df))) {
      p <- p %>%
        plotly::add_trace(
          data = plot_df[i, , drop = FALSE],
          x = ~size_kb,
          y = ~cluster_id,
          type = "bar",
          orientation = "h",
          name = as.character(plot_df$block_id[i]),
          hovertext = ~hover_txt,
          hoverinfo = "text",
          text = NULL,
          marker = list(
            color = plot_df$fill_col[i],
            line = list(width = 0.4, color = "rgba(80,80,80,0.6)")
          ),
          showlegend = FALSE
        )
    }
    
    p %>%
      plotly::layout(
        title = list(text = "Global LD blocks by cluster"),
        barmode = "stack",
        xaxis = list(
          title = "Total block span (kb)",
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = rev(cluster_levels)
        ),
        margin = list(l = 100, r = 30, t = 50, b = 60),
        hoverlabel = list(align = "left")
      )
  })
  
  output$ld_global_blocks_legend_ui <- renderUI({
    df <- ld_global_blocks_plot_df()
    
    if (!is.data.frame(df) || !nrow(df) || !"block_mean_ld" %in% names(df)) return(NULL)
    
    vals <- suppressWarnings(as.numeric(df$block_mean_ld))
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(NULL)
    
    vmin <- min(vals, na.rm = TRUE)
    vmax <- max(vals, na.rm = TRUE)
    
    tags$div(
      style = "margin-top:8px; font-size:12px; color:#444;",
      tags$span(style = "font-weight:600; margin-right:10px;", "Mean block LD"),
      tags$span(
        style = paste0(
          "display:inline-block; width:220px; height:12px; vertical-align:middle; ",
          "border:1px solid #999; border-radius:6px; margin-right:8px; ",
          "background: linear-gradient(to right, ",
          "#d9d9d9 0%, ",
          "#a6bddb 25%, ",
          "#3690c0 50%, ",
          "#756bb1 75%, ",
          "#54278f 100%);"
        )
      ),
      tags$span(formatC(vmin, format = "f", digits = 2)),
      tags$span(" – "),
      tags$span(formatC(vmax, format = "f", digits = 2))
    )
  })
  
  #############################################################################
  # ============================================================
  # BLOCK PRIORITY COMPONENTS · v2
  # ===========================================================
  
  prioritized_blocks_global_plot_df <- reactive({
    blk <- prioritized_block_df_v2()
    
    validate(
      need(is.data.frame(blk) && nrow(blk) > 0, "No prioritized LD blocks available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    count_semicolon_items <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      if (is.na(x) || !nzchar(x)) return(0L)
      
      vals <- unlist(strsplit(x, "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      
      as.integer(length(unique(vals)))
    }
    
    needed_chr <- c(
      "cluster_id",
      "block_id",
      "gwas_hits",
      "genes_in_block",
      "priority_class",
      "priority_class_relative"
    )
    
    needed_num <- c(
      "block_score",
      "block_score_from_hits",
      "top_gwas_hit_score",
      "other_gwas_hit_score",
      "block_gene_content_bonus",
      "block_bio_bonus",
      "n_gwas_hits",
      "n_apps_supported",
      "n_genes_in_block"
    )
    
    for (nm in needed_chr) {
      if (!nm %in% names(blk)) blk[[nm]] <- ""
    }
    
    for (nm in needed_num) {
      if (!nm %in% names(blk)) blk[[nm]] <- 0
    }
    
    blk_ord <- blk %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        chr_num = extract_chr_num(cluster_id),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$"))),
        block_num = suppressWarnings(as.integer(stringr::str_extract(block_id, "(?<=_B)\\d+$"))),
        
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        block_score_from_hits = dplyr::coalesce(as.numeric(block_score_from_hits), 0),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        other_gwas_hit_score = dplyr::coalesce(as.numeric(other_gwas_hit_score), 0),
        
        block_gene_content_bonus = dplyr::coalesce(
          as.numeric(block_gene_content_bonus),
          as.numeric(block_bio_bonus),
          0
        ),
        
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), ""),
        n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
        priority_class = dplyr::coalesce(as.character(priority_class), ""),
        priority_class_relative = dplyr::coalesce(as.character(priority_class_relative), "")
      ) %>%
      dplyr::mutate(
        n_genes_from_label = vapply(genes_in_block, count_semicolon_items, integer(1)),
        n_genes_in_block = dplyr::if_else(
          n_genes_in_block > 0L,
          n_genes_in_block,
          n_genes_from_label
        ),
        weighted_other_gwas_hit_score = 0.05 * other_gwas_hit_score,
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999),
        block_num = dplyr::coalesce(block_num, 999)
      ) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id, block_num, block_id)
    
    block_levels <- rev(unique(blk_ord$block_id))
    
    block_components_all <- blk_ord %>%
      dplyr::mutate(
        block_id = factor(block_id, levels = block_levels),
        block_id_chr = as.character(block_id),
        
        hover_top = paste0(
          "Block: ", block_id_chr,
          "<br>Cluster: ", cluster_id,
          "<br>Component: Top GWAS-hit support score",
          "<br>Value: ", round(top_gwas_hit_score, 2),
          "<br>Total block score: ", round(block_score, 2),
          "<br>Block score from hits: ", round(block_score_from_hits, 2),
          "<br>Other GWAS-hit support (raw): ", round(other_gwas_hit_score, 2),
          "<br>n GWAS hits: ", n_gwas_hits,
          "<br>n genes: ", n_genes_in_block,
          "<br>n apps: ", n_apps_supported,
          "<br>Priority class (ABS): ", priority_class,
          "<br>Priority class (REL): ", priority_class_relative
        ),
        
        hover_other = paste0(
          "Block: ", block_id_chr,
          "<br>Cluster: ", cluster_id,
          "<br>Component: 0.05 × other GWAS-hit support",
          "<br>Value: ", round(weighted_other_gwas_hit_score, 2),
          "<br>Raw other GWAS-hit support: ", round(other_gwas_hit_score, 2),
          "<br>Total block score: ", round(block_score, 2),
          "<br>Block score from hits: ", round(block_score_from_hits, 2),
          "<br>n GWAS hits: ", n_gwas_hits,
          "<br>n genes: ", n_genes_in_block,
          "<br>n apps: ", n_apps_supported,
          "<br>Priority class (ABS): ", priority_class,
          "<br>Priority class (REL): ", priority_class_relative
        ),
        
        hover_gene = paste0(
          "Block: ", block_id_chr,
          "<br>Cluster: ", cluster_id,
          "<br>Component: Gene-content bonus",
          "<br>Value: ", round(block_gene_content_bonus, 2),
          "<br>Total block score: ", round(block_score, 2),
          "<br>Block score from hits: ", round(block_score_from_hits, 2),
          "<br>n GWAS hits: ", n_gwas_hits,
          "<br>n genes: ", n_genes_in_block,
          "<br>n apps: ", n_apps_supported,
          "<br>Genes in block: ", ifelse(nzchar(genes_in_block), genes_in_block, "-"),
          "<br>Priority class (ABS): ", priority_class,
          "<br>Priority class (REL): ", priority_class_relative
        )
      ) %>%
      dplyr::transmute(
        block_id,
        block_id_chr,
        cluster_id,
        component_top = top_gwas_hit_score,
        component_other = weighted_other_gwas_hit_score,
        component_gene = block_gene_content_bonus,
        hover_top,
        hover_other,
        hover_gene
      ) %>%
      tidyr::pivot_longer(
        cols = c(component_top, component_other, component_gene),
        names_to = "component",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        component = dplyr::recode(
          component,
          component_top = "Top GWAS-hit support",
          component_other = "Other GWAS-hit support (x0.05)",
          component_gene = "Gene-content bonus"
        ),
        hover_txt = dplyr::case_when(
          component == "Top GWAS-hit support" ~ hover_top,
          component == "Other GWAS-hit support (x0.05)" ~ hover_other,
          TRUE ~ hover_gene
        )
      ) %>%
      dplyr::select(block_id, block_id_chr, cluster_id, component, value, hover_txt)
    
    block_components_all
  })
  
  output$prioritized_blocks_global_plot <- plotly::renderPlotly({
    df <- prioritized_blocks_global_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No block priority component data available.")
    )
    
    plotly::plot_ly(
      data = df,
      y = ~block_id,
      x = ~value,
      color = ~component,
      colors = c(
        "Top GWAS-hit support"           = "#1f78b4",
        "Other GWAS-hit support (x0.05)" = "#6baed6",
        "Gene-content bonus"             = "#33a02c"
      ),
      type = "bar",
      orientation = "h",
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis = list(
          title = "Block priority score",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Block",
          categoryorder = "array",
          categoryarray = levels(df$block_id),
          automargin = TRUE
        ),
        legend = list(
          title = list(text = "Score component"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 120, r = 20, t = 20, b = 60)
      )
  })
  
  observeEvent(input$show_block_priority_plot, {
    df <- prioritized_block_df_v2()
    req(is.data.frame(df), nrow(df) > 0)
    
    showModal(
      modalDialog(
        title = "LD block priority score components",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$p(
          "Each horizontal bar represents one prioritized LD block. ",
          "Bar segments show the contribution of the current canonical block-score components: ",
          "the top GWAS-hit support score in the block, the downweighted contribution of the remaining GWAS-hit support ",
          "(0.05 × other support), and the additional bonus for gene content in the block. "
        ),
        tags$p(
          "This block score is now aligned with the canonical prioritized gene score logic: ",
          "the strongest support unit dominates, while additional support contributes with a reduced weight."
        ),
        plotly::plotlyOutput("prioritized_blocks_global_plot", height = "650px")
      )
    )
  })
  
  #############################################################################
  # ============================================================
  # PRIORITIZED GENES BY CLUSTER · v2
  # base completa de clusters + overlay de gens prioritzats
  # ============================================================
  prioritized_genes_cluster_plot_df <- reactive({
    cl0 <- priority_source_clusters()
    pg  <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(cl0) && nrow(cl0) > 0, "No integrated clusters available."),
      need(is.data.frame(pg) && nrow(pg) > 0, "No prioritized gene table available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    collapse_unique_semicolon <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("")
      paste(sort(unique(x)), collapse = "; ")
    }
    
    cl_std <- cl0 %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        chr = as.character(chr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
      ) %>%
      dplyr::mutate(
        chr_num = extract_chr_num(cluster_id),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$"))),
        cluster_size_bp = dplyr::if_else(
          is.finite(start) & is.finite(end),
          end - start + 1,
          NA_real_
        ),
        cluster_size_kb = cluster_size_bp / 1000
      ) %>%
      dplyr::mutate(
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999),
        cluster_size_bp = dplyr::coalesce(cluster_size_bp, 0),
        cluster_size_kb = dplyr::coalesce(cluster_size_kb, 0)
      ) %>%
      dplyr::distinct(cluster_id, .keep_all = TRUE)
    
    gene_by_cluster <- pg %>%
      dplyr::mutate(
        gene = as.character(gene),
        cluster_ids = dplyr::coalesce(as.character(cluster_ids), ""),
        gene_score = suppressWarnings(as.numeric(gene_score))
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      tidyr::separate_rows(cluster_ids, sep = "\\s*;\\s*") %>%
      dplyr::mutate(
        cluster_id = trimws(cluster_ids)
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        n_genes = dplyr::n_distinct(gene),
        mean_gene_score = mean(gene_score, na.rm = TRUE),
        best_gene_score = max(gene_score, na.rm = TRUE),
        genes_list = collapse_unique_semicolon(gene),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        mean_gene_score = dplyr::if_else(is.finite(mean_gene_score), mean_gene_score, 0),
        best_gene_score = dplyr::if_else(is.finite(best_gene_score), best_gene_score, 0)
      )
    
    out <- cl_std %>%
      dplyr::left_join(gene_by_cluster, by = "cluster_id") %>%
      dplyr::mutate(
        n_genes = dplyr::coalesce(n_genes, 0L),
        mean_gene_score = dplyr::coalesce(mean_gene_score, 0),
        best_gene_score = dplyr::coalesce(best_gene_score, 0),
        genes_list = dplyr::coalesce(genes_list, "")
      )
    
    cluster_levels <- out %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
      dplyr::pull(cluster_id) %>%
      rev()
    
    vals <- out$mean_gene_score
    has_score <- any(is.finite(vals))
    
    if (has_score) {
      rng <- range(vals, na.rm = TRUE)
      
      if (!all(is.finite(rng))) {
        out$fill_col <- "#d9d9d9"
      } else if (rng[1] == rng[2]) {
        out$fill_col <- ifelse(out$n_genes > 0, "#993404", "#d9d9d9")
      } else {
        pal_fun <- scales::col_numeric(
          palette = c("#d9d9d9", "#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
          domain = rng,
          na.color = "#d9d9d9"
        )
        out$fill_col <- pal_fun(out$mean_gene_score)
        out$fill_col[out$n_genes == 0] <- "#d9d9d9"
      }
    } else {
      out$fill_col <- "#d9d9d9"
    }
    
    out %>%
      dplyr::mutate(
        cluster_id_chr = cluster_id,
        cluster_id = factor(cluster_id, levels = cluster_levels),
        hover_txt = paste0(
          "Cluster: ", cluster_id_chr,
          "<br>Genes in cluster: ", n_genes,
          "<br>Mean gene score: ", round(mean_gene_score, 2),
          "<br>Best gene score: ", round(best_gene_score, 2),
          "<br>Cluster size (kb): ", round(cluster_size_kb, 2),
          ifelse(
            nzchar(genes_list),
            paste0("<br>Genes: ", genes_list),
            ""
          )
        )
      )
  })
  
  # ============================================================
  # LEGEND · prioritized genes by cluster
  # ============================================================
  prioritized_genes_cluster_legend_ui <- function(df) {
    if (!is.data.frame(df) || !nrow(df) || !"mean_gene_score" %in% names(df)) {
      return(NULL)
    }
    
    vals <- suppressWarnings(as.numeric(df$mean_gene_score))
    vals <- vals[is.finite(vals)]
    
    if (!length(vals)) {
      return(
        tags$div(
          style = "margin-top:8px; font-size:12px; color:#444;",
          tags$span(style = "font-weight:600; margin-right:10px;", "Mean gene score"),
          tags$span(style = "color:#777;", "No score values available")
        )
      )
    }
    
    vmin <- min(vals, na.rm = TRUE)
    vmax <- max(vals, na.rm = TRUE)
    
    tags$div(
      style = "margin-top:8px; font-size:12px; color:#444;",
      tags$span(style = "font-weight:600; margin-right:10px;", "Mean gene score"),
      tags$span(
        style = paste0(
          "display:inline-block; width:220px; height:12px; vertical-align:middle; ",
          "border:1px solid #999; border-radius:6px; margin-right:8px; ",
          "background: linear-gradient(to right, ",
          "#d9d9d9 0%, ",
          "#fff7bc 20%, ",
          "#fec44f 40%, ",
          "#fe9929 60%, ",
          "#d95f0e 80%, ",
          "#993404 100%);"
        )
      ),
      tags$span(formatC(vmin, format = "f", digits = 2)),
      tags$span(" – "),
      tags$span(formatC(vmax, format = "f", digits = 2))
    )
  }
  
  output$prioritized_genes_cluster_legend <- renderUI({
    df <- prioritized_genes_cluster_plot_df()
    prioritized_genes_cluster_legend_ui(df)
  })
  
  # ============================================================
  # PLOT · prioritized genes by cluster
  # X = number of prioritized genes
  # color = mean gene score
  # ============================================================
  output$prioritized_genes_cluster_plot <- plotly::renderPlotly({
    df <- prioritized_genes_cluster_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized genes plot data available.")
    )
    
    plotly::plot_ly(
      data = df,
      source = "prioritized_genes_cluster_plot",
      y = ~cluster_id,
      x = ~n_genes,
      type = "bar",
      orientation = "h",
      customdata = ~cluster_id_chr,
      hovertext = ~hover_txt,
      hoverinfo = "text",
      marker = list(
        color = df$fill_col,
        line = list(color = "black", width = 0.5)
      ),
      showlegend = FALSE
    ) %>%
      plotly::layout(
        clickmode = "event+select",
        xaxis = list(
          title = "Number of prioritized genes",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Cluster",
          categoryorder = "array",
          categoryarray = levels(df$cluster_id),
          automargin = TRUE
        ),
        margin = list(l = 90, r = 20, t = 20, b = 60)
      )
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "prioritized_genes_cluster_plot"), {
    ed <- plotly::event_data("plotly_click", source = "prioritized_genes_cluster_plot")
    req(ed)
    
    cid <- ed$customdata[[1]] %||% ""
    req(nzchar(cid))
    
    session$sendCustomMessage(
      "filter_gene_table_by_cluster",
      list(cluster_id = cid)
    )
  })
  
# observeEvent(input$clear_prioritized_genes_table_filter, {
#   session$sendCustomMessage(
#     "filter_gene_table_by_cluster",
#     list(cluster_id = "")
#   )
# })
  # ============================================================
  
  consensus_cluster_app_plot_df <- reactive({
    cl <- clusters_consensus()
    ca <- candidates_consensus()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No consensus clusters available."),
      need(is.data.frame(ca) && nrow(ca) > 0, "No consensus candidates available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    all_apps <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
    
    cl_df <- cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        chr = as.character(chr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end)),
        cluster_size_kb = pmax((end - start + 1) / 1000, 0)
      ) %>%
      dplyr::distinct(cluster_id, .keep_all = TRUE)
    
    ca_df <- ca %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        source_app = tolower(trimws(as.character(source_app))),
        classe = tolower(trimws(as.character(classe))),
        chr = as.character(chr),
        position = suppressWarnings(as.integer(position)),
        rsid = dplyr::coalesce(dplyr::na_if(as.character(rsid), ""), NA_character_),
        id_hit = dplyr::coalesce(dplyr::na_if(as.character(id_hit), ""), NA_character_),
        hit_key = dplyr::coalesce(
          id_hit,
          rsid,
          dplyr::if_else(
            !is.na(chr) & nzchar(chr) & !is.na(position),
            paste0(chr, ":", position),
            NA_character_
          )
        ),
        source_app = dplyr::case_when(
          source_app %in% c("catalog") ~ "catalog",
          source_app %in% c("gtex") ~ "gtex",
          source_app %in% c("nonsyn", "dbnsfp") ~ "nonsyn",
          source_app %in% c("ewasdis", "disease") ~ "ewasdis",
          source_app %in% c("ewastum", "tumor") ~ "ewastum",
          TRUE ~ source_app
        )
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        source_app %in% all_apps,
        !is.na(hit_key), nzchar(hit_key),
        (
          (source_app == "catalog"  & classe == "catalog_hit")  |
            (source_app == "gtex"     & classe == "gtex_hit")     |
            (source_app == "nonsyn"   & classe == "nonsyn_hit")   |
            (source_app == "ewasdis"  & classe == "ewasdis_hit")  |
            (source_app == "ewastum"  & classe == "ewastum_hit")
        )
      ) %>%
      dplyr::distinct(cluster_id, source_app, hit_key, .keep_all = TRUE) %>%
      dplyr::count(cluster_id, source_app, name = "n_hits")
    
    plot_df <- ca_df %>%
      tidyr::complete(
        cluster_id,
        source_app = all_apps,
        fill = list(n_hits = 0L)
      ) %>%
      dplyr::left_join(cl_df, by = "cluster_id") %>%
      dplyr::filter(!is.na(cluster_size_kb)) %>%
      dplyr::group_by(cluster_id, chr, start, end, cluster_size_kb) %>%
      dplyr::mutate(
        total_hits = sum(n_hits, na.rm = TRUE),
        prop = dplyr::if_else(total_hits > 0, n_hits / total_hits, 0),
        seg_kb = cluster_size_kb * prop
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(total_hits > 0)
    
    ord <- plot_df %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::mutate(
        chr_num = extract_chr_num(cluster_id),
        cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
      ) %>%
      dplyr::mutate(
        chr_num = dplyr::coalesce(chr_num, 999),
        cluster_num = dplyr::coalesce(cluster_num, 999)
      ) %>%
      dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
      dplyr::pull(cluster_id) %>%
      rev()
    
    plot_df %>%
      dplyr::mutate(
        cluster_id_chr = cluster_id,
        cluster_id = factor(cluster_id, levels = ord),
        source_app = factor(source_app, levels = all_apps),
        hover_txt = paste0(
          "Cluster: ", cluster_id_chr,
          "<br>App: ", as.character(source_app),
          "<br>Cluster size (kb): ", round(cluster_size_kb, 2),
          "<br>Hits in app: ", n_hits,
          "<br>Total hits: ", total_hits,
          "<br>Proportion: ", sprintf("%.1f%%", 100 * prop),
          "<br>Segment size (kb): ", round(seg_kb, 2)
        )
      )
  })
  
  output$consensus_candidates_apps_plot <- plotly::renderPlotly({
    df <- consensus_cluster_app_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster composition data available.")
    )
    
    
    color_map <- c(
      catalog = "#e6ab02",
      gtex = "#f28e2b",
      nonsyn = "#e15759",
      ewasdis = "#b24745",
      ewastum = "#7b3f98"
    )
    app_order <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
    
    p <- plotly::plot_ly(source = "consensus_cluster_app_plot")
    
    for (app in app_order) {
      dsub <- df %>% dplyr::filter(as.character(source_app) == app)
      if (!nrow(dsub)) next
      
      p <- p %>%
        plotly::add_trace(
          data = dsub,
          y = ~cluster_id,
          x = ~seg_kb,
          type = "bar",
          orientation = "h",
          name = app,
          customdata = ~cluster_id_chr,
          hovertext = ~hover_txt,
          hoverinfo = "text",
          showlegend = FALSE,
          marker = list(
            color = unname(color_map[app]),
            line = list(color = "white", width = 0.6)
          )
        )
    }
    
    p %>%
      plotly::layout(
        barmode = "stack",
        clickmode = "event+select",
        showlegend = FALSE,
        xaxis = list(
          title = "Cluster size (kb)",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Cluster",
          categoryorder = "array",
          categoryarray = levels(df$cluster_id),
          automargin = TRUE
        ),
        margin = list(l = 90, r = 20, t = 20, b = 60)
      )
  })
  
  output$consensus_candidates_apps_legend <- renderUI({
    consensus_candidates_apps_legend_ui()
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "consensus_cluster_app_plot"), {
    ed <- plotly::event_data("plotly_click", source = "consensus_cluster_app_plot")
    req(ed)
    
    cid <- ed$customdata[[1]] %||% ""
    req(nzchar(cid))
    
    session$sendCustomMessage(
      "filter_candidates_consensus_table_by_cluster",
      list(cluster_id = cid)
    )
  })
  
  consensus_candidates_apps_table_df <- reactive({
    df <- consensus_cluster_app_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No consensus candidate app counts available.")
    )
    
    nm <- names(df)
    
    cluster_col <- if ("cluster_id" %in% nm) "cluster_id" else NA_character_
    
    app_col <- dplyr::case_when(
      "app" %in% nm ~ "app",
      "source_app" %in% nm ~ "source_app",
      "source" %in% nm ~ "source",
      "classe" %in% nm ~ "classe",
      TRUE ~ NA_character_
    )
    
    count_col <- dplyr::case_when(
      "n" %in% nm ~ "n",
      "count" %in% nm ~ "count",
      "value" %in% nm ~ "value",
      "hits" %in% nm ~ "hits",
      "n_hits" %in% nm ~ "n_hits",
      TRUE ~ NA_character_
    )
    
    validate(
      need(!is.na(cluster_col), paste("No cluster_id column found. Available columns:", paste(nm, collapse = ", "))),
      need(!is.na(app_col), paste("No app column found. Available columns:", paste(nm, collapse = ", "))),
      need(!is.na(count_col), paste("No count column found. Available columns:", paste(nm, collapse = ", ")))
    )
    
    out <- df %>%
      dplyr::transmute(
        cluster_id = as.character(.data[[cluster_col]]),
        app        = as.character(.data[[app_col]]),
        n          = suppressWarnings(as.integer(.data[[count_col]]))
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(app), nzchar(app)) %>%
      dplyr::group_by(cluster_id, app) %>%
      dplyr::summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from = app,
        values_from = n,
        values_fill = 0
      )
    
    wanted_apps <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
    for (aa in wanted_apps) {
      if (!aa %in% names(out)) out[[aa]] <- 0L
    }
    
    cand_all <- tryCatch(integrated_candidates_r(), error = function(e) NULL)
    
    gwas_counts <- if (is.data.frame(cand_all) && nrow(cand_all) > 0) {
      
      hit_col <- dplyr::case_when(
        "id_hit" %in% names(cand_all) ~ "id_hit",
        "rsid" %in% names(cand_all) ~ "rsid",
        TRUE ~ NA_character_
      )
      
      validate(
        need("cluster_id" %in% names(cand_all), "No cluster_id column found in integrated_candidates_from_rds()."),
        need("classe" %in% names(cand_all), "No classe column found in integrated_candidates_from_rds()."),
        need(!is.na(hit_col), "Neither id_hit nor rsid found in integrated_candidates_from_rds().")
      )
      
      gwas_only <- cand_all %>%
        dplyr::transmute(
          cluster_id = trimws(as.character(cluster_id)),
          classe2    = trimws(as.character(classe)),
          hit2       = trimws(as.character(.data[[hit_col]]))
        ) %>%
        dplyr::filter(
          classe2 == "GWAS",
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(hit2), nzchar(hit2)
        ) %>%
        dplyr::distinct(cluster_id, hit2)
      
      
      gwas_only %>%
        dplyr::count(cluster_id, name = "gwas")
      
    } else {
      tibble::tibble(cluster_id = character(), gwas = integer())
    }
    
    out <- out %>%
      dplyr::left_join(gwas_counts, by = "cluster_id")
    
    out$gwas[is.na(out$gwas)] <- 0L
    
    out %>%
      dplyr::select(cluster_id, gwas, dplyr::all_of(wanted_apps)) %>%
      dplyr::arrange(cluster_id)
    
  })
  
  output$consensus_candidates_apps_table_dt <- DT::renderDT({
    df <- consensus_candidates_apps_table_df()
    
    DT::datatable(
      df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = TRUE)
  
  ###################################################################################
  ############################# GENE NETWORKS #######################################
  ###################################################################################
  
  clean_top_genes_for_network <- function(df, top_n = 50, keep_noncoding = FALSE) {
    stopifnot(is.data.frame(df))
    
    out <- df %>%
      dplyr::filter(!is.na(gene), nzchar(trimws(gene))) %>%
      dplyr::mutate(
        gene = trimws(as.character(gene)),
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0)
      ) %>%
      # fora cadenes compostes / múltiples
      dplyr::filter(!grepl("[,;]", gene)) %>%
      dplyr::filter(!grepl("\\s-\\s", gene)) %>%
      dplyr::filter(!grepl("/", gene)) %>%
      dplyr::filter(!grepl("^LOC[0-9]+$", gene, ignore.case = TRUE))
    
    if (!isTRUE(keep_noncoding)) {
      out <- out %>%
        dplyr::filter(!grepl("^MIR[0-9]", gene, ignore.case = TRUE)) %>%
        dplyr::filter(!grepl("^LINC", gene, ignore.case = TRUE)) %>%
        dplyr::filter(!grepl("^RP[0-9]+-", gene, ignore.case = TRUE)) %>%
        dplyr::filter(!grepl("^AC[0-9]", gene, ignore.case = TRUE)) %>%
        dplyr::filter(!grepl("^AL[0-9]", gene, ignore.case = TRUE)) %>%
        dplyr::filter(!grepl("^AP[0-9]", gene, ignore.case = TRUE))
    }
    
    out %>%
      dplyr::arrange(dplyr::desc(gene_score), gene) %>%
      dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::slice_head(n = top_n)
  }
  
  top_genes_for_network <- reactive({
    df <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized genes available.")
    )
    
    df <- df %>%
      dplyr::mutate(
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0)
      ) %>%
      dplyr::filter(gene_score > 0)
    
    validate(
      need(nrow(df) > 0, "No prioritized genes with positive score available.")
    )
    
    clean_top_genes_for_network(df, top_n = 50, keep_noncoding = FALSE)
  })
  
  gene_pathway_network_obj <- reactive({
    tg <- top_genes_for_network()
    genes <- tg$gene
    
    validate(
      need(length(genes) >= 2, "Not enough clean prioritized genes for pathway network.")
    )
    
    entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    )
    
    entrez <- entrez[!is.na(entrez)]
    
    validate(
      need(length(entrez) >= 2, "Could not map enough genes to ENTREZID.")
    )
    
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = unname(entrez),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.2,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    
    validate(
      need(!is.null(ego), "GO enrichment failed."),
      need(nrow(as.data.frame(ego)) > 0, "No enriched GO pathways found.")
    )
    
    ego_df <- as.data.frame(ego) %>%
      dplyr::arrange(p.adjust, dplyr::desc(Count)) %>%
      dplyr::slice_head(n = 15)
    
    edges <- ego_df %>%
      dplyr::select(ID, Description, geneID, p.adjust, Count) %>%
      tidyr::separate_rows(geneID, sep = "/") %>%
      dplyr::transmute(
        from = geneID,
        to = Description,
        edge_type = "gene_pathway",
        p_adjust = p.adjust,
        pathway_count = Count
      )
    
    gene_nodes <- tg %>%
      dplyr::filter(gene %in% unique(edges$from)) %>%
      dplyr::transmute(
        id = gene,
        label = gene,
        node_type = "gene",
        score = gene_score,
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), "")
      )
    
    pathway_nodes <- ego_df %>%
      dplyr::transmute(
        id = Description,
        label = Description,
        node_type = "pathway",
        score = -log10(p.adjust),
        gene_link_mode = "",
        gene_evidence_types = "",
        gene_source_apps = "",
        apps_supported = ""
      )
    
    list(
      nodes = dplyr::bind_rows(gene_nodes, pathway_nodes),
      edges = edges
    )
  })
  
  gene_ppi_network_obj <- reactive({
    validate(
      need(
        requireNamespace("STRINGdb", quietly = TRUE),
        "STRINGdb is not available right now. Install it later when Bioconductor is reachable."
      )
    )
    
    tg <- top_genes_for_network()
    genes <- tg$gene
    
    validate(
      need(length(genes) >= 2, "Not enough clean prioritized genes for PPI network.")
    )
    
    string_db <- STRINGdb::STRINGdb$new(
      version = "12",
      species = 9606,
      score_threshold = 400,
      input_directory = ""
    )
    
    map_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
    
    mapped <- tryCatch(
      string_db$map(map_df, "gene", removeUnmappedRows = TRUE),
      error = function(e) NULL
    )
    
    validate(
      need(!is.null(mapped), "STRING mapping failed."),
      need(nrow(mapped) >= 2, "Too few genes mapped to STRING.")
    )
    
    ppi <- tryCatch(
      string_db$get_interactions(mapped$STRING_id),
      error = function(e) NULL
    )
    
    validate(
      need(!is.null(ppi), "STRING interaction retrieval failed."),
      need(nrow(ppi) > 0, "No PPI edges found.")
    )
    
    id2gene <- mapped %>%
      dplyr::select(STRING_id, gene)
    
    edges <- ppi %>%
      dplyr::left_join(id2gene, by = c("from" = "STRING_id")) %>%
      dplyr::rename(from_gene = gene) %>%
      dplyr::left_join(id2gene, by = c("to" = "STRING_id")) %>%
      dplyr::rename(to_gene = gene) %>%
      dplyr::filter(!is.na(from_gene), !is.na(to_gene), from_gene != to_gene) %>%
      dplyr::transmute(
        from = from_gene,
        to = to_gene
      ) %>%
      dplyr::distinct()
    
    validate(
      need(nrow(edges) > 0, "No mapped PPI edges available.")
    )
    
    nodes <- tg %>%
      dplyr::filter(gene %in% unique(c(edges$from, edges$to))) %>%
      dplyr::transmute(
        id = gene,
        label = gene,
        node_type = "gene",
        score = gene_score,
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), "")
      )
    
    list(
      nodes = nodes,
      edges = edges
    )
  })
  
  # ============================================================
  # HELPERS FOR LABEL DISPLAY
  # ============================================================
  
  get_label_subset <- function(df, mode = c("None", "Top", "All"), top_n = 10, order_col = NULL) {
    mode <- match.arg(mode)
    
    if (!is.data.frame(df) || !nrow(df)) {
      return(df[0, , drop = FALSE])
    }
    
    out <- df
    
    if (!is.null(order_col) && order_col %in% names(out)) {
      out <- out %>% dplyr::arrange(dplyr::desc(.data[[order_col]]))
    }
    
    if (identical(mode, "None")) {
      return(out[0, , drop = FALSE])
    }
    
    if (identical(mode, "All")) {
      return(out)
    }
    
    out %>% dplyr::slice_head(n = top_n)
  }
  
  # ============================================================
  # PATHWAY NETWORK PLOT
  # ============================================================
  
  output$gene_pathway_network <- plotly::renderPlotly({
    nw <- gene_pathway_network_obj()
    
    nodes <- nw$nodes
    edges <- nw$edges
    
    validate(
      need(is.data.frame(nodes) && nrow(nodes) > 0, "No pathway nodes available."),
      need(is.data.frame(edges) && nrow(edges) > 0, "No pathway edges available.")
    )
    
    g <- igraph::graph_from_data_frame(
      d = edges %>% dplyr::select(from, to),
      vertices = nodes,
      directed = FALSE
    )
    
    lay <- igraph::layout_with_fr(g)
    
    vdf <- data.frame(
      id = igraph::V(g)$name,
      x = lay[, 1],
      y = lay[, 2],
      label = igraph::V(g)$label,
      node_type = igraph::V(g)$node_type,
      score = igraph::V(g)$score,
      gene_link_mode = if ("gene_link_mode" %in% names(nodes)) igraph::V(g)$gene_link_mode else NA_character_,
      gene_evidence_types = if ("gene_evidence_types" %in% names(nodes)) igraph::V(g)$gene_evidence_types else NA_character_,
      gene_source_apps = if ("gene_source_apps" %in% names(nodes)) igraph::V(g)$gene_source_apps else NA_character_,
      apps_supported = if ("apps_supported" %in% names(nodes)) igraph::V(g)$apps_supported else NA_character_,
      stringsAsFactors = FALSE
    )
    
    vdf <- vdf %>%
      dplyr::mutate(
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        link_class = dplyr::case_when(
          node_type != "gene" ~ "pathway",
          grepl("physical", gene_link_mode) & grepl("effect", gene_link_mode) ~ "position+effect",
          grepl("physical", gene_link_mode) ~ "position",
          grepl("effect", gene_link_mode) ~ "effect",
          TRUE ~ "gene"
        ),
        marker_symbol = dplyr::case_when(
          link_class == "position" ~ "circle",
          link_class == "effect" ~ "diamond",
          link_class == "position+effect" ~ "circle-open",
          TRUE ~ "circle"
        ),
        hover_txt = dplyr::case_when(
          node_type == "gene" ~ paste0(
            label,
            "<br>gene_score = ", round(score, 2),
            ifelse(nzchar(gene_link_mode), paste0("<br>link mode = ", gene_link_mode), ""),
            ifelse(nzchar(gene_evidence_types), paste0("<br>evidence = ", gene_evidence_types), ""),
            ifelse(nzchar(gene_source_apps), paste0("<br>gene apps = ", gene_source_apps), ""),
            ifelse(nzchar(apps_supported), paste0("<br>hit-supported apps = ", apps_supported), "")
          ),
          TRUE ~ paste0(
            label,
            "<br>-log10(adj.p) = ", round(score, 2)
          )
        )
      )
    
    edf <- igraph::as_data_frame(g, what = "edges")
    edf$x <- vdf$x[match(edf$from, vdf$id)]
    edf$y <- vdf$y[match(edf$from, vdf$id)]
    edf$xend <- vdf$x[match(edf$to, vdf$id)]
    edf$yend <- vdf$y[match(edf$to, vdf$id)]
    
    gene_df <- vdf %>%
      dplyr::filter(node_type == "gene") %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(label_plot = label)
    
    path_df <- vdf %>%
      dplyr::filter(node_type == "pathway") %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(label_plot = stringr::str_trunc(label, 40))
    
    pathway_label_mode <- input$pathway_label_mode %||% "Top"
    n_labels <- input$network_top_labels %||% 10
    
    gene_lab_df <- get_label_subset(
      gene_df,
      mode = pathway_label_mode,
      top_n = n_labels,
      order_col = "score"
    )
    
    path_lab_df <- get_label_subset(
      path_df,
      mode = pathway_label_mode,
      top_n = n_labels,
      order_col = "score"
    )
    
    p <- plotly::plot_ly()
    
    for (i in seq_len(nrow(edf))) {
      p <- p %>%
        plotly::add_segments(
          x = edf$x[i], y = edf$y[i],
          xend = edf$xend[i], yend = edf$yend[i],
          inherit = FALSE,
          line = list(width = 0.6, color = "gray"),
          hoverinfo = "none",
          showlegend = FALSE
        )
    }
    
    gene_position_df <- gene_df %>% dplyr::filter(link_class == "position")
    gene_effect_df <- gene_df %>% dplyr::filter(link_class == "effect")
    gene_both_df <- gene_df %>% dplyr::filter(link_class == "position+effect")
    gene_other_df <- gene_df %>% dplyr::filter(!link_class %in% c("position", "effect", "position+effect"))
    
    if (nrow(gene_position_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_position_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: position-linked",
          marker = list(size = 12, symbol = "circle")
        )
    }
    
    if (nrow(gene_effect_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_effect_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: effect-linked",
          marker = list(size = 12, symbol = "diamond")
        )
    }
    
    if (nrow(gene_both_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_both_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: position + effect",
          marker = list(size = 12, symbol = "circle-open")
        )
    }
    
    if (nrow(gene_other_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_other_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes",
          marker = list(size = 12, symbol = "circle")
        )
    }
    
    p <- p %>%
      plotly::add_markers(
        data = path_df,
        x = ~x, y = ~y,
        text = ~hover_txt,
        hoverinfo = "text",
        name = "Pathways",
        marker = list(size = 18, symbol = "square")
      )
    
    if (nrow(gene_lab_df) > 0) {
      p <- p %>%
        plotly::add_text(
          data = gene_lab_df,
          x = ~x, y = ~y,
          text = ~label_plot,
          textposition = "top center",
          showlegend = FALSE,
          hoverinfo = "none"
        )
    }
    
    if (nrow(path_lab_df) > 0) {
      p <- p %>%
        plotly::add_text(
          data = path_lab_df,
          x = ~x, y = ~y,
          text = ~label_plot,
          textposition = "middle right",
          showlegend = FALSE,
          hoverinfo = "none"
        )
    }
    
    p %>%
      plotly::layout(
        title = "Prioritized gene–pathway network",
        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        legend = list(orientation = "h")
      )
  })
  
  # ============================================================
  # PPI NETWORK PLOT
  # ============================================================
  
  output$gene_ppi_network <- plotly::renderPlotly({
    validate(
      need(requireNamespace("STRINGdb", quietly = TRUE),
           "STRINGdb is not available right now. Install it later when Bioconductor is reachable.")
    )
    
    nw <- gene_ppi_network_obj()
    
    validate(
      need(is.data.frame(nw$nodes) && nrow(nw$nodes) > 0, "No PPI nodes available."),
      need(is.data.frame(nw$edges) && nrow(nw$edges) > 0, "No PPI edges available.")
    )
    
    g <- igraph::graph_from_data_frame(
      d = nw$edges,
      vertices = nw$nodes,
      directed = FALSE
    )
    
    lay <- igraph::layout_with_fr(g)
    
    vdf <- data.frame(
      id = igraph::V(g)$name,
      x = lay[, 1],
      y = lay[, 2],
      label = igraph::V(g)$label,
      score = igraph::V(g)$score,
      degree = igraph::degree(g),
      gene_link_mode = if ("gene_link_mode" %in% names(nw$nodes)) igraph::V(g)$gene_link_mode else NA_character_,
      gene_evidence_types = if ("gene_evidence_types" %in% names(nw$nodes)) igraph::V(g)$gene_evidence_types else NA_character_,
      gene_source_apps = if ("gene_source_apps" %in% names(nw$nodes)) igraph::V(g)$gene_source_apps else NA_character_,
      apps_supported = if ("apps_supported" %in% names(nw$nodes)) igraph::V(g)$apps_supported else NA_character_,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        gene_link_mode = dplyr::coalesce(as.character(gene_link_mode), ""),
        gene_evidence_types = dplyr::coalesce(as.character(gene_evidence_types), ""),
        gene_source_apps = dplyr::coalesce(as.character(gene_source_apps), ""),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        link_class = dplyr::case_when(
          grepl("physical", gene_link_mode) & grepl("effect", gene_link_mode) ~ "position+effect",
          grepl("physical", gene_link_mode) ~ "position",
          grepl("effect", gene_link_mode) ~ "effect",
          TRUE ~ "gene"
        ),
        marker_symbol = dplyr::case_when(
          link_class == "position" ~ "circle",
          link_class == "effect" ~ "diamond",
          link_class == "position+effect" ~ "circle-open",
          TRUE ~ "circle"
        ),
        hover_txt = paste0(
          label,
          "<br>gene_score = ", round(score, 2),
          "<br>degree = ", degree,
          ifelse(nzchar(gene_link_mode), paste0("<br>link mode = ", gene_link_mode), ""),
          ifelse(nzchar(gene_evidence_types), paste0("<br>evidence = ", gene_evidence_types), ""),
          ifelse(nzchar(gene_source_apps), paste0("<br>gene apps = ", gene_source_apps), ""),
          ifelse(nzchar(apps_supported), paste0("<br>hit-supported apps = ", apps_supported), "")
        )
      ) %>%
      dplyr::arrange(dplyr::desc(degree), dplyr::desc(score)) %>%
      dplyr::mutate(
        label_plot = label
      )
    
    edf <- igraph::as_data_frame(g, what = "edges")
    edf$x <- vdf$x[match(edf$from, vdf$id)]
    edf$y <- vdf$y[match(edf$from, vdf$id)]
    edf$xend <- vdf$x[match(edf$to, vdf$id)]
    edf$yend <- vdf$y[match(edf$to, vdf$id)]
    
    ppi_label_mode <- input$ppi_label_mode %||% "Top"
    n_labels <- input$network_top_labels %||% 10
    
    ppi_lab_df <- get_label_subset(
      vdf,
      mode = ppi_label_mode,
      top_n = n_labels,
      order_col = "degree"
    )
    
    p <- plotly::plot_ly()
    
    for (i in seq_len(nrow(edf))) {
      p <- p %>%
        plotly::add_segments(
          x = edf$x[i], y = edf$y[i],
          xend = edf$xend[i], yend = edf$yend[i],
          inherit = FALSE,
          line = list(width = 0.7, color = "gray"),
          hoverinfo = "none",
          showlegend = FALSE
        )
    }
    
    gene_position_df <- vdf %>% dplyr::filter(link_class == "position")
    gene_effect_df <- vdf %>% dplyr::filter(link_class == "effect")
    gene_both_df <- vdf %>% dplyr::filter(link_class == "position+effect")
    gene_other_df <- vdf %>% dplyr::filter(!link_class %in% c("position", "effect", "position+effect"))
    
    if (nrow(gene_position_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_position_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: position-linked",
          marker = list(
            symbol = "circle",
            size = ~pmax(10, 8 + degree)
          )
        )
    }
    
    if (nrow(gene_effect_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_effect_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: effect-linked",
          marker = list(
            symbol = "diamond",
            size = ~pmax(10, 8 + degree)
          )
        )
    }
    
    if (nrow(gene_both_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_both_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes: position + effect",
          marker = list(
            symbol = "circle-open",
            size = ~pmax(10, 8 + degree)
          )
        )
    }
    
    if (nrow(gene_other_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = gene_other_df,
          x = ~x, y = ~y,
          text = ~hover_txt,
          hoverinfo = "text",
          name = "Genes",
          marker = list(
            symbol = "circle",
            size = ~pmax(10, 8 + degree)
          )
        )
    }
    
    if (nrow(ppi_lab_df) > 0) {
      p <- p %>%
        plotly::add_text(
          data = ppi_lab_df,
          x = ~x, y = ~y,
          text = ~label_plot,
          textposition = "top center",
          showlegend = FALSE,
          hoverinfo = "none"
        )
    }
    
    p %>%
      plotly::layout(
        title = "Prioritized gene PPI network",
        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        legend = list(orientation = "h")
      )
  })
  
  # ============================================================
  # MODAL NETWORK PLOT
  # ============================================================
  
  observeEvent(input$show_gene_networks, {
    showModal(
      shiny::modalDialog(
        title = "Gene networks",
        size = "l",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        
        tags$style(HTML("
        .modal-dialog {
          width: 95vw !important;
          max-width: 95vw !important;
        }
        .modal-body {
          max-height: 85vh;
          overflow-y: auto;
        }
        .gene-net-sidepanel {
          background: #f8f9fa;
          border: 1px solid #e5e7eb;
          border-radius: 8px;
          padding: 14px 16px;
        }
      ")),
        
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::div(
              class = "gene-net-sidepanel",
              
              shiny::selectInput(
                "pathway_label_mode",
                "Pathway network labels",
                choices = c("None", "Top", "All"),
                selected = "Top"
              ),
              
              shiny::selectInput(
                "ppi_label_mode",
                "PPI labels",
                choices = c("None", "Top", "All"),
                selected = "Top"
              ),
              
              shiny::sliderInput(
                "network_top_labels",
                "Top labels shown",
                min = 5,
                max = 30,
                value = 10,
                step = 1
              ),
              
              shiny::tags$div(
                style = "font-size:12px; margin-top:14px; margin-bottom:10px;",
                shiny::HTML(
                  paste0(
                    "<b>Interpretation:</b> ",
                    "filled circles = genes (including position-linked genes); ",
                    "diamonds = effect-linked genes; ",
                    "open circles = genes with both position and effect support; ",
                    "squares = pathways. ",
                    "In PPI, larger nodes are more connected."
                  )
                )
              ),
              
              shiny::tags$div(
                style = "padding:8px 12px; background:#ffffff; border-radius:6px; font-size:13px; border:1px solid #e5e7eb;",
                
                shiny::tags$p(
                  "Networks are built from the current prioritized genes. Gene inclusion follows the canonical model, combining:"
                ),
                
                shiny::tags$ul(
                  shiny::tags$li(
                    shiny::tags$b("Position-linked genes: "),
                    "genes physically overlapping the cluster or located in the same LD block as a GWAS hit."
                  ),
                  shiny::tags$li(
                    shiny::tags$b("Effect-linked genes: "),
                    "genes associated through regulatory or epigenetic signals such as GTEx or EWAS, even if located outside the cluster interval."
                  )
                ),
                
                shiny::tags$p(
                  "Gene scores are derived from prioritized GWAS hits linked to each gene (MATCH or MARKER). ",
                  "Only genes with non-zero gene score contribute to network ranking."
                )
              )
            )
          ),
          
          shiny::column(
            width = 9,
            shiny::tabsetPanel(
              id = "gene_networks_tabs",
              
              shiny::tabPanel(
                "Gene–pathway",
                br(),
                shinycssloaders::withSpinner(
                  plotly::plotlyOutput("gene_pathway_network", height = "780px"),
                  type = 4
                )
              ),
              
              shiny::tabPanel(
                "PPI",
                br(),
                shinycssloaders::withSpinner(
                  plotly::plotlyOutput("gene_ppi_network", height = "780px"),
                  type = 4
                )
              )
            )
          )
        )
      )
    )
  })
  
  ############################# END GENE NETWORKS ##############################
  ##############################################################################
  
  # ============================================================

  gene_bridge_lookup_df <- reactive({
    sess_dir <- session_dir_r()   # 
    req(sess_dir)
    
    app_files <- c(
      catalog = "catalog_gene_bridge.rds",
      gtex    = "gtex_gene_bridge.rds",
      nonsyn  = "nonsyn_gene_bridge.rds",
      ewasdis = "ewasdis_gene_bridge.rds",
      ewastum = "ewastum_gene_bridge.rds"
    )
    
    lst <- lapply(names(app_files), function(app) {
      fp <- file.path(sess_dir, app_files[[app]])
      if (!file.exists(fp)) return(NULL)
      
      df <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (!is.data.frame(df) || !nrow(df)) return(NULL)
      
      df %>%
        dplyr::transmute(
          gene = trimws(as.character(gene)),
          cluster_id = trimws(as.character(cluster_id)),
          source_app = tolower(trimws(as.character(source_app))),
          start = suppressWarnings(as.numeric(start)),
          end   = suppressWarnings(as.numeric(end))
        ) %>%
        dplyr::filter(
          !is.na(gene), nzchar(gene),
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(source_app), nzchar(source_app)
        ) %>%
        dplyr::mutate(
          gene_mid = dplyr::if_else(
            is.finite(start) & is.finite(end),
            (start + end) / 2,
            NA_real_
          )
        ) %>%
        dplyr::distinct()
    })
    
    dplyr::bind_rows(lst)
  })
  
  nonsyn_gene_hit_bridge_df <- reactive({
    req(session_dir_r())
    
    fp <- file.path(session_dir_r(), "nonsyn_gene_hit_bridge.rds")
    
    if (!file.exists(fp)) {
      return(tibble::tibble(
        cluster_id = character(),
        source_app = character(),
        gene = character(),
        chr = integer(),
        position = integer(),
        id_hit = character(),
        hit_key = character()
      ))
    }
    
    df <- tryCatch(readRDS(fp), error = function(e) NULL)
    
    if (!is.data.frame(df) || !nrow(df)) {
      return(tibble::tibble(
        cluster_id = character(),
        source_app = character(),
        gene = character(),
        chr = integer(),
        position = integer(),
        id_hit = character(),
        hit_key = character()
      ))
    }
    
    df %>%
      dplyr::transmute(
        cluster_id = trimws(as.character(cluster_id)),
        source_app = tolower(trimws(as.character(source_app))),
        gene = trimws(as.character(gene)),
        chr = suppressWarnings(as.integer(chr)),
        position = suppressWarnings(as.integer(position)),
        id_hit = trimws(as.character(id_hit)),
        hit_key = trimws(as.character(hit_key))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(source_app), nzchar(source_app),
        !is.na(gene), nzchar(gene),
        !is.na(hit_key), nzchar(hit_key)
      ) %>%
      dplyr::distinct()
  })
  
  # ============================================================
  # GWAS HIT PRIORITY COMPONENTS · v2  (APP LEVEL)
  # ============================================================
  
  gwas_hit_priority_plot_df_v2 <- reactive({
    df <- gwas_hit_priority_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No GWAS hit priority data available.")
    )
    
    extract_chr_num <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      
      dplyr::case_when(
        x == "X" ~ 23,
        x == "Y" ~ 24,
        x %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(x))
      )
    }
    
    normalize_component_name <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x_low <- tolower(x)
      
      dplyr::case_when(
        x_low %in% c(
          "gwas_significance",
          "gwas_signif",
          "gwas_significance_bonus",
          "gwas_signif_bonus",
          "significance",
          "signif",
          "gwas_sig",
          "gwas_sig_bonus"
        ) ~ "gwas_significance",
        
        x_low %in% c("catalog") ~ "catalog",
        x_low %in% c("gtex") ~ "gtex",
        x_low %in% c("nonsyn", "dbnsfp", "non_syn", "non-syn") ~ "nonsyn",
        x_low %in% c("ewastum", "ewas_tumor", "tumor") ~ "ewastum",
        x_low %in% c("ewasdis", "ewas_disease", "disease") ~ "ewasdis",
        
        TRUE ~ x_low
      )
    }
    
    parse_named_numeric_pairs <- function(x) {
      x <- as.character(x %||% "")
      x <- trimws(x)
      
      if (!nzchar(x)) {
        return(tibble::tibble(
          source_app = character(0),
          value = numeric(0)
        ))
      }
      
      parts <- unlist(strsplit(x, "\\s*\\|\\s*"))
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      
      if (!length(parts)) {
        return(tibble::tibble(
          source_app = character(0),
          value = numeric(0)
        ))
      }
      
      out <- lapply(parts, function(part) {
        m <- stringr::str_match(
          part,
          "^([^=]+)=\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)$"
        )
        
        if (is.na(m[1, 1])) return(NULL)
        
        tibble::tibble(
          source_app = normalize_component_name(m[1, 2]),
          value = suppressWarnings(as.numeric(m[1, 3]))
        )
      })
      
      out <- out[!vapply(out, is.null, logical(1))]
      
      if (!length(out)) {
        return(tibble::tibble(
          source_app = character(0),
          value = numeric(0)
        ))
      }
      
      dplyr::bind_rows(out) %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(source_app) %>%
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
    }
    
    parse_named_state_pairs <- function(x) {
      x <- as.character(x %||% "")
      x <- trimws(x)
      
      if (!nzchar(x)) {
        return(tibble::tibble(
          source_app = character(0),
          support_state = character(0)
        ))
      }
      
      parts <- unlist(strsplit(x, "\\s*\\|\\s*"))
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      
      if (!length(parts)) {
        return(tibble::tibble(
          source_app = character(0),
          support_state = character(0)
        ))
      }
      
      out <- lapply(parts, function(part) {
        m <- stringr::str_match(part, "^([^:]+):\\s*(MATCH|MARKER|NOLINK|NOHIT)$")
        
        if (is.na(m[1, 1])) return(NULL)
        
        tibble::tibble(
          source_app = normalize_component_name(m[1, 2]),
          support_state = trimws(m[1, 3])
        )
      })
      
      out <- out[!vapply(out, is.null, logical(1))]
      
      if (!length(out)) {
        return(tibble::tibble(
          source_app = character(0),
          support_state = character(0)
        ))
      }
      
      dplyr::bind_rows(out) %>%
        dplyr::distinct(source_app, .keep_all = TRUE)
    }
    
    priority_col <- dplyr::case_when(
      "priority_class (ABS)" %in% names(df) ~ "priority_class (ABS)",
      "priority_class" %in% names(df) ~ "priority_class",
      TRUE ~ NA_character_
    )
    
    score_col <- dplyr::case_when(
      "score" %in% names(df) ~ "score",
      "gwas_hit_score" %in% names(df) ~ "gwas_hit_score",
      TRUE ~ NA_character_
    )
    
    validate(
      need(!is.na(priority_col), "Priority class column not found."),
      need(!is.na(score_col), "Score column not found."),
      need("score_breakdown" %in% names(df), "score_breakdown column not found."),
      need("support_signature" %in% names(df), "support_signature column not found.")
    )
    
    base_df <- df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        gwas_hit = as.character(gwas_hit),
        gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
        chr_num = extract_chr_num(cluster_id),
        total_score = dplyr::coalesce(as.numeric(.data[[score_col]]), 0),
        priority_class_abs = dplyr::coalesce(as.character(.data[[priority_col]]), ""),
        support_signature = dplyr::coalesce(as.character(support_signature), ""),
        score_breakdown = dplyr::coalesce(as.character(score_breakdown), ""),
        hit_label = paste0(cluster_id, " · ", gwas_hit)
      ) %>%
      dplyr::mutate(
        chr_num = dplyr::coalesce(chr_num, 999),
        gwas_pos = dplyr::coalesce(gwas_pos, Inf)
      )
    
    hit_levels <- base_df %>%
      dplyr::arrange(chr_num, gwas_pos, cluster_id, gwas_hit) %>%
      dplyr::pull(hit_label) %>%
      unique() %>%
      rev()
    
    plot_df <- purrr::pmap_dfr(
      list(
        base_df$cluster_id,
        base_df$gwas_hit,
        base_df$gwas_pos,
        base_df$total_score,
        base_df$priority_class_abs,
        base_df$support_signature,
        base_df$score_breakdown,
        base_df$hit_label
      ),
      function(cluster_id, gwas_hit, gwas_pos, total_score, priority_class_abs,
               support_signature, score_breakdown, hit_label) {
        
        breakdown_df <- parse_named_numeric_pairs(score_breakdown) %>%
          dplyr::rename(score_contribution = value)
        
        signature_df <- parse_named_state_pairs(support_signature)
        
        merged_df <- breakdown_df %>%
          dplyr::full_join(signature_df, by = "source_app") %>%
          dplyr::mutate(
            source_app = dplyr::coalesce(as.character(source_app), "unknown"),
            source_app = normalize_component_name(source_app),
            score_contribution = dplyr::coalesce(as.numeric(score_contribution), 0),
            support_state = dplyr::case_when(
              source_app == "gwas_significance" ~ "SIGNIF",
              TRUE ~ dplyr::coalesce(as.character(support_state), "UNKNOWN")
            )
          )
        
        if (!nrow(merged_df)) {
          return(tibble::tibble(
            hit_label = hit_label,
            cluster_id = cluster_id,
            gwas_hit = gwas_hit,
            gwas_pos = gwas_pos,
            total_score = total_score,
            plotted_sum = 0,
            score_diff = total_score,
            score_matches = abs(total_score) < 1e-8,
            priority_class = priority_class_abs,
            source_app = "unknown",
            support_state = "UNKNOWN",
            component_label = "No score components",
            score_contribution = 0
          ))
        }
        
        plotted_sum <- sum(merged_df$score_contribution, na.rm = TRUE)
        score_diff <- total_score - plotted_sum
        score_matches <- abs(total_score - plotted_sum) < 0.011
        
        merged_df %>%
          dplyr::mutate(
            hit_label = hit_label,
            cluster_id = cluster_id,
            gwas_hit = gwas_hit,
            gwas_pos = gwas_pos,
            total_score = total_score,
            plotted_sum = plotted_sum,
            score_diff = score_diff,
            score_matches = score_matches,
            priority_class = priority_class_abs,
            component_label = dplyr::case_when(
              source_app == "gwas_significance" ~ "GWAS significance",
              TRUE ~ paste0(source_app, " (", support_state, ")")
            )
          ) %>%
          dplyr::select(
            hit_label,
            cluster_id,
            gwas_hit,
            gwas_pos,
            total_score,
            plotted_sum,
            score_diff,
            score_matches,
            priority_class,
            source_app,
            support_state,
            component_label,
            score_contribution
          )
      }
    )
    
    app_levels <- c("gwas_significance", "catalog", "gtex", "nonsyn", "ewastum", "ewasdis", "unknown")
    state_levels <- c("SIGNIF", "MATCH", "MARKER", "NOLINK", "NOHIT", "UNKNOWN")
    
    plot_df %>%
      dplyr::mutate(
        source_app = factor(source_app, levels = app_levels),
        support_state = factor(support_state, levels = state_levels),
        component_label = factor(
          component_label,
          levels = c(
            "GWAS significance",
            "catalog (MATCH)", "catalog (MARKER)", "catalog (NOLINK)", "catalog (NOHIT)",
            "gtex (MATCH)", "gtex (MARKER)", "gtex (NOLINK)", "gtex (NOHIT)",
            "nonsyn (MATCH)", "nonsyn (MARKER)", "nonsyn (NOLINK)", "nonsyn (NOHIT)",
            "ewastum (MATCH)", "ewastum (MARKER)", "ewastum (NOLINK)", "ewastum (NOHIT)",
            "ewasdis (MATCH)", "ewasdis (MARKER)", "ewasdis (NOLINK)", "ewasdis (NOHIT)",
            "unknown (UNKNOWN)", "No score components"
          )
        ),
        hit_label = factor(hit_label, levels = hit_levels),
        hover_txt = dplyr::case_when(
          as.character(component_label) == "GWAS significance" ~ paste0(
            "Hit: ", as.character(hit_label),
            "<br>Component: GWAS significance",
            "<br>Contribution: ", round(score_contribution, 2),
            "<br>Total score: ", round(total_score, 2),
            "<br>Plotted sum: ", round(plotted_sum, 2),
            "<br>Priority class: ", priority_class,
            ifelse(
              score_matches,
              "",
              paste0(
                "<br><span style='color:#b22222;'><b>Warning:</b> plotted sum differs from score by ",
                round(score_diff, 4),
                "</span>"
              )
            )
          ),
          TRUE ~ paste0(
            "Hit: ", as.character(hit_label),
            "<br>App: ", as.character(source_app),
            "<br>Support: ", as.character(support_state),
            "<br>Contribution: ", round(score_contribution, 2),
            "<br>Total score: ", round(total_score, 2),
            "<br>Plotted sum: ", round(plotted_sum, 2),
            "<br>Priority class: ", priority_class,
            ifelse(
              score_matches,
              "",
              paste0(
                "<br><span style='color:#b22222;'><b>Warning:</b> plotted sum differs from score by ",
                round(score_diff, 4),
                "</span>"
              )
            )
          )
        )
      )
  })
  
  output$gwas_hit_priority_plot_popup_v2 <- plotly::renderPlotly({
    plot_df <- gwas_hit_priority_plot_df_v2()
    
    validate(
      need(is.data.frame(plot_df) && nrow(plot_df) > 0, "No plot data available.")
    )
    
    component_colors <- c(
      "GWAS significance" = "#222222",
      
      "catalog (MATCH)"  = "#1f78b4",
      "catalog (MARKER)" = "#a6cee3",
      "catalog (NOLINK)" = "#d9ecf7",
      "catalog (NOHIT)"  = "#edf5fb",
      
      "gtex (MATCH)"     = "#33a02c",
      "gtex (MARKER)"    = "#b2df8a",
      "gtex (NOLINK)"    = "#dff0c9",
      "gtex (NOHIT)"     = "#eef7e5",
      
      "nonsyn (MATCH)"   = "#6a3d9a",
      "nonsyn (MARKER)"  = "#cab2d6",
      "nonsyn (NOLINK)"  = "#e6dcef",
      "nonsyn (NOHIT)"   = "#f2edf7",
      
      "ewastum (MATCH)"  = "#ff7f00",
      "ewastum (MARKER)" = "#fdbf6f",
      "ewastum (NOLINK)" = "#fde3bf",
      "ewastum (NOHIT)"  = "#fff1df",
      
      "ewasdis (MATCH)"  = "#e31a1c",
      "ewasdis (MARKER)" = "#fb9a99",
      "ewasdis (NOLINK)" = "#f8c9c8",
      "ewasdis (NOHIT)"  = "#fde8e8",
      
      "unknown (UNKNOWN)" = "#bdbdbd",
      "No score components" = "#d9d9d9"
    )
    
    component_order <- c(
      "GWAS significance",
      
      "catalog (MATCH)", "catalog (MARKER)", "catalog (NOLINK)", "catalog (NOHIT)",
      "gtex (MATCH)", "gtex (MARKER)", "gtex (NOLINK)", "gtex (NOHIT)",
      "nonsyn (MATCH)", "nonsyn (MARKER)", "nonsyn (NOLINK)", "nonsyn (NOHIT)",
      "ewastum (MATCH)", "ewastum (MARKER)", "ewastum (NOLINK)", "ewastum (NOHIT)",
      "ewasdis (MATCH)", "ewasdis (MARKER)", "ewasdis (NOLINK)", "ewasdis (NOHIT)",
      "unknown (UNKNOWN)", "No score components"
    )
    
    plot_df <- plot_df %>%
      dplyr::mutate(
        component_label = as.character(component_label),
        component_label = ifelse(
          component_label %in% component_order,
          component_label,
          "unknown (UNKNOWN)"
        ),
        component_label = factor(component_label, levels = component_order)
      )
    
    max_score <- suppressWarnings(max(c(plot_df$plotted_sum, plot_df$total_score), na.rm = TRUE))
    if (!is.finite(max_score) || max_score <= 0) max_score <- 1
    
    plotly::plot_ly(
      data = plot_df,
      x = ~score_contribution,
      y = ~hit_label,
      color = ~component_label,
      colors = component_colors,
      type = "bar",
      orientation = "h",
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        barmode = "stack",
        title = list(
          text = "GWAS hit priority score components",
          y = 0.98
        ),
        xaxis = list(
          title = "Priority score",
          range = c(0, max_score * 1.05)
        ),
        yaxis = list(
          title = "",
          tickfont = list(size = 10)
        ),
        hovermode = "closest",
        showlegend = TRUE,
        legend = list(
          title = list(text = "Score components"),
          orientation = "h",
          x = 0,
          y = -0.18,
          font = list(size = 10)
        ),
        margin = list(l = 260, r = 30, t = 55, b = 120)
      )
  })
  
  observeEvent(input[["gwas_hit_priority-show_gwas_hit_priority_plot"]], {
    df <- gwas_hit_priority_df_v2()
    req(is.data.frame(df), nrow(df) > 0)
    
    showModal(
      modalDialog(
        title = "GWAS hit priority plot",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$p(
          "Each horizontal bar represents one prioritized GWAS hit. ",
          "Bar segments show the contribution of the current GWAS-hit priority components: ",
          "the app-specific evidence components (for example MATCH or MARKER support for each app) ",
          "together with the GWAS significance contribution of the hit itself. ",
          "The final GWAS-hit score is therefore driven by both cross-app evidence and the statistical strength of the GWAS association."
        ),
        shinycssloaders::withSpinner(plotly::plotlyOutput("gwas_hit_priority_plot_popup_v2", height = "680px"))
      )
    )
  })
  
  # ============================================================
  # GWAS HIT PRIORITY Manhattan hits

  gwas_hit_manhattan_plot_df_v2 <- reactive({
    df <- gwas_hit_priority_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No GWAS hit priority data available.")
    )
    
    extract_chr_label <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- sub("^chr", "", x, ignore.case = TRUE)
      x <- sub("_.*$", "", x)
      x <- toupper(x)
      x
    }
    
    extract_chr_num <- function(x) {
      lab <- extract_chr_label(x)
      
      dplyr::case_when(
        lab == "X" ~ 23,
        lab == "Y" ~ 24,
        lab %in% c("M", "MT") ~ 25,
        TRUE ~ suppressWarnings(as.numeric(lab))
      )
    }
    
    score_col <- dplyr::case_when(
      "score" %in% names(df) ~ "score",
      "gwas_hit_score" %in% names(df) ~ "gwas_hit_score",
      TRUE ~ NA_character_
    )
    
    abs_class_col <- dplyr::case_when(
      "priority_class (ABS)" %in% names(df) ~ "priority_class (ABS)",
      "priority_class" %in% names(df) ~ "priority_class",
      TRUE ~ NA_character_
    )
    
    rel_class_col <- dplyr::case_when(
      "priority_class (REL)" %in% names(df) ~ "priority_class (REL)",
      TRUE ~ NA_character_
    )
    
    validate(
      need(!is.na(score_col), "Score column not found."),
      need("cluster_id" %in% names(df), "cluster_id column not found."),
      need("gwas_hit" %in% names(df), "gwas_hit column not found."),
      need("gwas_pos" %in% names(df), "gwas_pos column not found.")
    )
    
    base_df <- df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        gwas_hit = as.character(gwas_hit),
        gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
        score = dplyr::coalesce(as.numeric(.data[[score_col]]), 0),
        chr_label = extract_chr_label(cluster_id),
        chr_num = extract_chr_num(cluster_id),
        priority_class_abs = if (!is.na(abs_class_col)) dplyr::coalesce(as.character(.data[[abs_class_col]]), "") else "",
        priority_class_rel = if (!is.na(rel_class_col)) dplyr::coalesce(as.character(.data[[rel_class_col]]), "") else "",
        evidence_summary = if ("evidence_summary" %in% names(df)) dplyr::coalesce(as.character(evidence_summary), "") else "",
        support_signature = if ("support_signature" %in% names(df)) dplyr::coalesce(as.character(support_signature), "") else "",
        score_breakdown = if ("score_breakdown" %in% names(df)) dplyr::coalesce(as.character(score_breakdown), "") else ""
      ) %>%
      dplyr::filter(!is.na(chr_num), !is.na(gwas_pos), is.finite(gwas_pos)) %>%
      dplyr::arrange(chr_num, gwas_pos, gwas_hit)
    
    validate(
      need(nrow(base_df) > 0, "No valid chromosome/position data available for Manhattan plot.")
    )
    
    chr_sizes <- base_df %>%
      dplyr::group_by(chr_num, chr_label) %>%
      dplyr::summarise(
        chr_len = max(gwas_pos, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(chr_num) %>%
      dplyr::mutate(
        offset = dplyr::lag(cumsum(chr_len), default = 0),
        center = offset + chr_len / 2
      )
    
    plot_df <- base_df %>%
      dplyr::left_join(
        chr_sizes %>% dplyr::select(chr_num, chr_label, offset, center),
        by = c("chr_num", "chr_label")
      ) %>%
      dplyr::mutate(
        pos_cum = gwas_pos + offset,
        point_group = ifelse(chr_num %% 2 == 0, "even", "odd"),
        hover_txt = paste0(
          "Hit: ", gwas_hit,
          "<br>Cluster: ", cluster_id,
          "<br>Chr: ", chr_label,
          "<br>Position: ", format(gwas_pos, big.mark = ",", scientific = FALSE),
          "<br>Score: ", round(score, 2),
          ifelse(nzchar(priority_class_abs), paste0("<br>Priority class (ABS): ", priority_class_abs), ""),
          ifelse(nzchar(priority_class_rel), paste0("<br>Priority class (REL): ", priority_class_rel), ""),
          ifelse(nzchar(evidence_summary), paste0("<br>Evidence: ", evidence_summary), ""),
          ifelse(nzchar(support_signature), paste0("<br>Support: ", support_signature), ""),
          ifelse(nzchar(score_breakdown), paste0("<br>Breakdown: ", score_breakdown), "")
        )
      )
    
    list(
      points = plot_df,
      chr_axis = chr_sizes
    )
  })
  
  output$gwas_hit_manhattan_plot_v2XXXXX <- plotly::renderPlotly({
    manh <- gwas_hit_manhattan_plot_df_v2()
    plot_df <- manh$points
    chr_axis <- manh$chr_axis
    
    validate(
      need(is.data.frame(plot_df) && nrow(plot_df) > 0, "No Manhattan plot data available."),
      need(is.data.frame(chr_axis) && nrow(chr_axis) > 0, "No chromosome axis data available.")
    )
    
    odd_df <- plot_df %>% dplyr::filter(point_group == "odd")
    even_df <- plot_df %>% dplyr::filter(point_group == "even")
    
    p <- plotly::plot_ly()
    
    if (nrow(odd_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = odd_df,
          x = ~pos_cum,
          y = ~score,
          text = ~hover_txt,
          hoverinfo = "text",
          marker = list(size = 7, color = "#1f78b4"),
          name = "Odd chromosomes",
          showlegend = FALSE
        )
    }
    
    if (nrow(even_df) > 0) {
      p <- p %>%
        plotly::add_markers(
          data = even_df,
          x = ~pos_cum,
          y = ~score,
          text = ~hover_txt,
          hoverinfo = "text",
          marker = list(size = 7, color = "#a6cee3"),
          name = "Even chromosomes",
          showlegend = FALSE
        )
    }
    
    sep_lines <- chr_axis %>%
      dplyr::mutate(xline = offset) %>%
      dplyr::filter(xline > 0)
    
    shapes <- lapply(seq_len(nrow(sep_lines)), function(i) {
      list(
        type = "line",
        x0 = sep_lines$xline[i],
        x1 = sep_lines$xline[i],
        y0 = 0,
        y1 = 1,
        xref = "x",
        yref = "paper",
        line = list(color = "rgba(120,120,120,0.35)", width = 1, dash = "dot")
      )
    })
    
    p %>%
      plotly::layout(
        title = list(
          text = "GWAS hit priority Manhattan plot",
          y = 0.98
        ),
        xaxis = list(
          title = "Chromosome / genomic position",
          tickmode = "array",
          tickvals = chr_axis$center,
          ticktext = chr_axis$chr_label
        ),
        yaxis = list(
          title = "Priority score",
          rangemode = "tozero"
        ),
        hovermode = "closest",
        margin = list(l = 70, r = 30, t = 55, b = 60),
        shapes = shapes
      )
  })
  
  output$gwas_hit_manhattan_plot_v2 <- plotly::renderPlotly({
    manh <- gwas_hit_manhattan_plot_df_v2()
    plot_df <- manh$points
    
    validate(
      need(is.data.frame(plot_df) && nrow(plot_df) > 0, "No GWAS hit plot data available.")
    )
    
    # ------------------------------------------------------------
    # Detect genomic position column safely
    # ------------------------------------------------------------
    pos_candidates <- c("position", "pos", "bp", "BP", "pos_cum")
    pos_col <- pos_candidates[pos_candidates %in% names(plot_df)][1]
    
    validate(
      need(!is.na(pos_col) && nzchar(pos_col), "No genomic position column available for GWAS hit plot."),
      need("chr" %in% names(plot_df), "Column 'chr' not found in GWAS hit plot data."),
      need("score" %in% names(plot_df), "Column 'score' not found in GWAS hit plot data.")
    )
    
    # ------------------------------------------------------------
    # Clean / harmonize
    # ------------------------------------------------------------
    df <- plot_df %>%
      dplyr::mutate(
        chr = as.character(.data$chr),
        chr_clean = gsub("^chr", "", .data$chr, ignore.case = TRUE),
        chr_num = dplyr::case_when(
          toupper(.data$chr_clean) == "X" ~ 23,
          toupper(.data$chr_clean) == "Y" ~ 24,
          toupper(.data$chr_clean) %in% c("M", "MT") ~ 25,
          TRUE ~ suppressWarnings(as.numeric(.data$chr_clean))
        ),
        chr_num = ifelse(is.na(.data$chr_num), 999, .data$chr_num),
        chr_label = dplyr::case_when(
          .data$chr_num == 23 ~ "X",
          .data$chr_num == 24 ~ "Y",
          .data$chr_num == 25 ~ "MT",
          .data$chr_num == 999 ~ .data$chr,
          TRUE ~ as.character(.data$chr_num)
        ),
        pos_chr = suppressWarnings(as.numeric(.data[[pos_col]])),
        score = suppressWarnings(as.numeric(.data$score))
      ) %>%
      dplyr::filter(
        !is.na(.data$chr_label),
        !is.na(.data$pos_chr),
        is.finite(.data$pos_chr),
        !is.na(.data$score),
        is.finite(.data$score)
      ) %>%
      dplyr::arrange(.data$chr_num, .data$pos_chr)
    
    validate(
      need(nrow(df) > 0, "No valid GWAS hit plot data available after filtering.")
    )
    
    # ------------------------------------------------------------
    # Ordered chromosome factor + alternating chromosome colors
    # ------------------------------------------------------------
    chr_info <- df %>%
      dplyr::distinct(.data$chr_label, .data$chr_num) %>%
      dplyr::arrange(.data$chr_num) %>%
      dplyr::mutate(
        chr_parity = ifelse(seq_len(dplyr::n()) %% 2 == 1, "odd_chr", "even_chr")
      )
    
    df <- df %>%
      dplyr::left_join(chr_info, by = c("chr_label", "chr_num")) %>%
      dplyr::mutate(
        chr_label = factor(.data$chr_label, levels = chr_info$chr_label)
      )
    
    # ------------------------------------------------------------
    # Shared Y range across all facets
    # ------------------------------------------------------------
    y_min <- min(df$score, na.rm = TRUE)
    y_max <- max(df$score, na.rm = TRUE)
    
    if (!is.finite(y_min)) y_min <- 0
    if (!is.finite(y_max)) y_max <- 1
    if (y_min == y_max) y_max <- y_min + 1
    
    # ------------------------------------------------------------
    # Plot: 1 single row, one facet per chromosome
    # ------------------------------------------------------------
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$pos_chr,
        y = .data$score,
        text = .data$hover_txt,
        color = .data$chr_parity
      )
    ) +
      ggplot2::geom_point(
        size = 2.4,
        alpha = 0.85
      ) +
      ggplot2::facet_wrap(
        ~ chr_label,
        scales = "free_x",
        nrow = 1
      ) +
      ggplot2::scale_color_manual(
        values = c(
          odd_chr = "#1f78b4",
          even_chr = "orange"
        )
      ) +
      ggplot2::scale_x_continuous(
        breaks = function(x) pretty(x, n = 3),
        labels = function(x) format(round(x), big.mark = ",", scientific = FALSE)
      ) +
      ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
      ggplot2::labs(
        title = "GWAS hit priority by chromosome",
        x = "Genomic position",
        y = "Priority score"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
        strip.text = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          size = 8
        ),
        legend.position = "none"
      )
    
    gg <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = FALSE)
    
    gg %>%
      plotly::layout(
        hovermode = "closest",
        margin = list(l = 75, r = 20, t = 60, b = 90),
        showlegend = FALSE
      )
  })
  
  
  observeEvent(input[["gwas_hit_priority-show_gwas_hit_manhattan_plot"]], {
    df <- gwas_hit_priority_df_v2()
    req(is.data.frame(df), nrow(df) > 0)
    
    showModal(
      modalDialog(
        title = "GWAS hit priority Manhattan plot",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$p(
          "Each point represents one prioritized GWAS hit. ",
          "The X axis shows chromosome and genomic position, and the Y axis shows the GWAS-hit priority score. ",
          "This provides a Manhattan-style view of the prioritized hits using the score reported in the priority table. ",
          "Hover over each point to inspect the hit, cluster, score, priority classes, support signature, and score breakdown."
        ),
        
        plotly::plotlyOutput("gwas_hit_manhattan_plot_v2", height = "650px")
      )
    )
  })
  
  # ============================================================
  # PRIORITY AUDIT · HELPERS
  # ============================================================
  collapse_unique_semicolon <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return("")
    
    vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
    vals <- trimws(vals)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) return("")
    
    paste(sort(unique(vals)), collapse = "; ")
  }
  
  count_semicolon_items <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    if (is.na(x) || !nzchar(x)) return(0L)
    
    vals <- unlist(strsplit(x, "\\s*;\\s*"))
    vals <- trimws(vals)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    
    as.integer(length(unique(vals)))
  }
  
  split_semicolon_unique <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character(0))
    
    vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
    vals <- trimws(vals)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    unique(vals)
  }
  
  safe_first_col <- function(df, candidates, default = NULL) {
    hit <- intersect(candidates, names(df))
    if (length(hit)) hit[1] else default
  }
  
  # ============================================================
  # CANONICAL GENE TABLE FOR AUDIT
  # ============================================================
  prioritized_gene_audit_base_df <- reactive({
    gdf <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(gdf) && nrow(gdf) > 0, "No prioritized genes available for audit.")
    )
    
    gene_col <- safe_first_col(gdf, c("gene", "gene_name"))
    score_col <- safe_first_col(gdf, c("gene_score", "priority_score", "score"))
    clusters_col <- safe_first_col(gdf, c("cluster_ids", "cluster_id"))
    support_col <- safe_first_col(gdf, c("support_signature", "apps_supported", "evidence_summary"), default = NULL)
    
    validate(
      need(!is.null(gene_col), "Gene column not available in prioritized_gene_df_v2()."),
      need(!is.null(score_col), "Gene score column not available in prioritized_gene_df_v2()."),
      need(!is.null(clusters_col), "cluster_ids / cluster_id column not available in prioritized_gene_df_v2().")
    )
    
    gdf %>%
      dplyr::transmute(
        gene = trimws(as.character(.data[[gene_col]])),
        gene_score = dplyr::coalesce(as.numeric(.data[[score_col]]), 0),
        cluster_ids_raw = trimws(as.character(.data[[clusters_col]])),
        gene_support_signature = if (!is.null(support_col)) {
          dplyr::coalesce(as.character(.data[[support_col]]), "")
        } else {
          ""
        }
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      tidyr::separate_rows(cluster_ids_raw, sep = "\\s*;\\s*") %>%
      dplyr::rename(cluster_id = cluster_ids_raw) %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id))
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::distinct(gene, cluster_id, .keep_all = TRUE)
  })
  
  # ============================================================
  # BLOCK ↔ GENE AUDIT
  # ============================================================
  block_gene_audit_df <- reactive({
    blk <- prioritized_block_df_v2()
    gen <- prioritized_gene_audit_base_df()
    
    validate(
      need(is.data.frame(blk) && nrow(blk) > 0, "No prioritized blocks available for audit."),
      need(is.data.frame(gen) && nrow(gen) > 0, "No prioritized genes available for audit.")
    )
    
    blk_std <- blk %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        block_score_from_hits = dplyr::coalesce(as.numeric(block_score_from_hits), 0),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), ""),
        n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
        apps_supported = dplyr::coalesce(as.character(apps_supported), ""),
        priority_class = dplyr::coalesce(as.character(priority_class), ""),
        priority_class_relative = dplyr::coalesce(as.character(priority_class_relative), "")
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(block_id), nzchar(block_id))
    
    block_gene_map <- blk_std %>%
      dplyr::select(cluster_id, block_id, block_score, block_score_from_hits, top_gwas_hit_score,
                    n_gwas_hits, genes_in_block, n_genes_in_block, apps_supported,
                    priority_class, priority_class_relative) %>%
      tidyr::separate_rows(genes_in_block, sep = "\\s*;\\s*") %>%
      dplyr::rename(gene = genes_in_block) %>%
      dplyr::mutate(
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::distinct(cluster_id, block_id, gene, .keep_all = TRUE)
    
    audit <- block_gene_map %>%
      dplyr::left_join(
        gen,
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::mutate(
        gene_score = dplyr::coalesce(as.numeric(gene_score), 0),
        gene_block_ratio = dplyr::if_else(
          is.finite(block_score) & block_score > 0,
          gene_score / block_score,
          NA_real_
        ),
        audit_flag = dplyr::case_when(
          block_score >= 20 & gene_score < 3 ~ "BLOCK_STRONG_GENE_WEAK",
          gene_score >= 15 & block_score < 10 ~ "GENE_STRONG_BLOCK_WEAK",
          n_genes_in_block == 0 ~ "BLOCK_WITHOUT_GENE_CONTENT",
          TRUE ~ "CONSISTENT"
        ),
        audit_comment = dplyr::case_when(
          audit_flag == "BLOCK_STRONG_GENE_WEAK" ~ "Strong block score but weak prioritized gene score inside the block.",
          audit_flag == "GENE_STRONG_BLOCK_WEAK" ~ "Strong prioritized gene score inside a relatively weak block.",
          audit_flag == "BLOCK_WITHOUT_GENE_CONTENT" ~ "Block has no annotated genes in genes_in_block.",
          TRUE ~ "Block and gene scores are broadly coherent."
        )
      )
    
    # Resum per bloc
    audit %>%
      dplyr::group_by(cluster_id, block_id) %>%
      dplyr::summarise(
        block_score = dplyr::first(block_score),
        block_score_from_hits = dplyr::first(block_score_from_hits),
        top_gwas_hit_score = dplyr::first(top_gwas_hit_score),
        n_gwas_hits = dplyr::first(n_gwas_hits),
        n_genes_in_block = dplyr::first(n_genes_in_block),
        apps_supported = dplyr::first(apps_supported),
        priority_class = dplyr::first(priority_class),
        priority_class_relative = dplyr::first(priority_class_relative),
        
        genes_in_block = paste(sort(unique(gene[!is.na(gene) & nzchar(gene)])), collapse = "; "),
        n_prioritized_genes_in_block = sum(is.finite(gene_score) & gene_score > 0, na.rm = TRUE),
        top_gene_score_in_block = suppressWarnings(max(gene_score, na.rm = TRUE)),
        sum_gene_score_in_block = round(sum(gene_score, na.rm = TRUE), 2),
        top_gene_in_block = gene[which.max(dplyr::coalesce(gene_score, 0))[1]],
        
        block_gene_ratio = dplyr::if_else(
          is.finite(block_score) & block_score > 0,
          suppressWarnings(max(gene_score, na.rm = TRUE)) / block_score,
          NA_real_
        ),
        
        audit_flag = dplyr::case_when(
          block_score >= 20 & suppressWarnings(max(gene_score, na.rm = TRUE)) < 3 ~ "BLOCK_DOMINANT_WITH_WEAK_GENES",
          block_score < 10 & suppressWarnings(max(gene_score, na.rm = TRUE)) >= 12 ~ "GENE_STRONG_IN_WEAK_BLOCK",
          n_genes_in_block == 0 ~ "TOP_BLOCK_WITHOUT_GENE_CONTENT",
          TRUE ~ "CONSISTENT"
        ),
        
        audit_comment = dplyr::case_when(
          audit_flag == "BLOCK_DOMINANT_WITH_WEAK_GENES" ~ "High-scoring block but no strong prioritized gene inside it.",
          audit_flag == "GENE_STRONG_IN_WEAK_BLOCK" ~ "A strong prioritized gene is placed inside a relatively weak block.",
          audit_flag == "TOP_BLOCK_WITHOUT_GENE_CONTENT" ~ "Block is prioritized but has no genes in genes_in_block.",
          TRUE ~ "Block-level and gene-level signals are broadly coherent."
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        top_gene_score_in_block = dplyr::if_else(is.finite(top_gene_score_in_block), top_gene_score_in_block, 0),
        top_gene_in_block = dplyr::coalesce(as.character(top_gene_in_block), "")
      ) %>%
      dplyr::arrange(dplyr::desc(block_score), cluster_id, block_id)
  })
  
  # ============================================================
  # CLUSTER ↔ BLOCK AUDIT
  # ============================================================
  cluster_block_audit_df <- reactive({
    cl <- prioritized_cluster_df_v2()
    blk <- prioritized_block_df_v2()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available for audit."),
      need(is.data.frame(blk) && nrow(blk) > 0, "No prioritized blocks available for audit.")
    )
    
    blk_sum <- blk %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        block_score_from_hits = dplyr::coalesce(as.numeric(block_score_from_hits), 0),
        n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        n_blocks = dplyr::n_distinct(block_id[!is.na(block_id) & nzchar(block_id)]),
        top_block_score_in_cluster = suppressWarnings(max(block_score, na.rm = TRUE)),
        sum_block_score_in_cluster = round(sum(block_score, na.rm = TRUE), 2),
        top_block_id = block_id[which.max(dplyr::coalesce(block_score, 0))[1]],
        blocks_in_cluster = paste(sort(unique(block_id[!is.na(block_id) & nzchar(block_id)])), collapse = "; "),
        genes_in_blocks = collapse_unique_semicolon(genes_in_block),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        n_unique_genes_in_blocks = vapply(genes_in_blocks, count_semicolon_items, integer(1)),
        top_block_score_in_cluster = dplyr::if_else(
          is.finite(top_block_score_in_cluster),
          top_block_score_in_cluster,
          0
        ),
        cluster_block_ratio = dplyr::if_else(
          is.finite(sum_block_score_in_cluster) & sum_block_score_in_cluster > 0,
          top_block_score_in_cluster / sum_block_score_in_cluster,
          NA_real_
        )
      )
    
    cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0),
        cluster_score_from_blocks = dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0),
        top_block_score = dplyr::coalesce(as.numeric(top_block_score), 0),
        other_block_score = dplyr::coalesce(as.numeric(other_block_score), 0),
        cluster_gene_bonus = dplyr::coalesce(as.numeric(cluster_gene_bonus), 0),
        cluster_fragmentation_bonus = dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        top_gwas_hit_score = dplyr::coalesce(as.numeric(top_gwas_hit_score), 0),
        n_apps_supported = dplyr::coalesce(as.integer(n_apps_supported), 0L),
        priority_class = dplyr::coalesce(as.character(priority_class), ""),
        priority_class_relative = dplyr::coalesce(as.character(priority_class_relative), "")
      ) %>%
      dplyr::left_join(blk_sum, by = "cluster_id") %>%
      dplyr::mutate(
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        top_block_score_in_cluster = dplyr::coalesce(as.numeric(top_block_score_in_cluster), 0),
        sum_block_score_in_cluster = dplyr::coalesce(as.numeric(sum_block_score_in_cluster), 0),
        top_block_id = dplyr::coalesce(as.character(top_block_id), ""),
        blocks_in_cluster = dplyr::coalesce(as.character(blocks_in_cluster), ""),
        genes_in_blocks = dplyr::coalesce(as.character(genes_in_blocks), ""),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        cluster_block_ratio = dplyr::if_else(
          is.finite(cluster_score) & cluster_score > 0,
          top_block_score_in_cluster / cluster_score,
          NA_real_
        ),
        audit_flag = dplyr::case_when(
          cluster_score >= 25 & top_block_score_in_cluster < 10 ~ "CLUSTER_STRONG_WITHOUT_STRONG_BLOCKS",
          is.finite(cluster_block_ratio) & cluster_block_ratio >= 0.80 ~ "CLUSTER_DOMINATED_BY_ONE_BLOCK",
          n_blocks >= 4 & top_block_score_in_cluster < 12 ~ "MULTI_BLOCK_DIFFUSE_SIGNAL",
          TRUE ~ "CONSISTENT"
        ),
        audit_comment = dplyr::case_when(
          audit_flag == "CLUSTER_STRONG_WITHOUT_STRONG_BLOCKS" ~ "High cluster score but no clearly strong supporting block.",
          audit_flag == "CLUSTER_DOMINATED_BY_ONE_BLOCK" ~ "Cluster score is largely driven by a single block.",
          audit_flag == "MULTI_BLOCK_DIFFUSE_SIGNAL" ~ "Cluster contains several blocks but no clearly dominant high-scoring block.",
          TRUE ~ "Cluster-level and block-level signals are broadly coherent."
        )
      ) %>%
      dplyr::arrange(dplyr::desc(cluster_score), cluster_id)
  })
  
  # ============================================================
  # GLOBAL AUDIT SUMMARY
  # ============================================================
  priority_audit_summary_df <- reactive({
    cb <- cluster_block_audit_df()
    bg <- block_gene_audit_df()
    gen <- prioritized_gene_audit_base_df()
    blk <- prioritized_block_df_v2()
    gha <- gene_gwas_hit_score_audit_df()
    
    validate(
      need(is.data.frame(cb) && nrow(cb) > 0, "No cluster-block audit available."),
      need(is.data.frame(bg) && nrow(bg) > 0, "No block-gene audit available.")
    )
    
    gene_best_block <- blk %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
      ) %>%
      tidyr::separate_rows(genes_in_block, sep = "\\s*;\\s*") %>%
      dplyr::rename(gene = genes_in_block) %>%
      dplyr::mutate(gene = trimws(as.character(gene))) %>%
      dplyr::filter(!is.na(gene), nzchar(gene)) %>%
      dplyr::group_by(cluster_id, gene) %>%
      dplyr::summarise(
        best_block_score = suppressWarnings(max(block_score, na.rm = TRUE)),
        best_block_id = block_id[which.max(dplyr::coalesce(block_score, 0))[1]],
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        best_block_score = dplyr::if_else(is.finite(best_block_score), best_block_score, 0),
        best_block_id = dplyr::coalesce(as.character(best_block_id), "")
      )
    
    gene_summary <- gen %>%
      dplyr::left_join(gene_best_block, by = c("cluster_id", "gene")) %>%
      dplyr::mutate(
        best_block_score = dplyr::coalesce(as.numeric(best_block_score), 0),
        best_block_id = dplyr::coalesce(as.character(best_block_id), ""),
        audit_flag = dplyr::case_when(
          gene_score >= 15 & best_block_score < 10 ~ "GENE_STRONG_IN_WEAK_BLOCK",
          gene_score < 3 & best_block_score >= 20 ~ "BLOCK_STRONG_GENE_WEAK",
          TRUE ~ "CONSISTENT"
        ),
        audit_comment = dplyr::case_when(
          audit_flag == "GENE_STRONG_IN_WEAK_BLOCK" ~ "Strong gene score but weak supporting block.",
          audit_flag == "BLOCK_STRONG_GENE_WEAK" ~ "Gene is weak relative to its best overlapping block.",
          TRUE ~ "Gene-level and block-level signals are broadly coherent."
        )
      ) %>%
      dplyr::transmute(
        entity_type = "gene",
        cluster_id,
        block_id = best_block_id,
        entity_id = gene,
        entity_score = round(gene_score, 2),
        supporting_score = round(best_block_score, 2),
        audit_flag,
        audit_comment
      )
    
    block_summary <- bg %>%
      dplyr::transmute(
        entity_type = "block",
        cluster_id,
        block_id,
        entity_id = block_id,
        entity_score = round(block_score, 2),
        supporting_score = round(top_gene_score_in_block, 2),
        audit_flag,
        audit_comment
      )
    
    cluster_summary <- cb %>%
      dplyr::transmute(
        entity_type = "cluster",
        cluster_id,
        block_id = top_block_id,
        entity_id = cluster_id,
        entity_score = round(cluster_score, 2),
        supporting_score = round(top_block_score_in_cluster, 2),
        audit_flag,
        audit_comment
      )
    
    hit_summary <- if (is.data.frame(gha) && nrow(gha) > 0) {
      
      find_col <- function(cands) {
        hit <- intersect(cands, names(gha))
        if (length(hit)) hit[1] else NULL
      }
      
      gene_col <- find_col(c("gene", "gene_name"))
      gwas_col <- find_col(c("gwas_hit", "id_hit", "matched_ids"))
      cluster_col <- find_col(c("cluster_id"))
      block_col <- find_col(c("block_id"))
      gwas_score_col <- find_col(c("gwas_hit_priority_score", "gwas_hit_score", "priority_score"))
      gene_component_col <- find_col(c("gene_score_component", "score_component"))
      app_col <- find_col(c("source_app", "score_app"))
      link_col <- find_col(c("link_state"))
      
      validate(
        need(!is.null(gene_col), "gene column not available in gene_gwas_hit_score_audit_df()."),
        need(!is.null(gwas_col), "gwas_hit column not available in gene_gwas_hit_score_audit_df()."),
        need(!is.null(cluster_col), "cluster_id column not available in gene_gwas_hit_score_audit_df()."),
        need(!is.null(gwas_score_col), paste0(
          "No usable GWAS-hit score column found in gene_gwas_hit_score_audit_df(). Available columns: ",
          paste(names(gha), collapse = ", ")
        )),
        need(!is.null(gene_component_col), paste0(
          "No usable gene-score component column found in gene_gwas_hit_score_audit_df(). Available columns: ",
          paste(names(gha), collapse = ", ")
        ))
      )
      
      gha %>%
        dplyr::transmute(
          cluster_id = as.character(.data[[cluster_col]]),
          block_id = if (!is.null(block_col)) as.character(.data[[block_col]]) else "",
          gene = as.character(.data[[gene_col]]),
          gwas_hit = as.character(.data[[gwas_col]]),
          gene_score_component = dplyr::coalesce(as.numeric(.data[[gene_component_col]]), 0),
          gwas_hit_priority_score = dplyr::coalesce(as.numeric(.data[[gwas_score_col]]), 0),
          source_app = if (!is.null(app_col)) as.character(.data[[app_col]]) else "",
          link_state = if (!is.null(link_col)) as.character(.data[[link_col]]) else ""
        ) %>%
        dplyr::mutate(
          cluster_id = trimws(cluster_id),
          block_id = trimws(block_id),
          gene = trimws(gene),
          gwas_hit = trimws(gwas_hit),
          source_app = trimws(source_app),
          link_state = trimws(link_state),
          hit_support_score = gene_score_component * gwas_hit_priority_score
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(gene), nzchar(gene),
          !is.na(gwas_hit), nzchar(gwas_hit)
        ) %>%
        dplyr::mutate(
          audit_flag = dplyr::case_when(
            hit_support_score >= 10 & gwas_hit_priority_score < 4 ~ "HIT_SUPPORT_HIGH_WITH_WEAK_GWAS",
            hit_support_score <= 1 & gwas_hit_priority_score >= 10 ~ "GWAS_STRONG_WITH_WEAK_HIT_SUPPORT",
            TRUE ~ "CONSISTENT"
          ),
          audit_comment = dplyr::case_when(
            audit_flag == "HIT_SUPPORT_HIGH_WITH_WEAK_GWAS" ~ "Hit-level support is high relative to the underlying GWAS-hit priority score.",
            audit_flag == "GWAS_STRONG_WITH_WEAK_HIT_SUPPORT" ~ "Strong GWAS-hit priority but weak contribution to the gene support score.",
            TRUE ~ "Gene-hit contribution is broadly coherent."
          ),
          entity_id = paste0(gene, " ← ", gwas_hit)
        ) %>%
        dplyr::transmute(
          entity_type = "hit",
          cluster_id,
          block_id,
          entity_id,
          entity_score = round(hit_support_score, 2),
          supporting_score = round(gwas_hit_priority_score, 2),
          audit_flag,
          audit_comment
        )
      
    } else {
      tibble::tibble(
        entity_type = character(),
        cluster_id = character(),
        block_id = character(),
        entity_id = character(),
        entity_score = numeric(),
        supporting_score = numeric(),
        audit_flag = character(),
        audit_comment = character()
      )
    }
    
    dplyr::bind_rows(
      cluster_summary,
      block_summary,
      gene_summary,
      hit_summary
    ) %>%
      dplyr::arrange(
        factor(entity_type, levels = c("cluster", "block", "gene", "hit")),
        dplyr::desc(entity_score),
        cluster_id,
        entity_id
      )
  })
  
  ### Global audit summary:
  
  output$priority_audit_summary_dt <- DT::renderDT({
    df <- priority_audit_summary_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No global priority audit summary available.")
    )
    
    show_df <- df %>%
      dplyr::mutate(
        entity_type = dplyr::coalesce(as.character(entity_type), ""),
        cluster_id = dplyr::coalesce(as.character(cluster_id), ""),
        block_id = dplyr::coalesce(as.character(block_id), ""),
        entity_id = dplyr::coalesce(as.character(entity_id), ""),
        entity_score = round(dplyr::coalesce(as.numeric(entity_score), 0), 2),
        supporting_score = round(dplyr::coalesce(as.numeric(supporting_score), 0), 2),
        audit_flag = dplyr::coalesce(as.character(audit_flag), ""),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        audit_flag_badge = make_priority_badge(audit_flag)
      ) %>%
      dplyr::select(
        entity_type,
        cluster_id,
        block_id,
        entity_id,
        entity_score,
        supporting_score,
        audit_flag = audit_flag_badge,
        audit_comment
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5")
      )
    )
  }, server = FALSE)
  
  # audir cluster - block
  
  output$cluster_block_audit_dt <- DT::renderDT({
    df <- cluster_block_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster-block audit available.")
    )
    
    show_df <- df %>%
      dplyr::mutate(
        cluster_score = round(dplyr::coalesce(as.numeric(cluster_score), 0), 2),
        cluster_score_from_blocks = round(dplyr::coalesce(as.numeric(cluster_score_from_blocks), 0), 2),
        top_block_score_in_cluster = round(dplyr::coalesce(as.numeric(top_block_score_in_cluster), 0), 2),
        sum_block_score_in_cluster = round(dplyr::coalesce(as.numeric(sum_block_score_in_cluster), 0), 2),
        cluster_gene_bonus = round(dplyr::coalesce(as.numeric(cluster_gene_bonus), 0), 2),
        cluster_fragmentation_bonus = round(dplyr::coalesce(as.numeric(cluster_fragmentation_bonus), 0), 2),
        cluster_block_ratio = round(dplyr::coalesce(as.numeric(cluster_block_ratio), 0), 2),
        genes_in_blocks = make_genecards_links(genes_in_blocks),
        audit_flag_badge = make_priority_badge(audit_flag)
      ) %>%
      dplyr::select(
        cluster_id,
        cluster_score,
        top_block_id,
        top_block_score_in_cluster,
        sum_block_score_in_cluster,
        n_blocks,
        n_unique_genes_in_blocks,
        cluster_gene_bonus,
        cluster_fragmentation_bonus,
        cluster_block_ratio,
        audit_flag = audit_flag_badge,
        audit_comment,
        genes_in_blocks
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5")
      )
    )
  }, server = FALSE)
  
  # audit block - gene
  
  output$block_gene_audit_dt <- DT::renderDT({
    df <- block_gene_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No block-gene audit available.")
    )
    
    show_df <- df %>%
      dplyr::mutate(
        block_score = round(dplyr::coalesce(as.numeric(block_score), 0), 2),
        block_score_from_hits = round(dplyr::coalesce(as.numeric(block_score_from_hits), 0), 2),
        top_gwas_hit_score = round(dplyr::coalesce(as.numeric(top_gwas_hit_score), 0), 2),
        top_gene_score_in_block = round(dplyr::coalesce(as.numeric(top_gene_score_in_block), 0), 2),
        sum_gene_score_in_block = round(dplyr::coalesce(as.numeric(sum_gene_score_in_block), 0), 2),
        block_gene_ratio = round(dplyr::coalesce(as.numeric(block_gene_ratio), 0), 2),
        genes_in_block = make_genecards_links(genes_in_block),
        top_gene_in_block = make_genecards_links(top_gene_in_block),
        audit_flag_badge = make_priority_badge(audit_flag)
      ) %>%
      dplyr::select(
        cluster_id,
        block_id,
        block_score,
        top_gwas_hit_score,
        n_gwas_hits,
        n_genes_in_block,
        n_prioritized_genes_in_block,
        top_gene_in_block,
        top_gene_score_in_block,
        sum_gene_score_in_block,
        block_gene_ratio,
        audit_flag = audit_flag_badge,
        audit_comment,
        genes_in_block
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5")
      )
    )
  }, server = FALSE)
  
  gene_hit_audit_df <- reactive({
    gha <- gene_gwas_hit_score_audit_df()
    
    validate(
      need(is.data.frame(gha) && nrow(gha) > 0, "No gene-hit audit data available.")
    )
    
    find_col <- function(cands) {
      hit <- intersect(cands, names(gha))
      if (length(hit)) hit[1] else NULL
    }
    
    gene_col <- find_col(c("gene", "gene_name"))
    gwas_col <- find_col(c("gwas_hit", "id_hit", "matched_ids"))
    cluster_col <- find_col(c("cluster_id"))
    block_col <- find_col(c("block_id"))
    gwas_score_col <- find_col(c(
      "gwas_hit_priority_score",
      "gwas_hit_score",
      "priority_score"
    ))
    gene_component_col <- find_col(c(
      "gene_score_component",
      "score_component"
    ))
    app_col <- find_col(c("source_app", "score_app"))
    link_col <- find_col(c("link_state"))
    
    validate(
      need(!is.null(gene_col), "gene column not available in gene_gwas_hit_score_audit_df()."),
      need(!is.null(gwas_col), "gwas_hit column not available in gene_gwas_hit_score_audit_df()."),
      need(!is.null(cluster_col), "cluster_id column not available in gene_gwas_hit_score_audit_df()."),
      need(!is.null(gwas_score_col), paste0(
        "No usable GWAS-hit score column found in gene_gwas_hit_score_audit_df(). Available columns: ",
        paste(names(gha), collapse = ", ")
      )),
      need(!is.null(gene_component_col), paste0(
        "No usable gene-score component column found in gene_gwas_hit_score_audit_df(). Available columns: ",
        paste(names(gha), collapse = ", ")
      ))
    )
    
    gha %>%
      dplyr::transmute(
        cluster_id = as.character(.data[[cluster_col]]),
        block_id = if (!is.null(block_col)) as.character(.data[[block_col]]) else "",
        gene = as.character(.data[[gene_col]]),
        gwas_hit = as.character(.data[[gwas_col]]),
        gene_score_component = dplyr::coalesce(as.numeric(.data[[gene_component_col]]), 0),
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(.data[[gwas_score_col]]), 0),
        source_app = if (!is.null(app_col)) as.character(.data[[app_col]]) else "",
        link_state = if (!is.null(link_col)) as.character(.data[[link_col]]) else ""
      ) %>%
      dplyr::mutate(
        cluster_id = trimws(cluster_id),
        block_id = trimws(block_id),
        gene = trimws(gene),
        gwas_hit = trimws(gwas_hit),
        source_app = trimws(source_app),
        link_state = trimws(link_state)
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(gene), nzchar(gene),
        !is.na(gwas_hit), nzchar(gwas_hit)
      ) %>%
      dplyr::mutate(
        hit_support_score = gene_score_component * gwas_hit_priority_score,
        hit_vs_gwas_ratio = dplyr::if_else(
          is.finite(gwas_hit_priority_score) & gwas_hit_priority_score > 0,
          hit_support_score / gwas_hit_priority_score,
          NA_real_
        ),
        audit_flag = dplyr::case_when(
          hit_support_score >= 10 & gwas_hit_priority_score < 4 ~ "HIT_SUPPORT_HIGH_WITH_WEAK_GWAS",
          hit_support_score <= 1 & gwas_hit_priority_score >= 10 ~ "GWAS_STRONG_WITH_WEAK_HIT_SUPPORT",
          link_state %in% c("NOHIT", "NOLINK") ~ "LINK_STATE_REVIEW",
          TRUE ~ "CONSISTENT"
        ),
        audit_comment = dplyr::case_when(
          audit_flag == "HIT_SUPPORT_HIGH_WITH_WEAK_GWAS" ~ "Hit-level support is high relative to the GWAS-hit priority score.",
          audit_flag == "GWAS_STRONG_WITH_WEAK_HIT_SUPPORT" ~ "Strong GWAS-hit priority but weak contribution to the gene support score.",
          audit_flag == "LINK_STATE_REVIEW" ~ "The gene-hit contribution should be reviewed in light of the link state.",
          TRUE ~ "Gene-hit contribution is broadly coherent."
        )
      ) %>%
      dplyr::arrange(
        dplyr::desc(hit_support_score),
        dplyr::desc(gwas_hit_priority_score),
        cluster_id,
        gene,
        gwas_hit
      )
  })
  
  output$gene_hit_audit_dt <- DT::renderDT({
    df <- gene_hit_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No gene-hit audit available.")
    )
    
    show_df <- df %>%
      dplyr::mutate(
        hit_support_score = round(dplyr::coalesce(as.numeric(hit_support_score), 0), 2),
        gwas_hit_priority_score = round(dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0), 2),
        hit_vs_gwas_ratio = round(dplyr::coalesce(as.numeric(hit_vs_gwas_ratio), 0), 2),
        source_app = dplyr::coalesce(as.character(source_app), ""),
        link_state = dplyr::coalesce(as.character(link_state), ""),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        audit_flag_badge = make_priority_badge(audit_flag),
        gene = make_genecards_links(gene)
      ) %>%
      dplyr::select(
        cluster_id,
        block_id,
        gene,
        gwas_hit,
        source_app,
        link_state,
        hit_support_score,
        gwas_hit_priority_score,
        hit_vs_gwas_ratio,
        audit_flag = audit_flag_badge,
        audit_comment
      )
    
    DT::datatable(
      show_df,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list("copyHtml5", "csvHtml5", "excelHtml5")
      )
    )
  }, server = FALSE)
  
  #--------------------------- Audit plots -------------------------------------
  
  observeEvent(input$show_priority_audit_plots, {
    showModal(
      modalDialog(
        title = "Priority audit plots",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$div(
          style = "font-size:13px; color:#444; margin-bottom:10px;",
          HTML(
            paste0(
              "These plots help interpret coherence between hierarchical priority levels using priority classes. ",
              "The Sankey hierarchy flow is available as a dedicated tab within <b>Priority audit</b>."
            )
          )
        ),
        
        tabsetPanel(
          id = "priority_audit_plot_tabs",
          
          tabPanel(
            title = "Cluster ↔ Block",
            tags$br(),
            fluidRow(
              column(
                4,
                radioButtons(
                  "cluster_block_audit_class_mode",
                  "Priority class mode",
                  choices = c("ABS", "REL"),
                  selected = "ABS",
                  inline = TRUE
                )
              )
            ),
            plotly::plotlyOutput("cluster_block_audit_class_plot", height = "650px")
          ),
          
          tabPanel(
            title = "Block ↔ Gene",
            tags$br(),
            fluidRow(
              column(
                4,
                radioButtons(
                  "block_gene_audit_class_mode",
                  "Priority class mode",
                  choices = c("ABS", "REL"),
                  selected = "ABS",
                  inline = TRUE
                )
              )
            ),
            plotly::plotlyOutput("block_gene_audit_class_plot", height = "650px")
          )
        )
      )
    )
  })
  
  cluster_block_audit_plot_df <- reactive({
    df <- cluster_block_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster-block audit data available.")
    )
    
    df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0),
        top_block_score_in_cluster = dplyr::coalesce(as.numeric(top_block_score_in_cluster), 0),
        sum_block_score_in_cluster = dplyr::coalesce(as.numeric(sum_block_score_in_cluster), 0),
        cluster_block_ratio = dplyr::coalesce(as.numeric(cluster_block_ratio), 0),
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        audit_flag = dplyr::coalesce(as.character(audit_flag), "CONSISTENT"),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        top_block_id = dplyr::coalesce(as.character(top_block_id), ""),
        hover_txt = paste0(
          "Cluster: ", cluster_id,
          "<br>Cluster score: ", round(cluster_score, 2),
          "<br>Top block score: ", round(top_block_score_in_cluster, 2),
          "<br>Sum block score: ", round(sum_block_score_in_cluster, 2),
          "<br>Cluster/block ratio: ", round(cluster_block_ratio, 2),
          "<br>n blocks: ", n_blocks,
          "<br>n genes in blocks: ", n_unique_genes_in_blocks,
          "<br>Top block: ", top_block_id,
          "<br>Audit flag: ", audit_flag,
          ifelse(nzchar(audit_comment), paste0("<br>Comment: ", audit_comment), "")
        )
      )
  })
  
  output$cluster_block_audit_scatter <- plotly::renderPlotly({
    df <- cluster_block_audit_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster-block audit plot data available.")
    )
    
    plotly::plot_ly(
      data = df,
      x = ~cluster_score,
      y = ~top_block_score_in_cluster,
      type = "scatter",
      mode = "markers",
      color = ~audit_flag,
      colors = c(
        "CONSISTENT" = "darkgreen",
        "CLUSTER_STRONG_WITHOUT_STRONG_BLOCKS" = "#6baed6",
        "CLUSTER_DOMINATED_BY_ONE_BLOCK" = "orange",
        "MULTI_BLOCK_DIFFUSE_SIGNAL" = "#f16913"
      ),
      size = ~n_blocks,
      sizes = c(8, 28),
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Cluster score",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Top block score in cluster",
          automargin = TRUE
        ),
        legend = list(
          title = list(text = "Audit flag"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 70, r = 20, t = 20, b = 60)
      )
  })
  
  block_gene_audit_plot_df <- reactive({
    df <- block_gene_audit_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No block-gene audit data available.")
    )
    
    df %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        top_gene_score_in_block = dplyr::coalesce(as.numeric(top_gene_score_in_block), 0),
        sum_gene_score_in_block = dplyr::coalesce(as.numeric(sum_gene_score_in_block), 0),
        block_gene_ratio = dplyr::coalesce(as.numeric(block_gene_ratio), 0),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
        n_prioritized_genes_in_block = dplyr::coalesce(as.integer(n_prioritized_genes_in_block), 0L),
        top_gene_in_block = dplyr::coalesce(as.character(top_gene_in_block), ""),
        audit_flag = dplyr::coalesce(as.character(audit_flag), "CONSISTENT"),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        hover_txt = paste0(
          "Block: ", block_id,
          "<br>Cluster: ", cluster_id,
          "<br>Block score: ", round(block_score, 2),
          "<br>Top gene score in block: ", round(top_gene_score_in_block, 2),
          "<br>Sum gene score in block: ", round(sum_gene_score_in_block, 2),
          "<br>Block/gene ratio: ", round(block_gene_ratio, 2),
          "<br>n GWAS hits: ", n_gwas_hits,
          "<br>n genes in block: ", n_genes_in_block,
          "<br>n prioritized genes in block: ", n_prioritized_genes_in_block,
          "<br>Top gene: ", top_gene_in_block,
          "<br>Audit flag: ", audit_flag,
          ifelse(nzchar(audit_comment), paste0("<br>Comment: ", audit_comment), "")
        )
      )
  })
  
  output$block_gene_audit_scatter <- plotly::renderPlotly({
    df <- block_gene_audit_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No block-gene audit plot data available.")
    )
    
    plotly::plot_ly(
      data = df,
      x = ~block_score,
      y = ~top_gene_score_in_block,
      type = "scatter",
      mode = "markers",
      color = ~audit_flag,
      colors = c(
        "CONSISTENT" = "darkgreen",
        "BLOCK_DOMINANT_WITH_WEAK_GENES" = "#6baed6",
        "GENE_STRONG_IN_WEAK_BLOCK" = "orange",
        "TOP_BLOCK_WITHOUT_GENE_CONTENT" = "#f16913"
      ),
      size = ~n_gwas_hits,
      sizes = c(8, 28),
      hovertext = ~hover_txt,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Block score",
          automargin = TRUE
        ),
        yaxis = list(
          title = "Top gene score in block",
          automargin = TRUE
        ),
        legend = list(
          title = list(text = "Audit flag"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 70, r = 20, t = 20, b = 60)
      )
  })
  
  priority_audit_sankey_data <- reactive({
    cl  <- prioritized_cluster_df_v2()
    blk <- prioritized_block_df_v2()
    gen <- prioritized_gene_audit_base_df()
    gha <- gene_hit_audit_df()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available for Sankey."),
      need(is.data.frame(blk) && nrow(blk) > 0, "No prioritized blocks available for Sankey."),
      need(is.data.frame(gen) && nrow(gen) > 0, "No prioritized genes available for Sankey."),
      need(is.data.frame(gha) && nrow(gha) > 0, "No gene-hit audit data available for Sankey.")
    )
    
    top_n_clusters_raw <- input$priority_audit_sankey_top_n_clusters %||% 15
    top_n_genes_raw <- input$priority_audit_sankey_top_n_genes_per_block %||% 4
    top_n_hits_raw <- input$priority_audit_sankey_top_n_hits_per_gene %||% 4
    
    top_n_clusters <- if (identical(top_n_clusters_raw, "All")) Inf else suppressWarnings(as.integer(top_n_clusters_raw))
    top_n_genes_per_block <- if (identical(top_n_genes_raw, "All")) Inf else suppressWarnings(as.integer(top_n_genes_raw))
    top_n_hits_per_gene <- if (identical(top_n_hits_raw, "All")) Inf else suppressWarnings(as.integer(top_n_hits_raw))
    
    if (!is.finite(top_n_clusters)) top_n_clusters <- nrow(cl)
    if (!is.finite(top_n_genes_per_block)) top_n_genes_per_block <- 1e9
    if (!is.finite(top_n_hits_per_gene)) top_n_hits_per_gene <- 1e9
    
    # ------------------------------------------------------------
    # 1) Top clusters
    # ------------------------------------------------------------
    cl_top <- cl %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0)
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::arrange(dplyr::desc(cluster_score), cluster_id)
    
    if (is.finite(top_n_clusters)) {
      cl_top <- cl_top %>% dplyr::slice_head(n = top_n_clusters)
    }
    
    validate(
      need(nrow(cl_top) > 0, "No clusters available after top-cluster filtering.")
    )
    
    # ------------------------------------------------------------
    # 2) Cluster -> Block
    # ------------------------------------------------------------
    block_edges <- blk %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
      ) %>%
      dplyr::semi_join(cl_top, by = "cluster_id") %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(block_id), nzchar(block_id),
        is.finite(block_score), block_score > 0
      ) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
    
    validate(
      need(nrow(block_edges) > 0, "No cluster → block relationships available for Sankey.")
    )
    
    # ------------------------------------------------------------
    # 3) Block -> Gene
    # ------------------------------------------------------------
    block_gene_candidates <- block_edges %>%
      dplyr::select(cluster_id, block_id, genes_in_block) %>%
      tidyr::separate_rows(genes_in_block, sep = "\\s*;\\s*") %>%
      dplyr::rename(gene = genes_in_block) %>%
      dplyr::mutate(
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene))
    
    block_gene_edges <- block_gene_candidates %>%
      dplyr::left_join(
        gen %>%
          dplyr::transmute(
            cluster_id = as.character(cluster_id),
            gene = as.character(gene),
            gene_score = dplyr::coalesce(as.numeric(gene_score), 0)
          ),
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::mutate(
        gene_score = dplyr::coalesce(gene_score, 0)
      ) %>%
      dplyr::group_by(cluster_id, block_id) %>%
      dplyr::arrange(dplyr::desc(gene_score), gene, .by_group = TRUE)
    
    if (is.finite(top_n_genes_per_block)) {
      block_gene_edges <- block_gene_edges %>% dplyr::slice_head(n = top_n_genes_per_block)
    }
    
    block_gene_edges <- block_gene_edges %>%
      dplyr::ungroup() %>%
      dplyr::filter(gene_score > 0) %>%
      dplyr::distinct(cluster_id, block_id, gene, .keep_all = TRUE)
    
    validate(
      need(nrow(block_gene_edges) > 0, "No block → gene relationships available for Sankey.")
    )
    
    # ------------------------------------------------------------
    # 4) Gene -> Hit
    # ------------------------------------------------------------
    gene_hit_edges <- gha %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = dplyr::coalesce(as.character(block_id), ""),
        gene = as.character(gene),
        gwas_hit = as.character(gwas_hit),
        hit_support_score = dplyr::coalesce(as.numeric(hit_support_score), 0),
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0)
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(gene), nzchar(gene),
        !is.na(gwas_hit), nzchar(gwas_hit),
        is.finite(hit_support_score), hit_support_score > 0
      ) %>%
      dplyr::semi_join(
        block_gene_edges %>% dplyr::select(cluster_id, gene) %>% dplyr::distinct(),
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::group_by(cluster_id, gene) %>%
      dplyr::arrange(
        dplyr::desc(hit_support_score),
        dplyr::desc(gwas_hit_priority_score),
        gwas_hit,
        .by_group = TRUE
      )
    
    if (is.finite(top_n_hits_per_gene)) {
      gene_hit_edges <- gene_hit_edges %>% dplyr::slice_head(n = top_n_hits_per_gene)
    }
    
    gene_hit_edges <- gene_hit_edges %>%
      dplyr::ungroup() %>%
      dplyr::distinct(cluster_id, gene, gwas_hit, .keep_all = TRUE)
    
    validate(
      need(nrow(gene_hit_edges) > 0, "No gene → hit relationships available for Sankey.")
    )
    
    # ------------------------------------------------------------
    # 5) Lookup block per cluster + gene
    # ------------------------------------------------------------
    gene_block_lookup <- block_gene_edges %>%
      dplyr::group_by(cluster_id, gene) %>%
      dplyr::summarise(
        block_id_final = dplyr::first(block_id),
        .groups = "drop"
      )
    
    # ------------------------------------------------------------
    # 6) Nodes (ordered hierarchically)
    # ------------------------------------------------------------
    cluster_nodes <- cl_top %>%
      dplyr::transmute(
        cluster_id = cluster_id,
        block_id = "",
        gene = "",
        gwas_hit = "",
        node_key = paste0("cluster::", cluster_id),
        label = paste0(cluster_id, " [C]"),
        node_type = "Cluster",
        node_name = cluster_id,
        node_score = round(cluster_score, 2),
        node_extra = paste0("Cluster score: ", round(cluster_score, 2))
      ) %>%
      dplyr::distinct(node_key, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(node_score), cluster_id)
    
    block_nodes <- block_edges %>%
      dplyr::transmute(
        cluster_id = cluster_id,
        block_id = block_id,
        gene = "",
        gwas_hit = "",
        node_key = paste0("block::", cluster_id, "::", block_id),
        label = paste0(block_id, " [B]"),
        node_type = "Block",
        node_name = block_id,
        node_score = round(block_score, 2),
        node_extra = paste0(
          "Cluster: ", cluster_id,
          "<br>Block score: ", round(block_score, 2)
        )
      ) %>%
      dplyr::distinct(node_key, .keep_all = TRUE) %>%
      dplyr::arrange(cluster_id, dplyr::desc(node_score), block_id)
    
    gene_nodes <- block_gene_edges %>%
      dplyr::transmute(
        cluster_id = cluster_id,
        block_id = block_id,
        gene = gene,
        gwas_hit = "",
        node_key = paste0("gene::", cluster_id, "::", block_id, "::", gene),
        label = paste0(gene, " [G]"),
        node_type = "Gene",
        node_name = gene,
        node_score = round(gene_score, 2),
        node_extra = paste0(
          "Cluster: ", cluster_id,
          "<br>Block: ", block_id,
          "<br>Gene score: ", round(gene_score, 2)
        )
      ) %>%
      dplyr::distinct(node_key, .keep_all = TRUE) %>%
      dplyr::arrange(cluster_id, block_id, dplyr::desc(node_score), gene)
    
    hit_nodes <- gene_hit_edges %>%
      dplyr::left_join(
        gene_block_lookup,
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::filter(!is.na(block_id_final), nzchar(block_id_final)) %>%
      dplyr::transmute(
        cluster_id = cluster_id,
        block_id = block_id_final,
        gene = gene,
        gwas_hit = gwas_hit,
        node_key = paste0("hit::", cluster_id, "::", block_id_final, "::", gene, "::", gwas_hit),
        label = paste0(gwas_hit, " [H]"),
        node_type = "GWAS hit",
        node_name = gwas_hit,
        node_score = round(hit_support_score, 2),
        node_extra = paste0(
          "Cluster: ", cluster_id,
          "<br>Block: ", block_id_final,
          "<br>Gene: ", gene,
          "<br>Hit support score: ", round(hit_support_score, 2),
          "<br>GWAS-hit priority score: ", round(gwas_hit_priority_score, 2)
        )
      ) %>%
      dplyr::distinct(node_key, .keep_all = TRUE) %>%
      dplyr::arrange(cluster_id, block_id, gene, dplyr::desc(node_score), gwas_hit)
    
    nodes <- dplyr::bind_rows(
      cluster_nodes,
      block_nodes,
      gene_nodes,
      hit_nodes
    ) %>%
      dplyr::distinct(node_key, .keep_all = TRUE) %>%
      dplyr::mutate(node_id = dplyr::row_number() - 1L)
    
    validate(
      need(nrow(nodes) > 0, "No Sankey nodes available after hierarchy assembly.")
    )
    
    # ------------------------------------------------------------
    # 7) Hierarchical fixed node layout
    # ------------------------------------------------------------
    
    scale_positions <- function(n, y_min, y_max) {
      if (n <= 0) return(numeric(0))
      if (n == 1) return((y_min + y_max) / 2)
      seq(y_min, y_max, length.out = n)
    }
    
    # x fixed by level
    nodes <- nodes %>%
      dplyr::mutate(
        x_pos = dplyr::case_when(
          node_type == "Cluster" ~ 0.05,
          node_type == "Block" ~ 0.32,
          node_type == "Gene" ~ 0.60,
          node_type == "GWAS hit" ~ 0.88,
          TRUE ~ 0.50
        ),
        y_pos = NA_real_
      )
    
    # ------------------------------------------------------------
    # 7.1 clusters
    # ------------------------------------------------------------
    cluster_order <- cluster_nodes %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::pull(cluster_id)
    
    cluster_y <- stats::setNames(
      scale_positions(length(cluster_order), 0.06, 0.94),
      cluster_order
    )
    
    nodes$y_pos[nodes$node_type == "Cluster"] <- cluster_y[nodes$cluster_id[nodes$node_type == "Cluster"]]
    
    # ------------------------------------------------------------
    # 7.2 blocks within each cluster band
    # ------------------------------------------------------------
    cluster_half_span <- 0.055
    
    for (cid in cluster_order) {
      y_c <- cluster_y[[cid]]
      
      blk_sub <- block_nodes %>%
        dplyr::filter(cluster_id == cid) %>%
        dplyr::arrange(dplyr::desc(node_score), block_id)
      
      if (nrow(blk_sub) == 0) next
      
      y_blk <- scale_positions(
        nrow(blk_sub),
        max(0.01, y_c - cluster_half_span),
        min(0.99, y_c + cluster_half_span)
      )
      
      idx <- match(blk_sub$node_key, nodes$node_key)
      nodes$y_pos[idx] <- y_blk
    }
    
    # helper: lookup y of a node by key
    node_y_lookup <- stats::setNames(nodes$y_pos, nodes$node_key)
    
    # ------------------------------------------------------------
    # 7.3 genes within each block band
    # ------------------------------------------------------------
    block_half_span <- 0.030
    
    block_order_df <- block_nodes %>%
      dplyr::select(cluster_id, block_id, node_key) %>%
      dplyr::distinct()
    
    for (i in seq_len(nrow(block_order_df))) {
      cid <- block_order_df$cluster_id[i]
      bid <- block_order_df$block_id[i]
      bkey <- block_order_df$node_key[i]
      y_b <- node_y_lookup[[bkey]]
      
      gene_sub <- gene_nodes %>%
        dplyr::filter(cluster_id == cid, block_id == bid) %>%
        dplyr::arrange(dplyr::desc(node_score), gene)
      
      if (nrow(gene_sub) == 0) next
      
      y_gene <- scale_positions(
        nrow(gene_sub),
        max(0.01, y_b - block_half_span),
        min(0.99, y_b + block_half_span)
      )
      
      idx <- match(gene_sub$node_key, nodes$node_key)
      nodes$y_pos[idx] <- y_gene
    }
    
    node_y_lookup <- stats::setNames(nodes$y_pos, nodes$node_key)
    
    # ------------------------------------------------------------
    # 7.4 hits within each gene band
    # ------------------------------------------------------------
    gene_half_span <- 0.018
    
    gene_order_df <- gene_nodes %>%
      dplyr::select(cluster_id, block_id, gene, node_key) %>%
      dplyr::distinct()
    
    for (i in seq_len(nrow(gene_order_df))) {
      cid <- gene_order_df$cluster_id[i]
      bid <- gene_order_df$block_id[i]
      gid <- gene_order_df$gene[i]
      gkey <- gene_order_df$node_key[i]
      y_g <- node_y_lookup[[gkey]]
      
      hit_sub <- hit_nodes %>%
        dplyr::filter(cluster_id == cid, block_id == bid, gene == gid) %>%
        dplyr::arrange(dplyr::desc(node_score), gwas_hit)
      
      if (nrow(hit_sub) == 0) next
      
      y_hit <- scale_positions(
        nrow(hit_sub),
        max(0.01, y_g - gene_half_span),
        min(0.99, y_g + gene_half_span)
      )
      
      idx <- match(hit_sub$node_key, nodes$node_key)
      nodes$y_pos[idx] <- y_hit
    }
    
    # ------------------------------------------------------------
    # 8) Links
    # ------------------------------------------------------------
    links_cb <- block_edges %>%
      dplyr::transmute(
        source_key = paste0("cluster::", cluster_id),
        target_key = paste0("block::", cluster_id, "::", block_id),
        value = pmax(block_score, 0.01),
        edge_type = "Cluster → Block",
        edge_info = paste0(
          "Cluster: ", cluster_id,
          "<br>Block: ", block_id,
          "<br>Flow value: ", round(block_score, 2)
        )
      )
    
    links_bg <- block_gene_edges %>%
      dplyr::transmute(
        source_key = paste0("block::", cluster_id, "::", block_id),
        target_key = paste0("gene::", cluster_id, "::", block_id, "::", gene),
        value = pmax(gene_score, 0.01),
        edge_type = "Block → Gene",
        edge_info = paste0(
          "Block: ", block_id,
          "<br>Gene: ", gene,
          "<br>Flow value: ", round(gene_score, 2)
        )
      )
    
    links_gh <- gene_hit_edges %>%
      dplyr::left_join(
        gene_block_lookup,
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::filter(!is.na(block_id_final), nzchar(block_id_final)) %>%
      dplyr::transmute(
        source_key = paste0("gene::", cluster_id, "::", block_id_final, "::", gene),
        target_key = paste0("hit::", cluster_id, "::", block_id_final, "::", gene, "::", gwas_hit),
        value = pmax(hit_support_score, 0.01),
        edge_type = "Gene → Hit",
        edge_info = paste0(
          "Gene: ", gene,
          "<br>GWAS hit: ", gwas_hit,
          "<br>Flow value: ", round(hit_support_score, 2)
        )
      )
    
    links <- dplyr::bind_rows(
      links_cb,
      links_bg,
      links_gh
    ) %>%
      dplyr::left_join(
        nodes %>% dplyr::select(source_key = node_key, source = node_id),
        by = "source_key"
      ) %>%
      dplyr::left_join(
        nodes %>% dplyr::select(target_key = node_key, target = node_id),
        by = "target_key"
      ) %>%
      dplyr::filter(
        is.finite(source),
        is.finite(target),
        is.finite(value),
        value > 0
      )
    
    validate(
      need(nrow(links) > 0, "No Sankey links available after hierarchy assembly.")
    )
    
    list(
      nodes = nodes,
      links = links
    )
  })
  
  output$priority_audit_sankey <- plotly::renderPlotly({
    sank <- priority_audit_sankey_data()
    
    nodes <- sank$nodes
    links <- sank$links
    
    validate(
      need(is.data.frame(nodes) && nrow(nodes) > 0, "No Sankey nodes available."),
      need(is.data.frame(links) && nrow(links) > 0, "No Sankey links available.")
    )
    
    node_labels <- nodes$label
    
    plotly::plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = node_labels,
        pad = 14,
        thickness = 14,
        line = list(color = "black", width = 0.3)
      ),
      link = list(
        source = links$source,
        target = links$target,
        value = links$value
      )
    ) %>%
      plotly::layout(
        font = list(size = 11),
        margin = list(l = 20, r = 20, t = 20, b = 20)
      )
  })
  
  output$priority_audit_sankey_full <- plotly::renderPlotly({
    sank <- priority_audit_sankey_data()
    
    nodes <- sank$nodes
    links <- sank$links
    
    validate(
      need(is.data.frame(nodes) && nrow(nodes) > 0, "No Sankey nodes available."),
      need(is.data.frame(links) && nrow(links) > 0, "No Sankey links available.")
    )
    
    node_colors <- dplyr::case_when(
      nodes$node_type == "Cluster" ~ "#1f78b4",
      nodes$node_type == "Block" ~ "#ff7f00",
      nodes$node_type == "Gene" ~ "#33a02c",
      nodes$node_type == "GWAS hit" ~ "#6a3d9a",
      TRUE ~ "#999999"
    )
    
    plotly::plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes$label,
        color = node_colors,
        x = nodes$x_pos,
        y = nodes$y_pos,
        pad = 14,
        thickness = 16,
        line = list(color = "black", width = 0.3)
      ),
      link = list(
        source = links$source,
        target = links$target,
        value = links$value
      )
    ) %>%
      plotly::layout(
        font = list(size = 11),
        margin = list(l = 20, r = 20, t = 20, b = 20)
      )
  })
  
  cluster_block_audit_class_plot_df <- reactive({
    cb <- cluster_block_audit_df()
    blk <- prioritized_block_df_v2()
    
    validate(
      need(is.data.frame(cb) && nrow(cb) > 0, "No cluster-block audit data available."),
      need(is.data.frame(blk) && nrow(blk) > 0, "No prioritized block data available.")
    )
    
    mode <- input$priority_audit_class_mode %||% "ABS"
    
    class_col_cluster <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    class_col_block   <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    
    cluster_base <- cb %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        cluster_score = dplyr::coalesce(as.numeric(cluster_score), 0),
        n_blocks = dplyr::coalesce(as.integer(n_blocks), 0L),
        n_unique_genes_in_blocks = dplyr::coalesce(as.integer(n_unique_genes_in_blocks), 0L),
        audit_flag = dplyr::coalesce(as.character(audit_flag), "CONSISTENT"),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        cluster_priority_class = dplyr::coalesce(as.character(.data[[class_col_cluster]]), "")
      )
    
    block_base <- blk %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
        block_priority_class = dplyr::coalesce(as.character(.data[[class_col_block]]), "")
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(block_id), nzchar(block_id))
    
    block_base %>%
      dplyr::left_join(cluster_base, by = "cluster_id") %>%
      dplyr::mutate(
        cluster_priority_class = factor(cluster_priority_class, levels = c("Low", "Medium", "High")),
        block_priority_class = factor(block_priority_class, levels = c("Low", "Medium", "High")),
        hover_txt = paste0(
          "Cluster: ", cluster_id,
          "<br>Block: ", block_id,
          "<br>Cluster class (", mode, "): ", cluster_priority_class,
          "<br>Block class (", mode, "): ", block_priority_class,
          "<br>Cluster score: ", round(cluster_score, 2),
          "<br>Block score: ", round(block_score, 2),
          "<br>Cluster n blocks: ", n_blocks,
          "<br>Block n GWAS hits: ", n_gwas_hits,
          "<br>Block n genes: ", n_genes_in_block,
          "<br>Audit flag: ", audit_flag,
          ifelse(nzchar(audit_comment), paste0("<br>Comment: ", audit_comment), "")
        )
      ) %>%
      dplyr::filter(!is.na(cluster_priority_class), !is.na(block_priority_class))
  })
  
  output$cluster_block_audit_class_plot <- plotly::renderPlotly({
    df <- cluster_block_audit_class_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No cluster-block class plot data available.")
    )
    
    class_to_num <- function(x) {
      dplyr::case_when(
        as.character(x) == "Low" ~ 1,
        as.character(x) == "Medium" ~ 2,
        as.character(x) == "High" ~ 3,
        TRUE ~ NA_real_
      )
    }
    
    set.seed(123)
    
    df <- df %>%
      dplyr::mutate(
        x_num = class_to_num(cluster_priority_class),
        y_num = class_to_num(block_priority_class),
        x_jit = x_num + stats::runif(dplyr::n(), -0.18, 0.18),
        y_jit = y_num + stats::runif(dplyr::n(), -0.18, 0.18)
      ) %>%
      dplyr::filter(is.finite(x_jit), is.finite(y_jit))
    
    plotly::plot_ly(
      data = df,
      x = ~x_jit,
      y = ~y_jit,
      type = "scatter",
      mode = "markers",
      color = ~audit_flag,
      colors = c(
        "CONSISTENT" = "darkgreen",
        "CLUSTER_STRONG_WITHOUT_STRONG_BLOCKS" = "#6baed6",
        "CLUSTER_DOMINATED_BY_ONE_BLOCK" = "orange",
        "MULTI_BLOCK_DIFFUSE_SIGNAL" = "#f16913"
      ),
      size = ~sqrt(pmax(n_gwas_hits, 1)),
      sizes = c(16, 42),
      hovertext = ~hover_txt,
      hoverinfo = "text",
      marker = list(opacity = 0.70)
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Cluster priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Block priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        legend = list(
          title = list(text = "Audit flag"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 80, r = 20, t = 20, b = 60)
      )
  })
  
  block_gene_audit_class_plot_df <- reactive({
    bg <- block_gene_audit_df()
    gen <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(bg) && nrow(bg) > 0, "No block-gene audit data available."),
      need(is.data.frame(gen) && nrow(gen) > 0, "No prioritized gene data available.")
    )
    
    mode <- input$priority_audit_class_mode %||% "ABS"
    
    class_col_block <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    gene_class_col  <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    
    gene_col <- intersect(c("gene", "gene_name"), names(gen))
    score_col <- intersect(c("gene_score", "priority_score", "score"), names(gen))
    clusters_col <- intersect(c("cluster_ids", "cluster_id"), names(gen))
    
    validate(
      need(length(gene_col) > 0, "Gene column not available in prioritized_gene_df_v2()."),
      need(length(score_col) > 0, "Gene score column not available in prioritized_gene_df_v2()."),
      need(length(clusters_col) > 0, "cluster_ids / cluster_id column not available in prioritized_gene_df_v2()."),
      need(gene_class_col %in% names(gen), paste0(mode, " gene priority class column not available in prioritized_gene_df_v2()."))
    )
    
    gene_col <- gene_col[1]
    score_col <- score_col[1]
    clusters_col <- clusters_col[1]
    
    gene_map <- gen %>%
      dplyr::transmute(
        gene = as.character(.data[[gene_col]]),
        gene_score = dplyr::coalesce(as.numeric(.data[[score_col]]), 0),
        cluster_ids_raw = as.character(.data[[clusters_col]]),
        gene_priority_class = dplyr::coalesce(as.character(.data[[gene_class_col]]), "")
      ) %>%
      tidyr::separate_rows(cluster_ids_raw, sep = "\\s*;\\s*") %>%
      dplyr::rename(cluster_id = cluster_ids_raw) %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(gene), nzchar(gene)) %>%
      dplyr::distinct(cluster_id, gene, .keep_all = TRUE)
    
    block_gene_map <- bg %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_score = dplyr::coalesce(as.numeric(block_score), 0),
        n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
        n_prioritized_genes_in_block = dplyr::coalesce(as.integer(n_prioritized_genes_in_block), 0L),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), ""),
        audit_flag = dplyr::coalesce(as.character(audit_flag), "CONSISTENT"),
        audit_comment = dplyr::coalesce(as.character(audit_comment), ""),
        block_priority_class = dplyr::coalesce(as.character(.data[[class_col_block]]), "")
      ) %>%
      tidyr::separate_rows(genes_in_block, sep = "\\s*;\\s*") %>%
      dplyr::rename(gene = genes_in_block) %>%
      dplyr::mutate(
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(!is.na(gene), nzchar(gene))
    
    block_gene_map %>%
      dplyr::left_join(
        gene_map,
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::mutate(
        block_priority_class = factor(block_priority_class, levels = c("Low", "Medium", "High")),
        gene_priority_class = factor(gene_priority_class, levels = c("Low", "Medium", "High")),
        hover_txt = paste0(
          "Block: ", block_id,
          "<br>Cluster: ", cluster_id,
          "<br>Gene: ", gene,
          "<br>Block class (", mode, "): ", block_priority_class,
          "<br>Gene class (", mode, "): ", gene_priority_class,
          "<br>Block score: ", round(block_score, 2),
          "<br>Gene score: ", round(gene_score, 2),
          "<br>n GWAS hits in block: ", n_gwas_hits,
          "<br>n prioritized genes in block: ", n_prioritized_genes_in_block,
          "<br>Audit flag: ", audit_flag,
          ifelse(nzchar(audit_comment), paste0("<br>Comment: ", audit_comment), "")
        )
      ) %>%
      dplyr::filter(!is.na(block_priority_class), !is.na(gene_priority_class))
  })
  
  output$block_gene_audit_class_plot <- plotly::renderPlotly({
    df <- block_gene_audit_class_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No block-gene class plot data available.")
    )
    
    class_to_num <- function(x) {
      dplyr::case_when(
        as.character(x) == "Low" ~ 1,
        as.character(x) == "Medium" ~ 2,
        as.character(x) == "High" ~ 3,
        TRUE ~ NA_real_
      )
    }
    
    set.seed(123)
    
    df <- df %>%
      dplyr::mutate(
        x_num = class_to_num(block_priority_class),
        y_num = class_to_num(gene_priority_class),
        x_jit = x_num + stats::runif(dplyr::n(), -0.18, 0.18),
        y_jit = y_num + stats::runif(dplyr::n(), -0.18, 0.18)
      ) %>%
      dplyr::filter(is.finite(x_jit), is.finite(y_jit))
    
    plotly::plot_ly(
      data = df,
      x = ~x_jit,
      y = ~y_jit,
      type = "scatter",
      mode = "markers",
      color = ~audit_flag,
      colors = c(
        "CONSISTENT" = "darkgreen",
        "BLOCK_DOMINANT_WITH_WEAK_GENES" = "#6baed6",
        "GENE_STRONG_IN_WEAK_BLOCK" = "orange",
        "TOP_BLOCK_WITHOUT_GENE_CONTENT" = "#f16913"
      ),
      size = ~sqrt(pmax(n_gwas_hits, 1)),
      sizes = c(16, 42),
      hovertext = ~hover_txt,
      hoverinfo = "text",
      marker = list(opacity = 0.70)
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Block priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Gene priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        legend = list(
          title = list(text = "Audit flag"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 80, r = 20, t = 20, b = 60)
      )
  })
  
  gene_hit_audit_class_plot_df <- reactive({
    gha <- gene_hit_audit_df()
    gen <- prioritized_gene_df_v2()
    gh  <- gwas_hit_priority_df_v2()
    
    validate(
      need(is.data.frame(gha) && nrow(gha) > 0, "No gene-hit audit data available."),
      need(is.data.frame(gen) && nrow(gen) > 0, "No prioritized gene data available."),
      need(is.data.frame(gh) && nrow(gh) > 0, "No prioritized GWAS-hit data available.")
    )
    
    mode <- input$priority_audit_class_mode %||% "ABS"
    
    gene_class_col <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    hit_class_col  <- if (identical(mode, "REL")) "priority_class_relative" else "priority_class"
    
    gene_col <- intersect(c("gene", "gene_name"), names(gen))
    score_col <- intersect(c("gene_score", "priority_score", "score"), names(gen))
    clusters_col <- intersect(c("cluster_ids", "cluster_id"), names(gen))
    
    validate(
      need(length(gene_col) > 0, "Gene column not available in prioritized_gene_df_v2()."),
      need(length(score_col) > 0, "Gene score column not available in prioritized_gene_df_v2()."),
      need(length(clusters_col) > 0, "cluster_ids / cluster_id column not available in prioritized_gene_df_v2()."),
      need(gene_class_col %in% names(gen), paste0(mode, " gene priority class column not available in prioritized_gene_df_v2().")),
      need(hit_class_col %in% names(gh), paste0(mode, " hit priority class column not available in gwas_hit_priority_df_v2()."))
    )
    
    gene_col <- gene_col[1]
    score_col <- score_col[1]
    clusters_col <- clusters_col[1]
    
    gene_map <- gen %>%
      dplyr::transmute(
        gene = as.character(.data[[gene_col]]),
        gene_score = dplyr::coalesce(as.numeric(.data[[score_col]]), 0),
        cluster_ids_raw = as.character(.data[[clusters_col]]),
        gene_priority_class = dplyr::coalesce(as.character(.data[[gene_class_col]]), "")
      ) %>%
      tidyr::separate_rows(cluster_ids_raw, sep = "\\s*;\\s*") %>%
      dplyr::rename(cluster_id = cluster_ids_raw) %>%
      dplyr::mutate(
        cluster_id = trimws(as.character(cluster_id)),
        gene = trimws(as.character(gene))
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(gene), nzchar(gene)) %>%
      dplyr::distinct(cluster_id, gene, .keep_all = TRUE)
    
    hit_map <- gh %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        gwas_hit = as.character(gwas_hit),
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0),
        hit_priority_class = dplyr::coalesce(as.character(.data[[hit_class_col]]), "")
      ) %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(gwas_hit), nzchar(gwas_hit)) %>%
      dplyr::distinct(cluster_id, gwas_hit, .keep_all = TRUE)
    
    gha %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = dplyr::coalesce(as.character(block_id), ""),
        gene = as.character(gene),
        gwas_hit = as.character(gwas_hit),
        hit_support_score = dplyr::coalesce(as.numeric(hit_support_score), 0),
        gwas_hit_priority_score = dplyr::coalesce(as.numeric(gwas_hit_priority_score), 0),
        source_app = dplyr::coalesce(as.character(source_app), ""),
        link_state = dplyr::coalesce(as.character(link_state), ""),
        audit_flag = dplyr::coalesce(as.character(audit_flag), "CONSISTENT"),
        audit_comment = dplyr::coalesce(as.character(audit_comment), "")
      ) %>%
      dplyr::left_join(
        gene_map %>% dplyr::select(cluster_id, gene, gene_score, gene_priority_class),
        by = c("cluster_id", "gene")
      ) %>%
      dplyr::left_join(
        hit_map %>% dplyr::select(cluster_id, gwas_hit, hit_priority_class),
        by = c("cluster_id", "gwas_hit")
      ) %>%
      dplyr::mutate(
        gene_priority_class = factor(gene_priority_class, levels = c("Low", "Medium", "High")),
        hit_priority_class = factor(hit_priority_class, levels = c("Low", "Medium", "High")),
        hover_txt = paste0(
          "Gene: ", gene,
          "<br>GWAS hit: ", gwas_hit,
          "<br>Cluster: ", cluster_id,
          ifelse(nzchar(block_id), paste0("<br>Block: ", block_id), ""),
          "<br>Gene class (", mode, "): ", gene_priority_class,
          "<br>Hit class (", mode, "): ", hit_priority_class,
          "<br>Gene score: ", round(gene_score, 2),
          "<br>Hit support score: ", round(hit_support_score, 2),
          "<br>GWAS-hit priority score: ", round(gwas_hit_priority_score, 2),
          ifelse(nzchar(source_app), paste0("<br>Source app: ", source_app), ""),
          ifelse(nzchar(link_state), paste0("<br>Link state: ", link_state), ""),
          "<br>Audit flag: ", audit_flag,
          ifelse(nzchar(audit_comment), paste0("<br>Comment: ", audit_comment), "")
        )
      ) %>%
      dplyr::filter(!is.na(gene_priority_class), !is.na(hit_priority_class))
  })
  
  output$gene_hit_audit_class_plot <- plotly::renderPlotly({
    df <- gene_hit_audit_class_plot_df()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No gene-hit class plot data available.")
    )
    
    class_to_num <- function(x) {
      dplyr::case_when(
        as.character(x) == "Low" ~ 1,
        as.character(x) == "Medium" ~ 2,
        as.character(x) == "High" ~ 3,
        TRUE ~ NA_real_
      )
    }
    
    set.seed(123)
    
    df <- df %>%
      dplyr::mutate(
        x_num = class_to_num(gene_priority_class),
        y_num = class_to_num(hit_priority_class),
        x_jit = x_num + stats::runif(dplyr::n(), -0.18, 0.18),
        y_jit = y_num + stats::runif(dplyr::n(), -0.18, 0.18)
      ) %>%
      dplyr::filter(is.finite(x_jit), is.finite(y_jit))
    
    plotly::plot_ly(
      data = df,
      x = ~x_jit,
      y = ~y_jit,
      type = "scatter",
      mode = "markers",
      color = ~audit_flag,
      colors = c(
        CONSISTENT = "darkgreen",
        HIT_SUPPORT_HIGH_WITH_WEAK_GWAS = "#6baed6",
        GWAS_STRONG_WITH_WEAK_HIT_SUPPORT =  "orange",
        LINK_STATE_REVIEW = "#f16913"
      ),
      size = ~sqrt(pmax(gwas_hit_priority_score, 1)),
      sizes = c(16, 42),
      hovertext = ~hover_txt,
      hoverinfo = "text",
      marker = list(opacity = 0.75)
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Gene priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "GWAS-hit priority class",
          tickmode = "array",
          tickvals = c(1, 2, 3),
          ticktext = c("Low", "Medium", "High"),
          range = c(0.5, 3.5),
          zeroline = FALSE
        ),
        legend = list(
          title = list(text = "Audit flag"),
          orientation = "h",
          x = 0,
          y = 1.08
        ),
        margin = list(l = 80, r = 20, t = 20, b = 60)
      )
  })
  
  priority_audit_hierarchy_module_server(
    "audit_table_temp",
    data_reactive = reactive({
      gene_hit_audit_df()
    })
  )
  
  #############################################################
  ################ GO / KEGG ENRICH ##########################
  #############################################################
  
  # ============================================================
  # SERVER — Integrator Enrichment (GO / KEGG)
  # Catalog/GTEx-style adaptation for prioritized genes
  # Input:
  #   - prioritized_gene_df_v2()
  #   - prioritized_cluster_df_v2()
  # Expected gene column:
  #   - gene
  # Optional cluster columns in gene table:
  #   - cluster_id
  #   - cluster_ids   (semicolon separated)
  # ============================================================
  
  pick_col <- function(df, candidates) {
    nm <- names(df)
    hit <- candidates[candidates %in% nm]
    if (length(hit)) hit[1] else NULL
  }
  
  split_semicolon_unique <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character(0))
    
    vals <- unlist(strsplit(x, "\\s*;\\s*"))
    vals <- trimws(vals)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    unique(vals)
  }
  
  # ------------------------------------------------------------
  # Base reactives
  # ------------------------------------------------------------
  priority_genes_all <- reactive({
    df <- prioritized_gene_df_v2()
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "No prioritized genes available."),
      need("gene" %in% names(df), "prioritized_gene_df_v2() must contain column 'gene'.")
    )
    
    df %>%
      dplyr::mutate(
        gene = as.character(.data$gene)
      ) %>%
      dplyr::filter(!is.na(.data$gene), nzchar(.data$gene))
  })
  
  priority_clusters_all <- reactive({
    cl <- prioritized_cluster_df_v2()
    
    validate(
      need(is.data.frame(cl) && nrow(cl) > 0, "No prioritized clusters available."),
      need("cluster_id" %in% names(cl), "prioritized_cluster_df_v2() must contain column 'cluster_id'."),
      need("chr" %in% names(cl), "prioritized_cluster_df_v2() must contain column 'chr'."),
      need("start" %in% names(cl), "prioritized_cluster_df_v2() must contain column 'start'."),
      need("end" %in% names(cl), "prioritized_cluster_df_v2() must contain column 'end'.")
    )
    
    cl %>%
      dplyr::mutate(
        cluster_id = as.character(.data$cluster_id),
        chr = as.character(.data$chr),
        start = suppressWarnings(as.numeric(.data$start)),
        end = suppressWarnings(as.numeric(.data$end))
      ) %>%
      dplyr::filter(!is.na(.data$cluster_id), nzchar(.data$cluster_id)) %>%
      dplyr::distinct(.data$cluster_id, .keep_all = TRUE)
  })
  
  gene_cluster_long <- reactive({
    df <- priority_genes_all()
    
    cid_col <- pick_col(df, c("cluster_ids", "cluster_id"))
    validate(
      need(!is.null(cid_col), "Prioritized gene table must contain 'cluster_id' or 'cluster_ids'.")
    )
    
    out <- df %>%
      dplyr::mutate(
        cluster_raw = as.character(.data[[cid_col]])
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        cluster_id_expanded = list(split_semicolon_unique(cluster_raw))
      ) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(cluster_id_expanded) %>%
      dplyr::rename(cluster_id = cluster_id_expanded) %>%
      dplyr::filter(!is.na(.data$cluster_id), nzchar(.data$cluster_id)) %>%
      dplyr::distinct(.data$gene, .data$cluster_id, .keep_all = TRUE)
    
    validate(
      need(nrow(out) > 0, "No gene ↔ cluster mapping could be built from prioritized genes.")
    )
    
    out
  })
  
  # ------------------------------------------------------------
  # UI cluster picker (same style as GTEx/Catalog)
  # ------------------------------------------------------------
  output$func_cluster_ui <- renderUI({
    cl <- priority_clusters_all()
    
    chr_choices <- sort(unique(suppressWarnings(as.integer(cl$chr))))
    chr_choices <- chr_choices[is.finite(chr_choices)]
    validate(need(length(chr_choices) > 0, "No valid chromosomes available."))
    
    chr_labels <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("func_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("func_cluster_id_ui")
    )
  })
  
  output$func_cluster_id_ui <- renderUI({
    cl <- priority_clusters_all()
    req(is.data.frame(cl), nrow(cl) > 0, input$func_chr)
    
    x <- cl %>%
      dplyr::filter(suppressWarnings(as.integer(.data$chr)) == suppressWarnings(as.integer(input$func_chr))) %>%
      dplyr::arrange(.data$start, .data$end)
    
    validate(need(nrow(x) > 0, "No clusters on this chromosome."))
    
    lab <- paste0(x$cluster_id, " (", x$start, "-", x$end, ")")
    choices <- stats::setNames(x$cluster_id, lab)
    
    selectInput("func_cluster_id", "Cluster", choices = choices, selected = x$cluster_id[1])
  })
  
  output$func_gene_ui <- renderUI({
    dt <- tryCatch(priority_genes_all(), error = function(e) NULL)
    if (!is.data.frame(dt) || !nrow(dt)) return(NULL)
    
    choices <- sort(unique(as.character(dt$gene)))
    choices <- choices[!is.na(choices) & nzchar(choices)]
    
    if (!length(choices)) return(helpText("No genes available in prioritized gene table."))
    
    selectInput("func_gene_id", "Gene", choices = choices)
  })
  
  # ------------------------------------------------------------
  # Genes for enrichment
  # ------------------------------------------------------------
  genes_all_symbols <- reactive({
    dt <- priority_genes_all()
    g <- unique(as.character(dt$gene))
    g <- g[!is.na(g) & nzchar(g)]
    g
  })
  
  genes_scope_symbols <- reactive({
    dt <- priority_genes_all()
    validate(
      need(is.data.frame(dt) && nrow(dt) > 0, "No prioritized genes loaded.")
    )
    
    if (identical(input$func_scope, "cluster")) {
      req(input$func_cluster_id)
      
      map_df <- gene_cluster_long() %>%
        dplyr::filter(.data$cluster_id == input$func_cluster_id)
      
      validate(
        need(nrow(map_df) > 0, "No prioritized genes found for the selected cluster.")
      )
      
      g <- unique(as.character(map_df$gene))
    } else {
      g <- unique(as.character(dt$gene))
    }
    
    g <- g[!is.na(g) & nzchar(g)]
    validate(need(length(g) > 0, "No valid genes found for the selected scope."))
    g
  })
  
  # ------------------------------------------------------------
  # SYMBOL -> ENTREZ mapping cache
  # ------------------------------------------------------------
  gene_map_cache_symbols <- reactiveVal(NULL)
  gene_map_cache_map     <- reactiveVal(NULL)
  
  gene_map_all <- reactive({
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
    syms <- genes_all_symbols()
    validate(need(length(syms) > 0, "No valid gene symbols found in prioritized genes."))
    
    prev_syms <- gene_map_cache_symbols()
    prev_map  <- gene_map_cache_map()
    
    if (!is.null(prev_syms) && !is.null(prev_map) && identical(sort(prev_syms), sort(syms))) {
      return(prev_map)
    }
    
    m <- tryCatch({
      suppressMessages(
        clusterProfiler::bitr(
          syms,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = org.Hs.eg.db
        )
      )
    }, error = function(e) NULL)
    
    validate(
      need(is.data.frame(m) && nrow(m) > 0, "SYMBOL → ENTREZ mapping failed (empty).")
    )
    
    gene_map_cache_symbols(sort(syms))
    gene_map_cache_map(m)
    m
  })
  
  entrez_scope <- reactive({
    m <- gene_map_all()
    syms <- genes_scope_symbols()
    
    validate(need(length(syms) > 0, "No valid genes for the selected scope."))
    
    ids <- unique(m$ENTREZID[m$SYMBOL %in% syms])
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    validate(
      need(length(ids) > 0, "No ENTREZIDs for the selected scope (mapping empty).")
    )
    
    ids
  })
  
  universe_entrez_dataset <- reactive({
    m <- gene_map_all()
    u <- unique(m$ENTREZID)
    u <- u[!is.na(u) & nzchar(u)]
    
    validate(
      need(length(u) > 0, "Universe (dataset) is empty after mapping.")
    )
    
    u
  })
  
  universe_entrez_for_scope <- reactive({
    scope <- input$func_scope %||% "global"
    bg    <- input$enrich_background %||% "dataset"
    
    if (identical(scope, "global")) return(NULL)
    
    if (identical(bg, "dataset")) return(universe_entrez_dataset())
    
    NULL
  })
  
  # ------------------------------------------------------------
  # Background note
  # ------------------------------------------------------------
  output$enrich_bg_note <- renderUI({
    bg  <- input$enrich_background %||% "dataset"
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    sc <- input$func_scope %||% "global"
    sc_txt <- if (identical(sc, "cluster")) {
      paste0("Foreground: prioritized genes in <b>", htmltools::htmlEscape(input$func_cluster_id %||% ""), "</b>.")
    } else {
      "Foreground: prioritized genes in <b>global</b> (all prioritized genes)."
    }
    
    safe_n <- function(expr) {
      tryCatch({
        v <- eval(expr)
        if (is.null(v)) return(NA_integer_)
        if (is.atomic(v) && length(v) == 1) return(as.integer(v))
        if (is.vector(v)) return(as.integer(length(v)))
        if (is.data.frame(v)) return(as.integer(nrow(v)))
        NA_integer_
      }, error = function(e) NA_integer_)
    }
    
    n_fg <- safe_n(quote(length(entrez_scope())))
    n_uni_dataset <- safe_n(quote(length(universe_entrez_dataset())))
    n_uni_orgdb <- safe_n(quote(length(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID"))))
    
    uni_txt <- function(default_name) {
      if (identical(bg, "dataset")) {
        if (is.finite(n_uni_dataset)) {
          paste0("<b>Background (universe):</b> dataset-mapped prioritized genes (N=", n_uni_dataset, ").")
        } else {
          "<b>Background (universe):</b> dataset-mapped prioritized genes."
        }
      } else {
        if (is.finite(n_uni_orgdb)) {
          paste0("<b>Background (universe):</b> ", default_name, " (N≈", n_uni_orgdb, ").")
        } else {
          paste0("<b>Background (universe):</b> ", default_name, ".")
        }
      }
    }
    
    fg_cnt_txt <- if (is.finite(n_fg)) paste0("<br><b>Foreground size:</b> n=", n_fg, " genes.") else ""
    gs_txt <- "<b>GSSize:</b> number of genes per term in the selected universe (filtered by minGSSize / maxGSSize)."
    
    msg <- if (identical(tab, "tab_enrich_go")) {
      paste0(
        "<b>GO enrichment</b><br>",
        sc_txt,
        fg_cnt_txt, "<br>",
        uni_txt("OrgDb annotated genes (org.Hs.eg.db)"), "<br>",
        "<b>Test unit:</b> genes (ENTREZID).<br>",
        gs_txt, "<br>",
        "<b>Method:</b> clusterProfiler::enrichGO (pAdjustMethod=BH; filtered by FDR cutoff)."
      )
    } else if (identical(tab, "tab_enrich_kegg")) {
      paste0(
        "<b>KEGG enrichment</b><br>",
        sc_txt,
        fg_cnt_txt, "<br>",
        uni_txt("KEGG default background (organism-specific, hsa)"), "<br>",
        "<b>Test unit:</b> genes (ENTREZID).<br>",
        gs_txt, "<br>",
        "<b>Method:</b> clusterProfiler::enrichKEGG (organism=hsa; pAdjustMethod=BH; filtered by FDR cutoff)."
      )
    } else {
      paste0("<b>Background:</b> ", if (identical(bg, "dataset")) "Dataset genes" else "Default background")
    }
    
    htmltools::HTML(paste0(
      "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
      msg, "</div>"
    ))
  })
  
  # ------------------------------------------------------------
  # Triggers
  # ------------------------------------------------------------
  go_trigger   <- reactiveVal(NULL)
  kegg_trigger <- reactiveVal(NULL)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    if (identical(tab, "tab_enrich_go")) {
      go_trigger((go_trigger() %||% 0L) + 1L)
    } else if (identical(tab, "tab_enrich_kegg")) {
      kegg_trigger((kegg_trigger() %||% 0L) + 1L)
    }
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # GO enrichment
  # ------------------------------------------------------------
  go_enrich_raw <- eventReactive(go_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
    withProgress(message = "Running GO enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (OrgDb)" else paste("Universe:", length(uni)))
      
      ontos <- input$go_ontos %||% c("BP", "CC", "MF")
      ontos <- intersect(c("BP", "CC", "MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one GO ontology (BP/CC/MF)."))
      
      minGS <- input$enrich_min_gs %||% 10
      maxGS <- input$enrich_max_gs %||% 500
      simplify_flag <- isTRUE(input$go_simplify)
      
      res <- purrr::map_dfr(ontos, function(ont) {
        incProgress(0.6 / length(ontos), detail = paste("Ontology:", ont))
        
        eg <- suppressMessages(
          clusterProfiler::enrichGO(
            gene          = gene,
            universe      = uni,
            OrgDb         = org.Hs.eg.db,
            keyType       = "ENTREZID",
            ont           = ont,
            pAdjustMethod = "BH",
            pvalueCutoff  = 1,
            qvalueCutoff  = 1,
            minGSSize     = minGS,
            maxGSSize     = maxGS,
            readable      = TRUE
          )
        )
        
        if (is.null(eg) || is.null(eg@result) || !nrow(eg@result)) return(NULL)
        
        if (simplify_flag && requireNamespace("clusterProfiler", quietly = TRUE)) {
          eg <- tryCatch(
            suppressMessages(clusterProfiler::simplify(eg, cutoff = 0.7, by = "p.adjust", select_fun = min)),
            error = function(e) eg
          )
        }
        
        df <- as.data.frame(eg@result, stringsAsFactors = FALSE)
        if (!nrow(df)) return(NULL)
        df$Ontology <- ont
        df
      })
      
      if (is.null(res) || !nrow(res)) return(tibble::tibble())
      res
    })
  }, ignoreInit = TRUE)
  
  go_enrich_tbl <- reactive({
    df <- go_enrich_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$go_topn %||% 10
    
    df2 <- df %>%
      dplyr::mutate(
        Ontology = as.character(Ontology),
        pvalue   = suppressWarnings(as.numeric(pvalue)),
        p.adjust = suppressWarnings(as.numeric(p.adjust))
      ) %>%
      dplyr::filter(is.finite(p.adjust), is.finite(pvalue)) %>%
      dplyr::filter(p.adjust <= pcut) %>%
      dplyr::group_by(Ontology) %>%
      dplyr::arrange(p.adjust, .by_group = TRUE) %>%
      dplyr::slice_head(n = topn) %>%
      dplyr::ungroup()
    
    if (!nrow(df2)) return(tibble::tibble())
    df2
  })
  
  output$go_table <- DT::renderDT({
    df <- go_enrich_tbl()
    
    if (!nrow(df)) {
      return(DT::datatable(
        data.frame(Message = "No GO terms passed the cutoff."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    out <- df %>%
      dplyr::select(Ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      dplyr::mutate(
        pvalue   = sprintf("%.3g", as.numeric(pvalue)),
        p.adjust = sprintf("%.3g", as.numeric(p.adjust))
      )
    
    DT::datatable(
      out,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  scope_label <- function() {
    sc <- input$func_scope %||% "global"
    if (identical(sc, "cluster")) {
      cid <- input$func_cluster_id %||% ""
      if (nzchar(cid)) paste0("cluster ", cid) else "cluster"
    } else {
      "global"
    }
  }
  
  go_top_df <- reactive({
    df <- go_enrich_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$go_topn %||% 10
    
    df0 <- df %>%
      dplyr::mutate(
        Ontology   = as.character(Ontology),
        pvalue     = suppressWarnings(as.numeric(pvalue)),
        p.adjust   = suppressWarnings(as.numeric(p.adjust)),
        term_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "."), Description),
        score      = -log10(pvalue)
      ) %>%
      dplyr::filter(
        Ontology %in% c("BP", "CC", "MF"),
        is.finite(pvalue), is.finite(p.adjust),
        nzchar(term_short)
      )
    
    if (!nrow(df0)) return(tibble::tibble())
    
    sig <- df0 %>%
      dplyr::filter(p.adjust <= pcut) %>%
      dplyr::group_by(Ontology) %>%
      dplyr::arrange(p.adjust, .by_group = TRUE) %>%
      dplyr::slice_head(n = topn) %>%
      dplyr::ungroup()
    
    if (nrow(sig) > 0) {
      attr(sig, "fallback") <- FALSE
      return(sig)
    }
    
    fb <- df0 %>%
      dplyr::group_by(Ontology) %>%
      dplyr::arrange(p.adjust, .by_group = TRUE) %>%
      dplyr::slice_head(n = topn) %>%
      dplyr::ungroup()
    
    attr(fb, "fallback") <- TRUE
    fb
  })
  
  output$go_bar <- plotly::renderPlotly({
    dat <- go_top_df()
    
    cat("\n[DBG] head(go_top_df()):\n")
    print(utils::head(dat, 10))
    cat("\n[DBG] names(go_top_df()):\n")
    print(names(dat))
    
    req(is.data.frame(dat), nrow(dat) > 0)
    
    fb   <- isTRUE(attr(dat, "fallback"))
    pcut <- input$enrich_pcut %||% 0.05
    
    ttl <- paste0("GO enrichment (", scope_label(), ")")
    subttl <- if (fb) {
      paste0("No terms at FDR ≤ ", pcut, " — showing top terms ranked by FDR; bar height = -log10(pvalue)")
    } else {
      paste0("FDR ≤ ", pcut, " (bar height = -log10(pvalue))")
    }
    
    dat <- dat %>%
      dplyr::mutate(
        Ontology = factor(as.character(Ontology), levels = c("BP", "CC", "MF")),
        p_adj    = suppressWarnings(as.numeric(p.adjust)),
        score    = suppressWarnings(as.numeric(score)),
        tooltip  = paste0(
          "<b>", htmltools::htmlEscape(Description), "</b>",
          "<br>Ontology: ", Ontology,
          "<br>-log10(pvalue): ", sprintf("%.3f", score),
          "<br>FDR: ", ifelse(is.finite(p_adj), formatC(p_adj, format = "e", digits = 2), "NA")
        )
      ) %>%
      dplyr::filter(
        !is.na(Ontology),
        is.finite(score),
        !is.na(term_short),
        nzchar(term_short)
      )
    
    req(nrow(dat) > 0)
    
    cols_go <- c(BP = "darkgreen", CC = "orange", MF = "darkblue")
    
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x    = term_short,
        y    = score,
        fill = Ontology,
        text = tooltip
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~Ontology, scales = "free_x", nrow = 1, drop = FALSE) +
      ggplot2::scale_fill_manual(values = cols_go, guide = "none", drop = FALSE) +
      ggplot2::labs(
        title = ttl,
        subtitle = subttl,
        x = NULL,
        y = "-log10(pvalue)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1)
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        margin = list(l = 70, r = 20, t = 80, b = 180),
        font   = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GO_barplot",
          width = 1800,
          height = 900,
          scale = 2
        )
      )
  })
  
  # ------------------------------------------------------------
  # KEGG enrichment
  # ------------------------------------------------------------
  kegg_enrich_raw <- eventReactive(kegg_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    
    withProgress(message = "Running KEGG enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (KEGG)" else paste("Universe:", length(uni)))
      
      minGS <- input$enrich_min_gs %||% 10
      maxGS <- input$enrich_max_gs %||% 500
      
      ek <- suppressMessages(
        clusterProfiler::enrichKEGG(
          gene          = gene,
          universe      = uni,
          organism      = "hsa",
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,
          qvalueCutoff  = 1,
          minGSSize     = minGS,
          maxGSSize     = maxGS
        )
      )
      
      if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(tibble::tibble())
      as.data.frame(ek@result, stringsAsFactors = FALSE)
    })
  }, ignoreInit = TRUE)
  
  kegg_enrich_tbl <- reactive({
    df <- kegg_enrich_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$enrich_kegg_top %||% 15
    
    df2 <- df %>%
      dplyr::mutate(
        pvalue   = suppressWarnings(as.numeric(pvalue)),
        p.adjust = suppressWarnings(as.numeric(p.adjust))
      ) %>%
      dplyr::filter(is.finite(p.adjust), is.finite(pvalue)) %>%
      dplyr::filter(p.adjust <= pcut) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = topn)
    
    if (!nrow(df2)) return(tibble::tibble())
    df2
  })
  
  output$kegg_table <- DT::renderDT({
    df <- kegg_enrich_tbl()
    
    if (!nrow(df)) {
      return(DT::datatable(
        data.frame(Message = "No KEGG pathways passed the cutoff."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    out <- df %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      dplyr::mutate(
        pvalue   = sprintf("%.3g", as.numeric(pvalue)),
        p.adjust = sprintf("%.3g", as.numeric(p.adjust))
      )
    
    DT::datatable(
      out,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  kegg_top_df <- reactive({
    df <- kegg_enrich_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$enrich_kegg_top %||% 15
    
    df0 <- df %>%
      dplyr::mutate(
        pvalue     = suppressWarnings(as.numeric(pvalue)),
        p.adjust   = suppressWarnings(as.numeric(p.adjust)),
        term_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "."), Description),
        score      = -log10(pvalue)
      ) %>%
      dplyr::filter(is.finite(pvalue), is.finite(p.adjust), nzchar(term_short))
    
    if (!nrow(df0)) return(tibble::tibble())
    
    sig <- df0 %>%
      dplyr::filter(p.adjust <= pcut) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = topn)
    
    if (nrow(sig) > 0) {
      attr(sig, "fallback") <- FALSE
      return(sig)
    }
    
    fb <- df0 %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = topn)
    
    attr(fb, "fallback") <- TRUE
    fb
  })
  
  output$kegg_bar <- plotly::renderPlotly({
    df <- kegg_top_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    fb   <- isTRUE(attr(df, "fallback"))
    pcut <- input$enrich_pcut %||% 0.05
    
    ttl <- paste0(
      "KEGG enrichment (", scope_label(), ")",
      if (fb) paste0(" — no pathways at FDR ≤ ", pcut, " (showing top by FDR)") else paste0(" — FDR ≤ ", pcut)
    )
    
    df <- df %>%
      dplyr::mutate(
        score = -log10(pvalue),
        padj = suppressWarnings(as.numeric(p.adjust)),
        tooltip = paste0(
          "<b>", htmltools::htmlEscape(Description), "</b>",
          "<br>-log10(pvalue): ", sprintf("%.3f", score),
          "<br>FDR: ", ifelse(is.finite(padj), formatC(padj, format = "e", digits = 2), "NA")
        )
      )
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = stats::reorder(term_short, score),
        y = score,
        fill = padj,
        text = tooltip
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(title = ttl, x = NULL, y = "-log10(pvalue)", fill = "FDR") +
      ggplot2::scale_fill_gradient(
        low = "red",
        high = "orange",
        na.value = "yellow",
        trans = "reverse"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        margin = list(l = 240, r = 20, t = 60, b = 50),
        font = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "KEGG_barplot",
          width = 1600,
          height = 900,
          scale = 2
        )
      )
  })
  
  #############################################################
  
}

shinyApp(ui, server)
