## app.R — GWAS Catalog Inspector (ADAPTED TO NonSyn: new/old + canonical clusters + Manhattan_combo + GO/KEGG + LD)
#/Volumes/DISK1TB/Inspector_app_slaves_ngroc/GItools/app/Catalog_inspector/app.R
options(shiny.maxRequestSize = 1024*1024^2)

library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(purrr)
library(tibble)
library(clusterProfiler)
library(igraph)
library(htmltools)
library(patchwork)
library(shinycssloaders)
library(shinyFiles)
library(shinyjs)

# Gene model (hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(IRanges)

# (opcional però recomanat per labels)
library(AnnotationDbi)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# ===========================
# Runtime logging (EARLY)
# ===========================
APP_KEY <- "catalog"  # <- canvia: catalog / nonsyn / ewastum / ewasdis

APP_LOG_DIR <- Sys.getenv("GITOOLS_LOG_DIR", "")
if (!nzchar(APP_LOG_DIR)) {
  APP_LOG_DIR <- file.path(dirname(getwd()), "_logs")  # -> GItools/app/_logs
}
dir.create(APP_LOG_DIR, showWarnings = FALSE, recursive = TRUE)

APP_LOG <- file.path(APP_LOG_DIR, paste0(APP_KEY, "_runtime.log"))

applog <- function(...) {
  cat(paste0(..., collapse = ""), "\n", file = APP_LOG, append = TRUE)
}

applog("[", APP_KEY, "] boot @ ", format(Sys.time()), " wd=", getwd(),
       " pid=", Sys.getpid(),
       " host=", getOption("shiny.host"),
       " port=", getOption("shiny.port"))

# Captura errors de Shiny al log
options(shiny.error = function() {
  applog("[", APP_KEY, "][ERROR] ", format(Sys.time()))
  applog(paste(capture.output(sys.calls()), collapse = "\n"))
})


# --- Fixar dplyr com a select/mutate/arrange oficial ---
select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
summarise  <- dplyr::summarise
rename     <- dplyr::rename

# --- GI shared state (portable via config.R) ---
cfg_file <- file.path("..", "..", "config.R")  # des de app/Catalog_inspector → GItools/config.R
stopifnot(file.exists(cfg_file))
source(cfg_file, local = TRUE)

cfg <- gi_cfg()

ld_file <- file.path(cfg$shared, "mod_ld.R")

cat("\n[LD-SOURCE] cfg$shared / gi_shared_root = ", cfg$shared %||% gi_shared_root, "\n", sep = "")
cat("[LD-SOURCE] ld_file = ", ld_file, "\n", sep = "")
cat("[LD-SOURCE] exists(ld_file) = ", file.exists(ld_file), "\n", sep = "")
cat("[LD-SOURCE] normalizePath(ld_file) = ", normalizePath(ld_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cl_file <- file.path(cfg$shared, "gi_clusters_canonical.R")

gi_state_file <- file.path(cfg$shared, "gi_state.R")
stopifnot(file.exists(gi_state_file))
source(gi_state_file, local = TRUE)

# Catalog as MASTER: write parameters for slaves
write_params_rds <- function(sid, input, min_hits_value = NULL, extra = list()) {
  stopifnot(nzchar(sid))
  p <- gi_state_paths(sid)$params
  
  par <- c(list(
    stamp          = as.integer(Sys.time()),
    cluster_method = input$cluster_method %||% NULL,
    hits_mode      = input$hits_mode %||% NULL,
    thr_type       = if (identical(input$cluster_method, "window")) "pthr" else "min_logp",
    thr_value      = if (identical(input$cluster_method, "window")) (input$pthr %||% NULL) else (input$min_logp %||% NULL),
    pthr           = input$pthr %||% NULL,
    min_logp       = input$min_logp %||% NULL,
    flank_bp       = input$flank %||% NULL,
    win_bp         = input$win_bp %||% NULL,
    step_bp        = input$step_bp %||% NULL,
    min_hits       = as.integer(min_hits_value %||% 3L)
  ), extra)
  
  saveRDS(par, p)
  invisible(TRUE)
}

stopifnot(file.exists(ld_file), file.exists(cl_file))
source(ld_file, local = TRUE)
source(cl_file, local = TRUE)

stopifnot(exists("ld_module_ui"), exists("ld_module_server"))

cat("[GI] shared =", cfg$shared, "\n")
cat("[GI] using gi_state.R at:", gi_state_file, "\n")
cat("[GI] gi_shared_root =", gi_shared_root, "\n")
# cat("[GI] gi_state_root  =", gi_state_root, "\n")
cat("[GI] gi_state_root  = ", gi_state_root(), "\n")



# =============================================================================
# Helpers generals
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

chr_map_plink19 <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT","M","MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

chr_label_plink <- function(chr_num) {
  out <- as.character(chr_num)
  out[out == "23"] <- "X"
  out[out == "24"] <- "Y"
  out[out == "26"] <- "MT"
  out
}

norm_chr_generic <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M","MT","m","mt","MTDNA","mtdna")] <- "MT"
  toupper(x)
}

parse_p_robust <- function(p){
  if (is.numeric(p)) return(as.numeric(p))
  p_chr <- trimws(as.character(p))
  p_chr <- gsub("\\s+", "", p_chr)
  needs_swap <- !grepl("\\.", p_chr) & grepl(",", p_chr)
  p_chr[needs_swap] <- gsub(",", ".", p_chr[needs_swap], fixed = TRUE)
  suppressWarnings(as.numeric(p_chr))
}

plotly_message <- function(msg) {
  plotly::plot_ly(
    x = 0, y = 0,
    type = "scatter", mode = "markers",
    marker = list(opacity = 0)
  ) %>%
    plotly::layout(
      title = list(text = paste0("<b>", msg, "</b>"), x = 0.5, xanchor = "center"),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

detect_rsid_col <- function(df) {
  cn <- names(df)
  patt <- "(^snp$|^snp\\.|^rsid$|^rsid\\.|marker|id)"
  hit <- grep(patt, cn, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) return(hit[1])
  return(NULL)
}

pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

pick_first <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[1] else NA_character_
}


make_mode_thr_tag <- function(cluster_method, pthr, min_logp) {
  mode_tag <- if (identical(cluster_method, "window")) "w" else "h"
  thr_val  <- if (identical(cluster_method, "window")) pthr else min_logp
  thr_txt  <- gsub("\\.", "p", sprintf("%.2f", thr_val %||% 0))
  list(mode_tag=mode_tag, thr_txt=thr_txt)
}

# =============================================================================
#  helpers LD module
# =============================================================================

# global/app.R (fora server, o al principi)
txdb <- NULL
if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
}

# =============================================================================
# helper Integrator
# =============================================================================

chr_map_integrator <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT", "M", "MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

norm_int <- function(x) {
  suppressWarnings(as.integer(readr::parse_number(as.character(x))))
}

sanitize_bridge <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  if ("chr" %in% names(df))   df$chr   <- chr_map_integrator(df$chr)
  if ("start" %in% names(df)) df$start <- norm_int(df$start)
  if ("end" %in% names(df))   df$end   <- norm_int(df$end)
  
  if ("start" %in% names(df) && "end" %in% names(df)) {
    s0 <- df$start
    e0 <- df$end
    df$start <- pmin(s0, e0, na.rm = FALSE)
    df$end   <- pmax(s0, e0, na.rm = FALSE)
  }
  
  df
}

sanitize_bridge_genes <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  if (!"gene" %in% names(df)) return(df)
  
  df %>%
    dplyr::mutate(
      gene = as.character(gene),
      gene = trimws(gene)
    ) %>%
    tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s+-\\s+") %>%
    dplyr::mutate(
      gene = trimws(gene)
    ) %>%
    dplyr::filter(
      !is.na(gene),
      nzchar(gene),
      tolower(gene) != "numeric",
      !grepl("^[0-9]+$", gene)
    ) %>%
    dplyr::distinct()
}

# =============================================================================
# Catalog helpers
# =============================================================================

catalog_prepare <- function(cat_df) {
  chr_col <- pick_col(cat_df, c("CHR","chr","chrom","CHROM","chromosome","CHR_ID"))
  pos_col <- pick_col(cat_df, c("POS","pos","BP","bp","position","POSITION","base_pair_location"))
  if (is.null(chr_col) || is.null(pos_col)) {
    stop("Catalog: no he pogut detectar CHR/POS.")
  }
  
  rs_col    <- pick_col(cat_df, c("rsid","RSID","SNPS","SNP","snp","variant_id","VARIANT_ID"))
  p_col     <- pick_col(cat_df, c("P_VALUE","p_value","P","p","PVAL","P-VALUE"))
  or_col    <- pick_col(cat_df, c("OR_BETA","OR","BETA","beta","ODDSRATIO","odds_ratio"))
  pm_col    <- pick_col(cat_df, c("PUBMEDID","PMID","PubMed","PUBMED_ID","PUBMED"))
  gene_col  <- pick_col(cat_df, c("MAPPED_GENE","mapped_gene","GENE","gene"))
  dis_col   <- pick_col(cat_df, c("DISEASE","Disease","disease"))
  trait_col <- pick_col(cat_df, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
  
  # ⚠️ clau: CHR com integer (PLINK) i POS com numeric, sense tocar OR/P massivament
  out <- cat_df %>%
    transmute(
      chr = chr_map_plink19(.data[[chr_col]]),
      pos = suppressWarnings(readr::parse_number(as.character(.data[[pos_col]]))),
      
      rsid_raw   = if (!is.null(rs_col))  .data[[rs_col]]  else NA,
      p_raw      = if (!is.null(p_col))   .data[[p_col]]   else NA,
      or_raw     = if (!is.null(or_col))  .data[[or_col]]  else NA,
      pmid_raw   = if (!is.null(pm_col))  .data[[pm_col]]  else NA,
      gene_raw   = if (!is.null(gene_col)) .data[[gene_col]] else NA,
      disease_raw = if (!is.null(dis_col))  .data[[dis_col]]  else NA,
      trait_raw   = if (!is.null(trait_col)) .data[[trait_col]] else NA
    ) %>%
    filter(!is.na(chr), !is.na(pos))
  
  out
}

extract_catalog_by_clusters <- function(clusters, cat_prepared) {
  if (is.null(clusters) || !nrow(clusters)) return(tibble::tibble())
  if (is.null(cat_prepared) || !nrow(cat_prepared)) return(tibble::tibble())
  
  out <- vector("list", nrow(clusters))
  
  for (i in seq_len(nrow(clusters))) {
    
    chr_i   <- suppressWarnings(as.integer(clusters$chr[i]))
    st      <- suppressWarnings(as.numeric(clusters$start[i]))
    en      <- suppressWarnings(as.numeric(clusters$end[i]))
    chr_lab <- chr_label_plink(chr_i)
    
    if (!is.finite(chr_i) || !is.finite(st) || !is.finite(en) || st > en) {
      out[[i]] <- NULL
      next
    }
    
    sub <- cat_prepared %>%
      dplyr::filter(chr == chr_i, pos >= st, pos <= en)
    
    if (!nrow(sub)) {
      out[[i]] <- NULL
      next
    }
    
    # --- cluster metadata per row (like NonSyn normalized output) ---
    sub$cluster_id    <- as.character(clusters$cluster_id[i])
    sub$cluster       <- if ("cluster_n" %in% names(clusters)) as.integer(clusters$cluster_n[i]) else NA_integer_
    sub$cluster_start <- as.integer(st)
    sub$cluster_end   <- as.integer(en)
    
    # --- standard display fields ---
    sub$CHR <- chr_lab
    sub$POS <- sub$pos
    
    # parse only on the subset
    sub2 <- sub %>%
      dplyr::mutate(
        rsid    = ifelse(is.na(rsid_raw), NA_character_, as.character(rsid_raw)),
        P_VALUE = ifelse(is.na(p_raw), NA_real_, parse_p_robust(p_raw)),
        OR_BETA = {
          x <- or_raw
          if (is.numeric(x)) as.numeric(x) else suppressWarnings(readr::parse_number(as.character(x)))
        },
        PUBMEDID     = ifelse(is.na(pmid_raw), NA_character_, as.character(pmid_raw)),
        MAPPED_GENE  = ifelse(is.na(gene_raw), NA_character_, as.character(gene_raw)),
        DISEASE      = ifelse(is.na(disease_raw), NA_character_, as.character(disease_raw)),
        MAPPED_TRAIT = ifelse(is.na(trait_raw), NA_character_, as.character(trait_raw))
      ) %>%
      dplyr::select(
        cluster_id, cluster, cluster_start, cluster_end,
        CHR, POS, rsid, P_VALUE, OR_BETA, PUBMEDID, MAPPED_GENE, DISEASE, MAPPED_TRAIT
      )
    
    out[[i]] <- sub2
  }
  
  hits <- dplyr::bind_rows(out)
  if (!nrow(hits)) return(tibble::tibble())
  
  hits %>%
    dplyr::mutate(hit_key = ifelse(!is.na(rsid) & nzchar(rsid), rsid, paste0(CHR, ":", POS))) %>%
    dplyr::distinct(cluster_id, hit_key, .keep_all = TRUE)
}



add_n_catalog <- function(clusters, hits) {
  if (is.null(clusters) || !nrow(clusters)) return(clusters)
  
  if (!"n_catalog" %in% names(clusters)) {
    clusters$n_catalog <- 0L
  } else {
    clusters$n_catalog <- dplyr::coalesce(as.integer(clusters$n_catalog), 0L)
  }
  
  if (is.null(hits) || !nrow(hits) || !"cluster_id" %in% names(hits)) {
    clusters$n_catalog <- 0L
    return(clusters)
  }
  
  counts <- hits %>% dplyr::count(cluster_id, name = "n_catalog_new")
  clusters %>%
    dplyr::left_join(counts, by = "cluster_id") %>%
    dplyr::mutate(n_catalog = dplyr::coalesce(as.integer(n_catalog_new), 0L)) %>%
    dplyr::select(-n_catalog_new)
}

infer_clusters_from_catalog_hits <- function(hits_df) {
  if (is.null(hits_df) || !is.data.frame(hits_df) || !nrow(hits_df)) return(tibble())
  
  validate(
    need("cluster_id" %in% names(hits_df), "Precomputed Catalog file must contain column: cluster_id."),
    need(any(c("POS","POS_std","position","BP","pos") %in% names(hits_df)), "Precomputed Catalog file must contain a position column (POS/POS_std/BP...)."),
    need(any(c("CHR","CHR_std","chr","chrom") %in% names(hits_df)), "Precomputed Catalog file must contain a chromosome column (CHR/CHR_std...).")
  )
  
  chr_col <- pick_col(hits_df, c("CHR_std","CHR","chr","chrom"))
  pos_col <- pick_col(hits_df, c("POS_std","POS","BP","pos","position","POSITION"))
  
  x <- hits_df %>%
    transmute(
      cluster_id = as.character(cluster_id),
      CHR_raw = .data[[chr_col]],
      POS_raw = .data[[pos_col]]
    ) %>%
    mutate(
      chrN = norm_chr_generic(CHR_raw),
      POSn = suppressWarnings(readr::parse_number(as.character(POS_raw)))
    ) %>%
    filter(!is.na(chrN), !is.na(POSn), nzchar(cluster_id))
  
  cl <- x %>%
    group_by(cluster_id, chrN) %>%
    summarise(start=min(POSn), end=max(POSn), .groups="drop") %>%
    mutate(chr = chr_map_plink19(chrN)) %>%
    filter(is.finite(chr)) %>%
    arrange(chr, start, end) %>%
    group_by(chr) %>%
    mutate(cluster_n = row_number()) %>%
    ungroup() %>%
    mutate(
      center = round((start+end)/2),
      n_snps = NA_integer_,
      top_snp = NA_character_,
      top_logp = NA_real_,
      cluster_size_kb = round((end-start)/1000, 2),
      n_catalog = 0L
    )
  
  cl
}


# =============================================================================
# Reference hg38 chr lengths (for Manhattan cum positions)
# =============================================================================

.ref_hg38 <- tibble::tibble(
  chr = c(as.character(1:22), "X", "Y", "MT"),
  len = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468,
    156040895, 57227415, 16569
  )
) %>%
  mutate(chr_cum = cumsum(len) - len)


# =============================================================================
# GO/KEGG helpers (NonSyn-style)
# =============================================================================

split_multi_unique <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(character(0))
  # separadors típics: ; , | espais múltiples
  y <- unlist(strsplit(x, "[;,|]"))
  y <- trimws(y)
  y <- y[nzchar(y)]
  unique(y)
}

# ============================================================
# GWAS significance bridge
# ============================================================

pick_col_first <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

normalize_chr_bridge <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT", "M", "MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

normalize_cluster_id_bridge <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^cluster_", "", x, ignore.case = TRUE)
  
  x <- ifelse(
    grepl("^[0-9XYM]+_[0-9]+$", x, ignore.case = TRUE),
    paste0("chr", x),
    x
  )
  
  x <- sub("^CHR", "chr", x, ignore.case = TRUE)
  x <- sub("^chr23_", "chrX_", x, ignore.case = TRUE)
  x <- sub("^chr24_", "chrY_", x, ignore.case = TRUE)
  x <- sub("^chr26_", "chrMT_", x, ignore.case = TRUE)
  
  x
}

build_gwas_significance_bridge <- function(gwas_df, clusters_df) {
  stopifnot(is.data.frame(gwas_df), is.data.frame(clusters_df))
  
  pick_col_first <- function(df, choices) {
    nm <- intersect(choices, names(df))
    if (length(nm)) nm[1] else NULL
  }
  
  chr_col <- pick_col_first(gwas_df, c("chr", "CHR", "chrom", "CHROM", "chromosome"))
  pos_col <- pick_col_first(gwas_df, c("position", "POS", "pos", "BP", "bp"))
  rs_col  <- pick_col_first(gwas_df, c("rsid", "RSID", "SNP", "snp", "marker", "MARKER", "id_hit", "ID"))
  p_col   <- pick_col_first(gwas_df, c("P", "p", "pval", "Pval", "PVAL", "p_value", "PVALUE"))
  
  if (is.null(chr_col) || is.null(pos_col) || is.null(p_col)) {
    stop("GWAS input must contain chr + position + P columns.")
  }
  
  gwas2 <- gwas_df %>%
    dplyr::transmute(
      chr = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[chr_col]])))),
      position = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
      rsid = if (!is.null(rs_col)) as.character(.data[[rs_col]]) else NA_character_,
      p_value = suppressWarnings(as.numeric(.data[[p_col]]))
    ) %>%
    dplyr::filter(
      is.finite(chr),
      is.finite(position),
      is.finite(p_value),
      p_value > 0
    ) %>%
    dplyr::mutate(
      logp = -log10(p_value)
    )
  
  cl2 <- clusters_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      start = suppressWarnings(as.integer(start)),
      end = suppressWarnings(as.integer(end))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(start), is.finite(end)
    )
  
  gwas2 %>%
    dplyr::left_join(
      cl2 %>% dplyr::rename(cluster_chr = chr, cluster_start = start, cluster_end = end),
      by = dplyr::join_by(chr == cluster_chr, position >= cluster_start, position <= cluster_end)
    ) %>%
    dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
    dplyr::select(cluster_id, chr, position, rsid, p_value, logp) %>%
    dplyr::distinct() %>%
    dplyr::arrange(cluster_id, position, dplyr::desc(logp))
}


# =============================================================================
# UI
# =============================================================================
# ===========================
# UI (Catalog) — UPDATED:
#   - Top-level tabs: Analysis / Enrichment / LD (same level)
#   - Enrichment includes:
#       (A) Gene/Disease/Trait enrichment
#       (B) Functional enrichment: GO / KEGG / GoSlim (GoSlim inside same tabset)
# ===========================

ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("📚 GWAS Catalog Inspector")
  ),
  id = "topnav",
  
  # -----------------------------
  # Global light-grey panels for plots/tables
  # -----------------------------
  tags$head(
    tags$style(HTML("
      .panel-lite{
        background:#f2f2f2;
        border:1px solid #e0e0e0;
        border-radius:10px;
        padding:12px;
        margin-bottom:12px;
      }
    "))
  ),
  
  #-----------------
  tags$head(tags$script(HTML("
  Shiny.addCustomMessageHandler('gwas_loading', function(x){
    var el = document.getElementById('gwas_loading_overlay');
    if(el) el.style.display = x ? 'block' : 'none';
  });
"))),
  #-----------------
  
  # ========================================================================
  # TAB 1 — ANALYSIS
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        # -----------------------------
        # Step 1: GWAS input
        # -----------------------------
        h3(
          "Step 1 · Load GWAS hits",
          style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
        ),
        fluidRow(
          column(4, actionButton("info_format", "ℹ️ format")),
          column(4, actionButton("perf_info", "ℹ️ Load")),
          column(4, actionButton(
            "reset_case",
            "Reset",
            icon = icon("rotate-left"),
            class = "btn-warning"
          ))
        ),
        
        fileInput("gwas_file", "GWAS p-value table (TSV/CSV)", accept = c(".tsv", ".txt", ".csv")),
        checkboxInput("gwas_header", "First row contains column names", TRUE),
        radioButtons(
          "gwas_sep",
          "Separator",
          c("Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"),
          selected = "\t"
        ),
        uiOutput("gwas_p_selector"),
        tags$hr(),
        
        # -----------------------------
        # Step 2: Clustering
        # -----------------------------
        h3("Step 2 · Clustering GWAS hits", style = "color:#1A4E8A; font-size:22px; font-weight:700;"),
        actionButton("info_01", "ℹ️ Clustering method", class = "btn btn-default"),
        
        radioButtons(
          "cluster_method", "Clustering method:",
          choices = c(
            "By hit intervals (thr + flank → merge)" = "window",
            "By hit density (min_logp + min_hits)"   = "hits"
          ),
          selected = "window"
        ),
        
        # Method: window
        conditionalPanel(
          condition = "input.cluster_method == 'window'",
          sliderInput(
            "pthr",
            "-log10(P) threshold",
            min = 2, max = 20, value = 8, step = 0.5
          ),
          numericInput(
            "flank",
            "Flank (+/- bp)",
            value = 10000, min = 0, max = 10000000, step = 1000
          ),
          numericInput(
            "min_hits_window",
            "Minimum GWAS hits per cluster",
            value = 3, min = 1, max = 1000, step = 1
          )
        ),
        
        # Method: hits (3 submodes)
        conditionalPanel(
          condition = "input.cluster_method == 'hits'",
          
          radioButtons(
            "hits_mode", "Hits mode:",
            choices = c(
              "1Mb hit-span (consecutive hits within 1Mb)" = "span1mb",
              "Tiled windows (non-overlapping)"            = "tiled",
              "Sliding windows (overlapping)"              = "sliding"
            ),
            selected = "span1mb"
          ),
          
          sliderInput(
            "min_logp",
            "-log10(P) threshold (hit significance)",
            min = 2, max = 20, value = 8, step = 0.5
          ),
          
          numericInput(
            "min_hits_hits",
            "Minimum GWAS hits per cluster/window",
            value = 3, min = 1, max = 1000, step = 1
          ),
          
          numericInput(
            "win_bp",
            "Window size (bp)",
            value = 1e6, min = 1e4, max = 5e7, step = 1e4
          ),
          
          conditionalPanel(
            condition = "input.hits_mode == 'sliding'",
            numericInput(
              "step_bp",
              "Step (bp)",
              value = 1e5, min = 1e3, max = 5e7, step = 1e3
            )
          ),
          
          conditionalPanel(
            condition = "input.hits_mode != 'sliding'",
            helpText("For 'tiled' and '1Mb hit-span', step is implicit.")
          )
        ),
        
        actionButton(
          "build_ranges",
          "➊ Generate intervals → merge → clusters",
          style = "background-color: #ffdd57; color: black; font-weight: bold;"
        ),
        h4("Preview selected_intervals.range (clusters)"),
        div(class = "panel-lite", verbatimTextOutput("ranges_preview")),
        tags$hr(),
        
        # -----------------------------
        # Step 3: Catalog extraction
        # -----------------------------
        h3(
          "Step 3 · Extract Catalog hits per cluster",
          style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
        ),
        textOutput("catalog_status"),
        actionButton(
          "run_catalog",
          "➋ Extract Catalog hits from ALL clusters",
          style = "background-color: #ffdd57; color: black; font-weight: bold;"
        ),
        textOutput("catalog_summary"),
        br(),
        div(class = "panel-lite", verbatimTextOutput("catalog_log")),
        br(),
        downloadButton("download_catalog", "⬇️ Catalog hits (CSV)"),
        downloadButton("download_catalog_rds", "⬇️ Catalog hits (RDS)"),
        tags$hr(),
        downloadButton("dl_candidates_zip", "⬇️ Download catalog candidates (ZIP)")
      ),
      
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            h4("Top: GWAS p-values · Bottom: GWAS Catalog hits inside clusters"),
            
            #--------
            tags$div(
              id = "gwas_loading_overlay",
              style = "display:none; position:fixed; inset:0; background:rgba(255,255,255,0.75); z-index:9999;",
              tags$div(
                style = "position:absolute; top:45%; left:50%; transform:translate(-50%,-50%); background:#fff; padding:18px 22px; border-radius:14px; box-shadow:0 8px 30px rgba(0,0,0,0.12); width:520px; max-width:92vw;",
                tags$div(style = "font-size:16px; font-weight:700; color:#1A4E8A; text-align:center;", "⏳ Loading GWAS…"),
                tags$div(
                  style = "margin-top:10px; font-size:14px; color:#444; line-height:1.35;",
                  HTML(paste0(
                    "<b>Performance note:</b> GItools is tuned for large input tables. ",
                    "Catalog Inspector apply an initial <b>P &lt; 0.05</b> filter and <b>smart downsampling</b> ",
                    "(e.g., target 50k variants while always retaining peaks at P ≤ 1e-6 and preserving genome-wide coverage via 10 Mb bins with a fixed seed) ",
                    "to keep plots and clustering fast and reproducible.<br><br>",
                    
                    "<b>What we are doing:</b><br>",
                    "1) <b>Initial filter</b>: keep only variants with <b>P &lt; 0.05</b> (pcut = 0.05).<br>",
                    "2) <b>Smart downsampling</b> to speed up plots and clustering:<br>",
                    "&nbsp;&nbsp;• target size: <b>50,000</b> SNPs (target_n = 50,000)<br>",
                    "&nbsp;&nbsp;• always keep strong signals: <b>P ≤ 1e-6</b> (peaks_keep_p = 1e-6)<br>",
                    "&nbsp;&nbsp;• preserve genome-wide coverage using <b>10 Mb</b> bins (bin_bp = 1e7)<br>",
                    "&nbsp;&nbsp;• reproducible sampling (seed = 123)<br><br>",
                    "<i>Large files may take a moment — thanks for your patience.</i>"
                  ))
                )
              )
            ),
            #--------
            
            div(class = "panel-lite", withSpinner(plotlyOutput("manhattan_combo", height = "700px"))),
            tags$hr(),
            helpText("Select a window on 'Catalog hits' to obtain a link to UCSC browser (see below)"),
            
            fluidRow(
              column(6, div(class = "panel-lite", uiOutput("debug_ucsc_state"))),
              column(2, h4("UCSC links to expanded region:")),
              column(
                4,
                div(
                  style = "margin-top:8px; background-color:#f2f2f2; padding:10px; border-radius:8px; width:100%;",
                  uiOutput("ucsc_link_gwas"),
                  uiOutput("ucsc_link_catalog")
                )
              )
            ),
            
            tags$hr(),
            h4("Clusters summary"),
            div(class = "panel-lite", withSpinner(DTOutput("cluster_dt")))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📚 Catalog hits</span>"),
            h4("GWAS Catalog hits assigned to clusters"),
            
            fluidRow(
              column(
                width = 6,
                radioButtons(
                  "catalog_scope", "Scope",
                  choices = c(
                    "Pick a cluster here"               = "pick",
                    "Selected cluster (Clusters table)" = "selected"
                  ),
                  selected = "pick",
                  inline = TRUE
                ),
                conditionalPanel(
                  condition = "input.catalog_scope == 'pick'",
                  uiOutput("catalog_cluster_pick_ui")
                )
              ),
              column(width = 6, uiOutput("catalog_gene_pick_ui"))
            ),
            
            tags$hr(),
            div(class = "panel-lite", verbatimTextOutput("catalog_summary")),
            div(class = "panel-lite", DTOutput("catalog_table"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔬 GWAS hits</span>"),
            h4("GWAS hits used for clustering"),
            div(class = "panel-lite", withSpinner(DTOutput("hits_tbl")))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Catalog visualizations</span>"),
            h4("Distribution of DISEASE"),
            div(class = "panel-lite", withSpinner(plotOutput("catalog_disease_bar", height = "380px"))),
            tags$hr(),
            
            h4("Top MAPPED_TRAIT terms"),
            div(class = "panel-lite", withSpinner(plotOutput("catalog_mapped_trait_bar", height = "380px"))),
            tags$hr(),
            
            h4("Gene ↔ Trait relationships & interpretation"),
            fluidRow(
              column(
                width = 7,
                div(class = "panel-lite", withSpinner(plotOutput("catalog_gene_trait_plot", height = "600px")))
              ),
              column(
                width = 5,
                h5("Interpretation summary"),
                div(class = "panel-lite", verbatimTextOutput("catalog_interpret"))
              )
            )
          )
        )
      )
    )
  ),
  
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 Table genes</span>"),
    
    radioButtons(
      "genes_scope", "Scope",
      choices = c("Global (clustered hits)" = "global", "Cluster" = "cluster"),
      selected = "global",
      inline = TRUE
    ),
    
    conditionalPanel(
      condition = "input.genes_scope == 'cluster'",
      uiOutput("genes_cluster_ui")
    ),
    
    tags$hr(),
    div(class = "panel-lite", withSpinner(DT::DTOutput("genes_table")))
  ),
  
  # ========================================================================
  # TAB 2 — ENRICHMENT
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>✳️ Enrichment</span>"),
    
    tabsetPanel(
      id = "enrich_top_tabs",
      div(class = "panel-lite", uiOutput("enrich_bg_note")),
      
      # ------------------------------------
      # (A) Disease/Trait Enrichment
      # ------------------------------------
      tabPanel(
        title = HTML("<span style='font-size:15px; font-weight:600;'>📊 Disease/Trait</span>"),
        br(),
        h3("Enrichment of GWAS Catalog terms", style = "color:#1A4E8A; font-weight:700; margin-top:10px;"),
        helpText("ORA DISEASE / TRAIT terms among Catalog hits."),
        
        fluidRow(
          column(
            width = 4,
            div(class = "panel-lite",
                selectInput(
                  "dt_scope", "Scope",
                  choices = c(
                    "Global (all clustered hits)"       = "global",
                    "Cluster (filter hits by cluster)"  = "pick",
                    "Gene (filter hits by gene)"        = "gene"
                  ),
                  selected = "global"
                ),
                
                conditionalPanel(
                  condition = "input.dt_scope === 'pick'",
                  uiOutput("dt_cluster_ui")
                ),
                
                conditionalPanel(
                  condition = "input.dt_scope === 'gene'",
                  uiOutput("dt_gene_ui")
                ),
                
                actionButton("run_dt_enrich", "Run enrichment", icon = icon("play"))
            )
          ),
          
          column(
            width = 8,
            div(class = "panel-lite", uiOutput("dt_scope_note"))
          )
        ),
        
        tabsetPanel(
          tabPanel(
            "Disease",
            fluidRow(
              column(6, div(class = "panel-lite", DTOutput("enrich_disease_table"))),
              column(6, div(class = "panel-lite", plotly::plotlyOutput("enrich_disease_plot", height = 700)))
            )
          ),
          tabPanel(
            "Trait",
            fluidRow(
              column(6, div(class = "panel-lite", DTOutput("enrich_trait_table"))),
              column(6, div(class = "panel-lite", plotly::plotlyOutput("enrich_trait_plot", height = 700)))
            )
          )
        )
      ),
      
      # ------------------------------------
      # (B) Functional enrichment
      # ------------------------------------
      tabPanel(
        title = HTML("<span style='font-size:15px; font-weight:600;'>🧬 GO/KEGG/GoSlim</span>"),
        value = "tab_enrich_func",
        
        sidebarLayout(
          sidebarPanel(
            width = 3,
            
            radioButtons(
              "func_scope", "Scope",
              choices = c("Global" = "global", "Cluster" = "cluster"),
              selected = "global"
            ),
            
            conditionalPanel(
              condition = "input.func_scope == 'cluster'",
              uiOutput("func_cluster_ui")
            ),
            
            tags$hr(),
            
            selectInput(
              "enrich_background", "Background",
              choices = c("Reference annotated genes" = "orgdb", "Dataset genes" = "dataset"),
              selected = "orgdb"
            ),
            
            numericInput(
              "enrich_pcut",
              "FDR cutoff",
              value = 0.05,
              min = 0,
              max = 1,
              step = 0.01
            ),
            
            conditionalPanel(
              condition = "input.enrich_tabs == 'tab_enrich_go' || input.enrich_tabs == 'tab_enrich_goslim'",
              checkboxGroupInput(
                "go_ontos",
                "GO ontologies",
                choices = c("BP", "CC", "MF"),
                selected = c("BP", "CC", "MF"),
                inline = TRUE
              ),
              numericInput("go_topn", "Top terms per ontology", value = 10, min = 1, max = 50, step = 1)
            ),
            
            conditionalPanel(
              condition = "input.enrich_tabs == 'tab_enrich_go'",
              checkboxInput("go_simplify", "Simplify GO terms (optional)", value = FALSE)
            ),
            
            conditionalPanel(
              condition = "input.enrich_tabs == 'tab_enrich_kegg'",
              numericInput("enrich_kegg_top", "Top KEGG pathways", value = 15, min = 1, max = 50, step = 1)
            ),
            
            tags$hr(),
            
            actionButton("info_12", "ℹ️ GSSize"),
            numericInput("enrich_min_gs", "minGSSize", value = 10, min = 1, step = 1),
            numericInput("enrich_max_gs", "maxGSSize", value = 500, min = 10, step = 10),
            
            tags$hr(),
            
            actionButton(
              "run_enrich",
              label = "Run enrichment",
              icon = icon("play"),
              class = "btn btn-primary"
            )
          ),
          
          mainPanel(
            width = 9,
            div(class = "panel-lite", uiOutput("enrich_bg_note")),
            
            tabsetPanel(
              id = "enrich_tabs",
              
              tabPanel(
                title = HTML("<span style='font-size:15px; font-weight:600;'>GO</span>"),
                value = "tab_enrich_go",
                div(class = "panel-lite", withSpinner(plotlyOutput("go_bar", height = 400))),
                tags$hr(),
                div(class = "panel-lite", withSpinner(DT::DTOutput("go_table")))
              ),
              
              tabPanel(
                title = HTML("<span style='font-size:15px; font-weight:600;'>KEGG</span>"),
                value = "tab_enrich_kegg",
                div(class = "panel-lite", withSpinner(plotOutput("kegg_bar", height = 400))),
                tags$hr(),
                div(class = "panel-lite", withSpinner(DT::DTOutput("kegg_table")))
              ),
              
              tabPanel(
                title = HTML("<span style='font-size:15px; font-weight:600;'>GoSlim</span>"),
                value = "tab_enrich_goslim",
                div(class = "panel-lite", withSpinner(plotlyOutput("goslim_bar", height = "450px"))),
                tags$hr(),
                div(class = "panel-lite", withSpinner(DT::DTOutput("goslim_table")))
              )
            )
          )
        )
      )
    )
  ),
  
  # ========================================================================
  # TAB 3 — LD
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 LD</span>"),
    value = "ld_tab",
    ld_module_ui("ld", app_tag = "catalog")
  )
)


# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  #############################################################
  # Catalog Inspector (inside server): shared helpers (portable)
  source(file.path(cfg$shared, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
  source(file.path(cfg$shared, "gi_clusters_canonical.R"), local = TRUE)
  
  ############## Catalog as MASTER >>>> save file to share
  observeEvent(list(input$gwas_file, input$gwas_p_col, input$gwas_sep, input$gwas_header), {
    df2 <- gwas_df()
    req(is.data.frame(df2), nrow(df2) > 0)
    
    sid <- gi_sid(session)
    p   <- gi_state_paths(sid)
    
    saveRDS(df2, p$gwas)
    #∫∫∫∫∫∫∫∫∫∫∫∫∫∫   
    st0 <- gi_read_state(sid)
    if (is.null(st0) || !is.list(st0)) st0 <- list()
    
    st_new <- modifyList(st0, list(
      has_gwas = TRUE,
      gwas_rds = p$gwas
    ))
    
    gi_write_state(sid, st_new)
  }, ignoreInit = TRUE)
  #∫∫∫∫∫∫∫∫∫∫∫∫∫∫∫∫
  
  
  #### ---------------------------------------------------------------------------
  # info modals
  #### ---------------------------------------------------------------------------
  
  observeEvent(input$info_format, {
    txt <- HTML('
<p style="text-align: justify;">
  <b>Input file format [(*) mandatory columns]</b>
</p>

<table style="border-collapse: collapse; width: 100%;" border="1">
  <thead>
    <tr>
      <th>CHR(*)</th>
      <th>SNP(*)</th>
      <th>BP(*)</th>
      <th>A1</th>
      <th>TEST</th>
      <th>NMISS</th>
      <th>OR</th>
      <th>STAT</th>
      <th>P(*)</th>
    </tr>
  </thead>

  <tbody>
    <tr><td>14</td><td>14:33967212:C:A</td><td>33967212</td><td>C</td><td>ADD</td><td>720</td><td>2.03</td><td>4.84</td><td>1.30E-06</td></tr>
    <tr><td>18</td><td>rs12955421</td><td>76002716</td><td>A</td><td>ADD</td><td>720</td><td>2.393</td><td>4.76</td><td>1.93E-06</td></tr>
    <tr><td>3</td><td>rs1145036</td><td>2328334</td><td>G</td><td>ADD</td><td>720</td><td>1.816</td><td>4.62</td><td>3.83E-06</td></tr>
  </tbody>
</table>
')
    showModal(modalDialog(
      title = NULL,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "m",
      tags$div(style="line-height:1.35;", txt)
    ))
  })
  
  observeEvent(input$perf_info, {
    showModal(modalDialog(
      title = "Performance note (large files)",
      easyClose = TRUE,
      size = "m",
      HTML(paste0(
        "<b>GItools is tuned for large input tables.</b> Some inspectors apply an initial <b>P &lt; 0.05</b> filter and ",
        "<b>smart downsampling</b> (e.g., target 50k variants while always retaining peaks at P ≤ 1e-6 and preserving genome-wide coverage ",
        "via 10 Mb bins with a fixed seed) to keep plots and clustering fast and reproducible.<br><br>",
        "<b>What we are doing:</b><br>",
        "1) <b>Initial filter</b>: keep only variants with <b>P &lt; 0.05</b> (pcut = 0.05).<br>",
        "2) <b>Smart downsampling</b>: target 50k; keep peaks at P ≤ 1e-6; 10 Mb bins; seed 123."
      )),
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$info_01, {
    txt <- HTML(
      "<b>Clustering methods — quick guide</b><br><br>",
      
      "<p>This app can generate GWAS <b>clusters</b> using two alternative strategies. ",
      "Both end up producing a table of clusters with <b>cluster_id</b>, <b>start</b>, and <b>end</b> (plus summary stats such as top SNP, top −log10(P), and EWAS bins when available).</p>",
      
      "<hr style='margin:10px 0;'>",
      
      "<b>1) By intervals (hits → flank → merge)</b><br>",
      "<p style='margin-top:6px;'>",
      "Starts from GWAS hits above the chosen threshold (or only the selected hits in the table). ",
      "For each hit, an interval is built as <b>[BP − flank, BP + flank]</b>. ",
      "Overlapping intervals are merged into continuous regions, which become the final clusters.",
      "</p>",
      
      "<ul>",
      "<li><b>Input used</b>: hits above threshold (optionally only selected rows).</li>",
      "<li><b>Main parameter</b>: <b>flank</b> (bp) controls how far each hit expands.</li>",
      "<li><b>Result</b>: merged genomic regions, each assigned a <b>cluster_id</b>.</li>",
      "</ul>",
      
      "<hr style='margin:10px 0;'>",
      
      "<b>2) By hit count (windows → merge)</b><br>",
      "<p style='margin-top:6px;'>",
      "First filters GWAS to SNPs with <b>−log10(P) ≥ min_logp</b>. ",
      "Then it looks for regions where significant hits are <b>dense</b>, by counting how many hits fall inside genomic windows. ",
      "Windows that pass <b>min_hits</b> are considered significant, and overlapping/adjacent significant windows are merged into final clusters. ",
      "Each final cluster is then summarised (start/end, number of hits, top SNP by −log10(P)).",
      "</p>",
      
      "<ul>",
      "<li><b>Input used</b>: all SNPs passing <b>min_logp</b>.</li>",
      "<li><b>Main parameters</b>: <b>min_logp</b> and <b>min_hits</b> (plus window settings below).</li>",
      "<li><b>Result</b>: clusters representing genomic regions enriched in significant hits.</li>",
      "</ul>",
      
      "<div style='margin-top:8px;'><b>Hit-count window modes</b></div>",
      "<p style='margin-top:6px;'>Choose how windows are defined before merging them into clusters:</p>",
      
      "<ul>",
      "<li><b>1 Mb hit-span</b>: within each chromosome, hits are ordered by position and grouped so that all hits in a group lie within a maximum span of <b>1,000,000 bp</b> from the <i>first</i> hit of that group. ",
      "These are <b>non-overlapping groups</b> (a hit belongs to only one group). ",
      "<br><i>Why clusters can be small:</i> the cluster boundaries are <b>min(BP)</b> to <b>max(BP)</b> of the hits in the group, so if hits are close together the cluster can be only a few kb even though the rule allows up to 1 Mb.</li>",
      
      "<li><b>Tiled windows</b>: the chromosome is split into <b>non-overlapping tiles</b> of size <b>win_bp</b> (step = win_bp). ",
      "Counts hits per tile; tiles passing <b>min_hits</b> are merged if they touch/overlap.</li>",
      
      "<li><b>Sliding windows</b>: windows of size <b>win_bp</b> move along the chromosome with step <b>step_bp</b> (so windows <b>overlap</b> when step_bp &lt; win_bp). ",
      "Counts hits per window; overlapping significant windows are merged into clusters. ",
      "<br><i>Typical behavior:</i> sliding is more sensitive to local density peaks and often finds narrower clusters than tiled.</li>",
      "</ul>",
      
      "<br><p style='margin:0;'><i>Tip:</i> Use <b>By intervals</b> for “hit-centered regions”. Use <b>By hit count</b> (especially <b>sliding</b>) when you want “hit-dense regions”.</p>"
    )
    
    showModal(modalDialog(
      title = NULL,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "m",
      tags$div(style = "line-height:1.35;", txt)
    ))
  })
  
  observeEvent(input$info_12, {
    txt <- HTML(
      "<b>GSSize (Gene Set Size) — quick guide</b><br><br>",
      "<p>GSSize is the number of genes in a GO/KEGG term gene set (depends on the background).</p>",
      "<ul>",
      "<li><b>minGSSize</b>: removes very small sets.</li><br>",
      "<li><b>maxGSSize</b>: removes very large generic sets.</li><br>",
      "<li>Only sets with <b>minGSSize ≤ GSSize ≤ maxGSSize</b> are tested.</li>",
      "</ul>"
    )
    
    showModal(modalDialog(
      title = "GSSize (Gene Set Size)",
      txt,
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  #≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  workdir <- file.path(tempdir(), "catalog_inspector")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  
  # --- state ---
  
  catalog_all  <- reactiveVal(NULL)
  
  # FINAL catalog hits path (unified new/old, NonSyn-style)
  catalog_final_path_csv <- reactiveVal(NULL)
  catalog_final_path_rds <- reactiveVal(NULL)
  
  
  # ---------------------------------------------------------------------------
  # Step 1: GWAS reader
  # ---------------------------------------------------------------------------
  
  rv <- reactiveValues(
    gwas = NULL,
    gwas_header_df = NULL,
    p_col_resolved = NULL,
    p_col_needs_user = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Lightweight preview/header read
  # ---------------------------------------------------------------------------
  
  gwas_preview <- reactive({
    req(input$gwas_file)
    
    tryCatch({
      readr::read_delim(
        input$gwas_file$datapath,
        delim = input$gwas_sep %||% "\t",
        col_names = isTRUE(input$gwas_header),
        show_col_types = FALSE,
        progress = FALSE,
        n_max = 200
      )
    }, error = function(e) {
      NULL
    })
  })
  
  observeEvent(
    list(input$gwas_file, input$gwas_sep, input$gwas_header),
    {
      req(input$gwas_file)
      
      rv$gwas <- NULL
      rv$gwas_header_df <- NULL
      rv$p_col_resolved <- NULL
      rv$p_col_needs_user <- FALSE
      
      f   <- input$gwas_file$datapath
      sep <- input$gwas_sep %||% "\t"
      hdr <- isTRUE(input$gwas_header)
      
      df0 <- tryCatch(
        data.table::fread(
          f,
          sep = sep,
          header = hdr,
          nrows = 0L,
          showProgress = FALSE,
          data.table = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(df0) || !is.data.frame(df0)) {
        showNotification("GWAS: cannot read file header.", type = "error")
        return()
      }
      
      rv$gwas_header_df <- df0
      
      auto_p <- pick_col(
        df0,
        c("P", "p", "PVAL", "pval", "P_VALUE", "p_value", "PVALUE", "pvalue")
      )
      
      if (!is.null(auto_p)) {
        rv$p_col_resolved <- auto_p
        rv$p_col_needs_user <- FALSE
      } else {
        rv$p_col_resolved <- NULL
        rv$p_col_needs_user <- TRUE
      }
    },
    ignoreInit = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Show selector only when auto-detection failed
  # ---------------------------------------------------------------------------
  
  output$gwas_p_selector <- renderUI({
    p <- gwas_preview()
    
    if (is.null(p) || !nrow(p)) {
      return(helpText("Could not read GWAS preview."))
    }
    
    if (!isTRUE(rv$p_col_needs_user)) {
      return(NULL)
    }
    
    cols <- names(p)
    
    tagList(
      div(
        style = paste(
          "margin-bottom:10px; padding:10px 12px; border-radius:8px;",
          "background:#fff3cd; color:#856404; border:1px solid #ffe69c;"
        ),
        HTML(
          "<b>P-value column not detected automatically.</b><br/>
         Please select the column containing raw p-values before continuing."
        )
      ),
      selectInput(
        "gwas_p_col",
        "Select p-value column",
        choices = cols,
        selected = cols[1]
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # Reactive that decides when the p column is truly resolved
  # ---------------------------------------------------------------------------
  
  gwas_p_col_final <- reactive({
    req(rv$gwas_header_df)
    
    hdr_names <- names(rv$gwas_header_df)
    
    if (isTRUE(rv$p_col_needs_user)) {
      req(input$gwas_p_col)
      validate(need(input$gwas_p_col %in% hdr_names, "Selected p-value column is not in file header."))
      return(input$gwas_p_col)
    }
    
    req(rv$p_col_resolved)
    rv$p_col_resolved
  })
  
  # ---------------------------------------------------------------------------
  # Full big-file read ONLY after p-value column is resolved
  # ---------------------------------------------------------------------------
  
  observeEvent(
    gwas_p_col_final(),
    {
      req(input$gwas_file)
      req(rv$gwas_header_df)
      
      p_col <- gwas_p_col_final()
      
      # Modal/loading starts ONLY here
      session$sendCustomMessage("gwas_loading", TRUE)
      on.exit(session$sendCustomMessage("gwas_loading", FALSE), add = TRUE)
      
      rv$gwas <- NULL
      
      f   <- input$gwas_file$datapath
      sep <- input$gwas_sep %||% "\t"
      hdr <- isTRUE(input$gwas_header)
      df0 <- rv$gwas_header_df
      
      chr_col <- pick_col(df0, c("CHR", "chr", "chrom", "CHROM", "chromosome"))
      bp_col  <- pick_col(df0, c("BP", "bp", "POS", "pos", "position", "POSITION"))
      snp_col <- pick_col(df0, c("SNP", "snp", "rsid", "RSID", "marker", "ID", "id"))
      
      if (is.null(chr_col)) {
        showNotification("GWAS: missing chromosome column.", type = "error")
        return()
      }
      
      if (is.null(bp_col)) {
        showNotification("GWAS: missing position column.", type = "error")
        return()
      }
      
      if (is.null(p_col)) {
        showNotification("GWAS: missing p-value column.", type = "error")
        return()
      }
      
      sel <- unique(c(chr_col, bp_col, p_col, snp_col))
      sel <- sel[!is.na(sel) & nzchar(sel)]
      
      dt <- tryCatch(
        data.table::fread(
          f,
          sep = sep,
          header = hdr,
          select = sel,
          showProgress = TRUE
        ),
        error = function(e) NULL
      )
      
      if (is.null(dt) || !nrow(dt)) {
        showNotification("GWAS: empty file or failed reading selected columns.", type = "error")
        return()
      }
      
      dt <- data.table::as.data.table(dt)
      
      rawP <- dt[[p_col]]
      BP   <- if (is.numeric(dt[[bp_col]])) {
        dt[[bp_col]]
      } else {
        suppressWarnings(readr::parse_number(as.character(dt[[bp_col]])))
      }
      
      Pval <- parse_p_robust(rawP)
      CHR  <- chr_map_plink19(dt[[chr_col]])
      
      snp <- if (!is.null(snp_col) && snp_col %in% names(dt)) {
        as.character(dt[[snp_col]])
      } else {
        paste0("chr", norm_chr_generic(dt[[chr_col]]), ":", BP)
      }
      
      out <- data.table::data.table(
        CHR  = as.integer(CHR),
        BP   = as.integer(BP),
        snp  = snp,
        Pval = as.numeric(Pval)
      )
      
      rm(dt, rawP, CHR, BP, Pval, snp)
      gc(FALSE)
      
      out <- out[
        is.finite(CHR) &
          is.finite(BP) &
          is.finite(Pval) &
          Pval > 0 &
          Pval <= 1
      ]
      
      if (!nrow(out)) {
        showNotification(
          "GWAS: after parsing, no valid rows remained. Check the selected p-value column.",
          type = "error",
          duration = 8
        )
        return()
      }
      
      # Initial filter + smart downsampling
      pcut         <- 0.05
      target_n     <- 50000L
      peaks_keep_p <- 1e-6
      bin_bp       <- 1e7
      seed         <- 123L
      
      out <- out[Pval < pcut]
      
      if (!nrow(out)) {
        showNotification(
          "GWAS: no variants passed the initial filter P < 0.05.",
          type = "warning",
          duration = 8
        )
        rv$gwas <- data.table::data.table(
          CHR = integer(),
          BP = integer(),
          snp = character(),
          Pval = numeric(),
          logp = numeric()
        )
        return()
      }
      
      if (nrow(out) > target_n) {
        set.seed(seed)
        
        peaks <- out[Pval <= peaks_keep_p]
        
        if (nrow(peaks) >= target_n) {
          data.table::setorder(peaks, Pval)
          out <- peaks[seq_len(target_n)]
        } else {
          rest <- out[Pval > peaks_keep_p]
          rest[, bin := paste0(CHR, "_", floor(BP / bin_bp))]
          
          remaining_n <- target_n - nrow(peaks)
          bin_sizes <- rest[, .N, by = bin]
          bin_sizes[, alloc := floor(remaining_n * (N / sum(N)))]
          
          if (remaining_n >= nrow(bin_sizes)) {
            bin_sizes[alloc == 0, alloc := 1L]
          }
          
          tot <- sum(bin_sizes$alloc)
          
          if (tot > remaining_n) {
            over <- tot - remaining_n
            data.table::setorder(bin_sizes, -alloc)
            i <- 1L
            while (over > 0L) {
              if (bin_sizes$alloc[i] > 0L) {
                bin_sizes$alloc[i] <- bin_sizes$alloc[i] - 1L
                over <- over - 1L
              }
              i <- i + 1L
              if (i > nrow(bin_sizes)) i <- 1L
            }
          } else if (tot < remaining_n) {
            under <- remaining_n - tot
            data.table::setorder(bin_sizes, -alloc)
            i <- 1L
            while (under > 0L) {
              bin_sizes$alloc[i] <- bin_sizes$alloc[i] + 1L
              under <- under - 1L
              i <- i + 1L
              if (i > nrow(bin_sizes)) i <- 1L
            }
          }
          
          picked <- rest[, {
            k <- bin_sizes[bin == .BY$bin, alloc][1]
            if (is.na(k) || k <= 0L) .SD[0]
            else if (.N <= k) .SD
            else .SD[sample.int(.N, k)]
          }, by = bin]
          
          picked[, bin := NULL]
          out <- data.table::rbindlist(list(peaks, picked), use.names = TRUE, fill = TRUE)
        }
      }
      
      out[, logp := -log10(pmax(Pval, .Machine$double.xmin))]
      
      rv$gwas <- out
    },
    ignoreInit = TRUE
  )
  
  gwas_df <- reactive({
    req(rv$gwas)
    rv$gwas
  })
  #==========================================================================
  
  min_hits_active <- reactive({
    if (identical(input$cluster_method, "window")) {
      input$min_hits_window %||% 3
    } else {
      input$min_hits_hits %||% 3
    }
  })
  
  #==========================================================================
  
  # --- IMPORTANT: only allow build_clusters in NEW mode (so canonical observer doesn't run in OLD) ---
  gwas_df_for_clustering <- reactive({
    gwas_df()
  })
  
  # --- Canonical clusters engine (same as GTEx) ---
  clusters_engine <- gi_clusters_canonical_init(
    session = session, input = input, output = output,
    gwas_df = gwas_df_for_clustering,
    build_btn_id   = "build_ranges",
    clusters_dt_id = "cluster_dt",
    hits_rows_id   = "hits_tbl_rows_selected",
    app_count_col  = "n_catalog"
  )
  
  # ---------------------------------------------------------------------------
  # SAVE GWAS (Catalog as MASTER → expose to slaves)
  # ---------------------------------------------------------------------------
  
  observeEvent(rv$gwas, {
    req(rv$gwas)
    
    sid <- gi_sid(session)
    if (is.null(sid) || !nzchar(sid)) return()
    
    gwas_df <- data.table::copy(rv$gwas)
    
    # --- compatibilitat: assegurar columna P ---
    if (!"P" %in% names(gwas_df) && "Pval" %in% names(gwas_df)) {
      gwas_df[, P := Pval]
    }
    
    # --- path ---
    p <- gi_state_paths(sid)
    gwas_path <- file.path(dirname(p$params), paste0("gwas_", sid, ".rds"))
    
    # --- save ---
    saveRDS(gwas_df, gwas_path)
    
    cat("[SAVE] GWAS saved:", gwas_path, "\n")
    
    # --- update state.json (CLAU!) ---
    st0 <- gi_read_state(sid)
    if (!is.list(st0)) st0 <- list()
    
    st_new <- modifyList(st0, list(
      gwas_rds  = gwas_path,
      stamp     = gi_bump_stamp(st0$stamp %||% 0),
      updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    
    gi_write_state(sid, st_new)
    
  }, ignoreInit = TRUE)
  
  # Catalog as MASTER, read params for slaves
  observeEvent(
    list(input$cluster_method, input$hits_mode,
         input$pthr, input$min_logp,
         input$flank, input$win_bp, input$step_bp,
         input$min_hits_window, input$min_hits_hits),
    {
      sid <- gi_sid(session)
      if (is.null(sid) || !nzchar(sid)) return()
      
      min_hits_value <- if (identical(input$cluster_method, "window")) {
        input$min_hits_window %||% 3
      } else {
        input$min_hits_hits %||% 3
      }
      
      write_params_rds(
        sid = sid,
        input = input,
        min_hits_value = as.integer(min_hits_active() %||% 3L)
      )
      
      p   <- gi_state_paths(sid)
      st0 <- gi_read_state(sid)
      if (!is.list(st0)) st0 <- list()
      
      st_new <- modifyList(st0, list(
        params_rds = p$params,
        stamp      = gi_bump_stamp(st0$stamp %||% 0),
        updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ))
      
      gi_write_state(sid, st_new)
    },
    ignoreInit = TRUE
  )
  
  
  # Replace Catalog's clusters_val with the canonical one
  clusters_val <- clusters_engine$clusters_cur
  gitools_deeplink_catalog(session, clusters_val)
  
  observeEvent(input$build_ranges, {
    
    # Reset downstream (Catalog outputs)
    catalog_final_path_csv(NULL)
    catalog_final_path_rds(NULL)
    
    # Read clusters computed by the canonical engine
    clusters <- clusters_val()
    validate(need(is.data.frame(clusters) && nrow(clusters) > 0, "No clusters generated."))
    
    # Ensure workdir exists
    validate(need(exists("workdir") && nzchar(workdir), "workdir is not defined."))
    if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
    
    write_range_file <- function(clusters_df) {
      cl <- as.data.frame(clusters_df)
      if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
      if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
      if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
      
      validate(need(all(c("chr","start","end","cluster_id") %in% names(cl)),
                    "clusters must contain chr/start/end/cluster_id."))
      
      lines <- sprintf(
        "%s %d %d %s",
        chr_label_plink(as.integer(cl$chr)),
        as.integer(cl$start),
        as.integer(cl$end),
        as.character(cl$cluster_id)
      )
      
      rng_path <- file.path(workdir, "selected_intervals.range")
      writeLines(lines, rng_path)
      list(lines = lines, path = rng_path)
    }
    
    sid <- gi_sid(session)
    p   <- gi_state_paths(sid)
    write_params_rds(
      sid = sid,
      input = input,
      min_hits_value = min_hits_active()
    )
    
    ok_save <- TRUE
    tryCatch({
      saveRDS(clusters, p$clus)
    }, error = function(e) {
      ok_save <<- FALSE
      cat("[GI][MASTER] could not save clusters rds:", conditionMessage(e), "\n")
    })
    validate(need(ok_save && file.exists(p$clus), "Failed to save clusters RDS (sync cannot proceed)."))
    
    st0 <- gi_read_state(sid)
    if (is.null(st0) || !is.list(st0)) st0 <- list()
    
    new_stamp <- gi_bump_stamp(st0$stamp %||% 0)
    
    st_new <- modifyList(st0, list(
      stamp        = new_stamp,
      sid          = sid,
      updated_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      has_gwas     = TRUE,
      gwas_rds     = p$gwas,
      has_clusters = TRUE,
      clusters_rds = p$clus,
      clus_rds     = p$clus,
      params_rds   = p$params,
      cluster_method = as.character(input$cluster_method %||% "window"),
      hits_mode      = as.character(input$hits_mode %||% NA_character_),
      thr_type  = if (identical(input$cluster_method %||% "window", "window")) "pthr" else "min_logp",
      thr_value = if (identical(input$cluster_method %||% "window", "window")) {
        as.numeric(input$pthr %||% NA_real_)
      } else {
        as.numeric(input$min_logp %||% NA_real_)
      },
      pthr     = as.numeric(input$pthr %||% NA_real_),
      min_logp = as.numeric(input$min_logp %||% NA_real_),
      flank_bp = as.integer(input$flank %||% 0L),
      flank    = as.integer(input$flank %||% 0L),
      min_hits = as.integer(min_hits_active() %||% NA_integer_),
      win_bp   = as.integer(input$win_bp %||% NA_integer_),
      step_bp  = as.integer(input$step_bp %||% NA_integer_)
    ))
    
    gi_write_state(sid, st_new)
    
    wf <- write_range_file(clusters)
    
    output$ranges_preview <- renderText({
      paste0(
        "thr=", input$pthr, " | flank=", input$flank, " bp\n",
        "Clusters: ", nrow(clusters), "\n\n",
        paste(head(wf$lines, 10), collapse = "\n")
      )
    })
    
    showNotification(
      "Clusters generated (canonical). Run Step 4 to assign Catalog hits.",
      type = "message", duration = 4
    )
    
  }, ignoreInit = TRUE, priority = -100)
  
  
  # ---------------------------------------------------------------------------
  # Manhattan helpers
  # ---------------------------------------------------------------------------
  
  dfp_manhattan <- reactive({
    df <- gwas_df(); req(nrow(df) > 0)
    # filter pval < 0.0 to a easely plot
    df <- df %>% dplyr::filter(Pval < 0.05)
    
    ref <- .ref_hg38
    df$chrN <- norm_chr_generic(df$CHR)
    
    df %>%
      inner_join(ref %>% select(chr, chr_cum), by = c("chrN" = "chr")) %>%
      arrange(chrN, BP) %>%
      mutate(BPcum = BP + chr_cum)
  })
  
  axis_df <- reactive({
    dfp <- dfp_manhattan()
    dfp %>% group_by(chrN) %>% summarise(center = mean(BPcum, na.rm = TRUE), .groups = "drop")
  })
  
  # hits_df <- reactive({
  #   df <- gwas_df(); req(nrow(df) > 0)
  #   req(input$pthr)
  #   df %>% filter(logp >= input$pthr) %>% arrange(desc(logp)) %>% select(CHR, BP, snp, p = Pval, logp)
  # })
  
  hits_df <- reactive({
    df <- gwas_df()
    req(nrow(df) > 0)
    
    cluster_method <- as.character(input$cluster_method %||% "")
    hits_mode      <- as.character(input$hits_mode %||% "")
    
    mode_resolved <- dplyr::case_when(
      cluster_method %in% c("window") ~ "window",
      cluster_method %in% c("hits") & hits_mode %in% c("span1mb", "tiled", "sliding") ~ hits_mode,
      cluster_method %in% c("span1mb", "tiled", "sliding") ~ cluster_method,
      TRUE ~ NA_character_
    )
    
    thr <- dplyr::case_when(
      mode_resolved == "window" ~ suppressWarnings(as.numeric(input$pthr)),
      mode_resolved %in% c("span1mb", "tiled", "sliding") ~ suppressWarnings(as.numeric(input$min_logp)),
      TRUE ~ NA_real_
    )
    
    req(!is.na(mode_resolved), is.finite(thr))
    
    out <- df %>%
      dplyr::filter(is.finite(logp), logp >= thr) %>%
      dplyr::arrange(dplyr::desc(logp)) %>%
      dplyr::select(CHR, BP, snp, p = Pval, logp)
    
    cat(
      "[HITS_DF][RESOLVED] cluster_method=", cluster_method,
      " | hits_mode=", hits_mode,
      " | mode_resolved=", mode_resolved,
      " | threshold=", thr,
      " | nrow=", nrow(out), "\n",
      sep = ""
    )
    
    out
  })
  
  
  hits_enriched <- reactive({
    df <- hits_df(); req(nrow(df) > 0)
    rscol <- detect_rsid_col(df) %||% "snp"
    if ("p" %in% names(df)) df$p <- formatC(df$p, format = "e", digits = 2)
    if ("logp" %in% names(df)) df$logp <- sprintf("%.2f", df$logp)
    df$dbSNP <- paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/", df[[rscol]], "' target='_blank'>", df[[rscol]], "</a>")
    df
  })
  
  output$hits_tbl <- renderDT({
    datatable(hits_enriched(), selection = "multiple", rownames = FALSE, escape = FALSE,
              extensions = "Buttons",
              options    = list(
                dom        = "Bfrtip",
                buttons = list(
                  list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                  list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                  list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                ),
                pageLength = 10,
                scrollX    = TRUE
              ))
  }, server = FALSE)
  
  # ---------------------------------------------------------------------------
  # Step 3: load catalog at startup
  # ---------------------------------------------------------------------------
  
  observe({
    path <- "www/gwas_catalog_simplified.rds"
    if (file.exists(path)) {
      df <- tryCatch(readRDS(path), error = function(e) NULL)
      if (!is.null(df) && is.data.frame(df)) {
        catalog_all(df)
        output$catalog_status <- renderText(sprintf("Catalog loaded: %s (rows: %d, cols: %d)", path, nrow(df), ncol(df)))
      } else {
        output$catalog_status <- renderText("Could not read www/gwas_catalog_simplified.rds.")
      }
    } else {
      output$catalog_status <- renderText("Catalog file not found: www/gwas_catalog_simplified.rds")
    }
  })
  
  # ---------------------------------------------------------------------------
  # Unified reactive: catalog hits dataset used everywhere (NonSyn-style)
  # - OLD mode: read uploaded CSV/RDS
  # - NEW mode: read FINAL paths (prefer RDS)
  # ---------------------------------------------------------------------------
  
  catalog_hits_df <- reactive({
    rds <- catalog_final_path_rds()
    csv <- catalog_final_path_csv()
    
    if ((is.null(rds) || !file.exists(rds)) && (is.null(csv) || !file.exists(csv))) {
      return(tibble::tibble())
    }
    
    df <- NULL
    
    if (!is.null(rds) && file.exists(rds)) {
      df <- tryCatch(readRDS(rds), error = function(e) NULL)
    } else if (!is.null(csv) && file.exists(csv)) {
      df <- tryCatch(utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    }
    
    validate(need(is.data.frame(df), "Could not read final Catalog hits file."))
    
    if (!"cluster_id" %in% names(df)) {
      return(tibble::tibble())
    }
    
    df
  })
  
  # canonical catalog_data for DT/plots
  catalog_data <- reactive({
    dt <- catalog_hits_df()
    # print(head(dt))
    if (is.null(dt) || !is.data.frame(dt)) tibble() else dt
  })
  
  
  # ---------------------------------------------------------------------------
  # Step 4: assign catalog hits (new mode only) + persist final
  # ---------------------------------------------------------------------------
  
  catalog_log_text <- reactiveVal("")
  
  append_catalog_log <- function(...) {
    old <- catalog_log_text()
    catalog_log_text(paste(c(old, paste(..., collapse = " ")), collapse = "\n"))
  }
  
  # logger generic per a integrator
  append_log <- append_catalog_log
  
  output$catalog_log <- renderText(catalog_log_text())
  
  
  observeEvent(input$run_catalog, {
    
    req(catalog_all())
    req(clusters_val())
    
    # reset log
    catalog_log_text("")
    append_catalog_log("[START] run_catalog")
    
    if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::disable("run_catalog")
    
    tryCatch({
      
      hits <- NULL
      cl2  <- NULL
      
      withProgress(message = "Assigning GWAS Catalog hits to clusters…", value = 0, {
        
        incProgress(0.05, detail = "Loading inputs")
        cl <- clusters_val()
        cat_raw <- catalog_all()
        append_catalog_log(sprintf("[1] Inputs loaded | clusters=%d | catalog=%d", nrow(cl), nrow(cat_raw)))
        
        incProgress(0.20, detail = "Preparing Catalog resource")
        cat_pre <- catalog_prepare(cat_raw)
        append_catalog_log("[2] Catalog prepared")
        
        # free memory early
        rm(cat_raw)
        gc()
        
        incProgress(0.45, detail = "Extracting hits per cluster")
        hits <- extract_catalog_by_clusters(cl, cat_pre)
        append_catalog_log(sprintf("[3] Hits extracted | hits=%d", nrow(hits)))
        
        append_catalog_log(paste0("[DEBUG] hits cols: ", paste(names(hits), collapse = ", ")))
        append_catalog_log(
          paste0(
            "[DEBUG] hits head:\n",
            paste(capture.output(print(utils::head(hits, 5))), collapse = "\n")
          )
        )
        
        rm(cat_pre)
        gc()
        
        incProgress(0.15, detail = "Updating n_catalog in clusters")
        cl2 <- add_n_catalog(cl, hits)
        clusters_val(cl2)
        append_catalog_log("[4] clusters updated (n_catalog)")
        
        output$catalog_summary <- renderText({
          sprintf("Catalog hits assigned: %d (clusters: %d)", nrow(hits), nrow(cl2))
        })
        
        incProgress(0.10, detail = "Saving outputs (CSV + RDS)")
        tg <- make_mode_thr_tag(input$cluster_method %||% "window", input$pthr, input$min_logp)
        stem <- paste0("catalog_hits_", tg$mode_tag, "_thr", tg$thr_txt, "_with_clusters")
        out_csv <- file.path(workdir, paste0(stem, ".csv"))
        out_rds <- file.path(workdir, paste0(stem, ".rds"))
        
        # ensure cluster columns exist
        if (!"cluster_id" %in% names(hits)) hits$cluster_id <- NA_character_
        
        saveRDS(hits, out_rds)
        utils::write.csv(hits, out_csv, row.names = FALSE)
        
        append_catalog_log(paste0("[DEBUG][RDS] out_rds path: ", out_rds))
        append_catalog_log(paste0("[DEBUG][RDS] file exists: ", file.exists(out_rds)))
        append_catalog_log(paste0("[DEBUG][RDS] hits nrow=", nrow(hits), " ncol=", ncol(hits)))
        append_catalog_log(paste0("[DEBUG][RDS] hits cols: ", paste(names(hits), collapse = ", ")))
        append_catalog_log(
          paste0(
            "[DEBUG][RDS] hits head:\n",
            paste(capture.output(print(utils::head(hits, 5))), collapse = "\n")
          )
        )
        
        catalog_final_path_rds(out_rds)
        catalog_final_path_csv(out_csv)
        
        append_catalog_log(paste0("[5] Saved: ", basename(out_csv), " | ", basename(out_rds)))
        
        ###### Start INTEGRATOR CODE ####################################################################
        # ------------------------------------------------------------
        # INTEGRATOR export (generic app bridges + LD inputs)
        # ------------------------------------------------------------
        tryCatch({
          
          # ============================================================
          # APP CONFIG  >>> ONLY CHANGE THIS BLOCK IN EACH APP
          # ============================================================
          app_slug      <- "catalog"
          app_hit_class <- "catalog_hit"
          
          # Generic logger expected in each app
          # Example: append_log <- append_catalog_log
          if (!exists("append_log", mode = "function", inherits = TRUE)) {
            stop("append_log() not found. Define a generic append_log() in this app.")
          }
          
          # Column mapping for THIS app
          col_gene    <- "MAPPED_GENE"
          col_trait   <- "MAPPED_TRAIT"
          col_disease <- "DISEASE"
          col_chr     <- "CHR"
          col_pos     <- "POS"
          col_id      <- "rsid"
          
          # Shared names used in all apps
          col_cluster_id    <- "cluster_id"
          col_cluster_start <- "cluster_start"
          col_cluster_end   <- "cluster_end"
          
          log_tag <- paste0("[INTEGRATOR][", app_slug, "]")
          ld_tag  <- paste0("[LD-INTEGRATOR][", app_slug, "]")
          
          `%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b
          
          # ============================================================
          # HELPERS
          # ============================================================
          fmt_threshold_for_path <- function(x) {
            x <- suppressWarnings(as.numeric(x))
            if (is.na(x)) return("NA")
            out <- format(x, scientific = FALSE, trim = TRUE)
            gsub("\\.", "p", out)
          }
          
          fmt_tag <- function(x) {
            x <- as.character(x %||% "NA")
            x <- trimws(x)
            x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
            if (!nzchar(x)) x <- "NA"
            x
          }
          
          get_catalog_cluster_method <- function(input) {
            cm <- as.character(input$cluster_method %||% NA_character_)
            hm <- as.character(input$hits_mode %||% NA_character_)
            
            if (identical(cm, "window")) {
              "window"
            } else if (identical(cm, "hits")) {
              if (identical(hm, "span1mb")) {
                "hits_span1mb"
              } else if (identical(hm, "tiled")) {
                "hits_tiled"
              } else if (identical(hm, "sliding")) {
                "hits_sliding"
              } else {
                "hits_unknown"
              }
            } else {
              "unknown"
            }
          }
          
          get_catalog_threshold <- function(input, cluster_method) {
            if (identical(cluster_method, "window")) {
              suppressWarnings(as.numeric(input$pthr))
            } else if (cluster_method %in% c("hits_span1mb", "hits_tiled", "hits_sliding", "hits_unknown")) {
              suppressWarnings(as.numeric(input$min_logp))
            } else {
              NA_real_
            }
          }
          
          get_catalog_gwas_name <- function(input) {
            nm <- tryCatch(as.character(input$gwas_file$name %||% NA_character_), error = function(e) NA_character_)
            if (!is.na(nm) && nzchar(trimws(nm))) return(trimws(nm))
            "unknown_gwas_source"
          }
          
          gi_get_integrator_session_dir <- function(gi_shared_root, cluster_method, threshold_used) {
            integrator_root <- file.path(gi_shared_root, "integrator_exports")
            dir.create(integrator_root, recursive = TRUE, showWarnings = FALSE)
            
            method_tag <- fmt_tag(cluster_method)
            thr_tag    <- fmt_threshold_for_path(threshold_used)
            
            session_id <- paste0(
              "clust_", method_tag,
              "__thr_", thr_tag
            )
            
            session_dir <- file.path(integrator_root, session_id)
            dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
            
            list(
              integrator_root = integrator_root,
              session_id = session_id,
              session_dir = session_dir,
              method_tag = method_tag,
              thr_tag = thr_tag
            )
          }
          
          gi_load_or_create_manifest <- function(session_dir, session_id, cluster_method, threshold_used, gwas_session_file) {
            manifest_path <- file.path(session_dir, "manifest.rds")
            
            if (file.exists(manifest_path)) {
              manifest <- tryCatch(readRDS(manifest_path), error = function(e) NULL)
              if (is.list(manifest)) {
                if (is.null(manifest$files)) manifest$files <- list()
                if (is.null(manifest$apps_present)) manifest$apps_present <- character(0)
                return(list(manifest = manifest, manifest_path = manifest_path))
              }
            }
            
            manifest <- list(
              session_id = session_id,
              created_at = as.character(Sys.time()),
              cluster_method = cluster_method,
              threshold_used = threshold_used,
              gwas_session_file = gwas_session_file,
              settings_key = paste0(
                "cluster_method=", cluster_method,
                ";threshold_used=", threshold_used
              ),
              apps_present = character(0),
              files = list()
            )
            
            list(manifest = manifest, manifest_path = manifest_path)
          }
          
          # ============================================================
          # SESSION / SETTINGS FOLDER
          # ============================================================
          current_cluster_method <- get_catalog_cluster_method(input)
          current_threshold      <- get_catalog_threshold(input, current_cluster_method)
          current_gwas_name      <- get_catalog_gwas_name(input)
          
          sess <- gi_get_integrator_session_dir(
            gi_shared_root = gi_shared_root,
            cluster_method = current_cluster_method,
            threshold_used = current_threshold
          )
          
          integrator_root <- sess$integrator_root
          session_id      <- sess$session_id
          session_dir     <- sess$session_dir
          method_tag      <- sess$method_tag
          thr_tag         <- sess$thr_tag
          
          manifest_obj <- gi_load_or_create_manifest(
            session_dir = session_dir,
            session_id = session_id,
            cluster_method = current_cluster_method,
            threshold_used = current_threshold,
            gwas_session_file = current_gwas_name
          )
          
          manifest      <- manifest_obj$manifest
          manifest_path <- manifest_obj$manifest_path
          
          # Validate existing manifest to avoid accidental mixing
          if (!identical(as.character(manifest$cluster_method %||% ""), as.character(current_cluster_method))) {
            stop("Existing manifest cluster_method does not match current cluster_method.")
          }
          
          if (!isTRUE(all.equal(
            suppressWarnings(as.numeric(manifest$threshold_used)),
            suppressWarnings(as.numeric(current_threshold))
          ))) {
            stop("Existing manifest threshold_used does not match current threshold_used.")
          }
          
          # Keep GWAS file recorded in manifest as a clean character vector
          manifest_gwas_files <- as.character(manifest$gwas_session_file %||% character(0))
          manifest_gwas_files <- trimws(manifest_gwas_files)
          manifest_gwas_files <- manifest_gwas_files[!is.na(manifest_gwas_files) & nzchar(manifest_gwas_files)]
          
          current_gwas_name_chr <- trimws(as.character(current_gwas_name %||% ""))
          if (nzchar(current_gwas_name_chr)) {
            manifest_gwas_files <- unique(c(manifest_gwas_files, current_gwas_name_chr))
          }
          
          manifest$gwas_session_file <- manifest_gwas_files
          
          gene_bridge_path     <- file.path(session_dir, paste0(app_slug, "_gene_bridge.rds"))
          term_bridge_path     <- file.path(session_dir, paste0(app_slug, "_term_bridge.rds"))
          clusters_master_path <- file.path(session_dir, paste0(app_slug, "_clusters_master.rds"))
          candidates_path      <- file.path(session_dir, paste0(app_slug, "_candidates.rds"))
          gwas_sig_bridge_path  <- file.path(session_dir, "gwas_significance_bridge.rds")
          
          append_log(paste0(log_tag, " integrator_root: ", integrator_root))
          append_log(paste0(log_tag, " cluster_method: ", current_cluster_method))
          append_log(paste0(log_tag, " threshold_used: ", current_threshold))
          append_log(paste0(log_tag, " gwas_session_file: ", paste(manifest$gwas_session_file, collapse = " | ")))
          append_log(paste0(log_tag, " session_id: ", session_id))
          append_log(paste0(log_tag, " session_dir: ", session_dir))
          append_log(paste0(log_tag, " hits nrow=", nrow(hits), " ncol=", ncol(hits)))
          append_log(paste0(log_tag, " hits cols: ", paste(names(hits), collapse = ", ")))
          
          # ============================================================
          # A) GENE BRIDGE
          # ============================================================
          gene_bridge <- tibble::tibble()
          
          if (col_gene %in% names(hits)) {
            gene_bridge <- hits %>%
              dplyr::transmute(
                gene = as.character(.data[[col_gene]]),
                source_app = app_slug,
                evidence_type = paste0(app_slug, "_gene"),
                cluster_id = as.character(.data[[col_cluster_id]]),
                chr = .data[[col_chr]],
                start = .data[[col_cluster_start]],
                end = .data[[col_cluster_end]]
              ) %>%
              sanitize_bridge() %>%
              sanitize_bridge_genes() %>%
              dplyr::distinct()
          } else {
            append_log(paste0(log_tag, " gene bridge skipped: column not found -> ", col_gene))
          }
          
          append_log(paste0(log_tag, " gene_bridge nrow=", nrow(gene_bridge)))
          append_log(
            paste0(
              log_tag, " gene_bridge head:\n",
              paste(capture.output(print(utils::head(gene_bridge, 5))), collapse = "\n")
            )
          )
          
          if (nrow(gene_bridge) > 0) {
            saveRDS(gene_bridge, gene_bridge_path)
            manifest$files[[paste0(app_slug, "_gene_bridge")]] <- basename(gene_bridge_path)
            append_log(paste0(log_tag, " saved ", basename(gene_bridge_path), " | exists=", file.exists(gene_bridge_path)))
          } else {
            append_log(paste0(log_tag, " ", basename(gene_bridge_path), " NOT saved (0 rows)"))
          }
          
          # ============================================================
          # B) TERM BRIDGE
          # ============================================================
          term_parts <- list()
          
          if (col_trait %in% names(hits)) {
            term_parts[[length(term_parts) + 1]] <- hits %>%
              dplyr::transmute(
                term = as.character(.data[[col_trait]]),
                term_type = "trait",
                source_app = app_slug,
                evidence_type = paste0(app_slug, "_trait"),
                cluster_id = as.character(.data[[col_cluster_id]]),
                chr = .data[[col_chr]],
                start = .data[[col_cluster_start]],
                end = .data[[col_cluster_end]]
              )
          } else {
            append_log(paste0(log_tag, " trait bridge skipped: column not found -> ", col_trait))
          }
          
          if (col_disease %in% names(hits)) {
            term_parts[[length(term_parts) + 1]] <- hits %>%
              dplyr::transmute(
                term = as.character(.data[[col_disease]]),
                term_type = "disease",
                source_app = app_slug,
                evidence_type = paste0(app_slug, "_disease"),
                cluster_id = as.character(.data[[col_cluster_id]]),
                chr = .data[[col_chr]],
                start = .data[[col_cluster_start]],
                end = .data[[col_cluster_end]]
              )
          } else {
            append_log(paste0(log_tag, " disease bridge skipped: column not found -> ", col_disease))
          }
          
          term_bridge <- if (length(term_parts)) {
            dplyr::bind_rows(term_parts) %>%
              sanitize_bridge() %>%
              dplyr::filter(!is.na(term), nzchar(term)) %>%
              dplyr::distinct()
          } else {
            tibble::tibble()
          }
          
          append_log(paste0(log_tag, " term_bridge nrow=", nrow(term_bridge)))
          append_log(
            paste0(
              log_tag, " term_bridge head:\n",
              paste(capture.output(print(utils::head(term_bridge, 5))), collapse = "\n")
            )
          )
          
          if (nrow(term_bridge) > 0) {
            saveRDS(term_bridge, term_bridge_path)
            manifest$files[[paste0(app_slug, "_term_bridge")]] <- basename(term_bridge_path)
            append_log(paste0(log_tag, " saved ", basename(term_bridge_path), " | exists=", file.exists(term_bridge_path)))
          } else {
            append_log(paste0(log_tag, " ", basename(term_bridge_path), " NOT saved (0 rows)"))
          }
          
          # ============================================================
          # C) LD INPUTS
          # ============================================================
          append_log(paste0(ld_tag, " dir: ", session_dir))
          append_log(paste0(ld_tag, " dir exists: ", dir.exists(session_dir)))
          
          # ----------------------------------------------------------
          # C1) CLUSTERS MASTER
          # ----------------------------------------------------------
          append_log(paste0(ld_tag, "[CLUSTERS] cl2 nrow=", nrow(cl2), " ncol=", ncol(cl2)))
          append_log(paste0(ld_tag, "[CLUSTERS] cl2 cols: ", paste(names(cl2), collapse = ", ")))
          append_log(
            paste0(
              ld_tag, "[CLUSTERS] cl2 head:\n",
              paste(capture.output(print(utils::head(cl2, 5))), collapse = "\n")
            )
          )
          
          chr_col <- intersect(c("chr","CHR","chrom","CHROM","chromosome"), names(cl2))
          st_col  <- intersect(c("start","START","cluster_start","start_bp","FROM","from","bp1"), names(cl2))
          en_col  <- intersect(c("end","END","cluster_end","end_bp","TO","to","bp2"), names(cl2))
          id_col  <- intersect(c("cluster_id","CLUSTER_ID","cluster","id"), names(cl2))
          key_col <- intersect(c("cluster_key","CLUSTER_KEY"), names(cl2))
          
          chr_col <- if (length(chr_col)) chr_col[1] else NULL
          st_col  <- if (length(st_col))  st_col[1]  else NULL
          en_col  <- if (length(en_col))  en_col[1]  else NULL
          id_col  <- if (length(id_col))  id_col[1]  else NULL
          key_col <- if (length(key_col)) key_col[1] else NULL
          
          append_log(paste0(ld_tag, "[CLUSTERS] chr_col=", chr_col %||% "NULL"))
          append_log(paste0(ld_tag, "[CLUSTERS] st_col=",  st_col  %||% "NULL"))
          append_log(paste0(ld_tag, "[CLUSTERS] en_col=",  en_col  %||% "NULL"))
          append_log(paste0(ld_tag, "[CLUSTERS] id_col=",  id_col  %||% "NULL"))
          append_log(paste0(ld_tag, "[CLUSTERS] key_col=", key_col %||% "NULL"))
          
          clusters_master <- tibble::tibble()
          
          if (is.null(chr_col) || is.null(st_col) || is.null(en_col) || is.null(id_col)) {
            append_log(paste0(ld_tag, "[CLUSTERS] ", basename(clusters_master_path), " NOT saved (required columns missing)"))
          } else {
            clusters_master <- cl2 %>%
              dplyr::transmute(
                cluster_id  = as.character(.data[[id_col]]),
                chr         = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[chr_col]])))),
                start_raw   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
                end_raw     = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]])))),
                cluster_key = if (!is.null(key_col)) as.character(.data[[key_col]]) else NA_character_
              ) %>%
              dplyr::mutate(
                cluster_id = trimws(cluster_id),
                start = pmin(start_raw, end_raw),
                end   = pmax(start_raw, end_raw),
                cluster_key = dplyr::if_else(
                  !is.na(cluster_key) & nzchar(cluster_key),
                  cluster_key,
                  paste0(cluster_id, "|chr", chr, ":", start, "-", end)
                )
              ) %>%
              dplyr::select(cluster_id, chr, start, end, cluster_key) %>%
              dplyr::filter(
                !is.na(cluster_id), nzchar(cluster_id),
                is.finite(chr), is.finite(start), is.finite(end)
              ) %>%
              dplyr::distinct()
            
            append_log(paste0(ld_tag, "[CLUSTERS] clusters_master nrow=", nrow(clusters_master), " ncol=", ncol(clusters_master)))
            append_log(paste0(ld_tag, "[CLUSTERS] clusters_master cols: ", paste(names(clusters_master), collapse = ", ")))
            append_log(
              paste0(
                ld_tag, "[CLUSTERS] clusters_master head:\n",
                paste(capture.output(print(utils::head(clusters_master, 5))), collapse = "\n")
              )
            )
            
            if (nrow(clusters_master) > 0) {
              saveRDS(clusters_master, clusters_master_path)
              manifest$files[[paste0(app_slug, "_clusters_master")]] <- basename(clusters_master_path)
              append_log(paste0(ld_tag, "[CLUSTERS] saved: ", clusters_master_path))
              append_log(paste0(ld_tag, "[CLUSTERS] exists after save: ", file.exists(clusters_master_path)))
            } else {
              append_log(paste0(ld_tag, "[CLUSTERS] ", basename(clusters_master_path), " NOT saved (0 rows)"))
            }
          }
          
          # ----------------------------------------------------------
          # C2) CANDIDATES FOR LD = GWAS + APP HIT
          # ----------------------------------------------------------
          get_gwas_hits_for_ld <- function() {
            if (exists("hits_df", mode = "function", inherits = TRUE)) {
              x <- tryCatch(hits_df(), error = function(e) NULL)
              if (is.data.frame(x) && nrow(x) > 0) {
                return(list(name = "hits_df()", data = x))
              }
            }
            
            list(name = NULL, data = NULL)
          }
          
          append_log(paste0(ld_tag, "[CANDIDATES] hits nrow=", nrow(hits), " ncol=", ncol(hits)))
          append_log(paste0(ld_tag, "[CANDIDATES] hits cols: ", paste(names(hits), collapse = ", ")))
          
          gwas_src   <- get_gwas_hits_for_ld()
          gwas_df_ld <- gwas_src$data
          
          append_log(paste0(ld_tag, "[CANDIDATES][GWAS] filtered source picked: ", gwas_src$name %||% "NULL"))
          
          gwas_candidates <- tibble::tibble()
          
          if (is.data.frame(gwas_df_ld) && nrow(gwas_df_ld)) {
            
            # hits_df() té: CHR, BP, snp, p, logp
            gwas_candidates0 <- gwas_df_ld %>%
              dplyr::transmute(
                cluster_id = NA_character_,
                chr        = suppressWarnings(as.integer(readr::parse_number(as.character(CHR)))),
                pos_ini    = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
                pos_end    = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
                id_hit     = as.character(snp),
                classe     = "GWAS"
              ) %>%
              dplyr::mutate(
                cluster_id = trimws(cluster_id),
                id_hit     = trimws(id_hit)
              ) %>%
              dplyr::filter(
                is.finite(chr), is.finite(pos_ini), is.finite(pos_end),
                !is.na(id_hit), nzchar(id_hit)
              ) %>%
              dplyr::distinct()
            
            gwas_candidates <- gwas_candidates0 %>%
              dplyr::left_join(
                clusters_master %>%
                  dplyr::rename(
                    cluster_chr   = chr,
                    cluster_start = start,
                    cluster_end   = end
                  ),
                by = dplyr::join_by(
                  chr == cluster_chr,
                  pos_ini >= cluster_start,
                  pos_ini <= cluster_end
                )
              ) %>%
              dplyr::mutate(
                cluster_id = cluster_id.y
              ) %>%
              dplyr::select(cluster_id, chr, pos_ini, pos_end, id_hit, classe) %>%
              dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
              dplyr::distinct()
          }
          
          append_log(paste0(ld_tag, "[CANDIDATES][GWAS] nrow=", nrow(gwas_candidates)))
          if (nrow(gwas_candidates) > 0) {
            append_log(
              paste0(
                ld_tag, "[CANDIDATES][GWAS] head:\n",
                paste(capture.output(print(utils::head(gwas_candidates, 10))), collapse = "\n")
              )
            )
          }
          
          app_candidates <- tibble::tibble()
          
          needed_hit_cols <- c(col_cluster_id, col_chr, col_pos, col_id)
          if (all(needed_hit_cols %in% names(hits))) {
            app_candidates <- hits %>%
              dplyr::transmute(
                cluster_id = as.character(.data[[col_cluster_id]]),
                chr        = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[col_chr]])))),
                pos_ini    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[col_pos]])))),
                pos_end    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[col_pos]])))),
                id_hit     = as.character(.data[[col_id]]),
                classe     = app_hit_class
              ) %>%
              dplyr::mutate(
                cluster_id = trimws(cluster_id),
                id_hit     = trimws(id_hit)
              ) %>%
              dplyr::filter(
                !is.na(cluster_id), nzchar(cluster_id),
                is.finite(chr), is.finite(pos_ini), is.finite(pos_end),
                !is.na(id_hit), nzchar(id_hit)
              ) %>%
              dplyr::distinct()
          } else {
            append_log(paste0(
              ld_tag, "[CANDIDATES][APP] skipped: missing columns -> ",
              paste(setdiff(needed_hit_cols, names(hits)), collapse = ", ")
            ))
          }
          
          append_log(paste0(ld_tag, "[CANDIDATES][APP] nrow=", nrow(app_candidates)))
          if (nrow(app_candidates) > 0) {
            append_log(
              paste0(
                ld_tag, "[CANDIDATES][APP] head:\n",
                paste(capture.output(print(utils::head(app_candidates, 10))), collapse = "\n")
              )
            )
          }
          
          candidates_ld <- dplyr::bind_rows(gwas_candidates, app_candidates) %>%
            dplyr::mutate(
              rsid     = id_hit,
              position = pos_ini
            ) %>%
            dplyr::select(cluster_id, chr, pos_ini, pos_end, id_hit, rsid, position, classe) %>%
            dplyr::filter(
              !is.na(cluster_id), nzchar(cluster_id),
              is.finite(chr), is.finite(pos_ini), is.finite(pos_end),
              !is.na(id_hit), nzchar(id_hit)
            ) %>%
            dplyr::distinct() %>%
            dplyr::arrange(chr, pos_ini, id_hit, classe)
          
          append_log(paste0(ld_tag, "[CANDIDATES] combined nrow=", nrow(candidates_ld), " ncol=", ncol(candidates_ld)))
          append_log(paste0(ld_tag, "[CANDIDATES] combined cols: ", paste(names(candidates_ld), collapse = ", ")))
          append_log(
            paste0(
              ld_tag, "[CANDIDATES] combined head:\n",
              paste(capture.output(print(utils::head(candidates_ld, 10))), collapse = "\n")
            )
          )
          
          if (nrow(candidates_ld) > 0) {
            append_log(paste0(
              ld_tag, "[CANDIDATES] classe counts:\n",
              paste(capture.output(print(table(candidates_ld$classe, useNA = "ifany"))), collapse = "\n")
            ))
          }
          
          if (nrow(candidates_ld) > 0) {
            saveRDS(candidates_ld, candidates_path)
            manifest$files[[paste0(app_slug, "_candidates")]] <- basename(candidates_path)
            append_log(paste0(ld_tag, "[CANDIDATES] saved: ", candidates_path))
            append_log(paste0(ld_tag, "[CANDIDATES] exists after save: ", file.exists(candidates_path)))
          } else {
            append_log(paste0(ld_tag, "[CANDIDATES] ", basename(candidates_path), " NOT saved (0 rows)"))
          }
          
          # ----------------------------------------------------------
          # C3) GWAS SIGNIFICANCE BRIDGE (ONLY FILTERED GWAS HITS)
          # ----------------------------------------------------------
          append_log(paste0(log_tag, " [GWAS SIGNIFICANCE BRIDGE] hits_df source requested"))
          
          gwas_hits_for_bridge <- tryCatch(
            hits_df(),
            error = function(e) {
              append_log(paste0(log_tag, " [GWAS SIGNIFICANCE BRIDGE][ERROR hits_df] ", conditionMessage(e)))
              tibble::tibble()
            }
          )
          
          append_log(paste0(log_tag, " [GWAS SIGNIFICANCE BRIDGE] hits_df nrow=", nrow(gwas_hits_for_bridge)))
          append_log(paste0(log_tag, " [GWAS SIGNIFICANCE BRIDGE] clusters_master nrow=", nrow(clusters_master)))
          
          gwas_sig_bridge <- tryCatch({
            
            if (!is.data.frame(gwas_hits_for_bridge) || !nrow(gwas_hits_for_bridge)) {
              tibble::tibble(
                cluster_id = character(),
                chr = integer(),
                position = integer(),
                rsid = character(),
                p_value = numeric(),
                logp = numeric()
              )
              
            } else {
              
              gwas_hits_std <- gwas_hits_for_bridge %>%
                dplyr::transmute(
                  chr      = suppressWarnings(as.integer(readr::parse_number(as.character(CHR)))),
                  position = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
                  rsid     = as.character(snp),
                  p_value  = suppressWarnings(as.numeric(p)),
                  logp     = suppressWarnings(as.numeric(logp))
                ) %>%
                dplyr::mutate(
                  rsid = trimws(rsid)
                ) %>%
                dplyr::filter(
                  is.finite(chr),
                  is.finite(position),
                  !is.na(rsid), nzchar(rsid)
                ) %>%
                dplyr::distinct()
              
              gwas_hits_std %>%
                dplyr::left_join(
                  clusters_master %>%
                    dplyr::rename(
                      cluster_chr   = chr,
                      cluster_start = start,
                      cluster_end   = end
                    ),
                  by = dplyr::join_by(
                    chr == cluster_chr,
                    position >= cluster_start,
                    position <= cluster_end
                  )
                ) %>%
                dplyr::transmute(
                  cluster_id = cluster_id,
                  chr = chr,
                  position = position,
                  rsid = rsid,
                  p_value = p_value,
                  logp = logp
                ) %>%
                dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
                dplyr::distinct()
            }
            
          }, error = function(e) {
            append_log(paste0(log_tag, " [GWAS SIGNIFICANCE BRIDGE][ERROR] ", conditionMessage(e)))
            tibble::tibble(
              cluster_id = character(),
              chr = integer(),
              position = integer(),
              rsid = character(),
              p_value = numeric(),
              logp = numeric()
            )
          })
          
          saveRDS(gwas_sig_bridge, gwas_sig_bridge_path)
          manifest$files[["gwas_significance_bridge"]] <- basename(gwas_sig_bridge_path)
          
          append_log(paste0(log_tag, " gwas_sig_bridge_path: ", gwas_sig_bridge_path))
          append_log(paste0(log_tag, " gwas_sig_bridge nrow=", nrow(gwas_sig_bridge)))
          
          if (nrow(gwas_sig_bridge) > 0) {
            append_log(
              paste0(
                log_tag, " gwas_sig_bridge head:\n",
                paste(capture.output(print(utils::head(gwas_sig_bridge, 10))), collapse = "\n")
              )
            )
          }
          
          cat(
            "[GWAS_SIGNIFICANCE_BRIDGE] saved: ",
            gwas_sig_bridge_path,
            " | nrow=", nrow(gwas_sig_bridge), "\n",
            sep = ""
          )
          
          # ============================================================
          # FINAL MANIFEST WRITE
          # ============================================================
          manifest$apps_present  <- sort(unique(c(manifest$apps_present, app_slug)))
          manifest$last_updated  <- as.character(Sys.time())
          
          saveRDS(manifest, manifest_path)
          
          append_log(paste0(log_tag, " manifest saved: ", manifest_path))
          append_log(
            paste0(
              log_tag, " manifest content:\n",
              paste(capture.output(str(manifest, max.level = 2)), collapse = "\n")
            )
          )
          
        }, error = function(e) {
          if (exists("append_log", mode = "function", inherits = TRUE)) {
            append_log(paste0("[LD-INTEGRATOR][ERROR] ", conditionMessage(e)))
          } else {
            message("[LD-INTEGRATOR][ERROR] ", conditionMessage(e))
          }
        })
        
        ####################### end integrator bridges #################################
        
        
        incProgress(0.05, detail = "Done")
        append_catalog_log("[END] run_catalog finished OK")
      })
      
      if (is.null(hits) || !nrow(hits)) {
        showNotification("No Catalog hits found in clusters.", type = "warning", duration = 6)
      } else {
        showNotification("Catalog hits assigned to clusters (n_catalog updated).", type = "message", duration = 4)
      }
      
    }, error = function(e) {
      
      catalog_final_path_csv(NULL)
      catalog_final_path_rds(NULL)
      
      cl_now <- clusters_val()
      if (!is.null(cl_now) && nrow(cl_now)) {
        cl_now$n_catalog <- 0L
        clusters_val(cl_now)
      }
      
      msg <- paste0("ERROR (run_catalog): ", conditionMessage(e))
      output$catalog_summary <- renderText(msg)
      
      append_catalog_log("[ERROR] ", msg)
      showNotification(conditionMessage(e), type = "error", duration = 12)
      
    }, finally = {
      
      if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::enable("run_catalog")
      
    })
    
  }, ignoreInit = TRUE)
  
  
  output$download_catalog <- downloadHandler(
    filename = function() {
      mode_tag <- if (identical(input$cluster_method, "window")) "w" else "h"
      thr_val  <- if (identical(input$cluster_method, "window")) input$pthr else input$min_logp
      thr_tag  <- gsub("\\.", "p", sprintf("%.2f", thr_val))
      paste0("catalog_hits_by_cluster_", format(Sys.Date(), "%Y%m%d"), "_", mode_tag, "_thr", thr_tag, ".csv")
    },
    content = function(file) {
      dt <- tryCatch(catalog_data(), error = function(e) NULL)
      if (is.null(dt) || !nrow(dt)) {
        writeLines("No Catalog hits extracted.", con = file)
        return()
      }
      
      # Ensure cluster_start/end exist (join from clusters_val if needed)
      cl <- tryCatch(clusters_val(), error = function(e) NULL)
      if (!"cluster_start" %in% names(dt) || !"cluster_end" %in% names(dt)) {
        if (is.data.frame(cl) && nrow(cl) && all(c("cluster_id","start","end") %in% names(cl))) {
          dt <- dt %>%
            dplyr::left_join(
              cl %>% dplyr::transmute(
                cluster_id,
                cluster_start = as.integer(start),
                cluster_end   = as.integer(end)
              ),
              by = "cluster_id"
            )
        } else {
          if (!"cluster_start" %in% names(dt)) dt$cluster_start <- NA_integer_
          if (!"cluster_end"   %in% names(dt)) dt$cluster_end   <- NA_integer_
        }
      }
      
      # Write CSV (not TSV)
      readr::write_csv(dt, file, na = "")
    }
  )
  
  output$download_catalog_rds <- downloadHandler(
    filename = function() {
      mode_tag <- if (identical(input$cluster_method, "window")) "w" else "h"
      thr_val  <- if (identical(input$cluster_method, "window")) input$pthr else input$min_logp
      thr_tag  <- gsub("\\.", "p", sprintf("%.2f", thr_val))
      paste0("catalog_hits_by_cluster_", format(Sys.Date(), "%Y%m%d"), "_", mode_tag, "_thr", thr_tag, ".rds")
    },
    content = function(file) {
      dt <- tryCatch(catalog_data(), error = function(e) NULL)
      validate(shiny::need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits extracted."))
      
      cl <- tryCatch(clusters_val(), error = function(e) NULL)
      if (!"cluster_start" %in% names(dt) || !"cluster_end" %in% names(dt)) {
        if (is.data.frame(cl) && nrow(cl) && all(c("cluster_id","start","end") %in% names(cl))) {
          dt <- dt %>%
            dplyr::left_join(
              cl %>% dplyr::transmute(
                cluster_id,
                cluster_start = as.integer(start),
                cluster_end   = as.integer(end)
              ),
              by = "cluster_id"
            )
        } else {
          if (!"cluster_start" %in% names(dt)) dt$cluster_start <- NA_integer_
          if (!"cluster_end"   %in% names(dt)) dt$cluster_end   <- NA_integer_
        }
      }
      
      saveRDS(dt, file)
    }
  )
  
  # ---------------------------------------------------------------------------
  # Cluster summary DT
  # ---------------------------------------------------------------------------
  
  clusters_for_ui <- reactive({
    cl <- tryCatch(clusters_val(), error = function(e) NULL)
    if (!is.data.frame(cl) || !nrow(cl)) return(cl)
    
    cl <- as.data.frame(cl)
    
    # aliases per tolerar esquemes (canonical / legacy)
    if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- suppressWarnings(as.integer(cl$CHR))
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- suppressWarnings(as.integer(cl$start_bp))
    if (!"end" %in% names(cl) && "end_bp" %in% names(cl)) cl$end <- suppressWarnings(as.integer(cl$end_bp))
    
    # cluster_id
    if (!"cluster_id" %in% names(cl) && "cluster_chr_n" %in% names(cl)) cl$cluster_id <- as.character(cl$cluster_chr_n)
    
    # ✅ aquest és el punt clau: cluster_n per al Catalog
    if (!"cluster_n" %in% names(cl)) {
      if ("cluster_chr" %in% names(cl)) cl$cluster_n <- suppressWarnings(as.integer(cl$cluster_chr))
      else if ("cluster" %in% names(cl)) cl$cluster_n <- suppressWarnings(as.integer(cl$cluster))
      else cl$cluster_n <- seq_len(nrow(cl))
    }
    
    # placeholders (perquè el DT no peti si encara no hi ha Step4)
    if (!"n_snps" %in% names(cl)) cl$n_snps <- NA_integer_
    if (!"top_snp" %in% names(cl)) cl$top_snp <- NA_character_
    if (!"top_logp" %in% names(cl)) cl$top_logp <- NA_real_
    if (!"cluster_size_kb" %in% names(cl)) cl$cluster_size_kb <- round((as.integer(cl$end) - as.integer(cl$start)) / 1000, 2)
    if (!"n_catalog" %in% names(cl)) cl$n_catalog <- 0L
    
    cl
  })
  
  
  output$cluster_dt <- renderDT({
    df <- clusters_for_ui()
    if (is.null(df) || !nrow(df)) {
      return(datatable(
        data.frame(Message="No clusters yet. NEW mode: click ➊ Generate intervals → merge → clusters (or OLD mode: load precomputed file)."),
        options = list(dom="t"),
        rownames = FALSE
      ))
    }
    
    df_show <- df %>%
      mutate(chr_label = chr_label_plink(as.integer(chr))) %>%
      select(
        cluster_chr_n = cluster_id,
        chr = chr_label,
        cluster = cluster_n,
        start, end,
        n_snps, top_snp, top_logp,
        cluster_size_kb,
        n_catalog
      )
    
    datatable(
      df_show,
      selection  = "single",
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX    = TRUE
      )
    ) %>%
      formatRound(columns = c("top_logp", "cluster_size_kb"), digits = 2)
  }, server = FALSE)
  
  
  selected_cluster <- reactive({
    idx <- input$cluster_dt_rows_selected
    df  <- clusters_for_ui()
    if (is.null(idx) || !length(idx) || is.null(df) || !nrow(df)) return(NULL)
    df[idx[1], , drop=FALSE]
  })
  
  
  selected_cluster <- reactive({
    idx <- input$cluster_dt_rows_selected
    df  <- clusters_val()
    if (is.null(idx) || !length(idx) || is.null(df) || !nrow(df)) return(NULL)
    df[idx, , drop=FALSE]
  })
  
  # ---------------------------------------------------------------------------
  # Catalog table + scope (Selected / Pick) + chr/cluster picker (ONLY clusters with hits)
  # + gene selector depending on the current cluster-filtered subset
  # ---------------------------------------------------------------------------
  
  # --- helper: parse chr from cluster_id (fallback) ---
  clusterid_to_chr <- function(x) {
    x <- as.character(x)
    m <- regmatches(x, regexpr("chr[^_]+", x, ignore.case = TRUE))
    m <- gsub("^chr", "", m, ignore.case = TRUE)
    m <- toupper(m)
    m[m == "X"]  <- "23"
    m[m == "Y"]  <- "24"
    m[m %in% c("M","MT")] <- "25"
    suppressWarnings(as.integer(m))
  }
  
  # --- helper: does a row-gene string contain any selected genes? ---
  row_has_any_gene <- function(gene_str, genes_sel) {
    if (is.null(genes_sel) || !length(genes_sel)) return(TRUE)
    if (is.na(gene_str) || !nzchar(gene_str)) return(FALSE)
    gs <- split_multi_unique(gene_str)
    any(gs %in% genes_sel)
  }
  
  # --- clusters with hits (from catalog_data), with chr (from clusters_val() if possible) ---
  catalog_clusters_with_hits <- reactive({
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits loaded."))
    validate(need("cluster_id" %in% names(dt), "Catalog hits missing cluster_id."))
    
    tab <- dt %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::count(cluster_id, sort = TRUE)
    
    validate(need(nrow(tab) > 0, "No clustered hits found (cluster_id empty/NA)."))
    
    cl <- tryCatch(clusters_val(), error = function(e) NULL)
    if (is.data.frame(cl) && nrow(cl) > 0 && all(c("cluster_id","chr") %in% names(cl))) {
      tab <- tab %>%
        dplyr::left_join(cl %>% dplyr::select(cluster_id, chr), by = "cluster_id") %>%
        dplyr::mutate(chr = suppressWarnings(as.integer(chr)))
    } else {
      tab <- tab %>% dplyr::mutate(chr = clusterid_to_chr(cluster_id))
    }
    
    tab %>%
      dplyr::filter(is.finite(chr)) %>%
      dplyr::arrange(chr, dplyr::desc(n))
  })
  
  # --- UI: pick a cluster here (chr -> cluster), only clusters with hits ---
  output$catalog_cluster_pick_ui <- renderUI({
    tab <- catalog_clusters_with_hits()
    validate(need(is.data.frame(tab) && nrow(tab) > 0, "No clusters with hits to pick."))
    
    chr_choices <- sort(unique(as.integer(tab$chr)))
    chr_labels  <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("catalog_pick_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("catalog_cluster_pick_id_ui")
    )
  })
  
  output$catalog_cluster_pick_id_ui <- renderUI({
    tab <- catalog_clusters_with_hits()
    req(input$catalog_pick_chr)
    
    x <- tab %>%
      dplyr::filter(as.integer(chr) == as.integer(input$catalog_pick_chr)) %>%
      dplyr::arrange(dplyr::desc(n), cluster_id)
    
    validate(need(nrow(x) > 0, "No clusters with hits for this chromosome."))
    
    choices <- stats::setNames(x$cluster_id, paste0(x$cluster_id, " (", x$n, " hits)"))
    
    selectInput(
      "catalog_pick_cluster_id", "Cluster (with hits)",
      choices  = choices,
      selected = x$cluster_id[1]
    )
  })
  
  # --- compute the cluster-filtered subset (NO gene filter here) ---
  catalog_cluster_subset <- reactive({
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits loaded."))
    validate(need("cluster_id" %in% names(dt), "Catalog hits missing cluster_id."))
    
    scope <- input$catalog_scope %||% "selected"
    
    if (identical(scope, "selected")) {
      cl <- selected_cluster()
      validate(need(is.data.frame(cl) && nrow(cl) > 0, "No cluster selected in Clusters table."))
      dt %>% dplyr::filter(cluster_id == cl$cluster_id[1])
      
    } else if (identical(scope, "pick")) {
      req(input$catalog_pick_cluster_id)
      dt %>% dplyr::filter(cluster_id == as.character(input$catalog_pick_cluster_id))
      
    } else {
      dt
    }
  })
  
  # --- UI: gene selector depends on current cluster-filtered subset ---
  output$catalog_gene_pick_ui <- renderUI({
    dt_sub <- catalog_cluster_subset()
    validate(need(is.data.frame(dt_sub) && nrow(dt_sub) > 0, "No hits in the chosen cluster."))
    
    gene_col <- pick_col(dt_sub, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    validate(need(!is.null(gene_col), "No gene column found (MAPPED_GENE/GENE)."))
    
    genes <- split_multi_unique(dt_sub[[gene_col]])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- sort(unique(genes))
    
    if (!length(genes)) return(NULL)
    
    selectizeInput(
      "catalog_gene_pick", "Gene(s)",
      choices  = genes,
      selected = NULL,
      multiple = TRUE,
      options  = list(placeholder = "Optional: filter by gene(s) within this cluster…")
    )
  })
  
  # --- final filtered table: cluster subset + optional gene filter ---
  catalog_table_filtered <- reactive({
    dt <- catalog_cluster_subset()
    if (!is.data.frame(dt) || !nrow(dt)) return(dt)
    
    gene_col  <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    genes_sel <- input$catalog_gene_pick
    
    if (!is.null(gene_col) && !is.null(genes_sel) && length(genes_sel) > 0) {
      keep <- vapply(dt[[gene_col]], row_has_any_gene, logical(1), genes_sel = genes_sel)
      dt <- dt[keep, , drop = FALSE]
    }
    print(head(dt))
    dt
  })
  
  # --- renderDT: keep ONLY requested columns + links for dbSNP / GeneCards / trait / PubMed ---
  output$catalog_table <- DT::renderDT({
    dt <- catalog_table_filtered()
    
    if (is.null(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message="No Catalog hits in this selection."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }
    
    tryCatch({
      chr_col   <- pick_col(dt, c("CHR","CHR_std","chr"))
      pos_col   <- pick_col(dt, c("POS","POS_std","BP","position","POSITION"))
      p_col     <- pick_col(dt, c("P_VALUE","p_value","PVAL","P","p","P-VALUE"))
      or_col    <- pick_col(dt, c("OR_BETA","OR","BETA","beta","ODDSRATIO","odds_ratio"))
      gene_col  <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
      dis_col   <- pick_col(dt, c("DISEASE","Disease","disease"))
      trait_col <- pick_col(dt, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
      rs_col    <- pick_col(dt, c("rsid","RSID","SNPS","SNP","snp","variant_id","VARIANT_ID"))
      pm_col    <- pick_col(dt, c("PUBMEDID","PMID","PubMed","PUBMED_ID","PUBMED"))
      
      disp <- dt %>%
        dplyr::transmute(
          cluster_chr_n = dplyr::coalesce(as.character(cluster_id), NA_character_),
          
          CHR = if (!is.null(chr_col)) as.character(.data[[chr_col]]) else NA_character_,
          POS = if (!is.null(pos_col)) suppressWarnings(readr::parse_number(as.character(.data[[pos_col]]))) else NA_real_,
          
          rsid_std = if (!is.null(rs_col)) as.character(.data[[rs_col]]) else NA_character_,
          P_std    = if (!is.null(p_col)) parse_p_robust(.data[[p_col]]) else NA_real_,
          OR_std   = if (!is.null(or_col)) suppressWarnings(readr::parse_number(as.character(.data[[or_col]]))) else NA_real_,
          
          gene_std  = if (!is.null(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
          dis_std   = if (!is.null(dis_col)) as.character(.data[[dis_col]]) else NA_character_,
          trait_std = if (!is.null(trait_col)) as.character(.data[[trait_col]]) else NA_character_,
          
          pmid_std  = if (!is.null(pm_col)) as.character(.data[[pm_col]]) else NA_character_
        ) %>%
        dplyr::mutate(
          # links
          dbSNP = ifelse(!is.na(rsid_std) & nzchar(rsid_std),
                         paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/", rsid_std,
                                "' target='_blank'>", rsid_std, "</a>"),
                         NA_character_),
          
          MAPPED_GENE = ifelse(!is.na(gene_std) & nzchar(gene_std),
                               paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                                      URLencode(gene_std, reserved = TRUE),
                                      "' target='_blank'>", gene_std, "</a>"),
                               NA_character_),
          
          MAPPED_TRAIT = ifelse(!is.na(trait_std) & nzchar(trait_std),
                                paste0("<a href='https://www.ebi.ac.uk/gwas/search?query=",
                                       URLencode(trait_std, reserved = TRUE),
                                       "' target='_blank'>", trait_std, "</a>"),
                                NA_character_),
          
          PubMed = ifelse(!is.na(pmid_std) & nzchar(pmid_std),
                          paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/",
                                 pmid_std, "/' target='_blank'>", pmid_std, "</a>"),
                          NA_character_),
          
          # formatting
          P = ifelse(!is.na(P_std), formatC(P_std, format = "e", digits = 2), NA_character_),
          OR_BETA = ifelse(!is.na(OR_std), sprintf("%.3f", OR_std), NA_character_)
        )
      
      show_df <- disp %>%
        dplyr::transmute(
          cluster_chr_n,
          CHR,
          POS,
          dbSNP,
          P,
          OR_BETA,
          MAPPED_GENE,
          DISEASE = dis_std,
          MAPPED_TRAIT,
          PubMed
        )
      
      DT::datatable(
        show_df,
        escape = FALSE,
        rownames = FALSE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons = list(
            list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
            list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
            list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
          ),
          pageLength = 10,
          scrollX    = TRUE
        )
      )
      
    }, error = function(e) {
      DT::datatable(
        data.frame(ERROR = paste("Catalog table render error:", conditionMessage(e))),
        options = list(dom = "t"), rownames = FALSE
      )
    })
  }, server = FALSE)
  
  # ---------------------------------------------------------------------------
  # Manhattan_combo (top GWAS + bottom Catalog)  + UCSC region via zoom (relayout)
  # ---------------------------------------------------------------------------
  
  #  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # -------- state: region selected from manhattan_combo zoom (BPcum range) -----
  combo_region <- reactiveVal(NULL)   # list(x0=, x1=) in BPcum units
  
  # Helper: robustly read x-range from plotly_relayout (xaxis or xaxis2)
  get_relayout_xrange <- function(d) {
    if (is.null(d) || !length(d)) return(NULL)
    
    # reset zoom: autorange true
    if (isTRUE(d[["xaxis.autorange"]]) || isTRUE(d[["xaxis2.autorange"]])) return(NULL)
    
    x0 <- d[["xaxis.range[0]"]] %||% NA
    x1 <- d[["xaxis.range[1]"]] %||% NA
    
    # subplot can send xaxis2
    if (!is.finite(suppressWarnings(as.numeric(x0))) || !is.finite(suppressWarnings(as.numeric(x1)))) {
      x0 <- d[["xaxis2.range[0]"]] %||% NA
      x1 <- d[["xaxis2.range[1]"]] %||% NA
    }
    
    x0 <- suppressWarnings(as.numeric(x0))
    x1 <- suppressWarnings(as.numeric(x1))
    if (!is.finite(x0) || !is.finite(x1)) return(NULL)
    
    if (x1 < x0) { tmp <- x0; x0 <- x1; x1 <- tmp }
    c(x0, x1)
  }
  
  # Helper: BPcum range -> UCSC region string using .ref_hg38
  # expects ref has: chr (numeric or character), chr_cum, len
  cumrange_to_ucsc_region <- function(x0, x1, ref) {
    if (!is.finite(x0) || !is.finite(x1)) return(NULL)
    if (x1 < x0) { tmp <- x0; x0 <- x1; x1 <- tmp }
    
    ref2 <- ref %>%
      dplyr::transmute(
        chr = suppressWarnings(as.integer(chr)),
        chr_cum = suppressWarnings(as.numeric(chr_cum)),
        len = suppressWarnings(as.numeric(len))
      ) %>%
      dplyr::filter(is.finite(chr), is.finite(chr_cum), is.finite(len), len > 0) %>%
      dplyr::arrange(chr)
    
    if (!nrow(ref2)) return(NULL)
    
    mid <- (x0 + x1) / 2
    idx <- which(mid >= ref2$chr_cum & mid <= (ref2$chr_cum + ref2$len))
    if (!length(idx)) {
      # fallback: nearest chromosome by distance to [chr_cum, chr_cum+len]
      dist <- pmin(abs(mid - ref2$chr_cum), abs(mid - (ref2$chr_cum + ref2$len)))
      idx <- which.min(dist)
    } else {
      idx <- idx[1]
    }
    
    chr_num <- ref2$chr[idx]
    chr_c0  <- ref2$chr_cum[idx]
    chr_end <- chr_c0 + ref2$len[idx]
    
    # clamp x0/x1 to this chromosome to avoid cross-chr weirdness
    xx0 <- max(x0, chr_c0)
    xx1 <- min(x1, chr_end)
    
    bp0 <- as.integer(floor(xx0 - chr_c0))
    bp1 <- as.integer(ceiling(xx1 - chr_c0))
    bp0 <- max(1L, bp0)
    bp1 <- max(bp0 + 1L, bp1)
    
    # label chr using your helper chr_label_plink()
    chr_lab <- chr_label_plink(as.integer(chr_num))
    paste0("chr", chr_lab, ":", bp0, "-", bp1)
  }
  
  # ---------------------------------------------------------------------------
  # Manhattan_combo plotly
  # ---------------------------------------------------------------------------
  
  
  output$manhattan_combo <- renderPlotly({
    
    src_combo <- "manhattan_combo"
    
    dfp    <- tryCatch(dfp_manhattan(), error = function(e) NULL)
    cat_dt <- tryCatch(catalog_data(),  error = function(e) NULL)
    
    if (is.null(dfp) || !is.data.frame(dfp) || !nrow(dfp)) {
      return(plotly_message("⚠️ GWAS table missing or incomplete."))
    }
    
    ref <- .ref_hg38
    ax  <- axis_df()
    axis_breaks <- ax$center
    axis_labels <- paste0("chr", ax$chrN)
    GENOME_END  <- max(ref$chr_cum + ref$len)
    
    # --- ref_map: chr_num -> chr_cum (robust) ---
    chr_key <- if ("chrN" %in% names(ref)) "chrN" else if ("chr" %in% names(ref)) "chr" else if ("CHR" %in% names(ref)) "CHR" else NA_character_
    ref_map <- NULL
    if (!is.na(chr_key) && "chr_cum" %in% names(ref)) {
      ref_map <- ref %>%
        dplyr::transmute(
          chr_num = suppressWarnings(as.integer(.data[[chr_key]])),
          chr_cum = as.numeric(.data$chr_cum)
        ) %>%
        dplyr::filter(is.finite(chr_num), is.finite(chr_cum))
    }
    
    thr_y <- if ((input$cluster_method %||% "window") == "window") (input$pthr %||% 7.3) else (input$min_logp %||% 6)
    
    # ============================
    # p1: GWAS (NO segments here)
    # ============================
    dfp <- dfp %>%
      dplyr::arrange(CHR, BP) %>%
      dplyr::mutate(
        rs_show = dplyr::if_else(!is.na(snp) & nzchar(snp), snp, paste0("chr", CHR, ":", BP)),
        chr_lab = chr_label_plink(as.integer(CHR)),
        col     = ifelse((as.integer(CHR) %% 2) == 0, "darkgreen", "#ff7f00")
      )
    
    dfp$tooltip <- paste0(
      "<b>GWAS hit</b>",
      "<br><b>rs:</b> ", dfp$rs_show,
      "<br><b>CHR:</b> ", dfp$chr_lab,
      "<br><b>BP:</b> ", dfp$BP,
      "<br><b>P:</b> ", signif(dfp$Pval, 3),
      "<br><b>-log10(P):</b> ", round(dfp$logp, 2)
    )
    
    p1 <- ggplot2::ggplot(dfp, ggplot2::aes(x = BPcum, y = logp, text = tooltip)) +
      ggplot2::geom_point(ggplot2::aes(color = col), size = 1)
    
    # threshold line only in "new" mode
    p1 <- p1 + ggplot2::geom_hline(yintercept = thr_y, linetype = "dashed")
    
    p1 <- p1 +
      ggplot2::scale_color_identity(guide = "none") +
      ggplot2::scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      ggplot2::labs(x = NULL, y = "-log10(P)") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    p1_pl <- plotly::ggplotly(p1, tooltip = "text", source = src_combo)
    
    # ============================
    # p2 base: eixos correctes sempre
    # ============================
    p2_base <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_blank() +
      ggplot2::scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::labs(x = "Genome", y = "-log10(P) [Catalog]") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    
    p2_pl <- plotly::ggplotly(p2_base, source = src_combo)
    
    # ============================
    # p2: Catalog hits (si n'hi ha)
    # ============================
    cat2 <- NULL
    y_max <- 1
    cat_ok <- !is.null(cat_dt) && is.data.frame(cat_dt) && nrow(cat_dt) > 0
    
    if (cat_ok) {
      
      pos_col <- pick_col(cat_dt, c("POS","POS_std","BP","position","POSITION"))
      p_col   <- pick_col(cat_dt, c("P_VALUE","p_value","P","p","PVAL","P-VALUE"))
      rs_col  <- pick_col(cat_dt, c("rsid","RSID","SNPS","SNP","snp"))
      
      if (!is.null(pos_col) && !is.null(p_col)) {
        
        cat2 <- cat_dt %>%
          dplyr::mutate(
            CHR0 = if ("CHR" %in% names(.)) as.character(CHR) else if ("CHR_std" %in% names(.)) as.character(CHR_std) else NA_character_,
            CHR0 = toupper(trimws(CHR0)),
            CHR0 = gsub("^CHR", "", CHR0, ignore.case = TRUE),
            CHR  = dplyr::case_when(
              CHR0 %in% as.character(1:22) ~ CHR0,
              CHR0 %in% c("23","X") ~ "X",
              CHR0 %in% c("24","Y") ~ "Y",
              CHR0 %in% c("MT","M") ~ "MT",
              TRUE ~ NA_character_
            ),
            POSn = suppressWarnings(readr::parse_number(as.character(.data[[pos_col]]))),
            Pn   = parse_p_robust(.data[[p_col]]),
            rsid = if (!is.null(rs_col)) as.character(.data[[rs_col]]) else NA_character_
          ) %>%
          dplyr::filter(!is.na(CHR), !is.na(POSn), !is.na(Pn), Pn > 0) %>%
          dplyr::inner_join(
            (ref %>% dplyr::mutate(chr = as.character(chr)) %>% dplyr::select(chr, chr_cum)),
            by = c("CHR" = "chr")
          ) %>%
          dplyr::mutate(
            BPcum = POSn + chr_cum,
            value = -log10(Pn)
          )
        
        # highlight selected cluster if available
        cl_sel <- selected_cluster()
        sel_id <- if (!is.null(cl_sel) && nrow(cl_sel) == 1) as.character(cl_sel$cluster_id[1]) else NA_character_
        
        cat2 <- cat2 %>%
          dplyr::mutate(
            highlight = if ("cluster_id" %in% names(.) && !is.na(sel_id)) {
              ifelse(as.character(cluster_id) == sel_id, "cluster", "other")
            } else "other"
          )
        
        cat2$tooltip <- paste0(
          "<b>Catalog hit</b>",
          "<br>rs: ", ifelse(!is.na(cat2$rsid) & nzchar(cat2$rsid), cat2$rsid, "(NA)"),
          "<br>CHR: ", cat2$CHR,
          "<br>POS: ", cat2$POSn,
          "<br>-log10(P): ", round(cat2$value, 2),
          if ("cluster_id" %in% names(cat2)) paste0("<br>cluster_id: ", as.character(cat2$cluster_id)) else ""
        )
        
        y_max <- max(cat2$value[is.finite(cat2$value)], na.rm = TRUE)
        if (!is.finite(y_max) || y_max <= 0) y_max <- 1
        
        p2 <- ggplot2::ggplot(cat2, ggplot2::aes(x = BPcum, y = value, text = tooltip)) +
          ggplot2::geom_point(ggplot2::aes(color = highlight), size = 1.2, alpha = 0.8) +
          ggplot2::scale_color_manual(values = c("other" = "grey60", "cluster" = "red"), guide = "none") +
          ggplot2::scale_x_continuous(
            limits = c(0, GENOME_END),
            breaks = axis_breaks,
            labels = axis_labels,
            expand = c(0, 0)
          ) +
          ggplot2::scale_y_continuous(limits = c(0, y_max * 1.25), expand = c(0, 0)) +
          ggplot2::labs(x = "Genome", y = "-log10(P) [Catalog]") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
        
        p2_pl <- plotly::ggplotly(p2, tooltip = "text", source = src_combo)
      }
    }
    
    # ============================
    # CLUSTER SEGMENT band (com GTEx)
    # ============================
    cl_df  <- clusters_for_ui()
    cl_sel <- selected_cluster()
    sel_id <- if (!is.null(cl_sel) && nrow(cl_sel) == 1) as.character(cl_sel$cluster_id[1]) else NA_character_
    
    if (is.data.frame(cl_df) && nrow(cl_df) && is.data.frame(ref_map) && nrow(ref_map)) {
      
      clseg <- as.data.frame(cl_df) %>%
        dplyr::mutate(
          chr   = if ("chr" %in% names(.) ) suppressWarnings(as.integer(.data$chr)) else if ("CHR" %in% names(.)) suppressWarnings(as.integer(.data$CHR)) else NA_integer_,
          start = if ("start" %in% names(.)) suppressWarnings(as.numeric(.data$start)) else if ("start_bp" %in% names(.)) suppressWarnings(as.numeric(.data$start_bp)) else NA_real_,
          end   = if ("end" %in% names(.)) suppressWarnings(as.numeric(.data$end)) else if ("end_bp" %in% names(.)) suppressWarnings(as.numeric(.data$end_bp)) else NA_real_,
          cluster_id = dplyr::coalesce(
            if ("cluster_id" %in% names(.)) as.character(.data$cluster_id) else NA_character_,
            if ("cluster_chr_n" %in% names(.)) as.character(.data$cluster_chr_n) else NA_character_
          )
        ) %>%
        dplyr::transmute(
          cluster_id = cluster_id,
          chr_num    = suppressWarnings(as.integer(chr)),
          start_i    = suppressWarnings(as.numeric(start)),
          end_i      = suppressWarnings(as.numeric(end)),
          is_sel     = if (!is.na(sel_id)) (cluster_id == sel_id) else FALSE
        ) %>%
        dplyr::filter(!is.na(cluster_id), nzchar(cluster_id),
                      is.finite(chr_num), is.finite(start_i), is.finite(end_i),
                      end_i >= start_i) %>%
        dplyr::distinct(cluster_id, chr_num, start_i, end_i, is_sel) %>%
        dplyr::left_join(ref_map, by = "chr_num") %>%
        dplyr::filter(is.finite(chr_cum)) %>%
        dplyr::transmute(
          cluster_id = cluster_id,
          x0   = pmax(0, start_i + chr_cum),
          x1   = pmin(GENOME_END, end_i + chr_cum),
          xmid = (x0 + x1) / 2,
          is_sel = is_sel,
          text = paste0(
            "Cluster: ", cluster_id,
            "<br>chr", chr_label_plink(chr_num), ":",
            format(start_i, scientific = FALSE), "-",
            format(end_i, scientific = FALSE)
          )
        ) %>%
        dplyr::filter(is.finite(x0), is.finite(x1), x1 >= x0) %>%
        dplyr::arrange(x0)
      
      if (is.data.frame(clseg) && nrow(clseg) > 0) {
        
        # Posició: si hi ha hits -> proporcional a y_max; si no -> dins [0..1]
        if (!is.null(cat2) && is.data.frame(cat2) && nrow(cat2)) {
          y_seg  <- y_max * 1.10   # dins del límit y_max*1.25
          y_tick <- y_max * 0.015
          y_txt  <- y_max * 1.18
        } else {
          y_seg  <- 0.90
          y_tick <- 0.02
          y_txt  <- 0.95
        }
        
        clseg$y_seg  <- y_seg
        clseg$y0tick <- y_seg - y_tick
        clseg$y1tick <- y_seg + y_tick
        clseg$y_txt  <- y_txt
        
        # segments + ticks
        p2_pl <- p2_pl %>%
          plotly::add_segments(
            data = clseg,
            x = ~x0, xend = ~x1,
            y = ~y_seg, yend = ~y_seg,
            inherit = FALSE,
            line = list(width = 3),
            hoverinfo = "text", text = ~text,
            showlegend = FALSE
          ) %>%
          plotly::add_segments(
            data = clseg,
            x = ~x0, xend = ~x0,
            y = ~y0tick, yend = ~y1tick,
            inherit = FALSE, line = list(width = 2),
            hoverinfo = "none", showlegend = FALSE
          ) %>%
          plotly::add_segments(
            data = clseg,
            x = ~x1, xend = ~x1,
            y = ~y0tick, yend = ~y1tick,
            inherit = FALSE, line = list(width = 2),
            hoverinfo = "none", showlegend = FALSE
          )
        
        # labels (més estable amb add_text que add_annotations dins subplot)
        p2_pl <- p2_pl %>%
          plotly::add_annotations(
            data = clseg,
            x = ~xmid, y = ~y_txt,
            text = ~cluster_id,
            xref = "x", yref = "y",
            showarrow = FALSE,
            textangle = 45,
            font = list(size = 10),
            xanchor = "center",
            yanchor = "bottom"
          )
        
        # highlight selected cluster (més gruixut)
        if ("is_sel" %in% names(clseg) && any(clseg$is_sel)) {
          seg2 <- clseg[clseg$is_sel, , drop = FALSE]
          p2_pl <- p2_pl %>% plotly::add_segments(
            data = seg2,
            x = ~x0, xend = ~x1,
            y = ~y_seg, yend = ~y_seg,
            inherit = FALSE,
            line = list(width = 6),
            hoverinfo = "text", text = ~text,
            showlegend = FALSE
          )
        }
      }
    }
    
    # ============================
    # combine + events
    # ============================
    out <- plotly::subplot(
      p1_pl, p2_pl,
      nrows = 2, shareX = TRUE,
      heights = c(0.55, 0.45),
      titleY = TRUE
    )
    
    out$x$source <- src_combo
    
    out <- plotly::event_register(out, "plotly_click")
    out <- plotly::event_register(out, "plotly_relayout")
    out <- plotly::event_register(out, "plotly_doubleclick")
    
    out %>%
      plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian")
      ) %>%
      plotly::layout(showlegend = FALSE) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "Combined_Manhattan_plot",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
    
  })
  
  
  
  # ---------------------------------------------------------------------------
  # UCSC links (custom tracks) — region driven by:
  #   (1) Manhattan_combo zoom (plotly_relayout) if present
  #   (2) else selected cluster
  # ---------------------------------------------------------------------------
  # ===========================
  # UCSC viewer (GWAS + Catalog)  # source="manhattan_combo"
  # ===========================
  
  safe_label_rsid <- function(df) {
    rscol <- if ("rsid" %in% names(df)) "rsid" else detect_rsid_col(df)
    if (!is.null(rscol) && rscol %in% names(df)) {
      out <- as.character(df[[rscol]])
      bad <- is.na(out) | out == ""
      out[bad] <- paste0("rs_unknown_", seq_len(sum(bad)))
      return(out)
    }
    paste0("rs_unknown_", seq_len(nrow(df)))
  }
  
  safe_label_catalog <- function(df) {
    cand <- intersect(c("SNPS","snp","rsid","RSID","variant_id","variant","ID","id","SNP_ID_CURRENT"), names(df))
    if (length(cand)) {
      out <- as.character(df[[cand[1]]])
      bad <- is.na(out) | out == ""
      out[bad] <- paste0("catalog_hit_", seq_len(sum(bad)))
      return(out)
    }
    paste0("catalog_hit_", seq_len(nrow(df)))
  }
  
  window_selected <- reactiveVal(NULL)
  ucsc_region     <- reactiveVal(NULL)
  
  # --- helper robust: extreu rang X d'un relayout (subplot: xaxis o xaxis2) ---
  extract_xrange_from_relayout <- function(ev) {
    
    if (is.null(ev)) return(NULL)
    
    if (is.data.frame(ev)) {
      if (nrow(ev) == 0) return(NULL)
      getv <- function(nm) ev[[nm]][1]
    } else if (is.list(ev)) {
      getv <- function(nm) ev[[nm]]
    } else {
      return(NULL)
    }
    
    is_autorange_true <- function(v) {
      if (is.null(v)) return(FALSE)
      if (length(v) >= 1) v <- v[1]
      identical(v, TRUE) || identical(v, 1L) || identical(v, 1) ||
        (is.character(v) && tolower(v) %in% c("true","t","1","yes","y"))
    }
    
    if (is_autorange_true(getv("xaxis.autorange")) ||
        is_autorange_true(getv("xaxis2.autorange"))) {
      return(NULL)
    }
    
    x0 <- getv("xaxis.range[0]")
    x1 <- getv("xaxis.range[1]")
    
    if (is.null(x0) || is.null(x1)) {
      x0 <- getv("xaxis2.range[0]")
      x1 <- getv("xaxis2.range[1]")
    }
    
    if ((is.null(x0) || is.null(x1)) && !is.null(getv("xaxis.range"))) {
      rr <- getv("xaxis.range")
      if (length(rr) >= 2) { x0 <- rr[1]; x1 <- rr[2] }
    }
    if ((is.null(x0) || is.null(x1)) && !is.null(getv("xaxis2.range"))) {
      rr <- getv("xaxis2.range")
      if (length(rr) >= 2) { x0 <- rr[1]; x1 <- rr[2] }
    }
    
    if (is.null(x0) || is.null(x1)) return(NULL)
    
    x0 <- suppressWarnings(as.numeric(x0))
    x1 <- suppressWarnings(as.numeric(x1))
    if (!is.finite(x0) || !is.finite(x1)) return(NULL)
    
    list(xmin = min(x0, x1), xmax = max(x0, x1))
  }
  
  # --- Convert BPcum -> chr + pos (hg38) using .ref_hg38 ---
  coord_from_bp_cum <- function(x, ref = .ref_hg38) {
    
    if (is.null(x) || length(x) == 0) return(list(chr = NULL, pos = NULL))
    x <- suppressWarnings(as.numeric(x[1]))
    if (!is.finite(x)) return(list(chr = NULL, pos = NULL))
    
    ref2 <- ref %>%
      dplyr::transmute(
        chr     = suppressWarnings(as.integer(chr)),
        chr_cum = suppressWarnings(as.numeric(chr_cum)),
        len     = suppressWarnings(as.numeric(len))
      ) %>%
      dplyr::filter(is.finite(chr), is.finite(chr_cum), is.finite(len), len > 0) %>%
      dplyr::arrange(chr_cum)
    
    if (!nrow(ref2)) return(list(chr = NULL, pos = NULL))
    
    idx <- which(x >= ref2$chr_cum & x <= (ref2$chr_cum + ref2$len))
    if (!length(idx)) return(list(chr = NULL, pos = NULL))
    i <- idx[1]
    
    pos <- as.integer(round(x - ref2$chr_cum[i]))
    if (!is.finite(pos)) return(list(chr = NULL, pos = NULL))
    if (pos < 1L) pos <- 1L
    if (pos > ref2$len[i]) pos <- as.integer(ref2$len[i])
    
    chr_lab <- chr_label_plink(as.integer(ref2$chr[i]))  # 23->X,24->Y,25->M
    list(chr = chr_lab, pos = pos)
  }
  
  # 1) escolta el relayout del combo correcte
  # 1) escolta el relayout del combo correcte
  observeEvent(plotly::event_data("plotly_relayout", source = "manhattan_combo"), {
    
    ev <- plotly::event_data("plotly_relayout", source = "manhattan_combo")
    if (is.null(ev) || !length(ev)) return()
    
    # ignora relayouts "brossa" d'inicialització
    if (all(names(ev) %in% c("autosize", "xaxis.autorange", "yaxis.autorange",
                             "xaxis2.autorange", "yaxis2.autorange"))) return()
    
    # extreu rang robust (xaxis o xaxis2, range[0]/range[1] o range vector)
    xr <- extract_xrange_from_relayout(ev)
    if (is.null(xr)) return()
    
    window_selected(xr)
    
    cat("[UCSC] relayout -> window_selected xmin=", xr$xmin, " xmax=", xr$xmax, "\n")
    flush.console()
    
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  
  # --- hits dins la finestra (per BPcum) ---
  get_window_hits_gwas <- reactive({
    win <- window_selected(); req(win)
    
    dfp <- dfp_manhattan()
    req(is.data.frame(dfp), nrow(dfp) > 0)
    
    dfp %>%
      dplyr::filter(is.finite(BPcum), BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  get_window_hits_catalog <- reactive({
    win <- window_selected()
    if (is.null(win) || !is.finite(win$xmin) || !is.finite(win$xmax)) {
      return(tibble::tibble())
    }
    
    # --- fallback chain (posa aquí la teva font "all catalog" si existeix) ---
    cat0 <- tryCatch(catalog_table_filtered(), error = function(e) NULL)
    
    if (is.null(cat0) || !is.data.frame(cat0) || nrow(cat0) == 0) {
      # si tens una taula "no filtrada", prova-la aquí:
      cat0 <- tryCatch(catalog_table(), error = function(e) NULL)       # <- canvia si no existeix
    }
    if (is.null(cat0) || !is.data.frame(cat0) || nrow(cat0) == 0) {
      cat0 <- tryCatch(catalog_data(), error = function(e) NULL)        # <- canvia si no existeix
    }
    
    if (is.null(cat0) || !is.data.frame(cat0) || nrow(cat0) == 0) {
      return(tibble::tibble())
    }
    
    # necessitem CHR i POS
    if (!all(c("CHR","POS") %in% names(cat0))) {
      return(tibble::tibble())
    }
    
    ref <- .ref_hg38
    if (is.null(ref) || !is.data.frame(ref) || !all(c("chr","chr_cum","len") %in% names(ref))) {
      return(tibble::tibble())
    }
    
    # Normalitza CHR a "1..22/X/Y/MT"
    chrN <- toupper(trimws(as.character(cat0$CHR)))
    chrN <- gsub("^chr", "", chrN, ignore.case = TRUE)
    chrN <- gsub("^CHR", "", chrN, ignore.case = TRUE)
    chrN[chrN %in% c("23","X")] <- "X"
    chrN[chrN %in% c("24","Y")] <- "Y"
    chrN[chrN %in% c("25","MT","M")] <- "MT"
    
    pos <- suppressWarnings(as.numeric(cat0$POS))
    
    # ref amb chr en el mateix format ("1","2",...,"X","Y","MT")
    ref2 <- ref %>%
      dplyr::transmute(
        chrN    = toupper(trimws(as.character(chr))),
        chr_cum = suppressWarnings(as.numeric(chr_cum)),
        len     = suppressWarnings(as.numeric(len))
      ) %>%
      dplyr::filter(!is.na(chrN), is.finite(chr_cum), is.finite(len), len > 0)
    
    out <- cat0 %>%
      dplyr::mutate(chrN = chrN, pos = pos) %>%
      dplyr::filter(!is.na(chrN), is.finite(pos)) %>%
      dplyr::inner_join(ref2 %>% dplyr::select(chrN, chr_cum, len), by = "chrN") %>%
      dplyr::mutate(BPcum = as.numeric(pos) + as.numeric(chr_cum)) %>%
      dplyr::filter(is.finite(BPcum), BPcum >= win$xmin, BPcum <= win$xmax)
    
    out
  })
  
  
  # --- construeix BED12 track (GWAS / Catalog) ---
  make_track_df_gwas <- function(df, R, G, B) {
    if (is.null(df) || !nrow(df)) return(tibble::tibble())
    if (!all(c("CHR","BP") %in% names(df))) {
      # si només tenim dfp_manhattan -> CHR/BP hi solen ser, però fem fallback
      validate(need(all(c("CHR","BP") %in% names(df)), "GWAS window hits need CHR and BP columns."))
    }
    
    chrom <- paste0("chr", chr_label_plink(as.integer(df$CHR)))
    start <- as.numeric(df$BP) - 1
    end   <- as.numeric(df$BP)
    lbl   <- safe_label_rsid(df)
    
    tibble::tibble(
      chrom       = chrom,
      start       = start,
      end         = end,
      name        = lbl,
      score       = 0L,
      strand      = "+",
      thickStart  = start,
      thickEnd    = end,
      itemRgb     = paste0(R, ",", G, ",", B),
      blockCount  = 1L,
      blockSizes  = end - start,
      blockStarts = 0L
    )
  }
  
  safe_label_catalog <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(character())
    }
    
    if ("rsid" %in% names(df)) {
      out <- as.character(df$rsid)
      bad <- is.na(out) | !nzchar(out)
      out[bad] <- paste0("catalog_hit_", seq_len(sum(bad)))
      return(out)
    }
    
    paste0("catalog_hit_", seq_len(nrow(df)))
  }
  
  make_track_df_catalog <- function(df, R = 227, G = 26, B = 28) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(tibble::tibble())
    
    # prefereix CHR/POS; si no hi són, prova chrN/pos
    if (all(c("CHR","POS") %in% names(df))) {
      chr_use <- as.character(df$CHR)
      pos_use <- suppressWarnings(as.integer(df$POS))
    } else {
      validate(need(all(c("chrN","pos") %in% names(df)), "Catalog hits need CHR/POS or chrN/pos."))
      chr_use <- as.character(df$chrN)
      pos_use <- suppressWarnings(as.integer(df$pos))
    }
    
    chr_use <- toupper(trimws(chr_use))
    chr_use <- gsub("^chr", "", chr_use, ignore.case = TRUE)
    chr_use <- gsub("^CHR", "", chr_use, ignore.case = TRUE)
    chr_use[chr_use %in% c("23","X")] <- "X"
    chr_use[chr_use %in% c("24","Y")] <- "Y"
    chr_use[chr_use %in% c("25","MT","M")] <- "MT"
    
    ok <- !is.na(chr_use) & is.finite(pos_use)
    if (!any(ok)) return(tibble::tibble())
    
    chrom <- paste0("chr", chr_use[ok])
    start <- pmax(0L, pos_use[ok] - 1L)
    end   <- pmax(start + 1L, pos_use[ok])
    
    lbl <- safe_label_catalog(df[ok, , drop = FALSE])
    
    tibble::tibble(
      chrom       = chrom,
      start       = start,
      end         = end,
      name        = lbl,
      score       = 0L,
      strand      = "+",
      thickStart  = start,
      thickEnd    = end,
      itemRgb     = paste0(R, ",", G, ",", B),
      blockCount  = 1L,
      blockSizes  = end - start,
      blockStarts = 0L
    )
  }
  
  
  clean_track <- function(df) {
    if (is.null(df) || !nrow(df)) return(tibble::tibble())
    df %>%
      dplyr::filter(!is.na(chrom), chrom != "chrNA",
                    is.finite(start), is.finite(end)) %>%
      dplyr::mutate(
        start = pmax(0L, as.integer(start)),
        end   = pmax(start + 1L, as.integer(end))
      ) %>%
      dplyr::distinct(chrom, start, end, name, .keep_all = TRUE)
  }
  
  # 2) Quan canvia la finestra: calcula region UCSC + tracks
  observeEvent(window_selected(), {
    win <- window_selected()
    if (is.null(win)) return()
    
    cat("[UCSC] window_selected xmin=", win$xmin, " xmax=", win$xmax, "\n"); flush.console()
    
    left  <- coord_from_bp_cum(win$xmin)
    right <- coord_from_bp_cum(win$xmax)
    
    if (is.null(left$chr) || is.null(right$chr) || left$chr != right$chr) {
      ucsc_region(NULL)
      session$userData$track_gwas_data    <- NULL
      session$userData$track_catalog_data <- NULL
      cat("[UCSC] window spans multiple chr -> reset tracks\n"); flush.console()
      return()
    }
    
    region <- sprintf("chr%s:%d-%d", left$chr, left$pos, right$pos)
    ucsc_region(region)
    cat("[UCSC] region=", region, "\n"); flush.console()
    
    gwas_hits <- tryCatch(get_window_hits_gwas(), error = function(e) {
      cat("[UCSC] get_window_hits_gwas ERROR:", conditionMessage(e), "\n"); flush.console()
      NULL
    })
    
    cat_hits <- tryCatch(get_window_hits_catalog(), error = function(e) {
      cat("[UCSC] get_window_hits_catalog ERROR class=", paste(class(e), collapse=","), "\n")
      cat("[UCSC] msg:", conditionMessage(e), "\n"); flush.console()
      NULL
    })
    
    cat("[UCSC] gwas window hits:", if (is.data.frame(gwas_hits)) nrow(gwas_hits) else NA, "\n"); flush.console()
    cat("[UCSC] catalog window hits:", if (is.data.frame(cat_hits)) nrow(cat_hits) else NA, "\n"); flush.console()
    
    # IMPORTANT: cat_hits ja porta CHR/POS originals, i també chrN/pos/BPcum afegits
    df_gwas <- clean_track(make_track_df_gwas(gwas_hits, 31, 120, 180))
    df_cat  <- clean_track(make_track_df_catalog(cat_hits, 227, 26, 28))
    
    session$userData$track_gwas_data    <- df_gwas
    session$userData$track_catalog_data <- df_cat
    
    cat("[UCSC] tracks: gwas=", nrow(df_gwas), " catalog=", nrow(df_cat), "\n"); flush.console()
    
  }, ignoreInit = TRUE)
  
  
  
  
  make_ucsc_track_text <- function(name, df, url_tpl) {
    header <- sprintf('track name="%s" description="%s" visibility=pack itemRgb=On url="%s"',
                      name, name, url_tpl)
    if (is.null(df) || !nrow(df)) return(header)
    cols <- c("chrom","start","end","name","score","strand",
              "thickStart","thickEnd","itemRgb",
              "blockCount","blockSizes","blockStarts")
    body <- apply(df[, cols], 1, paste, collapse = "\t")
    paste(c(header, body), collapse = "\n")
  }
  
  make_ucsc_url <- function(region, track_text) {
    base <- "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38"
    reg  <- paste0("&position=", utils::URLencode(region, reserved = TRUE))
    trk  <- paste0("&hgt.customText=", utils::URLencode(track_text, reserved = TRUE))
    paste0(base, reg, trk)
  }
  
  output$ucsc_link_gwas <- renderUI({
    region <- ucsc_region()
    df     <- session$userData$track_gwas_data
    req(!is.null(region), !is.null(df), nrow(df) > 0)
    
    trk <- make_ucsc_track_text("GWAS_hits", df, "https://www.ncbi.nlm.nih.gov/snp/$$")
    url <- make_ucsc_url(region, trk)
    
    tags$a(href = url, target = "_blank", "Open UCSC – GWAS hits (dbSNP links)")
  })
  
  output$ucsc_link_catalog <- renderUI({
    region <- ucsc_region()
    df     <- session$userData$track_catalog_data
    req(!is.null(region), !is.null(df), nrow(df) > 0)
    
    trk <- make_ucsc_track_text("Catalog_hits", df, "https://www.ebi.ac.uk/gwas/variants/$$")
    url <- make_ucsc_url(region, trk)
    
    tags$a(href = url, target = "_blank", "Open UCSC – Catalog hits (GWAS Catalog links)")
  })
  
  safe_nrow <- function(x) {
    if (is.null(x)) return(0L)
    if (is.data.frame(x)) return(nrow(x))
    if (inherits(x, "tbl")) return(nrow(x))
    0L
  }
  
  output$debug_ucsc_state <- renderUI({
    region <- ucsc_region() %||% "NULL"
    gwas_n <- safe_nrow(session$userData$track_gwas_data)
    cat_n  <- safe_nrow(session$userData$track_catalog_data)
    win    <- window_selected()
    
    tags$pre(
      style="background-color:#f6f6f6; border:1px solid #ddd; padding:10px; font-family: 'Courier New', monospace;",
      paste0(
        "window = ", if (is.null(win)) "NULL" else paste0(win$xmin, " .. ", win$xmax),
        "\nregion = ", region,
        "\nGWAS hits in track = ", gwas_n,
        "\nCatalog hits in track = ", cat_n
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # Catalog visualizations (fixed gene splitting + alternating bar colors)
  # ---------------------------------------------------------------------------
  
  # Helper: split per-row into a LIST (so unnest works)
  split_multi_list <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    lapply(x, function(s) {
      s <- trimws(s)
      if (!nzchar(s)) return(character(0))
      # separators commonly seen in GWAS Catalog mapped genes
      parts <- unlist(strsplit(s, "\\s*[,;|/]\\s*|\\s+\\+\\s+|\\s+-\\s+|\\s+and\\s+|\\s*&\\s*", perl = TRUE))
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      unique(parts)
    })
  }
  
  output$catalog_disease_bar <- renderPlot({
    dt <- catalog_data()
    if (is.null(dt) || !nrow(dt)) return(NULL)
    dis_col <- pick_col(dt, c("DISEASE","Disease","disease"))
    if (is.null(dis_col)) return(NULL)
    
    top <- dt %>%
      dplyr::filter(!is.na(.data[[dis_col]]), nzchar(.data[[dis_col]])) %>%
      dplyr::count(.data[[dis_col]], sort = TRUE) %>%
      dplyr::slice_head(n = 25) %>%
      dplyr::mutate(.alt = factor(dplyr::row_number() %% 2))
    
    ggplot2::ggplot(top, ggplot2::aes(x = n, y = stats::reorder(.data[[dis_col]], n), fill = .alt)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("0" = "orange", "1" = "darkgreen"), guide = "none") +
      ggplot2::labs(x = "Number of catalog hits", y = "Disease") +
      ggplot2::theme_minimal(base_size = 13)
  })
  
  output$catalog_mapped_trait_bar <- renderPlot({
    dt <- catalog_data()
    if (is.null(dt) || !nrow(dt)) return(NULL)
    trait_col <- pick_col(dt, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    if (is.null(trait_col)) return(NULL)
    
    top <- dt %>%
      dplyr::filter(!is.na(.data[[trait_col]]), nzchar(.data[[trait_col]])) %>%
      dplyr::count(.data[[trait_col]], sort = TRUE) %>%
      dplyr::slice_head(n = 25) %>%
      dplyr::mutate(.alt = factor(dplyr::row_number() %% 2))
    
    ggplot2::ggplot(top, ggplot2::aes(x = n, y = stats::reorder(.data[[trait_col]], n), fill = .alt)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("0" = "orange", "1" = "darkgreen"), guide = "none") +
      ggplot2::labs(x = "Number of catalog hits", y = "Mapped Trait") +
      ggplot2::theme_minimal(base_size = 13)
  })
  
  output$catalog_gene_trait_plot <- renderPlot({
    dt <- catalog_data()
    if (is.null(dt) || !nrow(dt)) return(NULL)
    
    gene_col  <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    trait_col <- pick_col(dt, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    if (is.null(gene_col) || is.null(trait_col)) return(NULL)
    
    edges <- dt %>%
      dplyr::filter(!is.na(.data[[gene_col]]), !is.na(.data[[trait_col]])) %>%
      dplyr::transmute(
        gene  = as.character(.data[[gene_col]]),
        trait = as.character(.data[[trait_col]])
      ) %>%
      dplyr::mutate(gene = split_multi_list(gene)) %>%      # <- LIST per row
      tidyr::unnest(cols = c(gene)) %>%
      dplyr::filter(nzchar(gene), nzchar(trait)) %>%
      dplyr::distinct()
    
    if (!nrow(edges)) return(NULL)
    
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    
    # Make it a bit more interpretable: gene vs trait colors
    genes <- unique(edges$gene)
    traits <- unique(edges$trait)
    V(g)$type <- ifelse(V(g)$name %in% genes, "gene", "trait")
    V(g)$color <- ifelse(V(g)$type == "gene", "#1f77b4", "#d62728")
    V(g)$frame.color <- "grey30"
    
    set.seed(1)
    plot(
      g,
      layout = igraph::layout_with_fr(g),
      vertex.size = 6,
      vertex.label.cex = 0.65,
      vertex.label.color = "black",
      edge.color = "grey70"
    )
  })
  
  output$catalog_interpret <- renderText({
    dt <- catalog_data()
    if (is.null(dt) || !nrow(dt)) return("No Catalog hits loaded.")
    
    gene_col  <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    trait_col <- pick_col(dt, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    dis_col   <- pick_col(dt, c("DISEASE","Disease","disease"))
    
    n_hits <- nrow(dt)
    
    n_clusters <- if ("cluster_id" %in% names(dt)) {
      length(unique(stats::na.omit(as.character(dt$cluster_id))))
    } else 0L
    
    n_genes <- if (!is.null(gene_col)) {
      length(split_multi_unique(dt[[gene_col]]))
    } else NA_integer_
    
    n_traits <- if (!is.null(trait_col)) {
      length(unique(stats::na.omit(as.character(dt[[trait_col]]))))
    } else NA_integer_
    
    n_dis <- if (!is.null(dis_col)) {
      length(unique(stats::na.omit(as.character(dt[[dis_col]]))))
    } else NA_integer_
    
    paste0(
      "Hits: ", n_hits, "\n",
      "Clusters: ", n_clusters, "\n",
      "Unique genes (approx): ", n_genes, "\n",
      "Unique traits: ", n_traits, "\n",
      "Unique diseases: ", n_dis, "\n"
    )
  })
  
  # ---------------------------------------------------------------------------
  # Disease / trait enrichment (ORA via Fisher) — scope + in-panel cluster picker + Run button
  # ---------------------------------------------------------------------------
  
  # ---- Catalog universe (reference) for enrichment ----
  catalog_universe_cache <- reactiveVal(NULL)
  
  catalog_universe_data <- reactive({
    x <- catalog_universe_cache()
    if (is.data.frame(x) && nrow(x) > 0) return(x)
    
    # robust paths (app dir + relative)
    app_dir <- tryCatch(normalizePath(".", winslash = "/", mustWork = TRUE), error = function(e) getwd())
    candidates <- c(
      file.path("www", "gwas_catalog_simplified.rds"),
      file.path(app_dir, "www", "gwas_catalog_simplified.rds")
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f) || !nzchar(f)) return(NULL)
    
    dt0 <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.data.frame(dt0) || !nrow(dt0)) return(NULL)
    
    catalog_universe_cache(dt0)
    dt0
  })
  
  # ---- trigger: run only on button click ----
  dt_enrich_trigger <- reactiveVal(0L)
  
  observeEvent(input$run_dt_enrich, {
    dt_enrich_trigger(dt_enrich_trigger() + 1L)
  }, ignoreInit = TRUE)
  
  # ---- cluster picker UI (like other panels) ----
  output$dt_cluster_ui <- renderUI({
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available."))
    
    chr_choices <- sort(unique(as.integer(cl$chr)))
    chr_labels  <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("dt_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("dt_cluster_id_ui")
    )
  })
  
  output$dt_cluster_id_ui <- renderUI({
    cl <- clusters_val(); req(is.data.frame(cl), nrow(cl) > 0)
    req(input$dt_chr)
    
    x <- cl %>%
      dplyr::filter(as.integer(chr) == as.integer(input$dt_chr)) %>%
      dplyr::arrange(start, end)
    
    validate(need(nrow(x) > 0, "No clusters for this chromosome."))
    selectInput("dt_cluster_id", "Cluster", choices = x$cluster_id, selected = x$cluster_id[1])
  })
  
  # ---- subset selector (dt2) according to scope ----
  # ---- subset selector (dt2) according to scope ----
  dt_subset_for_enrich <- reactive({
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits."))
    
    scope <- input$dt_scope %||% "pick"
    
    validate(need("cluster_id" %in% names(dt),
                  "Catalog hits missing cluster_id (cannot do cluster enrichment)."))
    
    if (identical(scope, "global")) {
      # GLOBAL = all hits that belong to ANY cluster (union of clusters)
      dt2 <- dt %>% dplyr::filter(!is.na(cluster_id), nzchar(cluster_id))
      validate(need(nrow(dt2) > 0, "No clustered hits found (cluster_id empty/NA everywhere)."))
      return(dt2)
    }
    
    if (identical(scope, "selected")) {
      cl <- selected_cluster()
      validate(need(is.data.frame(cl) && nrow(cl) > 0, "No cluster selected in Clusters table."))
      return(dt %>% dplyr::filter(cluster_id == cl$cluster_id[1]))
    }
    
    if (identical(scope, "pick")) {
      req(input$dt_cluster_id)
      return(dt %>% dplyr::filter(cluster_id == input$dt_cluster_id))
    }
    
    if (identical(scope, "gene")) {
      req(input$dt_gene)
      
      gene_col <- pick_col(dt, c(
        "MAPPED_GENE","mapped_gene",
        "REPORTED_GENE(S)","REPORTED_GENE","REPORTED_GENES",
        "GENE","gene","Gene","MAPPED_GENES"
      ))
      validate(need(!is.null(gene_col), "Gene column not found in current catalog."))
      
      # robust row-wise membership without false positives (wrap with ';')
      norm_gene_field <- function(x) {
        x <- toupper(trimws(as.character(x)))
        x[is.na(x)] <- ""
        x <- gsub("\\s*[,;|/]\\s*", ";", x)
        x <- gsub("\\s+and\\s+", ";", x, ignore.case = TRUE)
        x <- gsub("\\s*\\+\\s*", ";", x)
        x <- gsub("\\s*&\\s*",  ";", x)
        x <- gsub("\\s+", "", x)
        paste0(";", x, ";")
      }
      
      g <- toupper(trimws(as.character(input$dt_gene)))
      v <- norm_gene_field(dt[[gene_col]])
      keep <- grepl(paste0(";", g, ";"), v, fixed = TRUE)
      
      dt2 <- dt[keep, , drop = FALSE]
      validate(need(nrow(dt2) > 0, paste0("No hits found for gene: ", g)))
      return(dt2)
    }
    
    # fallback
    validate(need(FALSE, "Unknown dt_scope value."))
  })
  
  # ---- info note ----
  # ---- info note (scope + enrichment explanation) ----
  output$dt_scope_note <- renderUI({
    dt <- tryCatch(catalog_data(), error = function(e) NULL)
    if (!is.data.frame(dt) || !nrow(dt)) return(NULL)
    
    dt_uni <- catalog_universe_data()
    uni_ok <- is.data.frame(dt_uni) && nrow(dt_uni) > 0
    uni_n  <- if (uni_ok) nrow(dt_uni) else nrow(dt)
    
    scope <- input$dt_scope %||% "pick"
    msg <- NULL
    
    # -----------------------------
    # Scope-specific summary (subset/universe)
    # -----------------------------
    if (identical(scope, "global")) {
      if ("cluster_id" %in% names(dt)) {
        n_cl <- sum(!is.na(dt$cluster_id) & nzchar(dt$cluster_id))
        msg <- paste0(
          "<b>Global mode:</b> subset = <b>clustered hits</b> (", n_cl, "), universe = ",
          if (uni_ok) "<b>reference catalog</b>" else "<b>current catalog</b>",
          " (", uni_n, ")."
        )
        if (!uni_ok) {
          msg <- paste0(
            msg,
            "<hr style='margin:10px 0;'>",
            "<b>Enrichment (ORA):</b> For each term (Disease / Trait), we test whether it is <b>over-represented</b> in the subset vs the universe ",
            "using a one-sided Fisher/hypergeometric test (<i>alternative = greater</i>).",
            "<br>",
            "<b>Table columns:</b> ",
            "<code>in_subset</code> = <i>n/total subset</i>; ",
            "<code>n_universe</code> = <i>n/total universe</i>; ",
            "<code>OR</code> = odds ratio (with pseudocounts); ",
            "<code>p</code> = raw p-value; ",
            "<code>FDR</code> = BH-adjusted p-value.",
            "<br>",
            "<b>Plot:</b> bar length = <b>counts in subset</b>; color = <b>-log10(FDR)</b> (red → orange → yellow = more significant)."
          )
        }
      }
    }
    
    if (identical(scope, "pick")) {
      if (!is.null(input$dt_cluster_id) && nzchar(input$dt_cluster_id)) {
        # best-effort counts
        n_sub <- sum(dt$cluster_id == input$dt_cluster_id, na.rm = TRUE)
        msg <- paste0(
          "<b>Pick mode:</b> subset = hits in <b>", htmltools::htmlEscape(input$dt_cluster_id),
          "</b> (", n_sub, "), universe = ",
          if (uni_ok) "<b>reference catalog</b>" else "<b>current catalog</b>",
          " (", uni_n, ")."
        )
        if (!uni_ok) {
          msg <- paste0(
            msg,
            "<br><b>Warning:</b> missing <code>www/gwas_catalog_simplified.rds</code> → using current catalog as universe."
          )
        }
      } else {
        msg <- "<b>Pick mode:</b> choose a chromosome + cluster to define the subset."
      }
    }
    
    if (identical(scope, "selected")) {
      cl <- selected_cluster()
      if (is.null(cl) || !nrow(cl)) {
        msg <- "<b>Selected cluster mode:</b> No cluster selected in the <i>Clusters</i> table. Select one, or use <i>Pick a cluster here</i>."
      } else {
        cid <- as.character(cl$cluster_id[1])
        n_sub <- sum(dt$cluster_id == cid, na.rm = TRUE)
        msg <- paste0(
          "<b>Selected cluster mode:</b> subset = hits in <b>", htmltools::htmlEscape(cid),
          "</b> (", n_sub, "), universe = ",
          if (uni_ok) "<b>reference catalog</b>" else "<b>current catalog</b>",
          " (", uni_n, ")."
        )
        if (!uni_ok) {
          msg <- paste0(
            msg,
            "<br><b>Warning:</b> missing <code>www/gwas_catalog_simplified.rds</code> → using current catalog as universe."
          )
        }
      }
    }
    
    if (identical(scope, "gene")) {
      g <- input$dt_gene %||% ""
      if (nzchar(g)) {
        gene_hits <- 0L
        gene_col <- pick_col(dt, c(
          "MAPPED_GENE","mapped_gene",
          "REPORTED_GENE(S)","REPORTED_GENE","REPORTED_GENES",
          "GENE","gene","Gene","MAPPED_GENES"
        ))
        if (!is.null(gene_col)) {
          vv <- toupper(trimws(as.character(dt[[gene_col]])))
          vv[is.na(vv)] <- ""
          vv <- gsub("\\s*[,;|/]\\s*", ";", vv)
          vv <- gsub("\\s+and\\s+", ";", vv, ignore.case = TRUE)
          vv <- gsub("\\s*\\+\\s*", ";", vv)
          vv <- gsub("\\s*&\\s*",  ";", vv)
          vv <- gsub("\\s+", "", vv)
          vv <- paste0(";", vv, ";")
          gene_hits <- sum(grepl(paste0(";", toupper(g), ";"), vv, fixed = TRUE), na.rm = TRUE)
        }
        
        msg <- paste0(
          "<b>Gene mode:</b> subset = hits mapped/reported to <b>", htmltools::htmlEscape(g),
          "</b> (", gene_hits, "), universe = ",
          if (uni_ok) "<b>reference catalog</b>" else "<b>current catalog</b>",
          " (", uni_n, ")."
        )
        if (!uni_ok) {
          msg <- paste0(
            msg,
            "<br><b>Warning:</b> missing <code>www/gwas_catalog_simplified.rds</code> → using current catalog as universe."
          )
        }
      } else {
        msg <- "<b>Gene mode:</b> pick a gene to define the subset."
      }
    }
    
    # -----------------------------
    # Enrichment explanation (always appended when we have a scope msg)
    # -----------------------------
    if (!is.null(msg)) {
      msg <- paste0(
        msg,
        "<hr style='margin:10px 0;'>",
        "<b>Enrichment (ORA):</b> For each term (Disease / Trait), we test whether it is <b>over-represented</b> in the subset vs the universe ",
        "using a one-sided Fisher/hypergeometric test (<i>alternative = greater</i>).",
        "<br>",
        "<b>Table columns:</b> ",
        "<code>in_subset</code> = <i>n/total subset</i>; ",
        "<code>n_universe</code> = <i>n/total universe</i>; ",
        "<code>OR</code> = odds ratio (with pseudocounts); ",
        "<code>p</code> = raw p-value; ",
        "<code>FDR</code> = BH-adjusted p-value.",
        "<br>",
        "<b>Plot:</b> bar length = <b>counts in subset</b>; color = <b>-log10(FDR)</b> (red → orange → yellow = more significant)."
      )
    }
    
    if (!is.null(msg)) {
      return(HTML(paste0(
        "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
        msg,
        "</div>"
      )))
    }
    
    NULL
  })
  
  # -------------------------------
  # Formatting: 1 decimal; if <0.01 -> 3 decimals; never scientific
  # -------------------------------
  fmt_num_nosci <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(
      is.na(x),
      NA_character_,
      ifelse(
        abs(x) < 0.01,
        format(round(x, 3), nsmall = 3, scientific = FALSE, trim = TRUE),
        format(round(x, 1), nsmall = 1, scientific = FALSE, trim = TRUE)
      )
    )
  }
  
  # ---- Fisher ORA (FIX: split multi-terms + correct 2x2) ----
  
  fisher_enrich_terms_fast <- function(vec_terms, universe_terms,
                                       alternative = "greater",
                                       min_a = 2,        # min aparicions al subset
                                       min_Tt = 5,       # min aparicions a l'univers
                                       max_terms = 2000  # límit de termes a testar (després de filtres)
  ) {
    
    clean_vec <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      x
    }
    
    vec_terms      <- clean_vec(vec_terms)
    universe_terms <- clean_vec(universe_terms)
    
    if (!length(vec_terms) || !length(universe_terms)) return(tibble::tibble())
    
    n <- length(vec_terms)
    N <- length(universe_terms)
    
    hit_counts <- table(vec_terms)
    uni_counts <- table(universe_terms)
    
    terms <- intersect(names(hit_counts), names(uni_counts))
    if (!length(terms)) return(tibble::tibble())
    
    a  <- as.integer(hit_counts[terms])   # terme al subset
    Tt <- as.integer(uni_counts[terms])   # terme a l'univers
    
    # filtres per reduir càlcul
    keep <- (a >= min_a) & (Tt >= min_Tt) & (Tt <= N) & (a <= n)
    terms <- terms[keep]
    a  <- a[keep]
    Tt <- Tt[keep]
    
    if (!length(terms)) return(tibble::tibble())
    
    # (opcional) limita nombre de termes a testar, prioritzant els més freqüents al subset
    if (length(terms) > max_terms) {
      ord <- order(a, decreasing = TRUE)
      ord <- ord[seq_len(max_terms)]
      terms <- terms[ord]; a <- a[ord]; Tt <- Tt[ord]
    }
    
    b <- n - a
    c <- Tt - a
    d <- (N - n) - c
    
    # p-value (Fisher one-sided "greater") via hipergeomètrica
    if (identical(alternative, "greater")) {
      p <- stats::phyper(a - 1L, Tt, N - Tt, n, lower.tail = FALSE)
    } else if (identical(alternative, "less")) {
      p <- stats::phyper(a, Tt, N - Tt, n, lower.tail = TRUE)
    } else {
      # two-sided aproximat: fem servir fisher.test només si realment ho vols
      p <- rep(NA_real_, length(a))
    }
    
    # OR (aprox) amb pseudocomptes per evitar infinits
    or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    
    out <- tibble::tibble(
      term = terms, a = a, subset_n = n,
      uni_term = Tt, uni_N = N,
      or = or, p = p
    ) %>%
      dplyr::mutate(p_adj = stats::p.adjust(p, method = "BH")) %>%
      dplyr::arrange(p_adj, p)
    
    out
  }
  
  
  
  # ---- compute enrichment ONLY when triggered ----
  dt_enrich_res <- eventReactive(dt_enrich_trigger(), {
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits."))
    
    # Universe: prefer reference RDS; fallback to current dt if missing
    dt_uni <- catalog_universe_data()
    if (!is.data.frame(dt_uni) || !nrow(dt_uni)) dt_uni <- dt
    
    dt2 <- dt_subset_for_enrich()
    validate(need(is.data.frame(dt2) && nrow(dt2) > 0, "No hits in the chosen scope."))
    
    list(dt = dt, dt2 = dt2, dt_uni = dt_uni)
  }, ignoreInit = TRUE)
  
  
  # Shorting text
  short_label <- function(x, max_chars = 55, wrap_width = 18) {
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x)] <- ""
    x <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
    
    # wrap (si tens stringr; si no, queda tal qual)
    if (requireNamespace("stringr", quietly = TRUE)) {
      x <- stringr::str_wrap(x, width = wrap_width)
    }
    x
  }
  
  # ------------- Split genes
  
  split_genes_vec <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    
    # separadors habituals
    x <- gsub("\\s*[,;|/]\\s*", ";", x)     # , ; | /
    x <- gsub("\\s+and\\s+", ";", x, ignore.case = TRUE)
    x <- gsub("\\s*\\+\\s*", ";", x)
    x <- gsub("\\s*&\\s*",  ";", x)
    x <- gsub("\\s+-\\s+",  ";", x)         # "GENE1 - GENE2"
    
    parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[!is.na(parts) & nzchar(parts)]
    parts <- parts[parts != "-" & parts != "NA"]
    
    # NO unique() -> volem comptatges reals
    parts
  }
  
  # ---- gene picker UI (for scope == "gene") ----
  output$dt_gene_ui <- renderUI({
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits."))
    
    gene_col <- pick_col(dt, c(
      "MAPPED_GENE","mapped_gene",
      "REPORTED_GENE(S)","REPORTED_GENE","REPORTED_GENES",
      "GENE","gene","Gene","MAPPED_GENES"
    ))
    validate(need(!is.null(gene_col), "Gene column not found in current catalog."))
    
    genes <- split_genes_vec(dt[[gene_col]])
    genes <- sort(unique(toupper(genes)))
    validate(need(length(genes) > 0, "No genes found to select."))
    
    selectizeInput(
      "dt_gene", "Gene",
      choices = genes,
      selected = genes[1],
      multiple = FALSE,
      options = list(placeholder = "Type a gene symbol…")
    )
  })
  
  # -------------------------------
  # DISEASE table/plot (VERTICAL + zebra + proper universe column)
  # -------------------------------
  output$enrich_disease_table <- DT::renderDT({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    dis_col     <- pick_col(dt,     c("DISEASE","Disease","disease"))
    dis_col_uni <- pick_col(dt_uni, c("DISEASE","Disease","disease"))
    validate(need(!is.null(dis_col), "DISEASE column not found in current catalog."))
    validate(need(!is.null(dis_col_uni), "DISEASE column not found in reference universe."))
    
    tab <- fisher_enrich_terms_fast(dt2[[dis_col]], dt_uni[[dis_col]],
                                    alternative = "greater",
                                    min_a = 1, min_Tt = 1, max_terms = 2000)
    if (!nrow(tab)) {
      return(DT::datatable(data.frame(Message="No terms."), options=list(dom="t"), rownames=FALSE))
    }
    
    out <- tab %>%
      dplyr::transmute(
        term,
        in_subset  = paste0(a, "/", subset_n),
        n_universe = paste0(uni_term, "/", uni_N),
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format = "e", digits = 2),
        FDR = formatC(as.numeric(p_adj), format = "e", digits = 2)
      )
    
    DT::datatable(out, rownames=FALSE,               
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons = list(
                      list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                    ),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
  }, server = FALSE)
  
  output$enrich_disease_plot <- plotly::renderPlotly({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    dis_col     <- pick_col(dt,     c("DISEASE","Disease","disease"))
    dis_col_uni <- pick_col(dt_uni, c("DISEASE","Disease","disease"))
    validate(need(!is.null(dis_col), "DISEASE column not found."))
    validate(need(!is.null(dis_col_uni), "DISEASE column not found in universe."))
    
    tab <- fisher_enrich_terms_fast(
      dt2[[dis_col]], dt_uni[[dis_col]],
      alternative = "greater",
      min_a = 1, min_Tt = 1, max_terms = 2000
    )
    if (!nrow(tab)) return(NULL)
    
    top <- tab %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::mutate(
        counts     = a,
        is_zero    = is.finite(p_adj) & (p_adj == 0),
        logFDR_raw = suppressWarnings(-log10(p_adj)),
        # si FDR==0 o no-finit, assignem un valor alt per pintar "fort"
        logFDR     = dplyr::if_else(is.finite(logFDR_raw), logFDR_raw, max(logFDR_raw[is.finite(logFDR_raw)], na.rm = TRUE) + 1),
        term_short = short_label(term, max_chars = 55, wrap_width = 22),
        hover_txt  = paste0(
          "<b>", htmltools::htmlEscape(term), "</b>",
          "<br>in_subset: ", a, "/", subset_n,
          "<br>in_universe: ", uni_term, "/", uni_N,
          "<br>OR: ", sprintf("%.3f", as.numeric(or)),
          "<br>p: ", formatC(as.numeric(p), format = "e", digits = 2),
          "<br>FDR: ", ifelse(p_adj == 0, "0", formatC(as.numeric(p_adj), format = "e", digits = 2))
        )
      ) %>%
      dplyr::arrange(counts) %>%
      dplyr::mutate(term_short = factor(term_short, levels = term_short))
    
    
    
    p <- ggplot2::ggplot(top, ggplot2::aes(y = term_short, x = counts, text = hover_txt)) +
      
      # capa 1: no-zero → gradient
      ggplot2::geom_col(
        data = dplyr::filter(top, !is_zero),
        ggplot2::aes(fill = logFDR)
      ) +
      
      ggplot2::scale_fill_gradientn(
        colours = c("yellow", "orange", "red"),
        name = "-log10(FDR)"
      ) +
      
      # capa 2: zero → morat (no entra al gradient)
      ggplot2::geom_col(
        data = dplyr::filter(top, is_zero),
        fill = "#E34234"
      ) +
      
      ggplot2::labs(
        title = "Top enriched terms",
        y = NULL, x = "Counts in subset"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11),
        axis.text  = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 9),
        legend.text  = ggplot2::element_text(size = 8)
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        autosize = TRUE,
        margin = list(l = 210, r = 10, t = 35, b = 45),
        font   = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "enrich_disease_top10",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  # -------------------------------
  # TRAIT table/plot (FIX: uses trait_col; VERTICAL + zebra + proper universe column)
  # -------------------------------
  output$enrich_trait_table <- DT::renderDT({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    trait_col     <- pick_col(dt,     c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    trait_col_uni <- pick_col(dt_uni, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    validate(need(!is.null(trait_col), "MAPPED_TRAIT column not found in current catalog."))
    validate(need(!is.null(trait_col_uni), "MAPPED_TRAIT column not found in reference universe."))
    
    tab <- fisher_enrich_terms_fast(dt2[[trait_col]], dt_uni[[trait_col]],
                                    alternative = "greater",
                                    min_a = 1, min_Tt = 1, max_terms = 2000)
    if (!nrow(tab)) {
      return(DT::datatable(data.frame(Message="No terms."), options=list(dom="t"), rownames=FALSE))
    }
    
    out <- tab %>%
      dplyr::transmute(
        term,
        in_subset  = paste0(a, "/", subset_n),
        n_universe = paste0(uni_term, "/", uni_N),
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format = "e", digits = 2),
        FDR = formatC(as.numeric(p_adj), format = "e", digits = 2)
      )
    
    DT::datatable(out, rownames=FALSE,               
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons = list(
                      list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                    ),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
  }, server = FALSE)
  
  output$enrich_trait_plot <- plotly::renderPlotly({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    trait_col     <- pick_col(dt,     c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    trait_col_uni <- pick_col(dt_uni, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    validate(need(!is.null(trait_col), "MAPPED_TRAIT column not found."))
    validate(need(!is.null(trait_col_uni), "MAPPED_TRAIT column not found in universe."))
    
    tab <- fisher_enrich_terms_fast(
      dt2[[trait_col]], dt_uni[[trait_col]],
      alternative = "greater",
      min_a = 1, min_Tt = 1, max_terms = 2000
    )
    if (!nrow(tab)) return(NULL)
    
    top <- tab %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::mutate(
        counts     = a,
        is_zero    = is.finite(p_adj) & (p_adj == 0),
        logFDR_raw = suppressWarnings(-log10(p_adj)),
        # si FDR==0 o no-finit, assignem un valor alt per pintar "fort"
        logFDR     = dplyr::if_else(is.finite(logFDR_raw), logFDR_raw, max(logFDR_raw[is.finite(logFDR_raw)], na.rm = TRUE) + 1),
        term_short = short_label(term, max_chars = 55, wrap_width = 22),
        hover_txt  = paste0(
          "<b>", htmltools::htmlEscape(term), "</b>",
          "<br>in_subset: ", a, "/", subset_n,
          "<br>in_universe: ", uni_term, "/", uni_N,
          "<br>OR: ", sprintf("%.3f", as.numeric(or)),
          "<br>p: ", formatC(as.numeric(p), format = "e", digits = 2),
          "<br>FDR: ", ifelse(p_adj == 0, "0", formatC(as.numeric(p_adj), format = "e", digits = 2))
        )
      ) %>%
      dplyr::arrange(counts) %>%
      dplyr::mutate(term_short = factor(term_short, levels = term_short))
    
    
    
    p <- ggplot2::ggplot(top, ggplot2::aes(y = term_short, x = counts, text = hover_txt)) +
      
      # capa 1: no-zero → gradient
      ggplot2::geom_col(
        data = dplyr::filter(top, !is_zero),
        ggplot2::aes(fill = logFDR)
      ) +
      
      ggplot2::scale_fill_gradientn(
        colours = c("yellow", "orange", "red"),
        name = "-log10(FDR)"
      ) +
      
      # capa 2: zero → morat (no entra al gradient)
      ggplot2::geom_col(
        data = dplyr::filter(top, is_zero),
        fill = "#E34234"
      ) +
      
      ggplot2::labs(
        title = "Top enriched terms",
        y = NULL, x = "Counts in subset"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11),
        axis.text  = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 9),
        legend.text  = ggplot2::element_text(size = 8)
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        autosize = TRUE,
        margin = list(l = 210, r = 10, t = 35, b = 45),
        font   = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "enrich_trait_top10",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  
  # ---------------------------------------------------------------------------
  # GO / KEGG / GoSlim enrichment (NonSyn-style) driven by Catalog genes
  #   - Adds GoSlim tab support using:
  #       R/goslim_utils.R
  #       R/go_slim_generic.R   (expects GO_SLIM_GENERIC with GOID, ONTOLOGY, slim_term)
  #       R/goslim_generic.obo  (optional; not required if GO_SLIM_GENERIC exists)
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # UI helpers for cluster selector
  # ---------------------------------------------------------------------------
  
  output$func_cluster_ui <- renderUI({
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available."))
    
    chr_choices <- sort(unique(as.integer(cl$chr)))
    chr_labels  <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("func_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("func_cluster_id_ui")
    )
  })
  
  output$func_cluster_id_ui <- renderUI({
    cl <- clusters_val(); req(nrow(cl) > 0)
    req(input$func_chr)
    
    x <- cl %>% dplyr::filter(as.integer(chr) == as.integer(input$func_chr)) %>% dplyr::arrange(start, end)
    selectInput("func_cluster_id", "Cluster", choices = x$cluster_id, selected = x$cluster_id[1])
  })
  
  # ---------------------------------------------------------------------------
  # Gene symbols (global + scoped)
  # ---------------------------------------------------------------------------
  
  genes_all_symbols <- reactive({
    dt <- catalog_data()
    if (is.null(dt) || !nrow(dt)) return(character(0))
    gene_col <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    if (is.null(gene_col)) return(character(0))
    split_multi_unique(dt[[gene_col]])
  })
  
  genes_scope_symbols <- reactive({
    dt <- catalog_data()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No Catalog hits loaded."))
    
    gene_col <- pick_col(dt, c("MAPPED_GENE","mapped_gene","GENE","gene"))
    validate(need(!is.null(gene_col), "No MAPPED_GENE column found in Catalog hits."))
    
    if (identical(input$func_scope %||% "global", "cluster")) {
      validate(need("cluster_id" %in% names(dt), "Catalog hits missing cluster_id."))
      req(input$func_cluster_id)
      dt <- dt %>% dplyr::filter(cluster_id == input$func_cluster_id)
    }
    
    split_multi_unique(dt[[gene_col]])
  })
  
  # ---------------------------------------------------------------------------
  # SYMBOL -> ENTREZ mapping (cache)
  # ---------------------------------------------------------------------------
  
  gene_map_cache      <- reactiveVal(NULL)
  gene_map_cache_syms <- reactiveVal(NULL)
  
  gene_map_all <- reactive({
    syms <- genes_all_symbols()
    validate(need(length(syms) > 0, "No valid gene symbols found."))
    
    prev_syms <- gene_map_cache_syms()
    prev_map  <- gene_map_cache()
    
    if (!is.null(prev_syms) && identical(sort(prev_syms), sort(syms)) && !is.null(prev_map)) {
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
    
    validate(need(is.data.frame(m) && nrow(m) > 0, "Gene SYMBOL → ENTREZ mapping failed (empty)."))
    
    gene_map_cache_syms(sort(syms))
    gene_map_cache(m)
    m
  })
  
  entrez_scope <- reactive({
    m    <- gene_map_all()
    syms <- genes_scope_symbols()
    validate(need(length(syms) > 0, "No valid genes for the selected scope."))
    
    ids <- unique(m$ENTREZID[m$SYMBOL %in% syms])
    ids <- ids[!is.na(ids) & nzchar(ids)]
    validate(need(length(ids) > 0, "No ENTREZIDs for the selected scope (mapping empty)."))
    as.character(ids)
  })
  
  universe_entrez_dataset <- reactive({
    m <- gene_map_all()
    u <- unique(m$ENTREZID)
    u <- u[!is.na(u) & nzchar(u)]
    validate(need(length(u) > 0, "Universe (dataset) is empty after mapping."))
    as.character(u)
  })
  
  universe_entrez_for_scope <- reactive({
    scope <- input$func_scope %||% "global"
    bg    <- input$enrich_background %||% "dataset"
    
    # In GLOBAL mode: default background (OrgDb/KEGG default)
    if (identical(scope, "global")) return(NULL)
    
    # In CLUSTER mode: user decides
    if (identical(bg, "dataset")) return(universe_entrez_dataset())
    NULL
  })
  
  # ---------------------------------------------------------------------------
  # Background note (GO / KEGG / GoSlim)
  # ---------------------------------------------------------------------------
  
  output$enrich_bg_note <- renderUI({
    bg  <- input$enrich_background %||% "dataset"
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    bg_lbl <- if (identical(bg, "dataset")) "Dataset genes" else "Default background"
    
    # -----------------------------
    # Helpers (safe counts + scope)
    # -----------------------------
    sc <- input$func_scope %||% "global"
    sc_txt <- if (identical(sc, "cluster")) {
      paste0("Foreground: genes in <b>", htmltools::htmlEscape(input$func_cluster_id %||% ""), "</b>.")
    } else {
      "Foreground: genes in <b>global</b> (all extracted hits)."
    }
    
    safe_len <- function(expr) {
      tryCatch(as.integer(eval(expr)), error = function(e) NA_integer_)
    }
    
    # Foreground size (ENTREZ IDs used by GO/KEGG/GoSlim)
    n_fg <- safe_len(quote(length(entrez_scope())))
    
    # Dataset universe size (mapped dataset genes)
    n_uni_dataset <- safe_len(quote(length(universe_entrez_dataset())))
    
    # Default OrgDb universe size
    n_uni_orgdb <- tryCatch(
      as.integer(length(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID"))),
      error = function(e) NA_integer_
    )
    
    fg_cnt_txt <- if (is.finite(n_fg)) paste0("<br><b>Foreground size:</b> n=", n_fg, " genes.") else ""
    
    uni_txt <- function(default_name) {
      if (identical(bg, "dataset")) {
        if (is.finite(n_uni_dataset)) {
          paste0("<b>Background (universe):</b> dataset-mapped genes (N=", n_uni_dataset, ").")
        } else {
          "<b>Background (universe):</b> dataset-mapped genes."
        }
      } else {
        if (is.finite(n_uni_orgdb)) {
          paste0("<b>Background (universe):</b> ", default_name, " (N≈", n_uni_orgdb, ").")
        } else {
          paste0("<b>Background (universe):</b> ", default_name, ".")
        }
      }
    }
    
    gs_txt <- "<b>GSSize:</b> number of genes per term in the selected universe (filtered by minGSSize / maxGSSize)."
    
    # -----------------------------
    # Message per tab
    # -----------------------------
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
      
    } else if (identical(tab, "tab_enrich_goslim")) {
      
      tip <- "<br><br><b>Tip:</b> If you get no GoSlim terms, try setting <code>minGSSize</code> to <b>1–3</b> (GoSlim categories often have fewer hits)."
      
      paste0(
        "<b>GoSlim enrichment</b><br>",
        sc_txt,
        fg_cnt_txt, "<br>",
        if (identical(bg, "dataset")) {
          paste0("<b>Background (universe):</b> dataset-mapped genes",
                 if (is.finite(n_uni_dataset)) paste0(" (N=", n_uni_dataset, ").") else ".",
                 " (recommended in cluster mode).")
        } else {
          paste0("<b>Background (universe):</b> default OrgDb-based universe unless a dataset universe is provided.<br>",
                 "GoSlim categories come from the GO→GoSlim mapping (<i>goslim_generic</i>; TERM2GENE).")
        },
        "<br><b>Test unit:</b> genes (ENTREZID).<br>",
        gs_txt, "<br>",
        "<b>Method:</b> clusterProfiler::enricher with GoSlim TERM2GENE (GO ancestors mapped to <i>goslim_generic</i>).",
        tip
      )
      
    } else {
      paste0("<b>Background:</b> ", bg_lbl)
    }
    
    htmltools::HTML(paste0(
      "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
      msg, "</div>"
    ))
  })
  
  ##############################################################################
  # Load GOSlim
  # ------------------------------------------------------------
  # GoSlim helpers (FORCE global load, overwrite old defs)
  # ------------------------------------------------------------
  local({
    f1 <- file.path("R", "goslim_utils.R")
    f2 <- file.path("R", "go_slim_generic.R")
    f3 <- file.path("R", "goslim_generic.obo")
    
    if (file.exists(f1)) source(f1, local = FALSE)  # IMPORTANT (global)
    if (file.exists(f2)) source(f2, local = FALSE)  # IMPORTANT (global)
    
    if (!exists("GO_SLIM_GENERIC", inherits = TRUE) &&
        file.exists(f3) &&
        exists("load_goslim_generic_ids", inherits = TRUE)) {
      
      ids <- load_goslim_generic_ids(f3)
      
      GO_SLIM_GENERIC <- AnnotationDbi::select(
        GO.db::GO.db,
        keys    = ids,
        keytype = "GOID",
        columns = c("GOID", "TERM", "ONTOLOGY")
      ) |>
        dplyr::filter(!is.na(GOID), !is.na(TERM), !is.na(ONTOLOGY)) |>
        dplyr::transmute(GOID, ONTOLOGY, slim_term = TERM) |>
        dplyr::distinct(GOID, ONTOLOGY, slim_term)
      
      assign("GO_SLIM_GENERIC", GO_SLIM_GENERIC, envir = .GlobalEnv)
    }
  })
  
  
  # ---------------------------------------------------------------------------
  # Triggers (GO / KEGG / GoSlim) — avoid auto-run on load
  # ---------------------------------------------------------------------------
  
  go_trigger     <- reactiveVal(NULL)
  kegg_trigger   <- reactiveVal(NULL)
  goslim_trigger <- reactiveVal(NULL)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    if (identical(tab, "tab_enrich_go")) {
      go_trigger((go_trigger() %||% 0L) + 1L)
    } else if (identical(tab, "tab_enrich_kegg")) {
      kegg_trigger((kegg_trigger() %||% 0L) + 1L)
    } else if (identical(tab, "tab_enrich_goslim")) {
      goslim_trigger((goslim_trigger() %||% 0L) + 1L)
    }
  }, ignoreInit = TRUE)
  
  # ---------------------------------------------------------------------------
  # Shared scope label
  # ---------------------------------------------------------------------------
  
  scope_label <- function() {
    sc <- input$func_scope %||% "global"
    if (identical(sc, "cluster")) {
      cid <- input$func_cluster_id %||% ""
      if (nzchar(cid)) paste0("cluster ", cid) else "cluster"
    } else {
      "global"
    }
  }
  
  # ===========================================================================
  # GO enrichment (same as you had, kept)
  # ===========================================================================
  
  go_enrich_raw <- eventReactive(go_trigger(), {
    withProgress(message = "Running GO enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni  <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (OrgDb)" else paste("Universe:", length(uni)))
      
      ontos <- input$go_ontos
      if (is.null(ontos) || length(ontos) == 0 || all(is.na(ontos))) ontos <- c("BP","CC","MF")
      ontos <- intersect(c("BP","CC","MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one GO ontology (BP/CC/MF)."))
      
      minGS <- input$enrich_min_gs %||% 10
      maxGS <- input$enrich_max_gs %||% 500
      
      res <- purrr::map_dfr(ontos, function(ont) {
        incProgress(0.2 / length(ontos), detail = paste("Ontology:", ont))
        
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
        df <- as.data.frame(eg@result, stringsAsFactors = FALSE)
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
      dplyr::mutate(p.adjust = as.numeric(p.adjust), score = -log10(p.adjust)) %>%
      dplyr::filter(is.finite(p.adjust), p.adjust <= pcut) %>%
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
      return(DT::datatable(data.frame(Message="No GO terms passed the cutoff."),
                           options=list(dom="t"), rownames=FALSE))
    }
    
    out <- df %>%
      dplyr::select(Ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      dplyr::mutate(
        pvalue   = sprintf("%.3g", as.numeric(pvalue)),
        p.adjust = sprintf("%.3g", as.numeric(p.adjust))
      )
    
    DT::datatable(
      out, rownames = FALSE,
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX    = TRUE
      )
    )
  }, server = FALSE)
  
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
        term_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description),
        score      = -log10(pvalue)
      ) %>%
      dplyr::filter(Ontology %in% c("BP","CC","MF"), is.finite(pvalue), is.finite(p.adjust), nzchar(term_short))
    
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
  
  # ----------------------------
  # Helpers (format + wrap)
  # ----------------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  fmt_p <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(is.finite(x), sprintf("%.3g", x), NA_character_)
  }
  
  short_term <- function(x, max = 60) {
    x <- as.character(x)
    x <- stringr::str_squish(x)
    ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
  }
  
  # (opcional) si no tens shorten_term global
  if (!exists("shorten_term", mode = "function")) {
    shorten_term <- function(x, max = 55) {
      x <- as.character(x)
      x <- stringr::str_squish(x)
      x <- stringr::str_to_sentence(x)
      ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
    }
  }
  
  output$go_bar <- plotly::renderPlotly({
    
    dat0 <- go_top_df()
    
    if (!is.data.frame(dat0) || !nrow(dat0)) {
      return(plotly::plot_ly() %>% plotly::layout(title = "Run GO enrichment to see barplot."))
    }
    
    # --- Build Catalog-like input: term + n_genes + Ontology (+ p.adjust/pvalue) ---
    dat <- dat0 %>%
      dplyr::transmute(
        Ontology   = as.character(Ontology),
        go_term    = as.character(Description),
        term_short = if ("term_short" %in% names(dat0)) as.character(term_short) else as.character(Description),
        n_genes    = suppressWarnings(as.integer(Count)),
        p.adjust   = suppressWarnings(as.numeric(p.adjust)),
        pvalue     = suppressWarnings(as.numeric(pvalue))
      ) %>%
      dplyr::filter(!is.na(n_genes), n_genes > 0, !is.na(Ontology), nzchar(Ontology), nzchar(go_term))
    
    validate(need(nrow(dat) > 0, "GO: no data to plot yet."))
    
    top_n <- input$go_class_top %||% 15
    gap   <- 2
    
    # Helper: shorten/wrap terms
    if (!exists("shorten_term", mode = "function")) {
      shorten_term <- function(x, max = 55) {
        x <- as.character(x)
        x <- stringr::str_squish(x)
        x <- stringr::str_to_sentence(x)
        ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
      }
    }
    
    dat <- dat %>%
      dplyr::mutate(
        term_short = stringr::str_wrap(shorten_term(term_short, 55), width = 28)
      )
    
    ont_order <- c("BP","CC","MF")
    present   <- intersect(ont_order, unique(dat$Ontology))
    validate(need(length(present) > 0, "No BP/CC/MF data to plot."))
    
    dat <- dat %>%
      dplyr::filter(Ontology %in% present) %>%
      dplyr::mutate(Ontology = factor(Ontology, levels = ont_order))
    
    # Top-N per ontology: use n_genes (consistent with goslim_bar)
    dat <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::slice_max(order_by = n_genes, n = top_n, with_ties = FALSE) %>%
      dplyr::arrange(Ontology, dplyr::desc(n_genes)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    n_bp <- if ("BP" %in% present) max(dat$rank[dat$Ontology == "BP"], 0) else 0
    n_cc <- if ("CC" %in% present) max(dat$rank[dat$Ontology == "CC"], 0) else 0
    
    offsets <- c(BP = 0, CC = n_bp + gap, MF = n_bp + gap + n_cc + gap)
    dat$x <- dat$rank + offsets[as.character(dat$Ontology)]
    
    vlines <- c()
    if (all(c("BP","CC") %in% present)) vlines <- c(vlines, offsets["CC"] - gap/2)
    if (all(c("CC","MF") %in% present)) vlines <- c(vlines, offsets["MF"] - gap/2)
    
    centers <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::summarise(xmin=min(x), xmax=max(x), xmid=(xmin+xmax)/2, .groups="drop")
    
    ymax  <- max(dat$n_genes, na.rm = TRUE)
    ylab  <- -0.12 * ymax
    y_top <- ymax * 1.08
    
    sel_label <- if (exists("scope_label", mode = "function")) scope_label() else "global"
    pcut <- input$enrich_pcut %||% 0.05
    fb   <- isTRUE(attr(dat0, "fallback"))
    subttl <- if (fb) {
      paste0("No terms at FDR ≤ ", pcut, " — showing top terms ranked by FDR")
    } else {
      paste0("FDR ≤ ", pcut)
    }
    
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO term:</b> ", htmltools::htmlEscape(dat$go_term),
      "<br><b>Genes:</b> ", dat$n_genes,
      if ("p.adjust" %in% names(dat)) paste0("<br><b>p.adj:</b> ", fmt_p(dat$p.adjust)) else "",
      if ("pvalue"   %in% names(dat)) paste0("<br><b>p:</b> ", fmt_p(dat$pvalue)) else ""
    )
    
    p <- ggplot2::ggplot(dat, ggplot2::aes(x=x, y=n_genes, fill=Ontology, text=tooltip)) +
      ggplot2::geom_col(width=0.85, color="black", linewidth=0.25) +
      { if (length(vlines)) ggplot2::geom_vline(xintercept=vlines, linewidth=0.4) } +
      ggplot2::scale_x_continuous(breaks=dat$x, labels=dat$term_short, expand=c(0,0)) +
      ggplot2::coord_cartesian(ylim=c(ylab, ymax*1.20), clip="off") +
      { if ("BP" %in% present) ggplot2::annotate("text", x=centers$xmid[centers$Ontology=="BP"], y=y_top,
                                                 label="Biological Process", fontface="bold", size=3.6) } +
      { if ("CC" %in% present) ggplot2::annotate("text", x=centers$xmid[centers$Ontology=="CC"], y=y_top,
                                                 label="Cellular Component", fontface="bold", size=3.6) } +
      { if ("MF" %in% present) ggplot2::annotate("text", x=centers$xmid[centers$Ontology=="MF"], y=y_top,
                                                 label="Molecular Function", fontface="bold", size=3.6) } +
      ggplot2::labs(
        x=NULL, y="Number of genes",
        title=paste0("GO enrichment — ", sel_label),
        subtitle=subttl
      ) +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size=11, face="bold"),
        plot.subtitle= ggplot2::element_text(size=9),
        axis.text.x  = ggplot2::element_text(angle=90, size=7, vjust=0.5, hjust=1),
        legend.position="none",
        plot.margin  = grid::unit(c(28,10,35,10), "pt")
      ) +
      ggplot2::scale_fill_manual(values=c(BP="darkgreen", CC="orange", MF="darkblue"), guide="none")
    
    plotly::ggplotly(p, tooltip="text") %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "enrich_GO",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  
  # ===========================================================================
  # KEGG enrichment — FIXED mapping so FDR=1 still returns something (when possible)
  #   - DO NOT set keyType="kegg" here (that triggers "No gene can be mapped" with Entrez IDs)
  # ===========================================================================
  
  kegg_enrich_raw <- eventReactive(kegg_trigger(), {
    withProgress(message = "Running KEGG enrichment…", value = 0, {
      
      gene <- as.character(entrez_scope())
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni  <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (KEGG)" else paste("Universe:", length(uni)))
      
      minGS <- input$enrich_min_gs %||% 10
      maxGS <- input$enrich_max_gs %||% 500
      
      ek <- tryCatch({
        suppressMessages(
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
      }, error = function(e) NULL)
      
      if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(tibble::tibble())
      as.data.frame(ek@result, stringsAsFactors = FALSE)
    })
  }, ignoreInit = TRUE)
  
  kegg_top_df <- reactive({
    df <- kegg_enrich_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$enrich_kegg_top %||% 15
    
    df0 <- df %>%
      dplyr::mutate(
        pvalue     = suppressWarnings(as.numeric(pvalue)),
        p.adjust   = suppressWarnings(as.numeric(p.adjust)),
        term_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description),
        score      = -log10(pvalue)
      ) %>%
      dplyr::filter(is.finite(pvalue), is.finite(p.adjust), nzchar(term_short))
    
    if (!nrow(df0)) return(tibble::tibble())
    
    sig <- df0 %>% dplyr::filter(p.adjust <= pcut) %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = topn)
    if (nrow(sig) > 0) {
      attr(sig, "fallback") <- FALSE
      return(sig)
    }
    
    fb <- df0 %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = topn)
    attr(fb, "fallback") <- TRUE
    fb
  })
  
  output$kegg_bar <- renderPlot({
    df <- kegg_top_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    fb   <- isTRUE(attr(df, "fallback"))
    pcut <- input$enrich_pcut %||% 0.05
    
    ttl <- paste0(
      "KEGG enrichment (", scope_label(), ")",
      if (fb) paste0(" — No pathways at FDR ≤ ", pcut, " → showing top by FDR") else paste0(" — FDR ≤ ", pcut)
    )
    
    df <- df %>% dplyr::mutate(score = -log10(as.numeric(pvalue)))
    fill_id <- if ("ID" %in% names(df)) as.character(df$ID) else as.character(df$Description)
    df <- df %>% dplyr::mutate(fill_id = factor(fill_id, levels = fill_id))
    
    n   <- nlevels(df$fill_id)
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
    names(pal) <- levels(df$fill_id)
    
    ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Description, score), y = score, fill = fill_id)) +
      ggplot2::geom_col(color = "black", linewidth = 0.25) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::labs(title = ttl, x = NULL, y = "-log10(pvalue)") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }, height = function() {
    n <- tryCatch(nrow(kegg_top_df()), error = function(e) 0)
    max(420, 28 * n)
  })
  
  output$kegg_table <- DT::renderDT({
    df <- kegg_top_df()
    if (!nrow(df)) {
      return(DT::datatable(data.frame(Message="No KEGG pathways to display."),
                           options=list(dom="t"), rownames=FALSE))
    }
    
    out <- df %>%
      dplyr::mutate(
        KEGG = if ("ID" %in% names(df))
          sprintf("<a href='https://www.kegg.jp/pathway/%s' target='_blank'>%s</a>", ID, ID)
        else NA_character_
      ) %>%
      dplyr::transmute(
        KEGG,
        Description,
        GeneRatio,
        BgRatio,
        Count,
        pvalue   = sprintf("%.3g", as.numeric(pvalue)),
        p.adjust = sprintf("%.3g", as.numeric(p.adjust))
      )
    
    DT::datatable(
      out,
      escape     = FALSE,
      rownames   = FALSE,
      filter     = "top",
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 15,
        scrollX    = TRUE
      )
    )
  }, server = FALSE)
  
  ##############################################################################
  ###############################################################################
  # ==========================================================
  # GoSlim (OrgDb-based canonical) — Catalog genes
  #   - TERM2GENE / TERM2NAME / universe built from OrgDb + GO.db
  #   - Consistent with "OrgDb-based universe"
  #   - goslim_bar kept as CLASSIFICATION by gene counts
  # ==========================================================
  
  # ---- helper cache ----
  goslim_orgdb_cache <- reactiveVal(list())
  
  # ==========================================================
  # Build TERM2GENE / TERM2NAME / universe from OrgDb + GO.db
  # ==========================================================
  goslim_build_orgdb_sets_for_onto <- function(ontology = c("BP", "CC", "MF")) {
    ontology <- match.arg(ontology)
    
    validate(need(requireNamespace("AnnotationDbi", quietly = TRUE), "AnnotationDbi package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    validate(need(requireNamespace("GO.db", quietly = TRUE), "GO.db package is required."))
    validate(need(exists("GO_SLIM_GENERIC", inherits = TRUE), "GO_SLIM_GENERIC not loaded."))
    
    key <- paste0("orgdb_goslim_", ontology)
    cache <- goslim_orgdb_cache()
    if (!is.null(cache[[key]])) return(cache[[key]])
    
    slim_df <- GO_SLIM_GENERIC %>%
      dplyr::filter(ONTOLOGY == ontology) %>%
      dplyr::distinct(GOID, slim_term)
    
    slim_ids <- unique(as.character(slim_df$GOID))
    slim_ids <- slim_ids[!is.na(slim_ids) & nzchar(slim_ids)]
    validate(need(length(slim_ids) > 0, paste0("No GoSlim IDs loaded for ", ontology)))
    
    eg_keys <- tryCatch(
      AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"),
      error = function(e) NULL
    )
    eg_keys <- unique(as.character(eg_keys))
    eg_keys <- eg_keys[!is.na(eg_keys) & nzchar(eg_keys)]
    validate(need(length(eg_keys) > 0, "No ENTREZ keys found in org.Hs.eg.db"))
    
    org_map <- tryCatch({
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = eg_keys,
        keytype = "ENTREZID",
        columns = c("GO", "ONTOLOGY")
      )
    }, error = function(e) NULL)
    
    validate(need(is.data.frame(org_map) && nrow(org_map) > 0, "Could not retrieve GO mapping from org.Hs.eg.db"))
    validate(need(all(c("ENTREZID", "GO", "ONTOLOGY") %in% names(org_map)),
                  "OrgDb GO mapping did not return ENTREZID/GO/ONTOLOGY columns"))
    
    org_map <- org_map %>%
      dplyr::transmute(
        ENTREZID = as.character(ENTREZID),
        GOID     = as.character(GO),
        ONTOLOGY = as.character(ONTOLOGY)
      ) %>%
      dplyr::filter(!is.na(ENTREZID), nzchar(ENTREZID)) %>%
      dplyr::filter(!is.na(GOID), nzchar(GOID)) %>%
      dplyr::filter(ONTOLOGY == ontology) %>%
      dplyr::distinct(ENTREZID, GOID)
    
    validate(need(nrow(org_map) > 0, paste0("No OrgDb GO mapping found for ontology ", ontology)))
    
    anc_obj <- switch(
      ontology,
      BP = GO.db::GOBPANCESTOR,
      CC = GO.db::GOCCANCESTOR,
      MF = GO.db::GOMFANCESTOR
    )
    
    goids <- unique(org_map$GOID)
    
    go2slim <- lapply(goids, function(goid) {
      anc <- tryCatch(anc_obj[[goid]], error = function(e) NULL)
      anc <- unique(c(goid, as.character(anc %||% character(0))))
      intersect(anc, slim_ids)
    })
    names(go2slim) <- goids
    
    go2slim_df <- tibble::tibble(
      GOID = goids,
      slim_goid = I(go2slim)
    ) %>%
      tidyr::unnest(cols = c(slim_goid)) %>%
      dplyr::filter(!is.na(slim_goid), nzchar(slim_goid)) %>%
      dplyr::left_join(
        slim_df,
        by = c("slim_goid" = "GOID")
      ) %>%
      dplyr::filter(!is.na(slim_term), nzchar(slim_term)) %>%
      dplyr::distinct(GOID, slim_goid, slim_term)
    
    validate(need(nrow(go2slim_df) > 0, paste0("No GoSlim mapping produced for ", ontology)))
    
    term2gene <- org_map %>%
      dplyr::inner_join(
        go2slim_df,
        by = "GOID",
        relationship = "many-to-many"
      ) %>%
      dplyr::transmute(
        term = as.character(slim_goid),
        gene = as.character(ENTREZID)
      ) %>%
      dplyr::filter(!is.na(term), nzchar(term), !is.na(gene), nzchar(gene)) %>%
      dplyr::distinct()
    
    validate(need(nrow(term2gene) > 0, paste0("TERM2GENE empty for ", ontology)))
    
    term2name <- go2slim_df %>%
      dplyr::transmute(
        term = as.character(slim_goid),
        name = as.character(slim_term)
      ) %>%
      dplyr::filter(!is.na(term), nzchar(term), !is.na(name), nzchar(name)) %>%
      dplyr::distinct()
    
    universe <- unique(term2gene$gene)
    universe <- universe[!is.na(universe) & nzchar(universe)]
    validate(need(length(universe) > 0, paste0("Universe empty for ", ontology)))
    
    out <- list(
      term2gene = term2gene,
      term2name = term2name,
      universe  = universe
    )
    
    cache[[key]] <- out
    goslim_orgdb_cache(cache)
    out
  }
  
  
  goslim_orgdb_sets <- reactive({
    withProgress(message = "Preparing GO slim reference…", value = 0, {
      onts <- c("BP", "CC", "MF")
      out <- vector("list", length(onts))
      names(out) <- onts
      
      for (i in seq_along(onts)) {
        ont <- onts[i]
        incProgress(1 / length(onts), detail = paste("Reference:", ont))
        out[[ont]] <- goslim_build_orgdb_sets_for_onto(ont)
      }
      
      out
    })
  }) 
  
  # ------------------------------------------------------------
  # GoSlim ENRICHMENT (OrgDb-based canonical)
  # ------------------------------------------------------------
  goslim_enrich_raw <- eventReactive(goslim_trigger(), {
    validate(need(exists("GO_SLIM_GENERIC", inherits = TRUE), "GO_SLIM_GENERIC not loaded."))
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db required."))
    validate(need(requireNamespace("GO.db", quietly = TRUE), "GO.db required."))
    
    withProgress(message = "Running GO slim enrichment…", value = 0, {
      
      gene <- as.character(entrez_scope())
      gene <- unique(gene[!is.na(gene) & nzchar(gene)])
      validate(need(length(gene) > 0, "No genes for GoSlim (scope empty)."))
      incProgress(0.20, detail = paste("Genes:", length(gene)))
      
      ontos <- input$go_ontos
      if (is.null(ontos) || !length(ontos)) ontos <- c("BP", "CC", "MF")
      ontos <- intersect(c("BP", "CC", "MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one GO ontology (BP/CC/MF)."))
      
      minGS_in <- input$enrich_min_gs %||% 10
      maxGS    <- input$enrich_max_gs %||% 500
      
      gs_all <- goslim_orgdb_sets()
      incProgress(0.10, detail = "Reference ready")
      
      res <- purrr::map_dfr(ontos, function(ont) {
        incProgress(0.60 / length(ontos), detail = paste("Ontology:", ont))
        
        gs <- gs_all[[ont]]
        
        universe <- unique(gs$universe)
        gene_use <- intersect(unique(gene), universe)
        
        message(
          "[GoSlim][Catalog][", ont, "] ",
          "selected_genes=", length(gene),
          " | universe=", length(universe),
          " | gene_use=", length(gene_use),
          " | term2gene_rows=", nrow(gs$term2gene),
          " | n_terms=", dplyr::n_distinct(gs$term2gene$term)
        )
        
        if (!length(gene_use)) return(NULL)
        
        minGS <- max(1L, min(as.integer(minGS_in), length(gene_use)))
        
        eg <- suppressMessages(
          clusterProfiler::enricher(
            gene          = gene_use,
            universe      = universe,
            TERM2GENE     = gs$term2gene,
            TERM2NAME     = gs$term2name,
            pAdjustMethod = "BH",
            pvalueCutoff  = 1,
            qvalueCutoff  = 1,
            minGSSize     = minGS,
            maxGSSize     = maxGS
          )
        )
        
        if (is.null(eg) || is.null(eg@result) || !nrow(eg@result)) return(NULL)
        
        df <- as.data.frame(eg@result, stringsAsFactors = FALSE)
        if (!nrow(df)) return(NULL)
        df$Ontology <- ont
        df
      })
      
      if (!is.data.frame(res) || !nrow(res)) return(tibble::tibble())
      res
    })
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # GoSlim CLASSIFICATION for goslim_bar / goslim_table
  # Output: Ontology, slim_term, n_genes, n_terms
  # ------------------------------------------------------------
  
  scope_label_static <- reactive({
    mode <- input$go_mode %||% "total"
    
    if (identical(mode, "total")) return("Total")
    
    if (identical(mode, "chr")) {
      chr_sel <- as.character(input$go_chr %||% "")
      if (!nzchar(chr_sel)) return("Chromosome")
      return(paste0("chr", chr_sel))
    }
    
    cl_id <- as.character(input$go_cluster_id %||% "")
    if (!nzchar(cl_id)) return("Cluster")
    cl_id
  })
  
  go_slim_class_df_selected <- reactive({
    df <- goslim_enrich_tbl()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    out <- df %>%
      dplyr::transmute(
        Ontology  = as.character(Ontology),
        slim_term = as.character(Description),
        n_genes   = suppressWarnings(as.integer(Count)),
        n_terms   = NA_integer_
      ) %>%
      dplyr::filter(!is.na(n_genes), nzchar(slim_term)) %>%
      dplyr::arrange(Ontology, dplyr::desc(n_genes), slim_term)
    
    attr(out, "label") <- scope_label_static()
    out
  })
  
  
  
  # ==========================================================
  # GoSlim PLOT (classification by counts; titles clarified)
  # ==========================================================
  output$goslim_bar <- renderPlotly({
    
    dat <- go_slim_class_df_selected()
    validate(need(is.data.frame(dat) && nrow(dat) > 0, "GO slim: no data to plot yet."))
    
    top_n <- input$go_class_top %||% 15
    gap   <- 2
    
    validate(need(all(c("slim_term","n_genes","Ontology") %in% names(dat)),
                  "go_slim_class_df_selected() must return: slim_term, n_genes, Ontology"))
    
    if (!exists("shorten_term", mode = "function")) {
      shorten_term <- function(x, max = 55) {
        x <- as.character(x)
        x <- stringr::str_squish(x)
        x <- stringr::str_to_sentence(x)
        ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
      }
    }
    
    dat <- dat %>%
      dplyr::mutate(
        Ontology   = as.character(Ontology),
        term_short = stringr::str_wrap(shorten_term(slim_term, 55), width = 28)
      )
    
    ont_order <- c("BP","CC","MF")
    present   <- intersect(ont_order, unique(dat$Ontology))
    validate(need(length(present) > 0, "No BP/CC/MF data to plot."))
    
    dat <- dat %>%
      dplyr::filter(Ontology %in% present) %>%
      dplyr::mutate(Ontology = factor(Ontology, levels = ont_order))
    
    dat <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::slice_max(order_by = n_genes, n = top_n, with_ties = FALSE) %>%
      dplyr::arrange(Ontology, dplyr::desc(n_genes)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    n_bp <- if ("BP" %in% present) max(dat$rank[dat$Ontology == "BP"], 0) else 0
    n_cc <- if ("CC" %in% present) max(dat$rank[dat$Ontology == "CC"], 0) else 0
    
    offsets <- c(
      BP = 0,
      CC = n_bp + gap,
      MF = n_bp + gap + n_cc + gap
    )
    
    dat$x <- dat$rank + offsets[as.character(dat$Ontology)]
    
    vlines <- c()
    if (all(c("BP","CC") %in% present)) vlines <- c(vlines, offsets["CC"] - gap / 2)
    if (all(c("CC","MF") %in% present)) vlines <- c(vlines, offsets["MF"] - gap / 2)
    
    centers <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::summarise(
        xmin = min(x), xmax = max(x),
        xmid = (xmin + xmax) / 2,
        .groups = "drop"
      )
    
    ymax  <- max(dat$n_genes, na.rm = TRUE)
    ylab  <- -0.12 * ymax
    y_top <- ymax * 1.08
    
    sel_label <- attr(dat, "label") %||% scope_label_static()
    
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO slim term:</b> ", htmltools::htmlEscape(as.character(dat$slim_term)),
      "<br><b>Selected genes in term:</b> ", dat$n_genes
    )
    
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = x, y = n_genes, fill = Ontology, text = tooltip)
    ) +
      ggplot2::geom_col(width = 0.85, color = "black", linewidth = 0.25) +
      { if (length(vlines)) ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4) } +
      ggplot2::scale_x_continuous(
        breaks = dat$x,
        labels = dat$term_short,
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = c(ylab, ymax * 1.20), clip = "off") +
      { if ("BP" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology == "BP"], y = y_top,
                          label = "Biological Process", fontface = "bold", size = 3.6) } +
      { if ("CC" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology == "CC"], y = y_top,
                          label = "Cellular Component", fontface = "bold", size = 3.6) } +
      { if ("MF" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology == "MF"], y = y_top,
                          label = "Molecular Function", fontface = "bold", size = 3.6) } +
      ggplot2::labs(
        x = NULL,
        y = "Selected genes per GO slim term",
        title = paste0("GO slim classification by gene counts — ", sel_label)
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 11, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 9),
        axis.title.y = ggplot2::element_text(size = 9),
        axis.text.x  = ggplot2::element_text(angle = 90, size = 7, vjust = 0.5, hjust = 1),
        axis.text.y  = ggplot2::element_text(size = 8),
        legend.position = "none",
        plot.margin = grid::unit(c(28, 10, 35, 10), "pt")
      ) +
      ggplot2::scale_fill_manual(
        values = c(BP = "darkgreen", CC = "orange", MF = "darkblue"),
        guide = "none"
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "goslim_classification_counts",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  outputOptions(output, "goslim_bar", suspendWhenHidden = FALSE)
  
  # ==========================================================
  # ==========================================================
  
  goslim_enrich_tbl <- reactive({
    df <- goslim_enrich_raw()
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
    
    if (!nrow(df2)) {
      df2 <- df %>%
        dplyr::mutate(
          Ontology = as.character(Ontology),
          pvalue   = suppressWarnings(as.numeric(pvalue)),
          p.adjust = suppressWarnings(as.numeric(p.adjust))
        ) %>%
        dplyr::filter(is.finite(p.adjust), is.finite(pvalue)) %>%
        dplyr::group_by(Ontology) %>%
        dplyr::arrange(p.adjust, .by_group = TRUE) %>%
        dplyr::slice_head(n = topn) %>%
        dplyr::ungroup()
    }
    
    df2
  })
  
  # ==========================================================
  # GoSlim TABLE
  # ==========================================================
  output$goslim_table <- DT::renderDT({
    df <- goslim_enrich_tbl()
    
    if (!is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(
        data.frame(Message = "No GO slim terms passed the cutoff."),
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
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        )
      )
    )
  }, server = FALSE)
  
  
  
  ##############################################################################
  # ============================================================
  # CATALOG · Gene table (uses SAME inputs as enrichment)
  #   - genes_scope_symbols()  (scope)
  #   - genes_all_symbols()    (dataset)
  #   - gene_map_all()         (SYMBOL -> ENTREZ)
  # ============================================================
  # ---- cluster picker UI for genes table ----
  output$genes_cluster_ui <- renderUI({
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available."))
    
    chr_choices <- sort(unique(as.integer(cl$chr)))
    chr_labels  <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("genes_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("genes_cluster_id_ui")
    )
  })
  
  output$genes_cluster_id_ui <- renderUI({
    cl <- clusters_val(); req(is.data.frame(cl), nrow(cl) > 0)
    req(input$genes_chr)
    
    x <- cl %>%
      dplyr::filter(as.integer(chr) == as.integer(input$genes_chr)) %>%
      dplyr::arrange(start, end)
    
    validate(need(nrow(x) > 0, "No clusters for this chromosome."))
    selectInput("genes_cluster_id", "Cluster", choices = x$cluster_id, selected = x$cluster_id[1])
  })
  
  output$genes_table <- DT::renderDT({
    dt <- catalog_hits_df()
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No Catalog hits loaded."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # -----------------------------
    # Apply scope filter (global / cluster)
    # -----------------------------
    scope <- input$genes_scope %||% "global"
    
    if (!"cluster_id" %in% names(dt)) {
      # if cluster_id missing, only global makes sense (but still show table)
      scope <- "all"
    }
    
    caption_txt <- "Catalog genes"
    
    if (identical(scope, "global")) {
      # Global = clustered hits only (consistent with enrichment global)
      dt <- dt %>% dplyr::filter(!is.na(cluster_id), nzchar(cluster_id))
      caption_txt <- "Catalog genes (Global: clustered hits)"
    } else if (identical(scope, "cluster")) {
      req(input$genes_cluster_id)
      dt <- dt %>% dplyr::filter(cluster_id == input$genes_cluster_id)
      caption_txt <- paste0("Catalog genes (Cluster: ", input$genes_cluster_id, ")")
    } else {
      caption_txt <- "Catalog genes"
    }
    
    if (!nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No rows in this scope."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # -----------------------------
    # Helpers: pick first existing col
    # -----------------------------
    pick1 <- function(df, cand) {
      nm <- intersect(cand, names(df))
      if (length(nm)) nm[1] else NULL
    }
    
    # -----------------------------
    # Canonicalize required columns (create if missing)
    # -----------------------------
    
    # gene
    if (!"gene" %in% names(dt)) {
      gcol <- pick1(dt, c("MAPPED_GENE","mapped_gene","GENE","Gene","gene_name","nearest_gene"))
      if (!is.null(gcol)) dt$gene <- as.character(dt[[gcol]])
    } else {
      dt$gene <- as.character(dt$gene)
    }
    
    # SYMBOL
    if (!"SYMBOL" %in% names(dt)) {
      scol <- pick1(dt, c("symbol","SYMBOL","hgnc_symbol"))
      if (!is.null(scol)) {
        dt$SYMBOL <- as.character(dt[[scol]])
      } else if ("gene" %in% names(dt)) {
        dt$SYMBOL <- as.character(dt$gene)
        dt$SYMBOL <- sub("\\s*[;,|].*$", "", dt$SYMBOL)
      }
    } else {
      dt$SYMBOL <- as.character(dt$SYMBOL)
    }
    
    # ENTREZID
    if (!"ENTREZID" %in% names(dt)) {
      ecol <- pick1(dt, c("entrezid","ENTREZID","Entrez","ENTREZ","gene_id"))
      if (!is.null(ecol)) dt$ENTREZID <- as.character(dt[[ecol]])
    } else {
      dt$ENTREZID <- as.character(dt$ENTREZID)
    }
    
    # freq
    if (!"freq" %in% names(dt)) {
      fcol <- pick1(dt, c("freq","FREQ","frequency","FREQUENCY","allele_freq","AF"))
      if (!is.null(fcol)) dt$freq <- suppressWarnings(as.numeric(dt[[fcol]]))
    } else {
      dt$freq <- suppressWarnings(as.numeric(dt$freq))
    }
    
    # -----------------------------
    # Select columns
    # -----------------------------
    keep <- c("gene","SYMBOL","ENTREZID","PUBMEDID","DISEASE","MAPPED_TRAIT","freq")
    keep <- keep[keep %in% names(dt)]
    out  <- dt[, keep, drop = FALSE]
    
    # -----------------------------
    # Links + formatting
    # -----------------------------
    
    # GeneCards link on gene
    if ("gene" %in% names(out)) {
      g <- as.character(out$gene)
      ok <- !is.na(g) & nzchar(g)
      out$gene <- ifelse(
        ok,
        sprintf(
          "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank'>%s</a>",
          URLencode(g, reserved = TRUE),
          htmltools::htmlEscape(g)
        ),
        NA_character_
      )
    }
    
    # PubMed link
    if ("PUBMEDID" %in% names(out)) {
      pm <- as.character(out$PUBMEDID)
      ok <- !is.na(pm) & nzchar(pm)
      out$PUBMEDID <- ifelse(
        ok,
        sprintf(
          "<a href='https://pubmed.ncbi.nlm.nih.gov/%s/' target='_blank'>%s</a>",
          URLencode(pm, reserved = TRUE),
          htmltools::htmlEscape(pm)
        ),
        NA_character_
      )
    }
    
    # freq formatting
    if ("freq" %in% names(out)) {
      x <- suppressWarnings(as.numeric(out$freq))
      out$freq <- ifelse(is.finite(x), sprintf("%.3f", x), NA_character_)
    }
    
    DT::datatable(
      out,
      rownames = FALSE,
      escape   = FALSE,
      caption  = htmltools::tags$caption(
        style = "caption-side: top; text-align: left; font-weight: 600; font-size: 18px; padding: 6px 0;",
        caption_txt
      ),
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list( extend = "csv",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list( extend = "excel",exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        ),
        pageLength = 20,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  
  # --- DEBUG: print column names to console (catalog_hits_df + catalog_data) ---
  observeEvent(catalog_hits_df(), {
    dt <- catalog_hits_df()
    cat("\n====================\n")
    cat("[catalog_hits_df]  nrow=", nrow(dt), " ncol=", ncol(dt), "\n", sep = "")
    cat("[catalog_hits_df]  cols: ", paste(names(dt), collapse = ", "), "\n", sep = "")
    cat("====================\n")
  }, ignoreInit = FALSE)
  
  observeEvent(catalog_data(), {
    dt <- catalog_data()
    cat("\n====================\n")
    cat("[catalog_data]     nrow=", nrow(dt), " ncol=", ncol(dt), "\n", sep = "")
    cat("[catalog_data]     cols: ", paste(names(dt), collapse = ", "), "\n", sep = "")
    cat("====================\n")
  }, ignoreInit = FALSE)
  
  
  ##############################################################################
  # -------------------------
  # Cluster intervals source (must exist)
  #   -> assumes you already have cluster_intervals_src_ld() reactive defined
  # -------------------------
  # If you *don't* have it defined, keep your working version (the one that derives
  # from clusters_val() or catalog_data()).
  
  # -------------------------
  # Hits-of-interest (GWAS + Catalog) for the selected cluster
  # -------------------------
  
  cluster_hits_cl <- reactive({
    df_int  <- cluster_intervals_src_ld()
    chr_sel <- as_int_safe(input$ld_chr)
    cl_sel  <- as.character(input$ld_cluster_id %||% "")
    
    validate(need(is.finite(chr_sel) && nzchar(cl_sel), "Select chr + cluster."))
    
    one <- df_int %>% dplyr::filter(CHR == chr_sel, cluster_id == cl_sel)
    validate(need(nrow(one) == 1, "Cluster interval not found."))
    
    st <- as_int_safe(one$start_bp[1])
    en <- as_int_safe(one$end_bp[1])
    
    # GWAS hits (source = GWAS)
    g <- gwas_df()
    g2 <- g %>%
      dplyr::filter(as_int_safe(CHR) == chr_sel, BP >= st, BP <= en) %>%
      dplyr::transmute(
        SNP    = as.character(snp),
        source = "GWAS",
        CHR    = chr_sel,
        BP     = as_int_safe(BP)
      )
    
    # Catalog hits from normalized output (source = Catalog)
    cdt <- tryCatch(catalog_data(), error = function(e) NULL)
    c2 <- tibble::tibble()
    
    if (is.data.frame(cdt) && nrow(cdt) > 0 &&
        all(c("cluster_id","CHR","POS","rsid") %in% names(cdt))) {
      
      c2 <- cdt %>%
        dplyr::filter(cluster_id == cl_sel) %>%
        dplyr::transmute(
          SNP    = as.character(rsid),
          source = "Catalog",
          CHR    = chr_to_num_plink(CHR),
          BP     = as_int_safe(POS)
        ) %>%
        dplyr::filter(as_int_safe(CHR) == chr_sel, BP >= st, BP <= en, is.finite(BP), nzchar(SNP))
    }
    
    dplyr::bind_rows(g2, c2) %>%
      dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
      dplyr::distinct(SNP, .keep_all = TRUE) %>%
      dplyr::arrange(BP)
  })
  
  ############### save for LD 
  # =============================
  # Save cluster.csv + candidate.csv (Catalog Inspector)
  #  - clusters: clusters_val()
  #  - catalog hits: catalog_data()
  #  - output: workdir/exports/
  # =============================
  
  norm_chr_int <- function(x) {
    x <- toupper(as.character(x))
    x <- gsub("^CHR", "", x)
    x[x == "X"]  <- "23"
    x[x == "Y"]  <- "24"
    x[x %in% c("M","MT")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pick_first_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  ##############################################################################
  # ---- LD virtual inputs (session cache) ----
  rv_ld_inputs <- reactiveValues(
    dir = NULL,
    cluster_path = NULL,
    candidate_path = NULL,
    cluster_df = NULL,
    candidate_df = NULL
  )
  
  # folder de sessió per guardar CSVs (per debug / reutilització)
  ld_inputs_dir <- reactive({
    if (!is.null(rv_ld_inputs$dir) && dir.exists(rv_ld_inputs$dir)) return(rv_ld_inputs$dir)
    d <- file.path(tempdir(), "catalog_ld_inputs")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    rv_ld_inputs$dir <- d
    d
  })
  
  # helper: construeix els mateixos dataframes que al ZIP
  build_catalog_ld_inputs <- function() {
    
    # -----------------------------
    # Prepare clusters
    # -----------------------------
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available (build clusters first)."))
    
    cl <- as.data.frame(cl)
    if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
    
    validate(need(all(c("cluster_id","chr","start","end") %in% names(cl)),
                  "clusters_val() must contain cluster_id, chr, start, end (or compatible names)."))
    
    cl$chr   <- norm_chr_int(cl$chr)
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    validate(need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
    
    # -----------------------------
    # Prepare Catalog hits
    # -----------------------------
    cd <- catalog_data()
    validate(need(is.data.frame(cd) && nrow(cd) > 0, "No Catalog hits available (run Step 4 first)."))
    cd <- as.data.frame(cd)
    
    chr_col   <- pick_first_col(cd, c("chr","CHR","chrom","chromosome","CHR_ID"))
    start_col <- pick_first_col(cd, c("pos_ini","start","START","bp_start","POS","BP","pos"))
    end_col   <- pick_first_col(cd, c("pos_end","end","END","bp_end"))
    validate(need(!is.null(chr_col), "catalog_data() has no chromosome column (chr/CHR/...)."))
    validate(need(!is.null(start_col) || !is.null(end_col),
                  "catalog_data() has no position column (POS/BP/pos/start)."))
    if (is.null(end_col)) end_col <- start_col
    
    id_col <- pick_first_col(cd, c("id_hit","snp","SNPS","rsid","RSID","variant","ID","code","CODE"))
    if (is.null(id_col)) {
      id_col <- ".id_tmp"
      cd[[id_col]] <- paste0("Catalog_chr", cd[[chr_col]], ":", cd[[start_col]], "-", cd[[end_col]])
    }
    
    b <- data.frame(
      chr       = norm_chr_int(cd[[chr_col]]),
      pos_ini   = suppressWarnings(as.integer(cd[[start_col]])),
      pos_end   = suppressWarnings(as.integer(cd[[end_col]])),
      id_hit    = as.character(cd[[id_col]]),
      classe    = "Catalog_hit",
      stringsAsFactors = FALSE
    )
    b <- b[is.finite(b$chr) & is.finite(b$pos_ini) & is.finite(b$pos_end) & b$pos_end >= b$pos_ini, , drop = FALSE]
    validate(need(nrow(b) > 0, "All Catalog hits became invalid after coercion (chr/start/end)."))
    
    # -----------------------------
    # Keep clusters that have catalog hits
    # -----------------------------
    if ("cluster_id" %in% names(cd) && any(nzchar(as.character(cd$cluster_id)))) {
      keep_ids <- unique(as.character(cd$cluster_id))
      cl_sig <- cl[as.character(cl$cluster_id) %in% keep_ids, , drop = FALSE]
    } else {
      library(data.table)
      cl_dt <- data.table::as.data.table(cl)
      b_dt  <- data.table::as.data.table(b)
      ov <- cl_dt[b_dt, on = .(chr, start <= pos_end, end >= pos_ini), nomatch = 0L, allow.cartesian = TRUE]
      cl_sig <- unique(as.data.frame(ov[, .(cluster_id, chr, start, end)]))
    }
    validate(need(nrow(cl_sig) > 0, "No clusters matched Catalog hits (nothing to export)."))
    
    cluster_csv <- cl_sig[, c("cluster_id","chr","start","end")]
    names(cluster_csv) <- c("cluster_id","chr","cluster_start","cluster_end")
    
    # -----------------------------
    # Build candidate.csv (GWAS + Catalog)
    # -----------------------------
    h <- tryCatch(hits_df(), error = function(e) NULL)
    
    if (is.null(h) || !nrow(h)) {
      gwas_cand <- data.frame(chr=integer(), pos_ini=integer(), pos_end=integer(),
                              id_hit=character(), classe=character(), stringsAsFactors = FALSE)
    } else {
      validate(need(all(c("CHR","BP") %in% names(h)), "hits_df() must contain CHR and BP columns."))
      rsid <- if ("snp" %in% names(h)) as.character(h$snp) else NA_character_
      rsid[!nzchar(rsid)] <- paste0("chr", h$CHR, ":", h$BP)
      gwas_cand <- data.frame(
        chr     = norm_chr_int(h$CHR),
        pos_ini = as.integer(h$BP),
        pos_end = as.integer(h$BP),
        id_hit  = rsid,
        classe  = "GWAS",
        stringsAsFactors = FALSE
      )
    }
    
    candidate_csv <- rbind(
      gwas_cand,
      b[, c("chr","pos_ini","pos_end","id_hit","classe")]
    )
    
    list(cluster_csv = cluster_csv, candidate_csv = candidate_csv)
  }
  
  # escriu CSVs “virtuals” i deixa els DF a memòria per LD module
  refresh_catalog_ld_inputs <- function() {
    out <- build_catalog_ld_inputs()
    d   <- ld_inputs_dir()
    
    f1 <- file.path(d, "cluster_catalog.csv")
    f2 <- file.path(d, "candidate_catalog.csv")
    
    utils::write.csv(out$cluster_csv,  f1, row.names = FALSE, quote = FALSE)
    utils::write.csv(out$candidate_csv,f2, row.names = FALSE, quote = FALSE)
    
    rv_ld_inputs$cluster_path  <- f1
    rv_ld_inputs$candidate_path<- f2
    rv_ld_inputs$cluster_df    <- out$cluster_csv
    rv_ld_inputs$candidate_df  <- out$candidate_csv
  }
  
  observeEvent(list(clusters_val(), catalog_data(), hits_df()), {
    # només refresquem si hi ha dades mínimes
    ok1 <- is.data.frame(clusters_val()) && nrow(clusters_val()) > 0
    ok2 <- is.data.frame(catalog_data()) && nrow(catalog_data()) > 0
    if (ok1 && ok2) {
      tryCatch(refresh_catalog_ld_inputs(), error = function(e) NULL)
    }
  }, ignoreInit = TRUE)
  
  ##############################################################################
  ###########################Create cluster_candidate files#####################
  ##############################################################################
  build_ld_clusters_from_app <- function() {
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available. Build clusters first."))
    
    cl <- as.data.frame(cl)
    
    # Accepta variants de noms i les normalitza
    if (!"chr" %in% names(cl)   && "CHR" %in% names(cl))      cl$chr   <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end" %in% names(cl)   && "end_bp" %in% names(cl))   cl$end   <- cl$end_bp
    if (!"cluster_id" %in% names(cl) && "CLUSTER_ID" %in% names(cl)) cl$cluster_id <- cl$CLUSTER_ID
    
    validate(need(all(c("cluster_id","chr","start","end") %in% names(cl)),
                  paste0("clusters_val() must contain cluster_id, chr, start, end. Got: ",
                         paste(names(cl), collapse = ", "))))
    
    cl$chr   <- norm_chr_int(cl$chr)
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    validate(need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
    
    # IMPORTANT: els noms que vol el LD module
    out <- cl[, c("cluster_id","chr","start","end")]
    names(out) <- c("cluster_id","chr","cluster_start","cluster_end")
    out
  }
  
  build_ld_candidates_from_app <- function() {
    out <- data.frame(chr=integer(), pos_ini=integer(), pos_end=integer(),
                      id_hit=character(), classe=character(), stringsAsFactors = FALSE)
    
    # GWAS hits
    h <- tryCatch(hits_df(), error = function(e) NULL)
    if (is.data.frame(h) && nrow(h) > 0) {
      validate(need(all(c("CHR","BP") %in% names(h)),
                    paste0("hits_df() must contain CHR and BP. Got: ", paste(names(h), collapse = ", "))))
      
      rs <- if ("snp" %in% names(h)) as.character(h$snp)
      else if ("rsid" %in% names(h)) as.character(h$rsid)
      else if ("SNP" %in% names(h)) as.character(h$SNP)
      else NA_character_
      
      rs[is.na(rs) | !nzchar(rs)] <- paste0("chr", h$CHR, ":", h$BP)
      
      out <- rbind(out, data.frame(
        chr     = norm_chr_int(h$CHR),
        pos_ini = suppressWarnings(as.integer(h$BP)),
        pos_end = suppressWarnings(as.integer(h$BP)),
        id_hit  = rs,
        classe  = "GWAS",
        stringsAsFactors = FALSE
      ))
    }
    
    # Catalog hits
    cd <- tryCatch(catalog_data(), error = function(e) NULL)
    if (is.data.frame(cd) && nrow(cd) > 0) {
      
      chr_col   <- pick_first_col(cd, c("chr","CHR","chrom","chromosome","CHR_ID"))
      start_col <- pick_first_col(cd, c("pos_ini","start","START","bp_start","POS","BP","pos"))
      end_col   <- pick_first_col(cd, c("pos_end","end","END","bp_end"))
      if (is.null(end_col)) end_col <- start_col
      
      validate(need(!is.null(chr_col) && !is.null(start_col),
                    paste0("catalog_data() missing chr/pos. Got: ", paste(names(cd), collapse=", "))))
      
      id_col <- pick_first_col(cd, c("id_hit","snp","SNPS","rsid","RSID","variant","ID","code","CODE"))
      if (is.null(id_col)) {
        id_col <- ".id_tmp"
        cd[[id_col]] <- paste0("Catalog_chr", cd[[chr_col]], ":", cd[[start_col]], "-", cd[[end_col]])
      }
      
      out <- rbind(out, data.frame(
        chr     = norm_chr_int(cd[[chr_col]]),
        pos_ini = suppressWarnings(as.integer(cd[[start_col]])),
        pos_end = suppressWarnings(as.integer(cd[[end_col]])),
        id_hit  = as.character(cd[[id_col]]),
        classe  = "Catalog_hit",
        stringsAsFactors = FALSE
      ))
    }
    
    out <- out[is.finite(out$chr) & is.finite(out$pos_ini) & is.finite(out$pos_end) &
                 out$pos_end >= out$pos_ini & nzchar(out$id_hit), , drop = FALSE]
    
    validate(need(nrow(out) > 0, "No candidates available for LD."))
    out
  }
  
  output$dl_candidates_zip <- downloadHandler(
    filename = function() {
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      tag_txt <- tryCatch({
        tg <- make_mode_thr_tag(input$cluster_method %||% "window", input$pthr, input$min_logp)
        paste0(tg$mode_tag, "_thr", tg$thr_txt)
      }, error = function(e) "run")
      
      paste0("catalog_candidates_", tag_txt, "_", stamp, ".zip")
    },
    
    content = function(file) {
      
      # -----------------------------
      # Build CANONICAL clusters/candidates (same as LD module)
      # -----------------------------
      cluster_csv   <- build_ld_clusters_from_app()
      candidate_csv <- build_ld_candidates_from_app()
      
      validate(need(is.data.frame(cluster_csv)   && nrow(cluster_csv)   > 0, "No clusters available to export."))
      validate(need(is.data.frame(candidate_csv) && nrow(candidate_csv) > 0, "No candidates available to export."))
      
      # -----------------------------
      # OPTIONAL: keep only clusters that have candidate overlap
      # (matches your previous behavior: export only "clusters matched")
      # -----------------------------
      suppressWarnings({
        cl <- as.data.frame(cluster_csv)
        cd <- as.data.frame(candidate_csv)
      })
      
      # Overlap: cluster [start,end] intersects candidate [pos_ini,pos_end]
      # Candidate pos_end may be == pos_ini for GWAS
      if (nrow(cl) > 0 && nrow(cd) > 0) {
        # fast overlap with data.table if available; fallback to base otherwise
        if (requireNamespace("data.table", quietly = TRUE)) {
          cl_dt <- data.table::as.data.table(cl)
          cd_dt <- data.table::as.data.table(cd)
          
          ov <- cl_dt[cd_dt,
                      on = .(chr, cluster_start <= pos_end, cluster_end >= pos_ini),
                      nomatch = 0L, allow.cartesian = TRUE]
          
          cl_sig <- unique(as.data.frame(ov[, .(cluster_id, chr, cluster_start, cluster_end)]))
        } else {
          # base fallback
          keep <- logical(nrow(cl))
          for (i in seq_len(nrow(cl))) {
            keep[i] <- any(cd$chr == cl$chr[i] &
                             cd$pos_end >= cl$cluster_start[i] &
                             cd$pos_ini <= cl$cluster_end[i])
          }
          cl_sig <- cl[keep, , drop = FALSE]
        }
        
        validate(need(nrow(cl_sig) > 0, "No clusters matched candidates (nothing to export)."))
        cluster_csv <- cl_sig
      }
      
      # -----------------------------
      # Write to temp folder and zip
      # -----------------------------
      tmpdir <- tempfile("catalog_export_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      f1 <- file.path(tmpdir, "cluster_catalog.csv")
      f2 <- file.path(tmpdir, "candidate_catalog.csv")
      
      utils::write.csv(cluster_csv,   f1, row.names = FALSE, quote = FALSE)
      utils::write.csv(candidate_csv, f2, row.names = FALSE, quote = FALSE)
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      setwd(tmpdir)
      
      utils::zip(zipfile = file, files = c("cluster_catalog.csv", "candidate_catalog.csv"))
    }
  )
  
  
  #  ##############################################################################
  #  ############################ LD Module   #####################################
  #  ################# GWAS -> LD module canonical cluster input###################
  
  ld_mod <- ld_module_server(
    id = "ld",
    app_tag = "catalog",
    activate_r = reactive(identical(input$topnav, "ld_tab"))
  )
  
  ##############################################################################
  ##############################    RESET    ###################################
  ##############################################################################
  
  # dins el server del modul corresponent
  resetting <- reactiveVal(FALSE)
  
  observeEvent(input$reset_case, {
    req(!isTRUE(resetting()))
    resetting(TRUE)
    on.exit(resetting(FALSE), add = TRUE)
    
    safe <- function(expr) try(expr, silent = TRUE)
    
    # (A) Tanca modals
    safe(shiny::removeModal())
    
    # (B) Neteja seleccions DT (IMPORTANT en modul: usar session$ns)
    dt_unselect <- function(id_no_ns) {
      id <- session$ns(id_no_ns)
      proxy <- safe(DT::dataTableProxy(id, session = session))
      if (!inherits(proxy, "try-error")) {
        safe(DT::selectRows(proxy, NULL))
        safe(DT::selectColumns(proxy, NULL))
        safe(DT::selectCells(proxy, NULL))
      }
    }
    
    for (id in c(
      "hits_tbl", "cluster_dt", "catalog_table",
      "enrich_disease_table", "enrich_trait_table", "go_table", "kegg_table",
      "ld_pairs_dt"
    )) dt_unselect(id)
    
    # (C) Borra fitxers temporals (si tens workdir al scope del modul)
    wd <- safe({
      if (exists("workdir", inherits = TRUE)) {
        if (is.function(workdir)) workdir() else workdir
      } else NULL
    })
    
    if (!inherits(wd, "try-error") && is.character(wd) && nzchar(wd) && dir.exists(wd)) {
      for (f in c("selected_intervals.range", "selected_intervals.vcf", "selected_intervals.bed")) {
        ff <- file.path(wd, f)
        if (file.exists(ff)) safe(unlink(ff, force = TRUE))
      }
    }
    
    isolate({
      # (D) Reset del "core state" derivat (NO toquis catalog_all!)
      safe(clusters_val(NULL))
      safe(combo_region(NULL))
      safe(ucsc_region(NULL))
      
      # ✅ deixa catalog_all() intacte
      safe(catalog_final_path_csv(NULL))
      safe(catalog_final_path_rds(NULL))
      
      safe(gene_map_cache(NULL))
      safe(gene_map_cache_syms(NULL))
      
      safe(go_trigger(0L))
      safe(kegg_trigger(0L))
      
      safe(catalog_log_text(""))
      
      # (E) Reset LD state (si existeix)
      if (exists("ld_state", inherits = TRUE)) {
        safe({ ld_state$running      <- FALSE })
        safe({ ld_state$fl           <- NULL })
        safe({ ld_state$ld_pairs     <- NULL })
        safe({ ld_state$blocks_ij    <- NULL })
        safe({ ld_state$subset_prefix <- NULL })
        safe({ ld_state$tag          <- NULL })
        safe({ ld_state$log          <- "" })
      }
      
      # (F) Re-assegura que botons no queden disabled
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        safe(shinyjs::enable(session$ns("run_catalog")))
        safe(shinyjs::enable(session$ns("run_ld_cluster")))
        safe(shinyjs::enable(session$ns("run_ld_blocks")))
      }
    })
    
    # (G) Reset inputs (només si existeixen)
    upd <- function(id, fun) if (id %in% names(input)) safe(fun())
    
    
    upd("cluster_method",        function() shiny::updateRadioButtons(session, "cluster_method", selected = "window"))
    upd("pthr",                  function() shiny::updateSliderInput(session, "pthr", value = 5))
    upd("min_logp",              function() shiny::updateSliderInput(session, "min_logp", value = 5))
    upd("flank",                 function() shiny::updateNumericInput(session, "flank", value = 10000))
    upd("min_hits",              function() shiny::updateNumericInput(session, "min_hits", value = 3))
    upd("win_bp",                function() shiny::updateNumericInput(session, "win_bp", value = 1000000))
    upd("step_bp",               function() shiny::updateNumericInput(session, "step_bp", value = 250000))
    
    shiny::showNotification("Reset done. Build clusters again, then run 'Extract catalog hits'.",
                            type = "message", duration = 3)
  }, ignoreInit = TRUE)
  
  ##############################
  
}

shinyApp(ui, server)
