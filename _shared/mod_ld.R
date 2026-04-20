# R/mod_ld.R
# ------------------------------------------------------------------
# LD module
# - Reads shared RDS from gi_shared_root/integrator_exports
# - Uses:
#     * <app>_clusters_master.rds
#     * <app>_candidates.rds
# - Two workflows:
#     1) Plot / full LD matrix for selected cluster
#     2) Prioritization / seed + proxies for all GWAS hits in cluster
# ------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

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

chr_label_plink <- function(chr_num) {
  out <- as.character(chr_num)
  out[out == "23"] <- "X"
  out[out == "24"] <- "Y"
  out[out == "26"] <- "MT"
  out
}

# ============================================================
# Harmonització mínima i única per integrator exports
# ============================================================

harmonize_integrator_clusters <- function(df) {
  stopifnot(is.data.frame(df))
  
  need <- c("cluster_id", "chr", "start", "end", "cluster_key")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("clusters_master missing columns: ", paste(miss, collapse = ", "))
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  df$cluster_id  <- trimws(as.character(df$cluster_id))
  df$chr         <- chr_map_plink19(df$chr)
  df$start       <- suppressWarnings(as.integer(readr::parse_number(as.character(df$start))))
  df$end         <- suppressWarnings(as.integer(readr::parse_number(as.character(df$end))))
  df$cluster_key <- trimws(as.character(df$cluster_key))
  
  s0 <- pmin(df$start, df$end)
  e0 <- pmax(df$start, df$end)
  df$start <- s0
  df$end   <- e0
  
  idx_key <- is.na(df$cluster_key) | !nzchar(df$cluster_key)
  df$cluster_key[idx_key] <- paste0(
    df$cluster_id[idx_key], "|chr", chr_label_plink(df$chr[idx_key]), ":",
    df$start[idx_key], "-", df$end[idx_key]
  )
  
  df <- df[
    !is.na(df$cluster_id) & nzchar(df$cluster_id) &
      !is.na(df$chr) &
      !is.na(df$start) &
      !is.na(df$end),
    need,
    drop = FALSE
  ]
  
  df <- unique(df)
  rownames(df) <- NULL
  df
}

harmonize_integrator_candidates <- function(df) {
  stopifnot(is.data.frame(df))
  
  need <- c("cluster_id", "chr", "pos_ini", "pos_end", "id_hit", "rsid", "position", "classe")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("candidates missing columns: ", paste(miss, collapse = ", "))
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  df$cluster_id <- trimws(as.character(df$cluster_id))
  df$chr        <- chr_map_plink19(df$chr)
  df$pos_ini    <- suppressWarnings(as.integer(readr::parse_number(as.character(df$pos_ini))))
  df$pos_end    <- suppressWarnings(as.integer(readr::parse_number(as.character(df$pos_end))))
  df$id_hit     <- trimws(as.character(df$id_hit))
  df$rsid       <- trimws(as.character(df$rsid))
  df$position   <- suppressWarnings(as.integer(readr::parse_number(as.character(df$position))))
  df$classe     <- trimws(as.character(df$classe))
  
  # backfill de seguretat
  idx_rsid <- is.na(df$rsid) | !nzchar(df$rsid)
  df$rsid[idx_rsid] <- df$id_hit[idx_rsid]
  
  idx_pos <- is.na(df$position)
  df$position[idx_pos] <- df$pos_ini[idx_pos]
  
  idx_end <- is.na(df$pos_end)
  df$pos_end[idx_end] <- df$pos_ini[idx_end]
  
  s0 <- pmin(df$pos_ini, df$pos_end)
  e0 <- pmax(df$pos_ini, df$pos_end)
  df$pos_ini <- s0
  df$pos_end <- e0
  
  df <- df[
    !is.na(df$cluster_id) & nzchar(df$cluster_id) &
      !is.na(df$chr) &
      !is.na(df$pos_ini) &
      !is.na(df$pos_end) &
      !is.na(df$id_hit) & nzchar(df$id_hit),
    need,
    drop = FALSE
  ]
  
  df <- unique(df)
  rownames(df) <- NULL
  df
}

prepare_clusters <- function(df) {
  chr_col <- pick_col(df, c("chr", "CHR", "chrom", "CHROM", "chromosome"))
  st_col  <- pick_col(df, c("start", "START", "start_bp", "cluster_start", "FROM", "from", "bp1"))
  en_col  <- pick_col(df, c("end", "END", "end_bp", "cluster_end", "TO", "to", "bp2"))
  id_col  <- pick_col(df, c("cluster_id", "CLUSTER_ID", "cluster", "id"))
  key_col <- pick_col(df, c("cluster_key", "CLUSTER_KEY"))
  
  if (is.null(chr_col) || is.null(st_col) || is.null(en_col)) {
    stop("prepare_clusters: missing chr/start/end columns. Found: ", paste(names(df), collapse = ", "))
  }
  
  out <- df %>%
    dplyr::transmute(
      chr = chr_map_plink19(.data[[chr_col]]),
      start = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
      end   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]])))),
      cluster_id = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_,
      cluster_key = if (!is.null(key_col)) as.character(.data[[key_col]]) else NA_character_
    ) %>%
    dplyr::filter(is.finite(chr), is.finite(start), is.finite(end)) %>%
    dplyr::mutate(
      start0 = start,
      end0   = end,
      start  = pmin(start0, end0),
      end    = pmax(start0, end0),
      cluster_id = trimws(cluster_id),
      cluster_key = trimws(cluster_key)
    ) %>%
    dplyr::select(-start0, -end0) %>%
    dplyr::arrange(chr, start, end) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(cluster_n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cluster_id = ifelse(
        is.na(cluster_id) | !nzchar(cluster_id),
        paste0("cluster_chr", chr_label_plink(chr), "_", cluster_n),
        cluster_id
      ),
      cluster_key = dplyr::if_else(
        !is.na(cluster_key) & nzchar(cluster_key),
        cluster_key,
        paste0(cluster_id, "|chr", chr_label_plink(chr), ":", start, "-", end)
      ),
      size_kb = round((end - start) / 1000, 2)
    )
  
  out
}

prepare_candidates <- function(df) {
  rs_col  <- pick_col(df, c("rsid", "RSID", "snp", "SNP", "id", "ID", "marker", "MARKER", "id_hit", "ID_HIT"))
  chr_col <- pick_col(df, c("chr", "CHR", "chrom", "CHROM", "chromosome"))
  pos_col <- pick_col(df, c("position", "POS", "pos", "BP", "bp", "pos_ini", "pos_start", "start", "POS_INI"))
  end_col <- pick_col(df, c("pos_end", "POS_END", "end", "END"))
  cl_col  <- pick_col(df, c("classe", "class", "CLASS", "Classe", "type", "TYPE"))
  id_col  <- pick_col(df, c("cluster_id", "CLUSTER_ID", "cluster", "id"))
  
  if (is.null(rs_col) || is.null(chr_col) || is.null(pos_col)) {
    stop(
      "Candidates table must have id_hit/rsid + chr + position/pos_ini columns. Found: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  out <- df %>%
    dplyr::transmute(
      cluster_id = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_,
      rsid       = as.character(.data[[rs_col]]),
      chr        = chr_map_plink19(.data[[chr_col]]),
      position   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
      pos_end    = if (!is.null(end_col)) suppressWarnings(as.integer(readr::parse_number(as.character(.data[[end_col]])))) else NA_integer_,
      classe     = if (!is.null(cl_col)) as.character(.data[[cl_col]]) else NA_character_
    ) %>%
    dplyr::mutate(
      cluster_id = trimws(cluster_id),
      rsid       = trimws(rsid),
      classe     = trimws(classe),
      pos_end    = dplyr::coalesce(pos_end, position)
    ) %>%
    dplyr::filter(nzchar(rsid), is.finite(chr), is.finite(position)) %>%
    dplyr::distinct(cluster_id, rsid, chr, position, classe, .keep_all = TRUE)
  
  out
}

available_pops_dir <- function(ld_pops_dir) {
  d <- tryCatch(normalizePath(ld_pops_dir, winslash = "/", mustWork = FALSE),
                error = function(e) ld_pops_dir)
  if (is.null(d) || !nzchar(d) || !dir.exists(d)) return(character(0))
  
  ff <- list.files(d, pattern = "\\.[Tt][Xx][Tt]$", full.names = FALSE)
  pops <- sub("\\.[Tt][Xx][Tt]$", "", ff)
  pops <- pops[nzchar(pops)]
  sort(unique(pops))
}

read_keep_file_dir <- function(pop, ld_pops_dir) {
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
    file.path(d, paste0(pop, ".TXT")),
    file.path(d, paste0(pop, ".TxT"))
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f) || !nzchar(f)) {
    av <- available_pops_dir(d)
    stop("Missing keep file for pop='", pop, "' in ", d,
         ". Available: ", paste(av, collapse = ", "))
  }
  
  x <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 1) stop("Keep file must have at least 1 column (IID)")
  
  if (ncol(x) == 1) x <- cbind(x[, 1], x[, 1]) else x <- x[, 1:2]
  colnames(x) <- c("FID", "IID")
  
  out <- file.path(tempdir(), paste0("keep_", pop, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  write.table(x, out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  out
}

gi_ld_defaults <- function() {
  if (exists("gi_cfg", mode = "function")) {
    cfg <- gi_cfg()
    res <- cfg$resources %||% ""
    popdir <- cfg$pop_dir %||% ""
    shared <- cfg$shared %||% ""
  } else {
    res    <- Sys.getenv("GITOOLS_RESOURCES", unset = "")
    popdir <- Sys.getenv("GITOOLS_POP_DIR", unset = "")
    shared <- Sys.getenv("GITOOLS_SHARED", unset = "")
  }
  
  res    <- normalizePath(res, winslash = "/", mustWork = FALSE)
  popdir <- normalizePath(popdir, winslash = "/", mustWork = FALSE)
  shared <- normalizePath(shared, winslash = "/", mustWork = FALSE)
  
  plink <- Sys.getenv("GITOOLS_PLINK19", unset = "")
  if (!nzchar(plink) && nzchar(res)) {
    cand1 <- file.path(res, "software", "plink19", "plink")
    cand2 <- file.path(res, "software", "plink19", "plink19")
    cand3 <- file.path(res, "software", "plink19")
    plink <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else cand3
  }
  
  bfile <- Sys.getenv("GITOOLS_LD_BFILE", unset = "")
  if (!nzchar(bfile) && nzchar(res)) {
    bfile <- file.path(
      res, "LD_resources",
      "Merged_FULL_SET_hg38_hgdp.wgs_10000G3_ko07_MAF0.05"
    )
  }
  
  integrator_dir <- if (nzchar(shared)) file.path(shared, "integrator_exports") else ""
  
  list(
    plink = plink,
    bfile = bfile,
    popdir = popdir,
    resources = res,
    shared = shared,
    integrator_dir = integrator_dir
  )
}

read_shared_rds_safe <- function(path, label = "RDS", logger = NULL) {
  logf <- function(...) {
    if (is.function(logger)) logger(...)
  }
  
  if (is.null(path) || !nzchar(path)) {
    logf("[RDS][", label, "] empty path")
    return(NULL)
  }
  
  logf("[RDS][", label, "] path=", path)
  logf("[RDS][", label, "] exists=", file.exists(path))
  
  if (!file.exists(path)) return(NULL)
  
  x <- tryCatch(
    readRDS(path),
    error = function(e) {
      logf("[RDS][", label, "][ERROR] ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.data.frame(x)) {
    logf("[RDS][", label, "] nrow=", nrow(x), " ncol=", ncol(x))
    logf("[RDS][", label, "] cols=", paste(names(x), collapse = ", "))
    logf(
      "[RDS][", label, "] head:\n",
      paste(capture.output(print(utils::head(x, 5))), collapse = "\n")
    )
  } else {
    logf("[RDS][", label, "] object is not a data.frame")
  }
  
  x
}

run_plink <- function(args, out_prefix, plink_bin) {
  res <- tryCatch(
    system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) paste("system2 error:", conditionMessage(e))
  )
  st <- attr(res, "status")
  if (is.null(st)) st <- 0
  
  log_file <- paste0(out_prefix, ".log")
  log_tail <- character(0)
  if (file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE)
    log_tail <- utils::tail(log_lines, 40)
  }
  list(stdout = res, status = st, log_file = log_file, log_tail = log_tail)
}

read_plink_ld <- function(prefix) {
  f_gz <- paste0(prefix, ".ld.gz")
  f_tx <- paste0(prefix, ".ld")
  f <- if (file.exists(f_gz)) f_gz else if (file.exists(f_tx)) f_tx else NA_character_
  if (is.na(f)) stop("Missing LD file: expected ", f_gz, " or ", f_tx)
  
  con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con), add = TRUE)
  
  x <- utils::read.table(
    con, header = TRUE, sep = "",
    stringsAsFactors = FALSE, fill = TRUE, comment.char = ""
  )
  names(x) <- trimws(names(x))
  x
}

read_plink_blocks <- function(prefix) {
  f <- paste0(prefix, ".blocks.det")
  if (!file.exists(f)) stop("Missing blocks file: ", f)
  
  x <- tryCatch(
    read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = ""),
    error = function(e) NULL
  )
  
  if (is.null(x) || nrow(x) == 0 || ncol(x) <= 1) {
    x <- tryCatch(
      read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "", quote = "", comment.char = ""),
      error = function(e) NULL
    )
  }
  
  if (is.null(x) || nrow(x) == 0) {
    stop("Could not parse .blocks.det (bad delimiter or empty). File=", f)
  }
  x
}

blocks_det_to_ij <- function(blocks_raw, fl) {
  stopifnot(is.data.frame(blocks_raw), is.data.frame(fl))
  stopifnot(all(c("SNP", "BP") %in% names(fl)))
  
  fl <- fl %>%
    dplyr::arrange(BP, SNP) %>%
    dplyr::mutate(ix = dplyr::row_number())
  
  snp_to_ix <- stats::setNames(fl$ix, fl$SNP)
  
  nm <- names(blocks_raw)
  snps_col <- nm[match(TRUE, nm %in% c("SNPS", "Snps", "snps", "SNPs"))]
  bp1_col  <- nm[match(TRUE, nm %in% c("BP1", "bp1", "FROM_BP", "from_bp", "START_BP", "start_bp"))]
  bp2_col  <- nm[match(TRUE, nm %in% c("BP2", "bp2", "TO_BP", "to_bp", "END_BP", "end_bp"))]
  
  out <- NULL
  
  if (!is.na(snps_col)) {
    snps_str <- as.character(blocks_raw[[snps_col]])
    snps_str[is.na(snps_str)] <- ""
    
    split_snps <- function(s) {
      s <- trimws(s)
      if (!nzchar(s)) return(character(0))
      if (grepl("\\|", s, fixed = FALSE)) return(strsplit(s, "\\|", fixed = FALSE)[[1]])
      if (grepl(",", s, fixed = TRUE))   return(strsplit(s, ",", fixed = TRUE)[[1]])
      strsplit(s, "\\s+")[[1]]
    }
    
    lst <- lapply(snps_str, split_snps)
    first_snp <- vapply(lst, function(v) if (length(v)) v[1] else NA_character_, character(1))
    last_snp  <- vapply(lst, function(v) if (length(v)) v[length(v)] else NA_character_, character(1))
    
    i <- suppressWarnings(as.integer(snp_to_ix[first_snp]))
    j <- suppressWarnings(as.integer(snp_to_ix[last_snp]))
    
    out <- data.frame(i = i, j = j, stringsAsFactors = FALSE)
  }
  
  if (is.null(out) || all(is.na(out$i)) || all(is.na(out$j))) {
    if (is.na(bp1_col) || is.na(bp2_col)) {
      return(tibble::tibble(i = integer(0), j = integer(0)))
    }
    
    bp1 <- suppressWarnings(as.numeric(blocks_raw[[bp1_col]]))
    bp2 <- suppressWarnings(as.numeric(blocks_raw[[bp2_col]]))
    
    bp_to_ix <- function(bp) {
      if (!is.finite(bp)) return(NA_integer_)
      k <- which.min(abs(fl$BP - bp))
      fl$ix[k]
    }
    i <- vapply(bp1, bp_to_ix, integer(1))
    j <- vapply(bp2, bp_to_ix, integer(1))
    
    out <- data.frame(i = i, j = j, stringsAsFactors = FALSE)
  }
  
  out <- out %>%
    dplyr::mutate(
      i = suppressWarnings(as.integer(i)),
      j = suppressWarnings(as.integer(j))
    ) %>%
    dplyr::filter(is.finite(i), is.finite(j), i >= 1, j >= 1, j > i, j <= nrow(fl)) %>%
    dplyr::distinct(i, j, .keep_all = TRUE) %>%
    dplyr::arrange(i, j)
  
  tibble::as_tibble(out)
}

find_seed_snp <- function(hit_rsid, hit_pos, bim_df) {
  stopifnot(is.data.frame(bim_df), all(c("SNP", "BP") %in% names(bim_df)))
  
  hit_rsid <- trimws(as.character(hit_rsid %||% ""))
  hit_pos  <- suppressWarnings(as.integer(hit_pos))
  
  if (nzchar(hit_rsid) && hit_rsid %in% bim_df$SNP) {
    row <- bim_df[bim_df$SNP == hit_rsid, , drop = FALSE][1, , drop = FALSE]
    return(data.frame(
      seed_snp = as.character(row$SNP),
      seed_pos = as.integer(row$BP),
      seed_type = "exact",
      seed_dist_bp = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  if (!is.finite(hit_pos) || nrow(bim_df) == 0) {
    return(data.frame(
      seed_snp = NA_character_,
      seed_pos = NA_integer_,
      seed_type = "missing",
      seed_dist_bp = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  
  k <- which.min(abs(bim_df$BP - hit_pos))
  data.frame(
    seed_snp = as.character(bim_df$SNP[k]),
    seed_pos = as.integer(bim_df$BP[k]),
    seed_type = "nearest",
    seed_dist_bp = as.integer(abs(bim_df$BP[k] - hit_pos)),
    stringsAsFactors = FALSE
  )
}

compute_seed_ld <- function(subset_prefix, seed_snp, plink_bin, out_prefix, metric = "R2") {
  args <- c(
    "--bfile", subset_prefix,
    "--ld-snp", seed_snp,
    "--ld-window", "999999",
    "--ld-window-kb", "999999",
    "--ld-window-r2", "0"
  )
  
  if (identical(metric, "Dprime")) {
    args <- c(args, "--r2", "dprime", "gz")
  } else {
    args <- c(args, "--r2", "gz")
  }
  
  args <- c(args, "--out", out_prefix)
  run_plink(args, out_prefix, plink_bin)
}

build_hit_proxy_table <- function(hit_row, bim_df, subset_prefix, plink_bin, workdir, metric = "R2") {
  hit_rsid <- as.character(hit_row$rsid[1])
  hit_pos  <- as.integer(hit_row$position[1])
  cl_id    <- as.character(hit_row$cluster_id[1] %||% NA_character_)
  
  seed <- find_seed_snp(hit_rsid = hit_rsid, hit_pos = hit_pos, bim_df = bim_df)
  
  if (!is.finite(seed$seed_pos[1]) || is.na(seed$seed_snp[1]) || !nzchar(seed$seed_snp[1])) {
    return(NULL)
  }
  
  out_prefix <- file.path(
    workdir,
    paste0("seed_", gsub("[^A-Za-z0-9_\\-]", "_", seed$seed_snp[1]), "_", format(Sys.time(), "%H%M%S%OS3"))
  )
  
  rr <- compute_seed_ld(
    subset_prefix = subset_prefix,
    seed_snp = seed$seed_snp[1],
    plink_bin = plink_bin,
    out_prefix = out_prefix,
    metric = metric
  )
  if (is.null(rr$status) || rr$status != 0) return(NULL)
  
  ld <- read_plink_ld(out_prefix)
  val_col <- if (identical(metric, "Dprime")) pick_col(ld, c("Dprime", "D'", "DP", "DPRIME")) else pick_col(ld, c("R2", "r2"))
  snpA <- pick_col(ld, c("SNP_A", "SNP_A1", "SNP1"))
  snpB <- pick_col(ld, c("SNP_B", "SNP_B1", "SNP2"))
  
  if (is.null(val_col) || is.null(snpA) || is.null(snpB)) return(NULL)
  
  out <- ld %>%
    dplyr::transmute(
      SNP_A = as.character(.data[[snpA]]),
      SNP_B = as.character(.data[[snpB]]),
      ld_value = suppressWarnings(as.numeric(.data[[val_col]]))
    ) %>%
    dplyr::filter(is.finite(ld_value)) %>%
    dplyr::mutate(
      proxy_snp = dplyr::if_else(SNP_A == seed$seed_snp[1], SNP_B, SNP_A)
    ) %>%
    dplyr::left_join(
      bim_df %>% dplyr::transmute(proxy_snp = SNP, proxy_pos = BP),
      by = "proxy_snp"
    ) %>%
    dplyr::mutate(
      cluster_id    = cl_id,
      query_hit     = hit_rsid,
      query_pos     = hit_pos,
      seed_snp      = seed$seed_snp[1],
      seed_pos      = seed$seed_pos[1],
      seed_type     = seed$seed_type[1],
      seed_dist_bp  = seed$seed_dist_bp[1]
    ) %>%
    dplyr::select(
      cluster_id, query_hit, query_pos,
      seed_snp, seed_pos, seed_type, seed_dist_bp,
      proxy_snp, proxy_pos, ld_value
    ) %>%
    dplyr::distinct()
  
  out
}


# ------------------------------------------------------------------
# UI
# ------------------------------------------------------------------
ld_module_ui <- function(id, app_tag) {
  app_tag <- match.arg(app_tag, c("catalog", "gtex", "nonsyn", "ewastum", "ewasdis"))
  
  ns <- shiny::NS(id)
  ld_def <- gi_ld_defaults()
  
  shiny::tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        
        shiny::tags$details(
          shiny::tags$summary("Advanced LD paths"),
          shiny::wellPanel(
            shiny::textInput(ns("plink_bin"), "PLINK binary path", value = ld_def$plink %||% ""),
            shiny::textInput(ns("bfile_ref_prefix"), "Reference bfile prefix", value = ld_def$bfile %||% ""),
            shiny::textInput(ns("ld_pops_dir"), "POP keep-files directory (.txt)", value = ld_def$popdir %||% ""),
            shiny::textInput(ns("integrator_dir"), "Shared integrator_exports dir", value = ld_def$integrator_dir %||% ""),
            shiny::textInput(ns("clusters_rds_name"), "Clusters master RDS", value = paste0(app_tag, "_clusters_master.rds")),
            shiny::textInput(ns("candidates_rds_name"), "Candidates RDS", value = paste0(app_tag, "_candidates.rds"))
          )
        ),
        
        shiny::uiOutput(ns("ld_pop_ui")),
        shiny::tags$hr(),
        
        shiny::selectInput(ns("ld_chr"), "Chromosome", choices = character(0)),
        shiny::selectInput(ns("ld_cluster_id"), "Cluster", choices = character(0)),
        shiny::tags$hr(),
        
        shiny::radioButtons(
          ns("ld_metric"), "LD metric",
          choices = c("r²" = "R2", "D′" = "Dprime"),
          selected = "R2", inline = TRUE
        ),
        
        shiny::radioButtons(
          ns("ld_x_mode"), "X spacing",
          choices = c("Genomic (BP)" = "bp", "Equal spacing" = "equal"),
          selected = "bp"
        ),
        
        shiny::numericInput(ns("ld_max_snps_interval"), "Max SNPs in interval (plot safety)", value = 400, min = 50, step = 50),
        shiny::numericInput(ns("ld_proxy_r2_min"), "Min r² for priority proxies", value = 0.6, min = 0, max = 1, step = 0.05),
        
        shiny::tags$hr(),
        shiny::actionButton(ns("run_ld_cluster"), "Run LD plot", icon = shiny::icon("chart-area")),
        shiny::actionButton(ns("run_ld_prior"), "Run LD prioritization", icon = shiny::icon("list-check")),
        shiny::tags$hr(),
        
        shiny::checkboxInput(ns("ld_blocks_enable"), "Haploview blocks (plot mode)", value = FALSE),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("ld_blocks_enable")),
          shiny::actionButton(ns("run_ld_blocks"), "Run blocks (PLINK)", icon = shiny::icon("th")),
          shiny::tags$hr(),
          shiny::tags$div(
            style = "margin-left:10px; padding:8px; border-left:3px solid #ddd;",
            shiny::checkboxInput(ns("ld_blocks_no_small_max_span"), "no-small-max-span", value = TRUE),
            shiny::sliderInput(ns("ld_blocks_inform_frac"), "blocks-inform-frac", min = 0.10, max = 1.00, value = 0.60, step = 0.05),
            shiny::numericInput(ns("ld_blocks_max_kb"), "blocks-max-kb", value = 1000, min = 1, step = 50),
            shiny::numericInput(ns("ld_blocks_min_maf"), "blocks-min-maf", value = 0.05, min = 0, max = 0.50, step = 0.01),
            shiny::sliderInput(ns("ld_blocks_strong_highci"), "blocks-strong-highci", min = 0.84, max = 0.99, value = 0.90, step = 0.01),
            shiny::sliderInput(ns("ld_blocks_strong_lowci"), "blocks-strong-lowci", min = 0.10, max = 0.95, value = 0.55, step = 0.01)
          )
        ),
        
        shiny::tags$hr(),
        shiny::verbatimTextOutput(ns("ld_log"))
      ),
      
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel(
            "LD plot",
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("ld_plot"), height = "800px")),
            shiny::tags$hr(),
            shiny::h4("SNPs in interval"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_fl_dt"))),
            shiny::tags$hr(),
            shiny::h4("Top LD pairs"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_pairs_dt")))
          ),
          shiny::tabPanel(
            "Prioritization",
            shiny::h4("Seed summary"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_seed_dt"))),
            shiny::tags$hr(),
            shiny::h4("Proxy table"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_priority_dt")))
          ),
          shiny::tabPanel(
            "Debug shared RDS",
            shiny::h4("Clusters master"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_debug_clusters_dt"))),
            shiny::tags$hr(),
            shiny::h4("Candidates"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ld_debug_candidates_dt")))
          )
        )
      )
    )
  )
}

# ------------------------------------------------------------------
# SERVER
# ------------------------------------------------------------------
ld_module_server <- function(
    id,
    app_tag,
    activate_r = NULL,
    ...
) {
  app_tag <- match.arg(app_tag, c("catalog", "gtex", "nonsyn", "ewastum", "ewasdis"))
  
  moduleServer(id, function(input, output, session) {
    
    txdb <- NULL
    
    `%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b
    
    app_hit_class <- paste0(app_tag, "_hit")
    default_clusters_rds_name   <- paste0(app_tag, "_clusters_master.rds")
    default_candidates_rds_name <- paste0(app_tag, "_candidates.rds")
    
    normalize_classe <- function(x) {
      x <- trimws(as.character(x))
      dplyr::case_when(
        x %in% c("Catalog_hit", "catalog_hit") ~ "catalog_hit",
        x %in% c("GTEx_hit", "gtex_hit") ~ "gtex_hit",
        x %in% c("NonSyn_hit", "nonsyn_hit") ~ "nonsyn_hit",
        x %in% c("EWASTum_hit", "ewastum_hit") ~ "ewastum_hit",
        x %in% c("EWASDis_hit", "ewasdis_hit") ~ "ewasdis_hit",
        x %in% c("LD_hit", "ld_hit") ~ "ld_hit",
        x == "GWAS" ~ "GWAS",
        TRUE ~ x
      )
    }
    
    reassign_candidates_to_master_clusters <- function(ca, cl, logger = NULL) {
      logf <- function(...) {
        if (is.function(logger)) logger(...)
      }
      
      if (!is.data.frame(ca) || !nrow(ca)) return(ca)
      if (!is.data.frame(cl) || !nrow(cl)) return(ca)
      
      need_ca <- c("cluster_id", "chr", "position")
      need_cl <- c("cluster_id", "chr", "start", "end")
      
      if (!all(need_ca %in% names(ca))) {
        logf("[REMAP] candidates missing columns: ", paste(setdiff(need_ca, names(ca)), collapse = ", "))
        return(ca)
      }
      if (!all(need_cl %in% names(cl))) {
        logf("[REMAP] clusters missing columns: ", paste(setdiff(need_cl, names(cl)), collapse = ", "))
        return(ca)
      }
      
      ca2 <- ca %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          chr = suppressWarnings(as.integer(chr)),
          position = suppressWarnings(as.integer(position))
        )
      
      cl2 <- cl %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          chr = suppressWarnings(as.integer(chr)),
          start = suppressWarnings(as.integer(start)),
          end = suppressWarnings(as.integer(end))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          is.finite(chr), is.finite(start), is.finite(end)
        ) %>%
        dplyr::distinct(cluster_id, chr, start, end, .keep_all = TRUE)
      
      valid_ids <- unique(cl2$cluster_id)
      
      ca_exact <- ca2 %>%
        dplyr::mutate(cluster_id_valid = !is.na(cluster_id) & nzchar(cluster_id) & cluster_id %in% valid_ids)
      
      ca_keep <- ca_exact %>%
        dplyr::filter(cluster_id_valid) %>%
        dplyr::select(-cluster_id_valid)
      
      ca_remap <- ca_exact %>%
        dplyr::filter(!cluster_id_valid) %>%
        dplyr::select(-cluster_id_valid)
      
      if (nrow(ca_remap) > 0) {
        mapped <- ca_remap %>%
          dplyr::left_join(
            cl2 %>% dplyr::rename(cluster_chr = chr, cluster_start = start, cluster_end = end, cluster_id_master = cluster_id),
            by = dplyr::join_by(chr == cluster_chr, position >= cluster_start, position <= cluster_end)
          ) %>%
          dplyr::mutate(
            cluster_id = dplyr::coalesce(cluster_id_master, cluster_id)
          ) %>%
          dplyr::select(-cluster_id_master, -cluster_chr, -cluster_start, -cluster_end)
      } else {
        mapped <- ca_remap
      }
      
      out <- dplyr::bind_rows(ca_keep, mapped) %>%
        dplyr::mutate(cluster_id = trimws(as.character(cluster_id)))
      
      n_valid <- sum(!is.na(out$cluster_id) & nzchar(out$cluster_id) & out$cluster_id %in% valid_ids, na.rm = TRUE)
      logf("[REMAP] candidates total=", nrow(out), " | valid cluster_id after remap=", n_valid)
      
      out
    }
    
    log_buf <- character(0)
    ld_log  <- shiny::reactiveVal("")
    
    append_log <- function(...) {
      txt <- paste(..., collapse = "")
      cat(txt, "\n")
      log_buf <<- c(log_buf, txt)
      ld_log(paste(log_buf, collapse = "\n"))
    }
    
    output$ld_log <- shiny::renderText(ld_log())
    
    append_log("[INIT] ld_module_server started | app_tag=", app_tag, " | app_hit_class=", app_hit_class)
    
    is_active <- shiny::reactive({
      if (is.null(activate_r)) return(TRUE)
      isTRUE(activate_r())
    })
    
    app_hit_label <- app_hit_class
    
    observe({
      append_log("[ACTIVE] ", is_active())
    })
    
    workdir <- file.path(tempdir(), paste0("ldmod_", as.integer(Sys.time())))
    dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
    
    # --------------------------------------------------------------
    # Shared RDS paths
    # --------------------------------------------------------------
    clusters_rds_path <- shiny::reactive({
      if (!isTRUE(is_active())) return(NULL)
      dd <- input$integrator_dir %||% ""
      nm <- input$clusters_rds_name %||% default_clusters_rds_name
      file.path(dd, nm)
    })
    
    candidates_rds_path <- shiny::reactive({
      if (!isTRUE(is_active())) return(NULL)
      dd <- input$integrator_dir %||% ""
      nm <- input$candidates_rds_name %||% default_candidates_rds_name
      file.path(dd, nm)
    })
    
    # --------------------------------------------------------------
    # Shared data
    # --------------------------------------------------------------
    clusters_df <- shiny::reactive({
      if (!isTRUE(is_active())) return(NULL)
      
      f <- clusters_rds_path()
      append_log("[CLUSTERS] reading path: ", f)
      append_log("[CLUSTERS] exists: ", file.exists(f))
      
      x <- tryCatch(readRDS(f), error = function(e) {
        append_log("[CLUSTERS][readRDS ERROR] ", conditionMessage(e))
        NULL
      })
      
      if (is.null(x)) {
        append_log("[CLUSTERS] readRDS returned NULL")
        return(NULL)
      }
      
      append_log("[CLUSTERS] raw class: ", paste(class(x), collapse = ", "))
      
      if (!is.data.frame(x)) {
        append_log("[CLUSTERS] raw object is not a data.frame")
        return(NULL)
      }
      
      append_log("[CLUSTERS] raw nrow=", nrow(x), " ncol=", ncol(x))
      append_log("[CLUSTERS] raw cols: ", paste(names(x), collapse = ", "))
      append_log(
        paste0(
          "[CLUSTERS] raw head:\n",
          paste(capture.output(print(utils::head(x, 5))), collapse = "\n")
        )
      )
      
      x_h <- tryCatch(
        harmonize_integrator_clusters(as.data.frame(x)),
        error = function(e) {
          append_log("[CLUSTERS][harmonize ERROR] ", conditionMessage(e))
          NULL
        }
      )
      
      if (is.null(x_h)) {
        append_log("[CLUSTERS] harmonize_integrator_clusters returned NULL")
        return(NULL)
      }
      
      append_log("[CLUSTERS] harmonized nrow=", nrow(x_h), " ncol=", ncol(x_h))
      append_log("[CLUSTERS] harmonized cols: ", paste(names(x_h), collapse = ", "))
      
      out <- tryCatch(
        prepare_clusters(x_h),
        error = function(e) {
          append_log("[CLUSTERS][prepare ERROR] ", conditionMessage(e))
          NULL
        }
      )
      
      if (is.null(out)) {
        append_log("[CLUSTERS] prepare_clusters returned NULL")
        return(NULL)
      }
      
      append_log("[CLUSTERS] prepared nrow=", nrow(out), " ncol=", ncol(out))
      append_log("[CLUSTERS] prepared cols: ", paste(names(out), collapse = ", "))
      
      out
    })
    
    candidates_df <- shiny::reactive({
      if (!isTRUE(is_active())) return(NULL)
      
      f <- candidates_rds_path()
      append_log("[CANDIDATES] reading path: ", f)
      append_log("[CANDIDATES] exists: ", file.exists(f))
      
      x <- tryCatch(readRDS(f), error = function(e) {
        append_log("[CANDIDATES][readRDS ERROR] ", conditionMessage(e))
        NULL
      })
      
      if (is.null(x)) {
        append_log("[CANDIDATES] readRDS returned NULL")
        return(NULL)
      }
      
      append_log("[CANDIDATES] raw class: ", paste(class(x), collapse = ", "))
      
      if (!is.data.frame(x)) {
        append_log("[CANDIDATES] raw object is not a data.frame")
        return(NULL)
      }
      
      append_log("[CANDIDATES] raw nrow=", nrow(x), " ncol=", ncol(x))
      append_log("[CANDIDATES] raw cols: ", paste(names(x), collapse = ", "))
      append_log(
        paste0(
          "[CANDIDATES] raw head:\n",
          paste(capture.output(print(utils::head(x, 5))), collapse = "\n")
        )
      )
      
      x_h <- tryCatch(
        harmonize_integrator_candidates(as.data.frame(x)),
        error = function(e) {
          append_log("[CANDIDATES][harmonize ERROR] ", conditionMessage(e))
          NULL
        }
      )
      
      if (is.null(x_h)) {
        append_log("[CANDIDATES] harmonize_integrator_candidates returned NULL")
        return(NULL)
      }
      
      append_log("[CANDIDATES] harmonized nrow=", nrow(x_h), " ncol=", ncol(x_h))
      append_log("[CANDIDATES] harmonized cols: ", paste(names(x_h), collapse = ", "))
      
      out <- tryCatch(
        prepare_candidates(x_h),
        error = function(e) {
          append_log("[CANDIDATES][prepare ERROR] ", conditionMessage(e))
          NULL
        }
      )
      
      if (is.null(out)) {
        append_log("[CANDIDATES] prepare_candidates returned NULL")
        return(NULL)
      }
      
      if ("classe" %in% names(out)) {
        out <- out %>%
          dplyr::mutate(classe = normalize_classe(classe))
      }
      
      cl_master <- tryCatch(clusters_df(), error = function(e) NULL)
      
      if (is.data.frame(cl_master) && nrow(cl_master)) {
        out <- tryCatch(
          reassign_candidates_to_master_clusters(out, cl_master, logger = append_log),
          error = function(e) {
            append_log("[CANDIDATES][REMAP ERROR] ", conditionMessage(e))
            out
          }
        )
      }
      
      
      append_log("[CANDIDATES] prepared nrow=", nrow(out), " ncol=", ncol(out))
      append_log("[CANDIDATES] prepared cols: ", paste(names(out), collapse = ", "))
      
      out
    })
    
    observe({
      append_log("[DBG] integrator_dir=", input$integrator_dir %||% "")
      append_log("[DBG] clusters_rds_name=", input$clusters_rds_name %||% "")
      append_log("[DBG] candidates_rds_name=", input$candidates_rds_name %||% "")
      append_log("[DBG] clusters_rds_path=", clusters_rds_path())
      append_log("[DBG] candidates_rds_path=", candidates_rds_path())
      
      cl0 <- clusters_df()
      if (is.null(cl0)) {
        append_log("[DBG] clusters_df = NULL")
      } else {
        append_log("[DBG] clusters_df OK | nrow=", nrow(cl0), " | cols=", paste(names(cl0), collapse = ","))
      }
      
      ca0 <- candidates_df()
      if (is.null(ca0)) {
        append_log("[DBG] candidates_df = NULL")
      } else {
        append_log("[DBG] candidates_df OK | nrow=", nrow(ca0), " | cols=", paste(names(ca0), collapse = ","))
      }
    })
    
    # --------------------------------------------------------------
    clusters_available <- shiny::reactive({
      cl <- clusters_df()
      ca <- candidates_df()
      
      if (!isTRUE(is_active())) {
        append_log("[AVAILABLE] inactive")
        return(NULL)
      }
      
      if (!is.data.frame(cl) || !nrow(cl)) {
        append_log("[AVAILABLE] clusters_df NULL/empty")
        return(NULL)
      }
      
      if (!is.data.frame(ca) || !nrow(ca)) {
        append_log("[AVAILABLE] candidates_df NULL/empty")
        return(NULL)
      }
      
      if (!all(c("cluster_id", "classe") %in% names(ca))) {
        append_log("[AVAILABLE] candidates_df missing cluster_id/classe")
        return(NULL)
      }
      
      app_class <- app_hit_class
      
      append_log("[AVAILABLE] app_tag=", app_tag)
      append_log("[AVAILABLE] app_class=", app_class)
      
      ca2 <- ca %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          classe     = normalize_classe(classe)
        ) %>%
        dplyr::filter(!is.na(cluster_id), nzchar(cluster_id))
      
      cl_ids <- cl %>%
        dplyr::mutate(cluster_id = trimws(as.character(cluster_id))) %>%
        dplyr::pull(cluster_id) %>%
        unique()
      
      ca2 <- ca2 %>%
        dplyr::filter(cluster_id %in% cl_ids)
      
      cls_sum <- ca2 %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          has_gwas = any(classe == "GWAS", na.rm = TRUE),
          has_app  = any(classe == app_class, na.rm = TRUE),
          classes  = paste(sort(unique(classe)), collapse = "; "),
          .groups = "drop"
        )
      
      append_log(
        paste0(
          "[AVAILABLE] cls_sum:\n",
          paste(capture.output(print(cls_sum)), collapse = "\n")
        )
      )
      
      keep_ids <- cls_sum %>%
        dplyr::filter(has_gwas, has_app) %>%
        dplyr::pull(cluster_id) %>%
        unique()
      
      append_log("[AVAILABLE] keep_ids n=", length(keep_ids))
      append_log("[AVAILABLE] keep_ids=", paste(keep_ids, collapse = ", "))
      
      out <- cl %>%
        dplyr::mutate(cluster_id = as.character(cluster_id)) %>%
        dplyr::filter(cluster_id %in% keep_ids) %>%
        dplyr::arrange(chr, start, end)
      
      append_log("[AVAILABLE] final available clusters n=", nrow(out))
      
      if (nrow(out) > 0) {
        append_log(
          paste0(
            "[AVAILABLE] final head:\n",
            paste(capture.output(print(utils::head(out, 10))), collapse = "\n")
          )
        )
      }
      
      out
    })
    
    # --------------------------------------------------------------
    # Debug tablesmake_gene_track_plot
    # --------------------------------------------------------------
    output$ld_debug_clusters_dt <- DT::renderDT({
      cl <- tryCatch(clusters_df(), error = function(e) NULL)
      if (is.null(cl)) {
        return(DT::datatable(data.frame(Message = "Clusters RDS not available"), options = list(dom = "t"), rownames = FALSE))
      }
      DT::datatable(cl, extensions = "Buttons",
                    options = list(dom = "Bfrtip", buttons = c("copy", "csv", "excel"), pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    }, server = FALSE)
    
    output$ld_debug_candidates_dt <- DT::renderDT({
      ca <- tryCatch(candidates_df(), error = function(e) NULL)
      if (is.null(ca)) {
        return(DT::datatable(data.frame(Message = "Candidates RDS not available"), options = list(dom = "t"), rownames = FALSE))
      }
      DT::datatable(ca, extensions = "Buttons",
                    options = list(dom = "Bfrtip", buttons = c("copy", "csv", "excel"), pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    }, server = FALSE)
    
    # --------------------------------------------------------------
    # Cluster selector
    # --------------------------------------------------------------
    observe({
      if (!isTRUE(is_active())) return()
      
      cl <- clusters_available()
      
      if (is.null(cl) || !is.data.frame(cl) || !nrow(cl)) {
        append_log("[UI][CHR] clusters_available is NULL/empty")
        shiny::updateSelectInput(session, "ld_chr", choices = character(0), selected = character(0))
        return()
      }
      
      append_log("[UI][CHR] clusters_available n=", nrow(cl))
      append_log("[UI][CHR] unique chr raw: ", paste(unique(cl$chr), collapse = ", "))
      
      chrs <- suppressWarnings(as.integer(cl$chr))
      chrs <- sort(unique(chrs[is.finite(chrs)]))
      
      append_log("[UI][CHR] parsed chromosomes: ", paste(chrs, collapse = ", "))
      
      if (!length(chrs)) {
        append_log("[UI][CHR] no valid chromosomes found")
        shiny::updateSelectInput(session, "ld_chr", choices = character(0), selected = character(0))
        return()
      }
      
      chr_cur <- suppressWarnings(as.integer(input$ld_chr))
      chr_sel <- if (length(chr_cur) == 1 && is.finite(chr_cur) && chr_cur %in% chrs) chr_cur else chrs[1]
      
      append_log("[UI][CHR] selected chromosome: ", chr_sel)
      
      shiny::updateSelectInput(
        session,
        "ld_chr",
        choices = as.character(chrs),
        selected = as.character(chr_sel)
      )
    })
    
    observe({
      if (!isTRUE(is_active())) return()
      
      cl <- clusters_available()
      
      if (is.null(cl) || !is.data.frame(cl) || !nrow(cl)) {
        append_log("[UI][CLUSTER] clusters_available is NULL/empty")
        shiny::updateSelectInput(session, "ld_cluster_id", choices = character(0), selected = character(0))
        return()
      }
      
      chr_sel <- suppressWarnings(as.integer(input$ld_chr))
      append_log("[UI][CLUSTER] input ld_chr = ", input$ld_chr %||% "NULL")
      
      if (!is.finite(chr_sel)) {
        chrs <- suppressWarnings(as.integer(cl$chr))
        chrs <- sort(unique(chrs[is.finite(chrs)]))
        chr_sel <- if (length(chrs)) chrs[1] else NA_integer_
      }
      
      append_log("[UI][CLUSTER] using chr_sel = ", chr_sel)
      
      cl_chr <- cl %>%
        dplyr::filter(chr == chr_sel) %>%
        dplyr::arrange(start, end)
      
      append_log("[UI][CLUSTER] n clusters for chr ", chr_sel, " = ", nrow(cl_chr))
      
      if (!nrow(cl_chr)) {
        shiny::updateSelectInput(session, "ld_cluster_id", choices = character(0), selected = character(0))
        return()
      }
      
      labs <- paste0(
        cl_chr$cluster_id,
        " (chr", chr_label_plink(cl_chr$chr), ":", cl_chr$start, "-", cl_chr$end, ")",
        " (", cl_chr$size_kb, " kb)"
      )
      choices <- stats::setNames(as.character(cl_chr$cluster_key), labs)
      
      cur <- input$ld_cluster_id %||% ""
      if (!nzchar(cur) || !(cur %in% cl_chr$cluster_key)) cur <- as.character(cl_chr$cluster_key[1])
      
      append_log("[UI][CLUSTER] selected cluster key = ", cur)
      
      shiny::updateSelectInput(session, "ld_cluster_id", choices = choices, selected = cur)
    })
    
    selected_cluster <- shiny::reactive({
      cl <- clusters_available()
      
      if (is.null(cl) || !is.data.frame(cl) || !nrow(cl)) {
        append_log("[SELECTED_CLUSTER] clusters_available NULL/empty")
        return(NULL)
      }
      
      if (is.null(input$ld_chr) || !nzchar(input$ld_chr %||% "")) {
        append_log("[SELECTED_CLUSTER] ld_chr empty")
        return(NULL)
      }
      
      if (is.null(input$ld_cluster_id) || !nzchar(input$ld_cluster_id %||% "")) {
        append_log("[SELECTED_CLUSTER] ld_cluster_id empty")
        return(NULL)
      }
      
      one <- cl %>%
        dplyr::filter(
          suppressWarnings(as.integer(chr)) == suppressWarnings(as.integer(input$ld_chr)),
          cluster_key == as.character(input$ld_cluster_id)
        )
      
      append_log("[SELECTED_CLUSTER] matches found = ", nrow(one))
      
      if (!nrow(one)) return(NULL)
      if (nrow(one) > 1) one <- one %>% dplyr::arrange(start, end) %>% dplyr::slice(1)
      
      append_log("[SELECTED_CLUSTER] cluster_id = ", one$cluster_id[1])
      one
    })
    
    selected_candidates <- shiny::reactive({
      cl <- selected_cluster()
      ca <- candidates_df()
      
      if (!is.data.frame(cl) || !nrow(cl)) return(NULL)
      if (!is.data.frame(ca) || !nrow(ca)) return(NULL)
      
      chr_sel <- as.integer(cl$chr[1])
      st <- as.integer(cl$start[1])
      en <- as.integer(cl$end[1])
      cid <- as.character(cl$cluster_id[1])
      
      ca_sub <- ca %>%
        dplyr::filter(
          chr == chr_sel,
          position >= st,
          position <= en
        )
      
      if ("cluster_id" %in% names(ca_sub)) {
        ca_sub2 <- ca_sub %>%
          dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), cluster_id == cid)
        
        if (nrow(ca_sub2) > 0) ca_sub <- ca_sub2
      }
      
      ca_sub <- ca_sub %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          rsid = trimws(as.character(rsid)),
          classe = if ("classe" %in% names(.)) normalize_classe(classe) else classe
        ) %>%
        dplyr::distinct(cluster_id, rsid, chr, position, classe, .keep_all = TRUE) %>%
        dplyr::arrange(position, rsid, classe)
      
      append_log("[SELECTED_CANDIDATES] cluster_id=", cid, " | n=", nrow(ca_sub))
      
      if (nrow(ca_sub) > 0 && "classe" %in% names(ca_sub)) {
        cls_tab <- table(ca_sub$classe, useNA = "ifany")
        append_log(
          "[SELECTED_CANDIDATES] class counts: ",
          paste(names(cls_tab), as.integer(cls_tab), sep = "=", collapse = " | ")
        )
      }
      
      ca_sub
    })
    
    # --------------------------------------------------------------
    # Population UI
    # --------------------------------------------------------------
    output$ld_pop_ui <- shiny::renderUI({
      dd <- input$ld_pops_dir %||% ""
      d  <- tryCatch(normalizePath(dd, winslash = "/", mustWork = FALSE), error = function(e) dd)
      
      ff <- if (nzchar(d) && dir.exists(d)) list.files(d, pattern="\\.[Tt][Xx][Tt]$", full.names = FALSE) else character(0)
      pops <- sub("\\.[Tt][Xx][Tt]$", "", ff)
      
      if (length(pops) == 0) {
        return(shiny::tags$div(
          style="color:#b00020; font-weight:600;",
          paste0("No .txt keep files found in: ", dd),
          shiny::tags$br(),
          shiny::tags$small(paste0("normalizePath: ", d)),
          shiny::tags$br(),
          shiny::tags$small(paste0("dir.exists: ", dir.exists(d)))
        ))
      }
      
      sel <- if (!is.null(input$ld_pop) && input$ld_pop %in% pops) input$ld_pop else {
        if ("EUR" %in% pops) "EUR" else pops[1]
      }
      shiny::selectInput(session$ns("ld_pop"), "Population (keep file)", choices = sort(pops), selected = sel)
    })
    
    # --------------------------------------------------------------
    # State
    # --------------------------------------------------------------
    ld_state <- shiny::reactiveValues(
      fl = NULL,
      ld_pairs = NULL,
      blocks_ij = NULL,
      subset_prefix = NULL,
      subset_bim = NULL,
      tag = NULL,
      chr_sel = NULL,
      st = NULL,
      en = NULL,
      seeds = NULL,
      proxy_table = NULL,
      cand_tracks = NULL
    )
    
    # --------------------------------------------------------------
    # Helper: subset reference panel for selected cluster
    # --------------------------------------------------------------
    build_subset_for_selected_cluster <- function() {
      plink_bin <- input$plink_bin %||% ""
      bfile_ref <- input$bfile_ref_prefix %||% ""
      
      shiny::validate(shiny::need(file.exists(plink_bin), paste0("PLINK binary not found: ", plink_bin)))
      shiny::validate(shiny::need(nzchar(bfile_ref), "Reference bfile prefix is empty."))
      shiny::validate(shiny::need(
        file.exists(paste0(bfile_ref, ".bed")) &&
          file.exists(paste0(bfile_ref, ".bim")) &&
          file.exists(paste0(bfile_ref, ".fam")),
        "Reference bfile prefix invalid (missing .bed/.bim/.fam)."
      ))
      
      cl <- selected_cluster()
      shiny::validate(shiny::need(is.data.frame(cl) && nrow(cl) > 0, "No selected cluster available."))
      chr_sel <- suppressWarnings(as.integer(cl$chr[1]))
      
      st      <- suppressWarnings(as.integer(cl$start[1]))
      en      <- suppressWarnings(as.integer(cl$end[1]))
      cid     <- as.character(cl$cluster_id[1])
      
      shiny::validate(shiny::need(is.finite(chr_sel) && is.finite(st) && is.finite(en) && st < en,
                                  "Invalid cluster interval (chr/start/end)."))
      
      chr_lab <- chr_label_plink(chr_sel)
      
      keep_path <- tryCatch(read_keep_file_dir(input$ld_pop, input$ld_pops_dir), error = function(e) NULL)
      shiny::validate(shiny::need(!is.null(keep_path) && file.exists(keep_path), "Keep file not found/invalid."))
      
      tag <- paste0("ld_", cid, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
      subset_prefix <- file.path(workdir, paste0(tag, "_subset"))
      
      args_subset <- c(
        "--bfile", bfile_ref,
        "--chr", chr_lab,
        "--from-bp", st,
        "--to-bp", en,
        "--keep", keep_path,
        "--make-bed",
        "--out", subset_prefix
      )
      
      append_log("[LD] Running PLINK subset for cluster ", cid, " | chr", chr_lab, ":", st, "-", en)
      r1 <- run_plink(args_subset, subset_prefix, plink_bin)
      if (length(r1$stdout)) append_log(paste(r1$stdout, collapse = "\n"))
      if (is.null(r1$status) || r1$status != 0) {
        if (length(r1$log_tail)) append_log(paste(r1$log_tail, collapse = "\n"))
        stop("PLINK subset failed.")
      }
      
      bim_file <- paste0(subset_prefix, ".bim")
      shiny::validate(shiny::need(file.exists(bim_file), "Subset .bim not created."))
      
      bim <- utils::read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
      colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
      
      fl <- bim %>%
        dplyr::transmute(SNP = as.character(SNP), BP = suppressWarnings(as.integer(BP))) %>%
        dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
        dplyr::distinct(SNP, BP, .keep_all = TRUE) %>%
        dplyr::arrange(BP, SNP) %>%
        dplyr::mutate(ix = dplyr::row_number())
      
      shiny::validate(shiny::need(nrow(fl) > 1, "Not enough SNPs in interval to compute LD (need >= 2)."))
      
      list(
        subset_prefix = subset_prefix,
        fl = fl,
        bim = fl %>% dplyr::select(SNP, BP),
        chr_sel = chr_sel,
        st = st,
        en = en,
        cid = cid,
        tag = tag
      )
    }
    
    # --------------------------------------------------------------
    # Run LD plot mode
    # --------------------------------------------------------------
    make_gene_track_plot <- function(chr_sel, st, en, fl, x_limits2, x_mode) {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomicRanges", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) return(NULL)
      if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(NULL)
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
      if (!requireNamespace("S4Vectors", quietly = TRUE)) return(NULL)
      if (!requireNamespace("IRanges", quietly = TRUE)) return(NULL)
      if (!requireNamespace("dplyr", quietly = TRUE)) return(NULL)
      if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
      
      append_log("[GENE] building collapsed gene track | chr=", chr_sel, " start=", st, " end=", en)
      
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      chr_ucsc <- paste0("chr", chr_label_plink(chr_sel))
      
      region_gr <- GenomicRanges::GRanges(
        seqnames = chr_ucsc,
        ranges   = IRanges::IRanges(start = st, end = en)
      )
      
      # -----------------------------
      # genes in region
      # -----------------------------
      genes_gr <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
      GenomeInfoDb::seqlevelsStyle(genes_gr) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
      
      ov <- GenomicRanges::findOverlaps(genes_gr, region_gr, ignore.strand = TRUE)
      keep_idx <- unique(S4Vectors::queryHits(ov))
      
      if (!length(keep_idx)) {
        append_log("[GENE] no genes overlapping region")
        return(NULL)
      }
      
      gsub_gr <- genes_gr[keep_idx]
      
      gene_id <- S4Vectors::mcols(gsub_gr)$gene_id
      if (is.null(gene_id)) gene_id <- names(gsub_gr)
      gene_id <- as.character(gene_id)
      
      gene_df <- data.frame(
        gene_id = gene_id,
        gene_start = GenomicRanges::start(gsub_gr),
        gene_end   = GenomicRanges::end(gsub_gr),
        stringsAsFactors = FALSE
      )
      
      # -----------------------------
      # labels
      # -----------------------------
      mp <- tryCatch(
        suppressMessages(
          AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys    = unique(gene_df$gene_id),
            keytype = "ENTREZID",
            columns = c("SYMBOL")
          )
        ),
        error = function(e) NULL
      )
      
      if (is.data.frame(mp) && nrow(mp)) {
        mp <- mp %>%
          dplyr::mutate(gene_id = as.character(ENTREZID)) %>%
          dplyr::select(gene_id, SYMBOL) %>%
          dplyr::distinct(gene_id, .keep_all = TRUE)
        gene_df <- gene_df %>% dplyr::left_join(mp, by = "gene_id")
      } else {
        gene_df$SYMBOL <- gene_df$gene_id
      }
      
      gene_df <- gene_df %>%
        dplyr::mutate(label = dplyr::coalesce(SYMBOL, gene_id)) %>%
        dplyr::arrange(gene_start, gene_end) %>%
        dplyr::mutate(track_y = -seq_len(dplyr::n()))
      
      append_log("[GENE] genes overlapping region n=", nrow(gene_df))
      
      # -----------------------------
      # exons by gene (collapsed)
      # -----------------------------
      ex_all <- suppressMessages(GenomicFeatures::exons(txdb, columns = c("gene_id")))
      GenomeInfoDb::seqlevelsStyle(ex_all) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
      
      ov_ex <- GenomicRanges::findOverlaps(ex_all, region_gr, ignore.strand = TRUE)
      ex_sub <- ex_all[unique(S4Vectors::queryHits(ov_ex))]
      
      if (!length(ex_sub)) {
        append_log("[GENE] no exons overlapping region")
        return(NULL)
      }
      
      ex_gene_id <- S4Vectors::mcols(ex_sub)$gene_id
      if (is.null(ex_gene_id)) {
        append_log("[GENE] exon gene_id missing")
        return(NULL)
      }
      
      ex_df <- data.frame(
        gene_id = as.character(ex_gene_id),
        exon_start = GenomicRanges::start(ex_sub),
        exon_end   = GenomicRanges::end(ex_sub),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(gene_id %in% gene_df$gene_id)
      
      if (!nrow(ex_df)) {
        append_log("[GENE] exon_df empty after gene filter")
        return(NULL)
      }
      
      # merge overlapping exons within each gene
      collapse_exons_one_gene <- function(df_one) {
        df_one <- df_one[order(df_one$exon_start, df_one$exon_end), , drop = FALSE]
        starts <- df_one$exon_start
        ends   <- df_one$exon_end
        
        out_start <- integer(0)
        out_end   <- integer(0)
        
        cur_s <- starts[1]
        cur_e <- ends[1]
        
        if (length(starts) > 1) {
          for (k in 2:length(starts)) {
            if (starts[k] <= (cur_e + 1)) {
              cur_e <- max(cur_e, ends[k])
            } else {
              out_start <- c(out_start, cur_s)
              out_end   <- c(out_end, cur_e)
              cur_s <- starts[k]
              cur_e <- ends[k]
            }
          }
        }
        
        out_start <- c(out_start, cur_s)
        out_end   <- c(out_end, cur_e)
        
        data.frame(
          gene_id = df_one$gene_id[1],
          exon_start = out_start,
          exon_end = out_end,
          stringsAsFactors = FALSE
        )
      }
      
      ex_collapsed <- do.call(
        rbind,
        lapply(split(ex_df, ex_df$gene_id), collapse_exons_one_gene)
      )
      
      ex_collapsed <- ex_collapsed %>%
        dplyr::left_join(gene_df %>% dplyr::select(gene_id, label, track_y), by = "gene_id")
      
      # -----------------------------
      # coordinate mapper
      # -----------------------------
      bp_to_x <- function(bp_vec) {
        bp_vec <- suppressWarnings(as.numeric(bp_vec))
        if (!identical(x_mode, "equal")) return(bp_vec)
        stats::approx(
          x = fl$BP,
          y = fl$X,
          xout = bp_vec,
          rule = 2,
          ties = "ordered"
        )$y
      }
      
      gene_df <- gene_df %>%
        dplyr::mutate(
          x0 = bp_to_x(gene_start),
          x1 = bp_to_x(gene_end)
        )
      
      ex_collapsed <- ex_collapsed %>%
        dplyr::mutate(
          x0 = bp_to_x(exon_start),
          x1 = bp_to_x(exon_end)
        )
      
      # -----------------------------
      # plot
      # -----------------------------
      p_gene <- plotly::plot_ly()
      
      # gene backbone
      for (i in seq_len(nrow(gene_df))) {
        p_gene <- p_gene %>%
          plotly::add_segments(
            x = gene_df$x0[i], xend = gene_df$x1[i],
            y = gene_df$track_y[i], yend = gene_df$track_y[i],
            inherit = FALSE,
            hoverinfo = "text",
            text = paste0(
              "Gene: ", gene_df$label[i],
              "<br>Start: ", gene_df$gene_start[i],
              "<br>End: ", gene_df$gene_end[i]
            ),
            showlegend = FALSE,
            line = list(width = 1)
          ) %>%
          plotly::add_text(
            x = x_limits2[2],
            y = gene_df$track_y[i],
            text = gene_df$label[i],
            textposition = "middle left",
            inherit = FALSE,
            hoverinfo = "skip",
            showlegend = FALSE
          )
      }
      
      # exon boxes
      exon_half_height <- 0.28
      
      for (i in seq_len(nrow(ex_collapsed))) {
        y0 <- ex_collapsed$track_y[i] - exon_half_height
        y1 <- ex_collapsed$track_y[i] + exon_half_height
        
        p_gene <- p_gene %>%
          plotly::add_trace(
            type = "scatter",
            mode = "lines",
            x = c(ex_collapsed$x0[i], ex_collapsed$x1[i], ex_collapsed$x1[i], ex_collapsed$x0[i], ex_collapsed$x0[i]),
            y = c(y0, y0, y1, y1, y0),
            fill = "toself",
            inherit = FALSE,
            hoverinfo = "text",
            text = paste0(
              "Gene: ", ex_collapsed$label[i],
              "<br>Exon: ", ex_collapsed$exon_start[i], "-", ex_collapsed$exon_end[i]
            ),
            showlegend = FALSE,
            line = list(width = 1)
          )
      }
      
      p_gene %>%
        plotly::layout(
          xaxis = list(
            showticklabels = FALSE,
            title = "",
            range = x_limits2
          ),
          yaxis = list(
            showticklabels = FALSE,
            title = ""
          ),
          margin = list(l = 40, r = 120, t = 10, b = 0)
        )
    }
    
    # --------------------------------------------------------------
    # Run LD plot mode
    # --------------------------------------------------------------
    observeEvent(input$run_ld_cluster, {
      if (!isTRUE(is_active())) return()
      
      ld_log("")
      append_log("[LD-PLOT] Starting...")
      
      tryCatch({
        withProgress(message = "Building LD plot...", value = 0, {
          ld_state$fl <- NULL
          ld_state$ld_pairs <- NULL
          ld_state$blocks_ij <- NULL
          ld_state$subset_prefix <- NULL
          ld_state$subset_bim <- NULL
          ld_state$tag <- NULL
          ld_state$chr_sel <- NULL
          ld_state$st <- NULL
          ld_state$en <- NULL
          ld_state$cand_tracks <- NULL
          
          incProgress(0.20, detail = "Subsetting reference panel")
          ss <- build_subset_for_selected_cluster()
          
          fl <- ss$fl
          chr_sel <- ss$chr_sel
          st <- ss$st
          en <- ss$en
          tag <- ss$tag
          subset_prefix <- ss$subset_prefix
          
          cand_sub <- tryCatch(selected_candidates(), error = function(e) NULL)
          
          if (is.data.frame(cand_sub) && nrow(cand_sub) > 0 && "classe" %in% names(cand_sub)) {
            cls_tab0 <- table(cand_sub$classe, useNA = "ifany")
            append_log(
              "[LD-PLOT] cand_sub class counts: ",
              paste(names(cls_tab0), as.integer(cls_tab0), sep = "=", collapse = " | ")
            )
            append_log(
              paste0(
                "[LD-PLOT] cand_sub head:\n",
                paste(capture.output(print(utils::head(cand_sub, 10))), collapse = "\n")
              )
            )
          } else {
            append_log("[LD-PLOT] cand_sub is NULL/empty")
          }
          
          maxn <- suppressWarnings(as.integer(input$ld_max_snps_interval %||% 400))
          shiny::validate(shiny::need(is.finite(maxn) && maxn >= 2, "Invalid max SNPs limit."))
          
          if (nrow(fl) > maxn) {
            append_log("[LD-PLOT] Too many SNPs: ", nrow(fl), " > ", maxn, ". Auto-thinning...")
            
            keep_snps <- character(0)
            if (is.data.frame(cand_sub) && nrow(cand_sub) > 0) {
              keep_by_rsid <- intersect(fl$SNP, cand_sub$rsid)
              keep_by_pos  <- fl %>% dplyr::filter(BP %in% cand_sub$position) %>% dplyr::pull(SNP)
              keep_snps <- unique(c(keep_by_rsid, keep_by_pos))
            }
            
            shiny::validate(shiny::need(
              length(keep_snps) <= maxn,
              paste0("Candidates in interval (", length(keep_snps), ") exceed max limit (", maxn, "). Increase limit.")
            ))
            
            k <- maxn - length(keep_snps)
            rest <- fl %>% dplyr::filter(!(SNP %in% keep_snps)) %>% dplyr::arrange(BP, SNP)
            
            if (k > 0 && nrow(rest) > 0) {
              if (nrow(rest) <= k) {
                fill_snps <- rest$SNP
              } else {
                idx <- unique(pmax(1, pmin(nrow(rest), round(seq(1, nrow(rest), length.out = k)))))
                fill_snps <- rest$SNP[idx]
              }
            } else {
              fill_snps <- character(0)
            }
            
            snps_final <- unique(c(keep_snps, fill_snps))
            extract_file <- file.path(workdir, paste0(tag, "_extract_snps.txt"))
            writeLines(snps_final, extract_file)
            
            subset2_prefix <- file.path(workdir, paste0(tag, "_subset_thin"))
            args_thin <- c("--bfile", subset_prefix, "--extract", extract_file, "--make-bed", "--out", subset2_prefix)
            rthin <- run_plink(args_thin, subset2_prefix, input$plink_bin %||% "")
            if (length(rthin$stdout)) append_log(paste(rthin$stdout, collapse = "\n"))
            if (is.null(rthin$status) || rthin$status != 0) stop("PLINK thinning subset failed.")
            
            subset_prefix <- subset2_prefix
            
            bim2 <- utils::read.table(paste0(subset_prefix, ".bim"), header = FALSE, stringsAsFactors = FALSE)
            colnames(bim2) <- c("CHR","SNP","CM","BP","A1","A2")
            
            fl <- bim2 %>%
              dplyr::transmute(SNP = as.character(SNP), BP = suppressWarnings(as.integer(BP))) %>%
              dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
              dplyr::distinct(SNP, BP, .keep_all = TRUE) %>%
              dplyr::arrange(BP, SNP) %>%
              dplyr::mutate(ix = dplyr::row_number())
            
            append_log("[LD-PLOT] SNPs after thinning: ", nrow(fl))
          }
          
          incProgress(0.60, detail = "Computing full LD")
          ld_prefix <- file.path(workdir, paste0(tag, "_ld"))
          r2_args <- c(
            "--bfile", subset_prefix,
            "--ld-window", "999999",
            "--ld-window-kb", "999999",
            "--ld-window-r2", "0"
          )
          
          if (identical(input$ld_metric, "Dprime")) {
            r2_args <- c(r2_args, "--r2", "dprime", "gz")
          } else {
            r2_args <- c(r2_args, "--r2", "gz")
          }
          
          r2_args <- c(r2_args, "--out", ld_prefix)
          
          r2 <- run_plink(r2_args, ld_prefix, input$plink_bin %||% "")
          if (length(r2$stdout)) append_log(paste(r2$stdout, collapse = "\n"))
          if (is.null(r2$status) || r2$status != 0) stop("PLINK LD failed.")
          
          incProgress(0.80, detail = "Parsing LD")
          ld <- read_plink_ld(ld_prefix)
          snpA <- pick_col(ld, c("SNP_A","SNP_A1","SNP1"))
          snpB <- pick_col(ld, c("SNP_B","SNP_B1","SNP2"))
          vcol <- if (identical(input$ld_metric, "Dprime")) pick_col(ld, c("Dprime","D'","DP","DPRIME")) else pick_col(ld, c("R2","r2"))
          
          shiny::validate(shiny::need(!is.null(snpA) && !is.null(snpB) && !is.null(vcol), "LD file missing required columns."))
          
          fl_map <- fl %>% dplyr::select(SNP, BP, ix)
          
          ld2 <- ld %>%
            dplyr::transmute(
              SNP_A = as.character(.data[[snpA]]),
              SNP_B = as.character(.data[[snpB]]),
              value = suppressWarnings(as.numeric(.data[[vcol]]))
            ) %>%
            dplyr::filter(!is.na(SNP_A), !is.na(SNP_B), is.finite(value)) %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixA = ix, BPA = BP), by = c("SNP_A" = "SNP")) %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixB = ix, BPB = BP), by = c("SNP_B" = "SNP")) %>%
            dplyr::filter(is.finite(ixA), is.finite(ixB), ixA != ixB) %>%
            dplyr::mutate(i = pmax(ixA, ixB), j = pmin(ixA, ixB)) %>%
            dplyr::distinct(i, j, .keep_all = TRUE) %>%
            dplyr::filter(i > j)
          
          fl2 <- fl %>%
            dplyr::mutate(is_candidate = FALSE, classe = NA_character_) %>%
            dplyr::select(SNP, BP, ix, is_candidate, classe)
          
          cand_tracks <- tibble::tibble(
            track_id = character(),
            BP = integer(),
            classe = character(),
            label_id = character(),
            cluster_id = character()
          )
          
          if (is.data.frame(cand_sub) && nrow(cand_sub) > 0) {
            
            cand_tracks <- cand_sub %>%
              dplyr::mutate(
                cluster_id = as.character(cluster_id),
                BP         = suppressWarnings(as.integer(position)),
                classe     = normalize_classe(classe),
                label_id   = as.character(rsid),
                label_id   = trimws(label_id),
                label_id   = dplyr::if_else(
                  is.na(label_id) | !nzchar(label_id),
                  as.character(cluster_id),
                  label_id
                ),
                label_id   = trimws(label_id),
                label_id   = dplyr::if_else(!nzchar(label_id), paste0("hit_", dplyr::row_number()), label_id),
                track_id   = paste0(classe, "::", label_id, "::", BP)
              ) %>%
              dplyr::filter(
                !is.na(BP),
                is.finite(BP),
                BP >= st,
                BP <= en,
                !is.na(classe),
                nzchar(classe)
              ) %>%
              dplyr::distinct(track_id, .keep_all = TRUE) %>%
              dplyr::select(track_id, BP, classe, label_id, cluster_id) %>%
              dplyr::arrange(BP, classe, label_id)
            
            # Mantinc fl2 sense exigir match exacte amb SNP del panell.
            # Només marquem SNPs com a candidate si coincideixen exactament per posició.
            if (nrow(cand_tracks) > 0) {
              fl2 <- fl2 %>%
                dplyr::mutate(
                  is_candidate = BP %in% cand_tracks$BP
                ) %>%
                dplyr::left_join(
                  cand_tracks %>%
                    dplyr::group_by(BP) %>%
                    dplyr::summarise(
                      classe = paste(sort(unique(classe)), collapse = "; "),
                      .groups = "drop"
                    ),
                  by = "BP"
                ) %>%
                dplyr::mutate(classe = dplyr::coalesce(classe.y, classe.x)) %>%
                dplyr::select(SNP, BP, ix, is_candidate, classe)
            }
          }
          
          ld_state$cand_tracks <- cand_tracks
          
          append_log("[LD-PLOT] cand_tracks n=", nrow(cand_tracks))
          if (nrow(cand_tracks) > 0) {
            cls_tab2 <- table(cand_tracks$classe, useNA = "ifany")
            append_log(
              "[LD-PLOT] cand_tracks class counts: ",
              paste(names(cls_tab2), as.integer(cls_tab2), sep = "=", collapse = " | ")
            )
            append_log(
              paste0(
                "[LD-PLOT] cand_tracks head:\n",
                paste(capture.output(print(utils::head(cand_tracks, 10))), collapse = "\n")
              )
            )
          } else {
            append_log("[LD-PLOT] cand_tracks is empty")
          }
          
          ld_state$fl <- fl2
          ld_state$ld_pairs <- ld2
          ld_state$subset_prefix <- subset_prefix
          ld_state$subset_bim <- fl %>% dplyr::select(SNP, BP)
          ld_state$tag <- tag
          ld_state$chr_sel <- chr_sel
          ld_state$st <- st
          ld_state$en <- en
          
          incProgress(1, detail = "Done")
          append_log("[LD-PLOT] Done")
        })
      }, error = function(e) {
        append_log("[LD-PLOT][ERROR] ", conditionMessage(e))
        showNotification(paste("LD plot error:", conditionMessage(e)), type = "error", duration = NULL)
      })
    }, ignoreInit = TRUE)
    
    # --------------------------------------------------------------
    # Run blocks
    # --------------------------------------------------------------
    observeEvent(input$run_ld_blocks, {
      if (!isTRUE(is_active())) return()
      
      append_log("[BLOCKS] Starting...")
      
      tryCatch({
        withProgress(message = "Building LD blocks...", value = 0, {
          shiny::validate(shiny::need(!is.null(ld_state$subset_prefix), "Run LD plot first."))
          
          subset_prefix <- ld_state$subset_prefix
          fl <- ld_state$fl
          tag <- ld_state$tag %||% format(Sys.time(), "%Y%m%d_%H%M%S")
          
          blk_prefix <- file.path(workdir, paste0(tag, "_blocks"))
          blk_args <- c(
            "--bfile", subset_prefix,
            "--blocks", "no-pheno-req",
            if (isTRUE(input$ld_blocks_no_small_max_span)) "no-small-max-span" else NULL,
            "--blocks-inform-frac", as.character(input$ld_blocks_inform_frac %||% 0.60),
            "--blocks-max-kb", as.character(input$ld_blocks_max_kb %||% 1000),
            "--blocks-min-maf", as.character(input$ld_blocks_min_maf %||% 0.05),
            "--blocks-strong-highci", as.character(input$ld_blocks_strong_highci %||% 0.90),
            "--blocks-strong-lowci", as.character(input$ld_blocks_strong_lowci %||% 0.55),
            "--out", blk_prefix
          )
          blk_args <- blk_args[!is.null(blk_args) & nzchar(as.character(blk_args))]
          
          rb <- run_plink(blk_args, blk_prefix, input$plink_bin %||% "")
          if (length(rb$stdout)) append_log(paste(rb$stdout, collapse = "\n"))
          if (is.null(rb$status) || rb$status != 0) stop("PLINK --blocks failed.")
          
          br <- read_plink_blocks(blk_prefix)
          ij <- blocks_det_to_ij(br, fl %>% dplyr::select(SNP, BP))
          ld_state$blocks_ij <- ij
          
          append_log("[BLOCKS] Done | n=", nrow(ij))
        })
      }, error = function(e) {
        append_log("[BLOCKS][ERROR] ", conditionMessage(e))
        showNotification(paste("Blocks error:", conditionMessage(e)), type = "error", duration = NULL)
      })
    }, ignoreInit = TRUE)
    
    # --------------------------------------------------------------
    # Run prioritization
    # --------------------------------------------------------------
    observeEvent(input$run_ld_prior, {
      if (!isTRUE(is_active())) return()
      
      append_log("[LD-PRIOR] Starting...")
      
      tryCatch({
        withProgress(message = "Building LD prioritization...", value = 0, {
          
          ld_state$seeds <- NULL
          ld_state$proxy_table <- NULL
          
          incProgress(0.20, detail = "Subsetting reference panel")
          ss <- build_subset_for_selected_cluster()
          subset_prefix <- ss$subset_prefix
          bim_df <- ss$bim
          cid <- ss$cid
          
          hits_sub_all <- tryCatch(selected_candidates(), error = function(e) NULL)
          shiny::validate(shiny::need(
            is.data.frame(hits_sub_all) && nrow(hits_sub_all) > 0,
            "No candidates available in selected cluster."
          ))
          
          hits_sub <- hits_sub_all %>%
            dplyr::mutate(classe = normalize_classe(classe)) %>%
            dplyr::filter(classe == "GWAS")
          
          shiny::validate(shiny::need(
            nrow(hits_sub) > 0,
            "No GWAS hits available in selected cluster."
          ))
          
          append_log("[LD-PRIOR] cluster_id=", cid, " | GWAS query_hits=", nrow(hits_sub))
          
          incProgress(0.40, detail = "Finding seeds")
          seeds_df <- hits_sub %>%
            dplyr::rowwise() %>%
            dplyr::do({
              seed <- find_seed_snp(.$rsid, .$position, bim_df)
              dplyr::bind_cols(
                data.frame(
                  cluster_id = as.character(.$cluster_id),
                  query_hit = as.character(.$rsid),
                  query_pos = as.integer(.$position),
                  stringsAsFactors = FALSE
                ),
                seed
              )
            }) %>%
            dplyr::ungroup() %>%
            dplyr::distinct()
          
          append_log("[LD-PRIOR] seeds nrow=", nrow(seeds_df))
          
          incProgress(0.70, detail = "Computing proxies")
          proxy_list <- lapply(seq_len(nrow(hits_sub)), function(i) {
            hit_row <- hits_sub[i, , drop = FALSE]
            append_log("[LD-PRIOR] hit ", i, "/", nrow(hits_sub), " | ", hit_row$rsid[1])
            tryCatch(
              build_hit_proxy_table(
                hit_row = hit_row,
                bim_df = bim_df,
                subset_prefix = subset_prefix,
                plink_bin = input$plink_bin %||% "",
                workdir = workdir,
                metric = input$ld_metric %||% "R2"
              ),
              error = function(e) {
                append_log("[LD-PRIOR][WARN] ", hit_row$rsid[1], " -> ", conditionMessage(e))
                NULL
              }
            )
          })
          
          proxy_tbl <- dplyr::bind_rows(proxy_list)
          if (!is.data.frame(proxy_tbl) || nrow(proxy_tbl) == 0) {
            proxy_tbl <- tibble::tibble(
              cluster_id = character(0),
              query_hit = character(0),
              query_pos = integer(0),
              seed_snp = character(0),
              seed_pos = integer(0),
              seed_type = character(0),
              seed_dist_bp = integer(0),
              proxy_snp = character(0),
              proxy_pos = integer(0),
              ld_value = numeric(0)
            )
          }
          
          r2_min <- suppressWarnings(as.numeric(input$ld_proxy_r2_min %||% 0.6))
          proxy_tbl <- proxy_tbl %>%
            dplyr::filter(is.finite(ld_value)) %>%
            dplyr::filter(ld_value >= r2_min) %>%
            dplyr::arrange(query_hit, dplyr::desc(ld_value), proxy_pos, proxy_snp)
          
          append_log("[LD-PRIOR] proxy_table nrow=", nrow(proxy_tbl), " | r2_min=", r2_min)
          
          ld_state$subset_prefix <- subset_prefix
          ld_state$subset_bim <- bim_df
          ld_state$seeds <- seeds_df
          ld_state$proxy_table <- proxy_tbl
          
          incProgress(1, detail = "Done")
          append_log("[LD-PRIOR] Done")
        })
      }, error = function(e) {
        append_log("[LD-PRIOR][ERROR] ", conditionMessage(e))
        showNotification(paste("LD prioritization error:", conditionMessage(e)), type = "error", duration = NULL)
      })
    }, ignoreInit = TRUE)
    
    # --------------------------------------------------------------
    # Tables
    # --------------------------------------------------------------
    output$ld_fl_dt <- DT::renderDT({
      fl <- ld_state$fl
      if (is.null(fl) || !nrow(fl)) {
        return(DT::datatable(data.frame(Message = "Run LD plot to populate SNP list."), options = list(dom = "t"), rownames = FALSE))
      }
      
      show <- fl %>%
        dplyr::mutate(
          candidate = dplyr::if_else(isTRUE(is_candidate) | (is.logical(is_candidate) & is_candidate), "YES", ""),
          classe = dplyr::coalesce(as.character(classe), "")
        ) %>%
        dplyr::transmute(ix, SNP, BP, candidate, classe)
      
      DT::datatable(
        show,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          pageLength = 15,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_pairs_dt <- DT::renderDT({
      lp <- ld_state$ld_pairs
      fl <- ld_state$fl
      if (is.null(lp) || !nrow(lp) || is.null(fl) || !nrow(fl)) {
        return(DT::datatable(data.frame(Message = "Run LD plot to populate pairs."), options = list(dom = "t"), rownames = FALSE))
      }
      
      fl_small <- fl %>% dplyr::select(ix, SNP, BP, is_candidate, classe) %>% dplyr::distinct(ix, .keep_all = TRUE)
      
      top <- lp %>%
        dplyr::left_join(
          fl_small %>% dplyr::rename(i = ix, SNP_i = SNP, BP_i = BP, cand_i = is_candidate, class_i = classe),
          by = "i"
        ) %>%
        dplyr::left_join(
          fl_small %>% dplyr::rename(j = ix, SNP_j = SNP, BP_j = BP, cand_j = is_candidate, class_j = classe),
          by = "j"
        ) %>%
        dplyr::arrange(dplyr::desc(value))
      
      DT::datatable(
        top,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          pageLength = 12,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_seed_dt <- DT::renderDT({
      dd <- ld_state$seeds
      if (is.null(dd) || !nrow(dd)) {
        return(DT::datatable(data.frame(Message = "Run LD prioritization to populate seeds."), options = list(dom = "t"), rownames = FALSE))
      }
      
      DT::datatable(
        dd,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          pageLength = 12,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_priority_dt <- DT::renderDT({
      dd <- ld_state$proxy_table
      if (is.null(dd) || !nrow(dd)) {
        return(DT::datatable(data.frame(Message = "Run LD prioritization to populate proxy table."), options = list(dom = "t"), rownames = FALSE))
      }
      
      DT::datatable(
        dd,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          pageLength = 15,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    # --------------------------------------------------------------
    # Plot
    # --------------------------------------------------------------
    output$ld_plot <- plotly::renderPlotly({
      fl0   <- ld_state$fl
      ld    <- ld_state$ld_pairs
      cand0 <- ld_state$cand_tracks
      
      shiny::validate(
        shiny::need(is.data.frame(fl0) && nrow(fl0) >= 2, "Run LD plot first."),
        shiny::need(is.data.frame(ld) && nrow(ld) > 0, "No LD pairs available.")
      )
      
      metric <- input$ld_metric %||% "R2"
      x_mode <- input$ld_x_mode %||% "bp"
      
      shiny::validate(shiny::need(metric %in% c("R2","Dprime"), "ld_metric must be R2 or Dprime."))
      shiny::validate(shiny::need(x_mode %in% c("bp","equal"), "ld_x_mode must be bp or equal."))
      
      append_log("[LD-PLOT][RENDER] entered")
      
      fl <- fl0 %>%
        dplyr::mutate(
          SNP = as.character(SNP),
          BP  = suppressWarnings(as.integer(BP)),
          is_candidate = dplyr::coalesce(as.logical(is_candidate), FALSE),
          classe = dplyr::coalesce(as.character(classe), "")
        ) %>%
        dplyr::distinct(SNP, .keep_all = TRUE) %>%
        dplyr::arrange(BP, SNP)
      
      if (identical(x_mode, "equal")) {
        fl <- fl %>% dplyr::mutate(X = dplyr::row_number())
        x_title <- "Equal spacing"
      } else {
        fl <- fl %>% dplyr::mutate(X = BP)
        x_title <- "Genomic (BP)"
      }
      
      xmap   <- stats::setNames(fl$X, fl$SNP)
      bpmap  <- stats::setNames(fl$BP, fl$SNP)
      idxmap <- stats::setNames(seq_len(nrow(fl)), fl$SNP)
      
      x_limits <- range(fl$X, na.rm = TRUE)
      x_span   <- max(1, diff(x_limits))
      
      bp_span <- diff(range(fl$BP, na.rm = TRUE))
      if (!is.finite(bp_span) || bp_span <= 0) bp_span <- 1
      
      tri <- ld %>%
        dplyr::transmute(
          SNP_A = as.character(SNP_A),
          SNP_B = as.character(SNP_B),
          val   = suppressWarnings(as.numeric(value))
        ) %>%
        dplyr::filter(is.finite(val), SNP_A %in% fl$SNP, SNP_B %in% fl$SNP) %>%
        dplyr::mutate(
          X_A = unname(xmap[SNP_A]),
          X_B = unname(xmap[SNP_B]),
          i   = unname(idxmap[SNP_A]),
          j   = unname(idxmap[SNP_B]),
          bpA = unname(bpmap[SNP_A]),
          bpB = unname(bpmap[SNP_B]),
          x_mid = (X_A + X_B) / 2,
          y = if (identical(x_mode, "equal")) {
            -(abs(j - i) / max(1, nrow(fl) - 1))
          } else {
            -(abs(bpB - bpA) / bp_span)
          },
          hover = paste0(
            "SNP_A: ", SNP_A,
            "<br>SNP_B: ", SNP_B,
            "<br>", metric, ": ", signif(val, 3),
            "<br>BP_A: ", bpA,
            "<br>BP_B: ", bpB
          )
        ) %>%
        dplyr::filter(is.finite(x_mid), is.finite(y))
      
      shiny::validate(shiny::need(nrow(tri) > 0, "No LD pairs to plot."))
      
      cand <- if (is.data.frame(cand0) && nrow(cand0)) {
        
        bp_to_x_track <- function(bp_vec, fl_df, x_mode_cur) {
          bp_vec <- suppressWarnings(as.numeric(bp_vec))
          
          if (!identical(x_mode_cur, "equal")) {
            return(bp_vec)
          }
          
          stats::approx(
            x = fl_df$BP,
            y = fl_df$X,
            xout = bp_vec,
            rule = 2,
            ties = "ordered"
          )$y
        }
        
        cand0 %>%
          dplyr::mutate(
            BP = suppressWarnings(as.integer(BP)),
            classe = normalize_classe(classe),
            label_id = as.character(label_id),
            X = bp_to_x_track(BP, fl, x_mode)
          ) %>%
          dplyr::filter(is.finite(BP), is.finite(X)) %>%
          dplyr::mutate(
            hover = paste0(
              "ID: ", label_id,
              "<br>BP: ", BP,
              "<br>classe: ", classe
            )
          ) %>%
          dplyr::distinct(track_id, .keep_all = TRUE)
        
      } else {
        tibble::tibble()
      }
      
      append_log("[LD-PLOT][RENDER] cand nrow=", nrow(cand))
      if (nrow(cand) > 0) {
        cls_tab3 <- table(cand$classe, useNA = "ifany")
        append_log(
          "[LD-PLOT][RENDER] classes: ",
          paste(names(cls_tab3), as.integer(cls_tab3), sep = "=", collapse = " | ")
        )
      }
      
      cand_classes <- unique(cand$classe)
      cand_classes <- cand_classes[!is.na(cand_classes) & nzchar(cand_classes)]
      if (!length(cand_classes)) cand_classes <- "candidate"
      
      preferred_order <- c("GWAS", app_hit_class, "candidate")
      cand_classes <- unique(c(
        preferred_order[preferred_order %in% cand_classes],
        cand_classes[!(cand_classes %in% preferred_order)]
      ))
      
      track_base_y <- seq(0.04, by = 0.10, length.out = length(cand_classes))
      names(track_base_y) <- cand_classes
      
      if (nrow(cand) > 0) {
        cand <- cand %>%
          dplyr::mutate(
            track_y0 = unname(track_base_y[classe]),
            track_y1 = track_y0 + 0.06
          )
      }
      
      x_lab <- if (identical(x_mode, "equal")) {
        x_limits[2] + max(2, ceiling(0.08 * nrow(fl)))
      } else {
        x_limits[2] + 0.08 * x_span
      }
      x_limits2 <- c(x_limits[1], x_lab)
      
      cl <- selected_cluster()
      g0 <- NULL
      if (is.data.frame(cl) && nrow(cl)) {
        g0 <- tryCatch(
          make_gene_track_plot(
            chr_sel   = as.integer(cl$chr[1]),
            st        = as.integer(cl$start[1]),
            en        = as.integer(cl$end[1]),
            fl        = fl,
            x_limits2 = x_limits2,
            x_mode    = x_mode
          ),
          error = function(e) {
            append_log("[GENE][ERROR] ", conditionMessage(e))
            NULL
          }
        )
      }
      
      p <- plotly::plot_ly()
      
      # triangle
      p <- p %>%
        plotly::add_trace(
          data = tri,
          x = ~x_mid,
          y = ~y,
          type = "scatter",
          mode = "markers",
          text = ~hover,
          hoverinfo = "text",
          marker = list(
            symbol = "square",
            size = 9,
            color = tri$val,
            cmin = 0,
            cmax = 1,
            colorscale = list(
              c(0.0, "white"),
              c(0.2, "orange"),
              c(1.0, "red")
            ),
            colorbar = list(title = metric)
          ),
          showlegend = FALSE
        )
      
      # tracks
      if (nrow(cand) > 0) {
        for (cls in cand_classes) {
          dsub <- cand %>% dplyr::filter(classe == cls)
          if (!nrow(dsub)) next
          
          y0 <- track_base_y[[cls]]
          y1 <- y0 + 0.06
          
          p <- p %>%
            plotly::add_segments(
              x = x_limits[1],
              xend = x_limits[2],
              y = y0,
              yend = y0,
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE,
              line = list(width = 1)
            ) %>%
            plotly::add_segments(
              data = dsub,
              x = ~X,
              xend = ~X,
              y = y0,
              yend = y1,
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE,
              line = list(width = 1)
            ) %>%
            plotly::add_markers(
              data = dsub,
              x = ~X,
              y = ~track_y1,
              text = ~hover,
              hoverinfo = "text",
              name = cls,
              marker = list(size = 7),
              showlegend = FALSE
            ) %>%
            plotly::add_text(
              x = x_lab,
              y = y1,
              text = cls,
              textposition = "middle left",
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE
            )
        }
      }
      
      # blocks
      bd <- NULL
      if (is.data.frame(ld_state$blocks_ij) && nrow(ld_state$blocks_ij) > 0) {
        bd <- ld_state$blocks_ij %>%
          dplyr::mutate(
            i = suppressWarnings(as.integer(i)),
            j = suppressWarnings(as.integer(j))
          ) %>%
          dplyr::filter(is.finite(i), is.finite(j), i >= 1, j <= nrow(fl), j > i)
        
        if (nrow(bd) > 0) {
          bd <- bd %>%
            dplyr::mutate(
              xL  = fl$X[i],
              xR  = fl$X[j],
              xM  = (xL + xR) / 2,
              bpL = fl$BP[i],
              bpR = fl$BP[j]
            ) %>%
            dplyr::filter(is.finite(xL), is.finite(xR), is.finite(xM), is.finite(bpL), is.finite(bpR))
          
          if (nrow(bd) > 0) {
            y0_blk <- -0.01
            
            if (identical(x_mode, "equal")) {
              denom <- max(1, nrow(fl) - 1)
              bd <- bd %>% dplyr::mutate(yA = -(abs(j - i) / denom))
            } else {
              bd <- bd %>% dplyr::mutate(yA = -(abs(bpR - bpL) / bp_span))
            }
            
            bd <- bd %>% dplyr::mutate(yA = pmin(yA, y0_blk - 0.02))
            
            p <- p %>%
              plotly::add_segments(
                data = bd,
                x = ~xL, xend = ~xR,
                y = y0_blk, yend = y0_blk,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              ) %>%
              plotly::add_segments(
                data = bd,
                x = ~xL, xend = ~xM,
                y = y0_blk, yend = ~yA,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              ) %>%
              plotly::add_segments(
                data = bd,
                x = ~xR, xend = ~xM,
                y = y0_blk, yend = ~yA,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              )
          }
        }
      }
      
      ymax_tracks <- if (length(track_base_y)) max(track_base_y) + 0.12 else 0.25
      
      p <- p %>%
        plotly::layout(
          hovermode = "closest",
          xaxis = list(
            title = x_title,
            automargin = TRUE,
            range = x_limits2
          ),
          yaxis = list(
            title = "",
            range = c(-1.05, max(0.25, ymax_tracks))
          ),
          margin = list(l = 40, r = 120, t = 20, b = 60)
        )
      
      if (!is.null(g0)) {
        plotly::subplot(g0, p, nrows = 2, heights = c(0.20, 0.80), shareX = TRUE, titleX = TRUE) %>%
          plotly::layout(hovermode = "closest")
      } else {
        p
      }
    })
    
    # Optional return
    invisible(list(
      clusters_df = clusters_df,
      candidates_df = candidates_df,
      selected_cluster = selected_cluster,
      selected_candidates = selected_candidates,
      seed_table = shiny::reactive(ld_state$seeds),
      proxy_table = shiny::reactive(ld_state$proxy_table)
    ))
  })
}

############################
