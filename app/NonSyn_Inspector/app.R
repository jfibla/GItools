# app.R — NonSyn Inspector (GWAS + dbSNP intervals + dbNSFP non-synonymous)
# revised
# LD, GO, KEGG
# include clustering by hits
# annotation by .out file

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
library(fmsb)
library(shinycssloaders)
library(patchwork)

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
# Runtime logging (early)
# ===========================
APP_KEY <- "nonsyn"  # <- CANVIA: "catalog" / "ewasdis" / "ewastum" / "nonsyn"

LOG_DIR <- Sys.getenv("GITOOLS_LOG_DIR", "")
if (!nzchar(LOG_DIR)) {
  LOG_DIR <- file.path(dirname(getwd()), "_logs")  # -> GItools/app/_logs
}
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

APP_LOG <- file.path(LOG_DIR, paste0(APP_KEY, "_runtime.log"))

applog <- function(...) {
  cat(paste0(..., collapse = ""), "\n", file = APP_LOG, append = TRUE)
}

applog("[", APP_KEY, "] boot @ ", format(Sys.time()), " wd=", getwd(),
       " pid=", Sys.getpid(),
       " sid=", Sys.getenv("GITOOLS_SID", unset=""),
       " url_mode=", Sys.getenv("GITOOLS_URL_MODE", unset=""))

# IMPORTANT: captura el missatge REAL de l'error (no només calls)
options(shiny.error = function() {
  e <- geterrmessage()
  applog("[", APP_KEY, "][ERROR] shiny.error @ ", format(Sys.time()))
  applog("[", APP_KEY, "][ERROR] message: ", e)
  
  # traceback útil
  tb <- paste(utils::capture.output(traceback()), collapse = "\n")
  applog("[", APP_KEY, "][ERROR] traceback:\n", tb)
  
  # calls útils
  cl <- paste(utils::capture.output(sys.calls()), collapse = "\n")
  applog("[", APP_KEY, "][ERROR] calls:\n", cl)
})


# ============================
# Helpers (out of UI & SERVER)
# ============================

# --- Fix dplyr verbs ---
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
rename    <- dplyr::rename

# ---------------------------
# GItools portable config
# ---------------------------
cfg_file <- Sys.getenv("GITOOLS_CONFIG", unset = "")
if (nzchar(cfg_file) && file.exists(cfg_file)) {
  source(cfg_file, local = TRUE)
} else {
  # Busca config.R pujant directoris
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  found <- ""
  for (i in 1:10) {
    cand <- file.path(d, "config.R")
    if (file.exists(cand)) { found <- cand; break }
    d2 <- dirname(d)
    if (identical(d2, d)) break
    d <- d2
  }
  if (!nzchar(found)) stop("config.R not found. Set env GITOOLS_CONFIG or run inside GItools repo.")
  source(found, local = TRUE)
}

cfg <- gi_cfg()

ROOT   <- cfg$root
SHARED <- cfg$shared
RES    <- cfg$resources

# Shared canonical engine / state
source(file.path(SHARED, "gi_state.R"), local = TRUE)
source(file.path(SHARED, "gi_slave_canonical.R"), local = TRUE)
source(file.path(SHARED, "gi_clusters_canonical.R"), local = TRUE)

message("[NonSyn] sourced gi_slave_canonical from: ", file.path(cfg$shared, "gi_slave_canonical.R"))

PRED_UNIV_RDS <- file.path(gi_cfg()$shared, "predictor_universe.rds")
bg <- readRDS(PRED_UNIV_RDS)

# ---------------------------
# NonSyn defaults (portable)
# ---------------------------
gi_nonsyn_defaults <- function() {
  res <- gi_cfg()$resources
  
  # dbNSFP ALL
  dbnsfp_all <- Sys.getenv(
    "GITOOLS_DBNSFP_ALL",
    unset = if (nzchar(res)) file.path(res, "dbNSFP5_broad_hg38", "output_snps_biallelic_dbNSFP5_broad_hg38.ALL.out.gz") else ""
  )
  
  # PLINK 1.9 path (try common locations)
  plink <- Sys.getenv("GITOOLS_PLINK19", unset = "")
  if (!nzchar(plink) && nzchar(res)) {
    cand1 <- file.path(res, "software", "plink19", "plink")
    cand2 <- file.path(res, "software", "plink19", "plink19")
    cand3 <- file.path(res, "software", "plink19")  # legacy (single binary)
    plink <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else cand3
  }
  
  # LD reference bfile prefix
  bfile <- Sys.getenv(
    "GITOOLS_LD_BFILE",
    unset = if (nzchar(res)) file.path(res, "LD_resources", "Merged_FULL_SET_hg38_hgdp.wgs_10000G3_ko07_MAF0.05") else ""
  )
  
  # POP keep files dir
  popdir <- Sys.getenv(
    "GITOOLS_POP_DIR",
    unset = if (nzchar(res)) file.path(res, "LD_resources", "POP") else ""
  )
  
  list(res = res, dbnsfp_all = dbnsfp_all, plink = plink, bfile = bfile, popdir = popdir)
}

ns_def <- gi_nonsyn_defaults()

# -----------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

plotly_message <- function(msg) {
  plotly::plot_ly() %>%
    plotly::layout(
      title = list(
        text   = paste0("<b>", msg, "</b>"),
        x      = 0.5,
        xanchor = "center",
        font   = list(size = 12)
      ),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

# PLINK-like chr map (accepts "1..22,X,Y,MT", returns 1..22,23,24,26)
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


summarize_clusters_from_gwas <- function(merged_intervals, gwas_df) {
  if (is.null(merged_intervals) || !nrow(merged_intervals)) return(tibble())
  df <- gwas_df
  
  req_cols <- c("CHR","BP","snp","Pval","logp")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop(paste0("GWAS missing columns: ", paste(miss, collapse=", ")))
  
  cl <- merged_intervals %>%
    dplyr::arrange(chr, start, end) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(cluster_n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_id = paste0("cluster_chr", chr_label_plink(as.integer(chr)), "_", cluster_n))
  
  cl_sum <- purrr::pmap_dfr(
    list(cl$chr, cl$start, cl$end, cl$cluster_id, cl$cluster_n),
    function(chr, st, en, cid, cn) {
      sub <- df %>% dplyr::filter(as.integer(CHR) == as.integer(chr), BP >= st, BP <= en)
      if (!nrow(sub)) {
        tibble(chr=chr, cluster_n=cn, cluster_id=cid,
               start=st, end=en, center=round((st+en)/2),
               n_snps=0L, top_snp=NA_character_, top_logp=NA_real_,
               cluster_size_kb=round((en-st)/1000,2))
      } else {
        ii <- which.max(sub$logp)
        tibble(chr=chr, cluster_n=cn, cluster_id=cid,
               start=st, end=en, center=round(mean(sub$BP, na.rm=TRUE)),
               n_snps=nrow(sub), top_snp=as.character(sub$snp[ii]),
               top_logp=max(sub$logp, na.rm=TRUE),
               cluster_size_kb=round((en-st)/1000,2))
      }
    }
  )
  
  cl_sum %>% dplyr::arrange(chr, start, end)
}

# -----------------------------
# Helpers clustering (hit-density windows)
# -----------------------------
merge_significant_windows <- function(win_df, gap_bp = 0L) {
  gap_bp <- as.integer(gap_bp)
  if (is.na(gap_bp)) gap_bp <- 0L
  
  dt <- as.data.table(win_df)
  stopifnot(all(c("chr","start","end") %in% names(dt)))
  
  dt[, chr   := as.integer(chr)]
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt <- dt[!is.na(chr) & !is.na(start) & !is.na(end)]
  
  if (nrow(dt) == 0) return(as.data.frame(dt))
  
  setorder(dt, chr, start, end)
  
  merged <- dt[, {
    if (.N == 1L) {
      data.table(start = start[1], end = end[1])
    } else {
      out_s <- integer()
      out_e <- integer()
      
      cur_s <- start[1]
      cur_e <- end[1]
      
      for (i in 2L:.N) {
        s <- start[i]
        e <- end[i]
        if (is.na(s) || is.na(e) || is.na(cur_e)) next
        
        if (s <= (cur_e + gap_bp + 1L)) {
          cur_e <- max(cur_e, e, na.rm = TRUE)
        } else {
          out_s <- c(out_s, cur_s)
          out_e <- c(out_e, cur_e)
          cur_s <- s
          cur_e <- e
        }
      }
      out_s <- c(out_s, cur_s)
      out_e <- c(out_e, cur_e)
      
      data.table(start = out_s, end = out_e)
    }
  }, by = chr]
  
  merged[, center := as.integer(round((start + end) / 2))]
  merged[, cluster_size_kb := round((end - start) / 1000, 2)]
  as.data.frame(merged)
}

make_windows_dt <- function(chr_lengths, win_bp, step_bp) {
  win_bp  <- as.integer(win_bp)
  step_bp <- as.integer(step_bp)
  stopifnot(is.finite(win_bp), win_bp > 0L, is.finite(step_bp), step_bp > 0L)
  
  cl <- as.data.table(chr_lengths)
  stopifnot(all(c("chr","chr_len") %in% names(cl)))
  cl[, chr := as.integer(chr)]
  cl[, chr_len := as.integer(chr_len)]
  cl <- cl[!is.na(chr) & !is.na(chr_len) & chr_len > 0]
  
  wins <- cl[, {
    starts <- seq.int(1L, chr_len, by = step_bp)
    ends   <- pmin.int(starts + win_bp - 1L, chr_len)
    data.table(start = starts, end = ends)
  }, by = chr]
  
  wins[]
}

count_hits_in_windows <- function(hits, windows) {
  h <- as.data.table(hits)
  w <- as.data.table(windows)
  
  stopifnot(all(c("chr","pos") %in% names(h)))
  stopifnot(all(c("chr","start","end") %in% names(w)))
  
  h[, chr := as.integer(chr)]
  h[, pos := as.integer(pos)]
  h <- h[!is.na(chr) & !is.na(pos)]
  
  w[, chr := as.integer(chr)]
  w[, start := as.integer(start)]
  w[, end := as.integer(end)]
  w <- w[!is.na(chr) & !is.na(start) & !is.na(end)]
  setorder(w, chr, start, end)
  
  h[, `:=`(start = pos, end = pos)]
  
  setkey(w, chr, start, end)
  setkey(h, chr, start, end)
  
  ov <- foverlaps(h, w, nomatch = 0L)
  if (nrow(ov) == 0) {
    w[, n_hits := 0L]
    return(w)
  }
  
  cnt <- ov[, .(n_hits = .N), by = .(chr, start, end)]
  w <- merge(w, cnt, by = c("chr","start","end"), all.x = TRUE)
  w[is.na(n_hits), n_hits := 0L]
  w[]
}

add_cluster_stats_from_hits <- function(clusters_df, df_sig) {
  cl <- as.data.table(clusters_df)
  stopifnot(all(c("chr","start","end","cluster_id") %in% names(cl)))
  
  hits <- as.data.table(df_sig)
  stopifnot(all(c("CHR","BP","snp","logp") %in% names(hits)))
  
  cl[, `:=`(
    chr = as.integer(chr),
    start = as.integer(start),
    end = as.integer(end),
    cluster_id = as.character(cluster_id)
  )]
  hits[, `:=`(
    chr  = as.integer(CHR),
    pos  = as.integer(BP),
    snp  = as.character(snp),
    logp = suppressWarnings(as.numeric(logp))
  )]
  
  cl <- cl[is.finite(chr) & is.finite(start) & is.finite(end)]
  hits <- hits[is.finite(chr) & is.finite(pos)]
  
  cl_iv <- cl[, .(chr, start, end, cluster_id)]
  hits_iv <- hits[, .(chr, start = pos, end = pos, pos, snp, logp)]
  
  setkey(cl_iv, chr, start, end)
  setkey(hits_iv, chr, start, end)
  
  ov <- foverlaps(hits_iv, cl_iv, type = "within", nomatch = 0L)
  if ("i.cluster_id" %in% names(ov)) setnames(ov, "i.cluster_id", "cluster_id")
  
  cl[, `:=`(
    n_snps = 0L,
    center = as.integer(round((start + end) / 2)),
    top_snp = NA_character_,
    top_logp = NA_real_
  )]
  
  if (nrow(ov) > 0) {
    stats <- ov[, .(
      n_snps = .N,
      center = as.integer(round(mean(pos, na.rm = TRUE))),
      top_logp = {
        x <- logp
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
      },
      top_snp = {
        x <- logp
        if (all(is.na(x))) NA_character_ else snp[which.max(x)]
      }
    ), by = .(cluster_id)]
    
    cl[stats, on = .(cluster_id), `:=`(
      n_snps  = i.n_snps,
      center  = i.center,
      top_snp = i.top_snp,
      top_logp= i.top_logp
    )]
  }
  
  cl[, cluster_size_kb := round((end - start) / 1000, 2)]
  as.data.frame(cl)
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
##############################################

# hg38 chromosome lengths for cumulative positions
.chr_lengths_hg38 <- function() {
  lens <- c(
    `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555, `5`=181538259,
    `6`=170805979, `7`=159345973, `8`=145138636, `9`=138394717, `10`=133797422,
    `11`=135086622, `12`=133275309, `13`=114364328, `14`=107043718, `15`=101991189,
    `16`=90338345, `17`=83257441, `18`=80373285, `19`=58617616, `20`=64444167,
    `21`=46709983, `22`=50818468, `X`=156040895, `Y`=57227415, `MT`=16569
  )
  ord <- c(as.character(1:22), "X", "Y", "MT")
  df <- data.frame(chr = factor(names(lens), levels = ord), len = as.numeric(lens))
  df$chr_cum <- cumsum(df$len) - df$len
  df$center  <- df$chr_cum + df$len/2
  df
}
.ref_hg38 <- .chr_lengths_hg38()

# Robust (scalar) conversion from cumulative BP → (chr, pos) with clamping
coord_from_bp_cum <- function(bp_cum) {
  ref <- .ref_hg38
  bp_cum <- as.numeric(bp_cum)[1]
  
  idx <- max(which(ref$chr_cum <= bp_cum), na.rm = TRUE)
  if (!is.finite(idx) || length(idx) == 0) idx <- 1
  
  chr <- as.character(ref$chr[idx])
  pos <- round(bp_cum - ref$chr_cum[idx])
  
  # clamp to [1, chr_length]
  pos <- max(pos, 1)
  pos <- min(pos, ref$len[idx])
  
  list(chr = chr, pos = pos)
}

detect_rsid_col <- function(df) {
  cand <- c("rsid","RSID","snp","SNP","marker","ID","variant","rs")
  cand <- intersect(cand, names(df))
  if (length(cand) == 0) return(NULL)
  
  for (c in cand) {
    vals <- as.character(na.omit(df[[c]]))
    if (length(vals) == 0) next
    if (any(grepl("^rs[0-9]+$", vals))) return(c)
  }
  return(cand[1])
}

safe_label_rsid <- function(df) {
  rscol <- NULL
  if ("rsid" %in% names(df)) rscol <- "rsid"
  if (is.null(rscol)) rscol <- detect_rsid_col(df)
  
  if (!is.null(rscol)) {
    out <- as.character(df[[rscol]])
    bad <- is.na(out) | out == ""
    out[bad] <- paste0("rs_unknown_", seq_len(sum(bad)))
    return(out)
  }
  
  paste0("rs_unknown_", seq_len(nrow(df)))
}

# ---------- External tools helpers ----------
has_exec <- function(bin) {
  nzchar(Sys.which(bin))
}

run_cmd <- function(cmd) {
  # run system command, capture output
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) e$message)
  attr(out, "cmd") <- cmd
  out
}

nonsyn_metric_selector_ui <- function() {
  tagList(
    fluidRow(
      column(
        4,
        radioButtons(
          "nonsyn_scale",
          "Metric type",
          c("Score" = "score", "Rankscore" = "rankscore"),
          inline = TRUE
        )
      ),
      column(
        6,
        selectInput(
          "nonsyn_metric",
          "NonSyn metric",
          choices = NULL
        )
      ),
      column(
        2,
        actionButton(
          "info_10",
          label = "Metric info",
          icon  = icon("info-circle"),
          class = "btn btn-outline-secondary"
        )
      )
    ),
    tags$hr(style = "margin:4px 0;"),
    uiOutput("nonsyn_metric_warning")
  )
}


metric_has_data <- function(df, metric) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(FALSE)
  if (!metric %in% names(df)) return(FALSE)
  
  x <- df[[metric]]
  
  # if numeric
  if (is.numeric(x)) {
    return(any(is.finite(x), na.rm = TRUE))
  }
  
  # if character (as in dbNSFP: ".", "NA", etc.)
  x <- trimws(as.character(x))
  x[x %in% c(".", "", "NA", "NaN", "nan", "NULL")] <- NA_character_
  if (all(is.na(x))) return(FALSE)
  
  # european nomenclature: 1.234,56 -> 1234.56
  has_dot  <- grepl("\\.", x)
  has_coma <- grepl(",", x, fixed = TRUE)
  if (any(has_dot & has_coma, na.rm = TRUE)) {
    x2 <- gsub("\\.", "", x)                 # remove grouping dots
    x2 <- gsub(",", ".", x2, fixed = TRUE)   # decimal comma -> dot
    vals2 <- suppressWarnings(readr::parse_number(x2))
    if (any(is.finite(vals2), na.rm = TRUE)) return(TRUE)
  }
  
  # Parse robusto (punto decimal)
  vals <- suppressWarnings(
    readr::parse_number(x, locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
  )
  
  # Fallback: coma decimal simple
  if (!any(is.finite(vals), na.rm = TRUE) && any(grepl(",", x, fixed = TRUE), na.rm = TRUE)) {
    vals <- suppressWarnings(readr::parse_number(gsub(",", ".", x, fixed = TRUE)))
  }
  
  any(is.finite(vals), na.rm = TRUE)
}

make_mode_thr_tag <- function(method, pthr, min_logp) {
  method <- method %||% "window"
  mode_tag <- if (identical(method, "window")) "w" else "h"
  
  thr <- if (identical(method, "window")) pthr else min_logp
  thr <- suppressWarnings(as.numeric(thr))
  
  thr_txt <- if (is.na(thr)) "NA" else format(thr, scientific = FALSE, trim = TRUE)
  thr_txt <- gsub("\\s+", "", thr_txt)
  thr_txt <- gsub("\\.", "p", thr_txt)     # 8.5 -> 8p5
  thr_txt <- gsub("[^0-9A-Za-z_\\-]", "", thr_txt)
  
  list(mode_tag = mode_tag, thr_txt = thr_txt, thr = thr)
}

# -----------------------------
# Helpers clustering methods
# -----------------------------

make_windows_dt <- function(chr_lengths, win_bp, step_bp) {
  library(data.table)
  win_bp  <- as.integer(win_bp)
  step_bp <- as.integer(step_bp)
  stopifnot(is.finite(win_bp), win_bp > 0L, is.finite(step_bp), step_bp > 0L)
  
  cl <- as.data.table(chr_lengths)
  # espera: chr, chr_len
  stopifnot(all(c("chr","chr_len") %in% names(cl)))
  cl[, chr := as.integer(chr)]
  cl[, chr_len := as.integer(chr_len)]
  cl <- cl[!is.na(chr) & !is.na(chr_len) & chr_len > 0]
  
  wins <- cl[, {
    starts <- seq.int(1L, chr_len, by = step_bp)
    ends   <- pmin.int(starts + win_bp - 1L, chr_len)
    data.table(start = starts, end = ends)
  }, by = chr]
  
  wins[]
}

count_hits_in_windows <- function(hits, windows) {
  library(data.table)
  
  h <- as.data.table(hits)
  w <- as.data.table(windows)
  
  # hits: chr,pos
  stopifnot(all(c("chr","pos") %in% names(h)))
  stopifnot(all(c("chr","start","end") %in% names(w)))
  
  h[, chr := as.integer(chr)]
  h[, pos := as.integer(pos)]
  h <- h[!is.na(chr) & !is.na(pos)]
  
  w[, chr := as.integer(chr)]
  w[, start := as.integer(start)]
  w[, end := as.integer(end)]
  w <- w[!is.na(chr) & !is.na(start) & !is.na(end)]
  setorder(w, chr, start, end)
  
  # foverlaps requiere intervalos; para puntos usamos start=end=pos
  h[, `:=`(start = pos, end = pos)]
  
  setkey(w, chr, start, end)
  setkey(h, chr, start, end)
  
  ov <- foverlaps(h, w, nomatch = 0L)  # devuelve columnas i.* (hits) + start/end (ventana)
  if (nrow(ov) == 0) {
    w[, n_hits := 0L]
    return(w)
  }
  
  # cuenta hits por ventana
  cnt <- ov[, .(n_hits = .N), by = .(chr, start, end)]
  w <- merge(w, cnt, by = c("chr","start","end"), all.x = TRUE)
  w[is.na(n_hits), n_hits := 0L]
  w[]
}

# ------------------------------------------------------------
# Helper: 
# ------------------------------------------------------------
read_first_data_line_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    x <- readLines(con, n = 3, warn = FALSE)
  } else {
    x <- readLines(path, n = 3, warn = FALSE)
  }
  x <- x[!grepl("^\\s*$", x)]
  if (length(x) < 2) return(NA_character_)
  x[2]
}

detect_chr_prefix_from_file <- function(path) {
  # mai petis
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(FALSE)
  
  first_data <- tryCatch(read_first_data_line_any(path), error = function(e) NA_character_)
  if (is.na(first_data) || !nzchar(first_data)) return(FALSE)
  
  first_chr <- strsplit(first_data, "\t", fixed = TRUE)[[1]][1]
  isTRUE(grepl("^chr", tolower(first_chr)))
}


# -------------------------
# Helpers (normalization / harmonization / saving)
# -------------------------

sanitize_colnames <- function(df) {
  n <- names(df)
  n <- sub("^\ufeff", "", n)  # remove BOM if present
  n <- trimws(n)              # remove invisible spaces
  names(df) <- n
  df
}

harmonize_dbnsfp_loaded <- function(df) {
  # Identify metric columns (scores / rankscores / phred / raw)
  metric_cols <- grep("(_score$|rankscore$|_phred$|_raw$)", names(df), value = TRUE, ignore.case = TRUE)
  
  # Robust numeric coercion
  for (cc in metric_cols) {
    x <- df[[cc]]
    
    if (is.factor(x))  x <- as.character(x)
    if (is.logical(x)) x <- ifelse(is.na(x), NA_real_, as.numeric(x))
    
    if (is.character(x)) {
      # Use locale-aware parsing already defined in your helpers
      x_num <- suppressWarnings(parse_num_auto(x))
      df[[cc]] <- as.numeric(x_num)
    } else if (is.numeric(x)) {
      df[[cc]] <- x
    } else {
      df[[cc]] <- suppressWarnings(as.numeric(x))
    }
  }
  
  df
}

save_dbnsfp_final <- function(df, out_dir, stem = "dbnsfp_normalized") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(out_dir, paste0(stem, ".csv"))
  out_rds <- file.path(out_dir, paste0(stem, ".rds"))
  
  # Keep original column names (no automatic make.names)
  utils::write.csv(df, out_csv, row.names = FALSE)
  saveRDS(df, out_rds)
  
  list(csv = out_csv, rds = out_rds)
}

# ============================================================
# Robust numeric parsing (auto locale based on comma/dot usage)
# ============================================================

parse_num_auto <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c(".", "", "NA", "NaN", "nan", "NULL")] <- NA
  
  # If commas exist and dots do not, assume comma decimal (EU style)
  has_comma <- sum(grepl(",", x, fixed = TRUE), na.rm = TRUE)
  has_dot   <- sum(grepl("\\.", x), na.rm = TRUE)
  
  if (has_comma > 0 && has_dot == 0) {
    return(readr::parse_number(x, locale = readr::locale(decimal_mark = ",", grouping_mark = ".")))
  }
  
  # Default: dot decimal (comma as grouping)
  readr::parse_number(x, locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
}


# -------------
# Sources (paths robustos)
# -------------
app_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

source(file.path(app_dir, "R", "format_dbnsfp.R"), local = TRUE)
message("format_dbnsfp.R loaded? normalize_dbnsfp exists = ", exists("normalize_dbnsfp", mode = "function"))

# Evita matar l'app sencera si alguna cosa falla en deploy (deixa warning)
if (!exists("normalize_dbnsfp", mode = "function")) {
  warning("normalize_dbnsfp no està disponible (format_dbnsfp.R no s'ha carregat bé).")
}

source(file.path(app_dir, "R", "goslim_utils.R"))

goslim_path <- file.path(app_dir, "R", "goslim_generic.obo")
goslim_ids  <- load_goslim_generic_ids(goslim_path)
message("GO slim generic terms loaded: ", length(goslim_ids))

source(file.path(app_dir, "R", "go_slim_generic.R")) # donde creas GO_SLIM_GENERIC
source(file.path(app_dir, "R", "go_map.R"))

source(file.path(app_dir, "R", "clinical_metrics.R"))

cfg <- gi_cfg()
gi_shared_root <- cfg$shared

# --- LD module (portable via _shared) ---
ld_file <- file.path(gi_shared_root, "mod_ld.R")
stopifnot(file.exists(ld_file))
source(ld_file, local = TRUE)
stopifnot(exists("ld_module_ui"), exists("ld_module_server"))


st_file <- file.path(gi_shared_root, "gi_state.R")
stopifnot(file.exists(st_file))
source(st_file, local = TRUE)

# Helper: pick first existing column name
pick_first_col <- function(df, candidates) {
  candidates <- as.character(candidates)
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  NULL
}

##############################################################################
slave_master_settings_text <- function(state_obj) {
  if (is.null(state_obj) || !is.list(state_obj)) {
    return("Master settings: not available")
  }
  
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0) return("<empty>")
    x <- unlist(x, use.names = FALSE)
    x <- as.character(x)
    x <- x[!is.na(x)]
    if (!length(x)) return("<empty>")
    paste(x, collapse = " | ")
  }
  
  labels <- c(
    sid = "Session ID",
    updated_at = "Updated at",
    cluster_method = "Cluster method",
    hits_mode = "Hits mode",
    thr_value = "Threshold value",
    flank_bp = "Flank (bp)",
    min_hits = "Minimum hits",
    win_bp = "Window size (bp)",
    step_bp = "Step size (bp)",
    has_gwas = "GWAS loaded",
    has_clusters = "Clusters available"
  )
  
  present <- names(labels)[names(labels) %in% names(state_obj)]
  
  if (!length(present)) {
    return("Master settings: no known fields found")
  }
  
  paste(
    vapply(
      present,
      function(nm) paste0(labels[[nm]], ": ", fmt(state_obj[[nm]])),
      character(1)
    ),
    collapse = "\n"
  )
}
##############################################################################
############################# UI #############################################
##############################################################################
ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("🧬 NonSyn Inspector")
  ),
  id = "topnav",
  
  tags$head(
    tags$script(HTML("
(function(){
  function pushQS(){
    try {
      if (window.Shiny && Shiny.setInputValue) {
        Shiny.setInputValue('gi_qs', window.location.search || '', {priority:'event'});
        Shiny.setInputValue('gi_href', window.location.href || '', {priority:'event'});
        return true;
      }
    } catch(e) {}
    return false;
  }

  var tries = 0;
  var iv = setInterval(function(){
    tries++;
    if (pushQS() || tries > 50) clearInterval(iv);
  }, 100);

  window.addEventListener('pageshow', function(){ pushQS(); });
  document.addEventListener('visibilitychange', function(){
    if (!document.hidden) pushQS();
  });
})();
    ")),
    tags$style(HTML("
    /* ----------------------------------------------------
       Global light-gray background for plots + tables
       ---------------------------------------------------- */

    .plotly.html-widget {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    .shiny-plot-output {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    .dataTables_wrapper {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    table.dataTable,
    table.dataTable thead th,
    table.dataTable tbody td {
      background-color: transparent !important;
    }

    .shiny-spinner-output-container {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }
    "))
  ),
  
  # ========================================================================
  # TOP TAB 1: Analysis
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    
    sidebarLayout(
      # --------------------------------------------------------------------
      # LEFT SIDEBAR
      # --------------------------------------------------------------------
      sidebarPanel(
        width = 3,
        
        conditionalPanel(
          condition = "output.hub_mode == 'true'",
          div(
            style = "background:#fff3cd; border:1px solid #ffeeba; padding:10px; border-radius:10px; margin-bottom:10px;",
            HTML("<b>Hub mode</b><br>
         Inputs and clustering synchronized from the <b>Catalog Inspector</b>.<br>
         To change parameters (filtering, clustering, thresholds…), go back to the <b>Catalog Inspector</b> and re-sync.")
          ),
          div(
            class = "panelGrey",
            tags$b("Master session settings"),
            verbatimTextOutput("master_session_status")
          )
        ),
        
        # Hide ONLY Step 1 + Step 2 when everything comes from Hub
        conditionalPanel(
          condition = "output.hub_mode == 'false'",
          
          # -----------------------------
          # Step 1 · Load GWAS hits
          # -----------------------------
          h3("Step 1 · Load GWAS hits",
             style = "color:#1A4E8A; font-size:22px; font-weight:700;"),
          
          fluidRow(
            column(6, actionButton("info_00", "ℹ️ file format")),
            column(6, actionButton(
              "reset_case", "Reset",
              icon = icon("rotate-left"),
              class = "btn-warning"
            ))
          ),
          
          fileInput(
            "gwas_file",
            "GWAS p-value table (TSV/CSV)",
            accept = c(".tsv", ".txt", ".csv")
          ),
          
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
          # Step 2 · Clustering
          # -----------------------------
          h3("Step 2 · Clustering GWAS hits",
             style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
          actionButton("info_01", "ℹ️ Clustering method", class = "btn btn-default"),
          
          radioButtons(
            "cluster_method",
            "Clustering method:",
            choices  = c(
              "By hit intervals (thr + flank → merge)" = "window",
              "By hit density (min_logp + min_hits)"   = "hits"
            ),
            selected = "window"
          ),
          
          conditionalPanel(
            condition = "input.cluster_method == 'window'",
            sliderInput("pthr", "-log10(P) threshold", min = 2, max = 20, value = 8, step = 0.5),
            numericInput("flank", "Flank (+/- bp)", value = 10000, min = 0, max = 10000000, step = 1000)
          ),
          
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
              "min_hits",
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
              numericInput("step_bp", "Step (bp)", value = 1e5, min = 1e3, max = 5e7, step = 1e3)
            ),
            
            conditionalPanel(
              condition = "input.hits_mode != 'sliding'",
              helpText("For 'tiled' and '1Mb hit-span', step is implicit.")
            )
          ),
          
          actionButton(
            "build_ranges",
            "Generate intervals & clusters",
            style = "background-color:#ffdd57; font-weight:bold;"
          ),
          
          h5("Preview selected intervals range"),
          verbatimTextOutput("ranges_preview"),
          tags$hr()
        ),
        
        # ------------------------------------------------------------
        # ALWAYS VISIBLE: Step 3 + Step 4 (also in Hub mode)
        # ------------------------------------------------------------
        h3("Step 3 · Extract nonsyn variants in all clusters",
           style="color:#1A4E8A; font-weight:700;"),
        
        textInput(
          "dbnsfp_out_all",
          "Path to dbNSFP biallelic_nonsyn.out file",
          value = ns_def$dbnsfp_all
        ),
        
        actionButton(
          "run_nonsyn_clusters",
          "Extract nonsyn variants in all clusters",
          style = "background-color:#ffdd57; font-weight:bold;"
        ),
        
        h5("Extract log file"),
        verbatimTextOutput("dbnsfp_log"),
        tags$hr(),
        
        h3("Step 4 · Download files",
           style="color:#1A4E8A; font-weight:700;"),
        downloadButton("dl_dbnsfp_csv", "Download normalized CSV"),
        downloadButton("dl_dbnsfp_rds", "Download normalized RDS"),
        tags$hr(),
        downloadButton("dl_candidates_zip", "⬇️ Download NonSyn candidates (ZIP)")
      ),
      
      # --------------------------------------------------------------------
      # MAIN PANEL (Analysis outputs)
      # --------------------------------------------------------------------
      mainPanel(
        
        conditionalPanel(
          condition = "input.main_tabs == 'Manhattan' || input.main_tabs == 'Summary stats' || input.main_tabs == 'LD'",
          wellPanel(
            style = "margin-bottom:6px; padding:6px 10px;",
            nonsyn_metric_selector_ui()
          )
        ),
        
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            value = "Manhattan",
            
            fluidRow(
              column(
                12, 
                div(class = "panel-lite", shinycssloaders::withSpinner(plotlyOutput("manhattan_combo", height = "700px")))
              )
            ),
            
            tags$hr(),
            
            fluidRow(
              column(6, uiOutput("debug_ucsc_state")),
              column(2, h4("UCSC links to expanded region:")),
              column(
                4,
                div(
                  style = "margin-top:8px; background-color:#f2f2f2; padding:10px; border-radius:8px; width:100%;",
                  uiOutput("ucsc_link_gwas"),
                  uiOutput("ucsc_link_nonsyn")
                )
              )
            ),
            
            fluidRow(
              column(
                12,
                tags$hr(),
                h4("Clusters summary"),
                DTOutput("cluster_dt")
              )
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📄 NonSyn – dbNSFP</span>"),
            value = "NonSyn",
            uiOutput("hm_chr_selector2"),
            
            radioButtons(
              "clinical_class",
              "Functional class",
              c("Default","Pathogenicity","Conservation","Functional impact","LoF/Haploinsufficiency","Gene damage"),
              inline = TRUE
            ),
            
            radioButtons(
              "metric_kind",
              "Metric type",
              choices = c("All" = "all", "Score" = "score", "Rankscore" = "rankscore"),
              selected = "all",
              inline = TRUE
            ),
            
            DTOutput("nonsyn_table")
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔬 GWAS hits</span>"),
            value = "GWAS",
            uiOutput("gwas_tab_body")
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Summary stats</span>"),
            value = "Summary stats",
            
            fluidRow(
              column(
                12,
                shinycssloaders::withSpinner(plotlyOutput("nonsyn_box_combo_plotly", height = "350px")),
                fluidRow(
                  column(4, downloadButton("dl_nonsyn_chr_stats_csv",     "⬇️ Summary stats by chr CSV")),
                  column(4, downloadButton("dl_nonsyn_gene_stats_csv",    "⬇️ Summary stats by gene CSV")),
                  column(4, downloadButton("dl_nonsyn_cluster_stats_csv", "⬇️ Summary stats by cluster CSV"))
                )
              )
            ),
            
            tags$hr(),
            
            fluidRow(
              column(6, h4("Top genes barplot"), shinycssloaders::withSpinner(plotOutput("nonsyn_gene_bar", height = "300px"))),
              column(6, h4("Metric distribution"), shinycssloaders::withSpinner(plotOutput("nonsyn_metric_bar", height = "300px")))
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🕸️ Radar plots</span>"),
            value = "tab_radar",
            
            sidebarLayout(
              sidebarPanel(
                width = 2,
                
                radioButtons(
                  "radar_mode",
                  "Group by",
                  choices = c(
                    "Global"  = "global",
                    "Cluster" = "cluster",
                    "Gene"    = "gene"
                  ),
                  selected = "global"
                ),
                
                conditionalPanel(
                  condition = "input.radar_mode == 'cluster'",
                  selectInput("radar_chr_for_cluster", "Chromosome", choices = NULL),
                  selectInput("radar_cluster_id", "Cluster", choices = NULL)
                ),
                
                conditionalPanel(
                  condition = "input.radar_mode == 'gene'",
                  selectInput("radar_gene", "Gene (by cluster)", choices = NULL)
                )
              ),
              
              mainPanel(
                width = 10,
                actionButton("info_11", label = "Metric info", icon = icon("info-circle"), class = "btn btn-outline-secondary"),
                shinycssloaders::withSpinner(plotlyOutput("nonsyn_radar_plotly", height = "900px")),
                tags$hr(),
                DT::DTOutput("radar_selection_table")
              )
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧾 Function & Disease</span>"),
            value = "FuncDisease",
            DT::DTOutput("func_disease_tbl")
          )
        )
      )
    )
  ),
  
  # ========================================================================
  # TOP TAB 2: Enrichment (Predictors / Terms / GO / KEGG / GoSlim)
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>✳️ Enrichment</span>"),
    value = "TopEnrich",
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        radioButtons(
          "func_scope",
          "Scope",
          choices  = c("Global" = "global", "Cluster" = "cluster", "Gene" = "gene"),
          selected = "global"
        ),
        
        conditionalPanel(
          condition = "input.func_scope == 'cluster'",
          uiOutput("func_cluster_ui")
        ),
        
        conditionalPanel(
          condition = "input.func_scope == 'gene'",
          uiOutput("func_gene_ui")
        ),
        
        conditionalPanel(
          condition = "input.enrich_tabs != 'tab_enrich_pred' && input.enrich_tabs != 'tab_enrich_terms'",
          selectInput(
            "enrich_background", "Background",
            choices = c("Reference annotated genes" = "orgdb", "Dataset genes" = "dataset"),
            selected = "orgdb"
          ),
          numericInput("enrich_pcut", "FDR cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          tags$hr()
        ),
        
        conditionalPanel(
          condition = "input.enrich_tabs == 'tab_enrich_pred'",
          
          h4("Predictor enrichment"),
          checkboxInput("pred_split_pipe", "Split '|' (Aloft) into single labels", TRUE),
          numericInput("pred_topk", "Top items", value = 25, min = 5, max = 300, step = 5),
          numericInput("pred_min_a", "Min occurrences in foreground (min_a)", value = 2, min = 1, max = 1000, step = 1),
          tags$hr()
        ),
        
        conditionalPanel(
          condition = "input.enrich_tabs == 'tab_enrich_terms'",
          
          h4("NonSyn term enrichment"),
          uiOutput("nonsyn_term_field_ui"),
          helpText("Runs ORA on the selected dbNSFP field (unit = rows; multi-terms are split)."),
          tags$hr(),
          
          numericInput("nonsyn_terms_topk", "Preselect top terms (by foreground count)", value = 200, min = 20, max = 5000, step = 20),
          numericInput("nonsyn_terms_min_a", "Min occurrences in foreground (min_a)", value = 2, min = 1, max = 1000, step = 1),
          tags$hr()
        ),
        
        conditionalPanel(
          condition = "input.enrich_tabs == 'tab_enrich_go' || input.enrich_tabs == 'tab_enrich_goslim'",
          
          checkboxGroupInput(
            "go_ontos", "GO ontologies",
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
        
        conditionalPanel(
          condition = "input.enrich_tabs == 'tab_enrich_go' || input.enrich_tabs == 'tab_enrich_kegg' || input.enrich_tabs == 'tab_enrich_goslim'",
          actionButton("info_12", "ℹ GSSize"),
          numericInput("enrich_min_gs", "minGSSize", value = 10, min = 1, step = 1),
          numericInput("enrich_max_gs", "maxGSSize", value = 500, min = 10, step = 10),
          tags$hr()
        ),
        
        actionButton(
          "run_enrich",
          label = "Run enrichment",
          icon  = icon("play"),
          class = "btn btn-primary"
        ),
        tags$hr()
      ),
      
      mainPanel(
        width = 9,
        
        uiOutput("enrich_bg_note"),
        
        tabsetPanel(
          id = "enrich_tabs",
          selected = "tab_enrich_pred",
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>🧠 Predictors</span>"),
            value = "tab_enrich_pred",
            shinycssloaders::withSpinner(plotly::plotlyOutput("enrich_nonsyn_pred_plotly_radial", height = "520px")),
            
            tags$hr(),
            shinycssloaders::withSpinner(DT::DTOutput("enrich_nonsyn_pred_table")),
            tags$hr(),
            downloadButton("pred_download", "Download table")
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>🧾 Terms</span>"),
            value = "tab_enrich_terms",
            shinycssloaders::withSpinner(plotly::plotlyOutput("enrich_nonsyn_terms_plotly", height = "480px")),
            
            tags$hr(),
            shinycssloaders::withSpinner(DT::DTOutput("enrich_nonsyn_terms_table"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>GO</span>"),
            value = "tab_enrich_go",
            shinycssloaders::withSpinner(plotlyOutput("go_bar", height = 400)),
            tags$hr(),
            shinycssloaders::withSpinner(DT::DTOutput("go_table"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>KEGG</span>"),
            value = "tab_enrich_kegg",
            shinycssloaders::withSpinner(plotlyOutput("kegg_bar", height = 400)),
            tags$hr(),
            shinycssloaders::withSpinner(DT::DTOutput("kegg_table"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>GoSlim</span>"),
            value = "tab_enrich_goslim",
            shinycssloaders::withSpinner(plotlyOutput("goslim_bar", height = 450)),
            tags$hr(),
            shinycssloaders::withSpinner(DT::DTOutput("goslim_table"))
          )
        )
      )
    )
  ),
  
  # ========================================================================
  # TOP TAB 3: LD
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 LD</span>"),
    ld_module_ui("ld", app_tag = "nonsyn")
  )
)

##############################################################################
############################### SERVER #######################################
##############################################################################
server <- function(input, output, session) {
  
  source(file.path(SHARED, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
  
  # ===========================
  # 1) SLAVE sync (Hub -> app)
  # ===========================
  gi_sync <- gi_slave_canonical_init(session)
  
  # --- FLAGS per evitar flicker / resets mentre arriba sync del HUB ---
  rv_hub <- reactiveValues(
    hub_override = FALSE,
    applying     = FALSE
  )
  
  # Guardem GWAS rebut (si no el tens ja)
  if (!exists("gwas_shared_rv", inherits = TRUE)) {
    gwas_shared_rv <- reactiveVal(tibble::tibble())
  }
  
  observeEvent(gi_sync$gwas_shared(), {
    id <- paste0("hub_sync_", as.integer(Sys.time()))
    showNotification("🔄 HUB: sincronitzant GWAS…", type = "message", duration = NULL, id = id)
    
    df <- gi_sync$gwas_shared()
    req(is.data.frame(df), nrow(df) > 0)
    
    gwas_shared_rv(df)
    cat("[SLAVE] applied GWAS rows=", nrow(df), "\n")
    
    removeNotification(id)
    showNotification(sprintf("✅ HUB: GWAS sincronitzat (%d files).", nrow(df)),
                     type = "message", duration = 2)
  }, ignoreInit = FALSE)
  
  ################################ DETECT HUB MODE ##############################
  
  `%||%` <- function(a,b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b
  
  hub_mode_r <- reactive({
    qs <- session$clientData$url_search %||% ""
    q  <- tryCatch(shiny::parseQueryString(sub("^\\?", "", qs)), error = function(e) list())
    
    sid_from_url <- q$sid %||% q$gi_sid %||% q$SID %||% NULL
    has_sid_in_url <- !is.null(sid_from_url) && nzchar(as.character(sid_from_url))
    
    sid_ok <- FALSE
    if (exists("gi_sync") && is.list(gi_sync) && "sid" %in% names(gi_sync)) {
      s <- tryCatch(gi_sync$sid(), error = function(e) NULL)
      sid_ok <- !is.null(s) && nzchar(as.character(s))
    }
    
    hub_flag <- FALSE
    if (exists("gi_sync") && is.list(gi_sync) && "state_shared" %in% names(gi_sync)) {
      st <- tryCatch(gi_sync$state_shared(), error = function(e) NULL)
      hub_flag <- is.list(st) && isTRUE(st$from_hub %||% st$hub %||% st$sync %||% FALSE)
    }
    
    has_shared <- FALSE
    if (exists("gi_sync") && is.list(gi_sync) && "gwas_shared" %in% names(gi_sync)) {
      gw <- tryCatch(gi_sync$gwas_shared(), error = function(e) NULL)
      has_shared <- is.data.frame(gw) && nrow(gw) > 0
    }
    if (!has_shared && exists("gi_sync") && is.list(gi_sync) && "clusters_shared" %in% names(gi_sync)) {
      cl <- tryCatch(gi_sync$clusters_shared(), error = function(e) NULL)
      has_shared <- is.data.frame(cl) && nrow(cl) > 0
    }
    
    isTRUE(has_sid_in_url || sid_ok || hub_flag || has_shared)
  })
  
  master_session_status <- reactive({
    st <- NULL
    
    if (exists("gi_sync") && is.list(gi_sync) && "state_shared" %in% names(gi_sync)) {
      st <- tryCatch(gi_sync$state_shared(), error = function(e) NULL)
    }
    
    slave_master_settings_text(st)
  })
  
  
  output$master_session_status <- renderText({
    master_session_status()
  })
  
  output$hub_mode <- renderText({
    if (isTRUE(hub_mode_r())) "true" else "false"
  })
  outputOptions(output, "hub_mode", suspendWhenHidden = FALSE)
  
  output$dbg_hub_mode <- renderPrint({
    list(
      url_search = session$clientData$url_search %||% NA,
      hub_mode   = hub_mode_r()
    )
  })
  outputOptions(output, "dbg_hub_mode", suspendWhenHidden = FALSE)
  
  ######################################################################################
  ######################################################################################
  #  # ===========================
  #  # 2) CANONICAL cluster engine (local build_ranges)
  #  # ===========================
  #  gwas_df <- reactive({
  #    # ### CHANGE ### si a la teva app NO tens input$use_preloaded_*,
  #    # simplement retorna sempre gwas_shared_rv()
  #    if (!exists("use_preloaded_dbnsfp", where = names(input), inherits = FALSE)) {
  #      return(gwas_shared_rv())
  #    }
  #    if (!identical(input$use_preloaded_dbnsfp, "new")) return(tibble::tibble())
  #    gwas_shared_rv()
  #  })
  
  cl_engine <- gi_clusters_canonical_init(
    session = session, input = input, output = output,
    gwas_df = gwas_df,
    build_btn_id    = "build_ranges",   # ### CHANGE ### si el teu botó té un altre id
    clusters_dt_id  = "cluster_dt",     # ### CHANGE ### si la teva DT té un altre id
    hits_rows_id    = "hits_tbl_rows_selected",  # ### CHANGE ### si no existeix, posa el que toque o NULL si la funció ho accepta
    app_count_col   = "n_APP"           # ### CHANGE ### (veure taula més avall)
  )
  
  # Accessors canònics
  intervals_raw     <- cl_engine$intervals_raw
  clusters_cur      <- cl_engine$clusters_cur
  selected_cluster  <- cl_engine$selected_cluster
  
  
  # ===========================
  # 3) Hub clusters override (si arriben del master, manen)
  # ===========================
  hub_master_clusters <- reactiveVal(NULL)
  hub_master_locked   <- reactiveVal(FALSE)
  
  observeEvent(gi_sync$clusters_shared(), {
    cl <- gi_sync$clusters_shared()
    req(is.data.frame(cl), nrow(cl) > 0)
    
    cl <- standardize_cluster_ids(cl)
    
    if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
    
    cl$chr   <- suppressWarnings(as.integer(gsub("^chr", "", as.character(cl$chr), ignore.case = TRUE)))
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    req(nrow(cl) > 0)
    
    clusters_cur(cl)
    intervals_raw(
      dplyr::as_tibble(cl) %>%
        dplyr::transmute(
          chr = .data$chr, start = .data$start, end = .data$end,
          label = dplyr::coalesce(as.character(.data$cluster_id), as.character(.data$cluster_chr_n))
        )
    )
    
    hub_master_clusters(cl)
    hub_master_locked(TRUE)
    
    cat("[SLAVE] applied CLUSTERS rows=", nrow(cl), "\n")
    showNotification(sprintf("✅ HUB: CLUSTERS sincronitzats (%d).", nrow(cl)),
                     type = "message", duration = 2)
  }, ignoreInit = FALSE)
  
  
  # --- Watchdog: si el motor local trepitja/buida clusters_cur(), els restituïm del master ---
  observe({
    invalidateLater(600, session)
    
    if (!isTRUE(hub_master_locked())) return()
    cl_master <- hub_master_clusters()
    if (!is.data.frame(cl_master) || nrow(cl_master) == 0) return()
    
    cl_now <- tryCatch(clusters_cur(), error = function(e) NULL)
    
    if (!is.data.frame(cl_now) || nrow(cl_now) == 0) {
      clusters_cur(cl_master)
      intervals_raw(
        dplyr::as_tibble(cl_master) %>%
          dplyr::transmute(
            chr = .data$chr, start = .data$start, end = .data$end,
            label = dplyr::coalesce(as.character(.data$cluster_id), as.character(.data$cluster_chr_n))
          )
      )
      cat("[SLAVE] watchdog: restored master CLUSTERS rows=", nrow(cl_master), "\n")
    }
  })
  
  
  # ===========================
  # 4) clusters_src (una sola font per overlays/plots)
  # ===========================
  clusters_src <- reactive({
    cl <- tryCatch(clusters_cur(), error = function(e) NULL)
    if (is.data.frame(cl) && nrow(cl) > 0) return(cl)
    NULL
  })
  
  ####################################
  
  # ---- Global artifacts  ----
  vcf_out_path <- reactiveVal(NULL)
  dbsnp_hits_val <- reactiveVal(NULL)
  
  dbnsfp_final_path_csv <- reactiveVal(NULL)
  dbnsfp_final_path_rds <- reactiveVal(NULL)
  
  dbnsfp_norm_path_csv <- reactiveVal(NULL)
  dbnsfp_norm_path_rds <- reactiveVal(NULL)
  
  # ============================================================
  # Robust numeric parsing (auto locale based on comma/dot usage)
  # ============================================================
  
  parse_num_auto <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c(".", "", "NA", "NaN", "nan", "NULL")] <- NA
    
    # If commas exist and dots do not, assume comma decimal (EU style)
    has_comma <- sum(grepl(",", x, fixed = TRUE), na.rm = TRUE)
    has_dot   <- sum(grepl("\\.", x), na.rm = TRUE)
    
    if (has_comma > 0 && has_dot == 0) {
      return(readr::parse_number(x, locale = readr::locale(decimal_mark = ",", grouping_mark = ".")))
    }
    
    # Default: dot decimal (comma as grouping)
    readr::parse_number(x, locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
  }
  
  # ============================================================
  # Session workdir (temporary folder, removed on session end)
  # ============================================================
  
  workdir <- tempfile("nonsyn_inspector_")
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  
  session$onSessionEnded(function(...) {
    if (dir.exists(workdir)) unlink(workdir, recursive = TRUE, force = TRUE)
  })
  
  # ============================================================
  # Annotation state + logging buffer
  # ============================================================
  
  ann_path  <- reactiveVal(NULL)   # Path to generated .out
  ann_ready <- reactiveVal(FALSE)  # Flag for conditional UI
  
  dbnsfp_empty <- reactiveVal(FALSE) # message if final output is empty
  dbnsfp_msg   <- reactiveVal("")
  
  dbnsfp_log_text <- reactiveVal("")
  
  output$dbnsfp_log <- renderText({
    dbnsfp_log_text()
  })
  
  append_dbnsfp_log <- function(..., sep = "") {
    old <- dbnsfp_log_text()
    dbnsfp_log_text(paste(c(old, paste(..., collapse = sep)), collapse = sep))
  }
  
  # ============================================================
  # Step 1: Read GWAS table
  # ============================================================
  
  # -----------------------------
  # GWAS raw (shared-first)
  # -----------------------------
  dat_raw <- reactive({
    
    # 1) Preferència: GWAS compartit (master)
    df_shared <- gwas_shared_rv()
    if (is.data.frame(df_shared) && nrow(df_shared) > 0) {
      return(df_shared)
    }
    
    # 2) Fallback: fitxer upload
    req(input$gwas_file)
    
    readr::read_delim(
      file = input$gwas_file$datapath,
      delim = input$gwas_sep %||% "\t",
      col_names = isTRUE(input$gwas_header),
      locale = readr::locale(decimal_mark = ".", grouping_mark = ","),
      show_col_types = FALSE,
      na = c("", "NA")
    )
  })
  
  # -----------------------------
  # UI selector P column (shared-aware)
  # -----------------------------
  output$gwas_p_selector <- renderUI({
    
    df <- dat_raw()
    if (is.null(df) || !nrow(df)) return(helpText("GWAS table is empty."))
    
    cols <- names(df)
    
    # Si ve del master sovint ja tens Pval/logp; si ve d'upload, pot ser "P"/"pval"...
    guessed <- (
      grep("^p$|^pval$|p_value|pvalue|^pval_nominal$|^pval$|^Pval$", cols, ignore.case = TRUE, value = TRUE)[1]
      %||% cols[1]
    )
    
    selectInput(
      "gwas_col_p",
      "Select P-value column:",
      choices  = cols,
      selected = guessed
    )
  })
  
  # -----------------------------
  # Standardized GWAS dataframe used across the app (shared-aware)
  # -----------------------------
  gwas_df <- reactive({
    
    df <- dat_raw()
    validate(need(is.data.frame(df) && nrow(df) > 0, "GWAS table is empty."))
    
    # Detecta si ja està estandarditzat (ve del master)
    is_std <- all(c("CHR","BP","Pval") %in% names(df))
    
    if (is_std) {
      
      # Assegura camps mínims
      if (!"snp" %in% names(df)) df$snp <- NA_character_
      if (!"rsid" %in% names(df)) df$rsid <- df$snp
      
      out <- df %>%
        dplyr::mutate(
          CHR  = suppressWarnings(as.integer(CHR)),
          BP   = suppressWarnings(as.numeric(BP)),
          Pval = suppressWarnings(as.numeric(Pval)),
          logp = if ("logp" %in% names(df)) suppressWarnings(as.numeric(logp)) else -log10(Pval)
        ) %>%
        dplyr::filter(!is.na(CHR), !is.na(BP), !is.na(Pval), Pval > 0)
      
      return(out)
    }
    
    # --------- Upload mode (need P column selection) ----------
    req(input$gwas_col_p)
    
    colP <- input$gwas_col_p
    rawP <- df[[colP]]
    
    validate(
      need(!is.null(rawP), "Selected P column not found."),
      need(length(rawP) > 0, "Selected P column is empty.")
    )
    
    # Normalitza noms a lower-case
    names(df) <- tolower(names(df))
    
    validate(
      need(all(c("chr","bp") %in% names(df)), "GWAS table must contain columns: chr, bp")
    )
    
    if (!"snp" %in% names(df)) df$snp <- NA_character_
    
    BP   <- if (is.numeric(df$bp)) df$bp else readr::parse_number(as.character(df$bp))
    Pval <- parse_p_robust(rawP)
    CHR  <- chr_map_plink19(df$chr)
    
    out <- df %>%
      dplyr::mutate(
        CHR  = CHR,
        BP   = as.numeric(BP),
        Pval = as.numeric(Pval),
        logp = -log10(Pval)
      ) %>%
      dplyr::filter(!is.na(CHR), !is.na(BP), !is.na(Pval), Pval > 0)
    
    # rsid si es pot detectar
    rscol <- detect_rsid_col(out)
    if (!is.null(rscol) && rscol %in% names(out)) {
      out$rsid <- as.character(out[[rscol]])
    } else {
      out$rsid <- as.character(out$snp)
    }
    
    out
  })
  
  # ============================================================
  # Manhattan helpers
  # ============================================================
  
  dfp_manhattan <- reactive({
    df <- gwas_df()
    req(is.data.frame(df), nrow(df) > 0)
    # filter pval < 0.0 to a easely plot
    df <- df %>% dplyr::filter(Pval < 0.05)
    
    ref <- .ref_hg38
    
    # Convert numeric CHR to labels then normalize (1..22,X,Y,MT)
    df$chrN <- norm_chr_generic(chr_label_plink(df$CHR))
    
    df %>%
      inner_join(ref %>% dplyr::select(chr, chr_cum), by = c("chrN" = "chr")) %>%
      arrange(chrN, BP) %>%
      mutate(
        BPcum  = BP + chr_cum,
        colgrp = CHR %% 2
      )
  })
  
  axis_df <- reactive({
    dfp <- dfp_manhattan()
    dfp %>%
      group_by(chrN) %>%
      summarise(center = mean(BPcum, na.rm = TRUE), .groups = "drop")
  })
  
  # ============================================================
  # GWAS hits table (only meaningful in NEW mode)
  # ============================================================
  
  hits_df <- reactive({
    df <- gwas_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    cluster_method <- as.character(input$cluster_method %||% "")
    hits_mode      <- as.character(input$hits_mode %||% "")
    
    thr <- if (identical(cluster_method, "window")) {
      suppressWarnings(as.numeric(input$pthr))
    } else if (identical(cluster_method, "hits")) {
      suppressWarnings(as.numeric(input$min_logp))
    } else {
      NA_real_
    }
    
    req(is.finite(thr))
    
    # garantir rsid (si no existeix, usa snp)
    if (!"rsid" %in% names(df)) df$rsid <- NA_character_
    if (!"snp"  %in% names(df)) df$snp  <- NA_character_
    
    df$rsid <- as.character(df$rsid)
    df$snp  <- as.character(df$snp)
    df$rsid <- dplyr::coalesce(df$rsid, df$snp)
    
    df %>%
      dplyr::filter(is.finite(logp), .data$logp >= thr) %>%
      dplyr::arrange(dplyr::desc(.data$logp)) %>%
      dplyr::select(
        dplyr::all_of(c("CHR", "BP", "snp")),
        dplyr::any_of("rsid"),
        p = .data$Pval,
        .data$logp
      )
  })
  
  output$hits_tbl <- renderDT({
    h <- hits_df()
    
    if (is.null(h) || !nrow(h)) {
      return(DT::datatable(
        data.frame(Message = "No hits above threshold."),
        options = list(dom = "t")
      ))
    }
    
    h2 <- h %>%
      mutate(
        p    = formatC(p, format = "e", digits = 2),
        logp = sprintf("%.2f", as.numeric(logp)),
        dbSNP = ifelse(
          !is.na(rsid) & rsid != "" & grepl("^rs[0-9]+$", rsid),
          sprintf("<a href='https://www.ncbi.nlm.nih.gov/snp/%s' target='_blank'>%s</a>", rsid, rsid),
          ""
        )
      )
    
    DT::datatable(
      h2,
      selection = "multiple",
      rownames = FALSE,
      escape = FALSE,
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
  
  output$gwas_tab_body <- renderUI({
    tagList(
      h4("GWAS hits used for clustering"),
      DTOutput("hits_tbl")
    )
  })
  
  # ============================================================
  # STEP 2 · INTERVALS -> CLUSTERS (STABLE / CANONICAL)
  # ============================================================
  
  # ------------------------------------------------------------
  # Cluster view (for DT display + optional nonsyn counts)
  # ------------------------------------------------------------
  cluster_dt_view <- reactive({
    cl <- clusters_cur()   # <-- CANVI clau: ve del motor canonical (cl_engine)
    if (!is.data.frame(cl) || nrow(cl) == 0) return(tibble::tibble())
    
    cl <- as.data.frame(cl)
    cl <- standardize_cluster_ids(cl)
    
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
    if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
    
    out <- cl %>%
      dplyr::mutate(
        chr_num = suppressWarnings(as.integer(gsub("^chr", "", as.character(chr), ignore.case = TRUE))),
        start   = suppressWarnings(as.integer(start)),
        end     = suppressWarnings(as.integer(end))
      ) %>%
      dplyr::filter(is.finite(chr_num), is.finite(start), is.finite(end), end >= start) %>%
      dplyr::mutate(
        chr     = chr_label_plink(chr_num),
        size_kb = (end - start) / 1000
      )
    
    if ("cluster_size_kb" %in% names(out)) out$size_kb <- suppressWarnings(as.numeric(out$cluster_size_kb))
    
    ############
    final_rds <- tryCatch(dbnsfp_final_path_rds(), error = function(e) NULL)
    
    ns_ok <- is.character(final_rds) && length(final_rds) == 1 && nzchar(final_rds) && file.exists(final_rds)
    
    if (ns_ok) {
      ns <- readRDS(final_rds) %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          chr_num    = suppressWarnings(as.integer(readr::parse_number(as.character(chr)))),
          pos        = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
          ref        = as.character(ref),
          alt        = as.character(alt)
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          is.finite(chr_num), is.finite(pos)
        ) %>%
        dplyr::distinct(cluster_id, chr_num, pos, ref, alt)
      
      counts <- ns %>%
        dplyr::count(cluster_id, name = "n_nonsyn")
      
      out <- out %>%
        dplyr::left_join(counts, by = "cluster_id")
      
      out$n_nonsyn[is.na(out$n_nonsyn)] <- 0L
      
    } else {
      out$n_nonsyn <- NA_integer_
    }
    
    out
  })
  
  
  cluster_dt_view2 <- reactive(cluster_dt_view())
  
  # ------------------------------------------------------------
  # Outputs
  # NOTE: UI uses DTOutput("cluster_dt"), so keep this one as primary.
  # ------------------------------------------------------------
  
  output$cluster_dt <- DT::renderDT({
    
    dt <- cluster_dt_view2()
    
    # Avoid console spam in production (commented out)
    # cat("[DEBUG] cluster_dt_view2 =\n"); print(head(dt))
    
    if (is.null(dt) || !nrow(dt)) {
      return(DT::datatable(
        data.frame(Message = "No clusters yet. Click 'Generate intervals & clusters' in Step 2."),
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons = list(
            list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
            list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
            list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
          ),
          pageLength = 10,
          options = list(dom = "t"),
          scrollX    = TRUE
        )
        
      ))
    }
    
    dt2 <- dt %>%
      dplyr::select(cluster_id, chr, start, end, n_snps, n_nonsyn, top_snp, top_logp, size_kb) %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        chr        = as.character(chr),
        start      = as.integer(start),
        end        = as.integer(end),
        n_snps     = as.integer(n_snps),
        n_nonsyn   = as.integer(n_nonsyn),
        top_snp    = as.character(top_snp),
        top_logp   = as.numeric(top_logp),
        size_kb    = as.numeric(size_kb)
      )
    
    DT::datatable(dt2, rownames = FALSE, 
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
    )  %>%
      DT::formatRound(c("top_logp", "size_kb"), digits = 2)
  }, server = FALSE)
  
  
  # ------------------------------------------------------------
  # Helper for Step 3: annotate dbNSFP rows with cluster membership
  # ------------------------------------------------------------
  
  add_clusters_to_dbnsfp <- function(df, cl) {
    
    stopifnot(is.data.frame(df), is.data.frame(cl))
    validate(need("cluster_id" %in% names(cl), "Expected cluster_id in clusters (use cluster_dt_view2())."))
    
    # Prepare dbNSFP (expects columns 'chr' and 'BP' in df)
    df2 <- df %>%
      dplyr::mutate(
        .chr = norm_chr_generic(.data[["chr"]]),
        .pos = suppressWarnings(as.integer(.data[["BP"]]))
      )
    
    # Prepare clusters
    cl2 <- cl %>%
      dplyr::transmute(
        chr        = norm_chr_generic(.data[["chr"]]),
        start      = suppressWarnings(as.integer(.data[["start"]])),
        end        = suppressWarnings(as.integer(.data[["end"]])),
        cluster_id = as.character(.data[["cluster_id"]]),
        # --- FIX: cluster_chr_n sempre existeix ---
        cluster_chr_n = {
          if ("cluster_chr_n" %in% names(cl)) {
            suppressWarnings(as.integer(as.character(cl[["cluster_chr_n"]])))
          } else if ("cluster" %in% names(cl)) {
            suppressWarnings(as.integer(as.character(cl[["cluster"]])))
          } else {
            NA_integer_
          }
        }
      ) %>%
      dplyr::filter(!is.na(chr), is.finite(start), is.finite(end), end >= start, nzchar(cluster_id)) %>%
      dplyr::arrange(chr, start, end)
    
    # --- si cluster_chr_n continua sent NA, el generem per chr (1..n) ---
    if (all(is.na(cl2$cluster_chr_n))) {
      cl2 <- cl2 %>%
        dplyr::group_by(.data$chr) %>%
        dplyr::arrange(.data$start, .data$end, .by_group = TRUE) %>%
        dplyr::mutate(cluster_chr_n = dplyr::row_number()) %>%
        dplyr::ungroup()
    } else {
      # Omple NAs parcials per chr
      cl2 <- cl2 %>%
        dplyr::group_by(.data$chr) %>%
        dplyr::arrange(.data$start, .data$end, .by_group = TRUE) %>%
        dplyr::mutate(cluster_chr_n = dplyr::if_else(is.na(.data$cluster_chr_n),
                                                     dplyr::row_number(),
                                                     .data$cluster_chr_n)) %>%
        dplyr::ungroup()
    }
    
    n <- nrow(df2)
    out_cluster_id    <- rep(NA_character_, n)
    out_cluster_chr_n <- rep(NA_integer_,   n)
    out_cluster_start <- rep(NA_integer_,   n)
    out_cluster_end   <- rep(NA_integer_,   n)
    
    # Assign clusters per chromosome using findInterval (fast and stable)
    chrs <- unique(df2$.chr)
    for (cc in chrs) {
      
      idx_df <- which(!is.na(df2$.chr) & df2$.chr == cc & !is.na(df2$.pos))
      if (!length(idx_df)) next
      
      cl_cc <- cl2[cl2$chr == cc & !is.na(cl2$start) & !is.na(cl2$end), , drop = FALSE]
      if (!nrow(cl_cc)) next
      
      # IMPORTANT: starts ha d'estar ordenat per findInterval()
      cl_cc <- cl_cc[order(cl_cc$start, cl_cc$end), , drop = FALSE]
      
      starts <- cl_cc$start
      ends   <- cl_cc$end
      pos    <- df2$.pos[idx_df]
      
      ii <- findInterval(pos, starts)
      valid <- ii > 0 & pos <= ends[ii]
      if (!any(valid)) next
      
      sel_df <- idx_df[valid]
      sel_cl <- ii[valid]
      
      out_cluster_id[sel_df]    <- cl_cc$cluster_id[sel_cl]
      out_cluster_chr_n[sel_df] <- cl_cc$cluster_chr_n[sel_cl]
      out_cluster_start[sel_df] <- cl_cc$start[sel_cl]
      out_cluster_end[sel_df]   <- cl_cc$end[sel_cl]
    }
    
    df2$cluster_id    <- out_cluster_id
    df2$cluster_chr_n <- out_cluster_chr_n
    df2$cluster_start <- out_cluster_start
    df2$cluster_end   <- out_cluster_end
    df2$in_cluster    <- !is.na(df2$cluster_id)
    
    # Cleanup temp columns
    df2$.chr <- NULL
    df2$.pos <- NULL
    
    df2
  }
  
  
  
  # ---- Override GWAS reactive: prefer shared GWAS if present ----
  #  gwas_df_local <- gwas_df
  #  gwas_df <- reactive({
  #    df <- gwas_shared_rv()
  #    if (is.data.frame(df) && nrow(df) > 0) return(df)
  #    gwas_df_local()
  #  })
  
  
  # ============================================================
  # STEP NEW · Extract NONSYN variants in all clusters
  # (from a precomputed dbNSFP .out file)
  # ============================================================
  
  # Helper: robustly locate chr/pos columns in a dbNSFP output header
  guess_dbnsfp_chr_pos_cols <- function(cn) {
    cn0 <- cn
    cn  <- tolower(cn)
    
    chr_idx <- which(cn %in% c("#chr","chr","chrom","chromosome","chrm","chromosome_name"))
    pos_idx <- which(cn %in% c("pos(1-based)","pos","position","bp","start","hg38_pos","pos_1_based"))
    
    if (!length(chr_idx)) {
      chr_idx <- grep("(^#?chr$)|(^chrom$)|(^chromosome$)", cn)
    }
    if (!length(pos_idx)) {
      pos_idx <- grep("^pos\\(1-based\\)$|^pos$|^position$|^bp$", cn)
    }
    
    list(
      chr = if (length(chr_idx)) cn0[chr_idx[1]] else NULL,
      pos = if (length(pos_idx)) cn0[pos_idx[1]] else NULL
    )
  }
  
  
  ##--- helpers (keep near server) ---
  
  read_first_line_any <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      con <- gzfile(path, open = "rt")
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      readLines(con, n = 1, warn = FALSE)
    } else {
      readLines(path, n = 1, warn = FALSE)
    }
  }
  
  norm_chr_int_dbnsfp <- function(x) {
    x <- toupper(as.character(x))
    x <- sub("^CHR", "", x)
    x <- sub("^CHROM", "", x)
    x <- sub("^CHROMOSOME", "", x)
    x[x == "X"] <- "23"
    x[x == "Y"] <- "24"
    x[x %in% c("MT","M")] <- "26"
    suppressWarnings(as.integer(x))
  }
  
  save_dbnsfp_final_local <- function(df, out_dir, stem = "dbnsfp_normalized") {
    stopifnot(is.data.frame(df))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_csv <- file.path(out_dir, paste0(stem, ".csv"))
    out_rds <- file.path(out_dir, paste0(stem, ".rds"))
    readr::write_csv(df, out_csv)
    saveRDS(df, out_rds)
    list(csv = out_csv, rds = out_rds)
  }
  
  ##################
  
  observeEvent(input$run_nonsyn_clusters, {
    
    showNotification("START · Extract NonSyn variants in clusters", type = "message", duration = 2)
    
    safe_log <- function(...) {
      if (exists("append_dbnsfp_log", mode = "function")) {
        append_dbnsfp_log(...)
      } else {
        message(paste0(..., collapse = ""))
      }
    }
    
    safe_log("[STEP] run_nonsyn_clusters clicked\n")
    
    tryCatch({
      
      validate(need(exists("workdir") && nzchar(workdir), "workdir is not defined."))
      if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
      
      tabix_path <- Sys.which("tabix")
      safe_log("[INFO] tabix= ", tabix_path, " \n")
      validate(need(nzchar(tabix_path), "tabix not found in PATH (Sys.which('tabix') is empty)."))
      
      # --------------------------
      # Build regions for clusters (from clusters_cur)
      # --------------------------
      cl_cur <- clusters_cur()  # canònic
      validate(need(is.data.frame(cl_cur) && nrow(cl_cur) > 0, "No clusters available."))
      
      cl_cur <- as.data.frame(cl_cur)
      if (!"chr"   %in% names(cl_cur) && "CHR"      %in% names(cl_cur)) cl_cur$chr   <- cl_cur$CHR
      if (!"start" %in% names(cl_cur) && "start_bp" %in% names(cl_cur)) cl_cur$start <- cl_cur$start_bp
      if (!"end"   %in% names(cl_cur) && "end_bp"   %in% names(cl_cur)) cl_cur$end   <- cl_cur$end_bp
      
      cl_cur$chr_num <- suppressWarnings(as.integer(gsub("^chr", "", as.character(cl_cur$chr), ignore.case = TRUE)))
      cl_cur$start_i <- suppressWarnings(as.integer(cl_cur$start))
      cl_cur$end_i   <- suppressWarnings(as.integer(cl_cur$end))
      
      cl_cur <- cl_cur[
        is.finite(cl_cur$chr_num) & is.finite(cl_cur$start_i) & is.finite(cl_cur$end_i) &
          (cl_cur$end_i >= cl_cur$start_i),
        , drop = FALSE
      ]
      validate(need(nrow(cl_cur) > 0, "Clusters parsed but 0 valid regions."))
      
      # -------------------------
      # dbNSFP ALL file (bgzip + tabix)  [PATH HANDLING]
      # -------------------------
      dbnsfp_all_raw <- input$dbnsfp_out_all
      dbnsfp_all_gz  <- trimws(gsub('^"|"$', "", as.character(dbnsfp_all_raw)))
      
      safe_log("[INFO] dbnsfp_out_all(raw)= ", as.character(dbnsfp_all_raw), " \n")
      
      # Fix common typo: Inspector_app_slave -> Inspector_app_slaves
      dbnsfp_all_gz <- gsub("/Inspector_app_slave/", "/Inspector_app_slaves/", dbnsfp_all_gz, fixed = TRUE)
      
      # If empty, set default known path
      if (is.null(dbnsfp_all_gz) || is.na(dbnsfp_all_gz) || !nzchar(dbnsfp_all_gz)) {
        dbnsfp_all_gz <- ns_def$dbnsfp_all
      }
      
      # If user provided directory, pick an ALL file inside
      if (dir.exists(dbnsfp_all_gz)) {
        cand <- list.files(dbnsfp_all_gz, pattern = "\\.ALL\\.out\\.gz$", full.names = TRUE)
        if (!length(cand)) cand <- list.files(dbnsfp_all_gz, pattern = "\\.out\\.gz$", full.names = TRUE)
        validate(need(length(cand) > 0, "dbNSFP directory provided, but no .out.gz found inside."))
        cand <- cand[order(nchar(cand), decreasing = TRUE)]
        dbnsfp_all_gz <- cand[1]
        safe_log("[INFO] dbnsfp_all_from_dir= ", dbnsfp_all_gz, " \n")
      }
      
      safe_log("[INFO] file.exists(dbnsfp_all)= ", file.exists(dbnsfp_all_gz), " \n")
      validate(need(nzchar(dbnsfp_all_gz) && file.exists(dbnsfp_all_gz),
                    "dbNSFP .out.gz file not found (dbnsfp_out_all). Check path."))
      validate(need(isTRUE(grepl("\\.gz$", dbnsfp_all_gz, ignore.case = TRUE)),
                    "Provide a bgzip-compressed .gz file (tabix indexed)."))
      
      # --------------------------
      # Ensure tabix index exists (NA-safe)
      # --------------------------
      tbi1 <- paste0(dbnsfp_all_gz, ".tbi")
      tbi2 <- sub("\\.gz$", ".tbi", dbnsfp_all_gz)
      
      has_tbi1 <- isTRUE(file.exists(tbi1))
      has_tbi2 <- isTRUE(file.exists(tbi2))
      
      if (!has_tbi1 && !has_tbi2) {
        cat("[WARN] No tabix index found. Building index...\n")
        
        cmd_i <- sprintf('"%s" -f -s 1 -b 2 -e 2 "%s"', tabix_path, dbnsfp_all_gz)
        cat("[INFO] tabix index cmd:", cmd_i, "\n")
        rc <- system(cmd_i)
        
        has_tbi1 <- isTRUE(file.exists(tbi1))
        has_tbi2 <- isTRUE(file.exists(tbi2))
        validate(need(has_tbi1 || has_tbi2, "tabix indexing failed: .tbi was not created."))
      }
      
      # -------------------------
      # Helpers: read gz header / first data line
      # -------------------------
      read_first_line_gz <- function(path) {
        con <- gzfile(path, open = "rt")
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        readLines(con, n = 1, warn = FALSE)
      }
      
      read_first_data_line_gz <- function(path) {
        con <- gzfile(path, open = "rt")
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        x <- readLines(con, n = 3, warn = FALSE)
        x <- x[!grepl("^\\s*$", x)]
        if (length(x) < 2) return(NA_character_)
        x[2]
      }
      
      # Detect if dbNSFP chr column uses "chr" prefix
      first_data <- read_first_data_line_gz(dbnsfp_all_gz)
      validate(need(!is.na(first_data), "Could not read first data line from dbNSFP .gz."))
      first_chr <- strsplit(first_data, "\t", fixed = TRUE)[[1]][1]
      file_has_chr_prefix <- isTRUE(grepl("^chr", tolower(first_chr)))
      safe_log("[INFO] dbNSFP chr_prefix= ", file_has_chr_prefix, "  first_chr= ", first_chr, " \n")
      
      # -------------------------
      # Normalize clusters -> chr/start_bp/end_bp/cluster_id (from any source)
      # -------------------------
      get_clusters_any <- function() {
        cl <- NULL
        if (exists("clusters_cur", mode = "function")) {
          cl <- tryCatch(clusters_cur(), error = function(e) NULL)
        }
        if ((is.null(cl) || !is.data.frame(cl) || nrow(cl) == 0) && exists("clusters_src", mode = "function")) {
          cl <- tryCatch(clusters_src(), error = function(e) NULL)
        }
        if ((is.null(cl) || !is.data.frame(cl) || nrow(cl) == 0) && exists("clusters_rv", mode = "function")) {
          cl <- tryCatch(clusters_rv(), error = function(e) NULL)
        }
        if ((is.null(cl) || !is.data.frame(cl) || nrow(cl) == 0) && exists("rv", inherits = TRUE)) {
          cl <- tryCatch(rv$clusters, error = function(e) NULL)
        }
        if ((is.null(cl) || !is.data.frame(cl) || nrow(cl) == 0) && exists("cluster_dt_view2", mode = "function")) {
          cl <- tryCatch(cluster_dt_view2(), error = function(e) NULL)
        }
        cl
      }
      
      cl0 <- get_clusters_any()
      if (is.null(cl0) || !is.data.frame(cl0) || nrow(cl0) == 0) {
        showNotification("STOP · No clusters available. Run 'Generate intervals & clusters' first.", type = "error", duration = 10)
        safe_log("[STOP] clusters source empty\n")
        return()
      }
      
      cl0 <- as.data.frame(cl0)
      
      # --- CANONICAL: always create cluster_id ---
      if (!"cluster_id" %in% names(cl0)) cl0$cluster_id <- NA_character_
      
      if ("cluster_chr_n" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          as.character(cl0$cluster_chr_n))
      }
      
      if ("label" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          as.character(cl0$label))
      }
      if ("cluster" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          paste0("cluster_", as.character(cl0$cluster)))
      }
      if (all(is.na(cl0$cluster_id) | !nzchar(as.character(cl0$cluster_id)))) {
        cl0$cluster_id <- paste0("cluster_", seq_len(nrow(cl0)))
      }
      
      nm <- names(cl0)
      
      id_col <- dplyr::case_when(
        "cluster_id" %in% nm ~ "cluster_id",
        "label" %in% nm ~ "label",
        TRUE ~ NA_character_
      )
      
      chr_col <- dplyr::case_when(
        "chr" %in% nm ~ "chr",
        "CHR" %in% nm ~ "CHR",
        "chr_int" %in% nm ~ "chr_int",
        "chr_num" %in% nm ~ "chr_num",
        "chr_lbl" %in% nm ~ "chr_lbl",
        TRUE ~ NA_character_
      )
      
      st_col <- dplyr::case_when(
        "start_bp" %in% nm ~ "start_bp",
        "start" %in% nm ~ "start",
        "bp1" %in% nm ~ "bp1",
        "start_pos" %in% nm ~ "start_pos",
        TRUE ~ NA_character_
      )
      
      en_col <- dplyr::case_when(
        "end_bp" %in% nm ~ "end_bp",
        "end" %in% nm ~ "end",
        "bp2" %in% nm ~ "bp2",
        "end_pos" %in% nm ~ "end_pos",
        TRUE ~ NA_character_
      )
      
      id_col <- dplyr::case_when(
        "cluster_id" %in% nm ~ "cluster_id",
        "label" %in% nm ~ "label",
        "cluster_chr_n" %in% nm ~ "cluster_chr_n",
        TRUE ~ NA_character_
      )
      
      validate(
        need(!is.na(chr_col), "Cannot find chromosome column in clusters."),
        need(!is.na(st_col),  "Cannot find start column in clusters."),
        need(!is.na(en_col),  "Cannot find end column in clusters."),
        need(!is.na(id_col),  "Cannot find cluster id column in clusters.")
      )
      
      chr_raw <- as.character(cl0[[chr_col]])
      chr_num <- suppressWarnings(as.integer(gsub("^chr", "", chr_raw, ignore.case = TRUE)))
      
      cl_norm <- tibble::tibble(
        chr        = chr_num,
        start_bp   = suppressWarnings(as.integer(cl0[[st_col]])),
        end_bp     = suppressWarnings(as.integer(cl0[[en_col]])),
        cluster_id = as.character(cl0[[id_col]])
      ) %>%
        dplyr::filter(
          is.finite(.data$chr),
          is.finite(.data$start_bp),
          is.finite(.data$end_bp),
          .data$end_bp >= .data$start_bp,
          !is.na(.data$cluster_id), nzchar(.data$cluster_id)
        ) %>%
        dplyr::distinct(.data$cluster_id, .data$chr, .data$start_bp, .data$end_bp) %>%
        dplyr::arrange(.data$chr, .data$start_bp, .data$end_bp)
      
      validate(need(nrow(cl_norm) > 0, "After normalizing clusters, no valid clusters remained."))
      safe_log("[INFO] clusters_n= ", nrow(cl_norm), " \n")
      
      # -------------------------
      # cl_assign: compatibility for add_clusters_to_dbnsfp()
      # expects: chr, start, end, cluster (+ cluster_id ok)
      # -------------------------
      cl_assign <- cl_norm %>%
        dplyr::group_by(.data$chr) %>%
        dplyr::arrange(.data$start_bp, .data$end_bp, .by_group = TRUE) %>%
        dplyr::mutate(
          cluster = {
            x <- suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", .data$cluster_id)))
            ifelse(is.finite(x), x, dplyr::row_number())
          }
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          start = as.integer(.data$start_bp),
          end   = as.integer(.data$end_bp)
        )
      
      if (any(!is.finite(cl_assign$cluster))) {
        cl_assign <- cl_assign %>%
          dplyr::group_by(.data$chr) %>%
          dplyr::arrange(.data$start_bp, .data$end_bp, .by_group = TRUE) %>%
          dplyr::mutate(cluster = dplyr::row_number()) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(start = as.integer(.data$start_bp), end = as.integer(.data$end_bp))
      }
      
      # -------------------------
      # Prepare subset .out (unique header)
      # -------------------------
      subset_out <- file.path(workdir, "dbnsfp_subset_in_clusters.out")
      if (file.exists(subset_out)) file.remove(subset_out)
      
      header_line_raw <- read_first_line_gz(dbnsfp_all_gz)
      validate(need(length(header_line_raw) == 1 && nzchar(header_line_raw),
                    "Could not read header line from dbNSFP .gz."))
      
      header_names_raw <- strsplit(header_line_raw, "\t", fixed = TRUE)[[1]]
      validate(need(length(header_names_raw) >= 2, "Header is not tab-delimited."))
      
      header_names_u <- make.unique(header_names_raw)
      writeLines(paste(header_names_u, collapse = "\t"), subset_out)
      
      # -------------------------
      # Build tabix regions (dbNSFP)
      # -------------------------
      reg_chr <- if (isTRUE(file_has_chr_prefix)) paste0("chr", cl_assign$chr) else as.character(cl_assign$chr)
      regions <- paste0(reg_chr, ":", as.integer(cl_assign$start_bp), "-", as.integer(cl_assign$end_bp))
      regions <- regions[!is.na(regions) & nzchar(regions)]
      validate(need(length(regions) > 0, "No regions constructed from clusters."))
      
      safe_log("[INFO] regions_n= ", length(regions), " \n")
      showNotification(paste0("RUN · tabix querying ", length(regions), " regions..."), type = "message", duration = 4)
      
      # -------------------------
      # Run tabix (group by chromosome for speed)
      # -------------------------
      total_lines <- 0L
      chr_keys <- unique(reg_chr)
      
      withProgress(message = "Extracting dbNSFP rows in clusters (tabix)", value = 0, {
        nchr <- length(chr_keys)
        for (ii in seq_along(chr_keys)) {
          incProgress(0.9 / max(1, nchr), detail = chr_keys[ii])
          
          regs_chr <- regions[reg_chr == chr_keys[ii]]
          if (!length(regs_chr)) next
          
          out <- tryCatch(
            system2(tabix_path, args = c(dbnsfp_all_gz, regs_chr), stdout = TRUE, stderr = TRUE),
            error = function(e) character(0)
          )
          
          if (!length(out)) next
          out <- out[!grepl("^\\s*$", out)]
          out <- out[!grepl("^tabix:", out)]
          
          if (length(out)) {
            write(out, file = subset_out, append = TRUE)
            total_lines <- total_lines + length(out)
          }
        }
        incProgress(0.1, detail = "Done")
      })
      
      safe_log("[INFO] tabix_total_lines= ", total_lines, " \n")
      validate(need(total_lines > 0, "tabix returned 0 lines (no overlaps or wrong indexing/region format)."))
      
      # -------------------------
      # Normalize subset + assign clusters
      # -------------------------
      validate(need(exists("normalize_dbnsfp", mode = "function"),
                    "normalize_dbnsfp() not available (format_dbnsfp.R not loaded)."))
      
      res <- normalize_dbnsfp(in_file = subset_out, out_dir = workdir)
      validate(need(is.list(res) && !is.null(res$rds) && file.exists(res$rds),
                    "normalize_dbnsfp() did not produce an RDS."))
      
      df_norm <- readRDS(res$rds)
      validate(need(is.data.frame(df_norm) && nrow(df_norm) > 0, "Normalized subset is empty."))
      
      safe_log("\n====================\n")
      safe_log("[df_norm]  nrow=", nrow(df_norm), " ncol=", ncol(df_norm), "\n")
      safe_log("[df_norm]  cols: ", paste(names(df_norm), collapse = ", "), "\n")
      safe_log("====================\n")
      safe_log(
        paste0(
          "[df_norm] head:\n",
          paste(capture.output(print(utils::head(df_norm, 5))), collapse = "\n"),
          "\n"
        )
      )
      
      validate(need(exists("add_clusters_to_dbnsfp", mode = "function"),
                    "add_clusters_to_dbnsfp() not available."))
      
      df_final <- add_clusters_to_dbnsfp(df_norm, cl_assign)
      df_final <- df_final %>% dplyr::filter(.data$in_cluster)
      validate(need(nrow(df_final) > 0, "No normalized rows remained after cluster assignment."))
      
      safe_log("\n====================\n")
      safe_log("[df_final]  nrow=", nrow(df_final), " ncol=", ncol(df_final), "\n")
      safe_log("[df_final]  cols: ", paste(names(df_final), collapse = ", "), "\n")
      safe_log("====================\n")
      safe_log(
        paste0(
          "[df_final] head:\n",
          paste(capture.output(print(utils::head(df_final, 5))), collapse = "\n"),
          "\n"
        )
      )
      
      # -------------------------
      # Save final outputs
      # -------------------------
      tg <- make_mode_thr_tag(
        isolate(input$cluster_method %||% "window"),
        isolate(input$pthr),
        isolate(input$min_logp)
      )
      
      mode_tag <- as.character(tg$mode_tag %||% "run")
      thr_txt  <- as.character(tg$thr_txt  %||% "NA")
      if (length(mode_tag) != 1 || !nzchar(mode_tag)) mode_tag <- "run"
      if (length(thr_txt)  != 1 || !nzchar(thr_txt))  thr_txt  <- "NA"
      
      stem <- paste0("dbnsfp_in_clusters_", mode_tag, "_thr", thr_txt)
      
      save_dbnsfp_final_local <- function(df, out_dir, stem = "dbnsfp_normalized") {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        out_csv <- file.path(out_dir, paste0(stem, ".csv"))
        out_rds <- file.path(out_dir, paste0(stem, ".rds"))
        readr::write_csv(df, out_csv)
        saveRDS(df, out_rds)
        list(csv = out_csv, rds = out_rds)
      }
      
      final_paths <- save_dbnsfp_final_local(df_final, out_dir = workdir, stem = stem)
      
      if (exists("dbnsfp_final_path_rds", mode = "function") &&
          is.character(final_paths$rds) && length(final_paths$rds) == 1 && nzchar(final_paths$rds)) {
        dbnsfp_final_path_rds(final_paths$rds)
      }
      
      if (exists("dbnsfp_final_path_csv", mode = "function") &&
          is.character(final_paths$csv) && length(final_paths$csv) == 1 && nzchar(final_paths$csv)) {
        dbnsfp_final_path_csv(final_paths$csv)
      }
      
      if (exists("dbnsfp_norm_path_rds", mode = "function") &&
          is.character(res$rds) && length(res$rds) == 1 && nzchar(res$rds)) {
        dbnsfp_norm_path_rds(res$rds)
      }
      
      if (exists("dbnsfp_norm_path_csv", mode = "function") &&
          !is.null(res$csv) && is.character(res$csv) && length(res$csv) == 1 && nzchar(res$csv)) {
        dbnsfp_norm_path_csv(res$csv)
      }
      
      ###### Start INTEGRATOR CODE ####################################################################
      # ------------------------------------------------------------
      # INTEGRATOR export (NonSyn bridges)
      # ------------------------------------------------------------
      tryCatch({
        
        # ============================================================
        # APP CONFIG (NonSyn)
        # ============================================================
        app_slug      <- "nonsyn"
        app_hit_class <- "nonsyn_hit"
        
        # Main table in this app
        hits <- df_final
        
        # Generic logger fallback
        append_log <- if (exists("append_log", mode = "function", inherits = TRUE)) {
          get("append_log", mode = "function", inherits = TRUE)
        } else if (exists("safe_log", mode = "function", inherits = TRUE)) {
          function(x) safe_log(x, "\n")
        } else {
          function(x) cat(x, "\n", sep = "")
        }
        
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
        
        get_nonsyn_cluster_method <- function(input) {
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
        
        get_nonsyn_threshold <- function(input, cluster_method) {
          if (identical(cluster_method, "window")) {
            suppressWarnings(as.numeric(input$pthr))
          } else if (cluster_method %in% c("hits_span1mb", "hits_tiled", "hits_sliding", "hits_unknown")) {
            suppressWarnings(as.numeric(input$min_logp))
          } else {
            NA_real_
          }
        }
        
        get_nonsyn_gwas_name <- function(input) {
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
        current_cluster_method <- get_nonsyn_cluster_method(input)
        current_threshold      <- get_nonsyn_threshold(input, current_cluster_method)
        current_gwas_name      <- get_nonsyn_gwas_name(input)
        
        sess <- gi_get_integrator_session_dir(
          gi_shared_root = gi_shared_root,
          cluster_method = current_cluster_method,
          threshold_used = current_threshold
        )
        
        integrator_root <- sess$integrator_root
        session_id      <- sess$session_id
        session_dir     <- sess$session_dir
        
        manifest_obj <- gi_load_or_create_manifest(
          session_dir = session_dir,
          session_id = session_id,
          cluster_method = current_cluster_method,
          threshold_used = current_threshold,
          gwas_session_file = current_gwas_name
        )
        
        manifest      <- manifest_obj$manifest
        manifest_path <- manifest_obj$manifest_path
        
        if (!identical(as.character(manifest$cluster_method %||% ""), as.character(current_cluster_method))) {
          stop("Existing manifest cluster_method does not match current cluster_method.")
        }
        
        if (!isTRUE(all.equal(
          suppressWarnings(as.numeric(manifest$threshold_used)),
          suppressWarnings(as.numeric(current_threshold))
        ))) {
          stop("Existing manifest threshold_used does not match current threshold_used.")
        }
        
        manifest_gwas_files <- as.character(manifest$gwas_session_file %||% character(0))
        manifest_gwas_files <- trimws(manifest_gwas_files)
        manifest_gwas_files <- manifest_gwas_files[!is.na(manifest_gwas_files) & nzchar(manifest_gwas_files)]
        
        current_gwas_name_chr <- trimws(as.character(current_gwas_name %||% ""))
        if (nzchar(current_gwas_name_chr)) {
          manifest_gwas_files <- unique(c(manifest_gwas_files, current_gwas_name_chr))
        }
        
        manifest$gwas_session_file <- manifest_gwas_files
        
        gene_hit_bridge_path <- file.path(session_dir, paste0(app_slug, "_gene_hit_bridge.rds"))
        gene_bridge_path     <- file.path(session_dir, paste0(app_slug, "_gene_bridge.rds"))
        term_bridge_path     <- file.path(session_dir, paste0(app_slug, "_term_bridge.rds"))
        clusters_master_path <- file.path(session_dir, paste0(app_slug, "_clusters_master.rds"))
        candidates_path      <- file.path(session_dir, paste0(app_slug, "_candidates.rds"))
        
        append_log(paste0(log_tag, " integrator_root=", integrator_root))
        append_log(paste0(log_tag, " cluster_method=", current_cluster_method))
        append_log(paste0(log_tag, " threshold_used=", current_threshold))
        append_log(paste0(log_tag, " gwas_session_file=", paste(manifest$gwas_session_file, collapse = " | ")))
        append_log(paste0(log_tag, " session_id=", session_id))
        append_log(paste0(log_tag, " session_dir=", session_dir))
        append_log(paste0(log_tag, " hits nrow=", nrow(hits), " ncol=", ncol(hits)))
        append_log(paste0(log_tag, " hits cols: ", paste(names(hits), collapse = ", ")))
        
        # ============================================================
        # A) GENE BRIDGE
        # ============================================================
        needed_gene_cols <- c("genename", "cluster_id", "chr", "cluster_start", "cluster_end")
        
        gene_bridge <- if (all(needed_gene_cols %in% names(hits))) {
          hits %>%
            dplyr::transmute(
              gene = as.character(genename),
              source_app = app_slug,
              evidence_type = paste0(app_slug, "_gene"),
              cluster_id = as.character(cluster_id),
              chr = chr,
              start = cluster_start,
              end = cluster_end
            ) %>%
            sanitize_bridge() %>%
            sanitize_bridge_genes() %>%
            dplyr::filter(!is.na(gene), nzchar(gene)) %>%
            dplyr::distinct()
        } else {
          append_log(paste0(
            log_tag, " gene bridge skipped: missing columns -> ",
            paste(setdiff(needed_gene_cols, names(hits)), collapse = ", ")
          ))
          tibble::tibble()
        }
        
        append_log(paste0(log_tag, " gene_bridge nrow=", nrow(gene_bridge)))
        
        if (nrow(gene_bridge) > 0) {
          saveRDS(gene_bridge, gene_bridge_path)
          manifest$files[[paste0(app_slug, "_gene_bridge")]] <- basename(gene_bridge_path)
          append_log(paste0(
            log_tag, " saved ", basename(gene_bridge_path),
            " | exists=", file.exists(gene_bridge_path)
          ))
        } else {
          append_log(paste0(log_tag, " ", basename(gene_bridge_path), " NOT saved (0 rows)"))
        }
        
        # ============================================================
        # B) TERM BRIDGE
        # ============================================================
        # ============================================================
        # B) TERM BRIDGE
        # ============================================================
        
        term_bridge <- tibble::tibble(
          term = character(),
          term_type = character(),
          source_app = character(),
          evidence_type = character(),
          cluster_id = character(),
          chr = integer(),
          start = integer(),
          end = integer()
        )
        
        split_term_values <- function(x) {
          x <- as.character(x)
          x <- x[!is.na(x)]
          if (!length(x)) return(character())
          
          x <- trimws(x)
          x <- x[nzchar(x)]
          if (!length(x)) return(character())
          
          parts <- unlist(strsplit(x, "\\s*[;|]\\s*", perl = TRUE), use.names = FALSE)
          parts <- trimws(parts)
          parts <- parts[nzchar(parts)]
          unique(parts)
        }
        
        clean_term_value <- function(x, evidence_type) {
          x <- as.character(x)
          x <- trimws(x)
          x <- gsub("\\s+", " ", x, perl = TRUE)
          
          bad_vals <- c(
            "", ".", "-", "NA", "N/A", "na", "n/a",
            "none", "None", "not provided", "not_provided",
            "not_specified", "unknown", "Unknown", "unspecified", "Unspecified"
          )
          if (is.na(x) || x %in% bad_vals) return(NA_character_)
          
          x <- gsub("^[\"']+|[\"']+$", "", x, perl = TRUE)
          
          if (identical(evidence_type, "mim_disease")) {
            x <- gsub("\\[MIM\\s*:\\s*[0-9]+\\]", "", x, ignore.case = TRUE, perl = TRUE)
            x <- gsub("\\b(OMIM|MIM)\\s*[:#]*\\s*[0-9]+\\b", "", x, ignore.case = TRUE, perl = TRUE)
          }
          
          if (identical(evidence_type, "gwas_trait")) {
            # traiem PMIDs entre claudàtors
            x <- gsub("\\[[^\\]]*\\]", "", x, perl = TRUE)
          }
          
          x <- gsub("^[:;,_\\-\\s]+|[:;,_\\-\\s]+$", "", x, perl = TRUE)
          x <- trimws(x)
          if (!nzchar(x)) return(NA_character_)
          
          x_low <- tolower(x)
          generic_bad <- c(
            "disease", "diseases", "trait", "traits", "phenotype", "phenotypes",
            "clinical significance", "pathogenic", "likely pathogenic",
            "benign", "likely benign", "uncertain significance",
            "risk factor", "association", "not specified"
          )
          if (x_low %in% generic_bad) return(NA_character_)
          
          x
        }
        
        build_term_bridge_from_col <- function(df, col, term_type, evidence_type) {
          if (!(col %in% names(df))) {
            append_log(paste0(log_tag, " term source skipped: column not found -> ", col))
            return(NULL)
          }
          
          vals0 <- as.character(df[[col]])
          vals0 <- vals0[!is.na(vals0)]
          vals0 <- trimws(vals0)
          vals0 <- vals0[nzchar(vals0)]
          
          append_log(paste0(
            log_tag, " term source ", col,
            " | non_empty=", length(vals0),
            " | unique_non_empty=", length(unique(vals0))
          ))
          
          if (!length(vals0)) {
            append_log(paste0(log_tag, " built 0 rows from -> ", col))
            return(NULL)
          }
          
          out_list <- vector("list", nrow(df))
          
          for (i in seq_len(nrow(df))) {
            vals <- split_term_values(df[[col]][i])
            if (!length(vals)) next
            
            vals <- vapply(
              vals,
              FUN = clean_term_value,
              evidence_type = evidence_type,
              FUN.VALUE = character(1)
            )
            vals <- vals[!is.na(vals) & nzchar(vals)]
            vals <- unique(vals)
            if (!length(vals)) next
            
            out_list[[i]] <- tibble::tibble(
              term = vals,
              term_type = term_type,
              source_app = app_slug,
              evidence_type = evidence_type,
              cluster_id = as.character(df$cluster_id[i]),
              chr = as.integer(df$chr[i]),
              start = as.integer(df$cluster_start[i]),
              end = as.integer(df$cluster_end[i])
            )
          }
          
          out <- dplyr::bind_rows(out_list)
          
          if (!nrow(out)) {
            append_log(paste0(log_tag, " built 0 rows from -> ", col))
            return(NULL)
          }
          
          append_log(paste0(log_tag, " built ", nrow(out), " rows from -> ", col))
          out
        }
        
        tryCatch({
          
          append_log(paste0(log_tag, " hits nrow=", nrow(hits)))
          append_log(paste0(log_tag, " hits cols=", paste(names(hits), collapse = ", ")))
          
          # ------------------------------------------------
          # Fonts per a nonsyn shared terms
          # ------------------------------------------------
          # Mantinc MIM_disease i Trait_association_GWAS com a principals.
          # Disease_description i clinvar_trait es deixen fora del principal
          # perquè introdueixen massa soroll textual.
          term_specs <- list(
            list(col = "MIM_disease",            term_type = "disease", evidence_type = "mim_disease"),
            list(col = "Trait_association_GWAS", term_type = "trait",   evidence_type = "gwas_trait"),
            list(col = "Orphanet_disorder",      term_type = "disease", evidence_type = "orphanet_disorder"),
            list(col = "HPO_name",               term_type = "trait",   evidence_type = "hpo_name")
          )
          
          term_bridge_list <- lapply(term_specs, function(sp) {
            build_term_bridge_from_col(
              df = hits,
              col = sp$col,
              term_type = sp$term_type,
              evidence_type = sp$evidence_type
            )
          })
          
          term_bridge_list <- term_bridge_list[!vapply(term_bridge_list, is.null, logical(1))]
          
          if (length(term_bridge_list) > 0) {
            term_bridge <- dplyr::bind_rows(term_bridge_list) %>%
              sanitize_bridge() %>%
              dplyr::filter(!is.na(term), nzchar(term)) %>%
              dplyr::mutate(
                term = trimws(as.character(term)),
                term_type = as.character(term_type),
                source_app = as.character(source_app),
                evidence_type = as.character(evidence_type),
                cluster_id = as.character(cluster_id),
                chr = as.integer(chr),
                start = as.integer(start),
                end = as.integer(end)
              ) %>%
              dplyr::distinct()
          }
          
          append_log(paste0(log_tag, " term_bridge nrow=", nrow(term_bridge)))
          
          if (nrow(term_bridge) > 0) {
            append_log(
              paste0(
                log_tag, " term_bridge evidence counts: ",
                paste(names(table(term_bridge$evidence_type)), table(term_bridge$evidence_type), collapse = " | ")
              )
            )
          }
          
        }, error = function(e) {
          append_log(paste0(log_tag, " ERROR building term_bridge: ", conditionMessage(e)))
        })
        
        saveRDS(term_bridge, term_bridge_path)
        manifest$files[[paste0(app_slug, "_term_bridge")]] <- basename(term_bridge_path)
        append_log(paste0(
          log_tag, " saved ", if (nrow(term_bridge) > 0) "" else "empty ",
          basename(term_bridge_path),
          " | exists=", file.exists(term_bridge_path)
        ))
        
        # ============================================================
        # C) CLUSTERS MASTER
        # ============================================================
        append_log(paste0(ld_tag, " dir=", session_dir))
        append_log(paste0(ld_tag, " dir exists: ", dir.exists(session_dir)))
        
        clusters_master <- tibble::tibble()
        
        # Priority 1: use cl2 if available
        if (exists("cl2", inherits = TRUE) && is.data.frame(cl2) && nrow(cl2) > 0) {
          
          append_log(paste0(ld_tag, "[CLUSTERS] using cl2"))
          append_log(paste0(ld_tag, "[CLUSTERS] cl2 nrow=", nrow(cl2), " ncol=", ncol(cl2)))
          append_log(paste0(ld_tag, "[CLUSTERS] cl2 cols: ", paste(names(cl2), collapse = ", ")))
          
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
          
          if (!is.null(chr_col) && !is.null(st_col) && !is.null(en_col) && !is.null(id_col)) {
            clusters_master <- cl2 %>%
              dplyr::transmute(
                cluster_id = as.character(.data[[id_col]]),
                chr        = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[chr_col]])))),
                start_raw  = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
                end_raw    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]])))),
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
          } else {
            append_log(paste0(ld_tag, "[CLUSTERS] cl2 missing required columns"))
          }
          
        } else if (all(c("cluster_id", "chr", "cluster_start", "cluster_end") %in% names(hits))) {
          
          # Priority 2: derive from df_final
          append_log(paste0(ld_tag, "[CLUSTERS] deriving clusters_master from hits"))
          
          clusters_master <- hits %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              chr        = suppressWarnings(as.integer(readr::parse_number(as.character(chr)))),
              start_raw  = suppressWarnings(as.integer(readr::parse_number(as.character(cluster_start)))),
              end_raw    = suppressWarnings(as.integer(readr::parse_number(as.character(cluster_end))))
            ) %>%
            dplyr::mutate(
              cluster_id = trimws(cluster_id),
              start = pmin(start_raw, end_raw),
              end   = pmax(start_raw, end_raw),
              cluster_key = paste0(cluster_id, "|chr", chr, ":", start, "-", end)
            ) %>%
            dplyr::select(cluster_id, chr, start, end, cluster_key) %>%
            dplyr::filter(
              !is.na(cluster_id), nzchar(cluster_id),
              is.finite(chr), is.finite(start), is.finite(end)
            ) %>%
            dplyr::distinct()
        } else {
          append_log(paste0(ld_tag, "[CLUSTERS] no usable source found for clusters_master"))
        }
        
        append_log(paste0(ld_tag, "[CLUSTERS] clusters_master nrow=", nrow(clusters_master)))
        
        if (nrow(clusters_master) > 0) {
          saveRDS(clusters_master, clusters_master_path)
          manifest$files[[paste0(app_slug, "_clusters_master")]] <- basename(clusters_master_path)
          append_log(paste0(ld_tag, "[CLUSTERS] saved: ", clusters_master_path))
          append_log(paste0(ld_tag, "[CLUSTERS] exists after save: ", file.exists(clusters_master_path)))
        } else {
          append_log(paste0(ld_tag, "[CLUSTERS] ", basename(clusters_master_path), " NOT saved (0 rows)"))
        }
        
        # ============================================================
        # D) CANDIDATES FOR LD = GWAS + NONSYN HIT
        # ============================================================
        get_gwas_hits_for_ld <- function() {
          if (exists("hits_df", mode = "function", inherits = TRUE)) {
            x <- tryCatch(hits_df(), error = function(e) NULL)
            if (is.data.frame(x) && nrow(x) > 0) {
              return(list(name = "hits_df()", data = x))
            }
          }
          
          list(name = NULL, data = NULL)
        }
        
        gwas_src   <- get_gwas_hits_for_ld()
        gwas_df_ld <- gwas_src$data
        
        append_log(paste0(ld_tag, "[CANDIDATES][GWAS] filtered source picked: ", gwas_src$name %||% "NULL"))
        
        gwas_candidates <- tibble::tibble()
        
        if (is.data.frame(gwas_df_ld) && nrow(gwas_df_ld)) {
          
          # hits_df() té: CHR, BP, snp, rsid, p, logp
          gwas_candidates0 <- gwas_df_ld %>%
            dplyr::transmute(
              cluster_id = NA_character_,
              chr        = suppressWarnings(as.integer(readr::parse_number(as.character(CHR)))),
              pos_ini    = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
              pos_end    = suppressWarnings(as.integer(readr::parse_number(as.character(BP)))),
              id_hit     = dplyr::coalesce(as.character(rsid), as.character(snp)),
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
        
        # NonSyn hits
        # Construïm candidates a nivell de variant única per cluster,
        # no a nivell de fila anotacional de dbNSFP.
        
        id_hit_col_candidates <- c("rsid", "RSID", "SNP", "snp", "id", "ID")
        id_hit_col <- intersect(id_hit_col_candidates, names(hits))
        id_hit_col <- if (length(id_hit_col)) id_hit_col[1] else NULL
        
        pos_hit_col_candidates <- c("pos", "POS", "position", "Position", "bp", "BP", "pos_ini", "variant_pos")
        pos_hit_col <- intersect(pos_hit_col_candidates, names(hits))
        pos_hit_col <- if (length(pos_hit_col)) pos_hit_col[1] else NULL
        
        ref_col_candidates <- c("ref", "REF", "a1", "A1", "allele1")
        alt_col_candidates <- c("alt", "ALT", "a2", "A2", "allele2")
        
        ref_col <- intersect(ref_col_candidates, names(hits))
        alt_col <- intersect(alt_col_candidates, names(hits))
        
        ref_col <- if (length(ref_col)) ref_col[1] else NULL
        alt_col <- if (length(alt_col)) alt_col[1] else NULL
        
        nonsyn_candidates <- tibble::tibble()
        
        needed_hit_cols <- c("cluster_id", "chr")
        
        if (all(needed_hit_cols %in% names(hits)) && !is.null(pos_hit_col)) {
          
          nonsyn_candidates <- hits %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              chr        = suppressWarnings(as.integer(readr::parse_number(as.character(chr)))),
              pos_ini    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_hit_col]])))),
              pos_end    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_hit_col]])))),
              rsid_raw   = if (!is.null(id_hit_col)) as.character(.data[[id_hit_col]]) else NA_character_,
              ref_raw    = if (!is.null(ref_col)) as.character(.data[[ref_col]]) else NA_character_,
              alt_raw    = if (!is.null(alt_col)) as.character(.data[[alt_col]]) else NA_character_
            ) %>%
            dplyr::mutate(
              cluster_id = trimws(cluster_id),
              rsid_raw   = trimws(rsid_raw),
              ref_raw    = trimws(ref_raw),
              alt_raw    = trimws(alt_raw)
            ) %>%
            dplyr::filter(
              !is.na(cluster_id), nzchar(cluster_id),
              is.finite(chr), is.finite(pos_ini), is.finite(pos_end)
            ) %>%
            dplyr::mutate(
              variant_key = dplyr::case_when(
                !is.na(ref_raw) & nzchar(ref_raw) & !is.na(alt_raw) & nzchar(alt_raw) ~
                  paste0(cluster_id, "||chr", chr, ":", pos_ini, ":", ref_raw, ":", alt_raw),
                TRUE ~
                  paste0(cluster_id, "||chr", chr, ":", pos_ini)
              )
            ) %>%
            dplyr::distinct(variant_key, .keep_all = TRUE) %>%
            dplyr::mutate(
              id_hit = dplyr::case_when(
                !is.na(rsid_raw) & nzchar(rsid_raw) ~ rsid_raw,
                !is.na(ref_raw) & nzchar(ref_raw) & !is.na(alt_raw) & nzchar(alt_raw) ~
                  paste0("chr", chr, ":", pos_ini, ":", ref_raw, ":", alt_raw),
                TRUE ~
                  paste0("chr", chr, ":", pos_ini)
              ),
              classe = app_hit_class
            ) %>%
            dplyr::select(cluster_id, chr, pos_ini, pos_end, id_hit, classe)
          
        } else {
          append_log(paste0(
            ld_tag, "[CANDIDATES][APP] skipped: missing columns -> ",
            paste(
              c(
                setdiff(needed_hit_cols, names(hits)),
                if (is.null(pos_hit_col)) "variant_position_col" else NULL
              ),
              collapse = ", "
            )
          ))
        }
        
        append_log(paste0(ld_tag, "[CANDIDATES][APP] nrow=", nrow(nonsyn_candidates)))
        append_log(
          paste0(
            ld_tag, "[CANDIDATES][APP] cluster counts:\n",
            paste(
              capture.output(
                print(
                  nonsyn_candidates %>%
                    dplyr::count(cluster_id, name = "n_candidates_nonsyn")
                )
              ),
              collapse = "\n"
            )
          )
        )
        
        if (nrow(nonsyn_candidates) > 0) {
          append_log(
            paste0(
              ld_tag, "[CANDIDATES][APP] head:\n",
              paste(capture.output(print(utils::head(nonsyn_candidates, 10))), collapse = "\n")
            )
          )
        }
        
        candidates_ld <- dplyr::bind_rows(gwas_candidates, nonsyn_candidates) %>%
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
        
        # ============================================================
        # E) GENE-HIT BRIDGE (NonSyn hit -> gene)
        # ============================================================
        
        gene_hit_bridge <- tibble::tibble(
          cluster_id = character(),
          source_app = character(),
          gene = character(),
          chr = integer(),
          position = integer(),
          id_hit = character(),
          hit_key = character()
        )
        
        id_hit_col_candidates <- c("rsid", "RSID", "SNP", "snp", "id", "ID")
        id_hit_col <- intersect(id_hit_col_candidates, names(hits))
        id_hit_col <- if (length(id_hit_col)) id_hit_col[1] else NULL
        
        pos_hit_col_candidates <- c("pos", "POS", "position", "Position", "bp", "BP", "pos_ini", "variant_pos")
        pos_hit_col <- intersect(pos_hit_col_candidates, names(hits))
        pos_hit_col <- if (length(pos_hit_col)) pos_hit_col[1] else NULL
        
        ref_col_candidates <- c("ref", "REF", "a1", "A1", "allele1")
        alt_col_candidates <- c("alt", "ALT", "a2", "A2", "allele2")
        
        ref_col <- intersect(ref_col_candidates, names(hits))
        alt_col <- intersect(alt_col_candidates, names(hits))
        
        ref_col <- if (length(ref_col)) ref_col[1] else NULL
        alt_col <- if (length(alt_col)) alt_col[1] else NULL
        
        if (all(c("cluster_id", "chr", "genename") %in% names(hits)) && !is.null(pos_hit_col)) {
          
          gene_hit_bridge <- hits %>%
            dplyr::transmute(
              cluster_id = trimws(as.character(cluster_id)),
              source_app = app_slug,
              gene = trimws(as.character(genename)),
              chr = suppressWarnings(as.integer(readr::parse_number(as.character(chr)))),
              position = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_hit_col]])))),
              rsid_raw = if (!is.null(id_hit_col)) trimws(as.character(.data[[id_hit_col]])) else NA_character_,
              ref_raw = if (!is.null(ref_col)) trimws(as.character(.data[[ref_col]])) else NA_character_,
              alt_raw = if (!is.null(alt_col)) trimws(as.character(.data[[alt_col]])) else NA_character_
            ) %>%
            dplyr::filter(
              !is.na(cluster_id), nzchar(cluster_id),
              !is.na(gene), nzchar(gene),
              is.finite(chr), is.finite(position)
            ) %>%
            dplyr::mutate(
              hit_key = dplyr::case_when(
                !is.na(rsid_raw) & nzchar(rsid_raw) ~ rsid_raw,
                !is.na(ref_raw) & nzchar(ref_raw) & !is.na(alt_raw) & nzchar(alt_raw) ~
                  paste0("chr", chr, ":", position, ":", ref_raw, ":", alt_raw),
                TRUE ~
                  paste0("chr", chr, ":", position)
              ),
              id_hit = hit_key
            ) %>%
            dplyr::select(cluster_id, source_app, gene, chr, position, id_hit, hit_key) %>%
            dplyr::distinct()
        }
        
        append_log(paste0(ld_tag, "[GENE_HIT_BRIDGE] nrow=", nrow(gene_hit_bridge)))
        
        if (nrow(gene_hit_bridge) > 0) {
          saveRDS(gene_hit_bridge, gene_hit_bridge_path)
          manifest$files[[paste0(app_slug, "_gene_hit_bridge")]] <- basename(gene_hit_bridge_path)
          append_log(paste0(ld_tag, "[GENE_HIT_BRIDGE] saved: ", gene_hit_bridge_path))
          append_log(paste0(ld_tag, "[GENE_HIT_BRIDGE] exists after save: ", file.exists(gene_hit_bridge_path)))
        } else {
          append_log(paste0(ld_tag, "[GENE_HIT_BRIDGE] ", basename(gene_hit_bridge_path), " NOT saved (0 rows)"))
        }
        
        # ============================================================
        # FINAL MANIFEST WRITE
        # ============================================================
        manifest$apps_present <- sort(unique(c(manifest$apps_present, app_slug)))
        manifest$last_updated <- as.character(Sys.time())
        
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
        } else if (exists("safe_log", mode = "function", inherits = TRUE)) {
          safe_log("[LD-INTEGRATOR][ERROR] ", conditionMessage(e), "\n")
        } else {
          cat("[LD-INTEGRATOR][ERROR] ", conditionMessage(e), "\n", sep = "")
        }
      })
      
      ######################## END INTEGRATOR
      incProgress(0.05, detail = "Done")
      safe_log("[END] run_nonsyn_clusters finished OK\n")
      
      if (is.null(df_final) || !nrow(df_final)) {
        showNotification("No NonSyn variants found in clusters.", type = "warning", duration = 6)
      } else {
        showNotification("NonSyn variants assigned to clusters and final outputs saved.", type = "message", duration = 4)
      }
      
    }, error = function(e) {
      
      dbnsfp_final_path_csv(NULL)
      dbnsfp_final_path_rds(NULL)
      dbnsfp_norm_path_csv(NULL)
      dbnsfp_norm_path_rds(NULL)
      
      msg <- paste0("ERROR (run_nonsyn_clusters): ", conditionMessage(e))
      
      safe_log("[ERROR] ", msg, "\n")
      showNotification(conditionMessage(e), type = "error", duration = 12)
      
    }, finally = {
      
      # if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::enable("run_nonsyn_clusters")
      
    })
    
  }, ignoreInit = TRUE)
  ######################### END INTEGRATOR #########################################  
  # -------------------------
  # Single reactive: dbNSFP dataset used everywhere
  # - OLD mode: reads uploaded CSV/RDS
  # - NEW mode: reads FINAL paths (prefer RDS)
  # -------------------------
  dbnsfp_norm_df <- reactive({
    rds <- dbnsfp_final_path_rds()
    csv <- dbnsfp_final_path_csv()
    
    if ((is.null(rds) || !file.exists(rds)) && (is.null(csv) || !file.exists(csv))) {
      return(tibble::tibble())
    }
    
    df <- NULL
    
    if (!is.null(rds) && file.exists(rds)) {
      df <- readRDS(rds)
    } else {
      df <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    validate(need(is.data.frame(df), "Invalid final dbNSFP output."))
    
    df <- sanitize_colnames(df)
    df <- harmonize_dbnsfp_loaded(df)
    
    if (!"cluster_id" %in% names(df)) {
      return(tibble::tibble())
    }
    
    df
  })
  
  detect_gene_col <- function(df) {
    nm <- names(df)
    
    # el teu cas real
    if ("genename" %in% nm) return("genename")
    
    # altres candidats habituals
    cand <- c("gene","Gene","gene_name","Gene_name","symbol","SYMBOL",
              "hgnc_symbol","HGNC","Gene.refGene","Gene_refGene","GeneRefGene")
    hit <- cand[cand %in% nm]
    if (length(hit)) return(hit[1])
    
    # fallback regex
    low <- tolower(nm)
    ok <- grepl("gene|symbol", low) & !grepl("pred$", low)
    if (any(ok)) return(nm[which(ok)[1]])
    
    NULL
  }
  # ============================================================
  # Metric registry (score-only vs paired score/rankscore metrics)
  # ============================================================
  
  nonsyn_metric_type <- list(
    # Always SCORE only (no rankscore equivalent)
    score_only = c(
      "CADD_phred",
      "GERP++_RS", "GERP_91_mammals",
      "phyloP17way_primate", "phyloP100way_vertebrate",
      "phastCons17way_primate",
      "ExAC_pLI", "gnomAD_pLI", "ClinGen_Haploinsufficiency_Score",
      "GDI", "GDI_Phred", "LoFtool_score", "Gene_indispensability_score"
    ),
    
    # Base names that exist as *_score and *_rankscore pairs
    paired = c(
      "SIFT", "Polyphen2_HDIV", "Polyphen2_HVAR",
      "MutationAssessor", "PROVEAN", "REVEL",
      "MetaRNN", "MetaSVM", "MetaLR",
      "ClinPred", "PrimateAI", "AlphaMissense",
      "MutPred", "VEST4", "MPC",
      "DANN", "Eigen_raw_coding", "Eigen_PC_raw_coding",
      "bStatistic_converted",
      "CADD_raw"
    )
  )
  
  # Return available metrics in df for a given scale ("score" or "rankscore")
  get_metrics_by_scale <- function(df, scale) {
    all_names <- names(df)
    out <- character(0)
    
    # score-only metrics are shown only in score mode
    if (identical(scale, "score")) {
      out <- c(out, intersect(nonsyn_metric_type$score_only, all_names))
    }
    
    # paired metrics: pick *_score or *_rankscore depending on mode
    for (base in nonsyn_metric_type$paired) {
      score_name <- paste0(base, "_score")
      rank_name  <- paste0(base, "_rankscore")
      
      if (identical(scale, "score") && score_name %in% all_names) {
        out <- c(out, score_name)
      }
      if (identical(scale, "rankscore") && rank_name %in% all_names) {
        out <- c(out, rank_name)
      }
    }
    
    unique(out)
  }
  
  # ============================================================
  # Metric classes (used for table filtering + radar sets)
  # ============================================================
  
  nonsyn_metric_classes <- list(
    "Default" = c(
      # score
      "SIFT_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
      "MutationAssessor_score", "PROVEAN_score", "REVEL_score",
      "GERP_91_mammals", "phyloP17way_primate", "phastCons17way_primate",
      
      # rankscore (classic functional predictors)
      "SIFT_converted_rankscore",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_rankscore",
      "Polyphen2_HVAR_rankscore",
      "MutationAssessor_rankscore",
      "PROVEAN_converted_rankscore",      # note: converted
      "MutationTaster_converted_rankscore",
      "REVEL_rankscore"
    ),
    
    "Pathogenicity" = c(
      # score
      "CADD_phred", "CADD_raw",
      "REVEL_score", "MetaRNN_score", "MetaSVM_score", "MetaLR_score",
      "ClinPred_score", "PrimateAI_score", "AlphaMissense_score",
      "MutPred_score", "VEST4_score", "MPC_score",
      
      # rankscore (meta + ML/DL)
      "CADD_raw_rankscore",
      "REVEL_rankscore",
      "MetaRNN_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "ClinPred_rankscore",
      "PrimateAI_rankscore",
      "AlphaMissense_rankscore",
      "MutPred_rankscore",
      "VEST4_rankscore",
      "MPC_rankscore",
      
      # extra rankscores (meta)
      "DANN_rankscore",
      "fathmm_XF_coding_rankscore",
      "BayesDel_addAF_rankscore",
      "BayesDel_noAF_rankscore",
      "M_CAP_rankscore",
      "DEOGEN2_rankscore",
      "LIST_S2_rankscore",
      
      # extra rankscores (DL / modern)
      "ESM1b_rankscore",
      "MVP_rankscore",
      "gMVP_rankscore",
      "VARITY_R_rankscore",
      "VARITY_ER_rankscore",
      "VARITY_R_LOO_rankscore",
      "VARITY_ER_LOO_rankscore",
      "PHACTboost_rankscore",
      "MutFormer_rankscore",
      "MutScore_rankscore"
    ),
    
    "Conservation" = c(
      # score
      "GERP++_RS", "GERP_91_mammals",
      "phyloP17way_primate", "phyloP100way_vertebrate",
      "phastCons17way_primate",
      
      # rankscore (expanded)
      "GERP++_RS_rankscore",
      "GERP_91_mammals_rankscore",
      "phyloP17way_primate_rankscore",
      "phyloP100way_vertebrate_rankscore",
      "phyloP470way_mammalian_rankscore",
      "phastCons17way_primate_rankscore",
      "phastCons100way_vertebrate_rankscore",
      "phastCons470way_mammalian_rankscore"
    ),
    
    "Functional impact" = c(
      # score
      "CADD_raw", "DANN_score", "Eigen_raw_coding", "Eigen_PC_raw_coding",
      "bStatistic_converted",
      
      # rankscore
      "DANN_rankscore",
      "Eigen_raw_coding_rankscore",
      "Eigen_PC_raw_coding_rankscore",
      "bStatistic_converted_rankscore"
    ),
    
    "LoF/Haploinsufficiency" = c(
      "ExAC_pLI", "gnomAD_pLI", "ClinGen_Haploinsufficiency_Score"
    ),
    
    "Gene damage" = c(
      "GDI", "GDI_Phred", "LoFtool_score", "Gene_indispensability_score"
    )
  )
  
  # Build grouped choices for the selector, filtering by scale AND by real data presence
  build_nonsyn_metric_choices <- function(df, scale) {
    req(is.data.frame(df), nrow(df) > 0)
    
    scale <- scale %||% "score"
    out <- list()
    
    for (grp in names(nonsyn_metric_classes)) {
      mets <- nonsyn_metric_classes[[grp]]
      
      # Keep only metrics that match the selected scale
      if (identical(scale, "rankscore")) {
        mets <- mets[grepl("rankscore$", mets, ignore.case = TRUE)]
      } else {
        mets <- mets[!grepl("rankscore$", mets, ignore.case = TRUE)]
      }
      
      # Keep only columns that exist in the data
      mets <- intersect(mets, names(df))
      
      # Keep only metrics with at least some finite values (robust parsing)
      if (length(mets)) {
        keep <- vapply(mets, function(m) metric_has_data(df, m), logical(1))
        mets <- mets[keep]
      }
      
      if (length(mets)) out[[grp]] <- mets
    }
    
    out
  }
  
  # Convenience reactive: TRUE if current selected metric ends with "rankscore"
  is_rankscore <- reactive({
    met <- input$nonsyn_metric %||% ""
    isTRUE(nzchar(met) && grepl("rankscore$", met, ignore.case = TRUE))
  })
  
  # Update the metric selector whenever dbNSFP data or scale changes
  observeEvent(
    list(input$nonsyn_scale, dbnsfp_norm_df()),
    {
      df    <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
      scale <- input$nonsyn_scale %||% "score"
      
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        updateSelectInput(
          session, "nonsyn_metric",
          choices  = list("No metrics available (no dbNSFP data yet)" = character(0)),
          selected = NULL
        )
        return()
      }
      
      choices <- build_nonsyn_metric_choices(df, scale)
      
      if (!length(choices)) {
        updateSelectInput(
          session, "nonsyn_metric",
          choices  = list("No metrics with data" = character(0)),
          selected = NULL
        )
        return()
      }
      
      updateSelectInput(
        session, "nonsyn_metric",
        choices  = choices,
        selected = unlist(choices, use.names = FALSE)[1]
      )
    },
    ignoreInit = FALSE
  )
  
  
  # -------------------------
  # Cluster intervals derived from FINAL (useful for downstream summaries)
  # -------------------------
  cluster_intervals_from_final <- reactive({
    
    fin <- dbnsfp_norm_df()
    if (!is.data.frame(fin) || nrow(fin) == 0) return(tibble::tibble())
    
    req(all(c("cluster_id", "cluster_start", "cluster_end", "chr") %in% names(fin)))
    
    fin %>%
      dplyr::filter(!is.na(cluster_id), !is.na(cluster_start), !is.na(cluster_end)) %>%
      dplyr::mutate(
        CHR   = norm_chr_generic(.data$chr),
        start = suppressWarnings(as.integer(.data$cluster_start)),
        end   = suppressWarnings(as.integer(.data$cluster_end))
      ) %>%
      dplyr::filter(!is.na(CHR), !is.na(start), !is.na(end)) %>%
      dplyr::distinct(cluster_id, CHR, start, end) %>%
      dplyr::arrange(CHR, start)
  })
  
  # ---------------------------
  # Download FINAL normalized dbNSFP (with cluster_id)
  # ---------------------------
  
  output$dl_dbnsfp_csv <- downloadHandler(
    filename = function() {
      tg <- cluster_build_tag() %||% make_mode_thr_tag(input$cluster_method, input$pthr, input$min_logp)
      paste0(
        "dbnsfp_normalized_",
        tg$mode_tag, "_thr", tg$thr_txt,
        "_with_clusters_",
        format(Sys.Date(), "%Y%m%d"),
        ".csv"
      )
    },
    content = function(file) {
      
      # Prefer FINAL CSV; fallback to FINAL RDS -> write CSV on the fly
      src_csv <- dbnsfp_final_path_csv()
      src_rds <- dbnsfp_final_path_rds()
      
      if (!is.null(src_csv) && file.exists(src_csv) && file.info(src_csv)$size > 0) {
        file.copy(src_csv, file, overwrite = TRUE)
        return(invisible())
      }
      
      validate(need(
        !is.null(src_rds) && file.exists(src_rds) && file.info(src_rds)$size > 0,
        "Final normalized output is not available yet. Run Step 3 (annotation + normalization) first."
      ))
      
      df <- readRDS(src_rds)
      validate(need("cluster_id" %in% names(df), "Final output has no cluster_id column."))
      utils::write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$dl_dbnsfp_rds <- downloadHandler(
    filename = function() {
      tg <- cluster_build_tag() %||% make_mode_thr_tag(input$cluster_method, input$pthr, input$min_logp)
      paste0(
        "dbnsfp_normalized_",
        tg$mode_tag, "_thr", tg$thr_txt,
        "_with_clusters_",
        format(Sys.Date(), "%Y%m%d"),
        ".rds"
      )
    },
    content = function(file) {
      
      # Prefer FINAL RDS; fallback to FINAL CSV -> saveRDS on the fly
      src_rds <- dbnsfp_final_path_rds()
      src_csv <- dbnsfp_final_path_csv()
      
      if (!is.null(src_rds) && file.exists(src_rds) && file.info(src_rds)$size > 0) {
        file.copy(src_rds, file, overwrite = TRUE)
        return(invisible())
      }
      
      validate(need(
        !is.null(src_csv) && file.exists(src_csv) && file.info(src_csv)$size > 0,
        "Final normalized output is not available yet. Run Step 3 (annotation + normalization) first."
      ))
      
      df <- utils::read.csv(src_csv, check.names = FALSE, stringsAsFactors = FALSE)
      validate(need("cluster_id" %in% names(df), "Final output has no cluster_id column."))
      saveRDS(df, file)
    }
  )
  
  output$dl_candidates_zip <- downloadHandler(
    
    filename = function() {
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      tag_txt <- tryCatch({
        tg <- cluster_build_tag() %||% make_mode_thr_tag(input$cluster_method %||% "window", input$pthr, input$min_logp)
        paste0(tg$mode_tag, "_thr", tg$thr_txt)
      }, error = function(e) "run")
      
      paste0("nonsyn_candidates_", tag_txt, "_", stamp, ".zip")
    },
    
    content = function(file) {
      
      # -----------------------------
      # Helper: pick first existing column name
      # -----------------------------
      pick_first_col <- function(df, candidates) {
        candidates <- as.character(candidates)
        for (nm in candidates) {
          if (nm %in% names(df)) return(nm)
        }
        NULL
      }
      
      # -----------------------------
      # Prepare clusters (NonSyn naming)
      # -----------------------------
      cl <- cluster_dt_view2()
      validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available to export. Build clusters first."))
      
      cl <- as.data.frame(cl)
      
      # NonSyn cluster_dt_view2() should already have: cluster_id, chr, start, end
      validate(need(all(c("cluster_id", "chr", "start", "end") %in% names(cl)),
                    "cluster_dt_view2() must contain: cluster_id, chr, start, end."))
      
      cl$chr   <- chr_map_plink19(norm_chr_generic(cl$chr))
      cl$start <- suppressWarnings(as.integer(cl$start))
      cl$end   <- suppressWarnings(as.integer(cl$end))
      cl$cluster_id <- as.character(cl$cluster_id)
      
      cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start & nzchar(cl$cluster_id), , drop = FALSE]
      validate(need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
      
      # -----------------------------
      # Prepare NonSyn variants (dbNSFP) + ensure cluster assignment
      # -----------------------------
      ns <- dbnsfp_norm_df()
      validate(need(is.data.frame(ns) && nrow(ns) > 0,
                    "No NonSyn (dbNSFP) data available to export. Run 'Extract nonsyn variants in all clusters' first."))
      ns <- as.data.frame(ns)
      
      # Helper: pick first existing column
      pick_first_col <- function(df, candidates) {
        for (nm in candidates) if (nm %in% names(df)) return(nm)
        NULL
      }
      
      # Try to find cluster columns
      ns_clcol <- pick_first_col(ns, c("cluster_id", "cluster_chr_n", "cluster", "label"))
      
      # If cluster info missing (or empty), assign clusters now
      need_assign <- is.null(ns_clcol) || all(is.na(ns[[ns_clcol]]) | !nzchar(as.character(ns[[ns_clcol]])))
      if (need_assign) {
        validate(need(exists("add_clusters_to_dbnsfp", mode = "function"),
                      "add_clusters_to_dbnsfp() not found (needed to assign clusters for export)."))
        ns <- add_clusters_to_dbnsfp(ns, cluster_dt_view2())
        ns_clcol <- "cluster_id"
      }
      
      # If in_cluster exists, filter correctly (VECTORIZED)
      if ("in_cluster" %in% names(ns)) {
        keep <- !is.na(ns$in_cluster) & as.logical(ns$in_cluster)
        ns <- ns[keep, , drop = FALSE]
      }
      
      # If still no cluster_id-like values, stop
      validate(need(nrow(ns) > 0, "dbNSFP data has no rows assigned to clusters (nothing to export)."))
      validate(need(!is.null(ns_clcol) && ns_clcol %in% names(ns), "Could not find cluster column after assignment."))
      
      # Positions in dbNSFP
      ns_chr_col <- pick_first_col(ns, c("chr", "#chr", "CHR", "chrom", "chromosome"))
      ns_pos_col <- pick_first_col(ns, c("BP", "POS", "pos(1-based)", "pos", "position"))
      validate(need(!is.null(ns_chr_col), "dbNSFP data has no chromosome column (chr/#chr/CHR/...)."))
      validate(need(!is.null(ns_pos_col), "dbNSFP data has no position column (BP/POS/pos(1-based)/...)."))
      
      # ID for NonSyn variant
      ns_id_col <- pick_first_col(ns, c("rs_dbSNP", "rsid", "RSID", "variant_id", "variant", "ID", "id"))
      if (is.null(ns_id_col)) {
        ns_id_col <- ".id_tmp"
        ns[[ns_id_col]] <- paste0("chr", ns[[ns_chr_col]], ":", ns[[ns_pos_col]])
      }
      
      # Build b = point variants
      b <- data.frame(
        chr        = chr_map_plink19(norm_chr_generic(ns[[ns_chr_col]])),
        bin_start  = suppressWarnings(as.integer(ns[[ns_pos_col]])),
        bin_end    = suppressWarnings(as.integer(ns[[ns_pos_col]])),
        id_hit     = as.character(ns[[ns_id_col]]),
        cluster_id = as.character(ns[[ns_clcol]]),
        stringsAsFactors = FALSE
      )
      
      b <- b[is.finite(b$chr) & is.finite(b$bin_start) & is.finite(b$bin_end) &
               b$bin_end >= b$bin_start & nzchar(b$cluster_id), , drop = FALSE]
      
      validate(need(nrow(b) > 0, "All NonSyn variants became invalid after coercion (chr/pos/cluster_id)."))
      
      b <- data.frame(
        chr       = chr_map_plink19(norm_chr_generic(ns[[ns_chr_col]])),
        bin_start = suppressWarnings(as.integer(ns[[ns_pos_col]])),
        bin_end   = suppressWarnings(as.integer(ns[[ns_pos_col]])),
        id_hit    = as.character(ns[[ns_id_col]]),
        cluster_id = as.character(ns[[ns_clcol]]),
        stringsAsFactors = FALSE
      )
      
      b <- b[is.finite(b$chr) & is.finite(b$bin_start) & is.finite(b$bin_end) & b$bin_end >= b$bin_start & nzchar(b$cluster_id), , drop = FALSE]
      validate(need(nrow(b) > 0, "All NonSyn variants became invalid after coercion (chr/pos)."))
      
      # -----------------------------
      # Keep clusters that have NonSyn variants
      # -----------------------------
      keep_ids <- unique(as.character(b$cluster_id))
      cl_sig <- cl[as.character(cl$cluster_id) %in% keep_ids, , drop = FALSE]
      validate(need(nrow(cl_sig) > 0, "No clusters matched NonSyn variants (nothing to export)."))
      
      # -----------------------------
      # Build cluster CSV (export names)
      # -----------------------------
      cluster_csv <- cl_sig[, c("cluster_id", "chr", "start", "end")]
      names(cluster_csv) <- c("cluster_id", "chr", "cluster_start", "cluster_end")
      
      # -----------------------------
      # Build candidate CSV (GWAS + NonSyn)
      # -----------------------------
      h <- tryCatch(hits_df(), error = function(e) NULL)
      
      if (is.null(h) || !nrow(h)) {
        gwas_cand <- data.frame(
          chr = integer(), pos_ini = integer(), pos_end = integer(),
          id_hit = character(), classe = character(),
          stringsAsFactors = FALSE
        )
      } else {
        validate(need(all(c("CHR", "BP") %in% names(h)), "hits_df() must contain CHR and BP columns."))
        
        rsid <- if ("rsid" %in% names(h)) as.character(h$rsid) else if ("snp" %in% names(h)) as.character(h$snp) else NA_character_
        rsid[!nzchar(rsid)] <- paste0("chr", h$CHR, ":", h$BP)
        
        gwas_cand <- data.frame(
          chr     = chr_map_plink19(norm_chr_generic(h$CHR)),
          pos_ini = as.integer(h$BP),
          pos_end = as.integer(h$BP),
          id_hit  = rsid,
          classe  = "GWAS",
          stringsAsFactors = FALSE
        )
      }
      
      nonsyn_cand <- data.frame(
        chr     = as.integer(b$chr),
        pos_ini = as.integer(b$bin_start),
        pos_end = as.integer(b$bin_end),
        id_hit  = as.character(b$id_hit),
        classe  = "NonSyn_variant",
        stringsAsFactors = FALSE
      )
      
      candidate_csv <- rbind(gwas_cand, nonsyn_cand)
      
      # -----------------------------
      # Write to temp folder and zip
      # -----------------------------
      tmpdir <- tempfile("nonsyn_export_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      f1 <- file.path(tmpdir, "cluster_nonsyn.csv")
      f2 <- file.path(tmpdir, "candidate_nonsyn.csv")
      
      utils::write.csv(cluster_csv, f1, row.names = FALSE, quote = FALSE)
      utils::write.csv(candidate_csv, f2, row.names = FALSE, quote = FALSE)
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      setwd(tmpdir)
      
      utils::zip(zipfile = file, files = c("cluster_nonsyn.csv", "candidate_nonsyn.csv"))
    }
  )
  
  # ========================================================================
  # ============== Visualizing RESULTS and DATA ============================
  # ========================================================================
  
  # ------------------------------------------------------------------------
  # 0) Optional debug helper (safe: does nothing unless option is TRUE)
  # ------------------------------------------------------------------------
  .debug_on <- function() isTRUE(getOption("nonsyn.debug", FALSE))
  .debug_msg <- function(...) if (.debug_on()) message(...)
  
  # ========================================================================
  # 1) Helpers: robust numeric parsing for dbNSFP metrics
  # ========================================================================
  
  parse_metric_numeric <- function(x) {
    # numeric -> as is
    if (is.numeric(x)) return(x)
    
    # logical -> optional (TRUE/FALSE as 1/0)
    if (is.logical(x)) return(as.numeric(x))
    
    x <- as.character(x)
    x[x %in% c("", ".", "NA", "NaN")] <- NA_character_
    
    # If multi-valued ("0.12;0.34"), take first numeric-looking token
    x <- vapply(strsplit(x, ";", fixed = TRUE), FUN.VALUE = character(1), function(parts) {
      parts <- trimws(parts)
      parts <- parts[parts != "" & parts != "." & parts != "NA"]
      if (!length(parts)) return(NA_character_)
      
      i <- which(grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", parts))[1]
      if (is.na(i)) return(NA_character_)
      parts[i]
    })
    
    suppressWarnings(as.numeric(x))
  }
  
  metric_candidates <- function(nms, mode) {
    if (identical(mode, "rankscore")) {
      grep("rankscore$", nms, value = TRUE, ignore.case = TRUE)
    } else {
      # score mode
      c(
        grep("_score$", nms, value = TRUE, ignore.case = TRUE),
        intersect(c("CADD_raw", "CADD_phred"), nms)
      ) %>% unique()
    }
  }
  
  available_metrics <- function(df, mode = c("rankscore", "score")) {
    mode <- match.arg(mode)
    
    if (mode == "rankscore") {
      cand <- grep("rankscore$", names(df), value = TRUE, ignore.case = TRUE)
    } else {
      cand <- grep("(_score$|_phred$|_raw$)", names(df), value = TRUE, ignore.case = TRUE)
      cand <- setdiff(cand, grep("rankscore$", cand, value = TRUE, ignore.case = TRUE))
    }
    
    keep <- cand[vapply(cand, function(m) {
      x <- suppressWarnings(as.numeric(df[[m]]))
      any(is.finite(x))
    }, logical(1))]
    
    keep
  }
  
  # ========================================================================
  # 2) Metric status + warning UI (robust even if plot works)
  # ========================================================================
  
  nonsyn_metric_status <- reactive({
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    metric <- input$nonsyn_metric
    
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    if (is.null(metric) || !nzchar(metric) || !(metric %in% names(df))) return(NULL)
    
    x <- parse_metric_numeric(df[[metric]])
    n_ok <- sum(is.finite(x))
    
    list(metric = metric, n_ok = n_ok, n_total = length(x))
  })
  
  output$nonsyn_metric_warning <- renderUI({
    st <- nonsyn_metric_status()
    if (is.null(st)) return(NULL)
    
    if (st$n_ok == 0) {
      div(
        class = "alert alert-warning",
        style = "padding: 8px 12px; margin: 0;",
        shiny::icon("exclamation-triangle"),
        HTML(paste0(
          "&nbsp;<b>No numeric values found</b> for the selected metric: <code>",
          st$metric,
          "</code>. (All values are NA / non-numeric). Please choose another metric, ideally a <code>*rankscore</code> column."
        ))
      )
    } else {
      NULL
    }
  })
  
  # ========================================================================
  # 3) Build NonSyn dataframe for genome-wide plots (single source of truth)
  # ========================================================================
  
  nonsyn_manhattan_df <- reactive({
    
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(tibble::tibble())
    
    metric <- input$nonsyn_metric
    if (is.null(metric) || !nzchar(metric) || !(metric %in% names(df))) return(tibble::tibble())
    
    # Helper: return column if exists, otherwise NA vector
    col_or_na <- function(.df, col) {
      if (col %in% names(.df)) .df[[col]] else rep(NA, nrow(.df))
    }
    
    out <- df %>%
      dplyr::transmute(
        CHR  = norm_chr_generic(col_or_na(df, "chr")),
        POS  = suppressWarnings(as.integer(col_or_na(df, "BP"))),
        GENE = as.character(col_or_na(df, "genename")),
        rsid = as.character(col_or_na(df, "SNP")),
        value = parse_metric_numeric(.data[[metric]]),
        cluster_id = as.character(col_or_na(df, "cluster_id"))
      ) %>%
      dplyr::filter(!is.na(CHR), !is.na(POS), is.finite(value))
    
    
    if (!nrow(out)) return(tibble::tibble())
    
    ref <- .ref_hg38 %>% dplyr::mutate(chr = as.character(chr))
    
    out <- out %>%
      dplyr::left_join(ref %>% dplyr::select(chr, chr_cum), by = c("CHR" = "chr")) %>%
      dplyr::filter(!is.na(chr_cum)) %>%
      dplyr::mutate(BPcum = POS + chr_cum)
    
    if (!nrow(out)) return(tibble::tibble())
    
    # Ensure a usable ID
    out$rsid[is.na(out$rsid) | out$rsid %in% c("", ".", "NA")] <-
      paste0("chr", out$CHR, ":", out$POS)
    
    if (.debug_on()) {
      cat("[DEBUG] nonsyn_manhattan_df (head) =\n")
      print(head(out))
    }
    
    out
  })
  
  # Optional: debug availability of cluster_id in final output
  observe({
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    if (.debug_on() && is.data.frame(df)) {
      cat("[DEBUG] dbnsfp_norm_df cols has cluster_id =", "cluster_id" %in% names(df), "\n")
      if ("cluster_id" %in% names(df)) print(head(df$cluster_id))
    }
  })
  
  # ========================================================================
  # 4) Window selection from Plotly (shared by Manhattan + UCSC)
  # ========================================================================
  
  window_selected <- reactiveVal(NULL)
  ucsc_region     <- reactiveVal(NULL)
  
  observeEvent(plotly::event_data("plotly_relayout"), {
    ev <- plotly::event_data("plotly_relayout")
    req(ev)
    
    x0 <- ev[["xaxis.range[0]"]]
    x1 <- ev[["xaxis.range[1]"]]
    
    if (!is.null(x0) && !is.null(x1)) {
      window_selected(list(xmin = x0, xmax = x1))
    }
  })
  
  output$debug_window <- renderPrint({
    window_selected()
  })
  
  # Small debug log (safe)
  observe({
    win <- window_selected()
    ns  <- nonsyn_manhattan_df()
    
    .debug_msg(
      "Window = ",
      ifelse(is.null(win), "NULL", paste(round(win$xmin), round(win$xmax))),
      " | NonSyn total = ", nrow(ns),
      " | NonSyn in window = ",
      ifelse(is.null(win), NA,
             sum(ns$BPcum >= win$xmin & ns$BPcum <= win$xmax))
    )
  })
  
  get_window_hits_gwas <- reactive({
    win <- window_selected(); req(win)
    dfp <- dfp_manhattan()
    req(nrow(dfp) > 0)
    dfp %>% dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  get_window_hits_nonsyn <- reactive({
    win <- window_selected(); req(win)
    
    ns <- nonsyn_manhattan_df()
    if (!is.data.frame(ns) || nrow(ns) == 0) return(NULL)
    
    ns %>% dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  # ========================================================================
  # 5) Manhattan combo (GWAS on top, NonSyn metric on bottom)
  # ========================================================================
  
  output$manhattan_combo <- plotly::renderPlotly({
    
    src_combo <- "manhattan_combo"
    
    # -----------------------------
    # TOP: GWAS Manhattan (EWAS-style)
    # -----------------------------
    dfp <- tryCatch(dfp_manhattan(), error = function(e) NULL)
    if (is.null(dfp) || !is.data.frame(dfp) || !nrow(dfp)) {
      return(plotly_message("⚠️ GWAS table missing or incomplete."))
    }
    
    ref <- .ref_hg38
    ax  <- axis_df()
    
    axis_breaks <- ax$center
    axis_labels <- paste0("chr", ax$chrN)
    GENOME_END  <- max(ref$chr_cum + ref$len, na.rm = TRUE)
    
    thr_y <- if ((input$cluster_method %||% "window") == "window") (input$pthr %||% 5) else (input$min_logp %||% 6)
    
    dfp <- dfp %>%
      dplyr::arrange(CHR, BP) %>%
      dplyr::mutate(
        rs_show = ifelse(!is.na(snp) & nzchar(snp), snp, paste0("chr", CHR, ":", BP)),
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
    
    p1 <- ggplot(dfp, aes(x = BPcum, y = logp, text = tooltip)) +
      geom_point(aes(color = col), size = 1) +
      geom_hline(yintercept = thr_y, linetype = "dashed") +
      scale_color_identity(guide = "none") +
      scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      labs(x = NULL, y = "-log10(P)") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    p1_pl <- plotly::ggplotly(p1, tooltip = "text", source = src_combo)
    
    # -----------------------------------
    # BOTTOM: NonSyn metric (EWAS-style)
    # -----------------------------------
    metric <- input$nonsyn_metric
    ns2    <- tryCatch(nonsyn_manhattan_df(), error = function(e) NULL)
    
    p2_base <- ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
      geom_blank() +
      scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = "Genome", y = paste0(metric, " [NonSyn]")) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    p2_pl <- plotly::ggplotly(p2_base, source = src_combo)
    y_max <- 1
    
    if (is.data.frame(ns2) && nrow(ns2)) {
      
      # highlight selected cluster (keep your logic)
      cl_sel <- tryCatch(selected_cluster(), error = function(e) NULL)
      ns2$highlight <- "other"
      
      if (!is.null(cl_sel) && is.data.frame(cl_sel) && nrow(cl_sel) == 1 &&
          all(c("chr", "start", "end") %in% names(cl_sel))) {
        
        chr_sel_num <- suppressWarnings(as.integer(gsub("^chr", "", as.character(cl_sel$chr[1]), ignore.case = TRUE)))
        ns2_chr_num <- suppressWarnings(as.integer(gsub("^chr", "", as.character(ns2$CHR), ignore.case = TRUE)))
        
        ns2$highlight <- ifelse(
          is.finite(ns2_chr_num) & is.finite(chr_sel_num) &
            (ns2_chr_num == chr_sel_num) &
            (ns2$POS >= cl_sel$start[1]) & (ns2$POS <= cl_sel$end[1]),
          "cluster", "other"
        )
      }
      
      ns2$tooltip <- paste0(
        "<b>NonSyn</b>",
        "<br>rs: ", ns2$rsid,
        "<br>CHR: ", ns2$CHR,
        "<br>POS: ", ns2$POS,
        "<br>", metric, ": ", signif(ns2$value, 4),
        if ("cluster_id" %in% names(ns2)) paste0("<br>cluster_id: ", ns2$cluster_id) else ""
      )
      
      y_max <- max(ns2$value[is.finite(ns2$value)], na.rm = TRUE)
      if (!is.finite(y_max) || y_max <= 0) y_max <- 1
      
      ylab <- if (grepl("rankscore$", metric, ignore.case = TRUE)) paste0(metric, " (0–1)") else metric
      
      g2 <- ggplot(ns2, aes(x = BPcum, y = value, text = tooltip)) +
        geom_point(aes(color = highlight), size = 1.2, alpha = 0.8) +
        scale_color_manual(values = c("other" = "grey60", "cluster" = "red"), guide = "none") +
        scale_x_continuous(
          limits = c(0, GENOME_END),
          breaks = axis_breaks,
          labels = axis_labels,
          expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(0, y_max * 1.5), expand = c(0, 0)) +
        labs(x = "Genome", y = ylab) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      if (identical(metric, "CADD_phred")) {
        thr_cadd <- nonsyn_metric_thresholds[["CADD_phred"]]$value %||% NA_real_
        if (is.finite(thr_cadd)) {
          g2 <- g2 + geom_hline(yintercept = thr_cadd, linetype = "dashed", linewidth = 0.6, colour = "#1A4E8A")
        }
      }
      
      p2_pl <- plotly::ggplotly(g2, tooltip = "text", source = src_combo)
    }
    
    ns_ok <- is.data.frame(ns2) && nrow(ns2) > 0
    
    # ----------------------------------------------------------
    # Cluster segment band (FIXED: ref chr column robust)
    # ----------------------------------------------------------
    cl_raw <- NULL
    if (exists("clusters_src", mode = "function")) {
      cl_raw <- tryCatch(clusters_src(), error = function(e) NULL)
    }
    if (is.null(cl_raw) && exists("clusters_rv", mode = "function")) {
      cl_raw <- tryCatch(clusters_rv(), error = function(e) NULL)
    }
    if (is.null(cl_raw) && exists("rv", inherits = TRUE)) {
      cl_raw <- tryCatch(rv$clusters, error = function(e) NULL)
    }
    
    # Build robust ref mapping: chr_num + chr_cum
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
    
    if (is.data.frame(cl_raw) && nrow(cl_raw) && is.data.frame(ref_map) && nrow(ref_map)) {
      
      cl_tmp <- as.data.frame(cl_raw)
      if (!"chr" %in% names(cl_tmp) && "CHR" %in% names(cl_tmp)) cl_tmp$chr <- cl_tmp$CHR
      if (!"start" %in% names(cl_tmp) && "start_bp" %in% names(cl_tmp)) cl_tmp$start <- cl_tmp$start_bp
      if (!"end" %in% names(cl_tmp) && "end_bp" %in% names(cl_tmp)) cl_tmp$end <- cl_tmp$end_bp
      
      # --- ensure cluster_id exists (critical for overlay) ---
      if (!"cluster_id" %in% names(cl_tmp)) {
        if ("cluster_chr_n" %in% names(cl_tmp)) {
          cl_tmp$cluster_id <- as.character(cl_tmp$cluster_chr_n)
        } else if ("cluster" %in% names(cl_tmp)) {
          cl_tmp$cluster_id <- paste0("cluster_", as.character(cl_tmp$cluster))
        } else {
          # last resort: deterministic id
          cl_tmp$cluster_id <- paste0("cluster_", seq_len(nrow(cl_tmp)))
        }
      }
      
      if (all(c("cluster_id", "chr", "start", "end") %in% names(cl_tmp))) {
        
        clseg <- cl_tmp %>%
          dplyr::transmute(
            cluster_id = as.character(cluster_id),
            chr_num    = suppressWarnings(as.integer(gsub("^chr", "", as.character(chr), ignore.case = TRUE))),
            start_i    = suppressWarnings(as.numeric(start)),
            end_i      = suppressWarnings(as.numeric(end))
          ) %>%
          dplyr::filter(
            !is.na(cluster_id), nzchar(cluster_id),
            is.finite(chr_num), is.finite(start_i), is.finite(end_i),
            end_i >= start_i
          ) %>%
          dplyr::distinct(cluster_id, chr_num, start_i, end_i) %>%
          dplyr::left_join(ref_map, by = "chr_num") %>%
          dplyr::filter(is.finite(chr_cum)) %>%
          dplyr::transmute(
            cluster_id = cluster_id,
            x0   = pmax(0, start_i + chr_cum),
            x1   = pmin(GENOME_END, end_i + chr_cum),
            xmid = (x0 + x1) / 2,
            text = paste0(
              "Cluster: ", cluster_id,
              "<br>chr", chr_num, ":", format(start_i, scientific = FALSE),
              "-", format(end_i, scientific = FALSE)
            )
          ) %>%
          dplyr::filter(is.finite(x0), is.finite(x1), x1 >= x0) %>%
          dplyr::arrange(x0)
        
        if (nrow(clseg) > 0) {
          
          if (isTRUE(ns_ok)) {
            # Hi ha NonSyn punts: y-axis arriba fins y_max*1.5
            y_seg  <- y_max * 1.10   # (més segur que 1.32)
            y_tick <- y_max * 0.02
            y_txt  <- y_max * 1.18
          } else {
            # NO hi ha punts: p2_base té y = [0..1] → segments dins del rang
            y_seg  <- 0.90
            y_tick <- 0.02
            y_txt  <- 0.95
          }
          
          
          clseg$y_seg  <- y_seg
          clseg$y0tick <- y_seg - y_tick
          clseg$y1tick <- y_seg + y_tick
          clseg$y_txt  <- y_txt
          
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
            ) %>%
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
        }
      }
    }
    
    # -----------------------------
    # Combine
    # -----------------------------
    out <- plotly::subplot(
      p1_pl, p2_pl,
      nrows = 2, shareX = TRUE,
      heights = c(0.55, 0.45),
      titleY = TRUE
    )
    
    out$x$source <- src_combo
    out <- plotly::event_register(out, "plotly_relayout")
    
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
  
  
  
  # ========================================================================
  # 6) NonSyn table (filtered by chromosome) + external links
  # ========================================================================
  
  output$nonsyn_table <- DT::renderDT({
    
    df <- dbnsfp_norm_df()
    
    req(is.data.frame(df), nrow(df) > 0)
    
    # 1) Filter by chromosome (UI input: hm_chr)
    chr_col <- session$userData$chr_col %||% "chr"
    
    if (!is.null(input$hm_chr) && input$hm_chr != "all" && chr_col %in% names(df)) {
      
      norm_chr <- function(x) {
        x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
        x[x == "23"] <- "X"
        x[x == "24"] <- "Y"
        x[x %in% c("M","MT","m","mT","mtdna","MTDNA")] <- "MT"
        toupper(x)
      }
      
      df[[chr_col]] <- norm_chr(df[[chr_col]])
      df <- df[df[[chr_col]] == input$hm_chr, , drop = FALSE]
    }
    
    validate(need(nrow(df) > 0, "No variants match current filter."))
    
    # 2) Base columns + class columns
    base_cols <- c("SNP", "chr", "BP", "ref", "alt", "genename")
    
    classes <- nonsyn_metric_classes
    choice  <- input$clinical_class
    req(choice)
    
    class_pool <- classes[[choice]] %||% character(0)
    
    # Optional filter: score / rankscore
    kind <- input$metric_kind %||% "all"
    if (identical(kind, "score")) {
      class_pool <- class_pool[!grepl("rankscore$", class_pool, ignore.case = TRUE)]
    } else if (identical(kind, "rankscore")) {
      class_pool <- class_pool[grepl("rankscore$", class_pool, ignore.case = TRUE)]
    }
    
    class_cols <- intersect(class_pool, names(df))
    
    validate(
      need(length(class_cols) > 0,
           paste0("No columns available for class: ", choice,
                  " (metric type: ", kind, ")."))
    )
    
    cols <- unique(c(base_cols, class_cols))
    df <- df[, cols, drop = FALSE]
    
    # 3) Link columns
    if ("genename" %in% names(df)) {
      df$genename <- ifelse(
        is.na(df$genename) | df$genename == "",
        "",
        sprintf(
          "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank'>%s</a>",
          df$genename, df$genename
        )
      )
    }
    
    if ("SNP" %in% names(df)) {
      df$SNP <- ifelse(
        grepl("^rs", df$SNP, ignore.case = TRUE),
        sprintf(
          "<a href='https://www.ncbi.nlm.nih.gov/snp/%s' target='_blank'>%s</a>",
          df$SNP, df$SNP
        ),
        df$SNP
      )
    }
    
    if ("clinvar_id" %in% names(df)) {
      df$clinvar_id <- ifelse(
        !is.na(df$clinvar_id) & df$clinvar_id != "",
        sprintf(
          "<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/%s' target='_blank'>%s</a>",
          df$clinvar_id, df$clinvar_id
        ),
        ""
      )
    }
    
    if (all(c("chr", "BP") %in% names(df))) {
      df$UCSC <- sprintf(
        "<a href='https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr%s:%s-%s' target='_blank'>view</a>",
        df$chr, df$BP, df$BP
      )
    }
    
    # 4) Render
    DT::datatable(
      df,
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
  }, server = FALSE)
  
  # ========================================================================
  # 7) Summary tables + simple plots (gene bar + metric histogram)
  # ========================================================================
  
  # --- Text summary (optional) -> use verbatimTextOutput("nonsyn_interpret_text") in UI
  output$nonsyn_interpret_text <- renderText({
    
    dt <- nonsyn_manhattan_df()
    if (!is.data.frame(dt) || !nrow(dt)) return("No NonSyn variants available.")
    
    metric <- input$nonsyn_metric
    
    glob <- dt %>%
      dplyr::summarise(
        n       = n(),
        mean    = mean(value, na.rm = TRUE),
        median  = median(value, na.rm = TRUE),
        max     = max(value, na.rm = TRUE),
        sd      = sd(value, na.rm = TRUE),
        p90     = quantile(value, 0.90, na.rm = TRUE),
        p95     = quantile(value, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    
    txt <- paste0(
      "Clinical Summary stats\n",
      "====================\n",
      "Metric: ", metric, "\n\n",
      "GLOBAL SUMMARY\n",
      "--------------\n",
      sprintf(
        "n=%d | mean=%.3f | median=%.3f | max=%.3f | sd=%.3f | P90=%.3f | P95=%.3f\n\n",
        glob$n, glob$mean, glob$median, glob$max, glob$sd, glob$p90, glob$p95
      )
    )
    
    chr_sum <- dt %>%
      dplyr::group_by(CHR) %>%
      dplyr::summarise(
        n      = n(),
        mean   = mean(value, na.rm = TRUE),
        p90    = quantile(value, 0.90, na.rm = TRUE),
        p95    = quantile(value, 0.95, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        z_chr  = ifelse(glob$sd > 0, (mean - glob$mean) / glob$sd, NA_real_),
        fc_chr = mean / glob$mean
      ) %>%
      dplyr::arrange(CHR)
    
    txt <- paste0(
      txt,
      "BY CHROMOSOME\n",
      "-------------\n",
      paste(sprintf(
        "chr%s | n=%d | mean=%.3f | Z_chr=%.2f | FC_chr=%.2f | P90=%.3f | P95=%.3f",
        chr_sum$CHR, chr_sum$n, chr_sum$mean, chr_sum$z_chr, chr_sum$fc_chr, chr_sum$p90, chr_sum$p95
      ), collapse = "\n"),
      "\n\n"
    )
    
    cl <- selected_cluster()
    if (!is.null(cl) && nrow(cl) == 1) {
      chr_sel   <- chr_label_plink(as.integer(cl$chr[1]))
      start_sel <- cl$start[1]
      end_sel   <- cl$end[1]
      
      dt_cl <- dt %>% dplyr::filter(CHR == chr_sel, POS >= start_sel, POS <= end_sel)
      
      if (nrow(dt_cl)) {
        cl_sum <- dt_cl %>%
          dplyr::summarise(
            n      = n(),
            mean   = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            max    = max(value, na.rm = TRUE),
            p90    = quantile(value, 0.90, na.rm = TRUE),
            p95    = quantile(value, 0.95, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            z_cl  = ifelse(glob$sd > 0, (mean - glob$mean) / glob$sd, NA_real_),
            fc_cl = mean / glob$mean
          )
        
        txt <- paste0(
          txt,
          "SELECTED CLUSTER\n",
          "----------------\n",
          sprintf("Region: chr%s:%d-%d\n", chr_sel, start_sel, end_sel),
          sprintf(
            "n=%d | mean=%.3f | median=%.3f | max=%.3f | Z_cl=%.2f | FC_cl=%.2f | P90=%.3f | P95=%.3f\n\n",
            cl_sum$n, cl_sum$mean, cl_sum$median, cl_sum$max, cl_sum$z_cl, cl_sum$fc_cl, cl_sum$p90, cl_sum$p95
          )
        )
      } else {
        txt <- paste0(txt, "SELECTED CLUSTER\n----------------\nNo NonSyn variants in selected cluster.\n\n")
      }
    } else {
      txt <- paste0(txt, "SELECTED CLUSTER\n----------------\nNo cluster selected.\n\n")
    }
    
    paste0(
      txt,
      "NOTES\n-----\n",
      "Z_chr / Z_cl : (mean_group − mean_global) / sd_global\n",
      "FC_chr / FC_cl: mean_group / mean_global\n"
    )
  })
  
  # --- Main interpretation TABLE (this keeps the original output id)
  nonsyn_interpret_table <- reactive({
    dt <- nonsyn_manhattan_df() %>% dplyr::ungroup() %>% tibble::as_tibble()
    req(nrow(dt) > 0)
    
    metric <- input$nonsyn_metric
    
    glob <- dt %>%
      dplyr::summarise(
        scope  = "Global",
        id     = "All",
        metric = metric,
        n      = n(),
        mean   = mean(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        max    = max(value, na.rm = TRUE),
        sd     = sd(value, na.rm = TRUE),
        p90    = as.numeric(quantile(value, 0.90, na.rm = TRUE, names = FALSE)),
        p95    = as.numeric(quantile(value, 0.95, na.rm = TRUE, names = FALSE))
      )
    
    chr_tab <- dt %>%
      dplyr::group_by(CHR) %>%
      dplyr::summarise(
        scope  = "Chromosome",
        id     = paste0("chr", unique(CHR)[1]),
        metric = metric,
        n      = n(),
        mean   = mean(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        max    = max(value, na.rm = TRUE),
        sd     = sd(value, na.rm = TRUE),
        p90    = as.numeric(quantile(value, 0.90, na.rm = TRUE, names = FALSE)),
        p95    = as.numeric(quantile(value, 0.95, na.rm = TRUE, names = FALSE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        zscore  = (mean - glob$mean) / glob$sd,
        foldchg = mean / glob$mean
      )
    
    gene_tab <- dt %>%
      dplyr::group_by(GENE) %>%
      dplyr::summarise(
        scope  = "Gene",
        id     = unique(GENE)[1],
        metric = metric,
        n      = n(),
        mean   = mean(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        max    = max(value, na.rm = TRUE),
        sd     = sd(value, na.rm = TRUE),
        p90    = as.numeric(quantile(value, 0.90, na.rm = TRUE, names = FALSE)),
        p95    = as.numeric(quantile(value, 0.95, na.rm = TRUE, names = FALSE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        zscore  = (mean - glob$mean) / glob$sd,
        foldchg = mean / glob$mean
      )
    
    dplyr::bind_rows(
      glob %>% dplyr::mutate(zscore = NA_real_, foldchg = 1),
      chr_tab,
      gene_tab
    ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
  })
  
  output$nonsyn_interpret <- DT::renderDT({
    DT::datatable(
      nonsyn_interpret_table(),
      rownames = FALSE,
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 20,
        scrollX    = TRUE
      )
    )
  }, server = FALSE)
  
  output$nonsyn_gene_bar <- renderPlot({
    dt <- nonsyn_manhattan_df()
    if (is.null(dt) || !nrow(dt) || !"GENE" %in% names(dt)) return(NULL)
    
    df <- dt %>%
      dplyr::mutate(GENE = ifelse(is.na(GENE) | GENE == "", "NA", GENE)) %>%
      dplyr::count(GENE, sort = TRUE) %>%
      dplyr::slice_head(n = 15)
    
    ggplot(df, aes(x = reorder(GENE, n), y = n)) +
      geom_col() +
      coord_flip() +
      labs(x = "Gene", y = "NonSyn variants") +
      theme_minimal(base_size = 12)
  })
  
  output$nonsyn_metric_bar <- renderPlot({
    dt <- nonsyn_manhattan_df()
    req(nrow(dt) > 0)
    
    ggplot(dt, aes(x = value)) +
      geom_histogram(bins = 30, fill = "#1A4E8A", alpha = 0.8) +
      labs(
        x = input$nonsyn_metric,
        y = "Number of variants",
        title = "Distribution of selected NonSyn metric"
      ) +
      theme_minimal()
  })
  
  # ========================================================================
  # 8) Radar (rankscore) — clean single-source of truth (dbnsfp_norm_df)
  # ========================================================================
  
  # --- Metrics groups (rankscore)
  radar_metrics_rank <- list(
    "Functional_Classic" = c(
      "SIFT_converted_rankscore",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_rankscore",
      "Polyphen2_HVAR_rankscore",
      "PROVEAN_converted_rankscore",
      "MutationAssessor_rankscore",
      "MutationTaster_converted_rankscore"
    ),
    
    "Pathogenicity_Meta" = c(
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "MetaRNN_rankscore",
      "REVEL_rankscore",
      "ClinPred_rankscore",
      "VEST4_rankscore",
      "MutPred_rankscore",
      "CADD_raw_rankscore",
      "DANN_rankscore",
      "fathmm_XF_coding_rankscore",
      "BayesDel_addAF_rankscore",
      "BayesDel_noAF_rankscore",
      "M_CAP_rankscore",
      "DEOGEN2_rankscore",
      "LIST_S2_rankscore"
    ),
    
    "Pathogenicity_DL" = c(
      "PrimateAI_rankscore",
      "AlphaMissense_rankscore",
      "ESM1b_rankscore",
      "MVP_rankscore",
      "gMVP_rankscore",
      "VARITY_R_rankscore",
      "VARITY_ER_rankscore",
      "VARITY_R_LOO_rankscore",
      "VARITY_ER_LOO_rankscore",
      "PHACTboost_rankscore",
      "MutFormer_rankscore",
      "MutScore_rankscore"
    ),
    
    "Constraint_Intolerance" = c(
      "MPC_rankscore",
      "M_CAP_rankscore",
      "CADD_raw_rankscore",
      "bStatistic_converted_rankscore"
    ),
    
    "Conservation" = c(
      "GERP++_RS_rankscore",
      "GERP_91_mammals_rankscore",
      "phyloP17way_primate_rankscore",
      "phyloP100way_vertebrate_rankscore",
      "phyloP470way_mammalian_rankscore",
      "phastCons17way_primate_rankscore",
      "phastCons100way_vertebrate_rankscore",
      "phastCons470way_mammalian_rankscore"
    ),
    
    "Clinical_UsualSuspects" = c(
      "CADD_raw_rankscore",
      "REVEL_rankscore",
      "ClinPred_rankscore",
      "MetaRNN_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "PrimateAI_rankscore",
      "AlphaMissense_rankscore",
      "VEST4_rankscore",
      "MPC_rankscore"
    )
  )
  
  radar_metrics_all <- unique(unlist(radar_metrics_rank, use.names = FALSE))
  
  radar_chr_num <- function(x){
    x <- toupper(trimws(as.character(x)))
    x <- sub("^CHR", "", x)
    x[x == "X"] <- "23"
    x[x == "Y"] <- "24"
    x[x %in% c("MT","M")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pretty_metric_2lines <- function(x) {
    x <- as.character(x)
    x <- sub("_rankscore$", "", x, ignore.case = TRUE)
    x <- sub("_", "\n", x)
    x <- gsub("_", " ", x)
    x
  }
  
  radar_base_df <- reactive({
    df <- dbnsfp_norm_df()
    if (!is.data.frame(df) || nrow(df) == 0) return(tibble::tibble())
    
    # Core columns required
    if (!all(c("chr","BP","genename","cluster_id") %in% names(df))) return(tibble::tibble())
    
    mets <- intersect(radar_metrics_all, names(df))
    if (length(mets) < 3) return(tibble::tibble())
    
    df %>%
      dplyr::mutate(
        CHR        = radar_chr_num(.data[["chr"]]),
        POS        = suppressWarnings(as.integer(.data[["BP"]])),
        genename   = trimws(as.character(.data[["genename"]])),
        cluster_id = trimws(as.character(.data[["cluster_id"]]))
      ) %>%
      dplyr::filter(is.finite(CHR), is.finite(POS)) %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(mets), \(x) suppressWarnings(as.numeric(x)))
      )
  })
  
  # A) Fill chromosomes for cluster mode
  observeEvent(radar_base_df(), {
    df <- radar_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0 || !"CHR" %in% names(df)) {
      updateSelectInput(session, "radar_chr_for_cluster", choices = character(0), selected = NULL)
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    chr_choices <- sort(unique(suppressWarnings(as.integer(df$CHR))))
    chr_choices <- chr_choices[is.finite(chr_choices)]
    
    if (!length(chr_choices)) {
      updateSelectInput(session, "radar_chr_for_cluster", choices = character(0), selected = NULL)
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    chr_labels <- paste0("chr", chr_label_plink(chr_choices))
    chr_map    <- stats::setNames(as.character(chr_choices), chr_labels)
    
    sel_chr <- if (!is.null(input$radar_chr_for_cluster) &&
                   input$radar_chr_for_cluster %in% unname(chr_map)) {
      input$radar_chr_for_cluster
    } else {
      unname(chr_map)[1]
    }
    
    updateSelectInput(
      session,
      "radar_chr_for_cluster",
      choices = chr_map,
      selected = sel_chr
    )
  }, ignoreInit = FALSE)
  
  
  # B) Fill clusters for selected chromosome
  observeEvent(list(radar_base_df(), input$radar_chr_for_cluster), {
    df <- radar_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    chr_sel <- suppressWarnings(as.integer(input$radar_chr_for_cluster))
    if (!is.finite(chr_sel)) {
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    ids <- df %>%
      dplyr::filter(
        CHR == chr_sel,
        !is.na(cluster_id), nzchar(cluster_id)
      ) %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::arrange(cluster_id) %>%
      dplyr::pull(cluster_id)
    
    if (!length(ids)) {
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    sel_id <- if (!is.null(input$radar_cluster_id) &&
                  input$radar_cluster_id %in% ids) {
      input$radar_cluster_id
    } else {
      ids[1]
    }
    
    updateSelectInput(
      session,
      "radar_cluster_id",
      choices = ids,
      selected = sel_id
    )
  }, ignoreInit = FALSE)
  
  
  # C) Fill genes grouped by cluster (optgroup)
  observeEvent(radar_base_df(), {
    df <- radar_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "radar_gene", choices = character(0), selected = NULL)
      return()
    }
    
    dd <- df %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(genename),   nzchar(genename)
      ) %>%
      dplyr::distinct(cluster_id, genename) %>%
      dplyr::arrange(cluster_id, genename)
    
    if (!nrow(dd)) {
      updateSelectInput(session, "radar_gene", choices = character(0), selected = NULL)
      return()
    }
    
    choices_list <- dd %>%
      dplyr::group_split(cluster_id, .keep = TRUE) %>%
      purrr::set_names(purrr::map_chr(., ~ .x$cluster_id[1])) %>%
      purrr::map(function(x) {
        keys <- paste0(x$cluster_id, "||", x$genename)
        stats::setNames(keys, x$genename)
      })
    
    all_keys <- unlist(choices_list, use.names = FALSE)
    
    sel <- if (!is.null(input$radar_gene) && input$radar_gene %in% all_keys) {
      input$radar_gene
    } else {
      all_keys[1]
    }
    
    updateSelectInput(
      session,
      "radar_gene",
      choices = choices_list,
      selected = sel
    )
  }, ignoreInit = FALSE)
  
  
  # D) Aggregate row depending on mode (global / cluster / gene)
  radar_group_row <- reactive({
    df <- radar_base_df()
    if (!is.data.frame(df) || nrow(df) == 0) return(NULL)
    
    mets <- intersect(radar_metrics_all, names(df))
    if (length(mets) < 3) return(NULL)
    
    mode <- input$radar_mode %||% "cluster"
    
    if (mode == "global") {
      df_use <- df
      label  <- "Global"
      
    } else if (mode == "cluster") {
      chr_sel <- suppressWarnings(as.integer(input$radar_chr_for_cluster))
      cl      <- as.character(input$radar_cluster_id %||% "")
      
      if (!is.finite(chr_sel) || !nzchar(cl)) return(NULL)
      
      df_use <- df %>%
        dplyr::filter(
          CHR == chr_sel,
          cluster_id == cl
        )
      
      label <- paste0(cl, " (chr", chr_label_plink(chr_sel), ")")
      
    } else if (mode == "gene") {
      key <- as.character(input$radar_gene %||% "")
      if (!nzchar(key)) return(NULL)
      
      parts <- strsplit(key, "\\|\\|")[[1]]
      cluster_sel <- if (length(parts) >= 1) parts[1] else ""
      gene_sel    <- if (length(parts) >= 2) parts[2] else ""
      
      if (!nzchar(cluster_sel) || !nzchar(gene_sel)) return(NULL)
      
      df_use <- df %>%
        dplyr::filter(
          cluster_id == cluster_sel,
          genename == gene_sel
        )
      
      label <- paste0(gene_sel, " (", cluster_sel, ")")
      
    } else {
      return(NULL)
    }
    
    if (!is.data.frame(df_use) || nrow(df_use) == 0) return(NULL)
    
    row <- df_use %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(mets),
          ~ mean(.x, na.rm = TRUE)
        )
      )
    
    row$label <- label
    row
  })
  
  # E) Radar plot
  output$nonsyn_radar <- renderPlot({
    row <- radar_group_row()
    validate(need(is.data.frame(row) && nrow(row) == 1, "Radar: select a valid Cluster / Chromosome / Gene."))
    
    radar_long <- purrr::imap_dfr(radar_metrics_rank, function(mset, grp) {
      mset <- intersect(as.character(mset), names(row))
      if (!length(mset)) return(NULL)
      
      metric_lbl <- pretty_metric_2lines(mset)
      
      tibble::tibble(
        group  = grp,
        metric = factor(metric_lbl, levels = metric_lbl),
        value  = as.numeric(row[1, mset])
      )
    })
    
    radar_long <- radar_long %>%
      dplyr::filter(is.finite(value)) %>%
      dplyr::mutate(
        value = pmin(pmax(value, 0), 1),
        group = factor(group, levels = names(radar_metrics_rank))
      )
    
    validate(need(nrow(radar_long) >= 3, "Not enough rankscore metrics for radar plot (need >= 3)."))
    
    ggplot2::ggplot(radar_long, ggplot2::aes(x = metric, y = value, group = 1)) +
      ggplot2::geom_polygon(fill = "grey70", alpha = 0.18, colour = NA) +
      ggplot2::geom_path(colour = "grey35", linewidth = 0.7) +
      ggplot2::geom_point(ggplot2::aes(colour = value), size = 4) +
      ggplot2::coord_polar() +
      ggplot2::facet_wrap(~ group, scales = "free_x", ncol = 3) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_colour_gradient(
        limits = c(0, 1),
        low = "yellow",
        high = "red3",
        name = "Rankscore"
      ) +
      ggplot2::labs(
        title = paste0("Clinical radar – ", row$label),
        y = "Rankscore (0–1)",
        x = NULL
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 12),
        strip.text = ggplot2::element_text(face = "bold")
      )
  })
  
  output$nonsyn_radar_plotly <- plotly::renderPlotly({
    row <- radar_group_row()
    validate(need(
      is.data.frame(row) && nrow(row) == 1,
      "Radar: select a valid Cluster / Chromosome / Gene."
    ))
    
    # ----------------------------------------------------------
    # Build long table
    # ----------------------------------------------------------
    radar_long <- purrr::imap_dfr(radar_metrics_rank, function(mset, grp) {
      mset <- intersect(as.character(mset), names(row))
      if (!length(mset)) return(NULL)
      
      metric_lbl <- pretty_metric_2lines(mset)
      metric_lbl <- gsub("\n", "<br>", metric_lbl, fixed = TRUE)
      metric_lbl[is.na(metric_lbl) | !nzchar(metric_lbl)] <- mset[is.na(metric_lbl) | !nzchar(metric_lbl)]
      
      tibble::tibble(
        group      = as.character(grp),
        metric     = as.character(metric_lbl),
        metric_raw = as.character(mset),
        value      = as.numeric(row[1, mset])
      )
    })
    
    radar_long <- radar_long %>%
      dplyr::filter(
        is.finite(value),
        !is.na(group), nzchar(group),
        !is.na(metric), nzchar(metric)
      ) %>%
      dplyr::mutate(
        value = pmin(pmax(value, 0), 1)
      )
    
    validate(need(
      nrow(radar_long) >= 3,
      "Not enough rankscore metrics for radar plot (need >= 3)."
    ))
    
    groups_present <- unique(radar_long$group)
    validate(need(length(groups_present) > 0, "No radar groups available."))
    
    # ----------------------------------------------------------
    # Grid geometry
    # ----------------------------------------------------------
    n_grp <- length(groups_present)
    ncol  <- min(2, n_grp)
    nrow  <- ceiling(n_grp / ncol)
    
    x_gap       <- 0.05   # prova entre 0.04 i 0.12
    y_gap       <- 0.07   # prova entre 0.06 i 0.16
    
    cell_w <- (1 - x_gap * (ncol - 1)) / ncol
    cell_h <- (1 - y_gap * (nrow - 1)) / nrow
    
    inner_pad_x <- 0.03   # prova entre 0.02 i 0.10
    inner_pad_y <- 0.04   # prova entre 0.03 i 0.12
    
    # ----------------------------------------------------------
    # Base figure
    # ----------------------------------------------------------
    fig <- plotly::plot_ly()
    layout_args <- list()
    annotations <- list()
  
    
    for (i in seq_along(groups_present)) {
      grp <- groups_present[i]
      d <- radar_long %>% dplyr::filter(group == grp)
      
      validate(need(nrow(d) >= 3, paste("Not enough metrics in group:", grp)))
      
      # close polygon
      d2 <- dplyr::bind_rows(d, d[1, , drop = FALSE])
      
      # panel position
      col_i <- ((i - 1) %% ncol) + 1
      row_i <- ((i - 1) %/% ncol) + 1
      
      x0 <- (col_i - 1) * (cell_w + x_gap)
      x1 <- x0 + cell_w
      
      # row 1 = top
      y1 <- 1 - (row_i - 1) * (cell_h + y_gap)
      y0 <- y1 - cell_h
      
      polar_id <- if (i == 1) "polar" else paste0("polar", i)
      
      fig <- fig %>%
        plotly::add_trace(
          data = d2,
          type = "scatterpolar",
          mode = "lines+markers",
          subplot = polar_id,
          r = ~value,
          theta = ~metric,
          customdata = ~metric_raw,
          hovertemplate = paste0(
            "<b>", as.character(row$label[[1]]), "</b><br>",
            "Group: ", grp, "<br>",
            "Metric: %{customdata}<br>",
            "Rankscore: %{r:.3f}<extra></extra>"
          ),
          fill = "toself",
          fillcolor = "rgba(160,160,160,0.18)",
          line = list(color = "rgba(80,80,80,0.9)", width = 2),
          marker = list(
            size = 10,
            color = ~value,
            colorscale = list(
              list(0.00, "yellow"),
              list(0.50, "orange"),
              list(1.00, "red")
            ),
            cmin = 0,
            cmax = 1,
            showscale = TRUE,
            colorbar = list(
              title = list(text = "Rankscore", side = "right"),
              len = 0.35,
              thickness = 12,
              y = 0.5,
              yanchor = "middle",
              tickfont = list(size = 10)
            ),
            line = list(color = "rgba(60,60,60,0.8)", width = 1)
          ),
          showlegend = FALSE,
          inherit = FALSE
        )
      
      xp0 <- x0 + inner_pad_x
      xp1 <- x1 - inner_pad_x
      yp0 <- y0 + inner_pad_y
      yp1 <- y1 - inner_pad_y
      
      layout_args[[polar_id]] <- list(
        domain = list(x = c(xp0, xp1), y = c(yp0, yp1)),
        bgcolor = "white",
        radialaxis = list(
          range = c(0, 1),
          tickvals = c(0, 0.25, 0.5, 0.75, 1),
          ticktext = c("0", "0.25", "0.5", "0.75", "1"),
          angle = 90,
          gridcolor = "rgba(180,180,180,0.5)",
          linecolor = "rgba(180,180,180,0.8)"
        ),
        angularaxis = list(
          direction = "clockwise",
          rotation = 90,
          tickfont = list(size = 10),
          gridcolor = "rgba(220,220,220,0.35)",
          linecolor = "rgba(180,180,180,0.75)"
        )
      )
      
      annotations[[length(annotations) + 1]] <- list(
        text = paste0("<b>", grp, "</b>"),
        x = (x0 + x1) / 2,
        y = min(1, y1 + 0.04),
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
        font = list(size = 13)
      )
    }
    
    layout_args$title <- list(
      text = paste0("Clinical radar – ", as.character(row$label[[1]])),
      x = 0.02,
      xanchor = "left"
    )
    layout_args$paper_bgcolor <- "white"
    layout_args$plot_bgcolor  <- "white"
    layout_args$margin <- list(l = 20, r = 20, t = 80, b = 20)
    layout_args$annotations <- annotations
    
    fig <- do.call(plotly::layout, c(list(p = fig), layout_args))
    
    fig %>% plotly::config(displayModeBar = TRUE)
  })
  
  ##############################################################################
  
  
  
  # F) Radar selection table (summary of IDs + genes with links)
  radar_selected_df <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "Final output not available yet."))
    validate(need(all(c("cluster_id", "genename") %in% names(df)),
                  "Final output must contain: cluster_id, genename"))
    
    mode <- input$radar_mode %||% "cluster"
    chr_num <- suppressWarnings(as.integer(norm_chr_generic(df[["chr"]])))
    
    if (mode == "global") {
      df_use <- df
      label  <- "Global"
      
    } else if (mode == "cluster") {
      req(input$radar_chr_for_cluster, input$radar_cluster_id)
      
      chr_sel <- suppressWarnings(as.integer(input$radar_chr_for_cluster))
      cl_sel  <- as.character(input$radar_cluster_id %||% "")
      
      req(is.finite(chr_sel), nzchar(cl_sel))
      
      df_use <- df %>%
        dplyr::mutate(.CHR = chr_num) %>%
        dplyr::filter(
          .CHR == chr_sel,
          !is.na(cluster_id),
          cluster_id == cl_sel
        )
      
      label <- paste0(cl_sel, " (chr", chr_label_plink(chr_sel), ")")
      
    } else if (mode == "gene") {
      req(input$radar_gene)
      
      parts <- strsplit(as.character(input$radar_gene), "\\|\\|")[[1]]
      cluster_sel <- parts[1] %||% ""
      gene_sel    <- parts[2] %||% ""
      
      req(nzchar(cluster_sel), nzchar(gene_sel))
      
      df_use <- df %>%
        dplyr::filter(
          cluster_id == cluster_sel,
          genename == gene_sel
        )
      
      label <- paste0(gene_sel, " (", cluster_sel, ")")
      
    } else {
      return(NULL)
    }
    
    validate(need(nrow(df_use) > 0, "No variants for the current radar selection."))
    attr(df_use, "label") <- label
    df_use
  })
  
  output$radar_selection_table <- DT::renderDT({
    df_use <- radar_selected_df()
    req(is.data.frame(df_use))
    
    label <- attr(df_use, "label") %||% ""
    mode  <- input$radar_mode %||% "cluster"
    
    clean_unique <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x[x %in% c("", ".", "NA", "NaN", "nan")] <- NA_character_
      x <- unique(na.omit(x))
      x[nzchar(x)]
    }
    
    pick_variant_id_col <- function(df) {
      cand <- c("rsid", "SNP", "snp", "ID", "variant_id")
      cand[cand %in% names(df)][1] %||% NA_character_
    }
    
    make_dbsnp_url <- function(id) {
      id <- trimws(as.character(id))
      if (!nzchar(id)) return(NA_character_)
      if (grepl("^rs\\d+$", id, ignore.case = TRUE)) {
        paste0("https://www.ncbi.nlm.nih.gov/snp/", id)
      } else {
        paste0("https://www.ncbi.nlm.nih.gov/snp/?term=", utils::URLencode(id, reserved = TRUE))
      }
    }
    
    make_genecards_url <- function(g) {
      g <- trimws(as.character(g))
      if (!nzchar(g)) return(NA_character_)
      paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", utils::URLencode(g, reserved = TRUE))
    }
    
    linkify <- function(x, url_fun) {
      u <- url_fun(x)
      if (is.na(u) || !nzchar(u)) return(htmltools::htmlEscape(x))
      sprintf("<a href='%s' target='_blank'>%s</a>", u, htmltools::htmlEscape(x))
    }
    
    collapse_links <- function(vec, max_items = 250, url_fun) {
      vec <- clean_unique(vec)
      if (!length(vec)) return("")
      more <- 0L
      if (length(vec) > max_items) {
        more <- length(vec) - max_items
        vec  <- vec[seq_len(max_items)]
      }
      html <- paste(vapply(vec, linkify, character(1), url_fun = url_fun), collapse = ", ")
      if (more > 0) html <- paste0(html, sprintf(" … (+%d more)", more))
      html
    }
    
    id_col <- pick_variant_id_col(df_use)
    ids <- if (is.na(id_col)) {
      paste0("chr", norm_chr_generic(df_use[["chr"]]), ":", df_use[["BP"]])
    } else {
      df_use[[id_col]]
    }
    
    genes <- df_use[["genename"]]
    ids_u   <- clean_unique(ids)
    genes_u <- clean_unique(genes)
    
    out <- tibble::tibble(
      selection_type = mode,
      selection_id   = label,
      n_nonsyn_snps  = length(ids_u),
      n_genes        = length(genes_u),
      nonsyn_snps    = collapse_links(ids_u, 250, make_dbsnp_url),
      genes          = collapse_links(genes_u, 250, make_genecards_url)
    )
    
    DT::datatable(
      out,
      rownames = FALSE,
      escape = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy",  exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)),
          list(extend = "excel", exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE))
        ),
        pageLength = 1,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  # ========================================================================
  # 9) UCSC viewer (GWAS + NonSyn) — custom tracks (dbSNP clickable)
  #   - Només compta hits dins la finestra seleccionada (zoom/rang X)
  #   - Links només apareixen quan hi ha selecció vàlida i hits > 0
  #   - ADAPTAT a subplot + source="manhattan_combo"
  # ========================================================================
  
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
  
  safe_label_nonsyn <- function(df) {
    # Prioritza rsid/rsid_final si existeix
    if ("rsid" %in% names(df) && !all(is.na(df$rsid))) {
      out <- as.character(df$rsid)
      bad <- is.na(out) | out == ""
      out[bad] <- paste0("nonsyn_hit_", seq_len(sum(bad)))
      return(out)
    }
    if ("rsid_final" %in% names(df) && !all(is.na(df$rsid_final))) {
      out <- as.character(df$rsid_final)
      bad <- is.na(out) | out == ""
      out[bad] <- paste0("nonsyn_hit_", seq_len(sum(bad)))
      return(out)
    }
    paste0("nonsyn_hit_", seq_len(nrow(df)))
  }
  
  # Estat de selecció (igual que GTEx)
  window_selected <- reactiveVal(NULL)
  ucsc_region     <- reactiveVal(NULL)
  
  # --- helper robust: extreu rang X d'un relayout (accepta list o data.frame) ---
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
    
    # autorange pot ser TRUE/"TRUE"/1
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
    
    # subplot: pot venir per xaxis o xaxis2
    x0 <- getv("xaxis.range[0]")
    x1 <- getv("xaxis.range[1]")
    
    if (is.null(x0) || is.null(x1)) {
      x0 <- getv("xaxis2.range[0]")
      x1 <- getv("xaxis2.range[1]")
    }
    
    # alguns plotly envien "xaxis.range" com vector length 2
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
  
  # IMPORTANT: escolta el relayout del combo correcte
  observeEvent(plotly::event_data("plotly_relayout", source = "manhattan_combo"), {
    
    ev <- plotly::event_data("plotly_relayout", source = "manhattan_combo")
    xr <- extract_xrange_from_relayout(ev)
    
    if (is.null(xr)) {
      window_selected(NULL)
      ucsc_region(NULL)
      session$userData$track_gwas_data   <- NULL
      session$userData$track_nonsyn_data <- NULL
      return()
    }
    
    window_selected(xr)
    
  }, ignoreInit = TRUE)
  
  # -------------------------
  # Hits dins la finestra (GWAS)
  # -------------------------
  get_window_hits_gwas <- reactive({
    win <- window_selected()
    req(win)
    
    dfp <- dfp_manhattan()
    df  <- tryCatch(gwas_df(), error = function(e) NULL)
    
    req(is.data.frame(dfp), nrow(dfp) > 0)
    
    dfp2 <- dfp %>%
      dplyr::mutate(CHR = as.integer(CHR), BP = as.integer(BP))
    
    if (is.data.frame(df) && nrow(df) > 0) {
      df2 <- df %>%
        dplyr::mutate(CHR = as.integer(CHR), BP = as.integer(BP)) %>%
        dplyr::select(dplyr::any_of(c("CHR","BP","rsid","snp")))
      out <- dfp2 %>%
        dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax) %>%
        dplyr::left_join(df2, by = c("CHR","BP"))
      return(out)
    }
    
    # fallback: només dfp (sense join)
    dfp2 %>% dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  # -------------------------
  # Hits dins la finestra (NonSyn)
  # -------------------------
  get_window_hits_nonsyn <- reactive({
    win <- window_selected()
    req(win)
    
    ns <- nonsyn_manhattan_df()
    req(is.data.frame(ns), nrow(ns) > 0)
    
    # assegura BPcum
    if (!"BPcum" %in% names(ns)) {
      # si no ve BPcum, el reconstruïm amb ref hg38
      ref <- .ref_hg38
      chr_col <- if ("CHR" %in% names(ns)) "CHR" else if ("chr" %in% names(ns)) "chr" else NA_character_
      pos_col <- if ("POS" %in% names(ns)) "POS" else if ("pos" %in% names(ns)) "pos" else NA_character_
      validate(need(!is.na(chr_col) && !is.na(pos_col), "NonSyn table must have CHR/POS (or chr/pos) to compute BPcum."))
      
      ns <- ns %>%
        dplyr::mutate(
          chrN = norm_chr_generic(.data[[chr_col]]),
          pos1 = suppressWarnings(as.numeric(.data[[pos_col]]))
        ) %>%
        dplyr::inner_join(ref %>% dplyr::select(chr, chr_cum), by = c("chrN" = "chr")) %>%
        dplyr::mutate(BPcum = as.numeric(pos1) + as.numeric(chr_cum))
    }
    
    ns %>% dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  # -------------------------
  # Track builders
  # -------------------------
  make_track_df <- function(df, R, G, B, type = c("gwas","nonsyn")) {
    type <- match.arg(type)
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    
    if (type == "gwas") {
      if (!all(c("CHR","BP") %in% names(df))) return(tibble::tibble())
      chrom <- as.character(df$CHR)
      chrom <- gsub("^chr", "", chrom, ignore.case = TRUE)
      chrom <- paste0("chr", chr_label_plink(as.integer(chrom)))
      start <- suppressWarnings(as.numeric(df$BP) - 1)
      end   <- suppressWarnings(as.numeric(df$BP))
      lbl   <- safe_label_rsid(df)
      
    } else {
      # NonSyn: CHR + POS (o chr+pos)
      chr_col <- if ("CHR" %in% names(df)) "CHR" else if ("chr" %in% names(df)) "chr" else NA_character_
      pos_col <- if ("POS" %in% names(df)) "POS" else if ("pos" %in% names(df)) "pos" else NA_character_
      if (is.na(chr_col) || is.na(pos_col)) return(tibble::tibble())
      
      chr_use <- as.character(df[[chr_col]])
      chr_use <- gsub("^chr", "", chr_use, ignore.case = TRUE)
      chrom   <- paste0("chr", chr_use)
      
      start <- suppressWarnings(as.numeric(df[[pos_col]]) - 1)
      end   <- suppressWarnings(as.numeric(df[[pos_col]]))
      lbl   <- safe_label_nonsyn(df)
    }
    
    start[is.na(start)] <- 0
    end[is.na(end)]     <- start[is.na(end)] + 1
    end <- ifelse(end < start, start + 1, end)
    
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
      dplyr::filter(
        !is.na(chrom),
        chrom != "chrNA",
        !is.na(start),
        !is.na(end),
        start >= 0,
        end >= start
      ) %>%
      dplyr::mutate(
        chrom      = as.character(chrom),
        strand     = dplyr::if_else(strand %in% c("+", "-"), strand, "+"),
        thickStart = ifelse(is.na(thickStart), start, thickStart),
        thickEnd   = ifelse(is.na(thickEnd),   end,   thickEnd)
      )
  }
  
  # -------------------------
  # Quan hi ha finestra seleccionada: calcula regió UCSC + tracks
  # -------------------------
  observeEvent(window_selected(), {
    
    win <- window_selected()
    if (is.null(win)) return()
    
    left  <- coord_from_bp_cum(win$xmin)
    right <- coord_from_bp_cum(win$xmax)
    
    # UCSC region ha d’estar en un sol cromosoma
    if (is.null(left$chr) || is.null(right$chr) || left$chr != right$chr) {
      ucsc_region(NULL)
      session$userData$track_gwas_data   <- NULL
      session$userData$track_nonsyn_data <- NULL
      return()
    }
    
    region <- sprintf("chr%s:%d-%d", left$chr, left$pos, right$pos)
    ucsc_region(region)
    
    gwas_hits <- tryCatch(get_window_hits_gwas(),   error = function(e) NULL)
    ns_hits   <- tryCatch(get_window_hits_nonsyn(), error = function(e) NULL)
    
    if ((is.null(gwas_hits) || nrow(gwas_hits) == 0) &&
        (is.null(ns_hits)   || nrow(ns_hits)   == 0)) {
      session$userData$track_gwas_data   <- NULL
      session$userData$track_nonsyn_data <- NULL
      return()
    }
    
    df_gwas <- clean_track(make_track_df(gwas_hits,  31,120,180, type = "gwas"))
    df_ns   <- clean_track(make_track_df(ns_hits,   227, 26, 28, type = "nonsyn"))
    
    session$userData$track_gwas_data   <- df_gwas
    session$userData$track_nonsyn_data <- df_ns
    
  }, ignoreInit = TRUE)
  
  # -------------------------
  # UCSC text + URL + links
  # -------------------------
  make_ucsc_track_text <- function(name, df, url_base = "https://www.ncbi.nlm.nih.gov/snp/$$") {
    header <- sprintf(
      'track name="%s" description="%s" visibility=pack itemRgb=On url="%s"',
      name, name, url_base
    )
    if (is.null(df) || !nrow(df)) return(header)
    
    cols_bed12 <- c("chrom","start","end","name","score","strand",
                    "thickStart","thickEnd","itemRgb",
                    "blockCount","blockSizes","blockStarts")
    body <- apply(df[, cols_bed12], 1, paste, collapse = "\t")
    paste(c(header, body), collapse = "\n")
  }
  
  make_ucsc_url <- function(region, track_text) {
    base <- "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38"
    reg  <- paste0("&position=", utils::URLencode(region,  reserved = TRUE))
    trk  <- paste0("&hgt.customText=", utils::URLencode(track_text, reserved = TRUE))
    paste0(base, reg, trk)
  }
  
  output$ucsc_link_gwas <- renderUI({
    region <- ucsc_region()
    df     <- session$userData$track_gwas_data
    req(!is.null(region), !is.null(df), nrow(df) > 0)
    
    txt <- make_ucsc_track_text("GWAS_hits", df, url_base = "https://www.ncbi.nlm.nih.gov/snp/$$")
    url <- make_ucsc_url(region, txt)
    
    tags$a(href = url, target = "_blank", "Open UCSC – GWAS hits (dbSNP links)")
  })
  
  output$ucsc_link_nonsyn <- renderUI({
    region <- ucsc_region()
    df     <- session$userData$track_nonsyn_data
    req(!is.null(region), !is.null(df), nrow(df) > 0)
    
    txt <- make_ucsc_track_text("NonSyn_hits", df, url_base = "https://www.ncbi.nlm.nih.gov/snp/$$")
    url <- make_ucsc_url(region, txt)
    
    tags$a(href = url, target = "_blank", "Open UCSC – NonSyn hits (dbSNP links)")
  })
  
  output$debug_ucsc_state <- renderUI({
    region <- ucsc_region() %||% "NULL"
    gwas_n <- nrow(session$userData$track_gwas_data %||% tibble::tibble())
    ns_n   <- nrow(session$userData$track_nonsyn_data %||% tibble::tibble())
    
    tags$pre(
      style="background-color:#f6f6f6; border:1px solid #ddd; padding:10px; font-family:'Courier New', monospace;",
      paste0(
        "region = ", region,
        "\nGWAS hits in track = ", gwas_n,
        "\nNonSyn hits in track = ", ns_n,
        "\nwindow_selected = ", if (is.null(window_selected())) "NULL" else paste0(window_selected()$xmin, " .. ", window_selected()$xmax)
      )
    )
  })
  
  # Keep dbNSFP log rendering
  output$dbnsfp_log <- renderText(dbnsfp_log_text())
  
  # ========================================================================
  # 10) Metric thresholds + shared helpers (used by boxplot combo)
  # ========================================================================
  
  nonsyn_metric_thresholds <- list(
    
    # -------------------------
    # DEFAULT
    # -------------------------
    SIFT_score               = list(type = "lower",  value = 0.05),
    SIFT_converted_rankscore = list(type = "higher", value = 0.8),
    
    Polyphen2_HDIV_score     = list(type = "higher", value = 0.909),
    Polyphen2_HDIV_rankscore = list(type = "higher", value = 0.8),
    
    Polyphen2_HVAR_score     = list(type = "higher", value = 0.909),
    Polyphen2_HVAR_rankscore = list(type = "higher", value = 0.8),
    
    MutationAssessor_score     = list(type = "higher", value = 2.0),
    MutationAssessor_rankscore = list(type = "higher", value = 0.8),
    
    PROVEAN_score     = list(type = "lower",  value = -2.5),
    PROVEAN_rankscore = list(type = "higher", value = 0.8),
    
    REVEL_score     = list(type = "higher", value = 0.5),
    REVEL_rankscore = list(type = "higher", value = 0.8),
    
    # -------------------------
    # PATHOGENICITY
    # -------------------------
    CADD_phred         = list(type = "higher", value = 20),
    CADD_raw_rankscore = list(type = "higher", value = 0.8),
    
    MetaRNN_score     = list(type = "higher", value = 0.5),
    MetaRNN_rankscore = list(type = "higher", value = 0.8),
    
    MetaSVM_score     = list(type = "higher", value = 0),
    MetaSVM_rankscore = list(type = "higher", value = 0.8),
    
    MetaLR_score     = list(type = "higher", value = 0),
    MetaLR_rankscore = list(type = "higher", value = 0.8),
    
    ClinPred_score     = list(type = "higher", value = 0.5),
    ClinPred_rankscore = list(type = "higher", value = 0.8),
    
    PrimateAI_score     = list(type = "higher", value = 0.8),
    PrimateAI_rankscore = list(type = "higher", value = 0.8),
    
    AlphaMissense_score     = list(type = "higher", value = 0.56),
    AlphaMissense_rankscore = list(type = "higher", value = 0.8),
    
    MutPred_score     = list(type = "higher", value = 0.68),
    MutPred_rankscore = list(type = "higher", value = 0.8),
    
    VEST4_score     = list(type = "higher", value = 0.5),
    VEST4_rankscore = list(type = "higher", value = 0.8),
    
    MPC_score     = list(type = "higher", value = 2),
    MPC_rankscore = list(type = "higher", value = 0.8),
    
    # -------------------------
    # CONSERVATION
    # -------------------------
    `GERP++_RS`           = list(type = "higher", value = 2),
    GERP_91_mammals       = list(type = "higher", value = 2),
    `GERP++_RS_rankscore` = list(type = "higher", value = 0.8),
    GERP_91_mammals_rankscore = list(type = "higher", value = 0.8),
    
    phyloP17way_primate          = list(type = "higher", value = 2),
    phyloP17way_primate_rankscore = list(type = "higher", value = 0.8),
    
    phyloP100way_vertebrate = list(type = "higher", value = 2),
    
    phastCons17way_primate          = list(type = "higher", value = 0.8),
    phastCons17way_primate_rankscore = list(type = "higher", value = 0.8),
    
    # -------------------------
    # FUNCTIONAL IMPACT
    # -------------------------
    DANN_score     = list(type = "higher", value = 0.9),
    DANN_rankscore = list(type = "higher", value = 0.8),
    
    Eigen_raw_coding           = list(type = "higher", value = 2),
    Eigen_raw_coding_rankscore = list(type = "higher", value = 0.8),
    
    Eigen_PC_raw_coding           = list(type = "higher", value = 2),
    Eigen_PC_raw_coding_rankscore = list(type = "higher", value = 0.8),
    
    bStatistic_converted           = list(type = "higher", value = 1),
    bStatistic_converted_rankscore = list(type = "higher", value = 0.8),
    
    # -------------------------
    # LoF / HAPLOINSUFFICIENCY
    # -------------------------
    ExAC_pLI   = list(type = "higher", value = 0.9),
    gnomAD_pLI = list(type = "higher", value = 0.9),
    
    # -------------------------
    # GENE DAMAGE
    # -------------------------
    GDI                         = list(type = "higher", value = 13),
    GDI_Phred                   = list(type = "higher", value = 20),
    LoFtool_score               = list(type = "lower",  value = 0.1),
    Gene_indispensability_score = list(type = "higher", value = 0.7)
  )
  
  get_metric_threshold <- function(metric) {
    nonsyn_metric_thresholds[[metric]]
  }
  
  get_threshold_color <- function(thr) {
    if (thr$type == "higher") "red" else "darkgreen"
  }
  
  wilcox_vs_bg <- function(x, bg) {
    if (length(x) < 2 || length(bg) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(x, bg)$p.value)
  }
  
  # ========================================================================
  # 11) Summary boxplots (chr / gene / cluster) in one Plotly row
  # ========================================================================
  
  output$nonsyn_box_combo_plotly <- renderPlotly({
    
    dt <- nonsyn_manhattan_df()
    
    req(is.data.frame(dt), nrow(dt) > 0, input$nonsyn_metric)
    
    metric <- input$nonsyn_metric
    thr    <- get_metric_threshold(metric)
    
    # Shared Y range across panels
    yv <- dt$value
    yv <- yv[is.finite(yv)]
    req(length(yv) > 0)
    
    y_lim <- range(yv, na.rm = TRUE)
    pad <- diff(y_lim) * 0.05
    if (!is.finite(pad) || pad == 0) pad <- 1
    y_lim <- c(y_lim[1] - pad, y_lim[2] + pad)
    
    # Add threshold layer consistently
    add_thr_layer <- function(p, df, x_f) {
      if (is.null(thr) || is.null(thr$type) || is.null(thr$value)) return(p)
      
      col_thr <- get_threshold_color(thr)
      label_txt <- if (thr$type == "higher") {
        paste("Deleterious ≥", thr$value)
      } else {
        paste("Deleterious ≤", thr$value)
      }
      
      y_pos <- y_lim[1] + diff(y_lim) * 0.92
      x_pos <- mean(as.numeric(df[[x_f]]), na.rm = TRUE)
      label_df <- data.frame(x = x_pos, y = y_pos, label = label_txt)
      
      p +
        ggplot2::geom_hline(
          yintercept = thr$value,
          color = col_thr,
          linetype = "dashed",
          linewidth = 0.8
        ) +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          color = col_thr,
          size = 3.0,
          fontface = "bold",
          hjust = 0.5
        )
    }
    
    # 1) Chromosome boxplot
    chr_levels <- sort(unique(as.integer(dt$CHR)))
    
    dt_chr <- dt %>%
      dplyr::mutate(
        chr_f    = factor(CHR, levels = chr_levels),
        pos_id   = as.integer(chr_f),
        fill_grp = ifelse(pos_id %% 2 == 0, "even", "odd")
      )
    
    stats_chr <- dt_chr %>%
      dplyr::group_by(chr_f) %>%
      dplyr::summarise(
        p_wilcox = wilcox_vs_bg(value, dt_chr$value[dt_chr$CHR != unique(CHR)]),
        .groups = "drop"
      )
    
    dt_chr <- dt_chr %>% dplyr::left_join(stats_chr, by = "chr_f")
    
    p_chr <- ggplot2::ggplot(
      dt_chr,
      ggplot2::aes(
        x = chr_f, y = value, fill = fill_grp,
        text = paste0(
          "chr", CHR,
          "<br>value: ", signif(value, 4),
          "<br>Wilcox p: ", signif(p_wilcox, 3)
        )
      )
    ) +
      ggplot2::geom_boxplot(outlier.alpha = 0.25) +
      ggplot2::scale_fill_manual(values = c(even = "darkgreen", odd = "orange"), guide = "none") +
      ggplot2::coord_cartesian(ylim = y_lim) +
      # ggplot2::labs(x = "Chromosome", y = metric, title = "By chromosome") +
      ggplot2::labs(x = "Chromosome", y = metric) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 11, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 9),
        axis.title.y = ggplot2::element_text(size = 9),
        axis.text.x  = ggplot2::element_text(angle = 90, size = 7),
        axis.text.y  = ggplot2::element_text(size = 8),
        legend.position = "none"
      )
    
    p_chr <- add_thr_layer(p_chr, dt_chr, "chr_f")
    p_chr_pl <- plotly::ggplotly(p_chr, tooltip = "text")
    
    # 2) Gene boxplot (top 30)
    req("GENE" %in% names(dt))
    
    gene_levels <- dt %>%
      dplyr::count(GENE, sort = TRUE) %>%
      dplyr::slice_head(n = 30) %>%
      dplyr::pull(GENE)
    
    dt_gene <- dt %>%
      dplyr::filter(GENE %in% gene_levels) %>%
      dplyr::mutate(
        gene_f   = factor(GENE, levels = gene_levels),
        pos_id   = as.integer(gene_f),
        fill_grp = ifelse(pos_id %% 2 == 0, "even", "odd")
      )
    
    stats_gene <- dt_gene %>%
      dplyr::group_by(gene_f) %>%
      dplyr::summarise(
        p_wilcox = wilcox_vs_bg(value, dt_gene$value[dt_gene$GENE != unique(GENE)]),
        .groups = "drop"
      )
    
    dt_gene <- dt_gene %>% dplyr::left_join(stats_gene, by = "gene_f")
    
    p_gene <- ggplot2::ggplot(
      dt_gene,
      ggplot2::aes(
        x = gene_f, y = value, fill = fill_grp,
        text = paste0(
          "Gene: ", GENE,
          "<br>value: ", signif(value, 4),
          "<br>Wilcox p: ", signif(p_wilcox, 3)
        )
      )
    ) +
      ggplot2::geom_boxplot(outlier.alpha = 0.25) +
      ggplot2::scale_fill_manual(values = c(even = "darkgreen", odd = "orange"), guide = "none") +
      ggplot2::coord_cartesian(ylim = y_lim) +
      #  ggplot2::labs(x = "Gene (top 30)", y = metric, title = "By gene") +
      ggplot2::labs(x = "Gene (top 30)", y = metric) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 11, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 9),
        axis.title.y = ggplot2::element_text(size = 9),
        axis.text.x  = ggplot2::element_text(angle = 90, size = 7),
        axis.text.y  = ggplot2::element_text(size = 8),
        legend.position = "none"
      )
    
    p_gene <- add_thr_layer(p_gene, dt_gene, "gene_f")
    p_gene_pl <- plotly::ggplotly(p_gene, tooltip = "text")
    
    # 3) Cluster boxplot (ordered)
    cl_ord <- NULL
    
    if (!("cluster_id" %in% names(dt))) {
      
      p_cl_pl <- plotly::ggplotly(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = mean(y_lim), label = "No clusters available", size = 5) +
          ggplot2::xlim(0, 1) +
          ggplot2::coord_cartesian(ylim = y_lim) +
          ggplot2::theme_void(),
        tooltip = NULL
      )
      
    } else {
      
      dt_cl <- dt %>% dplyr::filter(!is.na(cluster_id), nzchar(cluster_id))
      
      if (nrow(dt_cl) == 0) {
        
        p_cl_pl <- plotly::ggplotly(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = mean(y_lim), label = "No points assigned to clusters", size = 5) +
            ggplot2::xlim(0, 1) +
            ggplot2::coord_cartesian(ylim = y_lim) +
            ggplot2::theme_void(),
          tooltip = NULL
        )
        
      } else {
        
        parse_chr_num <- function(x) {
          x <- as.character(x)
          x <- gsub("^chr", "", x, ignore.case = TRUE)
          x <- toupper(x)
          dplyr::case_when(
            x %in% c("X", "23") ~ 23L,
            x %in% c("Y", "24") ~ 24L,
            TRUE ~ suppressWarnings(as.integer(x))
          )
        }
        
        has_chr   <- "CHR" %in% names(dt_cl)
        has_cchn  <- "cluster_chr_n" %in% names(dt_cl)
        has_start <- "cluster_start" %in% names(dt_cl)
        has_end   <- "cluster_end" %in% names(dt_cl)
        
        cl_meta <- dt_cl %>%
          dplyr::mutate(
            CHR_raw = if (has_chr) as.character(.data[["CHR"]]) else NA_character_,
            cluster_chr_n_raw = if (has_cchn) as.character(.data[["cluster_chr_n"]]) else NA_character_,
            cluster_start = if (has_start) suppressWarnings(as.numeric(.data[["cluster_start"]])) else NA_real_,
            cluster_end   = if (has_end)   suppressWarnings(as.numeric(.data[["cluster_end"]]))   else NA_real_
          ) %>%
          dplyr::distinct(cluster_id, CHR_raw, cluster_chr_n_raw, cluster_start, cluster_end)
        
        cl_meta <- cl_meta %>%
          dplyr::mutate(
            id_chr = stringr::str_match(cluster_id, "^chr?([0-9]+|X|Y)")[, 2],
            id_cln = stringr::str_match(cluster_id, "(?:_|cl)([0-9]+)$")[, 2],
            
            chr_num = dplyr::coalesce(parse_chr_num(CHR_raw), parse_chr_num(id_chr)),
            chr_num = dplyr::if_else(is.na(chr_num), 999L, chr_num),
            
            cl_n = suppressWarnings(as.integer(cluster_chr_n_raw)),
            cl_n = dplyr::coalesce(cl_n, suppressWarnings(as.integer(id_cln))),
            cl_n = dplyr::coalesce(cl_n, dplyr::dense_rank(dplyr::coalesce(cluster_start, 0)))
          ) %>%
          dplyr::arrange(chr_num, cl_n, cluster_start, cluster_end, cluster_id)
        
        cl_ord <- cl_meta$cluster_id
        
        dt_cl <- dt_cl %>%
          dplyr::mutate(
            cluster_f = factor(cluster_id, levels = cl_ord),
            pos_id    = as.integer(cluster_f),
            fill_grp  = ifelse(pos_id %% 2 == 0, "even", "odd")
          )
        
        stats_cl <- dt_cl %>%
          dplyr::group_by(cluster_f) %>%
          dplyr::summarise(
            p_wilcox = wilcox_vs_bg(value, dt_cl$value[dt_cl$cluster_f != dplyr::first(cluster_f)]),
            .groups = "drop"
          )
        
        dt_cl <- dt_cl %>% dplyr::left_join(stats_cl, by = "cluster_f")
        
        has_range <- all(c("cluster_start", "cluster_end") %in% names(dt_cl))
        dt_cl <- dt_cl %>%
          dplyr::mutate(
            cluster_lbl = if (has_range) {
              paste0(cluster_id, " (", CHR, ":", cluster_start, "-", cluster_end, ")")
            } else {
              as.character(cluster_id)
            }
          )
        
        p_cl <- ggplot2::ggplot(
          dt_cl,
          ggplot2::aes(
            x = cluster_f, y = value, fill = fill_grp,
            text = paste0(
              "<b>", cluster_lbl, "</b>",
              "<br>value: ", signif(value, 4),
              "<br>Wilcox p: ", signif(p_wilcox, 3)
            )
          )
        ) +
          ggplot2::geom_boxplot(outlier.alpha = 0.25) +
          ggplot2::scale_fill_manual(values = c(even = "darkgreen", odd = "orange"), guide = "none") +
          ggplot2::coord_cartesian(ylim = y_lim) +
          #  ggplot2::labs(x = "Cluster", y = metric, title = "By cluster") +
          ggplot2::labs(x = "Cluster", y = metric) +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::theme(
            plot.title   = ggplot2::element_text(size = 11, face = "bold"),
            axis.title.x = ggplot2::element_text(size = 9),
            axis.title.y = ggplot2::element_text(size = 9),
            axis.text.x  = ggplot2::element_text(angle = 90, size = 7),
            axis.text.y  = ggplot2::element_text(size = 8),
            legend.position = "none"
          )
        
        p_cl <- add_thr_layer(p_cl, dt_cl, "cluster_f")
        
        p_cl_pl <- plotly::ggplotly(p_cl, tooltip = "text") %>%
          plotly::layout(xaxis = list(categoryorder = "array", categoryarray = cl_ord))
      }
    }
    
    global_title <- paste0("Summary stats for metric: ", input$nonsyn_metric)
    
    out <- plotly::subplot(
      p_chr_pl, p_gene_pl, p_cl_pl,
      nrows = 1,
      shareY = TRUE,
      titleX = TRUE,
      margin = 0.03
    ) %>%
      plotly::layout(
        title = list(
          text = paste0("<b>", htmltools::htmlEscape(global_title), "</b>"),
          x = 0.5, xanchor = "center",
          y = 0.995, yanchor = "top",
          font = list(size = 16),   # <-- AQUI: tamany del títol
          pad = list(t = 0, b = 0)          # menys “padding” del title
        ),
        showlegend = FALSE,
        margin = list(t = 65),              # ↓ abans 110/120; això evita “encongir” tant
        
        annotations = list(
          # posa els subtítols més a prop del plot (no tan amunt)
          list(text = "<b>By chromosome</b>", x = 0.16, y = 1.02, xref = "paper", yref = "paper",
               showarrow = FALSE, xanchor = "center", yanchor = "bottom"),
          list(text = "<b>By gene</b>",       x = 0.50, y = 1.02, xref = "paper", yref = "paper",
               showarrow = FALSE, xanchor = "center", yanchor = "bottom"),
          list(text = "<b>By cluster</b>",    x = 0.84, y = 1.02, xref = "paper", yref = "paper",
               showarrow = FALSE, xanchor = "center", yanchor = "bottom")
        )
      )
    
    if (!is.null(cl_ord)) {
      out <- out %>% plotly::layout(
        xaxis3 = list(categoryorder = "array", categoryarray = cl_ord)
      ) 
    }
    
    out %>% plotly::config(displaylogo = FALSE) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "summary_stats",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
    
  })
  
  # ========================================================================
  # 12) Stats tables (chr / gene / cluster) + downloads
  # ========================================================================
  
  nonsyn_chr_stats <- reactive({
    dt <- nonsyn_manhattan_df()
    req(is.data.frame(dt), nrow(dt) > 0)
    
    metric <- input$nonsyn_metric
    bg_values <- dt$value
    
    dt %>%
      dplyr::filter(!is.na(CHR), is.finite(value)) %>%
      dplyr::group_by(CHR) %>%
      dplyr::reframe(
        scope  = "Chromosome",
        id     = paste0("chr", unique(CHR)),
        metric = metric,
        n_snp  = n(),
        median_chr = median(value, na.rm = TRUE),
        median_bg  = median(bg_values, na.rm = TRUE),
        fold_change = median_chr / median_bg,
        wilcox_p = tryCatch(wilcox.test(value, bg_values, exact = FALSE)$p.value, error = function(e) NA_real_),
        ks_p     = tryCatch(ks.test(value, bg_values)$p.value, error = function(e) NA_real_)
      ) %>%
      dplyr::ungroup()
  })
  
  nonsyn_gene_stats <- reactive({
    dt <- nonsyn_manhattan_df()
    req(is.data.frame(dt), nrow(dt) > 0)
    
    metric <- input$nonsyn_metric
    bg_values <- dt$value
    
    dt %>%
      dplyr::filter(!is.na(GENE), is.finite(value)) %>%
      dplyr::group_by(GENE) %>%
      dplyr::reframe(
        scope  = "Gene",
        id     = unique(GENE),
        metric = metric,
        n_snp  = n(),
        median_gene = median(value, na.rm = TRUE),
        median_bg  = median(bg_values, na.rm = TRUE),
        fold_change = median_gene / median_bg,
        wilcox_p = tryCatch(wilcox.test(value, bg_values, exact = FALSE)$p.value, error = function(e) NA_real_),
        ks_p     = tryCatch(ks.test(value, bg_values)$p.value, error = function(e) NA_real_)
      ) %>%
      dplyr::ungroup()
  })
  
  nonsyn_cluster_stats <- reactive({
    dt0 <- nonsyn_manhattan_df()
    req(is.data.frame(dt0), nrow(dt0) > 0)
    req("cluster_id" %in% names(dt0))
    
    metric <- input$nonsyn_metric
    
    # GLOBAL background (abans de filtrar clusters)
    bg_values <- dt0$value
    bg_values <- bg_values[is.finite(bg_values)]
    
    dt <- dt0 %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(CHR), is.finite(value))
    req(nrow(dt) > 0)
    
    out <- dt %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::reframe(
        scope  = "Cluster",
        id     = unique(cluster_id),
        metric = metric,
        n_snp  = dplyr::n(),
        median_cluster = median(value, na.rm = TRUE),
        median_bg      = median(bg_values, na.rm = TRUE),
        fold_change    = median_cluster / median_bg,
        wilcox_p = tryCatch(wilcox.test(value, bg_values, exact = FALSE)$p.value, error = function(e) NA_real_),
        ks_p     = tryCatch(ks.test(value, bg_values)$p.value, error = function(e) NA_real_)
      ) %>%
      dplyr::ungroup()
    
    # Optional: stable ordering by chr + cluster number if available
    if (all(c("cluster_chr_n", "cluster_id") %in% names(dt))) {
      ord <- dt %>%
        dplyr::distinct(cluster_id, CHR, cluster_chr_n) %>%
        dplyr::mutate(
          CHR = suppressWarnings(as.integer(CHR)),
          cluster_chr_n = suppressWarnings(as.integer(cluster_chr_n))
        ) %>%
        dplyr::arrange(CHR, cluster_chr_n) %>%
        dplyr::pull(cluster_id)
      
      out <- out %>%
        dplyr::mutate(id = factor(id, levels = ord)) %>%
        dplyr::arrange(id) %>%
        dplyr::mutate(id = as.character(id))
    }
    
    out
  })
  
  output$nonsyn_chr_stats_tbl <- DT::renderDT({
    df <- nonsyn_chr_stats()
    req(nrow(df) > 0)
    
    DT::datatable(
      df %>%
        dplyr::mutate(
          dplyr::across(c(median_chr, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        ),
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
  }, server = FALSE)
  
  output$nonsyn_gene_stats_tbl <- DT::renderDT({
    df <- nonsyn_gene_stats()
    req(nrow(df) > 0)
    
    DT::datatable(
      df %>%
        dplyr::mutate(
          dplyr::across(c(median_gene, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        ),
      rownames = FALSE,
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
  
  output$nonsyn_cluster_stats_tbl <- DT::renderDT({
    df <- nonsyn_cluster_stats()
    req(nrow(df) > 0)
    
    DT::datatable(
      df %>%
        dplyr::mutate(
          dplyr::across(c(median_cluster, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        ),
      rownames = FALSE,
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
  
  output$dl_nonsyn_chr_stats_csv <- downloadHandler(
    filename = function() {
      tg <- make_mode_thr_tag(input$cluster_method, input$pthr, input$min_logp)
      paste0("nonsyn_chr_stats_", tg$mode_tag, "_thr", tg$thr_txt, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      df <- nonsyn_chr_stats()
      req(is.data.frame(df), nrow(df) > 0)
      
      df_out <- df %>%
        dplyr::mutate(
          dplyr::across(c(median_chr, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        )
      
      utils::write.csv(df_out, file, row.names = FALSE)
    }
  )
  
  output$dl_nonsyn_gene_stats_csv <- downloadHandler(
    filename = function() {
      tg <- make_mode_thr_tag(input$cluster_method, input$pthr, input$min_logp)
      paste0("nonsyn_gene_stats_", tg$mode_tag, "_thr", tg$thr_txt, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      df <- nonsyn_gene_stats()
      req(is.data.frame(df), nrow(df) > 0)
      
      df_out <- df %>%
        dplyr::mutate(
          dplyr::across(c(median_gene, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        )
      
      utils::write.csv(df_out, file, row.names = FALSE)
    }
  )
  
  output$dl_nonsyn_cluster_stats_csv <- downloadHandler(
    filename = function() {
      tg <- make_mode_thr_tag(input$cluster_method, input$pthr, input$min_logp)
      paste0("nonsyn_cluster_stats_", tg$mode_tag, "_thr", tg$thr_txt, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      df <- nonsyn_cluster_stats()
      req(is.data.frame(df), nrow(df) > 0)
      
      df_out <- df %>%
        dplyr::mutate(
          dplyr::across(c(median_cluster, median_bg, fold_change), ~ round(.x, 3)),
          dplyr::across(c(wilcox_p, ks_p), ~ signif(.x, 3))
        )
      
      utils::write.csv(df_out, file, row.names = FALSE)
    }
  )
  
  # =======================================================================
  # Helper predictors in enrich codes 
  # =======================================================================
  pred_item_meaning <- function(column, item) {
    col <- as.character(column); it <- as.character(item)
    
    # PolyPhen
    if (grepl("^Polyphen2_", col, ignore.case = TRUE)) {
      return(switch(it,
                    "B"="Benign",
                    "P"="Possibly damaging",
                    "D"="Probably damaging",
                    NA_character_
      ))
    }
    
    # MutationTaster
    if (grepl("^MutationTaster_", col, ignore.case = TRUE)) {
      return(switch(it,
                    "A"="Disease causing (automatic)",
                    "D"="Disease causing",
                    "N"="Polymorphism",
                    "P"="Polymorphism (automatic)",
                    NA_character_
      ))
    }
    
    # MutationAssessor
    if (grepl("^MutationAssessor_", col, ignore.case = TRUE)) {
      return(switch(it,
                    "H"="High impact",
                    "M"="Medium impact",
                    "L"="Low impact",
                    "N"="Neutral",
                    NA_character_
      ))
    }
    
    # PROVEAN
    if (grepl("^PROVEAN_", col, ignore.case = TRUE)) {
      return(switch(it, "D"="Damaging", "N"="Neutral", NA_character_))
    }
    
    # AlphaMissense
    if (grepl("^AlphaMissense_", col, ignore.case = TRUE)) {
      return(switch(it,
                    "P" ="Pathogenic",
                    "LP"="Likely pathogenic",
                    "A" ="Ambiguous",
                    "B" ="Benign",
                    "LB"="Likely benign",
                    NA_character_
      ))
    }
    
    # Binari per defecte (SIFT/Meta*/ClinPred/etc.)
    if (it %in% c("D","T")) {
      return(ifelse(it=="D","Damaging / deleterious","Tolerated"))
    }
    
    # Aloft (text)
    if (grepl("^Aloft_", col, ignore.case = TRUE)) {
      if (it %in% c("Tolerant","Recessive","Dominant")) return(it)
    }
    
    NA_character_
  }
  # ========================================================================
  # 13) Info modals (kept as-is, only relocated for readability)
  # ========================================================================
  
  observeEvent(input$info_00, {
    showModal(modalDialog(
      title = "Info: Hits to plot",
      HTML('
<p style="text-align: justify;">
  Input file format [(*) mandatory columns]
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
'),
      easyClose = TRUE,
      footer = NULL
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
  
  observeEvent(input$info_10, {
    showModal(modalDialog(
      title = "Non-synonymous variant metrics",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      HTML("... (your existing HTML unchanged) ...")
    ))
  })
  
  observeEvent(input$info_11, {
    txt <- HTML(
      "<b>Radar groups (rankscore) — quick guide</b><br><br>",
      "<ul>",
      "<li><b>Functional_Classic</b><br>Classic missense functional predictors.</li><br>",
      "<li><b>Pathogenicity_Meta</b><br>Meta-scores combining multiple features.</li><br>",
      "<li><b>Pathogenicity_DL</b><br>Deep-learning / language-model scores.</li><br>",
      "<li><b>Constraint_Intolerance</b><br>Population-based constraint evidence.</li><br>",
      "<li><b>Conservation</b><br>Phylogenetic conservation across species.</li><br>",
      "<li><b>Eigen_Unsup</b><br>Common practical clinical scores (overlap intended).</li>",
      "</ul>"
    )
    
    showModal(modalDialog(
      title = "Radar groups (rankscore)",
      txt,
      easyClose = TRUE,
      footer = modalButton("Close")
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
  
  # ========================================================================
  # ========================================================================
  # 14) Functional analysis (GO/KEGG) — left untouched in logic, de-duplicated UI
  # ========================================================================
  
  scope_df <- reactive({
    df <- dbnsfp_norm_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    if (identical(input$func_scope, "cluster")) {
      req(input$func_cluster_id)
      if ("cluster_id" %in% names(df)) {
        df <- df[dplyr::coalesce(as.character(df$cluster_id), "") == as.character(input$func_cluster_id), , drop = FALSE]
      }
    }
    df
  })
  
  empty_sel_df <- function() {
    tibble::tibble(CHR = integer(), POS = integer(), genename = character(), cluster_id = character())
  }
  
  split_semicol_unique <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x) & x != "." & x != "NA"]
    parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[nzchar(parts) & parts != "." & parts != "NA"]
    parts <- parts[grepl("^[A-Za-z0-9_.-]+$", parts)]
    unique(parts)
  }
  
  shorten_go <- function(x, max = 55) {
    x <- as.character(x)
    x <- gsub("^GO:[0-9]+\\s*;?\\s*", "", x)
    x <- gsub("\\bregulation of\\b", "reg. of", x)
    x <- gsub("\\bpositive\\b", "pos.", x)
    x <- gsub("\\bnegative\\b", "neg.", x)
    x <- gsub("\\bcellular\\b",  "cell.", x)
    x <- gsub("\\bprocess\\b",   "proc.", x)
    x <- gsub("\\bcomponent\\b", "comp.", x)
    x <- gsub("\\bactivity\\b",  "act.", x)
    x <- stringr::str_squish(x)
    ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
  }
  
  scope_label <- function() {
    if (identical(input$func_scope, "cluster")) {
      paste0("Cluster ", input$func_cluster_id %||% "")
    } else {
      "Global"
    }
  }
  
  cluster_gene_counts <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "No final output loaded."))
    validate(need(all(c("cluster_id","genename") %in% names(df)), "Final output must contain: cluster_id, genename"))
    
    df %>%
      dplyr::filter(!is.na(cluster_id), nzchar(as.character(cluster_id))) %>%
      dplyr::mutate(genename = as.character(genename)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(n_genes = length(split_semicol_unique(genename)), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n_genes), cluster_id)
  })
  
  # ---------- UI: cluster selector (for functional enrichment) ----------
  
  output$func_cluster_ui <- renderUI({
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    if (!is.data.frame(df) || !nrow(df)) return(NULL)
    
    # cluster col
    cc <- NULL
    for (cand in c("cluster_id","cluster_chr_n","cluster","cluster_chr")) {
      if (cand %in% names(df)) { cc <- cand; break }
    }
    validate(need(!is.null(cc), "No cluster column found in dbNSFP table."))
    
    # gene col (IMPORTANT)
    gene_col <- detect_gene_col(df)
    validate(need(!is.null(gene_col), "No gene column found (needed to count genes per cluster)."))
    
    df[[cc]] <- sub("^cluster_", "", as.character(df[[cc]]))
    
    dt <- data.table::as.data.table(df)
    byc <- dt[, .(n_genes = data.table::uniqueN(get(gene_col))), by = cc]
    data.table::setorder(byc, -n_genes)
    
    min_genes <- 1L  # si vols, ho fem input després
    ok <- byc[n_genes >= min_genes]
    
    validate(need(nrow(ok) > 0, paste0("No clusters with ≥ ", min_genes, " genes are available for enrichment.")))
    
    choices <- stats::setNames(
      ok[[cc]],
      paste0("cluster_", ok[[cc]], " (", ok$n_genes, " genes)")
    )
    
    selectInput("func_cluster_id", "Cluster", choices = choices, selected = ok[[cc]][1])
  })
  
  output$func_gene_ui <- renderUI({
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    if (!is.data.frame(df) || !nrow(df)) return(NULL)
    
    gene_col <- detect_gene_col(df)
    validate(need(!is.null(gene_col), "No gene column found in dbNSFP normalized table."))
    
    gs <- unique(na.omit(as.character(df[[gene_col]])))
    gs <- gs[nzchar(gs)]
    gs <- sort(gs)
    
    selectizeInput("func_gene", "Gene", choices = gs, selected = if (length(gs)) gs[1] else NULL,
                   options = list(placeholder = "Type gene…"))
  })
  
  
  # ---------- Genes (all + scope) ----------
  
  genes_all_symbols <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "No dbNSFP final output loaded."))
    validate(need("genename" %in% names(df), "Column 'genename' not found."))
    split_semicol_unique(df$genename)
  })
  
  genes_scope_symbols <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "No dbNSFP final output loaded."))
    validate(need("genename" %in% names(df), "Column 'genename' not found."))
    
    if (identical(input$func_scope, "cluster")) {
      validate(need("cluster_id" %in% names(df), "Column 'cluster_id' not found."))
      req(input$func_cluster_id)
      df <- df %>% dplyr::filter(.data[["cluster_id"]] == input$func_cluster_id)
    }
    
    split_semicol_unique(df$genename)
  })
  
  
  # ---------- Map ALL genes once (simple cache) ----------
  
  gene_map_cache      <- reactiveVal(NULL)
  gene_map_cache_syms <- reactiveVal(NULL)
  
  gene_map_all <- reactive({
    syms <- genes_all_symbols()
    validate(need(length(syms) > 0, "No valid gene symbols found."))
    
    prev_syms <- gene_map_cache_syms()
    prev_map  <- gene_map_cache()
    
    # Cache hit if the symbol list hasn't changed
    if (is.character(prev_syms) &&
        length(prev_syms) == length(syms) &&
        identical(sort(prev_syms), sort(syms))) {
      return(prev_map)
    }
    
    m <- suppressWarnings(
      clusterProfiler::bitr(
        syms,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Hs.eg.db
      )
    )
    validate(need(is.data.frame(m) && nrow(m) > 0, "Could not map any genes to ENTREZID."))
    
    gene_map_cache_syms(syms)
    gene_map_cache(m)
    m
  })
  
  
  # ---------- Mapping stats (debug/UX) ----------
  
  output$enrich_map_stats <- renderText({
    syms <- tryCatch(genes_all_symbols(), error = function(e) character(0))
    if (!length(syms)) return("Gene mapping: (no genes yet)")
    
    m <- tryCatch(gene_map_all(), error = function(e) NULL)
    if (!is.data.frame(m) || !nrow(m)) {
      return(paste0("Gene mapping: 0/", length(syms), " symbols mapped"))
    }
    
    mapped_syms <- length(unique(m$SYMBOL))
    pct <- round(100 * mapped_syms / max(1, length(syms)), 1)
    
    paste0(
      "Gene mapping\n",
      "- symbols in dataset: ", length(syms), "\n",
      "- mapped symbols:     ", mapped_syms, " (", pct, "%)\n",
      "- mapped pairs:       ", nrow(m)
    )
  })
  
  
  # ---------- Entrez (scope + universes) ----------
  
  entrez_scope <- reactive({
    m <- gene_map_all()
    syms <- genes_scope_symbols()
    validate(need(length(syms) > 0, "No valid genes for the selected scope."))
    
    ids <- unique(m$ENTREZID[m$SYMBOL %in% syms])
    ids <- ids[!is.na(ids) & nzchar(ids)]
    validate(need(length(ids) > 0, "No ENTREZIDs for the selected scope (mapping empty)."))
    ids
  })
  
  universe_entrez_dataset <- reactive({
    m <- gene_map_all()
    u <- unique(m$ENTREZID)
    u <- u[!is.na(u) & nzchar(u)]
    validate(need(length(u) > 0, "Universe (dataset) is empty after mapping."))
    u
  })
  
  universe_entrez_for_scope <- reactive({
    scope <- input$func_scope %||% "global"
    bg    <- input$enrich_background %||% "dataset"
    
    # In GLOBAL mode: use default background (OrgDb/KEGG), i.e. universe = NULL
    if (identical(scope, "global")) return(NULL)
    
    # In CLUSTER mode: user decides background
    if (identical(bg, "dataset")) return(universe_entrez_dataset())
    
    # bg == "orgdb" (or any other value): default OrgDb/KEGG background
    NULL
  })
  
  observeEvent(input$func_scope, {
    if (identical(input$func_scope, "global")) {
      updateRadioButtons(session, "enrich_background", selected = "orgdb")
    }
  }, ignoreInit = TRUE)
  
  
  # ============================================================
  # UI note: background explanation (GO vs KEGG)
  # ============================================================
  
  output$enrich_bg_note <- renderUI({
    bg  <- input$enrich_background
    if (is.null(bg) || !length(bg) || !nzchar(bg)) bg <- "dataset"
    
    tab <- input$enrich_tabs
    if (is.null(tab) || !length(tab) || !nzchar(tab)) tab <- "tab_enrich_terms"
    
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
    # TERMS (NonSyn dbNSFP) helpers
    # -----------------------------
    safe_int <- function(x) {
      x <- suppressWarnings(as.integer(x))
      if (!is.finite(x)) NA_integer_ else x
    }
    
    terms_info <- function() {
      # IMPORTANT: Terms enrichment NO usa input$enrich_background (usa universe_nonsyn_cache)
      fld <- input$nonsyn_term_field %||% ""
      
      # Foreground scope table (rows)
      dt2 <- tryCatch({
        x <- dt_enrich_res_nonsyn()
        x$dt2
      }, error = function(e) NULL)
      
      n_rows <- if (is.data.frame(dt2)) nrow(dt2) else NA_integer_
      
      # Universe table for selected field
      dir_bg <- "/Volumes/DISK1TB/Inspector_app_slaves_ngroc/GItools/_shared/universe_nonsyn_cache"
      tabU <- tryCatch({
        bg_all <- load_universe_nonsyn_cache(dir_bg)
        bg_key <- universe_key_for_field(fld)
        if (is.na(bg_key) || !(bg_key %in% names(bg_all))) return(NULL)
        bg_all[[bg_key]]
      }, error = function(e) NULL)
      
      n_uni_terms <- if (is.data.frame(tabU) && "term_id" %in% names(tabU)) length(unique(tabU$term_id)) else NA_integer_
      N_total <- if (is.data.frame(tabU) && "freq" %in% names(tabU)) suppressWarnings(sum(as.integer(tabU$freq), na.rm = TRUE)) else NA_integer_
      
      # Foreground parsed tokens (occurrences)
      n_tok <- NA_integer_
      n_uniq <- NA_integer_
      n_hit <- NA_integer_
      
      if (is.data.frame(dt2) && nrow(dt2) > 0 && nzchar(fld) && fld %in% names(dt2) && is.data.frame(tabU)) {
        toks <- tryCatch(
          fg_terms_to_term_id(fld, dt2[[fld]], bg_df = tabU),
          error = function(e) character(0)
        )
        toks <- as.character(toks)
        toks <- toks[!is.na(toks) & nzchar(trimws(toks))]
        
        n_tok  <- length(toks)
        n_uniq <- length(unique(toks))
        
        if (is.data.frame(tabU) && "term_id" %in% names(tabU)) {
          u <- unique(as.character(tabU$term_id))
          n_hit <- sum(unique(toks) %in% u)
        }
      }
      
      list(
        fld = fld,
        n_rows = n_rows,
        n_tok = n_tok,
        n_uniq = n_uniq,
        n_hit = n_hit,
        n_uni_terms = n_uni_terms,
        N_total = N_total
      )
    }
    
    # -----------------------------
    # -----------------------------
    predictors_info <- function() {
      
      # Foreground table (rows) en l’scope actual
      dt2 <- tryCatch({
        x <- dt_enrich_res_nonsyn()
        x$dt2
      }, error = function(e) NULL)
      
      n_rows <- if (is.data.frame(dt2)) nrow(dt2) else NA_integer_
      
      # Universe predictors
      univ_path <- tryCatch(PRED_UNIV_RDS, error = function(e) NA_character_)
      uni <- tryCatch({
        if (!is.character(univ_path) || !nzchar(univ_path) || !file.exists(univ_path)) return(NULL)
        readRDS(univ_path)
      }, error = function(e) NULL)
      
      n_uni_rows <- if (is.data.frame(uni)) nrow(uni) else NA_integer_
      n_uni_pred <- if (is.data.frame(uni) && "column" %in% names(uni)) length(unique(as.character(uni$column))) else NA_integer_
      
      # Foreground predictor columns present
      pred_cols_fg <- NA_integer_
      pred_cols_names <- character(0)
      if (is.data.frame(dt2)) {
        pred_cols_names <- names(dt2)[grepl("_pred$", names(dt2))]
        pred_cols_fg <- length(pred_cols_names)
      }
      
      # Parsed tokens (occurrences) for the current scope (rough)
      # (count items across all pred cols present; respects split_pipe)
      split_pipe <- isTRUE(input$pred_split_pipe)
      n_tok <- NA_integer_
      n_uniq <- NA_integer_
      
      if (is.data.frame(dt2) && pred_cols_fg > 0) {
        toks_all <- character(0)
        for (cc in pred_cols_names) {
          v <- as.character(dt2[[cc]])
          v <- v[!is.na(v)]
          v <- trimws(v)
          v <- v[nzchar(v)]
          v <- v[v != "."]
          if (!length(v)) next
          
          toks <- unlist(strsplit(v, ";", fixed = TRUE), use.names = FALSE)
          toks <- trimws(toks)
          toks <- toks[nzchar(toks)]
          toks <- toks[toks != "."]
          if (split_pipe && length(toks)) {
            toks <- unlist(strsplit(toks, "|", fixed = TRUE), use.names = FALSE)
            toks <- trimws(toks)
            toks <- toks[nzchar(toks)]
            toks <- toks[toks != "."]
          }
          toks_all <- c(toks_all, toks)
        }
        
        toks_all <- toks_all[!is.na(toks_all) & nzchar(toks_all)]
        n_tok  <- length(toks_all)
        n_uniq <- length(unique(toks_all))
      }
      
      # Universe total occurrences (sum N)
      N_total <- NA_integer_
      if (is.data.frame(uni) && "N" %in% names(uni)) {
        N_total <- suppressWarnings(sum(as.integer(uni$N), na.rm = TRUE))
      }
      
      list(
        n_rows = n_rows,
        pred_cols_fg = pred_cols_fg,
        n_tok = n_tok,
        n_uniq = n_uniq,
        univ_path = univ_path,
        n_uni_rows = n_uni_rows,
        n_uni_pred = n_uni_pred,
        N_total = N_total,
        split_pipe = split_pipe,
        topk = suppressWarnings(as.integer(input$pred_topk %||% NA_integer_)),
        min_a = suppressWarnings(as.integer(input$pred_min_a %||% NA_integer_))
      )
    }
    
    # -----------------------------
    # Message per tab
    # -----------------------------
    msg <- if (identical(tab, "tab_enrich_terms")) {
      
      ti <- terms_info()
      
      fld_txt <- if (nzchar(ti$fld)) htmltools::htmlEscape(ti$fld) else "(no field selected)"
      
      rows_txt <- if (is.finite(ti$n_rows)) paste0("<b>Foreground:</b> n=", ti$n_rows, " dbNSFP rows in current scope.<br>") else ""
      tok_txt  <- if (is.finite(ti$n_tok))  paste0("<b>Parsed occurrences:</b> n=", ti$n_tok, " (after splitting/cleaning).<br>") else ""
      uniq_txt <- if (is.finite(ti$n_uniq)) paste0("<b>Unique term_ids:</b> n=", ti$n_uniq, ".<br>") else ""
      hit_txt  <- if (is.finite(ti$n_hit))  paste0("<b>Universe matches:</b> n=", ti$n_hit, " unique term_ids found in the universe.<br>") else ""
      
      uni_terms_txt <- if (is.finite(ti$n_uni_terms)) paste0("<b>Background (universe):</b> ", ti$n_uni_terms, " unique terms") else "<b>Background (universe):</b> (unknown)"
      uni_total_txt <- if (is.finite(ti$N_total)) paste0(" (total occurrences N=", ti$N_total, ").") else "."
      
      paste0(
        "<b>NonSyn Terms enrichment</b><br>",
        "<b>Selected term field:</b> <code>", fld_txt, "</code><br>",
        rows_txt, tok_txt, uniq_txt, hit_txt,
        uni_terms_txt, uni_total_txt, "<br>",
        "<b>Test unit:</b> term occurrences (rows) — NOT unique genes.<br>",
        "<b>Method:</b> ORA with hypergeometric test on counts (a vs Tt; BH-FDR).<br>",
        "<b>Notes:</b> This panel ignores the global “Background” selector; it uses <code>universe_nonsyn_cache</code> for the selected field. Use <code>min_a</code> + <code>TopK</code> to control sparsity."
      )
      
    } else if (identical(tab, "tab_enrich_go")) {
      
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
      
    } else if (identical(tab, "tab_enrich_pred")) {
      
      pi <- predictors_info()
      
      rows_txt <- if (is.finite(pi$n_rows)) paste0("<b>Foreground:</b> n=", pi$n_rows, " dbNSFP rows in current scope.<br>") else ""
      cols_txt <- if (is.finite(pi$pred_cols_fg)) paste0("<b>Predictor columns in foreground:</b> n=", pi$pred_cols_fg, " (<code>*_pred</code>).<br>") else ""
      tok_txt  <- if (is.finite(pi$n_tok)) paste0("<b>Parsed occurrences:</b> n=", pi$n_tok, if (isTRUE(pi$split_pipe)) " (splitting <code>;</code> + <code>|</code>).<br>" else " (splitting <code>;</code> only).<br>") else ""
      uniq_txt <- if (is.finite(pi$n_uniq)) paste0("<b>Unique items:</b> n=", pi$n_uniq, ".<br>") else ""
      
      uni_path_txt <- if (is.character(pi$univ_path) && nzchar(pi$univ_path)) paste0("<code>", htmltools::htmlEscape(pi$univ_path), "</code>") else "(unknown path)"
      uni_txt2 <- paste0(
        "<b>Background (universe):</b> predictor item counts precomputed in ", uni_path_txt, ".<br>",
        if (is.finite(pi$n_uni_pred)) paste0("<b>Universe predictors:</b> n=", pi$n_uni_pred, ".<br>") else "",
        if (is.finite(pi$n_uni_rows)) paste0("<b>Universe rows:</b> n=", pi$n_uni_rows, " (predictor×item).<br>") else "",
        if (is.finite(pi$N_total)) paste0("<b>Universe total occurrences:</b> N=", pi$N_total, ".<br>") else ""
      )
      
      min_a_txt <- if (is.finite(pi$min_a)) paste0("<b>min_a:</b> ", pi$min_a, " (min occurrences in foreground per predictor×item).<br>") else ""
      topk_txt  <- if (is.finite(pi$topk))  paste0("<b>TopK:</b> ", pi$topk, " (preselect top by foreground count before ORA).<br>") else ""
      
      paste0(
        "<b>Predictor enrichment</b><br>",
        sc_txt, "<br>",
        rows_txt, cols_txt, tok_txt, uniq_txt,
        uni_txt2,
        "<b>Test unit:</b> predictor-item occurrences (rows/tokens) — NOT unique genes.<br>",
        "<b>Method:</b> ORA (Fisher exact / hypergeometric-style) per predictor×item with BH-FDR.<br>",
        min_a_txt, topk_txt,
        "<b>Notes:</b> This panel ignores the global “Background” selector (dataset/orgdb); it always uses the predictor universe RDS."
      ) 
    } else {
      paste0("<b>Background:</b> ", bg_lbl)
    }
    
    htmltools::HTML(paste0(
      "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
      msg, "</div>"
    ))
  })
  
  
  # ============================================================
  # EVENT: run enrichment (separate GO vs KEGG triggers)
  # ============================================================
  
  go_trigger     <- reactiveVal(0L)
  kegg_trigger   <- reactiveVal(0L)
  goslim_trigger <- reactiveVal(0L)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    if (identical(tab, "tab_enrich_go")) {
      go_trigger(go_trigger() + 1L)
    } else if (identical(tab, "tab_enrich_kegg")) {
      kegg_trigger(kegg_trigger() + 1L)
    } else if (identical(tab, "tab_enrich_goslim")) {
      goslim_trigger(goslim_trigger() + 1L)
    }
  }, ignoreInit = TRUE)
  
  
  
  # ---------- GO enrichment ----------
  
  go_enrich_raw <- eventReactive(go_trigger(), {
    
    withProgress(message = "Running GO enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni  <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (OrgDb)" else paste("Universe:", length(uni)))
      
      ontos <- input$go_ontos %||% c("BP", "CC", "MF")
      ontos <- intersect(c("BP","CC","MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one GO ontology (BP/CC/MF)."))
      
      minGS <- input$enrich_min_gs %||% 10
      maxGS <- input$enrich_max_gs %||% 500
      
      res_list <- lapply(seq_along(ontos), function(i) {
        ont <- ontos[i]
        incProgress(0.6/length(ontos), detail = paste("Ontology:", ont))
        
        eg <- clusterProfiler::enrichGO(
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
        
        if (is.null(eg) || is.null(eg@result) || !nrow(eg@result)) return(NULL)
        df <- as.data.frame(eg@result, stringsAsFactors = FALSE)
        df$Ontology <- ont
        df
      })
      
      out <- dplyr::bind_rows(res_list)
      validate(need(is.data.frame(out) && nrow(out) > 0, "No GO enrichment results."))
      
      incProgress(0.1, detail = "Done")
      out
    })
    
  }, ignoreInit = TRUE)
  
  
  
  go_top_df <- reactive({
    df <- go_enrich_raw()
    req(is.data.frame(df), nrow(df) > 0)
    
    pcut <- input$go_kegg_pcut %||% 0.05
    topn <- 10
    
    df_sig <- df %>%
      dplyr::filter(!is.na(p.adjust), p.adjust <= pcut)
    
    # If nothing passes FDR, fallback to top terms by FDR anyway
    using_fallback <- FALSE
    if (!nrow(df_sig)) {
      using_fallback <- TRUE
      df_sig <- df
    }
    
    out <- df_sig %>%
      dplyr::group_by(Ontology) %>%
      dplyr::slice_min(order_by = p.adjust, n = topn, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        term_short = stringr::str_wrap(shorten_go(Description, 55), 28),
        score      = -log10(pvalue + 1e-300),  # keep bars visible even if FDR ~ 1
        FDR        = p.adjust
      )
    
    attr(out, "fallback") <- using_fallback
    out
  })
  
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
    
    # fmt_p fallback (si no existeix a NonSyn)
    if (!exists("fmt_p", mode = "function")) {
      fmt_p <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
      }
    }
    
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO term:</b> ", htmltools::htmlEscape(dat$go_term),
      "<br><b>Genes:</b> ", dat$n_genes,
      "<br><b>p.adj:</b> ", fmt_p(dat$p.adjust),
      "<br><b>p:</b> ", fmt_p(dat$pvalue)
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
        plot.title    = ggplot2::element_text(size=11, face="bold"),
        plot.subtitle = ggplot2::element_text(size=9),
        axis.text.x   = ggplot2::element_text(angle=90, size=7, vjust=0.5, hjust=1),
        legend.position="none",
        plot.margin   = grid::unit(c(28,10,35,10), "pt")
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
  
  output$go_table <- DT::renderDT({
    df0 <- go_enrich_raw()
    req(is.data.frame(df0), nrow(df0) > 0)
    
    pcut <- input$go_kegg_pcut %||% 0.05
    
    df <- df0 %>%
      dplyr::mutate(
        FDR      = p.adjust,
        pass_FDR = dplyr::if_else(!is.na(p.adjust) & p.adjust <= pcut, "✓", "")
      ) %>%
      dplyr::arrange(p.adjust, pvalue) %>%
      dplyr::transmute(
        Ontology,
        pass_FDR,
        GOID        = ID,
        Description,
        Count,
        GeneRatio,
        BgRatio,
        pvalue,
        FDR,
        qvalue,
        Genes = geneID
      )
    
    DT::datatable(
      df,
      rownames   = FALSE,
      filter     = "top",
      extensions = c("Buttons"),
      options    = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        # FDR is column 9 in the displayed table (1-based index inside DataTables)
        order      = list(list(8, "asc"))
      )
    ) %>%
      DT::formatSignif(c("pvalue", "FDR", "qvalue"), digits = 3)
  }, server = FALSE)
  
  
  # ==========================
  # KEGG enrichment (robust)
  #   - tries: ncbi-geneid -> kegg (hsa:ID) -> internal KEGG.db
  #   - never hard-stops; returns tibble() when empty
  # ==========================
  
  run_kegg_try <- function(gene, universe = NULL, use_internal = FALSE, mode = c("ncbi", "kegg")) {
    mode <- match.arg(mode)
    
    # clusterProfiler expects character vectors
    gene <- unique(na.omit(as.character(gene)))
    gene <- gene[nzchar(gene)]
    if (!length(gene)) return(NULL)
    
    uni <- universe
    if (!is.null(uni)) {
      uni <- unique(na.omit(as.character(uni)))
      uni <- uni[nzchar(uni)]
      if (!length(uni)) uni <- NULL
    }
    
    # if mode = kegg -> prefix hsa:
    if (identical(mode, "kegg")) {
      gene <- ifelse(grepl("^hsa:", gene), gene, paste0("hsa:", gene))
      if (!is.null(uni)) uni <- ifelse(grepl("^hsa:", uni), uni, paste0("hsa:", uni))
    }
    
    minGS <- input$enrich_min_gs %||% 10
    maxGS <- input$enrich_max_gs %||% 500
    
    ek <- tryCatch({
      suppressMessages(clusterProfiler::enrichKEGG(
        gene              = gene,
        universe          = uni,
        organism          = "hsa",
        keyType           = if (identical(mode, "ncbi")) "ncbi-geneid" else "kegg",
        pvalueCutoff      = 1,
        pAdjustMethod     = "BH",
        qvalueCutoff      = 1,
        minGSSize         = minGS,
        maxGSSize         = maxGS,
        use_internal_data = use_internal
      ))
    }, error = function(e) NULL)
    
    if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(NULL)
    ek
  }
  
  kegg_raw <- eventReactive(kegg_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
    withProgress(message = "Running KEGG enrichment…", value = 0, {
      
      gene <- tryCatch(entrez_scope(), error = function(e) character(0))
      gene <- unique(na.omit(as.character(gene)))
      gene <- gene[nzchar(gene)]
      
      incProgress(0.25, detail = paste("Genes:", length(gene)))
      if (!length(gene)) return(NULL)
      
      uni <- tryCatch(universe_entrez_for_scope(), error = function(e) NULL)
      incProgress(0.10, detail = if (is.null(uni)) "Universe: default" else paste("Universe:", length(uni)))
      
      # 1) Standard: ncbi-geneid (online)
      ek <- run_kegg_try(gene, uni, use_internal = FALSE, mode = "ncbi")
      
      # 2) Retry using KEGG IDs format (hsa:xxxx)
      if (is.null(ek)) {
        incProgress(0.10, detail = "Retry: KEGG ID format (hsa:...)")
        ek <- run_kegg_try(gene, uni, use_internal = FALSE, mode = "kegg")
      }
      
      # 3) Retry using internal KEGG.db if available
      if (is.null(ek) && requireNamespace("KEGG.db", quietly = TRUE)) {
        incProgress(0.10, detail = "Retry: internal KEGG.db")
        ek <- run_kegg_try(gene, uni, use_internal = TRUE, mode = "ncbi")
        if (is.null(ek)) ek <- run_kegg_try(gene, uni, use_internal = TRUE, mode = "kegg")
      }
      
      if (is.null(ek)) return(NULL)
      
      incProgress(0.35, detail = "Setting readable…")
      suppressWarnings({
        ek <- clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      })
      
      incProgress(0.20, detail = "Done")
      ek
    })
  }, ignoreInit = TRUE)
  
  kegg_df_raw <- reactive({
    ek <- kegg_raw()
    if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(tibble::tibble())
    as.data.frame(ek@result, stringsAsFactors = FALSE)
  })
  
  kegg_top_df <- reactive({
    df <- kegg_df_raw()
    if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    
    # IMPORTANT: use your unified UI inputs
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$enrich_kegg_top %||% 15
    
    # ensure p.adjust
    if (!"p.adjust" %in% names(df) && "pvalue" %in% names(df)) {
      df$p.adjust <- p.adjust(df$pvalue, method = "BH")
    }
    
    df <- df %>% dplyr::mutate(
      pvalue   = suppressWarnings(as.numeric(pvalue)),
      p.adjust = suppressWarnings(as.numeric(p.adjust))
    ) %>% dplyr::filter(is.finite(pvalue), is.finite(p.adjust))
    
    df_fdr <- df %>% dplyr::filter(p.adjust <= pcut) %>% dplyr::arrange(p.adjust)
    
    if (!nrow(df_fdr)) {
      df_use <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = topn)
      attr(df_use, "fallback_mode") <- paste0("No KEGG with FDR ≤ ", pcut, " → showing top ", topn, " by FDR")
      return(df_use)
    }
    
    df_use <- df_fdr %>% dplyr::slice_head(n = topn)
    attr(df_use, "fallback_mode") <- ""
    df_use
  })
  
  output$kegg_bar <- plotly::renderPlotly({
    df0 <- kegg_top_df()
    
    if (!is.data.frame(df0) || !nrow(df0)) {
      p0 <- ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::labs(
          title = "KEGG enrichment",
          subtitle = "No KEGG pathways returned (mapping failed or none annotated)."
        )
      return(plotly::ggplotly(p0))
    }
    
    note <- attr(df0, "fallback_mode") %||% ""
    ttl  <- paste0(
      "KEGG enrichment (", scope_label(), ")",
      if (nzchar(note)) paste0(" — ", note) else paste0(" — FDR ≤ ", (input$enrich_pcut %||% 0.05))
    )
    
    # robust columns
    tcol <- pick_col(df0, c("Description","term","description","ID","Pathway","term_full","TERM"))
    if (is.null(tcol)) tcol <- "Description"
    
    ccol <- pick_col(df0, c("Count","count","N","n_genes","genes","n"))
    validate(need(!is.null(ccol), "KEGG: missing Count column (expected 'Count')."))
    
    fcol <- pick_col(df0, c("p.adjust","p_adj","FDR","padj","p_adjust"))
    if (is.null(fcol)) fcol <- NA_character_
    
    pcol <- pick_col(df0, c("pvalue","p.value","p_value","PValue","p"))
    if (is.null(pcol)) pcol <- NA_character_
    
    df <- df0 %>%
      dplyr::mutate(
        term_full = as.character(.data[[tcol]]),
        counts    = suppressWarnings(as.integer(.data[[ccol]])),
        pval      = if (!is.na(pcol)) suppressWarnings(as.numeric(.data[[pcol]])) else NA_real_,
        fdr       = if (!is.na(fcol)) suppressWarnings(as.numeric(.data[[fcol]])) else NA_real_,
        logFDR_raw = suppressWarnings(-log10(fdr)),
        is_zero    = is.finite(fdr) & (fdr == 0),
        logFDR     = dplyr::if_else(
          is.finite(logFDR_raw),
          logFDR_raw,
          max(logFDR_raw[is.finite(logFDR_raw)], na.rm = TRUE) + 1
        ),
        tooltip = paste0(
          "<b>", htmltools::htmlEscape(term_full), "</b>",
          "<br><b>Counts:</b> ", counts,
          dplyr::if_else(
            is.finite(pval),
            paste0("<br><b>p:</b> ", formatC(pval, format="e", digits=2)),
            ""
          ),
          dplyr::if_else(
            is.finite(fdr),
            paste0("<br><b>FDR:</b> ", formatC(fdr, format="e", digits=2)),
            "<br><b>FDR:</b> NA"
          ),
          "<br><b>-log10(FDR):</b> ",
          dplyr::if_else(is.finite(logFDR_raw), sprintf("%.2f", logFDR_raw), "Inf/NA")
        )
      ) %>%
      dplyr::filter(!is.na(counts), counts > 0)
    
    validate(need(nrow(df) > 0, "KEGG: no pathways with Count > 0."))
    
    # order: bigger counts first; tie-break by FDR if available
    df <- df %>%
      dplyr::arrange(dplyr::desc(counts), fdr) %>%
      dplyr::mutate(term_full = factor(term_full, levels = rev(unique(term_full)))) # for coord_flip ordering
    
    # dynamic height
    n <- nrow(df)
    h <- max(420, 28 * n)
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = term_full,
        y = counts,
        text = tooltip
      )
    ) +
      # non-zero gradient
      ggplot2::geom_col(
        data = dplyr::filter(df, !is_zero),
        ggplot2::aes(fill = logFDR),
        color = "black", linewidth = 0.25
      ) +
      ggplot2::scale_fill_gradientn(
        colours = c("red", "orange", "yellow"),
        name = "-log10(FDR)"
      ) +
      # FDR==0 -> vermilion
      ggplot2::geom_col(
        data = dplyr::filter(df, is_zero),
        fill = "#E34234",
        color = "black", linewidth = 0.25
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = ttl, x = NULL, y = "NonSyn hit counts by term") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        height = h,
        margin = list(l = 280, r = 40, t = 60, b = 40),
        font   = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "KEGG_barplot_counts",
          width = 1600,
          height = max(900, h),
          scale = 2
        )
      )
  })
  
  output$kegg_table <- DT::renderDT({
    df <- kegg_top_df()
    
    if (!is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(
        data.frame(Message = "No KEGG pathways returned (mapping failed or none annotated)."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    out <- df %>%
      dplyr::mutate(
        KEGG = if ("ID" %in% names(df)) sprintf(
          "<a href='https://www.kegg.jp/pathway/%s' target='_blank'>%s</a>", ID, ID
        ) else NA_character_,
        FDR  = if ("p.adjust" %in% names(df)) p.adjust else NA_real_
      ) %>%
      dplyr::transmute(
        KEGG,
        Description,
        GeneRatio,
        BgRatio,
        Count,
        pvalue,
        FDR,
        qvalue,
        Genes = geneID
      )
    
    DT::datatable(
      out,
      escape     = FALSE,
      rownames   = FALSE,
      filter     = "top",
      extensions = c("Buttons"),
      options    = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        )
      )
    ) %>%
      DT::formatSignif(c("pvalue", "FDR", "qvalue"), digits = 3)
  }, server = FALSE)
  
  ###############################################################################
  # GO selection helpers (chr / cluster / gene) - single source of truth via go_base_df()
  ###############################################################################
  
  # --- Base df used by the GO UI controls ---
  go_base_df <- reactive({
    df <- dbnsfp_norm_df()
    if (!is.data.frame(df) || nrow(df) == 0) return(tibble::tibble())
    
    req(all(c("chr","BP","genename","cluster_id") %in% names(df)))
    
    df %>%
      dplyr::mutate(
        CHR       = suppressWarnings(as.integer(norm_chr_generic(.data[["chr"]]))),
        POS       = suppressWarnings(as.integer(.data[["BP"]])),
        genename  = as.character(.data[["genename"]]),
        cluster_id = as.character(.data[["cluster_id"]])
      ) %>%
      dplyr::filter(!is.na(CHR), !is.na(POS))
  })
  
  # A) Fill chromosomes (go_chr / go_chr_for_cluster)
  observeEvent(go_base_df(), {
    df <- go_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0 || !"CHR" %in% names(df)) {
      updateSelectInput(session, "go_chr", choices = character(0), selected = NULL)
      updateSelectInput(session, "go_chr_for_cluster", choices = character(0), selected = NULL)
      return()
    }
    
    chr_choices <- sort(unique(df$CHR))
    chr_choices <- chr_choices[is.finite(chr_choices)]
    if (!length(chr_choices)) {
      updateSelectInput(session, "go_chr", choices = character(0), selected = NULL)
      updateSelectInput(session, "go_chr_for_cluster", choices = character(0), selected = NULL)
      return()
    }
    
    chr_labels <- paste0("chr", chr_label_plink(chr_choices))
    chr_map    <- stats::setNames(as.character(chr_choices), chr_labels)
    
    sel1 <- if (!is.null(input$go_chr) && input$go_chr %in% chr_map) input$go_chr else chr_map[[1]]
    sel2 <- if (!is.null(input$go_chr_for_cluster) && input$go_chr_for_cluster %in% chr_map) input$go_chr_for_cluster else chr_map[[1]]
    
    updateSelectInput(session, "go_chr", choices = chr_map, selected = sel1)
    updateSelectInput(session, "go_chr_for_cluster", choices = chr_map, selected = sel2)
  }, ignoreInit = FALSE)
  
  # B) Fill clusters by chromosome
  observeEvent(list(go_base_df(), input$go_chr_for_cluster), {
    df <- go_base_df()
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "go_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    # input can be NULL at init -> integer(0)
    chr_raw <- input$go_chr_for_cluster
    if (is.null(chr_raw) || length(chr_raw) == 0 || is.na(chr_raw) || !nzchar(as.character(chr_raw))) {
      updateSelectInput(session, "go_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    chr_sel <- suppressWarnings(as.integer(chr_raw))
    if (length(chr_sel) != 1L || is.na(chr_sel) || !is.finite(chr_sel)) {
      updateSelectInput(session, "go_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    ids <- df %>%
      dplyr::filter(CHR == chr_sel, !is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::arrange(cluster_id) %>%
      dplyr::pull(cluster_id)
    
    # If no clusters, keep empty and selected NULL
    if (!length(ids)) {
      updateSelectInput(session, "go_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    updateSelectInput(session, "go_cluster_id", choices = ids, selected = ids[[1]])
  }, ignoreInit = FALSE)
  
  
  # C) Fill genes grouped by cluster (optgroup)
  observeEvent(go_base_df(), {
    df <- go_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "go_gene", choices = character(0), selected = NULL)
      return()
    }
    
    req(all(c("cluster_id","genename") %in% names(df)))
    
    dd <- df %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id),
                    !is.na(genename),  nzchar(genename)) %>%
      dplyr::distinct(cluster_id, genename) %>%
      dplyr::arrange(cluster_id, genename)
    
    if (!nrow(dd)) {
      updateSelectInput(session, "go_gene", choices = character(0), selected = NULL)
      return()
    }
    
    choices_list <- dd %>%
      dplyr::group_split(cluster_id, .keep = TRUE) %>%
      purrr::set_names(purrr::map_chr(., ~ .x$cluster_id[1])) %>%
      purrr::map(function(x) {
        keys <- paste0(x$cluster_id, "||", x$genename)
        stats::setNames(keys, x$genename)
      })
    
    all_keys <- unlist(choices_list, use.names = FALSE)
    sel <- if (!is.null(input$go_gene) && input$go_gene %in% all_keys) input$go_gene else all_keys[1]
    
    updateSelectInput(session, "go_gene", choices = choices_list, selected = sel)
  }, ignoreInit = FALSE)
  
  
  # --- GO selection dataframe (total/chr/cluster) ---
  go_selected_df <- reactive({
    df <- dbnsfp_norm_df()
    if (!is.data.frame(df) || nrow(df) == 0) {
      out <- tibble::tibble()
      attr(out, "label") <- "No data"
      return(out)
    }
    
    mode <- input$go_mode %||% "total"
    
    if (mode == "total") {
      df_use <- df
      label  <- "Total"
      
    } else if (mode == "chr") {
      chr_sel <- suppressWarnings(as.integer(input$go_chr))
      if (!is.finite(chr_sel)) {
        out <- tibble::tibble()
        attr(out, "label") <- "Chromosome (not selected)"
        return(out)
      }
      
      df_use <- df %>%
        dplyr::mutate(.CHR = suppressWarnings(as.integer(norm_chr_generic(.data[["chr"]])))) %>%
        dplyr::filter(.CHR == chr_sel)
      
      label <- paste0("chr", chr_label_plink(chr_sel))
      
    } else { # cluster
      cl_id <- as.character(input$go_cluster_id %||% "")
      if (!nzchar(cl_id) || !"cluster_id" %in% names(df)) {
        out <- tibble::tibble()
        attr(out, "label") <- "Cluster (not selected)"
        return(out)
      }
      
      df_use <- df %>% dplyr::filter(cluster_id == cl_id)
      label  <- cl_id
    }
    
    if (!is.data.frame(df_use) || nrow(df_use) == 0) {
      out <- tibble::tibble()
      attr(out, "label") <- label
      return(out)
    }
    
    attr(df_use, "label") <- label
    df_use
  })
  
  # ==========================================================
  # GO slim classification (BP, CC, MF) GOSLIM_bar
  # ==========================================================
  
  go_slim_class_from_df <- function(df) {
    
    validate(need(is.data.frame(df) && nrow(df) > 0, "GO: no data for this selection."))
    
    go_cols <- intersect(
      c("GO_biological_process", "GO_cellular_component", "GO_molecular_function"),
      names(df)
    )
    if (length(go_cols) == 0) go_cols <- grep("^GO_", names(df), value = TRUE)
    
    validate(need(length(go_cols) > 0, "GO columns not found in dbNSFP final output."))
    
    long <- df %>%
      dplyr::select(dplyr::any_of(c("genename", "cluster_id", go_cols))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(go_cols), names_to = "go_col", values_to = "term_raw") %>%
      dplyr::filter(!is.na(term_raw), nzchar(term_raw)) %>%
      tidyr::separate_rows(term_raw, sep = ";", convert = FALSE) %>%
      dplyr::mutate(
        term_raw   = stringr::str_squish(term_raw),
        Ontology   = dplyr::recode(
          go_col,
          "GO_biological_process" = "BP",
          "GO_cellular_component" = "CC",
          "GO_molecular_function" = "MF",
          .default = NA_character_
        ),
        key_norm   = norm_go_term(term_raw),
        genename   = as.character(.data[["genename"]]),
        cluster_id = as.character(.data[["cluster_id"]])
      ) %>%
      dplyr::filter(!is.na(Ontology), nzchar(key_norm))
    
    validate(need(nrow(long) > 0, "No GO terms after parsing."))
    
    mapped_go <- long %>%
      dplyr::left_join(
        GO_TERM_LUT %>% dplyr::select(ONTOLOGY, key_norm, GOID),
        by = c("Ontology" = "ONTOLOGY", "key_norm" = "key_norm")
      ) %>%
      dplyr::filter(!is.na(GOID)) %>%
      dplyr::distinct(cluster_id, genename, Ontology, GOID)
    
    validate(need(nrow(mapped_go) > 0, "Could not map GO terms to GO IDs."))
    validate(need(exists("GO_SLIM_GENERIC"), "GO_SLIM_GENERIC not loaded."))
    validate(need(requireNamespace("GO.db", quietly = TRUE), "Package GO.db is required for GO ancestors."))
    
    slim_ids <- split(GO_SLIM_GENERIC$GOID, GO_SLIM_GENERIC$ONTOLOGY)
    
    get_ancestors <- function(goid, ont) {
      anc <- NULL
      if (identical(ont, "BP")) anc <- GO.db::GOBPANCESTOR[[goid]]
      if (identical(ont, "CC")) anc <- GO.db::GOCCANCESTOR[[goid]]
      if (identical(ont, "MF")) anc <- GO.db::GOMFANCESTOR[[goid]]
      unique(c(goid, anc))
    }
    
    uniq <- mapped_go %>% dplyr::distinct(Ontology, GOID)
    uniq$anc <- Map(get_ancestors, uniq$GOID, uniq$Ontology)
    uniq$slim_goid <- Map(
      function(anc, ont) intersect(anc, slim_ids[[ont]] %||% character(0)),
      uniq$anc, uniq$Ontology
    )
    
    uniq_slim <- uniq %>%
      dplyr::select(Ontology, GOID, slim_goid) %>%
      tidyr::unnest(cols = c(slim_goid))
    
    validate(need(nrow(uniq_slim) > 0, "No GO-slim matches found."))
    
    uniq_slim2 <- uniq_slim %>%
      dplyr::left_join(
        GO_SLIM_GENERIC %>% dplyr::select(GOID, ONTOLOGY, slim_term),
        by = c("slim_goid" = "GOID", "Ontology" = "ONTOLOGY")
      ) %>%
      dplyr::filter(!is.na(slim_term)) %>%
      dplyr::distinct(Ontology, GOID, slim_term)
    
    mapped_slim <- mapped_go %>%
      dplyr::inner_join(uniq_slim2, by = c("Ontology", "GOID"), relationship = "many-to-many") %>%
      dplyr::distinct(cluster_id, genename, Ontology, GOID, slim_term)
    
    out <- mapped_slim %>%
      dplyr::group_by(Ontology, slim_term) %>%
      dplyr::summarise(
        n_genes = dplyr::n_distinct(genename),
        n_terms = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(Ontology, dplyr::desc(n_genes))
    
    out
  }
  
  go_slim_class_df_selected <- reactive({
    df_sel <- go_selected_df()
    if (!is.data.frame(df_sel) || nrow(df_sel) == 0) return(tibble::tibble())
    go_slim_class_from_df(df_sel)
  })
  
  
  # ----------------------------
  # GO selection -> genes vector (used by downstream filters)
  # ----------------------------
  
  go_selected_genes <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "Final output not available yet."))
    
    mode <- input$go_mode %||% "total"
    
    if (mode == "cluster") {
      req(input$go_cluster_id)
      genes <- df %>% dplyr::filter(cluster_id == input$go_cluster_id) %>% dplyr::pull(genename)
      
    } else if (mode == "chr") {
      req(input$go_chr)
      chr_sel <- suppressWarnings(as.integer(input$go_chr))
      genes <- df %>%
        dplyr::mutate(.CHR = norm_chr_generic(.data[["chr"]])) %>%
        dplyr::filter(.CHR == chr_sel) %>%
        dplyr::pull(genename)
      
    } else if (mode == "gene") {
      req(input$go_gene)
      parts <- strsplit(as.character(input$go_gene), "\\|\\|")[[1]]
      cluster_sel <- parts[1] %||% ""
      gene_sel    <- parts[2] %||% ""
      req(nzchar(cluster_sel), nzchar(gene_sel))
      
      genes <- df %>%
        dplyr::filter(cluster_id == cluster_sel, genename == gene_sel) %>%
        dplyr::pull(genename)
      
    } else {
      genes <- df$genename
    }
    
    genes <- unique(na.omit(trimws(as.character(genes))))
    genes <- genes[nzchar(genes)]
    genes
  })
  
  go_filter_genes <- reactiveVal(NULL)
  
  observeEvent(go_selected_genes(), {
    gs <- go_selected_genes()
    gs <- unique(na.omit(trimws(as.character(gs))))
    gs <- gs[nzchar(gs)]
    go_filter_genes(if (length(gs)) gs else NULL)
  }, ignoreInit = FALSE)
  
  output$goslim_bar <- renderPlotly({
    
    dat <- goslim_enrich_tbl()
    validate(need(is.data.frame(dat) && nrow(dat) > 0, "GO slim: no enrichment data to plot yet."))
    
    top_n <- input$go_class_top %||% 15
    gap   <- 2
    
    validate(need(all(c("Description","Count","Ontology") %in% names(dat)),
                  "goslim_enrich_tbl() must return: Description, Count, Ontology"))
    
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
        slim_term  = as.character(Description),
        n_genes    = suppressWarnings(as.integer(Count)),
        term_short = stringr::str_wrap(shorten_term(Description, 55), width = 28)
      ) %>%
      dplyr::filter(!is.na(n_genes), nzchar(slim_term))
    
    ont_order <- c("BP","CC","MF")
    present   <- intersect(ont_order, unique(dat$Ontology))
    validate(need(length(present) > 0, "No BP/CC/MF data to plot."))
    
    dat <- dat %>%
      dplyr::filter(Ontology %in% present) %>%
      dplyr::mutate(Ontology = factor(Ontology, levels = ont_order))
    
    dat <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::arrange(p.adjust, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n) %>%
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
    
    sel_label <- tryCatch(scope_label(), error = function(e) "")
    
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO slim term:</b> ", htmltools::htmlEscape(as.character(dat$slim_term)),
      "<br><b>Count:</b> ", dat$n_genes,
      "<br><b>GeneRatio:</b> ", dat$GeneRatio,
      "<br><b>BgRatio:</b> ", dat$BgRatio,
      "<br><b>p:</b> ", sprintf("%.3g", as.numeric(dat$pvalue)),
      "<br><b>FDR:</b> ", sprintf("%.3g", as.numeric(dat$p.adjust))
    )
    
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = x, y = n_genes,
        fill = Ontology,
        text = tooltip
      )
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
        y = "Enriched genes per GO slim term",
        title = paste0("GO slim enrichment by gene counts — ", sel_label)
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
          filename = "GOSlim_enrichment_plot",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  output$goslim_barXXXXX <- renderPlotly({
    
    dat <- goslim_enrich_tbl()
    validate(need(is.data.frame(dat) && nrow(dat) > 0, "GO: no data to plot yet."))
    
    top_n <- input$go_class_top %||% 15
    gap   <- 2
    
    validate(need(all(c("slim_term","n_genes","Ontology") %in% names(dat)),
                  "go_slim_class_df_selected() must return: slim_term, n_genes, Ontology"))
    
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
        Ontology   = as.character(Ontology),
        term_short = stringr::str_wrap(shorten_term(slim_term, 55), width = 28)
      )
    
    # Enforce BP/CC/MF order
    ont_order <- c("BP","CC","MF")
    present   <- intersect(ont_order, unique(dat$Ontology))
    validate(need(length(present) > 0, "No BP/CC/MF data to plot."))
    
    dat <- dat %>%
      dplyr::filter(Ontology %in% present) %>%
      dplyr::mutate(Ontology = factor(Ontology, levels = ont_order))
    
    # Top-N per ontology
    dat <- dat %>%
      dplyr::group_by(Ontology) %>%
      dplyr::slice_max(order_by = n_genes, n = top_n, with_ties = FALSE) %>%
      dplyr::arrange(Ontology, dplyr::desc(n_genes)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    # Offsets per block
    n_bp <- if ("BP" %in% present) max(dat$rank[dat$Ontology=="BP"], 0) else 0
    n_cc <- if ("CC" %in% present) max(dat$rank[dat$Ontology=="CC"], 0) else 0
    
    offsets <- c(
      BP = 0,
      CC = n_bp + gap,
      MF = n_bp + gap + n_cc + gap
    )
    
    dat$x <- dat$rank + offsets[as.character(dat$Ontology)]
    
    # Separator lines between blocks
    vlines <- c()
    if (all(c("BP","CC") %in% present)) vlines <- c(vlines, offsets["CC"] - gap/2)
    if (all(c("CC","MF") %in% present)) vlines <- c(vlines, offsets["MF"] - gap/2)
    
    # Centers for top labels
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
    
    # Selection label (Total/chr/cluster)
    sel_label <- attr(go_selected_df(), "label") %||% ""
    
    # Tooltip summary
    # Tooltip summary
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO slim term:</b> ", htmltools::htmlEscape(as.character(dat$slim_term)),
      "<br><b>Selected genes in term:</b> ", dat$n_genes
    )
    
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = x, y = n_genes,
        fill = Ontology,
        text = tooltip
      )
    ) +
      ggplot2::geom_col(width = 0.85, color = "black", linewidth = 0.25) +
      { if (length(vlines)) ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4) } +
      ggplot2::scale_x_continuous(
        breaks = dat$x,
        labels = dat$term_short,
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = c(ylab, ymax * 1.20), clip = "off") +
      
      # Top class labels
      { if ("BP" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology=="BP"], y = y_top,
                          label = "Biological Process", fontface = "bold", size = 3.6) } +
      { if ("CC" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology=="CC"], y = y_top,
                          label = "Cellular Component", fontface = "bold", size = 3.6) } +
      { if ("MF" %in% present)
        ggplot2::annotate("text",
                          x = centers$xmid[centers$Ontology=="MF"], y = y_top,
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
          filename = "GOSlim_plot",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  output$goslim_barXXX <- renderPlotly({
    
    df <- goslim_enrich_tbl()
    validate(need(is.data.frame(df) && nrow(df) > 0, "GO slim: no enrichment data to plot yet."))
    
    top_n <- input$go_topn %||% 10
    
    df <- df %>%
      dplyr::mutate(
        Ontology    = as.character(Ontology),
        Description = as.character(Description),
        p.adjust    = suppressWarnings(as.numeric(p.adjust)),
        Count       = suppressWarnings(as.integer(Count))
      ) %>%
      dplyr::filter(is.finite(p.adjust), !is.na(Count), nzchar(Description)) %>%
      dplyr::mutate(
        term_short = stringr::str_wrap(Description, width = 35),
        tooltip = paste0(
          "<b>", htmltools::htmlEscape(Description), "</b>",
          "<br><b>Ontology:</b> ", Ontology,
          "<br><b>Count:</b> ", Count,
          "<br><b>GeneRatio:</b> ", GeneRatio,
          "<br><b>BgRatio:</b> ", BgRatio,
          "<br><b>p:</b> ", sprintf("%.3g", pvalue),
          "<br><b>FDR:</b> ", sprintf("%.3g", p.adjust)
        )
      )
    
    ont_order <- c("BP","CC","MF")
    df <- df %>%
      dplyr::filter(Ontology %in% ont_order) %>%
      dplyr::group_by(Ontology) %>%
      dplyr::arrange(p.adjust, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup()
    
    validate(need(nrow(df) > 0, "GO slim: no terms to plot."))
    
    df <- df %>%
      dplyr::arrange(Ontology, Count, p.adjust) %>%
      dplyr::mutate(term_short = factor(term_short, levels = unique(term_short)))
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = term_short,
        y = Count,
        fill = Ontology,
        text = tooltip
      )
    ) +
      ggplot2::geom_col(color = "black", linewidth = 0.25) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("GO slim enrichment by term counts — ", attr(go_selected_df(), "label") %||% ""),
        x = NULL,
        y = "Enriched genes per GO slim term"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      ) +
      ggplot2::scale_fill_manual(values = c(BP = "darkgreen", CC = "orange", MF = "darkblue"))
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "goslim_enrichment_counts",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  # ----------------------------
  # GO terms table
  # ----------------------------
  
  # Robust split for GO terms
  split_go_terms <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x %in% c("", ".", "NA", "NaN", "nan")] <- NA_character_
    if (all(is.na(x))) return(character(0))
    
    out <- unlist(strsplit(x, "[;|,]+", perl = TRUE), use.names = FALSE)
    out <- trimws(out)
    out <- out[nzchar(out)]
    unique(out)
  }
  
  collapse_go_terms_html <- function(x, max_items = 80) {
    x <- unique(na.omit(trimws(as.character(x))))
    x <- x[nzchar(x)]
    if (!length(x)) return("")
    if (length(x) > max_items) {
      x <- c(x[seq_len(max_items)], paste0("… (+", length(x) - max_items, " more)"))
    }
    paste(x, collapse = "<br/>")
  }
  
  make_amigo_gene_url <- function(g) {
    g <- trimws(as.character(g))
    if (!nzchar(g)) return(NA_character_)
    paste0(
      "https://amigo.geneontology.org/amigo/medial_search?q=",
      utils::URLencode(g, reserved = TRUE),
      "&searchtype=all"
    )
  }
  
  linkify_amigo_gene <- function(g) {
    g <- as.character(g)
    u <- make_amigo_gene_url(g)
    if (is.na(u) || !nzchar(u)) return(htmltools::htmlEscape(g))
    sprintf("<a href='%s' target='_blank'>%s</a>", u, htmltools::htmlEscape(g))
  }
  
  go_terms_table_df <- reactive({
    df_use <- go_selected_df()
    validate(need(is.data.frame(df_use) && nrow(df_use) > 0, "No data for the current GO selection."))
    
    need_cols <- c("chr","genename","cluster_id",
                   "GO_biological_process","GO_cellular_component","GO_molecular_function")
    validate(need(all(need_cols %in% names(df_use)),
                  "Missing GO columns in final output (GO_biological_process/GO_cellular_component/GO_molecular_function)."))
    
    out <- df_use %>%
      dplyr::transmute(
        chr        = paste0("chr", chr_label_plink(suppressWarnings(as.integer(norm_chr_generic(.data[["chr"]]))))),
        cluster_id = as.character(.data[["cluster_id"]]),
        genename   = as.character(.data[["genename"]]),
        GO_biological_process = .data[["GO_biological_process"]],
        GO_cellular_component = .data[["GO_cellular_component"]],
        GO_molecular_function = .data[["GO_molecular_function"]]
      ) %>%
      dplyr::filter(!is.na(genename), nzchar(genename)) %>%
      dplyr::group_by(chr, cluster_id, genename) %>%
      dplyr::summarise(
        GO_biological_process = collapse_go_terms_html(unlist(lapply(GO_biological_process, split_go_terms))),
        GO_cellular_component = collapse_go_terms_html(unlist(lapply(GO_cellular_component, split_go_terms))),
        GO_molecular_function = collapse_go_terms_html(unlist(lapply(GO_molecular_function, split_go_terms))),
        .groups = "drop"
      ) %>%
      dplyr::mutate(gene = vapply(genename, linkify_amigo_gene, character(1))) %>%
      dplyr::arrange(chr, cluster_id, genename) %>%
      dplyr::select(chr, cluster_id, gene,
                    GO_biological_process, GO_cellular_component, GO_molecular_function)
    
    out
  })
  
  output$go_terms_table <- DT::renderDT({
    out <- go_terms_table_df()
    validate(need(nrow(out) > 0, "No GO terms to show for the current selection."))
    
    DT::datatable(
      out,
      rownames = FALSE,
      escape   = FALSE,  # required for links and <br/>
      width    = "100%",
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX    = TRUE,
        autoWidth  = TRUE,
        columnDefs = list(
          list(targets = c(3,4,5), width = "380px")
        )
      )
    ) %>%
      DT::formatStyle(
        columns = c("GO_biological_process","GO_cellular_component","GO_molecular_function"),
        `white-space` = "normal"
      )
  }, server = FALSE)
  
  # ==========================================================
  # GoSlim enrichment (clusterProfiler::enricher) — NonSyn
  #   Creates:
  #     goslim_enrich_raw(), goslim_enrich_tbl(), output$goslim_table
  # ==========================================================
  
  # Helper: SYMBOL -> ENTREZ mapping as data.frame (safe)
  sym2ent_df_safe <- function(symbols) {
    symbols <- unique(na.omit(trimws(as.character(symbols))))
    symbols <- symbols[nzchar(symbols)]
    if (!length(symbols)) return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
    
    # Prefer existing mapping
    if (exists("gene_map_all", mode = "function")) {
      m <- tryCatch(gene_map_all(), error = function(e) NULL)
      if (is.data.frame(m) && all(c("SYMBOL", "ENTREZID") %in% names(m))) {
        out <- m[m$SYMBOL %in% symbols, c("SYMBOL", "ENTREZID"), drop = FALSE]
        out <- out[!is.na(out$ENTREZID) & nzchar(out$ENTREZID), ]
        out$SYMBOL   <- as.character(out$SYMBOL)
        out$ENTREZID <- as.character(out$ENTREZID)
        return(unique(out))
      }
    }
    
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
    mm <- tryCatch({
      suppressMessages(clusterProfiler::bitr(
        symbols,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Hs.eg.db
      ))
    }, error = function(e) NULL)
    
    if (!is.data.frame(mm) || !nrow(mm)) return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
    mm <- mm[, c("SYMBOL", "ENTREZID"), drop = FALSE]
    mm <- mm[!is.na(mm$ENTREZID) & nzchar(mm$ENTREZID), ]
    unique(mm)
  }
  
  # Helper: just ENTREZ vector
  entrez_from_symbols_safe <- function(symbols) {
    mm <- sym2ent_df_safe(symbols)
    ids <- unique(as.character(mm$ENTREZID))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    ids
  }
  
  
  # ==========================================================
  # GoSlim TERM2GENE / TERM2NAME / universe from OrgDb + GO.db
  #   ontology = BP / CC / MF
  #   Uses:
  #     - org.Hs.eg.db for ENTREZ <-> GO
  #     - GO.db for ancestors
  #     - GO_SLIM_GENERIC for slim terms
  # ==========================================================
  goslim_build_orgdb_sets_for_onto <- function(ontology = c("BP", "CC", "MF")) {
    ontology <- match.arg(ontology)
    
    validate(need(requireNamespace("AnnotationDbi", quietly = TRUE), "AnnotationDbi package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    validate(need(requireNamespace("GO.db", quietly = TRUE), "GO.db package is required."))
    validate(need(exists("GO_SLIM_GENERIC"), "GO_SLIM_GENERIC not loaded."))
    
    # --- GO slim IDs for this ontology ---
    slim_df <- GO_SLIM_GENERIC %>%
      dplyr::filter(ONTOLOGY == ontology) %>%
      dplyr::distinct(GOID, slim_term)
    
    slim_ids <- unique(as.character(slim_df$GOID))
    slim_ids <- slim_ids[!is.na(slim_ids) & nzchar(slim_ids)]
    
    validate(need(length(slim_ids) > 0, paste0("No GoSlim IDs loaded for ", ontology)))
    
    # --- Get ENTREZ -> GO mapping from OrgDb ---
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
    
    # --- Ancestor DB for this ontology ---
    anc_obj <- switch(
      ontology,
      BP = GO.db::GOBPANCESTOR,
      CC = GO.db::GOCCANCESTOR,
      MF = GO.db::GOMFANCESTOR
    )
    
    # --- GOID -> slim GOIDs via ancestors ---
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
    
    # --- TERM2GENE: slim_goid -> ENTREZID ---
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
    
    # --- TERM2NAME ---
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
    
    list(
      term2gene = term2gene,
      term2name = term2name,
      universe  = universe
    )
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
  
  # Trigger for GoSlim (assumed to exist)
  # goslim_trigger <- reactiveVal(0L)
  
  goslim_enrich_raw <- eventReactive(goslim_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    validate(need(requireNamespace("GO.db", quietly = TRUE), "GO.db package is required."))
    
    withProgress(message = "Running GO slim enrichment…", value = 0, {
      
      ontos <- input$go_ontos %||% c("BP", "CC", "MF")
      ontos <- intersect(c("BP", "CC", "MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one GO ontology (BP/CC/MF)."))
      
      df_sel <- go_selected_df()
      validate(need(is.data.frame(df_sel) && nrow(df_sel) > 0, "No data for current selection."))
      validate(need("genename" %in% names(df_sel), "go_selected_df() must include 'genename'."))
      
      symbols_sel <- unique(na.omit(trimws(as.character(df_sel$genename))))
      symbols_sel <- symbols_sel[nzchar(symbols_sel)]
      validate(need(length(symbols_sel) > 0, "No genes found for the current selection."))
      incProgress(0.15, detail = paste("Symbols:", length(symbols_sel)))
      
      gene_selected <- entrez_from_symbols_safe(symbols_sel)
      validate(need(length(gene_selected) > 0, "No ENTREZ IDs for selected genes (mapping empty)."))
      incProgress(0.10, detail = paste("ENTREZ:", length(gene_selected)))
      
      minGS_in <- input$enrich_min_gs %||% 10
      maxGS    <- input$enrich_max_gs %||% 500
      
      gs_all <- goslim_orgdb_sets()
      incProgress(0.10, detail = "Reference ready")
      
      res <- purrr::map_dfr(ontos, function(ont) {
        incProgress(0.50 / length(ontos), detail = paste("Ontology:", ont))
        
        gs <- gs_all[[ont]]
        
        universe <- unique(gs$universe)
        gene_use <- intersect(unique(gene_selected), universe)
        
        message(
          "[GoSlim][", ont, "] ",
          "selected_genes=", length(gene_selected),
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
      
      incProgress(0.05, detail = "Finalizing")
      
      if (!is.data.frame(res) || !nrow(res)) return(tibble::tibble())
      res
    })
  }, ignoreInit = TRUE)
  
  
  
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
    
    if (!nrow(df2)) return(tibble::tibble())
    df2
  })
  
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
  
  ############################################################################
  ###############################################################################
  # ============================================================
  # NonSyn · TERM enrichment (NEW) ✅ using universe_nonsyn_cache/*.rds
  #   Background: precomputed normalized tables (term_id,label,freq,src, ...)
  #   Foreground: dbnsfp_norm_df() rows (global or per-cluster)
  #   Unit: rows (NO unique); supports multi-term fields
  #   Plot: COUNTS in foreground (y=a)
  #
  # PATCH INCLUDED:
  # - Orphanet parser fixed (NO gsub(function(...)) bug; safe split + label->ORPHA id mapping)
  # - MIM parsing + fg_label_map support for missing universe labels
  # - Table shapes: uniprot / entrez / mim / interpro
  # - Plot avoids duplicated factor levels
  # ============================================================
  ###############################################################################
  
  # --------------------------------------------
  # Helpers (general)
  # --------------------------------------------
  `%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b
  
  .norm_key <- function(x) {
    x <- enc2utf8(as.character(x))
    x <- trimws(tolower(x))
    x <- gsub("[\u00a0]", " ", x) # NBSP
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)
    x
  }
  
  build_bg_termid_map <- function(bg_df) {
    stopifnot(is.data.frame(bg_df), nrow(bg_df) > 0, "term_id" %in% names(bg_df))
    
    key_cols <- intersect(
      c("term_id","label","label_in","term","symbol","name",
        "protein_name","gene_primary","accession","id","query_id"),
      names(bg_df)
    )
    
    mp <- new.env(parent = emptyenv())
    
    for (cc in key_cols) {
      v <- bg_df[[cc]]
      ok <- !is.na(v) & nzchar(trimws(as.character(v)))
      if (!any(ok)) next
      keys <- .norm_key(v[ok])
      tids <- as.character(bg_df$term_id[ok])
      
      for (i in seq_along(keys)) {
        k <- keys[i]
        if (!exists(k, envir = mp, inherits = FALSE)) mp[[k]] <- tids[i]
      }
    }
    
    mp
  }
  
  pick_col <- function(df, candidates) {
    nm <- names(df)
    hit <- candidates[candidates %in% nm]
    if (length(hit)) hit[1] else NULL
  }
  
  short_label <- function(x, max_chars = 45, wrap_width = 16) {
    x <- as.character(x)
    x <- stringr::str_squish(x)
    x <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
    stringr::str_wrap(x, width = wrap_width)
  }
  
  # --------------------------------------------
  # Split multi-term fields (IMPORTANT: diseases NO comma split)
  # --------------------------------------------
  split_terms_vec <- function(v, field = NULL) {
    v <- as.character(v)
    v <- v[!is.na(v)]
    v <- trimws(v)
    v <- v[nzchar(v)]
    if (!length(v)) return(character(0))
    
    fld <- tolower(as.character(field %||% ""))
    
    pat <- if (fld %in% c("mim_disease", "orphanet_disorder", "clinvar_trait")) {
      "[;|&]"   # <- NO comma
    } else {
      "[;|,&]"
    }
    
    toks <- unlist(strsplit(v, pat, perl = TRUE), use.names = FALSE)
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    if (!length(toks)) return(character(0))
    
    low <- tolower(toks)
    low <- gsub("[_\\-]+", " ", low)
    low <- gsub("\\s+", " ", low)
    low <- trimws(low)
    
    drop <- low %in% c("na","n/a","nan","null","none",".","-","_","") |
      grepl("^(not specified|no specified|not provided|unspecified|unknown|no information|not reported|not available)$", low) |
      grepl("^[[:punct:][:space:]]+$", toks)
    
    toks <- toks[!drop]
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    toks
  }
  
  # --------------------------------------------
  # MIM helpers
  # --------------------------------------------
  clean_mim_label <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- trimws(x)
    # remove leading "[MIM:123456]" or "MIM:123456"
    x <- gsub("^\\s*\\[?MIM:\\d+\\]?\\s*", "", x, perl = TRUE)
    # remove trailing "[dominant?]" etc (may be multiple; do twice)
    x <- gsub("\\s*\\[[^\\]]*\\]\\s*$", "", x, perl = TRUE)
    x <- gsub("\\s*\\[[^\\]]*\\]\\s*$", "", x, perl = TRUE)
    x <- gsub("\\s+", " ", x, perl = TRUE)
    trimws(x)
  }
  
  extract_mim_id_label_map <- function(raw_vec) {
    raw_vec <- as.character(raw_vec)
    raw_vec <- raw_vec[!is.na(raw_vec) & nzchar(trimws(raw_vec))]
    if (!length(raw_vec)) return(stats::setNames(character(0), character(0)))
    
    parts <- unlist(strsplit(raw_vec, ";", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    
    ids <- sub(".*\\[\\s*(MIM:\\d+)\\s*\\].*", "\\1", parts, perl = TRUE)
    ok  <- grepl("^MIM:\\d+$", ids)
    
    lab <- parts
    lab <- sub(".*\\[\\s*MIM:\\d+\\s*\\]\\s*", "", lab, perl = TRUE)
    lab <- sub("\\s*\\[[^\\]]*\\]\\s*$", "", lab, perl = TRUE)
    lab <- trimws(lab)
    
    ids <- ids[ok]
    lab <- lab[ok]
    
    out <- tapply(lab, ids, function(z) {
      z <- z[!is.na(z) & nzchar(z)]
      if (length(z)) z[1] else ""
    })
    
    out <- as.character(out)
    stats::setNames(out, names(out))
  }
  
  # --------------------------------------------
  # Orphanet helpers (FIXED)
  # --------------------------------------------
  
  # Replace semicolons inside parentheses iteratively (NO gsub(function(...))!)
  protect_semicolons_in_parens <- function(s) {
    s <- as.character(s)
    if (!length(s)) return(s)
    
    for (i in seq_along(s)) {
      x <- s[i]
      if (is.na(x) || !nzchar(x)) next
      while (grepl("\\([^)]*;[^)]*\\)", x, perl = TRUE)) {
        x <- sub("(\\([^)]*);([^)]*\\))", "\\1\u241E\\2", x, perl = TRUE)
      }
      s[i] <- x
    }
    s
  }
  
  split_orphanet_safe <- function(raw_vec) {
    raw <- as.character(raw_vec)
    raw <- raw[!is.na(raw)]
    raw <- trimws(raw)
    raw <- raw[nzchar(raw)]
    if (!length(raw)) return(character(0))
    
    raw2 <- protect_semicolons_in_parens(raw)
    parts <- unlist(strsplit(raw2, ";", fixed = TRUE), use.names = FALSE)
    parts <- gsub("\u241E", ";", parts, fixed = TRUE)
    
    parts <- trimws(parts)
    parts <- gsub("\\s+", " ", parts)
    parts <- parts[nzchar(parts)]
    
    parts <- sub("^Orphanet:\\s*", "", parts, ignore.case = TRUE)
    parts
  }
  
  norm_orpha_key <- function(z) {
    z <- enc2utf8(as.character(z))
    z <- trimws(z)
    z <- gsub("[\u00a0]", " ", z) # NBSP
    z <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", z)  # unicode dashes -> "-"
    z <- gsub("\\s+", " ", z)
    tolower(z)
  }
  
  map_orphanet_labels_to_ids <- function(labels, bg_df) {
    labels <- as.character(labels)
    labels[is.na(labels)] <- ""
    labels <- trimws(labels)
    
    if (!length(labels)) return(character(0))
    if (!is.data.frame(bg_df) || !nrow(bg_df)) return(rep(NA_character_, length(labels)))
    
    lab_col <- if ("label_in" %in% names(bg_df)) "label_in" else if ("label" %in% names(bg_df)) "label" else NULL
    if (is.null(lab_col) || !("orpha_id" %in% names(bg_df))) return(rep(NA_character_, length(labels)))
    
    kk <- norm_orpha_key(bg_df[[lab_col]])
    vv <- as.character(bg_df$orpha_id)
    
    ok <- !is.na(kk) & nzchar(kk) & !is.na(vv) & nzchar(vv)
    if (!any(ok)) return(rep(NA_character_, length(labels)))
    
    lut <- stats::setNames(paste0("ORPHA:", vv[ok]), kk[ok])
    
    lk <- norm_orpha_key(labels)
    mapped <- unname(lut[lk])
    mapped
  }
  
  # --------------------------------------------
  # Foreground terms -> term_id (patched)
  # --------------------------------------------
  fg_terms_to_term_id <- function(field, vec, bg_df = NULL) {
    fld <- tolower(as.character(field))
    fld_low <- tolower(as.character(field))
    
    # --- Use split_terms_vec() if you have it; otherwise fallback ---
    if (exists("split_terms_vec", mode = "function")) {
      toks <- split_terms_vec(vec, field = field)
    } else {
      v <- as.character(vec)
      v <- v[!is.na(v)]
      v <- trimws(v)
      v <- v[nzchar(v)]
      pat <- if (fld %in% c("mim_disease", "orphanet_disorder", "clinvar_trait")) "[;|&]" else "[;|,&]"
      toks <- unlist(strsplit(v, pat, perl = TRUE), use.names = FALSE)
      toks <- trimws(toks)
      toks <- toks[nzchar(toks)]
    }
    
    if (!length(toks)) return(character(0))
    
    # -----------------------------
    # Orphanet disorder
    # -----------------------------
    if (fld_low == "orphanet_disorder") {
      
      raw_chr <- as.character(vec)
      raw_chr <- raw_chr[!is.na(raw_chr)]
      raw_chr <- trimws(raw_chr)
      raw_chr <- raw_chr[nzchar(raw_chr)]
      if (!length(raw_chr)) return(character(0))
      
      # split ONLY by ';' safely (keeps per-row multi terms)
      parts <- split_orphanet_safe(raw_chr)
      parts <- parts[!is.na(parts)]
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      if (!length(parts)) return(character(0))
      
      # A) direct ORPHA:\d+ occurrences (KEEP DUPLICATES)
      ids_vec <- unlist(regmatches(parts, gregexpr("ORPHA:\\d+", parts, perl = TRUE)), use.names = FALSE)
      ids_vec <- trimws(ids_vec)
      ids_vec <- ids_vec[nzchar(ids_vec)]
      
      # B) labels (no ORPHA id inside) -> map to ORPHA:<id> (KEEP DUPLICATES)
      labels_only <- parts[!grepl("ORPHA:\\d+", parts, perl = TRUE)]
      labels_only <- trimws(labels_only)
      labels_only <- labels_only[nzchar(labels_only)]
      
      mapped_vec <- character(0)
      if (length(labels_only)) {
        mapped_vec <- map_orphanet_labels_to_ids(labels_only, bg_df)
        mapped_vec <- as.character(mapped_vec)
        mapped_vec <- mapped_vec[!is.na(mapped_vec) & nzchar(mapped_vec)]
      }
      
      out <- c(ids_vec, mapped_vec)
      out <- out[nzchar(out)]
      return(out)  # IMPORTANT: no unique()
    }
    
    # -----------------------------
    # MIM disease
    # -----------------------------
    if (fld_low == "mim_disease") {
      x <- as.character(vec)
      x <- x[!is.na(x)]
      x <- trimws(x)
      x <- x[nzchar(x)]
      if (!length(x)) return(character(0))
      
      # IMPORTANT: keep duplicates (no unique)
      mim <- unlist(regmatches(x, gregexpr("MIM:\\d+", x, perl = TRUE)), use.names = FALSE)
      mim <- trimws(mim)
      mim <- mim[nzchar(mim)]
      return(mim)
    }
    
    # -----------------------------
    # ClinVar trait (universe uses TERM:<normalized>)
    # foreground uses '|' between traits and commas inside trait names
    # -----------------------------
    if (fld_low == "clinvar_trait") {
      
      toks <- split_terms_vec(vec, field = "clinvar_trait")  # splits by ; | &
      toks <- as.character(toks)
      toks <- toks[!is.na(toks) & nzchar(trimws(toks))]
      if (!length(toks)) return(character(0))
      
      norm <- enc2utf8(toks)
      norm <- trimws(norm)
      
      norm <- gsub("[\u00a0]", " ", norm)  # NBSP
      norm <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", norm) # unicode dashes -> '-'
      norm <- gsub(",", "_", norm, fixed = TRUE)   # <<< clau: comes -> underscore
      norm <- gsub("\\s+", "_", norm)              # espais -> underscore
      norm <- gsub("[^0-9A-Za-z_\\-]+", "", norm)  # conserva _ i -
      norm <- gsub("_+", "_", norm)
      
      norm <- norm[nzchar(norm)]
      out  <- paste0("TERM:", norm)
      out  <- out[nzchar(out) & out != "TERM:"]
      
      # opcional: filtra a universe
      if (is.data.frame(bg_df) && "term_id" %in% names(bg_df)) {
        u <- unique(as.character(bg_df$term_id))
        out <- out[out %in% u]
      }
      
      # IMPORTANT: do NOT unique() -> we need true counts for 'a'
      return(out)
    }
    
    # -----------------------------
    # InterPro domain (your universe uses IPR_LABEL:<label>)
    # -----------------------------
    if (fld == "interpro_domain") {
      toks <- trimws(as.character(toks))
      toks <- toks[nzchar(toks)]
      return(paste0("IPR_LABEL:", toks))
    }
    
    # default
    toks
  }
  
  # --------------------------------------------
  # Universe cache loader
  # --------------------------------------------
  load_universe_nonsyn_cache <- function(dir_cache) {
    stopifnot(dir.exists(dir_cache))
    files <- list.files(dir_cache, pattern = "_norm\\.rds$", full.names = TRUE)
    validate(need(length(files) > 0, paste0("No *_norm.rds found in: ", dir_cache)))
    
    lst <- setNames(lapply(files, readRDS), basename(files))
    lst <- lapply(lst, function(df) {
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      validate(need(all(c("term_id","label","freq") %in% names(df)),
                    "Universe table missing required cols: term_id,label,freq"))
      df$term_id <- as.character(df$term_id)
      df$label   <- as.character(df$label)
      df$freq    <- suppressWarnings(as.integer(df$freq))
      df
    })
    lst
  }
  
  universe_key_for_field <- function(fld) {
    fld_low <- tolower(as.character(fld))
    if (fld_low == "interpro_domain")    return("tabU_interpro_domain_norm.rds")
    if (fld_low == "mim_disease")        return("tabU_mim_disease_norm.rds")
    if (fld_low == "orphanet_disorder")  return("tabU_orphanet_disorder_norm.rds")
    if (fld_low == "clinvar_trait")      return("tabU_clinvar_trait_norm.rds")
    NA_character_
  }
  
  # --------------------------------------------
  # ORA from counts (term_id + freq)
  # --------------------------------------------
  ora_from_counts_norm <- function(tabS, tabU, min_a = 1, min_Tt = 1, max_Tt = Inf, max_terms = 8000) {
    if (is.null(tabS) || length(tabS) == 0) return(data.table::data.table())
    if (!is.data.frame(tabU) || !nrow(tabU)) return(data.table::data.table())
    
    tabU <- as.data.frame(tabU, stringsAsFactors = FALSE)
    validate(need(all(c("term_id","freq") %in% names(tabU)), "Universe table must have term_id,freq"))
    
    tabU$term_id <- as.character(tabU$term_id)
    tabU$freq    <- suppressWarnings(as.integer(tabU$freq))
    
    terms <- intersect(names(tabS), tabU$term_id)
    if (!length(terms)) return(data.table::data.table())
    
    a  <- as.integer(tabS[terms])
    Tt <- as.integer(tabU$freq[match(terms, tabU$term_id)])
    
    keep <- which(a >= min_a & Tt >= min_Tt & Tt <= max_Tt)
    if (!length(keep)) return(data.table::data.table())
    
    terms <- terms[keep]; a <- a[keep]; Tt <- Tt[keep]
    
    if (length(terms) > max_terms) {
      ord <- order(Tt, decreasing = TRUE)[1:max_terms]
      terms <- terms[ord]; a <- a[ord]; Tt <- Tt[ord]
    }
    
    n <- sum(as.integer(tabS))
    N <- sum(as.integer(tabU$freq))
    
    p <- stats::phyper(a - 1L, Tt, N - Tt, n, lower.tail = FALSE)
    
    b <- n - a
    c <- Tt - a
    d <- (N - Tt) - b
    or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    
    out <- data.table::data.table(
      term_id = terms,
      a       = a,
      Tt      = Tt,
      or      = as.numeric(or),
      p       = as.numeric(p)
    )
    out[, p_adj := p.adjust(p, method = "BH")]
    out[, `:=`(n_total = n, N_total = N)]
    
    meta_cols <- setdiff(names(tabU), "freq")
    meta <- tabU[, meta_cols, drop = FALSE]
    meta <- meta[!duplicated(meta$term_id), , drop = FALSE]
    out <- merge(out, meta, by = "term_id", all.x = TRUE, sort = FALSE)
    
    data.table::setorder(out, p_adj, p)
    out
  }
  
  # --------------------------------------------
  # UI: selectable term fields
  # --------------------------------------------
  output$nonsyn_term_field_ui <- renderUI({
    df <- tryCatch(dbnsfp_norm_df(), error = function(e) NULL)
    if (!is.data.frame(df) || nrow(df) == 0) return(NULL)
    
    cn <- names(df)
    cn_low <- tolower(cn)
    
    pick_real <- function(want_name) {
      i <- match(tolower(want_name), cn_low)
      if (is.na(i)) return(NULL)
      cn[i]
    }
    
    want <- c(
      "Orphanet disorder" = "orphanet_disorder",
      "MIM disease"       = "mim_disease",
      "ClinVar trait"     = "clinvar_trait",
      "InterPro domain"   = "interpro_domain"
    )
    
    keep <- list()
    for (lab in names(want)) {
      real <- pick_real(unname(want[[lab]]))
      if (!is.null(real)) keep[[lab]] <- real
    }
    
    if (!length(keep)) return(helpText("No enrichable term columns found in dbNSFP table."))
    
    selectInput("nonsyn_term_field", "Term field", choices = keep, selected = unname(keep[[1]]))
  })
  
  # --------------------------------------------
  # Foreground bundle
  # --------------------------------------------
  dt_enrich_res_nonsyn <- reactive({
    dt_all <- dbnsfp_norm_df()
    validate(need(is.data.frame(dt_all) && nrow(dt_all) > 0,
                  "No NonSyn dbNSFP rows for current analysis yet (run extraction first)."))
    dt_all <- as.data.frame(dt_all)
    
    # -----------------------------
    # Detect cluster column
    # -----------------------------
    cc <- NULL
    for (cand in c("cluster_id","cluster_chr_n","cluster","cluster_chr")) {
      if (cand %in% names(dt_all)) { cc <- cand; break }
    }
    validate(need(!is.null(cc), "Foreground dbNSFP table has no cluster column (cluster_id/cluster_chr_n/cluster)."))
    
    dt_all[[cc]] <- sub("^cluster_", "", as.character(dt_all[[cc]]))
    
    dt2 <- dt_all
    
    # -----------------------------
    # Scope: cluster
    # -----------------------------
    if (identical(input$func_scope %||% "global", "cluster")) {
      req(input$func_cluster_id)
      cid <- sub("^cluster_", "", as.character(input$func_cluster_id))
      dt2 <- dt2[as.character(dt2[[cc]]) == cid, , drop = FALSE]
    }
    
    # -----------------------------
    # Scope: gene  (NEW)
    # -----------------------------
    if (identical(input$func_scope %||% "global", "gene")) {
      req(input$func_gene)
      
      gene_col <- detect_gene_col(dt2)
      validate(need(!is.null(gene_col),
                    paste0("Foreground dbNSFP table has no gene column. ",
                           "Expected e.g. genename/gene/symbol/Gene_name.")))
      
      g <- as.character(input$func_gene)
      dt2 <- dt2[as.character(dt2[[gene_col]]) == g, , drop = FALSE]
    }
    
    validate(need(nrow(dt2) > 0, "No NonSyn hits in this scope (empty after filtering)."))
    
    list(dt_all = dt_all, dt2 = dt2, cluster_col = cc)
  })
  
  # --------------------------------------------
  # Trigger: enrichment
  # --------------------------------------------
  nonsyn_terms_res <- eventReactive(input$run_enrich, {
    
    x <- dt_enrich_res_nonsyn()
    dt2 <- x$dt2
    
    fld <- input$nonsyn_term_field %||% ""
    validate(need(nzchar(fld), "Select a term field."))
    validate(need(fld %in% names(dt2), paste0("Selected field not present in foreground: ", fld)))
    
    dir_bg <- "/Volumes/DISK1TB/Inspector_app_slaves_ngroc/GItools/_shared/universe_nonsyn_cache"
    bg_all <- load_universe_nonsyn_cache(dir_bg)
    
    bg_key <- universe_key_for_field(fld)
    validate(need(!is.na(bg_key) && bg_key %in% names(bg_all),
                  paste0("No universe background available for field: ", fld)))
    
    tabU <- bg_all[[bg_key]]
    
    # Debug input examples
    rawv <- dt2[[fld]]
    rawv <- as.character(rawv)
    rawv <- rawv[!is.na(rawv) & nzchar(trimws(rawv)) & trimws(rawv) != "."]
    message("[ENRICH DBG] field=", fld,
            " | raw_nonempty_n=", length(rawv),
            " | raw_head=", paste(head(rawv, 3), collapse = " | "))
    
    raw_chr <- as.character(dt2[[fld]])
    raw_chr <- raw_chr[!is.na(raw_chr) & nzchar(trimws(raw_chr))]
    message("[ENRICH CHECK] field=", fld, " | raw_nonempty_n=", length(raw_chr))
    message("[ENRICH CHECK] raw examples:\n- ", paste(head(unique(raw_chr), 12), collapse = "\n- "))
    
    if (tolower(fld) == "mim_disease") {
      ids <- unique(unlist(regmatches(raw_chr, gregexpr("MIM:\\d+", raw_chr, perl = TRUE)), use.names = FALSE))
      message("[ENRICH CHECK] embedded MIM ids found n=", length(ids),
              " head=", paste(head(ids, 12), collapse = ", "))
    }
    
    if (tolower(fld) == "orphanet_disorder") {
      ids <- unique(unlist(regmatches(raw_chr, gregexpr("ORPHA:\\d+", raw_chr, perl = TRUE)), use.names = FALSE))
      message("[ENRICH CHECK] embedded ORPHA ids found n=", length(ids),
              " head=", paste(head(ids, 12), collapse = ", "))
    }
    
    fg_term_ids <- fg_terms_to_term_id(fld, dt2[[fld]], bg_df = tabU)
    
    # DEBUG: verify duplicates in foreground
    tabS_all <- table(fg_term_ids)
    message("[ENRICH DBG] field=", fld,
            " | fg_tokens=", length(fg_term_ids),
            " | fg_unique=", length(unique(fg_term_ids)),
            " | max_a=", if (length(tabS_all)) max(as.integer(tabS_all)) else NA_integer_,
            " | a_ge2=", sum(as.integer(tabS_all) >= 2),
            " | a_ge3=", sum(as.integer(tabS_all) >= 3))
    
    print(sort(tabS_all, decreasing = TRUE)[1:min(15, length(tabS_all))])
    
    if (tolower(fld) == "mim_disease") {
      tabS_all <- table(fg_term_ids)
      message("[MIM DBG] fg_tokens=", length(fg_term_ids),
              " | fg_unique=", length(unique(fg_term_ids)),
              " | max_a=", if (length(tabS_all)) max(as.integer(tabS_all)) else NA_integer_,
              " | a_ge2=", sum(as.integer(tabS_all) >= 2),
              " | a_ge3=", sum(as.integer(tabS_all) >= 3))
      print(sort(tabS_all, decreasing = TRUE)[1:min(15, length(tabS_all))])
    }
    
    # Map for MIM labels from foreground (to fill empty universe label)
    fg_label_map <- NULL
    if (tolower(fld) == "mim_disease") {
      fg_label_map <- extract_mim_id_label_map(dt2[[fld]])
    }
    
    # Debug intersect
    u <- unique(as.character(tabU$term_id))
    f <- unique(as.character(fg_term_ids))
    message("[ENRICH CHECK] fg unique=", length(f),
            " | universe unique=", length(u),
            " | intersect=", sum(f %in% u))
    if (sum(f %in% u) == 0) {
      message("[ENRICH CHECK] fg head:\n- ", paste(head(f, 20), collapse = "\n- "))
      message("[ENRICH CHECK] tabU head:\n- ", paste(head(u, 20), collapse = "\n- "))
    }
    
    message("[ENRICH DBG] fg_term_ids_n=", length(fg_term_ids),
            " | fg_head=", paste(head(unique(fg_term_ids), 6), collapse = " | "))
    message("[ENRICH DBG] tabU_term_id_head=",
            paste(head(unique(as.character(tabU$term_id)), 6), collapse = " | "))
    
    topk  <- input$nonsyn_terms_topk %||% 200
    min_a <- input$nonsyn_terms_min_a %||% 2
    topk  <- suppressWarnings(as.integer(topk)); if (!is.finite(topk) || topk < 20L) topk <- 200L
    min_a <- suppressWarnings(as.integer(min_a)); if (!is.finite(min_a) || min_a < 1L) min_a <- 1L
    
    minGS <- input$enrich_min_gs %||% 1
    maxGS <- input$enrich_max_gs %||% 500
    minGS <- suppressWarnings(as.integer(minGS)); if (!is.finite(minGS) || minGS < 1L) minGS <- 1L
    maxGS <- suppressWarnings(as.numeric(maxGS)); if (!is.finite(maxGS) || maxGS < minGS) maxGS <- 500000
    
    fld_low <- tolower(fld)
    trait_like <- fld_low %in% c("orphanet_disorder","mim_disease") || grepl("trait", fld_low) || grepl("disease", fld_low)
    max_Tt_use <- if (isTRUE(trait_like)) Inf else max(maxGS, 500000)
    
    withProgress(message = "Running NonSyn term enrichment…", value = 0, {
      
      incProgress(0.20, detail = "Foreground parsing…")
      
      fg_term_ids <- fg_terms_to_term_id(fld, dt2[[fld]], bg_df = tabU)
      
      if (!length(fg_term_ids)) {
        ex <- unique(na.omit(as.character(dt2[[fld]])))
        ex <- ex[nzchar(ex)]
        ex <- head(ex, 5)
        validate(need(FALSE, paste0(
          "No valid terms parsed in foreground for field: ", fld, ". ",
          "Likely all NA/empty OR parsing mismatch. ",
          "Examples: ", paste(ex, collapse = " | ")
        )))
      }
      
      tabS_all <- table(fg_term_ids)
      
      # AUTO min_a for sparse fields (diseases/traits often have a=1 only)
      if (tolower(fld) %in% c("mim_disease", "orphanet_disorder", "clinvar_trait")) {
        if (length(tabS_all) > 0 && max(as.integer(tabS_all)) < min_a) {
          message("[ENRICH] lowering min_a to 1 for sparse disease/trait terms (", fld, ")")
          min_a <- 1L
        }
      }
      
      tabS_keep <- tabS_all[as.integer(tabS_all) >= min_a]
      
      message("[ENRICH CHECK] tabS_keep n=", length(tabS_keep),
              " | min_a=", min_a,
              " | hits_in_universe=", sum(names(tabS_keep) %in% tabU$term_id))
      
      
      u_terms <- unique(as.character(tabU$term_id))
      hit_n   <- sum(names(tabS_keep) %in% u_terms)
      
      validate(need(hit_n > 0,
                    paste0("Foreground terms found for field ", fld,
                           ", but none match universe term_id. ",
                           "Tip: for Orphanet you likely need mapping name→ORPHA id via label_in.")))
      
      validate(need(length(tabS_all) > 0, "Foreground has 0 valid terms after cleaning (all NA / not specified)."))
      validate(need(length(tabS_keep) > 0, "After min_a filter, 0 terms remain. Lower min_a."))
      
      tabS_keep <- sort(tabS_keep, decreasing = TRUE)
      if (length(tabS_keep) > topk) tabS_keep <- tabS_keep[seq_len(topk)]
      
      incProgress(0.45, detail = "ORA…")
      tab <- ora_from_counts_norm(tabS_keep, tabU, min_a = min_a, min_Tt = minGS, max_Tt = max_Tt_use, max_terms = 8000)
      
      incProgress(0.35, detail = "Done.")
      list(tab = tab, field = fld, topk = topk, min_a = min_a, bg_key = bg_key, fg_label_map = fg_label_map)
    })
    
  }, ignoreInit = TRUE)
  
  #--------------------------------------------
  strip_orphanet_prefix <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- trimws(x)
    x <- sub("^Orphanet:\\s*", "", x, ignore.case = TRUE)
    trimws(x)
  }
  # --------------------------------------------
  # Outputs: table
  # --------------------------------------------
  output$enrich_nonsyn_terms_table <- DT::renderDT({
    r <- nonsyn_terms_res(); req(r)
    tab <- r$tab
    fld <- r$field
    fld_low <- tolower(as.character(fld))
    
    fg_label_map <- r$fg_label_map
    
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(DT::datatable(
        data.frame(Message = paste0("No enriched terms for field: ", fld)),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    n <- tab$n_total[1]
    N <- tab$N_total[1]
    
    tab <- tab %>%
      dplyr::mutate(
        in_subset_ratio = sprintf("%d/%d (%.2f%%)", a, n, 100 * a / max(1, n)),
        in_bg_ratio     = sprintf("%d/%d (%.2f%%)", Tt, N, 100 * Tt / max(1, N)),
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format = "e", digits = 2),
        FDR = formatC(as.numeric(p_adj), format = "e", digits = 2)
      )
    
    # dins output$enrich_nonsyn_terms_table ...
    
    if (fld_low == "mim_disease") {
      
      dis <- if ("label" %in% names(tab)) as.character(tab$label) else rep("", nrow(tab))
      dis[is.na(dis)] <- ""
      dis <- trimws(dis)
      
      if (!is.null(fg_label_map) && length(fg_label_map)) {
        fill <- unname(fg_label_map[as.character(tab$term_id)])
        fill[is.na(fill)] <- ""
        fill <- trimws(as.character(fill))
        dis <- ifelse(nzchar(dis), dis, fill)
      }
      
      # Fill missing disease names from OMIM dictionary (RDS)
      mp_omim <- omim_title_map()
      if (length(mp_omim)) {
        keys  <- mim_termid_to_num(as.character(tab$term_id))
        fill2 <- unname(mp_omim[keys])
        fill2[is.na(fill2)] <- ""
        fill2 <- trimws(fill2)
        dis <- ifelse(nzchar(dis), dis, fill2)
      }
      
      dis <- ifelse(nzchar(dis), dis, as.character(tab$term_id))
      dis <- clean_mim_label(dis)
      
      out <- tab %>%
        dplyr::transmute(
          `OMIM id` = as.character(.data$term_id),
          Disease   = dis,
          in_subset_ratio, in_bg_ratio, OR, p, FDR
        )
      
    } else if (fld_low == "orphanet_disorder") {
      
      dis0 <- dplyr::coalesce(as.character(tab$label), as.character(tab$label_in), "")
      dis0 <- strip_orphanet_prefix(dis0)
      dis0 <- ifelse(nzchar(dis0), dis0, as.character(tab$term_id))
      
      out <- tab %>%
        dplyr::mutate(Disorder = dis0) %>%
        dplyr::transmute(
          `Orphanet id` = as.character(.data$term_id),
          Disorder,
          in_subset_ratio, in_bg_ratio, OR, p, FDR
        )
      
    } else if (fld_low == "clinvar_trait") {
      
      # safe label: use label if present, else fallback to term_id
      trait_lab <- if ("label" %in% names(tab)) as.character(tab$label) else rep("", nrow(tab))
      trait_lab[is.na(trait_lab)] <- ""
      trait_lab <- trimws(trait_lab)
      trait_lab <- ifelse(nzchar(trait_lab), trait_lab, as.character(tab$term_id))
      
      out <- tab %>%
        dplyr::mutate(
          Trait_id = as.character(.data$term_id),
          Trait    = trait_lab
        ) %>%
        dplyr::transmute(
          Trait_id, Trait,
          in_subset_ratio, in_bg_ratio, OR, p, FDR
        )
      
    } else if (fld_low == "interpro_domain") {
      
      out <- tab %>%
        dplyr::transmute(
          `Domain name` = dplyr::coalesce(as.character(.data$label), ""),
          in_subset_ratio, in_bg_ratio, OR, p, FDR
        )
      
    } else {
      
      out <- tab %>%
        dplyr::transmute(
          term_id = as.character(.data$term_id),
          label   = dplyr::coalesce(as.character(.data$label), ""),
          in_subset_ratio, in_bg_ratio, OR, p, FDR
        )
    }
    
    DT::datatable(
      out, rownames = FALSE,
      extensions = "Buttons",
      options = list(
        extensions = "Buttons",
        dom = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 10,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  # --------------------------------------------
  # Outputs: plot (COUNTS)
  # --------------------------------------------
  
  short_label_1line <- function(x, max_chars = 55) {
    x <- as.character(x)
    x <- stringr::str_squish(x)               # treu salts/espais dobles
    x <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
    x
  }
  
  output$enrich_nonsyn_terms_plotly <- plotly::renderPlotly({
    r <- nonsyn_terms_res(); req(r)
    tab <- r$tab
    fld <- r$field
    fld_low <- tolower(as.character(fld))
    
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(plotly::plotly_empty(type = "bar") %>% plotly::layout(title = "No enriched terms"))
    }
    
    fg_label_map <- r$fg_label_map
    
    # 1) pick top rows first
    top <- tab %>%
      dplyr::filter(is.finite(p_adj), is.finite(p)) %>%
      dplyr::arrange(p_adj, p) %>%
      dplyr::slice_head(n = 15)
    
    # 2) build lab0 OUTSIDE dplyr NSE (fixes '... incorrect context')
    lab0 <- if ("label" %in% names(top)) as.character(top$label) else rep("", nrow(top))
    lab0[is.na(lab0)] <- ""
    lab0 <- ifelse(nzchar(trimws(lab0)), lab0, as.character(top$term_id))
    
    # MIM: fill missing labels then clean
    if (fld_low == "mim_disease") {
      
      # a) foreground-derived names (if any)
      if (!is.null(fg_label_map) && length(fg_label_map)) {
        fill <- unname(fg_label_map[as.character(top$term_id)])
        fill[is.na(fill)] <- ""
        fill <- trimws(as.character(fill))
        miss <- (!nzchar(trimws(lab0)) | lab0 == as.character(top$term_id)) & nzchar(fill)
        lab0[miss] <- fill[miss]
      }
      
      # b) OMIM dictionary (RDS)
      mp_omim <- omim_title_map()
      if (length(mp_omim)) {
        keys  <- mim_termid_to_num(as.character(top$term_id))
        fill2 <- unname(mp_omim[keys])
        fill2[is.na(fill2)] <- ""
        fill2 <- trimws(as.character(fill2))
        miss2 <- (!nzchar(trimws(lab0)) | lab0 == as.character(top$term_id)) & nzchar(fill2)
        lab0[miss2] <- fill2[miss2]
      }
      
      lab0 <- clean_mim_label(lab0)
    }
    
    # Orphanet: strip "Orphanet:"
    if (fld_low == "orphanet_disorder") {
      lab0 <- strip_orphanet_prefix(lab0)
    }
    
    # 3) one-line short labels for axis
    term_show_vec <- short_label_1line(lab0, max_chars = 55)
    
    # 4) attach plotting columns
    top <- top %>%
      dplyr::mutate(
        term_show = term_show_vec,
        term_key  = make.unique(as.character(term_show)),
        term_key  = factor(term_key, levels = rev(unique(term_key))),
        zebra     = factor(seq_len(dplyr::n()) %% 2),
        
        p_txt   = formatC(as.numeric(p),     format = "e", digits = 2),
        fdr_txt = formatC(as.numeric(p_adj), format = "e", digits = 2),
        or_txt  = sprintf("%.3f", as.numeric(or)),
        tt_txt  = as.integer(Tt),
        a_txt   = as.integer(a),
        
        tooltip = paste0(
          "<b>", as.character(term_show), "</b>",
          "<br><b>a</b> (fg count): ", a_txt,
          "<br><b>Tt</b> (bg count): ", tt_txt,
          "<br><b>OR</b>: ", or_txt,
          "<br><b>p</b>: ", p_txt,
          "<br><b>FDR</b>: ", fdr_txt
        )
      )
    
    top <- top %>%
      dplyr::mutate(
        neglogFDR = -log10(pmax(as.numeric(p_adj), 1e-300))
      )
    
    # 5) ggplot + ggplotly
    g <- ggplot2::ggplot(top, ggplot2::aes(x = term_key, y = a, fill = neglogFDR, text = tooltip)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_x_discrete(labels = stats::setNames(as.character(top$term_show), as.character(top$term_key))) +
      ggplot2::scale_fill_gradientn(
        colours = c("yellow", "orange", "red"),
        values  = c(0, 0.5, 1),
        name    = "−log10(FDR)"
      ) +
      ggplot2::labs(
        title = paste0("Top enriched terms (counts) — ", fld),
        x = NULL,
        y = paste0("Nonsynonym hit count in ", input$nonsyn_term_field, " term")
      ) +
      ggplot2::theme_minimal(base_size = 13)
    
    plotly::ggplotly(g, tooltip = "text", source = "terms_enrich") %>%
      plotly::layout(
        hoverlabel = list(align = "left"),
        yaxis = list(tickfont = list(size = 10)),
        margin = list(l = 240, r = 30, t = 60, b = 40)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "enrichment_nonsyn_terms",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  ###############################################################################
  # ============================================================
  # NonSyn · ENRICHMENT PREDICTORS (FULL BLOCK)
  # - Foreground: dt_enrich_res_nonsyn()$dt2 (global/cluster/gene scope)
  # - Background: predictor_universe.rds (column,item,N)
  # - ORA: Fisher exact (greater) per (predictor,item) + BH-FDR
  # - Outputs: DT + barplot (X=counts; fill=−log10(FDR) with yellow→orange→red)
  # - Adds: full predictor name + item meaning
  # ============================================================
  
  
  # -----------------------------
  # Helpers: predictor full name
  # -----------------------------
  pred_full_name <- function(column) {
    col <- as.character(column)
    
    # normalize common aliases you might have in df/universe
    if (col == "M_CAP_pred") col <- "M-CAP_pred"
    if (col == "LIST_S2_pred") col <- "LIST-S2_pred"
    if (col == "fathmm_XF_coding_pred") col <- "fathmm-XF_coding_pred"
    
    switch(col,
           "SIFT_pred"                  = "SIFT (Sorting Intolerant From Tolerant) prediction",
           "SIFT4G_pred"                = "SIFT4G prediction",
           "Polyphen2_HDIV_pred"        = "PolyPhen-2 HDIV (HumDiv) prediction",
           "Polyphen2_HVAR_pred"        = "PolyPhen-2 HVAR (HumVar) prediction",
           "MutationTaster_pred"        = "MutationTaster2 prediction",
           "MutationAssessor_pred"      = "MutationAssessor functional impact prediction",
           "PROVEAN_pred"               = "PROVEAN (Protein Variation Effect Analyzer) prediction",
           "LRT_pred"                   = "LRT (Likelihood Ratio Test) prediction",
           "MetaSVM_pred"               = "MetaSVM ensemble prediction",
           "MetaLR_pred"                = "MetaLR ensemble prediction",
           "MetaRNN_pred"               = "MetaRNN ensemble prediction",
           "M-CAP_pred"                 = "M-CAP (Mendelian Clinically Applicable Pathogenicity) prediction",
           "PrimateAI_pred"             = "PrimateAI prediction",
           "DEOGEN2_pred"               = "DEOGEN2 prediction",
           "BayesDel_addAF_pred"        = "BayesDel (with MaxAF) prediction",
           "BayesDel_noAF_pred"         = "BayesDel (without MaxAF) prediction",
           "ClinPred_pred"              = "ClinPred prediction",
           "LIST-S2_pred"               = "LIST-S2 prediction",
           "ESM1b_pred"                 = "ESM-1b (protein language model) prediction",
           "AlphaMissense_pred"         = "AlphaMissense prediction",
           "Aloft_pred"                 = "ALoFT (Annotation of Loss-of-Function Transcripts) classification",
           "fathmm-XF_coding_pred"      = "FATHMM-XF (coding) prediction",
           "fathmm-MKL_coding_pred"     = "FATHMM-MKL (coding) prediction",
           "FATHMM_pred"                = "FATHMM prediction",
           "HIPred"                     = "HIPred haploinsufficiency prediction",
           "Gene_indispensability_pred" = "Gene indispensability (essentiality) prediction",
           col
    )
  }
  
  # -----------------------------
  # Helpers: item meaning (dbNSFP 4.1c conventions + extensions used in your data)
  # -----------------------------
  pred_item_meaning <- function(column, item) {
    col <- as.character(column)
    it  <- as.character(item)
    
    # normalize common aliases you might have in df/universe
    if (col == "M_CAP_pred") col <- "M-CAP_pred"
    if (col == "LIST_S2_pred") col <- "LIST-S2_pred"
    if (col == "fathmm_XF_coding_pred") col <- "fathmm-XF_coding_pred"
    
    # SIFT / SIFT4G: D/T
    if (col %in% c("SIFT_pred","SIFT4G_pred")) {
      return(switch(it, "D"="Damaging", "T"="Tolerated", NA_character_))
    }
    
    # PolyPhen2: B/P/D
    if (grepl("^Polyphen2_", col, ignore.case = TRUE)) {
      return(switch(it,
                    "B"="Benign",
                    "P"="Possibly damaging",
                    "D"="Probably damaging",
                    NA_character_
      ))
    }
    
    # PROVEAN: D/N
    if (col == "PROVEAN_pred") {
      return(switch(it, "D"="Damaging", "N"="Neutral", NA_character_))
    }
    
    # LRT: D/N/U
    if (col == "LRT_pred") {
      return(switch(it, "D"="Deleterious", "N"="Neutral", "U"="Unknown", NA_character_))
    }
    
    # MutationTaster: A/D/N/P
    if (col == "MutationTaster_pred") {
      return(switch(it,
                    "A"="Disease causing (automatic)",
                    "D"="Disease causing",
                    "N"="Polymorphism",
                    "P"="Polymorphism (automatic)",
                    NA_character_
      ))
    }
    
    # MutationAssessor: H/M/L/N
    if (col == "MutationAssessor_pred") {
      return(switch(it,
                    "H"="High impact",
                    "M"="Medium impact",
                    "L"="Low impact",
                    "N"="Neutral",
                    NA_character_
      ))
    }
    
    # FATHMM: D/T
    if (col == "FATHMM_pred") {
      return(switch(it, "D"="Damaging", "T"="Tolerated", NA_character_))
    }
    
    # fathmm-MKL coding: D/N
    if (col == "fathmm-MKL_coding_pred") {
      return(switch(it, "D"="Damaging", "N"="Neutral", NA_character_))
    }
    
    # fathmm-XF coding: D/N
    if (col == "fathmm-XF_coding_pred") {
      return(switch(it, "D"="Damaging", "N"="Neutral", NA_character_))
    }
    
    # MetaSVM/MetaLR/MetaRNN/M-CAP/PrimateAI/DEOGEN2/BayesDel/ClinPred/LIST-S2/ESM1b: D/T (binary)
    if (col %in% c("MetaSVM_pred","MetaLR_pred","MetaRNN_pred",
                   "M-CAP_pred","PrimateAI_pred","DEOGEN2_pred",
                   "BayesDel_addAF_pred","BayesDel_noAF_pred",
                   "ClinPred_pred","LIST-S2_pred","ESM1b_pred")) {
      return(switch(it, "D"="Damaging / deleterious", "T"="Tolerated", NA_character_))
    }
    
    # AlphaMissense: P/A/B and sometimes LP/LB
    if (col == "AlphaMissense_pred") {
      return(switch(it,
                    "P" ="Pathogenic",
                    "LP"="Likely pathogenic",
                    "A" ="Ambiguous",
                    "B" ="Benign",
                    "LB"="Likely benign",
                    NA_character_
      ))
    }
    
    # Aloft: textual
    if (col == "Aloft_pred") {
      if (it %in% c("Tolerant","Recessive","Dominant")) return(it)
      return(NA_character_)
    }
    
    # HIPred: Y/N
    if (col == "HIPred") {
      return(switch(it, "Y"="Yes (haploinsufficient)", "N"="No", NA_character_))
    }
    
    # Gene indispensability: E/N
    if (col == "Gene_indispensability_pred") {
      return(switch(it, "E"="Essential", "N"="LoF-tolerant (non-essential)", NA_character_))
    }
    
    # fallback for generic D/T
    if (it %in% c("D","T")) return(ifelse(it == "D", "Damaging", "Tolerated"))
    
    NA_character_
  }
  
  # -----------------------------
  omim_title_map <- local({
    .cache <- NULL
    function(path = "/Volumes/DISK1TB/Inspector_app_slaves_ngroc/GItools/_shared/omim_mim_titles_map.rds") {
      if (!is.null(.cache)) return(.cache)
      if (!file.exists(path)) {
        .cache <- stats::setNames(character(0), character(0))
        return(.cache)
      }
      .cache <- readRDS(path)
      .cache
    }
  })
  
  mim_termid_to_num <- function(x) sub("^MIM:", "", as.character(x))
  
  # -----------------------------
  # Universe loader
  # -----------------------------
  load_predictor_universe <- function(rds) {
    validate(need(file.exists(rds), paste0("Missing predictor universe: ", rds)))
    u <- readRDS(rds)
    u <- data.table::as.data.table(u)
    validate(need(all(c("column","item","N") %in% names(u)), "Universe must have column,item,N"))
    u <- u[!is.na(item)]
    u <- u[item != column]  # remove header artefact counts (e.g. 22)
    u
  }
  
  #-------------------------------------------
  #-------- Calssifier -------------------
  #-------------------------------------------
  
  pred_item_class <- function(column, item) {
    col <- as.character(column)
    it  <- as.character(item)
    
    # normalize common aliases
    if (col == "M_CAP_pred") col <- "M-CAP_pred"
    if (col == "LIST_S2_pred") col <- "LIST-S2_pred"
    if (col == "fathmm_XF_coding_pred") col <- "fathmm-XF_coding_pred"
    
    # Binary D/T predictors: D=Bad, T=Good
    if (col %in% c("SIFT_pred","SIFT4G_pred","MetaSVM_pred","MetaLR_pred","MetaRNN_pred",
                   "M-CAP_pred","PrimateAI_pred","DEOGEN2_pred",
                   "BayesDel_addAF_pred","BayesDel_noAF_pred",
                   "ClinPred_pred","LIST-S2_pred","ESM1b_pred")) {
      if (it == "D") return("Bad")
      if (it == "T") return("Good")
      return("Unknown")
    }
    
    # PolyPhen2: B=Good, P=Neutral, D=Bad
    if (grepl("^Polyphen2_", col, ignore.case = TRUE)) {
      if (it == "B") return("Good")
      if (it == "P") return("Neutral")
      if (it == "D") return("Bad")
      return("Unknown")
    }
    
    # PROVEAN: N=Good, D=Bad
    if (col == "PROVEAN_pred") {
      if (it == "N") return("Good")
      if (it == "D") return("Bad")
      return("Unknown")
    }
    
    # MutationTaster: N/P = Good, A/D = Bad
    if (col == "MutationTaster_pred") {
      if (it %in% c("N","P")) return("Good")
      if (it %in% c("A","D")) return("Bad")
      return("Unknown")
    }
    
    # MutationAssessor: N/L = Good, M/H = Bad
    if (col == "MutationAssessor_pred") {
      if (it %in% c("N","L")) return("Good")
      if (it %in% c("M","H")) return("Bad")
      return("Unknown")
    }
    
    # fathmm-XF coding: N=Good, D=Bad
    if (col == "fathmm-XF_coding_pred") {
      if (it == "N") return("Good")
      if (it == "D") return("Bad")
      return("Unknown")
    }
    
    # AlphaMissense: B/LB=Good, A=Neutral, P/LP=Bad
    if (col == "AlphaMissense_pred") {
      if (it %in% c("B","LB")) return("Good")
      if (it == "A") return("Neutral")
      if (it %in% c("P","LP")) return("Bad")
      return("Unknown")
    }
    
    # Aloft: Tolerant=Good; Dominant/Recessive=Neutral (configurable) 
    if (col == "Aloft_pred") {
      if (it == "Tolerant") return("Good")
      if (it %in% c("Dominant","Recessive")) return("Neutral")
      return("Unknown")
    }
    
    # HIPred: N=Good, Y=Bad (haploinsufficient)
    if (col == "HIPred") {
      if (it == "N") return("Good")
      if (it == "Y") return("Bad")
      return("Unknown")
    }
    
    # Gene indispensability: N=Good (LoF-tolerant), E=Bad (essential)
    if (col == "Gene_indispensability_pred") {
      if (it == "N") return("Good")
      if (it == "E") return("Bad")
      return("Unknown")
    }
    
    "Unknown"
  }
  
  pred_item_sign <- function(column, item) {
    col <- as.character(column)
    it  <- as.character(item)
    
    # Normalitza aliases típics
    if (col == "M_CAP_pred") col <- "M-CAP_pred"
    if (col == "LIST_S2_pred") col <- "LIST-S2_pred"
    if (col == "fathmm_XF_coding_pred") col <- "fathmm-XF_coding_pred"
    
    # --- mappings explícits per predictor ---
    if (col %in% c("SIFT_pred","SIFT4G_pred","MetaSVM_pred","MetaLR_pred","MetaRNN_pred",
                   "M-CAP_pred","PrimateAI_pred","DEOGEN2_pred","BayesDel_addAF_pred",
                   "BayesDel_noAF_pred","ClinPred_pred","LIST-S2_pred","ESM1b_pred")) {
      return(ifelse(it == "D", "Negative", ifelse(it == "T", "Positive", "Unknown")))
    }
    
    if (grepl("^Polyphen2_", col, ignore.case = TRUE)) {
      return(ifelse(it %in% c("D","P"), "Negative", ifelse(it == "B", "Positive", "Unknown")))
    }
    
    if (col == "PROVEAN_pred") {
      return(ifelse(it == "D", "Negative", ifelse(it == "N", "Positive", "Unknown")))
    }
    
    if (col == "MutationTaster_pred") {
      return(ifelse(it %in% c("A","D"), "Negative", ifelse(it %in% c("N","P"), "Positive", "Unknown")))
    }
    
    if (col == "MutationAssessor_pred") {
      return(ifelse(it %in% c("H","M"), "Negative", ifelse(it %in% c("L","N"), "Positive", "Unknown")))
    }
    
    if (col %in% c("fathmm-XF_coding_pred","fathmm-MKL_coding_pred")) {
      return(ifelse(it == "D", "Negative", ifelse(it == "N", "Positive", "Unknown")))
    }
    
    if (col == "AlphaMissense_pred") {
      return(ifelse(it %in% c("P","LP"), "Negative", ifelse(it %in% c("B","LB"), "Positive", "Unknown")))
    }
    
    if (col == "Aloft_pred") {
      # interpretació habitual: Dominant/Recessive = LoF disease-causing (negatiu), Tolerant = positiu
      return(ifelse(it %in% c("Dominant","Recessive"), "Negative",
                    ifelse(it == "Tolerant", "Positive", "Unknown")))
    }
    
    if (col == "HIPred") {
      # Y = haploinsufficient (negatiu), N = positiu
      return(ifelse(it == "Y", "Negative", ifelse(it == "N", "Positive", "Unknown")))
    }
    
    if (col == "Gene_indispensability_pred") {
      # E = essential (interpretem com “negatiu” per tolerància a LoF), N = tolerat
      return(ifelse(it == "E", "Negative", ifelse(it == "N", "Positive", "Unknown")))
    }
    
    "Unknown"
  }
  # -----------------------------
  # Helpers: parse multi-transcript predictor fields
  # - split ';' always
  # - optionally split '|' (Aloft combos like Recessive|Recessive)
  # -----------------------------
  pred_split_tokens <- function(v, split_pipe = TRUE) {
    v <- as.character(v)
    v <- v[!is.na(v)]
    v <- trimws(v)
    v <- v[nzchar(v)]
    if (!length(v)) return(character(0))
    
    toks <- unlist(strsplit(v, ";", fixed = TRUE), use.names = FALSE)
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    toks <- toks[!(toks %in% c(".", "NA", ""))]
    
    if (split_pipe && length(toks)) {
      toks <- unlist(strsplit(toks, "|", fixed = TRUE), use.names = FALSE)
      toks <- trimws(toks)
      toks <- toks[nzchar(toks)]
      toks <- toks[!(toks %in% c(".", "NA", ""))]
    }
    
    toks
  }
  
  pred_count_items_fg <- function(dt, pred_cols, split_pipe = TRUE) {
    data.table::rbindlist(lapply(pred_cols, function(cc) {
      if (!cc %in% names(dt)) return(NULL)
      toks <- pred_split_tokens(dt[[cc]], split_pipe = split_pipe)
      if (!length(toks)) return(NULL)
      out <- data.table::data.table(item = toks)[, .N, by = item]
      out[, column := cc]
      data.table::setcolorder(out, c("column","item","N"))
      out
    }), use.names = TRUE, fill = TRUE)
  }
  
  # -----------------------------
  # ORA for predictors from counts using Fisher exact (greater)
  # Inputs:
  # - fg_counts: (column,item,N) where N is count in foreground
  # - bg_counts: (column,item,N) where N is count in background/universe
  # Output:
  # - columns: column,item,a,Tt,n_total,N_total,or,p,p_adj, plus helper cols
  # -----------------------------
  predictor_ora <- function(fg_counts, bg_counts, min_a = 1L) {
    
    if (!is.data.frame(fg_counts) || !nrow(fg_counts)) return(data.table::data.table())
    if (!is.data.frame(bg_counts) || !nrow(bg_counts)) return(data.table::data.table())
    
    fg <- data.table::as.data.table(fg_counts)
    bg <- data.table::as.data.table(bg_counts)
    
    fg <- fg[N >= min_a]
    if (!nrow(fg)) return(data.table::data.table())
    
    fg_tot <- fg[, .(n_total = sum(N)), by = column]
    bg_tot <- bg[, .(N_total = sum(N)), by = column]
    
    fg2 <- fg[, .(column, item, a = as.integer(N))]
    bg2 <- bg[, .(column, item, Tt = as.integer(N))]
    
    m <- merge(fg2, bg2, by = c("column","item"), all.x = TRUE)
    m[is.na(Tt), Tt := 0L]
    
    m <- merge(m, fg_tot, by = "column", all.x = TRUE)
    m <- merge(m, bg_tot, by = "column", all.x = TRUE)
    
    m[, `:=`(
      b = as.integer(pmax(n_total - a, 0)),
      c = as.integer(pmax(Tt - a, 0)),
      d = as.integer(pmax((N_total - Tt) - (n_total - a), 0))
    )]
    
    m[, p := mapply(function(a,b,c,d) {
      if ((a+b+c+d) == 0) return(NA_real_)
      stats::fisher.test(matrix(c(a,c,b,d), nrow = 2), alternative = "greater")$p.value
    }, a,b,c,d)]
    
    m[, or := mapply(function(a,b,c,d) {
      ft <- suppressWarnings(stats::fisher.test(matrix(c(a,c,b,d), nrow = 2), alternative = "greater"))
      unname(ft$estimate)
    }, a,b,c,d)]
    
    m[, p_adj := stats::p.adjust(p, method = "BH")]
    
    data.table::setorder(m, p_adj, p)
    m
  }
  
  # ============================================================
  # EVENT: run enrichment (only when Predictors tab is active)
  # NOTE: dt_enrich_res_nonsyn() already applies scope global/cluster/gene
  # ============================================================
  nonsyn_predictors_res <- eventReactive(input$run_enrich, {
    req(input$enrich_tabs == "tab_enrich_pred")
    
    x <- dt_enrich_res_nonsyn()
    dt2 <- x$dt2
    validate(need(is.data.frame(dt2) && nrow(dt2) > 0, "No NonSyn hits in this scope."))
    
    bg <- load_predictor_universe(PRED_UNIV_RDS)
    
    pred_cols_bg <- unique(as.character(bg$column))
    pred_cols_fg <- intersect(pred_cols_bg, names(dt2))
    validate(need(length(pred_cols_fg) > 0, "No *_pred columns found in foreground (dt2)."))
    
    split_pipe <- isTRUE(input$pred_split_pipe)
    
    topk  <- suppressWarnings(as.integer(input$pred_topk %||% 25))
    if (!is.finite(topk) || topk < 5L) topk <- 25L
    
    min_a <- suppressWarnings(as.integer(input$pred_min_a %||% 2))
    if (!is.finite(min_a) || min_a < 1L) min_a <- 1L
    
    withProgress(message = "Running predictor enrichment…", value = 0, {
      
      incProgress(0.30, detail = "Foreground parsing…")
      fg_counts <- pred_count_items_fg(dt2, pred_cols_fg, split_pipe = split_pipe)
      validate(need(is.data.frame(fg_counts) && nrow(fg_counts) > 0,
                    "No predictor items parsed in foreground (all NA/'.'?)."))
      
      incProgress(0.55, detail = "Fisher ORA…")
      res <- predictor_ora(fg_counts, bg, min_a = min_a)
      
      incProgress(0.15, detail = "Done.")
      list(tab = res, topk = topk, min_a = min_a, split_pipe = split_pipe)
    })
  }, ignoreInit = TRUE)
  
  # ============================================================
  # OUTPUT: TABLE (DT) with full names + item meanings
  # ============================================================
  output$enrich_nonsyn_pred_table <- DT::renderDT({
    r <- nonsyn_predictors_res(); req(r)
    tab <- r$tab
    
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(DT::datatable(
        data.frame(Message = "No enriched predictor items (try lowering min_a)."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    topk <- r$topk %||% 25
    
    df <- tab %>%
      dplyr::mutate(
        a_num   = suppressWarnings(as.numeric(as.character(a))),
        n_num   = suppressWarnings(as.numeric(as.character(n_total))),
        Tt_num  = suppressWarnings(as.numeric(as.character(Tt))),
        N_num   = suppressWarnings(as.numeric(as.character(N_total))),
        p_num   = suppressWarnings(as.numeric(as.character(p))),
        fdr_num = suppressWarnings(as.numeric(as.character(p_adj)))
      ) %>%
      dplyr::filter(is.finite(a_num), is.finite(Tt_num), is.finite(p_num), is.finite(fdr_num)) %>%
      dplyr::arrange(fdr_num, p_num) %>%
      dplyr::slice_head(n = topk) %>%
      dplyr::mutate(
        Predictor_full = vapply(column, pred_full_name, character(1)),
        Meaning        = mapply(pred_item_meaning, column, item),
        `In subset`    = sprintf("%d/%d (%.2f%%)", a_num, n_num, 100 * a_num / pmax(1, n_num)),
        `In bg`        = sprintf("%d/%d (%.2f%%)", Tt_num, N_num, 100 * Tt_num / pmax(1, N_num)),
        OR             = formatC(as.numeric(or), format = "f", digits = 3),
        p              = formatC(p_num,   format = "e", digits = 2),
        FDR            = formatC(fdr_num, format = "e", digits = 2)
      ) %>%
      dplyr::transmute(
        Predictor_code = as.character(column),
        Predictor_full = as.character(Predictor_full),
        Item           = as.character(item),
        Meaning        = as.character(Meaning),
        `In subset`, `In bg`,
        OR, p, FDR
      )
    
    DT::datatable(
      df, rownames = FALSE,
      extensions = "Buttons",
      options = list(
        extensions = "Buttons",
        dom = "Bfrtip",
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength = 15,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
  
  # ============================================================
  # OUTPUT: PLOT
  # - X visible = counts (a) via coord_flip()
  # - Fill = −log10(FDR) with yellow→orange→red
  # - Labels include predictor full name + item meaning
  # ============================================================
  
  pred_plot_df <- reactive({
    r <- nonsyn_predictors_res(); req(r)
    tab <- r$tab
    req(is.data.frame(tab), nrow(tab) > 0)
    
    tab %>%
      dplyr::mutate(
        a_num   = suppressWarnings(as.numeric(as.character(a))),
        fdr_num = suppressWarnings(as.numeric(as.character(p_adj)))
      ) %>%
      dplyr::filter(is.finite(a_num), is.finite(fdr_num)) %>%
      dplyr::arrange(fdr_num, dplyr::desc(a_num)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(
        neglogFDR = -log10(pmax(fdr_num, 1e-300)),
        cls = mapply(pred_item_class, column, item),
        cls = factor(cls, levels = c("Good","Neutral","Bad","Unknown")),
        Predictor_full = vapply(column, pred_full_name, character(1)),
        Meaning        = mapply(pred_item_meaning, column, item),
        term_show = dplyr::if_else(
          !is.na(Meaning) & nzchar(Meaning),
          paste0(Predictor_full, " (", column, ") : ", item, " (", Meaning, ")"),
          paste0(Predictor_full, " (", column, ") : ", item)
        ),
        row_id = dplyr::row_number()
      )
  }) 
  
  pred_click <- reactive({
    ed <- plotly::event_data("plotly_click", source = "pred_enrich")
    if (is.null(ed) || !nrow(ed)) return(NULL)
    if (is.null(ed$customdata)) return(NULL)
    as.data.frame(ed$customdata, stringsAsFactors = FALSE)
  })
  
  output$enrich_nonsyn_pred_click_table <- DT::renderDT({
    cd <- pred_click()
    if (is.null(cd)) {
      return(DT::datatable(
        data.frame(Message = "Click a bar to see details."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    cd$FDR <- formatC(as.numeric(cd$FDR), format = "e", digits = 2)
    cd$p   <- formatC(as.numeric(cd$p),   format = "e", digits = 2)
    cd$OR  <- formatC(as.numeric(cd$OR),  format = "f", digits = 3)
    cd$neglogFDR <- round(as.numeric(cd$neglogFDR), 3)
    
    DT::datatable(cd, 
                  rownames = FALSE, 
                  options = list(dom = "t", 
                                 extensions = "Buttons",
                                 scrollX = TRUE,
                                 buttons = list(
                                   list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                                   list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                                   list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                                 )
                  )
    )
  }, server = FALSE)
  
  
  
  output$enrich_nonsyn_pred_click_dbg <- renderPrint({
    ed <- plotly::event_data("plotly_click", source = "pred_enrich")
    if (is.null(ed) || !nrow(ed)) {
      cat("Click a bar to print raw plotly_click event_data().\n")
      return()
    }
    print(ed)
  })
  
  short_label_mid <- function(x, max_chars = 42) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("\\s*\\(.*?\\)\\s*", " ", x)  # treu parèntesis (acronims llargs)
    x <- gsub("\\s+prediction\\b", "", x, ignore.case = TRUE)
    x <- gsub("\\s+classification\\b", "", x, ignore.case = TRUE)
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)
    ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
  }
  
  output$enrich_nonsyn_pred_plotly <- plotly::renderPlotly({
    r <- nonsyn_predictors_res(); req(r)
    tab <- r$tab
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(plotly::plotly_empty(type = "bar") %>% plotly::layout(title = "No data"))
    }
    
    df <- tab %>%
      dplyr::mutate(
        a_num   = suppressWarnings(as.numeric(as.character(a))),
        fdr_num = suppressWarnings(as.numeric(as.character(p_adj)))
      ) %>%
      dplyr::filter(is.finite(a_num), is.finite(fdr_num)) %>%
      dplyr::arrange(fdr_num, dplyr::desc(a_num)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(
        neglogFDR = -log10(pmax(fdr_num, 1e-300)),
        cls = mapply(pred_item_class, column, item),
        cls = factor(cls, levels = c("Good","Neutral","Bad","Unknown")),
        Predictor_full = vapply(column, pred_full_name, character(1)),
        Meaning        = mapply(pred_item_meaning, column, item),
        term_show = dplyr::if_else(
          !is.na(Meaning) & nzchar(Meaning),
          paste0(Predictor_full, " (", column, ") : ", item, " (", Meaning, ")"),
          paste0(Predictor_full, " (", column, ") : ", item)
        )
      )
    
    # base colors per classe
    base_col <- c(Good="darkgreen", Neutral="gold", Bad="red", Unknown="grey70")
    df$base_color <- unname(base_col[as.character(df$cls)])
    df$base_color[is.na(df$base_color)] <- "grey70"
    
    # opacity per barra segons -logFDR (0.35..1.0)
    rng <- range(df$neglogFDR, finite = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
      df$opacity <- 1.0
    } else {
      z <- (df$neglogFDR - rng[1]) / (rng[2] - rng[1])
      df$opacity <- 0.35 + z * (1.0 - 0.35)
    }
    
    # ordre Y (més a dalt = més significatiu)
    df <- df %>% dplyr::arrange(fdr_num, dplyr::desc(a_num))
    df$term_show_full <- as.character(df$term_show)
    df$term_show_full <- factor(df$term_show_full, levels = rev(unique(df$term_show_full)))
    
    # etiqueta intermèdia per eix: "Nom curt (code) : item"
    pred_mid <- short_label_mid(df$Predictor_full, max_chars = 50)
    df$term_show_mid <- paste0(pred_mid, " (", df$column, ") : ", df$item)
    
    # tick mapping (yaxis)
    tickvals <- levels(df$term_show_full)
    ticktext <- df$term_show_mid[match(tickvals, df$term_show_full)]
    
    # customdata per click
    custom <- df %>%
      dplyr::transmute(
        Predictor_code = as.character(column),
        Predictor_full = as.character(Predictor_full),
        Item           = as.character(item),
        Meaning        = as.character(Meaning),
        Class          = as.character(cls),
        Count_a        = as.integer(a_num),
        FDR            = as.numeric(fdr_num),
        neglogFDR      = as.numeric(neglogFDR),
        OR             = as.numeric(or),
        p              = as.numeric(p)
      )
    
    # Tooltip (hover)
    txt <- paste0(
      "<b>", as.character(df$term_show), "</b>",
      "<br>Class: ", as.character(df$cls),
      "<br>Count (a): ", df$a_num,
      "<br>FDR: ", formatC(df$fdr_num, format="e", digits=2),
      "<br>-log10(FDR): ", round(df$neglogFDR, 3)
    )
    
    plotly::plot_ly(
      data = df,
      x = ~a_num,
      y = ~term_show_full,
      type = "bar",
      orientation = "h",
      source = "pred_enrich",
      text = txt,
      hoverinfo = "text",
      customdata = custom_list
    ) %>%
      plotly::style(
        textposition = "none",
        marker = list(
          color = df$base_color,
          opacity = df$opacity,
          line = list(width = 0)
        )
      ) %>%
      plotly::layout(
        title = "Top predictor items (X = counts; color = class; intensity = −log10(FDR))",
        xaxis = list(title = "Count of NonSynonym hits (a)"),
        yaxis = list(
          title = "",
          tickmode = "array",
          tickvals = tickvals,
          ticktext = ticktext
        ),
        margin = list(l = 340)
      )
  })
  
  
  # Helper: convert named color + alpha to RGBA hex string
  to_rgba <- function(col, alpha) {
    rgb <- grDevices::col2rgb(col)
    sprintf("rgba(%d,%d,%d,%.3f)", rgb[1,1], rgb[2,1], rgb[3,1], alpha)
  }
  
  output$enrich_nonsyn_pred_plotly_radial <- plotly::renderPlotly({
    r <- nonsyn_predictors_res(); req(r)
    tab <- r$tab
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(plotly::plotly_empty(type = "barpolar") %>% plotly::layout(title = "No data"))
    }
    
    df <- tab %>%
      dplyr::mutate(
        a_num   = suppressWarnings(as.numeric(as.character(a))),
        fdr_num = suppressWarnings(as.numeric(as.character(p_adj)))
      ) %>%
      dplyr::filter(is.finite(a_num), is.finite(fdr_num)) %>%
      dplyr::arrange(fdr_num, dplyr::desc(a_num)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(
        neglogFDR = -log10(pmax(fdr_num, 1e-300)),
        cls = mapply(pred_item_class, column, item),
        cls = factor(cls, levels = c("Good","Neutral","Bad","Unknown")),
        Predictor_full = vapply(column, pred_full_name, character(1)),
        Meaning        = mapply(pred_item_meaning, column, item),
        term_show_full = dplyr::if_else(
          !is.na(Meaning) & nzchar(Meaning),
          paste0(Predictor_full, " (", column, ") : ", item, " (", Meaning, ")"),
          paste0(Predictor_full, " (", column, ") : ", item)
        )
      )
    
    # etiqueta intermèdia per l’eix circular (theta)
    pred_mid <- short_label_mid(df$Predictor_full, max_chars = 34)
    df$theta_lab <- paste0(pred_mid, " (", df$column, ") : ", df$item)
    
    # colors base per classe
    base_col <- c(Good="darkgreen", Neutral="gold", Bad="red", Unknown="grey70")
    df$base_color <- unname(base_col[as.character(df$cls)])
    df$base_color[is.na(df$base_color)] <- "grey70"
    
    # alpha per barra segons -logFDR (0.35..1.0)
    rng <- range(df$neglogFDR, finite = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
      df$alpha <- 1.0
    } else {
      z <- (df$neglogFDR - rng[1]) / (rng[2] - rng[1])
      df$alpha <- 0.35 + z * (1.0 - 0.35)
    }
    
    # RGBA per barra (permet alpha variable per element)
    df$rgba <- mapply(to_rgba, df$base_color, df$alpha)
    
    # customdata per click
    custom <- df %>%
      dplyr::transmute(
        Predictor_code = as.character(column),
        Predictor_full = as.character(Predictor_full),
        Item           = as.character(item),
        Meaning        = as.character(Meaning),
        Class          = as.character(cls),
        Count_a        = as.integer(a_num),
        FDR            = as.numeric(fdr_num),
        neglogFDR      = as.numeric(neglogFDR),
        OR             = as.numeric(or),
        p              = as.numeric(p)
      )
    
    # IMPORTANT: convertir a llista per punt (longitud = nrow(df))
    custom_list <- lapply(seq_len(nrow(custom)), function(i) as.list(custom[i, , drop = FALSE]))
    
    # Tooltip (només hover/click; no es dibuixa dins la barra)
    txt <- paste0(
      "<b>", df$term_show_full, "</b>",
      "<br>Class: ", as.character(df$cls),
      "<br>Count (a): ", df$a_num,
      "<br>FDR: ", formatC(df$fdr_num, format="e", digits=2),
      "<br>-log10(FDR): ", round(df$neglogFDR, 3)
    )
    
    plotly::plot_ly(
      data = df,
      type = "barpolar",
      r = ~a_num,
      theta = ~theta_lab,
      text = txt,
      hoverinfo = "text",
      customdata = custom_list,
      source = "pred_enrich",
      marker = list(
        color = df$rgba,
        line  = list(color = "white", width = 1)
      )
    ) %>%
      plotly::layout(
        title = "Top 15 predictor items (radial) — r = Nonsyn hits count; color = class; intensity = −log10(FDR))",
        polar = list(
          radialaxis = list(title = "", ticksuffix = " "),
          angularaxis = list(direction = "clockwise")
        ),
        margin = list(l = 40, r = 40, t = 60, b = 40)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "enrichment_nonsyn_predictors",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  # ============================================================
  # DOWNLOAD: table
  # ============================================================
  output$pred_download <- downloadHandler(
    filename = function() paste0("enrichment_predictors_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) {
      r <- nonsyn_predictors_res(); req(r)
      tab <- r$tab
      if (!is.data.frame(tab) || !nrow(tab)) {
        data.table::fwrite(data.table::data.table(Message="No results"), file, sep = "\t")
        return()
      }
      
      df <- as.data.frame(tab, stringsAsFactors = FALSE)
      df$Predictor_full <- vapply(df$column, pred_full_name, character(1))
      df$Meaning        <- mapply(pred_item_meaning, df$column, df$item)
      
      data.table::fwrite(df, file, sep = "\t")
    }
  )
  
  ############################################################################
  # Func Disease Table
  ############################################################################
  
  func_disease_df <- reactive({
    
    path_rds <- dbnsfp_final_path_rds()
    path_csv <- dbnsfp_final_path_csv()
    
    df <- NULL
    if (!is.null(path_rds) && file.exists(path_rds) && file.info(path_rds)$size > 0) {
      df <- readRDS(path_rds)
    } else if (!is.null(path_csv) && file.exists(path_csv) && file.info(path_csv)$size > 0) {
      df <- utils::read.csv(path_csv, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    validate(
      need(is.data.frame(df) && nrow(df) > 0, "Final output not available yet (build + normalize + clusters first)."),
      need(all(c("SNP","genename","cluster_id","Function_description","Disease_description",
                 "Orphanet_disorder","MIM_disease","clinvar_trait",
                 "cluster_chr_n","cluster_start","cluster_end") %in% names(df)),
           "Missing required cluster/function/disease columns in final output.")
    )
    
    omim_col <- dplyr::case_when(
      "OMIM_id" %in% names(df) ~ "OMIM_id",
      "MIM_phenotype_id" %in% names(df) ~ "MIM_phenotype_id",
      TRUE ~ NA_character_
    )
    validate(need(!is.na(omim_col), "Missing OMIM column (expected OMIM_id or MIM_phenotype_id)."))
    
    make_ucsc_url <- function(chr_n, start, end, db = "hg38") {
      if (is.na(chr_n) || is.na(start) || is.na(end)) return(NA_character_)
      chr <- as.character(chr_n)
      if (!startsWith(chr, "chr")) chr <- paste0("chr", chr)
      pos <- paste0(chr, ":", as.integer(start), "-", as.integer(end))
      paste0(
        "https://genome.ucsc.edu/cgi-bin/hgTracks?db=", db,
        "&position=", utils::URLencode(pos, reserved = TRUE)
      )
    }
    
    make_genecards_url <- function(g) {
      if (is.null(g) || is.na(g)) return(NA_character_)
      s <- trimws(as.character(g))
      if (!nzchar(s) || s %in% c(".", "NA")) return(NA_character_)
      paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", utils::URLencode(s, reserved = TRUE))
    }
    
    # extreu el primer MIM:#### de MIM_disease (si existeix) i genera URL OMIM
    omim_from_mim_disease <- function(x) {
      s <- trimws(as.character(x %||% ""))
      if (!nzchar(s) || s %in% c(".", "NA")) return(NA_character_)
      id <- regmatches(s, regexpr("MIM:\\d+", s, perl = TRUE))
      if (!length(id) || !nzchar(id)) return(NA_character_)
      paste0("https://omim.org/entry/", gsub("^MIM:", "", id))
    }
    
    make_dbsnp_url <- function(snp) {
      s <- trimws(as.character(snp))
      if (!nzchar(s) || s %in% c(".", "NA")) return(NA_character_)
      # dbSNP: només per rsIDs (si tens chr:pos també es podria fer UCSC, però aquí fem dbSNP)
      if (grepl("^rs\\d+$", s, ignore.case = TRUE)) {
        paste0("https://www.ncbi.nlm.nih.gov/snp/", s)
      } else {
        NA_character_
      }
    }
    
    # ---- Multi-link helpers (MIM + Orphanet) ----
    
    strip_orphanet_prefix <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      x <- trimws(x)
      x <- sub("^Orphanet:\\s*", "", x, ignore.case = TRUE)
      x
    }
    
    # split by ';' but keep tokens
    split_semicol <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      x <- trimws(x)
      if (!nzchar(x)) return(character(0))
      v <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
      v <- trimws(v)
      v[nzchar(v)]
    }
    
    # Build OMIM links for a single MIM_disease cell (possibly multi ';')
    link_mim_cell <- function(x) {
      parts <- split_semicol(x)
      if (!length(parts)) return(as.character(x))
      
      out <- vapply(parts, function(p) {
        id <- regmatches(p, regexpr("MIM:\\d+", p, perl = TRUE))
        lab <- trimws(p)
        if (!length(id) || !nzchar(id)) return(lab)
        
        url <- paste0("https://omim.org/entry/", sub("^MIM:", "", id))
        sprintf("<a href='%s' target='_blank'>%s</a>", url, lab)
      }, character(1))
      
      paste(out, collapse = "<br>")
    }
    
    # ---- ClinVar helpers ----
    
    # ClinVar Variation link
    make_clinvar_variation_url <- function(variation_id, rsid) {
      vid <- trimws(as.character(variation_id %||% ""))
      rs  <- trimws(as.character(rsid %||% ""))
      
      vid_num <- gsub("[^0-9]", "", vid)
      if (nzchar(vid_num)) {
        q <- utils::URLencode(paste0('"', rs, '"'), reserved = TRUE)
        return(paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", vid_num, "/?term=", q))
      }
      
      # fallback: search by rs
      if (grepl("^rs\\d+$", rs, ignore.case = TRUE)) {
        q <- utils::URLencode(paste0('"', rs, '"'), reserved = TRUE)
        return(paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=", q))
      }
      
      NA_character_
    }
    
    # MedGen link (ClinVar trait)
    make_medgen_url <- function(medgen_id) {
      s <- trimws(as.character(medgen_id %||% ""))
      if (!nzchar(s)) return(NA_character_)
      # expect C########
      if (!grepl("^C\\d+$", s)) return(NA_character_)
      paste0("https://www.ncbi.nlm.nih.gov/medgen/", s)
    }
    
    # Split by ';' (your current convention)
    split_semicol <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      x <- trimws(x)
      if (!nzchar(x)) return(character(0))
      v <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
      v <- trimws(v)
      v[nzchar(v)]
    }
    
    # Multi-link for clinvar traits with aligned medgen ids (if provided)
    link_clinvar_trait_cell <- function(trait_cell, medgen_cell = NA_character_) {
      traits <- split_semicol(trait_cell)
      if (!length(traits)) return(as.character(trait_cell))
      
      ids <- split_semicol(medgen_cell)
      # align lengths (if 1 id only, recycle it)
      if (length(ids) == 1L && length(traits) > 1L) ids <- rep(ids, length(traits))
      if (length(ids) < length(traits)) ids <- c(ids, rep("", length(traits) - length(ids)))
      ids <- ids[seq_len(length(traits))]
      
      out <- mapply(function(tr, id) {
        tr <- trimws(tr)
        u <- make_medgen_url(id)
        if (is.na(u) || !nzchar(u)) return(tr)
        sprintf("<a href='%s' target='_blank'>%s</a>", u, tr)
      }, traits, ids, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      
      paste(out, collapse = "<br>")
    }
    # eliminar valors "escombreria"
    drop_bad_tokens <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      x <- trimws(x)
      x <- x[nzchar(x)]
      low <- tolower(x)
      bad <- low %in% c(".", "na", "n/a", "nan", "null", "none", "not_provided", "not provided", "not_specified", "not specified", "unspecified", "unknown")
      x <- x[!bad]
      x[nzchar(x)]
    }
    
    # Build Orphanet links for a single cell + optional id cell (both may be multi ';')
    link_orpha_cell <- function(name_cell, id_cell = NA_character_) {
      nm <- split_semicol(name_cell)
      if (!length(nm)) return(as.character(name_cell))
      
      ids <- split_semicol(id_cell)
      # if ids not provided, try to extract ORPHA:#### from name
      if (!length(ids)) {
        ids <- vapply(nm, function(p) {
          m <- regmatches(p, regexpr("ORPHA:\\d+", p, perl = TRUE))
          if (!length(m) || !nzchar(m)) "" else sub("^ORPHA:", "", m)
        }, character(1))
      } else {
        # if ids like "ORPHA:123" or "123" -> normalize to "123"
        ids <- sub("^ORPHA:", "", ids)
      }
      
      # align lengths
      if (length(ids) < length(nm)) ids <- c(ids, rep("", length(nm) - length(ids)))
      ids <- ids[seq_len(length(nm))]
      
      out <- mapply(function(p, id) {
        lab <- strip_orphanet_prefix(p)
        if (!nzchar(id)) return(lab)
        
        url <- paste0(
          "https://www.orpha.net/en/disease/detail/", id,
          "?name=", utils::URLencode(lab, reserved = TRUE),
          "&mode=name"
        )
        sprintf("<a href='%s' target='_blank'>%s</a>", url, lab)
      }, nm, ids, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      
      paste(out, collapse = "<br>")
    }
    
    out <- df %>%
      dplyr::transmute(
        SNP    = as.character(SNP),
        genename    = as.character(genename),
        cluster_id  = as.character(cluster_id),
        cluster_chr_n = as.character(cluster_chr_n),
        cluster_start = suppressWarnings(as.integer(cluster_start)),
        cluster_end   = suppressWarnings(as.integer(cluster_end)),
        Function_description = dplyr::na_if(trimws(as.character(Function_description)), ""),
        Disease_description  = dplyr::na_if(trimws(as.character(Disease_description)), ""),
        Orphanet_disorder  = dplyr::na_if(trimws(as.character(Orphanet_disorder)), ""),
        Orphanet_disorder_id = if ("Orphanet_disorder_id" %in% names(df)) as.character(Orphanet_disorder_id) else NA_character_,
        MIM_disease  = dplyr::na_if(trimws(as.character(MIM_disease)), ""),
        clinvar_trait  = dplyr::na_if(trimws(as.character(clinvar_trait)), ""),
        HPO_name = dplyr::na_if(trimws(as.character(HPO_name)), ""),
        clinvar_var_id   = if ("clinvar_id" %in% names(df)) dplyr::na_if(trimws(as.character(clinvar_id)), "") else NA_character_,
        clinvar_medgen_id = dplyr::case_when(
          "clinvar_MedGen_id" %in% names(df) ~ dplyr::na_if(trimws(as.character(clinvar_MedGen_id)), ""),
          # fallback if your clinvar_id is actually MedGen like C0950123
          "clinvar_id" %in% names(df) ~ dplyr::na_if(trimws(as.character(clinvar_id)), ""),
          TRUE ~ NA_character_
        ),
        OMIM_raw = .data[[omim_col]]
      ) %>%
      dplyr::filter(!is.na(cluster_id) & nzchar(cluster_id)) %>%
      
      # ---- aggregate SNPs into a single row per unique (cluster + gene + descriptions) ----
    dplyr::group_by(
      cluster_id, cluster_chr_n, cluster_start, cluster_end,
      genename,
      Function_description, Disease_description
    )%>%
      dplyr::summarise(
        
        # SNPs + links ClinVar (com ho tens) però ara agregant sobre totes les files del grup
        NonSyn_SNP = {
          dd <- dplyr::distinct(
            dplyr::tibble(SNP = SNP, clinvar_var_id = clinvar_var_id),
            SNP, clinvar_var_id
          )
          dd$SNP <- trimws(as.character(dd$SNP))
          dd <- dd[!is.na(dd$SNP) & nzchar(dd$SNP), , drop = FALSE]
          
          dd2 <- dd %>%
            dplyr::group_by(SNP) %>%
            dplyr::summarise(
              clinvar_var_id = {
                v <- trimws(as.character(clinvar_var_id))
                v[is.na(v)] <- ""
                v_num <- v[nzchar(gsub("[^0-9]", "", v))]
                if (length(v_num)) v_num[1] else (if (length(v) && nzchar(v[1])) v[1] else "")
              },
              .groups = "drop"
            )
          
          labs <- vapply(seq_len(nrow(dd2)), function(i) {
            rs  <- dd2$SNP[i]
            vid <- dd2$clinvar_var_id[i]
            url <- make_clinvar_variation_url(vid, rs)
            if (!is.na(url) && nzchar(url)) sprintf("<a href='%s' target='_blank'>%s</a>", url, rs) else rs
          }, character(1))
          
          paste(labs, collapse = "<br>")
        },
        
        n_NonSyn = {
          s <- unique(na.omit(trimws(as.character(SNP))))
          s <- s[nzchar(s)]
          length(s)
        },
        
        # Orphanet disorders (merge + clean)
        Orphanet_disorder = {
          v <- unlist(strsplit(paste0(na.omit(Orphanet_disorder), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          # optional: de-dup
          v <- unique(v)
          paste(v, collapse = ";")
        },
        
        Orphanet_disorder_id = {
          v <- unlist(strsplit(paste0(na.omit(Orphanet_disorder_id), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          v <- unique(sub("^ORPHA:", "", v))
          paste(v, collapse = ";")
        },
        
        # MIM disease (merge + clean)
        MIM_disease = {
          v <- unlist(strsplit(paste0(na.omit(MIM_disease), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          v <- unique(v)
          paste(v, collapse = ";")
        },
        
        # ClinVar trait + MedGen ids (merge + clean)
        clinvar_trait = {
          v <- unlist(strsplit(paste0(na.omit(clinvar_trait), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          v <- unique(v)
          paste(v, collapse = ";")
        },
        
        HPO_name = {
          v <- unlist(strsplit(paste0(na.omit(HPO_name), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          v <- trimws(v)
          v <- unique(v)
          paste(v, collapse = ";")
        },
        
        clinvar_medgen_id_agg = {
          v <- unlist(strsplit(paste0(na.omit(clinvar_medgen_id), collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
          v <- drop_bad_tokens(v)
          v <- unique(v)
          paste(v, collapse = ";")
        },
        
        .groups = "drop"
      ) %>%
      
      # ---- apply links to existing columns (no UCSC/OMIM columns) ----
    dplyr::mutate(
      # cluster_id → UCSC link (com ja ho tenies)
      cluster_id = {
        u <- mapply(make_ucsc_url, cluster_chr_n, cluster_start, cluster_end,
                    MoreArgs = list(db = "hg38"))
        cid <- as.character(cluster_id)
        ifelse(is.na(u) | !nzchar(u),
               cid,
               sprintf("<a href='%s' target='_blank'>%s</a>", u, cid))
      },
      
      # genename → GeneCards (com ja ho tenies)
      genename = {
        g <- vapply(genename, make_genecards_url, character(1))
        gn <- as.character(genename)
        ifelse(is.na(g), gn, sprintf("<a href='%s' target='_blank'>%s</a>", g, gn))
      },
      
      # MIM_disease → multiple OMIM links
      MIM_disease = vapply(MIM_disease, link_mim_cell, character(1)),
      # clinvar_trait
      clinvar_trait = mapply(link_clinvar_trait_cell, clinvar_trait, clinvar_medgen_id_agg),
      # Orphanet_disorder → multiple Orphanet links (needs ids if available)
      Orphanet_disorder = mapply(link_orpha_cell, Orphanet_disorder, Orphanet_disorder_id)
      
    ) %>%
      dplyr::select(
        cluster_id, genename,
        n_NonSyn, NonSyn_SNP,
        Function_description, Disease_description,
        Orphanet_disorder, MIM_disease, clinvar_trait,
        HPO_name
      )
    
    out
  })
  
  output$func_disease_tbl <- DT::renderDT({
    dat <- func_disease_df()
    req(nrow(dat) > 0)
    
    # ---- helper per text llarg pla ----
    mk_collapsible <- function(x, nchar_preview = 80) {
      x <- ifelse(is.na(x) | !nzchar(trimws(x)), "", as.character(x))
      
      vapply(x, function(txt) {
        txt_esc <- htmltools::htmlEscape(txt)
        
        if (nchar(txt) <= nchar_preview) {
          return(txt_esc)
        }
        
        preview <- paste0(substr(txt, 1, nchar_preview), "...")
        preview_esc <- htmltools::htmlEscape(preview)
        
        paste0(
          "<details>",
          "<summary style='cursor:pointer; color:#1A4E8A;'>", preview_esc, "</summary>",
          "<div style='white-space:normal; min-width:300px; max-width:600px; padding-top:6px;'>",
          txt_esc,
          "</div>",
          "</details>"
        )
      }, character(1))
    }
    
    # ---- helper per columnes HTML amb <br> (ex: NonSyn_SNP) ----
    # ---- helper per columnes HTML amb <br> (ex: NonSyn_SNP) ----
    mk_collapsible_br <- function(x, max_items_preview = 5) {
      x <- ifelse(is.na(x) | !nzchar(trimws(x)), "", as.character(x))
      
      vapply(x, function(txt) {
        parts <- unlist(strsplit(txt, "<br>", fixed = TRUE), use.names = FALSE)
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]
        
        if (length(parts) <= max_items_preview) {
          return(txt)
        }
        
        preview <- paste(parts[seq_len(max_items_preview)], collapse = "<br>")
        
        paste0(
          "<details>",
          "<summary style='cursor:pointer; color:#1A4E8A;'>",
          preview, "<br><span style='font-size:11px;color:#666;'>(", 
          length(parts) - max_items_preview, " more SNPs)</span>",
          "</summary>",
          "<div style='white-space:normal; min-width:220px; max-width:500px; padding-top:6px;'>",
          txt,
          "</div>",
          "</details>"
        )
      }, character(1))
    }
    
    # ---- convertir columnes llargues ----
    long_cols <- c(
      "Function_description",
      "Disease_description",
      "HPO_name"
    )
    
    for (cc in intersect(long_cols, names(dat))) {
      dat[[cc]] <- mk_collapsible(dat[[cc]], nchar_preview = 70)
    }
    
    # ---- NonSyn_SNP desplegable si hi ha > 5 SNPs ----
    if ("NonSyn_SNP" %in% names(dat)) {
      dat$NonSyn_SNP <- mk_collapsible_br(dat$NonSyn_SNP, max_items_preview = 5)
    }
    
    DT::datatable(
      dat,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      selection = "none",
      options = list(
        dom = "Bfrtip",
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          list(
            extend = "copy",
            exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)
          ),
          list(
            extend = "csv",
            exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)
          ),
          list(
            extend = "excel",
            exportOptions = list(modifier = list(page = "all"), stripHtml = TRUE)
          )
        ),
        columnDefs = list(
          list(
            targets = which(names(dat) %in% c(
              "NonSyn_SNP",
              "Function_description",
              "Disease_description",
              "HPO_name"
            )) - 1,
            className = "dt-left"
          )
        )
      ),
      callback = DT::JS("
      table.on('click', 'details, summary, details *', function(e) {
        e.stopPropagation();
      });
    ")
    )
  }, server = FALSE)
  
  #################################################################################  
  ########################### LD MODULE NonSyn #################################### 
  #################################################################################  
  
  build_ld_clusters_from_app_nonsyn <- function() {
    cl <- cluster_dt_view2()
    validate(need(is.data.frame(cl) && nrow(cl) > 0,
                  "No clusters available. Build clusters first."))
    
    cl <- as.data.frame(cl)
    
    validate(need(all(c("cluster_id", "chr", "start", "end") %in% names(cl)),
                  "cluster_dt_view2() must contain: cluster_id, chr, start, end."))
    
    cl$chr   <- chr_map_plink19(norm_chr_generic(cl$chr))
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    cl$cluster_id <- as.character(cl$cluster_id)
    
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) &
               cl$end >= cl$start & nzchar(cl$cluster_id), , drop = FALSE]
    validate(need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
    
    out <- cl[, c("cluster_id", "chr", "start", "end")]
    names(out) <- c("cluster_id", "chr", "cluster_start", "cluster_end")
    out
  }
  
  build_ld_candidates_from_app_nonsyn <- function() {
    
    # --- 1) GWAS hits (si n’hi ha al NonSyn) ---
    out <- data.frame(
      chr = integer(), pos_ini = integer(), pos_end = integer(),
      id_hit = character(), classe = character(),
      stringsAsFactors = FALSE
    )
    
    h <- tryCatch(hits_df(), error = function(e) NULL)
    if (is.data.frame(h) && nrow(h) > 0) {
      validate(need(all(c("CHR","BP") %in% names(h)), "hits_df() must contain CHR and BP columns."))
      
      rsid <- if ("rsid" %in% names(h)) as.character(h$rsid) else if ("snp" %in% names(h)) as.character(h$snp) else NA_character_
      rsid[is.na(rsid) | !nzchar(rsid)] <- paste0("chr", h$CHR, ":", h$BP)
      
      out <- rbind(out, data.frame(
        chr     = chr_map_plink19(norm_chr_generic(h$CHR)),
        pos_ini = suppressWarnings(as.integer(h$BP)),
        pos_end = suppressWarnings(as.integer(h$BP)),
        id_hit  = rsid,
        classe  = "GWAS",
        stringsAsFactors = FALSE
      ))
    }
    
    # --- 2) NonSyn variants (dbNSFP) ---
    ns <- dbnsfp_norm_df()
    validate(need(is.data.frame(ns) && nrow(ns) > 0,
                  "No NonSyn (dbNSFP) data available. Run extraction first."))
    
    ns <- as.data.frame(ns)
    
    # assegura cluster assignment (igual que al teu export)
    ns_clcol <- pick_first_col(ns, c("cluster_id", "cluster_chr_n", "cluster", "label"))
    
    need_assign <- is.null(ns_clcol) ||
      all(is.na(ns[[ns_clcol]]) | !nzchar(as.character(ns[[ns_clcol]])))
    
    if (need_assign) {
      validate(need(exists("add_clusters_to_dbnsfp", mode = "function"),
                    "add_clusters_to_dbnsfp() not found (needed to assign clusters)."))
      ns <- add_clusters_to_dbnsfp(ns, cluster_dt_view2())
      ns_clcol <- "cluster_id"
    }
    
    # filtra només els que estan dins cluster si tens in_cluster
    if ("in_cluster" %in% names(ns)) {
      keep <- !is.na(ns$in_cluster) & as.logical(ns$in_cluster)
      ns <- ns[keep, , drop = FALSE]
    }
    
    validate(need(nrow(ns) > 0, "dbNSFP data has no rows assigned to clusters (nothing for LD)."))
    validate(need(!is.null(ns_clcol) && ns_clcol %in% names(ns), "Could not find cluster column after assignment."))
    
    ns_chr_col <- pick_first_col(ns, c("chr", "#chr", "CHR", "chrom", "chromosome"))
    ns_pos_col <- pick_first_col(ns, c("BP", "POS", "pos(1-based)", "pos", "position"))
    validate(need(!is.null(ns_chr_col), "dbNSFP data has no chromosome column."))
    validate(need(!is.null(ns_pos_col), "dbNSFP data has no position column."))
    
    ns_id_col <- pick_first_col(ns, c("rs_dbSNP", "rsid", "RSID", "variant_id", "variant", "ID", "id"))
    if (is.null(ns_id_col)) {
      ns_id_col <- ".id_tmp"
      ns[[ns_id_col]] <- paste0("chr", ns[[ns_chr_col]], ":", ns[[ns_pos_col]])
    }
    
    nonsyn_cand <- data.frame(
      chr     = chr_map_plink19(norm_chr_generic(ns[[ns_chr_col]])),
      pos_ini = suppressWarnings(as.integer(ns[[ns_pos_col]])),
      pos_end = suppressWarnings(as.integer(ns[[ns_pos_col]])),
      id_hit  = as.character(ns[[ns_id_col]]),
      classe  = "NonSyn_variant",
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out, nonsyn_cand)
    
    out <- out[is.finite(out$chr) & is.finite(out$pos_ini) & is.finite(out$pos_end) &
                 out$pos_end >= out$pos_ini & nzchar(out$id_hit), , drop = FALSE]
    
    validate(need(nrow(out) > 0, "No candidates available for LD."))
    out
  }
  
  #  ##############################################################################
  #  ############################ LD Module   #####################################
  #  ################# GWAS -> LD module canonical cluster input###################
  
  ld_module_server(
    "ld",
    app_tag = "nonsyn",
    activate_r = reactive(TRUE)
  )
  
  #################################################################################  
  ##############################  RESET session   #################################  
  #################################################################################  
  
  observeEvent(input$reset_case, {
    
    safe <- function(expr) try(expr, silent = TRUE)
    
    removeModal()
    
    # (B) Neteja seleccions DT
    for (id in c("hits_tbl", "cluster_dt", "nonsyn_table", "go_table", "kegg_table",
                 "go_terms_table", "radar_selection_table",
                 "ld_cluster_fl_table", "func_disease_tbl")) {
      safe(DT::selectRows(DT::dataTableProxy(id), NULL))
    }
    
    # (C) Esborra selected_intervals.range si existeix
    if (exists("workdir", inherits = TRUE) && is.character(workdir) && dir.exists(workdir)) {
      f_rng <- file.path(workdir, "selected_intervals.range")
      if (file.exists(f_rng)) safe(unlink(f_rng, force = TRUE))
    }
    
    # ------------------------------------------------------------
    # (D) HUB payload: IMPORTANT
    # - NO esborrem GWAS (perquè no volem perdre el plot superior)
    # - SÍ esborrem CLUSTERS rebuts del Hub (això era el problema)
    # ------------------------------------------------------------
    if (exists("gi_sync", inherits = TRUE) && is.list(gi_sync)) {
      if (!is.null(gi_sync$clusters_shared) && is.function(gi_sync$clusters_shared)) {
        safe(gi_sync$clusters_shared(NULL))
      }
      # NO tocar:
      # if (!is.null(gi_sync$gwas_shared) && is.function(gi_sync$gwas_shared)) gi_sync$gwas_shared(NULL)
    }
    
    # (E) Reset del core clustering LOCAL (sense tocar GWAS)
    # Preferència: si tens motor canònic (gi_cl)
    if (exists("gi_cl", inherits = TRUE) && is.list(gi_cl)) {
      safe(gi_cl$clusters_cur(tibble::tibble()))
      safe(gi_cl$intervals_raw(tibble::tibble(chr=integer(), start=integer(), end=integer(), label=character())))
    } else {
      # fallback si uses reactiveVal directes
      if (exists("clusters_cur", inherits = TRUE) && is.function(clusters_cur))  safe(clusters_cur(NULL))
      if (exists("intervals_raw", inherits = TRUE) && is.function(intervals_raw)) safe(intervals_raw(NULL))
    }
    
    # IMPORTANT: NO cridar reset_intervals_clusters() si et mata el GWAS plot
    # if (exists("reset_intervals_clusters", mode = "function")) reset_intervals_clusters()
    
    # (F) Reset UCSC window + region
    if (exists("ucsc_region", inherits = TRUE) && is.function(ucsc_region)) safe(ucsc_region(NULL))
    if (exists("window_selected", inherits = TRUE) && is.function(window_selected)) safe(window_selected(NULL))
    
    # (G) Reset dbNSFP pipeline state
    if (exists("ann_path", inherits = TRUE) && is.function(ann_path)) safe(ann_path(NULL))
    if (exists("ann_ready", inherits = TRUE) && is.function(ann_ready)) safe(ann_ready(FALSE))
    if (exists("dbnsfp_log_text", inherits = TRUE) && is.function(dbnsfp_log_text)) safe(dbnsfp_log_text(""))
    if (exists("dbnsfp_norm_path_csv", inherits = TRUE) && is.function(dbnsfp_norm_path_csv)) safe(dbnsfp_norm_path_csv(NULL))
    if (exists("dbnsfp_norm_path_rds", inherits = TRUE) && is.function(dbnsfp_norm_path_rds)) safe(dbnsfp_norm_path_rds(NULL))
    if (exists("dbnsfp_final_path_csv", inherits = TRUE) && is.function(dbnsfp_final_path_csv)) safe(dbnsfp_final_path_csv(NULL))
    if (exists("dbnsfp_final_path_rds", inherits = TRUE) && is.function(dbnsfp_final_path_rds)) safe(dbnsfp_final_path_rds(NULL))
    
    # (H) Reset enrichment triggers
    if (exists("go_trigger", inherits = TRUE) && is.function(go_trigger)) safe(go_trigger(0L))
    if (exists("kegg_trigger", inherits = TRUE) && is.function(kegg_trigger)) safe(kegg_trigger(0L))
    
    # (I) Reset tracks UCSC
    # Manté GWAS track si el teu plot superior en depèn
    # session$userData$track_gwas_data <- tibble::tibble()   # <-- NO (per NonSyn)
    session$userData$track_nonsyn_data <- tibble::tibble()
    
    showNotification("Reset done. Ready for a new case.", type = "message", duration = 2)
    
  }, ignoreInit = TRUE)
  
  observeEvent(kegg_trigger(), {
    g <- tryCatch(entrez_scope(), error = function(e) character(0))
    g <- unique(na.omit(as.character(g))); g <- g[nzchar(g)]
    message("[KEGG] genes n=", length(g), " head=", paste(head(g, 6), collapse = ","))
  }, ignoreInit = TRUE)
  

  ##≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  
  
}  # end server

shinyApp(ui, server)

