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
merge_significant_windows <- function(win_df, gap_bp = 0L) {
  library(data.table)
  
  # gap segur
  gap_bp <- as.integer(gap_bp)
  if (is.na(gap_bp)) gap_bp <- 0L
  
  dt <- as.data.table(win_df)
  
  # assegura columnes
  stopifnot(all(c("chr","start","end") %in% names(dt)))
  
  # tipus + neteja NA
  dt[, chr   := as.integer(chr)]
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt <- dt[!is.na(chr) & !is.na(start) & !is.na(end)]
  
  # si no hi ha res, retorna buit
  if (nrow(dt) == 0) return(as.data.frame(dt))
  
  setorder(dt, chr, start, end)
  
  merged <- dt[, {
    # casos trivials
    if (.N == 1L) {
      data.table(start = start[1], end = end[1])
    } else {
      out_s <- integer()
      out_e <- integer()
      
      cur_s <- start[1]
      cur_e <- end[1]
      
      # només iterar si hi ha 2 o més files
      for (i in 2L:.N) {
        s <- start[i]
        e <- end[i]
        
        # si per algun motiu hi ha NA, salta
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
      # tanca últim
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

  // Try now + retries (Shiny may not be ready at DOMContentLoaded behind proxy)
  var tries = 0;
  var iv = setInterval(function(){
    tries++;
    if (pushQS() || tries > 50) clearInterval(iv); // ~5s max
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

    /* Plotly widgets (plotlyOutput) */
    .plotly.html-widget {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    /* Base R plots / ggplot outputs (plotOutput) */
    .shiny-plot-output {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    /* DT tables (DTOutput) container */
    .dataTables_wrapper {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }

    /* DT inner bits (optional: keeps header cohesive) */
    table.dataTable,
    table.dataTable thead th,
    table.dataTable tbody td {
      background-color: transparent !important;
    }

    /* If you use spinner wrappers, keep the same background */
    .shiny-spinner-output-container {
      background: #f2f2f2 !important;
      border-radius: 10px;
      padding: 10px;
    }
  "))
  ),
  
  # ========================================================================
  # Main tab: Analysis
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    
    sidebarLayout(
      # --------------------------------------------------------------------
      # LEFT SIDEBAR
      # --------------------------------------------------------------------
      sidebarPanel(
        width = 3,
        
        # -----------------------------
        # Analysis mode (new vs old)
        # -----------------------------
        h4(
          "Select analysis mode",
          style = "color:#1A4E8A; font-size:22px; font-weight:700;"
        ),
        
        radioButtons(
          "use_preloaded_dbnsfp",
          "",
          c(
            "🧪 Perform a new analysis (run steps 1–3)" = "new",
            "🧾 Analyze a precomputed dbNSFP output file (skip steps 2–3)" = "old"
          ),
          selected = "new"
        ),
        
        tags$hr(),
        
        # -----------------------------
        # Step 1: Load GWAS hits
        # -----------------------------
        h3(
          "Step 1 · Load GWAS hits",
          style = "color:#1A4E8A; font-size:22px; font-weight:700;"
        ),
        
        fluidRow(
          column(6, actionButton("info_00", "ℹ️ file format")),
          column(6, actionButton(
            "reset_case",
            "Reset",
            icon = icon("rotate-left"),
            class = "btn-warning"
          ))),
        
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
        # Old mode: load precomputed dbNSFP output
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_dbnsfp == 'old'",
          fileInput(
            "dbnsfp_preload",
            "Load dbNSFP output file (CSV or RDS)",
            accept = c(".csv", ".rds")
          )
        ),
        
        # -----------------------------
        # New mode: Steps 2–3
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_dbnsfp == 'new'",
          tagList(
            
            # ---- Step 2
            h3(
              "Step 2 · Clustering GWAS hits",
              style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
            ),
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
              sliderInput("pthr", "-log10(P) threshold", min = 2, max = 20, value = 5, step = 0.5),
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
                min = 2, max = 20, value = 5, step = 0.1
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
            
            tags$hr(),
            
            h3("Step 3 · Extract nonsyn variants in all clusters", style="color:#1A4E8A; font-weight:700;"),
            
            # Path al .out gran (dbNSFP ja anotat)
            textInput(
              "dbnsfp_out_all",
              "Path to dbNSFP biallelic_nonsyn.out file",
              value = ns_def$dbnsfp_all
              ),
            actionButton(
              "run_nonsyn_clusters",
              "Extract nonsyn variants in all clusters",
              #class = "btn-primary",
              style = "background-color:#ffdd57; font-weight:bold;"
            ),
            # tags$hr(),
            # verbatimTextOutput("ranges_preview"), # si ja el tens, deixa'l; si no, pots afegir-lo aquí
            h5("Extract log file"),
            verbatimTextOutput("dbnsfp_log"),
            tags$hr(),
            # uiOutput("dbnsfp_empty_msg"),
            h3("Step 4 · Download files", style="color:#1A4E8A; font-weight:700;"),
            downloadButton("dl_dbnsfp_csv", "Download normalized CSV"),
            downloadButton("dl_dbnsfp_rds", "Download normalized RDS"),
            tags$hr(),
            downloadButton("dl_candidates_zip", "⬇️ Download NonSyn candidates (ZIP)")
          )
        )
      ),
      
      # --------------------------------------------------------------------
      # MAIN PANEL (plots + tables)
      # --------------------------------------------------------------------
      mainPanel(
        # -----------------------------
        # Global selector (only for Manhattan + Summary stats)
        # -----------------------------
        conditionalPanel(
          condition = "input.main_tabs == 'Manhattan' || input.main_tabs == 'Summary stats'|| input.main_tabs == 'LD'",
          wellPanel(
            style = "margin-bottom:6px; padding:6px 10px;",  # antes 12px / 10px
            nonsyn_metric_selector_ui()
          )
        ),
        
        # -----------------------------
        # Main tabset
        # -----------------------------
        tabsetPanel(
          id = "main_tabs",
          
          # ============================================================
          # Manhattan
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            value = "Manhattan",
            
            fluidRow(
              column(
                12,
                withSpinner(plotlyOutput("manhattan_combo", height = "600px"))
              )
            ),
            
            tags$hr(),
            fluidRow(
              column(
                6,
                uiOutput("debug_ucsc_state")
              ),
              column(
                2,
                h4("UCSC links to expanded region:")),
              column(
                4, 
                div(
                  style = "margin-top:8px; background-color:#f2f2f2; padding:10px; border-radius:8px; width:100%;",
                  uiOutput("ucsc_link_gwas"),
                  uiOutput("ucsc_link_nonsyn")
                ),
              ),
            ),
            
            fluidRow(
              column(
                12,
                conditionalPanel(
                  condition = "input.use_preloaded_dbnsfp == 'new'",
                  tags$hr(),
                  h4("Clusters summary"),
                  DTOutput("cluster_dt")
                )
              )
            )
          )
          ,
          
          
          # ============================================================
          # NonSyn – dbNSFP
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📄 NonSyn – dbNSFP</span>"),
            value = "NonSyn",
            
            uiOutput("hm_chr_selector2"),
            
            radioButtons(
              "clinical_class",
              "Functional class",
              c(
                "Default",
                "Pathogenicity",
                "Conservation",
                "Functional impact",
                "LoF/Haploinsufficiency",
                "Gene damage"
              ),
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
          
          # ============================================================
          # GWAS hits
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔬 GWAS hits</span>"),
            value = "GWAS",
            uiOutput("gwas_tab_body")
          ),
          
          # ============================================================
          # Summary stats
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Summary stats</span>"),
            value = "Summary stats",
            
            fluidRow(
              column(
                12,
                withSpinner(plotlyOutput("nonsyn_box_combo_plotly", height = "350px")),
                fluidRow(
                  column(4, downloadButton("dl_nonsyn_chr_stats_csv",     "⬇️ Summary stats by chr CSV")),
                  column(4, downloadButton("dl_nonsyn_gene_stats_csv",    "⬇️ Summary stats by gene CSV")),
                  column(4, downloadButton("dl_nonsyn_cluster_stats_csv", "⬇️ Summary stats by cluster CSV"))
                )
              )
            ),
            
            tags$hr(),
            
            fluidRow(
              column(
                6,
                h4("Top genes barplot"),
                withSpinner(plotOutput("nonsyn_gene_bar", height = "300px"))
              ),
              column(
                6,
                h4("Metric distribution"),
                withSpinner(plotOutput("nonsyn_metric_bar", height = "300px"))
              )
            )
          ),
          
          
          # ============================================================
          # Radar plots
          # ============================================================
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
                    "Cluster" = "cluster",
                    "Chromosome" = "chr",
                    "Genes" = "gene"
                  ),
                  selected = "cluster"
                ),
                
                conditionalPanel(
                  condition = "input.radar_mode == 'cluster'",
                  selectInput("radar_chr_for_cluster", "Chromosome", choices = NULL),
                  selectInput("radar_cluster_id", "Cluster", choices = NULL)
                ),
                
                conditionalPanel(
                  condition = "input.radar_mode == 'chr'",
                  selectInput("radar_chr", "Chromosome", choices = NULL)
                ),
                
                conditionalPanel(
                  condition = "input.radar_mode == 'gene'",
                  selectInput("radar_gene", "Gene (by cluster)", choices = NULL)
                )
              ),
              
              mainPanel(
                width = 10,
                
                actionButton(
                  "info_11",
                  label = "Metric info",
                  icon  = icon("info-circle"),
                  class = "btn btn-outline-secondary"
                ),
                
                withSpinner(plotOutput("nonsyn_radar", height = "650px")),
                tags$hr(),
                DT::DTOutput("radar_selection_table")
              )
            )
          ),
          
          # ============================================================
          # GoSlim plots
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📄 GoSlim plots</span>"),
            value = "tab_go",
            
            fluidRow(
              column(
                width = 2,
                radioButtons(
                  "go_mode",
                  "Scope",
                  choices  = c("Total" = "total", "Chromosome" = "chr", "Cluster" = "cluster"),
                  selected = "total"
                ),
                
                conditionalPanel(
                  condition = "input.go_mode == 'cluster'",
                  selectInput("go_chr_for_cluster", "Chromosome", choices = NULL),
                  selectInput("go_cluster_id", "Cluster", choices = NULL)
                ),
                
                conditionalPanel(
                  condition = "input.go_mode == 'chr'",
                  selectInput("go_chr", "Chromosome", choices = NULL)
                )
              ),
              
              column(
                width = 10,
                withSpinner(plotlyOutput("go_class_plot", height = "650px")),
                tags$hr()
              )
            ),
            
            fluidRow(
              column(
                width = 12,
                div(
                  style = "width:100%; overflow-x:auto;",
                  DT::DTOutput("go_terms_table", width = "100%")
                )
              )
            )
          ),
          
          # ============================================================
          # Enrichment (GO/KEGG) - fixed duplicate input IDs
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧬 Enrichment</span>"),
            value = "tab_enrich",
            
            sidebarLayout(
              sidebarPanel(
                width = 3,
                
                radioButtons(
                  "func_scope",
                  "Scope",
                  choices  = c("Global" = "global", "Cluster" = "cluster"),
                  selected = "global"
                ),
                
                conditionalPanel(
                  condition = "input.func_scope == 'cluster'",
                  uiOutput("func_cluster_ui")
                ),
                
                tags$hr(),
                
                selectInput(
                  "enrich_background",
                  "Background",
                  choices  = c("Reference annotated genes" = "orgdb", "Dataset genes" = "dataset"),
                  selected = "orgdb"
                ),
                
                # -------------------------------
                # Shared p-value cutoff (single ID, used for GO + KEGG)
                # -------------------------------
                numericInput(
                  "enrich_pcut",
                  "FDR cutoff",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01
                ),
                
                # -------------------------------
                # GO-only controls (shown only when GO tab is active)
                # -------------------------------
                conditionalPanel(
                  condition = "input.enrich_tabs == 'tab_enrich_go'",
                  
                  checkboxGroupInput(
                    "go_ontos",
                    "GO ontologies",
                    choices  = c("BP", "CC", "MF"),
                    selected = c("BP", "CC", "MF"),
                    inline   = TRUE
                  ),
                  
                  numericInput("go_topn", "Top GO terms per ontology", value = 10, min = 1, max = 50, step = 1),
                  checkboxInput("go_simplify", "Simplify GO terms (optional)", value = FALSE)
                ),
                
                # -------------------------------
                # KEGG-only controls (shown only when KEGG tab is active)
                # -------------------------------
                conditionalPanel(
                  condition = "input.enrich_tabs == 'tab_enrich_kegg'",
                  numericInput("enrich_kegg_top", "Top KEGG pathways", value = 15, min = 1, max = 50, step = 1)
                ),
                
                tags$hr(),
                
                # Common controls (GO + KEGG)
                actionButton("info_12", "ℹ GSSize"),
                numericInput("enrich_min_gs", "minGSSize", value = 10, min = 1, step = 1),
                numericInput("enrich_max_gs", "maxGSSize", value = 500, min = 10, step = 10),
                
                tags$hr(),
                
                actionButton(
                  "run_enrich",
                  label = "Run enrichment",
                  icon  = icon("play"),
                  class = "btn btn-primary"
                )
              ),
              
              mainPanel(
                width = 9,
                uiOutput("enrich_bg_note"),
                
                tabsetPanel(
                  id = "enrich_tabs",
                  
                  tabPanel(
                    title = HTML("<span style='font-size:15px; font-weight:600;'>GO</span>"),
                    value = "tab_enrich_go",
                    withSpinner(plotOutput("go_bar", height = 400)),
                    tags$hr(),
                    withSpinner(DT::DTOutput("go_table"))
                  ),
                  
                  tabPanel(
                    title = HTML("<span style='font-size:15px; font-weight:600;'>KEGG</span>"),
                    value = "tab_enrich_kegg",
                    withSpinner(plotOutput("kegg_bar", height = 400)),
                    tags$hr(),
                    withSpinner(DT::DTOutput("kegg_table"))
                  )
                )
              )
            )
          ),
          # ============================================================
          # Function & Disease
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧾 Function & Disease</span>"),
            value = "FuncDisease",
            DT::DTOutput("func_disease_tbl")
          )
        )
      )
    )
  ),
  #  # ============================================================
  #  # LD
  #  # ============================================================
  
  # ==========================
  # TAB PRINCIPAL 2 (LD)
  # ==========================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 LD</span>"),
    ld_module_ui("ld")
  )
)



# ===========================
# SERVER
# ===========================
server <- function(input, output, session) {
  
  source(file.path(SHARED, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
  
  # ===========================
  # 1) SLAVE sync (Hub -> app)
  # ===========================
  gi_sync <- gi_slave_canonical_init(session)
  
  # Guardem GWAS rebut en un RV canònic (si no el tens ja)
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
  
  # ===========================
  # 2) CANONICAL cluster engine (local build_ranges)
  # ===========================
  gwas_df <- reactive({
    # només té sentit en NEW mode
    if (!identical(input$use_preloaded_dbnsfp, "new")) return(tibble::tibble())
    gwas_shared_rv()
  })
  
  cl_engine <- gi_clusters_canonical_init(
    session = session, input = input, output = output,
    gwas_df = gwas_df,
    build_btn_id    = "build_ranges",
    clusters_dt_id  = "cluster_dt",
    hits_rows_id    = "hits_tbl_rows_selected",
    app_count_col   = "n_nonsyn"   # columna comptador de l’app
  )
  
  # A partir d’aquí, aquests són ELS ACCESSORS CANÒNICS (per tota l’app)
  intervals_raw     <- cl_engine$intervals_raw
  clusters_cur      <- cl_engine$clusters_cur
  selected_cluster  <- cl_engine$selected_cluster
  
  # ===========================
  # 3) Hub clusters override (si arriben del master, manen)
  # ===========================
  observeEvent(gi_sync$clusters_shared(), {
    cl <- gi_sync$clusters_shared()
    req(is.data.frame(cl), nrow(cl) > 0)
    
    # normalitza ids (cluster_id/cluster_chr_n consistents)
    cl <- standardize_cluster_ids(cl)
    
    # normalitza camps mínims
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
        dplyr::transmute(chr = .data$chr, start = .data$start, end = .data$end, label = .data$cluster_id)
    )
    
    cat("[SLAVE] applied CLUSTERS rows=", nrow(cl), "\n")
  }, ignoreInit = FALSE)
  
  # ===========================
  # 4) clusters_src (una sola font per overlays/plots)
  # ===========================
  clusters_src <- reactive({
    cl <- tryCatch(clusters_cur(), error = function(e) NULL)
    if (is.data.frame(cl) && nrow(cl) > 0) return(cl)
    NULL
  })
  
  
  ################## end slave step 
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
    req(nrow(df) > 0)
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
    
    # threshold segons mètode
    thr <- if (identical(input$cluster_method, "window")) input$pthr else input$min_logp
    thr <- suppressWarnings(as.numeric(thr))
    if (!is.finite(thr)) thr <- -Inf
    
    # garantir rsid (si no existeix, usa snp)
    if (!"rsid" %in% names(df)) df$rsid <- NA_character_
    if (!"snp"  %in% names(df)) df$snp  <- NA_character_
    df$rsid <- as.character(df$rsid)
    df$snp  <- as.character(df$snp)
    df$rsid <- dplyr::coalesce(df$rsid, df$snp)
    
    df %>%
      dplyr::filter(.data$logp >= thr) %>%
      dplyr::arrange(dplyr::desc(.data$logp)) %>%
      dplyr::select(
        dplyr::all_of(c("CHR", "BP", "snp")),   # aquests han d'existir
        dplyr::any_of("rsid"),                  # tolerant
        p = .data$Pval,
        .data$logp
      )
  })
  
  
  output$hits_tbl <- renderDT({
    req(input$use_preloaded_dbnsfp %||% "new")
    req(identical(input$use_preloaded_dbnsfp, "new"))
    
    h <- hits_df()
    
    if (is.null(h) || !nrow(h)) {
      return(DT::datatable(data.frame(Message = "No hits above threshold."), 
                           options = list(dom = "t")))
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
      options    = list(
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE
      )
    )
  })
  
  output$gwas_tab_body <- renderUI({
    req(input$use_preloaded_dbnsfp)
    
    if (!identical(input$use_preloaded_dbnsfp, "new")) {
      return(
        div(
          style = "padding:20px; color:#777;",
          em("GWAS hits are not available when using precomputed dbNSFP results.")
        )
      )
    }
    
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
    
    ns <- tryCatch(nonsyn_manhattan_df(), error = function(e) NULL)
    ns_ok <- is.data.frame(ns) && nrow(ns) > 0 && all(c("CHR", "POS") %in% names(ns))
    
    if (ns_ok) {
      ns_chr <- suppressWarnings(as.integer(gsub("^chr", "", as.character(ns$CHR), ignore.case = TRUE)))
      ns_pos <- suppressWarnings(as.integer(ns$POS))
      
      out$n_nonsyn <- mapply(
        function(cc, s, e) sum(ns_chr == cc & ns_pos >= s & ns_pos <= e, na.rm = TRUE),
        out$chr_num, out$start, out$end
      )
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
          buttons    = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 10,
          options = list(dom = "t"),
          scrollX    = TRUE
        )
  
      ))
    }
    
    dt2 <- dt %>%
      dplyr::select(cluster_id, chr, start, end, n_snps, n_nonsyn, top_snp, top_logp, size_kb)
    
    DT::datatable(dt2, rownames = FALSE, 
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons    = c("copy", "csv", "excel", "pdf", "print"),
                    pageLength = 10,
                    scrollX    = TRUE
                  )
    )  %>%
      DT::formatRound(c("top_logp", "size_kb"), digits = 2)
  })
  
 
  # ------------------------------------------------------------
  # Canonical invalidation: if inputs change, old intervals/clusters are no longer valid
  # ------------------------------------------------------------
  
  observeEvent(
    list(input$cluster_method, input$pthr, input$flank, input$min_logp, input$min_snps),
    {
      if (!identical(input$use_preloaded_dbnsfp, "new")) return()
      cl_now <- tryCatch(clusters_cur(), error = function(e) NULL)
      if (!is.data.frame(cl_now) || nrow(cl_now) == 0) return()
      
      reset_intervals_clusters()
      vcf_out_path(NULL)
      dbsnp_hits_val(NULL)
      dbnsfp_final_path_csv(NULL)
      dbnsfp_final_path_rds(NULL)
      dbnsfp_norm_path_csv(NULL)
      dbnsfp_norm_path_rds(NULL)
    },
    ignoreInit = TRUE
  )
  
  # Optional: if switching clustering method, also reset immediately
  observeEvent(input$cluster_method, {
    if (!identical(input$use_preloaded_dbnsfp, "new")) return()
    reset_intervals_clusters()
  }, ignoreInit = TRUE)
  
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
  gwas_df_local <- gwas_df
  gwas_df <- reactive({
    df <- gwas_shared_rv()
    if (is.data.frame(df) && nrow(df) > 0) return(df)
    gwas_df_local()
  })
  
  
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
  
  # Helper: normalize chr to integer-like (1..23 etc.)
  norm_chr_int_dbnsfp <- function(x) {
    x <- toupper(as.character(x))
    x <- sub("^CHR", "", x)
    x <- sub("^CHROM", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHROMOSOME", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    
    x <- sub("^CHR", "", x) # defensive
    x <- sub("^CHROM", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- sub("^CHR", "", x)
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    x <- gsub("^CHR", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x <- gsub("^CHR", "", x)
    x <- gsub("^CHROM", "", x)
    
    x[x == "X"]  <- "23"
    x[x == "Y"]  <- "24"
    x[x %in% c("MT","M")] <- "26"
    
    suppressWarnings(as.integer(x))
  }
  
  # Helper: save final dataset paths (if you deleted the old helper)
  save_dbnsfp_final_local <- function(df, out_dir, stem = "dbnsfp_normalized") {
    stopifnot(is.data.frame(df))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_csv <- file.path(out_dir, paste0(stem, ".csv"))
    out_rds <- file.path(out_dir, paste0(stem, ".rds"))
    readr::write_csv(df, out_csv)
    saveRDS(df, out_rds)
    list(csv = out_csv, rds = out_rds)
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
      
      # -------------------------
      # Preconditions / blockers
      # -------------------------
      if (!identical(input$use_preloaded_dbnsfp, "new")) {
        showNotification("STOP · Mode must be 'new'", type = "error", duration = 8)
        safe_log("[STOP] use_preloaded_dbnsfp=", as.character(input$use_preloaded_dbnsfp), "\n")
        return()
      }
      
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
      
    # cl0 <- as.data.frame(cl0)
    # nm <- names(cl0)
      
      cl0 <- as.data.frame(cl0)
      
      # --- CANONICAL: mai assumeixis cluster_chr_n; crea cluster_id sempre ---
      # Evita warnings si algun origen NO té cluster_chr_n
      if (!"cluster_id" %in% names(cl0)) cl0$cluster_id <- NA_character_
      
      # Si hi ha cluster_chr_n, l'aprofitem PER OMPLIR cluster_id (sense referenciar-la si no hi és)
      if ("cluster_chr_n" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          as.character(cl0$cluster_chr_n))
      }
      
      # Fallbacks típics si encara està buit
      if ("label" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          as.character(cl0$label))
      }
      if ("cluster" %in% names(cl0)) {
        cl0$cluster_id <- dplyr::coalesce(as.character(cl0$cluster_id),
                                          paste0("cluster_", as.character(cl0$cluster)))
      }
      # Últim recurs: id determinista
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
      
      validate(need(exists("add_clusters_to_dbnsfp", mode = "function"),
                    "add_clusters_to_dbnsfp() not available."))
      
      df_final <- add_clusters_to_dbnsfp(df_norm, cl_assign)
      df_final <- df_final %>% dplyr::filter(.data$in_cluster)
      validate(need(nrow(df_final) > 0, "No normalized rows remained after cluster assignment."))
      
      # -------------------------
      # Save final outputs
      # -------------------------
      # -------------------------
      # Save final outputs (ROBUST)
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
      
      # --- setters NA/len-safe (EVITA "replacement has length zero") ---
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
      
   #   if (exists("dbnsfp_final_path_rds", mode = "function")) dbnsfp_final_path_rds(final_paths$rds)
   #   if (exists("dbnsfp_final_path_csv", mode = "function")) dbnsfp_final_path_csv(final_paths$csv)
      
   #  if (exists("dbnsfp_norm_path_rds", mode = "function")) dbnsfp_norm_path_rds(res$rds)
   #  if (exists("dbnsfp_norm_path_csv", mode = "function")) dbnsfp_norm_path_csv(res$csv)
      
      showNotification(
        paste0("OK · tabix lines=", total_lines, " | final rows=", nrow(df_final)),
        type = "message",
        duration = 8
      )
      
      safe_log("[OK] Finished. tabix lines=", total_lines, " final=", nrow(df_final), "\n")
      
    }, error = function(e) {
      msg <- conditionMessage(e)
      showNotification(paste0("ERROR · ", msg), type = "error", duration = 12)
      safe_log("[ERROR] ", msg, "\n")
    })
    
  }, ignoreInit = TRUE)
  
  
  ############################################################################  
  # -------------------------
  # Single reactive: dbNSFP dataset used everywhere
  # - OLD mode: reads uploaded CSV/RDS
  # - NEW mode: reads FINAL paths (prefer RDS)
  # -------------------------
  dbnsfp_norm_df <- reactive({
    
    # -------------------------
    # PRECOMPUTED MODE (old)
    # -------------------------
    if (identical(input$use_preloaded_dbnsfp, "old")) {
      
      req(input$dbnsfp_preload)
      path <- input$dbnsfp_preload$datapath
      ext  <- tolower(tools::file_ext(input$dbnsfp_preload$name))
      
      df <- tryCatch({
        if (ext == "rds") readRDS(path) else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      
      validate(
        need(is.data.frame(df), "Invalid dbNSFP file."),
        need(all(c("chr", "BP") %in% names(df)), "dbNSFP file must contain columns: chr, BP.")
      )
      
      df <- sanitize_colnames(df)
      df <- harmonize_dbnsfp_loaded(df)
      
      # Optional: expose "final" paths also in old mode (so downstream logic is unified)
      if (ext == "rds") {
        dbnsfp_final_path_rds(path)
        dbnsfp_final_path_csv(NULL)
      } else {
        dbnsfp_final_path_csv(path)
        dbnsfp_final_path_rds(NULL)
      }
      
      return(df)
    }
    
    # -------------------------
    # NEW MODE: read FINAL (prefer RDS)
    # -------------------------
    rds <- dbnsfp_final_path_rds()
    csv <- dbnsfp_final_path_csv()
    
    if ((is.null(rds) || !file.exists(rds)) && (is.null(csv) || !file.exists(csv))) {
      # Not ready yet
      return(tibble::tibble())
    }
    
    df <- NULL
    if (!is.null(rds) && file.exists(rds)) {
      df <- readRDS(rds)
    } else {
      df <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    df <- sanitize_colnames(df)
    df <- harmonize_dbnsfp_loaded(df)
    
    # In NEW mode we require cluster_id
    if (!"cluster_id" %in% names(df)) return(tibble::tibble())
    
    df
  })
  
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
      plotly::layout(showlegend = FALSE)
    
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE
      )
    )
  })
  
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 20,
        scrollX    = TRUE
      )
    )
  })
  
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
  
  # A) Fill chromosomes
  observeEvent(radar_base_df(), {
    df <- radar_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "radar_chr", choices = character(0), selected = NULL)
      updateSelectInput(session, "radar_chr_for_cluster", choices = character(0), selected = NULL)
      return()
    }
    
    chr_choices <- sort(unique(df$CHR))
    chr_choices <- chr_choices[is.finite(chr_choices)]
    if (!length(chr_choices)) {
      updateSelectInput(session, "radar_chr", choices = character(0), selected = NULL)
      updateSelectInput(session, "radar_chr_for_cluster", choices = character(0), selected = NULL)
      return()
    }
    
    chr_labels <- paste0("chr", chr_label_plink(chr_choices))
    chr_map    <- stats::setNames(as.character(chr_choices), chr_labels)
    
    sel_chr  <- if (!is.null(input$radar_chr) && input$radar_chr %in% chr_map) input$radar_chr else chr_map[[1]]
    sel_chr2 <- if (!is.null(input$radar_chr_for_cluster) && input$radar_chr_for_cluster %in% chr_map) input$radar_chr_for_cluster else chr_map[[1]]
    
    updateSelectInput(session, "radar_chr", choices = chr_map, selected = sel_chr)
    updateSelectInput(session, "radar_chr_for_cluster", choices = chr_map, selected = sel_chr2)
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
      dplyr::filter(CHR == chr_sel, !is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::arrange(cluster_id) %>%
      dplyr::pull(cluster_id)
    
    if (!length(ids)) {
      updateSelectInput(session, "radar_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    sel_id <- if (!is.null(input$radar_cluster_id) && input$radar_cluster_id %in% ids) input$radar_cluster_id else ids[1]
    updateSelectInput(session, "radar_cluster_id", choices = ids, selected = sel_id)
  }, ignoreInit = FALSE)
  
  # C) Fill genes grouped by cluster (optgroup)
  observeEvent(radar_base_df(), {
    df <- radar_base_df()
    
    if (!is.data.frame(df) || nrow(df) == 0) {
      updateSelectInput(session, "radar_gene", choices = character(0), selected = NULL)
      return()
    }
    
    dd <- df %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(genename), nzchar(genename)) %>%
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
    sel <- if (!is.null(input$radar_gene) && input$radar_gene %in% all_keys) input$radar_gene else all_keys[1]
    
    updateSelectInput(session, "radar_gene", choices = choices_list, selected = sel)
  }, ignoreInit = FALSE)
  
  # D) Aggregate row depending on mode (cluster / chr / gene)
  radar_group_row <- reactive({
    df <- radar_base_df()
    if (!is.data.frame(df) || nrow(df) == 0) return(NULL)
    
    mets <- intersect(radar_metrics_all, names(df))
    if (length(mets) < 3) return(NULL)
    
    mode <- input$radar_mode %||% "cluster"
    
    if (mode == "cluster") {
      cl <- as.character(input$radar_cluster_id %||% "")
      if (!nzchar(cl)) return(NULL)
      df_use <- df %>% dplyr::filter(cluster_id == cl)
      label  <- cl
      
    } else if (mode == "chr") {
      chr_sel <- suppressWarnings(as.integer(input$radar_chr))
      if (!is.finite(chr_sel)) return(NULL)
      df_use <- df %>% dplyr::filter(CHR == chr_sel)
      label  <- paste0("chr", chr_label_plink(chr_sel))
      
    } else { # gene: "cluster||gene"
      key <- as.character(input$radar_gene %||% "")
      if (!nzchar(key)) return(NULL)
      parts <- strsplit(key, "\\|\\|")[[1]]
      cluster_sel <- parts[1] %||% ""
      gene_sel    <- parts[2] %||% ""
      if (!nzchar(cluster_sel) || !nzchar(gene_sel)) return(NULL)
      
      df_use <- df %>% dplyr::filter(cluster_id == cluster_sel, genename == gene_sel)
      label  <- paste0(gene_sel, " (", cluster_sel, ")")
    }
    
    if (!is.data.frame(df_use) || nrow(df_use) == 0) return(NULL)
    
    row <- df_use %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(mets), \(x) mean(x, na.rm = TRUE)))
    
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
  
  # F) Radar selection table (summary of IDs + genes with links)
  radar_selected_df <- reactive({
    df <- dbnsfp_norm_df()
    validate(need(is.data.frame(df) && nrow(df) > 0, "Final output not available yet."))
    validate(need(all(c("cluster_id", "genename") %in% names(df)),
                  "Final output must contain: cluster_id, genename"))
    
    mode <- input$radar_mode %||% "cluster"
    chr_num <- suppressWarnings(as.integer(norm_chr_generic(df[["chr"]])))
    
    if (mode == "cluster") {
      req(input$radar_cluster_id)
      df_use <- df %>% dplyr::filter(!is.na(cluster_id), cluster_id == input$radar_cluster_id)
      label <- input$radar_cluster_id
      
    } else if (mode == "chr") {
      req(input$radar_chr)
      chr_sel <- suppressWarnings(as.integer(input$radar_chr))
      df_use <- df %>% dplyr::mutate(.CHR = chr_num) %>% dplyr::filter(.CHR == chr_sel)
      label <- paste0("chr", chr_label_plink(chr_sel))
      
    } else {
      req(input$radar_gene)
      parts <- strsplit(as.character(input$radar_gene), "\\|\\|")[[1]]
      cluster_sel <- parts[1] %||% ""
      gene_sel    <- parts[2] %||% ""
      req(nzchar(cluster_sel), nzchar(gene_sel))
      
      df_use <- df %>% dplyr::filter(cluster_id == cluster_sel, genename == gene_sel)
      label  <- paste0(gene_sel, " (", cluster_sel, ")")
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
      options    = list(
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        dom = "tip",
        pageLength = 1,
        scrollX    = TRUE
      )
    )
  })
  
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
    
    out %>% plotly::config(displaylogo = FALSE)
    
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
    dt <- nonsyn_manhattan_df()
    req(is.data.frame(dt), nrow(dt) > 0)
    req("cluster_id" %in% names(dt))
    
    metric <- input$nonsyn_metric
    
    dt <- dt %>%
      dplyr::filter(!is.na(cluster_id), nzchar(cluster_id), !is.na(CHR), is.finite(value))
    req(nrow(dt) > 0)
    
    bg_values <- dt$value
    
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE
      )
    )
  })
  
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 15,
        scrollX    = TRUE
      )
    )
  })
  
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 15,
        scrollX    = TRUE
      )
    )
  })
  
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
  
  # Cluster selector UI (single definition)
  output$func_cluster_ui <- renderUI({
    req(identical(input$func_scope, "cluster"))
    
    cc <- cluster_gene_counts()
    validate(need(nrow(cc) > 0, "No clusters available."))
    
    cc_ok <- cc %>% dplyr::filter(n_genes >= 5)
    validate(need(nrow(cc_ok) > 0, "No clusters with ≥ 5 genes available for enrichment."))
    
    choices <- stats::setNames(
      cc_ok$cluster_id,
      paste0(cc_ok$cluster_id, " (n=", cc_ok$n_genes, " genes)")
    )
    
    sel <- input$func_cluster_id
    if (is.null(sel) || !(sel %in% cc_ok$cluster_id)) sel <- cc_ok$cluster_id[1]
    
    selectInput("func_cluster_id", "Cluster (≥ 5 genes)", choices = choices, selected = sel)
  })
  
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
    req(identical(input$func_scope, "cluster"))
    
    cc <- cluster_gene_counts()
    validate(need(nrow(cc) > 0, "No clusters available."))
    
    cc_ok <- cc %>% dplyr::filter(n_genes >= 5)
    validate(need(nrow(cc_ok) > 0,
                  "No clusters with ≥ 5 genes are available for enrichment."))
    
    # Pretty label: cluster_id (n=XX genes)
    choices <- stats::setNames(
      cc_ok$cluster_id,
      paste0(cc_ok$cluster_id, " (n=", cc_ok$n_genes, " genes)")
    )
    
    # Keep selection if still valid
    sel <- input$func_cluster_id
    if (is.null(sel) || !(sel %in% cc_ok$cluster_id)) sel <- cc_ok$cluster_id[1]
    
    selectInput(
      "func_cluster_id",
      "Cluster (≥ 5 genes)",
      choices  = choices,
      selected = sel
    )
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
    bg  <- input$enrich_background %||% "dataset"
    tab <- input$enrich_tabs %||% "tab_enrich_go"  # tabsetPanel id
    
    bg_lbl <- if (identical(bg, "dataset")) "Dataset genes" else "Default background"
    
    msg <- if (identical(tab, "tab_enrich_go")) {
      
      if (identical(bg, "orgdb")) {
        paste0(
          "<b>Background:</b> ", bg_lbl, "<br>",
          "GO enrichment uses the default <i>OrgDb</i>-based background (org.Hs.eg.db; ontology-specific)."
        )
      } else {
        paste0(
          "<b>Background:</b> ", bg_lbl, "<br>",
          "GO enrichment uses your dataset-derived universe: genes present in the uploaded dataset (after SYMBOL→ENTREZ mapping)."
        )
      }
      
    } else if (identical(tab, "tab_enrich_kegg")) {
      
      if (identical(bg, "orgdb")) {
        paste0(
          "<b>Background:</b> ", bg_lbl, "<br>",
          "KEGG enrichment uses the default KEGG organism background (<i>hsa</i>). This background can differ from GO."
        )
      } else {
        paste0(
          "<b>Background:</b> ", bg_lbl, "<br>",
          "KEGG enrichment uses your dataset-derived universe: genes present in the uploaded dataset (after SYMBOL→ENTREZ mapping)."
        )
      }
      
    } else {
      paste0("<b>Background:</b> ", bg_lbl)
    }
    
    helpText(HTML(paste0("<span style='color:#555;'>", msg, "</span>")))
  })
  
  
  # ============================================================
  # EVENT: run enrichment (separate GO vs KEGG triggers)
  # ============================================================
  
  go_trigger   <- reactiveVal(0)
  kegg_trigger <- reactiveVal(0)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    if (identical(tab, "tab_enrich_go")) {
      go_trigger(go_trigger() + 1L)
    } else if (identical(tab, "tab_enrich_kegg")) {
      kegg_trigger(kegg_trigger() + 1L)
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
  
  output$go_bar <- renderPlot({
    dat <- go_top_df()
    req(is.data.frame(dat), nrow(dat) > 0)
    
    fb    <- isTRUE(attr(dat, "fallback"))
    ttl   <- paste0("GO enrichment (", scope_label(), ")")
    subttl <- if (fb) {
      paste0(
        "No terms at FDR ≤ ", input$go_kegg_pcut %||% 0.05,
        " — showing top terms ranked by FDR; bar height = -log10(pvalue)"
      )
    } else {
      paste0("FDR ≤ ", input$go_kegg_pcut %||% 0.05, " (bar height = -log10(pvalue))")
    }
    
    dat <- dat %>%
      dplyr::mutate(Ontology = factor(as.character(Ontology), levels = c("BP","CC","MF")))
    
    cols_go <- c(BP = "darkgreen", CC = "orange", MF = "darkblue")
    
    mk <- function(onto) {
      d <- dat %>% dplyr::filter(Ontology == onto)
      if (!nrow(d)) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = paste0(onto, " (none)")))
      }
      
      ggplot2::ggplot(
        d,
        ggplot2::aes(
          x    = stats::reorder(term_short, score),
          y    = score,
          fill = Ontology
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(title = onto, x = NULL, y = "-log10(pvalue)") +
        ggplot2::scale_fill_manual(values = cols_go, guide = "none", drop = FALSE) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    }
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
      (mk("BP") | mk("CC") | mk("MF")) +
        patchwork::plot_annotation(title = ttl, subtitle = subttl)
    } else {
      ggplot2::ggplot(
        dat,
        ggplot2::aes(
          x    = stats::reorder(term_short, score),
          y    = score,
          fill = Ontology
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~Ontology, scales = "free_y") +
        ggplot2::scale_fill_manual(values = cols_go, guide = "none", drop = FALSE) +
        ggplot2::labs(title = ttl, subtitle = subttl, x = NULL, y = "-log10(pvalue)") +
        ggplot2::theme_minimal(base_size = 12)
    }
  }, height = 450)
  
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
        buttons    = c("copy", "csv", "excel","pdf"),
        # FDR is column 9 in the displayed table (1-based index inside DataTables)
        order      = list(list(8, "asc"))
      )
    ) %>%
      DT::formatSignif(c("pvalue", "FDR", "qvalue"), digits = 3)
  })
  
  
  # ==========================
  # KEGG enrichment (with fallback when no terms pass FDR)
  # ==========================
  
  kegg_raw <- eventReactive(kegg_trigger(), {
    
    withProgress(message = "Running KEGG enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.3, detail = paste("Genes:", length(gene)))
      
      uni  <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (KEGG)" else paste("Universe:", length(uni)))
      
      ek <- clusterProfiler::enrichKEGG(
        gene         = gene,
        universe     = uni,
        organism     = "hsa",
        pvalueCutoff = 1
      )
      
      validate(need(!is.null(ek) && !is.null(ek@result) && nrow(ek@result) > 0,
                    "No KEGG enrichment results."))
      
      incProgress(0.4, detail = "Setting readable…")
      suppressWarnings({
        ek <- clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      })
      
      incProgress(0.2, detail = "Done")
      ek
    })
    
  }, ignoreInit = TRUE)
  
  
  kegg_top_df <- reactive({
    ek <- kegg_raw()
    validate(need(!is.null(ek) && nrow(ek@result) > 0, "No KEGG enrichment results."))
    
    pcut <- input$go_kegg_pcut %||% 0.05
    topn <- input$kegg_top %||% 15
    
    df <- as.data.frame(ek@result, stringsAsFactors = FALSE)
    
    # Ensure p.adjust exists
    if (!"p.adjust" %in% names(df) && "pvalue" %in% names(df)) {
      df$p.adjust <- p.adjust(df$pvalue, method = "BH")
    }
    
    df_fdr <- df %>%
      dplyr::filter(!is.na(p.adjust)) %>%
      dplyr::filter(p.adjust <= pcut) %>%
      dplyr::arrange(p.adjust)
    
    # Fallback if nothing passes FDR
    if (nrow(df_fdr) == 0) {
      if ("p.adjust" %in% names(df) && any(is.finite(df$p.adjust))) {
        df_use <- df %>%
          dplyr::filter(is.finite(p.adjust)) %>%
          dplyr::arrange(p.adjust) %>%
          dplyr::slice_head(n = topn)
        attr(df_use, "fallback_mode") <- paste0("No KEGG with FDR ≤ ", pcut, " → showing top ", topn, " by FDR")
      } else {
        df_use <- df %>%
          dplyr::filter(is.finite(pvalue)) %>%
          dplyr::arrange(pvalue) %>%
          dplyr::slice_head(n = topn)
        attr(df_use, "fallback_mode") <- paste0("No KEGG with FDR ≤ ", pcut, " → showing top ", topn, " by p-value")
      }
      validate(need(nrow(df_use) > 0, "No KEGG pathways to display (even by p-value)."))
      return(df_use)
    }
    
    df_use <- df_fdr %>% dplyr::slice_head(n = topn)
    attr(df_use, "fallback_mode") <- ""
    df_use
  })
  
  
  output$kegg_bar <- renderPlot({
    
    df <- kegg_top_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    note <- attr(df, "fallback_mode") %||% ""
    ttl  <- paste0(
      "KEGG enrichment (", scope_label(), ")",
      if (nzchar(note)) paste0(" — ", note) else paste0(" — FDR ≤ ", input$go_kegg_pcut %||% 0.05)
    )
    
    df <- df %>% dplyr::mutate(score = -log10(pvalue))
    ylab <- "-log10(p-value)"
    
    # Use a different fill per pathway (stable order)
    fill_id <- if ("ID" %in% names(df)) as.character(df$ID) else as.character(df$Description)
    df <- df %>% dplyr::mutate(fill_id = factor(fill_id, levels = fill_id))
    
    n   <- nlevels(df$fill_id)
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
    names(pal) <- levels(df$fill_id)
    
    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x    = stats::reorder(Description, score),
        y    = score,
        fill = fill_id
      )
    ) +
      ggplot2::geom_col(color = "black", linewidth = 0.25) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::labs(title = ttl, x = NULL, y = ylab) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
  }, height = function() {
    n <- tryCatch(nrow(kegg_top_df()), error = function(e) 0)
    max(420, 28 * n)
  })
  
  output$kegg_table <- DT::renderDT({
    df <- kegg_top_df()
    req(is.data.frame(df), nrow(df) > 0)
    
    out <- df %>%
      dplyr::mutate(
        KEGG = sprintf("<a href='https://www.kegg.jp/pathway/%s' target='_blank'>%s</a>", ID, ID),
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
        buttons    = c("copy", "csv", "excel","pdf")
      )
    ) %>%
      DT::formatSignif(c("pvalue", "FDR", "qvalue"), digits = 3)
  })
  
  
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
    
    chr_sel <- suppressWarnings(as.integer(input$go_chr_for_cluster))
    if (!is.finite(chr_sel)) {
      updateSelectInput(session, "go_cluster_id", choices = character(0), selected = NULL)
      return()
    }
    
    ids <- df %>%
      dplyr::filter(CHR == chr_sel, !is.na(cluster_id), nzchar(cluster_id)) %>%
      dplyr::distinct(cluster_id) %>%
      dplyr::arrange(cluster_id) %>%
      dplyr::pull(cluster_id)
    
    updateSelectInput(session, "go_cluster_id", choices = ids, selected = ids[1] %||% NULL)
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
  # ==========================================================
  # GO slim classification (BP, CC, MF)
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
  
  
  output$go_class_plot <- renderPlotly({
    
    dat <- go_slim_class_df_selected()
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
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO slim:</b> ", htmltools::htmlEscape(as.character(dat$slim_term)),
      "<br><b>Genes:</b> ", dat$n_genes
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
        y = "Number of genes",
        title = paste0("Gene Function Classification (GO slim generic) — ", sel_label)
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
    
    plotly::ggplotly(p, tooltip = "text")
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE,
        autoWidth  = TRUE,
        dom        = "tip",
        columnDefs = list(
          list(targets = c(3,4,5), width = "380px")
        )
      )
    ) %>%
      DT::formatStyle(
        columns = c("GO_biological_process","GO_cellular_component","GO_molecular_function"),
        `white-space` = "normal"
      )
  })
  
  
  output$go_debug <- renderPrint({
    syms <- tryCatch(genes_scope_symbols(), error = function(e) character(0))
    ids  <- tryCatch(entrez_scope(),        error = function(e) character(0))
    uni  <- tryCatch(universe_entrez_for_scope(), error = function(e) NULL)
    
    list(
      scope     = input$func_scope %||% "global",
      n_symbols = length(unique(syms)),
      n_entrez  = length(unique(ids)),
      universe  = if (is.null(uni)) "Default background" else length(unique(uni))
    )
  })
  
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
      need(all(c("genename","cluster_id","Function_description","Disease_description",
                 "cluster_chr_n","cluster_start","cluster_end") %in% names(df)),
           "Missing required cluster/function/disease columns in final output.")
    )
    
    # Detectar columna OMIM (acepta OMIM_id o MIM_phenotype_id)
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
    
    make_omim_url <- function(x) {
      if (is.null(x) || is.na(x)) return(NA_character_)
      s <- trimws(as.character(x))
      if (!nzchar(s) || s %in% c(".", "NA")) return(NA_character_)
      
      # Puede venir "616573;614123" o "616573" → consideramos el primero
      ids <- unlist(strsplit(s, "[,;|\\s]+"))
      ids <- ids[nzchar(ids)]
      if (!length(ids)) return(NA_character_)
      id1 <- gsub("[^0-9]", "", ids[1])
      if (!nzchar(id1)) return(NA_character_)
      
      paste0("https://omim.org/entry/", id1)
    }
    
    make_genecards_url <- function(g) {
      if (is.null(g) || is.na(g)) return(NA_character_)
      s <- trimws(as.character(g))
      if (!nzchar(s) || s %in% c(".", "NA")) return(NA_character_)
      paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", utils::URLencode(s, reserved = TRUE))
    }
    
    out <- df %>%
      dplyr::transmute(
        genename    = as.character(genename),
        cluster_id  = as.character(cluster_id),
        cluster_chr_n = as.character(cluster_chr_n),
        cluster_start = suppressWarnings(as.integer(cluster_start)),
        cluster_end   = suppressWarnings(as.integer(cluster_end)),
        Function_description = dplyr::na_if(trimws(as.character(Function_description)), ""),
        Disease_description  = dplyr::na_if(trimws(as.character(Disease_description)), ""),
        OMIM_raw = .data[[omim_col]]
      ) %>%
      dplyr::filter(!is.na(cluster_id) & nzchar(cluster_id)) %>%
      dplyr::distinct(
        cluster_id, cluster_chr_n, cluster_start, cluster_end, genename,
        Function_description, Disease_description, OMIM_raw
      ) %>%
      dplyr::mutate(
        UCSC = {
          u <- mapply(make_ucsc_url, cluster_chr_n, cluster_start, cluster_end,
                      MoreArgs = list(db = "hg38"))
          ifelse(is.na(u), "", sprintf("<a href='%s' target='_blank'>Open UCSC</a>", u))
        },
        OMIM = {
          o <- vapply(OMIM_raw, make_omim_url, character(1))
          ifelse(is.na(o), "", sprintf("<a href='%s' target='_blank'>Open OMIM</a>", o))
        },
        genename = {
          g <- vapply(genename, make_genecards_url, character(1))
          ifelse(is.na(g), genename, sprintf("<a href='%s' target='_blank'>%s</a>", g, genename))
        }
      ) %>%
      dplyr::select(UCSC, OMIM, genename, cluster_id, Function_description, Disease_description)
    
    
    out
  })
  
  output$func_disease_tbl <- DT::renderDT({
    dat <- func_disease_df()
    req(nrow(dat) > 0)
    
    DT::datatable(
      dat,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 15,
        scrollX    = TRUE
      )
    )
  })
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
  
  # -----------------------------
  # LD PANEL (NonSyn -> shared mod_ld.R)
  # -----------------------------
  clusters_r <- reactive({
    build_ld_clusters_from_app_nonsyn()
  })
  
  candidates_r <- reactive({
    build_ld_candidates_from_app_nonsyn()
  })
  
  ld_module_server(
    "ld",
    clusters_r   = clusters_r,
    candidates_r = candidates_r
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
  
  
  
}  # end server

shinyApp(ui, server)

