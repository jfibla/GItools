## app.R — GTEx eQTL Inspector (GWAS + GTEx clusters) GTEX Inspector App Patched
## Adaptat a l’esquema “canònic”:
## input GWAS -> 

## + n_gtex a Cluster summary i GTEx table amb cluster + cluster_chr_n
## + downloads amb nom incloent mode (w/h) i threshold

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

# Gene model (hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(IRanges)

# (opcional però recomanat per labels)
library(AnnotationDbi)
library(org.Hs.eg.db)

library("this.path")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# --- Fixar dplyr com a select/mutate/arrange oficial ---
select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
summarise  <- dplyr::summarise
rename     <- dplyr::rename
# ===========================
# ------------------------------------------------------------
# APP_DIR robust (no depèn de getwd)
# ------------------------------------------------------------
APP_DIR <- tryCatch({
  # quan s'executa amb runApp("path"), this.path funciona; si no, fallback
  if (requireNamespace("this.path", quietly = TRUE)) {
    normalizePath(dirname(this.path::this.path()), winslash = "/", mustWork = FALSE)
  } else {
    ""
  }
}, error = function(e) "")

# Fallback si this.path no està disponible:
if (!nzchar(APP_DIR)) {
  
  cand <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(cand, "R")) && dir.exists(file.path(cand, "www"))) APP_DIR <- cand
}

cat("[GTEX] APP_DIR=", APP_DIR, "\n")

# ===========================

# ===========================
# GTEx file logging (early)
# ===========================
GTEX_LOG_DIR <- Sys.getenv("GITOOLS_LOG_DIR", "")
if (!nzchar(GTEX_LOG_DIR)) {
  GTEX_LOG_DIR <- file.path(dirname(getwd()), "_logs")  # -> GItools/app/_logs
}
dir.create(GTEX_LOG_DIR, showWarnings = FALSE, recursive = TRUE)
GTEX_LOG <- file.path(GTEX_LOG_DIR, "gtex_runtime.log")

glog <- function(...) {
  cat(paste0(..., collapse = ""), "\n", file = GTEX_LOG, append = TRUE)
}

glog("[GTEx] boot @ ", format(Sys.time()), " wd=", getwd())


# --- GI shared state (PORTABLE root) ---

# 1) Guess GItools root from current working dir (runApp usually sets wd = app dir)
gi_root_guess <- tryCatch(
  normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE),
  error = function(e) ""
)

# 2) If <root>/config.R exists, source it (so gi_cfg() becomes available)
if (nzchar(gi_root_guess)) {
  cfg_file0 <- file.path(gi_root_guess, "config.R")
  if (file.exists(cfg_file0)) source(cfg_file0, local = TRUE)
}

# 3) Final root: prefer gi_cfg() if present; else fallback to guess

BASE <- Sys.getenv("GITOOLS_ROOT", "")
if (!nzchar(BASE)) {
  BASE <- if (exists("gi_cfg", mode="function")) gi_cfg()$root else gi_root_guess
}

BASE <- normalizePath(BASE, winslash="/", mustWork = FALSE)
if (!nzchar(BASE) || !dir.exists(BASE)) {
  glog("[GTEx] ERROR: BASE not found. wd=", getwd(), " gi_root_guess=", gi_root_guess, " env_GITOOLS_ROOT=", Sys.getenv("GITOOLS_ROOT",""))
  stop("GItools BASE not found (check HUB working dir/env vars).")
}

# 4) Load GI shared state (defines gi_shared_root / gi_state_root)
state_file <- file.path(BASE, "_shared", "gi_state.R")
stopifnot(file.exists(state_file))
source(state_file, local = TRUE)

cat("[GI] using gi_state.R at:", state_file, "\n")
cat("[GI] gi_shared_root =", gi_shared_root, "\n")
cat("[GI] gi_state_root  = ", gi_state_root(), "\n")


# --- LD module (portable via _shared) ---
ld_file <- file.path(gi_shared_root, "mod_ld.R")
stopifnot(file.exists(ld_file))
source(ld_file, local = TRUE)
stopifnot(exists("ld_module_ui"), exists("ld_module_server"))

# ===========================
# Helpers generals
# ===========================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

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

# detectar columna rsID generica
detect_rsid_col <- function(df) {
  cn <- names(df)
  patt <- "(^snp$|^snp\\.|^rsid$|^rsid\\.|marker|id|variant)"
  hit <- grep(patt, cn, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) return(hit[1])
  return(NULL)
}

# referència hg38 per cumulatiu
.chr_lengths_hg38 <- function() {
  lens <- c(
    `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555, `5`=181538259,
    `6`=170805979, `7`=159345973, `8`=145138636, `9`=138394717, `10`=133797422,
    `11`=135086622, `12`=133275309, `13`=114364328, `14`=107043718, `15`=101991189,
    `16`=90338345, `17`=83257441, `18`=80373285, `19`=58617616, `20`=64444167,
    `21`=46709983, `22`=50818468, `X`=156040895, `Y`=57227415, `MT`=16569
  )
  ord <- c(as.character(1:22),"X","Y","MT")
  df <- data.frame(chr = factor(names(lens), levels = ord), len = as.numeric(lens))
  df$chr_cum <- cumsum(df$len) - df$len
  df$center  <- df$chr_cum + df$len/2
  df
}
.ref_hg38 <- .chr_lengths_hg38()

coord_from_bp_cum <- function(bp_cum) {
  ref <- .ref_hg38
  idx <- which(ref$chr_cum <= bp_cum)
  if (length(idx) == 0) {
    chr <- as.character(ref$chr[1])
    pos <- round(bp_cum)
  } else {
    idx <- max(idx)
    chr <- as.character(ref$chr[idx])
    pos <- round(bp_cum - ref$chr_cum[idx])
  }
  pos <- max(pos, 1)
  list(chr = chr, pos = pos)
}

# downloads: tag mode + threshold
.make_mode_thr_tag <- function(input) {
  mode_tag <- if (identical(input$cluster_method, "window")) "w" else "h"
  thr_val <- if (identical(input$cluster_method, "window")) input$pthr else input$min_logp
  thr_tag <- gsub("\\.", "p", sprintf("%.2f", thr_val))
  list(mode_tag = mode_tag, thr_tag = thr_tag)
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

# -----------------------------
# Add per-cluster GWAS stats from df_sig (hits) to clusters
# clusters: data.frame with chr, start, end, cluster_id
# df_sig:  data.frame with CHR, BP, snp, logp
# -----------------------------
add_cluster_stats_from_hits <- function(clusters_df, df_sig) {
  library(data.table)
  
  cl <- as.data.table(clusters_df)
  stopifnot(all(c("chr","start","end") %in% names(cl)))
  
  # Accept either cluster_id (Catalog-like) OR cluster_chr_n (GTEx-like)
  if (!"cluster_id" %in% names(cl)) {
    if ("cluster_chr_n" %in% names(cl)) {
      cl[, cluster_id := as.character(cluster_chr_n)]
    } else {
      # last resort stable id
      cl[, cluster_id := paste0("cluster_", seq_len(.N))]
    }
  } else {
    cl[, cluster_id := as.character(cluster_id)]
  }
  
  hits <- as.data.table(df_sig)
  stopifnot(all(c("CHR","BP","logp") %in% names(hits)))
  if (!"snp" %in% names(hits)) hits[, snp := as.character(rsid %||% NA_character_)]
  
  # Types
  cl[, `:=`(
    chr   = as.integer(chr),
    start = as.integer(start),
    end   = as.integer(end)
  )]
  hits[, `:=`(
    chr  = as.integer(CHR),
    pos  = as.integer(BP),
    snp  = as.character(snp),
    logp = suppressWarnings(as.numeric(logp))
  )]
  
  cl <- cl[is.finite(chr) & is.finite(start) & is.finite(end)]
  hits <- hits[is.finite(chr) & is.finite(pos)]
  
  # Prepare intervals for overlaps
  cl_iv   <- cl[, .(chr, start, end, cluster_id)]
  hits_iv <- hits[, .(chr, start = pos, end = pos, pos, snp, logp)]
  
  setkey(cl_iv, chr, start, end)
  setkey(hits_iv, chr, start, end)
  
  ov <- foverlaps(hits_iv, cl_iv, type = "within", nomatch = 0L)
  
  # Default columns (in case of 0 overlaps)
  cl[, `:=`(
    n_snps   = 0L,
    center   = as.numeric((start + end) / 2),
    top_snp  = NA_character_,
    top_logp = NA_real_
  )]
  
  if (nrow(ov) > 0) {
    stats <- ov[, .(
      n_snps = .N,
      center = as.numeric(mean(pos, na.rm = TRUE)),
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
      n_snps   = i.n_snps,
      center   = i.center,
      top_snp  = i.top_snp,
      top_logp = i.top_logp
    )]
  }
  
  cl[, cluster_size_kb := round((end - start) / 1000, 2)]
  
  as.data.frame(cl)
}
# ------------------------------------------------------------
#### Helpers ORA + mapa gene_id→gene_name
# ------------------------------------------------------------
# ORA Fisher per files (eQTL rows)
# ------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

ora_fisher_rows <- function(universe_vec, foreground_vec) {
  U <- na.omit(as.character(universe_vec))
  F <- na.omit(as.character(foreground_vec))
  
  if (!length(F) || !length(U)) {
    return(data.table::data.table(
      term=character(), k=integer(), K=integer(), n=integer(), N=integer(),
      odds=numeric(), pvalue=numeric(), padj=numeric()
    ))
  }
  
  N <- length(U)
  n <- length(F)
  
  tabU <- table(U)
  tabF <- table(F)
  
  terms <- names(tabU)
  
  out <- data.table::rbindlist(lapply(terms, function(t) {
    K <- as.integer(tabU[[t]])
    k <- as.integer(tabF[[t]] %||% 0L)
    
    a <- k
    b <- n - k
    c <- K - k
    d <- (N - K) - (n - k)
    if (any(c(a,b,c,d) < 0)) return(NULL)
    
    ft <- stats::fisher.test(matrix(c(a,c,b,d), nrow = 2), alternative = "greater")
    data.table::data.table(
      term=t, k=a, K=K, n=n, N=N,
      odds=unname(ft$estimate), pvalue=ft$p.value
    )
  }), fill = TRUE)
  
  
  if (!nrow(out)) return(out)
  out <- data.table::data.table(term = terms, a = a, Tt = Tt, or = or, p = p)
  out[, p_adj := p.adjust(p, method = "BH")]
  out[, `:=`(n_total = n, N_total = N)]   # NEW
  data.table::setorder(out, p_adj, p)
  out
}
##############################################################################
# ------------------------------------------------------------
# Helpers label (Catalog-style)
# ------------------------------------------------------------
short_label <- function(x, max_chars = 30, wrap_width = 12) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
  stringr::str_wrap(x, width = wrap_width)
}

pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[1] else NULL
}
##############################################################################
################################################
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

##################################################
##############################################################################
############################### UI ###########################################
##############################################################################
# ===========================
# UI (FULL) — TOP NAV: Analysis / Enrichment / LD
# - Global grey background for plots + tables
# - Enrichment moved OUT of Analysis main_tabs (now its own navbar tab)
# ===========================

ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("🧠 GTEx eQTL Inspector")
  ),
  id = "topnav",
  
  # -----------------------------
  # Global QS push + global grey styling
  # -----------------------------
  tags$head(
    tags$script(HTML(
      "document.addEventListener('DOMContentLoaded', function(){
       if (window.Shiny) {
         Shiny.setInputValue('gi_qs', window.location.search || '', {priority: 'event'});
         Shiny.setInputValue('gi_href', window.location.href || '', {priority: 'event'});
       }
     });"
    )),
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

/* Plotly widgets */
.plotly.html-widget {
  background: #f2f2f2 !important;
  border-radius: 10px;
  padding: 10px;
}

/* Base R plots / ggplot outputs */
.shiny-plot-output {
  background: #f2f2f2 !important;
  border-radius: 10px;
  padding: 10px;
}

/* DT tables (wrapper container) */
.dataTables_wrapper {
  background: #f2f2f2 !important;
  border-radius: 10px;
  padding: 10px;
}

/* Keep DT inner bits cohesive (optional) */
table.dataTable,
table.dataTable thead th,
table.dataTable tbody td {
  background-color: transparent !important;
}

/* Spinner wrapper container */
.shiny-spinner-output-container {
  background: #f2f2f2 !important;
  border-radius: 10px;
  padding: 10px;
}

/* Optional: panel-lite shouldn't fight the global padding */
.panel-lite {
  background: transparent !important;
}
    "))
  ),
  
  # ========================================================================
  # NAV TAB 1: Analysis
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    
    sidebarLayout(
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
          # Step 1: GWAS input
          # -----------------------------
          h3(
            "Step 1 · Load GWAS hits",
            style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
          ),
          fluidRow(
            column(6, actionButton("info_00", "ℹ️ file format")),
            column(
              6,
              actionButton(
                "reset_case",
                "Reset",
                icon = icon("rotate-left"),
                class = "btn-warning"
              )
            )
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
          # Step 2: Clustering
          # -----------------------------
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
            "➊ Generate intervals + clusters",
            style = "background-color: #ffdd57; color: black; font-weight: bold;"
          ),
          h4("Preview clusters (derived from intervals)"),
          div(class = "panel-lite", verbatimTextOutput("ranges_preview")),
          tags$hr()
        ),
        
        # -----------------------------
        # Step 3: ALWAYS visible (even in Hub mode)
        # -----------------------------
        h3(
          "Step 3 · Extract eQTLs in all clusters",
          style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
        ),
        textOutput("gtex_status"),
        actionButton(
          "run_gtex",
          "➋ Extract GTEx eQTLs from ALL clusters",
          style = "background-color: #ffdd57; color: black; font-weight: bold;"
        ),
        textOutput("gtex_summary"),
        tags$hr(),
        downloadButton("dl_gtex_hits_csv", "Download GTEx hits (by cluster) CSV"),
        downloadButton("dl_gtex_hits_rds", "Download GTEx hits (by cluster) RDS"),
        tags$hr(),
        downloadButton("dl_candidates_zip", "⬇️ Download gtex candidates (ZIP)")
      ),
      
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            h4("Top: GWAS p-values · Bottom: GTEx eQTLs inside all clusters"),
            div(class = "panel-lite", shinycssloaders::withSpinner(plotlyOutput("manhattan_combo", height = "700px"))),
            tags$hr(),
            helpText("Select a window on 'GTEx hits' to obtain a link to UCSC browser (see below)"),
            
            fluidRow(
              column(6, div(class = "panel-lite", uiOutput("debug_ucsc_state"))),
              column(2, h4("UCSC links to expanded region:")),
              column(
                4,
                div(
                  style = "margin-top:8px; background-color:#f2f2f2; padding:10px; border-radius:8px; width:100%;",
                  uiOutput("ucsc_link_gwas"),
                  uiOutput("ucsc_link_gtex")
                )
              )
            ),
            
            h4("Clusters summary"),
            div(class = "panel-lite", DTOutput("cluster_dt"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔬 GTEX hits</span>"),
            h4("GTEx eQTLs inside all clusters"),
            div(class = "panel-lite", DTOutput("gtex_table"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔬 GWAS hits</span>"),
            h4("GWAS hits used for clustering (window mode)"),
            div(class = "panel-lite", DTOutput("hits_tbl"))
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🎨 GTEx visualizations</span>"),
            fluidRow(
              column(
                6,
                h4("GTEx hits by tissue families and slope sign"),
                div(class = "panel-lite", plotlyOutput("tissue_family_bar", height = "350px"))
              ),
              column(
                6,
                h4("Advanced dotplot (tissue × -log10(p))"),
                div(class = "panel-lite", plotlyOutput("tissue_heatmap", height = "350px"))
              )
            ),
            tags$hr(),
            fluidRow(
              column(
                8,
                h4("Functional heatmap: genes × tissue families"),
                div(class = "panel-lite", plotlyOutput("tissue_dotplot", height = "450px"))
              ),
              column(
                4,
                h4("Automatic interpretation"),
                div(class = "panel-lite", verbatimTextOutput("gtex_interpret"))
              )
            )
          )
        )
      )
    )
  ),
  
  # ========================================================================
  # NAV TAB 2: Enrichment (GO / KEGG / GoSlim)
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>✳️ Enrichment</span>"),
    value = "tab_enrich",
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        radioButtons(
          "func_scope", "Scope",
          choices  = c("Global" = "global", "Cluster" = "cluster", "Gene" = "gene"),
          selected = "global"
        ),
        conditionalPanel(
          condition = "input.func_scope == 'cluster'",
          selectInput("func_cluster_id", "Cluster", choices = NULL)
        ),
        conditionalPanel(
          condition = "input.func_scope == 'gene'",
          uiOutput("func_gene_ui")
        ),
        tags$hr(),
        
        selectInput(
          "enrich_background", "Background",
          choices  = c("Reference annotated genes" = "orgdb",
                       "Dataset genes"            = "dataset"),
          selected = "orgdb"
        ),
        
        numericInput(
          "enrich_pcut",
          "FDR cutoff",
          value = 0.05,
          min   = 0,
          max   = 1,
          step  = 0.01
        ),
        
        conditionalPanel(
          condition = "input.enrich_tabs == 'tab_enrich_go' || input.enrich_tabs == 'tab_enrich_goslim'",
          checkboxGroupInput(
            "go_ontos",
            "GO ontologies",
            choices  = c("BP", "CC", "MF"),
            selected = c("BP", "CC", "MF"),
            inline   = TRUE
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
          icon  = icon("play"),
          class = "btn btn-primary"
        )
      ),
      
      mainPanel(
        width = 9,
        
        div(class = "panel-lite", uiOutput("enrich_bg_note")),
        
        tabsetPanel(
          id = "enrich_tabs",
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>GTEx Tissue</span>"),
            value = "tab_enrich_gtex_terms",
            fluidRow(
              column(12, h4("Tissue enrichment"), shinycssloaders::withSpinner(plotOutput("enrich_tissue_plot", height=400)))
            ),
            tags$hr(),
            fluidRow(
              column(12, DT::DTOutput("enrich_tissue_table"))
            )
          ),
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>GO</span>"),
            value = "tab_enrich_go",
            div(class = "panel-lite",
                shinycssloaders::withSpinner(plotlyOutput("go_bar", height = 400))
            ),
            tags$hr(),
            div(class = "panel-lite",
                shinycssloaders::withSpinner(DT::DTOutput("go_table"))
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>KEGG</span>"),
            value = "tab_enrich_kegg",
            div(class = "panel-lite",
                shinycssloaders::withSpinner(plotlyOutput("kegg_bar", height = 400))
            ),
            tags$hr(),
            div(class = "panel-lite",
                shinycssloaders::withSpinner(DT::DTOutput("kegg_table"))
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:15px; font-weight:600;'>GoSlim</span>"),
            value = "tab_enrich_goslim",
            div(class = "panel-lite",
                shinycssloaders::withSpinner(plotlyOutput("goslim_bar", height = 450))
            ),
            tags$hr(),
            div(class = "panel-lite",
                shinycssloaders::withSpinner(DT::DTOutput("goslim_table"))
            )
          )
        )
      )
    )
  ),
  
  # ========================================================================
  # NAV TAB 3: LD
  # ========================================================================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 LD</span>"),
    ld_module_ui("ld", app_tag = "gtex")
  )
)


##############################################################################
############################### SERVER #######################################
##############################################################################

server <- function(input, output, session){
  
  # ------------------------------------------------------------
  # Canonical roots for SLAVE sync (MUST match MASTER)
  #   state : <GItools>/app/_state
  #   shared: <GItools>/_shared
  # ------------------------------------------------------------
  STATE_ROOT  <- normalizePath(file.path(BASE, "app", "_state"), winslash = "/", mustWork = FALSE)
  SHARED_ROOT <- normalizePath(file.path(BASE, "_shared"),       winslash = "/", mustWork = FALSE)
  
  options(gi_state_root  = STATE_ROOT)
  options(gi_shared_root = SHARED_ROOT)
  
  dir.create(STATE_ROOT,  recursive = TRUE, showWarnings = FALSE)
  
  cat("[GTEX][ROOTS] BASE      =", BASE, "\n")
  cat("[GTEX][ROOTS] STATE_ROOT =", STATE_ROOT, " exists=", dir.exists(STATE_ROOT), "\n")
  cat("[GTEX][ROOTS] SHARED_ROOT=", SHARED_ROOT," exists=", dir.exists(SHARED_ROOT), "\n")
  #-----------------------------------------------
  
  # SHARED <- gi_shared_root  # definit a gi_state.R
  SHARED <- getOption("gi_shared_root")
  # -----------------------------
  # Shared reactive holders (Hub OR local)
  # -----------------------------
  gwas_shared_rv <- shiny::reactiveVal(NULL)
  
  # -----------------------------
  # Canonical cluster engine (ALWAYS available)
  # -----------------------------
  # ------------------------------------------------------------
  # Resolve SHARED root robust (Hub OR standalone OR ngroc layout)
  # ------------------------------------------------------------
  SHARED <- NULL
  
  # 1) canonical: gi_shared_root (from gi_state.R) if exists
  if (exists("gi_shared_root", inherits = TRUE)) {
    SHARED <- tryCatch(get("gi_shared_root", inherits = TRUE), error = function(e) NULL)
  }
  
  # 2) option
  if (is.null(SHARED) || !is.character(SHARED) || !nzchar(SHARED)) {
    SHARED <- getOption("gi_shared_root", "")
  }
  
  # 3) env var
  if (!nzchar(SHARED)) {
    SHARED <- Sys.getenv("GITOOLS_SHARED", "")
  }
  
  # 4) fallback: relative to app dir (works if repo has _shared next to app/)
  if (!nzchar(SHARED)) {
    # getwd() when running app is typically the app folder
    cand <- normalizePath(file.path(getwd(), "..", "_shared"), winslash = "/", mustWork = FALSE)
    if (dir.exists(cand)) SHARED <- cand
  }
  
  # Final check + source canonical
  canon_file <- file.path(SHARED, "gi_clusters_canonical.R")
  validate(need(nzchar(SHARED) && file.exists(canon_file),
                paste0("Missing gi_clusters_canonical.R. Resolved SHARED='", SHARED, "'")))
  source(canon_file, local = TRUE)
  
  #-------------------------------
  # Stable proxy: never returns NULL after first valid GWAS
  # -----------------------------
  gwas_last_good <- shiny::reactiveVal(NULL)
  
  gwas_df_proxy <- shiny::reactive({
    df <- NULL
    
    # 1) prefer local gwas_df() si existeix i és bo
    if (exists("gwas_df", inherits = TRUE) && is.function(gwas_df)) {
      df2 <- tryCatch(gwas_df(), error = function(e) NULL)
      if (is.data.frame(df2) && nrow(df2) > 0) df <- df2
    }
    
    # 2) sinó, Hub GWAS
    if (is.null(df)) {
      df3 <- gwas_shared_rv()
      if (is.data.frame(df3) && nrow(df3) > 0) df <- df3
    }
    
    # 3) cache
    if (is.data.frame(df) && nrow(df) > 0) {
      gwas_last_good(df)
      return(df)
    }
    
    # 4) fallback: últim GWAS bo (evita flicker)
    gwas_last_good()
  })
  
  # -----------------------------
  source(file.path("R", "gi_slave_canonical_gtex.R"), local = TRUE)
  gi_sync <- gi_slave_canonical_init(session)
  
  observe({
    cat("[DBG] gi_sync names:", paste(names(gi_sync), collapse=", "), "\n")
    cat("[DBG] has gwas_shared:", "gwas_shared" %in% names(gi_sync), "\n")
    cat("[DBG] has clusters_shared:", "clusters_shared" %in% names(gi_sync), "\n")
  })
  
  observe({
    sid <- gi_sid(session)
    p   <- gi_state_paths(sid)
    cat("[GTEX][SYNC] sid=", sid, "\n")
    cat("[GTEX][SYNC] json=", p$json, " exists=", file.exists(p$json), "\n")
    cat("[GTEX][SYNC] gwas=", p$gwas, " exists=", file.exists(p$gwas), "\n")
    cat("[GTEX][SYNC] clus=", p$clus, " exists=", file.exists(p$clus), "\n")
    cat("[GTEX][SYNC] params=", p$params, " exists=", file.exists(p$params), "\n")
  })
  # ------------------------------------------------------------
  # HUB sync -> GWAS input (fills gwas_shared_rv)
  # ------------------------------------------------------------
  observeEvent(gi_sync$gwas_shared(), {
    
    df <- gi_sync$gwas_shared()
    req(is.data.frame(df), nrow(df) > 0)
    
    gwas_shared_rv(df)
    
    cat("[GTEx][HUB] GWAS synced rows=", nrow(df), " cols=", ncol(df), "\n")
    showNotification(sprintf("✅ HUB: GWAS sincronitzat (%d files).", nrow(df)),
                     type = "message", duration = 2)
    
  }, ignoreInit = FALSE)
  
  # -----------------------------
  # Step 1: GWAS reader + p selector (stable)
  # - Prioritza master (gwas_shared_rv)
  # - Fallback upload
  # -----------------------------
  
  dat_raw <- reactive({
    df_shared <- gwas_shared_rv()
    if (is.data.frame(df_shared) && nrow(df_shared) > 0) {
      return(df_shared)
    }
    
    req(input$gwas_file)
    readr::read_delim(
      file   = input$gwas_file$datapath,
      delim  = input$gwas_sep %||% "\t",
      col_names = isTRUE(input$gwas_header),
      locale = readr::locale(decimal_mark = ".", grouping_mark = ","),
      show_col_types = FALSE,
      na = c("", "NA")
    )
  })
  
  output$gwas_p_selector <- renderUI({
    df <- dat_raw()
    validate(need(is.data.frame(df) && nrow(df) > 0, "Could not read GWAS."))
    
    cols <- names(df)
    
    guess <- c(
      grep("^p$|pval|pvalue|p_value|p\\.value", cols, ignore.case = TRUE, value = TRUE),
      intersect(cols, c("Pval","PVAL","pval","pvalue","P","p","P_VALUE","p_value"))
    )
    guess <- guess[!is.na(guess) & nzchar(guess)]
    guess1 <- if (length(guess)) guess[1] else cols[1]
    
    selectInput(
      "gwas_col_p",
      "Select P-value column:",
      choices  = cols,
      selected = guess1
    )
  })
  
  gwas_df <- reactive({
    df_shared <- gwas_shared_rv()
    if (is.data.frame(df_shared) && nrow(df_shared) > 0) return(df_shared)
    
    df <- dat_raw()
    validate(need(is.data.frame(df) && nrow(df) > 0, "Empty GWAS table."))
    
    validate(need(!is.null(input$gwas_col_p) && nzchar(input$gwas_col_p), "Select a P-value column."))
    validate(need(input$gwas_col_p %in% names(df), "Selected P column not found in GWAS."))
    
    # 1) agafa rawP ABANS de tolower, amb el nom original seleccionat
    rawP <- df[[ input$gwas_col_p ]]
    validate(
      need(!is.null(rawP), "Selected P column is NULL."),
      need(length(rawP) == nrow(df), "Selected P column length mismatch.")
    )
    
    # 2) Ara sí: normalitza noms per detectar chr/bp/rsid
    df2 <- df
    names(df2) <- tolower(names(df2))
    
    validate(need(all(c("chr","bp") %in% names(df2)), "GWAS must contain columns chr and bp (case-insensitive)."))
    
    if (!"snp" %in% names(df2)) df2$snp <- NA_character_
    
    # rsid (si existeix)
    rscol <- detect_rsid_col(df2)
    if (!is.null(rscol) && rscol %in% names(df2)) {
      df2$rsid <- as.character(df2[[rscol]])
    } else {
      if ("snp" %in% names(df2) && any(grepl("^rs[0-9]+$", df2$snp))) {
        df2$rsid <- as.character(df2$snp)
      } else {
        df2$rsid <- NA_character_
      }
    }
    
    BP   <- if (is.numeric(df2$bp)) df2$bp else suppressWarnings(readr::parse_number(as.character(df2$bp)))
    Pval <- parse_p_robust(rawP)
    CHR  <- chr_map_plink19(df2$chr)
    
    out <- df2 %>%
      dplyr::mutate(
        CHR  = CHR,
        BP   = as.numeric(BP),
        Pval = as.numeric(Pval),
        logp = -log10(Pval)
      ) %>%
      dplyr::filter(is.finite(CHR), is.finite(BP), is.finite(Pval), Pval > 0)
    
    validate(need(nrow(out) > 0, "GWAS parsed but 0 valid rows (check columns/format)."))
    out
  })
  # ------------------------------------------------------------
  # GWAS STABLE (anti-flicker): fixa el primer GWAS vàlid
  # ------------------------------------------------------------
  gwas_df_stable_rv <- shiny::reactiveVal(NULL)
  
  observeEvent(gwas_df(), {
    df <- tryCatch(gwas_df(), error = function(e) NULL)
    if (!is.data.frame(df) || nrow(df) == 0) return()
    if (is.data.frame(gwas_df_stable_rv()) && nrow(gwas_df_stable_rv()) > 0) return()
    gwas_df_stable_rv(df)
    cat("[GTEx] GWAS_STABLE locked rows=", nrow(df), "\n")
  }, ignoreInit = FALSE)
  
  gwas_df_stable <- reactive({
    df <- gwas_df_stable_rv()
    if (is.data.frame(df) && nrow(df) > 0) return(df)
    # fallback mentre encara no s'ha “lockejat”
    df2 <- tryCatch(gwas_df(), error = function(e) NULL)
    if (is.data.frame(df2) && nrow(df2) > 0) return(df2)
    NULL
  })
  
  # ------------------------------------------------------------
  # Canonical clusters engine (local build_ranges)
  # - provides: intervals_raw(), clusters_cur(), selected_cluster()
  # ------------------------------------------------------------
  
  gi_clust <- gi_clusters_canonical_init(
    session = session, input = input, output = output,
    gwas_df = gwas_df_stable,          # <-- AIXÒ
    build_btn_id   = "build_ranges",
    clusters_dt_id = "cluster_dt",
    hits_rows_id   = "hits_tbl_rows_selected",
    app_count_col  = "n_gtex"
  )
  
  intervals_raw    <- gi_clust$intervals_raw
  clusters_cur     <- gi_clust$clusters_cur
  selected_cluster <- gi_clust$selected_cluster
  # ------------------------------------------------------------
  # --- stable latch for clusters (anti-flicker) ---
  clusters_stable_rv <- shiny::reactiveVal(NULL)
  
  observe({
    cl <- tryCatch(clusters_cur(), error = function(e) NULL)
    if (is.data.frame(cl) && nrow(cl) > 0) clusters_stable_rv(cl)
  })
  
  clusters_stable <- shiny::reactive({
    cl <- clusters_stable_rv()
    if (is.data.frame(cl) && nrow(cl) > 0) return(cl)
    # fallback: si encara no hi ha latch
    cl2 <- tryCatch(clusters_cur(), error = function(e) NULL)
    if (is.data.frame(cl2) && nrow(cl2) > 0) return(cl2)
    NULL
  })
  
  # ------------------------------------------------------------
  # CLUSTERS DT (ROBUST): shows errors instead of failing silently
  # ------------------------------------------------------------
  output$cluster_dt <- DT::renderDT({
    cl <- clusters_stable()
    
    if (!is.data.frame(cl) || nrow(cl) == 0) {
      return(DT::datatable(
        data.frame(Message = "No clusters yet (waiting Hub or click: Generate intervals + clusters)"),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    show_cols <- c("cluster_chr_n","cluster_id","chr","start","end","n_snps","top_snp","top_logp","n_gtex")
    keep <- intersect(show_cols, names(cl))
    cl2  <- cl[, keep, drop = FALSE]
    
    DT::datatable(
      cl2, rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom="Bfrtip", 
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength=10, scrollX=TRUE)
    ) %>% DT::formatRound(columns = intersect("top_logp", names(cl2)), digits = 2)
  }, server = FALSE)
  
  # ------------------------------------------------------------
  # ---- Preview text (Step 2) ----
  output$ranges_preview <- renderPrint({
    # cl <- clusters_cur()
    cl <- clusters_stable()
    if (!is.data.frame(cl) || nrow(cl) == 0) {
      cat("No clusters yet. Click: ➊ Generate intervals + clusters\n")
      return(invisible())
    }
    
    cat("Clusters generated:", nrow(cl), "\n\n")
    print(utils::head(cl, 10))
  })
  
  ###########################################
  # ------------------------------------------------------------
  # HUB sync -> apply CLUSTERS into canonical engine (FINAL)
  # - NO polling extra
  # - 1 sola font de veritat: gi_slave_canonical_init()
  # ------------------------------------------------------------
  
  # Helper: aplica clusters + intervals al motor canònic
  .apply_hub_clusters_to_engine <- function(cl) {
    if (!is.data.frame(cl) || !nrow(cl)) return(FALSE)
    cl <- as.data.frame(cl)
    
    # Normalitza camps mínims
    if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- cl$CHR
    
    if (!"start" %in% names(cl)) {
      if ("cluster_start" %in% names(cl)) cl$start <- cl$cluster_start
      if ("start_bp"      %in% names(cl)) cl$start <- cl$start_bp
    }
    if (!"end" %in% names(cl)) {
      if ("cluster_end" %in% names(cl)) cl$end <- cl$cluster_end
      if ("end_bp"      %in% names(cl)) cl$end <- cl$end_bp
    }
    
    if (!all(c("chr","start","end") %in% names(cl))) return(FALSE)
    
    cl$chr   <- suppressWarnings(as.integer(cl$chr))
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    if (!nrow(cl)) return(FALSE)
    
    # IDs canònics (necessaris per DT + highlights)
    if (!"cluster_chr_n" %in% names(cl) || all(!nzchar(as.character(cl$cluster_chr_n)))) {
      cl$cluster_chr_n <- paste0("chr", chr_label_plink(cl$chr), "_", ave(cl$chr, cl$chr, FUN = seq_along))
    } else {
      cl$cluster_chr_n <- sub("^cluster_", "", as.character(cl$cluster_chr_n))
    }
    
    if (!"cluster_id" %in% names(cl) || all(!nzchar(as.character(cl$cluster_id)))) {
      cl$cluster_id <- cl$cluster_chr_n
    } else {
      cl$cluster_id <- sub("^cluster_", "", as.character(cl$cluster_id))
    }
    
    if (!"cluster" %in% names(cl)) cl$cluster <- seq_len(nrow(cl))
    
    # intervals canònics per a bandes/plot inferior
    iv <- dplyr::transmute(
      cl,
      chr   = .data$chr,
      start = .data$start,
      end   = .data$end,
      label = .data$cluster_id
    )
    
    ok1 <- FALSE; ok2 <- FALSE
    
    # 1) Preferit: setters del motor canònic
    if (is.list(gi_clust) && !is.null(gi_clust$set_clusters) && is.function(gi_clust$set_clusters)) {
      gi_clust$set_clusters(cl); ok1 <- TRUE
    } else {
      # 2) Fallback: reactiveVal exposada (clusters_cur)
      try({ gi_clust$clusters_cur(cl); ok1 <- TRUE }, silent = TRUE)
    }
    
    if (is.list(gi_clust) && !is.null(gi_clust$set_intervals) && is.function(gi_clust$set_intervals)) {
      gi_clust$set_intervals(iv); ok2 <- TRUE
    } else {
      # 2) Fallback: reactiveVal exposada (intervals_raw)
      try({ gi_clust$intervals_raw(iv); ok2 <- TRUE }, silent = TRUE)
    }
    
    cat("[GTEx][HUB] applied CLUSTERS rows=", nrow(cl),
        " set_clusters=", ok1, " set_intervals=", ok2, "\n")
    
    ok1 && ok2
  }
  
  # IMPORTANT:
  # - NO fem polling manual de fitxers _state aquí
  # - ens fiem de gi_sync$clusters_shared() (ja fa polling del JSON del Hub)
  observeEvent(gi_sync$clusters_shared(), {
    cl <- gi_sync$clusters_shared()
    req(is.data.frame(cl), nrow(cl) > 0)
    
    .apply_hub_clusters_to_engine(cl)
    
    showNotification(sprintf("✅ HUB: CLUSTERS sincronitzats (%d).", nrow(cl)),
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
  
  
  # carpeta per aquesta sessió
  workdir <- file.path(tempdir(), paste0("gtex_inspector_", Sys.getpid()))
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  session$onSessionEnded(function(...) {
    if (dir.exists(workdir)) unlink(workdir, recursive = TRUE, force = TRUE)
  })
  
  # === Tissue family definitions (GTEx v10 robust) ===
  tissue_families <- list(
    Brain = c(
      "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
      "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex",
      "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus",
      "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
      "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra"
    ),
    Cardiovascular = c(
      "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
      "Heart_Atrial_Appendage", "Heart_Left_Ventricle"
    ),
    Immune = c("Whole_Blood", "Spleen"),
    Digestive = c(
      "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa",
      "Esophagus_Muscularis", "Stomach", "Colon_Sigmoid",
      "Colon_Transverse", "Small_Intestine_Terminal_Ileum", "Liver"
    ),
    Endocrine = c("Pancreas", "Pituitary", "Thyroid", "Adrenal_Gland"),
    Adipose = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum"),
    Muscle = c("Muscle_Skeletal"),
    Skin = c("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"),
    Lung = c("Lung"),
    Kidney = c("Kidney_Cortex"),
    Others = NULL
  )
  
  map_tissue_family <- function(tissue) {
    for (fam in names(tissue_families)) {
      tis_list <- tissue_families[[fam]]
      if (!is.null(tis_list) && tissue %in% tis_list) return(fam)
    }
    return("Others")
  }
  
  # ------------------------------------------------------------
  # GWAS for plots (HUB-first, robust)
  # ------------------------------------------------------------
  gwas_for_plot <- reactive({
    # 1) HUB (canonical) si existeix
    df <- NULL
    if (exists("gi_sync", inherits = TRUE) && is.list(gi_sync) &&
        "gwas_shared" %in% names(gi_sync) && is.function(gi_sync$gwas_shared)) {
      df <- tryCatch(gi_sync$gwas_shared(), error = function(e) NULL)
    }
    
    # 2) App-level cached GWAS (si el tens)
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      df <- tryCatch(gwas_shared_rv(), error = function(e) NULL)
    }
    
    # 3) Últim recurs: el teu gwas_df_stable() (si existeix)
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      if (exists("gwas_df_stable", inherits = TRUE) && is.function(gwas_df_stable)) {
        df <- tryCatch(gwas_df_stable(), error = function(e) NULL)
      }
    }
    
    validate(need(is.data.frame(df) && nrow(df) > 0,
                  "GWAS not available yet (Hub GWAS not loaded)."))
    
    # Normalitza tipus mínims
    df <- as.data.frame(df)
    validate(need(all(c("CHR","BP","Pval") %in% names(df)),
                  "GWAS must contain CHR, BP, Pval (from master)."))
    
    df$CHR  <- suppressWarnings(as.integer(df$CHR))
    df$BP   <- suppressWarnings(as.numeric(df$BP))
    df$Pval <- suppressWarnings(as.numeric(df$Pval))
    if (!"logp" %in% names(df)) df$logp <- -log10(df$Pval)
    if (!"snp"  %in% names(df)) df$snp  <- NA_character_
    
    df <- df[is.finite(df$CHR) & is.finite(df$BP) & is.finite(df$Pval) & df$Pval > 0, , drop = FALSE]
    validate(need(nrow(df) > 0, "GWAS arrived but became empty after coercion."))
    
    df
  })
  
  dfp_manhattan <- reactive({
    df <- gwas_for_plot()
    ref <- .ref_hg38
    
    # IMPORTANT: no filtris aquí si vols veure “sempre alguna cosa”
    # (si filtres a Pval < 0.05 i el master et passa un subset estrany, et quedes a 0)
    # Si vols filtrar, fes-ho amb validate fallback:
    df2 <- df
    # df2 <- df2 %>% dplyr::filter(Pval < 0.05)
    
    df2$chrN <- norm_chr_generic(df2$CHR)
    
    out <- df2 %>%
      dplyr::inner_join(ref %>% dplyr::select(chr, chr_cum), by = c("chrN" = "chr")) %>%
      dplyr::arrange(chrN, BP) %>%
      dplyr::mutate(BPcum = BP + chr_cum)
    
    validate(need(nrow(out) > 0, "GWAS available but dfp_manhattan became empty (chr join/filter)."))
    out
  })
  
  axis_df <- reactive({
    dfp <- dfp_manhattan()
    dfp %>%
      group_by(chrN) %>%
      summarise(center = mean(BPcum, na.rm = TRUE), .groups = "drop")
  })
  
  # ---------- Hits (GWAS significant hits for current mode) ----------
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
        p = Pval,
        logp
      )
  })
  
  # temporal
  observe({
    cluster_method <- as.character(input$cluster_method %||% "")
    hits_mode      <- as.character(input$hits_mode %||% "")
    
    thr <- if (identical(cluster_method, "window")) {
      suppressWarnings(as.numeric(input$pthr))
    } else if (identical(cluster_method, "hits")) {
      suppressWarnings(as.numeric(input$min_logp))
    } else {
      NA_real_
    }
    
    cat(
      "[GTEX HITS_DF][DEBUG] cluster_method=", cluster_method,
      " | hits_mode=", hits_mode,
      " | threshold=", thr, "\n",
      sep = ""
    )
    
    df <- try(hits_df(), silent = TRUE)
    
    if (!inherits(df, "try-error") && is.data.frame(df)) {
      cat("[GTEX HITS_DF][DEBUG] nrow=", nrow(df), " ncol=", ncol(df), "\n", sep = "")
      
      if (nrow(df) > 0 && "CHR" %in% names(df)) {
        chr_counts <- df %>%
          dplyr::count(CHR, name = "n_hits") %>%
          dplyr::arrange(suppressWarnings(as.integer(CHR)))
        
        cat("[GTEX HITS_DF][DEBUG] counts by chromosome:\n")
        print(chr_counts)
      }
    }
  })
  # end
  
  
  hits_enriched <- reactive({
    df <- hits_df()
    if (is.null(df) || !nrow(df)) return(df)
    
    if ("p" %in% names(df)) df$p <- formatC(df$p, format = "e", digits = 2)
    if ("logp" %in% names(df)) df$logp <- sprintf("%.2f", df$logp)
    
    rscol <- "rsid"
    df$dbSNP <- ifelse(
      !is.na(df[[rscol]]) & df[[rscol]] != "",
      paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/", df[[rscol]], "' target='_blank'>", df[[rscol]], "</a>"),
      NA
    )
    df
  })
  
  output$hits_tbl <- renderDT({
    df <- hits_enriched()
    if (is.null(df) || !nrow(df)) {
      return(datatable(data.frame(Message="No GWAS hits above threshold yet."), options=list(dom="t"), rownames=FALSE))
    }
    datatable(
      df %>% select(CHR, BP, snp, rsid, dbSNP, p, logp),
      selection = "multiple",
      rownames = FALSE,
      escape = FALSE,
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
  
  # ---------- Step 3: carregar GTEx (www/gtex_eqtl_light.rds) ----------
  gtex_all <- reactiveVal(NULL)
  
  observe({
    path <- "www/gtex_eqtl_light.rds"
    if (file.exists(path)) {
      df <- tryCatch(readRDS(path), error = function(e) NULL)
      if (!is.null(df) && is.data.frame(df)) {
        
        # Normalització mínima robusta (clau per evitar “taula buida” per mismatch chr)
        if ("chr" %in% names(df)) df$chr <- norm_chr_generic(df$chr)
        if ("pos" %in% names(df)) df$pos <- suppressWarnings(as.integer(df$pos))
        
        gtex_all(df)
        output$gtex_status <- renderText(
          sprintf("GTEx loaded: %s (rows: %d, cols: %d)",
                  path, nrow(df), ncol(df))
        )
      } else {
        output$gtex_status <- renderText("Could not read www/gtex_eqtl_light.rds.")
      }
    } else {
      output$gtex_status <- renderText("GTEx file not found: www/gtex_eqtl_light.rds")
    }
  })
  
  # ===========================
  # Step 4: extreure GTEx dins clusters + assignar cluster/cluster_chr_n a cada hit
  # ===========================
  
  gtex_hits_val <- reactiveVal(NULL)
  
  extract_gtex_from_clusters <- function(clusters_df, gtex_df) {
    if (is.null(clusters_df) || !is.data.frame(clusters_df) || !nrow(clusters_df)) {
      return(gtex_df[0, , drop = FALSE])
    }
    if (is.null(gtex_df) || !is.data.frame(gtex_df) || !nrow(gtex_df)) {
      return(gtex_df[0, , drop = FALSE])
    }
    
    shiny::validate(
      shiny::need(all(c("chr","pos") %in% colnames(gtex_df)),
                  "GTEx table must contain columns 'chr' and 'pos'.")
    )
    
    cl <- as.data.frame(clusters_df)
    
    # assegura chr/start/end
    if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end" %in% names(cl) && "end_bp" %in% names(cl)) cl$end <- cl$end_bp
    
    cl$chr   <- suppressWarnings(as.integer(cl$chr))
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    
    # cluster_chr_n / cluster_id (string) + normalitza prefix "cluster_"
    if (!"cluster_chr_n" %in% names(cl)) {
      if ("cluster_id" %in% names(cl)) cl$cluster_chr_n <- as.character(cl$cluster_id)
      else cl$cluster_chr_n <- NA_character_
    }
    if (!"cluster_id" %in% names(cl)) {
      cl$cluster_id <- cl$cluster_chr_n
    }
    
    cl$cluster_id    <- sub("^cluster_", "", as.character(cl$cluster_id))
    cl$cluster_chr_n <- sub("^cluster_", "", as.character(cl$cluster_chr_n))
    
    # cluster (num) si no existeix
    if (!"cluster" %in% names(cl)) {
      cl$cluster <- seq_len(nrow(cl))
    } else {
      cl$cluster <- suppressWarnings(as.integer(cl$cluster))
      if (any(!is.finite(cl$cluster))) {
        # si ve trencat, força seq
        cl$cluster <- seq_len(nrow(cl))
      }
    }
    
    # etiqueta chr per comparar amb gtex_df$chr (que és "1","2","X"? o "chr1"?)
    # tu estàs fent chr_label_plink() i després compares amb gtex_df$chr
    # deixem-ho igual però robust
    cl$chr_i <- vapply(cl$chr, function(x) chr_label_plink(as.integer(x)), character(1))
    
    out <- vector("list", length = 0)
    
    for (i in seq_len(nrow(cl))) {
      chr_i   <- cl$chr_i[i]
      start_i <- cl$start[i]
      end_i   <- cl$end[i]
      
      if (!is.finite(start_i) || !is.finite(end_i) || end_i < start_i) next
      if (is.na(chr_i) || !nzchar(chr_i)) next
      
      sub <- gtex_df %>%
        dplyr::filter(.data$chr == chr_i, .data$pos >= start_i, .data$pos <= end_i)
      
      if (nrow(sub) > 0) {
        # assignacions sempre de longitud 1
        sub$cluster       <- cl$cluster[i]
        sub$cluster_chr_n <- cl$cluster_chr_n[i] %||% cl$cluster_id[i] %||% paste0("chr", chr_i, "_", i)
        sub$cluster_start <- start_i
        sub$cluster_end   <- end_i
        out[[length(out) + 1]] <- sub
      }
    }
    
    if (!length(out)) return(gtex_df[0, , drop = FALSE])
    dplyr::bind_rows(out)
  }
  
  
 ######
  observeEvent(input$run_gtex, {
    
    tryCatch({
      
      req(gtex_all())
      
      # cl <- clusters_cur()
      cl <- clusters_stable()
      validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters. Generate intervals + clusters first."))
      
      # ------------------------------------------------------------
      # 1) Normalize GTEx chr/pos so extract works correctly
      # ------------------------------------------------------------
      gtex_norm <- gtex_all() %>%
        dplyr::mutate(
          chr = as.character(chr),
          chr = gsub("^CHR", "chr", chr),
          chr = gsub("^chr", "", chr),
          chr = toupper(chr),
          chr = dplyr::if_else(chr == "23", "X", chr),
          chr = dplyr::if_else(chr == "24", "Y", chr),
          pos = suppressWarnings(as.integer(pos))
        ) %>%
        dplyr::filter(!is.na(chr), !is.na(pos))
      
      cat("\n====================\n")
      cat("[gtex_norm]  nrow=", nrow(gtex_norm), " ncol=", ncol(gtex_norm), "\n", sep = "")
      cat("[gtex_norm]  cols: ", paste(names(gtex_norm), collapse = ", "), "\n", sep = "")
      cat("====================\n")
      if (nrow(gtex_norm) > 0) print(utils::head(gtex_norm, 5))
      
      # ------------------------------------------------------------
      # 2) Extract hits from clusters
      # ------------------------------------------------------------
      hits <- tryCatch(
        extract_gtex_from_clusters(cl, gtex_norm),
        error = function(e) {
          cat("[GTEx][run_gtex] extract_gtex_from_clusters ERROR: ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      
      if (is.null(hits) || !is.data.frame(hits)) {
        hits <- gtex_norm[0, , drop = FALSE]
      }
      
      cat("\n====================\n")
      cat("[gtex_hits_df]  nrow=", nrow(hits), " ncol=", ncol(hits), "\n", sep = "")
      cat("[gtex_hits_df]  cols: ", paste(names(hits), collapse = ", "), "\n", sep = "")
      cat("====================\n")
      if (nrow(hits) > 0) print(utils::head(hits, 5))
      
      gtex_hits_val(hits)
      
      # ------------------------------------------------------------
      # 3) Update n_gtex robustly
      # ------------------------------------------------------------
      cl_base <- cl %>% dplyr::select(-dplyr::any_of("n_gtex"))
      
      if (nrow(hits) > 0) {
        
        if ("cluster_chr_n" %in% names(hits)) {
          
          hits2 <- hits %>%
            dplyr::mutate(
              cluster_chr_n = as.character(cluster_chr_n),
              cluster_chr_n = sub("^cluster_", "", cluster_chr_n)
            )
          
          cl_key <- cl_base %>%
            dplyr::mutate(
              cluster_chr_n = dplyr::coalesce(
                as.character(.data$cluster_chr_n),
                as.character(.data$cluster_id)
              ),
              cluster_chr_n = sub("^cluster_", "", cluster_chr_n)
            )
          
          cnt <- hits2 %>%
            dplyr::count(cluster_chr_n, name = "n_gtex") %>%
            dplyr::mutate(n_gtex = as.integer(n_gtex))
          
          cl2 <- cl_key %>%
            dplyr::left_join(cnt, by = "cluster_chr_n") %>%
            dplyr::mutate(n_gtex = dplyr::coalesce(.data$n_gtex, 0L))
          
          clusters_cur(cl2)
          
        } else if ("cluster" %in% names(hits) && "cluster" %in% names(cl_base)) {
          
          cnt <- hits %>%
            dplyr::count(cluster, name = "n_gtex") %>%
            dplyr::mutate(n_gtex = as.integer(n_gtex))
          
          cl2 <- cl_base %>%
            dplyr::left_join(cnt, by = "cluster") %>%
            dplyr::mutate(n_gtex = dplyr::coalesce(.data$n_gtex, 0L))
          
          clusters_cur(cl2)
          
        } else {
          
          cl2 <- cl_base %>% dplyr::mutate(n_gtex = 0L)
          clusters_cur(cl2)
        }
        
      } else {
        cl2 <- cl_base %>% dplyr::mutate(n_gtex = 0L)
        clusters_cur(cl2)
      }
      
      ###### Start INTEGRATOR CODE ####################################################################
      tryCatch({
        
        app_slug          <- "gtex"
        app_hit_class     <- "gtex_hit"
        evidence_gene     <- "eqtl_gene"
        
        append_log <- if (exists("append_log", mode = "function", inherits = TRUE)) {
          get("append_log", mode = "function", inherits = TRUE)
        } else {
          function(...) cat(..., "\n", sep = "")
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
        
        get_gtex_cluster_method <- function(input) {
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
        
        get_gtex_threshold <- function(input, cluster_method) {
          if (identical(cluster_method, "window")) {
            suppressWarnings(as.numeric(input$pthr))
          } else if (cluster_method %in% c("hits_span1mb", "hits_tiled", "hits_sliding", "hits_unknown")) {
            suppressWarnings(as.numeric(input$min_logp))
          } else {
            NA_real_
          }
        }
        
        get_gtex_gwas_name <- function(input) {
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
        current_cluster_method <- get_gtex_cluster_method(input)
        current_threshold      <- get_gtex_threshold(input, current_cluster_method)
        current_gwas_name      <- get_gtex_gwas_name(input)
        
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
        
        if (!length(manifest_gwas_files)) {
          manifest$gwas_session_file <- current_gwas_name_chr
        } else if (nzchar(current_gwas_name_chr) && !(current_gwas_name_chr %in% manifest_gwas_files)) {
          manifest$gwas_session_file <- unique(c(manifest_gwas_files, current_gwas_name_chr))
        } else {
          manifest$gwas_session_file <- manifest_gwas_files
        }
        
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
        needed_gene_cols <- c("gene_name", "cluster_chr_n", "chr", "cluster_start", "cluster_end")
        has_gene_id_col <- "gene_id" %in% names(hits)
        
        gene_bridge <- if (all(needed_gene_cols %in% names(hits))) {
          
          gb0 <- hits %>%
            dplyr::transmute(
              gene = as.character(gene_name),
              gtex_gene_name = as.character(gene_name),
              source_app = app_slug,
              evidence_type = evidence_gene,
              cluster_id = as.character(cluster_chr_n),
              chr = chr,
              start = cluster_start,
              end = cluster_end
            )
          
          if (has_gene_id_col) {
            gb0 <- gb0 %>%
              dplyr::mutate(
                gtex_gene_id = as.character(hits$gene_id)
              )
          } else {
            gb0 <- gb0 %>%
              dplyr::mutate(
                gtex_gene_id = NA_character_
              )
          }
          
          gb0 %>%
            sanitize_bridge() %>%
            sanitize_bridge_genes() %>%
            dplyr::mutate(
              gtex_gene_name = dplyr::coalesce(as.character(gtex_gene_name), ""),
              gtex_gene_id = dplyr::coalesce(as.character(gtex_gene_id), "")
            ) %>%
            dplyr::filter(!is.na(gene), nzchar(gene)) %>%
            dplyr::distinct()
          
        } else {
          append_log(paste0(
            log_tag, " gene bridge skipped: missing columns -> ",
            paste(setdiff(needed_gene_cols, names(hits)), collapse = ", ")
          ))
          tibble::tibble()
        }
        
        append_log(paste0(log_tag, " gene_bridge_path=", gene_bridge_path))
        append_log(paste0(log_tag, " gene_bridge nrow=", nrow(gene_bridge), " ncol=", ncol(gene_bridge)))
        
        if (nrow(gene_bridge) > 0) {
          saveRDS(gene_bridge, gene_bridge_path)
          manifest$files[[paste0(app_slug, "_gene_bridge")]] <- basename(gene_bridge_path)
          append_log(paste0(
            log_tag, " saved ", basename(gene_bridge_path),
            " | exists=", file.exists(gene_bridge_path)
          ))
        } else {
          append_log(paste0(
            log_tag, " ", basename(gene_bridge_path),
            " NOT saved (0 rows)"
          ))
        }
        
        # ============================================================
        # B) TERM BRIDGE (empty for GTEx)
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
        
        append_log(paste0(log_tag, " term_bridge_path=", term_bridge_path))
        saveRDS(term_bridge, term_bridge_path)
        manifest$files[[paste0(app_slug, "_term_bridge")]] <- basename(term_bridge_path)
        append_log(paste0(
          log_tag, " saved empty ", basename(term_bridge_path),
          " | exists=", file.exists(term_bridge_path)
        ))
        
        # ============================================================
        # C) LD INPUTS
        # ============================================================
        append_log(paste0(ld_tag, " dir=", session_dir))
        append_log(paste0(ld_tag, " dir exists: ", dir.exists(session_dir)))
        
        # ----------------------------------------------------------
        # C1) CLUSTERS MASTER
        # ----------------------------------------------------------
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
          
          append_log(paste0(ld_tag, "[CLUSTERS] clusters_master nrow=", nrow(clusters_master), " ncol=", ncol(clusters_master)))
          
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
        
        id_hit_col_candidates <- c("rsid", "variant_id", "SNP", "snp", "id", "ID")
        id_hit_col <- intersect(id_hit_col_candidates, names(hits))
        id_hit_col <- if (length(id_hit_col)) id_hit_col[1] else NULL
        
        needed_hit_cols <- c("cluster_chr_n", "chr", "pos")
        
        gtex_candidates <- tibble::tibble()
        
        if (all(needed_hit_cols %in% names(hits)) && !is.null(id_hit_col)) {
          gtex_candidates <- hits %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_chr_n),
              chr        = suppressWarnings(as.integer(readr::parse_number(as.character(chr)))),
              pos_ini    = suppressWarnings(as.integer(readr::parse_number(as.character(pos)))),
              pos_end    = suppressWarnings(as.integer(readr::parse_number(as.character(pos)))),
              id_hit     = as.character(.data[[id_hit_col]]),
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
            paste(
              c(
                setdiff(needed_hit_cols, names(hits)),
                if (is.null(id_hit_col)) "id_hit_col" else NULL
              ),
              collapse = ", "
            )
          ))
        }
        
        append_log(paste0(ld_tag, "[CANDIDATES][APP] id_hit_col=", id_hit_col %||% "NULL"))
        append_log(paste0(ld_tag, "[CANDIDATES][APP] nrow=", nrow(gtex_candidates)))
        if (nrow(gtex_candidates) > 0) {
          append_log(
            paste0(
              ld_tag, "[CANDIDATES][APP] head:\n",
              paste(capture.output(print(utils::head(gtex_candidates, 10))), collapse = "\n")
            )
          )
        }
        
        candidates_ld <- dplyr::bind_rows(gwas_candidates, gtex_candidates) %>%
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
        
    # }, error = function(e) {
    #   if (exists("append_log", mode = "function", inherits = TRUE)) {
    #     append_log(paste0("[LD-INTEGRATOR][ERROR] ", conditionMessage(e)))
    #   } else {
    #     cat("[LD-INTEGRATOR][ERROR] ", conditionMessage(e), "\n", sep = "")
    #   }
    # })
        
      }, error = function(e) {
        cat("[LD-INTEGRATOR][ERROR] ", conditionMessage(e), "\n", sep = "")
      })
      
      ####################### end integrator bridges #################################
      
      
      incProgress(0.05, detail = "Done")
      cat("[END] run_gtex finished OK\n")
      
      if (is.null(hits) || !nrow(hits)) {
        showNotification("No GTEx eQTL hits found in clusters.", type = "warning", duration = 6)
      } else {
        showNotification("GTEx eQTL hits assigned to clusters and final outputs saved.", type = "message", duration = 4)
      }
      
    }, error = function(e) {
      
      gtex_final_path_csv(NULL)
      gtex_final_path_rds(NULL)
      gtex_norm_path_csv(NULL)
      gtex_norm_path_rds(NULL)
      
      msg <- paste0("ERROR (run_gtex): ", conditionMessage(e))
      
      output$gtex_summary <- renderText(msg)
      cat("[ERROR] ", msg, "\n", sep = "")
      showNotification(conditionMessage(e), type = "error", duration = 12)
      
    }, finally = {
      
      if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::enable("run_gtex")
      
    })
    
  }, ignoreInit = TRUE)
  ###  
  gtex_data <- reactive({
    dt <- gtex_hits_val()
    if (!is.data.frame(dt) || !nrow(dt)) return(dt)
    
    if (!"tissue" %in% names(dt)) dt$tissue <- NA_character_
    dt$tissue <- trimws(as.character(dt$tissue))
    
    if (!"gene_id" %in% names(dt)) dt$gene_id <- NA_character_
    dt$gene_id <- as.character(dt$gene_id)
    
    dt
  })
  
  # ===========================
  # GTEx table (cluster-aware)
  # ===========================
  
  gtex_table_filtered <- reactive({
    dt <- gtex_data()
    if (is.null(dt) || !nrow(dt)) return(dt)
    
    cl <- selected_cluster()
    if (is.null(cl) || !nrow(cl)) return(dt)
    
    filter(dt, cluster == cl$cluster[1])
  })
  
  gtex_table_all <- reactive({
    dt <- gtex_data()
    if (is.null(dt) || !nrow(dt)) return(dt)
    dt
  })
  
  output$gtex_table <- renderDT({
    
    dt <- gtex_table_filtered()
    
    if (is.null(dt) || !nrow(dt)) {
      return(datatable(
        data.frame(Message = "Run Step 4 to extract GTEx eQTLs."),
        options = list(dom = "t")
      ))
    }
    
    # --- assegura columnes per evitar errors quan faltin ---
    if (!"gene_name"  %in% names(dt)) dt$gene_name <- NA_character_
    if (!"gene_id"    %in% names(dt)) dt$gene_id   <- NA_character_
    if (!"variant_id" %in% names(dt)) dt$variant_id <- NA_character_
    if (!"ref"        %in% names(dt)) dt$ref <- NA_character_
    if (!"alt"        %in% names(dt)) dt$alt <- NA_character_
    if (!"pval"       %in% names(dt)) {
      if ("pval_nominal" %in% names(dt)) dt$pval <- dt$pval_nominal
      else if ("pvalue"  %in% names(dt)) dt$pval <- dt$pvalue
      else dt$pval <- NA_real_
    }
    if (!"slope" %in% names(dt)) dt$slope <- NA_real_
    
    dt <- dt %>%
      mutate(
        chr_clean = gsub("^chr", "", chr, ignore.case = TRUE),
        
        # si variant_id ja ve ple (GTEx), fem-lo servir; si no, el construïm
        gtex_variant = dplyr::if_else(
          !is.na(variant_id) & variant_id != "",
          variant_id,
          paste0("chr", chr_clean, "_", pos, "_", ref, "_", alt, "_b38")
        ),
        
        gtex_var_url = paste0("https://www.gtexportal.org/home/snp/", gtex_variant),
        
        # IMPORTANT: condició vectorial (NO &&)
        gtex_gene_url = dplyr::if_else(
          !is.na(gene_name) & gene_name != "",
          paste0("https://www.gtexportal.org/home/locusBrowserPage/", gene_name),
          NA_character_
        ),
        
        gtex_geneid_url = dplyr::if_else(
          !is.na(gene_id) & gene_id != "",
          paste0("https://www.gtexportal.org/home/gene/", gene_id),
          NA_character_
        ),
        
        Gene = dplyr::if_else(
          !is.na(gene_name) & gene_name != "",
          sprintf("<a href='%s' target='_blank'>%s</a>", gtex_gene_url, gene_name),
          gene_name
        ),
        
        Gene_ID = dplyr::if_else(
          !is.na(gene_id) & gene_id != "",
          sprintf("<a href='%s' target='_blank'>%s</a>", gtex_geneid_url, gene_id),
          gene_id
        ),
        
        Variant = sprintf("<a href='%s' target='_blank'>%s</a>", gtex_var_url, gtex_variant),
        
        pval_fmt  = dplyr::if_else(is.na(pval), NA_character_, formatC(pval, format = "e", digits = 2)),
        slope_fmt = dplyr::if_else(is.na(slope), NA_character_, sprintf("%.3f", slope))
      )
    
    datatable(
      dt %>% select(chr, pos, tissue, Gene, Gene_ID, Variant,
                    pval = pval_fmt, slope = slope_fmt, variant_id),
      escape = FALSE,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom="Bfrtip", 
        buttons = list(
          list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
          list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
        ),
        pageLength=10, scrollX=TRUE)
    )
  }, server = FALSE)
  
  # ===========================
  # Visualitzacions GTEx (cluster-aware)
  # ===========================
  
  gtex_vis_df <- reactive({
    dt <- gtex_table_filtered()
    if (is.null(dt) || !nrow(dt)) return(NULL)
    
    # yval
    if ("pval" %in% names(dt)) {
      dt$yval <- -log10(dt$pval)
    } else if ("pval_nominal" %in% names(dt)) {
      dt$yval <- -log10(dt$pval_nominal)
    } else if ("pvalue" %in% names(dt)) {
      dt$yval <- -log10(dt$pvalue)
    } else {
      num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
      num_cols <- setdiff(num_cols, c("pos","BPcum","chr_cum","cluster","cluster_start","cluster_end"))
      first <- num_cols[1]
      dt$yval <- if (!is.null(first)) dt[[first]] else NA_real_
    }
    
    # slope sign
    if ("slope" %in% names(dt)) {
      dt$slope_sign <- case_when(
        dt$slope > 0  ~ "positive",
        dt$slope < 0  ~ "negative",
        TRUE          ~ "zero"
      )
    } else {
      dt$slope_sign <- "unknown"
    }
    
    # tissue family
    if ("tissue" %in% names(dt)) {
      dt$tissue_family <- vapply(dt$tissue, map_tissue_family, FUN.VALUE = character(1))
    } else {
      dt$tissue_family <- "Unknown"
    }
    
    dt
  })
  
  output$tissue_family_bar <- plotly::renderPlotly({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"tissue_family" %in% names(dt)) return(NULL)
    
    df_bar <- dt %>%
      dplyr::mutate(
        slope_cat = dplyr::case_when(
          slope_sign == "positive" ~ "Positive slope",
          slope_sign == "negative" ~ "Negative slope",
          TRUE                     ~ "Other / unknown"
        )
      ) %>%
      dplyr::count(tissue_family, slope_cat, name = "n_hits") %>%
      dplyr::mutate(
        tissue_family = as.character(tissue_family),
        slope_cat     = as.character(slope_cat),
        key_click     = paste(tissue_family, slope_cat, sep = "||"),
        tooltip = paste0(
          "<b>Tissue family:</b> ", htmltools::htmlEscape(tissue_family),
          "<br><b>Slope sign:</b> ", htmltools::htmlEscape(slope_cat),
          "<br><b>eQTL hits:</b> ", n_hits
        )
      )
    
    fam_order <- df_bar %>%
      dplyr::group_by(tissue_family) %>%
      dplyr::summarise(total = sum(n_hits), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::pull(tissue_family)
    
    df_bar$tissue_family <- factor(df_bar$tissue_family, levels = fam_order)
    
    slope_levels <- c("Negative slope", "Positive slope", "Other / unknown")
    df_bar$slope_cat <- factor(df_bar$slope_cat, levels = slope_levels)
    
    cols <- c(
      "Negative slope"  = "#F8766D",
      "Positive slope"  = "#00BFC4",
      "Other / unknown" = "#999999"
    )
    
    plotly::plot_ly(source = "tissue_family_bar") %>%
      plotly::add_bars(
        data = df_bar %>% dplyr::filter(slope_cat == "Negative slope"),
        x = ~tissue_family, y = ~n_hits,
        name = "Negative slope",
        marker = list(color = cols[["Negative slope"]]),
        key = ~key_click,
        hovertext = ~tooltip, hoverinfo = "text",
        textposition = "none"
      ) %>%
      plotly::add_bars(
        data = df_bar %>% dplyr::filter(slope_cat == "Positive slope"),
        x = ~tissue_family, y = ~n_hits,
        name = "Positive slope",
        marker = list(color = cols[["Positive slope"]]),
        key = ~key_click,
        hovertext = ~tooltip, hoverinfo = "text",
        textposition = "none"
      ) %>%
      plotly::add_bars(
        data = df_bar %>% dplyr::filter(slope_cat == "Other / unknown"),
        x = ~tissue_family, y = ~n_hits,
        name = "Other / unknown",
        marker = list(color = cols[["Other / unknown"]]),
        key = ~key_click,
        hovertext = ~tooltip, hoverinfo = "text",
        textposition = "none"
      ) %>%
      plotly::layout(
        barmode = "stack",  # ✅ apilades
        xaxis = list(title = "Tissue family", tickangle = 45, automargin = TRUE),
        yaxis = list(title = "Number of eQTLs", automargin = TRUE),
        margin = list(l = 60, r = 20, t = 30, b = 90),
        font   = list(size = 11),
        legend = list(
          title = list(text = "Slope sign"),
          orientation = "h",
          x = 0, xanchor = "left",
          y = 1.12, yanchor = "top"
        )
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GTEx_tissue_family_bar",
          width = 1600, height = 900, scale = 2
        )
      )
  })
  
  output$tissue_dotplot <- plotly::renderPlotly({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"tissue" %in% names(dt)) return(NULL)
    
    dt2 <- dt %>%
      dplyr::arrange(dplyr::desc(yval)) %>%
      dplyr::slice_head(n = 200) %>%
      dplyr::mutate(
        size_val = if ("slope" %in% names(.)) abs(slope) else 1,
        tooltip = paste0(
          "<b>Tissue:</b> ", htmltools::htmlEscape(tissue),
          if ("tissue_family" %in% names(.)) paste0("<br><b>Family:</b> ", htmltools::htmlEscape(tissue_family)) else "",
          if ("gene_name" %in% names(.)) paste0("<br><b>Gene:</b> ", htmltools::htmlEscape(gene_name)) else "",
          if (all(c("chr","pos") %in% names(.))) paste0("<br><b>Pos:</b> chr", chr, ":", pos) else "",
          "<br><b>-log10(p):</b> ", sprintf("%.2f", yval),
          if ("slope" %in% names(.)) paste0("<br><b>Slope:</b> ", sprintf("%.3g", slope)) else "",
          if ("slope_sign" %in% names(.)) paste0("<br><b>Slope sign:</b> ", htmltools::htmlEscape(slope_sign)) else ""
        )
      )
    
    p <- ggplot2::ggplot(dt2, ggplot2::aes(
      x = tissue, y = yval,
      size = size_val,
      color = slope_sign,
      text = tooltip
    )) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Tissue", y = "-log10(p)", color = "Slope sign") +
      ggplot2::guides(
        size  = "none",
        color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "top",
        legend.justification = "left",
        legend.title = ggplot2::element_text(size = 9),
        legend.text  = ggplot2::element_text(size = 9)
      )
    
    plt <- plotly::ggplotly(p, tooltip = "text")
    
    # --- FIX: remove blank space by forcing all axes domains ---
    ax_names <- names(plt$x$layout)
    x_axes <- grep("^xaxis(\\d+)?$", ax_names, value = TRUE)
    y_axes <- grep("^yaxis(\\d+)?$", ax_names, value = TRUE)
    
    for (ax in x_axes) {
      if (is.null(plt$x$layout[[ax]])) plt$x$layout[[ax]] <- list()
      plt$x$layout[[ax]]$domain <- c(0, 1)
      plt$x$layout[[ax]]$automargin <- TRUE
    }
    for (ax in y_axes) {
      if (is.null(plt$x$layout[[ax]])) plt$x$layout[[ax]] <- list()
      plt$x$layout[[ax]]$domain <- c(0, 1)
      plt$x$layout[[ax]]$automargin <- TRUE
    }
    
    plt %>%
      plotly::layout(
        margin = list(l = 60, r = 20, t = 35, b = 10),
        font   = list(size = 11),
        legend = list(
          orientation = "h",
          x = 0, xanchor = "left",
          y = 1.02, yanchor = "bottom",
          tracegroupgap = 0,
          bgcolor = "rgba(255,255,255,0.7)"
        )
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GTEx_tissue_dotplot",
          width = 1600, height = 900, scale = 2
        )
      )
  })
  
  output$tissue_heatmap <- plotly::renderPlotly({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"gene_name" %in% names(dt)) return(NULL)
    
    # Aggregate: gene x tissue_family -> max -log10(p)
    df_hm <- dt %>%
      dplyr::filter(!is.na(gene_name), gene_name != "", !is.na(tissue_family), tissue_family != "") %>%
      dplyr::group_by(gene_name, tissue_family) %>%
      dplyr::summarise(max_yval = max(yval, na.rm = TRUE), .groups = "drop")
    
    if (!nrow(df_hm)) return(NULL)
    
    # Top genes by max signal
    top_genes <- df_hm %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(max_gene = max(max_yval, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(max_gene)) %>%
      dplyr::slice_head(n = 40) %>%
      dplyr::pull(gene_name)
    
    df_hm <- df_hm %>% dplyr::filter(gene_name %in% top_genes)
    
    # Order genes (top -> bottom) and tissue families (by total hits or alphabetic)
    gene_order <- df_hm %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(max_gene = max(max_yval, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(max_gene) %>%              # small bottom, big top (we'll reverse axis)
      dplyr::pull(gene_name)
    
    fam_order <- df_hm %>%
      dplyr::count(tissue_family, name = "n") %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(tissue_family)
    
    genes <- unique(gene_order)
    fams  <- unique(fam_order)
    
    # Build a complete grid (so missing combos show as NA)
    grid <- tidyr::expand_grid(
      gene_name = genes,
      tissue_family = fams
    ) %>%
      dplyr::left_join(df_hm, by = c("gene_name", "tissue_family"))
    
    # Matrix z (genes x families)
    zmat <- matrix(grid$max_yval, nrow = length(genes), ncol = length(fams), byrow = FALSE)
    
    # Hover text matrix
    txt <- matrix(
      paste0(
        "<b>Gene:</b> ", htmltools::htmlEscape(grid$gene_name),
        "<br><b>Tissue family:</b> ", htmltools::htmlEscape(grid$tissue_family),
        "<br><b>max -log10(p):</b> ", ifelse(is.finite(grid$max_yval), sprintf("%.2f", grid$max_yval), "NA")
      ),
      nrow = length(genes), ncol = length(fams), byrow = FALSE
    )
    
    # Plotly heatmap
    plotly::plot_ly(
      x = fams,
      y = genes,
      z = zmat,
      type = "heatmap",
      text = txt,
      hoverinfo = "text",
      colorscale = "Viridis",
      colorbar = list(
        title = "max -log10(p)",
        orientation = "v",      # vertical (default)
        x = 1.02, xanchor = "left",
        y = 0.5,  yanchor = "middle",
        len = 0.85,
        thickness = 14
      )
    ) %>%
      plotly::layout(
        title = list(text = "Advanced heatmap (gene × tissue family)"),
        xaxis = list(title = "Tissue family", tickangle = 45, automargin = TRUE),
        yaxis = list(title = "Gene", automargin = TRUE, autorange = "reversed"),
        margin = list(l = 140, r = 80, t = 30, b = 90),
        font = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GTEx_gene_tissuefamily_heatmap",
          width = 1600, height = 900, scale = 2
        )
      )
  })
  
  output$tissue_heatmapXXXX <- plotly::renderPlotly({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"gene_name" %in% names(dt)) return(NULL)
    
    df_hm <- dt %>%
      dplyr::group_by(gene_name, tissue_family) %>%
      dplyr::summarise(max_yval = max(yval, na.rm = TRUE), .groups = "drop")
    
    top_genes <- df_hm %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(max_gene = max(max_yval, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(max_gene)) %>%
      dplyr::slice_head(n = 40) %>%
      dplyr::pull(gene_name)
    
    df_hm <- df_hm %>%
      dplyr::filter(gene_name %in% top_genes) %>%
      dplyr::mutate(
        tooltip = paste0(
          "<b>Gene:</b> ", htmltools::htmlEscape(gene_name),
          "<br><b>Tissue family:</b> ", htmltools::htmlEscape(tissue_family),
          "<br><b>max -log10(p):</b> ", sprintf("%.2f", max_yval)
        )
      )
    
    p <- ggplot2::ggplot(df_hm, ggplot2::aes(
      x = tissue_family, y = gene_name, fill = max_yval, text = tooltip
    )) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "C", name = "max -log10(p)") +
      ggplot2::labs(x = "Tissue family", y = "Gene") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 9),
        legend.text  = ggplot2::element_text(size = 9)
      )
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        margin = list(l = 60, r = 20, t = 20, b = 160),
        font   = list(size = 11),
        legend = list(orientation = "h", x = 0, y = -0.25)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GTEx_gene_tissuefamily_heatmap",
          width = 1600, height = 900, scale = 2
        )
      )
  })
  output$tissue_family_barxxxx <- renderPlot({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"tissue_family" %in% names(dt)) return(NULL)
    
    df_bar <- dt %>%
      mutate(
        slope_cat = case_when(
          slope_sign == "positive" ~ "Positive slope",
          slope_sign == "negative" ~ "Negative slope",
          TRUE                     ~ "Other / unknown"
        )
      ) %>%
      count(tissue_family, slope_cat, name = "n_hits")
    
    ggplot(df_bar, aes(x = tissue_family, y = n_hits, fill = slope_cat)) +
      geom_col(position = "stack") +
      labs(x = "Tissue family", y = "Number of eQTLs", fill = "Slope sign") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tissue_dotplotxxxx <- renderPlot({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"tissue" %in% names(dt)) return(NULL)
    
    dt2 <- dt %>% arrange(desc(yval)) %>% slice_head(n = 200)
    
    ggplot(
      dt2,
      aes(
        x = tissue,
        y = yval,
        size = if ("slope" %in% names(dt2)) abs(slope) else 1,
        color = slope_sign
      )
    ) +
      geom_point(alpha = 0.7) +
      coord_flip() +
      labs(x = "Tissue", y = "-log10(p)", size = "|slope|", color = "Slope sign") +
      theme_minimal(base_size = 12)
  })
  
  output$tissue_heatmapxxxx <- renderPlot({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt) || !"gene_name" %in% names(dt)) return(NULL)
    
    df_hm <- dt %>%
      group_by(gene_name, tissue_family) %>%
      summarise(max_yval = max(yval, na.rm = TRUE), .groups = "drop")
    
    top_genes <- df_hm %>%
      group_by(gene_name) %>%
      summarise(max_gene = max(max_yval), .groups = "drop") %>%
      arrange(desc(max_gene)) %>%
      slice_head(n = 40) %>%
      pull(gene_name)
    
    df_hm <- df_hm %>% filter(gene_name %in% top_genes)
    
    ggplot(df_hm, aes(x = tissue_family, y = gene_name, fill = max_yval)) +
      geom_tile() +
      scale_fill_viridis_c(option = "C", name = "max -log10(p)") +
      labs(x = "Tissue family", y = "Gene") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$gtex_interpret <- renderText({
    dt <- gtex_vis_df()
    if (is.null(dt) || !nrow(dt)) {
      return("No GTEx eQTLs available yet. Run Step 4 and/or select a cluster.")
    }
    
    n_hits     <- nrow(dt)
    n_genes    <- length(unique(dt$gene_name %||% NA))
    n_tissues  <- length(unique(dt$tissue %||% NA))
    n_families <- length(unique(dt$tissue_family %||% NA))
    
    fam_summary <- dt %>%
      group_by(tissue_family) %>%
      summarise(n = n(), max_yval = max(yval, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(n)) %>%
      slice_head(n = 3)
    
    fam_text <- if (nrow(fam_summary)) {
      paste(sprintf("- %s: %d hits (max -log10(p) ≈ %.2f)",
                    fam_summary$tissue_family, fam_summary$n, fam_summary$max_yval),
            collapse = "\n")
    } else {
      "No tissue-family information."
    }
    
    slope_tab <- dt %>% count(slope_sign, name = "n") %>% arrange(desc(n))
    slope_text <- if (nrow(slope_tab)) paste(slope_tab$slope_sign, slope_tab$n, collapse = "; ") else "N/A"
    
    top_genes <- dt %>%
      filter(!is.na(gene_name) & gene_name != "") %>%
      group_by(gene_name, chr, pos) %>%
      summarise(max_yval = max(yval, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_yval)) %>%
      slice_head(n = 5)
    
    top5_text <- if (nrow(top_genes)) {
      paste(
        sprintf("%d. %s (chr %s:%d, max -log10(p) ≈ %.2f)",
                seq_len(nrow(top_genes)),
                top_genes$gene_name,
                top_genes$chr,
                top_genes$pos,
                top_genes$max_yval),
        collapse = "\n"
      )
    } else {
      "No gene name available."
    }
    
    paste0(
      "GTEx eQTLs summary (current view)\n",
      "--------------------------------\n",
      "Total eQTL hits: ", n_hits, "\n",
      "Distinct genes:  ", n_genes, "\n",
      "Tissues:         ", n_tissues, "\n",
      "Tissue families: ", n_families, "\n\n",
      "🔎 Top 5 most relevant genes (highest GTEx significance):\n",
      top5_text, "\n\n",
      "Most represented tissue families:\n",
      fam_text, "\n\n",
      "Slope sign distribution (positive/negative):\n",
      slope_text, "\n"
    )
  })
  
  # ===========================
  # Manhattan combinat (GWAS + GTEx)
  # ===========================
  
  get_xref <- function(p) {
    xs <- unique(na.omit(vapply(p$x$data, function(tr) tr$xaxis %||% NA_character_, character(1))))
    if (length(xs) == 0) "x" else xs[1]
  }
  
  output$manhattan_combo <- renderPlotly({
    src_combo <- "manhattan_combo"
    
    dfp     <- tryCatch(dfp_manhattan(), error = function(e) NULL)
    gtex_dt <- tryCatch(gtex_data(),     error = function(e) NULL)
    
    if (is.null(dfp) || !is.data.frame(dfp) || !nrow(dfp)) {
      return(plotly_message("⚠️ GWAS table missing or incomplete. Load CHR, BP and P-values in Step 1."))
    }
    
    ref <- .ref_hg38
    ax  <- axis_df()
    axis_breaks <- ax$center
    axis_labels <- paste0("chr", ax$chrN)
    GENOME_END  <- max(ref$chr_cum + ref$len)
    
    # --- ref_map correcte (no factor->integer) ---
    ref_map <- .ref_hg38 %>%
      dplyr::mutate(
        chr_chr = as.character(chr),
        chr_chr = gsub("^chr", "", chr_chr, ignore.case = TRUE),
        chr_num = dplyr::case_when(
          chr_chr %in% as.character(1:22) ~ suppressWarnings(as.integer(chr_chr)),
          chr_chr == "X"  ~ 23L,
          chr_chr == "Y"  ~ 24L,
          chr_chr == "MT" ~ 26L,
          TRUE ~ NA_integer_
        )
      ) %>%
      dplyr::select(chr_num, chr_cum) %>%
      dplyr::filter(is.finite(chr_num), is.finite(chr_cum))
    
    #-----------
    thr_y <- if ((input$cluster_method %||% "window") == "window") (input$pthr %||% 5) else (input$min_logp %||% 6)
    
    # ----------------------------
    # p1: GWAS (robust)
    # ----------------------------
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
    
    # ----------------------------
    # p2 base: eixos correctes sempre (evita autoscale)
    # ----------------------------
    p2_base <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_blank() +
      ggplot2::scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::labs(x = "Genome", y = "GTEx signal") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    
    
    p2_pl <- plotly::ggplotly(p2_base, source = src_combo)
    
    y_max <- 1
    
    # ----------------------------
    # p2: GTEx hits (si n'hi ha)
    # ----------------------------
    gtex_ok <- !is.null(gtex_dt) && is.data.frame(gtex_dt) && nrow(gtex_dt) > 0
    if (gtex_ok) {
      
      # columnes mínimes
      validate(need(all(c("chr","pos") %in% names(gtex_dt)), "GTEx data must have 'chr' and 'pos' columns."))
      
      g2 <- gtex_dt %>%
        dplyr::mutate(
          chrN = norm_chr_generic(chr)
        ) %>%
        dplyr::inner_join(ref %>% dplyr::select(chr, chr_cum), by = c("chrN" = "chr")) %>%
        dplyr::mutate(BPcum = as.numeric(pos) + as.numeric(chr_cum))
      
      # yval + label (mateixa lògica que ja feies)
      if ("pval" %in% names(g2)) {
        g2$yval <- -log10(g2$pval)
        ylab <- "-log10(p_GTEx)"
      } else if ("pval_nominal" %in% names(g2)) {
        g2$yval <- -log10(g2$pval_nominal)
        ylab <- "-log10(p_GTEx)"
      } else if ("pvalue" %in% names(g2)) {
        g2$yval <- -log10(g2$pvalue)
        ylab <- "-log10(p_GTEx)"
      } else if ("slope" %in% names(g2)) {
        g2$yval <- g2$slope
        ylab <- "Slope"
      } else {
        num_cols <- names(g2)[vapply(g2, is.numeric, logical(1))]
        num_cols <- setdiff(num_cols, c("pos","BPcum","chr_cum","cluster","cluster_start","cluster_end"))
        first <- num_cols[1] %||% "pos"
        g2$yval <- g2[[first]]
        ylab <- first
      }
      
      # highlight cluster seleccionat (si existeix columna cluster)
      cl_sel <- selected_cluster()
      if (!is.null(cl_sel) && nrow(cl_sel) == 1 && "cluster" %in% names(g2) && "cluster" %in% names(cl_sel)) {
        g2$highlight <- ifelse(g2$cluster == cl_sel$cluster[1], "cluster", "other")
      } else {
        g2$highlight <- "other"
      }
      
      g2$tooltip <- paste0(
        "<b>GTEx hit</b>",
        "<br><b>chr:</b> ", g2$chr,
        "<br><b>pos:</b> ", g2$pos,
        if ("variant_id" %in% names(g2)) paste0("<br><b>variant_id:</b> ", g2$variant_id) else "",
        if ("tissue" %in% names(g2))     paste0("<br><b>tissue:</b> ", g2$tissue) else "",
        if ("gene_name" %in% names(g2))  paste0("<br><b>gene:</b> ", g2$gene_name) else "",
        if ("cluster_chr_n" %in% names(g2)) paste0("<br><b>cluster:</b> ", g2$cluster_chr_n) else ""
      )
      
      y_max <- max(g2$yval[is.finite(g2$yval)], na.rm = TRUE)
      if (!is.finite(y_max) || y_max <= 0) y_max <- 1
      
      p2 <- ggplot2::ggplot(g2, ggplot2::aes(x = BPcum, y = yval, text = tooltip)) +
        ggplot2::geom_point(ggplot2::aes(color = highlight), size = 0.7) +
        ggplot2::scale_color_manual(values = c("other" = "grey50", "cluster" = "red"), guide = "none") +
        ggplot2::scale_x_continuous(
          limits = c(0, GENOME_END),
          breaks = axis_breaks,
          labels = axis_labels,
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(limits = c(min(0, min(g2$yval, na.rm = TRUE)), y_max * 1.25), expand = c(0, 0)) +
        ggplot2::labs(x = "Genome", y = ylab) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      
      p2_pl <- plotly::ggplotly(p2, tooltip = "text", source = src_combo)
    }
    
    # ----------------------------
    # CLUSTER SEGMENT band (sempre visible)
    #   GTEx: clusters_cur() té chr/start/end i cluster_chr_n
    # ----------------------------
    # --- Cluster segment band (ROBUST) ---
    
    ##### cl <- clusters_cur() 
    cl <- clusters_stable()
    if (is.data.frame(cl) && nrow(cl) && is.data.frame(ref_map) && nrow(ref_map)) {
      
      clseg <- as.data.frame(cl) %>%
        dplyr::transmute(
          cluster_id = dplyr::coalesce(
            as.character(.data$cluster_id),
            as.character(.data$cluster_chr_n)
          ),
          chr_num = suppressWarnings(as.integer(.data$chr)),
          start_i = suppressWarnings(as.numeric(.data$start)),
          end_i   = suppressWarnings(as.numeric(.data$end))
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
            "<br>chr", chr_label_plink(chr_num), ":",
            format(start_i, scientific = FALSE), "-",
            format(end_i, scientific = FALSE)
          )
        ) %>%
        dplyr::filter(is.finite(x0), is.finite(x1), x1 >= x0) %>%
        dplyr::arrange(x0)
      
      # --- Després d'afegir els punts GTEx a p2_pl (o abans, és igual) ---
      
      if (nrow(clseg) > 0) {
        
        # --- Posició dins el panell inferior ---
        # Si hi ha GTEx hits, usem y_max (i queda prop de la part alta del panell).
        # Si NO n'hi ha, p2_base té y=[0..1], així que ho fixem dins d'aquest rang.
        if (isTRUE(gtex_ok)) {
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
        
        # IMPORTANT: fem servir TRACES (add_segments) com al codi EWAS
        # perquè subplot() les reassigna bé al panell inferior.
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
            textangle = 45,          # <- 90º com vols
            font = list(size = 12),
            xanchor = "center",
            yanchor = "bottom"
          )
      }
    }
    
    out <- plotly::subplot(
      p1_pl, p2_pl,
      nrows = 2, shareX = TRUE,
      heights = c(0.55, 0.45),
      titleY = TRUE
    )
    
    # IMPORTANT: fixa source al widget FINAL (subplot sovint el "perd")
    out$x$source <- src_combo
    
    # registra events al widget FINAL
    #   out <- plotly::event_register(out, "plotly_click")
    out <- plotly::event_register(out, "plotly_relayout")
    
    out %>%
      plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d","lasso2d","hoverCompareCartesian")
      ) %>%
      plotly::layout(
        margin = list(t = 60),
        showlegend = FALSE,
        annotations = list(
          list(x=0.5, y=1.05, text="<b>GWAS Manhattan</b>", showarrow=FALSE, xref="paper", yref="paper"),
          list(x=0.5, y=0.47, text="<b>GTEx eQTLs inside clusters</b>", showarrow=FALSE, xref="paper", yref="paper")
        )
      )%>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "Combined_Mantattan_plot",
          width = 1400,
          height = 900,
          scale = 2
        )
      )
  })
  
  
  # ===========================
  # UCSC viewer (GWAS + GTEx)  [ADAPTAT A subplot + source="manhattan_combo"]
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
  
  safe_label_gtex <- function(df) {
    if ("variant_id" %in% names(df) && !all(is.na(df$variant_id))) {
      out <- as.character(df$variant_id)
      bad <- is.na(out) | out == ""
      out[bad] <- paste0("gtex_hit_", seq_len(sum(bad)))
      return(out)
    }
    paste0("gtex_hit_", seq_len(nrow(df)))
  }
  
  window_selected <- reactiveVal(NULL)
  ucsc_region     <- reactiveVal(NULL)
  
  # --- helper robust: extreu rang X d'un relayout (accepta list o data.frame) ---
  extract_xrange_from_relayout <- function(ev) {
    
    if (is.null(ev)) return(NULL)
    
    # event_data pot ser data.frame o llista
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
      session$userData$track_gwas_data <- NULL
      session$userData$track_gtex_data <- NULL
      return()
    }
    window_selected(xr)
  }, ignoreInit = TRUE)
  
  get_window_hits_gwas <- reactive({
    win <- window_selected(); req(win)
    
    dfp <- dfp_manhattan()
    df  <- gwas_df()
    req(is.data.frame(dfp), nrow(dfp) > 0)
    req(is.data.frame(df),  nrow(df)  > 0)
    
    dfp2 <- dfp %>% dplyr::mutate(
      CHR = suppressWarnings(as.integer(CHR)),
      BP  = suppressWarnings(as.integer(BP))
    )
    
    df2 <- df %>%
      dplyr::mutate(
        CHR = suppressWarnings(as.integer(CHR)),
        BP  = suppressWarnings(as.integer(BP))
      ) %>%
      dplyr::select(dplyr::any_of(c("CHR","BP","rsid","snp")))
    
    dfp2 %>%
      dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax) %>%
      dplyr::left_join(df2, by = c("CHR","BP"))
  })
  
  get_window_hits_gtex <- reactive({
    win <- window_selected(); req(win)
    gtex <- gtex_data()
    req(!is.null(gtex), nrow(gtex) > 0)
    
    ref <- .ref_hg38
    gtex2 <- gtex %>%
      dplyr::mutate(chrN = norm_chr_generic(chr)) %>%
      dplyr::inner_join(ref %>% dplyr::select(chr, chr_cum), by = c("chrN" = "chr")) %>%
      dplyr::mutate(BPcum = as.numeric(pos) + as.numeric(chr_cum))
    
    gtex2 %>% dplyr::filter(BPcum >= win$xmin, BPcum <= win$xmax)
  })
  
  make_track_df <- function(df, R, G, B, type = c("gwas", "gtex")) {
    type <- match.arg(type)
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    
    if (type == "gwas") {
      if (!"CHR" %in% names(df) || !"BP" %in% names(df)) return(tibble::tibble())
      chrom <- paste0("chr", chr_label_plink(as.integer(df$CHR)))
      start <- as.numeric(df$BP) - 1
      end   <- as.numeric(df$BP)
      lbl   <- safe_label_rsid(df)
    } else {
      # GTEx: usa chrN si existeix (ja normalitzat a 1..22/23..)
      if (!"pos" %in% names(df)) return(tibble::tibble())
      chr_use <- if ("chrN" %in% names(df)) df$chrN else df$chr
      chr_use <- as.character(chr_use)
      chr_use <- gsub("^chr", "", chr_use, ignore.case = TRUE)
      chrom   <- paste0("chr", chr_use)
      
      start <- as.numeric(df$pos) - 1
      end   <- as.numeric(df$pos)
      lbl   <- safe_label_gtex(df)
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
  
  observeEvent(window_selected(), {
    win <- window_selected()
    if (is.null(win)) return()
    
    left  <- coord_from_bp_cum(win$xmin)
    right <- coord_from_bp_cum(win$xmax)
    
    # IMPORTANT: UCSC region ha d’estar en un sol cromosoma
    if (is.null(left$chr) || is.null(right$chr) || left$chr != right$chr) {
      ucsc_region(NULL)
      session$userData$track_gwas_data <- NULL
      session$userData$track_gtex_data <- NULL
      return()
    }
    
    region <- sprintf("chr%s:%d-%d", left$chr, left$pos, right$pos)
    ucsc_region(region)
    
    gwas_hits <- try(get_window_hits_gwas(), silent = TRUE)
    gtex_hits <- try(get_window_hits_gtex(), silent = TRUE)
    if (inherits(gwas_hits, "try-error")) gwas_hits <- NULL
    if (inherits(gtex_hits, "try-error")) gtex_hits <- NULL
    
    if ((is.null(gwas_hits) || nrow(gwas_hits) == 0) &&
        (is.null(gtex_hits) || nrow(gtex_hits) == 0)) {
      session$userData$track_gwas_data <- NULL
      session$userData$track_gtex_data <- NULL
      return()
    }
    
    df_gwas <- clean_track(make_track_df(gwas_hits,  31,120,180, type = "gwas"))
    df_gtex <- clean_track(make_track_df(gtex_hits, 227, 26, 28, type = "gtex"))
    
    session$userData$track_gwas_data <- df_gwas
    session$userData$track_gtex_data <- df_gtex
  }, ignoreInit = TRUE)
  
  make_ucsc_track_text_gwas <- function(name, df) {
    header <- sprintf(
      'track name="%s" description="%s" visibility=pack itemRgb=On url="https://www.ncbi.nlm.nih.gov/snp/$$"',
      name, name
    )
    if (is.null(df) || !nrow(df)) return(header)
    
    cols_bed12 <- c("chrom","start","end","name","score","strand",
                    "thickStart","thickEnd","itemRgb",
                    "blockCount","blockSizes","blockStarts")
    body <- apply(df[, cols_bed12], 1, paste, collapse = "\t")
    paste(c(header, body), collapse = "\n")
  }
  
  make_ucsc_track_text_gtex <- function(name, df) {
    header <- sprintf(
      'track name="%s" description="%s" visibility=pack itemRgb=On url="https://www.gtexportal.org/home/snp/$$"',
      name, name
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
    
    txt <- make_ucsc_track_text_gwas("GWAS_hits", df)
    url <- make_ucsc_url(region, txt)
    
    tags$a(href = url, target = "_blank", "Open UCSC – GWAS hits (dbSNP links)")
  })
  
  output$ucsc_link_gtex <- renderUI({
    region <- ucsc_region()
    df     <- session$userData$track_gtex_data
    req(!is.null(region), !is.null(df), nrow(df) > 0)
    
    txt <- make_ucsc_track_text_gtex("GTEx_hits", df)
    url <- make_ucsc_url(region, txt)
    
    tags$a(href = url, target = "_blank", "Open UCSC – GTEx hits (clickable)")
  })
  
  output$debug_ucsc_state <- renderUI({
    region <- ucsc_region() %||% "NULL"
    gwas_n <- nrow(session$userData$track_gwas_data %||% tibble::tibble())
    gtex_n <- nrow(session$userData$track_gtex_data %||% tibble::tibble())
    tags$pre(
      style="background-color:#f6f6f6; border:1px solid #ddd; padding:10px; font-family: 'Courier New', monospace;",
      paste0(
        "region = ", region,
        "\nGWAS hits in track = ", gwas_n,
        "\nGTEx hits in track = ", gtex_n
      )
    )
  })
  
  # ===========================
  # Modal info input
  # ===========================
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
  
  # ===========================
  # Downloads (GTEx hits by cluster) amb nom mode+threshold
  # ===========================
  
  output$dl_gtex_hits_csv <- downloadHandler(
    filename = function() {
      tag <- .make_mode_thr_tag(input)
      paste0(
        "gtex_hits_by_cluster_",
        format(Sys.Date(), "%Y%m%d"),
        "_", tag$mode_tag,
        "_thr", tag$thr_tag,
        ".csv"
      )
    },
    content = function(file) {
      dt <- gtex_data()
      if (is.null(dt) || !nrow(dt)) {
        writeLines("No GTEx hits extracted.", con = file)
      } else {
        readr::write_csv(dt, file)
      }
    }
  )
  
  output$dl_gtex_hits_rds <- downloadHandler(
    filename = function() {
      tag <- .make_mode_thr_tag(input)
      paste0(
        "gtex_hits_by_cluster_",
        format(Sys.Date(), "%Y%m%d"),
        "_", tag$mode_tag,
        "_thr", tag$thr_tag,
        ".rds"
      )
    },
    content = function(file) {
      dt <- gtex_data()
      if (is.null(dt) || !nrow(dt)) {
        saveRDS(list(message = "No GTEx hits extracted."), file = file)
      } else {
        saveRDS(dt, file = file)
      }
    }
  )
  
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
  
  # ============================================================
  # SERVER — GTEx Enrichment (Catalog/NonSyn-style)  [ADAPTAT]
  # Manté denominacions GTEx: gene_id (ENSEMBL), gtex_table_filtered(), clusters_cur()
  # ============================================================
  
  pick_col <- function(df, candidates) {
    nm <- names(df)
    hit <- candidates[candidates %in% nm]
    if (length(hit)) hit[1] else NULL
  }
  
  strip_ens_version <- function(x) {
    x <- as.character(x)
    sub("\\..*$", "", x)
  }
  
  # ------------------------------------------------------------
  # UI cluster picker (com Catalog), però amb clusters_cur() i cluster_chr_n (GTEx)
  # ------------------------------------------------------------
  output$func_cluster_ui <- renderUI({
    # cl <- clusters_cur()
    cl <- clusters_stable()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available."))
    
    validate(need("chr" %in% names(cl), "clusters_cur() must contain column 'chr'."))
    validate(need("cluster_chr_n" %in% names(cl), "clusters_cur() must contain column 'cluster_chr_n'."))
    
    chr_choices <- sort(unique(as.integer(cl$chr)))
    chr_labels  <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("func_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("func_cluster_id_ui")
    )
  })
  
  output$func_cluster_id_ui <- renderUI({
    # cl <- clusters_cur()
    cl <- clusters_stable()
    req(is.data.frame(cl), nrow(cl) > 0, input$func_chr)
    
    x <- cl %>%
      dplyr::filter(as.integer(chr) == as.integer(input$func_chr)) %>%
      dplyr::arrange(start, end)
    
    validate(need(nrow(x) > 0, "No clusters on this chromosome."))
    
    lab <- paste0(x$cluster_chr_n, " (", x$start, "-", x$end, ")")
    choices <- stats::setNames(x$cluster_chr_n, lab)
    
    selectInput("func_cluster_id", "Cluster", choices = choices, selected = x$cluster_chr_n[1])
  })
  
  output$func_gene_ui <- renderUI({
    dt <- tryCatch(gtex_table_all(), error = function(e) NULL)
    if (!is.data.frame(dt) || !nrow(dt)) return(NULL)
    
    gene_id_col <- if ("gene_id" %in% names(dt)) "gene_id" else NULL
    if (is.null(gene_id_col)) return(helpText("No gene_id column available in GTEx hits."))
    
    # etiqueta si hi ha gene_name
    if ("gene_name" %in% names(dt)) {
      labs <- unique(paste0(dt$gene_name, " (", dt$gene_id, ")"))
      vals <- unique(dt$gene_id)
      # assegura mateixa longitud (si hi ha duplicats inconsistents, fem mapping pel primer)
      mp <- dt %>%
        dplyr::filter(!is.na(gene_id), nzchar(gene_id)) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::summarise(gene_name = dplyr::first(gene_name), .groups="drop")
      choices <- stats::setNames(mp$gene_id, paste0(mp$gene_name, " (", mp$gene_id, ")"))
      selectInput("func_gene_id", "Gene", choices = choices)
    } else {
      choices <- sort(unique(as.character(dt$gene_id)))
      selectInput("func_gene_id", "Gene", choices = choices)
    }
  })
  
  # ------------------------------------------------------------
  # Dades GTEx per enrichment
  #   - IMPORTANT: assumim que gtex_table_filtered() retorna GTEx hits (global o no)
  #   - Per CLUSTER mode filtrarem nosaltres per cluster_chr_n / cluster id column
  # ------------------------------------------------------------
  gtex_hits_all <- reactive({
    dt <- gtex_table_all()
    if (!is.data.frame(dt) || !nrow(dt)) return(NULL)
    dt
  })
  
  cluster_col_gtex <- reactive({
    dt <- gtex_hits_all()
    if (is.null(dt)) return(NULL)
    pick_col(dt, c("cluster_chr_n", "cluster_id", "cluster", "cluster_chr"))
  })
  
  genes_all_ensembl <- reactive({
    dt <- gtex_hits_all()
    if (is.null(dt) || !nrow(dt)) return(character(0))
    validate(need("gene_id" %in% names(dt), "GTEx hits missing 'gene_id' column."))
    g <- unique(strip_ens_version(dt$gene_id))
    g <- g[!is.na(g) & nzchar(g)]
    g
  })
  
  genes_scope_ensembl <- reactive({
    dt <- gtex_hits_all()
    validate(need(is.data.frame(dt) && nrow(dt) > 0, "No GTEx hits loaded."))
    
    validate(need("gene_id" %in% names(dt), "GTEx hits missing 'gene_id' column."))
    
    if (identical(input$func_scope, "cluster")) {
      cc <- cluster_col_gtex()
      validate(need(!is.null(cc), "GTEx hits table has no cluster column (cluster_chr_n/cluster_id/cluster/cluster_chr)."))
      req(input$func_cluster_id)
      dt <- dt %>% dplyr::filter(.data[[cc]] == input$func_cluster_id)
    }
    
    g <- unique(strip_ens_version(dt$gene_id))
    g <- g[!is.na(g) & nzchar(g)]
    g
  })
  
  # ------------------------------------------------------------
  # Mapatge ENSEMBL -> ENTREZ una sola vegada (cache)
  # ------------------------------------------------------------
  gene_map_cache_ens <- reactiveVal(NULL)
  gene_map_cache_map <- reactiveVal(NULL)
  
  gene_map_all <- reactive({
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
    ens <- genes_all_ensembl()
    validate(need(length(ens) > 0, "No valid ENSEMBL gene_ids found in GTEx hits."))
    
    prev_ens <- gene_map_cache_ens()
    prev_map <- gene_map_cache_map()
    
    if (!is.null(prev_ens) && !is.null(prev_map) && identical(sort(prev_ens), sort(ens))) {
      return(prev_map)
    }
    
    m <- tryCatch({
      suppressMessages(
        clusterProfiler::bitr(
          ens,
          fromType = "ENSEMBL",
          toType   = "ENTREZID",
          OrgDb    = org.Hs.eg.db
        )
      )
    }, error = function(e) NULL)
    
    validate(need(is.data.frame(m) && nrow(m) > 0, "ENSEMBL → ENTREZ mapping failed (empty)."))
    
    gene_map_cache_ens(sort(ens))
    gene_map_cache_map(m)
    m
  })
  
  entrez_scope <- reactive({
    m <- gene_map_all()
    ens <- genes_scope_ensembl()
    validate(need(length(ens) > 0, "No valid genes for the selected scope."))
    
    ids <- unique(m$ENTREZID[m$ENSEMBL %in% ens])
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
    
    # GLOBAL: background per defecte (OrgDb/KEGG)
    if (identical(scope, "global")) return(NULL)
    
    # CLUSTER: usuari decideix
    if (identical(bg, "dataset")) return(universe_entrez_dataset())
    
    NULL
  })
  
  output$enrich_bg_note <- renderUI({
    bg  <- input$enrich_background %||% "dataset"
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    
    # -----------------------------
    # Helpers (local)
    # -----------------------------
    bg_lbl <- if (identical(bg, "dataset")) "Dataset genes" else "Default background"
    
    sc <- input$func_scope %||% "global"
    sc_txt <- if (identical(sc, "cluster")) {
      paste0("Foreground: genes in <b>", htmltools::htmlEscape(input$func_cluster_id %||% ""), "</b>.")
    } else {
      "Foreground: genes in <b>global</b> (all extracted hits)."
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
    
    # foreground gene counts (ENTREZIDs used by GO/KEGG/GoSlim)
    n_fg <- safe_n(quote(length(entrez_scope())))
    
    # universe gene counts for "dataset" bg (mapped dataset universe)
    n_uni_dataset <- safe_n(quote(length(universe_entrez_dataset())))
    
    # universe gene counts for default backgrounds (OrgDb / KEGG)
    n_uni_orgdb <- safe_n(quote(length(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID"))))
    
    # choose N text
    uni_txt <- function(default_name) {
      if (identical(bg, "dataset")) {
        if (is.finite(n_uni_dataset)) {
          paste0("<b>Background (universe):</b> dataset-mapped genes (N=", n_uni_dataset, ").")
        } else {
          "<b>Background (universe):</b> dataset-mapped genes."
        }
      } else {
        # default background
        if (is.finite(n_uni_orgdb)) {
          paste0("<b>Background (universe):</b> ", default_name, " (N≈", n_uni_orgdb, ").")
        } else {
          paste0("<b>Background (universe):</b> ", default_name, ".")
        }
      }
    }
    
    fg_cnt_txt <- if (is.finite(n_fg)) paste0("<br><b>Foreground size:</b> n=", n_fg, " genes.") else ""
    
    # GSSize note (common)
    gs_txt <- "<b>GSSize:</b> number of genes per term in the selected universe (filtered by minGSSize / maxGSSize)."
    
    # -----------------------------
    # Build message per tab
    # -----------------------------
    msg <- if (identical(tab, "tab_enrich_go")) {
      
      # GO
      base_uni <- if (identical(bg, "orgdb")) "OrgDb annotated genes (org.Hs.eg.db)" else "dataset-mapped genes"
      tip <- "<br><b>Test unit:</b> genes (ENTREZID)."
      
      paste0(
        "<b>GO enrichment</b><br>",
        sc_txt,
        fg_cnt_txt, "<br>",
        uni_txt("OrgDb annotated genes (org.Hs.eg.db)"), "<br>",
        tip, "<br>",
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
          paste0("<b>Background (universe):</b> dataset-mapped genes (recommended in cluster mode).")
        } else {
          "GoSlim uses an OrgDb-derived GO→GoSlim TERM2GENE; universe is default unless a dataset universe is provided."
        },
        if (identical(bg, "dataset") && is.finite(n_uni_dataset)) paste0(" (N=", n_uni_dataset, ").") else "",
        "<br>",
        "<b>Test unit:</b> genes (ENTREZID).<br>",
        gs_txt, "<br>",
        "<b>Method:</b> clusterProfiler::enricher with GoSlim TERM2GENE (GO ancestors mapped to <i>goslim_generic</i>).",
        tip
      )
      
    } else if (identical(tab, "tab_enrich_gtex_terms")) {
      
      # (Keep your existing GTEx gene/tissue note unchanged)
      sc2 <- input$func_scope %||% "global"
      sc_txt2 <- if (identical(sc2, "cluster")) {
        paste0("Foreground: GTEx eQTL rows inside cluster <b>",
               htmltools::htmlEscape(input$func_cluster_id %||% ""),
               "</b>.")
      } else {
        "Foreground: GTEx eQTL rows inside <b>all clusters</b>."
      }
      
      n_uni <- tryCatch(nrow(gtex_all()), error = function(e) NA_integer_)
      n_fg2 <- tryCatch({ xx <- dt_enrich_res_gtex(); nrow(xx$dt2) }, error = function(e) NA_integer_)
      
      cnt_txt <- if (is.finite(n_uni) && is.finite(n_fg2)) {
        paste0("<br><b>Rows:</b> foreground n=", n_fg2, " · universe N=", n_uni, ".")
      } else {
        ""
      }
      
      paste0(
        "<b>GTEx Gene/Tissue enrichment</b><br>",
        sc_txt2, "<br>",
        "<b>Background (universe):</b> all GTEx eQTL rows from <code>www/gtex_eqtl_light.rds</code>.",
        cnt_txt, "<br>",
        "<b>Test unit:</b> eQTL rows (not unique genes/tissues).<br>",
        "<b>GSSize:</b> number of eQTL rows per term in the universe."
      )
      
    } else {
      paste0("<b>Background:</b> ", bg_lbl)
    }
    
    htmltools::HTML(paste0(
      "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
      msg, "</div>"
    ))
  })
  
  
  # ------------------------------------------------------------
  # Triggers (eviten execució en entrar al panell), com al Catalog
  # ------------------------------------------------------------
  go_trigger   <- reactiveVal(NULL)
  kegg_trigger <- reactiveVal(NULL)
  
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
  
  # ------------------------------------------------------------
  # GO enrichment (multi-ontology) — pvalueCutoff=1 i filtrem després per FDR cutoff
  # ------------------------------------------------------------
  go_enrich_raw <- eventReactive(go_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db package is required."))
    
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
    
    DT::datatable(out, rownames = FALSE,       
                  extensions = "Buttons",
                  options = list(
                    dom="Bfrtip", 
                    buttons = list(
                      list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                    ),
                    pageLength=10, scrollX=TRUE)
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
        Ontology %in% c("BP","CC","MF"),
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
        Ontology = factor(as.character(Ontology), levels = c("BP","CC","MF")),
        tooltip  = paste0(
          "<b>", htmltools::htmlEscape(Description), "</b>",
          "<br>Ontology: ", Ontology,
          "<br>-log10(pvalue): ", sprintf("%.3f", score),
          "<br>FDR: ", ifelse(is.finite(p_adj), formatC(p_adj, format="e", digits=2), "NA")
        )
      )
    
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
          fill = Ontology,
          text = tooltip
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(title = onto, x = NULL, y = "-log10(pvalue)") +
        ggplot2::scale_fill_manual(values = cols_go, guide = "none", drop = FALSE) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    }
    
    p <- NULL
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- (mk("BP") | mk("CC") | mk("MF")) +
        patchwork::plot_annotation(title = ttl, subtitle = subttl)
    } else {
      p <- ggplot2::ggplot(
        dat,
        ggplot2::aes(
          x    = stats::reorder(term_short, score),
          y    = score,
          fill = Ontology,
          text = tooltip
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~Ontology, scales = "free_y") +
        ggplot2::scale_fill_manual(values = cols_go, guide = "none", drop = FALSE) +
        ggplot2::labs(title = ttl, subtitle = subttl, x = NULL, y = "-log10(pvalue)") +
        ggplot2::theme_minimal(base_size = 12)
    }
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        margin = list(l = 210, r = 20, t = 60, b = 50),
        font   = list(size = 11)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "png",
          filename = "GO_barplot",
          width = 1600,
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
      
      uni  <- universe_entrez_for_scope()
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
    
    DT::datatable(out, rownames = FALSE,       
                  extensions = "Buttons",
                  options = list(
                    dom="Bfrtip", 
                    buttons = list(
                      list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                      list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                    ),
                    pageLength=10, scrollX=TRUE))
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
    
    cat("\n[DBG] head(kegg_top_df()):\n")
    print(utils::head(df, 10))
    cat("\n[DBG] names(kegg_top_df()):\n")
    print(names(df))
    
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
        padj  = suppressWarnings(as.numeric(p.adjust)),
        tooltip = paste0(
          "<b>", htmltools::htmlEscape(Description), "</b>",
          "<br>-log10(pvalue): ", sprintf("%.3f", score),
          "<br>FDR: ", ifelse(is.finite(padj), formatC(padj, format="e", digits=2), "NA")
        )
      )
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x    = stats::reorder(term_short, score),
        y    = score,
        fill = padj,
        text = tooltip
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(title = ttl, x = NULL, y = "-log10(pvalue)", fill = "FDR") +
      ggplot2::scale_fill_gradient(
        low  = "red",
        high = "orange",
        na.value = "yellow",
        trans = "reverse"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        margin = list(l = 240, r = 20, t = 60, b = 50),
        font   = list(size = 11)
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
  
  ############################################################################
  ############################################################################
  ############################## GENE / TISSUE ENRICHMENT ######################
  # Universe: gtex_all()  (www/gtex_eqtl_light.rds)
  # Foreground: gtex_table_all() == gtex_data() (hits inside clusters)
  # Scope: global / cluster (input$func_scope + input$func_cluster_id)
  # Triggered by: input$run_enrich when enrich_tabs == tab_enrich_gtex_terms
  ############################################################################
  
  # ------------------------------------------------------------
  # Helpers (Catalog-style)
  # ------------------------------------------------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  pick_col <- function(df, candidates) {
    nm <- names(df)
    hit <- candidates[candidates %in% nm]
    if (length(hit)) hit[1] else NULL
  }
  
  short_label <- function(x, max_chars = 35, wrap_width = 14) {
    x <- as.character(x)
    x <- stringr::str_squish(x)
    x <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "…"), x)
    stringr::str_wrap(x, width = wrap_width)
  }
  # ------------------------------------------------------------
  
  norm_tissue_id <- function(x) {
    x <- trimws(as.character(x))
    x <- toupper(x)
    # converteix qualsevol separador (espais, '-', '/', etc.) a '_'
    x <- gsub("[^A-Z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  # ------------------------------------------------------------
  # Fisher enrichment for terms (subset vs universe) on ROWS
  # (NO unique() used for the test)
  # Fast version: hypergeometric tail (equivalent to Fisher 1-sided "greater")
  # ------------------------------------------------------------
  fisher_enrich_terms_fast <- function(subset_terms, uni_terms,
                                       alternative = "greater",
                                       min_a = 1, min_Tt = 1, max_Tt = Inf,
                                       max_terms = 3000) {
    
    sub <- na.omit(as.character(subset_terms))
    uni <- na.omit(as.character(uni_terms))
    if (!length(sub) || !length(uni)) return(data.table::data.table())
    
    tabU <- table(uni)
    tabS <- table(sub)
    
    terms <- intersect(names(tabU), names(tabS))
    if (!length(terms)) return(data.table::data.table())
    
    Tt <- as.integer(tabU[terms])   # universe count per term (rows)
    a  <- as.integer(tabS[terms])   # subset count per term (rows)
    
    keep <- which(a >= min_a & Tt >= min_Tt & Tt <= max_Tt)
    if (!length(keep)) return(data.table::data.table())
    
    terms <- terms[keep]; a <- a[keep]; Tt <- Tt[keep]
    
    if (length(terms) > max_terms) {
      ord <- order(Tt, decreasing = TRUE)[1:max_terms]
      terms <- terms[ord]; a <- a[ord]; Tt <- Tt[ord]
    }
    
    n <- length(sub)  # subset rows
    N <- length(uni)  # universe rows
    
    p <- stats::phyper(a - 1L, Tt, N - Tt, n, lower.tail = FALSE)
    
    b <- n - a
    c <- Tt - a
    d <- (N - Tt) - b
    or <- rep(NA_real_, length(a))
    ok <- (b > 0 & c > 0 & d > 0)
    or[ok] <- (a[ok] * d[ok]) / (b[ok] * c[ok])
    
    out <- data.table::data.table(term = terms, a = a, Tt = Tt, or = or, p = p)
    out[, p_adj := p.adjust(p, method = "BH")]
    
    # NEW: totals for ratios
    out[, `:=`(n_total = n, N_total = N)]
    
    data.table::setorder(out, p_adj, p)
    out
  }
  
  # ------------------------------------------------------------
  # Bundle (Catalog-style): dt_uni, dt (all hits), dt2 (scope)
  # ------------------------------------------------------------
  dt_enrich_res_gtex <- reactive({
    dt_uni <- gtex_all()
    dt     <- gtex_table_all()   # == gtex_data()
    
    validate(need(is.data.frame(dt_uni) && nrow(dt_uni) > 0,
                  "GTEx universe not loaded (www/gtex_eqtl_light.rds)."))
    validate(need(is.data.frame(dt) && nrow(dt) > 0,
                  "No GTEx hits yet. Run Step 3/4 (Extract GTEx)."))
    
    dt_uni <- as.data.frame(dt_uni)
    dt     <- as.data.frame(dt)
    
    # validate columns exist
    validate(need("tissue" %in% names(dt),      "GTEx hits missing 'tissue' column."))
    validate(need("tissue" %in% names(dt_uni),  "GTEx universe missing 'tissue' column."))
    
    # normalize types
    dt$tissue      <- trimws(as.character(dt$tissue))
    dt_uni$tissue  <- trimws(as.character(dt_uni$tissue))
    
    # Normalize cluster id strings (Hub sometimes prepends cluster_)
    if ("cluster_chr_n" %in% names(dt)) {
      dt$cluster_chr_n <- sub("^cluster_", "", as.character(dt$cluster_chr_n))
    }
    cid_in <- sub("^cluster_", "", as.character(input$func_cluster_id %||% ""))
    
    # ---- scope ----
    dt2 <- dt
    
    if (identical(input$func_scope, "gene")) {
      req(input$func_gene_id)
      if ("gene_id" %in% names(dt2)) {
        dt2 <- dt2 %>% dplyr::filter(as.character(gene_id) == as.character(input$func_gene_id))
      } else {
        validate(need(FALSE, "GTEx hits missing gene_id column (cannot filter by gene)."))
      }
    }
    
    if (identical(input$func_scope, "cluster")) {
      req(nzchar(cid_in))
      if ("cluster_chr_n" %in% names(dt2) && any(nzchar(dt2$cluster_chr_n))) {
        dt2 <- dt2 %>% dplyr::filter(as.character(cluster_chr_n) == cid_in)
      } else if ("cluster_id" %in% names(dt2)) {
        dt2 <- dt2 %>% dplyr::filter(sub("^cluster_", "", as.character(cluster_id)) == cid_in)
      } else if ("cluster" %in% names(dt2)) {
        dt2 <- dt2 %>% dplyr::filter(as.integer(cluster) == suppressWarnings(as.integer(cid_in)))
      }
    }
    
    validate(need(is.data.frame(dt2) && nrow(dt2) > 0,
                  "No GTEx hits in this scope (cluster empty)."))
    
    list(dt = dt, dt2 = dt2, dt_uni = dt_uni)
  })
  
  # ------------------------------------------------------------
  # Trigger (same pattern as GO/KEGG/GoSlim)
  # ------------------------------------------------------------
  gtex_terms_trigger <- reactiveVal(NULL)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    if (identical(tab, "tab_enrich_gtex_terms")) {
      gtex_terms_trigger((gtex_terms_trigger() %||% 0L) + 1L)
    }
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # Compute ONCE per click (so DT doesn't "spin forever")
  # ------------------------------------------------------------
  norm_tissue_id <- function(x) {
    x <- trimws(as.character(x))
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)  # espais/guions/etc -> "_"
    x <- gsub("^_+|_+$", "", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  gtex_terms_res <- eventReactive(gtex_terms_trigger(), {
    x <- dt_enrich_res_gtex()
    
    dt2    <- x$dt2
    dt_uni <- x$dt_uni
    
    tissue_col     <- pick_col(dt2,    c("tissue"))
    tissue_col_uni <- pick_col(dt_uni, c("tissue"))
    
    
    validate(need(!is.null(tissue_col) && !is.null(tissue_col_uni),
                  "Tissue column not found (expected: tissue)."))
    
    
    minGS <- input$enrich_min_gs %||% 1
    maxGS <- input$enrich_max_gs %||% 500000
    
    withProgress(message = "Running GTEx tissue enrichment…", value = 0, {
      
      # ----------------------------
      # TISSUE (normalize + fallback to tissue_family if overlap==0)
      # ----------------------------
      incProgress(0.2, detail = "Counting tissues…")
      
      fg_tis_raw <- dt2[[tissue_col]]
      un_tis_raw <- dt_uni[[tissue_col_uni]]
      
      fg_id <- norm_tissue_id(fg_tis_raw)
      un_id <- norm_tissue_id(un_tis_raw)
      
      fg_ok <- fg_id[!is.na(fg_id) & nzchar(fg_id)]
      un_ok <- un_id[!is.na(un_id) & nzchar(un_id)]
      
      cat("[TISSUE] overlap_norm=",
          length(intersect(unique(fg_ok), unique(un_ok))),
          " example_fg=", paste(head(unique(fg_ok), 5), collapse = ","),
          " example_uni=", paste(head(unique(un_ok), 5), collapse = ","),
          "\n"
      )
      
      tab_tissue <- data.table::data.table()
      tissue_msg <- NULL
      
      if (!length(fg_ok)) {
        tissue_msg <- "Foreground has 0 valid tissue values (all NA/empty)."
      } else if (!length(un_ok)) {
        tissue_msg <- "Universe has 0 valid tissue values."
      } else {
        
        # check overlap after normalization
        ov <- intersect(unique(fg_ok), unique(un_ok))
        
        if (!length(ov)) {
          # ---- FALLBACK: use tissue families (Brain/Immune/...) on BOTH sides ----
          if (exists("map_tissue_family", mode = "function")) {
            
            fg_fam <- vapply(trimws(as.character(fg_tis_raw)), map_tissue_family, FUN.VALUE = character(1))
            un_fam <- vapply(trimws(as.character(un_tis_raw)), map_tissue_family, FUN.VALUE = character(1))
            
            fg_fam <- fg_fam[!is.na(fg_fam) & nzchar(fg_fam)]
            un_fam <- un_fam[!is.na(un_fam) & nzchar(un_fam)]
            
            tab_tissue <- fisher_enrich_terms_fast(
              fg_fam, un_fam,
              alternative = "greater",
              min_a  = 1,
              min_Tt = minGS,
              max_Tt = Inf,
              max_terms = 200
            )
            
            if (!nrow(tab_tissue)) {
              tissue_msg <- "No overlap for tissues; fallback to tissue families also returned 0 terms."
            } else {
              tab_tissue$label <- tab_tissue$term
              tissue_msg <- "No overlap for tissues; showing enrichment on tissue families instead."
            }
            
          } else {
            tissue_msg <- "Tissue overlap is 0 even after normalization. (No map_tissue_family() found for fallback.)"
          }
          
        } else {
          # ---- NORMAL PATH: enrich on normalized tissue IDs ----
          tab_tissue <- fisher_enrich_terms_fast(
            fg_ok, un_ok,
            alternative = "greater",
            min_a  = 1,
            min_Tt = minGS,
            max_Tt = Inf,
            max_terms = 5000
          )
          
          if (!nrow(tab_tissue)) {
            tissue_msg <- "Tissue enrichment returned 0 terms after normalization (unexpected)."
          } else {
            # label map: normalized id -> original universe label (first)
            u_lab <- data.frame(
              term  = un_id,
              label = trimws(as.character(un_tis_raw)),
              stringsAsFactors = FALSE
            )
            u_lab <- u_lab[!is.na(u_lab$term) & nzchar(u_lab$term), , drop = FALSE]
            u_lab <- u_lab[!duplicated(u_lab$term), , drop = FALSE]
            
            tab_tissue <- dplyr::left_join(tab_tissue, u_lab, by = "term")
            miss <- is.na(tab_tissue$label) | !nzchar(tab_tissue$label)
            tab_tissue$label[miss] <- tab_tissue$term[miss]
          }
        }
      }
      
      list(
        tab_tissue = tab_tissue,
        tissue_msg = tissue_msg
      )
    })
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # Tissue table + plot (render from computed result)
  # ------------------------------------------------------------
  output$enrich_tissue_table <- DT::renderDT({
    res <- gtex_terms_res(); req(res)
    tab <- res$tab_tissue
    msg <- res$tissue_msg %||% "No tissues pass filters."
    
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(DT::datatable(data.frame(Message = msg), options = list(dom="t"), rownames = FALSE))
    }
    
    # FIX: ensure label exists
    if (!"label" %in% names(tab)) tab$label <- tab$term
    
    # totals (agafem el primer valor disponible)
    n_tot <- if ("n_total" %in% names(tab)) suppressWarnings(as.integer(tab$n_total[1])) else NA_integer_
    N_tot <- if ("N_total" %in% names(tab)) suppressWarnings(as.integer(tab$N_total[1])) else NA_integer_
    
    out <- tab %>%
      dplyr::transmute(
        tissue = dplyr::coalesce(.data$label, .data$term),
        
        in_subset_ratio = if (is.finite(n_tot)) {
          sprintf("%d/%d (%.2f%%)", a, n_tot, 100 * a / max(1, n_tot))
        } else {
          as.character(a)
        },
        
        in_universe_ratio = if (is.finite(N_tot)) {
          sprintf("%d/%d (%.3f%%)", Tt, N_tot, 100 * Tt / max(1, N_tot))
        } else {
          as.character(Tt)
        },
        
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format="e", digits=2),
        FDR = formatC(as.numeric(p_adj), format="e", digits=2)
      )
    
    DT::datatable(out, rownames=FALSE,
                  extensions="Buttons",
                  options=list(dom="Bfrtip", 
                               buttons = list(
                                 list( extend = "copy", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                                 list( extend = "csv", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE)),
                                 list( extend = "excel", exportOptions = list(modifier = list(page = "all"),stripHtml = TRUE))
                               ),
                               pageLength=10, scrollX=TRUE))
  }, server = FALSE)
  
  
  output$enrich_tissue_plot <- renderPlot({
    res <- gtex_terms_res(); req(res)
    
    tab <- res$tab_tissue
    msg <- res$tissue_msg %||% "No tissues to plot."
    
    if (!is.data.frame(tab) || !nrow(tab)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 4) +
          ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1)
      )
    }
    
    if (!"label" %in% names(tab)) tab$label <- tab$term
    
    top <- tab %>%
      dplyr::filter(is.finite(p_adj), is.finite(a)) %>%
      dplyr::arrange(p_adj, p) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::mutate(
        neglogFDR = -log10(pmax(as.numeric(p_adj), 1e-300)),
        term_short = short_label(dplyr::coalesce(label, term), max_chars = 35, wrap_width = 14),
        term_short = factor(term_short, levels = rev(term_short))
      )
    
    ggplot2::ggplot(top, ggplot2::aes(x = term_short, y = a, fill = neglogFDR)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradientn(
        colours = c("yellow", "orange", "red"),
        values  = c(0, 0.5, 1),
        name    = "−log10(FDR)"
      ) +
      ggplot2::labs(
        title = "Top enriched tissues (GTEx)",
        x = NULL,
        y = "Count of GTEX hits by tissue"
      ) +
      ggplot2::theme_minimal(base_size = 13)
  })
  # ------------------------------------------------------------
  # Optional debug (runs only after click because depends on gtex_terms_res)
  # ------------------------------------------------------------
  output$enrich_gtex_debug <- renderPrint({
    res <- gtex_terms_res(); req(res)
    list(
      n_tissue_terms = nrow(res$tab_tissue),
      top_tissues    = head(res$tab_tissue$term, 10),
      tissue_msg     = res$tissue_msg
    )
  })
  
  ############################################################################
  ######## save candidates and hits
  ############################################################################
  
  norm_chr_int <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")] <- "23"
    x[x %in% c("Y","y")] <- "24"
    x[x %in% c("MT","Mt","mt","M","m")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pick_first_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  ################################################################################
  # ============================================================
  # GTEx → canonical CLUSTERS + CANDIDATES (shared by LD + ZIP)
  # ============================================================
  
  # ---- helpers (only if you don't already have them) ----
  
  norm_chr_int <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")] <- "23"
    x[x %in% c("Y","y")] <- "24"
    x[x %in% c("MT","Mt","mt","M","m")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pick_first_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  # ------------------------------------------------------------
  # A) CLUSTERS (canonical for LD module + export)
  #    returns: cluster_id, chr, cluster_start, cluster_end
  # ------------------------------------------------------------
  build_ld_clusters_from_gtex_app <- function() {
    
    # cl <- clusters_cur()
    cl <- clusters_stable()
    shiny::validate(shiny::need(is.data.frame(cl) && nrow(cl) > 0,
                                "No clusters available (run Step 2: build intervals + clusters)."))
    
    cl <- as.data.frame(cl)
    
    # tolerate variants
    if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
    
    shiny::validate(shiny::need(all(c("chr","start","end") %in% names(cl)),
                                "clusters_cur() must contain chr/start/end (or compatible names)."))
    
    # stable export id (DO NOT change GTEx internal naming)
    if ("cluster_chr_n" %in% names(cl) && any(nzchar(as.character(cl$cluster_chr_n)))) {
      cl$cluster_id <- as.character(cl$cluster_chr_n)     # e.g. chr1_3
    } else if ("cluster" %in% names(cl)) {
      cl$cluster_id <- paste0("cluster_", as.integer(cl$cluster))
    } else {
      cl$cluster_id <- paste0("cluster_", seq_len(nrow(cl)))
    }
    
    cl$chr   <- norm_chr_int(cl$chr)
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    shiny::validate(shiny::need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
    
    out <- cl[, c("cluster_id","chr","start","end")]
    names(out) <- c("cluster_id","chr","cluster_start","cluster_end")
    out
  }
  
  # ------------------------------------------------------------
  # B) CANDIDATES (canonical for LD module + export)
  #    returns: chr, pos_ini, pos_end, id_hit, classe
  # ------------------------------------------------------------
  build_ld_candidates_from_gtex_app <- function() {
    
    # ---------- GWAS candidates (hits_df) ----------
    h <- tryCatch(hits_df(), error = function(e) NULL)
    
    if (is.null(h) || !is.data.frame(h) || nrow(h) == 0) {
      gwas_cand <- data.frame(
        chr = integer(), pos_ini = integer(), pos_end = integer(),
        id_hit = character(), classe = character(),
        stringsAsFactors = FALSE
      )
    } else {
      shiny::validate(shiny::need(all(c("CHR","BP") %in% names(h)),
                                  "hits_df() must contain CHR and BP columns."))
      
      rsid <- if ("snp" %in% names(h)) as.character(h$snp) else NA_character_
      rsid[is.na(rsid) | !nzchar(rsid)] <- paste0("chr", h$CHR, ":", h$BP)
      
      gwas_cand <- data.frame(
        chr     = norm_chr_int(h$CHR),
        pos_ini = suppressWarnings(as.integer(h$BP)),
        pos_end = suppressWarnings(as.integer(h$BP)),
        id_hit  = rsid,
        classe  = "GWAS",
        stringsAsFactors = FALSE
      )
    }
    
    # ---------- GTEx candidates (gtex_data) ----------
    cd <- gtex_data()
    shiny::validate(shiny::need(is.data.frame(cd) && nrow(cd) > 0,
                                "No GTEx hits available (run Step 4: extract GTEx eQTLs)."))
    
    cd <- as.data.frame(cd)
    
    chr_col <- pick_first_col(cd, c("chr","CHR","chrom","chromosome"))
    shiny::validate(shiny::need(!is.null(chr_col), "gtex_data() has no chromosome column (chr/CHR/...)."))
    
    start_col <- pick_first_col(cd, c("start","pos","POS","BP","bp","position","variant_pos","pos_ini"))
    end_col   <- pick_first_col(cd, c("end","pos_end","bp_end"))
    shiny::validate(shiny::need(!is.null(start_col) || !is.null(end_col),
                                "gtex_data() has no position column (start/BP/POS/...)."))
    
    if (is.null(start_col)) start_col <- end_col
    if (is.null(end_col))   end_col   <- start_col
    
    id_col <- pick_first_col(cd, c("variant_id","variant","snp","rsid","RSID","ID","id_hit"))
    if (is.null(id_col)) {
      id_col <- ".id_tmp"
      cd[[id_col]] <- paste0("GTEx_chr", cd[[chr_col]], ":", cd[[start_col]], "-", cd[[end_col]])
    }
    
    gtex_cand <- data.frame(
      chr     = norm_chr_int(cd[[chr_col]]),
      pos_ini = suppressWarnings(as.integer(cd[[start_col]])),
      pos_end = suppressWarnings(as.integer(cd[[end_col]])),
      id_hit  = as.character(cd[[id_col]]),
      classe  = "GTEx_hit",
      stringsAsFactors = FALSE
    )
    
    # ---------- clean + combine ----------
    candidate_csv <- rbind(gwas_cand, gtex_cand)
    candidate_csv <- candidate_csv[
      is.finite(candidate_csv$chr) &
        is.finite(candidate_csv$pos_ini) &
        is.finite(candidate_csv$pos_end) &
        candidate_csv$pos_end >= candidate_csv$pos_ini &
        !is.na(candidate_csv$id_hit) & nzchar(candidate_csv$id_hit),
      , drop = FALSE
    ]
    
    shiny::validate(shiny::need(nrow(candidate_csv) > 0, "No candidates available for LD."))
    candidate_csv
  }
  
  # ============================================================
  # GTEx ZIP download (USES THE SAME BUILDERS)
  # ============================================================
  output$dl_candidates_zip <- downloadHandler(
    filename = function() {
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      tag_txt <- tryCatch({
        tg <- make_mode_thr_tag(input$cluster_method %||% "window", input$pthr, input$min_logp)
        paste0(tg$mode_tag, "_thr", tg$thr_txt)
      }, error = function(e) "run")
      
      paste0("gtex_candidates_", tag_txt, "_", stamp, ".zip")
    },
    
    content = function(file) {
      
      # Build canonical objects (single source of truth)
      cluster_csv   <- build_ld_clusters_from_gtex_app()
      candidate_csv <- build_ld_candidates_from_gtex_app()
      
      # Write to temp folder and zip
      tmpdir <- tempfile("gtex_export_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      f1 <- file.path(tmpdir, "cluster_gtex.csv")
      f2 <- file.path(tmpdir, "candidate_gtex.csv")
      
      utils::write.csv(cluster_csv,   f1, row.names = FALSE, quote = FALSE)
      utils::write.csv(candidate_csv, f2, row.names = FALSE, quote = FALSE)
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      setwd(tmpdir)
      
      utils::zip(zipfile = file, files = c("cluster_gtex.csv", "candidate_gtex.csv"))
    }
  )
  
  #  ##############################################################################
  #  ############################ LD Module   #####################################
  #  ################# GWAS -> LD module canonical cluster input###################
  
  ld_module_server(
    "ld",
    app_tag = "gtex",
    activate_r = reactive(TRUE)
  )
  
  ##############################################################################
  ##############################    RESET    ###################################
  ##############################################################################
  
  # -----------------------------
  # Reset session (GTEx Inspector)
  # -----------------------------
  observeEvent(input$reset_case, {
    # (A) Tanca modals si n'hi ha
    removeModal()
    
    # (B) Neteja seleccions de DT (si existeixen)
    for (id in c("hits_tbl", "cluster_dt", "gtex_table", "go_table", "kegg_table")) {
      try(DT::selectRows(DT::dataTableProxy(id), NULL), silent = TRUE)
    }
    
    # (C) Esborra fitxer temporal principal (selected_intervals.range)
    # workdir existeix al teu codi (per-session)
    if (exists("workdir", inherits = TRUE) && is.character(workdir) && dir.exists(workdir)) {
      f_rng <- file.path(workdir, "selected_intervals.range")
      if (file.exists(f_rng)) {
        try(unlink(f_rng, force = TRUE), silent = TRUE)
      }
    }
    
    # (D) Reset del "core state" real (reactiveVal del teu GTEx)
    intervals_raw(NULL)
    clusters_cur(NULL)
    ucsc_region(NULL)
    window_selected(NULL)
    
    # GTEx hits + caches
    gtex_hits_val(NULL)
    gene_map_cache_ens(NULL)
    gene_map_cache_map(NULL)
    
    # triggers d'enrichment
    go_trigger(0L)
    kegg_trigger(0L)
    
    # (E) Reset dels tracks UCSC (session$userData)
    # (els fas servir per construir els links/track hubs)
    session$userData$track_gwas_data <- tibble::tibble()
    session$userData$track_gtex_data <- tibble::tibble()
    
    # (F) Reset inputs a valors per defecte (segurs: si algun ID no existeix, no peta)
    safe <- function(expr) try(expr, silent = TRUE)
    
    safe(updateRadioButtons(session, "cluster_method", selected = "window"))
    safe(updateSliderInput(session, "pthr", value = 5))
    safe(updateSliderInput(session, "min_logp", value = 5))
    safe(updateNumericInput(session, "flank", value = 10000))
    safe(updateNumericInput(session, "min_hits", value = 3))
    safe(updateNumericInput(session, "win_bp", value = 1000000))
    safe(updateNumericInput(session, "step_bp", value = 250000))
    
    # (G) Feedback
    showNotification("Reset done. Ready for a new case.", type = "message", duration = 2)
  }, ignoreInit = TRUE)
  # ============================================================
  # ============================================================
  # ============================================================
  # GO SLIM (generic) enrichment — GTEx
  #  - gene set: ENTREZIDs (entrez_scope())
  #  - background: NULL (default) o dataset universe (universe_entrez_dataset())
  #  - GO annotations: org.Hs.eg.db (GOALL / ONTOLOGYALL)
  #  - Slim terms: goslim_generic.obo (GO slim generic)
  # ============================================================
  
  # --- 1) Sources required (put these files in your app under /R) ---
  #   R/goslim_utils.R        (load_goslim_generic_ids, map_go_to_goslim, shorten_go_label, etc.)
  #   R/go_slim_generic.R     (optional; loads GO_SLIM_GENERIC w/ TERM + ONTOLOGY from GO.db)
  #   R/goslim_generic.obo    (GO slim generic obo)
  #
  # If you already have them in a shared GItools folder, point source() there.
  
  observe({
    cl <- tryCatch(clusters_stable(), error = function(e) NULL)
    
    ids <- character(0)
    if (is.data.frame(cl) && nrow(cl)) {
      idcol <- if ("cluster_chr_n" %in% names(cl)) "cluster_chr_n" else
        if ("cluster_id" %in% names(cl)) "cluster_id" else names(cl)[1]
      ids <- unique(as.character(cl[[idcol]]))
    } else if (is.vector(cl) && length(cl)) {
      ids <- unique(as.character(cl))
    }
    
    ids <- sort(ids[!is.na(ids) & nzchar(ids)])
    updateSelectInput(session, "func_cluster_id", choices = ids, selected = ids[1] %||% "")
  })
  
  
  local({
    f1 <- file.path("R", "goslim_utils.R")
    f2 <- file.path("R", "go_slim_generic.R")
    f3 <- file.path("R", "goslim_generic.obo")
    
    if (file.exists(f1)) source(f1, local = TRUE)
    if (file.exists(f2)) source(f2, local = TRUE)
    
    # If GO_SLIM_GENERIC is not defined by go_slim_generic.R, define it minimally here:
    if (!exists("GO_SLIM_GENERIC", inherits = TRUE) && file.exists(f3) &&
        exists("load_goslim_generic_ids", inherits = TRUE)) {
      
      # Build a simple GOID->ONTOLOGY table from the IDs in the OBO
      ids <- load_goslim_generic_ids(f3)
      GO_SLIM_GENERIC <- AnnotationDbi::select(
        GO.db::GO.db,
        keys    = ids,
        keytype = "GOID",
        columns = c("GOID", "TERM", "ONTOLOGY")
      ) %>%
        dplyr::filter(!is.na(GOID), !is.na(TERM), !is.na(ONTOLOGY)) %>%
        dplyr::transmute(GOID, ONTOLOGY, slim_term = TERM) %>%
        dplyr::distinct(GOID, ONTOLOGY, slim_term)
      
      assign("GO_SLIM_GENERIC", GO_SLIM_GENERIC, envir = .GlobalEnv)
    }
  })
  
  # ------------------------------------------------------------
  # Trigger (com GO/KEGG)
  # ------------------------------------------------------------
  goslim_trigger <- reactiveVal(NULL)
  
  observeEvent(input$run_enrich, {
    tab <- input$enrich_tabs %||% "tab_enrich_go"
    if (identical(tab, "tab_enrich_goslim")) {
      goslim_trigger((goslim_trigger() %||% 0L) + 1L)
    }
  }, ignoreInit = TRUE)
  
  # ------------------------------------------------------------
  # Helpers: build GO slim TERM2GENE from org.Hs.eg.db annotations
  # ------------------------------------------------------------
  
  # Slim IDs per ontologia
  goslim_ids_by_onto <- reactive({
    validate(need(exists("GO_SLIM_GENERIC", inherits = TRUE), "GO_SLIM_GENERIC not loaded."))
    df <- get("GO_SLIM_GENERIC", inherits = TRUE)
    validate(need(is.data.frame(df) && nrow(df) > 0, "GO_SLIM_GENERIC is empty."))
    split(unique(df$GOID), as.character(df$ONTOLOGY))
  })
  
 

  # ------------------------------------------------------------
  # GO SLIM enrichment (multi-ontology) — enricher() + filtrem després per FDR
  # ------------------------------------------------------------
  
  
  
  goslim_orgdb_cache <- reactiveVal(list())
  
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
  
  goslim_enrich_raw <- eventReactive(goslim_trigger(), {
    validate(need(requireNamespace("clusterProfiler", quietly = TRUE), "clusterProfiler package is required."))
    validate(need(requireNamespace("dplyr", quietly = TRUE), "dplyr package is required."))
    
    withProgress(message = "Running GO slim enrichment…", value = 0, {
      
      gene <- entrez_scope()
      gene <- unique(as.character(gene))
      gene <- gene[!is.na(gene) & nzchar(gene)]
      validate(need(length(gene) > 0, "No valid genes for the selected scope."))
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      ontos <- input$go_ontos %||% c("BP", "CC", "MF")
      ontos <- intersect(c("BP", "CC", "MF"), ontos)
      validate(need(length(ontos) > 0, "Select at least one ontology (BP/CC/MF)."))
      
      minGS_in <- input$enrich_min_gs %||% 10
      maxGS    <- input$enrich_max_gs %||% 500
      
      gs_all <- goslim_orgdb_sets()
      incProgress(0.1, detail = "Reference ready")
      
      res <- purrr::map_dfr(ontos, function(ont) {
        incProgress(0.6 / length(ontos), detail = paste("Ontology:", ont))
        
        gs <- gs_all[[ont]]
        universe <- unique(gs$universe)
        gene_use <- intersect(unique(gene), universe)
        
        if (!length(gene_use)) return(NULL)
        
        minGS <- max(1L, min(as.integer(minGS_in), length(gene_use)))
        
        eg <- tryCatch({
          suppressMessages(
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
        }, error = function(e) NULL)
        
        if (is.null(eg) || is.null(eg@result) || !nrow(eg@result)) return(NULL)
        
        df <- as.data.frame(eg@result, stringsAsFactors = FALSE)
        if (!nrow(df)) return(NULL)
        df$Ontology <- ont
        df
      })
      
      if (is.null(res) || !nrow(res)) return(tibble::tibble())
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
  
  # ------------------------------------------------------------
  # Outputs: table + barplot (copiat del GO, canviant labels)
  # ------------------------------------------------------------
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
  
  
  output$goslim_bar <- plotly::renderPlotly({
    
    dat <- goslim_enrich_tbl()
    validate(need(is.data.frame(dat) && nrow(dat) > 0, "GoSlim: no enrichment data to plot yet. Run enrichment first."))
    
    top_n <- input$go_topn %||% 15
    gap   <- 2
    
    validate(need(all(c("Description","Count","Ontology","p.adjust","pvalue","GeneRatio","BgRatio") %in% names(dat)),
                  "goslim_enrich_tbl() must return: Description, Count, Ontology, p.adjust, pvalue, GeneRatio, BgRatio"))
    
    if (!exists("shorten_go_label", mode = "function")) {
      shorten_go_label <- function(x, max = 55) {
        x <- as.character(x)
        x <- stringr::str_squish(x)
        ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
      }
    }
    
    dat <- dat %>%
      dplyr::mutate(
        Ontology   = as.character(Ontology),
        slim_term  = as.character(Description),
        n_genes    = suppressWarnings(as.integer(Count)),
        p.adjust   = suppressWarnings(as.numeric(p.adjust)),
        pvalue     = suppressWarnings(as.numeric(pvalue)),
        term_short = stringr::str_wrap(shorten_go_label(slim_term, 55), width = 28)
      ) %>%
      dplyr::filter(Ontology %in% c("BP","CC","MF")) %>%
      dplyr::filter(is.finite(p.adjust), !is.na(n_genes), nzchar(slim_term))
    
    validate(need(nrow(dat) > 0, "GoSlim: no BP/CC/MF terms to plot."))
    
    ont_order <- c("BP","CC","MF")
    present   <- intersect(ont_order, unique(dat$Ontology))
    validate(need(length(present) > 0, "GoSlim: no BP/CC/MF data to plot."))
    
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
    if (all(c("BP","CC") %in% present)) vlines <- c(vlines, offsets["CC"] - gap/2)
    if (all(c("CC","MF") %in% present)) vlines <- c(vlines, offsets["MF"] - gap/2)
    
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
    
    ttl <- paste0("GO slim enrichment by gene counts — ", scope_label())
    
    dat$tooltip <- paste0(
      "<b>", dat$Ontology, "</b>",
      "<br><b>GO slim term:</b> ", htmltools::htmlEscape(dat$slim_term),
      "<br><b>Count:</b> ", dat$n_genes,
      "<br><b>GeneRatio:</b> ", dat$GeneRatio,
      "<br><b>BgRatio:</b> ", dat$BgRatio,
      "<br><b>p:</b> ", sprintf("%.3g", dat$pvalue),
      "<br><b>FDR:</b> ", sprintf("%.3g", dat$p.adjust)
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
        title = ttl
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 11, face = "bold"),
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
  
  
  output$goslim_debug <- renderPrint({
    tab <- input$enrich_tabs %||% NA_character_
    onto <- (input$go_ontos %||% c("BP","CC","MF"))[1]
    
    gene <- tryCatch(entrez_scope(), error = function(e) character(0))
    uni  <- tryCatch(universe_entrez_for_scope(), error = function(e) NULL)
    
    # TERM2GENE per una ontologia (BP per defecte)
    t2g <- tryCatch(goslim_term2gene_for_onto(onto, universe_entrez = uni), error = function(e) NULL)
    
    list(
      tab = tab,
      scope = input$func_scope %||% NA_character_,
      cluster = input$func_cluster_id %||% NA_character_,
      n_entrez_scope = length(gene),
      universe_is_null = is.null(uni),
      n_universe = if (is.null(uni)) NA_integer_ else length(uni),
      onto_tested = onto,
      t2g_is_null = is.null(t2g),
      n_term2gene = if (is.data.frame(t2g)) nrow(t2g) else NA_integer_,
      n_terms = if (is.data.frame(t2g)) length(unique(t2g$term)) else NA_integer_
    )
  })
  ##≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠

  ##≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  
}

shinyApp(ui, server)