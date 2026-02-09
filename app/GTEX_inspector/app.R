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

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# --- Fixar dplyr com a select/mutate/arrange oficial ---
select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
summarise  <- dplyr::summarise
rename     <- dplyr::rename

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
BASE <- if (exists("gi_cfg", mode = "function")) gi_cfg()$root else gi_root_guess
BASE <- normalizePath(BASE, winslash = "/", mustWork = TRUE)

# 4) Load GI shared state (defines gi_shared_root / gi_state_root)
state_file <- file.path(BASE, "_shared", "gi_state.R")
stopifnot(file.exists(state_file))
source(state_file, local = TRUE)

cat("[GI] using gi_state.R at:", state_file, "\n")
cat("[GI] gi_shared_root =", gi_shared_root, "\n")
# cat("[GI] gi_state_root  =", gi_state_root, "\n")
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

# ===========================
# UI  (FULL)
# ===========================

ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("🧠 GTEx eQTL Inspector")
  ),
  id = "topnav",
  
  # -----------------------------
  # Global light-grey panels for plots/tables
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
  "))
  ),
  
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    
    sidebarLayout(
      
      sidebarPanel(
        width = 3,
        
        # -----------------------------
        # Mode new/old
        # -----------------------------
        h4("Select new/old analysis", style = "color:#1A4E8A; font-size:22px; font-weight:700;"),
        radioButtons(
          "use_preloaded_catalog",
          label = NULL,
          choices = c(
            "🧪 Perform a new analysis (run steps 2–4)" = "new",
            "🧾 Analyze a precomputed hits file (skip steps 2–4)" = "old"
          ),
          selected = "new"
        ),
        tags$hr(),
        
        # -----------------------------
        # OLD MODE: load precomputed hits
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_catalog == 'old'",
          fileInput(
            "catalog_preload",
            "Load precomputed hits (CSV/RDS/TSV). Must include cluster_id.",
            accept = c(".csv", ".rds", ".tsv")
          ),
          tags$hr()
        ),
        
        # -----------------------------
        # Step 1: GWAS input (always visible)
        # -----------------------------
        h3(
          "Step 1 · Load GWAS hits",
          style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
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
        # NEW MODE: steps 2–4
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_catalog == 'new'",
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
              "➊ Generate intervals + clusters",
              style = "background-color: #ffdd57; color: black; font-weight: bold;"
            ),
            h4("Preview clusters (derived from intervals)"),
            div(class = "panel-lite", verbatimTextOutput("ranges_preview")),
            tags$hr(),
            
            # ---- Step 3
        #    h3(
        #      "Step 3 · Load GTEx eQTL resource",
        #      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"
        #    ),
        #    helpText("By default the app tries to load 'www/gtex_eqtl_light.rds' at startup."),
        #    textOutput("gtex_status"),
        #    tags$hr(),
            
            # ---- Step 3
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
          )
        )
      ),
      
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            h4("Top: GWAS p-values · Bottom: GTEx eQTLs inside all clusters"),
            div(class = "panel-lite", plotlyOutput("manhattan_combo", height = "700px")),
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
                div(class = "panel-lite", plotOutput("tissue_family_bar", height = "350px"))
              ),
              column(
                6,
                h4("Advanced dotplot (tissue × -log10(p))"),
                div(class = "panel-lite", plotOutput("tissue_dotplot", height = "350px"))
              )
            ),
            tags$hr(),
            fluidRow(
              column(
                7,
                h4("Functional heatmap: genes × tissue families"),
                div(class = "panel-lite", plotOutput("tissue_heatmap", height = "450px"))
              ),
              column(
                5,
                h4("Automatic interpretation"),
                div(class = "panel-lite", verbatimTextOutput("gtex_interpret"))
              )
            )
          ),
          
          # ============================================================
          # UI — GTEx Enrichment (Catalog/NonSyn-style)
          # ============================================================
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧬 Enrichment</span>"),
            value = "tab_enrich",
            
            sidebarLayout(
              sidebarPanel(
                width = 3,
                
                radioButtons(
                  "func_scope", "Scope",
                  choices  = c("Global" = "global", "Cluster" = "cluster"),
                  selected = "global"
                ),
                
                conditionalPanel(
                  condition = "input.func_scope == 'cluster'",
                  uiOutput("func_cluster_ui")
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
                    title = HTML("<span style='font-size:15px; font-weight:600;'>GO</span>"),
                    value = "tab_enrich_go",
                    div(class = "panel-lite",
                        shinycssloaders::withSpinner(plotOutput("go_bar", height = 400))
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
                        shinycssloaders::withSpinner(plotOutput("kegg_bar", height = 400))
                    ),
                    tags$hr(),
                    div(class = "panel-lite",
                        shinycssloaders::withSpinner(DT::DTOutput("kegg_table"))
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

server <- function(input, output, session){
  
# # ============================================================
# # CANONICAL shared objects (always exist)
# # ============================================================
# # These are used across the app (UI previews, plots, tables, LD, etc.)
# gwas_shared_rv <- shiny::reactiveVal(tibble::tibble())   # standardized GWAS table (may come from HUB/master)
# clusters_cur   <- shiny::reactiveVal(NULL)   # standardized clusters table (may come from HUB/master OR local clustering)
# intervals_raw  <- shiny::reactiveVal(NULL)   # canonical intervals derived from clusters (chr/start/end/label)
# 
# # Ensure 'rv' exists for legacy parts that still expect it
# if (!exists("rv", inherits = FALSE)) rv <- shiny::reactiveValues()
# 
# # ============================================================
# # 0) deeplink (after clusters_cur exists)
# # ============================================================
# source(file.path(gi_shared_root, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
# gitools_deeplink_gtex(session, clusters_cur)
# 
# source(file.path(gi_shared_root, "gi_slave_canonical.R"), local = TRUE)
# gi_sync <- gi_slave_canonical_init(session)
# 
# source(file.path(gi_shared_root, "gi_clusters_canonical.R"), local = TRUE)

  SHARED <- gi_shared_root  # definit a gi_state.R
  
  # -----------------------------
  # Shared reactive holders (Hub OR local)
  # -----------------------------
  gwas_shared_rv <- shiny::reactiveVal(NULL)
  
  # -----------------------------
  # Canonical cluster engine (ALWAYS available)
  # -----------------------------
  source(file.path(SHARED, "gi_clusters_canonical.R"), local = TRUE)
  
  # Proxy: prioritat al gwas_df() local si existeix; sinó al Hub
  gwas_df_proxy <- shiny::reactive({
    if (exists("gwas_df", inherits = TRUE) && is.function(gwas_df)) {
      df2 <- tryCatch(gwas_df(), error = function(e) NULL)
      if (is.data.frame(df2) && nrow(df2) > 0) return(df2)
    }
    df <- gwas_shared_rv()
    if (is.data.frame(df) && nrow(df) > 0) return(df)
    NULL
  })
  
  gi_cl <- gi_clusters_canonical_init(
    session, input, output,
    gwas_df = gwas_df_proxy,
    build_btn_id   = "build_ranges",
    clusters_dt_id = "cluster_dt",
    hits_rows_id   = "hits_tbl_rows_selected",
    app_count_col  = "n_gtex"   # <-- usa el nom real de la columna de recompte a GTEx
  )
  
  # Expose canonical names (la resta de l’app ha d’usar aquests)
  clusters_cur  <- gi_cl$clusters_cur
  intervals_raw <- gi_cl$intervals_raw
  
  # -----------------------------
  # Slave canonical sync (Hub -> this app)
  # -----------------------------
  source(file.path(SHARED, "gi_slave_canonical.R"), local = TRUE)
  gi_sync <- tryCatch(gi_slave_canonical_init(session), error = function(e) NULL)
  
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
  
  # tolerància: si ve del master pot tenir Pval o p-value original
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
# Canonical clusters engine (local build_ranges)
# - provides: intervals_raw(), clusters_cur(), selected_cluster()
# ------------------------------------------------------------
gi_clust <- gi_clusters_canonical_init(
  session = session, input = input, output = output,
  gwas_df = gwas_df,
  build_btn_id   = "build_ranges",
  clusters_dt_id = "cluster_dt",
  hits_rows_id   = "hits_tbl_rows_selected",
  app_count_col  = "n_gtex"
)

intervals_raw    <- gi_clust$intervals_raw
clusters_cur     <- gi_clust$clusters_cur
selected_cluster <- gi_clust$selected_cluster


# ---- Preview text (Step 2) ----
output$ranges_preview <- renderPrint({
  cl <- clusters_cur()
  if (!is.data.frame(cl) || nrow(cl) == 0) {
    cat("No clusters yet. Click: ➊ Generate intervals + clusters\n")
    return(invisible())
  }
  
  cat("Clusters generated:", nrow(cl), "\n\n")
  print(utils::head(cl, 10))
})

# ---- Clusters DT (Manhattan tab) ----
output$cluster_dt <- DT::renderDT({
  cl <- clusters_cur()
  
  if (!is.data.frame(cl) || nrow(cl) == 0) {
    return(DT::datatable(
      data.frame(Message = "No clusters yet. Click: ➊ Generate intervals + clusters"),
      options = list(dom = "t"),
      rownames = FALSE
    ))
  }
  
  # ordre i columnes “amables”
  show_cols <- c(
    "cluster_chr_n","cluster_id","chr","start","end",
    "n_snps","top_snp","top_logp",
    "n_gtex"
  )
  keep <- intersect(show_cols, names(cl))
  cl2 <- cl[, keep, drop = FALSE]
  
  DT::datatable(
    cl2,
    rownames  = FALSE,
    extensions = "Buttons",
    options    = list(
      dom        = "Bfrtip",
      buttons    = c("copy", "csv", "excel", "pdf", "print"),
      pageLength = 10,
      scrollX    = TRUE
    )
  ) %>%
    formatRound(columns = c("top_logp"), digits = 2)
})

# ------------------------------------------------------------
# Apply MASTER sync -> overwrite GWAS/clusters into canonical state
# (keeps local build_ranges available too)
# ------------------------------------------------------------
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


# Cluster assign from hub
observeEvent(gi_sync$clusters_shared(), {
  # avís curt i clar
  id <- paste0("hub_sync_", as.integer(Sys.time()))
  showNotification("🔄 HUB: sincronitzant CLUSTERS…", type = "message", duration = NULL, id = id)
  
  cl <- gi_sync$clusters_shared()
  req(is.data.frame(cl), nrow(cl) > 0)
  
  cl <- standardize_cluster_ids(cl)
  clusters_cur(cl)
  
  # intervals per passos downstream
  intervals_raw(cl |> dplyr::transmute(
    chr = chr, start = start, end = end,
    label = dplyr::coalesce(as.character(cluster_chr_n), as.character(cluster_id))
  ))
  
  cat("[SLAVE] applied CLUSTERS rows=", nrow(cl), "\n")
  
  removeNotification(id)
  showNotification(sprintf("✅ HUB: clusters sincronitzats (%d).", nrow(cl)),
                   type = "message", duration = 2)
  
}, ignoreInit = FALSE)





###########################################

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

dfp_manhattan <- reactive({
  df <- gwas_df(); req(nrow(df) > 0)
  # filter pval < 0.0 to a easely plot
  df <- df %>% dplyr::filter(Pval < 0.05)
  
  ref <- .ref_hg38
  df$chrN <- norm_chr_generic(df$CHR)
  
  df %>%
    inner_join(ref %>% select(chr, chr_cum), by = c("chrN" = "chr")) %>%
    arrange(chrN, BP) %>%
    mutate(
      BPcum = BP + chr_cum,
      colgrp = CHR %% 2
    )
})

axis_df <- reactive({
  dfp <- dfp_manhattan()
  dfp %>%
    group_by(chrN) %>%
    summarise(center = mean(BPcum, na.rm = TRUE), .groups = "drop")
})

# ---------- Hits (window mode) ----------
hits_df <- reactive({
  df <- gwas_df()
  req(is.data.frame(df), nrow(df) > 0)
  req(input$pthr)
  
  # garantir rsid (si no existeix, usa snp)
  if (!"rsid" %in% names(df)) df$rsid <- NA_character_
  if (!"snp"  %in% names(df)) df$snp  <- NA_character_
  df$rsid <- as.character(df$rsid)
  df$snp  <- as.character(df$snp)
  df$rsid <- dplyr::coalesce(df$rsid, df$snp)
  
  df %>%
    dplyr::filter(.data$logp >= as.numeric(input$pthr)) %>%
    dplyr::arrange(dplyr::desc(.data$logp)) %>%
    dplyr::select(
      dplyr::all_of(c("CHR", "BP", "snp")),
      dplyr::any_of("rsid"),
      p = Pval,
      logp
    )
})


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
      buttons    = c("copy", "csv", "excel", "pdf", "print"),
      pageLength = 10,
      scrollX    = TRUE
    )
  )
})

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


observeEvent(input$run_gtex, {
  
  req(gtex_all())
  
  cl <- clusters_cur()
  validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters. Generate intervals + clusters first."))
  
  # ------------------------------------------------------------
  # 1) Normalitza GTEx chr/pos perquè el filtre dins extract funcioni
  #    (la teva funció extract_gtex_from_clusters compara `chr` amb chr_label_plink())
  #    - gtex_all() pot venir com "chr1" -> passa a "1"
  #    - "chrX" -> "X" (si chr_label_plink retorna "X")
  # ------------------------------------------------------------
  gtex_norm <- gtex_all() %>%
    dplyr::mutate(
      chr = as.character(chr),
      chr = gsub("^CHR", "chr", chr),
      chr = gsub("^chr", "", chr),          # "1","X","23"...
      chr = toupper(chr),
      chr = dplyr::if_else(chr == "23", "X", chr),
      chr = dplyr::if_else(chr == "24", "Y", chr),
      pos = suppressWarnings(as.integer(pos))
    ) %>%
    dplyr::filter(!is.na(chr), !is.na(pos))
  
  # ------------------------------------------------------------
  # 2) Extract: usa la teva funció, però ara ja li passen chr coherent
  # ------------------------------------------------------------
  hits <- tryCatch(
    extract_gtex_from_clusters(cl, gtex_norm),
    error = function(e) {
      cat("[GTEx][run_gtex] extract_gtex_from_clusters ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(hits) || !is.data.frame(hits)) {
    hits <- gtex_norm[0, , drop = FALSE]
  }
  
  gtex_hits_val(hits)
  
  # ------------------------------------------------------------
  # 3) Actualitza n_gtex (robust)
  #    IMPORTANT: no comptis per `cluster` perquè pot fallar quan:
  #      - hits buit
  #      - al HUB els noms venen com "cluster_chr1_1" vs "chr1_1"
  #    Fem servir cluster_chr_n quan existeix, i normalitzem prefix.
  # ------------------------------------------------------------
  cl_base <- cl %>% dplyr::select(-dplyr::any_of("n_gtex"))
  
  if (nrow(hits) > 0) {
    
    # assegura columna clau:
    # preferim cluster_chr_n si existeix; si no, fem fallback a cluster (numeric)
    if ("cluster_chr_n" %in% names(hits)) {
      
      hits2 <- hits %>%
        dplyr::mutate(
          cluster_chr_n = as.character(cluster_chr_n),
          cluster_chr_n = sub("^cluster_", "", cluster_chr_n)  # HUB fix
        )
      
      # cl també normalitzat
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
      
      # retorna a l'estructura original (no forces canvis de columnes)
      # però mantén n_gtex
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
      
      # no tenim cap clau fiable -> posa 0 i no trenquis res
      clusters_cur(cl_base %>% dplyr::mutate(n_gtex = 0L))
    }
    
  } else {
    clusters_cur(cl_base %>% dplyr::mutate(n_gtex = 0L))
  }
  
  output$gtex_summary <- renderText({
    sprintf("GTEx eQTLs inside all clusters: %d hits (clusters: %d)",
            nrow(hits), nrow(cl))
  })
  
  if (!nrow(hits)) {
    showNotification("No GTEx eQTLs found inside the current clusters.", type = "warning")
  } else {
    showNotification("GTEx eQTLs extracted successfully.", type = "message")
  }
  
}, ignoreInit = TRUE)



###  
gtex_data <- reactive({
  gtex_hits_val()
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
    options = list(pageLength = 15, scrollX = TRUE)
  )
})


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

output$tissue_family_bar <- renderPlot({
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

output$tissue_dotplot <- renderPlot({
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

output$tissue_heatmap <- renderPlot({
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
  
  # --- ref_map (chr_num -> chr_cum) per poder projectar coordenades a BPcum ---
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
  
  # Show threshold line only in "new" mode (or any mode NOT "old")
  if (!identical(input$use_preloaded_catalog, "old")) {
    p1 <- p1 + ggplot2::geom_hline(yintercept = thr_y, linetype = "dashed")
  }
  
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
  
  cat("[DBG] clusters_cur n=", if (is.data.frame(clusters_cur())) nrow(clusters_cur()) else -1, "\n")
  
  cl <- clusters_cur()
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
  out <- plotly::event_register(out, "plotly_click")
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

# --- helper: extreu rang X d'un relayout (subplot: xaxis o xaxis2) ---
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
  req(nrow(dfp) > 0, nrow(df) > 0)
  
  dfp2 <- dfp %>% dplyr::mutate(CHR = as.integer(CHR), BP = as.integer(BP))
  df2  <- df  %>% dplyr::mutate(CHR = as.integer(CHR), BP = as.integer(BP)) %>%
    dplyr::select(CHR, BP, rsid, snp)
  
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
  cl <- clusters_cur()
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
  cl <- clusters_cur()
  req(is.data.frame(cl), nrow(cl) > 0, input$func_chr)
  
  x <- cl %>%
    dplyr::filter(as.integer(chr) == as.integer(input$func_chr)) %>%
    dplyr::arrange(start, end)
  
  validate(need(nrow(x) > 0, "No clusters on this chromosome."))
  
  lab <- paste0(x$cluster_chr_n, " (", x$start, "-", x$end, ")")
  choices <- stats::setNames(x$cluster_chr_n, lab)
  
  selectInput("func_cluster_id", "Cluster", choices = choices, selected = x$cluster_chr_n[1])
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
  
  bg_lbl <- if (identical(bg, "dataset")) "Dataset genes" else "Default background"
  
  msg <- if (identical(tab, "tab_enrich_go")) {
    if (identical(bg, "orgdb")) {
      paste0("<b>Background:</b> ", bg_lbl, "<br>",
             "GO enrichment uses the default <i>OrgDb</i>-based background (org.Hs.eg.db).")
    } else {
      paste0("<b>Background:</b> ", bg_lbl, "<br>",
             "GO enrichment uses only genes present in the GTEx hits dataset as background (cluster mode).")
    }
  } else {
    if (identical(bg, "orgdb")) {
      paste0("<b>Background:</b> ", bg_lbl, "<br>",
             "KEGG enrichment uses the default KEGG background (organism-specific).")
    } else {
      paste0("<b>Background:</b> ", bg_lbl, "<br>",
             "KEGG enrichment uses dataset genes as background (cluster mode).")
    }
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
  
  DT::datatable(out, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
})

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

output$go_bar <- renderPlot({
  dat <- go_top_df()
  req(is.data.frame(dat), nrow(dat) > 0)
  
  fb   <- isTRUE(attr(dat, "fallback"))
  pcut <- input$enrich_pcut %||% 0.05
  
  ttl <- paste0("GO enrichment (", scope_label(), ")")
  subttl <- if (fb) {
    paste0("No terms at FDR ≤ ", pcut, " — showing top terms ranked by FDR; bar height = -log10(pvalue)")
  } else {
    paste0("FDR ≤ ", pcut, " (bar height = -log10(pvalue))")
  }
  
  dat <- dat %>% dplyr::mutate(Ontology = factor(as.character(Ontology), levels = c("BP","CC","MF")))
  cols_go <- c(BP = "darkgreen", CC = "orange", MF = "darkblue")
  
  mk <- function(onto) {
    d <- dat %>% dplyr::filter(Ontology == onto)
    if (!nrow(d)) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = paste0(onto, " (none)")))
    
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
  
  DT::datatable(out, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
})

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

output$kegg_bar <- renderPlot({
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
      padj  = suppressWarnings(as.numeric(p.adjust))
    )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = stats::reorder(term_short, score),
      y    = score,
      fill = padj
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = ttl, x = NULL, y = "-log10(pvalue)", fill = "FDR") +
    ggplot2::scale_fill_gradient(
      low  = "red",   # més significatiu (FDR petit)
      high = "orange",   # menys significatiu
      na.value = "yellow",
      trans = "reverse"   # perquè low (color fort) correspongui a FDR petit
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}, height = function() {
  n <- tryCatch(nrow(kegg_top_df()), error = function(e) 0)
  max(420, 28 * n)
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
  
  cl <- clusters_cur()
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

##############################################################################
############################ LD Module  GTEX #################################
##############################################################################
# -----------------------------
# LD PANEL
# -----------------------------
clusters_r <- reactive({
  build_ld_clusters_from_gtex_app()
})

candidates_r <- reactive({
  build_ld_candidates_from_gtex_app()
})

ld_module_server(
  "ld",
  clusters_r   = clusters_r,
  candidates_r = candidates_r
)

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
  
  safe(updateRadioButtons(session, "use_preloaded_catalog", selected = "new"))
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



}

shinyApp(ui, server)