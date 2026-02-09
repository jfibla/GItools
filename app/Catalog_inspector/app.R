## app.R — GWAS Catalog Inspector (ADAPTED TO NonSyn: new/old + canonical clusters + Manhattan_combo + GO/KEGG + LD)
#~/Library/CloudStorage/Dropbox/Catalog_inspector/app.R
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
cl_file <- file.path(cfg$shared, "gi_clusters_canonical.R")

gi_state_file <- file.path(cfg$shared, "gi_state.R")
stopifnot(file.exists(gi_state_file))
source(gi_state_file, local = TRUE)

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



# =============================================================================
# UI
# =============================================================================
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
  
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis</span>"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        # -----------------------------
        # Mode new/old (NonSyn-style)
        # -----------------------------
        h4("Select new/old analysis", style = "color:#1A4E8A; font-size:22px; font-weight:700;"),
        radioButtons(
          "use_preloaded_catalog",
          label = NULL,
          choices = c(
            "🧪 Perform a new analysis (run steps 2–4)" = "new",
            "🧾 Analyze a precomputed Catalog hits file (skip steps 2–4)" = "old"
          ),
          selected = "new"
        ),
        tags$hr(),
        
        # -----------------------------
        # Step 1: GWAS input
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
        
        fileInput("gwas_file", "GWAS p-value table (TSV/CSV)", accept = c(".tsv",".txt",".csv")),
        checkboxInput("gwas_header", "First row contains column names", TRUE),
        radioButtons("gwas_sep", "Separator", c("Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"), selected = "\t"),
        uiOutput("gwas_p_selector"),
        tags$hr(),
        
        # -----------------------------
        # OLD MODE: load precomputed Catalog hits
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_catalog == 'old'",
          fileInput(
            "catalog_preload",
            "Load precomputed Catalog hits (CSV or RDS). Must include cluster_id.",
            accept = c(".csv", ".rds", ".tsv")
          )
        ),
        
        # -----------------------------
        # NEW MODE: steps 2–4
        # -----------------------------
        conditionalPanel(
          condition = "input.use_preloaded_catalog == 'new'",
          tagList(
            h3("Step 2 · Clustering GWAS hits", style="color:#1A4E8A; font-size:22px; font-weight:700;"),
            actionButton("info_01", "ℹ️ Clustering method", class = "btn btn-default"),
            
            radioButtons(
              "cluster_method", "Clustering method:",
              choices  = c(
                "By hit intervals (thr + flank → merge)" = "window",
                "By hit density (min_logp + min_hits)"   = "hits"
              ),
              selected = "window"
            ),
            
            # -----------------------------
            # Method: window
            # -----------------------------
            conditionalPanel(
              condition = "input.cluster_method == 'window'",
              sliderInput("pthr", "-log10(P) threshold", min = 2, max = 20, value = 5, step = 0.5),
              numericInput("flank", "Flank (+/- bp)", value = 10000, min = 0, max = 10000000, step = 1000)
            ),
            
            # -----------------------------
            # Method: hits (3 submodes)
            # -----------------------------
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
              "➊ Generate intervals → merge → clusters",
              style = "background-color: #ffdd57; color: black; font-weight: bold;"
            ),
            h4("Preview selected_intervals.range (clusters)"),
            div(class = "panel-lite", verbatimTextOutput("ranges_preview")),
            tags$hr(),
 
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
            
          )
        )
      ),
      
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            h4("Top: GWAS p-values · Bottom: GWAS Catalog hits inside clusters"),
            div(class = "panel-lite", withSpinner(plotlyOutput("manhattan_combo", height = "700px"))),
            tags$hr(),
            
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
                  inline   = TRUE
                ),
                conditionalPanel(
                  condition = "input.catalog_scope == 'pick'",
                  uiOutput("catalog_cluster_pick_ui")
                )
              ),
              
              column(
                width = 6,
                uiOutput("catalog_gene_pick_ui")
              )
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
          ),
          
          ###################################
          # DISEASE/TRAIT Enrichment
          ##################################
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Disease/Trait Enrichment</span>"),
            br(),
            h3("Enrichment of GWAS Catalog terms", style="color:#1A4E8A; font-weight:700; margin-top:10px;"),
            helpText("ORA of DISEASE / MAPPED_TRAIT terms among Catalog hits."),
            
            radioButtons(
              "dt_scope", "Scope",
              choices = c(
                "Global (all Catalog hits) — not meaningful" = "global",
                "Selected cluster (from Clusters panel)"     = "selected",
                "Pick a cluster here"                        = "pick"
              ),
              selected = "global",
              inline = TRUE
            ),
            
            conditionalPanel(
              condition = "input.dt_scope == 'pick'",
              uiOutput("dt_cluster_ui")
            ),
            
            uiOutput("dt_scope_note"),
            br(),
            
            actionButton("run_dt_enrich", "Run enrichment", icon = icon("play")),
            br(), br(),
            
            tabsetPanel(
              tabPanel(
                HTML("<b>DISEASE Enrichment</b>"),
                br(),
                fluidRow(
                  column(6, div(class = "panel-lite", withSpinner(DTOutput("enrich_disease_table")))),
                  column(6, div(class = "panel-lite", withSpinner(plotOutput("enrich_disease_plot", height = "450px"))))
                )
              ),
              tabPanel(
                HTML("<b>MAPPED_TRAIT Enrichment</b>"),
                br(),
                fluidRow(
                  column(6, div(class = "panel-lite", withSpinner(DTOutput("enrich_trait_table")))),
                  column(6, div(class = "panel-lite", withSpinner(plotOutput("enrich_trait_plot", height = "450px"))))
                )
              )
            )
          ),
          
          # -----------------------------
          # GO/KEGG enrichment (NonSyn-style)
          # -----------------------------
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧬 GO/KEGG Enrichment</span>"),
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
                  choices  = c("Reference annotated genes" = "orgdb", "Dataset genes" = "dataset"),
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
                    div(class = "panel-lite", withSpinner(plotOutput("go_bar", height = 400))),
                    tags$hr(),
                    div(class = "panel-lite", withSpinner(DT::DTOutput("go_table")))
                  ),
                  
                  tabPanel(
                    title = HTML("<span style='font-size:15px; font-weight:600;'>KEGG</span>"),
                    value = "tab_enrich_kegg",
                    div(class = "panel-lite", withSpinner(plotOutput("kegg_bar", height = 400))),
                    tags$hr(),
                    div(class = "panel-lite", withSpinner(DT::DTOutput("kegg_table")))
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


# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  #############################################################
  # Catalog Inspector (inside server): shared helpers (portable)
  source(file.path(cfg$shared, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
  source(file.path(cfg$shared, "gi_clusters_canonical.R"), local = TRUE)
  
  ##############save file to share
  observeEvent(list(input$gwas_file, input$gwas_p_col, input$gwas_sep, input$gwas_header), {
    df2 <- gwas_df()
    req(is.data.frame(df2), nrow(df2) > 0)
    
    sid <- gi_sid(session)
    p   <- gi_state_paths(sid)
    
    saveRDS(df2, p$gwas)
    
    gi_write_state(sid, list(
      has_gwas  = TRUE,
      gwas_rds  = p$gwas
    ))
  }, ignoreInit = TRUE)
  
  # INFO panels
  
  #### ---------------------------------------------------------------------------
  # info modals
  #### ---------------------------------------------------------------------------
  
  observeEvent(input$info_00, {
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
    }, error = function(e) NULL)
  })
  
  output$gwas_p_selector <- renderUI({
    p <- gwas_preview()
    if (is.null(p) || !nrow(p)) return(helpText("Could not read GWAS (preview)."))
    cols <- names(p)
    selectInput(
      "gwas_p_col",
      "Select p-value column",
      choices = cols,
      selected = if ("P" %in% cols) "P" else cols[1]
    )
  })
  
  
  gwas_df <- reactive({
    req(input$gwas_file)
    df <- readr::read_delim(
      input$gwas_file$datapath,
      delim = input$gwas_sep %||% "\t",
      col_names = isTRUE(input$gwas_header),
      show_col_types = FALSE,
      progress = FALSE
    )
    validate(need(is.data.frame(df) && nrow(df) > 0, "Empty GWAS file."))
    
    # detect / standardize columns
    chr_col <- pick_col(df, c("CHR","chr","chrom","CHROM","chromosome"))
    bp_col  <- pick_col(df, c("BP","bp","POS","pos","position","POSITION"))
    snp_col <- pick_col(df, c("SNP","snp","rsid","RSID","marker","ID","id"))
    p_col   <- input$gwas_p_col %||% pick_col(df, c("P","p","PVAL","pval","P_VALUE","p_value"))
    
    validate(
      need(!is.null(chr_col), "GWAS: missing chromosome column."),
      need(!is.null(bp_col),  "GWAS: missing position column."),
      need(!is.null(p_col),   "GWAS: missing p-value column.")
    )
    
    rawP <- df[[p_col]]
    BP   <- if (is.numeric(df[[bp_col]])) df[[bp_col]] else suppressWarnings(readr::parse_number(as.character(df[[bp_col]])))
    Pval <- parse_p_robust(rawP)
    CHR  <- chr_map_plink19(df[[chr_col]])
    snp  <- if (!is.null(snp_col)) as.character(df[[snp_col]]) else paste0("chr", norm_chr_generic(df[[chr_col]]), ":", BP)
    
    df2 <- tibble::tibble(
      CHR = CHR,
      BP  = as.numeric(BP),
      snp = snp,
      Pval = as.numeric(Pval),
      logp = -log10(as.numeric(Pval))
    ) %>%
      filter(!is.na(CHR), !is.na(BP), !is.na(Pval), Pval > 0)
    
    df2
  })
  
  
  # --- IMPORTANT: only allow build_clusters in NEW mode (so canonical observer doesn't run in OLD) ---
  gwas_df_for_clustering <- reactive({
    validate(need(identical(input$use_preloaded_catalog, "new"),
                  "Switch to NEW mode to build clusters."))
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
  
  # Replace Catalog's clusters_val with the canonical one
  clusters_val <- clusters_engine$clusters_cur
  gitools_deeplink_catalog(session, clusters_val)
  
  observeEvent(input$build_ranges, {
    validate(need(identical(input$use_preloaded_catalog, "new"),
                  "Switch to NEW mode to build clusters."))
    
    # Reset downstream (Catalog path) — keep your current behavior
    catalog_final_path_csv(NULL)
    catalog_final_path_rds(NULL)
    
    # Read clusters computed by the canonical engine
    clusters <- clusters_val()
    validate(need(is.data.frame(clusters) && nrow(clusters) > 0, "No clusters generated."))
    
    # Ensure workdir exists (keep your behavior)
    validate(need(exists("workdir") && nzchar(workdir), "workdir is not defined."))
    if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
    
    # ---- helper to write .range (reuse your exact function body if you prefer) ----
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
    
    # -------------------------
    # MASTER: write hub state
    # -------------------------
    sid <- gi_sid(session)
    p   <- gi_state_paths(sid)
    
    # Save clusters and update state stamp (keep exactly your master logic)
    tryCatch({
      saveRDS(clusters, p$clus)
    }, error = function(e) {
      cat("[GI][MASTER] could not save clusters rds:", conditionMessage(e), "\n")
    })
    
    prev <- gi_read_state(sid)
    
    st <- list(
      stamp          = gi_bump_stamp(prev$stamp %||% 0),
      sid            = sid,
      updated_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      has_gwas       = TRUE,
      has_clusters   = TRUE,
      
      cluster_method = as.character(input$cluster_method %||% "window"),
      thr_type       = if (identical(input$cluster_method %||% "window", "window")) "logp" else "min_logp",
      thr_value      = if (identical(input$cluster_method %||% "window", "window")) {
        as.numeric(input$pthr %||% NA_real_)
      } else {
        as.numeric(input$min_logp %||% NA_real_)
      },
      
      pthr           = as.numeric(input$pthr %||% NA_real_),
      min_logp       = as.numeric(input$min_logp %||% NA_real_),
      flank          = as.integer(input$flank %||% 0L),
      
      flank_bp       = as.integer(input$flank %||% 0L),
      hits_mode      = as.character(input$hits_mode %||% NA_character_),
      min_hits       = as.integer(input$min_hits %||% NA_integer_),
      win_bp         = as.integer(input$win_bp %||% NA_integer_),
      step_bp        = as.integer(input$step_bp %||% NA_integer_)
    )
    
    gi_write_state(sid, st)
    
    cat("[GI][MASTER] wrote:\n",
        "  ", p$json, "\n",
        "  ", p$gwas, "\n",
        "  ", p$clus, "\n")
    
    # -------------------------
    # UI preview
    # -------------------------
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
  
  hits_df <- reactive({
    df <- gwas_df(); req(nrow(df) > 0)
    req(input$pthr)
    df %>% filter(logp >= input$pthr) %>% arrange(desc(logp)) %>% select(CHR, BP, snp, p = Pval, logp)
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
                buttons    = c("copy", "csv", "excel", "pdf", "print"),
                pageLength = 10,
                scrollX    = TRUE
              ))
  })
  
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
    
    # PRECOMPUTED MODE (old)
    if (identical(input$use_preloaded_catalog, "old")) {
      
      req(input$catalog_preload)
      path <- input$catalog_preload$datapath
      ext  <- tolower(tools::file_ext(input$catalog_preload$name))
      
      df <- tryCatch({
        if (ext == "rds") readRDS(path) else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      
      validate(need(is.data.frame(df), "Invalid precomputed Catalog file."))
      validate(need("cluster_id" %in% names(df), "Precomputed Catalog file must contain 'cluster_id'."))
      
      if (ext == "rds") {
        catalog_final_path_rds(path)
        catalog_final_path_csv(NULL)
      } else {
        catalog_final_path_csv(path)
        catalog_final_path_rds(NULL)
      }
      
      return(df)
    }
    
    # NEW MODE: read FINAL (prefer RDS)
    rds <- catalog_final_path_rds()
    csv <- catalog_final_path_csv()
    
    if ((is.null(rds) || !file.exists(rds)) && (is.null(csv) || !file.exists(csv))) {
      return(tibble::tibble())
    }
    
    df <- NULL
    if (!is.null(rds) && file.exists(rds)) {
      df <- readRDS(rds)
    } else {
      df <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    if (!"cluster_id" %in% names(df)) return(tibble::tibble())
    df
  })
  
  # canonical catalog_data for DT/plots
  catalog_data <- reactive({
    dt <- catalog_hits_df()
    # print(head(dt))
    if (is.null(dt) || !is.data.frame(dt)) tibble() else dt
  })
  
  # ---------------------------------------------------------------------------
  # OLD MODE: infer clusters from preloaded file (so Manhattan/LD/Enrich work)
  # ---------------------------------------------------------------------------
  
  observeEvent(catalog_hits_df(), {
    if (!identical(input$use_preloaded_catalog, "old")) return()
    dt <- catalog_hits_df()
    if (!is.data.frame(dt) || !nrow(dt)) return()
    
    cl <- infer_clusters_from_catalog_hits(dt)
    if (is.data.frame(cl) && nrow(cl)) {
      
      # Canonicalize cluster_id format (remove "cluster_" prefix, enforce chrN_n)
      cl <- standardize_cluster_ids(cl)
      
      cl <- add_n_catalog(cl, dt)
      clusters_val(cl)
    }
  }, ignoreInit = TRUE)
  
  
  
  # ---------------------------------------------------------------------------
  # Step 4: assign catalog hits (new mode only) + persist final
  # ---------------------------------------------------------------------------
  
  catalog_log_text <- reactiveVal("")
  append_catalog_log <- function(...) {
    old <- catalog_log_text()
    catalog_log_text(paste(c(old, paste(..., collapse = " ")), collapse = "\n"))
  }
  output$catalog_log <- renderText(catalog_log_text())
  
  observeEvent(input$run_catalog, {
    
    validate(need(identical(input$use_preloaded_catalog, "new"),
                  "Switch to NEW mode to run Catalog assignment."))
    
    req(catalog_all())
    req(clusters_val())
    
    # reset log
    catalog_log_text("")
    append_catalog_log("[START] run_catalog")
    
    # opcional: desactivar botó mentre corre (si tens shinyjs::useShinyjs() al UI)
    if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::disable("run_catalog")
    
    tryCatch({
      
      withProgress(message = "Assigning GWAS Catalog hits to clusters…", value = 0, {
        
        incProgress(0.05, detail = "Loading inputs")
        cl <- clusters_val()
        cat_raw <- catalog_all()
        append_catalog_log(sprintf("[1] Inputs loaded | clusters=%d | catalog=%d", nrow(cl), nrow(cat_raw)))
        
        incProgress(0.20, detail = "Preparing Catalog resource")
        cat_pre <- catalog_prepare(cat_raw)
        append_catalog_log("[2] Catalog prepared")
        
        # free memory early
        rm(cat_raw); gc()
        
        incProgress(0.45, detail = "Extracting hits per cluster")
        hits <- extract_catalog_by_clusters(cl, cat_pre)
        append_catalog_log(sprintf("[3] Hits extracted | hits=%d", nrow(hits)))
        
        rm(cat_pre); gc()
        
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
        
        catalog_final_path_rds(out_rds)
        catalog_final_path_csv(out_csv)
        
        append_catalog_log(paste0("[5] Saved: ", basename(out_csv), " | ", basename(out_rds)))
        
        incProgress(0.05, detail = "Done")
        append_catalog_log("[END] run_catalog finished OK")
      })
      
      if (!nrow(hits)) {
        showNotification("No Catalog hits found in clusters.", type = "warning", duration = 6)
      } else {
        showNotification("Catalog hits assigned to clusters (n_catalog updated).", type = "message", duration = 4)
      }
      
    }, error = function(e) {
      
      catalog_final_path_csv(NULL)
      catalog_final_path_rds(NULL)
      
      cl <- clusters_val()
      if (!is.null(cl) && nrow(cl)) {
        cl$n_catalog <- 0L
        clusters_val(cl)
      }
      
      msg <- paste0("ERROR (run_catalog): ", conditionMessage(e))
      output$catalog_summary <- renderText(msg)
      
      append_catalog_log("[ERROR]", msg)
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
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE
      )
    ) %>%
      formatRound(columns = c("top_logp", "cluster_size_kb"), digits = 2)
  })
  
  
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
          buttons    = c("copy", "csv", "excel", "pdf", "print"),
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
  })
  
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
      plotly::layout(showlegend = FALSE)
    
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
  observeEvent(plotly::event_data("plotly_relayout", source = "manhattan_combo"), {
    ev <- plotly::event_data("plotly_relayout", source = "manhattan_combo")
    if (is.null(ev) || !length(ev)) return()
    
    # Ignora relayouts d'inicialització / autosize
    if (all(names(ev) %in% c("autosize", "xaxis.autorange", "yaxis.autorange"))) return()
    if (!any(grepl("^xaxis\\.range\\[", names(ev)))) return()  # només quan hi ha rangs explícits
    
    # ---- el teu codi actual de captura de rangs ----
    # xmin <- ev[["xaxis.range[0]"]]
    # xmax <- ev[["xaxis.range[1]"]]
    # ... actualitza region_selected ...
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
    
    # scope == "pick"
    req(input$dt_cluster_id)
    dt %>% dplyr::filter(cluster_id == input$dt_cluster_id)
  })
  
  # ---- info note ----
  output$dt_scope_note <- renderUI({
    dt <- tryCatch(catalog_data(), error = function(e) NULL)
    if (!is.data.frame(dt) || !nrow(dt)) return(NULL)
    
    dt_uni <- catalog_universe_data()
    uni_ok <- is.data.frame(dt_uni) && nrow(dt_uni) > 0
    uni_n  <- if (uni_ok) nrow(dt_uni) else nrow(dt)
    
    scope <- input$dt_scope %||% "pick"
    msg <- NULL
    
    if (identical(scope, "global")) {
      if ("cluster_id" %in% names(dt)) {
        n_cl <- sum(!is.na(dt$cluster_id) & nzchar(dt$cluster_id))
        msg <- paste0(
          "<b>Global mode:</b> subset = clustered hits (", n_cl, "), universe = ",
          if (uni_ok) "reference catalog" else "current catalog",
          " (", uni_n, ")."
        )
        if (!uni_ok) {
          msg <- paste0(
            msg,
            "<br><b>Warning:</b> missing www/gwas_catalog_simplified.rds → using current catalog as universe."
          )
        }
      }
    }
    
    if (identical(scope, "selected")) {
      cl <- selected_cluster()
      if (is.null(cl) || !nrow(cl)) {
        msg <- "<b>Note:</b> No cluster selected in the Clusters table. Select one, or use <i>Pick a cluster here</i>."
      }
    }
    
    if (!is.null(msg)) {
      return(HTML(paste0(
        "<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>",
        msg, "</div>"
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
                                    min_a = 2, min_Tt = 5, max_terms = 2000)
    if (!nrow(tab)) {
      return(DT::datatable(data.frame(Message="No terms."), options=list(dom="t"), rownames=FALSE))
    }
    
    out <- tab %>%
      dplyr::transmute(
        term,
        in_subset = a,
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format="e", digits=2),
        FDR = formatC(as.numeric(p_adj), format="e", digits=2)
      )
    
    DT::datatable(out, rownames=FALSE,               
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons    = c("copy", "csv", "excel", "pdf", "print"),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
  })
  
  output$enrich_disease_plot <- renderPlot({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    dis_col     <- pick_col(dt,     c("DISEASE","Disease","disease"))
    dis_col_uni <- pick_col(dt_uni, c("DISEASE","Disease","disease"))
    validate(need(!is.null(dis_col), "DISEASE column not found."))
    validate(need(!is.null(dis_col_uni), "DISEASE column not found in universe."))
    
    tab <- fisher_enrich_terms_fast(dt2[[dis_col]], dt_uni[[dis_col]],
                                    alternative = "greater",
                                    min_a = 2, min_Tt = 5, max_terms = 2000)
    if (!nrow(tab)) return(NULL)
    
    top <- tab %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::mutate(score = -log10(p_adj)) %>%
      dplyr::arrange(score) %>%
      dplyr::mutate(
        term_short = short_label(term, max_chars = 55, wrap_width = 18),
        term_short = factor(term_short, levels = term_short),
        zebra = factor(seq_len(dplyr::n()) %% 2)
      )
    
    ggplot2::ggplot(top, ggplot2::aes(x = term_short, y = score, fill = zebra)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("0"="orange","1"="darkgreen"), guide = "none") +
      ggplot2::labs(title = "Top enriched terms", x = NULL, y = "-log10(FDR)") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
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
                                    min_a = 2, min_Tt = 5, max_terms = 2000)
    if (!nrow(tab)) {
      return(DT::datatable(data.frame(Message="No terms."), options=list(dom="t"), rownames=FALSE))
    }
    
    out <- tab %>%
      dplyr::transmute(
        term,
        in_subset = a,
        OR  = sprintf("%.3f", as.numeric(or)),
        p   = formatC(as.numeric(p),     format="e", digits=2),
        FDR = formatC(as.numeric(p_adj), format="e", digits=2)
      )
    
    DT::datatable(out, rownames=FALSE,               
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons    = c("copy", "csv", "excel", "pdf", "print"),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
  })
  
  output$enrich_trait_plot <- renderPlot({
    x <- dt_enrich_res(); req(x)
    dt  <- x$dt
    dt2 <- x$dt2
    dt_uni <- x$dt_uni
    
    trait_col     <- pick_col(dt,     c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    trait_col_uni <- pick_col(dt_uni, c("MAPPED_TRAIT","mapped_trait","TRAIT","trait"))
    validate(need(!is.null(trait_col), "MAPPED_TRAIT column not found."))
    validate(need(!is.null(trait_col_uni), "MAPPED_TRAIT column not found in universe."))
    
    tab <- fisher_enrich_terms_fast(dt2[[trait_col]], dt_uni[[trait_col]],
                                    alternative = "greater",
                                    min_a = 2, min_Tt = 5, max_terms = 2000)
    if (!nrow(tab)) return(NULL)
    
    top <- tab %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::mutate(score = -log10(p_adj)) %>%
      dplyr::arrange(score) %>%
      dplyr::mutate(
        term_short = short_label(term, max_chars = 55, wrap_width = 18),
        term_short = factor(term_short, levels = term_short),
        zebra = factor(seq_len(dplyr::n()) %% 2)
      )
    
    ggplot2::ggplot(top, ggplot2::aes(x = term_short, y = score, fill = zebra)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("0"="orange","1"="darkgreen"), guide = "none") +
      ggplot2::labs(title = "Top enriched terms", x = NULL, y = "-log10(FDR)") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
    
  })
  
  # ---------------------------------------------------------------------------
  # GO/KEGG enrichment (NonSyn-style) driven by Catalog genes
  # ---------------------------------------------------------------------------
  
  output$func_cluster_ui <- renderUI({
    cl <- clusters_val()
    validate(need(is.data.frame(cl) && nrow(cl) > 0, "No clusters available."))
    
    chr_choices <- sort(unique(as.integer(cl$chr)))
    chr_labels <- chr_label_plink(chr_choices)
    names(chr_choices) <- paste0("chr", chr_labels)
    
    tagList(
      selectInput("func_chr", "Chromosome", choices = chr_choices, selected = chr_choices[1]),
      uiOutput("func_cluster_id_ui")
    )
  })
  
  output$func_cluster_id_ui <- renderUI({
    cl <- clusters_val(); req(nrow(cl) > 0)
    req(input$func_chr)
    x <- cl %>% filter(as.integer(chr) == as.integer(input$func_chr)) %>% arrange(start, end)
    selectInput("func_cluster_id", "Cluster", choices = x$cluster_id, selected = x$cluster_id[1])
  })
  
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
    
    if (identical(input$func_scope, "cluster")) {
      validate(need("cluster_id" %in% names(dt), "Catalog hits missing cluster_id."))
      req(input$func_cluster_id)
      dt <- dt %>% filter(cluster_id == input$func_cluster_id)
    }
    
    split_multi_unique(dt[[gene_col]])
  })
  
  # Map genes once (cache)
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
    
    # In GLOBAL mode: default background
    if (identical(scope, "global")) return(NULL)
    
    # In CLUSTER mode: user decides background
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
               "GO enrichment uses the default <i>OrgDb</i>-based background (org.Hs.eg.db; ontology-specific).")
      } else {
        paste0("<b>Background:</b> ", bg_lbl, "<br>",
               "GO enrichment uses only genes present in your Catalog dataset as background (cluster mode).")
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
    
    HTML(paste0("<div style='background:#f6f6f6;border:1px solid #ddd;padding:10px;border-radius:8px;'>", msg, "</div>"))
  })
  
  # ---- triggers: IMPORTANT (NULL evita execució al carregar panell) ----
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
  
  go_enrich_raw <- eventReactive(go_trigger(), {
    withProgress(message = "Running GO enrichment…", value = 0, {
      
      gene <- entrez_scope()
      incProgress(0.2, detail = paste("Genes:", length(gene)))
      
      uni  <- universe_entrez_for_scope()
      incProgress(0.1, detail = if (is.null(uni)) "Universe: default (OrgDb)" else paste("Universe:", length(uni)))
      
      ontos <- input$go_ontos
      if (is.null(ontos) || length(ontos) == 0 || all(is.na(ontos))) {
        ontos <- c("BP","CC","MF")
      }
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
    if (!is.data.frame(df) || !nrow(df)) return(tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$go_topn %||% 10
    
    df2 <- df %>%
      mutate(p.adjust = as.numeric(p.adjust),
             score = -log10(p.adjust)) %>%
      filter(is.finite(p.adjust), p.adjust <= pcut) %>%
      group_by(Ontology) %>%
      arrange(p.adjust, .by_group = TRUE) %>%
      slice_head(n = topn) %>%
      ungroup()
    
    if (!nrow(df2)) return(tibble())
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
        pvalue    = sprintf("%.3f", as.numeric(pvalue)),
        p.adjust  = sprintf("%.3f", as.numeric(p.adjust))
      )
    
    DT::datatable(out, rownames = FALSE,
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons    = c("copy", "csv", "excel", "pdf", "print"),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
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
        Ontology    = as.character(Ontology),
        pvalue      = suppressWarnings(as.numeric(pvalue)),
        p.adjust    = suppressWarnings(as.numeric(p.adjust)),
        term_short  = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description),
        score       = -log10(pvalue)   # IMPORTANT: bar height = -log10(pvalue)
      ) %>%
      dplyr::filter(
        Ontology %in% c("BP","CC","MF"),
        is.finite(pvalue), is.finite(p.adjust),
        nzchar(term_short)
      )
    
    if (!nrow(df0)) return(tibble::tibble())
    
    # Primer intent: FDR <= pcut
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
    
    # Fallback: top terms per FDR
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
    
    fb     <- isTRUE(attr(dat, "fallback"))
    ttl    <- paste0("GO enrichment (", scope_label(), ")")
    pcut   <- input$enrich_pcut %||% 0.05
    
    subttl <- if (fb) {
      paste0(
        "No terms at FDR ≤ ", pcut,
        " — showing top terms ranked by FDR; bar height = -log10(pvalue)"
      )
    } else {
      paste0("FDR ≤ ", pcut, " (bar height = -log10(pvalue))")
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
  
  
  kegg_enrich_raw <- eventReactive(kegg_trigger(), {
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
          keyType       = "kegg",
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
    if (!is.data.frame(df) || !nrow(df)) return(tibble())
    
    pcut <- input$enrich_pcut %||% 0.05
    topn <- input$enrich_kegg_top %||% 15
    
    df2 <- df %>%
      mutate(p.adjust = as.numeric(p.adjust),
             score = -log10(p.adjust)) %>%
      filter(is.finite(p.adjust), p.adjust <= pcut) %>%
      arrange(p.adjust) %>%
      slice_head(n = topn)
    
    if (!nrow(df2)) return(tibble())
    df2
  })
  
  output$kegg_table <- DT::renderDT({
    df <- kegg_enrich_tbl()
    if (!nrow(df)) {
      return(DT::datatable(data.frame(Message="No KEGG pathways passed the cutoff."),
                           options=list(dom="t"), rownames=FALSE))
    }
    
    out <- df %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      dplyr::mutate(
        pvalue    = sprintf("%.3f", as.numeric(pvalue)),
        p.adjust  = sprintf("%.3f", as.numeric(p.adjust))
      )
    
    DT::datatable(out, rownames = FALSE,
                  extensions = "Buttons",
                  options    = list(
                    dom        = "Bfrtip",
                    buttons    = c("copy", "csv", "excel", "pdf", "print"),
                    pageLength = 10,
                    scrollX    = TRUE
                  ))
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
        term_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description),
        score      = -log10(pvalue)   # IMPORTANT: bar height = -log10(pvalue)
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
  
  
  ##############################################################################
  ############################ LD Module   #####################################
  ################# GWAS -> LD module canonical cluster input###################
  # -----------------------------
  # LD PANEL
  # -----------------------------
  clusters_r <- reactive({
    build_ld_clusters_from_app()
  })
  
  candidates_r <- reactive({
    build_ld_candidates_from_app()
  })
  
  ld_module_server(
    "ld",
    clusters_r   = clusters_r,
    candidates_r = candidates_r
    # keep_dir = file.path(cfg$resources, "POP")  # opcional
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
    
    upd("use_preloaded_catalog", function() shiny::updateRadioButtons(session, "use_preloaded_catalog", selected = "new"))
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
  
  ##≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  ##≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  
  
}

shinyApp(ui, server)
