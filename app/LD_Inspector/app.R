# ==============================
# app.R (LD_Inspector) - PORTABLE + feeds _shared/mod_ld.R
# ==============================

options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)
library(readr)
library(dplyr)
library(DT)
library(tibble)

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

# FIX: assegura que tags és el de shiny (llista), no una funció
tags <- shiny::tags
HTML <- shiny::HTML
div  <- shiny::div

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------
# Portable config (GItools)
# -----------------------------
cfg_file <- normalizePath(file.path("..", "..", "config.R"), winslash = "/", mustWork = FALSE)
if (!file.exists(cfg_file)) {
  # fallback (si obres l'app des del root)
  cfg_file <- normalizePath("config.R", winslash = "/", mustWork = FALSE)
}
stopifnot(file.exists(cfg_file))
source(cfg_file, local = TRUE)

cfg <- gi_cfg()
gi_root        <- cfg$root
gi_shared_root <- cfg$shared
gi_res_root    <- cfg$resources
gi_pop_dir     <- cfg$pop_dir

# ------------------------------------------------------------------
# IMPORTANT: força POP dir per al mòdul LD (si no, queda buit i peta)
# ------------------------------------------------------------------
Sys.setenv(
  GITOOLS_ROOT      = gi_root,
  GITOOLS_SHARED    = gi_shared_root,
  GITOOLS_RESOURCES = gi_res_root,
  GITOOLS_POP_DIR   = gi_pop_dir
)

# Debug ràpid (opcional, però útil)
cat("[LD_Inspector] gi_pop_dir =", gi_pop_dir, "\n")
cat("[LD_Inspector] dir.exists(gi_pop_dir) =", dir.exists(gi_pop_dir), "\n")


# -----------------------------
# Shared modules
# -----------------------------
dl_file <- file.path(gi_shared_root, "GItools_local_deeplinks_ALL_IN_ONE.R")
if (file.exists(dl_file)) source(dl_file, local = TRUE)

stopifnot(nzchar(gi_pop_dir))
stopifnot(dir.exists(gi_pop_dir))
stopifnot(length(list.files(gi_pop_dir, pattern = "\\.txt$", full.names = TRUE)) > 0)

ld_file <- file.path(gi_shared_root, "mod_ld.R")
stopifnot(file.exists(ld_file))

# IMPORTANT: posa el mòdul al globalenv perquè ld_module_ui/ld_module_server
# siguin visibles sense problemes dins server()
source(ld_file, local = globalenv())

stopifnot(exists("ld_module_ui"), exists("ld_module_server"))

cat("gi_pop_dir =", gi_pop_dir, "\n")
# -----------------------------
# Robust file reader (auto delim)
# -----------------------------
guess_delim_by_ext <- function(filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext == "csv") return(",")
  if (ext == "tsv") return("\t")
  return(NULL)
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

# -----------------------------
# Helpers (schema normalization)
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
  x[x %in% c("MT","M","MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}


# -----------------------------
# UI
# -----------------------------
ui <- navbarPage(
  title = div(style = "font-weight:700; font-size:22px; color:#1A4E8A;", HTML("🔗 LD Inspector")),
  id = "topnav",
  tags$head(
    tags$style(HTML("
      .panelGrey{ background:#f2f2f2; padding:14px; border-radius:10px; }
      .panelGrey.padDT .dataTables_wrapper{ background:#f2f2f2; padding:10px; border-radius:8px; }
      .panelGrey.padDT table.dataTable{ background:transparent !important; }
      .panelGrey.padDT table.dataTable thead th{ background:transparent !important; }
      .panelGrey.padDT{ min-height:420px; overflow:auto; }
    "))
  ),
  tabPanel(
    title = div(style = "font-weight:700; font-size:22px; color:#1A4E8A;", HTML("🧾 Input files")),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        h4("1) Load clusters", style = "color:#1A4E8A; font-weight:700;"),
        fileInput(
          "clusters_file",
          "Clusters (CSV/TSV/TXT) — multiple files",
          accept = c(".csv",".tsv",".txt"),
          multiple = TRUE
        ),
        checkboxInput("clusters_header", "Header", TRUE),
        radioButtons(
          "clusters_sep", "Separator",
          choices  = c("Auto" = "auto", "Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"),
          selected = "auto"
        ),
        tags$hr(),
        
        h4("2) Load candidate SNPs", style = "color:#1A4E8A; font-weight:700;"),
        fileInput(
          "candidates_file",
          "Candidates (CSV/TSV/TXT) — multiple files",
          accept = c(".csv",".tsv",".txt"),
          multiple = TRUE
        ),
        checkboxInput("candidates_header", "Header", TRUE),
        radioButtons(
          "candidates_sep", "Separator",
          choices  = c("Auto" = "auto", "Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"),
          selected = "auto"
        ),
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = div(style = "font-weight:700; font-size:22px; color:#1A4E8A;", HTML("🧾 Tables")),
          h3("Cluster table"),  
          div(class="panelGrey",
          DTOutput("clusters_dt")),
          h3("Candidate table"), 
          div(class="panelGrey",
          DTOutput("candidates_dt"))),
          tabPanel(
            "Notes",
            tags$div(
              style = "background:#f6f6f6;border:1px solid #ddd;padding:12px;border-radius:10px;",
              HTML(paste0(
                "<b>Clusters</b> accepts: chr + (cluster_start/cluster_end) or start/end or start_bp/end_bp. ",
                "Optional: cluster_id.<br>",
                "<b>Candidates</b> accepts: chr + (pos_ini or position/POS/BP) + (id_hit or rsid).<br><br>",
                "The module expects (and here we give it) canonical schema: ",
                "<code>clusters: chr,start,end,cluster_id</code> and <code>candidates: rsid,chr,position,classe</code>."
              ))
            )
          )
        )
      )
    )
  ),
  
  tabPanel(
    title = div(style = "font-weight:700; font-size:22px; color:#1A4E8A;", HTML("🧩 LD plot")),
    ld_module_ui("ld")
  )
)

# -----------------------------
# SERVER
# -----------------------------
server <- function(input, output, session) {
  
  # --- RAW tables (merge multiple files) ---
  clusters_raw <- reactive({
    req(input$clusters_file)
    files <- input$clusters_file
    req(nrow(files) >= 1)
    
    user_delim <- input$clusters_sep %||% "auto"
    if (identical(user_delim, "auto") || identical(user_delim, "")) user_delim <- NULL
    
    lst <- lapply(seq_len(nrow(files)), function(i) {
      df <- read_any_delim(
        path       = files$datapath[i],
        name       = files$name[i],
        header     = isTRUE(input$clusters_header),
        user_delim = user_delim
      )
      validate(need(is.data.frame(df) && nrow(df) > 0, paste0("Empty clusters file: ", files$name[i])))
      df
    })
    
    dplyr::bind_rows(lst)
  })
  
  candidates_raw <- reactive({
    req(input$candidates_file)
    files <- input$candidates_file
    req(nrow(files) >= 1)
    
    user_delim <- input$candidates_sep %||% "auto"
    if (identical(user_delim, "auto") || identical(user_delim, "")) user_delim <- NULL
    
    lst <- lapply(seq_len(nrow(files)), function(i) {
      df <- read_any_delim(
        path       = files$datapath[i],
        name       = files$name[i],
        header     = isTRUE(input$candidates_header),
        user_delim = user_delim
      )
      validate(need(is.data.frame(df) && nrow(df) > 0, paste0("Empty candidates file: ", files$name[i])))
      df
    })
    
    dplyr::bind_rows(lst)
  })
  
  # --- NORMALIZED tables for mod_ld.R (KEY FIX) ---
  clusters_for_ld <- reactive({
    df <- clusters_raw()
    validate(need(is.data.frame(df) && nrow(df) > 0, "No clusters loaded."))
    
    chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
    st_col  <- pick_col(df, c("cluster_start","start","START","start_bp","FROM","from","bp1"))
    en_col  <- pick_col(df, c("cluster_end","end","END","end_bp","TO","to","bp2"))
    id_col  <- pick_col(df, c("cluster_id","CLUSTER_ID","cluster","id","ID","label","cluster_chr_n"))
    
    validate(need(!is.null(chr_col) && !is.null(st_col) && !is.null(en_col),
                  "Clusters must have chr + start + end (e.g. chr,cluster_start,cluster_end)."))
    
    out <- df %>%
      dplyr::transmute(
        chr = chr_map_plink19(.data[[chr_col]]),
        start = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
        end   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]])))),
        cluster_id = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_
      ) %>%
      dplyr::filter(is.finite(chr), is.finite(start), is.finite(end)) %>%
      dplyr::mutate(
        start = pmin(start, end),
        end   = pmax(start, end),
        cluster_id = trimws(cluster_id),
        cluster_id = ifelse(is.na(cluster_id) | !nzchar(cluster_id),
                            paste0("cluster_chr", chr, "_", dplyr::row_number()),
                            cluster_id)
      ) %>%
      dplyr::arrange(chr, start, end)
    
    validate(need(nrow(out) > 0, "Clusters normalized to 0 rows (check chr/start/end columns)."))
    out
  })
  
  candidates_for_ld <- reactive({
    df <- candidates_raw()
    validate(need(is.data.frame(df) && nrow(df) > 0, "No candidates loaded."))
    
    chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
    pos_col <- pick_col(df, c("pos_ini","POS","pos","position","BP","bp","start","pos_start"))
    id_col  <- pick_col(df, c("id_hit","rsid","RSID","SNP","snp","id","ID","marker","MARKER"))
    cl_col  <- pick_col(df, c("classe","class","CLASS","type","TYPE"))
    
    validate(need(!is.null(chr_col) && !is.null(pos_col) && !is.null(id_col),
                  "Candidates must have chr + position + id (e.g. chr,pos_ini,id_hit)."))
    
    out <- df %>%
      dplyr::transmute(
        chr = chr_map_plink19(.data[[chr_col]]),
        position = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
        rsid = trimws(as.character(.data[[id_col]])),
        classe = if (!is.null(cl_col)) trimws(as.character(.data[[cl_col]])) else NA_character_
      ) %>%
      dplyr::filter(is.finite(chr), is.finite(position), !is.na(rsid), nzchar(rsid)) %>%
      dplyr::distinct(chr, position, rsid, .keep_all = TRUE) %>%
      dplyr::arrange(chr, position)
    
    validate(need(nrow(out) > 0, "Candidates normalized to 0 rows (check chr/pos/id columns)."))
    out
  })
  
  # --- previews ---
  output$clusters_dt <- DT::renderDT({
    req(clusters_for_ld())
    DT::datatable(
      head(clusters_for_ld(), 200),
         extensions = "Buttons",
    options    = list(
      dom         = "Bfrtip",
      buttons     = c("copy", "csv", "excel", "pdf", "print"),
      scrollX     = TRUE,
      pageLength  = 10
    )
  )
  })
  
  output$candidates_dt <- DT::renderDT({
    req(candidates_for_ld())
    DT::datatable(
      head(candidates_for_ld(), 200),
      extensions = "Buttons",
      options    = list(
        dom         = "Bfrtip",
        buttons     = c("copy", "csv", "excel", "pdf", "print"),
        scrollX     = TRUE,
        pageLength  = 10
      )
    )
  })
  
  # --- FEED the shared LD module ---
  # (Això ha d'omplir el selector de cromosoma dins del mòdul)
  ld_module_server(
    id = "ld",
    clusters_r = clusters_for_ld,
    candidates_r = candidates_for_ld,
    default_pops_dir = gi_pop_dir,
    app_tag = "ld_inspector"
  )
  
}

shinyApp(ui, server)
