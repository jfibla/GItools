# R/mod_ld.R  ---------------------------------------------------------------
# LD module (copied from LD_Inspector; ONLY input source changed)
# - clusters_r(): data.frame with chr/start/end (canonical or compatible)
# - candidates_r(): data.frame with rsid/chr/position (canonical or compatible)
# - Paths via UI textInputs: plink_bin, bfile_ref_prefix, ld_pops_dir
# --------------------------------------------------------------------------

# deps used in LD_Inspector
# library(shiny)
# library(dplyr)
# library(readr)
# library(ggplot2)
# library(plotly)
# library(scales)
# library(DT)
# library(tibble)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(GenomeInfoDb)
# library(S4Vectors)
# library(IRanges)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
# library(AnnotationDbi)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

chr_map_plink19 <- function(x){
  x <- toupper(as.character(x))
  x <- trimws(x)
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

run_plink <- function(args, out_prefix, plink_bin) {
  res <- tryCatch(
    system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) paste("system2 error:", conditionMessage(e))
  )
  st <- attr(res, "status"); if (is.null(st)) st <- 0
  
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
  stopifnot(all(c("SNP","BP") %in% names(fl)))
  
  fl <- fl %>%
    dplyr::arrange(BP, SNP) %>%
    dplyr::mutate(ix = dplyr::row_number())
  
  snp_to_ix <- stats::setNames(fl$ix, fl$SNP)
  
  nm <- names(blocks_raw)
  snps_col <- nm[match(TRUE, nm %in% c("SNPS","Snps","snps","SNPs"))]
  bp1_col  <- nm[match(TRUE, nm %in% c("BP1","bp1","FROM_BP","from_bp","START_BP","start_bp"))]
  bp2_col  <- nm[match(TRUE, nm %in% c("BP2","bp2","TO_BP","to_bp","END_BP","end_bp"))]
  
  out <- NULL
  
  # Plan A: SNPS column
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
  
  # Plan B: BP1/BP2
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

# --- prepare inputs (same helpers as LD_Inspector, but fed by reactives) -----
prepare_clusters <- function(df) {
  chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
  st_col  <- pick_col(df, c("start","START","start_bp","cluster_start","FROM","from","bp1"))
  en_col  <- pick_col(df, c("end","END","end_bp","cluster_end","TO","to","bp2"))
  id_col  <- pick_col(df, c("cluster_id","CLUSTER_ID","cluster","id"))
  
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
      cluster_id = trimws(cluster_id)
    ) %>%
    dplyr::arrange(chr, start, end) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(cluster_n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cluster_id = ifelse(is.na(cluster_id) | !nzchar(cluster_id),
                          paste0("cluster_chr", chr_label_plink(chr), "_", cluster_n),
                          cluster_id),
      size_kb = round((end - start) / 1000, 2),
      cluster_key = paste0(cluster_id, "|chr", chr_label_plink(chr), ":", start, "-", end)
    )
  
  out
}

prepare_candidates <- function(df) {
  rs_col  <- pick_col(df, c("rsid","RSID","snp","SNP","id","ID","marker","MARKER","id_hit","ID_HIT"))
  chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
  pos_col <- pick_col(df, c("position","POS","pos","BP","bp","pos_ini","pos_start","start","POS_INI"))
  cl_col  <- pick_col(df, c("classe","class","CLASS","Classe","type","TYPE"))

  if (is.null(rs_col) || is.null(chr_col) || is.null(pos_col)) {
    stop("Candidates table must have rsid/id_hit + chr + position/pos_ini columns. Found: ",
         paste(names(df), collapse = ", "))
  }

  out <- df %>%
    dplyr::transmute(
      rsid     = as.character(.data[[rs_col]]),
      chr      = chr_map_plink19(.data[[chr_col]]),
      position = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
      classe   = if (!is.null(cl_col)) as.character(.data[[cl_col]]) else NA_character_
    ) %>%
    dplyr::mutate(
      rsid   = trimws(.data$rsid),
      classe = trimws(.data$classe)
    ) %>%
    dplyr::filter(nzchar(.data$rsid), is.finite(.data$chr), is.finite(.data$position)) %>%
    dplyr::distinct(.data$rsid, .data$chr, .data$position, .keep_all = TRUE)

  out
}



# keep files now come from input$ld_pops_dir (NOT www/POP)
# keep files now come from input$ld_pops_dir (NOT www/POP)

available_pops_dir <- function(ld_pops_dir) {
  d <- tryCatch(normalizePath(ld_pops_dir, winslash = "/", mustWork = FALSE),
                error = function(e) ld_pops_dir)
  if (is.null(d) || !nzchar(d) || !dir.exists(d)) return(character(0))
  
  # accept .txt / .TXT / .TxT ...
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
  
  # try pop.txt then pop.TXT variants (in case FS is case-sensitive)
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
  
  if (ncol(x) == 1) x <- cbind(x[,1], x[,1]) else x <- x[,1:2]
  colnames(x) <- c("FID","IID")
  
  out <- file.path(tempdir(), paste0("keep_", pop, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  write.table(x, out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  out
}

gi_ld_defaults <- function() {
  
  # Base resources dir (portable)
  if (exists("gi_cfg", mode = "function")) {
    cfg <- gi_cfg()
    res <- cfg$resources %||% ""
    popdir <- cfg$pop_dir %||% ""
  } else {
    res <- Sys.getenv("GITOOLS_RESOURCES", unset = "")
    popdir <- Sys.getenv("GITOOLS_POP_DIR", unset = "")
  }
  
  res <- normalizePath(res, winslash = "/", mustWork = FALSE)
  popdir <- normalizePath(popdir, winslash = "/", mustWork = FALSE)
  
  # plink executable (user can override)
  plink <- Sys.getenv("GITOOLS_PLINK19", unset = "")
  if (!nzchar(plink) && nzchar(res)) {
    cand1 <- file.path(res, "software", "plink19", "plink")
    cand2 <- file.path(res, "software", "plink19", "plink19")
    cand3 <- file.path(res, "software", "plink19")
    plink <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else cand3
  }
  
  # default LD bfile prefix (user can override)
  bfile <- Sys.getenv("GITOOLS_LD_BFILE", unset = "")
  if (!nzchar(bfile) && nzchar(res)) {
    bfile <- file.path(res, "LD_resources",
                       "Merged_FULL_SET_hg38_hgdp.wgs_10000G3_ko07_MAF0.05")
  }
  
  list(plink = plink, bfile = bfile, popdir = popdir, resources = res)
}



# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------
ld_module_ui <- function(id) {
  ns <- shiny::NS(id)
  ld_def <- gi_ld_defaults()
  
  shiny::tags$pre(
    style="font-size:11px; color:#666;",
    paste0(
      "DBG popdir: ", ld_def$popdir, "\n",
      "DBG exists: ", dir.exists(ld_def$popdir)
    )
  )
  
  shiny::tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        
        shiny::wellPanel(
          shiny::h4("Inputs (paths)"),
          
          shiny::textInput(
            ns("plink_bin"),
            "PLINK binary path",
            value = ld_def$plink %||% ""
          ),
          
          shiny::textInput(
            ns("bfile_ref_prefix"),
            "Reference bfile prefix (no .bed/.bim/.fam)",
            value = ld_def$bfile %||% ""
          ),
          
          shiny::textInput(
            ns("ld_pops_dir"),
            "POP keep-files directory (.txt)",
            value = ld_def$popdir %||% ""
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
        
        shiny::numericInput(ns("ld_max_snps_interval"), "Max SNPs in interval (safety)", value = 400, min = 50, step = 50),
        shiny::tags$hr(),
        shiny::actionButton(ns("run_ld_cluster"), "Run LD (PLINK)", icon = shiny::icon("play")),
        shiny::tags$hr(),
        
        shiny::checkboxInput(ns("ld_blocks_enable"), "Haploview blocks (PLINK --blocks)", value = FALSE),
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
            shiny::sliderInput(ns("ld_blocks_strong_highci"), "blocks-strong-highci (≥0.84 in PLINK 1.9)", min = 0.84, max = 0.99, value = 0.90, step = 0.01),
            shiny::sliderInput(ns("ld_blocks_strong_lowci"), "blocks-strong-lowci", min = 0.10, max = 0.95, value = 0.55, step = 0.01)
          ),
         # shiny::actionButton(ns("run_ld_blocks"), "Run blocks (PLINK)", icon = shiny::icon("th"))
        ),
        
        shiny::tags$hr(),
        shiny::verbatimTextOutput(ns("ld_log"))
      ),
      
      shiny::mainPanel(
        width = 9,
        shinycssloaders::withSpinner(plotly::plotlyOutput(ns("ld_plot"), height = "800px")),
        shiny::tags$hr(),
        shiny::h4("SNPs in interval"),
        shinycssloaders::withSpinner(DT::DTOutput(ns("ld_fl_dt"))),
        shiny::tags$hr(),
        shiny::h4("Top LD pairs"),
        shinycssloaders::withSpinner(DT::DTOutput(ns("ld_pairs_dt")))
      )
    )
  )
}



# --------------------------------------------------------------------------
# SERVER
# --------------------------------------------------------------------------
ld_module_server <- function(id, clusters_r, candidates_r, app_tag = NULL, default_pops_dir = NULL, ...) {
  moduleServer(id, function(input, output, session) {

 #   # If caller provided a default POP dir, set it once (won't override user edits)
 #   if (!is.null(default_pops_dir) && nzchar(default_pops_dir)) {
 #     shiny::observeEvent(TRUE, {
 #       cur <- shiny::isolate(input$ld_pops_dir %||% "")
 #       if (!nzchar(cur) || identical(cur, gi_ld_defaults()$popdir)) {
 #         shiny::updateTextInput(session, "ld_pops_dir", value = default_pops_dir)
 #       }
 #     }, once = TRUE)
 #   }
    
    # If caller provided a default POP dir, set it once (won't override user edits)
    if (!is.null(default_pops_dir) && nzchar(default_pops_dir)) {
      shiny::observeEvent(TRUE, {
        cur <- shiny::isolate(input$ld_pops_dir %||% "")
        # Si és buit o no existeix, aplica default
        if (!nzchar(cur) || !dir.exists(cur)) {
          shiny::updateTextInput(session, "ld_pops_dir", value = default_pops_dir)
        }
      }, once = TRUE)
    }

    
    stopifnot(is.function(clusters_r))
    stopifnot(is.function(candidates_r))
    
    append_log <- function(...) {
      txt <- paste0(..., collapse = "")
      cat(txt, "\n")
      isolate({
        old <- input$.__dummy__ %||% ""
      })
    }
    
    # DEBUG immediat
    observe({
      c0 <- tryCatch(candidates_r(), error = function(e) e)
      if (inherits(c0, "error")) {
        cat("[LD][DBG] candidates_r() ERROR: ", conditionMessage(c0), "\n")
      } else {
        cat("[LD][DBG] candidates_r() nrow=", nrow(c0), " | cols=", paste(names(c0), collapse=","), "\n")
        if (nrow(c0) > 0) print(utils::head(c0, 3))
      }
    })
    
    observe({
      c0 <- tryCatch(clusters_r(), error = function(e) e)
      if (inherits(c0, "error")) {
        cat("[LD][DBG] clusters_r() ERROR: ", conditionMessage(c0), "\n")
      } else {
        cat("[LD][DBG] clusters_r() nrow=", nrow(c0), " | cols=", paste(names(c0), collapse=","), "\n")
        if (nrow(c0) > 0) print(utils::head(c0, 3))
      }
    })
    
    # TxDb (same as LD_Inspector)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    # log buffer
    ld_log <- shiny::reactiveVal("")
    append_log <- function(...) {
      old <- ld_log()
      ld_log(paste(c(old, paste(..., collapse = "\n")), collapse = "\n"))
    }
    output$ld_log <- shiny::renderText(ld_log())
    
    workdir <- file.path(tempdir(), paste0("ldmod_", as.integer(Sys.time())))
    dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
    
    # --- canonicalized inputs from app ---
    clusters_df <- shiny::reactive({
      cl0 <- clusters_r()
      shiny::validate(shiny::need(is.data.frame(cl0) && nrow(cl0) > 0, "No clusters available for LD."))
      cl0 <- as.data.frame(cl0)

      # If already standardized, just coerce types / fill missing helpers
      if (all(c("chr","start","end") %in% names(cl0))) {
        cl0$chr   <- suppressWarnings(as.integer(chr_map_plink19(cl0$chr)))
        if (!"cluster_id" %in% names(cl0)) cl0$cluster_id <- NA_character_
        cl0$start <- suppressWarnings(as.integer(readr::parse_number(as.character(cl0$start))))
        cl0$end   <- suppressWarnings(as.integer(readr::parse_number(as.character(cl0$end))))
        cl0 <- cl0 %>% dplyr::filter(is.finite(.data$chr), is.finite(.data$start), is.finite(.data$end))
        cl0$start <- pmin(cl0$start, cl0$end)
        cl0$end   <- pmax(cl0$start, cl0$end)

        if (!"cluster_n" %in% names(cl0)) {
          cl0 <- cl0 %>% dplyr::arrange(.data$chr, .data$start, .data$end) %>%
            dplyr::group_by(.data$chr) %>% dplyr::mutate(cluster_n = dplyr::row_number()) %>% dplyr::ungroup()
        }
        if (!"cluster_id" %in% names(cl0) || all(is.na(cl0$cluster_id) | !nzchar(trimws(as.character(cl0$cluster_id))))) {
          cl0$cluster_id <- paste0("cluster_chr", chr_label_plink(cl0$chr), "_", cl0$cluster_n)
        } else {
          cl0$cluster_id <- trimws(as.character(cl0$cluster_id))
          cl0$cluster_id[!nzchar(cl0$cluster_id)] <- paste0("cluster_chr", chr_label_plink(cl0$chr[!nzchar(cl0$cluster_id)]), "_", cl0$cluster_n[!nzchar(cl0$cluster_id)])
        }

        if (!"cluster_key" %in% names(cl0)) {
          cl0$cluster_key <- paste0(cl0$cluster_id, "|chr", chr_label_plink(cl0$chr), ":", cl0$start, "-", cl0$end)
        }
        if (!"size_kb" %in% names(cl0)) cl0$size_kb <- round((cl0$end - cl0$start)/1000, 2)

        return(cl0)
      }

      # Otherwise, try to parse with flexible column detection
      prepare_clusters(cl0)
    })
    
    candidates_df <- shiny::reactive({
      ca0 <- candidates_r()
      shiny::validate(shiny::need(is.data.frame(ca0) && nrow(ca0) > 0, "No candidates available for LD."))
      ca0 <- as.data.frame(ca0)

      # If already standardized, just coerce types
      if (all(c("rsid","chr","position") %in% names(ca0))) {
        ca0$rsid     <- trimws(as.character(ca0$rsid))
        ca0$chr      <- suppressWarnings(as.integer(chr_map_plink19(ca0$chr)))
        ca0$position <- suppressWarnings(as.integer(readr::parse_number(as.character(ca0$position))))
        if (!"classe" %in% names(ca0)) ca0$classe <- NA_character_
        ca0$classe <- trimws(as.character(ca0$classe))
        ca0 <- ca0 %>% dplyr::filter(nzchar(.data$rsid), is.finite(.data$chr), is.finite(.data$position)) %>%
          dplyr::distinct(.data$rsid, .data$chr, .data$position, .keep_all = TRUE)
        return(ca0)
      }

      # Otherwise, try to parse with flexible column detection
      prepare_candidates(ca0)
    })
    
    # --- Refresh chromosome + cluster selectors when clusters become available ---
    # Quan canvien els clusters: omple chr i cluster_id SENSE dependre de candidates
#   shiny::observeEvent(clusters_df(), {
#     cl <- clusters_df()
#     shiny::req(is.data.frame(cl), nrow(cl) > 0)
#     
#     # Assegura columnes canòniques (per si prepare_clusters encara no les ha deixat bé)
#     if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- cl$CHR
#     if (!"cluster_id" %in% names(cl) && "cluster" %in% names(cl)) cl$cluster_id <- as.character(cl$cluster)
#     
#     # Si encara falta alguna columna, no facis res (evites crash)
#     shiny::req(all(c("chr","cluster_id") %in% names(cl)))
#     
#     # chr disponibles, però robust (sempre numeric + sense NA)
#     chrs <- suppressWarnings(as.integer(cl$chr))
#     chrs <- sort(unique(chrs[is.finite(chrs)]))
#     if (length(chrs) == 0) return()
#     
#     # selecció chr robusta (mai NA a l'if)
#     chr_cur <- suppressWarnings(as.integer(input$ld_chr))
#     ok_chr  <- (length(chr_cur) == 1) && is.finite(chr_cur) && (chr_cur %in% chrs)
#     chr_sel <- if (ok_chr) chr_cur else chrs[1]
#     
#     shiny::updateSelectInput(
#       session,
#       "ld_chr",
#       choices  = as.character(chrs),
#       selected = as.character(chr_sel)
#     )
#     
#     # Clusters per aquest chr (robust: compara sempre com integer)
#     cl_chr_num <- suppressWarnings(as.integer(cl$chr))
#     cl_chr <- cl[is.finite(cl_chr_num) & (cl_chr_num == chr_sel), , drop = FALSE]
#     
#     ids <- unique(as.character(cl_chr$cluster_id))
#     ids <- ids[!is.na(ids) & nzchar(ids)]
#     
#     if (length(ids) == 0) {
#       shiny::updateSelectInput(session, "ld_cluster_id",
#                                choices = character(0),
#                                selected = character(0))
#       return()
#     }
#     
#     # selecció cluster robusta
#     id_cur <- as.character(input$ld_cluster_id %||% "")
#     ok_id  <- nzchar(id_cur) && (id_cur %in% ids)
#     id_sel <- if (ok_id) id_cur else ids[1]
#     
#     shiny::updateSelectInput(
#       session,
#       "ld_cluster_id",
#       choices  = ids,
#       selected = id_sel
#     )
#     
#   }, ignoreInit = TRUE)
    
    # --- Refresh chromosome selector when clusters become available ---
    shiny::observeEvent(clusters_df(), {
      cl <- clusters_df()
      shiny::req(is.data.frame(cl), nrow(cl) > 0)
      
      if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- cl$CHR
      shiny::req("chr" %in% names(cl))
      
      chrs <- suppressWarnings(as.integer(cl$chr))
      chrs <- sort(unique(chrs[is.finite(chrs)]))
      if (length(chrs) == 0) return()
      
      chr_cur <- suppressWarnings(as.integer(input$ld_chr))
      ok_chr  <- (length(chr_cur) == 1) && is.finite(chr_cur) && (chr_cur %in% chrs)
      chr_sel <- if (ok_chr) chr_cur else chrs[1]
      
      shiny::updateSelectInput(
        session,
        "ld_chr",
        choices  = as.character(chrs),
        selected = as.character(chr_sel)
      )
    }, ignoreInit = FALSE)
    
    
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
    
 
    # cluster selector
    shiny::observeEvent(list(clusters_df(), input$ld_chr), {
      cl <- clusters_df()
      chr_sel <- suppressWarnings(as.integer(input$ld_chr))
      if (length(chr_sel) != 1 || is.na(chr_sel) || !is.finite(chr_sel)) {
        chr_sel <- sort(unique(suppressWarnings(as.integer(cl$chr))))[1]
      }
      
      cl_chr <- cl %>%
        dplyr::filter(.data$chr == chr_sel) %>%
        dplyr::arrange(start, end)
      
      shiny::validate(shiny::need(nrow(cl_chr) > 0, "No clusters for selected chromosome."))
      
      labs <- paste0(
        cl_chr$cluster_id,
        " (chr", chr_label_plink(cl_chr$chr), ":", cl_chr$start, "-", cl_chr$end, ")",
        " (", cl_chr$size_kb, " kb)"
      )
      choices <- stats::setNames(as.character(cl_chr$cluster_key), labs)
      
      cur <- input$ld_cluster_id
      cur <- if (length(cur) == 1) as.character(cur) else ""
      if (is.na(cur) || !nzchar(cur) || !(cur %in% cl_chr$cluster_key)) {
        cur <- as.character(cl_chr$cluster_key[1])
      }
      shiny::updateSelectInput(session, "ld_cluster_id", choices = choices, selected = cur)
    }, ignoreInit = FALSE)
    
#   selected_cluster <- shiny::reactive({
#     shiny::req(input$ld_chr, input$ld_cluster_id)
#     
#     cl <- clusters_df()
#     shiny::validate(shiny::need(is.data.frame(cl) && nrow(cl) > 0, "No clusters loaded."))
#     
#     chr_sel <- suppressWarnings(as.integer(input$ld_chr))
#     key_sel <- as.character(input$ld_cluster_id)
#     
#     one <- cl %>%
#       dplyr::filter(
#         suppressWarnings(as.integer(.data$chr)) == chr_sel,
#         as.character(.data$cluster_key) == key_sel
#       )
#     
#     shiny::validate(shiny::need(nrow(one) >= 1, "Cluster interval not found."))
#     
#     if (nrow(one) > 1) {
#       one <- one %>% dplyr::arrange(.data$start, .data$end) %>% dplyr::slice(1)
#     }
#     
#     one
#   })
 
    selected_cluster <- shiny::reactive({
      shiny::req(input$ld_chr, input$ld_cluster_id)
      cl <- clusters_df()
      
      one <- cl %>%
        dplyr::filter(
          suppressWarnings(as.integer(.data$chr)) == suppressWarnings(as.integer(input$ld_chr)),
          .data$cluster_key == as.character(input$ld_cluster_id)
        )
      
      shiny::validate(shiny::need(nrow(one) >= 1, "Cluster interval not found."))
      if (nrow(one) > 1) one <- one %>% dplyr::arrange(.data$start, .data$end) %>% dplyr::slice(1)
      one
    })
    
    
    # LD state
    ld_state <- shiny::reactiveValues(
      fl = NULL,
      ld_pairs = NULL,
      blocks_ij = NULL,
      subset_prefix = NULL,
      tag = NULL,
      chr_sel = NULL,
      st = NULL,
      en = NULL
    )
    
    # -----------------------------
    # Build LD (PLINK subset + LD)
    # -----------------------------
    shiny::observeEvent(input$run_ld_cluster, {
      ld_log("")
      append_log("[LD] Starting…")
      
      tryCatch({
        shiny::withProgress(message = "Building LD…", value = 0, {
          
          ld_state$fl <- NULL
          ld_state$ld_pairs <- NULL
          ld_state$blocks_ij <- NULL
          ld_state$subset_prefix <- NULL
          ld_state$tag <- NULL
          ld_state$chr_sel <- NULL
          ld_state$st <- NULL
          ld_state$en <- NULL
          
          plink_bin <- input$plink_bin %||% ""
          bfile_ref <- input$bfile_ref_prefix %||% ""
          
          shiny::validate(shiny::need(file.exists(plink_bin), paste0("PLINK binary not found: ", plink_bin)))
          shiny::validate(shiny::need(nzchar(bfile_ref), "Reference bfile prefix is empty."))
          shiny::validate(shiny::need(file.exists(paste0(bfile_ref, ".bed")) &&
                                        file.exists(paste0(bfile_ref, ".bim")) &&
                                        file.exists(paste0(bfile_ref, ".fam")),
                                      "Reference bfile prefix invalid (missing .bed/.bim/.fam)."))
          
          cl <- selected_cluster()
          chr_sel <- suppressWarnings(as.integer(cl$chr[1]))
          st      <- suppressWarnings(as.integer(cl$start[1]))
          en      <- suppressWarnings(as.integer(cl$end[1]))
          shiny::validate(shiny::need(is.finite(chr_sel) && is.finite(st) && is.finite(en) && st < en,
                                      "Invalid cluster interval (chr/start/end)."))
          
          chr_lab <- chr_label_plink(chr_sel)
          append_log(sprintf("[LD] Cluster: %s | chr%s:%d-%d", cl$cluster_id[1], chr_lab, st, en))
          
          shiny::incProgress(0.10, detail = "Resolving keep file…")
          keep_path <- tryCatch(read_keep_file_dir(input$ld_pop, input$ld_pops_dir), error = function(e) NULL)
          shiny::validate(shiny::need(!is.null(keep_path) && file.exists(keep_path), "Keep file not found/invalid."))
          append_log(sprintf("[LD] keep_dir=%s", input$ld_pops_dir %||% ""))
          append_log(sprintf("[LD] pop=%s", input$ld_pop %||% ""))
          append_log(sprintf("[LD] Using keep file: %s", keep_path))
          
          # 1) subset
          shiny::incProgress(0.20, detail = "Subsetting interval…")
          tag <- paste0("ld_", cl$cluster_id[1], "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
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
          
          append_log("[LD] Running PLINK subset…")
          r1 <- run_plink(args_subset, subset_prefix, plink_bin)
          if (length(r1$stdout)) append_log(paste(r1$stdout, collapse = "\n"))
          if (is.null(r1$status) || r1$status != 0) {
            if (length(r1$log_tail)) append_log(paste(r1$log_tail, collapse = "\n"))
            shiny::showNotification("PLINK subset failed. See log.", type = "error", duration = NULL)
            return()
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
          
          # thinning (same as LD_Inspector)
          maxn <- suppressWarnings(as.integer(input$ld_max_snps_interval %||% 400))
          shiny::validate(shiny::need(is.finite(maxn) && maxn >= 2, "Invalid max SNPs limit."))
          
          if (nrow(fl) > maxn) {
            append_log(sprintf("[LD] Too many SNPs: %d > %d. Auto-thinning…", nrow(fl), maxn))
            
            cand_sub <- NULL
            cand <- tryCatch(candidates_df(), error = function(e) NULL)
            if (is.data.frame(cand) && nrow(cand) > 0) {
              cand_sub <- cand %>% dplyr::filter(chr == chr_sel, position >= st, position <= en)
            }
            
            keep_snps <- character(0)
            if (is.data.frame(cand_sub) && nrow(cand_sub) > 0) {
              keep_by_rsid <- intersect(fl$SNP, cand_sub$rsid)
              keep_by_pos  <- fl %>% dplyr::filter(BP %in% cand_sub$position) %>% dplyr::pull(SNP)
              keep_snps <- unique(c(keep_by_rsid, keep_by_pos))
            }
            
            shiny::validate(shiny::need(length(keep_snps) <= maxn,
                                        paste0("Candidates in interval (", length(keep_snps),
                                               ") exceed max limit (", maxn, "). Increase limit.")))
            
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
            snps_final <- snps_final[!is.na(snps_final) & nzchar(snps_final)]
            shiny::validate(shiny::need(length(snps_final) >= 2, "After thinning, <2 SNPs remain."))
            
            extract_file <- file.path(workdir, paste0(tag, "_extract_snps.txt"))
            writeLines(snps_final, extract_file)
            
            subset2_prefix <- file.path(workdir, paste0(tag, "_subset_thin"))
            
            append_log(sprintf("[LD] Creating thinned subset with %d SNPs…", length(snps_final)))
            args_thin <- c("--bfile", subset_prefix, "--extract", extract_file, "--make-bed", "--out", subset2_prefix)
            rthin <- run_plink(args_thin, subset2_prefix, plink_bin)
            if (length(rthin$stdout)) append_log(paste(rthin$stdout, collapse = "\n"))
            if (is.null(rthin$status) || rthin$status != 0) {
              if (length(rthin$log_tail)) append_log(paste(rthin$log_tail, collapse = "\n"))
              shiny::showNotification("PLINK thinning subset failed. See log.", type = "error", duration = NULL)
              return()
            }
            
            subset_prefix <- subset2_prefix
            
            bim_file2 <- paste0(subset_prefix, ".bim")
            shiny::validate(shiny::need(file.exists(bim_file2), "Thinned subset .bim not created."))
            
            bim2 <- utils::read.table(bim_file2, header = FALSE, stringsAsFactors = FALSE)
            colnames(bim2) <- c("CHR","SNP","CM","BP","A1","A2")
            
            fl <- bim2 %>%
              dplyr::transmute(SNP = as.character(SNP), BP = suppressWarnings(as.integer(BP))) %>%
              dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
              dplyr::distinct(SNP, BP, .keep_all = TRUE) %>%
              dplyr::arrange(BP, SNP) %>%
              dplyr::mutate(ix = dplyr::row_number())
            
            append_log(sprintf("[LD] SNPs after thinning: %d", nrow(fl)))
          } else {
            append_log(sprintf("[LD] SNPs in interval: %d", nrow(fl)))
          }
          
          # 2) LD
          shiny::incProgress(0.45, detail = "Computing LD…")
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
          
          append_log("[LD] Running PLINK LD…")
          r2 <- run_plink(r2_args, ld_prefix, plink_bin)
          if (length(r2$stdout)) append_log(paste(r2$stdout, collapse = "\n"))
          if (is.null(r2$status) || r2$status != 0) {
            if (length(r2$log_tail)) append_log(paste(r2$log_tail, collapse = "\n"))
            shiny::showNotification("PLINK LD failed. See log.", type = "error", duration = NULL)
            return()
          }
          
          shiny::incProgress(0.60, detail = "Parsing LD output…")
          ld <- tryCatch(read_plink_ld(ld_prefix), error = function(e) NULL)
          shiny::validate(shiny::need(is.data.frame(ld) && nrow(ld) > 0, "LD output empty/unreadable."))
          
          snpA <- pick_col(ld, c("SNP_A","SNP_A1","SNP1"))
          snpB <- pick_col(ld, c("SNP_B","SNP_B1","SNP2"))
          shiny::validate(shiny::need(!is.null(snpA) && !is.null(snpB), "LD file missing SNP_A / SNP_B."))
          
          if (identical(input$ld_metric, "Dprime")) {
            vcol <- pick_col(ld, c("Dprime","D'","DP","DPRIME"))
          } else {
            vcol <- pick_col(ld, c("R2","r2"))
          }
          shiny::validate(shiny::need(!is.null(vcol), "LD file missing requested metric column."))
          
          ld2 <- ld %>%
            dplyr::transmute(
              SNP_A = as.character(.data[[snpA]]),
              SNP_B = as.character(.data[[snpB]]),
              value = suppressWarnings(as.numeric(.data[[vcol]]))
            ) %>%
            dplyr::filter(!is.na(SNP_A), !is.na(SNP_B), is.finite(value))
          
          fl_map <- fl %>% dplyr::select(SNP, BP, ix)
          
          ld2 <- ld2 %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixA = ix, BPA = BP), by = c("SNP_A" = "SNP")) %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixB = ix, BPB = BP), by = c("SNP_B" = "SNP")) %>%
            dplyr::filter(is.finite(ixA), is.finite(ixB), ixA != ixB) %>%
            dplyr::mutate(i = pmax(ixA, ixB), j = pmin(ixA, ixB)) %>%
            dplyr::distinct(i, j, .keep_all = TRUE) %>%
            dplyr::filter(i > j)
          
          shiny::validate(shiny::need(nrow(ld2) > 0, "No LD pairs after filtering."))
          ld_state$ld_pairs <- ld2
          
          # 3) candidate annotation (same as LD_Inspector)
          shiny::incProgress(0.85, detail = "Annotating candidates…")
          cand <- tryCatch(candidates_df(), error = function(e) NULL)
          
          if (!is.data.frame(cand) || nrow(cand) == 0) {
            append_log("[LD] Candidates empty. Skipping annotation.")
            fl2 <- fl %>% dplyr::mutate(is_candidate = FALSE, classe = NA_character_) %>%
              dplyr::select(SNP, BP, ix, is_candidate, classe)
          } else {
            cand_sub <- cand %>% dplyr::filter(chr == chr_sel, position >= st, position <= en)
            
            fl2 <- fl %>%
              dplyr::left_join(
                cand_sub %>%
                  dplyr::select(rsid, position, classe) %>%
                  dplyr::rename(cand_pos = position, cand_class = classe),
                by = c("SNP" = "rsid")
              ) %>%
              dplyr::mutate(is_candidate = !is.na(cand_pos) | (BP %in% cand_sub$position))
            
            if (nrow(cand_sub) > 0) {
              pos_class <- cand_sub %>%
                dplyr::select(position, classe) %>%
                dplyr::distinct() %>%
                dplyr::rename(BP = position, classe_pos = classe)
              
              fl2 <- fl2 %>%
                dplyr::left_join(pos_class, by = "BP") %>%
                dplyr::mutate(classe = dplyr::coalesce(cand_class, classe_pos)) %>%
                dplyr::select(SNP, BP, ix, is_candidate, classe)
            } else {
              fl2 <- fl2 %>%
                dplyr::mutate(classe = cand_class) %>%
                dplyr::select(SNP, BP, ix, is_candidate, classe)
            }
          }
          
          ld_state$fl <- fl2
          append_log(sprintf("[LD] Candidate SNPs in interval: %d", sum(fl2$is_candidate, na.rm = TRUE)))
          
          # save context for blocks
          ld_state$subset_prefix <- subset_prefix
          ld_state$tag <- tag
          ld_state$chr_sel <- chr_sel
          ld_state$st <- st
          ld_state$en <- en
          
          shiny::incProgress(1, detail = "Done.")
          shiny::showNotification("LD computed. (Blocks can be run separately.)", type = "message", duration = 3)
        })
      }, error = function(e) {
        append_log(paste0("[LD] ERROR: ", (e$message %||% "") %||% conditionMessage(e)))
        shiny::showNotification(paste("Build LD error:", e$message), type = "error", duration = NULL)
      })
    }, ignoreInit = TRUE)
    
    # -----------------------------
    # Build BLOCKS (PLINK --blocks)
    # IMPORTANT: plink1.9 modifiers are tokens AFTER --blocks (NO --blocks-no-pheno-req)
    # -----------------------------
    shiny::observeEvent(input$run_ld_blocks, {
      append_log("[BLOCKS] Starting…")
      
      tryCatch({
        shiny::withProgress(message = "Building LD blocks…", value = 0, {
          
          shiny::incProgress(0.10, detail = "Checking prerequisites…")
          
          plink_bin <- input$plink_bin %||% ""
          shiny::validate(shiny::need(file.exists(plink_bin), paste0("PLINK binary not found: ", plink_bin)))
          
          subset_prefix <- ld_state$subset_prefix
          shiny::validate(shiny::need(!is.null(subset_prefix) && nzchar(subset_prefix), "Run LD first (no subset available)."))
          shiny::validate(shiny::need(file.exists(paste0(subset_prefix, ".bed")) &&
                                        file.exists(paste0(subset_prefix, ".bim")) &&
                                        file.exists(paste0(subset_prefix, ".fam")),
                                      "LD subset files not found. Run LD again."))
          
          fl <- ld_state$fl
          shiny::validate(shiny::need(is.data.frame(fl) && nrow(fl) >= 2, "Run LD first (SNP list missing)."))
          
          ld_state$blocks_ij <- NULL
          
          shiny::incProgress(0.30, detail = "Running PLINK --blocks…")
          
          tag <- ld_state$tag %||% format(Sys.time(), "%Y%m%d_%H%M%S")
          blk_prefix <- file.path(workdir, paste0(tag, "_blocks"))
          
          blk_args <- c(
            "--bfile", subset_prefix,
            "--blocks", "no-pheno-req",
            if (isTRUE(input$ld_blocks_no_small_max_span)) "no-small-max-span" else NULL,
            "--blocks-inform-frac",  as.character(input$ld_blocks_inform_frac %||% 0.60),
            "--blocks-max-kb",       as.character(input$ld_blocks_max_kb %||% 1000),
            "--blocks-min-maf",      as.character(input$ld_blocks_min_maf %||% 0.05),
            "--blocks-strong-highci", as.character(input$ld_blocks_strong_highci %||% 0.90),
            "--blocks-strong-lowci",  as.character(input$ld_blocks_strong_lowci %||% 0.55),
            "--out", blk_prefix
          )
          blk_args <- blk_args[!is.null(blk_args) & nzchar(as.character(blk_args))]
          
          append_log("[BLOCKS] Running PLINK --blocks…")
          rb <- run_plink(blk_args, blk_prefix, plink_bin)
          if (length(rb$stdout)) append_log(paste(rb$stdout, collapse = "\n"))
          
          if (is.null(rb$status) || rb$status != 0) {
            append_log("[BLOCKS] --blocks failed.")
            if (length(rb$log_tail)) append_log(paste(rb$log_tail, collapse = "\n"))
            shiny::showNotification("PLINK --blocks failed. See log.", type = "error", duration = NULL)
            return()
          }
          
          shiny::incProgress(0.65, detail = "Parsing blocks…")
          
          br <- tryCatch(read_plink_blocks(blk_prefix), error = function(e) NULL)
          shiny::validate(shiny::need(is.data.frame(br) && nrow(br) > 0, "Blocks file empty/unreadable."))
          
          ij <- tryCatch(blocks_det_to_ij(br, fl %>% dplyr::select(SNP, BP)), error = function(e) NULL)
          shiny::validate(shiny::need(is.data.frame(ij) && nrow(ij) > 0, "Blocks parsed but no valid i/j intervals."))
          
          ld_state$blocks_ij <- ij
          append_log(sprintf("[BLOCKS] Blocks: %d", nrow(ij)))
          
          shiny::incProgress(1, detail = "Done.")
          shiny::showNotification("Blocks computed.", type = "message", duration = 3)
        })
      }, error = function(e) {
        append_log(paste0("[BLOCKS] ERROR: ", e$message))
        shiny::showNotification(paste("Build blocks error:", e$message), type = "error", duration = NULL)
      })
    }, ignoreInit = TRUE)
    
    # tables
    output$ld_fl_dt <- DT::renderDT({
      fl <- ld_state$fl
      if (is.null(fl) || !nrow(fl)) {
        return(DT::datatable(
          data.frame(Message="Run LD to populate SNP list."),
          options = list(dom="t"),
          rownames = FALSE
        ))
      }
      
      show <- fl %>%
        dplyr::mutate(
          candidate = dplyr::if_else(isTRUE(is_candidate) | (is.logical(is_candidate) & is_candidate), "YES", ""),
          classe    = dplyr::coalesce(as.character(classe), "")
        ) %>%
        dplyr::mutate(classe = dplyr::if_else(nzchar(classe), classe, "")) %>%
        dplyr::transmute(ix, SNP, BP, candidate, classe)
      
      DT::datatable(
        show,
        rownames   = FALSE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons    = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 15,
          scrollX    = TRUE
        )
      )
    })
    
    
    output$ld_pairs_dt <- DT::renderDT({
      lp <- ld_state$ld_pairs
      fl <- ld_state$fl
      if (is.null(lp) || !nrow(lp) || is.null(fl) || !nrow(fl)) {
        return(DT::datatable(
          data.frame(Message="Run LD to populate pairs."),
          options = list(dom="t"),
          rownames = FALSE
        ))
      }
      
      fl_small <- fl %>% dplyr::select(ix, SNP, BP, is_candidate, classe) %>%
        dplyr::distinct(ix, .keep_all = TRUE)
      
      lp2 <- lp %>% dplyr::distinct(i, j, .keep_all = TRUE)
      
      top <- lp2 %>%
        dplyr::left_join(
          fl_small %>% dplyr::rename(i = ix, SNP_i = SNP, BP_i = BP, cand_i = is_candidate, class_i = classe),
          by = "i", relationship = "many-to-one"
        ) %>%
        dplyr::left_join(
          fl_small %>% dplyr::rename(j = ix, SNP_j = SNP, BP_j = BP, cand_j = is_candidate, class_j = classe),
          by = "j", relationship = "many-to-one"
        )
      
      DT::datatable(
        top,
        rownames   = FALSE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons    = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 12,
          scrollX    = TRUE
        )
      )
    })
    
    # ----------------------------------------------------------------------
    # PLOT (copied from LD_Inspector; uses findOverlaps; NO subsetByOverlaps)
    # ----------------------------------------------------------------------
    output$ld_plot <- plotly::renderPlotly({
      
      shiny::req(requireNamespace("ggplot2", quietly = TRUE))
      shiny::req(requireNamespace("dplyr", quietly = TRUE))
      shiny::req(requireNamespace("scales", quietly = TRUE))
      shiny::req(requireNamespace("plotly", quietly = TRUE))
      
      g0 <- NULL
      fl0 <- ld_state$fl
      ld  <- ld_state$ld_pairs
      
      shiny::validate(
        shiny::need(is.data.frame(fl0) && nrow(fl0) >= 2, "Run LD first (cluster + population)."),
        shiny::need(is.data.frame(ld)  && nrow(ld)  >  0, "No LD pairs available (run LD first).")
      )
      
      metric <- input$ld_metric %||% "R2"
      x_mode <- input$ld_x_mode %||% "bp"
      shiny::validate(shiny::need(metric %in% c("R2","Dprime"), "ld_metric must be R2 or Dprime."))
      shiny::validate(shiny::need(x_mode %in% c("bp","equal"), "ld_x_mode must be bp or equal."))
      
      cl <- tryCatch(selected_cluster(), error = function(e) NULL)
      chr_sel <- if (!is.null(cl)) suppressWarnings(as.integer(cl$chr[1]))   else NA_integer_
      st      <- if (!is.null(cl)) suppressWarnings(as.integer(cl$start[1])) else NA_integer_
      en      <- if (!is.null(cl)) suppressWarnings(as.integer(cl$end[1]))   else NA_integer_
      
      bp_to_int <- function(x) suppressWarnings(as.integer(x))
      
      # A) fl net + X
      fl <- fl0 %>%
        dplyr::mutate(
          SNP = as.character(SNP),
          BP  = bp_to_int(BP),
          is_candidate = dplyr::coalesce(as.logical(is_candidate), FALSE),
          classe = if ("classe" %in% names(.)) as.character(classe) else ""
        ) %>%
        dplyr::mutate(
          classe = trimws(classe),
          classe = dplyr::na_if(classe, "")
        ) %>%
        dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
        dplyr::distinct(SNP, .keep_all = TRUE) %>%
        dplyr::arrange(BP, SNP)
      
      shiny::validate(shiny::need(nrow(fl) >= 2, "Subset (fl) has <2 SNPs after cleaning."))
      
      if (identical(x_mode, "equal")) {
        fl <- fl %>% dplyr::mutate(X = dplyr::row_number())
      } else {
        fl <- fl %>% dplyr::mutate(X = BP)
      }
      
      x_limits <- range(fl$X, na.rm = TRUE)
      x_span   <- max(1, x_limits[2] - x_limits[1])
      
      xmap <- stats::setNames(fl$X, fl$SNP)
      
      bp_to_x <- function(bp_vec) {
        bp_vec <- suppressWarnings(as.numeric(bp_vec))
        if (!identical(x_mode, "equal")) return(bp_vec)
        stats::approx(x = fl$BP, y = fl$X, xout = bp_vec, rule = 2, ties = "ordered")$y
      }
      
      # B) candidates in interval (outside LD set)
      cand_all <- NULL
      if (is.finite(chr_sel) && is.finite(st) && is.finite(en)) {
        cand <- tryCatch(candidates_df(), error = function(e) NULL)
        if (is.data.frame(cand) && nrow(cand) > 0) {
          cand_all <- cand %>%
            dplyr::mutate(
              rsid   = as.character(rsid),
              chr    = suppressWarnings(as.integer(chr)),
              BP     = bp_to_int(position),
              classe = if ("classe" %in% names(.)) as.character(classe) else ""
            ) %>%
            dplyr::mutate(
              classe = trimws(classe),
              classe = dplyr::na_if(classe, "")
            ) %>%
            dplyr::filter(chr == chr_sel, is.finite(BP), BP >= st, BP <= en, !is.na(rsid), nzchar(rsid)) %>%
            dplyr::distinct(rsid, .keep_all = TRUE)
        }
      }
      
      # B2) enrich fl classe from candidates_df (SNP==rsid)
      if (is.data.frame(cand_all) && nrow(cand_all) > 0) {
        cand_map <- cand_all %>%
          dplyr::transmute(SNP = rsid, classe_cand = dplyr::coalesce(classe, "candidate")) %>%
          dplyr::distinct(SNP, .keep_all = TRUE)
        
        fl <- fl %>%
          dplyr::left_join(cand_map, by = "SNP") %>%
          dplyr::mutate(classe = dplyr::coalesce(classe, classe_cand)) %>%
          dplyr::select(-classe_cand) %>%
          dplyr::mutate(
            classe = trimws(classe),
            classe = dplyr::na_if(classe, "")
          )
      }
      
      # C) classe2 + palette
      fl <- fl %>% dplyr::mutate(classe2 = dplyr::if_else(is_candidate, dplyr::coalesce(classe, "candidate"), "Other"))
      
      # D) ticks + hover
      y_base       <- 0
      y_tick_other <- 0.18
      y_tick_hit   <- 0.28
      
      ticks <- fl %>%
        dplyr::mutate(
          tick_h = dplyr::if_else(is_candidate, y_tick_hit, y_tick_other),
          hover  = paste0("SNP: ", SNP,
                          "<br>BP: ", BP,
                          "<br>candidate: ", is_candidate,
                          ifelse(is_candidate, paste0("<br>classe: ", classe2), ""))
        )
      
      labs_out <- NULL
      if (is.data.frame(cand_all) && nrow(cand_all) > 0) {
        labs_out <- cand_all %>%
          dplyr::filter(!(rsid %in% fl$SNP)) %>%
          dplyr::mutate(SNP = rsid, X = bp_to_x(BP), grp = dplyr::coalesce(classe, "candidate")) %>%
          dplyr::filter(is.finite(X)) %>%
          dplyr::transmute(SNP, BP, X, grp = as.character(grp)) %>%
          dplyr::arrange(X)
      }
      
      ticks_out <- NULL
      if (is.data.frame(labs_out) && nrow(labs_out) > 0) {
        ticks_out <- labs_out %>%
          dplyr::transmute(
            X = X,
            tick_h = y_tick_hit,
            classe2 = grp,
            hover = paste0("Candidate (not in LD set): ", SNP, "<br>BP: ", BP, "<br>classe: ", grp)
          )
      }
      
      cand_classes_all <- sort(unique(c(
        fl %>% dplyr::filter(is_candidate) %>% dplyr::pull(classe2),
        if (is.data.frame(labs_out) && nrow(labs_out) > 0) labs_out$grp else character(0)
      )))
      if (length(cand_classes_all) == 0) cand_classes_all <- "candidate"
      
      pal <- scales::hue_pal()(max(1, length(cand_classes_all)))
      class_cols <- stats::setNames(pal, cand_classes_all)
      cols_all <- c(class_cols, Other = "darkgrey")
      
      x_title <- if (identical(x_mode, "equal")) "Equal spacing" else "Genomic (BP)"
      x_breaks_fun <- scales::pretty_breaks(n = 4)
      x_labels_fun <- if (identical(x_mode, "equal")) {
        function(z) z
      } else {
        scales::label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01)
      }
      
      # label overlap helper (same)
      assign_label_rows_greedy <- function(x, min_dx, n_rows = 6) {
        x <- as.numeric(x)
        last_x <- rep(-Inf, n_rows)
        row_id <- integer(length(x))
        for (i in seq_along(x)) {
          placed <- FALSE
          for (r in seq_len(n_rows)) {
            if (!is.finite(last_x[r]) || (x[i] - last_x[r] >= min_dx)) {
              row_id[i] <- r
              last_x[r] <- x[i]
              placed <- TRUE
              break
            }
          }
          if (!placed) {
            r <- which.min(last_x)
            row_id[i] <- r
            last_x[r] <- x[i]
          }
        }
        row_id
      }
      
      n_rows_labels <- 6
      min_dx <- if (identical(x_mode, "equal")) 1 else x_span * 0.015
      y0_lab <- 0.42
      dy_lab <- 0.075
      
      add_y_lab_by_grp <- function(df) {
        if (!is.data.frame(df) || nrow(df) == 0) return(df)
        df %>%
          dplyr::group_by(grp) %>%
          dplyr::arrange(X, .by_group = TRUE) %>%
          dplyr::group_modify(~{
            rid <- assign_label_rows_greedy(.x$X, min_dx = min_dx, n_rows = n_rows_labels)
            .x$y_lab <- y0_lab + (rid - 1) * dy_lab
            .x
          }) %>%
          dplyr::ungroup()
      }
      
      labs_in <- fl %>%
        dplyr::filter(is_candidate) %>%
        dplyr::transmute(SNP = SNP, BP = BP, X = X, grp = as.character(classe2)) %>%
        dplyr::arrange(grp, X)
      
      labs_in  <- add_y_lab_by_grp(labs_in)
      labs_out <- add_y_lab_by_grp(labs_out)
      
      df_in <- fl %>%
        dplyr::filter(is_candidate) %>%
        dplyr::transmute(
          X = X, grp = as.character(classe2), tick_h = y_tick_hit, kind = "in",
          hover = paste0("Candidate: ", SNP, "<br>BP: ", BP, "<br>classe: ", classe2, "<br>in LD set: YES")
        )
      
      df_out <- NULL
      if (!is.null(ticks_out) && nrow(ticks_out) > 0) {
        df_out <- ticks_out %>%
          dplyr::transmute(X = X, grp = as.character(classe2), tick_h = tick_h, kind = "out", hover = hover)
      }
      
      df_track <- dplyr::bind_rows(df_in, df_out) %>%
        dplyr::filter(!is.na(grp), nzchar(grp), is.finite(X)) %>%
        dplyr::mutate(grp = factor(grp, levels = cand_classes_all),
                      lty = ifelse(kind == "out", "dashed", "solid"))
      
      shiny::validate(shiny::need(nrow(df_track) > 0, "No candidate tracks to plot."))
      
      max_lab_chars <- max(nchar(cand_classes_all), na.rm = TRUE)
      x_pad <- if (identical(x_mode, "equal")) {
        max(3, ceiling(max_lab_chars * 0.9))
      } else {
        x_span * min(0.20, max(0.06, max_lab_chars / 80))
      }
      x_limits2 <- c(x_limits[1], x_limits[2] + x_pad)
      
      base_df <- data.frame(
        grp = factor(cand_classes_all, levels = cand_classes_all),
        x1 = x_limits2[1], x2 = x_limits2[2], y = y_base
      )
      
      p_tracks <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = base_df,
          ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
          inherit.aes = FALSE,
          linewidth = 1.15
        ) +
        ggplot2::geom_segment(
          data = df_track,
          ggplot2::aes(x = X, xend = X, y = y_base, yend = y_base + tick_h, color = grp, linetype = lty),
          linewidth = 0.40, alpha = 0.90
        ) +
        ggplot2::geom_point(
          data = df_track,
          ggplot2::aes(x = X, y = y_base + tick_h, color = grp, text = hover),
          alpha = 0.01, size = 2
        ) +
        ggplot2::facet_grid(rows = ggplot2::vars(grp), switch = "y") +
        ggplot2::scale_x_continuous(limits = x_limits2, breaks = x_breaks_fun, labels = x_labels_fun, expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(-0.05, 0.90), clip = "off") +
        ggplot2::scale_color_manual(values = cols_all, guide = "none") +
        ggplot2::scale_linetype_identity() +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          strip.text.y.left = ggplot2::element_blank(),
          strip.background  = ggplot2::element_blank(),
          strip.placement   = "outside",
          axis.title.y = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 11),
          axis.text.x  = ggplot2::element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = ggplot2::element_line(),
          plot.margin  = ggplot2::margin(t = 6, r = 18, b = 6, l = 6)
        ) +
        ggplot2::labs(x = x_title)
      
      g1 <- plotly::ggplotly(p_tracks, tooltip = "text", source = "ld") %>%
        plotly::layout(
          margin = list(l = 80, r = 140, t = 10, b = 55),
          hovermode = "closest",
          xaxis = list(range = x_limits2, tickangle = 45, ticklen = 4, tickwidth = 1, ticks = "outside", automargin = TRUE)
        )
      
      # remove y titles created by facets + add colored facet labels as annotations (same logic)
      lay_names <- names(g1$x$layout)
      ykeys <- lay_names[grepl("^yaxis(\\d+)?$", lay_names)]
      for (yk in ykeys) {
        g1$x$layout[[yk]]$title <- list(text = "")
        g1$x$layout[[yk]]$showticklabels <- FALSE
        g1$x$layout[[yk]]$ticks <- ""
        g1$x$layout[[yk]]$zeroline <- FALSE
      }
      
      # remove strip annotations if present
      if (!is.null(g1$x$layout$annotations) && length(g1$x$layout$annotations) > 0) {
        ann <- g1$x$layout$annotations
        ann_keep <- Filter(function(a) {
          tx <- as.character(a$text %||% "")
          !(tx %in% cand_classes_all)
        }, ann)
        g1 <- g1 %>% plotly::layout(annotations = ann_keep)
      }
      
      get_dom <- function(fig, k) {
        d <- fig$x$layout[[k]]$domain
        if (length(d) != 2) return(c(NA_real_, NA_real_))
        c(as.numeric(d[[1]]), as.numeric(d[[2]]))
      }
      yinfo <- data.frame(
        key = ykeys,
        d0  = vapply(ykeys, function(k) get_dom(g1, k)[1], numeric(1)),
        d1  = vapply(ykeys, function(k) get_dom(g1, k)[2], numeric(1)),
        stringsAsFactors = FALSE
      )
      yinfo <- yinfo[is.finite(yinfo$d0) & is.finite(yinfo$d1), , drop = FALSE]
      yinfo <- yinfo[order((yinfo$d0 + yinfo$d1)/2, decreasing = TRUE), , drop = FALSE]
      
      n_map <- min(length(cand_classes_all), nrow(yinfo))
      axis_map <- data.frame(grp = cand_classes_all[seq_len(n_map)],
                             yaxis_key = yinfo$key[seq_len(n_map)],
                             stringsAsFactors = FALSE)
      
      ykey_to_yref <- function(ykey) sub("^yaxis", "y", ykey)
      ykey_to_xref <- function(fig, ykey) {
        anc <- fig$x$layout[[ykey]]$anchor
        if (is.null(anc) || !nzchar(as.character(anc))) "x" else as.character(anc)
      }
      
      ylim0 <- c(-0.05, 0.90)
      y_name_data <- y_base + 0.06
      frac_name <- (y_name_data - ylim0[1]) / (ylim0[2] - ylim0[1])
      frac_name <- max(0.02, min(0.98, frac_name))
      
      name_anns <- lapply(seq_len(nrow(axis_map)), function(i) {
        grp_i <- axis_map$grp[i]
        yk <- axis_map$yaxis_key[i]
        dom <- get_dom(g1, yk)
        y_paper <- dom[1] + frac_name * (dom[2] - dom[1])
        list(
          x = 0.992, xref = "paper",
          y = y_paper, yref = "paper",
          text = grp_i,
          showarrow = FALSE,
          xanchor = "right",
          yanchor = "bottom",
          font = list(size = 12, color = cols_all[[grp_i]] %||% "black")
        )
      })
      
      build_snp_label_anns <- function(df, alpha = 1, italic = FALSE, per_grp_max = 25, textangle = 10) {
        if (!is.data.frame(df) || nrow(df) == 0) return(list())
        df2 <- df %>%
          dplyr::filter(is.finite(X), is.finite(y_lab), !is.na(grp), nzchar(grp)) %>%
          dplyr::group_by(grp) %>% dplyr::slice_head(n = per_grp_max) %>% dplyr::ungroup()
        
        out <- list()
        for (k in seq_len(nrow(df2))) {
          grp_k <- as.character(df2$grp[k])
          map_k <- axis_map[axis_map$grp == grp_k, , drop = FALSE]
          if (nrow(map_k) != 1) next
          yk <- map_k$yaxis_key
          xref_k <- ykey_to_xref(g1, yk)
          yref_k <- ykey_to_yref(yk)
          
          col_k <- cols_all[[grp_k]] %||% "black"
          if (alpha < 1) {
            rgb <- grDevices::col2rgb(col_k)
            col_k <- sprintf("rgba(%d,%d,%d,%g)", rgb[1,1], rgb[2,1], rgb[3,1], alpha)
          }
          
          txt <- df2$SNP[k]
          if (italic) txt <- paste0("<i>", txt, "</i>")
          
          out[[length(out) + 1]] <- list(
            x = df2$X[k], y = df2$y_lab[k],
            xref = xref_k, yref = yref_k,
            text = txt,
            showarrow = FALSE,
            textangle = textangle,
            xanchor = "left",
            yanchor = "middle",
            font = list(size = 10, color = col_k)
          )
        }
        out
      }
      
      anns_in  <- build_snp_label_anns(labs_in,  alpha = 1.0,  italic = FALSE, per_grp_max = 30, textangle = 10)
      anns_out <- build_snp_label_anns(labs_out, alpha = 0.75, italic = TRUE,  per_grp_max = 20, textangle = 10)
      
      old_anns <- g1$x$layout$annotations %||% list()
      g1 <- g1 %>% plotly::layout(annotations = c(old_anns, name_anns, anns_in, anns_out))
      
      # E) LD triangle (NO TOCAR)
      idxmap <- stats::setNames(seq_len(nrow(fl)), fl$SNP)
      bpmap  <- stats::setNames(fl$BP, fl$SNP)
      
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
          y = if (identical(x_mode, "equal")) -(abs(j - i) / max(1, nrow(fl) - 1)) else -(abs(bpB - bpA) / bp_span),
          hover = paste0("SNP_A: ", SNP_A,
                         "<br>SNP_B: ", SNP_B,
                         "<br>", metric, ": ", signif(val, 3),
                         "<br>BP_A: ", bpA,
                         "<br>BP_B: ", bpB)
        ) %>%
        dplyr::filter(is.finite(x_mid), is.finite(y), X_A <= X_B) %>%
        dplyr::select(x_mid, y, val, hover)
      
      shiny::validate(shiny::need(nrow(tri) > 0, "No LD pairs to plot (maybe all NA?)."))
      
      sq_size <- { n <- nrow(fl); if (n <= 60) 3.2 else if (n <= 120) 2.6 else 2.2 }
      
      p2 <- ggplot2::ggplot(tri, ggplot2::aes(x = x_mid, y = y, fill = val, text = hover)) +
        ggplot2::geom_point(shape = 22, size = sq_size, stroke = 0) +
        ggplot2::scale_fill_gradientn(
          colours = c("white", "orange", "red"),
          values  = scales::rescale(c(0, 0.2, 1)),
          limits  = c(0, 1),
          name    = metric
        ) +
        ggplot2::scale_x_continuous(limits = x_limits2, expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(-1.02, 0.02), clip = "off") +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right",
          plot.margin = ggplot2::margin(t = 18, r = 6, b = 6, l = 6)
        )
      
      # F) overlay blocks (optional)
      ij0 <- ld_state$blocks_ij
      if (is.data.frame(ij0) && nrow(ij0) > 0 && all(c("i","j") %in% names(ij0))) {
        
        # --- ENSURE fl has the columns used by block overlay ---
        # fl must have BP and an X coordinate (axis position)
        if (!"BP" %in% names(fl)) {
          # fallback si algun dia el df ve amb POS en lloc de BP
          if ("POS" %in% names(fl)) fl$BP <- fl$POS
        }
        if (!"X" %in% names(fl) || all(!is.finite(suppressWarnings(as.numeric(fl$X))))) {
          fl$X <- if (identical(x_mode, "equal")) seq_len(nrow(fl)) else fl$BP
        }
        # bp_span ha d'existir i ser >0 per evitar yA = -Inf/NA
        if (!exists("bp_span", inherits = TRUE) || !is.finite(bp_span) || bp_span <= 0) {
          bp_span <- max(1, max(fl$BP, na.rm = TRUE) - min(fl$BP, na.rm = TRUE))
        }
        # -------------------------------------------------------
        
        bd <- ij0 %>%
          dplyr::mutate(i = suppressWarnings(as.integer(i)),
                        j = suppressWarnings(as.integer(j))) %>%
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
            
            y0 <- -0.01
            if (identical(x_mode, "equal")) {
              denom <- max(1, nrow(fl) - 1)
              bd <- bd %>% dplyr::mutate(yA = -(abs(j - i) / denom))
            } else {
              bd <- bd %>% dplyr::mutate(yA = -(abs(bpR - bpL) / bp_span))
            }
            bd <- bd %>% dplyr::mutate(yA = pmin(yA, y0 - 0.02))
            
            p2 <- p2 +
              ggplot2::geom_segment(data = bd, ggplot2::aes(x = xL, xend = xR, y = y0, yend = y0),
                                    inherit.aes = FALSE, linewidth = 0.7, colour = "black", alpha = 0.9) +
              ggplot2::geom_segment(data = bd, ggplot2::aes(x = xL, xend = xM, y = y0, yend = yA),
                                    inherit.aes = FALSE, linewidth = 0.7, colour = "black", alpha = 0.9) +
              ggplot2::geom_segment(data = bd, ggplot2::aes(x = xR, xend = xM, y = y0, yend = yA),
                                    inherit.aes = FALSE, linewidth = 0.7, colour = "black", alpha = 0.9) +
              ggplot2::coord_cartesian(ylim = c(-1.02, 0.05), clip = "off")
          }
        }
      }
      
      g2 <- plotly::ggplotly(p2, tooltip = "text", source = "ld") %>%
        plotly::layout(margin = list(l = 40, r = 10, t = 10, b = 10))
      
      # --- FIX: overlay blocks as plotly segments (always on top) ---
      if (exists("bd") && is.data.frame(bd) && nrow(bd) > 0 &&
          all(c("xL","xR","xM","yA") %in% names(bd))) {
        
        y0 <- -0.01
        
        g2 <- g2 %>%
          plotly::add_segments(
            x = bd$xL, xend = bd$xR,
            y = rep(y0, nrow(bd)), yend = rep(y0, nrow(bd)),
            inherit = FALSE, showlegend = FALSE, hoverinfo = "skip",
            line = list(width = 2)
          ) %>%
          plotly::add_segments(
            x = bd$xL, xend = bd$xM,
            y = rep(y0, nrow(bd)), yend = bd$yA,
            inherit = FALSE, showlegend = FALSE, hoverinfo = "skip",
            line = list(width = 2)
          ) %>%
          plotly::add_segments(
            x = bd$xR, xend = bd$xM,
            y = rep(y0, nrow(bd)), yend = bd$yA,
            inherit = FALSE, showlegend = FALSE, hoverinfo = "skip",
            line = list(width = 2)
          )
      }
      
      ############# G) gene model (findOverlaps; no subsetByOverlaps)
      # ------------------------------------------------------------
      # G) Gene model (TOP) richer: backbone + introns + exons + labels + hover
      #    (copiat del LD_Inspector)
      # ------------------------------------------------------------
      
      # TxDb (hg38) com a LD_Inspector
      txdb <- NULL
      if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      }
      
      # helper: greedy lane assignment to avoid overlaps
      assign_lanes_greedy <- function(starts, ends) {
        starts <- as.numeric(starts)
        ends   <- as.numeric(ends)
        o <- order(starts, ends, na.last = TRUE)
        lane_end <- numeric(0)
        lane <- rep(NA_integer_, length(starts))
        
        for (ii in o) {
          s <- starts[ii]; e <- ends[ii]
          if (!is.finite(s) || !is.finite(e)) next
          
          placed <- FALSE
          if (length(lane_end) > 0) {
            for (k in seq_along(lane_end)) {
              if (s > lane_end[k]) { # non-overlapping
                lane[ii] <- k
                lane_end[k] <- e
                placed <- TRUE
                break
              }
            }
          }
          if (!placed) {
            lane_end <- c(lane_end, e)
            lane[ii] <- length(lane_end)
          }
        }
        lane
      }
      
      # helper: SYMBOL/GENENAME mapping (ENTREZ or ENSEMBL)
      map_gene_names <- function(gene_ids) {
        
        out <- data.frame(
          gene_id = as.character(gene_ids),
          stringsAsFactors = FALSE
        )
        
        if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
            !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
          out$SYMBOL   <- NA_character_
          out$GENENAME <- NA_character_
          out$gene_name <- out$gene_id
          return(out)
        }
        
        ids <- unique(out$gene_id)
        ids <- ids[!is.na(ids) & nzchar(ids)]
        
        if (length(ids) == 0) {
          out$SYMBOL   <- NA_character_
          out$GENENAME <- NA_character_
          out$gene_name <- out$gene_id
          return(out)
        }
        
        is_ensg    <- any(grepl("^ENSG", ids, ignore.case = TRUE))
        is_numeric <- all(grepl("^[0-9]+$", ids))
        
        if (is_ensg) {
          
          keys <- unique(sub("\\..*$", "", ids))  # remove ENSG version suffix
          
          mp <- tryCatch(
            suppressMessages(
              AnnotationDbi::select(
                org.Hs.eg.db::org.Hs.eg.db,
                keys    = keys,
                keytype = "ENSEMBL",
                columns = c("SYMBOL", "GENENAME")
              )
            ),
            error = function(e) NULL
          )
          
          if (is.data.frame(mp) && nrow(mp) > 0) {
            mp <- mp %>%
              dplyr::mutate(gene_id_clean = as.character(ENSEMBL)) %>%
              dplyr::select(gene_id_clean, SYMBOL, GENENAME) %>%
              dplyr::arrange(gene_id_clean, dplyr::desc(!is.na(SYMBOL) & nzchar(SYMBOL))) %>%
              dplyr::distinct(gene_id_clean, .keep_all = TRUE)
            
            out <- out %>%
              dplyr::mutate(gene_id_clean = sub("\\..*$", "", gene_id)) %>%
              dplyr::left_join(mp, by = "gene_id_clean") %>%
              dplyr::select(-gene_id_clean)
          }
          
        } else if (is_numeric) {
          
          keys <- unique(ids)
          
          mp <- tryCatch(
            suppressMessages(
              AnnotationDbi::select(
                org.Hs.eg.db::org.Hs.eg.db,
                keys    = keys,
                keytype = "ENTREZID",
                columns = c("SYMBOL", "GENENAME")
              )
            ),
            error = function(e) NULL
          )
          
          if (is.data.frame(mp) && nrow(mp) > 0) {
            mp <- mp %>%
              dplyr::mutate(gene_id = as.character(ENTREZID)) %>%
              dplyr::select(gene_id, SYMBOL, GENENAME) %>%
              dplyr::arrange(gene_id, dplyr::desc(!is.na(SYMBOL) & nzchar(SYMBOL))) %>%
              dplyr::distinct(gene_id, .keep_all = TRUE)
            
            out <- out %>%
              dplyr::left_join(mp, by = "gene_id")
          }
          
        } else {
          out$SYMBOL   <- NA_character_
          out$GENENAME <- NA_character_
        }
        
        if (!"SYMBOL" %in% names(out))   out$SYMBOL <- NA_character_
        if (!"GENENAME" %in% names(out)) out$GENENAME <- NA_character_
        
        out <- out %>%
          dplyr::mutate(
            gene_name = dplyr::coalesce(SYMBOL, gene_id),
            gene_name = dplyr::if_else(is.na(gene_name) | !nzchar(gene_name), gene_id, gene_name)
          )
        
        out
      }
      
      # Funció que construeix el plotly del gene-track (retorna NULL si no hi ha gens/txdb)
      make_gene_track_plot <- function(chr_sel, st, en, bp_to_x, x_limits2) {
        
        if (is.null(txdb)) return(NULL)
        if (!requireNamespace("GenomicRanges", quietly = TRUE)) return(NULL)
        if (!requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
        if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) return(NULL)
        if (!requireNamespace("S4Vectors", quietly = TRUE)) return(NULL)
        
        chr_ucsc <- paste0("chr", chr_label_plink(chr_sel))
        
        region_gr <- GenomicRanges::GRanges(
          seqnames = chr_ucsc,
          ranges   = IRanges::IRanges(start = st, end = en)
        )
        
        genes_gr <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
        GenomeInfoDb::seqlevelsStyle(genes_gr) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
        
        ov <- GenomicRanges::findOverlaps(genes_gr, region_gr, ignore.strand = TRUE)
        keep_idx <- unique(S4Vectors::queryHits(ov))
        if (length(keep_idx) == 0) return(NULL)
        
        gsub_gr <- genes_gr[keep_idx]
        
        gene_id <- S4Vectors::mcols(gsub_gr)$gene_id
        if (is.null(gene_id)) gene_id <- names(gsub_gr)
        gene_id <- as.character(gene_id)
        
        gene_df <- data.frame(
          gene_id = gene_id,
          start   = GenomicRanges::start(gsub_gr),
          end     = GenomicRanges::end(gsub_gr),
          strand  = as.character(GenomicRanges::strand(gsub_gr)),
          stringsAsFactors = FALSE
        )
        gene_df <- gene_df[is.finite(gene_df$start) & is.finite(gene_df$end), , drop = FALSE]
        gene_df <- gene_df[order(gene_df$start, gene_df$end), , drop = FALSE]
        if (nrow(gene_df) == 0) return(NULL)
        
        nm <- map_gene_names(gene_df$gene_id)
        gene_df <- gene_df %>%
          dplyr::left_join(nm %>% dplyr::select(gene_id, SYMBOL, GENENAME, gene_name), by = "gene_id") %>%
          dplyr::mutate(
            gene_name = dplyr::coalesce(gene_name, gene_id),
            SYMBOL    = dplyr::coalesce(SYMBOL, gene_name),
            strand    = dplyr::coalesce(strand, "?")
          )
        
        gene_df$lane <- assign_lanes_greedy(gene_df$start, gene_df$end)
        gene_df$y    <- -as.numeric(gene_df$lane)
        
        # Exons by gene (robust)
        ex_by_gene <- tryCatch(
          suppressMessages(GenomicFeatures::exonsBy(txdb, by = "gene", use.names = FALSE)),
          error = function(e) NULL
        )
        
        ex_df <- NULL
        if (!is.null(ex_by_gene) && length(ex_by_gene) > 0) {
          
          ex_gr_all <- suppressWarnings(unlist(ex_by_gene, use.names = FALSE))
          if (!is.null(ex_gr_all) && length(ex_gr_all) > 0) {
            
            GenomeInfoDb::seqlevelsStyle(ex_gr_all) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
            
            ov_ex_reg <- GenomicRanges::findOverlaps(ex_gr_all, region_gr, ignore.strand = TRUE)
            if (length(ov_ex_reg) > 0) {
              ex_gr_all <- ex_gr_all[unique(S4Vectors::queryHits(ov_ex_reg))]
            } else {
              ex_gr_all <- ex_gr_all[0]
            }
            
            if (length(ex_gr_all) > 0) {
              ov_ex_gene <- GenomicRanges::findOverlaps(ex_gr_all, gsub_gr, ignore.strand = TRUE)
              if (length(ov_ex_gene) > 0) {
                ex_q <- S4Vectors::queryHits(ov_ex_gene)
                ex_s <- S4Vectors::subjectHits(ov_ex_gene)
                
                ex_df <- data.frame(
                  gene_id = gene_id[ex_s],
                  start   = GenomicRanges::start(ex_gr_all)[ex_q],
                  end     = GenomicRanges::end(ex_gr_all)[ex_q],
                  stringsAsFactors = FALSE
                ) %>%
                  dplyr::distinct(gene_id, start, end) %>%
                  dplyr::filter(is.finite(start), is.finite(end)) %>%
                  dplyr::arrange(gene_id, start, end)
              }
            }
          }
        }
        
        if (is.data.frame(ex_df) && nrow(ex_df) > 0) {
          ex_df <- ex_df %>%
            dplyr::left_join(gene_df %>% dplyr::select(gene_id, y, gene_name, SYMBOL, strand), by = "gene_id") %>%
            dplyr::filter(is.finite(y))
        }
        
        # Introns
        intr_df <- NULL
        if (is.data.frame(ex_df) && nrow(ex_df) > 0) {
          intr_list <- lapply(split(ex_df, ex_df$gene_id), function(dd) {
            dd <- dd[order(dd$start, dd$end), , drop = FALSE]
            if (nrow(dd) < 2) return(NULL)
            data.frame(
              gene_id = dd$gene_id[1],
              x0_bp   = dd$end[-nrow(dd)],
              x1_bp   = dd$start[-1],
              y       = dd$y[1],
              gene_name = dd$gene_name[1],
              SYMBOL    = dd$SYMBOL[1],
              strand    = dd$strand[1],
              stringsAsFactors = FALSE
            )
          })
          intr_df <- do.call(rbind, intr_list)
          if (!is.null(intr_df) && nrow(intr_df) > 0) {
            intr_df <- intr_df %>% dplyr::filter(is.finite(x0_bp), is.finite(x1_bp), x1_bp > x0_bp)
          } else {
            intr_df <- NULL
          }
        }
        
        # map to X axis
        gene_df$x0 <- bp_to_x(gene_df$start)
        gene_df$x1 <- bp_to_x(gene_df$end)
        
        if (is.data.frame(ex_df) && nrow(ex_df) > 0) {
          ex_df$x0 <- bp_to_x(ex_df$start)
          ex_df$x1 <- bp_to_x(ex_df$end)
        }
        
        if (!is.null(intr_df) && nrow(intr_df) > 0) {
          intr_df$x0 <- bp_to_x(intr_df$x0_bp)
          intr_df$x1 <- bp_to_x(intr_df$x1_bp)
        }
        
        # hover
        gene_df <- gene_df %>%
          dplyr::mutate(
            hover_gene = paste0(
              "Gene: ", dplyr::coalesce(gene_name, gene_id),
              "<br>SYMBOL: ", dplyr::coalesce(SYMBOL, ""),
              "<br>ID: ", gene_id,
              ifelse(!is.na(GENENAME) & nzchar(GENENAME), paste0("<br>Name: ", GENENAME), ""),
              "<br>Strand: ", strand,
              "<br>Range: ", chr_ucsc, ":", start, "-", end
            )
          )
        
        if (is.data.frame(ex_df) && nrow(ex_df) > 0) {
          ex_df <- ex_df %>%
            dplyr::mutate(
              hover_exon = paste0(
                "Gene: ", dplyr::coalesce(gene_name, gene_id),
                "<br>Feature: exon",
                "<br>Range: ", chr_ucsc, ":", start, "-", end
              )
            )
        }
        
        if (!is.null(intr_df) && nrow(intr_df) > 0) {
          intr_df <- intr_df %>%
            dplyr::mutate(
              hover_intron = paste0(
                "Gene: ", dplyr::coalesce(gene_name, gene_id),
                "<br>Feature: intron",
                "<br>Range: ", chr_ucsc, ":", x0_bp, "-", x1_bp
              )
            )
        }
        
        hover_pts <- dplyr::bind_rows(
          gene_df %>% dplyr::transmute(x = (x0 + x1)/2, y = y, text = hover_gene),
          if (is.data.frame(ex_df) && nrow(ex_df) > 0) ex_df %>% dplyr::transmute(x = (x0 + x1)/2, y = y, text = hover_exon) else NULL,
          if (!is.null(intr_df) && nrow(intr_df) > 0) intr_df %>% dplyr::transmute(x = (x0 + x1)/2, y = y, text = hover_intron) else NULL
        ) %>%
          dplyr::filter(is.finite(x), is.finite(y), !is.na(text), nzchar(text))
        
        label_max <- 35
        gene_labels <- gene_df %>%
          dplyr::arrange(y) %>%
          dplyr::slice_head(n = label_max) %>%
          dplyr::transmute(
            x_lab = x_limits2[2],
            y     = y + 0.35,
            label = gene_name,
            hover = hover_gene
          )
        
        p_gene <- ggplot2::ggplot() +
          ggplot2::geom_segment(
            data = gene_df,
            ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
            linewidth = 0.8, alpha = 0.95
          ) +
          { if (!is.null(intr_df) && nrow(intr_df) > 0)
            ggplot2::geom_segment(
              data = intr_df,
              ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
              linewidth = 0.55, alpha = 0.90
            ) else NULL } +
          { if (is.data.frame(ex_df) && nrow(ex_df) > 0)
            ggplot2::geom_rect(
              data = ex_df,
              ggplot2::aes(xmin = x0, xmax = x1, ymin = y - 0.22, ymax = y + 0.22),
              fill  = "orange",
              color = "darkgrey",
              alpha = 0.95
            ) else NULL } +
          ggplot2::geom_point(
            data = hover_pts,
            ggplot2::aes(x = x, y = y, text = text),
            alpha = 0, size = 2, inherit.aes = FALSE
          ) +
          ggplot2::scale_x_continuous(limits = x_limits2, expand = c(0, 0)) +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::theme_void() +
          ggplot2::theme(plot.margin = ggplot2::margin(t = 6, r = 6, b = 0, l = 6))
        
        g_gene <- suppressWarnings(plotly::ggplotly(p_gene, tooltip = "text", source = "ld")) %>%
          plotly::layout(
            hovermode = "closest",
            showlegend = FALSE,
            margin = list(l = 80, r = 10, t = 8, b = 0),
            xaxis = list(showticklabels = FALSE, title = "", range = x_limits2, ticks = ""),
            yaxis = list(showticklabels = FALSE, title = "")
          )
        
        old_anns <- g_gene$x$layout$annotations %||% list()
        
        lab_anns <- lapply(seq_len(nrow(gene_labels)), function(i) {
          list(
            x = gene_labels$x_lab[i],
            y = gene_labels$y[i],
            xref = "x",
            yref = "y",
            text = gene_labels$label[i],
            showarrow = FALSE,
            xanchor = "right",
            yanchor = "middle",
            font = list(size = 11, color = "black")
          )
        })
        
        g_gene <- g_gene %>%
          plotly::add_markers(
            data = gene_labels,
            x = ~x_lab, y = ~y,
            text = ~hover,
            hoverinfo = "text",
            marker = list(opacity = 0, size = 12),
            inherit = FALSE,
            showlegend = FALSE
          ) %>%
          plotly::layout(annotations = c(old_anns, lab_anns))
        
        g_gene
      }
      
      ######## end gene model
      # --- g0 = Gene model (TOP), com LD_Inspector ---
      g0 <- NULL
      if (is.finite(chr_sel) && is.finite(st) && is.finite(en)) {
        g0 <- make_gene_track_plot(
          chr_sel   = chr_sel,
          st        = st,
          en        = en,
          bp_to_x   = bp_to_x,
          x_limits2 = x_limits2
        )
      }
      
      # combine
      if (!is.null(g0)) {
        plotly::subplot(g0, g1, g2, nrows = 3, heights = c(0.10, 0.35, 0.55), shareX = TRUE, titleX = TRUE) %>%
          plotly::layout(hovermode = "closest")
      } else {
        plotly::subplot(g1, g2, nrows = 2, heights = c(0.35, 0.65), shareX = TRUE, titleX = TRUE) %>%
          plotly::layout(hovermode = "closest")
      }
    })
  })
}
