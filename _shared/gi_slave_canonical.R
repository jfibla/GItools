# R/gi_slave_canonical.R
# ------------------------------------------------------------
# GItools CANONICAL SLAVE SYNC (common to ALL apps)
# - Reads SID from input$gi_qs or url_search
# - Polls state JSON
# - Loads standardized GWAS + CLUSTERS from RDS
# - Applies UI params (best effort)
# - Exposes gwas_shared_r(), clusters_shared_r()
# ------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

gi_sanitize_sid <- function(s) {
  s <- as.character(s %||% "")
  s <- gsub("[^0-9A-Za-z_\\-]", "", s)
  if (!nzchar(s)) NULL else s
}

# IMPORTANT: adapt these roots to your project
# Portable defaults: use options if set, else globals from gi_state.R, else cfg, else infer.
gi_get_shared_root <- function() {
  x <- getOption("gi_shared_root", "")
  if (nzchar(x)) return(x)
  if (exists("gi_shared_root")) return(gi_shared_root)
  if (exists("gi_cfg", mode = "function")) return(gi_cfg()$shared)
  this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "")
  if (nzchar(this_dir)) return(normalizePath(this_dir, winslash = "/", mustWork = FALSE))
  ""
}

gi_get_state_root <- function() {
  x <- getOption("gi_state_root", "")
  if (nzchar(x)) return(x)
  if (exists("gi_state_root")) return(gi_state_root)
  if (exists("gi_cfg", mode = "function")) return(file.path(gi_cfg()$root, "app", "_state"))
  this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "")
  if (nzchar(this_dir)) {
    root_guess <- normalizePath(file.path(this_dir, ".."), winslash = "/", mustWork = FALSE)
    return(file.path(root_guess, "app", "_state"))
  }
  ""
}

# and then wherever you used getOption(..., "/Volumes/...") replace with:
shared_root <- gi_get_shared_root()
state_root  <- gi_get_state_root()


gi_state_paths <- function(sid) {
 # root <- gi_state_root()
  root <- if (exists("gi_state_root") && is.function(gi_state_root)) gi_state_root() else gi_state_root
  
  list(
    json = file.path(root, paste0("state_", sid, ".json")),
    gwas = file.path(root, paste0("gwas_", sid, ".rds")),
    clus = file.path(root, paste0("clusters_", sid, ".rds"))
  )
}

gi_read_state <- function(sid) {
  p <- gi_state_paths(sid)$json
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::fromJSON(p), error = function(e) NULL)
}

gi_std_gwas <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  
  # tolerate various incoming formats
  nms <- names(df)
  
  # ensure required names exist
  if (!"CHR" %in% nms && "chr" %in% nms) df$CHR <- df$chr
  if (!"BP"  %in% nms && "bp"  %in% nms) df$BP  <- df$bp
  if (!"Pval" %in% nms && "pval" %in% nms) df$Pval <- df$pval
  if (!"logp" %in% nms && "LOGP" %in% nms) df$logp <- df$LOGP
  
  if (!"snp" %in% nms && "SNP" %in% nms) df$snp <- df$SNP
  if (!"rsid" %in% nms && "RSID" %in% nms) df$rsid <- df$RSID
  
  df$CHR  <- suppressWarnings(as.integer(df$CHR))
  df$BP   <- suppressWarnings(as.integer(df$BP))
  df$Pval <- suppressWarnings(as.numeric(df$Pval))
  if (!"logp" %in% names(df)) df$logp <- -log10(df$Pval)
  
  if (!"snp" %in% names(df))  df$snp  <- NA_character_
  if (!"rsid" %in% names(df)) df$rsid <- NA_character_
  
  df <- df[is.finite(df$CHR) & is.finite(df$BP) & is.finite(df$Pval) & df$Pval > 0, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  
  # fallback snp display
  df$snp <- ifelse(is.na(df$snp) | !nzchar(df$snp), paste0("chr", df$CHR, ":", df$BP), as.character(df$snp))
  df$rsid <- as.character(df$rsid)
  
  df[, intersect(c("CHR","BP","snp","rsid","Pval","logp","BPcum"), names(df)), drop = FALSE]
}

gi_std_clusters <- function(cl0) {
  if (!is.data.frame(cl0) || !nrow(cl0)) return(NULL)
  cl <- as.data.frame(cl0)
  
  if (!"chr" %in% names(cl) && "CHR" %in% names(cl)) cl$chr <- cl$CHR
  if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
  if (!"end" %in% names(cl) && "end_bp" %in% names(cl)) cl$end <- cl$end_bp
  
  cl$chr   <- suppressWarnings(as.integer(cl$chr))
  cl$start <- suppressWarnings(as.integer(cl$start))
  cl$end   <- suppressWarnings(as.integer(cl$end))
  cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
  if (!nrow(cl)) return(NULL)
  
  cl <- cl[order(cl$chr, cl$start, cl$end), , drop = FALSE]
  
  cl$cluster_chr <- ave(cl$chr, cl$chr, FUN = seq_along)
  cl$cluster_id  <- if ("cluster_id" %in% names(cl)) as.character(cl$cluster_id) else paste0("chr", cl$chr, "_", cl$cluster_chr)
  cl$cluster_chr_n <- paste0("chr", cl$chr, "_", cl$cluster_chr)
  
  cl$center <- as.integer(round((cl$start + cl$end) / 2))
  cl$cluster_size_kb <- round((cl$end - cl$start + 1) / 1000, 2)
  
  # keep only canonical cols + any extras
  keep <- unique(c(
    "chr","start","end",
    "cluster_id","cluster_chr","cluster_chr_n",
    "center","cluster_size_kb",
    "n_snps","top_snp","top_logp"
  ))
  cl[, intersect(keep, names(cl)), drop = FALSE]
}

gi_slave_canonical_init <- function(session) {
  
  # -------------------------------------------------------------------
  # Helpers (local)
  # -------------------------------------------------------------------
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b
  }
  
  safe_update <- function(expr) {
    tryCatch(force(expr), error = function(e) NULL)
    invisible(TRUE)
  }
  
  # -------------------------------------------------------------------
  # SID (shared session identifier)
  # -------------------------------------------------------------------
  sid_rv <- shiny::reactiveVal(NULL)
  
  # 1) Preferred: input$gi_qs (ngrok/caddy). IGNORE empty.
  shiny::observeEvent(session$input$gi_qs, {
    qs <- as.character(session$input$gi_qs %||% "")
    if (!nzchar(qs)) return()
    q  <- shiny::parseQueryString(sub("^\\?", "", qs))
    s0 <- gi_sanitize_sid(q$sid %||% "")
    if (!is.null(s0)) sid_rv(s0)
  }, ignoreInit = FALSE)
  
  # 2) Fallback: url_search (local). IGNORE empty.
  shiny::observe({
    if (!is.null(sid_rv())) return()
    qs <- as.character(session$clientData$url_search %||% "")
    if (!nzchar(qs)) return()
    q  <- shiny::parseQueryString(sub("^\\?", "", qs))
    s0 <- gi_sanitize_sid(q$sid %||% "")
    if (!is.null(s0)) sid_rv(s0)
  })
  
  # Debug SID when it becomes available
  shiny::observeEvent(sid_rv(), {
    s <- sid_rv()
    if (is.null(s)) {
      cat("[SLAVE] sid = <NULL> (waiting)\n")
      return()
    }
    p <- gi_state_paths(s)
    cat("[SLAVE] sid =", s, "\n")
    cat("[SLAVE] json =", p$json, "exists=", file.exists(p$json), "\n")
    cat("[SLAVE] gwas =", p$gwas, "exists=", file.exists(p$gwas), "\n")
    cat("[SLAVE] clus =", p$clus, "exists=", file.exists(p$clus), "\n")
  }, ignoreInit = FALSE)
  
  # -------------------------------------------------------------------
  # Detect app name + decide if we skip UI update sync
  # (GTEx/EWAS often hang if we push update*Input too early)
  # -------------------------------------------------------------------
  app_wd   <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  app_name <- basename(app_wd)
  
  # More robust than regex $ (handles trailing slash, etc.)
  skip_ui_updates <- app_name %in% c("GTEX_inspector", "EWAS_disease", "EWAS_cancer")
  
  cat("[SLAVE] app_name=", app_name,
      " app_wd=", app_wd,
      " skip_ui_updates=", skip_ui_updates, "\n")
  
  # -------------------------------------------------------------------
  # Poll state JSON (mtime)
  # -------------------------------------------------------------------
  gi_state <- shiny::reactivePoll(
    intervalMillis = 4000,
    session = session,
    checkFunc = function() {
      s <- sid_rv()
      if (is.null(s)) return(0)
      pj <- gi_state_paths(s)$json
      if (!file.exists(pj)) return(0)
      as.numeric(file.info(pj)$mtime)
    },
    valueFunc = function() {
      s <- sid_rv()
      if (is.null(s)) return(NULL)
      gi_read_state(s)
    }
  )
  
  # -------------------------------------------------------------------
  # Shared RVs exported to apps
  # -------------------------------------------------------------------
  gi_last_stamp   <- shiny::reactiveVal(NA_integer_)
  gwas_shared_rv  <- shiny::reactiveVal(NULL)
  clus_shared_rv  <- shiny::reactiveVal(NULL)
  
  # -------------------------------------------------------------------
  # Apply state when JSON changes
  # -------------------------------------------------------------------
  shiny::observeEvent(gi_state(), {
    
    s <- sid_rv()
    if (is.null(s)) return()
    
    st <- gi_state()
    if (is.null(st) || !is.list(st)) return()
    
    stamp <- suppressWarnings(as.integer(st$stamp %||% 0L))
    if (!is.na(gi_last_stamp()) && identical(stamp, gi_last_stamp())) return()
    gi_last_stamp(stamp)
    
    p <- gi_state_paths(s)
    cat("[SLAVE] sync stamp=", stamp, " sid=", s, "\n")
    
    # ---- 1) GWAS ----
    if (file.exists(p$gwas)) {
      df0 <- tryCatch(readRDS(p$gwas), error = function(e) {
        cat("[SLAVE] readRDS(gwas) ERROR:", conditionMessage(e), "\n")
        NULL
      })
      df <- tryCatch(gi_std_gwas(df0), error = function(e) df0)
      if (is.data.frame(df) && nrow(df) > 0) {
        gwas_shared_rv(df)
        cat("[SLAVE] GWAS synced rows=", nrow(df), "\n")
      } else {
        cat("[SLAVE] GWAS file read but empty/invalid\n")
      }
    } else {
      cat("[SLAVE] GWAS file not found:", p$gwas, "\n")
    }
    
    # ---- 2) UI params (best-effort) ----
    # IMPORTANT: optional skip for apps that hang on update*Input
    if (!isTRUE(skip_ui_updates)) {
      
      safe_update(if ("cluster_method" %in% names(session$input) && !is.null(st$cluster_method))
        shiny::updateRadioButtons(session, "cluster_method", selected = st$cluster_method))
      
      safe_update(if ("hits_mode" %in% names(session$input) && !is.null(st$hits_mode))
        shiny::updateRadioButtons(session, "hits_mode", selected = st$hits_mode))
      
      safe_update(if ("pthr" %in% names(session$input) &&
                      identical(st$thr_type, "logp") &&
                      is.finite(as.numeric(st$thr_value)))
        shiny::updateSliderInput(session, "pthr", value = as.numeric(st$thr_value)))
      
      safe_update(if ("min_logp" %in% names(session$input) &&
                      identical(st$thr_type, "min_logp") &&
                      is.finite(as.numeric(st$thr_value)))
        shiny::updateSliderInput(session, "min_logp", value = as.numeric(st$thr_value)))
      
      safe_update(if ("flank" %in% names(session$input) &&
                      is.finite(as.numeric(st$flank_bp %||% NA)))
        shiny::updateNumericInput(session, "flank", value = as.integer(st$flank_bp)))
      
      safe_update(if ("win_bp" %in% names(session$input) &&
                      is.finite(as.numeric(st$win_bp %||% NA)))
        shiny::updateNumericInput(session, "win_bp", value = as.integer(st$win_bp)))
      
      safe_update(if ("step_bp" %in% names(session$input) &&
                      is.finite(as.numeric(st$step_bp %||% NA)))
        shiny::updateNumericInput(session, "step_bp", value = as.integer(st$step_bp)))
      
      safe_update(if ("min_hits" %in% names(session$input) &&
                      is.finite(as.numeric(st$min_hits %||% NA)))
        shiny::updateNumericInput(session, "min_hits", value = as.integer(st$min_hits)))
    }
    
    # ---- 3) CLUSTERS ----
    if (file.exists(p$clus)) {
      cl0 <- tryCatch(readRDS(p$clus), error = function(e) {
        cat("[SLAVE] readRDS(clus) ERROR:", conditionMessage(e), "\n")
        NULL
      })
      cl <- tryCatch(gi_std_clusters(cl0), error = function(e) cl0)
      if (is.data.frame(cl) && nrow(cl) > 0) {
        clus_shared_rv(cl)
        cat("[SLAVE] CLUSTERS synced rows=", nrow(cl), "\n")
      } else {
        cat("[SLAVE] Clusters file read but empty/invalid\n")
      }
    } else {
      cat("[SLAVE] Clusters file not found:", p$clus, "\n")
    }
    
  }, ignoreInit = FALSE)
  
  # -------------------------------------------------------------------
  # Return "handle" used by apps
  # -------------------------------------------------------------------
  list(
    sid            = sid_rv,
    gi_state       = gi_state,
    gwas_shared    = gwas_shared_rv,
    clusters_shared = clus_shared_rv
  )
}

