# GItools_local_deeplinks_ALL_IN_ONE.R
# -------------------------------------------------------------------------
# Copy/paste file. Includes:
#   (A) start_gitools_local(): launcher to run all apps locally (fixed ports)
#   (B) gi_build_url(): helper to build deep-link URLs
#   (C) Deeplink handlers (drop-in blocks) for each app:
#       - Catalog Inspector
#       - GTEX Inspector
#       - NonSyn Inspector
#       - EWAStum Inspector (EWAS_cancer)
#       - EWASDis Inspector (EWAS_disease)
#       - LD Inspector
#
# IMPORTANT:
# - The “handlers” are written as functions you call INSIDE each app's server:
#     gitools_deeplink_<app>(session, <clusters reactive>, ...)
# - You still need to paste ONE LINE per app to call the right handler.
# - This is the cleanest “single file” approach: one shared file + one-line hook.
# -------------------------------------------------------------------------

#`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# -------------------------------------------------------------------------
# (A) Launcher: run all apps locally on fixed ports
# -------------------------------------------------------------------------
# --- portable app directories (based on GItools root) ---
gi_app_dirs <- function() {
  
  # Prefer gi_cfg() (portable)
  if (exists("gi_cfg", mode = "function")) {
    root <- gi_cfg()$root
    base <- file.path(root, "app")
    
  } else {
    
    # fallback: infer from this file in <root>/_shared/
    this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
    this_dir  <- if (nzchar(this_file)) dirname(this_file) else ""
    root <- if (nzchar(this_dir)) file.path(this_dir, "..") else getwd()
    
    # try to source <root>/config.R so gi_cfg() becomes available
    cfg_file <- normalizePath(file.path(root, "config.R"), winslash = "/", mustWork = FALSE)
    if (file.exists(cfg_file)) source(cfg_file, local = TRUE)
    
    # if gi_cfg now exists, use it; else keep the inferred root
    if (exists("gi_cfg", mode = "function")) {
      root <- gi_cfg()$root
    }
    
    base <- file.path(root, "app")
  }
  
  list(
    catalog = file.path(base, "Catalog_inspector"),
    gtex    = file.path(base, "GTEX_inspector"),
    nonsyn  = file.path(base, "NonSyn_Inspector"),
    ewastum = file.path(base, "EWAS_cancer"),
    ewasdis = file.path(base, "EWAS_disease"),
    ld      = file.path(base, "LD_Inspector")
  )
}

start_gitools_local <- function(
    open  = c("nonsyn","catalog","gtex","ewastum","ewasdis","ld"),
    ports = c(catalog=7001, gtex=7002, nonsyn=7003, ewastum=7004, ewasdis=7005, ld=7006),
    paths = NULL,
    start_all = FALSE
) {
  open <- match.arg(open)
  if (is.null(paths)) paths <- gi_app_dirs()
  
  run_one <- function(app_dir, port) {
    stopifnot(dir.exists(app_dir))
    cmd <- sprintf("shiny::runApp('%s', port=%d, host='127.0.0.1', launch.browser=FALSE)", app_dir, port)
    system2("R", c("-e", shQuote(cmd)), wait = FALSE)
    message(sprintf("[OK] %s -> http://127.0.0.1:%d/", basename(app_dir), port))
  }
  
  if (isTRUE(start_all)) {
    for (nm in names(paths)) run_one(paths[[nm]], ports[[nm]])
  } else {
    run_one(paths[[open]], ports[[open]])
  }
  
  utils::browseURL(sprintf("http://127.0.0.1:%d/", ports[[open]]))
  invisible(TRUE)
}


# -------------------------------------------------------------------------
# (B) Deep-link builder (use in Hub or anywhere)
# -------------------------------------------------------------------------
gi_build_url <- function(base, params) {
  qs <- vapply(names(params), function(k) {
    v <- params[[k]]
    if (is.null(v) || !nzchar(as.character(v))) return(NA_character_)
    paste0(k, "=", utils::URLencode(as.character(v), reserved = TRUE))
  }, FUN.VALUE = character(1))
  qs <- qs[!is.na(qs)]
  paste0(sub("/+$", "/", base), "?", paste(qs, collapse = "&"))
}

gi_local_urls <- function(ports = c(catalog=7001, gtex=7002, nonsyn=7003, ewastum=7004, ewasdis=7005, ld=7006)) {
  list(
    catalog = sprintf("http://127.0.0.1:%d/", ports[["catalog"]]),
    gtex    = sprintf("http://127.0.0.1:%d/", ports[["gtex"]]),
    nonsyn  = sprintf("http://127.0.0.1:%d/", ports[["nonsyn"]]),
    ewastum = sprintf("http://127.0.0.1:%d/", ports[["ewastum"]]),
    ewasdis = sprintf("http://127.0.0.1:%d/", ports[["ewasdis"]]),
    ld      = sprintf("http://127.0.0.1:%d/", ports[["ld"]])
  )
}

# -------------------------------------------------------------------------
# (C) Shared internal helpers (used by handlers)
# -------------------------------------------------------------------------
gi_read_qs <- function(session) {
  shiny::parseQueryString(session$clientData$url_search %||% "")
}

gi_select_dt_row <- function(dt_id, row_index) {
  if (!length(row_index) || is.na(row_index) || row_index < 1) return(invisible(FALSE))
  DT::selectRows(DT::dataTableProxy(dt_id), row_index)
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# (D) Deeplink handlers per app (call inside each app's server)
# -------------------------------------------------------------------------

# -------------------------
# Catalog Inspector handler
# -------------------------
# Expects:
# - clusters_val() -> data.frame with 'cluster_id' (preferred)
# - DT outputId: "cluster_dt"
# - Optional inputs: "catalog_pick_chr", "catalog_pick_cluster_id"
gitools_deeplink_catalog <- function(session, clusters_val) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), clusters_val()), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- clusters_val()
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    # Select DT row by cluster_id
    if (!is.null(qs$cluster) && "cluster_id" %in% names(cl)) {
      idx <- which(as.character(cl$cluster_id) == as.character(qs$cluster))
      if (length(idx) == 1) gi_select_dt_row("cluster_dt", idx[1])
    }
    
    # Optional: update pickers if they exist
    if (!is.null(qs$chr) && "catalog_pick_chr" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "catalog_pick_chr", selected = as.character(as.integer(qs$chr))))
    }
    if (!is.null(qs$cluster) && "catalog_pick_cluster_id" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "catalog_pick_cluster_id", selected = as.character(qs$cluster)))
    }
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# ----------------------
# GTEX Inspector handler
# ----------------------
# Expects:
# - clusters_cur() -> data.frame with cluster_chr_n and/or cluster
# - DT outputId: "cluster_dt"
# - Optional inputs: "func_chr", "func_cluster_id"
gitools_deeplink_gtex <- function(session, clusters_cur) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), clusters_cur()), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- clusters_cur()
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    if (!is.null(qs$cluster)) {
      idx <- integer(0)
      if ("cluster_chr_n" %in% names(cl)) idx <- which(as.character(cl$cluster_chr_n) == as.character(qs$cluster))
      if (!length(idx) && "cluster" %in% names(cl)) idx <- which(as.character(cl$cluster) == as.character(qs$cluster))
      if (length(idx) == 1) gi_select_dt_row("cluster_dt", idx[1])
    }
    
    if (!is.null(qs$chr) && "func_chr" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "func_chr", selected = as.character(as.integer(qs$chr))))
    }
    if (!is.null(qs$cluster) && "func_cluster_id" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "func_cluster_id", selected = as.character(qs$cluster)))
    }
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# -----------------------
# NonSyn Inspector handler
# -----------------------
# Expects:
# - cluster_dt_view2() -> data.frame with 'cluster_id'
# - DT outputId: "cluster_dt"
# - Optional LD inputs: "ld_chr", "ld_cluster_id"
gitools_deeplink_nonsyn <- function(session, cluster_dt_view2) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), cluster_dt_view2()), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- cluster_dt_view2()
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    if (!is.null(qs$cluster) && "cluster_id" %in% names(cl)) {
      idx <- which(as.character(cl$cluster_id) == as.character(qs$cluster))
      if (length(idx) == 1) gi_select_dt_row("cluster_dt", idx[1])
    }
    
    # Optional: sync internal LD module selectors
    if (!is.null(qs$chr) && "ld_chr" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "ld_chr", selected = as.character(as.integer(qs$chr))))
    }
    if (!is.null(qs$cluster) && "ld_cluster_id" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "ld_cluster_id", selected = as.character(qs$cluster)))
    }
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# -------------------------
# EWAStum (EWAS_cancer) handler
# -------------------------
# Expects:
# - rv$clusters data.frame with 'cluster_id'
# - DT outputId: "cluster_dt"
# - Optional: qs$cancer -> stored to rv$cancer_sel if present
gitools_deeplink_ewastum <- function(session, rv) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), rv$clusters), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- rv$clusters
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    if (!is.null(qs$cluster) && "cluster_id" %in% names(cl)) {
      idx <- which(as.character(cl$cluster_id) == as.character(qs$cluster))
      if (length(idx) == 1) gi_select_dt_row("cluster_dt", idx[1])
    }
    
    if (!is.null(qs$cancer)) rv$cancer_sel <- as.character(qs$cancer)
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# -------------------------
# EWASDis (EWAS_disease) handler
# -------------------------
# Expects:
# - rv$clusters data.frame with 'cluster_id'
# - DT outputId: "cluster_dt"
gitools_deeplink_ewasdis <- function(session, rv) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), rv$clusters), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- rv$clusters
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    if (!is.null(qs$cluster) && "cluster_id" %in% names(cl)) {
      idx <- which(as.character(cl$cluster_id) == as.character(qs$cluster))
      if (length(idx) == 1) gi_select_dt_row("cluster_dt", idx[1])
    }
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# -------------------
# LD Inspector handler
# -------------------
# Expects:
# - clusters_df() -> data.frame with cluster_key + cluster_id (+ optional start/end)
# - Inputs: ld_chr (SelectInput), ld_cluster_id (SelectInput; expects cluster_key)
gitools_deeplink_ld <- function(session, clusters_df) {
  qs_applied <- reactiveVal(FALSE)
  qs_store   <- reactiveVal(NULL)
  
  observe({
    qs <- gi_read_qs(session)
    if (length(qs)) qs_store(qs)
  })
  
  observeEvent(list(qs_store(), clusters_df()), {
    if (qs_applied()) return()
    qs <- qs_store(); cl <- clusters_df()
    req(length(qs) > 0, is.data.frame(cl), nrow(cl) > 0)
    
    if (!is.null(qs$chr) && "ld_chr" %in% names(session$input)) {
      suppressWarnings(updateSelectInput(session, "ld_chr", selected = as.character(as.integer(qs$chr))))
    }
    
    if (!is.null(qs$cluster) && "cluster_key" %in% names(cl)) {
      cid <- as.character(qs$cluster)
      
      idx <- integer(0)
      if ("cluster_id" %in% names(cl)) idx <- which(as.character(cl$cluster_id) == cid)
      if (!length(idx)) idx <- which(as.character(cl$cluster_key) == cid)
      
      if (!length(idx) && !is.null(qs$start) && !is.null(qs$end)) {
        st <- suppressWarnings(as.integer(qs$start))
        en <- suppressWarnings(as.integer(qs$end))
        if (all(c("start","end") %in% names(cl))) idx <- which(as.integer(cl$start) == st & as.integer(cl$end) == en)
      }
      
      if (length(idx) >= 1 && "ld_cluster_id" %in% names(session$input)) {
        key <- as.character(cl$cluster_key[idx[1]])
        suppressWarnings(updateSelectInput(session, "ld_cluster_id", selected = key))
      }
    }
    
    qs_applied(TRUE)
  }, ignoreInit = TRUE)
}

# -------------------------------------------------------------------------
# (E) ONE-LINE HOOKS (paste into each app's server)
# -------------------------------------------------------------------------
# Catalog Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_catalog(session, clusters_val)
#
# GTEX Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_gtex(session, clusters_cur)
#
# NonSyn Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_nonsyn(session, cluster_dt_view2)
#
# EWAStum Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_ewastum(session, rv)
#
# EWASDis Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_ewasdis(session, rv)
#
# LD Inspector (inside server):
#   source("GItools_local_deeplinks_ALL_IN_ONE.R", local=TRUE)
#   gitools_deeplink_ld(session, clusters_df)
#
# -------------------------------------------------------------------------
# (F) Quick test URLs (after running start_gitools_local())
# -------------------------------------------------------------------------
# Catalog: http://127.0.0.1:7001/?chr=5&cluster=chr5_1
# GTEx:    http://127.0.0.1:7002/?cluster=chr5_1
# NonSyn:  http://127.0.0.1:7003/?cluster=chr5_1
# LD:      http://127.0.0.1:7006/?chr=5&cluster=chr5_1
# -------------------------------------------------------------------------
