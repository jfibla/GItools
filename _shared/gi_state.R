# GItools shared state per-session (sid)
# Portable layout (recommended):
#   - <GItools>/app/_state/state_<sid>.json
#   - <GItools>/app/_state/gwas_<sid>.rds
#   - <GItools>/app/_state/clusters_<sid>.rds
#
# Roots come from GItools/config.R via gi_cfg() (preferred),
# with fallbacks based on this file location.

# ---- GItools state root (ALWAYS defined) ----
gi_state_root <- function() {
  opt <- getOption("gi_state_root", NULL)
  if (is.character(opt) && length(opt) && nzchar(opt[1])) {
    return(normalizePath(opt[1], winslash="/", mustWork=FALSE))
  }
  env <- Sys.getenv("GITOOLS_STATE_ROOT", "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash="/", mustWork=FALSE))
  }
  normalizePath(file.path(dirname(getwd()), "_state"), winslash="/", mustWork=FALSE)
}

# simple %||% (if missing)
`%||%` <- get0("%||%", ifnotfound = function(a, b) {
  if (is.null(a) || length(a) == 0) b else a
})

# ---- Resolve GItools roots (portable) ----
gi_resolve_roots <- function() {
  
  # 1) Preferred: config.R already sourced by app, gi_cfg() exists
  if (exists("gi_cfg", mode = "function")) {
    cfg <- gi_cfg()
    root   <- cfg$root
    shared <- cfg$shared
    
    # State directory:
    # - Allow override via env var (useful in servers / Docker)
    # - Default: <root>/app/_state (matches your current repo layout)
    st <- Sys.getenv("GITOOLS_STATE", unset = file.path(root, "app", "_state"))
    st <- normalizePath(st, winslash = "/", mustWork = FALSE)
    
    return(list(root = root, shared = shared, state = st))
  }
  
  # 2) Fallback: infer root from this file location (<root>/_shared/gi_state.R)
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  this_dir  <- tryCatch(dirname(this_file), error = function(e) NULL)
  
  if (!is.null(this_dir) && nzchar(this_dir)) {
    
    # Assume: this_dir == <root>/_shared
    root_guess <- normalizePath(file.path(this_dir, ".."), winslash = "/", mustWork = FALSE)
    shared_guess <- normalizePath(this_dir, winslash = "/", mustWork = FALSE)
    
    # Try to source <root>/config.R if present
    cfg_file <- file.path(root_guess, "config.R")
    cfg_file <- normalizePath(cfg_file, winslash = "/", mustWork = FALSE)
    
    if (file.exists(cfg_file)) {
      source(cfg_file, local = TRUE)
      if (exists("gi_cfg", mode = "function")) {
        cfg <- gi_cfg()
        root   <- cfg$root
        shared <- cfg$shared
        st <- Sys.getenv("GITOOLS_STATE", unset = file.path(root, "app", "_state"))
        st <- normalizePath(st, winslash = "/", mustWork = FALSE)
        return(list(root = root, shared = shared, state = st))
      }
    }
    
    # If config.R not found or gi_cfg still missing, use guesses
    st_guess <- Sys.getenv("GITOOLS_STATE", unset = file.path(root_guess, "app", "_state"))
    st_guess <- normalizePath(st_guess, winslash = "/", mustWork = FALSE)
    
    return(list(root = root_guess, shared = shared_guess, state = st_guess))
  }
  
  # 3) Last resort: current working dir as root
  root_guess <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  shared_guess <- normalizePath(file.path(root_guess, "_shared"), winslash = "/", mustWork = FALSE)
  st_guess <- normalizePath(file.path(root_guess, "app", "_state"), winslash = "/", mustWork = FALSE)
  
  list(root = root_guess, shared = shared_guess, state = st_guess)
}

# ---- roots (globals expected by existing code) ----
gi_roots <- gi_resolve_roots()

gi_base_root   <- gi_roots$root
gi_shared_root <- gi_roots$shared
gi_state_root  <- gi_roots$state


gi_ensure_dirs <- function() {
  if (!dir.exists(gi_shared_root)) dir.create(gi_shared_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(gi_state_root))  dir.create(gi_state_root,  recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

gi_sanitize_sid <- function(x) {
  x <- as.character(x %||% "")
  x <- trimws(x)
  if (!nzchar(x)) return("default")
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  substr(x, 1, 64)
}

gi_sid <- function(session) {
  # Try multiple sources because behind reverse proxies url_search can be empty.
  qs <- NULL
  
  # 1) normal: url_search like "?sid=..."
  try(qs <- session$clientData$url_search, silent = TRUE)
  
  # 2) fallback: full url (session$clientData$url) -> extract query
  if (is.null(qs) || !nzchar(as.character(qs))) {
    u <- NULL
    try(u <- session$clientData$url, silent = TRUE)
    u <- as.character(u %||% "")
    if (nzchar(u) && grepl("\\?", u)) {
      qs <- sub("^[^?]*\\?", "?", u)
    }
  }
  
  # 3) last resort: nothing
  qs <- as.character(qs %||% "")
  qp <- list()
  if (nzchar(qs)) {
    qp <- tryCatch(shiny::parseQueryString(qs), error = function(e) list())
  }
  
  gi_sanitize_sid(qp[["sid"]] %||% "default")
}

# ---- GItools state root (ALWAYS defined) ----
gi_state_root <- function() {
  
  # 1) explicit option (preferred)
  opt <- getOption("gi_state_root", NULL)
  if (is.character(opt) && length(opt) && nzchar(opt[1])) {
    return(normalizePath(opt[1], winslash = "/", mustWork = FALSE))
  }
  
  # 2) env var (portable override)
  env <- Sys.getenv("GITOOLS_STATE_ROOT", "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  
  # 3) portable default: app/<ThisApp>  ->  app/_state
  normalizePath(file.path(dirname(getwd()), "_state"), winslash = "/", mustWork = FALSE)
}



#gi_state_paths <- function(sid) {
#  gi_ensure_dirs()
#  sid <- gi_sanitize_sid(sid)
#  list(
#    json = file.path(gi_state_root, paste0("state_", sid, ".json")),
#    gwas = file.path(gi_state_root, paste0("gwas_", sid, ".rds")),
#    clus = file.path(gi_state_root, paste0("clusters_", sid, ".rds"))
#  )
#}

gi_state_paths <- function(sid) {
  root <- gi_state_root()
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  
  list(
    root = root,
    json = file.path(root, paste0("state_", sid, ".json")),
    gwas = file.path(root, paste0("gwas_", sid, ".rds")),
    clus = file.path(root, paste0("clusters_", sid, ".rds"))
  )
}

gi_read_state <- function(sid) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  p <- gi_state_paths(sid)$json
  if (!file.exists(p)) return(NULL)
  txt <- tryCatch(readLines(p, warn = FALSE), error = function(e) NULL)
  if (is.null(txt) || !length(txt)) return(NULL)
  tryCatch(
    jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = TRUE),
    error = function(e) NULL
  )
}

gi_write_state <- function(sid, st) {
  stopifnot(requireNamespace("jsonlite", quietly = TRUE))
  p <- gi_state_paths(sid)$json
  jsonlite::write_json(st, p, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(p)
}

gi_bump_stamp <- function(prev) {
  if (is.null(prev) || !is.finite(as.numeric(prev))) return(1L)
  as.integer(prev) + 1L
}

# helper: build URLs with sid, preserving other params if you want
gi_url_with_sid <- function(base_url, sid) {
  sid <- gi_sanitize_sid(sid)
  if (grepl("\\?", base_url)) paste0(base_url, "&sid=", sid) else paste0(base_url, "?sid=", sid)
}

# OPTIONAL helpers (if you want to use them in master/slaves)
gi_write_gwas <- function(sid, gwas_df) {
  p <- gi_state_paths(sid)$gwas
  tryCatch(saveRDS(gwas_df, p), error = function(e) NULL)
  invisible(p)
}

gi_read_gwas <- function(sid) {
  p <- gi_state_paths(sid)$gwas
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
}

gi_write_clusters <- function(sid, clusters_df) {
  p <- gi_state_paths(sid)$clus
  tryCatch(saveRDS(clusters_df, p), error = function(e) NULL)
  invisible(p)
}

gi_read_clusters <- function(sid) {
  p <- gi_state_paths(sid)$clus
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
}
