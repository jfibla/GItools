#!/usr/bin/env Rscript

# start_ALL_local.R
#
# Start GItools Hub + all apps in background (portable, GitHub-friendly)
# - logs to <repo>/app/_logs
# - auto-detect repo root by finding config.R + app/GItools_Hub
# - optional: --ngrok to start tunnels and write ngrok_urls.json
#
# Usage (from repo root):
#   Rscript --vanilla scripts/start_ALL_local.R
#   Rscript --vanilla scripts/start_ALL_local.R --ngrok
#
# Tip: you can also run from anywhere if the script can find the repo root.
# Fallback:
#   export GITOOLS_REPO=/path/to/GItools

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b
has_cmd <- function(cmd) nzchar(Sys.which(cmd))

# -----------------------------
# Script path helpers (robust)
# -----------------------------
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[1]), winslash="/", mustWork = TRUE))
  NA_character_
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIR  <- if (is.na(SCRIPT_PATH)) normalizePath(getwd(), winslash="/", mustWork=TRUE) else dirname(SCRIPT_PATH)

# -----------------------------
# Auto-detect GItools repo root
# -----------------------------
detect_repo_root <- function(seed_dir) {
  seed_dir <- normalizePath(seed_dir, winslash = "/", mustWork = FALSE)
  
  # Walk up (seed -> parents) and test for repo signature
  up <- seed_dir
  cand <- character(0)
  for (i in 0:6) {
    if (nzchar(up)) cand <- c(cand, up)
    parent <- dirname(up)
    if (identical(parent, up)) break
    up <- parent
  }
  
  # Also accept explicit env var (highest priority if set)
  env <- Sys.getenv("GITOOLS_REPO", unset = "")
  if (nzchar(env)) cand <- c(normalizePath(env, winslash="/", mustWork=FALSE), cand)
  
  cand <- unique(cand)
  ok <- vapply(cand, function(p) {
    file.exists(file.path(p, "config.R")) &&
      dir.exists(file.path(p, "app", "GItools_Hub"))
  }, logical(1))
  
  if (!any(ok)) return(NA_character_)
  cand[which(ok)[1]]
}

REPO_ROOT <- detect_repo_root(SCRIPT_DIR)
if (is.na(REPO_ROOT)) {
  stop(paste0(
    "Could not detect GItools repo root.\n",
    "Expected: <repo>/config.R and <repo>/app/GItools_Hub\n",
    "Tip: set env var GITOOLS_REPO=/path/to/GItools\n",
    "Current seed: ", SCRIPT_DIR
  ))
}

APPS_DIR <- normalizePath(file.path(REPO_ROOT, "app"), winslash="/", mustWork=TRUE)
HUB_DIR  <- normalizePath(file.path(APPS_DIR, "GItools_Hub"), winslash="/", mustWork=TRUE)

message("[start] REPO_ROOT = ", REPO_ROOT)
message("[start] LOG_DIR   = ", file.path(REPO_ROOT, "app", "_logs"))
message("[start] APPS_DIR   = ", APPS_DIR)
message("[start] HUB_DIR    = ", HUB_DIR)

# -----------------------------
# Load portable config + shared
# -----------------------------
source(file.path(REPO_ROOT, "config.R"), local = TRUE)

# IMPORTANT: gi_cfg() / gi_find_root() uses getwd()
old_wd <- getwd()
setwd(REPO_ROOT)
on.exit(setwd(old_wd), add = TRUE)

cfg <- gi_cfg()

source(file.path(cfg$shared, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
source(file.path(cfg$shared, "gi_state.R"), local = TRUE)

# -----------------------------
# Logs (same as Hub)
# -----------------------------
LOG_DIR <- file.path(APPS_DIR, "_logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

HUB_LOG <- file.path(LOG_DIR, "hub.log")
hub_log <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", file = HUB_LOG, append = TRUE)
}

# -----------------------------
# Console tee -> file (portable)
# -----------------------------
START_LOG <- file.path(LOG_DIR, "start_all_console.log")
zz_out <- file(START_LOG, open = "at")
sink(zz_out, type="output", split=TRUE)

on.exit({
  try(sink(type="output"), silent = TRUE)
  try(close(zz_out), silent = TRUE)
}, add = TRUE)

logi <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n")
}
logi("[start] console tee -> ", START_LOG)

# -----------------------------
# Rscript binary
# -----------------------------
RSCRIPT <- file.path(R.home("bin"), "Rscript")
if (!file.exists(RSCRIPT)) RSCRIPT <- Sys.which("Rscript")
if (!nzchar(RSCRIPT) || !file.exists(RSCRIPT)) stop("Rscript not found.")

# -----------------------------
# Ports (match Hub)
# -----------------------------
HUB_PORT      <- 7101L
APP_PORT_BASE <- 7200L

apps <- list(
  catalog = list(dir=file.path(APPS_DIR, "Catalog_inspector"), port=APP_PORT_BASE + 1L),
  gtex    = list(dir=file.path(APPS_DIR, "GTEX_inspector"),    port=APP_PORT_BASE + 2L),
  nonsyn  = list(dir=file.path(APPS_DIR, "NonSyn_Inspector"),  port=APP_PORT_BASE + 3L),
  ewastum = list(dir=file.path(APPS_DIR, "EWAS_cancer"),       port=APP_PORT_BASE + 4L),
  ewasdis = list(dir=file.path(APPS_DIR, "EWAS_disease"),      port=APP_PORT_BASE + 5L),
  ld      = list(dir=file.path(APPS_DIR, "LD_Inspector"),      port=APP_PORT_BASE + 6L),
  hub     = list(dir=HUB_DIR,                                 port=HUB_PORT)
)

start_order <- c("catalog","gtex","nonsyn","ewasdis","ewastum","ld","hub")

# -----------------------------
# Port / HTTP helpers
# -----------------------------
find_listen_pid <- function(port) {
  if (!has_cmd("lsof")) return(NA_integer_)
  port <- as.integer(port)
  cmd <- sprintf("lsof -nP -iTCP:%d -sTCP:LISTEN -t 2>/dev/null | head -n 1", port)
  out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error = function(e) character(0)))
  pid <- suppressWarnings(as.integer(out[1]))
  if (!is.finite(pid)) NA_integer_ else pid
}
is_port_listening <- function(port) is.finite(find_listen_pid(port))

http_status_fast <- function(url) {
  if (!has_cmd("curl")) return(NA_integer_)
  args <- c("-s","-o","/dev/null","-w","%{http_code}","--connect-timeout","0.25","--max-time","0.8", url)
  out <- suppressWarnings(tryCatch(system2("curl", args, stdout = TRUE, stderr = TRUE),
                                   error = function(e) ""))
  code <- suppressWarnings(as.integer(out[1] %||% ""))
  if (!is.finite(code)) NA_integer_ else code
}

wait_for_listen <- function(port, timeout_sec = 40) {
  t0 <- Sys.time()
  repeat {
    pid <- find_listen_pid(port)
    if (is.finite(pid)) return(pid)
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) >= timeout_sec) break
    Sys.sleep(0.2)
  }
  NA_integer_
}

wait_for_http <- function(url, timeout_sec = 60) {
  t0 <- Sys.time()
  repeat {
    code <- http_status_fast(url)
    if (is.finite(code) && code >= 200 && code < 500) return(code)
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) >= timeout_sec) break
    Sys.sleep(0.5)
  }
  NA_integer_
}

# -----------------------------
# Runner-based background start
# -----------------------------
start_app_bg <- function(key, app_dir, port, log_dir = LOG_DIR, extra_env = list()) {
  stopifnot(dir.exists(app_dir))
  port <- as.integer(port)
  
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_log <- file.path(log_dir, sprintf("gitools_%d.out.log", port))
  err_log <- file.path(log_dir, sprintf("gitools_%d.err.log", port))
  runner  <- file.path(log_dir, sprintf("runner_%s_%d.R", key, port))
  
  if (!file.exists(out_log)) file.create(out_log)
  if (!file.exists(err_log)) file.create(err_log)
  
  app_dir_q <- paste0('"', gsub("\\\\", "/", app_dir), '"')
  log_dir_s <- gsub("\\\\", "/", log_dir)
  
  env_lines <- c(
    sprintf('Sys.setenv(GITOOLS_LOG_DIR = "%s")', log_dir_s),
    sprintf('Sys.setenv(SHINY_PORT = "%d")', port),
    'Sys.setenv(SHINY_HOST = "127.0.0.1")',
    sprintf('Sys.setenv(GITOOLS_REPO = "%s")', gsub('"', '\\"', normalizePath(REPO_ROOT, winslash="/", mustWork=TRUE)))
  )
  
  if (length(extra_env)) {
    env_lines <- c(env_lines, vapply(names(extra_env), function(k) {
      v <- extra_env[[k]]
      sprintf('Sys.setenv(%s = "%s")', k, gsub('"', '\\"', as.character(v)))
    }, ""))
  }
  
  runner_txt <- c(
    "options(shiny.port = NULL)",
    "options(shiny.host = NULL)",
    env_lines,
    sprintf("setwd(%s)", app_dir_q),
    sprintf('cat("[runner] starting %s on port %d\\n")', key, port),
    sprintf('shiny::runApp(%s, port=%d, host="127.0.0.1", launch.browser=FALSE)', app_dir_q, port)
  )
  writeLines(runner_txt, runner)
  
  system2(RSCRIPT, args = c("--vanilla", runner), wait = FALSE, stdout = out_log, stderr = err_log)
  
  pid <- NA_integer_
  for (i in 1:250) {
    pid <- find_listen_pid(port)
    if (is.finite(pid)) break
    Sys.sleep(0.1)
  }
  
  list(pid = pid, out_log = out_log, err_log = err_log, runner = runner)
}

# -----------------------------
# Optional NGROK
# -----------------------------
use_ngrok <- "--ngrok" %in% commandArgs(trailingOnly = TRUE)

start_ngrok_tunnels <- function(ports, log_dir = LOG_DIR) {
  if (!has_cmd("ngrok")) {
    message("[ngrok] ngrok not found in PATH -> skipping")
    return(list(ok = FALSE, urls = list()))
  }
  
  ngrok_log <- file.path(log_dir, "ngrok.log")
  if (!file.exists(ngrok_log)) file.create(ngrok_log)
  
  if (has_cmd("pkill")) suppressWarnings(system2("pkill", c("-f", "ngrok http"), stdout=NULL, stderr=NULL))
  
  urls <- list()
  for (p in ports) {
    p <- as.integer(p)
    system2("ngrok",
            args = c("http", as.character(p), "--log", "stdout"),
            wait = FALSE,
            stdout = ngrok_log,
            stderr = ngrok_log
    )
    Sys.sleep(0.25)
  }
  
  if (!has_cmd("curl")) return(list(ok = TRUE, urls = list()))
  api <- "http://127.0.0.1:4040/api/tunnels"
  j <- suppressWarnings(tryCatch(system2("curl", c("-s", api), stdout=TRUE, stderr=TRUE), error=function(e) ""))
  txt <- paste(j, collapse = "")
  if (!nzchar(txt)) return(list(ok = TRUE, urls = list()))
  
  addr  <- regmatches(txt, gregexpr('"addr":"[^"]+"', txt))[[1]]
  pub   <- regmatches(txt, gregexpr('"public_url":"[^"]+"', txt))[[1]]
  addr  <- sub('^"addr":"', "", sub('"$', "", addr))
  pub   <- sub('^"public_url":"', "", sub('"$', "", pub))
  
  for (i in seq_along(addr)) {
    a <- addr[i]
    prt <- suppressWarnings(as.integer(sub(".*:(\\d+)$", "\\1", a)))
    if (is.finite(prt) && i <= length(pub)) urls[[as.character(prt)]] <- pub[i]
  }
  
  list(ok = TRUE, urls = urls)
}

write_ngrok_urls_file <- function(urls_by_port, file) {
  ports <- names(urls_by_port)
  kv <- vapply(ports, function(p) {
    sprintf('  "%s": "%s"', p, gsub('"', '\\"', urls_by_port[[p]]))
  }, "")
  txt <- paste0("{\n", paste(kv, collapse = ",\n"), "\n}\n")
  writeLines(txt, file)
  invisible(TRUE)
}

# -----------------------------
# Start sequence
# -----------------------------
message("=== GItools start (portable) ===")
message("REPO_ROOT: ", REPO_ROOT)
message("LOG_DIR  : ", LOG_DIR)
message("Ports    : hub=", HUB_PORT, " | apps=", APP_PORT_BASE+1L, "..", APP_PORT_BASE+6L)
message("")

ngrok_urls_file <- file.path(LOG_DIR, "ngrok_urls.json")
extra_env <- list()

if (isTRUE(use_ngrok)) {
  message("[ngrok] starting tunnels...")
  ports <- vapply(apps, `[[`, 0L, "port")
  ng <- start_ngrok_tunnels(ports, LOG_DIR)
  if (isTRUE(ng$ok) && length(ng$urls)) {
    write_ngrok_urls_file(ng$urls, ngrok_urls_file)
    message("[ngrok] wrote: ", ngrok_urls_file)
    extra_env <- list(
      GITOOLS_URL_MODE = "ngrok",
      GITOOLS_NGROK_URLS_FILE = ngrok_urls_file
    )
    message("[ngrok] public URLs (by port):")
    for (p in names(ng$urls)) message("  ", p, " -> ", ng$urls[[p]])
  } else {
    message("[ngrok] tunnels started but could not read public URLs from 4040 API (yet).")
    message("[ngrok] check ", file.path(LOG_DIR, "ngrok.log"))
  }
  message("")
}

res <- list()

for (nm in start_order) {
  info <- apps[[nm]]
  stopifnot(!is.null(info$dir), !is.null(info$port))
  
  url <- sprintf("http://127.0.0.1:%d/", as.integer(info$port))
  
  if (!dir.exists(info$dir)) {
    message(sprintf("[SKIP] %-8s DIR NOT FOUND: %s", nm, info$dir))
    res[[nm]] <- list(app = nm, port = info$port, pid = NA_integer_, http = NA_integer_,
                      url = url, status = "DIR_MISSING",
                      err_log = NA_character_)
    next
  }
  
  if (is_port_listening(info$port)) {
    pid0  <- find_listen_pid(info$port)
    code0 <- http_status_fast(url)
    message(sprintf("[OK]   %-8s already LISTEN pid=%s (http=%s)  %s",
                    nm, pid0, ifelse(is.finite(code0), code0, "NA"), url))
    res[[nm]] <- list(app = nm, port = info$port, pid = pid0, http = code0,
                      url = url, status = "ALREADY_RUNNING",
                      err_log = file.path(LOG_DIR, sprintf("gitools_%d.err.log", info$port)))
    next
  }
  
  logs <- start_app_bg(nm, info$dir, info$port, LOG_DIR, extra_env = extra_env)
  pid <- wait_for_listen(info$port, timeout_sec = 40)
  
  if (!is.finite(pid)) {
    message(sprintf("[WARN] %-8s did not LISTEN within 40s (port %d). Check: %s",
                    nm, info$port, logs$err_log))
    res[[nm]] <- list(app = nm, port = info$port, pid = NA_integer_, http = NA_integer_,
                      url = url, status = "NO_LISTEN",
                      err_log = logs$err_log)
    next
  }
  
  code <- wait_for_http(url, timeout_sec = if (nm == "hub") 30 else 60)
  status <- if (is.finite(code) && code >= 200 && code < 500) "UP" else "LISTEN_ONLY"
  
  message(sprintf("[OK]   %-8s pid=%s  http=%s  status=%s  %s",
                  nm, pid, ifelse(is.finite(code), code, "NA"), status, url))
  
  res[[nm]] <- list(app = nm, port = info$port, pid = pid, http = code,
                    url = url, status = status,
                    err_log = logs$err_log)
}

to_df <- function(res) {
  do.call(rbind, lapply(res, function(x) {
    data.frame(
      app = x$app,
      port = x$port,
      pid = ifelse(is.finite(x$pid), as.character(x$pid), ""),
      http_code = ifelse(is.finite(x$http), x$http, 0L),
      status = x$status,
      url = x$url,
      err_log = x$err_log %||% "",
      stringsAsFactors = FALSE
    )
  }))
}

df <- to_df(res)
message("\n=== SUMMARY ===")
print(df, row.names = FALSE)

hub_url <- sprintf("http://127.0.0.1:%d/", apps$hub$port)
cat_url <- sprintf("http://127.0.0.1:%d/", apps$catalog$port)

hub_pid  <- find_listen_pid(apps$hub$port)

if (is.finite(hub_pid)) {
  message("\nOpening Hub (LISTEN pid=", hub_pid, "): ", hub_url)
  if (interactive()) utils::browseURL(hub_url)
} else {
  message("\nHub is DOWN (no listen on ", apps$hub$port, "). Opening Catalog instead: ", cat_url)
  if (interactive()) utils::browseURL(cat_url)
}

if (isTRUE(use_ngrok)) {
  message("\n[ngrok] If your Hub supports ngrok URLs, it can read: ", ngrok_urls_file)
}