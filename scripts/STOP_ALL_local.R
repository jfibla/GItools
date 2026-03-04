#!/usr/bin/env Rscript

# STOP_ALL_local.R
# Stop GItools Hub + apps started locally (portable, GitHub-friendly)
#
# Usage:
#   Rscript --vanilla scripts/STOP_ALL_local.R
#   Rscript --vanilla scripts/STOP_ALL_local.R --kill-ngrok
#   Rscript --vanilla scripts/STOP_ALL_local.R --clean-logs
#   Rscript --vanilla scripts/STOP_ALL_local.R --ports=7101,7201,7202,7203,7204,7205,7206
#
# Tip:
#   export GITOOLS_REPO=/path/to/GItools   (helps locate LOG_DIR for --clean-logs)

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b
has_cmd <- function(cmd) nzchar(Sys.which(cmd))

args <- commandArgs(trailingOnly = TRUE)

default_ports <- c(7101L, 7201L, 7202L, 7203L, 7204L, 7205L, 7206L)

ports_arg <- grep("^--ports=", args, value = TRUE)
ports <- default_ports
if (length(ports_arg)) {
  s <- sub("^--ports=", "", ports_arg[1])
  s <- gsub("[[:space:]]+", "", s)
  v <- unlist(strsplit(s, ",", fixed = TRUE))
  v <- suppressWarnings(as.integer(v))
  v <- v[is.finite(v)]
  if (length(v)) ports <- v
}

kill_ngrok <- "--kill-ngrok" %in% args
clean_logs <- "--clean-logs" %in% args

# -----------------------------
# Repo/log dir best-effort
# -----------------------------
detect_repo_root <- function(seed_dir) {
  seed_dir <- normalizePath(seed_dir, winslash = "/", mustWork = FALSE)
  
  # env var first
  env <- Sys.getenv("GITOOLS_REPO", unset = "")
  cand <- character(0)
  if (nzchar(env)) cand <- c(cand, normalizePath(env, winslash="/", mustWork=FALSE))
  
  up <- seed_dir
  for (i in 0:6) {
    cand <- c(cand, up)
    parent <- dirname(up)
    if (identical(parent, up)) break
    up <- parent
  }
  
  cand <- unique(cand)
  ok <- vapply(cand, function(p) {
    file.exists(file.path(p, "config.R")) &&
      dir.exists(file.path(p, "app", "GItools_Hub"))
  }, logical(1))
  
  if (!any(ok)) return(NA_character_)
  cand[which(ok)[1]]
}

REPO_ROOT <- detect_repo_root(getwd())
LOG_DIR <- if (!is.na(REPO_ROOT)) file.path(REPO_ROOT, "app", "_logs") else NA_character_

# -----------------------------
# Helpers: kill by port
# -----------------------------
find_listen_pids <- function(port) {
  if (!has_cmd("lsof")) return(integer(0))
  port <- as.integer(port)
  cmd <- sprintf("lsof -nP -iTCP:%d -sTCP:LISTEN -t 2>/dev/null", port)
  out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error = function(e) character(0)))
  pids <- suppressWarnings(as.integer(out))
  pids <- pids[is.finite(pids)]
  unique(pids)
}

kill_pids <- function(pids, sig = "-TERM") {
  if (!length(pids)) return(invisible(FALSE))
  for (pid in pids) {
    suppressWarnings(system2("kill", c(sig, as.character(pid)), stdout = FALSE, stderr = FALSE))
  }
  invisible(TRUE)
}

# -----------------------------
# Stop sequence
# -----------------------------
cat("=== GItools STOP (local) ===\n")
cat("Ports:", paste(ports, collapse = ", "), "\n")
if (!is.na(REPO_ROOT)) cat("REPO_ROOT:", REPO_ROOT, "\n")
if (!is.na(LOG_DIR))   cat("LOG_DIR  :", LOG_DIR, "\n")
cat("\n")

if (!has_cmd("lsof")) {
  cat("[ERROR] lsof not found. Install it or stop processes manually.\n")
  quit(status = 1)
}

killed <- data.frame(port = integer(0), pid = integer(0), action = character(0), stringsAsFactors = FALSE)

# First try TERM
for (p in ports) {
  pids <- find_listen_pids(p)
  if (!length(pids)) {
    cat(sprintf("[OK]  port %-5d : free\n", p))
    next
  }
  cat(sprintf("[STOP] port %-5d : pids=%s (TERM)\n", p, paste(pids, collapse = ",")))
  kill_pids(pids, sig = "-TERM")
  killed <- rbind(killed, data.frame(port = p, pid = pids, action = "TERM", stringsAsFactors = FALSE))
}

# Wait, then KILL if still listening
Sys.sleep(0.6)
for (p in ports) {
  pids2 <- find_listen_pids(p)
  if (!length(pids2)) next
  cat(sprintf("[STOP] port %-5d : still listening pids=%s (KILL)\n", p, paste(pids2, collapse = ",")))
  kill_pids(pids2, sig = "-KILL")
  killed <- rbind(killed, data.frame(port = p, pid = pids2, action = "KILL", stringsAsFactors = FALSE))
}

# Optional: ngrok
if (isTRUE(kill_ngrok)) {
  if (has_cmd("pkill")) {
    cat("\n[ngrok] pkill -f 'ngrok http'\n")
    suppressWarnings(system2("pkill", c("-f", "ngrok http"), stdout = FALSE, stderr = FALSE))
  } else {
    cat("\n[ngrok] pkill not found -> skip\n")
  }
}

# Optional: clean logs
if (isTRUE(clean_logs)) {
  if (!is.na(LOG_DIR) && dir.exists(LOG_DIR)) {
    cat("\n[logs] cleaning:", LOG_DIR, "\n")
    pat <- "^(gitools_|runner_).*(\\.log|\\.R)$|^(ngrok\\.log|ngrok_urls\\.json|start_all_console\\.log)$"
    files <- list.files(LOG_DIR, full.names = TRUE)
    todel <- files[grepl(pat, basename(files))]
    if (length(todel)) {
      unlink(todel, recursive = FALSE, force = TRUE)
      cat("[logs] removed:", length(todel), "files\n")
    } else {
      cat("[logs] nothing to remove\n")
    }
  } else {
    cat("\n[logs] LOG_DIR not found -> skip\n")
  }
}

cat("\n=== SUMMARY ===\n")
if (nrow(killed)) print(killed, row.names = FALSE) else cat("No listening processes found on specified ports.\n")

cat("\n=== VERIFY ===\n")
for (p in ports) {
  pids <- find_listen_pids(p)
  if (!length(pids)) cat(sprintf("port %-5d : free\n", p))
  else               cat(sprintf("port %-5d : STILL LISTENING pids=%s\n", p, paste(pids, collapse = ",")))
}
cat("\nDone.\n")