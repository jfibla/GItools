# /Volumes/DISK1TB/Inspector_app_slaves_ngroc/GItools/app/GItools_Hub/app.R
#
# to kill all processes and remove logs:
#  pkill -f "shiny::runApp"
#  rm -rf /Volumes/DISK1TB/Inspector_app_slaves/_logs
#  mkdir -p /Volumes/DISK1TB/Inspector_app_slaves/_logs
#
# to run hub:
# shiny::runApp("/Volumes/DISK1TB/Inspector_app_slaves_github/GItools/app/GItools_Hub",
#              port=7101, host="127.0.0.1", launch.browser=TRUE)

# sequencia d'arrancada:'
# pkill -f "shiny::runApp"; pkill caddy; pkill ngrok
# Rscript --vanilla /Volumes/DISK1TB/Inspector_app_slaves_ngroc/start_ALL_local.R
# caddy run --config /Volumes/DISK1TB/Inspector_app_slaves_ngroc/Caddyfile
# ngrok http 8080

library(shiny)
library(DT)

RSCRIPT <- file.path(R.home("bin"), "Rscript")
if (!file.exists(RSCRIPT)) {
  # fallback per alguns entorns
  RSCRIPT <- Sys.which("Rscript")
}

# --- Portable config (repo-root + _shared) ---
source(file.path("..", "..", "config.R"), local = TRUE)  # des de GItools/app/GItools_Hub
cfg <- gi_cfg()

# Shared
source(file.path(cfg$shared, "GItools_local_deeplinks_ALL_IN_ONE.R"), local = TRUE)
source(file.path(cfg$shared, "gi_state.R"), local = TRUE)

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b

# Repo layout REAL: app/<AppName>
HUB_DIR  <- normalizePath(getwd(), winslash="/", mustWork=TRUE)         # .../GItools/app/GItools_Hub
APPS_DIR <- normalizePath(dirname(HUB_DIR), winslash="/", mustWork=TRUE) # .../GItools/app

message("[hub] HUB_DIR  = ", HUB_DIR)
message("[hub] APPS_DIR = ", APPS_DIR)
message("[hub] APPS_DIR contents: ", paste(list.files(APPS_DIR), collapse=", "))

#############################

gi_public_base_url <- function() {
  x <- Sys.getenv("GI_PUBLIC_BASE_URL", unset = "")
  x <- trimws(x)
  x <- sub("/+$", "", x)
  x
}

gi_app_url <- function(app_slug, local_port) {
  base <- gi_public_base_url()
  
  if (nzchar(base)) {
    return(paste0(base, "/", app_slug))
  }
  
  paste0("http://127.0.0.1:", local_port)
}

#############################

# ---- Ports (portable) ----
APP_PORT_BASE <- 7200 

apps <- list(
  hub        = list(name="🕸️ GItools Hub",            dir=HUB_DIR,                                  port=7101),
  catalog    = list(name="📚 Catalog Inspector",      dir=file.path(APPS_DIR, "Catalog_inspector"), port=APP_PORT_BASE + 1),
  gtex       = list(name="🧠 GTEx Inspector",         dir=file.path(APPS_DIR, "GTEX_inspector"),    port=APP_PORT_BASE + 2),
  nonsyn     = list(name="🧬 NonSyn Inspector",       dir=file.path(APPS_DIR, "NonSyn_Inspector"),  port=APP_PORT_BASE + 3),
  ewasdis    = list(name="🧫 EWASdis Inspector",      dir=file.path(APPS_DIR, "EWAS_disease"),      port=APP_PORT_BASE + 5),
  ewastum    = list(name="🧪 EWAStum Inspector",      dir=file.path(APPS_DIR, "EWAS_cancer"),       port=APP_PORT_BASE + 4),
  integrator = list(name="🧩 Integrator Inspector",   dir=file.path(APPS_DIR, "Integrator_Inspector"), port=APP_PORT_BASE + 6)
)


# Logs (en temp por defecto, siempre escribible)
LOG_DIR <- file.path(APPS_DIR, "_logs")   # -> GItools/app/_logs
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

HUB_LOG <- file.path(LOG_DIR, "hub_actions.log")
if (!file.exists(HUB_LOG)) file.create(HUB_LOG)

hub_log <- function(...) {
  f <- file.path(LOG_DIR, "hub.log")
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", file = f, append = TRUE)
}

options(shiny.error = function() {
  hub_log("[hub][ERROR] ", format(Sys.time()))
  hub_log(paste(capture.output(sys.calls()), collapse = "\n"))
})

hub_log("[hub] boot ", format(Sys.time()))
hub_log("[hub] LOG_DIR = ", LOG_DIR)
hub_log("[hub] APPS_DIR = ", APPS_DIR)
hub_log("[hub] RSCRIPT = ", RSCRIPT, " exists=", file.exists(RSCRIPT))

#if (!isTRUE(ok_logdir)) {
#  LOG_DIR <- file.path(tempdir(), "gitools_logs")
#  dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
#}

PID_FILE <- file.path(LOG_DIR, "_pids_gitools.rds")
# RSCRIPT  <- file.path(R.home("bin"), "Rscript")

busy_css <- "
.gitools-spinner {
  display:inline-block; width:14px; height:14px; margin-right:8px;
  border:2px solid #999; border-top-color: transparent; border-radius:50%;
  animation: spin 0.8s linear infinite;
}
@keyframes spin { to { transform: rotate(360deg); } }
.gitools-banner {
  padding:10px 12px; border-radius:8px; margin:10px 0;
  background:#f2f6ff; border:1px solid #d6e3ff; color:#163b72;
  font-weight:600;
}
"

btn_style <- "width:100%; margin-bottom:8px; font-size:16px; padding:12px 14px;"

has_cmd <- function(cmd) nzchar(Sys.which(cmd))

# --- single-port lsof (quiet) ---
find_listen_pid <- function(port) {
  if (!has_cmd("lsof")) return(NA_integer_)
  port <- as.integer(port)
  
  cmd <- sprintf("lsof -nP -iTCP:%d -sTCP:LISTEN -t 2>/dev/null | head -n 1", port)
  out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error = function(e) character(0)))
  
  pid <- suppressWarnings(as.integer(out[1]))
  if (!is.finite(pid)) NA_integer_ else pid
}

is_port_listening <- function(port) is.finite(find_listen_pid(port))

# --- fast local curl status (quiet + short timeouts) ---

http_status_fast <- function(url) {
  if (!has_cmd("curl")) return(NA_integer_)
  
  args <- c(
    "-s", "-o", "/dev/null",
    "-w", "%{http_code}",
    "--connect-timeout", "1.5",
    "--max-time", "3",
    url
  )
  out <- suppressWarnings(tryCatch(system2("curl", args, stdout = TRUE, stderr = TRUE),
                                   error = function(e) ""))
  code <- suppressWarnings(as.integer(out[1] %||% ""))
  if (!is.finite(code)) NA_integer_ else code
}


read_pids <- function() {
  if (file.exists(PID_FILE)) {
    x <- tryCatch(readRDS(PID_FILE), error=function(e) list())
    if (is.list(x)) return(x)
  }
  list()
}

write_pids <- function(x) {
  tryCatch(saveRDS(x, PID_FILE), error=function(e) NULL)
  invisible(TRUE)
}

start_app_bg <- function(key, app_dir, port, log_dir = NULL) {
  
  cat(sprintf("[hub] start_app_bg ENTER key=%s port=%d\n", key, port),
      file = HUB_LOG, append = TRUE)
  
  stopifnot(dir.exists(app_dir))
  port <- as.integer(port)
  
  # ✅ agafa LOG_DIR “live” (no enganxat al default)
  if (is.null(log_dir) || !nzchar(log_dir)) log_dir <- LOG_DIR
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_log <- file.path(log_dir, sprintf("gitools_%d.out.log", port))
  err_log <- file.path(log_dir, sprintf("gitools_%d.err.log", port))
  runner  <- file.path(log_dir, sprintf("runner_%s_%d.R", key, port))
  
  # ✅ crea fitxers (encara que estiguin buits)
  if (!file.exists(out_log)) file.create(out_log)
  if (!file.exists(err_log)) file.create(err_log)
  
  app_dir_q <- paste0('"', gsub("\\\\", "/", app_dir), '"')
  
  runner_txt <- c(
    "Sys.unsetenv('SHINY_PORT')",
    "Sys.unsetenv('SHINY_HOST')",
    "Sys.unsetenv('PORT')",
    "options(shiny.port = NULL)",
    "options(shiny.host = NULL)",
    sprintf('Sys.setenv(GITOOLS_LOG_DIR = "%s")', gsub("\\\\", "/", log_dir)),
    sprintf('Sys.setenv(TMPDIR = "%s")', gsub("\\\\","/", file.path(log_dir, "tmp"))),
    'dir.create(Sys.getenv("TMPDIR"), showWarnings=FALSE, recursive=TRUE)',
    sprintf("setwd(%s)", app_dir_q),
    sprintf('shiny::runApp(%s, port=%d, host="127.0.0.1", launch.browser=FALSE)', app_dir_q, port)
  )
  writeLines(runner_txt, runner)
  
  cat(sprintf("[hub] wrote runner %s\n", runner), file = HUB_LOG, append = TRUE)
  
  if (!file.exists(RSCRIPT)) stop(sprintf("Rscript not found at: %s", RSCRIPT))
  
  system2(
    RSCRIPT,
    args   = c("--vanilla", runner),
    wait   = FALSE,
    stdout = out_log,
    stderr = err_log
  )
  
  pid <- NA_integer_
  for (i in 1:250) { # 25s max
    pid <- find_listen_pid(port)
    if (is.finite(pid)) break
    Sys.sleep(0.1)
  }
  
  list(pid = pid, out_log = out_log, err_log = err_log, runner = runner)
}



kill_by_port <- function(port) {
  port <- as.integer(port)
  
  pids <- suppressWarnings(tryCatch(
    system2(
      "lsof",
      c("-tiTCP", paste0(":", port), "-sTCP:LISTEN"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0)
  ))
  
  pids <- unique(trimws(pids))
  pids <- pids[nzchar(pids)]
  
  if (!length(pids)) return(FALSE)
  
  for (pid in pids) {
    suppressWarnings(system2("kill", c("-9", pid), stdout = NULL, stderr = NULL))
  }
  
  Sys.sleep(0.5)
  !is_port_listening(port)
}

open_many_js <- function(url_vec) {
  js <- paste(vapply(url_vec, function(u) sprintf("window.open('%s','_blank');", u), ""),
              collapse = "")
  sprintf("<script>%s</script>", js)
}



ui <- fluidPage(
  title =   HTML("🕸️ GItools Hub"),
  tags$head(
    tags$style(HTML(busy_css)),
    tags$script(HTML("
      (function(){
        function setBanner(show){
          var el = document.getElementById('hubBanner');
          if(!el) return;
          el.style.display = show ? 'block' : 'none';
        }

        // Visible from first paint
        document.addEventListener('DOMContentLoaded', function(){
          setBanner(true);
        });

        // Server can hide/show
        if(window.Shiny){
          Shiny.addCustomMessageHandler('hub_banner', function(msg){
            setBanner(!!(msg && msg.show));
          });
        }
      })();
    "))
  ),
  
  # Banner (ONLY one)
  tags$div(
    id = "hubBanner",
    class = "gitools-banner",
    tags$span(class="gitools-spinner"),
    "Loading GItools Hub… Please wait while links are prepared."
  ),
  
  tags$h1("Genomic Inspector tools (GItools Hub - Web Launcher)"),
  tags$hr(),
  
  fluidRow(
    column(
      4,
    #  tags$h4("Open apps"),
      actionButton("info_00", "ℹ️ GItools Hub help"),
     tags$hr(),
      uiOutput("open_buttons_ui"),
    tags$hr(),
    downloadButton("dl_example_zip", "⬇️ Download example files (zip)", style="width:50%; font-size:14px; padding:12px 14px;"),

      tags$hr(),
      tags$div(style="margin:6px 0; font-size:14px; color:#555;",
               "Session SID: ", tags$code(textOutput("sid_txt", inline = TRUE)))
    ),
    
    column(
      8,
      tags$h4("Status"),
      
      fluidRow(
        column(4, actionButton("refresh_status", "Refresh status", class="btn btn-info", style="width:100%;")),
        column(4, actionButton("start_all", "Start all apps", class="btn btn-success", style="width:100%;")),
        column(4, actionButton("kill_all", "Kill all apps", class="btn btn-danger", style="width:100%;"))
      ),
      
      tags$div(style="margin-top:10px;"),
      DTOutput("status_table"),
      
      tags$hr(),
      tags$h4("Kill / Restart (one app)"),
      fluidRow(
        column(4, selectInput("manage_app", "App",
                              choices  = setdiff(names(apps), "hub"),
                              selected = "catalog")),
        column(
          8,
          actionButton("start_one", "Start", class="btn btn-success"),
          actionButton("kill_one", "Kill", class="btn btn-danger"),
          actionButton("restart_one", "Restart", class="btn btn-warning")
        )
      ),
      verbatimTextOutput("manage_out"),
      
      tags$hr(),
      tags$h4("Diagnostics"),
      verbatimTextOutput("diag_out"),
      tags$small(paste0("LOG_DIR: ", LOG_DIR))
    )
  )
)

server <- function(input, output, session) {
  
  # ---- Banner control: wait for BOTH open_buttons_ui + first status ----
  banner_hidden  <- reactiveVal(FALSE)
  buttons_ready  <- reactiveVal(FALSE)
  status_ready   <- reactiveVal(FALSE)
  banner_started <- Sys.time()
  
  BANNER_MIN_SECS <- 2.5   # <- ajusta: 2.0 / 3.0 segons el que vulguis
  
  observe({
    req(isTRUE(buttons_ready()))
    req(isTRUE(status_ready()))
    
    elapsed <- as.numeric(difftime(Sys.time(), banner_started, units = "secs"))
    wait_left <- max(0, BANNER_MIN_SECS - elapsed)
    
    if (wait_left > 0) {
      invalidateLater(ceiling(wait_left * 1000), session)
      return()
    }
    
    if (!isTRUE(banner_hidden())) {
      session$sendCustomMessage("hub_banner", list(show = FALSE))
      banner_hidden(TRUE)
    }
  })
  
  # --- SID builder ---
  sid <- reactiveVal(NULL)
  observeEvent(TRUE, {
    q <- shiny::parseQueryString(session$clientData$url_search %||% "")
    sid0 <- as.character(q$sid %||% "")
    if (!nzchar(sid0)) sid0 <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample(1000:9999, 1))
    sid(sid0)
  }, once = TRUE)
  

  # ------------------------------------------------------------
  # NGROK-aware URL builder (falls back to localhost)
  # Expects a JSON file written by start_ALL_local.R --ngrok:
  #   { "7101":"https://....", "7201":"https://....", ... }
  # Env vars:
  #   GITOOLS_URL_MODE="ngrok" (optional)
  #   GITOOLS_NGROK_URLS_FILE="/path/to/ngrok_urls.json"
  # ------------------------------------------------------------
  
  read_ngrok_urls <- function() {
    f <- Sys.getenv("GITOOLS_NGROK_URLS_FILE", unset = "")
    if (!nzchar(f) || !file.exists(f)) return(list())
    txt <- tryCatch(paste(readLines(f, warn = FALSE), collapse = "\n"), error = function(e) "")
    if (!nzchar(txt)) return(list())
    
    # Minimal JSON parser (no jsonlite dependency)
    # Extract "port":"url" pairs
    m <- gregexpr('"[0-9]+"\\s*:\\s*"[^"]+"', txt, perl = TRUE)
    hits <- regmatches(txt, m)[[1]]
    if (!length(hits)) return(list())
    
    out <- list()
    for (h in hits) {
      # h like: "7201": "https://xxxx.ngrok-free.app"
      port <- sub('^"([0-9]+)".*$', "\\1", h)
      url  <- sub('^"[0-9]+"\\s*:\\s*"([^"]+)".*$', "\\1", h)
      if (nzchar(port) && nzchar(url)) out[[port]] <- url
    }
    out
  }
  
  base_url_for_port <- function(port) {
    sprintf("http://127.0.0.1:%d", as.integer(port))
  }
  
  url_for_app_key <- function(key) {
    local_port <- apps[[key]]$port
    
    route <- switch(
      key,
      catalog    = "catalog",
      gtex       = "gtex",
      nonsyn     = "nonsyn",
      ewasdis    = "ewasdis",
      ewastum    = "ewastum",
      integrator = "integrator",
      hub        = "",
      key
    )
    
    base <- gi_public_base_url()
    
    if (nzchar(base)) {
      if (nzchar(route)) {
        return(paste0(base, "/", route))
      } else {
        return(base)
      }
    }
    
    paste0("http://127.0.0.1:", local_port)
  }
  
  urls_ui <- reactive({
    s <- sid()
    
    list(
      catalog    = paste0(url_for_app_key("catalog"),    "?sid=", s),
      gtex       = paste0(url_for_app_key("gtex"),       "?sid=", s),
      nonsyn     = paste0(url_for_app_key("nonsyn"),     "?sid=", s),
      ewasdis    = paste0(url_for_app_key("ewasdis"),    "?sid=", s),
      ewastum    = paste0(url_for_app_key("ewastum"),    "?sid=", s),
      integrator = paste0(url_for_app_key("integrator"), "?sid=", s)
    )
  })
  
  output$sid_txt <- renderText(sid() %||% "")
  
  output$open_buttons_ui <- renderUI({
    req(sid())
    
    u <- urls_ui()
    
    if (!isTRUE(buttons_ready())) {
      buttons_ready(TRUE)
    }
    
    tagList(
      tags$a(
        "📚 Catalog (MASTER)",
        href = u$catalog,
        target = "_blank",
        class = "btn btn-primary",
        style = btn_style
      ),
      tags$a(
        "🧠 GTEx (SLAVE)",
        href = u$gtex,
        target = "_blank",
        class = "btn btn-success",
        style = btn_style
      ),
      tags$a(
        "🧬 NonSyn (SLAVE)",
        href = u$nonsyn,
        target = "_blank",
        class = "btn btn-success",
        style = btn_style
      ),
      tags$a(
        "🧫 EWASdis (SLAVE)",
        href = u$ewasdis,
        target = "_blank",
        class = "btn btn-success",
        style = btn_style
      ),
      tags$a(
        "🧪 EWAStum (SLAVE)",
        href = u$ewastum,
        target = "_blank",
        class = "btn btn-success",
        style = btn_style
      ),
      tags$a(
        "🧩 Integrator",
        href = u$integrator,
        target = "_blank",
        class = "btn btn-warning",
        style = btn_style
      )
    )
  })
  
###############################
  rv <- reactiveValues(
    pids = read_pids(),
    last_status = NULL,
    manage_msg = "",
    starting = list()
  )
  
  output$diag_out <- renderText({
    paste0(
      "Rscript used: ", RSCRIPT, " (exists=", file.exists(RSCRIPT), ")\n",
      "curl: ", Sys.which("curl"), "\n",
      "lsof: ", Sys.which("lsof"), "\n",
      "kill: ", Sys.which("kill"), "\n",
      "URL_MODE: ", Sys.getenv("GITOOLS_URL_MODE"), "\n",
      "NGROK_URLS_FILE: ", Sys.getenv("GITOOLS_NGROK_URLS_FILE"), " (exists=",
      file.exists(Sys.getenv("GITOOLS_NGROK_URLS_FILE", "")), ")\n"
    )
  })
  
  compute_status <- function() {
    rows <- lapply(names(apps), function(key) {
      info <- apps[[key]]
      url <- paste0(url_for_app_key(key), "/")
      
      pid <- find_listen_pid(info$port)
      port_up <- is.finite(pid)
      
      # si no escolta, no fem curl
      code <- if (port_up) http_status_fast(url) else NA_integer_
      http_up <- is.finite(code) && code >= 200 && code < 500
      
      ts <- rv$starting[[key]] %||% NA
      is_recent_start <- !is.na(ts) && difftime(Sys.time(), ts, units = "secs") < 120
      
      status_txt <- if (port_up) {
        "UP"
      } else if (is_recent_start) {
        "STARTING"
      } else {
        "DOWN"
      }
      
      if (http_up || port_up) {
        rv$starting[[key]] <- NULL
      }
      
      data.frame(
        app = key,
        url = url,
        port = as.integer(info$port),
        http_code = ifelse(is.na(code), 0L, as.integer(code)),
        status = status_txt,
        pid = ifelse(is.na(pid), "", as.character(pid)),
        pid_alive = is.finite(pid),
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, rows)
  }
  
  status_df_for_display <- function(df) {
    if (is.null(df) || !is.data.frame(df)) {
      return(data.frame(
        app = character(0),
        port = integer(0),
        http_code = integer(0),
        status = character(0),
        pid = character(0),
        pid_alive = logical(0),
        url_link = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    if (!nrow(df)) {
      return(data.frame(
        app = character(0),
        port = integer(0),
        http_code = integer(0),
        status = character(0),
        pid = character(0),
        pid_alive = logical(0),
        url_link = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    if (!"url" %in% names(df)) {
      df$url <- ""
    }
    
    df$url_link <- vapply(df$url, function(u) {
      if (is.na(u) || !nzchar(u)) return("")
      as.character(htmltools::tags$a(href = u, target = "_blank", u))
    }, character(1))
    
    df$url <- NULL
    
    keep <- c("app", "port", "http_code", "status", "pid", "pid_alive", "url_link")
    keep <- intersect(keep, names(df))
    df[, keep, drop = FALSE]
  }
  
  refresh_status_table <- function() {
    rv$last_status <- tryCatch(
      compute_status(),
      error = function(e) {
        cat("[HUB][status] compute_status failed:", conditionMessage(e), "\n")
        data.frame(
          app = names(apps),
          port = vapply(apps, function(x) as.integer(x$port), integer(1)),
          http_code = rep.int(0L, length(apps)),
          status = rep.int("DOWN", length(apps)),
          pid = rep.int("", length(apps)),
          pid_alive = rep.int(FALSE, length(apps)),
          url_link = rep("", length(apps)),
          stringsAsFactors = FALSE
        )
      }
    )
    
    if (!isTRUE(status_ready())) status_ready(TRUE)
  }
  
  observeEvent(TRUE, {
    refresh_status_table()
  }, once = TRUE)
  
#  observe({
#    invalidateLater(60000, session)
#    refresh_status_table()
#  })
  
  observeEvent(input$refresh_status, {
    refresh_status_table()
  }, ignoreInit = TRUE)
  
  output$status_table <- DT::renderDT({
    req(!is.null(rv$last_status))
    
    df <- status_df_for_display(rv$last_status)
    
    DT::datatable(
      df,
      rownames = FALSE,
      selection = "none",
      escape = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "tip"
      )
    ) %>%
      DT::formatStyle(
        "status",
        target = "row",
        backgroundColor = DT::styleEqual(
          c("UP", "STARTING", "DOWN"),
          c("#e7f7ee", "#fff7e6", "#fdecea")
        )
      )
  }, server = FALSE)
  
##############################
  observeEvent(input$open_all_slaves, {
    u <- urls_ui()
    insertUI(
      selector = "body",
      where = "beforeEnd",
      ui = HTML(open_many_js(c(
        u$gtex,
        u$nonsyn,
        u$ewasdis,
        u$ewastum
      ))),
      immediate = TRUE
    )
  }, ignoreInit = TRUE)
##############################  
  observeEvent(input$start_one, {
    key <- input$manage_app
    info <- apps[[key]]
    
    showNotification(paste0("Starting ", key, "…"), type="message", duration=2)
    
    if (!dir.exists(info$dir)) {
      rv$manage_msg <- sprintf("[%s] DIR NOT FOUND: %s", key, info$dir)
      refresh_status_table()
      return()
    }
    
    if (is_port_listening(info$port)) {
      rv$manage_msg <- sprintf("[%s] already running (port %d).", key, info$port)
      refresh_status_table()
      return()
    }
    
    rv$starting[[key]] <- Sys.time()
    
    res <- tryCatch(start_app_bg(key, info$dir, info$port), error=function(e) e)
    if (inherits(res, "error")) {
      rv$manage_msg <- sprintf("[%s] START ERROR: %s", key, conditionMessage(res))
      showNotification(paste0(key, ": start error (check logs)"), type="error", duration=5)
    } else {
      rv$pids[[key]] <- res$pid
      write_pids(reactiveValuesToList(rv)$pids)
      rv$manage_msg <- sprintf("[%s] launched.\nerr_log=%s\npid(listen)=%s",
                               key, res$err_log, res$pid)
      if (is.na(res$pid)) {
        showNotification(paste0(key, " still STARTING. Check err_log if it stays DOWN."),
                         type="warning", duration=6)
      }
    }
    
    refresh_status_table()
  }, ignoreInit = TRUE)
  
  observeEvent(input$kill_one, {
    key <- input$manage_app
    info <- apps[[key]]
    
    ok <- kill_by_port(info$port)
    rv$starting[[key]] <- NULL
    rv$manage_msg <- sprintf("[%s] kill port=%d -> %s", key, info$port, ok)
    
    refresh_status_table()
  }, ignoreInit = TRUE)
  
  observeEvent(input$restart_one, {
    key <- input$manage_app
    info <- apps[[key]]
    
    kill_by_port(info$port)
    Sys.sleep(0.4)
    
    rv$starting[[key]] <- Sys.time()
    
    res <- tryCatch(start_app_bg(key, info$dir, info$port), error=function(e) e)
    if (inherits(res, "error")) {
      rv$manage_msg <- sprintf("[%s] RESTART ERROR: %s", key, conditionMessage(res))
      showNotification(paste0(key, ": restart error (check logs)"), type="error", duration=6)
    } else {
      rv$pids[[key]] <- res$pid
      write_pids(reactiveValuesToList(rv)$pids)
      rv$manage_msg <- sprintf("[%s] restarted.\nerr_log=%s\npid(listen)=%s",
                               key, res$err_log, res$pid)
      if (is.na(res$pid)) {
        showNotification(paste0(key, " still STARTING. Check err_log if it stays DOWN."),
                         type="warning", duration=6)
      }
    }
    
    refresh_status_table()
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$start_all, {
    hub_log("[hub] start_all clicked value=", input$start_all, " @ ", format(Sys.time()))
    
    cat(sprintf("[hub] start_all clicked @ %s\n", format(Sys.time())),
        file = HUB_LOG, append = TRUE)
    
    msgs <- character(0)
    
    withProgress(message = "Starting apps…", value = 0, {
      keys <- setdiff(names(apps), "hub")
      n <- length(keys)
      
      for (i in seq_along(keys)) {
        key <- keys[i]
        info <- apps[[key]]
        
        incProgress(1/n, detail = paste0(key, " (port ", info$port, ")"))
        
        if (!dir.exists(info$dir)) {
          msgs <- c(msgs, sprintf("[%s] DIR NOT FOUND: %s", key, info$dir))
          next
        }
        
        if (is_port_listening(info$port)) next
        
        rv$starting[[key]] <- Sys.time()
        
        cat(sprintf("[hub] start_all -> %s port=%d dir=%s exists=%s\n",
                    key, info$port, info$dir, dir.exists(info$dir)),
            file = HUB_LOG, append = TRUE)
        
        res <- tryCatch(start_app_bg(key, info$dir, info$port), error=function(e) e)
        if (inherits(res, "error")) {
          msgs <- c(msgs, sprintf("[%s] START ERROR: %s", key, conditionMessage(res)))
        } else {
          rv$pids[[key]] <- res$pid
          msgs <- c(msgs, sprintf("[%s] launched (pid listen=%s)", key, res$pid))
        }
        
        Sys.sleep(0.2)
      }
    })
    
    write_pids(reactiveValuesToList(rv)$pids)
    rv$manage_msg <- paste(c("[ALL] start_all done.", msgs), collapse = "\n")
    refresh_status_table()
    
    showNotification("Start all finished. Check Status table.", type="message", duration=3)
  }, ignoreInit = TRUE)
  
  observeEvent(input$kill_all, {
    
    keys <- setdiff(names(apps), "hub")  # <-- IMPORTANT: no matar el hub
    
    for (key in keys) {
      info <- apps[[key]]
      kill_by_port(info$port)
      rv$starting[[key]] <- NULL
    }
    
    rv$manage_msg  <- "[ALL] kill_all executed (hub preserved)."
    refresh_status_table()
    
    showNotification("Killed all slave apps (Hub stays running).", type="message", duration=3)
  }, ignoreInit = TRUE)
  
  output$manage_out <- renderText(rv$manage_msg)
  
  # =============================================================================
  # Hub info modal (input$info_00) — GItools Quick Guide (EN)
  # 
  # =============================================================================
  
  observeEvent(input$info_00, {
    
    info_ui <- tags$div(
      style = "max-height:70vh; overflow-y:auto; padding-right:8px; line-height:1.35;",
      
      tags$h3("GItools Hub — Quick Guide", style = "margin-top:0;"),
      
      tags$p(
        tags$b("Goal:"),
        " Define the shared genomic regions once in ", tags$b("Catalog Inspector"),
        " (master) and reuse them across the other GItools apps. ",
        tags$b("Integrator Inspector"),
        " acts as the cross-app integration layer, including LD evidence, once the complementary inspectors have generated their exports."
      ),
      
      tags$div(
        style = "text-align:center; margin:10px 0 14px 0;",
        tags$img(
          src = "gi_tools_ideogram.png",
          style = "max-width:100%; width:700px; height:auto; border:1px solid #e5e5e5; border-radius:10px;"
        )
      ),
      tags$hr(),
      
      tags$h4("1) Start with Catalog Inspector (Master)"),
      tags$ol(
        tags$li("Open ", tags$b("Catalog Inspector"), " from the Hub."),
        tags$li(tags$b("Load Catalog/GWAS data"), " (your selected dataset)."),
        tags$li("Choose how you want to define the shared regions:"),
        tags$ul(
          tags$li(
            tags$b("Standard mode:"),
            " build clusters from GWAS hits using the clustering controls."
          ),
          tags$li(
            tags$b("User-defined mode:"),
            " upload a file with predefined genomic intervals (cluster_id / chr / start / end)."
          )
        ),
        tags$li(
          "If using ", tags$b("user-defined intervals"),
          ", the uploaded intervals become the canonical clusters, and the selected GWAS threshold is used to assign significant GWAS hits to those intervals."
        ),
        tags$li("Explore hits using the Manhattan plot and summary tables."),
        tags$li("Adjust thresholds or clustering parameters as needed."),
        tags$li("Click ", tags$b("Build clusters"), "."),
        tags$li("Confirm that clusters appear in the cluster table and plots.")
      ),
      tags$p(
        tags$i(
          "Once clusters are defined, the rest of the apps inspect the same canonical cluster IDs and genomic intervals, regardless of whether they were generated from GWAS clustering or uploaded as user-defined regions."
        )
      ),
      
      tags$hr(),
      
      tags$h4("2) Open the complementary inspectors"),
      tags$ul(
        tags$li(
          "After defining clusters in ", tags$b("Catalog Inspector"),
          ", open the other inspectors from the Hub."
        ),
        tags$li(
          "Each inspector loads its own reference dataset and ",
          tags$b("reuses the same canonical cluster framework"), "."
        ),
        tags$li(
          "The goal is to inspect ", tags$b("the same loci across complementary biological evidence layers"),
          " rather than rebuilding each analysis independently."
        ),
        tags$li(
          "These inspectors generate the harmonized evidence exports that will later be consolidated by ",
          tags$b("Integrator Inspector"), "."
        )
      ),
      
      tags$hr(),
      
      tags$h4("3) Integrate with Integrator Inspector"),
      tags$ul(
        tags$li(
          "Once the complementary inspectors have been opened and their evidence has been generated, open ",
          tags$b("Integrator Inspector"), "."
        ),
        tags$li(
          tags$b("Integrator Inspector"),
          " is the central cross-app layer that consolidates evidence from the other inspectors."
        ),
        tags$li(
          "It reads the harmonized exports generated across apps and combines them into shared summaries for the ",
          tags$b("same canonical clusters"), "."
        ),
        tags$li(
          "This includes integrated interpretation of ",
          tags$b("Catalog, GTEx, NonSyn, EWAS, and LD evidence"),
          "."
        ),
        tags$li(
          "Use it to review shared genes, shared traits/diseases, consensus loci, prioritized candidates, and LD support in one place."
        )
      ),
      
      tags$hr(),
      
      tags$h4("What each app is for"),
      tags$div(
        
        tags$h5("• Catalog Inspector"),
        tags$ul(
          tags$li("Explore GWAS/Catalog associations."),
          tags$li("Create canonical clusters shared across GItools."),
          tags$li("Support both GWAS-derived clustering and uploaded user-defined cluster intervals."),
          tags$li("Define the genomic backbone for downstream inspection.")
        ),
        
        tags$h5("• GTEx Inspector"),
        tags$ul(
          tags$li("Inspect GTEx eQTL signals within the same clusters."),
          tags$li("Compare tissue-specific regulatory evidence."),
          tags$li("Highlight genes with expression-based support.")
        ),
        
        tags$h5("• NonSyn Inspector"),
        tags$ul(
          tags$li("Annotate clusters with coding and non-synonymous variants."),
          tags$li("Review functional impact predictions and gene context."),
          tags$li("Generate candidate gene/variant evidence for integration.")
        ),
        
        tags$h5("• EWAS Disease Inspector"),
        tags$ul(
          tags$li("Evaluate disease vs control methylation patterns per cluster."),
          tags$li("Summarize EWAS hits by disease and genomic region."),
          tags$li(
            tags$b("Click actions:"),
            " disease → CpG violins for that disease; ",
            "probe → probe across diseases; ",
            "hyper/hypo → probe+disease detail (case/control)."
          )
        ),
        
        tags$h5("• EWAS Tumor Inspector"),
        tags$ul(
          tags$li("Evaluate tumor vs normal methylation patterns per cluster."),
          tags$li("Summarize EWAS hits by cancer and genomic region."),
          tags$li("Use the same cluster framework to compare cancer-related methylation evidence.")
        ),
        
        tags$h5("• Integrator Inspector"),
        tags$ul(
          tags$li("Collect harmonized evidence exported by the other apps."),
          tags$li("Combine gene-level and trait/disease-level bridge tables across inspectors."),
          tags$li("Provide a unified view of shared genes, shared traits, consensus loci, and prioritized candidates."),
          tags$li("Include ", tags$b("LD evidence"), " as one component of the broader integration workflow.")
        )
      ),
      
      tags$hr(),
      
      tags$h4("Integrator Inspector and LD evidence"),
      tags$ul(
        tags$li(
          tags$b("Integrator Inspector"),
          " is the central destination for combined interpretation."
        ),
        tags$li(
          "Instead of using a standalone LD destination, LD is now handled inside ",
          tags$b("Integrator Inspector"), " as the ", tags$b("LD evidence"), " section."
        ),
        tags$li(
          "This makes it possible to interpret LD support together with Catalog, GTEx, NonSyn, and EWAS evidence."
        ),
        tags$li(
          "Each inspector exports harmonized bridge tables to a shared location, and Integrator Inspector reads those exports to build the integrated summaries."
        )
      ),
      
      tags$hr(),
      
      tags$h4("Recommended fast workflow"),
      tags$ol(
        tags$li(
          tags$b("Hub → Catalog Inspector:"),
          " load data and define clusters either from GWAS-based clustering or from uploaded user-defined intervals."
        ),
        tags$li(
          tags$b("Hub → open the other inspectors:"),
          " review complementary evidence for the same canonical regions."
        ),
        tags$li(
          tags$b("Hub → Integrator Inspector:"),
          " consolidate all available evidence, including ", tags$b("LD evidence"),
          ", shared genes, shared traits/diseases, consensus loci, and prioritized candidates."
        )
      ),
      
      tags$hr(),
      
      tags$h4("Practical notes"),
      tags$ul(
        tags$li(
          tags$b("Clusters are the common language."),
          " If an inspector looks empty, first confirm that clusters were defined in Catalog Inspector."
        ),
        tags$li(
          "If you rebuild clusters or upload a different user-defined interval set in Catalog, reload the other apps so they work with the updated cluster set."
        ),
        tags$li(
          "Cluster IDs (cluster_id) should remain canonical and consistent across all apps."
        ),
        tags$li(
          "In user-defined mode, uploaded intervals define the regions, while the selected GWAS threshold is used only to assign GWAS hits to those intervals."
        ),
        tags$li(
          "Integrator summaries depend on the bridge export files written by each app. If integrated tables are empty, check that the corresponding bridge files are being generated correctly."
        ),
        tags$li(
          "For LD evidence, if results are empty, verify PLINK paths, chromosome naming, population selection, and that the region contains enough variants."
        )
      )
    )
    
    showModal(modalDialog(
      title = tags$span("ℹ️ GItools Hub — Help"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      info_ui
    ))
    
  }, ignoreInit = TRUE)
  
  output$dl_example_zip <- downloadHandler(
    filename = function() {
      paste0("gitools_example_files_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      src <- file.path(APPS_DIR, "GItools_Hub", "www", "example_files.zip")
      if (!file.exists(src)) stop("Example file not found: ", src)
      file.copy(src, file, overwrite = TRUE)
    }
  )
  
}

shinyApp(ui, server)
