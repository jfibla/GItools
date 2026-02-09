# /Volumes/DISK1TB/Inspector_app_slaves/GItools_Hub/app.R
#
# to kill all processes and remove logs:
#  pkill -f "shiny::runApp"
#  rm -rf /Volumes/DISK1TB/Inspector_app_slaves/_logs
#  mkdir -p /Volumes/DISK1TB/Inspector_app_slaves/_logs
#
# to run hub:
# shiny::runApp("/Volumes/DISK1TB/Inspector_app_slaves_github/GItools/app/GItools_Hub",
#              port=7101, host="127.0.0.1", launch.browser=TRUE)

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

# ---- Ports (portable) ----
APP_PORT_BASE <- 7200 

apps <- list(
  catalog = list(name="📚 Catalog Inspector", dir=file.path(APPS_DIR, "Catalog_inspector"), port=APP_PORT_BASE + 1),
  gtex    = list(name="🧠 GTEx Inspector",    dir=file.path(APPS_DIR, "GTEX_inspector"),    port=APP_PORT_BASE + 2),
  nonsyn  = list(name="🧬 NonSyn Inspector",  dir=file.path(APPS_DIR, "NonSyn_Inspector"),  port=APP_PORT_BASE + 3),
  ewastum = list(name="🧪 EWAStum Inspector", dir=file.path(APPS_DIR, "EWAS_cancer"),       port=APP_PORT_BASE + 4),
  ewasdis = list(name="🧫 EWASDis Inspector", dir=file.path(APPS_DIR, "EWAS_disease"),      port=APP_PORT_BASE + 5),
  ld      = list(name="🔗 LD Inspector",      dir=file.path(APPS_DIR, "LD_Inspector"),      port=APP_PORT_BASE + 6)
)


# Logs (en temp por defecto, siempre escribible)
LOG_DIR <- file.path(APPS_DIR, "_logs")   # -> GItools/app/_logs
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

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
RSCRIPT  <- file.path(R.home("bin"), "Rscript")

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
http_statusXXXX <- function(url) {
  if (!has_cmd("curl")) return(NA_integer_)
  
  args <- c(
    "-s", "-o", "/dev/null",
    "-w", "%{http_code}",
    "--connect-timeout", "3",
    "--max-time", "5",
    url
  )
  out <- suppressWarnings(tryCatch(system2("curl", args, stdout = TRUE, stderr = TRUE),
                                   error = function(e) ""))
  
  code <- suppressWarnings(as.integer(as.character(out[1] %||% "")))
  if (!is.finite(code)) NA_integer_ else code
}

http_status_fast <- function(url) {
  if (!has_cmd("curl")) return(NA_integer_)
  
  args <- c(
    "-s", "-o", "/dev/null",
    "-w", "%{http_code}",
    "--connect-timeout", "0.25",
    "--max-time", "0.8",
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
    "options(shiny.port = NULL)",
    "options(shiny.host = NULL)",
    sprintf('Sys.setenv(GITOOLS_LOG_DIR = "%s")', gsub("\\\\", "/", log_dir)),  # opcional
    sprintf("setwd(%s)", app_dir_q),
    sprintf('cat("[runner] starting %s on port %d\\n")', key, port),
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
  pid <- find_listen_pid(port)
  if (!is.finite(pid)) return(FALSE)
  suppressWarnings(system2("kill", c("-15", as.character(pid)), stdout = NULL, stderr = NULL))
  TRUE
}

open_many_js <- function(url_vec) {
  js <- paste(vapply(url_vec, function(u) sprintf("window.open('%s','_blank');", u), ""),
              collapse = "")
  sprintf("<script>%s</script>", js)
}

ui <- fluidPage(
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
        column(4, selectInput("manage_app", "App", choices = names(apps), selected = "catalog")),
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
  
  # --- SID builder ---
  sid <- reactiveVal(NULL)
  observeEvent(TRUE, {
    q <- shiny::parseQueryString(session$clientData$url_search %||% "")
    sid0 <- as.character(q$sid %||% "")
    if (!nzchar(sid0)) sid0 <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample(1000:9999, 1))
    sid(sid0)
  }, once = TRUE)
  
#  urls_ui <- reactive({
#    s <- sid()
#    list(
#      catalog = paste0("/catalog/?sid=", s),
#      gtex    = paste0("/gtex/?sid=", s),
#      nonsyn  = paste0("/nonsyn/?sid=", s),
#      ewastum = paste0("/ewastum/?sid=", s),
#      ewasdis = paste0("/ewasdis/?sid=", s),
#      ld      = paste0("/ld/?sid=", s)
#    )
#  })
  
  urls_ui <- reactive({
    s <- sid()
    
    list(
      catalog = sprintf("http://127.0.0.1:%d/?sid=%s", apps$catalog$port, s),
      gtex    = sprintf("http://127.0.0.1:%d/?sid=%s", apps$gtex$port,    s),
      nonsyn  = sprintf("http://127.0.0.1:%d/?sid=%s", apps$nonsyn$port,  s),
      ewastum = sprintf("http://127.0.0.1:%d/?sid=%s", apps$ewastum$port, s),
      ewasdis = sprintf("http://127.0.0.1:%d/?sid=%s", apps$ewasdis$port, s),
      ld      = sprintf("http://127.0.0.1:%d/?sid=%s", apps$ld$port,      s)
    )
  })
  
  output$sid_txt <- renderText(sid() %||% "")
  
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
      "kill: ", Sys.which("kill"), "\n"
    )
  })
  
  # ✅ This controls the banner closing exactly when open_buttons_ui is rendered
  banner_hidden <- reactiveVal(FALSE)
  
  output$open_buttons_ui <- renderUI({
    u <- urls_ui()
    
    # mark ready, but DO NOT hide banner here
    if (!isTRUE(buttons_ready())) buttons_ready(TRUE)
    
    tagList(
      tags$a("📚 Catalog (MASTER)", href=u$catalog, target="_blank", class="btn btn-primary", style=btn_style),
      tags$a("🧠 GTEx (SLAVE)",     href=u$gtex,    target="_blank", class="btn btn-primary", style=btn_style),
      tags$a("🧬 NonSyn (SLAVE)",   href=u$nonsyn,  target="_blank", class="btn btn-primary", style=btn_style),
      tags$a("🧪 EWASDis (SLAVE)",  href=u$ewasdis, target="_blank", class="btn btn-success", style=btn_style),
      tags$a("🧫 EWAStum (SLAVE)",  href=u$ewastum, target="_blank", class="btn btn-success", style=btn_style),
      tags$a("🔗 LD",              href=u$ld,      target="_blank", class="btn btn-warning", style=btn_style)
    )
  })
  # --- Status computation (light + stable) ---
  compute_status <- function() {
    rows <- lapply(names(apps), function(key) {
      info <- apps[[key]]
      url  <- sprintf("http://127.0.0.1:%d/", info$port)
      
      pid <- find_listen_pid(info$port)
      port_up <- is.finite(pid)
      
      # ✅ si NO escolta, NO fem curl (evita bloqueig)
      code <- if (port_up) http_status_fast(url) else NA_integer_
      
      http_up <- is.finite(code) && code >= 200 && code < 500
      
      ts <- rv$starting[[key]] %||% NA
      is_recent_start <- !is.na(ts) && difftime(Sys.time(), ts, units="secs") < 120
      
      status_txt <- if (http_up) {
        "UP"
      } else if (port_up) {
        "UP"        # listening => UP (estable)
      } else if (is_recent_start) {
        "STARTING"
      } else {
        "DOWN"
      }
      
      if (http_up || port_up) rv$starting[[key]] <- NULL
      
      data.frame(
        app = key,
        url = url,
        port = info$port,
        http_code = ifelse(is.na(code), 0L, code),
        status = status_txt,
        pid = ifelse(is.na(pid), "", as.character(pid)),
        pid_alive = is.finite(pid),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
  
  
  output$status_table <- renderDT({
    df <- rv$last_status
    if (is.null(df)) {
      df <- tryCatch(compute_status(), error = function(e) {
        data.frame(
          app = names(apps),
          url = vapply(apps, function(x) sprintf("http://127.0.0.1:%d/", x$port), ""),
          port = vapply(apps, `[[`, 0L, "port"),
          http_code = 0L,
          status = paste0("ERR: ", conditionMessage(e)),
          pid = "",
          pid_alive = FALSE,
          stringsAsFactors = FALSE
        )
      })
    }
    
    
    datatable(df, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatStyle(
        "status", target="row",
        backgroundColor = styleEqual(
          c("UP","STARTING","DOWN"),
          c("#e7f7ee","#fff7e6","#fdecea")
        )
      )
  })
  
  # light polling (avoid heavy loops over ngrok)
  observe({
    invalidateLater(12000, session)
    
   # rv$last_status <- compute_status()
    rv$last_status <- tryCatch(compute_status(), error = function(e) rv$last_status)
    
    # mark first status computed
    if (!isTRUE(status_ready())) status_ready(TRUE)
  })
  
  observe({
    # Close banner only once, when both UI buttons + first status are ready
    if (isTRUE(banner_hidden())) return()
    req(isTRUE(buttons_ready()), isTRUE(status_ready()))
    
    elapsed <- as.numeric(difftime(Sys.time(), banner_started, units = "secs"))
    if (elapsed < BANNER_MIN_SECS) {
      invalidateLater(as.integer((BANNER_MIN_SECS - elapsed) * 1000), session)
      return()
    }
    
    session$sendCustomMessage("hub_banner", list(show = FALSE))
    banner_hidden(TRUE)
  })
  
  
  observeEvent(input$refresh_status, {
    rv$last_status <- compute_status()
  }, ignoreInit = TRUE)
  
  observeEvent(input$open_4x, {
    u <- urls_ui()
    insertUI(selector="body", where="beforeEnd", ui=HTML(open_many_js(c(u$gtex, u$nonsyn,u$ewasdis, u$ewastum))), immediate=TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$open_gtex_nonsyn, {
    u <- urls_ui()
    insertUI(selector="body", where="beforeEnd", ui=HTML(open_many_js(c(u$gtex, u$nonsyn))), immediate=TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$open_catalog_ld, {
    u <- urls_ui()
    insertUI(selector="body", where="beforeEnd", ui=HTML(open_many_js(c(u$catalog, u$ld))), immediate=TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$open_ewas_both, {
    u <- urls_ui()
    insertUI(selector="body", where="beforeEnd", ui=HTML(open_many_js(c(u$ewasdis, u$ewastum))), immediate=TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$start_one, {
    key <- input$manage_app
    info <- apps[[key]]
    
    showNotification(paste0("Starting ", key, "…"), type="message", duration=2)
    
    if (!dir.exists(info$dir)) {
      rv$manage_msg <- sprintf("[%s] DIR NOT FOUND: %s", key, info$dir)
      rv$last_status <- compute_status()
      return()
    }
    
    if (is_port_listening(info$port)) {
      rv$manage_msg <- sprintf("[%s] already running (port %d).", key, info$port)
      rv$last_status <- compute_status()
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
    
    rv$last_status <- compute_status()
  }, ignoreInit = TRUE)
  
  observeEvent(input$kill_one, {
    key <- input$manage_app
    info <- apps[[key]]
    
    ok <- kill_by_port(info$port)
    rv$starting[[key]] <- NULL
    rv$manage_msg <- sprintf("[%s] kill port=%d -> %s", key, info$port, ok)
    
    rv$last_status <- compute_status()
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
    
    rv$last_status <- compute_status()
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$start_all, {
    hub_log("[hub] start_all clicked value=", input$start_all, " @ ", format(Sys.time()))
    
    cat(sprintf("[hub] start_all clicked @ %s\n", format(Sys.time())),
        file = HUB_LOG, append = TRUE)
    
    msgs <- character(0)
    
    withProgress(message = "Starting apps…", value = 0, {
      keys <- names(apps)
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
    rv$last_status <- compute_status()
    
    showNotification("Start all finished. Check Status table.", type="message", duration=3)
  }, ignoreInit = TRUE)
  
  observeEvent(input$kill_all, {
    for (key in names(apps)) {
      info <- apps[[key]]
      kill_by_port(info$port)
      rv$starting[[key]] <- NULL
    }
    rv$manage_msg <- "[ALL] kill_all executed."
    rv$last_status <- compute_status()
  }, ignoreInit = TRUE)
  
  output$manage_out <- renderText(rv$manage_msg)
  
  # =============================================================================
  # Hub info modal (input$info_00) — GItools Quick Guide (EN)
  # Includes: per-app LD module + standalone LD Inspector using candidates.zip
  # =============================================================================
  
  observeEvent(input$info_00, {
    
    info_ui <- tags$div(
      style = "max-height:70vh; overflow-y:auto; padding-right:8px; line-height:1.35;",
      
      tags$h3("GItools Hub — Quick Guide", style="margin-top:0;"),
      
      tags$p(
        tags$b("Goal:"),
        "Build clusters once in ", tags$b("Catalog Inspector"),
        " (master) and reuse them across the other GItools apps (slaves)."
      ),
      
      tags$hr(),
      
      tags$h4("1) Start with Catalog Inspector (Master)"),
      tags$ol(
        tags$li("Open ", tags$b("Catalog Inspector"), " from the Hub."),
        tags$li(tags$b("Load Catalog/GWAS data"), " (your selected dataset)."),
        tags$li("Explore hits (Manhattan plot / tables)."),
        tags$li("Set parameters (thresholds, method, etc.)."),
        tags$li("Click ", tags$b("Build clusters"), "."),
        tags$li("Confirm clusters appear in the cluster table and plots.")
      ),
      tags$p(tags$i(
        "Once clusters are built, the rest of the apps can import the same cluster IDs and genomic intervals."
      )),
      
      tags$hr(),
      
      tags$h4("2) Open the other apps (Slaves)"),
      tags$ul(
        tags$li(
          "When you open another app from the Hub, it loads its own reference dataset and ",
          tags$b("imports the clusters built in Catalog Inspector"), "."
        ),
        tags$li(
          "You usually do not need to rebuild clusters in each app — the idea is that all apps inspect ",
          tags$b("the same regions"), "."
        )
      ),
      
      tags$hr(),
      
      tags$h4("What each app is for"),
      tags$div(
        
        tags$h5("• Catalog Inspector"),
        tags$ul(
          tags$li("Explore GWAS/Catalog associations."),
          tags$li("Create canonical clusters shared across GItools."),
          tags$li("Export clusters/candidates for downstream steps.")
        ),
        
        tags$h5("• GTEx Inspector"),
        tags$ul(
          tags$li("Inspect GTEx eQTL signals within the same clusters."),
          tags$li("Compare tissue-specific regulatory evidence."),
          tags$li("Prioritize candidate genes based on eQTL support.")
        ),
        
        tags$h5("• NonSyn Inspector (dbNSFP)"),
        tags$ul(
          tags$li("Annotate clusters with coding / non-synonymous variants."),
          tags$li("Review functional impact predictions and gene context."),
          tags$li("Generate candidate variant lists per cluster.")
        ),
        
        tags$h5("• EWAS Disease Inspector"),
        tags$ul(
          tags$li("Evaluate disease vs control methylation patterns per cluster."),
          tags$li("Build EWAS hits per cluster and disease."),
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
          tags$li("Build EWAS hits per cluster and cancer."),
          tags$li("Same click-to-plot workflow as Disease, using cancer labels.")
        )
      ),
      
      tags$hr(),
      
      tags$h4("LD (Linkage Disequilibrium) options"),
      
      tags$h5("A) Built-in LD module inside each app"),
      tags$ul(
        tags$li(
          "Each GItools app includes an embedded ", tags$b("LD module"),
          " to visualize LD structure within a selected cluster region (e.g., r² or D')."
        ),
        tags$li(
          "Typical flow: select a ", tags$b("cluster"), " → choose a ", tags$b("population"),
          " and reference genotype panel (PLINK) → run LD → view LD triangle and (optionally) LD blocks."
        ),
        tags$li(
          "Use it to refine candidate regions, identify LD blocks, and interpret whether multiple hits may reflect the same signal."
        )
      ),
      
      tags$h5("B) Standalone LD Inspector (candidates.zip)"),
      tags$ul(
        tags$li(
          "There is also an ", tags$b("independent LD Inspector"),
          " that runs outside the individual apps."
        ),
        tags$li(
          "From any app, download the ", tags$b("candidates.zip"),
          " export (clusters + candidate hits) and load it into the LD Inspector."
        ),
        tags$li(
          "This standalone tool can visualize LD structure within the selected cluster region ",
          tags$b("for all candidate hits together"),
          " (useful for cross-app candidate consolidation)."
        )
      ),
      
      tags$hr(),
      
      tags$h4("Recommended fast workflow"),
      tags$ol(
        tags$li(tags$b("Hub → Catalog Inspector:"), " load data → build clusters."),
        tags$li(tags$b("Hub → open any other app:"), " it loads its dataset + imports clusters."),
        tags$li(
          "Inspect each cluster across evidence layers: GWAS → eQTL → coding impact → EWAS → ",
          tags$b("LD (built-in or standalone)"), "."
        )
      ),
      
      tags$hr(),
      
      tags$h4("Practical notes"),
      tags$ul(
        tags$li(
          tags$b("Clusters are the common language."),
          " If a slave app looks empty, confirm you built clusters in Catalog and the slave app has loaded them."
        ),
        tags$li(
          "If you rebuild clusters in Catalog, reload the other apps so they use the updated cluster set."
        ),
        tags$li(
          "Cluster IDs (cluster_id) should remain canonical/consistent across all apps."
        ),
        tags$li(
          "For LD: if results are empty, verify PLINK paths, chromosome naming, population selection, and that the region contains enough variants."
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
  
  
  
}

shinyApp(ui, server)
