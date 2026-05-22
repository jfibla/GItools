`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) a else b

gi_nav_base_url <- function(port) {
  sprintf("http://127.0.0.1:%d", as.integer(port))
}

gi_nav_urls <- function(sid) {
  sid <- as.character(sid %||% "")
  
  list(
    hub        = paste0(gi_nav_base_url(7101), "/?sid=", sid),
    catalog    = paste0(gi_nav_base_url(7201), "/?sid=", sid),
    gtex       = paste0(gi_nav_base_url(7202), "/?sid=", sid),
    nonsyn     = paste0(gi_nav_base_url(7203), "/?sid=", sid),
    ewastum    = paste0(gi_nav_base_url(7204), "/?sid=", sid),
    ewasdis    = paste0(gi_nav_base_url(7205), "/?sid=", sid),
    integrator = paste0(gi_nav_base_url(7206), "/?sid=", sid)
  )
}

gi_top_nav_ui <- function(sid, active = NULL) {
  u <- gi_nav_urls(sid)
  
  nav_btn <- function(label, href, key, active_key = NULL) {
    cls <- if (identical(key, active_key)) {
      "background:#1A4E8A; color:white; border:1px solid #1A4E8A;"
    } else {
      "background:#f7f7f7; color:#222; border:1px solid #d9d9d9;"
    }
    
    tags$a(
      label,
      href = href,
      style = paste(
        "display:inline-block;",
        "padding:6px 10px;",
        "margin-right:6px;",
        "margin-bottom:6px;",
        "border-radius:8px;",
        "text-decoration:none;",
        cls
      )
    )
  }
  
  tags$div(
    style = paste(
      "padding:8px 10px;",
      "margin-bottom:12px;",
      "background:#fbfbfb;",
      "border:1px solid #e6e6e6;",
      "border-radius:10px;"
    ),
    nav_btn("🧭 Hub",        u$hub,        "hub",        active),
    nav_btn("📚 Catalog",    u$catalog,    "catalog",    active),
    nav_btn("🧠 GTEx",       u$gtex,       "gtex",       active),
    nav_btn("🧬 NonSyn",     u$nonsyn,     "nonsyn",     active),
    nav_btn("🧫 EWASDis",    u$ewasdis,    "ewasdis",    active),
    nav_btn("🧪 EWAStum",    u$ewastum,    "ewastum",    active),
    nav_btn("🧩 Integrator", u$integrator, "integrator", active)
  )
}

################################################################################