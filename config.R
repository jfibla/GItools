# GItools/config.R
# Config portable para todas las apps

gi_find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(d, "config.R")) && dir.exists(file.path(d, "_shared"))) return(d)
    d2 <- dirname(d)
    if (identical(d2, d)) break
    d <- d2
  }
  stop("GItools root not found. Put config.R + _shared at repo root, or set env GITOOLS_ROOT.")
}

gi_cfg <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    
    root <- Sys.getenv("GITOOLS_ROOT", unset = "")
    if (!nzchar(root)) root <- gi_find_root(getwd())
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    
    shared <- Sys.getenv("GITOOLS_SHARED", unset = file.path(root, "_shared"))
    shared <- normalizePath(shared, winslash = "/", mustWork = TRUE)
    
    # Recursos grandes (NO van en git)
    # ✅ IMPORTANT: en el teu repo estan sota root/app/Inspector_resources
    resources <- Sys.getenv("GITOOLS_RESOURCES", unset = file.path(root, "app", "Inspector_resources"))
    resources <- normalizePath(resources, winslash = "/", mustWork = FALSE)
    
    # POP keep-files directory
    pop_dir <- Sys.getenv("GITOOLS_POP_DIR", unset = "")
    if (!nzchar(pop_dir)) {
      # ✅ IMPORTANT: en el teu layout real és LD_resources/POP
      pop_dir <- if (nzchar(resources)) file.path(resources, "LD_resources", "POP") else ""
    }
    pop_dir <- normalizePath(pop_dir, winslash = "/", mustWork = FALSE)
    
    cache <<- list(
      root      = root,
      shared    = shared,
      resources = resources,
      pop_dir   = pop_dir
    )
    cache
  }
})
