# R/gi_slave_canonical_gtex.R
# ------------------------------------------------------------
# GItools CANONICAL SLAVE SYNC (GTEx slave)
# - Reads SID from input$gi_qs, url_search, url_hash
# - Polls state JSON (mtime)
# - Loads GWAS + CLUSTERS RDS
# - Applies UI params (best effort)
# - Exposes: gwas_shared(), clusters_shared(), state_shared(), sid()
# ------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

gi_sanitize_sid <- function(s) {
  s <- as.character(s %||% "")
  s <- gsub("[^0-9A-Za-z_\\-]", "", s)
  if (!nzchar(s)) NULL else s
}

sid_ok <- function(x) {
  if (is.null(x)) return(FALSE)
  x <- as.character(x)
  if (length(x) != 1) return(FALSE)
  x <- trimws(x)
  isTRUE(nzchar(x)) && !is.na(x)
}

# Prefer options() roots (portable + hub)
gi_get_state_root <- function() {
  x <- getOption("gi_state_root", "")
  if (nzchar(x)) return(x)
  # fallback very conservative
  ""
}

gi_state_paths <- function(sid) {
  root <- gi_get_state_root()
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  
  list(
    json = file.path(root, paste0("state_",    sid, ".json")),
    gwas = file.path(root, paste0("gwas_",     sid, ".rds")),
    clus = file.path(root, paste0("clusters_", sid, ".rds")),
    params = file.path(root, paste0("params_",   sid, ".rds"))  # <-- NEW
  )
}

gi_read_state <- function(sid) {
  p <- gi_state_paths(sid)$json
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::fromJSON(p), error = function(e) NULL)
}

gi_std_gwas <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  nms <- names(df)
  
  if (!"CHR"  %in% nms && "chr"  %in% nms) df$CHR  <- df$chr
  if (!"BP"   %in% nms && "bp"   %in% nms) df$BP   <- df$bp
  if (!"Pval" %in% nms && "pval" %in% nms) df$Pval <- df$pval
  if (!"logp" %in% nms && "LOGP" %in% nms) df$logp <- df$LOGP
  if (!"snp"  %in% nms && "SNP"  %in% nms) df$snp  <- df$SNP
  if (!"rsid" %in% nms && "RSID" %in% nms) df$rsid <- df$RSID
  
  df$CHR  <- suppressWarnings(as.integer(df$CHR))
  df$BP   <- suppressWarnings(as.integer(df$BP))
  df$Pval <- suppressWarnings(as.numeric(df$Pval))
  if (!"logp" %in% names(df)) df$logp <- -log10(df$Pval)
  
  if (!"snp"  %in% names(df)) df$snp  <- NA_character_
  if (!"rsid" %in% names(df)) df$rsid <- NA_character_
  
  df <- df[is.finite(df$CHR) & is.finite(df$BP) & is.finite(df$Pval) & df$Pval > 0, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  
  df$snp  <- ifelse(is.na(df$snp)  | !nzchar(df$snp),  paste0("chr", df$CHR, ":", df$BP), as.character(df$snp))
  df$rsid <- as.character(df$rsid)
  
  df[, intersect(c("CHR","BP","snp","rsid","Pval","logp","BPcum"), names(df)), drop = FALSE]
}

gi_std_clusters <- function(cl0) {
  if (!is.data.frame(cl0) || !nrow(cl0)) return(NULL)
  cl <- as.data.frame(cl0)
  
  if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
  if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
  if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
  
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
  
  keep <- unique(c(
    "chr","start","end",
    "cluster_id","cluster_chr","cluster_chr_n",
    "center","cluster_size_kb",
    "n_snps","top_snp","top_logp",
    "n_gtex"
  ))
  cl[, intersect(keep, names(cl)), drop = FALSE]
}

# ------------------------------------------------------------
# MAIN INIT
# ------------------------------------------------------------
gi_slave_canonical_init <- function(session) {
  
  safe_update <- function(expr) {
    tryCatch(force(expr), error = function(e) NULL)
    invisible(TRUE)
  }
  
  sid_rv <- shiny::reactiveVal(NULL)
  
  # 1) Preferred: input$gi_qs (Hub/ngrok/caddy)
  shiny::observeEvent(session$input$gi_qs, {
    qs <- as.character(session$input$gi_qs %||% "")
    if (!nzchar(qs)) return()
    q  <- shiny::parseQueryString(sub("^\\?", "", qs))
    s0 <- gi_sanitize_sid(q$sid %||% "")
    if (!is.null(s0)) sid_rv(s0)
  }, ignoreInit = FALSE)
  
  # 2) Fallback: url_search (?sid=...)
  shiny::observe({
    if (!is.null(sid_rv())) return()
    qs <- as.character(session$clientData$url_search %||% "")
    if (!nzchar(qs)) return()
    q  <- shiny::parseQueryString(sub("^\\?", "", qs))
    s0 <- gi_sanitize_sid(q$sid %||% "")
    if (!is.null(s0)) sid_rv(s0)
  })
  
  # 3) Extra fallback: url_hash (#?sid=...)
  shiny::observe({
    if (!is.null(sid_rv())) return()
    hs <- as.character(session$clientData$url_hash %||% "")
    if (!nzchar(hs)) return()
    hs <- sub("^#", "", hs)
    q  <- shiny::parseQueryString(sub("^\\?", "", hs))
    s0 <- gi_sanitize_sid(q$sid %||% "")
    if (!is.null(s0)) sid_rv(s0)
  })
  
  # Debug when sid becomes available
  shiny::observeEvent(sid_rv(), {
    s <- sid_rv()
    if (!sid_ok(s)) {
      cat("[SLAVE] sid not ready yet:", s, "\n")
      return()
    }
    p <- gi_state_paths(s)
    cat("[SLAVE] sid =", s, "\n")
    cat("[SLAVE] json =", p$json, "exists=", file.exists(p$json), "\n")
    cat("[SLAVE] gwas =", p$gwas, "exists=", file.exists(p$gwas), "\n")
    cat("[SLAVE] clus =", p$clus, "exists=", file.exists(p$clus), "\n")
  }, ignoreInit = TRUE)
  
  # Poll JSON mtime
  gi_state <- shiny::reactivePoll(
    intervalMillis = 2000,
    session = session,
    checkFunc = function() {
      s <- sid_rv()
      if (!sid_ok(s)) return(0)
      pj <- gi_state_paths(s)$json
      if (!file.exists(pj)) return(0)
      as.numeric(file.info(pj)$mtime)
    },
    valueFunc = function() {
      s <- sid_rv()
      if (!sid_ok(s)) return(NULL)
      gi_read_state(s)
    }
  )
  
  gi_last_stamp  <- shiny::reactiveVal(NA_integer_)
  gi_state_rv    <- shiny::reactiveVal(NULL)
  gwas_shared_rv <- shiny::reactiveVal(NULL)
  clus_shared_rv <- shiny::reactiveVal(NULL)
  
  shiny::observeEvent(gi_state(), {
    
    s <- sid_rv()
    if (is.null(s)) return()
    
    st <- gi_state()
    if (is.null(st) || !is.list(st)) return()
    
    # Sempre definim st2 abans de fer-lo servir
    st2 <- st
    
    # --- paths
    p <- gi_state_paths(s)
    
    # --- stamp (del JSON; és el trigger de poll)
    stamp <- suppressWarnings(as.integer(st$stamp %||% 0L))
    if (!is.na(gi_last_stamp()) && identical(stamp, gi_last_stamp())) return()
    gi_last_stamp(stamp)
    
    cat("[SLAVE] sync stamp=", stamp, " sid=", s, "\n")
    
    # ---- 0) PARAMS RDS (optional)
    par <- NULL
    if (!is.null(p$params) && file.exists(p$params)) {
      par <- tryCatch(readRDS(p$params), error = function(e) {
        cat("[SLAVE] readRDS(params) ERROR:", conditionMessage(e), "\n")
        NULL
      })
      if (is.list(par)) {
        cat("[SLAVE] PARAMS synced keys=", paste(names(par), collapse=","), "\n")
        st2 <- modifyList(st2, par)  # params tenen prioritat
      }
    }
    
    # Normalize legacy thr_type values (compat)
    if (!is.null(st2$thr_type)) {
      st2$thr_type <- trimws(as.character(st2$thr_type))
      if (identical(st2$thr_type, "logp")) st2$thr_type <- "pthr"
    }
    
    # ✅ ara sí: guardem l'estat “merged” (st2) que farà servir l’app
    gi_state_rv(st2)
    
    # ---- 1) GWAS
    if (file.exists(p$gwas)) {
      
      df0 <- tryCatch(readRDS(p$gwas), error = function(e) {
        cat("[SLAVE] readRDS(gwas) ERROR:", conditionMessage(e), "\n")
        NULL
      })
      
      df_std <- tryCatch(gi_std_gwas(df0), error = function(e) NULL)
      
      df <- NULL
      if (is.data.frame(df_std) && nrow(df_std) > 0) {
        df <- df_std
      } else if (is.data.frame(df0) && nrow(df0) > 0) {
        df <- df0
      } else if (is.list(df0)) {
        cand <- df0[["gwas"]] %||% df0[["data"]] %||% df0[[1]]
        if (is.data.frame(cand) && nrow(cand) > 0) df <- cand
      }
      
      if (is.data.frame(df) && nrow(df) > 0) {
        gwas_shared_rv(df)
        cat("[SLAVE] GWAS synced rows=", nrow(df),
            " cols=", paste(names(df), collapse=","), "\n")
      } else {
        cat("[SLAVE] GWAS file read but empty/invalid | class(df0)=",
            paste(class(df0), collapse=","), "\n")
      }
      
    } else {
      cat("[SLAVE] GWAS file not found:", p$gwas, "\n")
    }
    
    # ---- helpers: update slider OR numeric, depèn de cada app
    update_num_like <- function(id, value) {
      if (!(id %in% names(session$input))) return(invisible(FALSE))
      v <- suppressWarnings(as.numeric(value))
      if (!is.finite(v)) return(invisible(FALSE))
      
      x <- session$input[[id]]
      if (is.numeric(x) && length(x) == 1) {
        # pot ser slider o numeric; provem slider primer
        ok <- tryCatch({ shiny::updateSliderInput(session, id, value = v); TRUE },
                       error = function(e) FALSE)
        if (!ok) {
          tryCatch({ shiny::updateNumericInput(session, id, value = v); TRUE },
                   error = function(e) FALSE)
        }
      } else {
        tryCatch({ shiny::updateNumericInput(session, id, value = v); TRUE },
                 error = function(e) FALSE)
      }
      invisible(TRUE)
    }
    
    do_ui_updates <- function() {
      
      safe_update(if ("cluster_method" %in% names(session$input) && !is.null(st2$cluster_method))
        shiny::updateRadioButtons(session, "cluster_method", selected = as.character(st2$cluster_method)))
      
      safe_update(if ("hits_mode" %in% names(session$input) && !is.null(st2$hits_mode))
        shiny::updateRadioButtons(session, "hits_mode", selected = as.character(st2$hits_mode)))
      
      # Threshold sync (thr_type = "pthr" | "min_logp")
      safe_update({
        tt <- trimws(as.character(st2$thr_type %||% ""))
        tv <- st2$thr_value %||% NA_real_
        
        if (identical(tt, "pthr")) {
          update_num_like("pthr", tv)
        } else if (identical(tt, "min_logp")) {
          update_num_like("min_logp", tv)
        }
      })
      
      safe_update(if ("flank" %in% names(session$input) && is.finite(as.numeric(st2$flank_bp %||% NA)))
        shiny::updateNumericInput(session, "flank", value = as.integer(st2$flank_bp)))
      
      safe_update(if ("win_bp" %in% names(session$input) && is.finite(as.numeric(st2$win_bp %||% NA)))
        shiny::updateNumericInput(session, "win_bp", value = as.integer(st2$win_bp)))
      
      safe_update(if ("step_bp" %in% names(session$input) && is.finite(as.numeric(st2$step_bp %||% NA)))
        shiny::updateNumericInput(session, "step_bp", value = as.integer(st2$step_bp)))
      
      safe_update(if ("min_hits" %in% names(session$input) && is.finite(as.numeric(st2$min_hits %||% NA)))
        shiny::updateNumericInput(session, "min_hits", value = as.integer(st2$min_hits)))
    }
    
    do_ui_updates()
    
    # ---- 3) CLUSTERS
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
  
  
  list(
    sid            = sid_rv,
    gi_state       = gi_state,
    state_shared   = gi_state_rv,
    gwas_shared    = gwas_shared_rv,
    clusters_shared= clus_shared_rv
  )
}
