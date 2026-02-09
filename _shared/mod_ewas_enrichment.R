# R/mod_ewas_enrichment.R  (MAX-NA-SAFE / SINGLE TABLE+PLOT / AUTO-RECALC + INFO MODALS)
# =============================================================================
# Weights modes kept ONLY:
#   - |delta|
#   - -log10(FDR)
#   - |delta| * -log10(FDR)
# (plus Unweighted as fallback)
# =============================================================================
# /Volumes/DISK1TB/Inspector_app_slaves/EWAS_disease/R/mod_ewas_enrichment.R

library(shiny)
library(DT)
library(plotly)
library(data.table)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# NA-safe string + numeric helpers
# -----------------------------------------------------------------------------

safe_chr <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- gsub("\xC2\xA0", " ", x, useBytes = TRUE)
  x <- enc2utf8(x)
  x
}

clean_item <- function(x) {
  x <- safe_chr(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  up <- toupper(x)
  x[up %in% c("NA", "N/A", "NULL", "NONE", ".", "-", "?", "NAN")] <- NA_character_
  x[x == ""] <- NA_character_
  x
}

safe_num <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  xx <- safe_chr(x)
  xx <- trimws(xx)
  xx[xx == ""] <- NA_character_
  xx <- gsub(",", ".", xx, fixed = TRUE)
  xx <- gsub("[^0-9eE+\\-\\.]", "", xx)
  suppressWarnings(as.numeric(xx))
}

safe_int <- function(x) suppressWarnings(as.integer(round(safe_num(x))))

norm_chr <- function(x) {
  x <- clean_item(x)
  x[!is.na(x) & !startsWith(x, "chr")] <- paste0("chr", x[!is.na(x) & !startsWith(x, "chr")])
  x
}

# -----------------------------------------------------------------------------
# Robust event table builder (unit = unique (bin_id, item))
# -----------------------------------------------------------------------------

pick_first_existing <- function(nms, candidates) {
  candidates <- candidates[candidates %in% nms]
  if (length(candidates)) candidates[1] else NULL
}

as_events_tbl <- function(dt, item_candidates = c("disease","cancer")) {
  validate(need(is.data.frame(dt) || is.data.table(dt), "events table is not a data.frame/data.table"))
  dt <- as.data.table(dt)
  
  item_col <- pick_first_existing(names(dt), item_candidates)
  validate(need(!is.null(item_col), "Observed/ref table must have 'disease' or 'cancer' column."))
  
  req_bin <- c("chr","bin_start","bin_end")
  miss_bin <- setdiff(req_bin, names(dt))
  validate(need(!length(miss_bin), paste("Missing bin columns:", paste(miss_bin, collapse = ", "))))
  validate(need("padj" %in% names(dt), "Missing column: padj"))
  validate(need("delta" %in% names(dt), "Missing column: delta"))
  
  dt[, chr := norm_chr(get("chr"))]
  dt[, bin_start := safe_int(get("bin_start"))]
  dt[, bin_end   := safe_int(get("bin_end"))]
  dt[, padj      := safe_num(get("padj"))]
  dt[, delta     := safe_num(get("delta"))]
  dt[, item      := clean_item(get(item_col))]
  
  dt[, bin_id := ifelse(!is.na(chr) & is.finite(bin_start) & is.finite(bin_end),
                        paste0(chr, ":", bin_start, "-", bin_end),
                        NA_character_)]
  
  dt <- dt[!is.na(bin_id) & nzchar(bin_id) & !is.na(item) & nzchar(item)]
  
  keep <- c("item","chr","bin_start","bin_end","bin_id","padj","delta")
  if ("cluster_id" %in% names(dt)) {
    dt[, cluster_id := clean_item(cluster_id)]
    keep <- c(keep, "cluster_id")
  }
  dt <- dt[, ..keep]
  
  setkey(dt, bin_id, item)
  unique(dt)
}

apply_direction_filter <- function(dt, direction = c("all","hyper","hypo")) {
  direction <- match.arg(direction)
  if (identical(direction, "all")) return(dt)
  dt <- dt[is.finite(delta)]
  if (identical(direction, "hyper")) dt <- dt[delta > 0]
  if (identical(direction, "hypo"))  dt <- dt[delta < 0]
  dt
}

# -----------------------------------------------------------------------------
# Weights (ONLY allowed modes)
# -----------------------------------------------------------------------------

make_weights <- function(padj, delta,
                         mode = c("none","neglog10","abs_delta","absdelta_neglog10"),
                         w_cap = 1e6,
                         fdr_floor = 1e-300) {
  
  mode <- match.arg(mode)
  
  if (identical(mode, "none")) return(rep(1, length(padj)))
  
  f <- safe_num(padj)
  f[!is.finite(f) | f <= 0] <- NA_real_
  f <- pmax(f, fdr_floor)
  
  d <- abs(safe_num(delta))
  d[!is.finite(d)] <- NA_real_
  
  w <- switch(
    mode,
    "neglog10"          = -log10(f),
    "abs_delta"         = d,
    "absdelta_neglog10" = d * (-log10(f))
  )
  
  w[!is.finite(w)] <- NA_real_
  w <- pmin(w, w_cap)
  w[!is.finite(w) | w <= 0] <- NA_real_
  w
}

cap_by_quantile <- function(w, q = 0.995) {
  w2 <- w
  ok <- is.finite(w2)
  if (!any(ok)) return(w2)
  qq <- suppressWarnings(as.numeric(q))
  if (!is.finite(qq) || qq <= 0 || qq >= 1) return(w2)
  cap <- as.numeric(stats::quantile(w2[ok], probs = qq, na.rm = TRUE, names = FALSE))
  if (is.finite(cap)) w2[ok] <- pmin(w2[ok], cap)
  w2
}

# -----------------------------------------------------------------------------
# Formatting for DT (safe)
# -----------------------------------------------------------------------------

format_num_3_or_sci <- function(x) {
  x <- safe_num(x)
  out <- rep("", length(x))
  ok <- is.finite(x)
  small <- ok & (abs(x) < 1e-3) & (x != 0)
  norm  <- ok & !small
  out[small] <- format(x[small], scientific = TRUE, digits = 3)
  out[norm]  <- format(round(x[norm], 3), nsmall = 3, scientific = FALSE, trim = TRUE)
  out[!ok] <- ""
  out
}

format_int_plain <- function(x) {
  x <- safe_num(x)
  out <- rep("", length(x))
  ok <- is.finite(x)
  out[ok] <- as.character(as.integer(round(x[ok])))
  out[!ok] <- ""
  out
}

# -----------------------------------------------------------------------------
# Core enrichment (NA-safe)
# -----------------------------------------------------------------------------

compute_ewas_enrichment <- function(obs_events,
                                    ref_events,
                                    scope = c("global","cluster"),
                                    cluster_id = NULL,
                                    direction = c("all","hyper","hypo"),
                                    weight_mode = c("none","neglog10","abs_delta","absdelta_neglog10"),
                                    winsorize = TRUE,
                                    wins_q = 0.995,
                                    bin_normalize = FALSE,
                                    w_cap = 1e6,
                                    fdr_floor = 1e-300) {
  
  scope <- match.arg(scope)
  direction <- match.arg(direction)
  weight_mode <- match.arg(weight_mode)
  
  obs <- as_events_tbl(obs_events)
  ref <- as_events_tbl(ref_events)
  
  validate(need(nrow(obs) > 0, "Observed events empty after cleaning (NA/empty removed)."))
  validate(need(nrow(ref) > 0, "Reference events empty after cleaning (NA/empty removed)."))
  
  if (identical(scope, "cluster")) {
    cid <- clean_item(cluster_id %||% "")
    validate(need(length(cid) == 1 && nzchar(cid) && !is.na(cid), "Select a cluster."))
    validate(need("cluster_id" %in% names(obs), "Observed events missing 'cluster_id'."))
    obs <- obs[cluster_id == cid]
    validate(need(nrow(obs) > 0, "No observed events in selected cluster (after cleaning)."))
  }
  
  obs <- apply_direction_filter(obs, direction)
  ref <- apply_direction_filter(ref, direction)
  validate(need(nrow(obs) > 0, "No observed events after direction filter."))
  validate(need(nrow(ref) > 0, "No reference events after direction filter."))
  
  ref_items <- unique(ref$item)
  obs <- obs[item %in% ref_items]
  validate(need(nrow(obs) > 0, "No observed items match reference items (after cleaning)."))
  
  # ---- COUNTS ----
  N_ref <- nrow(ref)
  ref_cnt <- ref[, .(n_ref = .N), by = item]
  obs_cnt <- obs[, .(obs_bins = .N), by = item]
  N_obs <- sum(obs_cnt$obs_bins)
  
  cnt <- merge(obs_cnt, ref_cnt, by = "item", all.x = TRUE, all.y = FALSE)
  cnt <- cnt[!is.na(n_ref) & is.finite(n_ref) & n_ref > 0]
  cnt[, p_ref := n_ref / N_ref]
  cnt[, expected := N_obs * p_ref]
  cnt[, enrich_ratio := obs_bins / pmax(expected, 1e-12)]
  cnt[, pval := phyper(obs_bins - 1, m = n_ref, n = N_ref - n_ref, k = N_obs, lower.tail = FALSE)]
  cnt[, fdr := p.adjust(pval, method = "BH")]
  setorder(cnt, fdr, -enrich_ratio)
  
  # ---- WEIGHTS (ONLY selected modes) ----
  ref[, w := make_weights(padj, delta, weight_mode, w_cap = w_cap, fdr_floor = fdr_floor)]
  obs[, w := make_weights(padj, delta, weight_mode, w_cap = w_cap, fdr_floor = fdr_floor)]
  ref <- ref[is.finite(w)]
  obs <- obs[is.finite(w)]
  validate(need(nrow(ref) > 0 && nrow(obs) > 0, "No finite weights available after cleaning."))
  
  if (!identical(weight_mode, "none") && isTRUE(winsorize)) {
    ref[, w := cap_by_quantile(w, q = wins_q)]
    obs[, w := cap_by_quantile(w, q = wins_q)]
  }
  
  if (!identical(weight_mode, "none") && isTRUE(bin_normalize)) {
    ref[, w := w / sum(w, na.rm = TRUE), by = bin_id]
    obs[, w := w / sum(w, na.rm = TRUE), by = bin_id]
    ref <- ref[is.finite(w)]
    obs <- obs[is.finite(w)]
    validate(need(nrow(ref) > 0 && nrow(obs) > 0, "No finite weights after per-bin normalization."))
  }
  
  ref_w <- ref[, .(w_ref = sum(w, na.rm = TRUE), n_ref_w = .N), by = item]
  tot_w_ref <- sum(ref_w$w_ref, na.rm = TRUE)
  ref_w <- ref_w[is.finite(w_ref) & w_ref > 0]
  ref_w[, p_w_ref := w_ref / pmax(tot_w_ref, 1e-12)]
  
  obs_w <- obs[, .(obs_w = sum(w, na.rm = TRUE), obs_bins_w = .N), by = item]
  tot_w_obs <- sum(obs_w$obs_w, na.rm = TRUE)
  
  wtbl <- merge(obs_w, ref_w[, .(item, w_ref, p_w_ref, n_ref_w)], by = "item", all.x = TRUE)
  wtbl <- wtbl[!is.na(p_w_ref) & is.finite(p_w_ref) & p_w_ref > 0]
  wtbl[, expected_w := p_w_ref * tot_w_obs]
  wtbl[, enrich_ratio_w := obs_w / pmax(expected_w, 1e-12)]
  
  tot_w2_obs <- sum(obs$w^2, na.rm = TRUE)
  wtbl[, var_w := p_w_ref * (1 - p_w_ref) * pmax(tot_w2_obs, 1e-12)]
  wtbl[, z_w := (obs_w - expected_w) / sqrt(var_w)]
  wtbl[, pval_w := pnorm(z_w, lower.tail = FALSE)]
  wtbl[, fdr_w := p.adjust(pval_w, method = "BH")]
  setorder(wtbl, fdr_w, -enrich_ratio_w)
  
  list(counts = cnt[], weights = wtbl[])
}

# -----------------------------------------------------------------------------
# Module UI (single table + plot) + INFO buttons
# -----------------------------------------------------------------------------

# UI: add grey background + padding + rounded corners to the LEFT sidebar column
# (wrap the controls in a div with inline style OR a CSS class)

mod_ewas_enrich_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      
      # =========================
      # LEFT SIDEBAR (GREY BG)
      # =========================
      column(
        width = 3,
        
        div(
          class = "panelGrey",
          style = "background:#F3F4F6; border:1px solid #E5E7EB; border-radius:10px; padding:12px;",
          
          selectInput(
            ns("view"), "View",
            choices = c("Counts" = "counts", "Weights" = "weights"),
            selected = "counts"
          ),
          
          tags$hr(style = "margin:10px 0;"),
          
          radioButtons(
            ns("scope"), "Scope",
            choices = c("Global" = "global", "Cluster" = "cluster"),
            selected = "global"
          ),
          uiOutput(ns("cluster_ui")),
          
          selectInput(
            ns("direction"), "Direction",
            choices = c("All" = "all", "Hyper (delta>0)" = "hyper", "Hypo (delta<0)" = "hypo"),
            selected = "all"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'weights'", ns("view")),
            selectInput(
              ns("weight_mode"), "Weights",
              choices = c(
                "|delta|"               = "abs_delta",
                "-log10(FDR)"           = "neglog10",
                "|delta| * -log10(FDR)" = "absdelta_neglog10",
                "Unweighted (w=1)"      = "none"
              ),
              selected = "abs_delta"
            ),
            checkboxInput(ns("winsorize"), "Winsorize weights (recommended)", value = TRUE),
            numericInput(ns("wins_q"), "Winsorize quantile", value = 0.995, min = 0.90, max = 0.9999, step = 0.001),
            checkboxInput(ns("bin_normalize"), "Per-bin normalization (each bin sums to 1)", value = FALSE),
            numericInput(ns("w_cap"), "Weight cap (hard)", value = 1e6, min = 1, step = 1),
            numericInput(ns("fdr_floor"), "FDR floor", value = 1e-10, min = 0, step = 0.0)
          )
        )
      ),
      
      # =========================
      # RIGHT MAIN AREA
      # =========================
      column(
        width = 9,
        
        # --- PANEL 1: PLOT ---
        div(
          class = "panelGrey",
          div(class = "panelTitle",
              HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Enrichment plot</span>")
          ),
          plotlyOutput(ns("plot"), height = "440px")
        ),
        
        tags$div(style = "height:10px;"),
        
        # --- PANEL 2: TABLE ---
        div(
          class = "panelGrey",
          div(class = "panelTitle",
              HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📋 Enrichment table</span>")
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'counts'", ns("view")),
            actionLink(ns("info_enrich_counts"), "ℹ️ info")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'weights'", ns("view")),
            actionLink(ns("info_enrich_weights"), "ℹ️ info")
          ),
          
          DTOutput(ns("tbl"))
        )
      )
    )
  )
}


# -----------------------------------------------------------------------------
# Module server + INFO MODALS
# -----------------------------------------------------------------------------

mod_ewas_enrich_server <- function(id,
                                   obs_events,          # reactive observed detail table
                                   app_mode,            # reactive: "disease" or "tumor"
                                   ref_paths,           # list(disease="...rds", tumor="...rds")
                                   top_n = 20) {
  
  moduleServer(id, function(input, output, session) {
    
    # ----------------------------
    # INFO MODALS (English)
    # ----------------------------
    
    observeEvent(input$info_enrich_weights, {
      txt <- HTML("
      <h4 style='margin-top:0;'>Weighted enrichment table — column definitions</h4>
      <ul>
        <li><b>item</b>: Disease/cancer name being tested.</li>
        <li><b>obs_bins_w</b>: Number of observed <i>events</i> for this item (integer). An event is one unique (bin_id, item) pair in the observed dataset.</li>
        <li><b>obs_w</b>: Total observed <i>weight</i> for this item (sum of event weights <code>w</code> in the observed dataset, after winsorization/per-bin normalization if enabled).</li>
        <li><b>expected_w</b>: Expected total weight under the genome reference: <code>expected_w = p_w_ref × total_obs_w</code>.</li>
        <li><b>enrich_ratio_w</b>: Weighted enrichment ratio: <code>obs_w / expected_w</code>. Values &gt; 1 indicate enrichment.</li>
        <li><b>z_w</b>: Z-score (normal approximation): <code>(obs_w − expected_w) / sqrt(var_w)</code>, where <code>var_w ≈ p_w_ref × (1 − p_w_ref) × Σ(w²)</code> over observed events.</li>
        <li><b>n_ref_w</b>: Number of reference events for this item (integer). Reference is genome-wide (bin_id, item) pairs.</li>
        <li><b>pval_w</b>: One-sided p-value for enrichment (upper tail): <code>P(Z ≥ z_w)</code>.</li>
        <li><b>fdr_w</b>: Benjamini–Hochberg FDR across all items.</li>
      </ul>
      ")
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        footer = modalButton("Close"),
        size = "m",
        tags$div(style = "line-height:1.35;", txt)
      ))
    })
    
    observeEvent(input$info_enrich_counts, {
      txt <- HTML("
      <h4 style='margin-top:0;'>Count enrichment table — column definitions</h4>
      <ul>
        <li><b>item</b>: Disease/cancer name being tested.</li>
        <li><b>obs_bins</b>: Number of observed <i>events</i> for this item (integer). An event is one unique (bin_id, item) pair in the observed dataset.</li>
        <li><b>expected</b>: Expected number of observed events under the reference: <code>expected = N_obs × p_ref</code>.</li>
        <li><b>enrich_ratio</b>: Count enrichment ratio: <code>obs_bins / expected</code>. Values &gt; 1 indicate enrichment.</li>
        <li><b>n_ref</b>: Number of reference events for this item (integer).</li>
        <li><b>p_ref</b>: Reference probability: <code>p_ref = n_ref / N_ref</code>.</li>
        <li><b>pval</b>: One-sided hypergeometric p-value (upper tail) for enrichment.</li>
        <li><b>fdr</b>: Benjamini–Hochberg FDR across all items.</li>
      </ul>
      ")
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        footer = modalButton("Close"),
        size = "m",
        tags$div(style = "line-height:1.35;", txt)
      ))
    })
    
    # ----------------------------
    # Reference cache
    # ----------------------------
    
    ref_cache <- reactiveVal(NULL)
    
    observeEvent(app_mode(), {
      mode <- clean_item(app_mode() %||% "disease")
      validate(need(!is.na(mode) && nzchar(mode), "app_mode() is empty/NA."))
      path <- ref_paths[[mode]]
      validate(need(!is.null(path) && nzchar(path), "Reference path not set for this app_mode."))
      validate(need(file.exists(path), paste("Reference RDS not found:", path)))
      ref_cache(readRDS(path))
    }, ignoreInit = FALSE)
    
    # ----------------------------
    # Cluster selector UI
    # ----------------------------
    
    output$cluster_ui <- renderUI({
      req(input$scope)
      if (!identical(input$scope, "cluster")) return(NULL)
      
      dt <- obs_events()
      validate(need(is.data.frame(dt) && nrow(dt) > 0, "No observed EWAS events available."))
      validate(need("cluster_id" %in% names(dt), "Observed table has no 'cluster_id'."))
      
      ids <- clean_item(unique(dt$cluster_id))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      ids <- sort(unique(ids))
      validate(need(length(ids) > 0, "No valid cluster_id values found (after cleaning)."))
      
      cur <- clean_item(isolate(input$cluster_id %||% ""))
      sel <- if (!is.na(cur) && nzchar(cur) && cur %in% ids) cur else ids[1]
      
      selectInput(session$ns("cluster_id"), "Cluster", choices = ids, selected = sel)
    })
    
    # ----------------------------
    # Auto recompute
    # ----------------------------
    
    res <- reactive({
      ref_dt <- ref_cache()
      validate(need(is.data.frame(ref_dt) && nrow(ref_dt) > 0, "Reference not loaded or empty."))
      
      obs_dt <- obs_events()
      validate(need(is.data.frame(obs_dt) && nrow(obs_dt) > 0, "Observed table is empty."))
      
      scope <- input$scope %||% "global"
      cid   <- if (identical(scope, "cluster")) clean_item(input$cluster_id %||% "") else NULL
      dir   <- input$direction %||% "all"
      
      wm    <- input$weight_mode %||% "absdelta_neglog10"
      wins  <- isTRUE(input$winsorize)
      qv    <- suppressWarnings(as.numeric(input$wins_q %||% 0.995))
      bn    <- isTRUE(input$bin_normalize)
      wc    <- suppressWarnings(as.numeric(input$w_cap %||% 1e6))
      ff    <- suppressWarnings(as.numeric(input$fdr_floor %||% 1e-300))
      
      compute_ewas_enrichment(
        obs_events    = obs_dt,
        ref_events    = ref_dt,
        scope         = scope,
        cluster_id    = cid,
        direction     = dir,
        weight_mode   = wm,
        winsorize     = wins,
        wins_q        = qv,
        bin_normalize = bn,
        w_cap         = wc,
        fdr_floor     = ff
      )
    })
    
    # ----------------------------
    # Current view
    # ----------------------------
    
    cur_tbl <- reactive({
      r <- res()
      view <- input$view %||% "counts"
      
      if (identical(view, "weights")) {
        df <- as.data.frame(r$weights)
        validate(need(nrow(df) > 0, "No weighted enrichment results."))
        
        show <- c("item","obs_bins_w","obs_w","expected_w","enrich_ratio_w","z_w","n_ref_w","pval_w","fdr_w")
        show <- show[show %in% names(df)]
        df <- df[, show, drop = FALSE]
        
        disp <- df
        if ("obs_bins_w" %in% names(disp)) disp$obs_bins_w <- format_int_plain(disp$obs_bins_w)
        if ("n_ref_w" %in% names(disp))    disp$n_ref_w    <- format_int_plain(disp$n_ref_w)
        
        num_cols <- setdiff(names(disp), c("item","obs_bins_w","n_ref_w"))
        for (cc in num_cols) disp[[cc]] <- format_num_3_or_sci(df[[cc]])
        
        return(list(type = "weights", raw = df, disp = disp))
      }
      
      df <- as.data.frame(r$counts)
      validate(need(nrow(df) > 0, "No count enrichment results."))
      
      show <- c("item","obs_bins","expected","enrich_ratio","n_ref","p_ref","pval","fdr")
      show <- show[show %in% names(df)]
      df <- df[, show, drop = FALSE]
      
      disp <- df
      if ("n_ref" %in% names(disp))    disp$n_ref    <- format_int_plain(disp$n_ref)
      if ("obs_bins" %in% names(disp)) disp$obs_bins <- format_int_plain(disp$obs_bins)
      
      num_cols <- setdiff(names(disp), c("item","n_ref","obs_bins"))
      for (cc in num_cols) disp[[cc]] <- format_num_3_or_sci(df[[cc]])
      
      list(type = "counts", raw = df, disp = disp)
    })
    
    output$tbl <- renderDT({
      x <- cur_tbl()
      
      DT::datatable(
        x$disp,
        rownames   = FALSE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons    = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 15,
          scrollX    = TRUE
        )
      )
    })
    
    output$plot <- renderPlotly({
      x <- cur_tbl()
      df <- x$raw
      validate(need(is.data.frame(df) && nrow(df) > 0, ""))
      
      if (identical(x$type, "weights")) {
        df <- df[order(df$fdr_w, -df$enrich_ratio_w), , drop = FALSE]
        df <- head(df, top_n)
        df$item <- factor(df$item, levels = rev(df$item))
        
        return(
          plot_ly(
            data = df,
            x = ~enrich_ratio_w,
            y = ~item,
            type = "bar",
            marker = list(
              color = ~fdr_w,
              colorscale = "Viridis",
              reversescale = TRUE,
              colorbar = list(title = "FDR (w)")
            ),
            customdata = ~cbind(obs_bins_w, obs_w, expected_w, fdr_w),
            hovertemplate = paste(
              "<b>%{y}</b><br>",
              "Enrich ratio (w): %{x:.3f}<br>",
              "Obs bins: %{customdata[0]}<br>",
              "Obs weight: %{customdata[1]:.3g}<br>",
              "Expected weight: %{customdata[2]:.3g}<br>",
              "FDR(w): %{customdata[3]:.3g}<extra></extra>"
            )
          ) %>%
            layout(
              xaxis = list(title = "Observed / Expected (weights)"),
              yaxis = list(title = "", automargin = TRUE),
              margin = list(l = 210, r = 50, t = 10, b = 40)
            )
        )
      }
      
      df <- df[order(df$fdr, -df$enrich_ratio), , drop = FALSE]
      df <- head(df, top_n)
      df$item <- factor(df$item, levels = rev(df$item))
      
      plot_ly(
        data = df,
        x = ~enrich_ratio,
        y = ~item,
        type = "bar",
        marker = list(
          color = ~fdr,
          colorscale = "Viridis",
          reversescale = TRUE,
          colorbar = list(title = "FDR")
        ),
        customdata = ~cbind(obs_bins, expected, fdr),
        hovertemplate = paste(
          "<b>%{y}</b><br>",
          "Enrich ratio: %{x:.3f}<br>",
          "Obs bins: %{customdata[0]}<br>",
          "Expected: %{customdata[1]:.3f}<br>",
          "FDR: %{customdata[2]:.3g}<extra></extra>"
        )
      ) %>%
        layout(
          xaxis = list(title = "Observed / Expected (counts)"),
          yaxis = list(title = "", automargin = TRUE),
          margin = list(l = 210, r = 50, t = 10, b = 40)
        )
    })
    
    invisible(NULL)
  })
}
