# _shared/gi_clusters_canonical.R
# ------------------------------------------------------------
# GItools CANONICAL CLUSTER ENGINE (common to ALL apps)
# >>>>>>>>>> Generate clusters from LOCAL app demand
#
# Responsibilities:
#  - Given a GWAS dataframe (CHR,BP,Pval,logp,snp[,rsid]), compute:
#      * intervals_raw(): raw intervals (per-hit) and/or per-cluster intervals
#      * clusters_cur(): canonical cluster table
#      * selected_cluster(): selected row from DT clusters table
#  - Supports two clustering methods:
#      1) window: +/- flank around hits with logp >= pthr, then merge overlaps
#      2) hits  : density-based windows (span1mb / tiled / sliding)
#
# NOTE:
#  - This file does NOT know anything about GTEx/NonSyn/EWAS specifics.
#  - Each app can:
#      * use build_ranges (local clustering) via gi_clusters_canonical_init()
#      * OR override clusters_cur()/intervals_raw() with clusters from the Master (gi_slave_canonical.R)
# ------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# -----------------------------
# Canonical clustering helpers
# -----------------------------

merge_overlapping_intervals <- function(df_int) {
  if (is.null(df_int) || !is.data.frame(df_int) || nrow(df_int) == 0) {
    return(tibble::tibble(chr=integer(), start=integer(), end=integer()))
  }
  
  x <- as.data.frame(df_int)
  
  # accept either (start,end) OR (bp1,bp2)
  if (!all(c("start","end") %in% names(x))) {
    if (all(c("bp1","bp2") %in% names(x))) {
      x$start <- pmin(as.integer(x$bp1), as.integer(x$bp2))
      x$end   <- pmax(as.integer(x$bp1), as.integer(x$bp2))
    } else {
      stop("merge_overlapping_intervals(): need columns (chr,start,end) or (chr,bp1,bp2).")
    }
  }
  
  x <- x %>%
    dplyr::transmute(
      chr   = suppressWarnings(as.integer(.data$chr)),
      start = suppressWarnings(as.integer(.data$start)),
      end   = suppressWarnings(as.integer(.data$end))
    ) %>%
    dplyr::filter(is.finite(.data$chr), is.finite(.data$start), is.finite(.data$end), .data$end >= .data$start) %>%
    dplyr::arrange(.data$chr, .data$start, .data$end)
  
  if (nrow(x) == 0) return(tibble::tibble(chr=integer(), start=integer(), end=integer()))
  
  out <- vector("list", length = 0)
  
  for (cc in unique(x$chr)) {
    xi <- x[x$chr == cc, , drop = FALSE]
    if (nrow(xi) == 0) next
    
    cur_s <- xi$start[1]
    cur_e <- xi$end[1]
    
    if (nrow(xi) > 1) {
      for (i in 2:nrow(xi)) {
        if (xi$start[i] <= cur_e) {
          cur_e <- max(cur_e, xi$end[i])
        } else {
          out[[length(out)+1]] <- tibble::tibble(chr=cc, start=cur_s, end=cur_e)
          cur_s <- xi$start[i]
          cur_e <- xi$end[i]
        }
      }
    }
    out[[length(out)+1]] <- tibble::tibble(chr=cc, start=cur_s, end=cur_e)
  }
  
  dplyr::bind_rows(out) %>% dplyr::arrange(.data$chr, .data$start, .data$end)
}

summarize_clusters_from_hits <- function(merged_intervals, gwas_df) {
  if (is.null(merged_intervals) || !is.data.frame(merged_intervals) || nrow(merged_intervals) == 0) {
    return(tibble::tibble())
  }
  df <- gwas_df
  req_cols <- c("CHR","BP","snp","logp")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop(paste0("summarize_clusters_from_hits(): GWAS missing: ", paste(miss, collapse=", ")))
  
  cl <- merged_intervals %>%
    dplyr::transmute(
      chr   = suppressWarnings(as.integer(.data$chr)),
      start = suppressWarnings(as.integer(.data$start)),
      end   = suppressWarnings(as.integer(.data$end))
    ) %>%
    dplyr::filter(is.finite(.data$chr), is.finite(.data$start), is.finite(.data$end), .data$end >= .data$start) %>%
    dplyr::arrange(.data$chr, .data$start, .data$end) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(cluster_chr = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cluster_chr_n = paste0("chr", chr_label_plink(.data$chr), "_", .data$cluster_chr),
      cluster_id    = .data$cluster_chr_n
    ) %>%
    dplyr::mutate(cluster = dplyr::row_number())
  
  # compute stats per interval
  cl_sum <- purrr::pmap_dfr(
    list(cl$chr, cl$start, cl$end, cl$cluster_id, cl$cluster, cl$cluster_chr, cl$cluster_chr_n),
    function(chr, st, en, cid, cnum, cchr, cchrn) {
      sub <- df %>% dplyr::filter(as.integer(.data$CHR) == as.integer(chr), .data$BP >= st, .data$BP <= en)
      if (!nrow(sub)) {
        tibble::tibble(
          chr = chr, start = st, end = en,
          center = as.integer(round((st + en)/2)),
          n_snps = 0L,
          top_snp = NA_character_,
          top_logp = NA_real_,
          cluster_size_kb = round((en - st) / 1000, 2),
          cluster = cnum, cluster_chr = cchr, cluster_chr_n = cchrn, cluster_id = cid
        )
      } else {
        ii <- which.max(sub$logp)
        tibble::tibble(
          chr = chr, start = st, end = en,
          center = as.integer(round(mean(sub$BP, na.rm = TRUE))),
          n_snps = nrow(sub),
          top_snp = as.character(sub$snp[ii]),
          top_logp = max(sub$logp, na.rm = TRUE),
          cluster_size_kb = round((en - st) / 1000, 2),
          cluster = cnum, cluster_chr = cchr, cluster_chr_n = cchrn, cluster_id = cid
        )
      }
    }
  )
  
  cl_sum %>% dplyr::arrange(.data$chr, .data$start, .data$end)
}

# -----------------------------
# Helpers clustering (hit-density windows)
# -----------------------------

merge_significant_windows <- function(win_df, gap_bp = 0L) {
  gap_bp <- suppressWarnings(as.integer(gap_bp)); if (!is.finite(gap_bp)) gap_bp <- 0L
  dt <- data.table::as.data.table(win_df)
  stopifnot(all(c("chr","start","end") %in% names(dt)))
  
  dt[, chr := as.integer(chr)]
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  dt <- dt[is.finite(chr) & is.finite(start) & is.finite(end)]
  if (nrow(dt) == 0) return(as.data.frame(dt))
  
  data.table::setorder(dt, chr, start, end)
  
  merged <- dt[, {
    if (.N == 1L) {
      data.table::data.table(start = start[1], end = end[1])
    } else {
      out_s <- integer(); out_e <- integer()
      cur_s <- start[1]; cur_e <- end[1]
      for (i in 2L:.N) {
        s <- start[i]; e <- end[i]
        if (!is.finite(s) || !is.finite(e) || !is.finite(cur_e)) next
        if (s <= (cur_e + gap_bp + 1L)) cur_e <- max(cur_e, e, na.rm = TRUE)
        else { out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e); cur_s <- s; cur_e <- e }
      }
      out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e)
      data.table::data.table(start = out_s, end = out_e)
    }
  }, by = chr]
  
  merged[, center := as.integer(round((start + end) / 2))]
  merged[, cluster_size_kb := round((end - start) / 1000, 2)]
  as.data.frame(merged)
}

make_windows_dt <- function(chr_lengths, win_bp, step_bp) {
  win_bp  <- suppressWarnings(as.integer(win_bp))
  step_bp <- suppressWarnings(as.integer(step_bp))
  stopifnot(is.finite(win_bp), win_bp > 0L, is.finite(step_bp), step_bp > 0L)
  
  cl <- data.table::as.data.table(chr_lengths)
  stopifnot(all(c("chr","chr_len") %in% names(cl)))
  cl[, chr := as.integer(chr)]
  cl[, chr_len := as.integer(chr_len)]
  cl <- cl[is.finite(chr) & is.finite(chr_len) & chr_len > 0L]
  
  wins <- cl[, {
    starts <- seq.int(1L, chr_len, by = step_bp)
    ends   <- pmin.int(starts + win_bp - 1L, chr_len)
    data.table::data.table(start = starts, end = ends)
  }, by = chr]
  
  wins[]
}

count_hits_in_windows <- function(hits, windows) {
  h <- data.table::as.data.table(hits)
  w <- data.table::as.data.table(windows)
  stopifnot(all(c("chr","pos") %in% names(h)))
  stopifnot(all(c("chr","start","end") %in% names(w)))
  
  h[, chr := as.integer(chr)]
  h[, pos := as.integer(pos)]
  h <- h[is.finite(chr) & is.finite(pos)]
  
  w[, chr := as.integer(chr)]
  w[, start := as.integer(start)]
  w[, end := as.integer(end)]
  w <- w[is.finite(chr) & is.finite(start) & is.finite(end)]
  data.table::setorder(w, chr, start, end)
  
  # point intervals
  h[, `:=`(start = pos, end = pos)]
  
  data.table::setkey(w, chr, start, end)
  data.table::setkey(h, chr, start, end)
  
  ov <- data.table::foverlaps(h, w, nomatch = 0L)
  if (nrow(ov) == 0) {
    w[, n_hits := 0L]
    return(as.data.frame(w))
  }
  
  cnt <- ov[, .(n_hits = .N), by = .(chr, start, end)]
  w <- merge(w, cnt, by = c("chr","start","end"), all.x = TRUE)
  w[is.na(n_hits), n_hits := 0L]
  as.data.frame(w[])
}

add_cluster_stats_from_hits <- function(clusters_df, df_sig) {
  cl <- data.table::as.data.table(clusters_df)
  stopifnot(all(c("chr","start","end","cluster_id") %in% names(cl)))
  
  hits <- data.table::as.data.table(df_sig)
  stopifnot(all(c("CHR","BP","snp","logp") %in% names(hits)))
  
  cl[, `:=`(
    chr = as.integer(chr),
    start = as.integer(start),
    end = as.integer(end),
    cluster_id = as.character(cluster_id)
  )]
  
  hits[, `:=`(
    chr  = as.integer(CHR),
    pos  = as.integer(BP),
    snp  = as.character(snp),
    logp = suppressWarnings(as.numeric(logp))
  )]
  
  cl <- cl[is.finite(chr) & is.finite(start) & is.finite(end)]
  hits <- hits[is.finite(chr) & is.finite(pos)]
  
  cl_iv   <- cl[, .(chr, start, end, cluster_id)]
  hits_iv <- hits[, .(chr, start = pos, end = pos, pos, snp, logp)]
  
  data.table::setkey(cl_iv, chr, start, end)
  data.table::setkey(hits_iv, chr, start, end)
  
  ov <- data.table::foverlaps(hits_iv, cl_iv, type = "within", nomatch = 0L)
  if ("i.cluster_id" %in% names(ov)) data.table::setnames(ov, "i.cluster_id", "cluster_id")
  
  cl[, `:=`(
    n_snps = 0L,
    center = as.integer(round((start + end) / 2)),
    top_snp = NA_character_,
    top_logp = NA_real_
  )]
  
  if (nrow(ov) > 0) {
    stats <- ov[, .(
      n_snps = .N,
      center = as.integer(round(mean(pos, na.rm = TRUE))),
      top_logp = { x <- logp; if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) },
      top_snp  = { x <- logp; if (all(is.na(x))) NA_character_ else snp[which.max(x)] }
    ), by = .(cluster_id)]
    
    cl[stats, on = .(cluster_id), `:=`(
      n_snps  = i.n_snps,
      center  = i.center,
      top_snp = i.top_snp,
      top_logp= i.top_logp
    )]
  }
  
  cl[, cluster_size_kb := round((end - start) / 1000, 2)]
  as.data.frame(cl)
}


# helper normalize cluster names
standardize_cluster_idsXXXgtex <- function(cl) {
  if (!is.data.frame(cl) || !nrow(cl)) return(cl)
  
  cl2 <- cl
  
  # assegura strings
  if ("cluster_chr_n" %in% names(cl2)) cl2$cluster_chr_n <- as.character(cl2$cluster_chr_n)
  if ("cluster_id"    %in% names(cl2)) cl2$cluster_id    <- as.character(cl2$cluster_id)
  
  # 1) Treu prefix "cluster_" (p.ex. "cluster_chr1_2" -> "chr1_2")
  if ("cluster_chr_n" %in% names(cl2)) cl2$cluster_chr_n <- sub("^cluster_", "", cl2$cluster_chr_n)
  if ("cluster_id"    %in% names(cl2)) cl2$cluster_id    <- sub("^cluster_", "", cl2$cluster_id)
  
  # 2) Si no hi ha cluster_chr_n (o queda buit/NA), reconstrueix "chrN_n"
  missing_chr_n <- (!"cluster_chr_n" %in% names(cl2)) ||
    all(is.na(cl2$cluster_chr_n) | !nzchar(cl2$cluster_chr_n))
  
  if (missing_chr_n) {
    if (all(c("chr","cluster_chr") %in% names(cl2))) {
      cl2$cluster_chr_n <- paste0(
        "chr", chr_label_plink(as.integer(cl2$chr)),
        "_", as.integer(cl2$cluster_chr)
      )
    } else if (all(c("chr","cluster") %in% names(cl2))) {
      cl2$cluster_chr_n <- paste0(
        "chr", chr_label_plink(as.integer(cl2$chr)),
        "_", as.integer(cl2$cluster)
      )
    }
  }
  
  # 3) Força cluster_id coherent amb cluster_chr_n (mateix format a tot arreu)
  if ("cluster_chr_n" %in% names(cl2)) {
    cl2$cluster_id <- as.character(cl2$cluster_chr_n)
  }
  
  cl2
}

standardize_cluster_ids <- function(cl) {
  if (!is.data.frame(cl) || !nrow(cl)) return(cl)
  
  cl2 <- cl
  
  # Normalitza noms a caràcter si existeixen
  if ("cluster_chr_n" %in% names(cl2)) {
    cl2$cluster_chr_n <- as.character(cl2$cluster_chr_n)
    cl2$cluster_chr_n <- sub("^cluster_", "", cl2$cluster_chr_n)
  }
  if ("cluster_id" %in% names(cl2)) {
    cl2$cluster_id <- as.character(cl2$cluster_id)
    cl2$cluster_id <- sub("^cluster_", "", cl2$cluster_id)
  }
  
  # Crea cluster_id canònic SEMPRE
  if (!"cluster_id" %in% names(cl2)) cl2$cluster_id <- NA_character_
  
  # Si tenim cluster_chr_n, l’usem per omplir cluster_id
  if ("cluster_chr_n" %in% names(cl2)) {
    cl2$cluster_id[is.na(cl2$cluster_id) | !nzchar(cl2$cluster_id)] <- cl2$cluster_chr_n[is.na(cl2$cluster_id) | !nzchar(cl2$cluster_id)]
  }
  
  # Fallbacks
  if ("label" %in% names(cl2)) {
    idx <- is.na(cl2$cluster_id) | !nzchar(cl2$cluster_id)
    cl2$cluster_id[idx] <- as.character(cl2$label[idx])
  }
  if ("cluster" %in% names(cl2)) {
    idx <- is.na(cl2$cluster_id) | !nzchar(cl2$cluster_id)
    cl2$cluster_id[idx] <- paste0("cluster_", as.character(cl2$cluster[idx]))
  }
  
  # Últim recurs: id determinista
  idx <- is.na(cl2$cluster_id) | !nzchar(cl2$cluster_id)
  if (any(idx)) cl2$cluster_id[idx] <- paste0("cluster_", seq_len(sum(idx)))
  
  cl2
}


# -----------------------------
# Main init: wire build_ranges
# -----------------------------

gi_clusters_canonical_init <- function(session, input, output,
                                       gwas_df,
                                       build_btn_id = "build_ranges",
                                       clusters_dt_id = "cluster_dt",
                                       hits_rows_id = "hits_tbl_rows_selected",
                                       app_count_col = "n_app") {
  
  intervals_raw <- shiny::reactiveVal(tibble::tibble(chr=integer(), start=integer(), end=integer(), label=character()))
  clusters_cur  <- shiny::reactiveVal(tibble::tibble())
  
  safe_num <- function(x, default = NA_real_) {
    y <- suppressWarnings(as.numeric(x))
    if (!is.finite(y)) default else y
  }
  safe_int <- function(x, default = NA_integer_) {
    y <- suppressWarnings(as.integer(x))
    if (!is.finite(y)) default else y
  }
  
  selected_cluster <- shiny::reactive({
    cl <- clusters_cur()
    idx <- input[[paste0(clusters_dt_id, "_rows_selected")]] %||% input$cluster_dt_rows_selected
    if (!is.data.frame(cl) || !nrow(cl) || is.null(idx) || !length(idx)) return(NULL)
    cl[idx[1], , drop = FALSE]
  })
  
  # reset on method change
  shiny::observeEvent(input$cluster_method, {
    intervals_raw(tibble::tibble(chr=integer(), start=integer(), end=integer(), label=character()))
    clusters_cur(tibble::tibble())
    if ("ranges_preview" %in% names(output)) {
      output$ranges_preview <- shiny::renderText("No clusters generated yet.")
    }
  }, ignoreInit = TRUE)
  
  # BUILD
  shiny::observeEvent(input[[build_btn_id]], {
    
    df <- gwas_df()
    shiny::validate(shiny::need(is.data.frame(df) && nrow(df) > 0, "GWAS is empty."))
    shiny::validate(shiny::need(all(c("CHR","BP","snp","Pval","logp") %in% names(df)),
                                "GWAS must have CHR,BP,snp,Pval,logp."))
    
    method <- input$cluster_method %||% "window"
    
    # reset current
    clusters_cur(tibble::tibble())
    intervals_raw(tibble::tibble(chr=integer(), start=integer(), end=integer(), label=character()))
    
    # -----------------------------
    # METHOD 1: window (+/- flank around hits above pthr)
    # -----------------------------
    if (identical(method, "window")) {
      
      thr   <- safe_num(input$pthr, default = -Inf)
      flank <- safe_int(input$flank, default = 0L)
      
      df_sig <- df %>% dplyr::filter(is.finite(.data$logp), .data$logp >= thr) %>% dplyr::arrange(.data$CHR, .data$BP)
      shiny::validate(shiny::need(nrow(df_sig) > 0, "No GWAS hits above threshold."))
      
      # optional: keep only selected hit rows (if exists)
      sel <- NULL
      if (!is.null(input[[hits_rows_id]])) sel <- input[[hits_rows_id]]
      picks <- if (length(sel) > 0) df_sig[sel, , drop = FALSE] else df_sig
      
      raw_int <- tibble::tibble(
        chr   = as.integer(picks$CHR),
        start = pmax(1L, as.integer(picks$BP) - flank),
        end   = as.integer(picks$BP) + flank,
        label = ifelse(!is.na(picks$snp) & nzchar(picks$snp), as.character(picks$snp),
                       paste0("chr", chr_label_plink(as.integer(picks$CHR)), ":", picks$BP))
      ) %>%
        dplyr::arrange(.data$chr, .data$start, .data$end) %>%
        dplyr::mutate(label = make.unique(.data$label, sep = "_"))
      
      intervals_raw(raw_int)
      
      merged <- merge_overlapping_intervals(raw_int)
      cl <- summarize_clusters_from_hits(merged, df_sig)
      
      if (nrow(cl)) {
        cl[[app_count_col]] <- 0L
        clusters_cur(cl)
        intervals_raw(cl %>% dplyr::transmute(chr = .data$chr, start = .data$start, end = .data$end, label = .data$cluster_id))
      }
      
      if ("ranges_preview" %in% names(output)) {
        output$ranges_preview <- shiny::renderText({
          paste0("method=window | thr=", thr, " | flank=", flank, " bp\n",
                 "raw intervals=", nrow(raw_int), " | clusters=", nrow(cl))
        })
      }
      
      return(invisible(TRUE))
    }
    
    # -----------------------------
    # METHOD 2: hits (density-based)
    # hits_mode: span1mb / tiled / sliding
    # -----------------------------
    if (method %in% c("hit","hits")) {
      
      hits_mode <- input$hits_mode %||% "span1mb"
      min_logp  <- safe_num(input$min_logp, default = NA_real_)
      min_hits  <- safe_int(input$min_hits %||% input$min_snps, default = 3L)
      win_bp    <- safe_int(input$win_bp,  default = 1e6L)
      step_bp   <- safe_int(input$step_bp, default = floor(win_bp / 10))
      
      shiny::validate(shiny::need(is.finite(min_logp), "min_logp is missing."))
      shiny::validate(shiny::need(is.finite(min_hits) && min_hits >= 1, "min_hits must be >= 1."))
      shiny::validate(shiny::need(is.finite(win_bp) && win_bp > 0, "win_bp must be > 0."))
      
      df_sig <- df %>% dplyr::filter(is.finite(.data$logp), .data$logp >= min_logp) %>% dplyr::arrange(.data$CHR, .data$BP)
      shiny::validate(shiny::need(nrow(df_sig) > 0, "No GWAS hits above min_logp."))
      
      # --- span1mb: cluster by distance (generalized to win_bp)
      if (identical(hits_mode, "span1mb")) {
        
        clustered <- df_sig %>%
          dplyr::group_by(.data$CHR) %>%
          dplyr::group_modify(~{
            x <- .x
            n <- nrow(x)
            if (n == 0) return(x)
            cid <- integer(n)
            start_idx <- 1
            cur <- 1
            for (i in seq_len(n)) {
              if (x$BP[i] - x$BP[start_idx] > win_bp) {
                cur <- cur + 1
                start_idx <- i
              }
              cid[i] <- cur
            }
            x$cluster_chr <- cid
            x
          }) %>%
          dplyr::ungroup()
        
        cl <- clustered %>%
          dplyr::group_by(.data$CHR, .data$cluster_chr) %>%
          dplyr::summarise(
            chr = as.integer(dplyr::first(.data$CHR)),
            start = min(.data$BP),
            end   = max(.data$BP),
            center = as.integer(round(mean(.data$BP))),
            n_snps = dplyr::n(),
            top_snp  = .data$snp[which.max(.data$logp)],
            top_logp = max(.data$logp),
            cluster_size_kb = round((max(.data$BP) - min(.data$BP)) / 1000, 2),
            .groups = "drop"
          ) %>%
          dplyr::filter(.data$n_snps >= min_hits, .data$top_logp >= min_logp) %>%
          dplyr::arrange(.data$chr, .data$start, .data$end) %>%
          dplyr::mutate(cluster = dplyr::row_number()) %>%
          dplyr::group_by(.data$chr) %>%
          dplyr::mutate(cluster_chr = dplyr::row_number()) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            cluster_chr_n = paste0("chr", chr_label_plink(.data$chr), "_", .data$cluster_chr),
            cluster_id    = .data$cluster_chr_n
          )
        
        if (nrow(cl)) {
          cl[[app_count_col]] <- 0L
          clusters_cur(cl)
          intervals_raw(cl %>% dplyr::transmute(chr = .data$chr, start = .data$start, end = .data$end, label = .data$cluster_id))
        }
        
        return(invisible(TRUE))
      }
      
      # tiled: step = win
      if (identical(hits_mode, "tiled")) step_bp <- win_bp
      shiny::validate(shiny::need(is.finite(step_bp) && step_bp > 0, "step_bp must be > 0."))
      
      # chr lengths (max BP from full GWAS)
      chr_lengths <- df %>%
        dplyr::group_by(.data$CHR) %>%
        dplyr::summarise(chr_len = max(.data$BP, na.rm = TRUE), .groups = "drop") %>%
        dplyr::transmute(chr = as.integer(.data$CHR), chr_len = as.integer(.data$chr_len))
      
      hits <- df_sig %>% dplyr::transmute(chr = as.integer(.data$CHR), pos = as.integer(.data$BP))
      
      windows  <- make_windows_dt(chr_lengths, win_bp = win_bp, step_bp = step_bp)
      win_hits <- count_hits_in_windows(hits, windows)
      
      sig_win <- win_hits[win_hits$n_hits >= min_hits, , drop = FALSE]
      shiny::validate(shiny::need(nrow(sig_win) > 0,
                                  paste0("No windows pass min_hits=", min_hits,
                                         " (mode=", hits_mode, ", win_bp=", win_bp, ", step_bp=", step_bp, ").")))
      
      merged <- merge_significant_windows(sig_win[, c("chr","start","end")], gap_bp = 0L)
      
      cl <- as.data.frame(merged) %>%
        dplyr::arrange(.data$chr, .data$start, .data$end) %>%
        dplyr::group_by(.data$chr) %>%
        dplyr::mutate(cluster_chr = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(cluster = dplyr::row_number()) %>%
        dplyr::mutate(
          cluster_chr_n = paste0("chr", chr_label_plink(.data$chr), "_", .data$cluster_chr),
          cluster_id    = .data$cluster_chr_n
        )
      
      cl <- add_cluster_stats_from_hits(cl, df_sig)
      
      if (nrow(cl)) {
        cl[[app_count_col]] <- 0L
        clusters_cur(cl)
        intervals_raw(cl %>% dplyr::transmute(chr = .data$chr, start = .data$start, end = .data$end, label = .data$cluster_id))
      }
      
      return(invisible(TRUE))
    }
    
  }, ignoreInit = TRUE)
  
  list(
    intervals_raw = intervals_raw,
    clusters_cur  = clusters_cur,
    selected_cluster = selected_cluster
  )
}
