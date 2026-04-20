# R/mod_gwas_hit_priority.R
# ------------------------------------------------------------------
# GWAS hit priority module for GItools Integrator
# ADDITIVE MODULE:
# - Does NOT modify LD/block canonical pipeline
# - Consumes already integrated/canonical objects
# - Builds a GWAS-hit-centric prioritization table
# ------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
}

safe_int0 <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x[is.na(x)] <- 0L
  x
}

safe_num0 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

normalize_cluster_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^cluster_", "", x, ignore.case = TRUE)
  
  x <- ifelse(
    grepl("^[0-9XYM]+_[0-9]+$", x, ignore.case = TRUE),
    paste0("chr", x),
    x
  )
  
  x <- sub("^CHR", "chr", x, ignore.case = TRUE)
  x <- sub("^chr23_", "chrX_", x, ignore.case = TRUE)
  x <- sub("^chr24_", "chrY_", x, ignore.case = TRUE)
  x <- sub("^chr26_", "chrMT_", x, ignore.case = TRUE)
  
  x
}

normalize_classe <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("Catalog_hit", "catalog_hit") ~ "catalog_hit",
    x %in% c("GTEx_hit", "gtex_hit") ~ "gtex_hit",
    x %in% c("NonSyn_hit", "nonsyn_hit") ~ "nonsyn_hit",
    x %in% c("EWASTum_hit", "ewastum_hit") ~ "ewastum_hit",
    x %in% c("EWASDis_hit", "ewasdis_hit") ~ "ewasdis_hit",
    x %in% c("LD_hit", "ld_hit") ~ "ld_hit",
    x == "GWAS" ~ "GWAS",
    TRUE ~ x
  )
}

collapse_chr_unique <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- unique(x)
  paste(sort(x), collapse = "; ")
}

extract_chr_num <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("_.*$", "", x)
  x <- toupper(x)
  
  dplyr::case_when(
    x == "X" ~ 23,
    x == "Y" ~ 24,
    x %in% c("M", "MT") ~ 25,
    TRUE ~ suppressWarnings(as.numeric(x))
  )
}


 # ===========================================================
build_gwas_hit_match_raw_audit <- function(candidates_df) {
  apps_ref <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
  
  if (!is.data.frame(candidates_df) || !nrow(candidates_df)) {
    return(tibble::tibble())
  }
  
  x <- candidates_df %>%
    dplyr::mutate(
      cluster_id = normalize_cluster_id(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      pos = suppressWarnings(as.integer(position)),
      rsid = trimws(as.character(rsid)),
      id_hit = trimws(as.character(id_hit)),
      classe = normalize_classe(classe),
      source_app = dplyr::case_when(
        classe == "catalog_hit" ~ "catalog",
        classe == "gtex_hit" ~ "gtex",
        classe == "nonsyn_hit" ~ "nonsyn",
        classe == "ewasdis_hit" ~ "ewasdis",
        classe == "ewastum_hit" ~ "ewastum",
        classe == "GWAS" ~ "gwas",
        TRUE ~ NA_character_
      )
    )
  
  gwas_hits <- x %>%
    dplyr::filter(
      classe == "GWAS",
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(pos)
    ) %>%
    dplyr::mutate(
      gwas_hit = dplyr::coalesce(dplyr::na_if(rsid, ""), dplyr::na_if(id_hit, "")),
      gwas_hit = dplyr::if_else(
        !is.na(gwas_hit) & nzchar(gwas_hit),
        gwas_hit,
        paste0("pos_", pos)
      )
    ) %>%
    dplyr::transmute(
      cluster_id,
      chr,
      gwas_pos = pos,
      gwas_hit
    ) %>%
    dplyr::distinct()
  
  app_hits <- x %>%
    dplyr::filter(
      !is.na(source_app), source_app %in% apps_ref,
      !is.na(cluster_id), nzchar(cluster_id)
    ) %>%
    dplyr::transmute(
      cluster_id,
      source_app,
      app_pos = pos,
      app_rsid = dplyr::na_if(rsid, ""),
      app_id_hit = dplyr::na_if(id_hit, "")
    )
  
  gwas_hits %>%
    tidyr::crossing(source_app = apps_ref) %>%
    dplyr::left_join(
      app_hits,
      by = c("cluster_id", "source_app"),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      exact_pos_match = is.finite(app_pos) & is.finite(gwas_pos) & (app_pos == gwas_pos),
      exact_id_match =
        (!is.na(app_rsid) & nzchar(app_rsid) & app_rsid == gwas_hit) |
        (!is.na(app_id_hit) & nzchar(app_id_hit) & app_id_hit == gwas_hit),
      verified_match = exact_pos_match | exact_id_match
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      !is.na(gwas_hit), nzchar(gwas_hit)
    ) %>%
    dplyr::transmute(
      cluster_id,
      chr,
      gwas_hit,
      gwas_pos,
      source_app,
      app_rsid = dplyr::coalesce(app_rsid, ""),
      app_id_hit = dplyr::coalesce(app_id_hit, ""),
      app_pos,
      exact_pos_match = as.integer(exact_pos_match),
      exact_id_match = as.integer(exact_id_match),
      verified_match = as.integer(verified_match),
      raw_state = dplyr::if_else(verified_match == 1L, "MATCH", "NO_MATCH")
    ) %>%
    dplyr::arrange(cluster_id, gwas_pos, gwas_hit, source_app, dplyr::desc(verified_match))
}

add_gwas_hit_marker_state <- function(raw_audit_df, block_hits_df, candidates_df) {
  apps_ref <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
  
  if (!is.data.frame(raw_audit_df) || !nrow(raw_audit_df)) {
    return(tibble::tibble())
  }
  
  # ------------------------------------------------------------
  # 1) Presència real d'app al cluster
  #    NOMÉS a partir de *_hit, mai de files GWAS repetides
  # ------------------------------------------------------------
  cluster_app_presence <- cluster_app_presence_from_candidates(candidates_df)
  
  # ------------------------------------------------------------
  # 2) Parelles GWAS hit × app amb relació MARKER estricta
  #    Definició:
  #    - mateix cluster_id
  #    - mateix block_id
  #    - hit real de l'app dins aquell bloc
  # ------------------------------------------------------------
  marker_pairs <- if (is.data.frame(block_hits_df) && nrow(block_hits_df)) {
    bh <- block_hits_df %>%
      dplyr::mutate(
        cluster_id = normalize_cluster_id(cluster_id),
        block_id = trimws(as.character(block_id)),
        gwas_hit = trimws(as.character(lead_snp)),
        gwas_pos = suppressWarnings(as.numeric(lead_pos)),
        source_app = tolower(trimws(as.character(source_app))),
        source_app = dplyr::case_when(
          source_app %in% c("catalog") ~ "catalog",
          source_app %in% c("gtex") ~ "gtex",
          source_app %in% c("nonsyn", "dbnsfp", "non_syn", "non-syn") ~ "nonsyn",
          source_app %in% c("ewasdis", "ewas_disease", "disease") ~ "ewasdis",
          source_app %in% c("ewastum", "ewas_tumor", "tumor") ~ "ewastum",
          TRUE ~ source_app
        ),
        classe = if ("classe" %in% names(.)) normalize_classe(classe) else NA_character_
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(block_id), nzchar(block_id),
        !is.na(gwas_hit), nzchar(gwas_hit),
        is.finite(gwas_pos),
        !is.na(source_app), source_app %in% apps_ref
      )
    
    app_real_hits_in_block <- bh %>%
      dplyr::filter(
        !is.na(classe),
        (
          (source_app == "catalog"  & classe == "catalog_hit")  |
            (source_app == "gtex"     & classe == "gtex_hit")     |
            (source_app == "nonsyn"   & classe == "nonsyn_hit")   |
            (source_app == "ewasdis"  & classe == "ewasdis_hit")  |
            (source_app == "ewastum"  & classe == "ewastum_hit")
        )
      ) %>%
      dplyr::distinct(cluster_id, block_id, source_app) %>%
      dplyr::mutate(app_has_real_hit_in_block = TRUE)
    
    bh %>%
      dplyr::distinct(cluster_id, block_id, gwas_hit, gwas_pos, source_app) %>%
      dplyr::left_join(
        app_real_hits_in_block,
        by = c("cluster_id", "block_id", "source_app")
      ) %>%
      dplyr::filter(app_has_real_hit_in_block) %>%
      dplyr::distinct(cluster_id, gwas_hit, gwas_pos, source_app) %>%
      dplyr::mutate(has_marker_relation = TRUE)
  } else {
    tibble::tibble(
      cluster_id = character(),
      gwas_hit = character(),
      gwas_pos = numeric(),
      source_app = character(),
      has_marker_relation = logical()
    )
  }
  
  # ------------------------------------------------------------
  # 3) Estat final
  #    MATCH > MARKER > NOLINK > NOHIT
  # ------------------------------------------------------------
  raw_audit_df %>%
    dplyr::mutate(
      cluster_id = normalize_cluster_id(cluster_id),
      gwas_hit = trimws(as.character(gwas_hit)),
      gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
      source_app = tolower(trimws(as.character(source_app))),
      source_app = dplyr::case_when(
        source_app %in% c("catalog") ~ "catalog",
        source_app %in% c("gtex") ~ "gtex",
        source_app %in% c("nonsyn", "dbnsfp", "non_syn", "non-syn") ~ "nonsyn",
        source_app %in% c("ewasdis", "ewas_disease", "disease") ~ "ewasdis",
        source_app %in% c("ewastum", "ewas_tumor", "tumor") ~ "ewastum",
        TRUE ~ source_app
      )
    ) %>%
    dplyr::left_join(
      cluster_app_presence,
      by = c("cluster_id", "source_app")
    ) %>%
    dplyr::left_join(
      marker_pairs,
      by = c("cluster_id", "gwas_hit", "gwas_pos", "source_app")
    ) %>%
    dplyr::mutate(
      cluster_has_app = dplyr::coalesce(cluster_has_app, FALSE),
      has_marker_relation = dplyr::coalesce(has_marker_relation, FALSE),
      link_state = dplyr::case_when(
        verified_match == 1L ~ "MATCH",
        has_marker_relation ~ "MARKER",
        cluster_has_app ~ "NOLINK",
        TRUE ~ "NOHIT"
      )
    )
}


gwas_hit_app_heatmap_df <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  
  hit_order <- df %>%
    dplyr::distinct(cluster_id, gwas_hit, gwas_pos, hit_label) %>%
    dplyr::mutate(
      chr_num = extract_chr_num(cluster_id),
      cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
    ) %>%
    dplyr::mutate(
      chr_num = dplyr::coalesce(chr_num, 999),
      cluster_num = dplyr::coalesce(cluster_num, 999)
    ) %>%
    dplyr::arrange(chr_num, cluster_num, gwas_pos, gwas_hit) %>%
    dplyr::pull(hit_label)
  
  df %>%
    dplyr::mutate(
      source_app = factor(source_app, levels = c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")),
      link_state = factor(link_state, levels = c("MATCH", "MARKER", "NOLINK", "NOHIT")),
      hit_label = factor(hit_label, levels = rev(unique(hit_order)))
    )
}

summarize_gwas_hit_app_audit <- function(audit_df) {
  if (!is.data.frame(audit_df) || !nrow(audit_df)) {
    return(tibble::tibble())
  }
  
  audit_df %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      gwas_hit = as.character(gwas_hit),
      gwas_pos = suppressWarnings(as.integer(gwas_pos)),
      source_app = as.character(source_app),
      link_state = dplyr::coalesce(as.character(link_state), "NOHIT"),
      verified_match = dplyr::coalesce(as.integer(verified_match), 0L),
      app_rsid = dplyr::coalesce(as.character(app_rsid), ""),
      app_id_hit = dplyr::coalesce(as.character(app_id_hit), ""),
      app_pos = suppressWarnings(as.integer(app_pos)),
      state_rank = dplyr::case_when(
        link_state == "MATCH" ~ 3L,
        link_state == "MARKER" ~ 2L,
        link_state == "NOLINK" ~ 1L,
        link_state == "NOHIT" ~ 0L,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::group_by(cluster_id, chr, gwas_hit, gwas_pos, source_app) %>%
    dplyr::summarise(
      best_rank = max(state_rank, na.rm = TRUE),
      n_rows = dplyr::n(),
      n_verified = sum(verified_match == 1L, na.rm = TRUE),
      matched_ids = paste(
        sort(unique(c(
          app_rsid[verified_match == 1L & !is.na(app_rsid) & nzchar(app_rsid)],
          app_id_hit[verified_match == 1L & !is.na(app_id_hit) & nzchar(app_id_hit)]
        ))),
        collapse = "; "
      ),
      matched_pos = paste(
        sort(unique(app_pos[verified_match == 1L & is.finite(app_pos)])),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      link_state = dplyr::case_when(
        best_rank >= 3L ~ "MATCH",
        best_rank == 2L ~ "MARKER",
        best_rank == 1L ~ "NOLINK",
        TRUE ~ "NOHIT"
      ),
      hit_label = paste0(cluster_id, " · ", gwas_hit)
    ) %>%
    dplyr::select(
      cluster_id, chr, gwas_hit, gwas_pos, source_app,
      link_state, n_rows, n_verified, matched_ids, matched_pos, hit_label
    ) %>%
    dplyr::arrange(cluster_id, gwas_pos, gwas_hit, source_app)
}

# ---------- helper per generar links

make_matched_ids_links <- function(matched_ids, source_app, link_state, chr, gwas_pos) {
  ids <- trimws(as.character(matched_ids))
  app <- trimws(tolower(as.character(source_app)))
  state <- trimws(as.character(link_state))
  chrv <- as.character(chr)
  posv <- suppressWarnings(as.integer(gwas_pos))
  
  out <- purrr::pmap_chr(
    list(ids, app, state, chrv, posv),
    function(id_str, app_i, state_i, chr_i, pos_i) {
      
      # MARKER -> UCSC +/- 1000 bp al voltant del GWAS hit
      if (!is.na(state_i) && state_i == "MARKER" &&
          !is.na(chr_i) && nzchar(chr_i) && is.finite(pos_i)) {
        
        chr_clean <- gsub("^chr", "", chr_i, ignore.case = TRUE)
        start_pos <- max(1L, pos_i - 10000L)
        end_pos   <- pos_i + 10000L
        
        ucsc_url <- paste0(
          "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr",
          utils::URLencode(chr_clean, reserved = TRUE),
          "%3A",
          start_pos,
          "-",
          end_pos
        )
        
        return(paste0(
          "<a href='", ucsc_url, "' target='_blank'>UCSC±10kb</a>"
        ))
      }
      
      # Resta de casos: links per matched_ids segons app
      if (is.na(id_str) || !nzchar(id_str)) return("")
      
      vals <- unlist(strsplit(id_str, "\\s*;\\s*"))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return("")
      
      links <- vapply(vals, function(v) {
        lab <- htmltools::htmlEscape(v)
        
        href <- switch(
          app_i,
          "catalog" = paste0("https://www.ebi.ac.uk/gwas/variants/", utils::URLencode(v, reserved = TRUE)),
          "gtex"    = paste0("https://gtexportal.org/home/snp/", utils::URLencode(v, reserved = TRUE)),
          "nonsyn"  = paste0("https://www.ncbi.nlm.nih.gov/snp/", utils::URLencode(v, reserved = TRUE)),
          "ewasdis" = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=", utils::URLencode(v, reserved = TRUE)),
          "ewastum" = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=", utils::URLencode(v, reserved = TRUE)),
          NULL
        )
        
        if (is.null(href) || !nzchar(href)) {
          lab
        } else {
          paste0("<a href='", href, "' target='_blank'>", lab, "</a>")
        }
      }, character(1))
      
      paste(links, collapse = "; ")
    }
  )
  
  out
}

# helper links generics

make_simple_links <- function(x, base_url, sep = "\\s*;\\s*") {
  vals <- as.character(x)
  
  vapply(vals, function(v) {
    v <- trimws(v)
    if (is.na(v) || !nzchar(v)) return("")
    
    parts <- unlist(strsplit(v, sep))
    parts <- trimws(parts)
    parts <- parts[!is.na(parts) & nzchar(parts)]
    if (!length(parts)) return("")
    
    links <- vapply(parts, function(p) {
      lab <- htmltools::htmlEscape(p)
      href <- paste0(base_url, utils::URLencode(p, reserved = TRUE))
      paste0("<a href='", href, "' target='_blank'>", lab, "</a>")
    }, character(1))
    
    paste(links, collapse = "; ")
  }, character(1))
}

make_dbsnp_links <- function(x) {
  make_simple_links(x, "https://www.ncbi.nlm.nih.gov/snp/")
}

make_genecards_links <- function(x) {
  make_simple_links(x, "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
}

make_gwascatalog_term_links <- function(x) {
  vals <- as.character(x)
  
  vapply(vals, function(v) {
    v <- trimws(v)
    if (is.na(v) || !nzchar(v)) return("")
    
    lab <- htmltools::htmlEscape(v)
    href <- paste0(
      "https://www.ebi.ac.uk/gwas/search?query=",
      utils::URLencode(v, reserved = TRUE)
    )
    
    paste0("<a href='", href, "' target='_blank'>", lab, "</a>")
  }, character(1))
}

make_omim_term_links <- function(x) {
  vals <- as.character(x)
  
  vapply(vals, function(v) {
    v <- trimws(v)
    if (is.na(v) || !nzchar(v)) return("")
    
    lab <- htmltools::htmlEscape(v)
    href <- paste0(
      "https://www.omim.org/search?search=",
      utils::URLencode(v, reserved = TRUE)
    )
    
    paste0("<a href='", href, "' target='_blank'>", lab, "</a>")
  }, character(1))
}

# helper reusable per priority class
make_priority_badge <- function(x) {
  vals <- as.character(x)
  
  vapply(vals, function(v) {
    v <- trimws(v)
    if (is.na(v) || !nzchar(v)) return("")
    
    bg <- dplyr::case_when(
      v == "High"   ~ "#993404",
      v == "Medium" ~ "#fe9929",
      v == "Low"    ~ "#fff3cd",
      TRUE          ~ "#eeeeee"
    )
    
    fg <- dplyr::case_when(
      v == "High" ~ "#ffffff",
      TRUE        ~ "#000000"
    )
    
    paste0(
      "<span style='display:inline-block; padding:3px 10px; border-radius:999px; font-weight:600; font-size:12px; background:",
      bg,
      "; color:",
      fg,
      ";'>",
      htmltools::htmlEscape(v),
      "</span>"
    )
  }, character(1))
}

# export sempre de tota la taula, no només la pàgina o el filtre actiu
dt_buttons_full <- list(
  list(extend = "copy",  exportOptions = list(modifier = list(page = "all", search = "none"))),
  list(extend = "csv",   exportOptions = list(modifier = list(page = "all", search = "none"))),
  list(extend = "excel", exportOptions = list(modifier = list(page = "all", search = "none")))
)

# ============================================================
# UI
# ============================================================
gwas_hit_priority_module_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::tabsetPanel(
      
      # ==========================================================
      # TAB 1: Priority components
      # ==========================================================
      shiny::tabPanel(
        title = "GWAS hit scores",
        
        shiny::fluidPage(
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::tags$div(
                class = "panelGrey",
                style = "margin-bottom:12px;",
                shiny::tags$div(
                  class = "smallNote",
                  shiny::HTML(
                    paste0(
                      "<b>GWAS hit priority:</b> ",
                      "Each GWAS hit is ranked using cross-app evidence integrated at the variant level. ",
                      "The table reports both its overall score and the detailed structure of supporting evidence across apps. ",
                      "For each hit, the <b>support signature</b> indicates, per app, whether the evidence corresponds to an exact ",
                      "<b>MATCH</b> (same variant) or a <b>MARKER</b> (linked through the same LD block).",
                      
                      "<br><br>",
                      
                      "<b>GWAS hit priority score:</b> ",
                      "The score is computed as the sum of two components, as shown in the <b>score breakdown</b>. ",
                      "First, each app contributes a weight depending on both its evidence type and its support state ",
                      "(typically higher for MATCH than for MARKER). ",
                      "Second, the hit receives an additional <b>GWAS significance contribution</b>, derived from its own statistical significance ",
                      "(higher for more significant GWAS hits). ",
                      "Therefore, the final score reflects both the diversity and quality of cross-app support and the intrinsic strength of the GWAS association itself.",
                      
                      "<br><br>",
                      
                      "<b>Score components:</b> ",
                      "The score can be decomposed into individual contributions shown in the plot as stacked segments. ",
                      "These include the per-app evidence components (for example <i>catalog (MATCH)</i> or <i>gtex (MARKER)</i>) ",
                      "together with the <b>GWAS significance</b> component. ",
                      "The total bar length equals the GWAS-hit priority score, and differences between hits arise from both ",
                      "the number and type of supporting apps and the statistical significance of the GWAS hit.",
                      
                      "<br><br>",
                      
                      "<b>Priority classes:</b> ",
                      "Two complementary classifications are shown.<br>",
                      "<b>Absolute priority class</b> uses fixed score thresholds and supports comparison across sessions.<br>",
                      "<b>Relative priority class</b> uses tertiles computed within the current result table and indicates ",
                      "whether a hit falls in the lower, middle, or upper part of the current score distribution."
                    )
                  )
                )
              ),
       
              shiny::tags$br(),
              
              shiny::div(
                style = "margin-bottom:10px;",
                shiny::actionButton(
                  ns("show_gwas_hit_priority_plot"),
                  "📊 Show GWAS hit priority components",
                  class = "btn-primary"
                ),
                shiny::actionButton(
                  ns("show_gwas_hit_manhattan_plot"),
                  "GWAS hit Manhattan plot"
                )
              )
            ),
            
            shiny::column(
              6,
              shiny::h3("GWAS hits by cluster"),
              shiny::div(
                class = "panelGrey",
                style = "margin-top:12px;",
                withSpinner(plotly::plotlyOutput(ns("gwas_hit_priority_cluster_plot"), height = "300px")),
                uiOutput(ns("gwas_hit_cluster_legend")),
              )
            )
          ),
          
          shiny::tags$br(),
          
          shiny::fluidRow(
            shiny::column(
              12,
              uiOutput(ns("priority_gwas_warning")),
              shiny::h3("GWAS Hit priority scores"),
              withSpinner(DT::DTOutput(ns("gwas_hit_priority_dt")))
            )
          )
        )
      ),
      
      # ==========================================================
      # TAB 2: GWAS x app
      # ==========================================================
      shiny::tabPanel(
        title = "GWAS hits × app",
        
        shiny::fluidPage(
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::tags$div(
                class = "panelGrey",
                style = "margin-bottom:12px;",
                shiny::tags$div(
                  class = "smallNote",
                  shiny::HTML(
                    paste0(
                      "<b>GWAS hit × app relationships:</b> ",
                      "For each GWAS hit and each app, this table summarizes the current relationship state between the hit and the app evidence. ",
                      "These states describe whether the hit is directly supported, indirectly linked, present without linkage, or absent for that app within the cluster."
                    )
                  )
                )
              ),
              
              shiny::tags$div(
                class = "panelGrey",
                style = "margin-top:12px;",
                shiny::tags$div(
                  class = "smallNote",
                  shiny::HTML(
                    paste0(
                      "<b>Interpretation:</b> ",
                      "<br><b>MATCH</b>: the GWAS hit is directly supported by that app.",
                      "<br><b>MARKER</b>: the GWAS hit is linked to that app through the same LD-block context.",
                      "<br><b>NOLINK</b>: the app has evidence in the cluster, but that GWAS hit is not linked to it.",
                      "<br><b>NOHIT</b>: the app has no evidence in that cluster for the evaluated hit."
                    )
                  )
                )
              ),
              
              shiny::tags$br(),
              
              shiny::div(
                style = "margin-bottom:10px;",
                shiny::actionButton(
                  ns("show_gwas_hit_app_heatmap"),
                  "🧮 Show GWAS hit × app heatmap",
                  class = "btn-primary"
                )
              ),
            ),
            
            shiny::column(
              6,
              shiny::h3("GWAS × app by cluster"),
              shiny::div(
                class = "panelGrey",
                style = "margin-top:12px;",
                withSpinner(plotly::plotlyOutput(ns("gwas_hit_match_summary_cluster_plot"), height = "300px")),
                uiOutput(ns("gwas_hit_match_summary_cluster_legend")),
              )
            )
          ),
          
          shiny::tags$br(),
          
          shiny::fluidRow(
            shiny::column(
              12,
              shiny::h3("GWAS × app summary"),
              withSpinner(DT::DTOutput(ns("gwas_hit_match_summary_dt")))
            )
          )
        )
      )
    )
  )
}

# ============================================================
# SERVER
# ============================================================
gwas_hit_priority_module_server <- function(
    id,
    gene_bridge_lookup_df,
    nonsyn_gene_hit_bridge_df,
    integrated_candidates_r,
    block_hits_global_r,
    clusters_consensus_r,
    canonical_prioritized_cluster_genes_r = NULL,
    gene_gwas_hit_score_audit_r = NULL,
    gwas_hit_priority_df_v2_r = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {
    
    ##############################################
    ca_base_r <- shiny::reactive({
      x <- tryCatch(integrated_candidates_r(), error = function(e) NULL)
      if (is.data.frame(x)) x else tibble::tibble()
    })
    
    bh_base_r <- shiny::reactive({
      x <- tryCatch(block_hits_global_r(), error = function(e) NULL)
      if (is.data.frame(x)) x else tibble::tibble()
    })
    
    cl_base_r <- shiny::reactive({
      x <- tryCatch(clusters_consensus_r(), error = function(e) NULL)
      if (is.data.frame(x)) x else tibble::tibble()
    })
    
    gene_bridge_base_r <- shiny::reactive({
      x <- tryCatch(gene_bridge_lookup_df(), error = function(e) NULL)
      if (is.data.frame(x)) x else tibble::tibble()
    })
    
    nonsyn_gene_hit_bridge_base_r <- shiny::reactive({
      x <- tryCatch(nonsyn_gene_hit_bridge_df(), error = function(e) NULL)
      if (is.data.frame(x)) x else tibble::tibble()
    })
    
    gene_bridge_by_cluster_app_r <- shiny::reactive({
      gl <- gene_bridge_base_r()
      
      if (!is.data.frame(gl) || !nrow(gl)) {
        return(tibble::tibble(
          cluster_id = character(),
          source_app = character(),
          gene = character(),
          gene_mid = numeric()
        ))
      }
      
      out <- gl %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          source_app = tolower(trimws(as.character(source_app))),
          gene = trimws(as.character(gene))
        )
      
      if (!"gene_mid" %in% names(out)) {
        if (all(c("start", "end") %in% names(out))) {
          out <- out %>%
            dplyr::mutate(
              gene_mid = (
                suppressWarnings(as.numeric(start)) +
                  suppressWarnings(as.numeric(end))
              ) / 2
            )
        } else {
          out$gene_mid <- NA_real_
        }
      } else {
        out <- out %>%
          dplyr::mutate(gene_mid = suppressWarnings(as.numeric(gene_mid)))
      }
      
      out %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(source_app), nzchar(source_app),
          !is.na(gene), nzchar(gene)
        ) %>%
        dplyr::distinct(cluster_id, source_app, gene, gene_mid)
    })
    
    
#############################################################################
    
    ns <- session$ns
    js_id <- gsub("[^A-Za-z0-9_]", "_", ns(""))
    
    
    shiny::insertUI(
      selector = "body",
      where = "beforeEnd",
      ui = shiny::tags$script(shiny::HTML(sprintf("
    if (!window.__gwasHitFilterHandler_%s) {
      window.__gwasHitFilterHandler_%s = true;
      
      Shiny.addCustomMessageHandler('%s', function(message) {
        var cid = (message && message.cluster_id) ? message.cluster_id : '';
        
        if (!window.gwasHitPriorityTable_%s || window.gwasHitPriorityClusterCol_%s === undefined) {
          console.log('[filter_gwashit_table_by_cluster] table not ready');
          return;
        }
        
        var tbl = window.gwasHitPriorityTable_%s;
        var col = window.gwasHitPriorityClusterCol_%s;
        
        tbl.search('');
        tbl.columns().search('');
        
        if (cid !== '') {
          var esc = $.fn.dataTable.util.escapeRegex(cid);
          tbl.column(col).search('^' + esc + '$', true, false).draw();
        } else {
          tbl.draw();
        }
        
        console.log('[filter_gwashit_table_by_cluster] filtered to', cid);
      });
    }
  ", js_id, js_id, ns("filter_gwashit_table_by_cluster"),
                                                  js_id, js_id, js_id, js_id)))
    )
    
   
    ##### Barplot 
    gwas_hit_priority_cluster_plot_df <- reactive({
      hits_df <- gwas_hit_priority_df_v2_r()
      cl_df   <- clusters_consensus_r()
      
      validate(
        need(is.data.frame(cl_df) && nrow(cl_df) > 0, "No cluster information available."),
        need("cluster_id" %in% names(cl_df), "clusters_consensus_r() must contain cluster_id.")
      )
      
      extract_chr_num <- function(x) {
        x <- as.character(x)
        x <- trimws(x)
        x <- sub("^chr", "", x, ignore.case = TRUE)
        x <- sub("_.*$", "", x)
        x <- toupper(x)
        
        dplyr::case_when(
          x == "X" ~ 23,
          x == "Y" ~ 24,
          x %in% c("M", "MT") ~ 25,
          TRUE ~ suppressWarnings(as.numeric(x))
        )
      }
      
      # ------------------------------------------------------------
      # 1) Base de TOTS els clusters
      # ------------------------------------------------------------
      all_clusters <- cl_df %>%
        dplyr::mutate(
          cluster_id = as.character(cluster_id)
        ) %>%
        dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
        dplyr::distinct(cluster_id) %>%
        dplyr::mutate(
          chr_num = extract_chr_num(cluster_id),
          cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
        ) %>%
        dplyr::mutate(
          chr_num = dplyr::coalesce(chr_num, 999),
          cluster_num = dplyr::coalesce(cluster_num, 999)
        )
      
      # ------------------------------------------------------------
      # 2) Resum dels hits prioritzats
      # ------------------------------------------------------------
      hits_summary <- if (is.data.frame(hits_df) && nrow(hits_df) > 0) {
        hits_df %>%
          dplyr::mutate(
            cluster_id = as.character(cluster_id),
            gwas_hit_score = suppressWarnings(as.numeric(gwas_hit_score))
          ) %>%
          dplyr::filter(!is.na(cluster_id), nzchar(cluster_id)) %>%
          dplyr::group_by(cluster_id) %>%
          dplyr::summarise(
            n_hits = dplyr::n(),
            mean_hit_score = if (all(is.na(gwas_hit_score))) 0 else mean(gwas_hit_score, na.rm = TRUE),
            best_hit_score = if (all(is.na(gwas_hit_score))) 0 else max(gwas_hit_score, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        tibble::tibble(
          cluster_id = character(0),
          n_hits = numeric(0),
          mean_hit_score = numeric(0),
          best_hit_score = numeric(0)
        )
      }
      
      # ------------------------------------------------------------
      # 3) Join -> els clusters sense score queden a 0
      # ------------------------------------------------------------
      out <- all_clusters %>%
        dplyr::left_join(hits_summary, by = "cluster_id") %>%
        dplyr::mutate(
          n_hits = dplyr::coalesce(as.numeric(n_hits), 0),
          mean_hit_score = dplyr::coalesce(as.numeric(mean_hit_score), 0),
          best_hit_score = dplyr::coalesce(as.numeric(best_hit_score), 0)
        )
      
      cluster_levels <- out %>%
        dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
        dplyr::pull(cluster_id) %>%
        rev()
      
      # ------------------------------------------------------------
      # 4) Colors
      # ------------------------------------------------------------
      has_positive_score <- any(is.finite(out$mean_hit_score) & out$mean_hit_score > 0)
      
      if (has_positive_score) {
        rng <- range(out$mean_hit_score[out$mean_hit_score > 0], na.rm = TRUE)
        
        if (!all(is.finite(rng))) {
          out$fill_col <- ifelse(out$mean_hit_score > 0, "#d9ead3", "#d9d9d9")
        } else if (rng[1] == rng[2]) {
          out$fill_col <- ifelse(out$mean_hit_score > 0, "#6aaed6", "#d9d9d9")
        } else {
          pal_fun <- scales::col_numeric(
            palette = c("#d9d9d9", "#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
            domain = rng,
            na.color = "#d9d9d9"
          )
          
          out$fill_col <- ifelse(
            out$mean_hit_score > 0,
            pal_fun(out$mean_hit_score),
            "#d9d9d9"
          )
        }
      } else {
        out$fill_col <- "#d9d9d9"
      }
      
      out %>%
        dplyr::mutate(
          cluster_id_chr = cluster_id,
          cluster_id = factor(cluster_id, levels = cluster_levels),
          hover_txt = paste0(
            "Cluster: ", cluster_id_chr,
            "<br>Prioritized hits: ", n_hits,
            "<br>Mean hit score: ", round(mean_hit_score, 2),
            "<br>Best hit score: ", round(best_hit_score, 2)
          )
        )
    })
    
    output$gwas_hit_priority_cluster_plot <- plotly::renderPlotly({
      df <- gwas_hit_priority_cluster_plot_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS hit cluster plot data available.")
      )
      
      plotly::plot_ly(
        data = df,
        source = session$ns("gwas_hit_priority_cluster_plot"),
        y = ~cluster_id,
        x = ~n_hits,
        type = "bar",
        orientation = "h",
        customdata = ~cluster_id_chr,
        hovertext = ~hover_txt,
        hoverinfo = "text",
        marker = list(
          color = df$fill_col,
          line = list(color = "black", width = 0.5)
        ),
        showlegend = FALSE
      ) %>%
        plotly::layout(
          clickmode = "event+select",
          xaxis = list(
            title = "Number of prioritized GWAS hits",
            automargin = TRUE
          ),
          yaxis = list(
            title = "Cluster",
            categoryorder = "array",
            categoryarray = levels(df$cluster_id),
            automargin = TRUE
          ),
          margin = list(l = 90, r = 20, t = 20, b = 60)
        )
    })
    
    output$gwas_hit_cluster_legend <- shiny::renderUI({
      df <- gwas_hit_priority_cluster_plot_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "")
      )
      
      vals <- suppressWarnings(as.numeric(df$mean_hit_score))
      vals <- vals[is.finite(vals)]
      
      if (!length(vals)) return(NULL)
      
      vmin <- min(vals, na.rm = TRUE)
      vmax <- max(vals, na.rm = TRUE)
      
      shiny::div(
        style = "margin-top:8px; font-size:12px; color:#444;",
        shiny::tags$span(style = "font-weight:600; margin-right:10px;", "Mean hit score"),
        shiny::tags$span(
          style = paste0(
            "display:inline-block; width:220px; height:12px; vertical-align:middle; ",
            "border:1px solid #999; border-radius:6px; margin-right:8px; ",
            "background: linear-gradient(to right, ",
            "#d9d9d9 0%, ",
            "#fff7bc 20%, ",
            "#fec44f 40%, ",
            "#fe9929 60%, ",
            "#d95f0e 80%, ",
            "#993404 100%);"
          )
        ),
        shiny::tags$span(formatC(vmin, format = "f", digits = 2)),
        shiny::tags$span(" – "),
        shiny::tags$span(formatC(vmax, format = "f", digits = 2))
      )
    })
    
    
    observeEvent(plotly::event_data("plotly_click", source = session$ns("gwas_hit_priority_cluster_plot")), {
      ed <- plotly::event_data("plotly_click", source = session$ns("gwas_hit_priority_cluster_plot"))
      req(ed)
      
      cid <- ed$customdata[[1]] %||% ""
      req(nzchar(cid))
      
      session$sendCustomMessage(
        session$ns("filter_gwashit_table_by_cluster"),
        list(cluster_id = cid)
      )
    })
    
#   observeEvent(input$clear_gwashit_table_filter, {
#     session$sendCustomMessage(
#       session$ns("filter_gwashit_table_by_cluster"),
#       list(cluster_id = "")
#     )
#   })
    

 # ==================== Heatmap & tables ==================================
  
    gwas_hit_match_audit_summary_with_gene_cache <- shiny::reactiveVal(NULL)
    
    observeEvent(
      list(
        gwas_hit_match_audit_summary_df(),
        gene_bridge_by_cluster_app_r(),
        nonsyn_gene_hit_bridge_base_r()
      ),
      {
        gwas_hit_match_audit_summary_with_gene_cache(NULL)
      },
      ignoreInit = TRUE
    )
    
    gwas_hit_match_audit_summary_with_gene_df_XXXX <- reactive({
      cached <- gwas_hit_match_audit_summary_with_gene_cache()
      if (is.data.frame(cached)) {
        return(cached)
      }
      
      sm <- gwas_hit_match_audit_summary_df()
      gl <- gene_bridge_by_cluster_app_r()
      nhb <- nonsyn_gene_hit_bridge_base_r()
      
      if (!is.data.frame(sm) || !nrow(sm)) {
        return(tibble::tibble())
      }
      
      # base
      base <- sm %>%
        dplyr::mutate(
          .row_id = dplyr::row_number(),
          cluster_id = trimws(as.character(cluster_id)),
          source_app = tolower(trimws(as.character(source_app))),
          gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
          link_state = trimws(as.character(link_state)),
          matched_ids = dplyr::coalesce(as.character(matched_ids), ""),
          nearest_gene = "",
          nearest_gene_dist_bp = NA_real_
        )
      
      # ------------------------------------------------------------
      # 1) NONSYN MATCH  -> matched_ids -> hit_key -> gene
      # ------------------------------------------------------------
      nonsyn_match_rows <- base %>%
        dplyr::filter(source_app == "nonsyn", link_state == "MATCH")
      
      nonsyn_match_ann <- if (nrow(nonsyn_match_rows) &&
                              is.data.frame(nhb) && nrow(nhb) &&
                              all(c("cluster_id", "source_app", "hit_key", "gene") %in% names(nhb))) {
        
        nhb2 <- nhb %>%
          dplyr::mutate(
            cluster_id = trimws(as.character(cluster_id)),
            source_app = tolower(trimws(as.character(source_app))),
            hit_key = trimws(as.character(hit_key)),
            gene = trimws(as.character(gene))
          ) %>%
          dplyr::filter(
            source_app == "nonsyn",
            !is.na(cluster_id), nzchar(cluster_id),
            !is.na(hit_key), nzchar(hit_key),
            !is.na(gene), nzchar(gene)
          ) %>%
          dplyr::distinct(cluster_id, source_app, hit_key, gene)
        
        nonsyn_match_rows %>%
          dplyr::select(.row_id, cluster_id, source_app, matched_ids) %>%
          tidyr::separate_rows(matched_ids, sep = "\\s*;\\s*") %>%
          dplyr::mutate(
            hit_key = trimws(matched_ids)
          ) %>%
          dplyr::filter(!is.na(hit_key), nzchar(hit_key)) %>%
          dplyr::left_join(
            nhb2,
            by = c("cluster_id", "source_app", "hit_key"),
            relationship = "many-to-many"
          ) %>%
          dplyr::filter(!is.na(gene), nzchar(gene)) %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            nearest_gene = paste(sort(unique(gene)), collapse = "; "),
            nearest_gene_dist_bp = 0,
            .groups = "drop"
          )
      } else {
        tibble::tibble(
          .row_id = integer(),
          nearest_gene = character(),
          nearest_gene_dist_bp = numeric()
        )
      }
      
      # ------------------------------------------------------------
      # 2) MATCH / MARKER (excepte nonsyn MATCH)
      #    nearest gene by cluster_id + source_app
      # ------------------------------------------------------------
      other_rows <- base %>%
        dplyr::filter(
          link_state %in% c("MATCH", "MARKER"),
          !(source_app == "nonsyn" & link_state == "MATCH")
        )
      
      other_ann <- if (nrow(other_rows) &&
                       is.data.frame(gl) && nrow(gl) &&
                       all(c("cluster_id", "source_app", "gene") %in% names(gl))) {
        
        gl2 <- gl %>%
          dplyr::mutate(
            cluster_id = trimws(as.character(cluster_id)),
            source_app = tolower(trimws(as.character(source_app))),
            gene = trimws(as.character(gene)),
            gene_mid = suppressWarnings(as.numeric(gene_mid))
          ) %>%
          dplyr::filter(
            !is.na(cluster_id), nzchar(cluster_id),
            !is.na(source_app), nzchar(source_app),
            !is.na(gene), nzchar(gene)
          ) %>%
          dplyr::distinct(cluster_id, source_app, gene, gene_mid)
        
        joined <- other_rows %>%
          dplyr::select(.row_id, cluster_id, source_app, gwas_pos) %>%
          dplyr::left_join(
            gl2,
            by = c("cluster_id", "source_app"),
            relationship = "many-to-many"
          )
        
        # A) files amb gene_mid usable -> nearest by distance
        ann_dist <- joined %>%
          dplyr::filter(is.finite(gwas_pos), is.finite(gene_mid)) %>%
          dplyr::mutate(
            dist_bp = abs(gene_mid - gwas_pos)
          ) %>%
          dplyr::arrange(.row_id, dist_bp, gene) %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            nearest_gene = dplyr::first(gene),
            nearest_gene_dist_bp = dplyr::first(dist_bp),
            .groups = "drop"
          )
        
        # B) files sense gene_mid usable -> collapse genes
        ann_nodist <- joined %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            any_dist = any(is.finite(gwas_pos) & is.finite(gene_mid)),
            genes_txt = paste(sort(unique(gene[!is.na(gene) & nzchar(gene)])), collapse = "; "),
            .groups = "drop"
          ) %>%
          dplyr::filter(!any_dist, nzchar(genes_txt)) %>%
          dplyr::transmute(
            .row_id,
            nearest_gene = genes_txt,
            nearest_gene_dist_bp = NA_real_
          )
        
        dplyr::bind_rows(ann_dist, ann_nodist)
      } else {
        tibble::tibble(
          .row_id = integer(),
          nearest_gene = character(),
          nearest_gene_dist_bp = numeric()
        )
      }
      
      # ------------------------------------------------------------
      # 3) Merge final
      # ------------------------------------------------------------
      ann <- dplyr::bind_rows(nonsyn_match_ann, other_ann) %>%
        dplyr::group_by(.row_id) %>%
        dplyr::summarise(
          nearest_gene = dplyr::first(nearest_gene),
          nearest_gene_dist_bp = dplyr::first(nearest_gene_dist_bp),
          .groups = "drop"
        )
      
      out <- base %>%
        dplyr::select(-nearest_gene, -nearest_gene_dist_bp) %>%
        dplyr::left_join(ann, by = ".row_id") %>%
        dplyr::mutate(
          nearest_gene = dplyr::coalesce(nearest_gene, ""),
          nearest_gene_dist_bp = suppressWarnings(as.numeric(nearest_gene_dist_bp))
        ) %>%
        dplyr::select(-.row_id)
      
      gwas_hit_match_audit_summary_with_gene_cache(out)
      out
    })
    
    ##############################
    gwas_hit_gene_support_long_df <- reactive({
      df <- gwas_hit_match_audit_summary_with_gene_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS × app gene support data available."),
        need(all(c(
          "cluster_id", "chr", "gwas_hit", "gwas_pos",
          "source_app", "link_state", "nearest_gene",
          "nearest_gene_dist_bp", "matched_ids", "matched_pos",
          "n_rows", "n_verified"
        ) %in% names(df)), "Required columns missing in GWAS × app gene support source.")
      )
      
      app_levels <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
      state_levels <- c("MATCH", "MARKER")
      
      out <- df %>%
        dplyr::mutate(
          cluster_id = trimws(as.character(cluster_id)),
          chr = as.character(chr),
          gwas_hit = trimws(as.character(gwas_hit)),
          gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
          source_app = tolower(trimws(as.character(source_app))),
          link_state = trimws(as.character(link_state)),
          nearest_gene = dplyr::coalesce(as.character(nearest_gene), ""),
          nearest_gene_dist_bp = suppressWarnings(as.numeric(nearest_gene_dist_bp)),
          matched_ids = dplyr::coalesce(as.character(matched_ids), ""),
          matched_pos = dplyr::coalesce(as.character(matched_pos), ""),
          n_rows = dplyr::coalesce(as.integer(n_rows), 0L),
          n_verified = dplyr::coalesce(as.integer(n_verified), 0L)
        ) %>%
        dplyr::filter(
          link_state %in% state_levels,
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(gwas_hit), nzchar(gwas_hit),
          !is.na(source_app), nzchar(source_app),
          !is.na(nearest_gene), nzchar(nearest_gene)
        ) %>%
        tidyr::separate_rows(nearest_gene, sep = "\\s*;\\s*") %>%
        dplyr::mutate(
          gene = trimws(nearest_gene),
          source_app = factor(source_app, levels = app_levels),
          link_state = factor(link_state, levels = state_levels)
        ) %>%
        dplyr::filter(!is.na(gene), nzchar(gene)) %>%
        dplyr::transmute(
          cluster_id,
          chr,
          gwas_hit,
          gwas_pos,
          source_app = as.character(source_app),
          link_state = as.character(link_state),
          gene,
          gene_dist_bp = nearest_gene_dist_bp,
          matched_ids,
          matched_pos,
          n_rows,
          n_verified
        ) %>%
        dplyr::distinct(
          cluster_id, gwas_hit, source_app, link_state, gene,
          .keep_all = TRUE
        ) %>%
        dplyr::arrange(
          cluster_id,
          gwas_pos,
          gwas_hit,
          factor(source_app, levels = app_levels),
          factor(link_state, levels = state_levels),
          gene
        )
      
      out
    })
    ##############################################################
    
    gwas_hit_match_raw_audit_df <- shiny::reactive({
      ca <- ca_base_r()
      build_gwas_hit_match_raw_audit(ca)
    })
    
    gwas_hit_match_audit_df <- shiny::reactive({
      raw <- gwas_hit_match_raw_audit_df()
      ca  <- ca_base_r()
      bh  <- bh_base_r()
      
      add_gwas_hit_marker_state(
        raw_audit_df = raw,
        block_hits_df = bh,
        candidates_df = ca
      )
    })
    
    gwas_hit_match_audit_summary_df_XXXX <- shiny::reactive({
      au <- gwas_hit_match_audit_df()
      summarize_gwas_hit_app_audit(au)
    })
    
    output$gwas_hit_match_summary_dt_XXXX <- DT::renderDT({
     # df <- gwas_hit_match_audit_summary_df()
      df <- gwas_hit_match_audit_summary_with_gene_df()
      
      shiny::validate(
        shiny::need(is.data.frame(df) && nrow(df) > 0, "No summarized GWAS hit × app table available.")
      )
      
      show_df <- df %>%
        dplyr::transmute(
          cluster_id,
          chr,
          gwas_hit,
          gwas_pos,
          source_app,
          link_state_raw = link_state,
          n_rows = dplyr::coalesce(as.integer(n_rows), 0L),
          n_verified_raw = dplyr::coalesce(as.integer(n_verified), 0L),
          matched_ids_raw = dplyr::coalesce(as.character(matched_ids), ""),
          matched_pos = dplyr::coalesce(as.character(matched_pos), ""),
          nearest_gene = dplyr::coalesce(as.character(nearest_gene), ""),
          nearest_gene_dist_kb = suppressWarnings(as.numeric(nearest_gene_dist_bp) / 1000)
        ) %>%
        dplyr::mutate(
          matched_ids = make_matched_ids_links(
            matched_ids = matched_ids_raw,
            source_app = source_app,
            link_state = link_state_raw,
            chr = chr,
            gwas_pos = gwas_pos
          ),
          nearest_gene = dplyr::if_else(
            nzchar(nearest_gene),
            make_genecards_links(nearest_gene),
            ""
          ),
          nearest_gene_dist_kb = dplyr::if_else(
            is.finite(nearest_gene_dist_kb),
            round(nearest_gene_dist_kb, 2),
            NA_real_
          ),
          link_state = dplyr::case_when(
            link_state_raw == "MATCH" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#993404; color:white; font-weight:700;'>MATCH</span>",
            link_state_raw == "MARKER" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#fe9929; color:black; font-weight:700;'>MARKER</span>",
            link_state_raw == "NOLINK" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#fff3cd; color:black; font-weight:700;'>NOLINK</span>",
            TRUE ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#f2f2f2; color:black; font-weight:700;'>NOHIT</span>"
          ),
          n_verified = dplyr::case_when(
            n_verified_raw > 0 ~ paste0(
              "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#d9ead3; color:#1b5e20; font-weight:700;'>",
              n_verified_raw, "</span>"
            ),
            TRUE ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#f4cccc; color:#7f0000; font-weight:700;'>0</span>"
          )
        ) %>%
        dplyr::select(
          cluster_id, chr, gwas_hit, gwas_pos, source_app, link_state,
          n_rows, n_verified, matched_ids, matched_pos, nearest_gene, nearest_gene_dist_kb
        )
      
      DT::datatable(
        show_df,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        escape = FALSE,
        extensions = "Buttons",
        callback = DT::JS(sprintf("
    var tbl = table;
    window.gwasxappSummaryTable = tbl;
    window.gwasxappSummaryClusterCol = %d;
    
    tbl.columns.adjust();
    
    $(window).on('resize.gwasxappSummaryTable', function() {
      tbl.columns.adjust();
    });
    
    $(document).on('shown.bs.tab.gwasxappSummaryTable', 'button[data-bs-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function() {
      $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
    });
  ", which(names(show_df) == "cluster_id") - 1L)),
        options = list(
          scrollX = TRUE,
          pageLength = 15,
          autoWidth = FALSE,
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          )
        )
      ) %>% DT::formatRound("nearest_gene_dist_kb", 2)
    }, server = FALSE)

    output$gwas_hit_app_heatmap_popup <- plotly::renderPlotly({
      df <- gwas_hit_app_link_selected_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS hit × app data available for the selected cluster(s).")
      )
      
      plot_df <- gwas_hit_app_heatmap_df(df)
      
      validate(
        need(is.data.frame(plot_df) && nrow(plot_df) > 0, "No plot data available.")
      )
      
      app_levels <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
      
      state_to_num <- c(
        "NOHIT" = 0,
        "NOLINK" = 1,
        "MARKER" = 2,
        "MATCH" = 3
      )
      
      plot_df <- plot_df %>%
        dplyr::mutate(
          source_app = factor(as.character(source_app), levels = app_levels),
          hit_label = as.character(hit_label),
          link_state = factor(as.character(link_state), levels = c("MATCH", "MARKER", "NOLINK", "NOHIT")),
          z_value = unname(state_to_num[as.character(link_state)]),
          hover_txt = paste0(
            "<b>GWAS hit:</b> ", htmltools::htmlEscape(hit_label),
            "<br><b>App:</b> ", as.character(source_app),
            "<br><b>State:</b> ", as.character(link_state)
          )
        )
      
      hit_order <- rev(unique(plot_df$hit_label))
      
      p <- plotly::plot_ly(
        data = plot_df,
        x = ~source_app,
        y = ~hit_label,
        z = ~z_value,
        type = "heatmap",
        text = ~hover_txt,
        hoverinfo = "text",
        colorscale = list(
          c(0.00, "#f2f2f2"),
          c(0.24, "#f2f2f2"),
          c(0.25, "#fff3cd"),
          c(0.49, "#fff3cd"),
          c(0.50, "#fe9929"),
          c(0.74, "#fe9929"),
          c(0.75, "#993404"),
          c(1.00, "#993404")
        ),
        zmin = 0,
        zmax = 3,
        showscale = FALSE,
        xgap = 2,
        ygap = 2
      )
      
      p <- plotly::layout(
        p,
        xaxis = list(
          title = "",
          side = "top",
          tickangle = 0,
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = app_levels
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = hit_order
        ),
        margin = list(l = 220, r = 20, t = 40, b = 40),
        plot_bgcolor = "white",
        paper_bgcolor = "white"
      )
      
      p
    })
    
    output$gwas_hit_heatmap_legend_ui <- renderUI({
      tags$div(
        style = "font-size:12px; color:#444; line-height:1.5; padding-top:4px;",
        tags$div(tags$b("Heatmap legend")),
        tags$div(HTML("<span style='color:#993404; font-weight:700;'>MATCH</span> = exact match with app hit")),
        tags$div(HTML("<span style='color:#fe9929; font-weight:700;'>MARKER</span> = same block as app hit")),
        tags$div(HTML("<span style='color:#b8860b; font-weight:700;'>NOLINK</span> = app hit present in cluster but not linked")),
        tags$div(HTML("<span style='color:#808080; font-weight:700;'>NOHIT</span> = no app hit in cluster"))
      )
    })
    
    observeEvent(input$show_gwas_hit_app_heatmap, {
      showModal(
        modalDialog(
          title = "GWAS hit × app map",
          size = "l",
          easyClose = TRUE,
          footer = modalButton("Close"),
          tags$p(
            "Select one cluster to inspect the relationship between GWAS hits and apps."
          ),
          shiny::fluidRow(
            shiny::column(
              width = 5,
              uiOutput(session$ns("gwas_hit_heatmap_cluster_ui"))
            ),
            shiny::column(
              width = 7,
              uiOutput(session$ns("gwas_hit_heatmap_legend_ui"))
            )
          ),
          tags$br(),
          uiOutput(session$ns("gwas_hit_app_heatmap_ui"))
        )
      )
    })
    
  # ============================================================
    
    output$gwas_hit_heatmap_cluster_ui <- renderUI({
      df <- gwas_hit_match_audit_summary_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "")
      )
      
      clusters <- df %>%
        dplyr::distinct(cluster_id) %>%
        dplyr::mutate(
          chr_num = extract_chr_num(cluster_id),
          cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
        ) %>%
        dplyr::mutate(
          chr_num = dplyr::coalesce(chr_num, 999),
          cluster_num = dplyr::coalesce(cluster_num, 999)
        ) %>%
        dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
        dplyr::pull(cluster_id)
      
      shiny::selectizeInput(
        session$ns("gwas_hit_heatmap_cluster"),
        "Cluster(s)",
        choices = c("All", clusters),
        selected = "All",
        multiple = TRUE,
        options = list(
          plugins = list("remove_button"),
          placeholder = "Select one or more clusters"
        )
      )
    })
    
    gwas_hit_app_link_selected_df <- reactive({
      df <- gwas_hit_match_audit_summary_df()
      req(is.data.frame(df), nrow(df) > 0)
      req(input$gwas_hit_heatmap_cluster)
      
      sel <- input$gwas_hit_heatmap_cluster
      sel <- as.character(sel)
      sel <- sel[!is.na(sel) & nzchar(sel)]
      
      if (!length(sel) || "All" %in% sel) {
        return(df)
      }
      
      df %>%
        dplyr::filter(cluster_id %in% sel)
    })
    
    output$gwas_hit_app_heatmap_ui <- renderUI({
      df <- try(gwas_hit_app_link_selected_df(), silent = TRUE)
      
      if (inherits(df, "try-error") || !is.data.frame(df) || !nrow(df)) {
        return(tags$div("No GWAS hit × app data available for the selected cluster(s)."))
      }
      
      plot_df <- try(gwas_hit_app_heatmap_df(df), silent = TRUE)
      
      if (inherits(plot_df, "try-error") || !is.data.frame(plot_df) || !nrow(plot_df)) {
        return(tags$div("No GWAS hit × app plot data available for the selected cluster(s)."))
      }
      
      n_hits <- dplyr::n_distinct(plot_df$hit_label)
      h <- max(350, min(1600, 20 * n_hits))
      
      plotly::plotlyOutput(
        session$ns("gwas_hit_app_heatmap_popup"),
        height = paste0(h, "px")
      )
    })
  
  # ============================================================
    #----------------------------------------------------------------------------
    gwas_hit_match_summary_cluster_plot_df <- reactive({
      df <- gwas_hit_match_audit_summary_df()
      cl <- cl_base_r()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS × app summary data available."),
        need(all(c("cluster_id", "link_state") %in% names(df)), "Required columns missing in GWAS × app summary."),
        need(is.data.frame(cl) && nrow(cl) > 0, "No cluster consensus table available.")
      )
      
      cl2 <- cl %>%
        dplyr::mutate(
          cluster_id = as.character(cluster_id)
        )
      
      if (!"cluster_size_kb" %in% names(cl2)) {
        start_col <- intersect(c("start", "pos_ini", "cluster_start"), names(cl2))[1]
        end_col   <- intersect(c("end", "pos_end", "cluster_end"), names(cl2))[1]
        size_col  <- intersect(c("size_kb", "cluster_size_kb"), names(cl2))[1]
        
        if (!is.na(size_col) && nzchar(size_col)) {
          cl2$cluster_size_kb <- suppressWarnings(as.numeric(cl2[[size_col]]))
        } else if (!is.na(start_col) && nzchar(start_col) && !is.na(end_col) && nzchar(end_col)) {
          cl2$cluster_size_kb <- (
            suppressWarnings(as.numeric(cl2[[end_col]])) -
              suppressWarnings(as.numeric(cl2[[start_col]])) + 1
          ) / 1000
        } else {
          cl2$cluster_size_kb <- NA_real_
        }
      }
      
      cluster_sizes <- cl2 %>%
        dplyr::transmute(
          cluster_id,
          cluster_size_kb = pmax(dplyr::coalesce(as.numeric(cluster_size_kb), 0), 0)
        ) %>%
        dplyr::distinct(cluster_id, .keep_all = TRUE)
      
      cluster_order <- df %>%
        dplyr::distinct(cluster_id) %>%
        dplyr::mutate(
          chr_num = extract_chr_num(cluster_id),
          cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
        ) %>%
        dplyr::mutate(
          chr_num = dplyr::coalesce(chr_num, 999),
          cluster_num = dplyr::coalesce(cluster_num, 999)
        ) %>%
        dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
        dplyr::pull(cluster_id) %>%
        rev()
      
      df %>%
        dplyr::mutate(
          link_state = factor(link_state, levels = c("MATCH", "MARKER", "NOLINK", "NOHIT"))
        ) %>%
        dplyr::count(cluster_id, link_state, name = "n") %>%
        tidyr::complete(
          cluster_id,
          link_state = factor(c("MATCH", "MARKER", "NOLINK", "NOHIT"),
                              levels = c("MATCH", "MARKER", "NOLINK", "NOHIT")),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::mutate(
          total_n = sum(n, na.rm = TRUE),
          prop = dplyr::if_else(total_n > 0, n / total_n, 0)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(cluster_sizes, by = "cluster_id") %>%
        dplyr::mutate(
          cluster_size_kb = dplyr::coalesce(cluster_size_kb, 0),
          x_value = cluster_size_kb * prop,
          cluster_id_chr = cluster_id,
          cluster_id = factor(cluster_id, levels = cluster_order),
          hover_txt = paste0(
            "Cluster: ", cluster_id_chr,
            "<br>State: ", as.character(link_state),
            "<br>Count: ", n,
            "<br>Proportion: ", sprintf("%.1f%%", 100 * prop),
            "<br>Cluster size (kb): ", round(cluster_size_kb, 2),
            "<br>Segment size (kb): ", round(x_value, 2)
          )
        )
    })
    
    output$gwas_hit_match_summary_cluster_plot_XXXXX <- plotly::renderPlotly({
      df <- gwas_hit_match_summary_cluster_plot_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS × app cluster plot data available.")
      )
      
      color_map <- c(
        MATCH  = "#993404",
        MARKER = "#fe9929",
        NOLINK = "#fff3cd",
        NOHIT  = "#f2f2f2"
      )
      
      state_order <- c("MATCH", "MARKER", "NOLINK", "NOHIT")
      
      p <- plotly::plot_ly(source = session$ns("gwasxapp_cluster_plot"))
      
      for (st in state_order) {
        dsub <- df %>% dplyr::filter(as.character(link_state) == st)
        if (!nrow(dsub)) next
        
        p <- p %>%
          plotly::add_trace(
            data = dsub,
            y = ~cluster_id,
            x = ~x_value,
            type = "bar",
            orientation = "h",
            name = st,
            customdata = ~cluster_id_chr,
            hovertext = ~hover_txt,
            hoverinfo = "text",
            showlegend = FALSE,
            marker = list(
              color = unname(color_map[st]),
              line = list(color = "black", width = 0.4)
            )
          )
      }
      
      p %>%
        plotly::layout(
          barmode = "stack",
          clickmode = "event+select",
          xaxis = list(
            title = "Cluster size (kb)",
            automargin = TRUE
          ),
          yaxis = list(
            title = "Cluster",
            categoryorder = "array",
            categoryarray = levels(df$cluster_id),
            automargin = TRUE
          ),
          margin = list(l = 90, r = 20, t = 20, b = 60)
        )
    })
    
    output$gwas_hit_match_summary_cluster_legend <- renderUI({
      shiny::tags$div(
        style = "margin-top:8px; font-size:12px; color:#444;",
        shiny::HTML(
          paste0(
            "<b>Composition:</b> ",
            "<span style='display:inline-block; width:12px; height:12px; background:#993404; border:1px solid black; border-radius:50%; margin-left:8px;'></span> MATCH ",
            "<span style='display:inline-block; width:12px; height:12px; background:#fe9929; border:1px solid black; border-radius:50%; margin-left:12px;'></span> MARKER ",
            "<span style='display:inline-block; width:12px; height:12px; background:#fff3cd; border:1px solid black; border-radius:50%; margin-left:12px;'></span> NOLINK ",
            "<span style='display:inline-block; width:12px; height:12px; background:#f2f2f2; border:1px solid black; border-radius:50%; margin-left:12px;'></span> NOHIT"
          )
        )
      )
    })
    
    observeEvent(plotly::event_data("plotly_click", source = session$ns("gwasxapp_cluster_plot")), {
      ed <- plotly::event_data("plotly_click", source = session$ns("gwasxapp_cluster_plot"))
      req(ed)
      
      cid <- ed$customdata[[1]]
      req(!is.null(cid), nzchar(cid))
      
      session$sendCustomMessage("filter_gwasxapp_table_by_cluster", list(cluster_id = cid))
    })
    

  # ============================================================
    # ============================================================
    # GWAS hit × app summary: shared cached bundle
    # One computation -> one cache -> reused by DT + plot
    # ============================================================
    
    gwas_hit_match_bundle_cache <- reactiveVal(list(
      summary = tibble::tibble(),
      summary_with_gene = tibble::tibble(),
      plot_df = tibble::tibble()
    ))
    
    build_gwas_hit_match_summary_with_gene <- function(sm, gl, nhb) {
      if (!is.data.frame(sm) || !nrow(sm)) {
        return(tibble::tibble())
      }
      
      base <- sm %>%
        dplyr::mutate(
          .row_id = dplyr::row_number(),
          cluster_id = trimws(as.character(cluster_id)),
          source_app = tolower(trimws(as.character(source_app))),
          gwas_pos = suppressWarnings(as.numeric(gwas_pos)),
          link_state = trimws(as.character(link_state)),
          matched_ids = dplyr::coalesce(as.character(matched_ids), ""),
          nearest_gene = "",
          nearest_gene_dist_bp = NA_real_
        )
      
      # ------------------------------------------------------------
      # 1) NONSYN MATCH -> matched_ids -> hit_key -> gene
      # ------------------------------------------------------------
      nonsyn_match_rows <- base %>%
        dplyr::filter(source_app == "nonsyn", link_state == "MATCH")
      
      nonsyn_match_ann <- if (
        nrow(nonsyn_match_rows) &&
        is.data.frame(nhb) && nrow(nhb) &&
        all(c("cluster_id", "source_app", "hit_key", "gene") %in% names(nhb))
      ) {
        nhb2 <- nhb %>%
          dplyr::mutate(
            cluster_id = trimws(as.character(cluster_id)),
            source_app = tolower(trimws(as.character(source_app))),
            hit_key = trimws(as.character(hit_key)),
            gene = trimws(as.character(gene))
          ) %>%
          dplyr::filter(
            !is.na(cluster_id), nzchar(cluster_id),
            !is.na(source_app), nzchar(source_app),
            !is.na(hit_key), nzchar(hit_key),
            !is.na(gene), nzchar(gene),
            source_app == "nonsyn"
          ) %>%
          dplyr::distinct(cluster_id, source_app, hit_key, gene)
        
        nonsyn_match_rows %>%
          dplyr::select(.row_id, cluster_id, source_app, matched_ids) %>%
          tidyr::separate_rows(matched_ids, sep = "\\s*;\\s*") %>%
          dplyr::mutate(hit_key = trimws(matched_ids)) %>%
          dplyr::filter(!is.na(hit_key), nzchar(hit_key)) %>%
          dplyr::left_join(
            nhb2,
            by = c("cluster_id", "source_app", "hit_key"),
            relationship = "many-to-many"
          ) %>%
          dplyr::filter(!is.na(gene), nzchar(gene)) %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            nearest_gene = paste(sort(unique(gene)), collapse = "; "),
            nearest_gene_dist_bp = 0,
            .groups = "drop"
          )
      } else {
        tibble::tibble(
          .row_id = integer(),
          nearest_gene = character(),
          nearest_gene_dist_bp = numeric()
        )
      }
      
      # ------------------------------------------------------------
      # 2) MATCH / MARKER (except nonsyn MATCH)
      #    nearest gene by cluster_id + source_app
      # ------------------------------------------------------------
      other_rows <- base %>%
        dplyr::filter(
          link_state %in% c("MATCH", "MARKER"),
          !(source_app == "nonsyn" & link_state == "MATCH")
        )
      
      other_ann <- if (
        nrow(other_rows) &&
        is.data.frame(gl) && nrow(gl) &&
        all(c("cluster_id", "source_app", "gene") %in% names(gl))
      ) {
        gl2 <- gl %>%
          dplyr::mutate(
            cluster_id = trimws(as.character(cluster_id)),
            source_app = tolower(trimws(as.character(source_app))),
            gene = trimws(as.character(gene)),
            gene_mid = suppressWarnings(as.numeric(gene_mid))
          ) %>%
          dplyr::filter(
            !is.na(cluster_id), nzchar(cluster_id),
            !is.na(source_app), nzchar(source_app),
            !is.na(gene), nzchar(gene)
          ) %>%
          dplyr::distinct(cluster_id, source_app, gene, gene_mid)
        
        joined <- other_rows %>%
          dplyr::select(.row_id, cluster_id, source_app, gwas_pos) %>%
          dplyr::left_join(
            gl2,
            by = c("cluster_id", "source_app"),
            relationship = "many-to-many"
          )
        
        ann_dist <- joined %>%
          dplyr::filter(is.finite(gwas_pos), is.finite(gene_mid)) %>%
          dplyr::mutate(dist_bp = abs(gene_mid - gwas_pos)) %>%
          dplyr::arrange(.row_id, dist_bp, gene) %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            nearest_gene = dplyr::first(gene),
            nearest_gene_dist_bp = dplyr::first(dist_bp),
            .groups = "drop"
          )
        
        ann_nodist <- joined %>%
          dplyr::group_by(.row_id) %>%
          dplyr::summarise(
            any_dist = any(is.finite(gwas_pos) & is.finite(gene_mid)),
            genes_txt = paste(sort(unique(gene[!is.na(gene) & nzchar(gene)])), collapse = "; "),
            .groups = "drop"
          ) %>%
          dplyr::filter(!any_dist, nzchar(genes_txt)) %>%
          dplyr::transmute(
            .row_id,
            nearest_gene = genes_txt,
            nearest_gene_dist_bp = NA_real_
          )
        
        dplyr::bind_rows(ann_dist, ann_nodist)
      } else {
        tibble::tibble(
          .row_id = integer(),
          nearest_gene = character(),
          nearest_gene_dist_bp = numeric()
        )
      }
      
      ann <- dplyr::bind_rows(nonsyn_match_ann, other_ann) %>%
        dplyr::group_by(.row_id) %>%
        dplyr::summarise(
          nearest_gene = dplyr::first(nearest_gene),
          nearest_gene_dist_bp = dplyr::first(nearest_gene_dist_bp),
          .groups = "drop"
        )
      
      base %>%
        dplyr::select(-nearest_gene, -nearest_gene_dist_bp) %>%
        dplyr::left_join(ann, by = ".row_id") %>%
        dplyr::mutate(
          nearest_gene = dplyr::coalesce(nearest_gene, ""),
          nearest_gene_dist_bp = suppressWarnings(as.numeric(nearest_gene_dist_bp))
        ) %>%
        dplyr::select(-.row_id)
    }
    
    build_gwas_hit_match_summary_cluster_plot_df <- function(df, cl) {
      if (!is.data.frame(df) || !nrow(df)) {
        return(tibble::tibble())
      }
      
      if (!is.data.frame(cl) || !nrow(cl)) {
        return(tibble::tibble())
      }
      
      if (!all(c("cluster_id", "link_state") %in% names(df))) {
        return(tibble::tibble())
      }
      
      cl2 <- cl %>%
        dplyr::mutate(
          cluster_id = as.character(cluster_id)
        )
      
      if (!"cluster_size_kb" %in% names(cl2)) {
        start_col <- intersect(c("start", "pos_ini", "cluster_start"), names(cl2))[1]
        end_col   <- intersect(c("end", "pos_end", "cluster_end"), names(cl2))[1]
        size_col  <- intersect(c("size_kb", "cluster_size_kb"), names(cl2))[1]
        
        if (!is.na(size_col) && nzchar(size_col)) {
          cl2$cluster_size_kb <- suppressWarnings(as.numeric(cl2[[size_col]]))
        } else if (!is.na(start_col) && nzchar(start_col) && !is.na(end_col) && nzchar(end_col)) {
          cl2$cluster_size_kb <- (
            suppressWarnings(as.numeric(cl2[[end_col]])) -
              suppressWarnings(as.numeric(cl2[[start_col]])) + 1
          ) / 1000
        } else {
          cl2$cluster_size_kb <- NA_real_
        }
      }
      
      cluster_sizes <- cl2 %>%
        dplyr::transmute(
          cluster_id,
          cluster_size_kb = pmax(dplyr::coalesce(as.numeric(cluster_size_kb), 0), 0)
        ) %>%
        dplyr::distinct(cluster_id, .keep_all = TRUE)
      
      cluster_order <- df %>%
        dplyr::distinct(cluster_id) %>%
        dplyr::mutate(
          chr_num = extract_chr_num(cluster_id),
          cluster_num = suppressWarnings(as.integer(stringr::str_extract(cluster_id, "(?<=_)\\d+$")))
        ) %>%
        dplyr::mutate(
          chr_num = dplyr::coalesce(chr_num, 999),
          cluster_num = dplyr::coalesce(cluster_num, 999)
        ) %>%
        dplyr::arrange(chr_num, cluster_num, cluster_id) %>%
        dplyr::pull(cluster_id) %>%
        rev()
      
      df %>%
        dplyr::mutate(
          cluster_id = as.character(cluster_id),
          link_state = factor(link_state, levels = c("MATCH", "MARKER", "NOLINK", "NOHIT"))
        ) %>%
        dplyr::count(cluster_id, link_state, name = "n") %>%
        tidyr::complete(
          cluster_id,
          link_state = factor(
            c("MATCH", "MARKER", "NOLINK", "NOHIT"),
            levels = c("MATCH", "MARKER", "NOLINK", "NOHIT")
          ),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::mutate(
          total_n = sum(n, na.rm = TRUE),
          prop = dplyr::if_else(total_n > 0, n / total_n, 0)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(cluster_sizes, by = "cluster_id") %>%
        dplyr::mutate(
          cluster_size_kb = dplyr::coalesce(cluster_size_kb, 0),
          x_value = cluster_size_kb * prop,
          cluster_id_chr = cluster_id,
          cluster_id = factor(cluster_id, levels = cluster_order),
          hover_txt = paste0(
            "Cluster: ", cluster_id_chr,
            "<br>State: ", as.character(link_state),
            "<br>Count: ", n,
            "<br>Proportion: ", sprintf("%.1f%%", 100 * prop),
            "<br>Cluster size (kb): ", round(cluster_size_kb, 2),
            "<br>Segment size (kb): ", round(x_value, 2)
          )
        )
    }
    
    # ------------------------------------------------------------
    # Single shared computation
    # ------------------------------------------------------------
    observeEvent(
      list(
        gwas_hit_match_audit_df(),
        gene_bridge_by_cluster_app_r(),
        nonsyn_gene_hit_bridge_base_r(),
        cl_base_r()
      ),
      {
        au  <- gwas_hit_match_audit_df()
        gl  <- gene_bridge_by_cluster_app_r()
        nhb <- nonsyn_gene_hit_bridge_base_r()
        cl  <- cl_base_r()
        
        sm <- if (is.data.frame(au) && nrow(au)) {
          summarize_gwas_hit_app_audit(au)
        } else {
          tibble::tibble()
        }
        
        sm_with_gene <- build_gwas_hit_match_summary_with_gene(
          sm = sm,
          gl = gl,
          nhb = nhb
        )
        
        plot_df <- build_gwas_hit_match_summary_cluster_plot_df(
          df = sm_with_gene,
          cl = cl
        )
        
        gwas_hit_match_bundle_cache(list(
          summary = sm,
          summary_with_gene = sm_with_gene,
          plot_df = plot_df
        ))
      },
      ignoreInit = FALSE
    )
    
    gwas_hit_match_audit_summary_df <- reactive({
      gwas_hit_match_bundle_cache()$summary
    })
    
    gwas_hit_match_audit_summary_with_gene_df <- reactive({
      gwas_hit_match_bundle_cache()$summary_with_gene
    })
    
    gwas_hit_match_summary_cluster_plot_df <- reactive({
      gwas_hit_match_bundle_cache()$plot_df
    })
    
    # ============================================================
    # DT
    # ============================================================
    output$gwas_hit_match_summary_dt <- DT::renderDT({
      df <- gwas_hit_match_audit_summary_with_gene_df()
      
      shiny::validate(
        shiny::need(is.data.frame(df) && nrow(df) > 0, "No summarized GWAS hit × app table available.")
      )
      
      show_df <- df %>%
        dplyr::transmute(
          cluster_id,
          chr,
          gwas_hit,
          gwas_pos,
          source_app,
          link_state_raw = link_state,
          n_rows = dplyr::coalesce(as.integer(n_rows), 0L),
          n_verified_raw = dplyr::coalesce(as.integer(n_verified), 0L),
          matched_ids_raw = dplyr::coalesce(as.character(matched_ids), ""),
          matched_pos = dplyr::coalesce(as.character(matched_pos), ""),
          nearest_gene = dplyr::coalesce(as.character(nearest_gene), ""),
          nearest_gene_dist_kb = suppressWarnings(as.numeric(nearest_gene_dist_bp) / 1000)
        ) %>%
        dplyr::mutate(
          matched_ids = make_matched_ids_links(
            matched_ids = matched_ids_raw,
            source_app = source_app,
            link_state = link_state_raw,
            chr = chr,
            gwas_pos = gwas_pos
          ),
          nearest_gene = dplyr::if_else(
            nzchar(nearest_gene),
            make_genecards_links(nearest_gene),
            ""
          ),
          nearest_gene_dist_kb = dplyr::if_else(
            is.finite(nearest_gene_dist_kb),
            round(nearest_gene_dist_kb, 2),
            NA_real_
          ),
          link_state = dplyr::case_when(
            link_state_raw == "MATCH" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#993404; color:white; font-weight:700;'>MATCH</span>",
            link_state_raw == "MARKER" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#fe9929; color:black; font-weight:700;'>MARKER</span>",
            link_state_raw == "NOLINK" ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#fff3cd; color:black; font-weight:700;'>NOLINK</span>",
            TRUE ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#f2f2f2; color:black; font-weight:700;'>NOHIT</span>"
          ),
          n_verified = dplyr::case_when(
            n_verified_raw > 0 ~ paste0(
              "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#d9ead3; color:#1b5e20; font-weight:700;'>",
              n_verified_raw, "</span>"
            ),
            TRUE ~ "<span style='display:inline-block; padding:3px 10px; border-radius:999px; background:#f4cccc; color:#7f0000; font-weight:700;'>0</span>"
          )
        ) %>%
        dplyr::select(
          cluster_id, chr, gwas_hit, gwas_pos, source_app, link_state,
          n_rows, n_verified, matched_ids, matched_pos, nearest_gene, nearest_gene_dist_kb
        )
      
      DT::datatable(
        show_df,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        escape = FALSE,
        extensions = "Buttons",
        callback = DT::JS(sprintf("
      var tbl = table;
      window.gwasxappSummaryTable = tbl;
      window.gwasxappSummaryClusterCol = %d;
      
      tbl.columns.adjust();
      
      $(window).off('resize.gwasxappSummaryTable').on('resize.gwasxappSummaryTable', function() {
        tbl.columns.adjust();
      });
      
      $(document).off('shown.bs.tab.gwasxappSummaryTable').on(
        'shown.bs.tab.gwasxappSummaryTable',
        'button[data-bs-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]',
        function() {
          $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
        }
      );
    ", which(names(show_df) == "cluster_id") - 1L)),
        options = list(
          scrollX = TRUE,
          pageLength = 15,
          autoWidth = FALSE,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "copyHtml5",
              exportOptions = list(columns = ":visible")
            ),
            list(
              extend = "csvHtml5",
              exportOptions = list(columns = ":visible")
            ),
            list(
              extend = "excelHtml5",
              exportOptions = list(columns = ":visible")
            )
          )
        )
      ) %>%
        DT::formatRound("nearest_gene_dist_kb", 2)
    }, server = FALSE)
    
    # ============================================================
    # Plot
    # ============================================================
    output$gwas_hit_match_summary_cluster_plot <- plotly::renderPlotly({
      df <- gwas_hit_match_summary_cluster_plot_df()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No GWAS × app cluster plot data available.")
      )
      
      color_map <- c(
        MATCH  = "#993404",
        MARKER = "#fe9929",
        NOLINK = "#fff3cd",
        NOHIT  = "#f2f2f2"
      )
      
      state_order <- c("MATCH", "MARKER", "NOLINK", "NOHIT")
      
      p <- plotly::plot_ly(source = session$ns("gwasxapp_cluster_plot"))
      
      for (st in state_order) {
        dsub <- df %>% dplyr::filter(as.character(link_state) == st)
        if (!nrow(dsub)) next
        
        p <- p %>%
          plotly::add_trace(
            data = dsub,
            y = ~cluster_id,
            x = ~x_value,
            type = "bar",
            orientation = "h",
            name = st,
            customdata = ~cluster_id_chr,
            hovertext = ~hover_txt,
            hoverinfo = "text",
            showlegend = FALSE,
            marker = list(
              color = unname(color_map[st]),
              line = list(color = "black", width = 0.4)
            )
          )
      }
      
      p %>%
        plotly::layout(
          barmode = "stack",
          clickmode = "event+select",
          xaxis = list(
            title = "Cluster size (kb)",
            automargin = TRUE
          ),
          yaxis = list(
            title = "Cluster",
            categoryorder = "array",
            categoryarray = levels(df$cluster_id),
            automargin = TRUE
          ),
          margin = list(l = 90, r = 20, t = 20, b = 60)
        )
    })
    
    # ============================================================
    # Expose shared reactives
    # ============================================================
    invisible(list(
      gwas_hit_match_audit_summary_df = gwas_hit_match_audit_summary_df,
      gwas_hit_match_audit_summary_with_gene_df = gwas_hit_match_audit_summary_with_gene_df,
      gwas_hit_gene_support_long_df = gwas_hit_gene_support_long_df
    ))
  })
  # ============================================================
    
    
#    invisible(list(
#      gwas_hit_match_audit_summary_df = gwas_hit_match_audit_summary_df,
#      gwas_hit_match_audit_summary_with_gene_df = gwas_hit_match_audit_summary_with_gene_df,
#      gwas_hit_gene_support_long_df = gwas_hit_gene_support_long_df
#    ))
#  })
  
  # end module_server
  # ============================================================
  
}