# ============================================================
# module_audit_plots.R
# Audit hierarchy module: selectors + filtered audit table
# with clickable links for genes and rsIDs
# ============================================================

priority_audit_hierarchy_module_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      column(
        3,
        selectInput(
          ns("cluster"),
          "Cluster",
          choices = character(0),
          selected = NULL
        )
      ),
      column(
        3,
        selectInput(
          ns("block"),
          "Block",
          choices = character(0),
          selected = NULL
        )
      ),
      column(
        3,
        selectInput(
          ns("gene"),
          "Gene",
          choices = character(0),
          selected = NULL
        )
      ),
      column(
        3,
        selectInput(
          ns("level"),
          "Audit level",
          choices = c("All", "Cluster", "Block", "Gene", "Hit"),
          selected = "All"
        )
      )
    ),
    
    fluidRow(
      column(
        12,
        div(
          style = "text-align: right; margin-top: 4px; margin-bottom: 10px;",
          actionButton(
            ns("reset_filters"),
            label = "Reset filters",
            icon = icon("rotate-left"),
            class = "btn-default"
          )
        )
      )
    ),
    
    fluidRow(
      column(
        12,
        div(
          class = "panelGrey",
          tags$div(
            class = "smallNote",
            HTML(
              paste0(
                "<b>Audit hierarchy explorer:</b> ",
                "Use the selectors to filter the audit table by cluster, block, gene, and audit level. ",
                "Selectors are dynamically coordinated, so you can start from any field and refine the view in any order."
              )
            )
          )
        )
      )
    ),
    
    tags$br(),
    
    fluidRow(
      column(
        12,
        div(
          class = "panelGrey padDT",
          DT::DTOutput(ns("audit_table"))
        )
      )
    )
  )
}

priority_audit_hierarchy_module_server <- function(id, data_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    split_semicolon_unique <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return(character(0))
      
      vals <- unlist(strsplit(paste(x, collapse = "; "), "\\s*;\\s*"))
      vals <- trimws(vals)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      unique(vals)
    }
    
    safe_chr <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      trimws(x)
    }
    
    `%||%` <- function(x, y) {
      if (is.null(x) || length(x) == 0) y else x
    }
    
    is_rsid <- function(x) {
      x <- safe_chr(x)
      grepl("^rs[0-9]+$", x, ignore.case = TRUE)
    }
    
    gene_to_genecards_link <- function(x) {
      x <- safe_chr(x)
      ifelse(
        nzchar(x),
        paste0(
          '<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',
          utils::URLencode(x, reserved = TRUE),
          '" target="_blank">', x, "</a>"
        ),
        ""
      )
    }
    
    rsid_to_app_link <- function(rsid, source_app = "") {
      rsid <- safe_chr(rsid)
      source_app <- tolower(safe_chr(source_app))
      
      out <- rsid
      
      idx <- is_rsid(rsid)
      if (!any(idx)) return(out)
      
      url <- rep("", length(rsid))
      
      url[idx & source_app == "gtex"] <- paste0(
        "https://gtexportal.org/home/snp/",
        rsid[idx & source_app == "gtex"]
      )
      
      url[idx & source_app == "catalog"] <- paste0(
        "https://www.ebi.ac.uk/gwas/variants/",
        rsid[idx & source_app == "catalog"]
      )
      
      url[idx & source_app == "nonsyn"] <- paste0(
        "https://www.ncbi.nlm.nih.gov/snp/",
        rsid[idx & source_app == "nonsyn"]
      )
      
      url[idx & !nzchar(url)] <- paste0(
        "https://www.ncbi.nlm.nih.gov/snp/",
        rsid[idx & !nzchar(url)]
      )
      
      out[idx] <- paste0(
        '<a href="', url[idx], '" target="_blank">', rsid[idx], "</a>"
      )
      
      out
    }
    
    infer_source_app_from_row <- function(df) {
      n <- nrow(df)
      out <- rep("", n)
      
      if ("source_app" %in% names(df)) {
        out <- safe_chr(df$source_app)
      }
      
      if ("score_app" %in% names(df)) {
        fill_idx <- !nzchar(out)
        out[fill_idx] <- safe_chr(df$score_app)[fill_idx]
      }
      
      out
    }
    
    round_numeric_df <- function(df, digits = 2) {
      num_cols <- vapply(df, is.numeric, logical(1))
      df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
      df
    }
    
    # ------------------------------------------------------------
    # Canonical input
    # Expected flexible columns if present:
    # entity_type / cluster_id / block_id / entity_id / gene / gwas_hit / source_app
    # ------------------------------------------------------------
    audit_df <- reactive({
      df <- data_reactive()
      
      validate(
        need(is.data.frame(df), "Audit module input is not a data frame.")
      )
      
      if (!nrow(df)) return(df)
      
      nm <- names(df)
      
      if (!"entity_type" %in% nm) {
        if ("level" %in% nm) {
          df$entity_type <- safe_chr(df$level)
        } else {
          df$entity_type <- ""
        }
      }
      
      if (!"cluster_id" %in% nm) df$cluster_id <- ""
      if (!"block_id" %in% nm) df$block_id <- ""
      if (!"entity_id" %in% nm) df$entity_id <- ""
      if (!"gene" %in% nm) df$gene <- ""
      if (!"gwas_hit" %in% nm) df$gwas_hit <- ""
      if (!"source_app" %in% nm) df$source_app <- ""
      
      df %>%
        dplyr::mutate(
          entity_type = safe_chr(.data$entity_type),
          cluster_id = safe_chr(.data$cluster_id),
          block_id = safe_chr(.data$block_id),
          entity_id = safe_chr(.data$entity_id),
          gene = safe_chr(.data$gene),
          gwas_hit = safe_chr(.data$gwas_hit),
          source_app = safe_chr(.data$source_app)
        )
    })
    
    # ------------------------------------------------------------
    # Reset filters
    # ------------------------------------------------------------
    observeEvent(input$reset_filters, {
      updateSelectInput(session, "cluster", selected = "")
      updateSelectInput(session, "block", selected = "")
      updateSelectInput(session, "gene", selected = "")
      updateSelectInput(session, "level", selected = "All")
    })
    
    # ------------------------------------------------------------
    # Generic filter builder
    # exclude = selector que estem actualitzant, per tal que
    # aquest selector es calculi a partir dels altres filtres actius
    # i es pugui triar qualsevol camp en qualsevol ordre.
    # ------------------------------------------------------------
    filtered_for_selector <- function(exclude = NULL) {
      df <- audit_df()
      
      if (!is.data.frame(df) || !nrow(df)) return(df)
      
      cid <- input$cluster %||% ""
      bid <- input$block %||% ""
      gid <- input$gene %||% ""
      lvl <- input$level %||% "All"
      
      out <- df
      
      if (!identical(exclude, "cluster") && nzchar(cid)) {
        out <- out %>% dplyr::filter(.data$cluster_id == cid)
      }
      
      if (!identical(exclude, "block") && nzchar(bid)) {
        out <- out %>% dplyr::filter(.data$block_id == bid)
      }
      
      if (!identical(exclude, "gene") && nzchar(gid)) {
        out <- out %>%
          dplyr::filter(
            .data$gene == gid |
              .data$entity_id == gid
          )
      }
      
      if (!identical(exclude, "level") && !identical(lvl, "All")) {
        out <- out %>%
          dplyr::filter(tolower(.data$entity_type) == tolower(lvl))
      }
      
      out
    }
    
    # ------------------------------------------------------------
    # Selector choices: coordinated in all directions
    # ------------------------------------------------------------
    observe({
      df <- filtered_for_selector(exclude = "cluster")
      req(is.data.frame(df))
      
      cluster_choices <- sort(unique(df$cluster_id[nzchar(df$cluster_id)]))
      selected_now <- isolate(input$cluster %||% "")
      selected_final <- if (selected_now %in% cluster_choices) selected_now else ""
      
      updateSelectInput(
        session,
        "cluster",
        choices = c("All" = "", cluster_choices),
        selected = selected_final
      )
    })
    
    observe({
      df <- filtered_for_selector(exclude = "block")
      req(is.data.frame(df))
      
      block_choices <- sort(unique(df$block_id[nzchar(df$block_id)]))
      selected_now <- isolate(input$block %||% "")
      selected_final <- if (selected_now %in% block_choices) selected_now else ""
      
      updateSelectInput(
        session,
        "block",
        choices = c("All" = "", block_choices),
        selected = selected_final
      )
    })
    
    observe({
      df <- filtered_for_selector(exclude = "gene")
      req(is.data.frame(df))
      
      gene_pool <- c(
        df$gene[nzchar(df$gene)],
        df$entity_id[tolower(df$entity_type) %in% c("gene") & nzchar(df$entity_id)]
      )
      
      gene_choices <- sort(unique(trimws(as.character(gene_pool))))
      gene_choices <- gene_choices[nzchar(gene_choices)]
      
      selected_now <- isolate(input$gene %||% "")
      selected_final <- if (selected_now %in% gene_choices) selected_now else ""
      
      updateSelectInput(
        session,
        "gene",
        choices = c("All" = "", gene_choices),
        selected = selected_final
      )
    })
    
    # ------------------------------------------------------------
    # Final filtered table
    # ------------------------------------------------------------
    filtered_audit_df <- reactive({
      df <- audit_df()
      
      req(is.data.frame(df))
      if (!nrow(df)) return(df)
      
      cid <- input$cluster %||% ""
      bid <- input$block %||% ""
      gid <- input$gene %||% ""
      lvl <- input$level %||% "All"
      
      out <- df
      
      if (nzchar(cid)) {
        out <- out %>% dplyr::filter(.data$cluster_id == cid)
      }
      
      if (nzchar(bid)) {
        out <- out %>% dplyr::filter(.data$block_id == bid)
      }
      
      if (nzchar(gid)) {
        out <- out %>%
          dplyr::filter(
            .data$gene == gid |
              .data$entity_id == gid
          )
      }
      
      if (!identical(lvl, "All")) {
        out <- out %>%
          dplyr::filter(tolower(.data$entity_type) == tolower(lvl))
      }
      
      out
    })
    
    # ------------------------------------------------------------
    # Output table
    # ------------------------------------------------------------
    output$audit_table <- DT::renderDT({
      df <- filtered_audit_df()
      
      validate(
        need(is.data.frame(df), "No audit data available."),
        need(nrow(df) > 0, "No rows match the current selector combination.")
      )
      
      show_df <- df
      show_df <- round_numeric_df(show_df, digits = 2)
      
      inferred_source_app <- infer_source_app_from_row(show_df)
      
      # --------------------------------------------------------
      # clickable gene links
      # --------------------------------------------------------
      if ("gene" %in% names(show_df)) {
        gene_idx <- nzchar(show_df$gene)
        show_df$gene[gene_idx] <- gene_to_genecards_link(show_df$gene[gene_idx])
      }
      
      if (all(c("entity_type", "entity_id") %in% names(show_df))) {
        gene_entity_idx <- tolower(safe_chr(show_df$entity_type)) == "gene" & nzchar(show_df$entity_id)
        show_df$entity_id[gene_entity_idx] <- gene_to_genecards_link(show_df$entity_id[gene_entity_idx])
      }
      
      # --------------------------------------------------------
      # clickable rsid links
      # --------------------------------------------------------
      if ("gwas_hit" %in% names(show_df)) {
        rs_idx <- is_rsid(show_df$gwas_hit)
        show_df$gwas_hit[rs_idx] <- rsid_to_app_link(
          show_df$gwas_hit[rs_idx],
          inferred_source_app[rs_idx]
        )
      }
      
      if ("entity_id" %in% names(show_df)) {
        rs_entity_idx <- is_rsid(show_df$entity_id)
        show_df$entity_id[rs_entity_idx] <- rsid_to_app_link(
          show_df$entity_id[rs_entity_idx],
          inferred_source_app[rs_entity_idx]
        )
      }
      
      DT::datatable(
        show_df,
        rownames = FALSE,
        escape = FALSE,
        filter = "top",
        extensions = "Buttons",
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copyHtml5", exportOptions = list(columns = ":all")),
            list(extend = "csvHtml5", exportOptions = list(columns = ":all")),
            list(extend = "excelHtml5", exportOptions = list(columns = ":all"))
          )
        )
      )
    }, server = FALSE)
  })
}