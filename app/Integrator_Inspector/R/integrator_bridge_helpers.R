# R/integrator_bridge_helpers.R


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

chr_map_integrator <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT", "M", "MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

norm_int <- function(x) {
  suppressWarnings(as.integer(readr::parse_number(as.character(x))))
}

sanitize_bridge <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  if ("chr" %in% names(df))   df$chr   <- chr_map_integrator(df$chr)
  if ("start" %in% names(df)) df$start <- norm_int(df$start)
  if ("end" %in% names(df))   df$end   <- norm_int(df$end)
  
  if ("start" %in% names(df) && "end" %in% names(df)) {
    s0 <- df$start
    e0 <- df$end
    df$start <- pmin(s0, e0, na.rm = FALSE)
    df$end   <- pmax(s0, e0, na.rm = FALSE)
  }
  
  df
}

split_bridge_gene_tokens <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x)] <- ""
  
  vals <- unlist(strsplit(x, "\\s*;\\s*|\\s+-\\s+"))
  vals <- trimws(vals)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  vals <- vals[tolower(vals) != "numeric"]
  vals <- vals[!grepl("^[0-9]+$", vals)]
  
  unique(vals)
}

sanitize_bridge_genes <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  if (!"gene" %in% names(df)) return(df)
  
  df %>%
    dplyr::mutate(
      gene = as.character(gene),
      gene = trimws(gene)
    ) %>%
    tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s+-\\s+") %>%
    dplyr::mutate(
      gene = trimws(gene)
    ) %>%
    dplyr::filter(
      !is.na(gene),
      nzchar(gene),
      tolower(gene) != "numeric",
      !grepl("^[0-9]+$", gene)
    ) %>%
    dplyr::distinct()
}

all_gene_bridges <- reactive({
  idir <- selected_session_dir()
  
  dplyr::bind_rows(
    read_bridge_rds(file.path(idir, "catalog_gene_bridge.rds"), empty_gene_bridge()),
    read_bridge_rds(file.path(idir, "gtex_gene_bridge.rds"), empty_gene_bridge()),
    read_bridge_rds(file.path(idir, "nonsyn_gene_bridge.rds"), empty_gene_bridge()),
    read_bridge_rds(file.path(idir, "ewasdis_gene_bridge.rds"), empty_gene_bridge()),
    read_bridge_rds(file.path(idir, "ewastum_gene_bridge.rds"), empty_gene_bridge()),
    read_bridge_rds(file.path(idir, "ld_gene_bridge.rds"), empty_gene_bridge())
  ) %>%
    sanitize_bridge() %>%
    sanitize_bridge_genes() %>%
    dplyr::distinct()
})

as_gene_bridge <- function(df, source_app, evidence_type) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  out <- df %>%
    dplyr::mutate(
      gene = as.character(gene %||% NA),
      source_app = as.character(source_app),
      evidence_type = as.character(evidence_type),
      cluster_id = as.character(cluster_id %||% NA),
      chr = chr,
      start = start,
      end = end
    )
  
  # Assegurem columnes extres útils per GTEx sense trencar altres apps
  if (!"gtex_gene_name" %in% names(out)) {
    out$gtex_gene_name <- NA_character_
  } else {
    out$gtex_gene_name <- as.character(out$gtex_gene_name)
  }
  
  if (!"gtex_gene_id" %in% names(out)) {
    out$gtex_gene_id <- NA_character_
  } else {
    out$gtex_gene_id <- as.character(out$gtex_gene_id)
  }
  
  out <- out %>%
    sanitize_bridge() %>%
    dplyr::mutate(
      gene = trimws(as.character(gene))
    ) %>%
    tidyr::separate_rows(gene, sep = "\\s*;\\s*|\\s+-\\s+") %>%
    dplyr::mutate(
      gene = trimws(gene),
      gtex_gene_name = dplyr::coalesce(as.character(gtex_gene_name), ""),
      gtex_gene_id = dplyr::coalesce(as.character(gtex_gene_id), "")
    ) %>%
    dplyr::filter(
      !is.na(gene),
      nzchar(gene),
      tolower(gene) != "numeric",
      !grepl("^[0-9]+$", gene)
    ) %>%
    dplyr::distinct()
  
  out
}

as_term_bridge <- function(df, term_col = c("trait", "disease"), source_app, evidence_type) {
  term_col <- match.arg(term_col)
  
  out <- df %>%
    dplyr::transmute(
      term = as.character(.data[[term_col]]),
      term_type = term_col,
      source_app = as.character(source_app),
      evidence_type = as.character(evidence_type),
      cluster_id = as.character(cluster_id %||% NA),
      chr = chr,
      start = start,
      end = end
    ) %>%
    sanitize_bridge() %>%
    dplyr::filter(!is.na(term), nzchar(term))
  
  out
}