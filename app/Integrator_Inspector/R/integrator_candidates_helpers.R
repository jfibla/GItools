#R/integrator_candidates_helpers.R

prepare_candidates_canonical <- function(candidates_df) {
  if (!is.data.frame(candidates_df) || !nrow(candidates_df)) {
    return(tibble::tibble())
  }
  
  out <- candidates_df %>%
    dplyr::mutate(
      cluster_id = normalize_cluster_id(cluster_id),
      chr = as.character(chr),
      position = suppressWarnings(as.numeric(position)),
      rsid = trimws(as.character(rsid)),
      id_hit = trimws(as.character(id_hit)),
      classe = normalize_classe(classe),
      source_app_raw = if ("source_app" %in% names(.)) trimws(as.character(source_app)) else NA_character_,
      source_app_canonical = dplyr::case_when(
        classe == "catalog_hit" ~ "catalog",
        classe == "gtex_hit" ~ "gtex",
        classe == "nonsyn_hit" ~ "nonsyn",
        classe == "ewasdis_hit" ~ "ewasdis",
        classe == "ewastum_hit" ~ "ewastum",
        classe == "GWAS" ~ "gwas",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      !is.na(chr), nzchar(chr),
      is.finite(position),
      !is.na(classe), nzchar(classe)
    )
  
  gwas_df <- out %>%
    dplyr::filter(classe == "GWAS") %>%
    dplyr::group_by(cluster_id, chr, position, rsid, id_hit, classe, source_app_canonical) %>%
    dplyr::summarise(
      source_app_raw = paste(sort(unique(source_app_raw[!is.na(source_app_raw) & nzchar(source_app_raw)])), collapse = "; "),
      .groups = "drop"
    )
  
  app_hits_df <- out %>%
    dplyr::filter(classe != "GWAS") %>%
    dplyr::distinct(
      cluster_id, chr, position, rsid, id_hit, classe,
      source_app_raw, source_app_canonical,
      .keep_all = TRUE
    )
  
  dplyr::bind_rows(gwas_df, app_hits_df) %>%
    dplyr::arrange(cluster_id, chr, position, classe, id_hit)
}

cluster_app_presence_from_candidates <- function(candidates_df) {
  apps_ref <- c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum")
  
  cand <- prepare_candidates_canonical(candidates_df)
  
  if (!is.data.frame(cand) || !nrow(cand)) {
    return(tibble::tibble(
      cluster_id = character(),
      source_app = character(),
      cluster_has_app = logical()
    ))
  }
  
  cand %>%
    dplyr::filter(
      classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit"),
      source_app_canonical %in% apps_ref
    ) %>%
    dplyr::transmute(
      cluster_id = cluster_id,
      source_app = source_app_canonical
    ) %>%
    dplyr::distinct(cluster_id, source_app) %>%
    dplyr::mutate(cluster_has_app = TRUE)
}