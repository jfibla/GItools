# go_map.R
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(GO.db)

detect_go_cols <- function(df) {
  nms <- names(df)
  nms_norm <- tolower(gsub("\\.", "_", nms))
  
  wanted <- c(
    GO_biological_process = "go_biological_process",
    GO_cellular_component = "go_cellular_component",
    GO_molecular_function = "go_molecular_function"
  )
  
  out <- setNames(rep(NA_character_, length(wanted)), names(wanted))
  for (k in names(wanted)) {
    idx <- which(nms_norm == wanted[[k]])
    if (length(idx)) out[[k]] <- nms[idx[1]]  # <- nombre REAL en df
  }
  
  out
}


norm_go_term <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(tolower(x))
  x
}

make_go_term_lut <- function() {
  
  ids <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
  
  # 1) Tabla principal GOID -> TERM/ONTOLOGY (SIN SYNONYM)
  tab <- AnnotationDbi::select(
    GO.db::GO.db,
    keys    = ids,
    keytype = "GOID",
    columns = c("GOID", "TERM", "ONTOLOGY")
  ) %>%
    dplyr::filter(!is.na(ONTOLOGY))
  
  # 2) Sinónimos desde GOSYNONYM (si existe)
  syn_tab <- NULL
  if (exists("GOSYNONYM", where = asNamespace("GO.db"), inherits = FALSE)) {
    syn_tab <- AnnotationDbi::toTable(GO.db::GOSYNONYM)
    
    # Normalmente viene como: go_id / synonym
    if ("go_id" %in% names(syn_tab))    syn_tab <- dplyr::rename(syn_tab, GOID = go_id)
    if ("synonym" %in% names(syn_tab))  syn_tab <- dplyr::rename(syn_tab, SYNONYM = synonym)
    
    # Por si acaso, filtra columnas esperadas
    syn_tab <- syn_tab %>%
      dplyr::select(dplyr::any_of(c("GOID", "SYNONYM"))) %>%
      dplyr::filter(!is.na(GOID), !is.na(SYNONYM), nzchar(SYNONYM))
  }
  
  # 3) Crear “keys” por TERM y por SYNONYM
  lut <- tab %>%
    dplyr::transmute(GOID, ONTOLOGY, key = TERM)
  
  if (!is.null(syn_tab) && nrow(syn_tab) > 0) {
    # Une ONTOLOGY desde 'tab' y añade SYNONYM como keys
    syn2 <- syn_tab %>%
      dplyr::inner_join(tab %>% dplyr::select(GOID, ONTOLOGY), by = "GOID") %>%
      dplyr::transmute(GOID, ONTOLOGY, key = SYNONYM)
    
    lut <- dplyr::bind_rows(lut, syn2)
  }
  
  lut %>%
    dplyr::mutate(key_norm = norm_go_term(key)) %>%
    dplyr::filter(nzchar(key_norm)) %>%
    dplyr::distinct(ONTOLOGY, key_norm, GOID)
}

GO_TERM_LUT <- make_go_term_lut()
message("GO term LUT loaded: ", nrow(GO_TERM_LUT), " rows")
