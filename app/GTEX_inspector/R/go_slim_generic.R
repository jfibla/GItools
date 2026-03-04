# R/go_slim_generic.R
library(AnnotationDbi)
library(GO.db)
library(dplyr)

load_goslim_generic <- function(obo_path = file.path("R", "goslim_generic.obo")) {
  
  stopifnot(file.exists(obo_path))
  
  lines <- readLines(obo_path, warn = FALSE)
  ids <- grep("^id:\\s*GO:\\d+", lines, value = TRUE)
  ids <- unique(trimws(sub("^id:\\s*", "", ids)))
  
  tab <- AnnotationDbi::select(
    GO.db::GO.db,
    keys    = ids,
    keytype = "GOID",
    columns = c("GOID", "TERM", "ONTOLOGY")
  )
  
  tab <- tab %>%
    dplyr::filter(!is.na(GOID), !is.na(TERM), !is.na(ONTOLOGY)) %>%
    dplyr::transmute(
      GOID,
      ONTOLOGY,
      slim_term = TERM
    ) %>%
    dplyr::distinct(GOID, ONTOLOGY, slim_term)
  
  tab
}

GO_SLIM_GENERIC <- load_goslim_generic()
message("GO slim generic terms loaded: ", nrow(GO_SLIM_GENERIC))
