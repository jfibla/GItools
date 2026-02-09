# R/goslim_utils.R

# Lee GO slim generic desde un .obo:
# - Si es el archivo "goslim_generic.obo" típico, contiene solo los términos slim.
# - Si en cambio metes un go-basic.obo, también funciona si lleva líneas "subset: goslim_generic".
load_goslim_generic_ids <- function(obo_path) {
  stopifnot(file.exists(obo_path))
  x <- readLines(obo_path, warn = FALSE)
  
  is_term <- grepl("^\\[Term\\]$", x)
  term_starts <- which(is_term)
  
  # Caso 1: el .obo es un GO completo con subset tags
  has_subset <- any(grepl("^subset:\\s*goslim_generic\\b", x))
  if (has_subset) {
    # Extraer IDs de términos que tengan subset: goslim_generic
    ids <- character(0)
    for (i in seq_along(term_starts)) {
      s <- term_starts[i]
      e <- if (i < length(term_starts)) term_starts[i + 1] - 1 else length(x)
      block <- x[s:e]
      if (any(grepl("^subset:\\s*goslim_generic\\b", block))) {
        id_line <- block[grepl("^id:\\s*GO:\\d{7}\\b", block)]
        if (length(id_line)) {
          ids <- c(ids, sub("^id:\\s*", "", id_line[1]))
        }
      }
    }
    return(unique(ids))
  }
  
  # Caso 2: el .obo es "goslim_generic.obo" (solo slim terms)
  ids <- x[grepl("^id:\\s*GO:\\d{7}\\b", x)]
  ids <- sub("^id:\\s*", "", ids)
  unique(ids)
}

# Extrae GO IDs aunque vengan con texto, p.e. "GO:0008150; biological_process"
extract_go_ids <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  ids <- regmatches(x, gregexpr("GO:\\d{7}", x, perl = TRUE))
  # devuelve lista (una celda puede tener varios IDs)
  ids
}

# Nombre bonito del GO ID (si existe en GO.db)
go_term_name <- function(go_id) {
  if (!requireNamespace("GO.db", quietly = TRUE)) return(go_id)
  tt <- tryCatch(GO.db::GOTERM[[go_id]], error = function(e) NULL)
  if (is.null(tt)) return(go_id)
  nm <- tryCatch(GO.db::Term(tt), error = function(e) go_id)
  ifelse(is.na(nm) | !nzchar(nm), go_id, nm)
}

# Map GO IDs -> GO slim generic (puede devolver 0, 1 o varios slim IDs por GO ID)
map_go_to_goslim <- function(go_ids, ontology = c("BP","CC","MF"), goslim_ids) {
  ontology <- match.arg(ontology)
  stopifnot(length(goslim_ids) > 0)
  
  anc_db <- switch(
    ontology,
    BP = GO.db::GOBPANCESTOR,
    CC = GO.db::GOCCANCESTOR,
    MF = GO.db::GOMFANCESTOR
  )
  
  # mget devuelve lista de ancestros por GO
  anc <- mget(go_ids, anc_db, ifnotfound = list(character(0)))
  
  out <- lapply(seq_along(go_ids), function(i) {
    id <- go_ids[i]
    a  <- unique(c(id, anc[[i]]))
    slim <- intersect(a, goslim_ids)
    slim
  })
  out
}

# Shortener para etiquetas
shorten_go_label <- function(x, max = 45) {
  x <- as.character(x)
  x <- gsub("\\bregulation of\\b", "reg. of", x)
  x <- gsub("\\bpositive\\b", "pos.", x)
  x <- gsub("\\bnegative\\b", "neg.", x)
  x <- gsub("\\bcellular\\b",  "cell.", x)
  x <- gsub("\\bprocess\\b",   "proc.", x)
  x <- gsub("\\bcomponent\\b", "comp.", x)
  x <- gsub("\\bactivity\\b",  "act.", x)
  x <- stringr::str_squish(x)
  ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
}
