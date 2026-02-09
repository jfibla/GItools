# app.R — EWAS Inspector Disease (GWAS clusters + EWAS significant disease bins)
# Canonical path:
# Step 1) Load GWAS → Step 2) Build clusters → Step 3) Load EWAS resources → Step 4) Compute EWAS bins for ALL clusters
# Then: Manhattan + Cluster summary + Per-disease drill-down (violin/density/CpG delta plots) + UCSC custom tracks
# /Volumes/DISK1TB/Inspector_app_slaves_github/GItools/app/EWAS_disease/app.R

options(shiny.maxRequestSize = 1024*1024^2)

library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(purrr)
library(tibble)
library(data.table)
library(ggrepel)
library(shinyjs)
library(shinycssloaders)
library(htmltools)

# Gene model (hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(IRanges)

# (opcional però recomanat per labels)
library(AnnotationDbi)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ============================================================
# GItools portable config
# ============================================================

cfg_file <- normalizePath(file.path("..", "..", "config.R"), winslash = "/", mustWork = FALSE)
if (!file.exists(cfg_file)) {
  # fallback (si obres l'app des del root)
  cfg_file <- normalizePath("config.R", winslash = "/", mustWork = FALSE)
}
stopifnot(file.exists(cfg_file))
source(cfg_file, local = TRUE)

cfg <- gi_cfg()

options(gi_state_root = file.path(cfg$root, "app", "_state"))

gi_root        <- cfg$root
gi_shared_root <- cfg$shared
gi_res_root    <- cfg$resources
gi_pop_dir     <- cfg$pop_dir
#----------------------------------------------------------------------------
# --- EWAS reference paths (portable) ---
ref_paths <- list(
  slim_dir = file.path(gi_res_root, "EWAS_disease", "slim_disease_by_chr"),
  disease  = file.path(gi_res_root, "EWAS_disease", "ewas_detail_dis_genome.rds")
)

# (opcional però molt útil) avisos si falta res
if (!dir.exists(ref_paths$slim_dir)) {
  warning("EWAS slim_dir not found: ", ref_paths$slim_dir)
}
if (!file.exists(ref_paths$disease)) {
  warning("EWAS disease reference not found: ", ref_paths$disease)
}

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
gi_path <- function(...) {
  if (!nzchar(gi_res_root)) return("")
  file.path(gi_res_root, ...)
}
gi_file <- function(...) {
  p <- gi_path(...)
  if (nzchar(p) && file.exists(p)) p else ""
}
gi_dir <- function(...) {
  p <- gi_path(...)
  if (nzchar(p) && dir.exists(p)) p else ""
}


#----------------------------------------------------------------------------
# --- Canonical engines / modules from _shared ---
source(file.path(gi_shared_root, "gi_clusters_canonical.R"), local = TRUE)

ld_file <- file.path(gi_shared_root, "mod_ld.R")
stopifnot(file.exists(ld_file))
source(ld_file, local = TRUE)
stopifnot(exists("ld_module_ui"), exists("ld_module_server"))

# --- EWAS enrichment module (portable via _shared) ---
enrich_file <- file.path(gi_shared_root, "mod_ewas_enrichment.R")
stopifnot(file.exists(enrich_file))
source(enrich_file, local = TRUE)

stopifnot(exists("mod_ewas_enrichment_ui") || exists("mod_ewas_enrich_ui"))
# -----------------------------
# Safe load: coord_hg19 (optional)
# -----------------------------
coord_hg19 <- NULL
coord_hg19_path <- file.path("www", "coord_hg19.rds")
if (file.exists(coord_hg19_path)) {
  coord_hg19 <- readRDS(coord_hg19_path)
  coord_hg19 <- as.data.frame(coord_hg19, stringsAsFactors = FALSE)
}

# -----------------------------
# Helpers generals
# -----------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

chr_map_plink19 <- function(x){
  x <- toupper(as.character(x))
  x <- trimws(x)
  x <- sub("^CHR", "", x)
  x <- sub("^CHROM", "", x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT","M","MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}

chr_label_plink <- function(chr_num) {
  out <- as.character(chr_num)
  out[out == "23"] <- "X"
  out[out == "24"] <- "Y"
  out[out == "26"] <- "MT"
  out
}

norm_chr_generic <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M","MT","m","mt","MTDNA","mtdna")] <- "MT"
  toupper(x)
}

parse_p_robust <- function(p){
  if (is.numeric(p)) return(as.numeric(p))
  p_chr <- trimws(as.character(p))
  p_chr <- gsub("\\s+", "", p_chr)
  needs_swap <- !grepl("\\.", p_chr) & grepl(",", p_chr)
  p_chr[needs_swap] <- gsub(",", ".", p_chr[needs_swap], fixed = TRUE)
  suppressWarnings(as.numeric(p_chr))
}

plotly_message <- function(msg) {
  plotly::plot_ly(
    x = 0, y = 0,
    type = "scatter", mode = "markers",
    marker = list(opacity = 0)
  ) %>%
    plotly::layout(
      title = list(text = paste0("<b>", msg, "</b>"), x = 0.5, xanchor = "center"),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

################################################################################
# -----------------------------
# Reference hg38 chr lengths
# -----------------------------
.ref_hg38 <- tibble::tibble(
  chr = c(as.character(1:22), "X", "Y", "MT"),
  len = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468,
    156040895, 57227415, 16569
  )
) %>%
  dplyr::mutate(
    chrN = dplyr::case_when(
      chr == "X"  ~ 23L,
      chr == "Y"  ~ 24L,
      chr == "MT" ~ 25L,
      TRUE ~ suppressWarnings(as.integer(chr))
    ),
    chr_cum = cumsum(len) - len
  )

# -----------------------------
# Python extractor (slim subset)
# -----------------------------
extract_rows_python_slim <- function(txt_path, probes, out_path,
                                     report_every = 20000L,
                                     log_path = tempfile(fileext = ".log")) {
  stopifnot(length(probes) > 0)
  stopifnot(file.exists(txt_path))
  
  py <- tempfile(fileext = ".py")
  probes_file <- tempfile(fileext = ".txt")
  writeLines(unique(probes), probes_file)
  
  writeLines(c(
    "import sys, gzip",
    "",
    "txt_path = sys.argv[1]",
    "probes_path = sys.argv[2]",
    "out_path = sys.argv[3]",
    "report_every = int(sys.argv[4])",
    "",
    "def open_any(p):",
    "    if p.endswith('.gz'):",
    "        return gzip.open(p, 'rt', encoding='utf-8', errors='replace')",
    "    return open(p, 'r', encoding='utf-8', errors='replace')",
    "",
    "with open(probes_path, 'r', encoding='utf-8') as f:",
    "    probes = set([line.strip().replace('\\r','') for line in f if line.strip()])",
    "need = len(probes)",
    "found = 0",
    "scanned = 0",
    "",
    "with open_any(txt_path) as fin, open(out_path, 'w', encoding='utf-8') as fout:",
    "    header = fin.readline()",
    "    if not header:",
    "        raise SystemExit('Empty input')",
    "    fout.write(header)",
    "",
    "    for line in fin:",
    "        scanned += 1",
    "        if not line.startswith('cg'):",
    "            continue",
    "        t = line.find('\\t')",
    "        key = (line[:t] if t != -1 else line).strip().replace('\\r','')",
    "        if key in probes:",
    "            fout.write(line)",
    "            probes.remove(key)",
    "            found += 1",
    "            if found >= need:",
    "                break",
    "        if scanned % report_every == 0:",
    "            print(f\"scanned={scanned} found={found}/{need}\", file=sys.stderr)",
    "",
    "print(f\"DONE scanned={scanned} found={found}/{need}\", file=sys.stderr)"
  ), py)
  
  system2("python3",
          c(py, txt_path, probes_file, out_path, report_every),
          stdout = NULL, stderr = log_path)
  
  list(out_path = out_path, log_path = log_path)
}

# -----------------------------
# EWAS helpers (bins + diseases)
# -----------------------------
as_chr_disease <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^CHR", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "M"] <- "MT"
  paste0("chr", x)
}

chr_for_ref <- function(chr_disease) {
  x <- toupper(gsub("^CHR", "", chr_disease, ignore.case = TRUE))
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x
}

probes_in_interval <- function(coord, chr_disease, st, en) {
  chr_disease <- as_chr_disease(chr_disease)
  coord %>%
    dplyr::filter(chr == chr_disease, pos >= st, pos <= en) %>%
    dplyr::pull(probe) %>%
    unique()
}


sig_bins_for_disease <- function(
    sub_full, coord, meta, disease_name,
    chr_disease, st, en,
    bin_size = 200000L,
    min_cpg  = 10L,
    min_n    = 10L,
    test_m   = c("wilcox","ttest"),
    use_adj_as_ctl = TRUE
){
  test_m <- match.arg(test_m)
  
  meta_c <- meta[disease == disease_name &
                   sample_type %in% c("disease tissue","control","adjacent normal"),
                 .(sample_id, disease, sample_type)]
  meta_c <- unique(meta_c, by = "sample_id")
  
  dis_ids <- meta_c[sample_type == "disease tissue", sample_id]
  ctl_ids <- meta_c[sample_type == "control",       sample_id]
  adj_ids <- meta_c[sample_type == "adjacent normal", sample_id]
  
  if (length(ctl_ids) == 0 && isTRUE(use_adj_as_ctl) && length(adj_ids) > 0) {
    ctl_ids <- adj_ids
  }
  
  if (length(dis_ids) < min_n || length(ctl_ids) < min_n) return(NULL)
  
  sel <- unique(c(dis_ids, ctl_ids))
  
  # ------------------------------------------------------------
  # READ subset robustly (probe_id column may be misnamed sample_id)
  # ------------------------------------------------------------
  if (is.character(sub_full) && length(sub_full) == 1L) {
    
    # read header only
    hdr <- names(data.table::fread(sub_full, nrows = 0, check.names = FALSE))
    if (!length(hdr)) return(NULL)
    
    idcol <- hdr[1]
    # if the file still uses "sample_id" for probes, treat it as probe_id logically
    if (identical(idcol, "sample_id")) idcol <- "probe_id"
    
    # GSM columns available in file
    avail <- setdiff(hdr, hdr[1])
    sel2  <- intersect(sel, avail)
    if (length(sel2) < (2L * min_n)) return(NULL)  # rough guard; avoids useless reads
    
    cols_to_read <- c(hdr[1], sel2)   # IMPORTANT: the physical column name in file is hdr[1]
    dt <- data.table::fread(
      sub_full,
      select = cols_to_read,
      na.strings = c("NA","<NA>","NaN",""),
      check.names = FALSE
    )
    
    # normalize id column name in memory
    if (names(dt)[1] == "sample_id") data.table::setnames(dt, 1, "probe_id")
    
  } else if (data.table::is.data.table(sub_full) || is.data.frame(sub_full)) {
    
    dt <- data.table::as.data.table(sub_full)
    if (names(dt)[1] == "sample_id") data.table::setnames(dt, 1, "probe_id")
    
    # keep only probe_id + selected samples that exist
    sel2 <- intersect(sel, names(dt))
    if (length(sel2) < (2L * min_n)) return(NULL)
    dt <- dt[, c("probe_id", sel2), with = FALSE]
    
  } else {
    stop("sub_full must be a file path or a data.frame/data.table.")
  }
  
  # Now we use probe_id consistently
  if (!"probe_id" %in% names(dt)) return(NULL)
  
  # keep only CpG-like probes
  beta_dt <- dt[grepl("^cg\\d{8}$", probe_id)]
  if (nrow(beta_dt) == 0) return(NULL)
  
  # coerce numeric for selected columns (only those present)
  sel_present <- setdiff(names(beta_dt), "probe_id")
  beta_dt[, (sel_present) := lapply(.SD, as.numeric), .SDcols = sel_present]
  
  coord2 <- coord[chr == chr_disease & pos >= st & pos <= en]
  pos <- coord2$pos[match(beta_dt$probe_id, coord2$probe)]
  keep <- is.finite(pos)
  
  beta_dt <- beta_dt[keep]
  pos <- as.integer(pos[keep])
  if (nrow(beta_dt) == 0) return(NULL)
  
  bin_size <- as.integer(bin_size)
  st0 <- as.integer(st)
  
  bin_start <- st0 + ((pos - st0) %/% bin_size) * bin_size
  bin_mid   <- bin_start + bin_size/2
  
  idx_list <- split(seq_along(bin_mid), bin_mid)
  idx_list <- idx_list[vapply(idx_list, length, 1L) >= min_cpg]
  if (!length(idx_list)) return(NULL)
  
  out_list <- lapply(names(idx_list), function(bm) {
    ii <- idx_list[[bm]]
    
    # matrix over selected sample columns
    mat <- as.matrix(beta_dt[ii, ..sel_present])
    mb <- colMeans(mat, na.rm = TRUE)
    
    x <- mb[names(mb) %in% dis_ids]
    y <- mb[names(mb) %in% ctl_ids]
    if (length(x) < min_n || length(y) < min_n) return(NULL)
    
    p <- if (test_m == "ttest") {
      tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
    } else {
      tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    }
    
    data.table::data.table(
      disease   = disease_name,
      chr       = chr_disease,
      bin_mid   = as.numeric(bm),
      bin_start = as.integer(as.numeric(bm) - bin_size/2),
      bin_end   = as.integer(as.numeric(bm) + bin_size/2),
      n_cpg     = length(ii),
      n_dis     = length(x),
      n_ctl     = length(y),
      p         = as.numeric(p),
      delta     = mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE)
    )
  })
  
  res <- data.table::rbindlist(out_list, fill = TRUE)
  if (!nrow(res)) return(NULL)
  
  res[, padj := p.adjust(p, method = "BH")]
  res
}


compute_bins_anydisease_for_cluster <- function(
    sub_full, coord, meta,
    chr_disease, st, en,
    alpha = 0.05,
    bin_size = 200000L,
    min_n = 10L,
    min_cpg = 10L,
    test_m = c("wilcox","ttest"),
    use_adj_as_ctl = TRUE
) {
  
  test_m <- match.arg(test_m)
  
  # --- Robust: allow sub_full as a path or as a data.table/data.frame ---
  if (is.character(sub_full) && length(sub_full) == 1) {
    if (!file.exists(sub_full)) stop("sub_full path does not exist: ", sub_full)
    
    # Ensure normalized header (just in case older cached subset exists)
    hdr <- names(data.table::fread(sub_full, nrows = 0))
    if (length(hdr) && identical(hdr[1], "sample_id")) {
      dt0 <- data.table::fread(sub_full)
      data.table::setDT(dt0)
      data.table::setnames(dt0, 1, "probe_id")
      data.table::fwrite(dt0, sub_full, sep = "\t")
    }
    
  } else if (is.data.frame(sub_full) || data.table::is.data.table(sub_full)) {
    
    sub_full <- data.table::as.data.table(sub_full)
    if (names(sub_full)[1] == "sample_id") {
      data.table::setnames(sub_full, 1, "probe_id")
    }
    
  } else {
    stop("sub_full must be a file path or a data.frame/data.table.")
  }
  
  diseases <- sort(unique(as.character(meta$disease)))
  diseases <- diseases[!is.na(diseases) & nzchar(diseases)]
  
  all_res <- lapply(diseases, function(dd) {
    sig_bins_for_disease(
      sub_full = sub_full,
      coord = coord,
      meta = meta,
      disease_name = dd,
      chr_disease = chr_disease,
      st = st, en = en,
      bin_size = bin_size,
      min_cpg = min_cpg,
      min_n = min_n,
      test_m = test_m,
      use_adj_as_ctl = use_adj_as_ctl
    )
  })
  
  dt <- data.table::rbindlist(all_res, fill = TRUE)
  if (!nrow(dt)) return(list(detail = dt, bins = data.table::data.table()))
  
  dt_sig <- dt[is.finite(padj) & padj < alpha]
  if (!nrow(dt_sig)) return(list(detail = dt, bins = data.table::data.table()))
  
  bins <- dt_sig[, .(
    n_diseases   = data.table::uniqueN(disease),
    diseases     = paste(sort(unique(disease)), collapse = ", "),
    best_padj    = min(padj, na.rm = TRUE),
    best_disease = disease[which.min(padj)],
    best_delta   = delta[which.min(padj)]
  ), by = .(chr, bin_start, bin_end, bin_mid)]
  
  data.table::setorder(bins, best_padj)
  list(detail = dt_sig, bins = bins)
}

bins_to_manhattan_df <- function(bins_dt) {
  if (!is.data.frame(bins_dt) || !nrow(bins_dt)) return(tibble())
  
  ref <- .ref_hg38
  bins_dt2 <- as.data.table(bins_dt)
  
  bins_dt2[, CHRref := chr_for_ref(chr)]
  bins_dt2 <- bins_dt2[CHRref %in% ref$chr]
  
  ref2 <- as.data.table(ref)[, .(CHRref = as.character(chr), chr_cum = as.numeric(chr_cum))]
  bins_dt2 <- merge(bins_dt2, ref2, by = "CHRref", all.x = TRUE)
  bins_dt2 <- bins_dt2[is.finite(chr_cum)]
  
  bins_dt2[, BPcum := as.numeric(bin_mid) + chr_cum]
  bins_dt2[, y := -log10(as.numeric(best_padj))]
  if (!"cluster_id" %in% names(bins_dt2)) bins_dt2[, cluster_id := NA_character_]
  as_tibble(bins_dt2)
}

# ----------------------------
# Relayout helpers
# ----------------------------
get_relayout_xrange <- function(d) {
  if (is.null(d) || length(d) == 0) return(NULL)
  auto <- d[["xaxis.autorange"]] %||% d[["xaxis.autorange[0]"]] %||% NULL
  if (!is.null(auto) && isTRUE(as.logical(auto))) return(NULL)
  
  x0 <- d[["xaxis.range[0]"]] %||% d[["xaxis.range0"]] %||% NULL
  x1 <- d[["xaxis.range[1]"]] %||% d[["xaxis.range1"]] %||% NULL
  
  if ((is.null(x0) || is.null(x1)) &&
      !is.null(d[["xaxis"]]) &&
      !is.null(d[["xaxis"]][["range"]]) &&
      length(d[["xaxis"]][["range"]]) >= 2) {
    x0 <- d[["xaxis"]][["range"]][[1]]
    x1 <- d[["xaxis"]][["range"]][[2]]
  }
  
  if (is.null(x0) || is.null(x1)) return(NULL)
  
  x0 <- suppressWarnings(as.numeric(x0))
  x1 <- suppressWarnings(as.numeric(x1))
  if (!is.finite(x0) || !is.finite(x1)) return(NULL)
  
  if (x1 < x0) { tmp <- x0; x0 <- x1; x1 <- tmp }
  c(x0, x1)
}

cumrange_to_ucsc_region <- function(x0, x1, ref = .ref_hg38) {
  x0 <- suppressWarnings(as.numeric(x0))
  x1 <- suppressWarnings(as.numeric(x1))
  if (!is.finite(x0) || !is.finite(x1)) return(NULL)
  if (x1 < x0) { tmp <- x0; x0 <- x1; x1 <- tmp }
  
  ref2 <- ref %>%
    dplyr::transmute(
      chr     = toupper(as.character(chr)),
      chr_cum = suppressWarnings(as.numeric(chr_cum)),
      len     = suppressWarnings(as.numeric(len))
    ) %>%
    dplyr::filter(is.finite(chr_cum), is.finite(len), len > 0) %>%
    dplyr::arrange(chr_cum)
  
  if (!nrow(ref2)) return(NULL)
  
  GENOME_END <- max(ref2$chr_cum + ref2$len, na.rm = TRUE)
  x0 <- max(0, min(x0, GENOME_END))
  x1 <- max(0, min(x1, GENOME_END))
  
  xm <- (x0 + x1) / 2
  idx <- which(xm >= ref2$chr_cum & xm <= (ref2$chr_cum + ref2$len))
  if (!length(idx)) return(NULL)
  i <- idx[1]
  
  chr <- ref2$chr[i]
  st  <- floor(x0 - ref2$chr_cum[i])
  en  <- ceiling(x1 - ref2$chr_cum[i])
  
  st <- max(1L, as.integer(st))
  en <- min(as.integer(ref2$len[i]), as.integer(en))
  if (en < st) { tmp <- st; st <- en; en <- tmp }
  
  paste0("chr", chr, ":", format(st, scientific = FALSE), "-", format(en, scientific = FALSE))
}


# -----------------------------
# UI
# -----------------------------
ui <- navbarPage(
  title = div(style="font-weight:700; font-size:22px; color:#1A4E8A;", HTML("🧫 EWASdis Inspector (Disease)")),
  id = "topnav",
  header = tagList(
    tags$head(
      tags$script(HTML("
(function(){
  function pushQS(){
    try {
      if (window.Shiny && Shiny.setInputValue) {
        Shiny.setInputValue('gi_qs', window.location.search || '', {priority:'event'});
        Shiny.setInputValue('gi_href', window.location.href || '', {priority:'event'});
        return true;
      }
    } catch(e) {}
    return false;
  }

  // Try now + retries (Shiny may not be ready at DOMContentLoaded behind proxy)
  var tries = 0;
  var iv = setInterval(function(){
    tries++;
    if (pushQS() || tries > 50) clearInterval(iv); // ~5s max
  }, 100);

  window.addEventListener('pageshow', function(){ pushQS(); });
  document.addEventListener('visibilitychange', function(){
    if (!document.hidden) pushQS();
  });
})();
  "))
    ),
    tags$style(HTML("
      .diseasePick{
        color:#1A4E8A;
        text-decoration:underline;
        cursor:pointer;
      }
      #tbl_bin_diseases table.dataTable tbody tr.selected,
      #tbl_bin_diseases table.dataTable tbody td.selected {
        background-color: transparent !important;
      }
      #tbl_bin_diseases table.dataTable tbody td:focus {
        outline: none !important;
      }
      .panelGrey{
        background:#f2f3f5;
        padding:12px;
        border-radius:10px;
        border:1px solid rgba(0,0,0,0.06);
        margin:12px 0;
      }
    ")),
    
    tags$head(
      tags$style(HTML("
/* --- Popup diseases --- */
#giDiseasePopup{
  position: absolute;
  z-index: 99999;
  background: #ffffff;
  border: 1px solid rgba(0,0,0,0.15);
  border-radius: 12px;
  box-shadow: 0 10px 30px rgba(0,0,0,0.15);
  padding: 10px 12px;
  min-width: 260px;
  max-width: 360px;
  max-height: 320px;
  overflow: auto;
  font-size: 13px;
}

#giDiseasePopup .gi-title{
  font-weight: 700;
  color: #1A4E8A;
  margin-bottom: 8px;
}

#giDiseasePopup .gi-close{
  float: right;
  cursor: pointer;
  font-weight: 700;
  padding: 0 6px;
  border-radius: 8px;
}

#giDiseasePopup .gi-close:hover{
  background: #f3f6fb;
}

#giDiseasePopup .gi-item{
  margin: 4px 0;
}

#giDiseasePopup .diseasePick{
  cursor: pointer;
  color: #1A4E8A;
}

#giDiseasePopup .diseasePick:hover{
  text-decoration: underline;
}

/* --- Clickable pill in table cell --- */
.gi-openList{
  display: inline-block;
  padding: 4px 10px;
  border-radius: 10px;
  background: #f3f6fb;
  color: #1A4E8A;
  font-weight: 600;
  cursor: pointer;
}
.gi-openList:hover{
  background: #e7eefb;
}
"))
    ),
    shinyjs::useShinyjs()
  ),
  
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Disease analysis</span>"),
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        h3("Step 1 · Load GWAS hits", style="color:#1A4E8A; font-weight:700;"),
        fluidRow(
          column(6, actionButton("info_00", "ℹ️ file format")),
          column(6, actionButton(
            "reset_case",
            "Reset",
            icon = icon("rotate-left"),
            class = "btn-warning"
          ))),
        fileInput("gwas_file", "GWAS p-value table (TSV/CSV)", accept = c(".tsv",".txt",".csv")),
        checkboxInput("gwas_header", "First row contains column names", TRUE),
        radioButtons("gwas_sep", "Separator",
                     c("Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"),
                     selected = "\t"),
        uiOutput("gwas_p_selector"),
        
        tags$hr(),
        
        h3("Step 2 · Clustering GWAS hits", style="color:#1A4E8A; font-weight:700;"),
        actionButton("info_01", "ℹ️ Clustering method", class = "btn btn-default"),
        
        radioButtons(
          "cluster_method", "Clustering method:",
          choices  = c("By hit intervals (thr + flank → merge)" = "window",
                       "By hit density (min_logp + min_hits)"   = "hits"),
          selected = "window"
        ),
        
        conditionalPanel(
          condition = "input.cluster_method == 'window'",
          sliderInput("pthr", "-log10(P) threshold", min = 2, max = 20, value = 5, step = 0.5),
          numericInput("flank", "Flank (+/- bp)", value = 10000, min = 0, max = 10000000, step = 1000)
        ),
        
        conditionalPanel(
          condition = "input.cluster_method == 'hits'",
          
          radioButtons(
            "hits_mode", "Hits mode:",
            choices = c("1Mb hit-span (consecutive hits within win_bp)" = "span1mb",
                        "Tiled windows (non-overlapping)"              = "tiled",
                        "Sliding windows (overlapping)"                = "sliding"),
            selected = "span1mb"
          ),
          
          sliderInput("min_logp", "-log10(P) threshold (hit significance)",
                      min = 2, max = 20, value = 5, step = 0.1),
          
          numericInput("min_hits", "Minimum GWAS hits per cluster/window",
                       value = 3, min = 1, max = 1000, step = 1),
          
          numericInput("win_bp", "Window size (bp)", value = 1e6, min = 1e4, max = 5e7, step = 1e4),
          
          conditionalPanel(
            condition = "input.hits_mode == 'sliding'",
            numericInput("step_bp", "Step (bp)", value = 1e5, min = 1e3, max = 5e7, step = 1e3)
          ),
          
          conditionalPanel(
            condition = "input.hits_mode != 'sliding'",
            helpText("For 'tiled' and '1Mb hit-span', step is implicit.")
          )
        ),
        
        actionButton(
          "build_ranges", "➊ intervals → merge → clusters",
          style="background-color:#ffdd57; color:black; font-weight:bold;"
        ),
        verbatimTextOutput("ranges_preview"),
        
        tags$hr(),
        
        h3("Step 3 · EWAS resources", style="color:#1A4E8A; font-weight:700;"),
        textInput("slim_dir", "slim_disease_by_chr directory",
                  value = ref_paths$slim_dir),
        textInput("meta_rds", "meta_slim_disease.rds path",
                  value = "www/meta_slim_disease.rds"),
        
        tags$hr(),
        
        h3("Step 4 · EWAS bin-disease hits", style="color:#1A4E8A; font-weight:700;"),
        actionButton("run_ewas_all","➋ Compute EWAS bins for ALL clusters",
                     style="background-color:#ffdd57; color:black; font-weight:bold;"),
        
        tags$hr(),
        downloadButton("dl_candidates_zip", "⬇️ Download catalog candidates (ZIP)"),
        tags$hr(),
        strong("Log"),
        verbatimTextOutput("run_log", placeholder = TRUE)
      ),
      
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📈 Manhattan</span>"),
            value = "tab_manhattan",
            
            h4("Top: GWAS p-values · Bottom: EWAS significant bins (BH)"),
            
            div(class="panelGrey",
                withSpinner(plotlyOutput("manhattan_combo", height = "720px"))
            ),
            
            helpText("Click diamond symbols to inspect bins; zoom to update UCSC region links."),
            
            fluidRow(
              column(6, uiOutput("debug_ucsc_state")),
              column(2, h4("UCSC links to region:")),
              column(4,
                     div(
                       class = "panelGrey",
                       style = "margin-top:8px; width:100%;",
                       uiOutput("ucsc_link_gwas"),
                       uiOutput("ucsc_link_cpg"),
                       uiOutput("ucsc_link_ewasbins"),
                       helpText("Click a link to open a UCSC window.")
                     )
              )
            ),
            
            tags$hr(),
            h4("Clusters summary"),
            div(class="panelGrey",
                withSpinner(DTOutput("cluster_dt"))
            ),
            helpText("Then go to 'Cluster details' for per-disease drill-down.")
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔎 Cluster details</span>"),
            value = "tab_cluster_detail",
            
            sidebarLayout(
              sidebarPanel(
                width = 2,
                h4("Per-disease bin test controls"),
                numericInput("alpha_any", "Alpha (BH) per disease", value = 0.05, min = 1e-6, max = 0.5, step = 0.01),
                numericInput("min_n_any", "Min samples/group (per disease)", value = 10, min = 3, step = 1),
                numericInput("min_cpg_any", "Min CpGs/bin", value = 10, min = 1, step = 1),
                numericInput("bin_size", "Bin size (bp)", value = 5000, min = 1000, step = 500),
                selectInput("bin_test_any", "Test per disease", choices = c("Wilcoxon"="wilcox","t-test"="ttest"), selected = "wilcox"),
                checkboxInput("use_adj_as_ctl", "If no 'control', use 'adjacent normal' as control", value = TRUE)
              ),
              
              mainPanel(
                width = 10,
                h4("Summary for selected cluster"),
                helpText("Select a cluster on 'Clusters summary' and click the button below."),
                actionButton("build_summary_sel_cluster", "Build summary (selected cluster)", icon = icon("table")),
                tags$hr(),
                
                div(class="panelGrey",
                    uiOutput("bin_disease_cluster_filter_ui"),
                    withSpinner(DTOutput("tbl_bin_diseases"))
                ),
                helpText("Click a disease name to render plots."),
                tags$hr(),
                
                fluidRow(
                  column(
                    8,
                    h4("Mean beta per bin (disease vs control)"),
                    div(class="panelGrey",
                        withSpinner(plotlyOutput("p_violin_sel", height = "853px"))
                    )
                  ),
                  column(
                    4,
                    h4("Density plot (disease vs control)"),
                    div(class="panelGrey",
                        withSpinner(plotOutput("p_beta_dist", height = "388px"))
                    )
                  ),
                  column(
                    4,
                    h4("Delta-beta histogram"),
                    div(class="panelGrey",
                        withSpinner(plotOutput("p_hist", height = "388px"))
                    )
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧬 GWAS hits</span>"),
            value = "tab_hits",
            div(class="panelGrey",
                withSpinner(DTOutput("hits_tbl"))
            )
          ),
          
          tabPanel(
            title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧮 EWAS hits</span>"),
            value = "tab_ewas_hits",
            actionButton("build_ewas_hits", "Build EWAS hits (CpG table)"),
            div(class="panelGrey",
                fluidRow(
                  column(3, uiOutput("ewas_cluster_filter_ui")),
                  column(3, uiOutput("ewas_disease_filter_ui")),
                  column(3, uiOutput("ewas_probe_filter_ui")),
                  column(3, textInput("ewas_free_search", "Search (any)", value = ""))
                ),
                DTOutput("tbl_ewas_hits")
            ),
            helpText("Click on 'disease', 'probe' or 'hyper_hypo' to expand plots.")
          )
        )
      )
    )
  ),
  
  # ✅ TAB PRINCIPAL 2 (mateix nivell)
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧩 LD</span>"),
    # opcional: si vols que LD tingui sidebar + main:
    # sidebarLayout(sidebarPanel(...), mainPanel(ld_module_ui("ld")))
    ld_module_ui("ld")
  ),
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>✳️ Enrichment</span>"),
    mod_ewas_enrich_ui("enrich")
  )
)


# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  # -----------------------------
  # Locate _shared (Hub + standalone)
  # -----------------------------
  # Portable: _shared comes from GItools/config.R
  SHARED <- gi_shared_root  # (already defined from gi_cfg())
  stopifnot(dir.exists(SHARED))
  
  # -----------------------------
  # Shared core (state/helpers used across GItools)
  # - defines things like: .ref_hg38, make_mode_thr_tag(), etc.
  # -----------------------------
  if (file.exists(file.path(SHARED, "gi_state.R"))) {
    source(file.path(SHARED, "gi_state.R"), local = TRUE)
  }
  
  # -----------------------------
  # Shared reactive holders (filled by Hub OR local upload)
  # -----------------------------
  gwas_shared_rv <- shiny::reactiveVal(NULL)
  
  # -----------------------------
  # Canonical cluster engine (ALWAYS available)
  # - standalone: user clicks Build ranges
  # - hub: we inject clusters_shared into it
  # -----------------------------
  source(file.path(SHARED, "gi_clusters_canonical.R"), local = TRUE)
  
  # NOTE: gwas_df is defined later in this server() in the original file.
  # We therefore create a thin proxy here, and re-bind it right after gwas_df exists.
  gwas_df_proxy <- shiny::reactive({
    # 1) If later code defines gwas_df() (upload/local), prefer it.
    if (exists("gwas_df", inherits = TRUE) && is.function(gwas_df)) {
      df2 <- tryCatch(gwas_df(), error = function(e) NULL)
      if (is.data.frame(df2) && nrow(df2) > 0) return(df2)
    }
    # 2) Otherwise, use Hub-shared GWAS if present
    df <- gwas_shared_rv()
    if (is.data.frame(df) && nrow(df) > 0) return(df)
    NULL
  })
  
  gi_cl <- gi_clusters_canonical_init(
    session, input, output,
    gwas_df = gwas_df_proxy,
    build_btn_id   = "build_ranges",
    clusters_dt_id = "cluster_dt",
    hits_rows_id   = "hits_tbl_rows_selected",
    app_count_col  = "n_EWAS_bins"
  )
  
  # expose canonical names (used by deeplinks + rest of app)
  clusters_cur  <- gi_cl$clusters_cur
  intervals_raw <- gi_cl$intervals_raw
  
  # -----------------------------
  # Slave canonical sync (Hub -> this app)
  # -----------------------------
  source(file.path(SHARED, "gi_slave_canonical.R"), local = TRUE)
  gi_sync <- tryCatch(gi_slave_canonical_init(session), error = function(e) NULL)
  
  # robust accessors (avoid '$ on closure' at startup)
  get_sync_fn <- function(obj, nm) {
    if (is.null(obj)) return(NULL)
    if (is.list(obj) && !is.null(obj[[nm]])) return(obj[[nm]])
    if (is.environment(obj) && exists(nm, envir = obj, inherits = FALSE)) return(get(nm, envir = obj))
    NULL
  }
  
  gwas_shared_fn     <- get_sync_fn(gi_sync, "gwas_shared")
  clusters_shared_fn <- get_sync_fn(gi_sync, "clusters_shared")
  
  if (!is.function(gwas_shared_fn))     gwas_shared_fn     <- shiny::reactiveVal(NULL)
  if (!is.function(clusters_shared_fn)) clusters_shared_fn <- shiny::reactiveVal(NULL)
  
  # Apply Hub GWAS (if present)
  shiny::observeEvent(gwas_shared_fn(), {
    df <- gwas_shared_fn()
    req(is.data.frame(df), nrow(df) > 0)
    gwas_shared_rv(df)
    cat("[SLAVE] applied GWAS rows=", nrow(df), "\n")
  }, ignoreInit = FALSE)
  
  # Apply Hub CLUSTERS (if present)
  shiny::observeEvent(clusters_shared_fn(), {
    cl <- clusters_shared_fn()
    req(is.data.frame(cl), nrow(cl) > 0)
    
    if (exists("standardize_cluster_ids", mode = "function")) {
      cl <- standardize_cluster_ids(cl)
    }
    
    # push into canonical engine
    try(gi_cl$clusters_cur(cl), silent = TRUE)
    try(gi_cl$intervals_raw(cl |> dplyr::transmute(chr = .data$chr, start = .data$start, end = .data$end, label = .data$cluster_id)), silent = TRUE)
    
    cat("[SLAVE] applied CLUSTERS rows=", nrow(cl), "\n")
  }, ignoreInit = FALSE)
  # -----------------------------
  # Deeplinks
  # Disabled in EWAS canonical apps because the deeplink helper registers observers
  # that assume a different object shape and can crash startup (closure-not-subsettable).
  # Hub ↔ slave communication is handled exclusively by gi_slave_canonical.R.
  # -----------------------------
  
  # -----------------------------
  # HUB -> EWAS : GWAS sync feedback
  # -----------------------------
  observeEvent(gi_sync$gwas_shared(), {
    
    id <- paste0("hub_sync_gwas_", as.integer(Sys.time()))
    showNotification("🔄 HUB: sincronitzant GWAS…", type = "message", duration = NULL, id = id)
    
    df <- gi_sync$gwas_shared()
    req(is.data.frame(df), nrow(df) > 0)
    
    # guarda GWAS rebut (igual que GTEx/NonSyn)
    gwas_shared_rv(df)
    
    cat("[SLAVE] applied GWAS rows=", nrow(df), "\n")
    
    removeNotification(id)
    showNotification(sprintf("✅ HUB: GWAS sincronitzat (%d files).", nrow(df)),
                     type = "message", duration = 2)
    
  }, ignoreInit = FALSE)
  
  
  ################## end slave step
  
  # UCSC region state (MUST be defined before observers that use it)
  ucsc_region  <- reactiveVal(NULL)
  combo_region <- reactiveVal(NULL)
  
  # -----------------------------
  # State
  # -----------------------------
  rv <- reactiveValues(
    log_txt = "",
    workdir = NULL,
    
    meta = NULL,
    coord = NULL,
    
    clusters = NULL,
    
    summary_cluster = NULL,
    summary_cluster_ids = character(0),
    
    ewas_sub_full = NULL,
    ewas_chr = NULL,
    ewas_st = NULL,
    ewas_en = NULL,
    
    ewas_bins_all   = NULL,
    ewas_detail_all = NULL,
    
    ewas_bins   = NULL,
    ewas_detail = NULL,
    
    disease_sel = NULL,
    
    res = NULL,
    beta_mat = NULL,
    meta_c = NULL,
    
    ewas_subset_cache = new.env(parent = emptyenv())
  )
  
  # ---- subset cache (per-session, NOT reactive) ----
  ewas_subset_cache <- new.env(parent = emptyenv())
  
  # (opcional) netejar al final de sessió
  session$onSessionEnded(function() {
    rm(list = ls(envir = ewas_subset_cache, all.names = TRUE), envir = ewas_subset_cache)
  })
  
  rv$workdir <- file.path(tempdir(), "ewas_inspector_disease")
  dir.create(isolate(rv$workdir), showWarnings = FALSE, recursive = TRUE)
  
  append_log <- function(...) {
    rv$log_txt <- paste(c(rv$log_txt, paste(..., collapse = " ")), collapse = "\n")
  }
  output$run_log <- renderText(rv$log_txt)
  
  notify_err <- function(msg) {
    append_log("[ERROR]", msg)
    showNotification(msg, type = "error", duration = 8)
  }
  notify_info <- function(msg) {
    append_log("[INFO]", msg)
    showNotification(msg, type = "message", duration = 4)
  }
  
  
  # -----------------------------
  # Clickable disease tags JS
  # -----------------------------
  observe({
    shinyjs::runjs("
      $(document).off('click', '.diseasePick');
      $(document).on('click', '.diseasePick', function(){
        var dd = $(this).data('disease');
        var cl = $(this).data('cluster');
        Shiny.setInputValue('disease_pick', {disease: dd, cluster_id: cl, nonce: Math.random()}, {priority: 'event'});
      });
    ")
  })
  
  # -----------------------------
  # Load coord_hg38 at startup (www/coord_hg38.rds)
  # -----------------------------
  observeEvent(TRUE, {
    coord_path <- file.path("www", "coord_hg38.rds")
    
    if (!file.exists(coord_path)) {
      append_log("[coord] ERROR: file not found: ", coord_path)
      return()
    }
    
    coord <- readRDS(coord_path)
    coord <- data.table::as.data.table(coord)
    
    stopifnot(all(c("probe","chr","pos") %in% names(coord)))
    
    coord[, probe := as.character(probe)]
    coord[, chr   := as_chr_disease(chr)]
    coord[, pos   := as.integer(pos)]
    
    rv$coord <- coord
    append_log("[coord] loaded at startup from www/coord_hg38.rds | rows=", nrow(coord))
  }, once = TRUE)
  
  # -----------------------------
  # Step 3 resources: meta
  # -----------------------------
  
  observeEvent(input$meta_rds, {
    
    req(input$meta_rds)
    
    # Soporta: fileInput() (lista con $datapath) o textInput() (string)
    meta_path <- if (is.list(input$meta_rds) && !is.null(input$meta_rds$datapath)) {
      input$meta_rds$datapath
    } else {
      as.character(input$meta_rds)
    }
    
    req(nzchar(meta_path))
    req(file.exists(meta_path))
    
    meta <- readRDS(meta_path)
    data.table::setDT(meta)
    
    need_cols <- c("sample_id", "disease", "sample_type")
    miss <- setdiff(need_cols, names(meta))
    validate(shiny::need(length(miss) == 0,
                         paste("meta is missing columns:", paste(miss, collapse = ", "))))
    
    meta[, sample_id   := as.character(sample_id)]
    meta[, disease     := as.character(disease)]
    meta[, sample_type := as.character(sample_type)]
    
    rv$meta <- meta
    
    append_log("[meta] loaded file=", basename(meta_path),
               " rows=", nrow(meta),
               " diseases=", data.table::uniqueN(meta$disease))
    
  }, ignoreInit = FALSE)
  
  
  # -----------------------------
  # Step 1: GWAS reader + p selector
  # -----------------------------
  gwas_preview <- reactive({
    req(input$gwas_file)
    tryCatch({
      readr::read_delim(
        input$gwas_file$datapath,
        delim = input$gwas_sep %||% "\t",
        col_names = isTRUE(input$gwas_header),
        show_col_types = FALSE,
        progress = FALSE,
        n_max = 200
      )
    }, error = function(e) NULL)
  })
  
  output$gwas_p_selector <- renderUI({
    p <- gwas_preview()
    if (is.null(p) || !nrow(p)) return(helpText("Could not read GWAS (preview)."))
    
    cols <- names(p)
    
    # Si ve del master, sovint ja tens Pval/logp; fem selector tolerant
    guess <- intersect(cols, c("Pval","PVAL","pval","pvalue","P","p","P_VALUE","p_value"))
    
    selectInput(
      "gwas_col_p",
      "Select p-value column",
      choices  = cols,
      selected = if (length(guess)) guess[1] else cols[1]
    )
  })
  
  gwas_df <- reactive({
    req(input$gwas_file)
    
    df <- readr::read_delim(
      input$gwas_file$datapath,
      delim = input$gwas_sep %||% "\t",
      col_names = isTRUE(input$gwas_header),
      show_col_types = FALSE,
      progress = FALSE
    )
    validate(need(is.data.frame(df) && nrow(df) > 0, "Empty GWAS file."))
    
    chr_col <- pick_col(df, c("CHR","chr","chrom","CHROM","chromosome"))
    bp_col  <- pick_col(df, c("BP","bp","POS","pos","position","POSITION"))
    snp_col <- pick_col(df, c("SNP","snp","rsid","RSID","marker","ID","id"))
    
    p_col_in <- input$gwas_col_p
    if (is.null(p_col_in) || length(p_col_in) != 1 || !nzchar(p_col_in) || !(p_col_in %in% names(df))) {
      p_col_in <- pick_col(df, c("P","p","PVAL","pval","P_VALUE","p_value"))
    }
    
    validate(
      need(!is.null(chr_col), "GWAS: missing chromosome column."),
      need(!is.null(bp_col),  "GWAS: missing position column."),
      need(!is.null(p_col_in), "GWAS: missing p-value column.")
    )
    
    BP   <- if (is.numeric(df[[bp_col]])) df[[bp_col]] else suppressWarnings(readr::parse_number(as.character(df[[bp_col]])))
    Pval <- parse_p_robust(df[[p_col_in]])
    CHR  <- chr_map_plink19(df[[chr_col]])
    snp  <- if (!is.null(snp_col)) as.character(df[[snp_col]]) else paste0("chr", norm_chr_generic(df[[chr_col]]), ":", BP)
    
    out <- tibble(
      CHR = CHR,
      BP  = as.numeric(BP),
      snp = snp,
      Pval = as.numeric(Pval)
    ) %>%
      dplyr::filter(is.finite(CHR), is.finite(BP), is.finite(Pval), Pval > 0) %>%
      dplyr::mutate(logp = -log10(Pval))
    
    validate(need(nrow(out) > 0, "GWAS parsed but resulted in 0 valid rows (check columns/format)."))
    out
  })
  
  # -----------------------------
  # Canonical cluster init (local build_ranges)
  # -----------------------------
  
  # keep placeholders in sync with canonical engine (local clustering)
  shiny::observeEvent(gi_cl$clusters_cur(), {
    cl <- gi_cl$clusters_cur()
    if (is.data.frame(cl) && nrow(cl)) {
      clusters_cur(cl)
      if (exists("rv", inherits = TRUE)) {
        try({ if (inherits(rv, "reactivevalues")) rv$clusters <- cl }, silent = TRUE)
      }
    }
  }, ignoreInit = FALSE)
  
  shiny::observeEvent(gi_cl$intervals_raw(), {
    iv <- gi_cl$intervals_raw()
    if (is.data.frame(iv)) intervals_raw(iv)
  }, ignoreInit = FALSE)
  
  dfp_manhattan <- reactive({
    df <- gwas_df()
    # filter pval < 0.0 to a easely plot
    df <- df %>% dplyr::filter(Pval < 0.05)
    ref <- .ref_hg38
    
    df %>%
      dplyr::mutate(
        CHR = suppressWarnings(as.integer(CHR)),
        BP  = as.numeric(BP)
      ) %>%
      dplyr::filter(is.finite(CHR), is.finite(BP)) %>%
      dplyr::left_join(
        ref %>% dplyr::transmute(CHR = as.integer(chrN), chr_cum = chr_cum),
        by = "CHR"
      ) %>%
      dplyr::filter(is.finite(chr_cum)) %>%
      dplyr::mutate(BPcum = BP + chr_cum)
  })
  
  axis_df <- reactive({
    ref <- .ref_hg38
    ref %>%
      dplyr::transmute(
        chrN   = as.integer(chrN),
        center = chr_cum + (len / 2)
      )
  })
  
  hits_df <- reactive({
    df <- gwas_df()
    if ((input$cluster_method %||% "window") == "window") {
      req(input$pthr)
      df %>% dplyr::filter(logp >= input$pthr) %>%
        dplyr::arrange(dplyr::desc(logp)) %>%
        dplyr::select(CHR, BP, snp, p = Pval, logp)
    } else {
      req(input$min_logp)
      df %>% dplyr::filter(logp >= input$min_logp) %>%
        dplyr::arrange(dplyr::desc(logp)) %>%
        dplyr::select(CHR, BP, snp, p = Pval, logp)
    }
  })
  
  output$hits_tbl <- renderDT({
    
    h <- tryCatch(hits_df(), error = function(e) NULL)
    
    if (is.null(h) || !nrow(h)) {
      return(datatable(
        data.frame(Message = "No hits above threshold."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # -----------------------------
    # Assignar cluster a cada SNP
    # -----------------------------
    if (is.data.frame(rv$clusters) && nrow(rv$clusters) > 0) {
      
      # data.table no-equi join (ràpid i robust)
      h_dt  <- as.data.table(h)
      cl_dt <- as.data.table(rv$clusters)
      
      # assegura noms esperats (hits: CHR/BP; clusters: chr/start/end/cluster_id)
      if ("CHR" %in% names(h_dt)) setnames(h_dt, "CHR", "chr")
      if ("BP"  %in% names(h_dt)) setnames(h_dt, "BP",  "pos")
      
      # si clusters té CHR/BP en comptes de chr/start/end, adapta-ho
      if ("CHR" %in% names(cl_dt)) setnames(cl_dt, "CHR", "chr")
      if ("start_bp" %in% names(cl_dt)) setnames(cl_dt, "start_bp", "start")
      if ("end_bp"   %in% names(cl_dt)) setnames(cl_dt, "end_bp",   "end")
      
      # evita problemes de tipus
      h_dt[, chr := as.integer(chr)]
      h_dt[, pos := as.integer(pos)]
      cl_dt[, chr   := as.integer(chr)]
      cl_dt[, start := as.integer(start)]
      cl_dt[, end   := as.integer(end)]
      
      # crea columna cluster_id si no existeix
      if (!"cluster_id" %in% names(h_dt)) h_dt[, cluster_id := NA_character_]
      
      # join per interval
      # (cada SNP hauria de caure en 0 o 1 cluster; si en caigués en >1, data.table replicaria files)
      h_dt <- cl_dt[h_dt, on = .(chr, start <= pos, end >= pos)]
      # ara h_dt té columnes de cl_dt prefixades (i.) o directes segons el cas; unifica:
      if ("cluster_id" %in% names(h_dt)) {
        # OK
      } else if ("i.cluster_id" %in% names(h_dt)) {
        setnames(h_dt, "i.cluster_id", "cluster_id")
      }
      
      h <- as.data.frame(h_dt)
    } else {
      # si encara no hi ha clusters, igualment crea la columna perquè l'UI sigui consistent (opcional)
      if (!"cluster_id" %in% names(h)) h$cluster_id <- NA_character_
    }
    
    # -----------------------------
    # Neteja columnes + ordre
    # -----------------------------
    
    # Treu columnes no desitjades si existeixen
    drop_cols <- c("n_snps", "cluster_n", "cluster_chr", "center", "i.cluster_id", "top_snp", "top_logp","start","end","cluster_chr_n")
    h <- h[, setdiff(names(h), drop_cols), drop = FALSE]
    
    # Assegura que p i logp vagin just després de snp (i cluster_id després)
    first_cols <- intersect(c("cluster_id","snp", "p", "logp","cluster_size_kb","n_EWAS_bins"), names(h))
    h <- h[, c(first_cols, setdiff(names(h), first_cols)), drop = FALSE]
    
    # -----------------------------
    # Link dbSNP a 'snp'
    # -----------------------------
    if ("snp" %in% names(h)) {
      snp_txt <- as.character(h$snp)
      h$snp <- sprintf(
        '<a href="https://www.ncbi.nlm.nih.gov/snp/%s" target="_blank" rel="noopener">%s</a>',
        snp_txt, snp_txt
      )
    }
    
    # (opcional) mou cluster_id al costat de snp si vols
    if (all(c("cluster_id","snp") %in% names(h))) {
      h <- h[, c("snp", "cluster_id", setdiff(names(h), c("snp","cluster_id")))]
    }
    
    # -----------------------------
    # DT + formats
    # -----------------------------
    dt <- datatable(
      h,
      selection = "multiple",
      rownames  = FALSE,
      escape    = FALSE,  # IMPORTANT per renderitzar l'<a href=...>
      extensions = "Buttons",
      options   = list(
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE,
        columnDefs = list(
          # Format condicional per 'p': científic si < 1e-3, sinó 3 decimals
          list(
            targets = which(names(h) == "p") - 1,
            render  = JS(
              "function(data, type, row, meta){",
              "  if(data === null || data === undefined || data === '') return data;",
              "  var x = parseFloat(data);",
              "  if(isNaN(x)) return data;",
              "  if(type === 'display'){",
              "    if(x < 1e-3) return x.toExponential(2);",
              "    return x.toFixed(3);",
              "  }",
              "  return x;",  # per ordenar/filtrar
              "}"
            )
          )
        )
      )
    )
    
    # logp a 3 decimals (si existeix)
    if ("logp" %in% names(h)) {
      dt <- DT::formatRound(dt, columns = "logp", digits = 3)
    }
    
    dt
  })
  # -----------------------------
  # Step 2: build clusters (from gi_clusters_canonical_init())
  # -----------------------------
  
  # -----------------------------
  
  # Cluster summary table
  # -----------------------------
  clusters_show <- reactive({
    cl <- rv$clusters
    req(is.data.frame(cl), nrow(cl) > 0)
    
    cl <- as.data.frame(cl)
    
    if (!"cluster_size_kb" %in% names(cl) && all(c("start","end") %in% names(cl))) {
      cl$cluster_size_kb <- round((as.numeric(cl$end) - as.numeric(cl$start)) / 1000, 2)
    }
    if (!"n_snps" %in% names(cl)) cl$n_snps <- NA_integer_
    if (!"center" %in% names(cl) && all(c("start","end") %in% names(cl))) {
      cl$center <- as.integer(round((as.numeric(cl$start) + as.numeric(cl$end)) / 2))
    }
    if (!"top_snp" %in% names(cl)) cl$top_snp <- NA_character_
    if (!"top_logp" %in% names(cl)) cl$top_logp <- NA_real_
    if (!"n_EWAS_bins" %in% names(cl)) cl$n_EWAS_bins <- 0L
    
    dplyr::as_tibble(cl) %>%
      dplyr::transmute(
        cluster_id,
        chr,
        start,
        end,
        cluster_size_kb = as.numeric(cluster_size_kb),
        n_snps          = as.integer(n_snps),
        center          = as.integer(center),
        top_snp,
        top_logp        = as.numeric(top_logp),
        n_EWAS_bins     = as.integer(n_EWAS_bins)
      )
  })
  
  output$cluster_dt <- renderDT({
    df_show <- tryCatch(clusters_show(), error = function(e) NULL)
    if (is.null(df_show) || !nrow(df_show)) {
      return(datatable(
        data.frame(Message="No clusters yet. Run Step 2 (➊)."),
        options = list(dom="t"), rownames = FALSE
      ))
    }
    
    datatable(
      dplyr::select(df_show, -center),  # omit center
      selection = "single",
      rownames = FALSE,
      extensions = "Buttons",
      options    = list(
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel", "pdf", "print"),
        pageLength = 10,
        scrollX    = TRUE
      )
    ) %>%
      formatRound(columns = c("top_logp","cluster_size_kb"), digits = 2)
  })
  
  # -----------------------------
  # Subset extraction helper (per cluster)
  # -----------------------------
  extract_subset_if_needed <- function(chr_disease, st, en) {
    
    slim_chr_path <- file.path(input$slim_dir, paste0("slim_ds_", chr_disease, ".tsv.gz"))
    if (!file.exists(slim_chr_path)) {
      append_log("[EWAS] missing slim file: ", slim_chr_path)
      return(NULL)
    }
    
    probes_int <- probes_in_interval(rv$coord, chr_disease, st, en)
    if (!length(probes_int)) return(NULL)
    
    sub_full <- file.path(rv$workdir, paste0("subset__", chr_disease, "_", st, "_", en, ".tsv"))
    
    # If already exists, still ensure it is normalized
    if (file.exists(sub_full)) {
      
      # --- normalize header if needed (sample_id -> probe_id) ---
      hdr <- names(data.table::fread(sub_full, nrows = 0))
      if (length(hdr) && identical(hdr[1], "sample_id")) {
        dt0 <- data.table::fread(sub_full)
        data.table::setDT(dt0)
        data.table::setnames(dt0, 1, "probe_id")
        data.table::fwrite(dt0, sub_full, sep = "\t")
      }
      
      return(sub_full)
    }
    
    ex <- extract_rows_python_slim(
      txt_path = slim_chr_path,
      probes   = probes_int,
      out_path = sub_full,
      report_every = 200000L
    )
    
    if (file.exists(ex$log_path)) {
      append_log(paste(readLines(ex$log_path), collapse = "\n"))
    }
    
    if (!file.exists(sub_full)) return(NULL)
    
    # --- normalize header: first column is probe id but file calls it "sample_id" ---
    hdr <- names(data.table::fread(sub_full, nrows = 0))
    if (length(hdr) && identical(hdr[1], "sample_id")) {
      dt0 <- data.table::fread(sub_full)
      data.table::setDT(dt0)
      data.table::setnames(dt0, 1, "probe_id")
      data.table::fwrite(dt0, sub_full, sep = "\t")
    }
    
    sub_full
  }
  
  
  # ---- Override GWAS reactive: prefer shared GWAS if present ----
  gwas_df_local <- gwas_df
  gwas_df <- reactive({
    df <- gwas_shared_rv()
    if (is.data.frame(df) && nrow(df) > 0) return(df)
    gwas_df_local()
  })
  
  
  # -----------------------------
  # Step 4: EWAS bins for ALL clusters
  # -----------------------------
  observeEvent(input$run_ewas_all, {
    
    append_log("[DEBUG] run_ewas_all clicked @", format(Sys.time(), "%H:%M:%S"))
    
    tryCatch({
      
      cl <- rv$clusters
      if (!is.data.frame(cl) || !nrow(cl)) { notify_err("No clusters yet. Run Step 2 (➊)."); return() }
      if (is.null(rv$meta)  || !nrow(rv$meta))  { notify_err("meta_slim_ds not loaded."); return() }
      if (is.null(rv$coord) || !nrow(rv$coord)) { notify_err("coord_hg38 not loaded."); return() }
      if (!dir.exists(input$slim_dir)) { notify_err("slim_ds_by_chr directory not found."); return() }
      if (!nzchar(Sys.which("python3"))) { notify_err("python3 not found on PATH."); return() }
      
      rv$ewas_bins_all   <- data.table()
      rv$ewas_detail_all <- data.table()
      
      alpha   <- suppressWarnings(as.numeric(input$alpha_any))
      bin_sz  <- suppressWarnings(as.integer(input$bin_size))
      min_n   <- suppressWarnings(as.integer(input$min_n_any))
      min_cpg <- suppressWarnings(as.integer(input$min_cpg_any))
      test_m  <- as.character(input$bin_test_any)
      use_adj <- isTRUE(input$use_adj_as_ctl)
      
      validate(
        need(!is.na(alpha)   && alpha > 0 && alpha < 1, "alpha_any must be a number in (0,1)."),
        need(!is.na(bin_sz)  && bin_sz > 0,            "bin_size must be a positive integer."),
        need(!is.na(min_n)   && min_n >= 2,            "min_n_any must be an integer >= 2."),
        need(!is.na(min_cpg) && min_cpg >= 1,          "min_cpg_any must be an integer >= 1."),
        need(nzchar(test_m),                         "bin_test_any is empty.")
      )
      
      append_log(
        "[EWAS-ALL] start clusters=", nrow(cl),
        " alpha=", alpha, " bin_size=", bin_sz,
        " min_n=", min_n, " min_cpg=", min_cpg, " test=", test_m
      )
      
      withProgress(message = "Computing EWAS bins for ALL clusters…", value = 0, {
        
        ncl <- nrow(cl)
        inc <- 1 / max(1, ncl)
        
        all_bins_list   <- vector("list", ncl)
        all_detail_list <- vector("list", ncl)
        
        for (i in seq_len(ncl)) {
          
          cl0 <- cl[i, , drop = FALSE]
          
          # --- robust column access ---
          cid <- if ("cluster_id" %in% names(cl0)) as.character(cl0$cluster_id[1]) else paste0("cluster_", i)
          
          # chr can be numeric (5) or "chr5"
          chr_raw <- NULL
          if ("chr" %in% names(cl0)) chr_raw <- cl0$chr[1]
          if (is.null(chr_raw) && "cluster_chr" %in% names(cl0)) chr_raw <- cl0$cluster_chr[1]
          if (is.null(chr_raw) && "cluster_chr_n" %in% names(cl0)) chr_raw <- cl0$cluster_chr_n[1]
          
          chr_num <- suppressWarnings(as.integer(gsub("^chr", "", as.character(chr_raw))))
          validate(need(!is.na(chr_num), paste0("Cluster ", cid, ": cannot parse chr value: ", chr_raw)))
          
          chr_disease <- paste0("chr", chr_label_plink(chr_num))
          
          st <- if ("start" %in% names(cl0)) suppressWarnings(as.integer(cl0$start[1])) else NA_integer_
          en <- if ("end"   %in% names(cl0)) suppressWarnings(as.integer(cl0$end[1]))   else NA_integer_
          validate(need(!is.na(st) && !is.na(en) && st < en, paste0("Cluster ", cid, ": invalid start/end.")))
          
          incProgress(inc, detail = paste0(cid, " (", i, "/", ncl, ")"))
          append_log("[EWAS-ALL] ", cid, " ", chr_disease, ":", st, "-", en)
          
          sub_full <- extract_subset_if_needed(chr_disease, st, en)
          if (is.null(sub_full) || (is.data.frame(sub_full) && !nrow(sub_full))) {
            append_log("[EWAS-ALL] skip (no probes/subset): ", cid)
            next
          }
          
          out <- compute_bins_anydisease_for_cluster(
            sub_full = sub_full,
            coord    = rv$coord,
            meta     = rv$meta,
            chr_disease   = chr_disease, st = st, en = en,
            alpha    = alpha,
            bin_size = bin_sz,
            min_n    = min_n,
            min_cpg  = min_cpg,
            test_m   = test_m,
            use_adj_as_ctl = use_adj
          )
          
          # --- guard: out must be a list ---
          if (!is.list(out)) {
            append_log("[EWAS-ALL] skip (compute returned non-list): ", cid)
            next
          }
          
          if (is.data.frame(out$bins) && nrow(out$bins)) {
            b <- as.data.table(out$bins)
            b[, `:=`(cluster_id = cid, cluster_chr = chr_disease, cluster_start = st, cluster_end = en)]
            all_bins_list[[i]] <- b
          }
          
          if (is.data.frame(out$detail) && nrow(out$detail)) {
            d <- as.data.table(out$detail)
            d[, `:=`(cluster_id = cid, cluster_chr = chr_disease, cluster_start = st, cluster_end = en)]
            all_detail_list[[i]] <- d
          }
        }
        
        # rbind only non-null entries
        all_bins_list   <- Filter(Negate(is.null), all_bins_list)
        all_detail_list <- Filter(Negate(is.null), all_detail_list)
        
        rv$ewas_bins_all   <- if (length(all_bins_list))   rbindlist(all_bins_list, fill = TRUE)   else data.table()
        rv$ewas_detail_all <- if (length(all_detail_list)) rbindlist(all_detail_list, fill = TRUE) else data.table()
      })
      
      counts <- data.table(cluster_id = character(), n_EWAS_bins = integer())
      if (is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all)) {
        counts <- as.data.table(rv$ewas_bins_all)[, .(n_EWAS_bins = .N), by = .(cluster_id)]
      }
      
      cldt <- as.data.table(rv$clusters)
      cldt[, cluster_id := as.character(cluster_id)]
      
      if ("n_EWAS_bins" %in% names(cldt)) cldt[, n_EWAS_bins := NULL]
      cldt <- merge(cldt, counts, by = "cluster_id", all.x = TRUE)
      cldt[is.na(n_EWAS_bins), n_EWAS_bins := 0L]
      cldt[, n_EWAS_bins := as.integer(n_EWAS_bins)]
      
      # order safely
      # order safely (data.table: setorder wants column names, not expressions)
      if (all(c("chr","start","end") %in% names(cldt))) {
        cldt[, chr   := suppressWarnings(as.integer(chr))]
        cldt[, start := suppressWarnings(as.integer(start))]
        cldt[, end   := suppressWarnings(as.integer(end))]
        
        cldt <- cldt[is.finite(chr) & is.finite(start) & is.finite(end)]
        data.table::setorder(cldt, chr, start, end)
      }
      
      rv$clusters <- as.data.frame(cldt)
      
      nb <- if (!is.data.frame(rv$ewas_bins_all) || !nrow(rv$ewas_bins_all)) 0L else nrow(rv$ewas_bins_all)
      
      if (nb == 0) {
        notify_err("No significant EWAS bins across ALL clusters (try higher alpha / lower min_n / lower min_cpg).")
      } else {
        notify_info(paste0("EWAS ALL clusters computed: ", nb, " significant bins."))
      }
      
    }, error = function(e) {
      notify_err(paste0("run_ewas_all failed: ", conditionMessage(e)))
    })
    
  }, ignoreInit = TRUE)
  
  
  # -----------------------------
  # Build summary for selected cluster (frozen)
  # -----------------------------
  observeEvent(input$build_summary_sel_cluster, {
    
    if (!is.data.frame(rv$clusters) || nrow(rv$clusters) == 0) {
      notify_err("No clusters yet. Run Step 2 (➊) first.")
      return()
    }
    
    if (!is.data.frame(rv$ewas_bins_all) || nrow(rv$ewas_bins_all) == 0) {
      notify_err("No EWAS bins computed. Run Step 4 (➋) first.")
      return()
    }
    
    cl <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
    
    idx <- input$cluster_dt_rows_selected
    if (length(idx) == 0) {
      if ("n_EWAS_bins" %in% names(cl)) idx <- which(as.integer(cl$n_EWAS_bins) > 0L)[1]
      if (length(idx) == 0 || is.na(idx)) idx <- 1
      notify_info(paste0("No row selected in Clusters table. Using cluster row #", idx, "."))
    }
    
    if (!is.finite(idx[1]) || idx[1] < 1 || idx[1] > nrow(cl)) {
      notify_err("Selected cluster row is invalid.")
      return()
    }
    
    cl0 <- cl[idx[1], , drop = FALSE]
    if (!"cluster_id" %in% names(cl0) || !nzchar(as.character(cl0$cluster_id[1]))) {
      notify_err("Internal error: cluster_id missing in rv$clusters.")
      return()
    }
    
    cid <- as.character(cl0$cluster_id[1])
    
    rv$summary_cluster     <- cl0
    rv$summary_cluster_ids <- cid
    
    # load bins/detail for this cluster
    rv$ewas_bins   <- NULL
    rv$ewas_detail <- NULL
    
    b <- as.data.table(rv$ewas_bins_all)[as.character(cluster_id) == cid]
    if (nrow(b)) rv$ewas_bins <- b
    
    if (is.data.frame(rv$ewas_detail_all) && nrow(rv$ewas_detail_all)) {
      d <- as.data.table(rv$ewas_detail_all)[as.character(cluster_id) == cid]
      if (nrow(d)) rv$ewas_detail <- d
    }
    
    rv$ewas_chr <- paste0("chr", chr_label_plink(as.integer(cl0$chr[1])))
    rv$ewas_st  <- as.integer(cl0$start[1])
    rv$ewas_en  <- as.integer(cl0$end[1])
    
    reg <- paste0(rv$ewas_chr, ":", rv$ewas_st, "-", rv$ewas_en)
    ucsc_region(reg)
    
    proxy <- DT::dataTableProxy("cluster_dt")
    DT::selectRows(proxy, idx[1])
    
    if (!is.data.frame(rv$ewas_bins) || nrow(rv$ewas_bins) == 0) {
      notify_err(paste0("Selected cluster has 0 significant EWAS bins: ", cid))
      append_log("[summary] cluster=", cid, " bins=0")
      return()
    }
    
    notify_info(paste0("Summary ready for ", cid, " (", nrow(rv$ewas_bins), " bins). Click a disease in the table."))
    append_log("[summary] cluster=", cid, " bins=", nrow(rv$ewas_bins))
    
  }, ignoreInit = TRUE)
  
  # -----------------------------
  # Table: bins + clickable diseases (for frozen summary cluster)
  # -----------------------------
  current_bins_for_selected <- reactive({
    cl0 <- rv$summary_cluster
    if (is.null(cl0) || !nrow(cl0)) return(NULL)
    
    if (is.data.frame(rv$ewas_bins) && nrow(rv$ewas_bins)) {
      return(rv$ewas_bins)
    }
    
    if (is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all)) {
      cid <- as.character(cl0$cluster_id[1])
      b <- as.data.frame(rv$ewas_bins_all, stringsAsFactors = FALSE)
      b <- b[as.character(b$cluster_id) == cid, , drop = FALSE]
      if (nrow(b)) return(b)
    }
    
    NULL
  })
  
  
  current_bins_for_summary <- reactive({
    validate(need(is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all) > 0,
                  "No ewas_bins_all loaded. Run ALL clusters (➋) first."))
    
    b <- as.data.frame(rv$ewas_bins_all, stringsAsFactors = FALSE)
    if (!nrow(b)) return(NULL)
    
    # assegura tipus
    b$cluster_id <- as.character(b$cluster_id)
    
    # IMPORTANT: per al filtre de clusters, només considerem bins amb diseases no buida
    if (!"diseases" %in% names(b)) return(NULL)
    
    b$diseases <- as.character(b$diseases)
    b <- b[!is.na(b$diseases) & nzchar(trimws(b$diseases)), , drop = FALSE]
    if (!nrow(b)) return(NULL)
    
    b
  })
  
  output$bin_disease_cluster_filter_ui <- renderUI({
    bins <- current_bins_for_summary()
    
    choices <- if (is.null(bins) || !nrow(bins)) character(0) else
      sort(unique(as.character(bins$cluster_id)))
    
    selectizeInput(
      "bin_disease_cluster_filter", "Filter cluster",
      choices  = c("All" = "", choices),
      selected = "",
      multiple = FALSE,
      options  = list(placeholder = "Type to search cluster_id.", allowEmptyOption = TRUE)
    )
  })
  
  output$tbl_bin_diseases <- renderDT({
    
    bins <- current_bins_for_summary()
    
    if (is.null(bins) || !nrow(bins)) {
      return(DT::datatable(
        data.frame(Message = "No Summary yet. Run ALL clusters (➋) first."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # filtre per cluster (si n’hi ha)
    cl_f <- as.character(input$bin_disease_cluster_filter %||% "")
    if (nzchar(cl_f)) {
      bins <- bins[as.character(bins$cluster_id) == cl_f, , drop = FALSE]
    }
    validate(need(nrow(bins) > 0, "No bins for the selected cluster filter."))
    
    # -----------------------------
    # HTML cell builder
    # -----------------------------
    mk_dis_html <- function(s, cluster_id) {
      dd <- trimws(unlist(strsplit(as.character(s), ",")))
      dd <- dd[nzchar(dd)]
      if (!length(dd)) return("")
      
      cl2 <- htmltools::htmlEscape(as.character(cluster_id))
      
      diseases_raw <- paste(dd, collapse = "\n")
      diseases_enc <- utils::URLencode(diseases_raw, reserved = TRUE)
      
      paste0(
        "<span class='gi-openList diseaseOpen' ",
        "data-cluster='", cl2, "' ",
        "data-diseases-enc='", diseases_enc, "'>",
        "View diseases (", length(dd), ")",
        "</span>"
      )
    }
    
    show_df <- as.data.frame(bins, stringsAsFactors = FALSE)
    
    # HTML column
    show_df$diseases_click <- mapply(
      mk_dis_html,
      show_df$diseases,
      show_df$cluster_id,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )
    
    # numeric formatting
    show_df$best_padj <- suppressWarnings(as.numeric(show_df$best_padj))
    show_df$best_log10FDR <- suppressWarnings(round(-log10(show_df$best_padj), 2))
    
    # keep structure; diseases column is HTML
    show_df <- show_df %>%
      dplyr::transmute(
        cluster_id,
        chr,
        bin_start,
        bin_end,
        n_diseases,
        best_disease,
        best_log10FDR,
        diseases = diseases_click
      )
    
    # --- JS formatter for exports: replace "View diseases" with the REAL list ---
    js_export_body <- DT::JS("
function(data, row, column, node){
  var $sp = $(node).find('span.diseaseOpen');
  if($sp.length){
    var enc = $sp.attr('data-diseases-enc') || '';
    var raw = '';
    try { raw = decodeURIComponent(enc); } catch(e){ raw = enc; }

    // export as a single cell: disease1; disease2; ...
    return raw.split('\\n')
              .map(function(x){ return x.trim(); })
              .filter(function(x){ return x.length > 0; })
              .join('; ');
  }
  return $(node).text().trim();
}
")
    
    DT::datatable(
      show_df,
      escape = FALSE,
      rownames = FALSE,
      selection = "none",
      extensions = "Buttons",
      options = list(
        dom        = "Bfrtip",
        pageLength = 15,
        scrollX    = TRUE,
        buttons    = list(
          list(extend = "copy",  exportOptions = list(columns=":visible", format=list(body = js_export_body))),
          list(extend = "csv",   exportOptions = list(columns=":visible", format=list(body = js_export_body))),
          list(extend = "excel", exportOptions = list(columns=":visible", format=list(body = js_export_body))),
          list(extend = "print", exportOptions = list(columns=":visible", format=list(body = js_export_body)))
        )
      ),
      callback = DT::JS("
function giCloseDiseasePopup(){
  $('#giDiseasePopup').remove();
  $(document).off('mousedown.giDiseasePopup');
  $(document).off('keydown.giDiseasePopup');
}

function giEscHtml(s){
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

table.off('click', 'span.diseaseOpen');
table.on('click', 'span.diseaseOpen', function(e){
  e.preventDefault(); e.stopPropagation();
  giCloseDiseasePopup();

  var disEnc = $(this).attr('data-diseases-enc') || '';
  var cl = $(this).attr('data-cluster') || '';

  // decode raw list (newline separated)
  var disRaw = '';
  try { disRaw = decodeURIComponent(disEnc); } catch(err){ disRaw = disEnc; }
  var diseases = disRaw.split('\\n').map(function(x){ return x.trim(); }).filter(function(x){ return x.length > 0; });

  var rect = this.getBoundingClientRect();
  var left = rect.left + window.scrollX;
  var top  = rect.bottom + window.scrollY + 6;

  var html = \"<div id='giDiseasePopup'>\" +
             \"<div class='gi-title'>Diseases in cluster \" + giEscHtml(cl) +
             \"<span class='gi-close' id='giDiseasePopupClose'>✕</span></div>\";

  if(diseases.length === 0){
    html += \"<div style='color:#777;'>No diseases.</div>\";
  } else {
    for(var i=0;i<diseases.length;i++){
      var ddRaw = diseases[i];                  // RAW
      var ddLbl = giEscHtml(ddRaw);             // display
      var ddEnc = encodeURIComponent(ddRaw);    // safe attribute
      html += \"<div class='gi-item'>\" +
              \"<span class='diseasePick' data-disease-enc='\" + ddEnc + \"' data-cluster='\" + giEscHtml(cl) + \"'>\" +
              ddLbl +
              \"</span></div>\";
    }
  }
  html += \"</div>\";

  $('body').append(html);

  var $p = $('#giDiseasePopup');
  $p.css({ left: left + 'px', top: top + 'px' });

  $('#giDiseasePopupClose').on('click', function(){ giCloseDiseasePopup(); });

  $(document).on('mousedown.giDiseasePopup', function(ev){
    if($(ev.target).closest('#giDiseasePopup').length === 0){
      giCloseDiseasePopup();
    }
  });

  $(document).on('keydown.giDiseasePopup', function(ev){
    if(ev.key === 'Escape'){ giCloseDiseasePopup(); }
  });
});

$(document).off('click', '#giDiseasePopup .diseasePick');
$(document).on('click', '#giDiseasePopup .diseasePick', function(e){
  e.preventDefault(); e.stopPropagation();

  var ddEnc = $(this).attr('data-disease-enc') || '';
  var cl = $(this).attr('data-cluster') || '';

  var ddRaw = '';
  try { ddRaw = decodeURIComponent(ddEnc); } catch(err){ ddRaw = ddEnc; }

  // close first, then notify Shiny
  giCloseDiseasePopup();

  setTimeout(function(){
    Shiny.setInputValue('disease_pick', {disease: ddRaw, cluster_id: cl}, {priority: 'event'});
  }, 0);
});
")
    ) %>%
      DT::formatStyle("n_diseases", fontWeight = "bold")
    
  })
  
  output$tbl_bin_diseasesXXXX <- renderDT({
    
    bins <- current_bins_for_summary()
    
    if (is.null(bins) || !nrow(bins)) {
      return(DT::datatable(
        data.frame(Message = "No Summary yet. Run ALL clusters (➋) first."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # filtre per cluster (si n’hi ha)
    cl_f <- as.character(input$bin_disease_cluster_filter %||% "")
    if (nzchar(cl_f)) {
      bins <- bins[as.character(bins$cluster_id) == cl_f, , drop = FALSE]
    }
    validate(need(nrow(bins) > 0, "No bins for the selected cluster filter."))
    
    # -----------------------------
    # HTML cell builder
    # - clickable element in the column
    # - stores disease list URL-encoded (newline-separated)
    # -----------------------------
    mk_dis_html <- function(s, cluster_id) {
      dd <- trimws(unlist(strsplit(as.character(s), ",")))
      dd <- dd[nzchar(dd)]
      if (!length(dd)) return("")
      
      cl2 <- htmltools::htmlEscape(as.character(cluster_id))
      
      diseases_raw <- paste(dd, collapse = "\n")
      diseases_enc <- utils::URLencode(diseases_raw, reserved = TRUE)
      
      paste0(
        "<span class='gi-openList diseaseOpen' ",
        "data-cluster='", cl2, "' ",
        "data-diseases-enc='", diseases_enc, "'>",
        "View diseases (", length(dd), ")",
        "</span>"
      )
    }
    
    show_df <- as.data.frame(bins, stringsAsFactors = FALSE)
    
    # HTML column
    show_df$diseases_click <- mapply(
      mk_dis_html,
      show_df$diseases,
      show_df$cluster_id,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )
    
    # numeric formatting
    show_df$best_padj <- suppressWarnings(as.numeric(show_df$best_padj))
    show_df$best_log10FDR <- suppressWarnings(round(-log10(show_df$best_padj), 2))
    
    # keep structure; diseases column is HTML
    show_df <- show_df %>%
      dplyr::transmute(
        cluster_id,
        chr,
        bin_start,
        bin_end,
        n_diseases,
        best_disease,
        best_log10FDR,
        diseases = diseases_click
      )
    
    DT::datatable(
      show_df,
      escape = FALSE,
      rownames = FALSE,
      selection = "none",
      options = list(pageLength = 15, scrollX = TRUE),
      callback = DT::JS("
function giCloseDiseasePopup(){
  $('#giDiseasePopup').remove();
  $(document).off('mousedown.giDiseasePopup');
  $(document).off('keydown.giDiseasePopup');
}

function giEscHtml(s){
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

table.off('click', 'span.diseaseOpen');
table.on('click', 'span.diseaseOpen', function(e){
  e.preventDefault(); e.stopPropagation();
  giCloseDiseasePopup();

  var disEnc = $(this).attr('data-diseases-enc') || '';
  var cl = $(this).attr('data-cluster') || '';

  // decode raw list (newline separated)
  var disRaw = '';
  try { disRaw = decodeURIComponent(disEnc); } catch(err){ disRaw = disEnc; }
  var diseases = disRaw.split('\\n').map(function(x){ return x.trim(); }).filter(function(x){ return x.length > 0; });

  var rect = this.getBoundingClientRect();
  var left = rect.left + window.scrollX;
  var top  = rect.bottom + window.scrollY + 6;

  var html = \"<div id='giDiseasePopup'>\" +
             \"<div class='gi-title'>Diseases in cluster \" + giEscHtml(cl) +
             \"<span class='gi-close' id='giDiseasePopupClose'>✕</span></div>\";

  if(diseases.length === 0){
    html += \"<div style='color:#777;'>No diseases.</div>\";
  } else {
    for(var i=0;i<diseases.length;i++){
      var ddRaw = diseases[i];                  // RAW
      var ddLbl = giEscHtml(ddRaw);             // display
      var ddEnc = encodeURIComponent(ddRaw);    // safe attribute
      html += \"<div class='gi-item'>\" +
              \"<span class='diseasePick' data-disease-enc='\" + ddEnc + \"' data-cluster='\" + giEscHtml(cl) + \"'>\" +
              ddLbl +
              \"</span></div>\";
    }
  }
  html += \"</div>\";

  $('body').append(html);

  var $p = $('#giDiseasePopup');
  $p.css({ left: left + 'px', top: top + 'px' });

  $('#giDiseasePopupClose').on('click', function(){ giCloseDiseasePopup(); });

  $(document).on('mousedown.giDiseasePopup', function(ev){
    if($(ev.target).closest('#giDiseasePopup').length === 0){
      giCloseDiseasePopup();
    }
  });

  $(document).on('keydown.giDiseasePopup', function(ev){
    if(ev.key === 'Escape'){ giCloseDiseasePopup(); }
  });
});

$(document).off('click', '#giDiseasePopup .diseasePick');
$(document).on('click', '#giDiseasePopup .diseasePick', function(e){
  e.preventDefault(); e.stopPropagation();

  var ddEnc = $(this).attr('data-disease-enc') || '';
  var cl = $(this).attr('data-cluster') || '';

  var ddRaw = '';
  try { ddRaw = decodeURIComponent(ddEnc); } catch(err){ ddRaw = ddEnc; }

  // close first, then notify Shiny (avoids redraw/DOM conflicts)
  giCloseDiseasePopup();

  setTimeout(function(){
    Shiny.setInputValue('disease_pick', {disease: ddRaw, cluster_id: cl}, {priority: 'event'});
  }, 0);
});
")
    ) %>%
      DT::formatStyle("n_diseases", fontWeight = "bold")
    
  })
  
  
  
  # -----------------------------
  # Ensure selected subset exists for frozen cluster
  # -----------------------------
  ensure_selected_subset <- function() {
    cl0 <- rv$summary_cluster
    validate(need(is.data.frame(cl0) && nrow(cl0) == 1, "Click 'Build summary (selected cluster)' first."))
    validate(need(is.data.table(rv$meta) && nrow(rv$meta) > 0, "meta_slim_ds not loaded."))
    validate(need(is.data.table(rv$coord) && nrow(rv$coord) > 0, "coord_hg38 not loaded."))
    validate(need(dir.exists(input$slim_dir), "slim_ds_by_chr directory not found."))
    validate(need(nzchar(Sys.which("python3")), "python3 not found on PATH."))
    
    chr_disease <- paste0("chr", chr_label_plink(as.integer(cl0$chr[1])))
    st <- as.integer(cl0$start[1])
    en <- as.integer(cl0$end[1])
    
    sub_full <- extract_subset_if_needed(chr_disease, st, en)
    validate(need(!is.null(sub_full) && file.exists(sub_full), "Subset extraction failed / no probes in interval."))
    
    rv$ewas_sub_full <- sub_full
    rv$ewas_chr <- chr_disease
    rv$ewas_st  <- st
    rv$ewas_en  <- en
    TRUE
  }
  
  # read columns at ewasdis
  read_subset_agnosticXXXXXXXXXXX <- function(path, sel_cols) {
    stopifnot(file.exists(path))
    # read only header to detect first column name
    hdr <- data.table::fread(path, nrows = 0, check.names = FALSE)
    idcol <- names(hdr)[1]
    
    cols_to_read <- unique(c(idcol, sel_cols))
    
    dt <- data.table::fread(
      path,
      select = cols_to_read,
      na.strings = c("NA","<NA>","NaN",""),
      check.names = FALSE
    )
    
    # normalize first column to sample_id internally
    if (idcol != "sample_id") data.table::setnames(dt, idcol, "sample_id")
    dt
  }
  # =========================
  # Helpers
  # =========================
  
  must_path1 <- function(p, tag="path") {
    if (is.null(p) || length(p) != 1L || !is.character(p) || is.na(p) || !nzchar(p)) {
      stop(sprintf("[%s] invalid path: %s", tag, paste(capture.output(str(p)), collapse=" ")), call. = FALSE)
    }
    p
  }
  
  read_subset_agnostic <- function(path, sel_cols) {
    path <- must_path1(path, "subset_file")
    if (!file.exists(path)) stop("[subset_file] not found: ", path, call. = FALSE)
    
    hdr <- data.table::fread(path, nrows = 0, check.names = FALSE)
    idcol <- names(hdr)[1]
    
    cols_to_read <- unique(c(idcol, sel_cols))
    
    dt <- data.table::fread(
      path,
      select = cols_to_read,
      na.strings = c("NA","<NA>","NaN",""),
      check.names = FALSE
    )
    
    if (idcol != "sample_id") data.table::setnames(dt, idcol, "sample_id")
    dt
  }
  
  
  # -----------------------------
  # Per-disease CpG delta computation
  # -----------------------------
  run_for_disease <- function(disease_name) {
    req(rv$meta, rv$coord)
    ensure_selected_subset()
    validate(need(file.exists(rv$ewas_sub_full), "Subset file not found (select cluster again)."))
    
    meta <- rv$meta
    chr_disease <- rv$ewas_chr; st <- rv$ewas_st; en <- rv$ewas_en
    
    meta_c <- meta[disease == disease_name &
                     sample_type %in% c("disease tissue", "control", "adjacent normal"),
                   .(sample_id, disease, sample_type)]
    meta_c <- unique(meta_c, by = "sample_id")
    validate(need(nrow(meta_c) > 0, "No samples found for this disease."))
    
    dis_ids <- meta_c[sample_type == "disease tissue", sample_id]
    ctl_ids <- meta_c[sample_type == "control", sample_id]
    adj_ids <- meta_c[sample_type == "adjacent normal", sample_id]
    
    if (length(ctl_ids) == 0 && isTRUE(input$use_adj_as_ctl) && length(adj_ids) > 0) {
      ctl_ids <- adj_ids
    }
    
    validate(need(length(dis_ids) > 0, "No disease samples."))
    validate(need(length(ctl_ids) > 0, "No control samples (and adjacent normal not used/available)."))
    
    sel <- unique(c(dis_ids, ctl_ids))
    cols_to_read <- c("sample_id", sel)
    
    dt <- data.table::fread(rv$ewas_sub_full,
                            select = cols_to_read,
                            na.strings = c("NA","<NA>","NaN",""),
                            check.names = FALSE)
    
    beta_dt <- dt[grepl("^cg\\d{8}$", sample_id)]
    validate(need(nrow(beta_dt) > 0, "Subset contains 0 CpG rows."))
    
    beta_dt[, (sel) := lapply(.SD, as.numeric), .SDcols = sel]
    beta_mat <- as.matrix(beta_dt[, ..sel])
    rownames(beta_mat) <- beta_dt$sample_id
    
    coord2 <- rv$coord[chr == chr_disease & pos >= st & pos <= en]
    pos <- coord2$pos[match(rownames(beta_mat), coord2$probe)]
    
    dis <- intersect(colnames(beta_mat), dis_ids)
    ctl <- intersect(colnames(beta_mat), ctl_ids)
    
    delta <- rowMeans(beta_mat[, dis, drop = FALSE], na.rm = TRUE) -
      rowMeans(beta_mat[, ctl, drop = FALSE], na.rm = TRUE)
    
    res <- data.table(
      probe = names(delta),
      delta_beta = as.numeric(delta),
      hyper_hypo = fifelse(delta > 0, "hyper", fifelse(delta < 0, "hypo", "no_change")),
      chr = chr_disease,
      pos = as.integer(pos)
    )
    
    rv$res <- res
    rv$beta_mat <- beta_mat
    rv$meta_c <- meta_c
  }
  
  # Click disease from bins table
  
  observeEvent(input$disease_pick, {
    req(input$disease_pick$disease, input$disease_pick$cluster_id)
    
    disease     <- as.character(input$disease_pick$disease)
    cluster_id <- as.character(input$disease_pick$cluster_id)
    
    append_log(paste0("[disease_pick] ", disease, " | cluster=", cluster_id))
    
    # Si per qualsevol motiu el click ve d'un altre cluster, actualitza el "frozen cluster"
    if (is.null(rv$summary_cluster) || !nrow(rv$summary_cluster) ||
        as.character(rv$summary_cluster$cluster_id[1]) != cluster_id) {
      
      cl_all <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
      one <- cl_all[as.character(cl_all$cluster_id) == cluster_id, , drop = FALSE]
      validate(need(nrow(one) == 1, "Cluster not found for this disease click."))
      
      rv$summary_cluster     <- one
      rv$summary_cluster_ids <- cluster_id
      
      rv$ewas_chr <- paste0("chr", chr_label_plink(as.integer(one$chr[1])))
      rv$ewas_st  <- as.integer(one$start[1])
      rv$ewas_en  <- as.integer(one$end[1])
    }
    
    rv$disease_sel <- as.character(input$disease_pick$disease)
    run_for_disease(disease)
    
  }, ignoreInit = TRUE)
  
  
  
  # -----------------------------
  # CpG plots
  # -----------------------------
  
  output$p_hist <- renderPlot({
    req(rv$res)
    
    ggplot(as.data.frame(rv$res), aes(x = delta_beta, fill = hyper_hypo)) +
      geom_histogram(bins = 40, alpha = 0.85) +
      geom_vline(xintercept = 0) +
      scale_fill_manual(
        values = c("hyper" = "red", "hypo" = "darkblue"),
        breaks = c("hyper", "hypo"),
        drop = FALSE
      ) +
      labs(x = "Δbeta (disease - control)", y = "CpGs") +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center"
      )
  })
  
  
  output$p_beta_dist <- renderPlot({
    req(rv$beta_mat, rv$meta_c)
    
    meta <- rv$meta_c
    dis_ids <- meta[sample_type == "disease tissue", as.character(sample_id)]
    ctl_ids <- meta[sample_type %in% c("control","adjacent normal"), as.character(sample_id)]
    
    dis_ids <- intersect(colnames(rv$beta_mat), dis_ids)
    ctl_ids <- intersect(colnames(rv$beta_mat), ctl_ids)
    validate(need(length(dis_ids) > 0 && length(ctl_ids) > 0, "No disease/control samples."))
    
    maxS <- 600L
    if (length(dis_ids) > maxS) dis_ids <- sample(dis_ids, maxS)
    if (length(ctl_ids) > maxS) ctl_ids <- sample(ctl_ids, maxS)
    sel <- unique(c(dis_ids, ctl_ids))
    
    mat <- rv$beta_mat[, sel, drop = FALSE]
    long <- as.data.table(as.table(mat))
    setnames(long, c("probe","sample_id","bval"))
    long[, bval := as.numeric(bval)]
    long <- long[is.finite(bval) & bval >= 0 & bval <= 1]
    long[, grp := fifelse(sample_id %in% dis_ids, "disease", "control")]
    
    max_points <- 200000L
    if (nrow(long) > max_points) long <- long[sample(.N, max_points)]
    
    ggplot(long, aes(x = bval, fill = grp, color = grp)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 60, alpha = 0.30, position = "identity",
        linewidth = 0
      ) +
      geom_density(alpha = 0.15) +
      scale_fill_manual(values = c("disease" = "#ff7f00", "control" = "darkgreen")) +
      scale_color_manual(values = c("disease" = "#ff7f00", "control" = "darkgreen")) +
      labs(
        x = "Beta", y = "Density",
        title = paste0(rv$disease_sel %||% "", " — disease vs control")
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center"
      )
  })
  
  
  # =========================
  # p_violin_sel (FOCUSED BIN)
  # - Violin: mean beta per sample for ONE focused bin (control vs disease)
  # - Scatter: CpG Δbeta inside that SAME bin interval
  # - Shared X axis in bp (start..end) for BOTH subplot panels (xaxis + xaxis2)
  # - Split violins (left/right) with width proportional to bin_size (no compression)
  # =========================
  output$p_violin_sel <- plotly::renderPlotly({
    
    req(rv$meta, rv$coord, rv$disease_sel)
    ensure_selected_subset()
    
    validate(need(file.exists(rv$ewas_sub_full), "Subset file not found."))
    
    # -----------------------------
    # Helpers: normalize + resolve disease name vs meta
    # -----------------------------
    norm_key_disease <- function(x) {
      x <- enc2utf8(trimws(as.character(x)))
      x <- tolower(x)
      x <- gsub("&apos;|&#39;|&#x27;|&#x2019;|&#8217;", "'", x, ignore.case = TRUE)
      x <- gsub("[\u2018\u2019\u02BC\u0060\u00B4]", "'", x, perl = TRUE)
      x <- gsub("['’`´]", "", x)
      x <- gsub("[^a-z0-9]+", " ", x)
      x <- gsub("\\s+", " ", x)
      trimws(x)
    }
    
    resolve_disease_in_meta <- function(disease_sel, meta_disease_vec) {
      sel_raw <- as.character(disease_sel)
      sel_key <- norm_key_disease(sel_raw)
      
      meta_u <- unique(as.character(meta_disease_vec))
      meta_u <- meta_u[!is.na(meta_u) & nzchar(trimws(meta_u))]
      meta_k <- norm_key_disease(meta_u)
      
      hit <- which(meta_k == sel_key)
      if (length(hit)) return(meta_u[hit[1]])
      
      hit <- which(meta_k == paste0(sel_key, " disease"))
      if (length(hit)) return(meta_u[hit[1]])
      
      hit <- which(startsWith(meta_k, sel_key))
      if (length(hit)) {
        cand <- meta_u[hit]
        cand <- cand[order(nchar(cand))]
        return(cand[1])
      }
      
      hit <- which(grepl(paste0("\\b", sel_key, "\\b"), meta_k))
      if (length(hit)) {
        cand <- meta_u[hit]
        cand <- cand[order(nchar(cand))]
        return(cand[1])
      }
      
      d <- adist(sel_key, meta_k)
      j <- which.min(d)
      if (length(j) == 1 && is.finite(d[j]) && d[j] <= max(2, floor(nchar(sel_key) * 0.25))) {
        return(meta_u[j])
      }
      
      sel_raw
    }
    
    # -----------------------------
    # Selected disease (raw -> resolved -> key)
    # -----------------------------
    disease_sel <- as.character(rv$disease_sel %||% input$disease_sel %||% input$disease %||% "")
    validate(need(nzchar(disease_sel), "No disease selected."))
    
    meta0 <- data.table::as.data.table(rv$meta)
    meta0[, disease_key := norm_key_disease(disease)]
    
    disease_name_meta <- resolve_disease_in_meta(disease_sel, meta0$disease)
    sel_key <- norm_key_disease(disease_name_meta)
    
    # -----------------------------
    # Params
    # -----------------------------
    bin_size <- max(1000L, as.integer(input$bin_size))
    alpha    <- as.numeric(input$alpha_any)
    min_n    <- as.integer(input$min_n_any)
    min_cpg  <- as.integer(input$min_cpg_any)
    test_m   <- as.character(input$bin_test_any)
    use_adj  <- isTRUE(input$use_adj_as_ctl)
    
    test_name <- if (identical(test_m, "ttest")) "t-test" else "Wilcoxon"
    
    chr_disease <- rv$ewas_chr
    st0    <- as.integer(rv$ewas_st)
    en0    <- as.integer(rv$ewas_en)
    
    validate(need(nzchar(chr_disease) && is.finite(st0) && is.finite(en0) && st0 < en0, "Invalid chr/interval."))
    
    # -----------------------------
    # 1) Read subset (CpG x samples)
    # -----------------------------
    dt_all <- data.table::fread(
      rv$ewas_sub_full,
      na.strings  = c("NA","<NA>","NaN",""),
      check.names = FALSE,
      showProgress = FALSE
    )
    
    beta_dt <- dt_all[grepl("^cg\\d{8}$", sample_id)]
    validate(need(nrow(beta_dt) > 0, "Subset has 0 CpG rows."))
    
    samp_cols <- setdiff(names(beta_dt), "sample_id")
    validate(need(length(samp_cols) > 0, "Subset has 0 sample columns."))
    
    beta_dt[, (samp_cols) := lapply(.SD, suppressWarnings(as.numeric)), .SDcols = samp_cols]
    beta_mat_all <- as.matrix(beta_dt[, ..samp_cols])
    rownames(beta_mat_all) <- beta_dt$sample_id
    
    # -----------------------------
    # 2) CpG positions + bins (within current interval st0..en0)
    # -----------------------------
    coord2 <- rv$coord[chr == chr_disease & pos >= st0 & pos <= en0]
    validate(need(nrow(coord2) > 0, "No coord rows for this chr/interval."))
    
    pos_vec <- coord2$pos[match(rownames(beta_mat_all), coord2$probe)]
    keep <- is.finite(pos_vec)
    
    beta_mat_all <- beta_mat_all[keep, , drop = FALSE]
    pos_vec <- as.integer(pos_vec[keep])
    validate(need(nrow(beta_mat_all) > 0, "No CpGs with valid hg38 position in this interval."))
    
    bin_start <- st0 + ((pos_vec - st0) %/% bin_size) * bin_size
    bin_mid   <- bin_start + bin_size/2
    
    # -----------------------------
    # 3) Mean beta per sample per bin (sb)
    # -----------------------------
    idx_list <- split(seq_along(bin_mid), bin_mid)
    idx_list <- idx_list[vapply(idx_list, length, 1L) >= min_cpg]
    validate(need(length(idx_list) > 0, "No bins pass min CpGs/bin."))
    
    sb_list <- lapply(names(idx_list), function(bm) {
      ii <- idx_list[[bm]]
      mb <- colMeans(beta_mat_all[ii, , drop = FALSE], na.rm = TRUE)
      data.table::data.table(
        bin_mid   = as.numeric(bm),
        sample_id = colnames(beta_mat_all),
        mean_beta = as.numeric(mb),
        n_cpg     = length(ii)
      )
    })
    sb <- data.table::rbindlist(sb_list, use.names = TRUE, fill = TRUE)
    sb <- sb[is.finite(mean_beta) & mean_beta >= 0 & mean_beta <= 1]
    validate(need(nrow(sb) > 0, "No valid mean_beta values in bins."))
    
    # -----------------------------
    # 4) Samples: disease vs control (robust mapping)
    # -----------------------------
    stype <- tolower(trimws(as.character(meta0$sample_type)))
    
    is_case <- stype %in% tolower(c("disease tissue")) |
      grepl("\\b(case|patient|affected|disease)\\b", stype)
    
    is_control <- stype %in% tolower(c("control","adjacent normal","normal","healthy")) |
      grepl("\\b(control|healthy|normal|unaffected|adjacent)\\b", stype)
    
    meta0[, grp := data.table::fifelse(is_case, "disease",
                                       data.table::fifelse(is_control, "control", NA_character_))]
    
    meta_c <- meta0[disease_key == sel_key & !is.na(grp),
                    .(sample_id = as.character(sample_id), grp)]
    meta_c <- unique(meta_c, by = "sample_id")
    
    validate(need(nrow(meta_c) > 0,
                  paste0("No usable case/control samples for: ", disease_name_meta)))
    
    dis_ids <- meta_c[grp == "disease", sample_id]
    ctl_ids <- meta_c[grp == "control", sample_id]
    
    if (length(ctl_ids) == 0 && isTRUE(use_adj)) {
      ctl_ids <- meta0[disease_key == sel_key & grepl("adjacent", tolower(sample_type)),
                       as.character(sample_id)]
    }
    
    validate(need(length(dis_ids) > 0, "No disease/case samples found."))
    validate(need(length(ctl_ids) > 0, "No control/normal samples found."))
    
    # -----------------------------
    # 5) Bin tests + BH (pdt)
    # -----------------------------
    sb_c <- sb[sample_id %in% c(dis_ids, ctl_ids)]
    validate(need(nrow(sb_c) > 0, "No beta values after filtering disease/control samples."))
    
    sb_c[, grp := data.table::fifelse(sample_id %in% dis_ids, "disease", "control")]
    
    pdt <- sb_c[, .(
      n_dis = sum(grp == "disease"),
      n_ctl = sum(grp == "control"),
      p = {
        x <- mean_beta[grp == "disease"]
        y <- mean_beta[grp == "control"]
        if (length(x) < min_n || length(y) < min_n) NA_real_ else {
          if (identical(test_m, "ttest")) {
            tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
          } else {
            tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
          }
        }
      }
    ), by = .(bin_mid)]
    pdt[, padj := p.adjust(p, method = "BH")]
    
    padj_min <- suppressWarnings(min(pdt$padj, na.rm = TRUE))
    if (!is.finite(padj_min)) padj_min <- NA_real_
    
    i_min <- suppressWarnings(which.min(pdt$padj))
    bin_min <- if (length(i_min) && is.finite(pdt$padj[i_min])) pdt$bin_mid[i_min] else NA_real_
    
    padj_txt <- if (is.na(padj_min)) "NA" else formatC(padj_min, format = "e", digits = 2)
    
    sig_sel <- unique(pdt[is.finite(padj) & padj < alpha, bin_mid])
    validate(need(length(sig_sel) > 0, paste0("No significant bins at BH<", alpha)))
    
    # -----------------------------
    # 5b) Focus bin (robust): choose ONE bin to plot
    # - priority: rv$bin_mid_sel (if you store it elsewhere) else bin_min else first sig
    # -----------------------------
    bin_focus <- suppressWarnings(as.numeric(rv$bin_mid_sel %||% bin_min))
    if (!is.finite(bin_focus) && length(sig_sel) > 0) bin_focus <- suppressWarnings(as.numeric(sig_sel[1]))
    if (!is.finite(bin_focus) && nrow(pdt) > 0) bin_focus <- suppressWarnings(as.numeric(pdt$bin_mid[1]))
    validate(need(is.finite(bin_focus), "Cannot determine focused bin."))
    
    bin_focus_start <- as.integer(round(bin_focus - bin_size/2))
    bin_focus_end   <- as.integer(round(bin_focus + bin_size/2))
    bin_focus_start <- max(bin_focus_start, st0)
    bin_focus_end   <- min(bin_focus_end, en0)
    
    # -----------------------------
    # 6) Violin data: ONLY focused bin
    # -----------------------------
    sb_plot <- sb_c[bin_mid == bin_focus]
    validate(need(nrow(sb_plot) > 0, "No rows for focused bin."))
    
    # factor order
    sb_plot[, grp := factor(grp, levels = c("control","disease"))]
    
    # hover
    sb_plot[, hover := paste0(
      "bin_mid=", format(bin_mid, scientific = FALSE),
      "<br>grp=", grp,
      "<br>mean_beta=", signif(mean_beta, 4),
      "<br>n_cpg=", n_cpg
    )]
    
    # violin width in bp (CRITICAL when x is bp)
    w_bin <- 0.90 * bin_size
    
    p_violin <- plotly::plot_ly() %>%
      plotly::add_trace(
        data = sb_plot[grp == "control"],
        type = "violin",
        x = ~bin_mid, y = ~mean_beta,
        name = "control",
        text = ~hover, hoverinfo = "text",
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        points = "outliers",
        legendgroup = "control",
        offsetgroup = "control",
        width = w_bin,
        fillcolor = "darkgreen",
        line = list(color = "black"),
        side = "negative"
      ) %>%
      plotly::add_trace(
        data = sb_plot[grp == "disease"],
        type = "violin",
        x = ~bin_mid, y = ~mean_beta,
        name = "disease",
        text = ~hover, hoverinfo = "text",
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        points = "outliers",
        legendgroup = "disease",
        offsetgroup = "disease",
        width = w_bin,
        fillcolor = "orange",
        line = list(color = "black"),
        side = "positive"
      ) %>%
      plotly::layout(
        yaxis = list(title = "Mean beta per sample (bin)"),
        violinmode = "overlay"  # best for split violins
      )
    
    # -----------------------------
    # 7) CpG-level Δbeta within focused bin
    # -----------------------------
    dis_cols <- intersect(colnames(beta_mat_all), dis_ids)
    ctl_cols <- intersect(colnames(beta_mat_all), ctl_ids)
    validate(need(length(dis_cols) >= 2 && length(ctl_cols) >= 2, "Not enough sample columns in subset for CpG scatter."))
    
    mean_dis <- rowMeans(beta_mat_all[, dis_cols, drop = FALSE], na.rm = TRUE)
    mean_ctl <- rowMeans(beta_mat_all[, ctl_cols, drop = FALSE], na.rm = TRUE)
    
    res_cpg <- data.table::data.table(
      probe      = rownames(beta_mat_all),
      pos        = as.integer(pos_vec),
      delta_beta = as.numeric(mean_dis - mean_ctl)
    )
    
    # assign CpG to bin (same scheme)
    res_cpg[, bin_start := st0 + ((pos - st0) %/% bin_size) * bin_size]
    res_cpg[, bin_mid   := bin_start + bin_size/2]
    
    # ONLY focused bin
    res_cpg <- res_cpg[bin_mid == bin_focus]
    validate(need(nrow(res_cpg) > 0, "No CpGs found inside focused bin."))
    
    res_cpg <- res_cpg[is.finite(pos) & is.finite(delta_beta)]
    
    res_cpg[, hyper_hypo := data.table::fifelse(
      delta_beta > 0, "hyper",
      data.table::fifelse(delta_beta < 0, "hypo", "no_change")
    )]
    
    res_cpg[, hover2 := paste0(
      "probe=", probe,
      "<br>pos=", pos,
      "<br>Δbeta=", signif(delta_beta, 4),
      "<br>", hyper_hypo
    )]
    
    # Top CpGs labels
    topN <- 20L
    top  <- res_cpg[order(-abs(delta_beta))][1:min(topN, .N)]
    
    # Scatter: x = pos (bp) aligns naturally to start..end
    p_scatter <- plotly::plot_ly(
      data = res_cpg,
      type = "scatter",
      mode = "markers",
      x = ~pos,
      y = ~delta_beta,
      color = ~hyper_hypo,
      colors = c("hypo" = "darkblue", "hyper" = "red", "no_change" = "gray70"),
      text = ~hover2, hoverinfo = "text",
      marker = list(size = 5, opacity = 0.75),
      showlegend = FALSE
    ) %>%
      # y=0 baseline
      plotly::add_lines(
        x = c(bin_focus_start, bin_focus_end), y = c(0, 0),
        inherit = FALSE,
        line = list(width = 1),
        showlegend = FALSE,
        hoverinfo = "skip"
      ) %>%
      # top labels
      plotly::add_trace(
        data = top,
        type = "scatter",
        mode = "markers+text",
        x = ~pos,
        y = ~delta_beta,
        text = ~probe,
        textposition = "top center",
        marker = list(size = 7),
        showlegend = FALSE
      ) %>%
      plotly::layout(
        yaxis = list(title = "Δbeta (disease − control)")
      )
    
    # -----------------------------
    # 8) Shared X axis settings (apply to BOTH xaxis and xaxis2)
    # -----------------------------
    xax <- list(
      title = paste0("Focused bin interval (hg38) — bin_size=", bin_size, " bp"),
      range = c(bin_focus_start, bin_focus_end),
      tickmode = "array",
      tickvals = c(bin_focus_start, bin_focus, bin_focus_end),
      ticktext = c("start", "mid", "end"),
      tickangle = 0
    )
    
    # -----------------------------
    # 9) Subplot: shared X, aligned axes
    # -----------------------------
    plotly::subplot(
      p_violin, p_scatter,
      nrows = 2,
      shareX = TRUE,
      heights = c(0.58, 0.42),
      titleX = TRUE
    ) %>%
      plotly::layout(
        title = list(
          text = paste0(
            disease_name_meta, " — ", chr_disease, ":", bin_focus_start, "-", bin_focus_end,
            "<br>| BH<", alpha,
            " | min BH padj=", padj_txt,
            " | Test: ", test_name
          ),
          font = list(size = 12),
          y = 0.95, yanchor = "top"
        ),
        margin = list(t = 90),
        legend = list(
          orientation = "h",
          x = 0, xanchor = "left",
          y = 0.92, yanchor = "top"
        ),
        xaxis  = xax,
        xaxis2 = xax,
        yaxis  = list(title = "Mean beta per sample (bin)"),
        yaxis2 = list(title = "Δbeta (disease − control)")
      )
    
  })
  
  # -----------------------------
  # Manhattan_combo (GWAS + EWAS bins) + cluster segments
  # -----------------------------
  
  output$manhattan_combo <- renderPlotly({
    src_combo <- "manhattan_combo"
    
    dfp <- tryCatch(dfp_manhattan(), error = function(e) NULL)
    if (is.null(dfp) || !nrow(dfp)) return(plotly_message("⚠️ GWAS table missing or incomplete."))
    
    ref <- .ref_hg38
    ax  <- axis_df()
    axis_breaks <- ax$center
    axis_labels <- paste0("chr", ax$chrN)
    GENOME_END  <- max(ref$chr_cum + ref$len)
    
    thr_y <- if ((input$cluster_method %||% "window") == "window") (input$pthr %||% 5) else (input$min_logp %||% 6)
    
    dfp <- dfp %>%
      dplyr::arrange(CHR, BP) %>%
      dplyr::mutate(
        rs_show = ifelse(!is.na(snp) & nzchar(snp), snp, paste0("chr", CHR, ":", BP)),
        chr_lab = chr_label_plink(as.integer(CHR)),
        col     = ifelse((as.integer(CHR) %% 2) == 0, "darkgreen", "#ff7f00")
      )
    
    dfp$tooltip <- paste0(
      "<b>GWAS hit</b>",
      "<br><b>rs:</b> ", dfp$rs_show,
      "<br><b>CHR:</b> ", dfp$chr_lab,
      "<br><b>BP:</b> ", dfp$BP,
      "<br><b>P:</b> ", signif(dfp$Pval, 3),
      "<br><b>-log10(P):</b> ", round(dfp$logp, 2)
    )
    
    p1 <- ggplot(dfp, aes(x = BPcum, y = logp, text = tooltip)) +
      geom_point(aes(color = col), size = 1) +
      geom_hline(yintercept = thr_y, linetype = "dashed") +
      scale_color_identity(guide = "none") +
      scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      labs(x = NULL, y = "-log10(P)") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    p1_pl <- ggplotly(p1, tooltip = "text", source = src_combo)
    
    # ----------------------------
    # p2: EWAS bins (prefer ALL clusters if available, else selected)
    #   FIX: sempre fixa rang X (0..GENOME_END) perquè els segments de clusters
    #        NO quedin fora de rang quan no hi ha bins.
    # ----------------------------
    
    p2_base <- ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
      geom_blank() +
      scale_x_continuous(
        limits = c(0, GENOME_END),
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = "Genome", y = "-log10(FDR) [EWAS bins]") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    p2_pl <- ggplotly(p2_base, source = src_combo)
    
    bins <- NULL
    mode <- NULL
    if (is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all)) {
      bins <- rv$ewas_bins_all
      mode <- "ALL"
    } else if (is.data.frame(rv$ewas_bins) && nrow(rv$ewas_bins)) {
      bins <- rv$ewas_bins
      mode <- "SEL"
    }
    
    df2 <- NULL
    
    # y_max = top of bins if present; else default
    has_bins <- is.data.frame(bins) && nrow(bins)
    y_max <- 1
    
    if (has_bins) {
      df2 <- bins_to_manhattan_df(bins)
      
      if (is.data.frame(df2) && nrow(df2)) {
        df2 <- df2 %>%
          dplyr::mutate(
            cluster_id = if ("cluster_id" %in% names(df2)) as.character(cluster_id) else NA_character_,
            customdata = Map(list, chr, bin_start, bin_end),
            tooltip = paste0(
              "<b>EWAS bin (", mode, ")</b>",
              ifelse(!is.na(cluster_id) & nzchar(cluster_id), paste0("<br><b>cluster:</b> ", cluster_id), ""),
              "<br><b>chr:</b> ", chr,
              "<br><b>bin:</b> ", bin_start, "-", bin_end,
              "<br><b>n_disease:</b> ", n_diseases,
              "<br><b>best_disease:</b> ", best_disease,
              "<br><b>-log10(FDR):</b> ", round(y, 2)
            )
          )
        
        y_max <- max(df2$y[is.finite(df2$y)], na.rm = TRUE)
        if (!is.finite(y_max) || y_max <= 0) y_max <- 1
        
        # Placeholder; we'll set final y-range after adding cluster band
        p2 <- ggplot(df2, aes(x = BPcum, y = y, text = tooltip, key = cluster_id, customdata = customdata)) +
          geom_point(
            aes(fill = y),
            shape = 23,
            size  = 6,
            stroke = 0.3,
            color = "black",
            alpha = 0.90
          ) +
          scale_fill_gradient(
            low  = "yellow",
            high = "red",
            name = "-log10(FDR)"
          ) +
          scale_x_continuous(
            limits = c(0, GENOME_END),
            breaks = axis_breaks,
            labels = axis_labels,
            expand = c(0, 0)
          ) +
          # provisional; final after segments
          scale_y_continuous(limits = c(0, y_max * 1.5), expand = c(0, 0)) +
          labs(x = "Genome", y = "-log10(FDR) [EWAS bins]") +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        p2_pl <- ggplotly(p2, tooltip = "text", source = src_combo)
      }
    }
    
    # ----------------------------
    # CLUSTER SEGMENT BAND on p2 (TOP when no bins, adaptive when bins)
    # ----------------------------
    y_axis_max <- if (is.finite(y_max) && y_max > 0) (y_max * 1.5) else 1.5
    
    if (is.data.frame(rv$clusters) && nrow(rv$clusters)) {
      
      clseg <- as.data.frame(rv$clusters) %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          chr_num    = suppressWarnings(as.integer(chr)),
          start_i    = suppressWarnings(as.numeric(start)),
          end_i      = suppressWarnings(as.numeric(end))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          is.finite(chr_num), is.finite(start_i), is.finite(end_i),
          end_i >= start_i
        ) %>%
        dplyr::distinct(cluster_id, chr_num, start_i, end_i) %>%
        dplyr::left_join(
          ref %>% dplyr::transmute(chr_num = as.integer(chrN), chr_cum = as.numeric(chr_cum)),
          by = "chr_num"
        ) %>%
        dplyr::filter(is.finite(chr_cum)) %>%
        dplyr::transmute(
          cluster_id = cluster_id,
          x0   = pmax(0, start_i + chr_cum),
          x1   = pmin(GENOME_END, end_i + chr_cum),
          xmid = (x0 + x1) / 2,
          text = paste0(
            "Cluster: ", cluster_id,
            "<br>chr", chr_num, ":", format(start_i, scientific = FALSE),
            "-", format(end_i, scientific = FALSE)
          )
        ) %>%
        dplyr::filter(is.finite(x0), is.finite(x1), x1 >= x0) %>%
        dplyr::arrange(x0)
      
      if (nrow(clseg) > 0) {
        
        # y_ref: top of bins if present; else a safe top based on thresholds
        thr_guess <- suppressWarnings(as.numeric(thr_y))
        if (!is.finite(thr_guess) || thr_guess <= 0) thr_guess <- 8
        
        y_ref <- if (isTRUE(has_bins) && is.finite(y_max) && y_max > 0) y_max else max(8, thr_guess, 10)
        
        bump   <- max(0.6, 0.06 * y_ref)
        y_seg  <- y_ref + bump
        y_tick <- max(0.15, 0.02 * y_ref)
        y_txt  <- y_seg + max(0.4, 0.05 * y_ref)
        
        clseg$y_seg  <- y_seg
        clseg$y0tick <- y_seg - y_tick
        clseg$y1tick <- y_seg + y_tick
        clseg$y_txt  <- y_txt
        
        # ensure y-axis includes segment + text band
        y_axis_max <- max(y_axis_max, y_txt + max(0.4, 0.05 * y_ref))
        
        p2_pl <- p2_pl %>%
          plotly::add_segments(
            data = clseg,
            x = ~x0, xend = ~x1,
            y = ~y_seg, yend = ~y_seg,
            key = ~cluster_id,
            inherit = FALSE,
            line = list(width = 3),
            hoverinfo = "text", text = ~text,
            showlegend = FALSE
          ) %>%
          plotly::add_segments(
            data = clseg,
            x = ~x0, xend = ~x0,
            y = ~y0tick, yend = ~y1tick,
            inherit = FALSE, line = list(width = 2),
            hoverinfo = "none", showlegend = FALSE
          ) %>%
          plotly::add_segments(
            data = clseg,
            x = ~x1, xend = ~x1,
            y = ~y0tick, yend = ~y1tick,
            inherit = FALSE, line = list(width = 2),
            hoverinfo = "none", showlegend = FALSE
          ) %>%
          plotly::add_annotations(
            data = clseg,
            x = ~xmid, y = ~y_txt,
            text = ~cluster_id,
            xref = "x", yref = "y",
            showarrow = FALSE,
            textangle = 45,
            font = list(size = 10),
            xanchor = "center",
            yanchor = "bottom"
          )
      }
    }
    
    # enforce final y-range AFTER adding segments (works both with/without bins)
    p2_pl <- p2_pl %>%
      plotly::layout(yaxis = list(range = c(0, y_axis_max)))
    
    out <- plotly::subplot(
      p1_pl, p2_pl,
      nrows = 2, shareX = TRUE,
      heights = c(0.55, 0.45),
      titleY = TRUE
    )
    
    out$x$source <- src_combo
    out <- plotly::event_register(out, "plotly_click")
    out <- plotly::event_register(out, "plotly_relayout")
    
    out <- out %>%
      plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d","lasso2d","hoverCompareCartesian")
      ) %>%
      plotly::layout(showlegend = FALSE)
    
    out
  })
  
  # -----------------------------
  # Update UCSC region on zoom
  # -----------------------------
  observeEvent(plotly::event_data("plotly_relayout", source = "manhattan_combo"), {
    d <- plotly::event_data("plotly_relayout", source = "manhattan_combo")
    xr <- get_relayout_xrange(d)
    if (is.null(xr)) { combo_region(NULL); return() }
    
    max_width <- 50e6
    if (!is.finite(xr[1]) || !is.finite(xr[2])) return()
    if ((xr[2] - xr[1]) > max_width) return()
    
    combo_region(list(x0 = xr[1], x1 = xr[2]))
    reg <- cumrange_to_ucsc_region(xr[1], xr[2], .ref_hg38)
    if (!is.null(reg)) ucsc_region(reg)
  }, ignoreInit = TRUE)
  
  # ----------------------------
  # Track builders for UCSC
  # ----------------------------
  parse_ucsc_region <- function(region) {
    if (is.null(region) || !nzchar(region)) return(NULL)
    m <- regexec("^chr([^:]+):([0-9,]+)-([0-9,]+)$", region)
    r <- regmatches(region, m)[[1]]
    if (length(r) != 4) return(NULL)
    chr <- toupper(r[2])
    st  <- suppressWarnings(as.integer(gsub(",", "", r[3])))
    en  <- suppressWarnings(as.integer(gsub(",", "", r[4])))
    if (!is.finite(st) || !is.finite(en)) return(NULL)
    if (en < st) { tmp <- st; st <- en; en <- tmp }
    list(chr = chr, start = st, end = en)
  }
  
  clean_track <- function(df) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(tibble::tibble())
    df %>%
      dplyr::mutate(
        start = pmax(0L, as.integer(start)),
        end   = pmax(start + 1L, as.integer(end))
      ) %>%
      dplyr::distinct(chrom, start, end, name, .keep_all = TRUE)
  }
  
  cap_track_rows <- function(df, max_rows = 5000L) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)
    max_rows <- as.integer(max_rows)
    if (!is.finite(max_rows) || max_rows <= 0) return(df[0, , drop = FALSE])
    if (nrow(df) <= max_rows) return(df)
    df[seq_len(max_rows), , drop = FALSE]
  }
  
  normalize_hits_for_track <- function(hits, ref) {
    if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) return(tibble::tibble())
    
    chr_col  <- pick_col(hits, c("CHR","chr","CHR_std","chrom","chromosome"))
    bp_col   <- pick_col(hits, c("BP","POS","POS_std","position","POSITION","bp","pos"))
    bpc_col  <- pick_col(hits, c("BPcum","pos_cum","bp_cum","poscum","BP_CUM","POS_CUM"))
    name_col <- pick_col(hits, c("snp","SNP","rsid","RSID","marker","ID","name"))
    
    ref2 <- ref %>%
      dplyr::transmute(
        chr     = toupper(as.character(chr)),
        chr_cum = suppressWarnings(as.numeric(chr_cum)),
        len     = suppressWarnings(as.numeric(len))
      ) %>%
      dplyr::filter(!is.na(chr), is.finite(chr_cum), is.finite(len), len > 0) %>%
      dplyr::arrange(chr_cum)
    
    out <- hits %>%
      dplyr::transmute(
        CHR_raw = if (!is.null(chr_col)) as.character(.data[[chr_col]]) else NA_character_,
        BP_raw  = if (!is.null(bp_col))  suppressWarnings(readr::parse_number(as.character(.data[[bp_col]]))) else NA_real_,
        BPcum   = if (!is.null(bpc_col)) suppressWarnings(as.numeric(.data[[bpc_col]])) else NA_real_,
        name    = if (!is.null(name_col)) as.character(.data[[name_col]]) else NA_character_
      )
    
    if (all(!is.finite(out$BP_raw)) && any(is.finite(out$BPcum))) {
      out <- out %>%
        dplyr::mutate(
          chr_idx = vapply(BPcum, function(v) {
            if (!is.finite(v)) return(NA_integer_)
            idx <- which(v >= ref2$chr_cum & v <= (ref2$chr_cum + ref2$len))
            if (length(idx)) idx[1] else NA_integer_
          }, integer(1)),
          CHR = ifelse(is.na(chr_idx), NA_character_, ref2$chr[chr_idx]),
          BP  = ifelse(is.na(chr_idx), NA_real_, BPcum - ref2$chr_cum[chr_idx])
        )
    } else {
      out <- out %>%
        dplyr::mutate(
          CHR = toupper(trimws(CHR_raw)),
          CHR = gsub("^CHR", "", CHR, ignore.case = TRUE),
          CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
          CHR = dplyr::case_when(
            CHR %in% as.character(1:22) ~ CHR,
            CHR %in% c("23","X") ~ "X",
            CHR %in% c("24","Y") ~ "Y",
            CHR %in% c("MT","M","25") ~ "MT",
            TRUE ~ NA_character_
          ),
          BP = BP_raw
        )
    }
    
    out %>%
      dplyr::mutate(
        BP   = suppressWarnings(as.integer(BP)),
        name = dplyr::if_else(!is.na(name) & nzchar(name), name, paste0("chr", CHR, ":", BP))
      ) %>%
      dplyr::filter(!is.na(CHR), is.finite(BP), !is.na(name), nzchar(name)) %>%
      dplyr::distinct(CHR, BP, name, .keep_all = TRUE)
  }
  
  make_track_df_hits <- function(hits, rgb_r, rgb_g, rgb_b, region_info = NULL, ref = .ref_hg38) {
    nh <- normalize_hits_for_track(hits, ref)
    if (!nrow(nh)) return(tibble::tibble())
    
    if (!is.null(region_info)) {
      chr_need <- toupper(region_info$chr)
      nh <- nh %>% dplyr::filter(CHR == chr_need, BP >= region_info$start, BP <= region_info$end)
      if (!nrow(nh)) return(tibble::tibble())
    }
    
    nh %>%
      dplyr::transmute(
        chrom = paste0("chr", CHR),
        start = as.integer(BP) - 1L,
        end   = as.integer(BP),
        name  = as.character(name),
        score = 0,
        strand = "+",
        thickStart = start,
        thickEnd   = end,
        itemRgb = paste0(rgb_r, ",", rgb_g, ",", rgb_b),
        blockCount = 1,
        blockSizes = 1,
        blockStarts = 0
      ) %>%
      dplyr::filter(is.finite(start), is.finite(end))
  }
  
  make_track_df_cpg <- function(coord, rgb_r, rgb_g, rgb_b,
                                region_info = NULL,
                                max_rows = 2000L) {
    
    if (is.null(region_info)) return(tibble::tibble())
    if (is.null(coord) || !is.data.frame(coord) || !nrow(coord)) return(tibble::tibble())
    if (!all(c("probe","chr","pos") %in% names(coord))) return(tibble::tibble())
    
    dt <- coord %>%
      dplyr::transmute(
        probe = as.character(probe),
        chr0  = as.character(chr),
        pos0  = suppressWarnings(as.integer(pos))
      ) %>%
      dplyr::filter(!is.na(probe), nzchar(probe), is.finite(pos0)) %>%
      dplyr::mutate(
        CHR = toupper(trimws(chr0)),
        CHR = gsub("^CHR", "", CHR, ignore.case = TRUE),
        CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
        CHR = dplyr::case_when(
          CHR %in% as.character(1:22) ~ CHR,
          CHR %in% c("23","X") ~ "X",
          CHR %in% c("24","Y") ~ "Y",
          CHR %in% c("MT","M","25") ~ "MT",
          TRUE ~ NA_character_
        ),
        BP = pos0
      ) %>%
      dplyr::filter(!is.na(CHR), is.finite(BP))
    
    chr_need <- toupper(region_info$chr)
    dt <- dt %>% dplyr::filter(CHR == chr_need, BP >= region_info$start, BP <= region_info$end)
    if (!nrow(dt)) return(tibble::tibble())
    
    if (nrow(dt) > max_rows) dt <- dt[order(dt$BP), ][unique(round(seq(1, nrow(dt), length.out = max_rows))), , drop = FALSE]
    
    dt %>%
      dplyr::transmute(
        chrom = paste0("chr", CHR),
        start = as.integer(BP) - 1L,
        end   = as.integer(BP),
        name  = probe,
        score = 0,
        strand = "+",
        thickStart = start,
        thickEnd   = end,
        itemRgb = paste0(rgb_r, ",", rgb_g, ",", rgb_b),
        blockCount = 1,
        blockSizes = 1,
        blockStarts = 0
      ) %>%
      dplyr::filter(is.finite(start), is.finite(end)) %>%
      dplyr::distinct(chrom, start, end, name, .keep_all = TRUE)
  }
  
  make_track_df_bins <- function(bins, rgb_r, rgb_g, rgb_b,
                                 region_info = NULL,
                                 coord_hg38 = NULL,
                                 coord_hg19 = NULL) {
    
    if (is.null(bins) || !is.data.frame(bins) || !nrow(bins)) return(tibble::tibble())
    bins <- as.data.frame(bins, stringsAsFactors = FALSE)
    
    chr_col <- pick_col(bins, c("chr","CHR","chrom"))
    st_col  <- pick_col(bins, c("bin_start","start","BIN_START","start_bp"))
    en_col  <- pick_col(bins, c("bin_end","end","BIN_END","end_bp"))
    if (is.null(chr_col) || is.null(st_col) || is.null(en_col)) return(tibble::tibble())
    
    dt <- bins %>%
      dplyr::transmute(
        chr0 = as.character(.data[[chr_col]]),
        st0  = suppressWarnings(as.integer(.data[[st_col]])),
        en0  = suppressWarnings(as.integer(.data[[en_col]]))
      ) %>%
      dplyr::filter(!is.na(chr0), is.finite(st0), is.finite(en0)) %>%
      dplyr::mutate(
        CHR = toupper(trimws(chr0)),
        CHR = gsub("^CHR", "", CHR, ignore.case = TRUE),
        CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
        CHR = dplyr::case_when(
          CHR %in% as.character(1:22) ~ CHR,
          CHR %in% c("23","X") ~ "X",
          CHR %in% c("24","Y") ~ "Y",
          CHR %in% c("MT","M","25") ~ "MT",
          TRUE ~ NA_character_
        ),
        start = pmin(st0, en0),
        end   = pmax(st0, en0)
      ) %>%
      dplyr::filter(!is.na(CHR), is.finite(start), is.finite(end), end > start)
    
    if (!is.null(region_info)) {
      chr_need <- toupper(region_info$chr)
      dt <- dt %>% dplyr::filter(CHR == chr_need, end >= region_info$start, start <= region_info$end)
    }
    if (!nrow(dt)) return(tibble::tibble())
    
    dt %>%
      dplyr::transmute(
        chrom = paste0("chr", CHR),
        start = as.integer(start),
        end   = as.integer(end),
        name  = paste0("bin:", start, "-", end),
        score = 0,
        strand = "+",
        thickStart = start,
        thickEnd   = end,
        itemRgb = paste0(rgb_r, ",", rgb_g, ",", rgb_b),
        blockCount = 1,
        blockSizes = pmax(1L, end - start),
        blockStarts = 0
      ) %>%
      dplyr::distinct(chrom, start, end, name, .keep_all = TRUE)
  }
  
  # Update session track data filtered to ucsc_region()
  observe({
    region <- ucsc_region()
    rinfo  <- parse_ucsc_region(region)
    
    if (is.null(rinfo)) {
      session$userData$track_gwas_data     <- tibble::tibble()
      session$userData$track_cpg_data      <- tibble::tibble()
      session$userData$track_ewasbin_data  <- tibble::tibble()
      return()
    }
    
    gwas_hits <- tryCatch(hits_df(), error = function(e) NULL)
    if (is.null(gwas_hits) || !is.data.frame(gwas_hits) || !nrow(gwas_hits)) {
      gwas_hits <- tryCatch(gwas_df(), error = function(e) NULL)
    }
    
    bins <- NULL
    if (is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all)) bins <- rv$ewas_bins_all
    
    df_gwas <- clean_track(make_track_df_hits(gwas_hits, 31,120,180, region_info = rinfo, ref = .ref_hg38))
    df_bins <- clean_track(make_track_df_bins(bins, 227,26,28, region_info = rinfo, coord_hg38 = rv$coord, coord_hg19 = coord_hg19))
    df_cpg  <- clean_track(make_track_df_cpg(rv$coord, 128,0,128, region_info = rinfo, max_rows = 2000L))
    
    df_cpg  <- cap_track_rows(df_cpg,  max_rows = 5000L)
    df_bins <- cap_track_rows(df_bins, max_rows = 8000L)
    df_gwas <- cap_track_rows(df_gwas, max_rows = 8000L)
    
    session$userData$track_gwas_data      <- df_gwas
    session$userData$track_cpg_data       <- df_cpg
    session$userData$track_ewasbin_data   <- df_bins
  })
  
  make_ucsc_track_text <- function(name, df, url_tpl) {
    name2 <- gsub("[^A-Za-z0-9_\\-\\.]", "_", as.character(name))
    header <- paste(
      "track",
      paste0("name=", name2),
      paste0("description=", name2),
      "visibility=pack",
      "itemRgb=On",
      paste0("url=", url_tpl)
    )
    
    if (is.null(df) || !nrow(df)) return(header)
    
    cols <- c(
      "chrom","start","end","name","score","strand",
      "thickStart","thickEnd","itemRgb",
      "blockCount","blockSizes","blockStarts"
    )
    df2 <- as.data.frame(df, stringsAsFactors = FALSE)
    miss <- setdiff(cols, names(df2))
    if (length(miss)) return(header)
    
    df2 <- df2[, cols, drop = FALSE]
    df2$chrom <- gsub("[\t\r\n ]+", "", df2$chrom)
    df2$name  <- gsub("[\t\r\n]", " ", df2$name)
    
    body <- apply(df2, 1, paste, collapse = "\t")
    paste(c(header, body), collapse = "\n")
  }
  
  make_ucsc_url <- function(region, track_text) {
    base <- "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38"
    reg_enc <- utils::URLencode(region, reserved = TRUE)
    trk_enc <- utils::URLencode(track_text, reserved = TRUE)
    trk_enc <- gsub("\\+", "%20", trk_enc)
    paste0(base, "&position=", reg_enc, "&hgt.customText=", trk_enc)
  }
  
  # --- Add these helpers ONCE (near the top of app.R / before server) --------
  
  as_df_or_empty <- function(x) {
    if (is.null(x)) return(tibble::tibble())
    if (is.data.frame(x)) return(x)
    tibble::tibble()
  }
  
  safe_nrow <- function(x) {
    if (is.null(x)) return(0L)
    if (is.data.frame(x)) return(nrow(x))
    0L
  }
  
  # --- Replace your EWAS disease UCSC blocks with this (copy/paste) ----------
  
  output$ucsc_link_gwas <- renderUI({
    region <- ucsc_region()
    req(region)
    
    df  <- as_df_or_empty(session$userData$track_gwas_data)
    trk <- make_ucsc_track_text("GWAS_hits", df, "https://www.ncbi.nlm.nih.gov/snp/$$")
    url <- make_ucsc_url(region, trk)
    
    tags$a(href = url, target = "_blank", "Open UCSC – GWAS hits (dbSNP links)")
  })
  
  output$ucsc_link_cpg <- renderUI({
    region <- ucsc_region()
    req(region)
    
    df  <- as_df_or_empty(session$userData$track_cpg_data)
    trk <- make_ucsc_track_text("CpG_probes", df, "https://ngdc.cncb.ac.cn/ewas/datahub/probe/$$")
    url <- make_ucsc_url(region, trk)
    
    tags$a(href = url, target = "_blank", "Open UCSC – CpGs in interval (EWAS Datahub)")
  })
  
  output$ucsc_link_ewasbins <- renderUI({
    region <- ucsc_region()
    req(region)
    
    df  <- as_df_or_empty(session$userData$track_ewasbin_data)
    trk <- make_ucsc_track_text("EWAS_sig_bins", df, "")
    url <- make_ucsc_url(region, trk)
    
    tags$a(href = url, target = "_blank", "Open UCSC – EWAS bins")
  })
  
  output$debug_ucsc_state <- renderUI({
    region <- ucsc_region() %||% "NULL"
    
    gwas_n <- safe_nrow(session$userData$track_gwas_data)
    cpg_n  <- safe_nrow(session$userData$track_cpg_data)
    bin_n  <- safe_nrow(session$userData$track_ewasbin_data)
    
    combo <- combo_region()
    
    tags$pre(
      style="background-color:#f6f6f6; border:1px solid #ddd; padding:10px; font-family: 'Courier New', monospace;",
      paste0(
        "region = ", region,
        "\nGWAS hits in track = ", gwas_n,
        "\nCpGs in track = ", cpg_n,
        "\nEWAS bins in track = ", bin_n,
        "\ncombo_region = ", if (is.null(combo)) "NULL" else paste0(combo$x0, " - ", combo$x1)
      )
    )
  })
  
  ########### EWAS hits table (DISEASE version: case/control) ####################
  # ============================================================================
  # EWAS HITS TABLE (CpG-level) — Build from ALL clusters with sig bins
  # Triggered ONLY by input$build_ewas_hits
  # DISEASE MODE: grp = "case" vs "control" (NO tumor/adjacent normal)
  # ============================================================================
  #
  # Assumptions (already exist in your app):
  # - rv$clusters        : data.frame with cluster_id, chr, start, end
  # - rv$ewas_bins_all   : data.frame with cluster_id, bin_start, bin_end, and "diseases" (CSV) (or similar name)
  # - rv$meta            : data.table/data.frame with sample_id, disease, grp ("case"/"control")  [1:1 no repetition already applied]
  # - rv$coord           : data.table with chr ("chr1"...), pos, probe
  # - extract_subset_if_needed(chr_disease, st, en) exists and returns subset file path (.tsv or .tsv.gz)
  # - chr_label_plink() exists
  # - rv$ewas_subset_cache is an env (create it if missing)
  # - append_log() exists
  #
  # Outputs produced here:
  # - rv$ewas_hits_all   : CpG hits table for ALL clusters
  # - output$tbl_ewas_hits : DT table with clickable disease/probe/hyper_hypo
  
  # -----------------------------
  # Helper: from bins_df -> list(disease -> unique(bin_mid))
  # bins_df needs: bin_start, bin_end, diseases (CSV string)
  # -----------------------------
  sig_bins_by_disease_from_bins <- function(bins_df) {
    if (is.null(bins_df) || !nrow(bins_df)) return(list())
    
    bins_df <- as.data.frame(bins_df, stringsAsFactors = FALSE)
    
    bins_df$bin_start <- suppressWarnings(as.numeric(bins_df$bin_start))
    bins_df$bin_end   <- suppressWarnings(as.numeric(bins_df$bin_end))
    bins_df$bin_mid   <- (bins_df$bin_start + bins_df$bin_end) / 2
    
    out <- list()
    for (i in seq_len(nrow(bins_df))) {
      ds <- trimws(unlist(strsplit(as.character(bins_df$diseases[i]), ",")))
      ds <- ds[nzchar(ds)]
      if (!length(ds)) next
      bm <- as.numeric(bins_df$bin_mid[i])
      for (d in ds) out[[d]] <- unique(c(out[[d]], bm))
    }
    out
  }
  
  # -----------------------------
  # Helper: extract/read subset for ONE cluster (cached)
  # returns subset file path
  # -----------------------------
  
  get_subset_for_cluster <- function(cluster_id, chr_disease, st, en,
                                     cache_env = ewas_subset_cache) {
    key <- as.character(cluster_id)
    
    if (exists(key, envir = cache_env, inherits = FALSE)) {
      p <- get(key, envir = cache_env, inherits = FALSE)
      if (!is.null(p) && file.exists(p)) return(p)
    }
    
    sub_full <- extract_subset_if_needed(chr_disease, st, en)
    validate(need(!is.null(sub_full) && file.exists(sub_full), "Subset extraction failed."))
    
    assign(key, sub_full, envir = cache_env)
    sub_full
  }
  
  # -----------------------------
  # Helper: read ONE probe row from subset (fast via awk)
  # -----------------------------
  read_one_probe_row <- function(sub_full, probe, cols_keep = NULL) {
    validate(need(file.exists(sub_full), "Subset file not found."))
    
    is_gz <- grepl("\\.gz$", sub_full, ignore.case = TRUE)
    catbin <- if (is_gz) {
      gz <- Sys.which("gzip")
      validate(need(nzchar(gz), "gzip not found on PATH (needed for .gz subset)."))
      paste0(shQuote(gz), " -cd ", shQuote(sub_full))
    } else {
      paste0("cat ", shQuote(sub_full))
    }
    
    awkcmd <- paste0(
      "awk -F'\t' -v p=", shQuote(probe), " 'NR==1{print; next} $1==p{print; exit}'"
    )
    
    cmd <- paste(catbin, awkcmd, sep = " | ")
    
    dt <- data.table::fread(
      cmd = cmd,
      na.strings = c("NA", "<NA>", "NaN", ""),
      check.names = FALSE,
      showProgress = FALSE
    )
    
    validate(need(nrow(dt) == 1, "Probe not found in subset interval."))
    
    if (!is.null(cols_keep)) {
      cols_keep <- intersect(cols_keep, names(dt))
      validate(need(length(cols_keep) >= 2, "Not enough columns found for this CpG."))
      dt <- dt[, ..cols_keep]
    }
    
    dt
  }
  
  # -----------------------------
  # Helper: compute p + BH(FDR) for a set of probes (same subset)
  # -----------------------------
  compute_p_fdr_for_hitset <- function(sub_full, probes, case_ids, ctl_ids,
                                       test_m = c("wilcox", "ttest")) {
    test_m <- match.arg(test_m)
    
    probes <- unique(as.character(probes))
    probes <- probes[nzchar(probes)]
    if (!length(probes)) {
      return(data.table::data.table(probe = character(), p = numeric(), fdr = numeric()))
    }
    
    sel <- unique(c(case_ids, ctl_ids))
    cols_to_read <- c("sample_id", sel)
    
    dt <- data.table::fread(
      sub_full,
      select = cols_to_read,
      na.strings = c("NA", "<NA>", "NaN", ""),
      check.names = FALSE,
      showProgress = FALSE
    )
    
    dt <- dt[sample_id %in% probes]
    if (!nrow(dt)) {
      return(data.table::data.table(probe = character(), p = numeric(), fdr = numeric()))
    }
    
    dt[, (sel) := lapply(.SD, suppressWarnings(as.numeric)), .SDcols = sel]
    
    pvals <- apply(dt[, ..sel], 1, function(v) {
      x <- v[sel %in% case_ids]
      y <- v[sel %in% ctl_ids]
      
      x <- x[is.finite(x) & x >= 0 & x <= 1]
      y <- y[is.finite(y) & y >= 0 & y <= 1]
      if (length(x) < 2 || length(y) < 2) return(NA_real_)
      
      if (test_m == "ttest") {
        tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
      } else {
        tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      }
    })
    
    out <- data.table::data.table(
      probe = as.character(dt$sample_id),
      p     = as.numeric(pvals)
    )
    out[, fdr := p.adjust(p, method = "BH")]
    out
  }
  
  # -----------------------------
  # Core: build EWAS CpG hits for ONE cluster (from its bins)
  # DISEASE MODE: loop diseases, case/control from meta (grp)
  # -----------------------------
  build_ewas_hits_one_clusterXXXXXXXX <- function(cl_row, bins_cluster) {
    req(rv$meta, rv$coord)
    
    cl_row <- as.data.frame(cl_row, stringsAsFactors = FALSE)
    validate(need(nrow(cl_row) == 1, "Internal: cluster row must be length 1."))
    
    cluster_id <- as.character(cl_row$cluster_id[1])
    chr_disease     <- paste0("chr", chr_label_plink(as.integer(cl_row$chr[1])))
    st0        <- as.integer(cl_row$start[1])
    en0        <- as.integer(cl_row$end[1])
    
    # list disease -> sig bin_mid
    sig_bins_list <- sig_bins_by_disease_from_bins(bins_cluster)
    if (!length(sig_bins_list)) return(NULL)
    
    # subset file for this cluster (cached)
    sub_full <- get_subset_for_cluster(cluster_id, chr_disease, st0, en0)
    validate(need(file.exists(sub_full), "Subset file not found after extraction."))
    
    # coord interval
    coord2 <- rv$coord[chr == chr_disease & pos >= st0 & pos <= en0, .(probe, pos)]
    if (!nrow(coord2)) return(NULL)
    
    bin_size <- as.integer(input$bin_size %||% 5000L)
    bin_size <- max(1000L, bin_size)
    
    meta <- as.data.table(rv$meta)
    meta[, sample_id := as.character(sample_id)]
    meta[, disease   := as.character(disease)]
    meta[, grp       := as.character(grp)]
    
    # For efficiency: read subset ONCE with union of sample cols across diseases present here
    dis_names <- names(sig_bins_list)
    
    dis_info <- lapply(dis_names, function(d) {
      case_ids <- meta[disease == d & grp == "case", unique(sample_id)]
      ctl_ids  <- meta[disease == d & grp == "control", unique(sample_id)]
      list(disease = d, case_ids = case_ids, ctl_ids = ctl_ids)
    })
    # keep only diseases with enough samples
    dis_info <- Filter(function(x) length(x$case_ids) >= 2 && length(x$ctl_ids) >= 2, dis_info)
    if (!length(dis_info)) return(NULL)
    
    sel_all <- unique(unlist(lapply(dis_info, function(x) c(x$case_ids, x$ctl_ids))))
 
    
    dt <- read_subset_agnostic(rv$ewas_sub_full, sel_all)
    
    beta_dt <- dt[grepl("^cg\\d{8}$", sample_id)]
    if (!nrow(beta_dt)) return(NULL)
    
    beta_dt[, (sel_all) := lapply(.SD, suppressWarnings(as.numeric)), .SDcols = sel_all]
    beta_mat <- as.matrix(beta_dt[, ..sel_all])
    rownames(beta_mat) <- beta_dt$sample_id
    
    # positions
    pos <- coord2$pos[match(rownames(beta_mat), coord2$probe)]
    keep <- is.finite(pos)
    if (!any(keep)) return(NULL)
    
    beta_mat <- beta_mat[keep, , drop = FALSE]
    pos      <- as.integer(pos[keep])
    
    # assign CpGs to bins (same scheme as violin)
    bin_start <- st0 + ((pos - st0) %/% bin_size) * bin_size
    bin_end   <- bin_start + bin_size
    bin_mid   <- bin_start + bin_size / 2
    
    hits_all <- list()
    
    for (info in dis_info) {
      d <- info$disease
      sig_bins <- sig_bins_list[[d]]
      if (!length(sig_bins)) next
      
      case_ids <- info$case_ids
      ctl_ids  <- info$ctl_ids
      
      case_cols <- intersect(colnames(beta_mat), case_ids)
      ctl_cols  <- intersect(colnames(beta_mat), ctl_ids)
      if (length(case_cols) < 2 || length(ctl_cols) < 2) next
      
      mean_case <- rowMeans(beta_mat[, case_cols, drop = FALSE], na.rm = TRUE)
      mean_ctl  <- rowMeans(beta_mat[, ctl_cols,  drop = FALSE], na.rm = TRUE)
      delta     <- mean_case - mean_ctl
      
      keep_bin <- bin_mid %in% sig_bins
      if (!any(keep_bin)) next
      
      res <- data.table::data.table(
        cluster_id = cluster_id,
        chr        = chr_disease,
        start      = st0,
        end        = en0,
        bin_start  = as.integer(bin_start[keep_bin]),
        bin_end    = as.integer(bin_end[keep_bin]),
        bin_mid    = as.numeric(bin_mid[keep_bin]),
        pos        = as.integer(pos[keep_bin]),
        probe      = rownames(beta_mat)[keep_bin],
        disease    = as.character(d),
        n_case     = as.integer(length(case_cols)),
        n_ctl      = as.integer(length(ctl_cols)),
        mean_beta_case    = as.numeric(mean_case[keep_bin]),
        mean_beta_control = as.numeric(mean_ctl[keep_bin]),
        delta_beta = as.numeric(delta[keep_bin]),
        hyper_hypo = data.table::fifelse(delta[keep_bin] > 0, "hyper",
                                         data.table::fifelse(delta[keep_bin] < 0, "hypo", "no_change"))
      )
      
      hits_all[[d]] <- res
    }
    
    if (!length(hits_all)) return(NULL)
    data.table::rbindlist(hits_all, use.names = TRUE, fill = TRUE)
  }
  
  build_ewas_hits_one_cluster <- function(cl_row, bins_cluster) {
    req(rv$meta, rv$coord)
    
    cl_row <- as.data.frame(cl_row, stringsAsFactors = FALSE)
    validate(need(nrow(cl_row) == 1, "Internal: cluster row must be length 1."))
    
    cluster_id  <- as.character(cl_row$cluster_id[1])
    chr_disease <- paste0("chr", chr_label_plink(as.integer(cl_row$chr[1])))
    st0         <- as.integer(cl_row$start[1])
    en0         <- as.integer(cl_row$end[1])
    
    # list disease -> sig bin_mid
    sig_bins_list <- sig_bins_by_disease_from_bins(bins_cluster)
    if (!length(sig_bins_list)) return(NULL)
    
    # subset file for this cluster (cached)
    sub_full <- get_subset_for_cluster(cluster_id, chr_disease, st0, en0)
    
    # --- DEBUG + robust guards (AQUI) ---
    cat("[DEBUG] build_ewas_hits_one_cluster:",
        "cluster_id=", cluster_id,
        " chr=", chr_disease,
        " start=", st0,
        " end=", en0, "\n")
    cat("[DEBUG] sub_full =", paste0(capture.output(str(sub_full)), collapse=" "), "\n")
    
    # (si tens must_path1() global, usa-la; si no, substitueix per validate/stop)
    sub_full <- must_path1(sub_full, "subset_file")
    if (!file.exists(sub_full)) stop("[subset_file] not found: ", sub_full, call. = FALSE)
    
    # coord interval
    coord2 <- rv$coord[chr == chr_disease & pos >= st0 & pos <= en0, .(probe, pos)]
    if (!nrow(coord2)) return(NULL)
    
    bin_size <- as.integer(input$bin_size %||% 5000L)
    bin_size <- max(1000L, bin_size)
    
    meta <- as.data.table(rv$meta)
    meta[, sample_id := as.character(sample_id)]
    meta[, disease   := as.character(disease)]
    meta[, grp       := as.character(grp)]
    
    # For efficiency: read subset ONCE with union of sample cols across diseases present here
    dis_names <- names(sig_bins_list)
    
    dis_info <- lapply(dis_names, function(d) {
      case_ids <- meta[disease == d & grp == "case", unique(sample_id)]
      ctl_ids  <- meta[disease == d & grp == "control", unique(sample_id)]
      list(disease = d, case_ids = case_ids, ctl_ids = ctl_ids)
    })
    # keep only diseases with enough samples
    dis_info <- Filter(function(x) length(x$case_ids) >= 2 && length(x$ctl_ids) >= 2, dis_info)
    if (!length(dis_info)) return(NULL)
    
    sel_all <- unique(unlist(lapply(dis_info, function(x) c(x$case_ids, x$ctl_ids))))
    sel_all <- sel_all[nzchar(sel_all)]
    if (!length(sel_all)) return(NULL)
    
    # >>> IMPORTANT FIX: use sub_full (NOT rv$ewas_sub_full) <<<
    dt <- read_subset_agnostic(sub_full, sel_all)
    
    beta_dt <- dt[grepl("^cg\\d{8}$", sample_id)]
    if (!nrow(beta_dt)) return(NULL)
    
    beta_dt[, (sel_all) := lapply(.SD, suppressWarnings(as.numeric)), .SDcols = sel_all]
    beta_mat <- as.matrix(beta_dt[, ..sel_all])
    rownames(beta_mat) <- beta_dt$sample_id
    
    # positions
    pos <- coord2$pos[match(rownames(beta_mat), coord2$probe)]
    keep <- is.finite(pos)
    if (!any(keep)) return(NULL)
    
    beta_mat <- beta_mat[keep, , drop = FALSE]
    pos      <- as.integer(pos[keep])
    
    # assign CpGs to bins (same scheme as violin)
    bin_start <- st0 + ((pos - st0) %/% bin_size) * bin_size
    bin_end   <- bin_start + bin_size
    bin_mid   <- bin_start + bin_size / 2
    
    hits_all <- list()
    
    for (info in dis_info) {
      d <- info$disease
      sig_bins <- sig_bins_list[[d]]
      if (!length(sig_bins)) next
      
      case_ids <- info$case_ids
      ctl_ids  <- info$ctl_ids
      
      case_cols <- intersect(colnames(beta_mat), case_ids)
      ctl_cols  <- intersect(colnames(beta_mat), ctl_ids)
      if (length(case_cols) < 2 || length(ctl_cols) < 2) next
      
      mean_case <- rowMeans(beta_mat[, case_cols, drop = FALSE], na.rm = TRUE)
      mean_ctl  <- rowMeans(beta_mat[, ctl_cols,  drop = FALSE], na.rm = TRUE)
      delta     <- mean_case - mean_ctl
      
      keep_bin <- bin_mid %in% sig_bins
      if (!any(keep_bin)) next
      
      res <- data.table::data.table(
        cluster_id = cluster_id,
        chr        = chr_disease,
        start      = st0,
        end        = en0,
        bin_start  = as.integer(bin_start[keep_bin]),
        bin_end    = as.integer(bin_end[keep_bin]),
        bin_mid    = as.numeric(bin_mid[keep_bin]),
        pos        = as.integer(pos[keep_bin]),
        probe      = rownames(beta_mat)[keep_bin],
        disease    = as.character(d),
        n_case     = as.integer(length(case_cols)),
        n_ctl      = as.integer(length(ctl_cols)),
        mean_beta_case    = as.numeric(mean_case[keep_bin]),
        mean_beta_control = as.numeric(mean_ctl[keep_bin]),
        delta_beta = as.numeric(delta[keep_bin]),
        hyper_hypo = data.table::fifelse(
          delta[keep_bin] > 0, "hyper",
          data.table::fifelse(delta[keep_bin] < 0, "hypo", "no_change")
        )
      )
      
      hits_all[[d]] <- res
    }
    
    if (!length(hits_all)) return(NULL)
    data.table::rbindlist(hits_all, use.names = TRUE, fill = TRUE)
  }
  
  # -----------------------------
  # Build EWAS hits for ALL clusters that have bins with diseases
  # Uses rv$clusters + rv$ewas_bins_all (NOT selected_cluster())
  # -----------------------------
  build_ewas_hits_all_clusters <- function() {
    req(rv$clusters, rv$ewas_bins_all, rv$meta, rv$coord)
    
    cl <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
    bins_all <- as.data.frame(rv$ewas_bins_all, stringsAsFactors = FALSE)
    
    validate(need(nrow(cl) > 0, "No clusters available."))
    validate(need(nrow(bins_all) > 0, "No EWAS bins available. Run ALL clusters first."))
    
    validate(need("cluster_id" %in% names(bins_all),
                  paste0("bins_all missing 'cluster_id'. Available: ", paste(names(bins_all), collapse = ", "))))
    
    # Standardize the column holding CSV diseases
    if (!"diseases" %in% names(bins_all)) {
      alt <- intersect(c("disease", "sig_disease", "sig_diseases"), names(bins_all))
      if (length(alt) >= 1) {
        names(bins_all)[names(bins_all) == alt[1]] <- "diseases"
      } else {
        bins_all$diseases <- rep(NA_character_, nrow(bins_all))
      }
    }
    
    bins_all$diseases <- as.character(bins_all$diseases)
    bins_all <- bins_all[!is.na(bins_all$diseases) & nzchar(trimws(bins_all$diseases)), , drop = FALSE]
    validate(need(nrow(bins_all) > 0, "No bins have diseases assigned."))
    
    cids <- sort(unique(as.character(bins_all$cluster_id)))
    validate(need(length(cids) > 0, "No clusters with significant bins/diseases."))
    
    cl <- cl[as.character(cl$cluster_id) %in% cids, , drop = FALSE]
    validate(need(nrow(cl) > 0, "Clusters list does not match bins cluster_id."))
    
    out <- list()
    for (i in seq_len(nrow(cl))) {
      cid <- as.character(cl$cluster_id[i])
      bcl <- bins_all[as.character(bins_all$cluster_id) == cid, , drop = FALSE]
      if (!nrow(bcl)) next
      
      ht <- build_ewas_hits_one_cluster(cl[i, , drop = FALSE], bcl)
      if (!is.null(ht) && nrow(ht)) out[[cid]] <- ht
    }
    
    if (!length(out)) return(NULL)
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }
  
  # -----------------------------
  # Trigger: Build EWAS hits (ALL clusters)
  # -----------------------------
  observeEvent(input$build_ewas_hits, {
    append_log("[EWAS hits] building CpG table for ALL clusters with sig bins (DISEASE: case/control)...")
    
    withProgress(message = "Building EWAS hits (CpG-level)…", value = 0, {
      incProgress(0.10, detail = "Checking inputs…")
      incProgress(0.35, detail = "Scanning clusters & computing deltas…")
      
      ht <- build_ewas_hits_all_clusters()
      
      incProgress(0.90, detail = "Preparing table…")
      rv$ewas_hits_all <- ht
      
      incProgress(1)
    })
  }, ignoreInit = TRUE)
  
  # -----------------------------
  # Table renderer (DISEASE version)
  # - click disease/probe/hyper_hypo
  # -----------------------------
  
  output$ewas_cluster_filter_ui <- renderUI({
    hits <- rv$ewas_hits_all
    choices <- if (is.null(hits) || !nrow(hits)) character(0) else sort(unique(as.character(hits$cluster_id)))
    selectizeInput(
      "ewas_cluster_filter", "Filter cluster",
      choices = c("All" = "", choices),
      selected = "", multiple = FALSE,
      options = list(placeholder = "Type to search cluster_id...", allowEmptyOption = TRUE)
    )
  })
  
  output$ewas_disease_filter_ui <- renderUI({
    hits <- rv$ewas_hits_all
    choices <- if (is.null(hits) || !nrow(hits)) character(0) else sort(unique(as.character(hits$disease)))
    selectizeInput(
      "ewas_disease_filter", "Filter disease",
      choices = c("All" = "", choices),
      selected = "", multiple = FALSE,
      options = list(placeholder = "Type to search disease...", allowEmptyOption = TRUE)
    )
  })
  
  output$ewas_probe_filter_ui <- renderUI({
    hits <- rv$ewas_hits_all
    choices <- if (is.null(hits) || !nrow(hits)) character(0) else sort(unique(as.character(hits$probe)))
    selectizeInput(
      "ewas_probe_filter", "Filter probe",
      choices = c("All" = "", choices),
      selected = "", multiple = FALSE,
      options = list(placeholder = "Type to search probe...", allowEmptyOption = TRUE)
    )
  })
  
  
  # =========================
  # Violin design (shared)
  # =========================
  ewas_violin_scales <- function() {
    list(
      ggplot2::scale_fill_manual(values = c(
        "case"   = "#F39C12",  # cases (taronja)
        "control" = "darkgreen"   # controls (verd fosc)
      )),
      ggplot2::guides(fill = "none")
    )
  }
  
  ewas_violin_theme <- function() {
    ggplot2::theme_bw() +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 9),
        axis.text.x = ggplot2::element_text(angle = 0),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  
  output$tbl_ewas_hits <- DT::renderDT({
    
    hits <- rv$ewas_hits_all
    
    if (is.null(hits) || !nrow(hits)) {
      return(DT::datatable(
        data.frame(Message = "No EWAS hits yet. Click 'Build EWAS hits'."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # -----------------------------
    # Filters
    # -----------------------------
    cl_f <- as.character(input$ewas_cluster_filter %||% "")
    d_f  <- as.character(input$ewas_disease_filter %||% "")
    p_f  <- as.character(input$ewas_probe_filter   %||% "")
    q_f  <- as.character(input$ewas_free_search    %||% "")
    
    hits2 <- hits
    if (nzchar(cl_f)) hits2 <- hits2[as.character(hits2$cluster_id) == cl_f, ]
    if (nzchar(d_f))  hits2 <- hits2[as.character(hits2$disease)    == d_f,  ]
    if (nzchar(p_f))  hits2 <- hits2[as.character(hits2$probe)      == p_f,  ]
    
    if (nzchar(trimws(q_f))) {
      qq  <- tolower(trimws(q_f))
      hay <- tolower(paste(hits2$cluster_id, hits2$disease, hits2$probe, hits2$chr, hits2$pos))
      hits2 <- hits2[grepl(qq, hay, fixed = TRUE), ]
    }
    
    hits <- hits2
    
    if (is.null(hits) || !nrow(hits)) {
      return(DT::datatable(
        data.frame(Message = "No rows after filters."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    df <- as.data.frame(hits, stringsAsFactors = FALSE)
    
    # -----------------------------
    # Raw columns (hidden) to keep exact values
    # -----------------------------
    df$probe_raw   <- as.character(df$probe)
    df$disease_raw <- as.character(df$disease)
    df$hh_raw      <- as.character(df$hyper_hypo)
    
    # IMPORTANT: URL-encode disease so it is safe inside HTML attributes (handles ' and unicode)
    df$disease_enc <- vapply(df$disease_raw, function(x) utils::URLencode(x, reserved = TRUE), character(1))
    
    # -----------------------------
    # Clickable PROBE (unchanged pattern)
    # -----------------------------
    df$probe_html <- mapply(function(probe, cluster_id) {
      paste0(
        "<span class='probePick' ",
        "data-probe='", htmltools::htmlEscape(probe), "' ",
        "data-cluster='", htmltools::htmlEscape(cluster_id), "'>",
        htmltools::htmlEscape(probe),
        "</span>"
      )
    }, df$probe_raw, df$cluster_id, USE.NAMES = FALSE)
    
    # -----------------------------
    # Clickable DISEASE (SAFE: uses data-disease-enc)
    # -----------------------------
    df$disease_html <- mapply(function(disease, cluster_id) {
      
      dis <- htmltools::htmlEscape(as.character(disease))
      cl  <- htmltools::htmlEscape(as.character(cluster_id))
      
      paste0(
        "<span class=\"diseasePickHit\" ",
        "data-disease=\"", dis, "\" ",
        "data-cluster=\"", cl, "\">",
        dis,
        "</span>"
      )
    }, df$disease_raw, df$cluster_id, USE.NAMES = FALSE)
    
    # -----------------------------
    # Clickable hyper/hypo (SAFE: uses data-disease-enc too)
    # -----------------------------
    df$hyper_hypo_html <- mapply(function(hh, probe, disease_enc, cluster_id) {
      paste0(
        "<span class='hhPick' ",
        "data-probe='", htmltools::htmlEscape(probe), "' ",
        "data-disease-enc='", disease_enc, "' ",
        "data-cluster='", htmltools::htmlEscape(cluster_id), "'>",
        htmltools::htmlEscape(hh),
        "</span>"
      )
    }, df$hh_raw, df$probe_raw, df$disease_enc, df$cluster_id, USE.NAMES = FALSE)
    
    
    # --- DISEASE clickable (IMPORTANT: to open violins modal) ---
    df$disease_html <- mapply(function(disease, cluster_id) {
      paste0(
        "<span class='diseasePickHit' ",
        "data-disease='", htmltools::htmlEscape(disease), "' ",
        "data-cluster='", htmltools::htmlEscape(cluster_id), "'>",
        htmltools::htmlEscape(disease),
        "</span>"
      )
    }, df$disease_raw, df$cluster_id, USE.NAMES = FALSE)
    
    # --- hyper/hypo clickable ---
    df$hyper_hypo_html <- mapply(function(hh, probe, disease, cluster_id) {
      paste0(
        "<span class='hhPick' ",
        "data-probe='",   htmltools::htmlEscape(probe),   "' ",
        "data-disease='", htmltools::htmlEscape(disease), "' ",
        "data-cluster='", htmltools::htmlEscape(cluster_id), "'>",
        htmltools::htmlEscape(hh),
        "</span>"
      )
    }, df$hh_raw, df$probe_raw, df$disease_raw, df$cluster_id, USE.NAMES = FALSE)
    
    # numèrics (com tenies)
    df$mean_beta_case    <- round(as.numeric(df$mean_beta_case), 4)
    df$mean_beta_control <- round(as.numeric(df$mean_beta_control), 4)
    df$delta_beta        <- round(as.numeric(df$delta_beta), 4)
    
    show_df <- df %>%
      dplyr::transmute(
        cluster_id,
        chr,
        
        disease_raw,          # hidden
        disease = disease_html,
        
        probe_raw,            # hidden
        probe = probe_html,
        
        pos,
        mean_beta_case,
        mean_beta_control,
        delta_beta,
        
        hh_raw,               # hidden
        hyper_hypo = hyper_hypo_html,
        
        n_case,
        n_ctl
      )
    
    col_probe_raw   <- which(names(show_df) == "probe_raw") - 1
    col_hh_raw      <- which(names(show_df) == "hh_raw") - 1
    col_disease_raw <- which(names(show_df) == "disease_raw") - 1
    
    DT::datatable(
      show_df,
      escape = FALSE,
      rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        columnDefs = list(
          list(targets = c(col_probe_raw, col_hh_raw, col_disease_raw), visible = FALSE)
        )
      ),
 callback = DT::JS(sprintf("
  // column indices (0-based) for hidden raw columns
  var COL_PROBE_RAW   = %d;
  var COL_DISEASE_RAW = %d;
  var COL_HH_RAW      = %d;

  // avoid duplicates on re-render
  table.off('click', 'span.hhPick');
  table.off('click', 'span.probePick');
  table.off('click', 'span.diseasePickHit');

  // helper: get row data array from clicked element
  function giRowData(el){
    var tr = $(el).closest('tr');
    return table.row(tr).data(); // array of cell values
  }

  // DISEASE click -> disease_hit_pick (modal violins)
  table.on('click', 'span.diseasePickHit', function(e){
    e.preventDefault(); e.stopPropagation();

    var row = giRowData(this);
    if(!row) return;

    // disease from hidden raw column (NOT from data-* attribute)
    var disease = row[COL_DISEASE_RAW];
    var cluster = $(this).data('cluster'); // cluster_id is safe

    Shiny.setInputValue('disease_hit_pick',
      {disease: disease, cluster_id: cluster},
      {priority: 'event'}
    );
  });

  // PROBE click -> probe_pick
  table.on('click', 'span.probePick', function(e){
    e.preventDefault(); e.stopPropagation();

    var row = giRowData(this);
    if(!row) return;

    var probe   = row[COL_PROBE_RAW];       // from hidden raw column
    var cluster = $(this).data('cluster');

    Shiny.setInputValue('probe_pick',
      {probe: probe, cluster_id: cluster},
      {priority: 'event'}
    );
  });

  // hyper/hypo click -> hh_pick
  table.on('click', 'span.hhPick', function(e){
    e.preventDefault(); e.stopPropagation();

    var row = giRowData(this);
    if(!row) return;

    var probe   = row[COL_PROBE_RAW];       // raw probe
    var disease = row[COL_DISEASE_RAW];     // raw disease (works with ')
    var cluster = $(this).data('cluster');

    Shiny.setInputValue('hh_pick',
      {probe: probe, disease: disease, cluster_id: cluster},
      {priority: 'event'}
    );
  });

", col_probe_raw, col_disease_raw, col_hh_raw))
      
    ) %>%
      DT::formatStyle("probe",      cursor = "pointer", textDecoration = "underline", fontWeight = "bold") %>%
      DT::formatStyle("disease",    cursor = "pointer", textDecoration = "underline", fontWeight = "bold") %>%
      DT::formatStyle("hyper_hypo", cursor = "pointer", textDecoration = "underline", fontWeight = "bold") %>%
      DT::formatStyle(
        "hh_raw",
        target = "row",
        backgroundColor = DT::styleEqual(c("hyper", "hypo"), c("#ffd6d6", "#d9fdd9"))
      )
  })
  
  
  # end table ewas hits
  
  norm_quote <- function(x) {
    x <- enc2utf8(trimws(as.character(x)))
    x <- gsub("&apos;|&#39;|&#x27;|&#x2019;|&#8217;", "'", x, ignore.case = TRUE)
    x <- gsub("[\u2018\u2019\u02BC\u0060\u00B4]", "'", x, perl = TRUE)
    x
  }
  
  # ============== Downloaders =================================================
  
  # --- download buffers (ONE TIME) ---
  rv$dl_hit_hh_txt      <- reactiveVal(NULL)
  rv$dl_hit_disease_txt <- reactiveVal(NULL)
  rv$dl_hit_probe_txt   <- reactiveVal(NULL)
  
  output$dl_hit_hh <- downloadHandler(
    filename = function() paste0("ewasdis_hit_hh_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      txt <- rv$dl_hit_hh_txt()
      validate(need(!is.null(txt) && nzchar(txt), "Click a hyper/hypo item first to generate the text file."))
      writeLines(txt, file, useBytes = TRUE)
    }
  )
  
  output$dl_hit_disease <- downloadHandler(
    filename = function() paste0("ewasdis_hit_disease_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      txt <- rv$dl_hit_disease_txt()
      validate(need(!is.null(txt) && nzchar(txt), "Click a disease item first to generate the text file."))
      writeLines(txt, file, useBytes = TRUE)
    }
  )
  
  output$dl_hit_probe <- downloadHandler(
    filename = function() paste0("ewasdis_hit_probe_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      txt <- rv$dl_hit_probe_txt()
      validate(need(!is.null(txt) && nzchar(txt), "Click a probe item first to generate the text file."))
      writeLines(txt, file, useBytes = TRUE)
    }
  )
  
 
  # =============================================================================
  # hh_pick -> modal violin for ONE CpG within (cluster + disease), with p + FDR
  # + Download TXT
  # =============================================================================
  rv$cpg_violin_dt      <- NULL
  rv$cpg_violin_probe   <- NULL
  rv$cpg_violin_disease <- NULL
  rv$cpg_violin_cluster <- NULL
  rv$cpg_violin_test    <- NULL
  rv$cpg_violin_p       <- NULL
  rv$cpg_violin_fdr     <- NULL
  
  observeEvent(input$hh_pick, {
    req(input$hh_pick$probe, input$hh_pick$disease, input$hh_pick$cluster_id)
    req(rv$meta)
    validate(need(is.data.frame(rv$clusters) && nrow(rv$clusters) > 0, "No clusters loaded."))
    
    probe       <- as.character(input$hh_pick$probe)
    disease_sel <- as.character(input$hh_pick$disease)
    cluster_id  <- as.character(input$hh_pick$cluster_id)
    
    # 1) cluster chr/start/end
    cl0 <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
    cl0 <- cl0[as.character(cl0$cluster_id) == cluster_id, , drop = FALSE]
    validate(need(nrow(cl0) == 1, "Cluster not found for this hit."))
    
    chr_disease <- paste0("chr", chr_label_plink(as.integer(cl0$chr[1])))
    st <- as.integer(cl0$start[1])
    en <- as.integer(cl0$end[1])
    
    # 2) Subset (cache)
    sub_full <- get_subset_for_cluster(cluster_id, chr_disease, st, en)
    validate(need(!is.null(sub_full) && length(sub_full) == 1 && nzchar(sub_full), "Subset path is NULL/empty."))
    validate(need(file.exists(sub_full), "Subset file not found."))
    
    # 3) samples case/control for this disease (meta grp)
    meta <- as.data.table(rv$meta)
    meta[, sample_id := as.character(sample_id)]
    meta[, disease   := as.character(disease)]
    meta[, grp       := as.character(grp)]
    
    meta_d <- unique(
      meta[disease == disease_sel & grp %in% c("case", "control"), .(sample_id, grp)],
      by = "sample_id"
    )
    validate(need(nrow(meta_d) > 0, "No samples for this disease in meta."))
    
    case_ids <- meta_d[grp == "case",    as.character(sample_id)]
    ctl_ids  <- meta_d[grp == "control", as.character(sample_id)]
    case_ids <- case_ids[nzchar(case_ids)]
    ctl_ids  <- ctl_ids[nzchar(ctl_ids)]
    validate(need(length(case_ids) >= 2, "Not enough case samples."))
    validate(need(length(ctl_ids)  >= 2, "Not enough control samples."))
    
    sel <- unique(c(case_ids, ctl_ids))
    
    # 4) Read ONE CpG row only
    cols_keep <- c("sample_id", sel)
    row_dt <- read_one_probe_row(sub_full, probe, cols_keep = cols_keep)
    validate(need(is.data.frame(row_dt) && nrow(row_dt) == 1, "CpG not found in subset."))
    
    sel2 <- intersect(sel, names(row_dt))
    validate(need(length(sel2) >= 4, "Not enough sample columns found for this CpG."))
    
    b <- unlist(as.list(row_dt[1, ..sel2]), use.names = TRUE)
    b <- suppressWarnings(as.numeric(b))
    names(b) <- sel2
    
    long <- data.table::data.table(
      sample_id = sel2,
      beta = as.numeric(b),
      grp  = data.table::fifelse(sel2 %in% case_ids, "case", "control")
    )
    long <- long[is.finite(beta) & beta >= 0 & beta <= 1]
    long[, grp := factor(grp, levels = c("control", "case"))]
    
    validate(need(
      nrow(long[grp == "case"]) >= 2 && nrow(long[grp == "control"]) >= 2,
      "Not enough beta values to draw violins."
    ))
    
    # 5) p-value + BH(FDR) within (cluster + disease) hitset
    test_m <- as.character(input$bin_test_any %||% "wilcox")
    
    x <- long[grp == "case", beta]
    y <- long[grp == "control", beta]
    p_one <- if (test_m == "ttest") {
      tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
    } else {
      tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    }
    
    fdr_one <- NA_real_
    hits <- rv$ewas_hits_all
    if (is.data.frame(hits) && nrow(hits)) {
      pr_col <- if ("probe_raw" %in% names(hits)) "probe_raw" else "probe"
      probes_hitset <- unique(as.character(hits[[pr_col]][hits$cluster_id == cluster_id & hits$disease == disease_sel]))
      probes_hitset <- probes_hitset[nzchar(probes_hitset)]
      if (length(probes_hitset) > 0 && length(probes_hitset) <= 2000) {
        pfdr <- compute_p_fdr_for_hitset(sub_full, probes_hitset, case_ids, ctl_ids, test_m = test_m)
        if (is.data.frame(pfdr) && nrow(pfdr)) fdr_one <- pfdr$fdr[match(probe, pfdr$probe)]
      }
    }
    
    # save state for plot
    rv$cpg_violin_dt      <- long
    rv$cpg_violin_probe   <- probe
    rv$cpg_violin_disease <- disease_sel
    rv$cpg_violin_cluster <- cluster_id
    rv$cpg_violin_test    <- test_m
    rv$cpg_violin_p       <- p_one
    rv$cpg_violin_fdr     <- fdr_one
    
    # --- build TXT for download ---
    one_df <- data.frame(
      cluster_id = cluster_id,
      disease = disease_sel,
      probe = probe,
      chr = chr_disease,
      start = st,
      end = en,
      test = test_m,
      p = p_one,
      fdr = fdr_one,
      n_case = sum(long$grp == "case"),
      n_ctl  = sum(long$grp == "control"),
      mean_case = mean(long$beta[long$grp == "case"], na.rm = TRUE),
      mean_ctl  = mean(long$beta[long$grp == "control"], na.rm = TRUE),
      delta = mean(long$beta[long$grp == "case"], na.rm = TRUE) - mean(long$beta[long$grp == "control"], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    hdr <- c(
      "TYPE\tEWASDIS_HH",
      paste0("CLUSTER_ID\t", cluster_id),
      paste0("DISEASE\t", disease_sel),
      paste0("PROBE\t", probe),
      paste0("CHR\t", chr_disease),
      paste0("START\t", st),
      paste0("END\t", en),
      paste0("TEST\t", test_m),
      paste0("DATE\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ""
    )
    
    tab_one  <- capture.output(write.table(one_df, sep = "\t", row.names = FALSE, quote = FALSE))
    tab_long <- capture.output(write.table(as.data.frame(long), sep = "\t", row.names = FALSE, quote = FALSE))
    
    rv$dl_hit_hh_txt(paste(c(
      hdr,
      "ONE_ROW_STATS",
      tab_one,
      "",
      "LONG_BETAS",
      tab_long,
      ""
    ), collapse = "\n"))
    
    showModal(modalDialog(
      title = paste0("CpG ", probe, " — ", disease_sel, " (Cluster ", cluster_id, ")"),
      size = "l",
      easyClose = TRUE,
      plotOutput("p_violin_cpg", height = "420px"),
      footer = tagList(
        downloadButton("dl_hit_hh", "Download hyper/hypo stats (txt)"),
        modalButton("Close")
      )
    ))
  }, ignoreInit = TRUE)
  
  # =============================================================================
  # disease_hit_pick -> modal “all CpG violins” for (cluster + disease)
  # (EWAS DIS) + Download TXT
  # =============================================================================
  rv$disease_violin_long    <- NULL
  rv$disease_violin_stats   <- NULL
  rv$disease_violin_disease <- NULL
  rv$disease_violin_cluster <- NULL
  rv$disease_violin_test    <- NULL
  rv$disease_violin_total_n <- NULL
  rv$disease_violin_top_n   <- NULL
  
  observeEvent(input$disease_hit_pick, {
    req(input$disease_hit_pick$disease, input$disease_hit_pick$cluster_id)
    req(rv$meta)
    validate(need(is.data.frame(rv$clusters) && nrow(rv$clusters) > 0, "No clusters loaded."))
    validate(need(is.data.frame(rv$ewas_hits_all) && nrow(rv$ewas_hits_all) > 0, "No EWAS hits loaded."))
    
    disease_sel <- as.character(input$disease_hit_pick$disease)
    cluster_id  <- as.character(input$disease_hit_pick$cluster_id)
    
    # 1) cluster chr/start/end
    cl0 <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
    cl0 <- cl0[as.character(cl0$cluster_id) == cluster_id, , drop = FALSE]
    validate(need(nrow(cl0) == 1, "Cluster not found."))
    
    chr_disease <- paste0("chr", chr_label_plink(as.integer(cl0$chr[1])))
    st <- as.integer(cl0$start[1])
    en <- as.integer(cl0$end[1])
    
    # 2) subset (cache)
    sub_full <- get_subset_for_cluster(cluster_id, chr_disease, st, en)
    validate(need(!is.null(sub_full) && length(sub_full) == 1 && nzchar(sub_full), "Subset path is NULL/empty."))
    validate(need(file.exists(sub_full), "Subset file not found."))
    
    # 3) probes d’aquest (cluster, disease)
    hits <- as.data.frame(rv$ewas_hits_all, stringsAsFactors = FALSE)
    pr_col <- if ("probe_raw" %in% names(hits)) "probe_raw" else "probe"
    
    probes <- unique(as.character(hits[[pr_col]][hits$cluster_id == cluster_id & hits$disease == disease_sel]))
    probes <- probes[nzchar(probes)]
    validate(need(length(probes) > 0, "No CpG hits for this disease in this cluster."))
    
    max_probes <- 60L
    if (length(probes) > max_probes) probes <- probes[1:max_probes]
    
    # 4) case/control ids pel disease (EWASdis meta: grp)
    meta <- as.data.table(rv$meta)
    meta[, sample_id := as.character(sample_id)]
    meta[, disease   := as.character(disease)]
    meta[, grp       := as.character(grp)]
    
    meta_d <- unique(meta[disease == disease_sel & grp %in% c("case","control"), .(sample_id, grp)], by="sample_id")
    validate(need(nrow(meta_d) > 0, "No samples for this disease in meta."))
    
    case_ids <- meta_d[grp == "case",    as.character(sample_id)]
    ctl_ids  <- meta_d[grp == "control", as.character(sample_id)]
    case_ids <- case_ids[nzchar(case_ids)]
    ctl_ids  <- ctl_ids[nzchar(ctl_ids)]
    validate(need(length(case_ids) >= 2, "Not enough case samples."))
    validate(need(length(ctl_ids)  >= 2, "Not enough control samples."))
    
    sel <- unique(c(case_ids, ctl_ids))
    cols_to_read <- c("sample_id", sel)
    
    dt <- data.table::fread(
      sub_full,
      select = cols_to_read,
      na.strings = c("NA","<NA>","NaN",""),
      check.names = FALSE,
      showProgress = FALSE
    )
    validate(need(nrow(dt) > 0, "No CpGs found in subset for this cluster interval."))
    
    data.table::setnames(dt, "sample_id", "probe")
    
    dt <- dt[probe %in% probes]
    validate(need(nrow(dt) > 0, "None of the CpG hits were found inside the subset file."))
    
    test_m <- as.character(input$bin_test_any %||% "wilcox")
    
    # stats per CpG (p + BH/FDR dins aquest set)
    stats_all <- dt[, {
      x <- as.numeric(unlist(.SD[, ..case_ids], use.names = FALSE))
      y <- as.numeric(unlist(.SD[, ..ctl_ids],  use.names = FALSE))
      x <- x[is.finite(x) & x >= 0 & x <= 1]
      y <- y[is.finite(y) & y >= 0 & y <= 1]
      
      p <- if (length(x) >= 2 && length(y) >= 2) {
        if (test_m == "ttest") {
          tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
        } else {
          tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
        }
      } else NA_real_
      
      .(p = p)
    }, by = .(probe)]
    
    stats_all <- stats_all[is.finite(p)]
    validate(need(nrow(stats_all) > 0, "No valid p-values for these CpGs (after filtering)."))
    
    stats_all[, fdr := p.adjust(p, method = "BH")]
    stats_all <- stats_all[is.finite(fdr)]
    validate(need(nrow(stats_all) > 0, "No valid FDR values for these CpGs (after filtering)."))
    
    stats_all[, probe_lab := paste0(
      probe,
      "\n p=", format.pval(p, digits = 2, eps = 1e-12),
      " | FDR=", format.pval(fdr, digits = 2, eps = 1e-12)
    )]
    
    max_plot <- 16L
    data.table::setorder(stats_all, fdr, p)
    stats_plot <- stats_all[1:min(max_plot, .N)]
    probes_plot <- stats_plot$probe
    validate(need(length(probes_plot) > 0, "No CpGs available for plotting."))
    
    dt_plot <- dt[probe %in% probes_plot, c("probe", sel), with = FALSE]
    validate(need(nrow(dt_plot) > 0, "No data rows for plotting after selecting top CpGs."))
    
    long <- data.table::melt(
      dt_plot,
      id.vars = "probe",
      variable.name = "sample_id",
      value.name = "beta"
    )
    long[, beta := as.numeric(beta)]
    long <- long[is.finite(beta) & beta >= 0 & beta <= 1]
    validate(need(nrow(long) > 0, "No finite beta values to plot."))
    
    long[, grp := data.table::fifelse(sample_id %in% case_ids, "case", "control")]
    long[, grp := factor(grp, levels = c("control","case"))]
    
    long <- merge(long, stats_plot[, .(probe, probe_lab)], by = "probe", all.x = TRUE)
    levs <- stats_plot$probe_lab
    long[, probe_lab := factor(probe_lab, levels = levs)]
    
    # save
    rv$disease_violin_long    <- long
    rv$disease_violin_stats   <- stats_all
    rv$disease_violin_test    <- test_m
    rv$disease_violin_total_n <- nrow(stats_all)
    rv$disease_violin_top_n   <- nrow(stats_plot)
    rv$disease_violin_disease <- disease_sel
    rv$disease_violin_cluster <- cluster_id
    
    # --- build TXT for download ---
    hdr <- c(
      "TYPE\tEWASDIS_DISEASE_ALLCPG",
      paste0("CLUSTER_ID\t", cluster_id),
      paste0("DISEASE\t", disease_sel),
      paste0("CHR\t", chr_disease),
      paste0("START\t", st),
      paste0("END\t", en),
      paste0("TEST\t", test_m),
      paste0("N_PROBES_TOTAL\t", nrow(stats_all)),
      paste0("TOP_N_PLOTTED\t", nrow(stats_plot)),
      paste0("DATE\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ""
    )
    
    tab_stats <- capture.output(write.table(as.data.frame(stats_all), sep="\t", row.names=FALSE, quote=FALSE))
    tab_long  <- capture.output(write.table(as.data.frame(long),      sep="\t", row.names=FALSE, quote=FALSE))
    
    rv$dl_hit_disease_txt(paste(c(
      hdr,
      "STATS_TABLE",
      tab_stats,
      "",
      "LONG_TABLE_TOPN",
      tab_long,
      ""
    ), collapse="\n"))
    
    showModal(modalDialog(
      title = paste0("Disease: ", disease_sel, " (Cluster ", cluster_id, ")"),
      size = "l",
      easyClose = TRUE,
      plotOutput("p_violin_disease_all", height = "520px"),
      footer = tagList(
        downloadButton("dl_hit_disease", "Download disease stats (txt)"),
        modalButton("Close")
      )
    ))
  }, ignoreInit = TRUE)
  
  # =============================================================================
  # probe_pick -> modal: CpG distributions across DISEASES in this cluster
  # (EWAS DIS terminology; FDR across diseases for this CpG) + Download TXT
  # =============================================================================
  rv$probe_cluster_long    <- NULL
  rv$probe_cluster_stats   <- NULL
  rv$probe_cluster_probe   <- NULL
  rv$probe_cluster_cluster <- NULL
  rv$probe_cluster_test    <- NULL
  rv$probe_cluster_total_n <- NULL
  rv$probe_cluster_top_n   <- NULL
  
  observeEvent(input$probe_pick, {
    req(input$probe_pick$probe, input$probe_pick$cluster_id)
    req(rv$meta)
    validate(need(is.data.frame(rv$clusters) && nrow(rv$clusters) > 0, "No clusters loaded."))
    
    probe      <- as.character(input$probe_pick$probe)
    cluster_id <- as.character(input$probe_pick$cluster_id)
    
    cl0 <- as.data.frame(rv$clusters, stringsAsFactors = FALSE)
    cl0 <- cl0[as.character(cl0$cluster_id) == cluster_id, , drop = FALSE]
    validate(need(nrow(cl0) == 1, "Cluster not found."))
    
    chr_disease <- paste0("chr", chr_label_plink(as.integer(cl0$chr[1])))
    st <- as.integer(cl0$start[1])
    en <- as.integer(cl0$end[1])
    
    sub_full <- get_subset_for_cluster(cluster_id, chr_disease, st, en)
    validate(need(!is.null(sub_full) && length(sub_full) == 1 && nzchar(sub_full), "Subset path is NULL/empty."))
    validate(need(file.exists(sub_full), "Subset file not found."))
    
    meta <- as.data.table(rv$meta)
    meta[, sample_id := as.character(sample_id)]
    meta[, disease   := as.character(disease)]
    meta[, grp       := as.character(grp)]
    
    test_m <- as.character(input$bin_test_any %||% "wilcox")
    
    diseases_all <- sort(unique(as.character(meta$disease)))
    diseases_all <- diseases_all[nzchar(diseases_all)]
    validate(need(length(diseases_all) > 0, "No diseases found in meta."))
    
    maxS <- 250L
    disease_info <- lapply(diseases_all, function(di) {
      meta_d <- unique(meta[disease == di & grp %in% c("case","control"), .(sample_id, grp)], by="sample_id")
      if (!nrow(meta_d)) return(NULL)
      
      case_ids <- sort(as.character(meta_d[grp == "case", sample_id]))
      ctl_ids  <- sort(as.character(meta_d[grp == "control", sample_id]))
      case_ids <- case_ids[nzchar(case_ids)]
      ctl_ids  <- ctl_ids[nzchar(ctl_ids)]
      if (length(case_ids) < 2 || length(ctl_ids) < 2) return(NULL)
      
      if (length(case_ids) > maxS) case_ids <- case_ids[1:maxS]
      if (length(ctl_ids)  > maxS) ctl_ids  <- ctl_ids[1:maxS]
      
      list(disease = di, case_ids = case_ids, ctl_ids = ctl_ids)
    })
    disease_info <- Filter(Negate(is.null), disease_info)
    validate(need(length(disease_info) > 0, "No diseases have enough case/control samples."))
    
    row_dt <- read_one_probe_row(sub_full, probe, cols_keep = NULL)
    validate(need(is.data.frame(row_dt) && nrow(row_dt) == 1, "CpG not found in subset."))
    
    all_cols <- setdiff(names(row_dt), "sample_id")
    validate(need(length(all_cols) > 0, "Subset row has no sample columns."))
    
    b <- suppressWarnings(as.numeric(unlist(as.list(row_dt[1, ..all_cols]))))
    names(b) <- all_cols
    
    stats_all <- data.table::rbindlist(lapply(disease_info, function(di) {
      d0   <- di$disease
      case <- intersect(names(b), di$case_ids)
      ctl  <- intersect(names(b), di$ctl_ids)
      if (length(case) < 2 || length(ctl) < 2) return(NULL)
      
      x <- b[case]; y <- b[ctl]
      x <- x[is.finite(x) & x >= 0 & x <= 1]
      y <- y[is.finite(y) & y >= 0 & y <= 1]
      if (length(x) < 2 || length(y) < 2) return(NULL)
      
      p <- if (test_m == "ttest") {
        tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
      } else {
        tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      }
      
      data.table::data.table(disease = d0, n_case = length(x), n_ctl = length(y), p = p)
    }), fill = TRUE)
    
    validate(need(is.data.frame(stats_all) && nrow(stats_all) > 0, "No diseases produced valid tests for this CpG."))
    
    stats_all[, fdr := p.adjust(p, method = "BH")]
    stats_all <- stats_all[is.finite(p) & is.finite(fdr)]
    validate(need(nrow(stats_all) > 0, "No valid p/FDR after filtering."))
    
    stats_all[, disease_lab := paste0(
      disease,
      "\n p=", format.pval(p, digits = 2, eps = 1e-12),
      " | FDR=", format.pval(fdr, digits = 2, eps = 1e-12)
    )]
    
    max_diseases <- 16L
    data.table::setorder(stats_all, fdr, p)
    stats_plot <- stats_all[1:min(max_diseases, .N)]
    
    plot_set <- stats_plot$disease
    info_plot <- Filter(function(x) x$disease %in% plot_set, disease_info)
    
    long_list <- lapply(info_plot, function(di) {
      d0   <- di$disease
      case <- intersect(names(b), di$case_ids)
      ctl  <- intersect(names(b), di$ctl_ids)
      if (length(case) < 2 || length(ctl) < 2) return(NULL)
      
      dt1 <- data.table::data.table(
        disease = d0,
        sample_id = c(case, ctl),
        grp = c(rep("case", length(case)), rep("control", length(ctl)))
      )
      dt1[, beta := as.numeric(b[sample_id])]
      dt1 <- dt1[is.finite(beta) & beta >= 0 & beta <= 1]
      if (!nrow(dt1)) return(NULL)
      dt1
    })
    
    long <- data.table::rbindlist(Filter(Negate(is.null), long_list), use.names = TRUE, fill = TRUE)
    validate(need(is.data.frame(long) && nrow(long) > 0, "No valid beta values to plot."))
    
    long <- merge(long, stats_plot[, .(disease, disease_lab)], by = "disease", all.x = TRUE)
    levs <- stats_plot$disease_lab
    long[, disease_lab := factor(disease_lab, levels = levs)]
    long[, grp := factor(grp, levels = c("control","case"))]
    
    rv$probe_cluster_long    <- long
    rv$probe_cluster_stats   <- stats_all
    rv$probe_cluster_probe   <- probe
    rv$probe_cluster_cluster <- cluster_id
    rv$probe_cluster_test    <- test_m
    rv$probe_cluster_total_n <- nrow(stats_all)
    rv$probe_cluster_top_n   <- nrow(stats_plot)
    
    # --- build TXT for download ---
    hdr <- c(
      "TYPE\tEWASDIS_PROBE_ACROSS_DISEASES",
      paste0("CLUSTER_ID\t", cluster_id),
      paste0("PROBE\t", probe),
      paste0("CHR\t", chr_disease),
      paste0("START\t", st),
      paste0("END\t", en),
      paste0("TEST\t", test_m),
      paste0("N_DISEASES_TOTAL\t", nrow(stats_all)),
      paste0("TOP_N_PLOTTED\t", nrow(stats_plot)),
      paste0("DATE\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ""
    )
    
    tab_stats <- capture.output(write.table(as.data.frame(stats_all), sep="\t", row.names=FALSE, quote=FALSE))
    tab_long  <- capture.output(write.table(as.data.frame(long),      sep="\t", row.names=FALSE, quote=FALSE))
    
    rv$dl_hit_probe_txt(paste(c(
      hdr,
      "STATS_TABLE",
      tab_stats,
      "",
      "LONG_TABLE_TOPN",
      tab_long,
      ""
    ), collapse="\n"))
    
    nd <- length(unique(long$disease))
    h_px <- max(520L, min(2000L, 120L + 40L * nd))
    
    showModal(modalDialog(
      title = paste0("CpG ", probe, " — across diseases (Cluster ", cluster_id, ")"),
      size = "l",
      easyClose = TRUE,
      plotOutput("p_violin_probe_cluster", height = paste0(h_px, "px")),
      footer = tagList(
        downloadButton("dl_hit_probe", "Download probe stats (txt)"),
        modalButton("Close")
      )
    ))
  }, ignoreInit = TRUE)
  
  
  # ----------------- Modal Plots
  # plot after click on cpg column
  output$p_violin_cpg <- renderPlot({
    req(rv$cpg_violin_dt)
    df <- rv$cpg_violin_dt
    
    p_one   <- rv$cpg_violin_p
    fdr_one <- rv$cpg_violin_fdr
    test_m  <- rv$cpg_violin_test %||% "wilcox"
    
    lab_p   <- if (is.finite(p_one)) format(p_one, digits = 3, scientific = TRUE) else "NA"
    lab_fdr <- if (is.finite(fdr_one)) format(fdr_one, digits = 3, scientific = TRUE) else "NA"
    
    ggplot(df, aes(x = grp, y = beta, fill = grp)) +
      geom_violin(trim = TRUE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.5) +
      scale_fill_manual(values = c("control" = "darkgreen", "case" = "#ff7f00")) +
      labs(
        title = paste0("CpG ", rv$cpg_violin_probe, " — ", rv$cpg_violin_disease,
                       " (Cluster ", rv$cpg_violin_cluster, ")"),
        subtitle = paste0("Test: ", test_m, " | p=", lab_p, " | BH(FDR)=", lab_fdr,
                          " (within cluster+disease hitset)"),
        x = NULL, y = "Beta"
      ) +
      theme_minimal()
  })
  
  # plot at click on disease column
  output$p_violin_disease_all <- renderPlot({
    req(rv$disease_violin_long)
    
    long   <- rv$disease_violin_long
    test_m <- rv$disease_violin_test %||% (input$bin_test_any %||% "wilcox")
    
    long[, grp := factor(grp, levels = c("control","case"))]
    ncol_fac <- 4L
    
    ggplot2::ggplot(long, ggplot2::aes(x = grp, y = beta, fill = grp)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      scale_fill_manual(values = c("control" = "darkgreen", "case" = "#ff7f00")) +
      ggplot2::facet_wrap(~ probe_lab, ncol = ncol_fac, scales = "free_y") +
      ewas_violin_scales() +
      ggplot2::labs(
        x = NULL,
        y = "Beta",
        subtitle = paste0(
          "Test: ", test_m,
          " | BH(FDR) computed across ALL CpGs in the cluster (N=",
          rv$disease_violin_total_n %||% NA, "). ",
          "Showing Top ", rv$disease_violin_top_n %||% NA, " by FDR."
        )
      ) +
      ewas_violin_theme()
  })
  
  
  # plot on click at hyper_hypo column
  output$p_violin_probe_cluster <- renderPlot({
    req(rv$probe_cluster_long)
    
    long   <- rv$probe_cluster_long
    test_m <- rv$probe_cluster_test %||% (input$bin_test_any %||% "wilcox")
    
    if (!"disease_lab" %in% names(long)) {
      long$disease_lab <- as.character(long$disease)
    }
    
    long[, grp := factor(grp, levels = c("control","case"))]
    ncol_fac <- 4L
    
    ggplot2::ggplot(long, ggplot2::aes(x = grp, y = beta, fill = grp)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      scale_fill_manual(values = c("control" = "darkgreen", "case" = "#ff7f00")) +
      ggplot2::facet_wrap(~ disease_lab, ncol = ncol_fac, scales = "free_y") +
      ewas_violin_scales() +
      ggplot2::labs(
        x = NULL,
        y = "Beta",
        subtitle = paste0(
          "Test: ", test_m,
          " | BH(FDR) across all eligible diseases (N=",
          rv$probe_cluster_total_n %||% NA, "). Showing Top ",
          rv$probe_cluster_top_n %||% NA, " by FDR."
        )
      ) +
      ewas_violin_theme()
  })
  
  
  outputOptions(output, "p_hist",    suspendWhenHidden = FALSE)
  outputOptions(output, "p_beta_dist", suspendWhenHidden = FALSE)
  outputOptions(output, "p_violin_sel", suspendWhenHidden = FALSE)  # opcional
  
  
  ############################################################################
  ######## EWAS DISEASE — export clusters + EWAS bins (ZIP)
  ############################################################################
  
  norm_chr_int <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")] <- "23"
    x[x %in% c("Y","y")] <- "24"
    x[x %in% c("MT","Mt","mt","M","m")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pick_first_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  # assegura bin_start/bin_end a partir de bin_mid si cal
  ensure_bin_start_end <- function(bins_dt, bin_size) {
    b <- data.table::as.data.table(bins_dt)
    
    if (!"bin_start" %in% names(b) && "start" %in% names(b)) b[, bin_start := as.integer(start)]
    if (!"bin_end"   %in% names(b) && "end"   %in% names(b)) b[, bin_end   := as.integer(end)]
    
    if (!("bin_start" %in% names(b) && "bin_end" %in% names(b))) {
      validate(need("bin_mid" %in% names(b), "EWAS bins must have bin_start/bin_end or bin_mid."))
      bs <- as.integer(bin_size)
      half <- as.integer(floor(bs/2))
      b[, bin_start := as.integer(round(bin_mid)) - half]
      b[, bin_end   := as.integer(bin_start) + bs - 1L]
    }
    
    b[, bin_start := as.integer(bin_start)]
    b[, bin_end   := as.integer(bin_end)]
    b
  }
  
  norm_chr_int <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")] <- "23"
    x[x %in% c("Y","y")] <- "24"
    x[x %in% c("MT","Mt","mt","M","m")] <- "25"
    suppressWarnings(as.integer(x))
  }
  
  pick_first_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  # Assegura bin_start/bin_end si només tens "bin" (index) o si tens start sense end.
  # ADAPTA aquesta funció si ja en tens una a la teva app.
  ensure_bin_start_end <- function(b, bin_size = 5000L) {
    b <- data.table::as.data.table(b)
    
    # si ja hi són, ok
    if ("bin_start" %in% names(b) && "bin_end" %in% names(b)) return(b)
    
    # possibles camps
    stc <- pick_first_col(b, c("bin_start","start","start_bp","pos_ini","from","begin"))
    enc <- pick_first_col(b, c("bin_end","end","end_bp","pos_end","to","stop"))
    
    # si tens start però no end: end = start + bin_size - 1
    if (!is.null(stc) && is.null(enc)) {
      b[, bin_start := suppressWarnings(as.integer(get(stc)))]
      b[, bin_end   := bin_start + as.integer(bin_size) - 1L]
      return(b)
    }
    
    # si tens end però no start: start = end - bin_size + 1
    if (is.null(stc) && !is.null(enc)) {
      b[, bin_end   := suppressWarnings(as.integer(get(enc)))]
      b[, bin_start := pmax(1L, bin_end - as.integer(bin_size) + 1L)]
      return(b)
    }
    
    # si tens els dos però amb altres noms
    if (!is.null(stc) && !is.null(enc)) {
      b[, bin_start := suppressWarnings(as.integer(get(stc)))]
      b[, bin_end   := suppressWarnings(as.integer(get(enc)))]
      return(b)
    }
    
    # fallback: si tens "bin" com a índex de bin (0/1-based), calcula start/end
    binc <- pick_first_col(b, c("bin","bin_id","bin_index"))
    if (!is.null(binc)) {
      bi <- suppressWarnings(as.integer(b[[binc]]))
      bi[is.na(bi)] <- 0L
      # assumeixo bin index 0-based si hi ha zeros
      zero_based <- any(bi == 0L, na.rm = TRUE)
      if (zero_based) {
        b[, bin_start := bi * as.integer(bin_size) + 1L]
      } else {
        b[, bin_start := (bi - 1L) * as.integer(bin_size) + 1L]
      }
      b[, bin_end := bin_start + as.integer(bin_size) - 1L]
      return(b)
    }
    
    stop("ensure_bin_start_end(): cannot derive bin_start/bin_end (missing position/bin columns).")
  }
  
  build_ld_clusters_from_ewas_app <- function(rv) {
    
    cl <- rv$clusters
    shiny::validate(shiny::need(is.data.frame(cl) && nrow(cl) > 0,
                                "No clusters available for LD (rv$clusters empty)."))
    
    cl <- as.data.frame(cl)
    
    # tolerate variants
    if (!"chr"   %in% names(cl) && "CHR"      %in% names(cl)) cl$chr   <- cl$CHR
    if (!"start" %in% names(cl) && "start_bp" %in% names(cl)) cl$start <- cl$start_bp
    if (!"end"   %in% names(cl) && "end_bp"   %in% names(cl)) cl$end   <- cl$end_bp
    
    shiny::validate(shiny::need(all(c("chr","start","end") %in% names(cl)),
                                "rv$clusters must contain chr/start/end (or compatible names)."))
    
    # stable cluster_id
    if ("cluster_chr_n" %in% names(cl) && any(nzchar(as.character(cl$cluster_chr_n)))) {
      cl$cluster_id <- as.character(cl$cluster_chr_n)   # ex: cluster_chr2_1
    } else if ("cluster_id" %in% names(cl) && any(nzchar(as.character(cl$cluster_id)))) {
      cl$cluster_id <- as.character(cl$cluster_id)
    } else if ("cluster" %in% names(cl)) {
      cl$cluster_id <- paste0("cluster_", as.integer(cl$cluster))
    } else {
      cl$cluster_id <- paste0("cluster_", seq_len(nrow(cl)))
    }
    
    cl$chr   <- norm_chr_int(cl$chr)
    cl$start <- suppressWarnings(as.integer(cl$start))
    cl$end   <- suppressWarnings(as.integer(cl$end))
    
    cl <- cl[is.finite(cl$chr) & is.finite(cl$start) & is.finite(cl$end) & cl$end >= cl$start, , drop = FALSE]
    shiny::validate(shiny::need(nrow(cl) > 0, "All clusters became invalid after coercion (chr/start/end)."))
    
    out <- cl[, c("cluster_id","chr","start","end")]
    names(out) <- c("cluster_id","chr","cluster_start","cluster_end")
    out
  }
  
  build_ld_candidates_from_ewas_app <- function(rv) {
    
    # -----------------------------
    # GWAS candidates (opcional)
    # -----------------------------
    h <- tryCatch(hits_df(), error = function(e) NULL)
    
    if (is.null(h) || !nrow(h)) {
      gwas_cand <- data.frame(
        chr=integer(), pos_ini=integer(), pos_end=integer(),
        id_hit=character(), classe=character(), stringsAsFactors = FALSE
      )
    } else {
      h <- as.data.frame(h)
      shiny::validate(shiny::need(all(c("CHR","BP") %in% names(h)),
                                  "hits_df() must contain CHR and BP columns."))
      
      id_col <- pick_first_col(h, c("snp","rsid","RSID","SNP","ID","id_hit","variant_id"))
      if (is.null(id_col)) {
        id_col <- ".id_tmp"
        h[[id_col]] <- paste0("chr", h$CHR, ":", h$BP)
      }
      
      rsid <- as.character(h[[id_col]])
      rsid[is.na(rsid) | !nzchar(rsid)] <- paste0("chr", h$CHR, ":", h$BP)
      
      gwas_cand <- data.frame(
        chr     = norm_chr_int(h$CHR),
        pos_ini = suppressWarnings(as.integer(h$BP)),
        pos_end = suppressWarnings(as.integer(h$BP)),
        id_hit  = rsid,
        classe  = "GWAS_hit",
        stringsAsFactors = FALSE
      )
      gwas_cand <- gwas_cand[is.finite(gwas_cand$chr) & is.finite(gwas_cand$pos_ini), , drop = FALSE]
    }
    
    # -----------------------------
    # EWAS bins (SIGNIFICANTS)
    # source-of-truth: rv$ewas_bins_all
    # -----------------------------
    shiny::validate(shiny::need(is.data.frame(rv$ewas_bins_all) && nrow(rv$ewas_bins_all) > 0,
                                "No EWAS bins available for LD/export. Run 'EWAS-ALL' first."))
    
    bin_sz <- as.integer(input$bin_size %||% 5000L)
    b <- data.table::as.data.table(rv$ewas_bins_all)
    b <- ensure_bin_start_end(b, bin_size = bin_sz)
    
    shiny::validate(shiny::need("cluster_id" %in% names(b), "rv$ewas_bins_all must contain cluster_id."))
    b[, cluster_id := as.character(cluster_id)]
    
    b[, chr_int := {
      if ("cluster_chr" %in% names(b))      norm_chr_int(cluster_chr)
      else if ("chr" %in% names(b))         norm_chr_int(chr)
      else NA_integer_
    }]
    
    b <- b[is.finite(chr_int) & is.finite(bin_start) & is.finite(bin_end) & bin_end >= bin_start]
    shiny::validate(shiny::need(nrow(b) > 0, "All EWAS bins became invalid after coercion (chr/bin_start/bin_end)."))
    
    ewas_cand <- data.frame(
      chr     = as.integer(b$chr_int),
      pos_ini = as.integer(b$bin_start),
      pos_end = as.integer(b$bin_end),
      id_hit  = paste0("EWAS_bin_", b$cluster_id, "_chr", b$chr_int, ":", b$bin_start, "-", b$bin_end),
      classe  = "EWAS_dis_bin",
      stringsAsFactors = FALSE
    )
    ewas_cand <- ewas_cand[is.finite(ewas_cand$chr) & is.finite(ewas_cand$pos_ini) & is.finite(ewas_cand$pos_end) &
                             ewas_cand$pos_end >= ewas_cand$pos_ini & nzchar(ewas_cand$id_hit), , drop = FALSE]
    
    out <- rbind(gwas_cand, ewas_cand)
    out <- out[is.finite(out$chr) & is.finite(out$pos_ini) & is.finite(out$pos_end) &
                 out$pos_end >= out$pos_ini & nzchar(out$id_hit), , drop = FALSE]
    
    shiny::validate(shiny::need(nrow(out) > 0, "No candidates available for LD/export."))
    out
  }
  
  output$dl_candidates_zip <- downloadHandler(
    filename = function() {
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      tag_txt <- tryCatch({
        tg <- make_mode_thr_tag(input$cluster_method %||% "window", input$pthr, input$min_logp)
        paste0(tg$mode_tag, "_thr", tg$thr_txt)
      }, error = function(e) "run")
      
      paste0("ewas_disease_candidates_", tag_txt, "_", stamp, ".zip")
    },
    
    content = function(file) {
      
      # 1) clusters canonical
      cl_can <- build_ld_clusters_from_ewas_app(rv)   # cluster_id chr cluster_start cluster_end
      
      # 2) candidates canonical
      cand_can <- build_ld_candidates_from_ewas_app(rv)  # chr pos_ini pos_end id_hit classe
      
      # 3) keep only clusters that truly have EWAS bins (from rv$ewas_bins_all)
      b <- data.table::as.data.table(rv$ewas_bins_all)
      shiny::validate(shiny::need("cluster_id" %in% names(b), "rv$ewas_bins_all must contain cluster_id."))
      keep_ids <- unique(as.character(b$cluster_id))
      keep_ids <- keep_ids[!is.na(keep_ids) & nzchar(keep_ids)]
      
      cl_sig <- cl_can[as.character(cl_can$cluster_id) %in% keep_ids, , drop = FALSE]
      shiny::validate(shiny::need(nrow(cl_sig) > 0, "No clusters matched EWAS bins (nothing to export)."))
      
      # 4) write + zip
      tmpdir <- tempfile("ewas_export_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      f1 <- file.path(tmpdir, "cluster_ewas_dis.csv")
      f2 <- file.path(tmpdir, "candidate_ewas_dis.csv")
      
      utils::write.csv(cl_sig,   f1, row.names = FALSE, quote = FALSE)
      utils::write.csv(cand_can, f2, row.names = FALSE, quote = FALSE)
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      setwd(tmpdir)
      
      utils::zip(zipfile = file, files = c("cluster_ewas_dis.csv", "candidate_ewas_dis.csv"))
    }
  )
  
  ##############################################################################
  ############################ LD Module   #####################################
  ################# EWAS -> LD module canonical cluster input###################
  
  clusters_r <- reactive({
    build_ld_clusters_from_ewas_app(rv)
  })
  
  candidates_r <- reactive({
    build_ld_candidates_from_ewas_app(rv)
  })
  
  ld_module_server("ld", clusters_r = clusters_r, candidates_r = candidates_r)
  
  
  # -----------------------------
  # reset session
  # -----------------------------
  observeEvent(input$reset_case, {
    # --- (A) Esborra fitxers temporals (subset) si existeixen ---
    if (!is.null(rv$ewas_sub_full) && nzchar(rv$ewas_sub_full) && file.exists(rv$ewas_sub_full)) {
      try(unlink(rv$ewas_sub_full, force = TRUE), silent = TRUE)
    }
    
    # --- (B) Reset caches/envs ---
    rv$ewas_subset_cache <- new.env(parent = emptyenv())
    
    # --- (C) Reset “core” state (plots/taules) ---
    rv$clusters        <- NULL
    rv$ewas_bins_all   <- NULL
    rv$ewas_hits_all   <- NULL
    rv$res             <- NULL
    
    # selections
    rv$disease_sel     <- NULL
    rv$cancer_sel      <- NULL
    
    # interval / subset selection
    rv$ewas_chr        <- NULL
    rv$ewas_st         <- NULL
    rv$ewas_en         <- NULL
    rv$ewas_sub_full   <- NULL
    
    # matrices/long tables used in extra plots
    rv$beta_mat        <- NULL
    rv$meta_c          <- NULL
    
    # modals / violins helpers
    rv$cpg_violin_dt   <- NULL
    rv$cpg_violin_probe <- NULL
    rv$cpg_violin_cancer <- NULL
    rv$cpg_violin_cluster <- NULL
    rv$cpg_violin_p    <- NULL
    rv$cpg_violin_fdr  <- NULL
    
    rv$cancer_violin_long    <- NULL
    rv$cancer_violin_stats   <- NULL
    rv$cancer_violin_cancer  <- NULL
    rv$cancer_violin_cluster <- NULL
    
    rv$probe_cluster_long    <- NULL
    rv$probe_cluster_stats   <- NULL
    rv$probe_cluster_probe   <- NULL
    rv$probe_cluster_cluster <- NULL
    rv$probe_cluster_test    <- NULL
    
    # si tens logs
    # rv$log_text <- ""
    
    # --- (D) Reset inputs (opcional, però molt útil) ---
    # Exemple (canvia IDs pels teus)
    # updateSelectInput(session, "disease_sel", selected = "")
    # updateSelectInput(session, "cancer_sel", selected = "")
    # updateNumericInput(session, "alpha_any", value = 0.05)
    # updateNumericInput(session, "bin_size", value = 5000)
    
    # --- (E) Tanca modals si n'hi ha ---
    removeModal()
    
    # --- (F) Feedback ---
    # showNotification("Reset done. Ready for a new case.", type = "message")
  }, ignoreInit = TRUE)
  
  # -----------------------------
  # Info modals
  # -----------------------------
  observeEvent(input$info_00, {
    txt <- HTML('
<p style="text-align: justify;">
  <b>Input file format [(*) mandatory columns]</b>
</p>

<table style="border-collapse: collapse; width: 100%;" border="1">
  <thead>
    <tr>
      <th>CHR(*)</th>
      <th>SNP(*)</th>
      <th>BP(*)</th>
      <th>A1</th>
      <th>TEST</th>
      <th>NMISS</th>
      <th>OR</th>
      <th>STAT</th>
      <th>P(*)</th>
    </tr>
  </thead>

  <tbody>
    <tr><td>14</td><td>14:33967212:C:A</td><td>33967212</td><td>C</td><td>ADD</td><td>720</td><td>2.03</td><td>4.84</td><td>1.30E-06</td></tr>
    <tr><td>18</td><td>rs12955421</td><td>76002716</td><td>A</td><td>ADD</td><td>720</td><td>2.393</td><td>4.76</td><td>1.93E-06</td></tr>
    <tr><td>3</td><td>rs1145036</td><td>2328334</td><td>G</td><td>ADD</td><td>720</td><td>1.816</td><td>4.62</td><td>3.83E-06</td></tr>
  </tbody>
</table>
')
    showModal(modalDialog(
      title = NULL,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "m",
      tags$div(style="line-height:1.35;", txt)
    ))
  })
  
  observeEvent(input$info_01, {
    txt <- HTML(
      "<b>Clustering methods — quick guide</b><br><br>",
      "<p>This app can generate GWAS <b>clusters</b> using two alternative strategies.</p>",
      "<hr style='margin:10px 0;'>",
      "<b>1) By intervals (hits → flank → merge)</b><br>",
      "<p style='margin-top:6px;'>Filter GWAS hits above threshold, build hit-centered intervals [BP±flank], then merge overlaps.</p>",
      "<hr style='margin:10px 0;'>",
      "<b>2) By hit density (windows → merge)</b><br>",
      "<p style='margin-top:6px;'>Filter hits by min_logp, count in windows, keep windows with ≥ min_hits, merge into clusters.</p>"
    )
    
    showModal(modalDialog(
      title = NULL,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "m",
      tags$div(style = "line-height:1.35;", txt)
    ))
  })
  
  #### ---------------------------------------------------------------------------
  disease_path <- gi_file("EWAS_disease", "ewas_detail_dis_genome.rds")
  
  mod_ewas_enrich_server(
    "enrich",
    obs_events = reactive(rv$ewas_detail_all),
    app_mode   = reactive("disease"),
    ref_paths  = ref_paths
  )
  #### ---------------------------------------------------------------------------
  
}

shinyApp(ui, server)
