# ============================================================
# EWAS → nearby genes (hg38) for enrichment
# - Input: bins table (chr + bin_start + bin_end)
# - Output:
#   (1) detailed table: bin ↔ genes (ENTREZID + SYMBOL)
#   (2) gene vector for enrichment (unique ENTREZIDs or SYMBOLs)
# Requires: TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicRanges, IRanges,
#           AnnotationDbi, org.Hs.eg.db, data.table (optional but recommended)
# ============================================================

# -----------------------------
# Helpers: normalize chromosome labels
# -----------------------------
norm_chr_hg38 <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^CHR", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "M"] <- "MT"
  paste0("chr", x)
}

# -----------------------------
# Build GRanges from bins table
# Expected columns (any case accepted):
#   chr, bin_start, bin_end
# Optional: cluster_id / disease / best_padj / etc (kept if present)
# -----------------------------
bins_to_granges <- function(bins_df) {
  if (is.null(bins_df) || !is.data.frame(bins_df) || !nrow(bins_df)) return(NULL)
  
  pick_col <- function(df, candidates) {
    nm <- intersect(candidates, names(df))
    if (length(nm)) nm[1] else NULL
  }
  
  chr_col <- pick_col(bins_df, c("chr","CHR","chrom","chromosome"))
  st_col  <- pick_col(bins_df, c("bin_start","start","BIN_START","start_bp","start1"))
  en_col  <- pick_col(bins_df, c("bin_end","end","BIN_END","end_bp","end1"))
  if (is.null(chr_col) || is.null(st_col) || is.null(en_col)) return(NULL)
  
  chr <- norm_chr_hg38(bins_df[[chr_col]])
  st  <- suppressWarnings(as.integer(bins_df[[st_col]]))
  en  <- suppressWarnings(as.integer(bins_df[[en_col]]))
  
  ok <- !is.na(chr) & is.finite(st) & is.finite(en)
  if (!any(ok)) return(NULL)
  
  st2 <- pmin(st[ok], en[ok])
  en2 <- pmax(st[ok], en[ok])
  chr2 <- chr[ok]
  
  # Build an ID per bin row to keep join stable
  bin_row_id <- seq_len(sum(ok))
  
  gr <- GenomicRanges::GRanges(
    seqnames = chr2,
    ranges   = IRanges::IRanges(start = st2, end = en2)
  )
  GenomicRanges::mcols(gr)$bin_row_id <- bin_row_id
  
  # Keep original columns (subsetted to ok rows)
  meta <- bins_df[ok, , drop = FALSE]
  meta$bin_row_id <- bin_row_id
  attr(gr, "bins_meta") <- meta
  
  gr
}

# -----------------------------
# Main: annotate bins with nearby genes
# flank_bp: genes within +/- flank_bp of bin interval (0 = overlap only)
# gene_model: "gene" uses TxDb genes() ranges; "transcript" uses transcripts() (usually not needed)
# return_ids: "ENTREZID" (recommended) or "SYMBOL"
# -----------------------------
ewas_bins_nearby_genes <- function(
    bins_df,
    txdb,
    flank_bp   = 50000L,
    return_ids = c("ENTREZID","SYMBOL"),
    keep_all_bins = TRUE
) {
  return_ids <- match.arg(return_ids)
  
  gr_bins <- bins_to_granges(bins_df)
  if (is.null(gr_bins)) {
    return(list(detail = data.frame(), genes = character()))
  }
  
  bins_meta <- attr(gr_bins, "bins_meta")
  if (is.null(bins_meta)) bins_meta <- bins_df
  
  flank_bp <- as.integer(flank_bp)
  if (!is.finite(flank_bp) || flank_bp < 0) flank_bp <- 0L
  
  # Expand bins by flank (if requested)
  gr_q <- gr_bins
  if (flank_bp > 0L) {
    gr_q <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr_bins),
      ranges   = IRanges::resize(GenomicRanges::ranges(gr_bins),
                                 width = GenomicRanges::width(gr_bins) + 2L * flank_bp,
                                 fix = "center")
    )
    GenomicRanges::mcols(gr_q)$bin_row_id <- GenomicRanges::mcols(gr_bins)$bin_row_id
  }
  
  # Gene model from TxDb (hg38 knownGene => ENTREZID as gene_id)
  gr_genes <- GenomicFeatures::genes(txdb)
  # Ensure seqlevels style compatible (chr-prefixed)
  GenomeInfoDb::seqlevelsStyle(gr_genes) <- GenomeInfoDb::seqlevelsStyle(gr_q)[1]
  GenomeInfoDb::seqlevelsStyle(gr_q)     <- GenomeInfoDb::seqlevelsStyle(gr_genes)[1]
  
  # Overlaps: bins (query) ↔ genes (subject)
  hits <- GenomicRanges::findOverlaps(gr_q, gr_genes, ignore.strand = TRUE)
  
  if (length(hits) == 0) {
    # Optionally return empty but keep bins
    detail0 <- if (keep_all_bins && is.data.frame(bins_meta) && nrow(bins_meta)) {
      bins_meta$ENTREZID <- NA_character_
      bins_meta$SYMBOL   <- NA_character_
      bins_meta$dist_to_gene_bp <- NA_integer_
      bins_meta
    } else data.frame()
    return(list(detail = detail0, genes = character()))
  }
  
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)
  
  # Bin row id
  bin_row_id <- GenomicRanges::mcols(gr_q)$bin_row_id[q_idx]
  
  # ENTREZID comes from TxDb knownGene gene_id
  entrez <- as.character(GenomicRanges::mcols(gr_genes)$gene_id[s_idx])
  entrez[!nzchar(entrez)] <- NA_character_
  
  # Compute distance bin ↔ gene (0 if overlapping)
  # Use ORIGINAL (non-flanked) bin coordinates for distance reporting
  gr_bins_map <- gr_bins[match(bin_row_id, GenomicRanges::mcols(gr_bins)$bin_row_id)]
  dist_bp <- IRanges::distance(GenomicRanges::ranges(gr_bins_map),
                               GenomicRanges::ranges(gr_genes[s_idx]))
  dist_bp <- suppressWarnings(as.integer(dist_bp))
  
  # SYMBOL mapping (org.Hs.eg.db)
  sym <- rep(NA_character_, length(entrez))
  ok_entrez <- !is.na(entrez) & nzchar(entrez)
  if (any(ok_entrez)) {
    sym_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys     = unique(entrez[ok_entrez]),
      column   = "SYMBOL",
      keytype  = "ENTREZID",
      multiVals = "first"
    )
    sym[ok_entrez] <- unname(sym_map[entrez[ok_entrez]])
  }
  
  # Build detail table: bins_meta (by bin_row_id) + ENTREZID/SYMBOL + distance
  detail <- bins_meta[match(bin_row_id, bins_meta$bin_row_id), , drop = FALSE]
  detail$ENTREZID <- entrez
  detail$SYMBOL   <- sym
  detail$dist_to_gene_bp <- dist_bp
  
  # De-duplicate (same bin ↔ same gene)
  detail <- detail[!duplicated(detail[, c("bin_row_id","ENTREZID","SYMBOL")]), , drop = FALSE]
  
  # If keep_all_bins, ensure bins with no hit appear with NA gene
  if (keep_all_bins) {
    all_ids <- unique(bins_meta$bin_row_id)
    hit_ids <- unique(detail$bin_row_id)
    miss_ids <- setdiff(all_ids, hit_ids)
    if (length(miss_ids)) {
      miss <- bins_meta[bins_meta$bin_row_id %in% miss_ids, , drop = FALSE]
      miss$ENTREZID <- NA_character_
      miss$SYMBOL   <- NA_character_
      miss$dist_to_gene_bp <- NA_integer_
      detail <- rbind(detail, miss)
    }
  }
  
  # Final genes vector for enrichment
  genes_vec <- if (return_ids == "ENTREZID") detail$ENTREZID else detail$SYMBOL
  genes_vec <- unique(as.character(genes_vec))
  genes_vec <- genes_vec[!is.na(genes_vec) & nzchar(genes_vec)]
  
  list(detail = detail, genes = genes_vec)
}

# -----------------------------
# Convenience wrapper:
# from ewas_bins_all (or filtered bins) -> genes for enrichment
# Example:
#   out <- ewas_genes_for_enrichment(rv$ewas_bins_all, txdb, flank_bp=50000, id="ENTREZID")
#   genes <- out$genes
# -----------------------------
ewas_genes_for_enrichment <- function(
    bins_df,
    txdb,
    flank_bp = 50000L,
    id = c("ENTREZID","SYMBOL")
) {
  id <- match.arg(id)
  ewas_bins_nearby_genes(
    bins_df   = bins_df,
    txdb      = txdb,
    flank_bp  = flank_bp,
    return_ids = id,
    keep_all_bins = FALSE
  )
}

# ============================================================
# GItools wrapper expected by the app:
# gi_ewas_genes_nearby(bins_df, flank_bp) -> data.frame with gene/symbol
# ============================================================
gi_ewas_genes_nearby <- function(bins_df, flank_bp = 25000L, id = c("SYMBOL","ENTREZID"), txdb = NULL) {
  id <- match.arg(id)
  
  # txdb default (hg38 knownGene)
  if (is.null(txdb)) {
    if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
      stop("TxDb.Hsapiens.UCSC.hg38.knownGene not available.")
    }
  }
  
  out <- ewas_bins_nearby_genes(
    bins_df       = bins_df,
    txdb          = txdb,
    flank_bp      = as.integer(flank_bp),
    return_ids    = if (id == "ENTREZID") "ENTREZID" else "SYMBOL",
    keep_all_bins = TRUE
  )
  
  df <- out$detail
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  
  # Column expected by the module
  if (id == "ENTREZID") {
    df$gene <- as.character(df$ENTREZID)
  } else {
    df$symbol <- as.character(df$SYMBOL)
    df$gene   <- df$symbol
  }
  
  df
}
