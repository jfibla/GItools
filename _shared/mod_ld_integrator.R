# R/mod_ld_integrator.R
# ------------------------------------------------------------------
# LD module for GItools Integrator
# - Reads all shared RDS from gi_shared_root/integrator_exports
# - Uses:
#     * <app>_clusters_master.rds
#     * <app>_candidates.rds
# - Apps:
#     catalog, gtex, nonsyn, ewastum, ewasdis
# - Workflows:
#     1) Plot / full LD matrix for selected integrated cluster
#     2) Prioritization / seed + proxies for GWAS hits in selected cluster
# ------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
}

pick_col <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm)) nm[1] else NULL
}

chr_map_plink19 <- function(x){
  x <- toupper(as.character(x))
  x <- trimws(x)
  x <- sub("^CHR", "", x)
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

harmonize_integrator_clusters <- function(df, source_app = NA_character_) {
  stopifnot(is.data.frame(df))
  
  need <- c("cluster_id", "chr", "start", "end", "cluster_key")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("clusters_master missing columns: ", paste(miss, collapse = ", "))
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  df$cluster_id      <- normalize_cluster_id(df$cluster_id)
  df$cluster_id_raw  <- as.character(df$cluster_id)
  df$chr             <- chr_map_plink19(df$chr)
  df$start           <- suppressWarnings(as.integer(readr::parse_number(as.character(df$start))))
  df$end             <- suppressWarnings(as.integer(readr::parse_number(as.character(df$end))))
  df$cluster_key     <- trimws(as.character(df$cluster_key))
  df$source_app      <- as.character(source_app)
  
  s0 <- pmin(df$start, df$end)
  e0 <- pmax(df$start, df$end)
  df$start <- s0
  df$end   <- e0
  
  idx_key <- is.na(df$cluster_key) | !nzchar(df$cluster_key)
  df$cluster_key[idx_key] <- paste0(
    df$cluster_id[idx_key], "|chr", chr_label_plink(df$chr[idx_key]), ":",
    df$start[idx_key], "-", df$end[idx_key]
  )
  
  df <- df[
    !is.na(df$cluster_id) & nzchar(df$cluster_id) &
      !is.na(df$chr) &
      !is.na(df$start) &
      !is.na(df$end),
    c("source_app", "cluster_id", "chr", "start", "end", "cluster_key"),
    drop = FALSE
  ]
  
  df <- unique(df)
  rownames(df) <- NULL
  df
}

harmonize_integrator_candidates <- function(df, source_app = NA_character_) {
  stopifnot(is.data.frame(df))
  
  need <- c("cluster_id", "chr", "pos_ini", "pos_end", "id_hit", "rsid", "position", "classe")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("candidates missing columns: ", paste(miss, collapse = ", "))
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  df$source_app <- as.character(source_app)
  df$cluster_id <- normalize_cluster_id(df$cluster_id)
  df$chr        <- chr_map_plink19(df$chr)
  df$pos_ini    <- suppressWarnings(as.integer(readr::parse_number(as.character(df$pos_ini))))
  df$pos_end    <- suppressWarnings(as.integer(readr::parse_number(as.character(df$pos_end))))
  df$id_hit     <- trimws(as.character(df$id_hit))
  df$rsid       <- trimws(as.character(df$rsid))
  df$position   <- suppressWarnings(as.integer(readr::parse_number(as.character(df$position))))
  df$classe     <- normalize_classe(df$classe)
  
  df$rsid[df$rsid == "" | is.na(df$rsid)] <- df$id_hit[df$rsid == "" | is.na(df$rsid)]
  df$position[is.na(df$position)] <- df$pos_ini[is.na(df$position)]
  df$pos_end[is.na(df$pos_end)] <- df$pos_ini[is.na(df$pos_end)]
  
  df <- df[
    !is.na(df$cluster_id) & nzchar(df$cluster_id) &
      !is.na(df$chr) &
      !is.na(df$pos_ini) &
      !is.na(df$pos_end) &
      !is.na(df$id_hit) & nzchar(df$id_hit),
    c("source_app", "cluster_id", "chr", "pos_ini", "pos_end", "id_hit", "rsid", "position", "classe"),
    drop = FALSE
  ]
  
  df <- unique(df)
  rownames(df) <- NULL
  df
}

prepare_clusters <- function(df) {
  chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
  st_col  <- pick_col(df, c("start","START","start_bp","cluster_start","FROM","from","bp1"))
  en_col  <- pick_col(df, c("end","END","end_bp","cluster_end","TO","to","bp2"))
  id_col  <- pick_col(df, c("cluster_id","CLUSTER_ID","cluster","id"))
  key_col <- pick_col(df, c("cluster_key","CLUSTER_KEY"))
  
  out <- df %>%
    dplyr::transmute(
      chr = chr_map_plink19(.data[[chr_col]]),
      start = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[st_col]])))),
      end   = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[en_col]])))),
      cluster_id = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_,
      cluster_key = if (!is.null(key_col)) as.character(.data[[key_col]]) else NA_character_,
      source_apps = if ("source_apps" %in% names(df)) as.character(source_apps) else NA_character_
    ) %>%
    dplyr::filter(is.finite(chr), is.finite(start), is.finite(end)) %>%
    dplyr::mutate(
      start0 = start,
      end0   = end,
      start  = pmin(start0, end0),
      end    = pmax(start0, end0),
      cluster_id = normalize_cluster_id(cluster_id)
    ) %>%
    dplyr::select(-start0, -end0) %>%
    dplyr::arrange(chr, start, end) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(cluster_n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cluster_id = ifelse(
        is.na(cluster_id) | !nzchar(cluster_id),
        paste0("chr", chr_label_plink(chr), "_", cluster_n),
        cluster_id
      ),
      cluster_key = dplyr::if_else(
        !is.na(cluster_key) & nzchar(cluster_key),
        cluster_key,
        paste0(cluster_id, "|chr", chr_label_plink(chr), ":", start, "-", end)
      ),
      size_kb = round((end - start) / 1000, 2)
    )
  
  out
}

prepare_candidates <- function(df) {
  rs_col  <- pick_col(df, c("rsid","RSID","snp","SNP","id","ID","marker","MARKER","id_hit","ID_HIT"))
  chr_col <- pick_col(df, c("chr","CHR","chrom","CHROM","chromosome"))
  pos_col <- pick_col(df, c("position","POS","pos","BP","bp","pos_ini","pos_start","start","POS_INI"))
  end_col <- pick_col(df, c("pos_end","POS_END","end","END"))
  cl_col  <- pick_col(df, c("classe","class","CLASS","Classe","type","TYPE"))
  id_col  <- pick_col(df, c("cluster_id","CLUSTER_ID","cluster","id"))
  
  if (is.null(rs_col) || is.null(chr_col) || is.null(pos_col)) {
    stop("Candidates table must have id_hit/rsid + chr + position/pos_ini columns. Found: ",
         paste(names(df), collapse = ", "))
  }
  
  out <- df %>%
    dplyr::transmute(
      source_app  = if ("source_app" %in% names(df)) as.character(source_app) else NA_character_,
      cluster_id  = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_,
      id_hit      = if ("id_hit" %in% names(df)) as.character(id_hit) else as.character(.data[[rs_col]]),
      rsid        = as.character(.data[[rs_col]]),
      chr         = chr_map_plink19(.data[[chr_col]]),
      position    = suppressWarnings(as.integer(readr::parse_number(as.character(.data[[pos_col]])))),
      pos_end     = if (!is.null(end_col)) suppressWarnings(as.integer(readr::parse_number(as.character(.data[[end_col]])))) else NA_integer_,
      classe      = if (!is.null(cl_col)) as.character(.data[[cl_col]]) else NA_character_
    ) %>%
    dplyr::mutate(
      cluster_id = normalize_cluster_id(cluster_id),
      id_hit     = trimws(as.character(id_hit)),
      rsid       = trimws(as.character(rsid)),
      classe     = normalize_classe(classe),
      pos_end    = dplyr::coalesce(pos_end, position),
      rsid       = dplyr::if_else(is.na(rsid) | !nzchar(rsid), id_hit, rsid)
    ) %>%
    dplyr::filter(nzchar(id_hit), is.finite(chr), is.finite(position)) %>%
    dplyr::distinct(source_app, cluster_id, id_hit, rsid, chr, position, classe, .keep_all = TRUE)
  
  out
}

available_pops_dir <- function(ld_pops_dir) {
  d <- tryCatch(normalizePath(ld_pops_dir, winslash = "/", mustWork = FALSE),
                error = function(e) ld_pops_dir)
  if (is.null(d) || !nzchar(d) || !dir.exists(d)) return(character(0))
  
  ff <- list.files(d, pattern = "\\.[Tt][Xx][Tt]$", full.names = FALSE)
  pops <- sub("\\.[Tt][Xx][Tt]$", "", ff)
  pops <- pops[nzchar(pops)]
  sort(unique(pops))
}

read_keep_file_dir <- function(pop, ld_pops_dir) {
  d <- tryCatch(normalizePath(ld_pops_dir, winslash = "/", mustWork = FALSE),
                error = function(e) ld_pops_dir)
  if (is.null(d) || !nzchar(d) || !dir.exists(d)) {
    stop("POP keep-files directory not found: ", ld_pops_dir)
  }
  
  pop <- as.character(pop %||% "")
  if (!nzchar(pop)) {
    av <- available_pops_dir(d)
    stop("Population not selected. Available: ", paste(av, collapse = ", "))
  }
  
  cand <- c(
    file.path(d, paste0(pop, ".txt")),
    file.path(d, paste0(pop, ".TXT")),
    file.path(d, paste0(pop, ".TxT"))
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f) || !nzchar(f)) {
    av <- available_pops_dir(d)
    stop("Missing keep file for pop='", pop, "' in ", d,
         ". Available: ", paste(av, collapse = ", "))
  }
  
  x <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 1) stop("Keep file must have at least 1 column (IID)")
  
  if (ncol(x) == 1) x <- cbind(x[,1], x[,1]) else x <- x[,1:2]
  colnames(x) <- c("FID","IID")
  
  out <- file.path(tempdir(), paste0("keep_", pop, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  write.table(x, out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  out
}

gi_ld_defaults <- function() {
  if (exists("gi_cfg", mode = "function")) {
    cfg <- gi_cfg()
    res <- cfg$resources %||% ""
    popdir <- cfg$pop_dir %||% ""
    shared <- cfg$shared %||% ""
  } else {
    res    <- Sys.getenv("GITOOLS_RESOURCES", unset = "")
    popdir <- Sys.getenv("GITOOLS_POP_DIR", unset = "")
    shared <- Sys.getenv("GITOOLS_SHARED", unset = "")
  }
  
  res    <- normalizePath(res, winslash = "/", mustWork = FALSE)
  popdir <- normalizePath(popdir, winslash = "/", mustWork = FALSE)
  shared <- normalizePath(shared, winslash = "/", mustWork = FALSE)
  
  plink <- Sys.getenv("GITOOLS_PLINK19", unset = "")
  if (!nzchar(plink) && nzchar(res)) {
    cand1 <- file.path(res, "software", "plink19", "plink")
    cand2 <- file.path(res, "software", "plink19", "plink19")
    cand3 <- file.path(res, "software", "plink19")
    plink <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else cand3
  }
  
  bfile <- Sys.getenv("GITOOLS_LD_BFILE", unset = "")
  if (!nzchar(bfile) && nzchar(res)) {
    bfile <- file.path(res, "LD_resources",
                       "Merged_FULL_SET_hg38_hgdp.wgs_10000G3_ko07_MAF0.05")
  }
  
  list(
    plink = plink,
    bfile = bfile,
    popdir = popdir,
    resources = res,
    shared = shared
  )
}

run_plink <- function(args, out_prefix, plink_bin) {
  res <- tryCatch(
    system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) paste("system2 error:", conditionMessage(e))
  )
  st <- attr(res, "status")
  if (is.null(st)) st <- 0
  
  log_file <- paste0(out_prefix, ".log")
  log_tail <- character(0)
  if (file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE)
    log_tail <- utils::tail(log_lines, 40)
  }
  list(stdout = res, status = st, log_file = log_file, log_tail = log_tail)
}

read_plink_ld <- function(prefix) {
  f_gz <- paste0(prefix, ".ld.gz")
  f_tx <- paste0(prefix, ".ld")
  f <- if (file.exists(f_gz)) f_gz else if (file.exists(f_tx)) f_tx else NA_character_
  if (is.na(f)) stop("Missing LD file: expected ", f_gz, " or ", f_tx)
  
  con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con), add = TRUE)
  
  x <- utils::read.table(
    con, header = TRUE, sep = "",
    stringsAsFactors = FALSE, fill = TRUE, comment.char = ""
  )
  names(x) <- trimws(names(x))
  x
}

read_plink_blocks <- function(prefix) {
  f <- paste0(prefix, ".blocks.det")
  if (!file.exists(f)) stop("Missing blocks file: ", f)
  
  x <- tryCatch(
    read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = ""),
    error = function(e) NULL
  )
  
  if (is.null(x) || nrow(x) == 0 || ncol(x) <= 1) {
    x <- tryCatch(
      read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "", quote = "", comment.char = ""),
      error = function(e) NULL
    )
  }
  
  if (is.null(x) || nrow(x) == 0) {
    stop("Could not parse .blocks.det (bad delimiter or empty). File=", f)
  }
  x
}

blocks_det_to_ij <- function(blocks_raw, fl) {
  stopifnot(is.data.frame(blocks_raw), is.data.frame(fl))
  stopifnot(all(c("SNP","BP") %in% names(fl)))
  
  fl <- fl %>%
    dplyr::arrange(BP, SNP) %>%
    dplyr::mutate(ix = dplyr::row_number())
  
  snp_to_ix <- stats::setNames(fl$ix, fl$SNP)
  
  nm <- names(blocks_raw)
  snps_col <- nm[match(TRUE, nm %in% c("SNPS","Snps","snps","SNPs"))]
  bp1_col  <- nm[match(TRUE, nm %in% c("BP1","bp1","FROM_BP","from_bp","START_BP","start_bp"))]
  bp2_col  <- nm[match(TRUE, nm %in% c("BP2","bp2","TO_BP","to_bp","END_BP","end_bp"))]
  
  out <- NULL
  
  if (!is.na(snps_col)) {
    snps_str <- as.character(blocks_raw[[snps_col]])
    snps_str[is.na(snps_str)] <- ""
    
    split_snps <- function(s) {
      s <- trimws(s)
      if (!nzchar(s)) return(character(0))
      if (grepl("\\|", s, fixed = FALSE)) return(strsplit(s, "\\|", fixed = FALSE)[[1]])
      if (grepl(",", s, fixed = TRUE)) return(strsplit(s, ",", fixed = TRUE)[[1]])
      strsplit(s, "\\s+")[[1]]
    }
    
    lst <- lapply(snps_str, split_snps)
    first_snp <- vapply(lst, function(v) if (length(v)) v[1] else NA_character_, character(1))
    last_snp  <- vapply(lst, function(v) if (length(v)) v[length(v)] else NA_character_, character(1))
    
    i <- suppressWarnings(as.integer(snp_to_ix[first_snp]))
    j <- suppressWarnings(as.integer(snp_to_ix[last_snp]))
    
    out <- data.frame(i = i, j = j, stringsAsFactors = FALSE)
  }
  
  if (is.null(out) || all(is.na(out$i)) || all(is.na(out$j))) {
    if (is.na(bp1_col) || is.na(bp2_col)) {
      return(tibble::tibble(i = integer(0), j = integer(0)))
    }
    
    bp1 <- suppressWarnings(as.numeric(blocks_raw[[bp1_col]]))
    bp2 <- suppressWarnings(as.numeric(blocks_raw[[bp2_col]]))
    
    bp_to_ix <- function(bp) {
      if (!is.finite(bp)) return(NA_integer_)
      k <- which.min(abs(fl$BP - bp))
      fl$ix[k]
    }
    i <- vapply(bp1, bp_to_ix, integer(1))
    j <- vapply(bp2, bp_to_ix, integer(1))
    
    out <- data.frame(i = i, j = j, stringsAsFactors = FALSE)
  }
  
  out <- out %>%
    dplyr::mutate(
      i = suppressWarnings(as.integer(i)),
      j = suppressWarnings(as.integer(j))
    ) %>%
    dplyr::filter(is.finite(i), is.finite(j), i >= 1, j >= 1, j > i, j <= nrow(fl)) %>%
    dplyr::distinct(i, j, .keep_all = TRUE) %>%
    dplyr::arrange(i, j)
  
  tibble::as_tibble(out)
}

find_seed_snp <- function(hit_rsid, hit_pos, bim_df) {
  stopifnot(is.data.frame(bim_df), all(c("SNP","BP") %in% names(bim_df)))
  
  hit_rsid <- trimws(as.character(hit_rsid %||% ""))
  hit_pos  <- suppressWarnings(as.integer(hit_pos))
  
  if (nzchar(hit_rsid) && hit_rsid %in% bim_df$SNP) {
    row <- bim_df[bim_df$SNP == hit_rsid, , drop = FALSE][1, , drop = FALSE]
    return(data.frame(
      seed_snp = as.character(row$SNP),
      seed_pos = as.integer(row$BP),
      seed_type = "exact",
      seed_dist_bp = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  if (!is.finite(hit_pos) || nrow(bim_df) == 0) {
    return(data.frame(
      seed_snp = NA_character_,
      seed_pos = NA_integer_,
      seed_type = "missing",
      seed_dist_bp = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  
  k <- which.min(abs(bim_df$BP - hit_pos))
  data.frame(
    seed_snp = as.character(bim_df$SNP[k]),
    seed_pos = as.integer(bim_df$BP[k]),
    seed_type = "nearest",
    seed_dist_bp = as.integer(abs(bim_df$BP[k] - hit_pos)),
    stringsAsFactors = FALSE
  )
}

compute_seed_ld <- function(subset_prefix, seed_snp, plink_bin, out_prefix, metric = "R2") {
  args <- c(
    "--bfile", subset_prefix,
    "--ld-snp", seed_snp,
    "--ld-window", "999999",
    "--ld-window-kb", "999999",
    "--ld-window-r2", "0"
  )
  
  if (identical(metric, "Dprime")) {
    args <- c(args, "--r2", "dprime", "gz")
  } else {
    args <- c(args, "--r2", "gz")
  }
  
  args <- c(args, "--out", out_prefix)
  run_plink(args, out_prefix, plink_bin)
}

build_hit_proxy_table <- function(hit_row, bim_df, subset_prefix, plink_bin, workdir, metric = "R2") {
  hit_rsid <- as.character(hit_row$rsid[1])
  hit_pos  <- as.integer(hit_row$position[1])
  cl_id    <- as.character(hit_row$cluster_id[1] %||% NA_character_)
  
  seed <- find_seed_snp(hit_rsid = hit_rsid, hit_pos = hit_pos, bim_df = bim_df)
  
  if (!is.finite(seed$seed_pos[1]) || is.na(seed$seed_snp[1]) || !nzchar(seed$seed_snp[1])) {
    return(NULL)
  }
  
  out_prefix <- file.path(
    workdir,
    paste0("seed_", gsub("[^A-Za-z0-9_\\-]", "_", seed$seed_snp[1]), "_", format(Sys.time(), "%H%M%S%OS3"))
  )
  
  rr <- compute_seed_ld(
    subset_prefix = subset_prefix,
    seed_snp = seed$seed_snp[1],
    plink_bin = plink_bin,
    out_prefix = out_prefix,
    metric = metric
  )
  if (is.null(rr$status) || rr$status != 0) return(NULL)
  
  ld <- read_plink_ld(out_prefix)
  val_col <- if (identical(metric, "Dprime")) pick_col(ld, c("Dprime","D'","DP","DPRIME")) else pick_col(ld, c("R2","r2"))
  snpA <- pick_col(ld, c("SNP_A","SNP_A1","SNP1"))
  snpB <- pick_col(ld, c("SNP_B","SNP_B1","SNP2"))
  
  if (is.null(val_col) || is.null(snpA) || is.null(snpB)) return(NULL)
  
  out <- ld %>%
    dplyr::transmute(
      SNP_A = as.character(.data[[snpA]]),
      SNP_B = as.character(.data[[snpB]]),
      ld_value = suppressWarnings(as.numeric(.data[[val_col]]))
    ) %>%
    dplyr::filter(is.finite(ld_value)) %>%
    dplyr::mutate(
      proxy_snp = dplyr::if_else(SNP_A == seed$seed_snp[1], SNP_B, SNP_A)
    ) %>%
    dplyr::left_join(
      bim_df %>% dplyr::transmute(proxy_snp = SNP, proxy_pos = BP),
      by = "proxy_snp"
    ) %>%
    dplyr::mutate(
      cluster_id    = cl_id,
      query_hit     = hit_rsid,
      query_pos     = hit_pos,
      seed_snp      = seed$seed_snp[1],
      seed_pos      = seed$seed_pos[1],
      seed_type     = seed$seed_type[1],
      seed_dist_bp  = seed$seed_dist_bp[1]
    ) %>%
    dplyr::select(
      cluster_id, query_hit, query_pos,
      seed_snp, seed_pos, seed_type, seed_dist_bp,
      proxy_snp, proxy_pos, ld_value
    ) %>%
    dplyr::distinct()
  
  out
}

# detail-level helper for per-lead / per-block support tables
# not used as canonical source for block priority summary
add_block_marker_support <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(df)
  }
  
  if (!"lead_pos" %in% names(df)) df$lead_pos <- NA_integer_
  
  ext_sum <- df %>%
    dplyr::mutate(
      source_app = as.character(source_app),
      classe = as.character(classe),
      is_external_hit = !is.na(classe) & classe %in% c(
        "catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit"
      )
    ) %>%
    dplyr::group_by(cluster_id, lead_snp, lead_pos, block_id) %>%
    dplyr::summarise(
      has_external_support = any(is_external_hit, na.rm = TRUE),
      n_external_hits = sum(is_external_hit, na.rm = TRUE),
      n_external_apps = dplyr::n_distinct(source_app[is_external_hit & !is.na(source_app) & nzchar(source_app)]),
      external_apps = paste(sort(unique(stats::na.omit(source_app[is_external_hit & nzchar(source_app)]))), collapse = "; "),
      external_classes = paste(sort(unique(stats::na.omit(classe[is_external_hit & nzchar(classe)]))), collapse = "; "),
      marker_status = dplyr::if_else(any(is_external_hit, na.rm = TRUE), "MARKER", "NO_EXTERNAL_SUPPORT"),
      .groups = "drop"
    )
  
  df %>%
    dplyr::left_join(
      ext_sum,
      by = c("cluster_id", "lead_snp", "lead_pos", "block_id")
    ) %>%
    dplyr::mutate(
      has_external_support = dplyr::coalesce(has_external_support, FALSE),
      n_external_hits = dplyr::coalesce(as.integer(n_external_hits), 0L),
      n_external_apps = dplyr::coalesce(as.integer(n_external_apps), 0L),
      external_apps = dplyr::coalesce(external_apps, ""),
      external_classes = dplyr::coalesce(external_classes, ""),
      marker_status = dplyr::coalesce(marker_status, "NO_EXTERNAL_SUPPORT")
    )
}

summarize_hits_by_block <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble(
      cluster_id = character(),
      lead_snp = character(),
      lead_pos = integer(),
      block_id = character(),
      block_label = character(),
      block_start = integer(),
      block_end = integer(),
      block_size_bp = integer(),
      chr = integer(),
      gwas_hits = character(),
      catalog_hits = character(),
      gtex_hits = character(),
      nonsyn_hits = character(),
      ewasdis_hits = character(),
      ewastum_hits = character(),
      catalog_pos = character(),
      gtex_pos = character(),
      nonsyn_pos = character(),
      ewasdis_pos = character(),
      ewastum_pos = character(),
      n_gwas_hits = integer(),
      n_catalog_hits = integer(),
      n_gtex_hits = integer(),
      n_nonsyn_hits = integer(),
      n_ewasdis_hits = integer(),
      n_ewastum_hits = integer(),
      has_external_support = logical(),
      n_external_apps = integer(),
      external_apps = character(),
      marker_status = character()
    ))
  }
  
  df2 <- df
  
  if (!("block_label" %in% names(df2))) {
    df2$block_label <- as.character(df2$block_id)
  }
  
  if (!("chr" %in% names(df2))) {
    df2$chr <- NA_integer_
  }
  
  if (!("source_app" %in% names(df2))) {
    df2$source_app <- NA_character_
  }
  
  if (!("classe" %in% names(df2))) {
    df2$classe <- NA_character_
  }
  
  if (!("hit_rsid" %in% names(df2))) {
    df2$hit_rsid <- NA_character_
  }
  
  if (!("hit_id" %in% names(df2))) {
    df2$hit_id <- NA_character_
  }
  
  if (!("hit_pos" %in% names(df2))) {
    df2$hit_pos <- NA_integer_
  }
  
  df2 <- df2 %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      lead_snp = as.character(lead_snp),
      block_id = as.character(block_id),
      block_label = as.character(block_label),
      chr = suppressWarnings(as.integer(chr)),
      source_app = as.character(source_app),
      classe = as.character(classe),
      hit_rsid = as.character(hit_rsid),
      hit_id = as.character(hit_id),
      hit_pos = suppressWarnings(as.integer(hit_pos)),
      hit_show = dplyr::coalesce(hit_rsid, hit_id)
    )
  
  collapse_hits <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- x[nzchar(x)]
    paste(sort(x), collapse = "; ")
  }
  
  collapse_pos <- function(x) {
    x <- unique(stats::na.omit(suppressWarnings(as.integer(x))))
    x <- x[is.finite(x)]
    paste(sort(x), collapse = "; ")
  }
  
  df2 %>%
    dplyr::group_by(
      cluster_id, lead_snp, lead_pos,
      block_id, block_label, block_start, block_end, block_size_bp, chr
    ) %>%
    dplyr::summarise(
      gwas_hits    = collapse_hits(hit_show[classe == "GWAS"]),
      catalog_hits = collapse_hits(hit_show[classe == "catalog_hit"]),
      gtex_hits    = collapse_hits(hit_show[classe == "gtex_hit"]),
      nonsyn_hits  = collapse_hits(hit_show[classe == "nonsyn_hit"]),
      ewasdis_hits = collapse_hits(hit_show[classe == "ewasdis_hit"]),
      ewastum_hits = collapse_hits(hit_show[classe == "ewastum_hit"]),
      
      catalog_pos = collapse_pos(hit_pos[classe == "catalog_hit"]),
      gtex_pos    = collapse_pos(hit_pos[classe == "gtex_hit"]),
      nonsyn_pos  = collapse_pos(hit_pos[classe == "nonsyn_hit"]),
      ewasdis_pos = collapse_pos(hit_pos[classe == "ewasdis_hit"]),
      ewastum_pos = collapse_pos(hit_pos[classe == "ewastum_hit"]),
      
      n_gwas_hits    = dplyr::n_distinct(hit_show[classe == "GWAS" & !is.na(hit_show) & nzchar(hit_show)]),
      n_catalog_hits = dplyr::n_distinct(hit_show[classe == "catalog_hit" & !is.na(hit_show) & nzchar(hit_show)]),
      n_gtex_hits    = dplyr::n_distinct(hit_show[classe == "gtex_hit" & !is.na(hit_show) & nzchar(hit_show)]),
      n_nonsyn_hits  = dplyr::n_distinct(hit_show[classe == "nonsyn_hit" & !is.na(hit_show) & nzchar(hit_show)]),
      n_ewasdis_hits = dplyr::n_distinct(hit_show[classe == "ewasdis_hit" & !is.na(hit_show) & nzchar(hit_show)]),
      n_ewastum_hits = dplyr::n_distinct(hit_show[classe == "ewastum_hit" & !is.na(hit_show) & nzchar(hit_show)]),
      
      has_external_support = any(
        classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit"),
        na.rm = TRUE
      ),
      
      n_external_apps = dplyr::n_distinct(
        source_app[
          classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit") &
            !is.na(source_app) & nzchar(source_app)
        ]
      ),
      
      external_apps = paste(
        sort(unique(stats::na.omit(
          source_app[
            classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit") &
              !is.na(source_app) & nzchar(source_app)
          ]
        ))),
        collapse = "; "
      ),
      
      marker_status = dplyr::if_else(
        any(
          classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit"),
          na.rm = TRUE
        ),
        "MARKER",
        "NO_EXTERNAL_SUPPORT"
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gwas_hits = dplyr::coalesce(gwas_hits, ""),
      catalog_hits = dplyr::coalesce(catalog_hits, ""),
      gtex_hits = dplyr::coalesce(gtex_hits, ""),
      nonsyn_hits = dplyr::coalesce(nonsyn_hits, ""),
      ewasdis_hits = dplyr::coalesce(ewasdis_hits, ""),
      ewastum_hits = dplyr::coalesce(ewastum_hits, ""),
      catalog_pos = dplyr::coalesce(catalog_pos, ""),
      gtex_pos = dplyr::coalesce(gtex_pos, ""),
      nonsyn_pos = dplyr::coalesce(nonsyn_pos, ""),
      ewasdis_pos = dplyr::coalesce(ewasdis_pos, ""),
      ewastum_pos = dplyr::coalesce(ewastum_pos, ""),
      external_apps = dplyr::coalesce(external_apps, ""),
      n_gwas_hits = dplyr::coalesce(as.integer(n_gwas_hits), 0L),
      n_catalog_hits = dplyr::coalesce(as.integer(n_catalog_hits), 0L),
      n_gtex_hits = dplyr::coalesce(as.integer(n_gtex_hits), 0L),
      n_nonsyn_hits = dplyr::coalesce(as.integer(n_nonsyn_hits), 0L),
      n_ewasdis_hits = dplyr::coalesce(as.integer(n_ewasdis_hits), 0L),
      n_ewastum_hits = dplyr::coalesce(as.integer(n_ewastum_hits), 0L),
      n_external_apps = dplyr::coalesce(as.integer(n_external_apps), 0L),
      has_external_support = dplyr::coalesce(has_external_support, FALSE)
    ) %>%
    dplyr::arrange(
      dplyr::desc(has_external_support),
      dplyr::desc(n_external_apps),
      lead_pos,
      block_start
    )
}



html_escape_basic <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

split_semicolon <- function(x) {
  x <- as.character(x %||% "")
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))
  out <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  out <- trimws(out)
  out[nzchar(out)]
}

make_ucsc_link <- function(chr, pos, flank = 1000L, label = NULL) {
  chr <- suppressWarnings(as.integer(chr))
  pos <- suppressWarnings(as.integer(pos))
  if (!is.finite(chr) || !is.finite(pos)) return(html_escape_basic(label %||% ""))
  
  start <- max(1L, pos - flank)
  end   <- pos + flank
  chr_lab <- paste0("chr", chr_label_plink(chr))
  url <- paste0(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=",
    utils::URLencode(paste0(chr_lab, ":", start, "-", end), reserved = TRUE)
  )
  lab <- html_escape_basic(label %||% paste0(chr_lab, ":", pos))
  paste0("<a href='", url, "' target='_blank'>", lab, "</a>")
}

make_catalog_link <- function(x) {
  x2 <- trimws(as.character(x))
  if (!nzchar(x2)) return("")
  url <- paste0(
    "https://www.ebi.ac.uk/gwas/search?query=",
    utils::URLencode(x2, reserved = TRUE)
  )
  paste0("<a href='", url, "' target='_blank'>", html_escape_basic(x2), "</a>")
}

make_gtex_link <- function(x) {
  x2 <- trimws(as.character(x))
  if (!nzchar(x2)) return("")
  url <- paste0(
    "https://gtexportal.org/home/snp/",
    utils::URLencode(x2, reserved = TRUE)
  )
  paste0("<a href='", url, "' target='_blank'>", html_escape_basic(x2), "</a>")
}

make_dbsnp_link <- function(x) {
  x2 <- trimws(as.character(x))
  if (!nzchar(x2)) return("")
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/snp/",
    utils::URLencode(x2, reserved = TRUE)
  )
  paste0("<a href='", url, "' target='_blank'>", html_escape_basic(x2), "</a>")
}

collapse_link_list_html <- function(links, max_visible = 5L, summary_label = NULL) {
  links <- as.character(links)
  links <- links[nzchar(links)]
  
  if (!length(links)) return("")
  
  if (length(links) <= max_visible) {
    return(paste(links, collapse = "; "))
  }
  
  first_part <- paste(links[seq_len(max_visible)], collapse = "; ")
  rest_part  <- paste(links[(max_visible + 1):length(links)], collapse = "; ")
  
  summ <- summary_label %||% paste0(length(links), " elements")
  
  paste0(
    "<details>",
    "<summary>", html_escape_basic(summ), "</summary>",
    "<div style='margin-top:6px;'>", first_part,
    if (nzchar(rest_part)) paste0("; ", rest_part) else "",
    "</div>",
    "</details>"
  )
}

collapse_plain_list_html <- function(x, max_visible = 5L, summary_label = NULL) {
  vals <- split_semicolon(x)
  
  if (!length(vals)) return("")
  
  if (length(vals) <= max_visible) {
    return(html_escape_basic(paste(vals, collapse = "; ")))
  }
  
  first_part <- html_escape_basic(paste(vals[seq_len(max_visible)], collapse = "; "))
  rest_part  <- html_escape_basic(paste(vals[(max_visible + 1):length(vals)], collapse = "; "))
  
  summ <- summary_label %||% paste0(length(vals), " elements")
  
  paste0(
    "<details>",
    "<summary>", html_escape_basic(summ), "</summary>",
    "<div style='margin-top:6px;'>", first_part,
    if (nzchar(rest_part)) paste0("; ", rest_part) else "",
    "</div>",
    "</details>"
  )
}

format_block_overlap_grouped_for_dt <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  
  out <- df %>%
    dplyr::mutate(
      lead_snps = vapply(
        lead_snps,
        function(z) collapse_plain_list_html(z, max_visible = 5L, summary_label = paste0(length(split_semicolon(z)), " lead SNPs")),
        character(1)
      ),
      lead_positions = vapply(
        lead_positions,
        function(z) collapse_plain_list_html(z, max_visible = 5L, summary_label = paste0(length(split_semicolon(z)), " positions")),
        character(1)
      ),
      gwas_hits = vapply(
        gwas_hits,
        function(z) collapse_plain_list_html(z, max_visible = 5L, summary_label = paste0(length(split_semicolon(z)), " GWAS hits")),
        character(1)
      ),
      external_apps = vapply(
        external_apps,
        function(z) collapse_plain_list_html(z, max_visible = 5L, summary_label = paste0(length(split_semicolon(z)), " apps")),
        character(1)
      )
    )
  
  out <- format_block_overlap_summary_for_dt(out)
  out
}


build_link_column_html <- function(hits_str, type, chr = NA_integer_, pos_str = "") {
  hits <- split_semicolon(hits_str)
  posv <- suppressWarnings(as.integer(split_semicolon(pos_str)))
  
  if (!length(hits)) return("")
  
  links <- switch(
    type,
    catalog = vapply(hits, make_catalog_link, character(1)),
    gtex    = vapply(hits, make_gtex_link, character(1)),
    nonsyn  = vapply(hits, make_dbsnp_link, character(1)),
    ewasdis = vapply(seq_along(hits), function(i) {
      pp <- if (length(posv) >= i) posv[i] else NA_integer_
      make_ucsc_link(chr = chr, pos = pp, flank = 1000L, label = hits[i])
    }, character(1)),
    ewastum = vapply(seq_along(hits), function(i) {
      pp <- if (length(posv) >= i) posv[i] else NA_integer_
      make_ucsc_link(chr = chr, pos = pp, flank = 1000L, label = hits[i])
    }, character(1)),
    hits
  )
  
  collapse_link_list_html(
    links = links,
    max_visible = 5L,
    summary_label = paste0(length(links), " elements")
  )
}

format_block_overlap_summary_for_dt <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  
  x <- df
  
  if (!"chr" %in% names(x)) x$chr <- NA_integer_
  
  if (!"catalog_hits" %in% names(x)) x$catalog_hits <- ""
  if (!"gtex_hits" %in% names(x)) x$gtex_hits <- ""
  if (!"nonsyn_hits" %in% names(x)) x$nonsyn_hits <- ""
  if (!"ewasdis_hits" %in% names(x)) x$ewasdis_hits <- ""
  if (!"ewastum_hits" %in% names(x)) x$ewastum_hits <- ""
  
  if (!"catalog_pos" %in% names(x)) x$catalog_pos <- ""
  if (!"gtex_pos" %in% names(x)) x$gtex_pos <- ""
  if (!"nonsyn_pos" %in% names(x)) x$nonsyn_pos <- ""
  if (!"ewasdis_pos" %in% names(x)) x$ewasdis_pos <- ""
  if (!"ewastum_pos" %in% names(x)) x$ewastum_pos <- ""
  
  x %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      catalog_hits = build_link_column_html(
        hits_str = catalog_hits,
        type = "catalog",
        chr = chr,
        pos_str = catalog_pos
      ),
      gtex_hits = build_link_column_html(
        hits_str = gtex_hits,
        type = "gtex",
        chr = chr,
        pos_str = gtex_pos
      ),
      nonsyn_hits = build_link_column_html(
        hits_str = nonsyn_hits,
        type = "nonsyn",
        chr = chr,
        pos_str = nonsyn_pos
      ),
      ewasdis_hits = build_link_column_html(
        hits_str = ewasdis_hits,
        type = "ewasdis",
        chr = chr,
        pos_str = ewasdis_pos
      ),
      ewastum_hits = build_link_column_html(
        hits_str = ewastum_hits,
        type = "ewastum",
        chr = chr,
        pos_str = ewastum_pos
      )
    ) %>%
    dplyr::ungroup()
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

build_block_summary_from_hits <- function(block_hits_df) {
  if (!is.data.frame(block_hits_df) || !nrow(block_hits_df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      block_label = character(),
      block_start = numeric(),
      block_end = numeric(),
      block_size_bp = numeric(),
      n_block_hits = integer(),
      n_catalog_block_hits = integer(),
      n_gtex_block_hits = integer(),
      n_nonsyn_block_hits = integer(),
      n_ewasdis_block_hits = integer(),
      n_ewastum_block_hits = integer(),
      n_block_apps = integer(),
      n_block_genes = integer(),
      genes_in_block = character()
    ))
  }
  
  df <- block_hits_df
  
  if (!"cluster_id" %in% names(df)) df$cluster_id <- NA_character_
  if (!"block_id" %in% names(df)) df$block_id <- NA_character_
  if (!"block_label" %in% names(df)) df$block_label <- df$block_id
  if (!"block_start" %in% names(df)) df$block_start <- NA_real_
  if (!"block_end" %in% names(df)) df$block_end <- NA_real_
  if (!"block_size_bp" %in% names(df)) df$block_size_bp <- NA_real_
  if (!"classe" %in% names(df)) df$classe <- NA_character_
  if (!"hit_rsid" %in% names(df)) df$hit_rsid <- NA_character_
  if (!"hit_id" %in% names(df)) df$hit_id <- NA_character_
  if (!"hit_pos" %in% names(df)) df$hit_pos <- NA_real_
  
  df %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id   = as.character(block_id),
      block_label = as.character(block_label),
      block_start = suppressWarnings(as.numeric(block_start)),
      block_end   = suppressWarnings(as.numeric(block_end)),
      block_size_bp = suppressWarnings(as.numeric(block_size_bp)),
      classe = as.character(classe),
      hit_rsid = as.character(hit_rsid),
      hit_id   = as.character(hit_id),
      hit_pos  = suppressWarnings(as.integer(hit_pos)),
      hit_key  = dplyr::coalesce(
        dplyr::na_if(hit_rsid, ""),
        dplyr::na_if(hit_id, ""),
        ifelse(is.finite(hit_pos), paste0("pos_", hit_pos), NA_character_)
      )
    ) %>%
    dplyr::filter(!is.na(block_id), nzchar(block_id)) %>%
    dplyr::distinct(
      cluster_id, block_id, classe, hit_key,
      .keep_all = TRUE
    ) %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      block_label   = dplyr::first(block_label),
      block_start   = suppressWarnings(min(block_start, na.rm = TRUE)),
      block_end     = suppressWarnings(max(block_end, na.rm = TRUE)),
      block_size_bp = suppressWarnings(max(block_size_bp, na.rm = TRUE)),
      
      n_block_hits = dplyr::n(),
      n_catalog_block_hits = sum(classe == "catalog_hit", na.rm = TRUE),
      n_gtex_block_hits    = sum(classe == "gtex_hit",    na.rm = TRUE),
      n_nonsyn_block_hits  = sum(classe == "nonsyn_hit",  na.rm = TRUE),
      n_ewasdis_block_hits = sum(classe == "ewasdis_hit", na.rm = TRUE),
      n_ewastum_block_hits = sum(classe == "ewastum_hit", na.rm = TRUE),
      
      n_block_apps =
        as.integer(any(classe == "catalog_hit", na.rm = TRUE)) +
        as.integer(any(classe == "gtex_hit",    na.rm = TRUE)) +
        as.integer(any(classe == "nonsyn_hit",  na.rm = TRUE)) +
        as.integer(any(classe == "ewasdis_hit", na.rm = TRUE)) +
        as.integer(any(classe == "ewastum_hit", na.rm = TRUE)),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      block_start   = ifelse(is.infinite(block_start), NA_real_, block_start),
      block_end     = ifelse(is.infinite(block_end), NA_real_, block_end),
      block_size_bp = ifelse(is.infinite(block_size_bp), NA_real_, block_size_bp),
      n_block_hits = safe_int0(n_block_hits),
      n_catalog_block_hits = safe_int0(n_catalog_block_hits),
      n_gtex_block_hits = safe_int0(n_gtex_block_hits),
      n_nonsyn_block_hits = safe_int0(n_nonsyn_block_hits),
      n_ewasdis_block_hits = safe_int0(n_ewasdis_block_hits),
      n_ewastum_block_hits = safe_int0(n_ewastum_block_hits),
      n_block_apps = safe_int0(n_block_apps)
    )
}

build_block_gene_overlap <- function(block_df) {
  if (!is.data.frame(block_df) || !nrow(block_df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_genes_in_block = integer(),
      genes_in_block = character()
    ))
  }
  
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) ||
      !requireNamespace("GenomicFeatures", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("GenomeInfoDb", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    return(
      block_df %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          n_genes_in_block = 0L,
          genes_in_block = ""
        )
    )
  }
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  genes_gr <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
  
  if (!"chr" %in% names(block_df)) {
    block_df$chr <- NA_integer_
  }
  
  out_list <- lapply(seq_len(nrow(block_df)), function(i) {
    chr_sel <- suppressWarnings(as.integer(block_df$chr[i]))
    st <- suppressWarnings(as.integer(block_df$block_start[i]))
    en <- suppressWarnings(as.integer(block_df$block_end[i]))
    
    if (!is.finite(chr_sel) || !is.finite(st) || !is.finite(en) || st > en) {
      return(data.frame(
        cluster_id = as.character(block_df$cluster_id[i]),
        block_id = as.character(block_df$block_id[i]),
        n_genes_in_block = 0L,
        genes_in_block = "",
        stringsAsFactors = FALSE
      ))
    }
    
    chr_ucsc <- paste0("chr", chr_label_plink(chr_sel))
    block_gr <- GenomicRanges::GRanges(
      seqnames = chr_ucsc,
      ranges = IRanges::IRanges(start = st, end = en)
    )
    
    ov <- GenomicRanges::findOverlaps(genes_gr, block_gr, ignore.strand = TRUE)
    keep_idx <- unique(S4Vectors::queryHits(ov))
    
    if (!length(keep_idx)) {
      return(data.frame(
        cluster_id = as.character(block_df$cluster_id[i]),
        block_id = as.character(block_df$block_id[i]),
        n_genes_in_block = 0L,
        genes_in_block = "",
        stringsAsFactors = FALSE
      ))
    }
    
    gsub_gr <- genes_gr[keep_idx]
    gene_id <- as.character(S4Vectors::mcols(gsub_gr)$gene_id)
    
    mp <- tryCatch(
      suppressMessages(
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = unique(gene_id),
          keytype = "ENTREZID",
          columns = c("SYMBOL")
        )
      ),
      error = function(e) NULL
    )
    
    if (is.data.frame(mp) && nrow(mp)) {
      mp <- mp %>%
        dplyr::mutate(gene_id = as.character(ENTREZID)) %>%
        dplyr::select(gene_id, SYMBOL) %>%
        dplyr::distinct(gene_id, .keep_all = TRUE)
      
      labs <- data.frame(gene_id = gene_id, stringsAsFactors = FALSE) %>%
        dplyr::left_join(mp, by = "gene_id") %>%
        dplyr::mutate(
          label = dplyr::if_else(
            !is.na(SYMBOL) & nzchar(SYMBOL),
            SYMBOL,
            NA_character_
          )
        ) %>%
        dplyr::pull(label)
    } else {
      labs <- gene_id
    }
    
    labs <- sort(unique(labs[!is.na(labs) & nzchar(labs)]))
    
    data.frame(
      cluster_id = as.character(block_df$cluster_id[i]),
      block_id = as.character(block_df$block_id[i]),
      n_genes_in_block = length(labs),
      genes_in_block = paste(labs, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out_list)
}

read_manifest_safe_mod <- function(session_dir) {
  mf <- file.path(session_dir, "manifest.rds")
  if (!file.exists(mf)) return(NULL)
  x <- tryCatch(readRDS(mf), error = function(e) NULL)
  if (!is.list(x)) return(NULL)
  x
}
############# Helpers commons
make_block_canonical_table_display <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  x <- df
  
  if (!"cluster_id" %in% names(x)) x$cluster_id <- NA_character_
  if (!"block_id" %in% names(x)) x$block_id <- NA_character_
  if (!"block_label" %in% names(x)) x$block_label <- NA_character_
  
  if (!"block_start" %in% names(x) && "start" %in% names(x)) x$block_start <- x$start
  if (!"block_end" %in% names(x) && "end" %in% names(x)) x$block_end <- x$end
  if (!"block_start" %in% names(x)) x$block_start <- NA_real_
  if (!"block_end" %in% names(x)) x$block_end <- NA_real_
  
  if (!"block_size_bp" %in% names(x) && "size_bp" %in% names(x)) x$block_size_bp <- x$size_bp
  if (!"block_size_bp" %in% names(x)) x$block_size_bp <- NA_real_
  
  if (!"block_size_kb" %in% names(x) && "size_kb" %in% names(x)) x$block_size_kb <- x$size_kb
  if (!"block_size_kb" %in% names(x)) x$block_size_kb <- NA_real_
  
  if (!"block_support_apps" %in% names(x) && "support_apps" %in% names(x)) x$block_support_apps <- x$support_apps
  if (!"block_support_apps" %in% names(x) && "n_block_apps" %in% names(x)) x$block_support_apps <- x$n_block_apps
  if (!"block_support_apps" %in% names(x)) x$block_support_apps <- 0L
  
  if (!"block_catalog_hits" %in% names(x) && "catalog_hits" %in% names(x)) x$block_catalog_hits <- x$catalog_hits
  if (!"block_catalog_hits" %in% names(x) && "n_catalog_block_hits" %in% names(x)) x$block_catalog_hits <- x$n_catalog_block_hits
  if (!"block_catalog_hits" %in% names(x)) x$block_catalog_hits <- 0L
  
  if (!"block_gtex_hits" %in% names(x) && "gtex_hits" %in% names(x)) x$block_gtex_hits <- x$gtex_hits
  if (!"block_gtex_hits" %in% names(x) && "n_gtex_block_hits" %in% names(x)) x$block_gtex_hits <- x$n_gtex_block_hits
  if (!"block_gtex_hits" %in% names(x)) x$block_gtex_hits <- 0L
  
  if (!"block_nonsyn_hits" %in% names(x) && "nonsyn_hits" %in% names(x)) x$block_nonsyn_hits <- x$nonsyn_hits
  if (!"block_nonsyn_hits" %in% names(x) && "n_nonsyn_block_hits" %in% names(x)) x$block_nonsyn_hits <- x$n_nonsyn_block_hits
  if (!"block_nonsyn_hits" %in% names(x)) x$block_nonsyn_hits <- 0L
  
  if (!"block_ewasdis_hits" %in% names(x) && "ewasdis_hits" %in% names(x)) x$block_ewasdis_hits <- x$ewasdis_hits
  if (!"block_ewasdis_hits" %in% names(x) && "n_ewasdis_block_hits" %in% names(x)) x$block_ewasdis_hits <- x$n_ewasdis_block_hits
  if (!"block_ewasdis_hits" %in% names(x)) x$block_ewasdis_hits <- 0L
  
  if (!"block_ewastum_hits" %in% names(x) && "ewastum_hits" %in% names(x)) x$block_ewastum_hits <- x$ewastum_hits
  if (!"block_ewastum_hits" %in% names(x) && "n_ewastum_block_hits" %in% names(x)) x$block_ewastum_hits <- x$n_ewastum_block_hits
  if (!"block_ewastum_hits" %in% names(x)) x$block_ewastum_hits <- 0L
  
  if (!"gwas_sig_hits" %in% names(x) && "gwas_hits" %in% names(x)) x$gwas_sig_hits <- x$gwas_hits
  if (!"gwas_sig_hits" %in% names(x)) x$gwas_sig_hits <- ""
  if (!"max_gwas_logp" %in% names(x) && "gwas_top_logp" %in% names(x)) x$max_gwas_logp <- x$gwas_top_logp
  if (!"max_gwas_logp" %in% names(x)) x$max_gwas_logp <- NA_real_
  
  if (!"n_ld_proxy_hits" %in% names(x)) x$n_ld_proxy_hits <- 0L
  if (!"block_max_ld_value" %in% names(x) && "block_max_ld" %in% names(x)) x$block_max_ld_value <- x$block_max_ld
  if (!"block_mean_ld_value" %in% names(x) && "block_mean_ld" %in% names(x)) x$block_mean_ld_value <- x$block_mean_ld
  if (!"block_max_ld_value" %in% names(x)) x$block_max_ld_value <- NA_real_
  if (!"block_mean_ld_value" %in% names(x)) x$block_mean_ld_value <- NA_real_
  if (!"ld_proxy_snps" %in% names(x)) x$ld_proxy_snps <- ""
  
  if (!"genes_in_block" %in% names(x) && "genes_name" %in% names(x)) x$genes_in_block <- x$genes_name
  if (!"genes_in_block" %in% names(x)) x$genes_in_block <- ""
  
  x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = as.character(block_label),
      start = safe_num0(block_start),
      end = safe_num0(block_end),
      size_kb = dplyr::if_else(
        is.finite(suppressWarnings(as.numeric(block_size_kb))) &
          suppressWarnings(as.numeric(block_size_kb)) > 0,
        suppressWarnings(as.numeric(block_size_kb)),
        pmax(safe_num0(block_size_bp) / 1000, 1)
      ),
      support_apps = safe_int0(block_support_apps),
      catalog_hits = safe_int0(block_catalog_hits),
      gtex_hits = safe_int0(block_gtex_hits),
      nonsyn_hits = safe_int0(block_nonsyn_hits),
      ewasdis_hits = safe_int0(block_ewasdis_hits),
      ewastum_hits = safe_int0(block_ewastum_hits),
      gwas_hits = dplyr::coalesce(as.character(gwas_sig_hits), ""),
      gwas_top_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      gwas_min_p = dplyr::if_else(
        is.finite(gwas_top_logp),
        10^(-gwas_top_logp),
        NA_real_
      ),
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), ""),
      genes_name = dplyr::coalesce(as.character(genes_in_block), "")
    ) %>%
    dplyr::select(
      cluster_id,
      block_id,
      block_label,
      start,
      end,
      size_kb,
      support_apps,
      catalog_hits,
      gtex_hits,
      nonsyn_hits,
      ewasdis_hits,
      ewastum_hits,
      gwas_hits,
      gwas_top_logp,
      gwas_min_p,
      n_ld_proxy_hits,
      block_max_ld,
      block_mean_ld,
      ld_proxy_snps,
      genes_name
    )
}

build_subset_for_cluster_global <- function(
    cluster_row,
    bfile_ref,
    keep_path,
    plink_bin,
    workdir,
    tag = NULL
) {
  chr_sel <- suppressWarnings(as.integer(cluster_row$chr[1]))
  st      <- suppressWarnings(as.integer(cluster_row$start[1]))
  en      <- suppressWarnings(as.integer(cluster_row$end[1]))
  cid     <- as.character(cluster_row$cluster_id[1])
  
  if (!is.finite(chr_sel) || !is.finite(st) || !is.finite(en) || st >= en) {
    stop("Invalid cluster interval for cluster_id=", cid)
  }
  
  chr_lab <- chr_label_plink(chr_sel)
  tag <- tag %||% paste0(
    "ldglobal_",
    gsub("[^A-Za-z0-9_\\-]", "_", cid),
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
  
  subset_prefix <- file.path(workdir, paste0(tag, "_subset"))
  
  args_subset <- c(
    "--bfile", bfile_ref,
    "--chr", chr_lab,
    "--from-bp", st,
    "--to-bp", en,
    "--keep", keep_path,
    "--make-bed",
    "--out", subset_prefix
  )
  
  r1 <- run_plink(args_subset, subset_prefix, plink_bin)
  if (is.null(r1$status) || r1$status != 0) {
    if (length(r1$stdout)) message(paste(r1$stdout, collapse = "\n"))
    stop("PLINK subset failed for cluster_id=", cid)
  }
  
  bim_file <- paste0(subset_prefix, ".bim")
  if (!file.exists(bim_file)) {
    stop("Subset .bim not created for cluster_id=", cid)
  }
  
  bim <- utils::read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
  
  fl <- bim %>%
    dplyr::transmute(
      SNP = as.character(SNP),
      BP = suppressWarnings(as.integer(BP))
    ) %>%
    dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
    dplyr::distinct(SNP, BP, .keep_all = TRUE) %>%
    dplyr::arrange(BP, SNP) %>%
    dplyr::mutate(ix = dplyr::row_number())
  
  list(
    subset_prefix = subset_prefix,
    fl = fl,
    bim = fl %>% dplyr::select(SNP, BP),
    chr_sel = chr_sel,
    st = st,
    en = en,
    cid = cid,
    tag = tag
  )
}

annotate_block_hits_for_cluster <- function(
    cluster_row,
    candidates_df,
    block_df,
    hits_sub,
    proxy_tbl = NULL
) {
  empty_out <- tibble::tibble(
    cluster_id = character(),
    chr = integer(),
    lead_snp = character(),
    lead_pos = integer(),
    block_id = character(),
    block_label = character(),
    block_start = integer(),
    block_end = integer(),
    block_size_bp = integer(),
    source_app = character(),
    classe = character(),
    hit_id = character(),
    hit_rsid = character(),
    hit_pos = integer(),
    relation_to_lead = character()
  )
  
  if (!is.data.frame(block_df) || !nrow(block_df)) return(empty_out)
  if (!is.data.frame(hits_sub) || !nrow(hits_sub)) return(empty_out)
  
  chr_sel <- suppressWarnings(as.integer(cluster_row$chr[1]))
  cid <- as.character(cluster_row$cluster_id[1])
  
  cand_sub <- candidates_df %>%
    dplyr::mutate(
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      classe = normalize_classe(classe),
      source_app = if ("source_app" %in% names(.)) as.character(source_app) else NA_character_,
      hit_id = dplyr::coalesce(
        if ("id_hit" %in% names(.)) as.character(id_hit) else NA_character_,
        as.character(rsid)
      )
    ) %>%
    dplyr::filter(chr == chr_sel, is.finite(position))
  
  out_list <- lapply(seq_len(nrow(hits_sub)), function(k) {
    lead_rsid <- as.character(hits_sub$rsid[k])
    lead_pos  <- suppressWarnings(as.integer(hits_sub$position[k]))
    
    blk_one <- block_df %>%
      dplyr::filter(block_start <= lead_pos, block_end >= lead_pos)
    
    if (!nrow(blk_one)) return(NULL)
    
    proxy_set <- character(0)
    if (is.data.frame(proxy_tbl) && nrow(proxy_tbl)) {
      proxy_set <- proxy_tbl %>%
        dplyr::filter(query_hit == lead_rsid) %>%
        dplyr::pull(proxy_snp) %>%
        unique() %>%
        as.character()
    }
    
    dplyr::bind_rows(lapply(seq_len(nrow(blk_one)), function(bi) {
      bb <- blk_one[bi, , drop = FALSE]
      
      hh <- cand_sub %>%
        dplyr::filter(position >= bb$block_start[1], position <= bb$block_end[1])
      
      if (!nrow(hh)) return(NULL)
      
      hh %>%
        dplyr::mutate(
          cluster_id = cid,
          lead_snp = .env$lead_rsid,
          lead_pos = .env$lead_pos,
          block_id = bb$block_id[1],
          block_label = if ("block_label" %in% names(bb)) bb$block_label[1] else bb$block_id[1],
          block_start = bb$block_start[1],
          block_end = bb$block_end[1],
          block_size_bp = bb$block_size_bp[1],
          hit_rsid = rsid,
          hit_pos = position,
          relation_to_lead = dplyr::case_when(
            (!is.na(hit_rsid) & hit_rsid == lead_snp) | (!is.na(hit_pos) & hit_pos == lead_pos) ~ "DIRECT",
            !is.na(hit_rsid) & hit_rsid %in% proxy_set ~ "LD",
            TRUE ~ "BLOCK"
          )
        ) %>%
        dplyr::select(
          cluster_id, chr, lead_snp, lead_pos,
          block_id, block_label, block_start, block_end, block_size_bp,
          source_app, classe, hit_id, hit_rsid, hit_pos, relation_to_lead
        ) %>%
        dplyr::distinct()
    }))
  })
  
  out <- dplyr::bind_rows(out_list)
  
  if (!is.data.frame(out) || !nrow(out)) {
    return(empty_out)
  }
  
  if (!"lead_pos" %in% names(out)) out$lead_pos <- NA_integer_
  if (!"block_start" %in% names(out)) out$block_start <- NA_integer_
  if (!"hit_pos" %in% names(out)) out$hit_pos <- NA_integer_
  if (!"classe" %in% names(out)) out$classe <- NA_character_
  if (!"hit_rsid" %in% names(out)) out$hit_rsid <- NA_character_
  
  out <- out %>%
    dplyr::arrange(
      dplyr::coalesce(lead_pos, hit_pos),
      block_start,
      hit_pos,
      classe,
      hit_rsid
    )

  add_block_marker_support(out)
  
}


compute_ld_bundle_common <- function(
    cluster_row,
    candidates_df,
    bfile_ref,
    keep_path,
    plink_bin,
    workdir,
    pop,
    ld_metric = "R2",
    r2_min = 0.6,
    max_snps_interval = 400,
    compute_blocks = TRUE,
    append_log = NULL
) {
  ss <- build_subset_for_cluster_global(
    cluster_row = cluster_row,
    bfile_ref = bfile_ref,
    keep_path = keep_path,
    plink_bin = plink_bin,
    workdir = workdir
  )
  
  cid <- ss$cid
  fl <- ss$fl
  subset_prefix <- ss$subset_prefix
  bim_df <- ss$bim
  
  cand_sub <- select_candidates_for_cluster_common(candidates_df, cluster_row)
  
  thin_res <- thin_ld_subset_common(
    subset_prefix = subset_prefix,
    fl = fl,
    cand_sub = cand_sub,
    max_snps = max_snps_interval,
    workdir = workdir,
    tag = ss$tag,
    plink_bin = plink_bin,
    append_log = append_log
  )
  
  subset_prefix <- thin_res$subset_prefix
  fl <- thin_res$fl
  bim_df <- fl %>% dplyr::select(SNP, BP)
  
  hits_sub <- cand_sub %>%
    dplyr::filter(classe == "GWAS") %>%
    dplyr::distinct(cluster_id, rsid, position, .keep_all = TRUE) %>%
    dplyr::arrange(position, rsid)
  
  seeds_df <- if (nrow(hits_sub) > 0) {
    hits_sub %>%
      dplyr::rowwise() %>%
      dplyr::do({
        seed <- find_seed_snp(.$rsid, .$position, bim_df)
        dplyr::bind_cols(
          data.frame(
            cluster_id = as.character(.$cluster_id),
            query_hit = as.character(.$rsid),
            query_pos = as.integer(.$position),
            stringsAsFactors = FALSE
          ),
          seed
        )
      }) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()
  } else {
    tibble::tibble(
      cluster_id = character(),
      query_hit = character(),
      query_pos = integer(),
      seed_snp = character(),
      seed_pos = integer(),
      seed_type = character(),
      seed_dist_bp = integer()
    )
  }
  
  proxy_list <- if (nrow(hits_sub) > 0) {
    lapply(seq_len(nrow(hits_sub)), function(i) {
      hit_row <- hits_sub[i, , drop = FALSE]
      tryCatch(
        build_hit_proxy_table(
          hit_row = hit_row,
          bim_df = bim_df,
          subset_prefix = subset_prefix,
          plink_bin = plink_bin,
          workdir = workdir,
          metric = ld_metric
        ),
        error = function(e) NULL
      )
    })
  } else {
    list()
  }
  
  proxy_tbl <- dplyr::bind_rows(proxy_list)
  if (!is.data.frame(proxy_tbl) || nrow(proxy_tbl) == 0) {
    proxy_tbl <- tibble::tibble(
      cluster_id = character(),
      query_hit = character(),
      query_pos = integer(),
      seed_snp = character(),
      seed_pos = integer(),
      seed_type = character(),
      seed_dist_bp = integer(),
      proxy_snp = character(),
      proxy_pos = integer(),
      ld_value = numeric()
    )
  }
  
  proxy_tbl <- proxy_tbl %>%
    dplyr::filter(is.finite(ld_value), ld_value >= r2_min) %>%
    dplyr::arrange(query_hit, dplyr::desc(ld_value), proxy_pos, proxy_snp)
  
  blocks_raw <- NULL
  blocks_ij <- tibble::tibble(i = integer(), j = integer())
  block_ranges <- tibble::tibble()
  block_hits <- tibble::tibble()
  block_genes <- tibble::tibble()
  block_summary <- tibble::tibble()
  
  if (isTRUE(compute_blocks)) {
    blk_prefix <- file.path(workdir, paste0(ss$tag, "_blocks"))
    blk_args <- c(
      "--bfile", subset_prefix,
      "--blocks", "no-pheno-req",
      "--out", blk_prefix
    )
    
    rb <- run_plink(blk_args, blk_prefix, plink_bin)
    if (is.null(rb$status) || rb$status != 0) {
      stop("PLINK --blocks failed.")
    }
    
    blocks_raw <- read_plink_blocks(blk_prefix)
    blocks_ij <- blocks_det_to_ij(blocks_raw, fl %>% dplyr::select(SNP, BP))
    
    block_ranges <- blocks_ij_to_ranges_common(
      blocks_ij = blocks_ij,
      fl = fl,
      cluster_id = cid,
      chr = cluster_row$chr[1]
    )
    
    block_hits <- annotate_block_hits_for_cluster(
      cluster_row = cluster_row,
      candidates_df = cand_sub,
      block_df = block_ranges,
      hits_sub = hits_sub,
      proxy_tbl = proxy_tbl
    )

    
    block_genes <- build_block_gene_overlap(block_ranges)
    
    if (!is.data.frame(block_genes)) {
      block_genes <- tibble::tibble()
    }
    
    if (!"cluster_id" %in% names(block_genes)) block_genes$cluster_id <- NA_character_
    if (!"block_id" %in% names(block_genes)) block_genes$block_id <- NA_character_
    if (!"n_genes_in_block" %in% names(block_genes)) block_genes$n_genes_in_block <- 0L
    if (!"genes_in_block" %in% names(block_genes)) block_genes$genes_in_block <- ""
    
    block_genes <- block_genes %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        n_genes_in_block = safe_int0(n_genes_in_block),
        genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
      ) %>%
      dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
    

    bs_hits <- build_block_summary_from_hits(block_hits)
    
    if (!is.data.frame(bs_hits)) {
      bs_hits <- tibble::tibble()
    }
    
    bs_hits <- bs_hits %>%
      dplyr::select(
        -dplyr::any_of(c("n_genes_in_block", "genes_in_block"))
      )
    
    
    if (is.data.frame(block_ranges) && nrow(block_ranges)) {
      
      # provisional block summary stored in bundle details;
      # canonical scoring/display summary is built later with block_overlap_summary_df()
      block_summary <- block_ranges %>%
        dplyr::select(cluster_id, chr, block_id, block_label, block_start, block_end, block_size_bp) %>%
        dplyr::left_join(
          bs_hits,
          by = c("cluster_id", "block_id", "block_label", "block_start", "block_end", "block_size_bp")
        ) %>%
        dplyr::left_join(block_genes, by = c("cluster_id", "block_id")) %>%
        dplyr::mutate(
          n_genes_in_block = dplyr::coalesce(as.integer(n_genes_in_block), 0L),
          genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
        )
    }
    
    
  }
  
  list(
    cluster_id = cid,
    chr = ss$chr_sel,
    start = ss$st,
    end = ss$en,
    population = as.character(pop),
    ld_metric = as.character(ld_metric),
    subset_prefix = subset_prefix,
    fl = fl,
    bim = bim_df,
    candidates = cand_sub,
    gwas_hits = hits_sub,
    seeds = seeds_df,
    proxies = proxy_tbl,
    blocks = blocks_raw,
    blocks_ij = blocks_ij,
    block_ranges = block_ranges,
    block_hits = block_hits,
    block_genes = block_genes,
    block_summary = block_summary
  )
}

# funcio integradora:
thin_ld_subset_common <- function(
    subset_prefix,
    fl,
    cand_sub,
    max_snps,
    workdir,
    tag,
    plink_bin,
    append_log = NULL
) {
  stopifnot(is.data.frame(fl), all(c("SNP", "BP") %in% names(fl)))
  
  if (!is.finite(max_snps) || max_snps < 2 || nrow(fl) <= max_snps) {
    return(list(
      subset_prefix = subset_prefix,
      fl = fl %>% dplyr::arrange(BP, SNP) %>% dplyr::mutate(ix = dplyr::row_number())
    ))
  }
  
  if (is.function(append_log)) {
    append_log("[LD-THIN] Too many SNPs: ", nrow(fl), " > ", max_snps, ". Auto-thinning...")
  }
  
  keep_bp <- integer(0)
  
  if (is.data.frame(cand_sub) && nrow(cand_sub) > 0) {
    cand_keep <- cand_sub %>%
      dplyr::mutate(
        position = suppressWarnings(as.integer(position)),
        classe = as.character(classe),
        keep_rank = dplyr::case_when(
          classe == "GWAS" ~ 1L,
          classe == "catalog_hit" ~ 2L,
          classe == "gtex_hit" ~ 3L,
          classe == "nonsyn_hit" ~ 4L,
          classe == "ewasdis_hit" ~ 5L,
          classe == "ewastum_hit" ~ 6L,
          TRUE ~ 99L
        )
      ) %>%
      dplyr::filter(is.finite(position)) %>%
      dplyr::distinct(position, .keep_all = TRUE) %>%
      dplyr::arrange(keep_rank, position)
    
    if (nrow(cand_keep) > max_snps) {
      idx <- unique(round(seq(1, nrow(cand_keep), length.out = max_snps)))
      cand_keep <- cand_keep[idx, , drop = FALSE] %>%
        dplyr::arrange(keep_rank, position)
      
      if (is.function(append_log)) {
        append_log("[LD-THIN] Candidate SNPs exceed max limit. Keeping prioritised subset n=", nrow(cand_keep))
      }
    }
    
    keep_bp <- unique(as.integer(cand_keep$position))
  }
  
  keep_snps <- fl %>%
    dplyr::filter(BP %in% keep_bp) %>%
    dplyr::pull(SNP) %>%
    unique()
  
  k <- max_snps - length(keep_snps)
  
  rest <- fl %>%
    dplyr::filter(!(SNP %in% keep_snps)) %>%
    dplyr::arrange(BP, SNP)
  
  if (k > 0 && nrow(rest) > 0) {
    if (nrow(rest) <= k) {
      fill_snps <- rest$SNP
    } else {
      idx <- unique(pmax(1, pmin(nrow(rest), round(seq(1, nrow(rest), length.out = k)))))
      fill_snps <- rest$SNP[idx]
    }
  } else {
    fill_snps <- character(0)
  }
  
  snps_final <- unique(c(keep_snps, fill_snps))
  if (length(snps_final) > max_snps) {
    snps_final <- snps_final[seq_len(max_snps)]
  }
  
  extract_file <- file.path(workdir, paste0(tag, "_extract_snps.txt"))
  writeLines(snps_final, extract_file)
  
  subset2_prefix <- file.path(workdir, paste0(tag, "_subset_thin"))
  args_thin <- c(
    "--bfile", subset_prefix,
    "--extract", extract_file,
    "--make-bed",
    "--out", subset2_prefix
  )
  
  rthin <- run_plink(args_thin, subset2_prefix, plink_bin)
  if (is.null(rthin$status) || rthin$status != 0) {
    stop("PLINK thinning subset failed.")
  }
  
  bim2 <- utils::read.table(paste0(subset2_prefix, ".bim"), header = FALSE, stringsAsFactors = FALSE)
  colnames(bim2) <- c("CHR","SNP","CM","BP","A1","A2")
  
  fl2 <- bim2 %>%
    dplyr::transmute(SNP = as.character(SNP), BP = suppressWarnings(as.integer(BP))) %>%
    dplyr::filter(!is.na(SNP), nzchar(SNP), is.finite(BP)) %>%
    dplyr::distinct(SNP, BP, .keep_all = TRUE) %>%
    dplyr::arrange(BP, SNP) %>%
    dplyr::mutate(ix = dplyr::row_number())
  
  if (is.function(append_log)) {
    append_log("[LD-THIN] SNPs after thinning: ", nrow(fl2))
  }
  
  list(
    subset_prefix = subset2_prefix,
    fl = fl2
  )
}


select_candidates_for_cluster_common <- function(candidates_df, cluster_row) {
  stopifnot(is.data.frame(candidates_df), is.data.frame(cluster_row), nrow(cluster_row) > 0)
  
  cid <- as.character(cluster_row$cluster_id[1])
  chr_sel <- suppressWarnings(as.integer(cluster_row$chr[1]))
  st <- suppressWarnings(as.integer(cluster_row$start[1]))
  en <- suppressWarnings(as.integer(cluster_row$end[1]))
  
  ca_sub <- candidates_df %>%
    dplyr::mutate(
      cluster_id = normalize_cluster_id(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      classe = normalize_classe(classe)
    ) %>%
    dplyr::filter(
      chr == chr_sel,
      position >= st,
      position <= en
    )
  
  ca_sub2 <- ca_sub %>% dplyr::filter(cluster_id == cid)
  if (nrow(ca_sub2) > 0) ca_sub <- ca_sub2
  
  ca_sub %>%
    dplyr::distinct(source_app, cluster_id, id_hit, rsid, chr, position, classe, .keep_all = TRUE) %>%
    dplyr::arrange(position, classe, id_hit)
}

blocks_ij_to_ranges_common <- function(blocks_ij, fl, cluster_id, chr) {
  if (!is.data.frame(blocks_ij) || !nrow(blocks_ij)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      block_id = character(),
      block_label = character(),
      i = integer(),
      j = integer(),
      block_start = integer(),
      block_end = integer(),
      block_size_bp = integer()
    ))
  }
  
  blocks_ij %>%
    dplyr::mutate(
      i = suppressWarnings(as.integer(i)),
      j = suppressWarnings(as.integer(j))
    ) %>%
    dplyr::filter(
      is.finite(i), is.finite(j),
      i >= 1, j >= 1,
      i <= nrow(fl), j <= nrow(fl),
      j > i
    ) %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      block_start = as.integer(fl$BP[i]),
      block_end   = as.integer(fl$BP[j]),
      block_size_bp = as.integer(block_end - block_start),
      block_label = paste0("B", dplyr::row_number()),
      block_id = paste0(cluster_id, "_", block_label)
    ) %>%
    dplyr::select(cluster_id, chr, block_id, block_label, i, j, block_start, block_end, block_size_bp)
}


# final funcio integradora Global Cluster LD

rebuild_block_hits_from_components <- function(cl, ca, hits, blk, proxy_tbl = NULL) {
  empty_hits <- tibble::tibble(
    cluster_id = character(),
    chr = integer(),
    lead_snp = character(),
    lead_pos = integer(),
    block_id = character(),
    block_label = character(),
    block_start = integer(),
    block_end = integer(),
    block_size_bp = integer(),
    source_app = character(),
    classe = character(),
    hit_id = character(),
    hit_rsid = character(),
    hit_pos = integer(),
    relation_to_lead = character()
  )
  
  if (!is.data.frame(cl) || !nrow(cl) ||
      !is.data.frame(ca) || !nrow(ca) ||
      !is.data.frame(hits) || !nrow(hits) ||
      !is.data.frame(blk) || !nrow(blk)) {
    return(empty_hits)
  }
  
  cid <- as.character(cl$cluster_id[1])
  chr_sel <- suppressWarnings(as.integer(cl$chr[1]))
  
  ca2 <- ca %>%
    dplyr::mutate(
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      classe = as.character(classe),
      source_app = if ("source_app" %in% names(.)) as.character(source_app) else NA_character_,
      hit_id = dplyr::coalesce(
        if ("id_hit" %in% names(.)) as.character(id_hit) else NA_character_,
        as.character(rsid)
      )
    ) %>%
    dplyr::filter(chr == chr_sel, is.finite(position))
  
  out_list <- lapply(seq_len(nrow(hits)), function(k) {
    lead_rsid <- as.character(hits$rsid[k])
    lead_pos  <- suppressWarnings(as.integer(hits$position[k]))
    
    blk_one <- blk %>%
      dplyr::filter(block_start <= lead_pos, block_end >= lead_pos)
    
    if (!nrow(blk_one)) return(NULL)
    
    proxy_set <- character(0)
    if (is.data.frame(proxy_tbl) && nrow(proxy_tbl)) {
      proxy_set <- proxy_tbl %>%
        dplyr::filter(query_hit == lead_rsid) %>%
        dplyr::pull(proxy_snp) %>%
        unique() %>%
        as.character()
    }
    
    dplyr::bind_rows(lapply(seq_len(nrow(blk_one)), function(bi) {
      bb <- blk_one[bi, , drop = FALSE]
      
      hh <- ca2 %>%
        dplyr::filter(
          position >= bb$block_start[1],
          position <= bb$block_end[1],
          classe %in% c(
            "GWAS",
            "catalog_hit",
            "gtex_hit",
            "nonsyn_hit",
            "ewasdis_hit",
            "ewastum_hit"
          )
        )
      if (!nrow(hh)) return(NULL)
      
      hh %>%
        dplyr::mutate(
          cluster_id = cid,
          lead_snp = lead_rsid,
          lead_pos = lead_pos,
          block_id = bb$block_id[1],
          block_label = bb$block_label[1],
          block_start = bb$block_start[1],
          block_end = bb$block_end[1],
          block_size_bp = bb$block_size_bp[1],
          hit_rsid = rsid,
          hit_pos = position,
          relation_to_lead = dplyr::case_when(
            (!is.na(hit_rsid) & hit_rsid == lead_snp) | (!is.na(hit_pos) & hit_pos == lead_pos) ~ "DIRECT",
            !is.na(hit_rsid) & hit_rsid %in% proxy_set ~ "LD",
            TRUE ~ "BLOCK"
          )
        ) %>%
        dplyr::select(
          cluster_id, chr, lead_snp, lead_pos,
          block_id, block_label, block_start, block_end, block_size_bp,
          source_app, classe, hit_id, hit_rsid, hit_pos, relation_to_lead
        ) %>%
        dplyr::distinct()
    }))
  })
  
  out <- dplyr::bind_rows(out_list)
  if (!is.data.frame(out) || !nrow(out)) return(empty_hits)
  
  out %>%
    dplyr::arrange(lead_pos, block_start, hit_pos, classe, hit_rsid)
}

# ============================================================
# Block table display helpers
# ============================================================

safe_block_label_from_id <- function(block_id, cluster_id = NA_character_) {
  block_id <- as.character(block_id)
  cluster_id <- as.character(cluster_id)
  
  out <- block_id
  
  ok <- !is.na(block_id) & nzchar(block_id)
  
  # cas preferit: chr19_1_B47 -> B47
  out[ok] <- sub("^.*_(B[0-9]+)$", "\\1", block_id[ok])
  
  # fallback: chr19_1_block_47 -> B47
  idx_block <- ok & grepl("_block_[0-9]+$", block_id)
  out[idx_block] <- sub("^.*_block_([0-9]+)$", "B\\1", block_id[idx_block])
  
  # fallback final: si no s'ha pogut transformar, deixa block_id
  out
}

make_block_ld_table_display <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  x <- df
  
  if (!"gwas_sig_hits" %in% names(x)) x$gwas_sig_hits <- ""
  if (!"ld_proxy_snps" %in% names(x)) x$ld_proxy_snps <- ""
  if (!"genes_in_block" %in% names(x)) x$genes_in_block <- ""
  
  x %>%
    dplyr::mutate(
      start = safe_num0(block_start),
      end = safe_num0(block_end),
      size_kb = safe_num0(block_size_kb),
      support_any = safe_int0(block_support_any),
      support_apps = safe_int0(block_support_apps),
      support_hits = safe_int0(block_support_hits),
      hit_density = safe_num0(block_hit_density),
      catalog_hits = safe_int0(block_catalog_hits),
      gtex_hits = safe_int0(block_gtex_hits),
      nonsyn_hits = safe_int0(block_nonsyn_hits),
      ewasdis_hits = safe_int0(block_ewasdis_hits),
      ewastum_hits = safe_int0(block_ewastum_hits),
      gwas_top_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      gwas_mean_logp = suppressWarnings(as.numeric(mean_gwas_logp)),
      gwas_min_p = dplyr::if_else(
        is.finite(gwas_top_logp),
        10^(-gwas_top_logp),
        NA_real_
      ),
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld = suppressWarnings(as.numeric(block_mean_ld_value)),
      n_genes = safe_int0(n_genes_in_block),
      genes_name = dplyr::coalesce(as.character(genes_in_block), ""),
      gwas_sig_hits = dplyr::coalesce(as.character(gwas_sig_hits), ""),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
    ) %>%
    dplyr::select(
      cluster_id,
      block_id,
      block_label,
      start,
      end,
      size_kb,
      support_any,
      support_apps,
      support_hits,
      hit_density,
      catalog_hits,
      gtex_hits,
      nonsyn_hits,
      ewasdis_hits,
      ewastum_hits,
      gwas_top_logp,
      gwas_mean_logp,
      gwas_min_p,
      n_ld_proxy_hits,
      block_max_ld,
      block_mean_ld,
      n_genes,
      genes_name,
      gwas_sig_hits,
      ld_proxy_snps
    )
}

make_block_ld_merged_table_display <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  x <- df
  
  if (!"cluster_id" %in% names(x)) x$cluster_id <- NA_character_
  if (!"block_id" %in% names(x)) x$block_id <- NA_character_
  if (!"block_label" %in% names(x)) x$block_label <- NA_character_
  
  if (!"block_start" %in% names(x)) x$block_start <- NA_real_
  if (!"block_end" %in% names(x)) x$block_end <- NA_real_
  if (!"block_size_bp" %in% names(x)) x$block_size_bp <- NA_real_
  if (!"block_size_kb" %in% names(x)) x$block_size_kb <- NA_real_
  
  if (!"block_support_any" %in% names(x)) x$block_support_any <- 0L
  if (!"block_support_apps" %in% names(x)) x$block_support_apps <- 0L
  if (!"block_support_hits" %in% names(x)) x$block_support_hits <- 0L
  if (!"block_hit_density" %in% names(x)) x$block_hit_density <- 0
  
  if (!"block_catalog_hits" %in% names(x)) x$block_catalog_hits <- 0L
  if (!"block_gtex_hits" %in% names(x)) x$block_gtex_hits <- 0L
  if (!"block_nonsyn_hits" %in% names(x)) x$block_nonsyn_hits <- 0L
  if (!"block_ewasdis_hits" %in% names(x)) x$block_ewasdis_hits <- 0L
  if (!"block_ewastum_hits" %in% names(x)) x$block_ewastum_hits <- 0L
  
  if (!"n_gwas_sig_hits" %in% names(x)) x$n_gwas_sig_hits <- 0L
  if (!"max_gwas_logp" %in% names(x)) x$max_gwas_logp <- NA_real_
  if (!"mean_gwas_logp" %in% names(x)) x$mean_gwas_logp <- NA_real_
  if (!"gwas_sig_hits" %in% names(x)) x$gwas_sig_hits <- ""
  
  if (!"n_ld_proxy_hits" %in% names(x)) x$n_ld_proxy_hits <- 0L
  if (!"block_max_ld_value" %in% names(x)) x$block_max_ld_value <- NA_real_
  if (!"block_mean_ld_value" %in% names(x)) x$block_mean_ld_value <- NA_real_
  if (!"ld_proxy_snps" %in% names(x)) x$ld_proxy_snps <- ""
  
  if (!"n_genes_in_block" %in% names(x)) x$n_genes_in_block <- 0L
  if (!"genes_in_block" %in% names(x)) x$genes_in_block <- ""
  
  if (!"n_lead_snps" %in% names(x)) x$n_lead_snps <- 0L
  if (!"lead_snps" %in% names(x)) x$lead_snps <- ""
  if (!"lead_positions" %in% names(x)) x$lead_positions <- ""
  
  if (!"gwas_hits" %in% names(x)) x$gwas_hits <- ""
  if (!"catalog_hits" %in% names(x)) x$catalog_hits <- ""
  if (!"gtex_hits" %in% names(x)) x$gtex_hits <- ""
  if (!"nonsyn_hits" %in% names(x)) x$nonsyn_hits <- ""
  if (!"ewasdis_hits" %in% names(x)) x$ewasdis_hits <- ""
  if (!"ewastum_hits" %in% names(x)) x$ewastum_hits <- ""
  
  if (!"n_external_apps" %in% names(x)) x$n_external_apps <- 0L
  if (!"external_apps" %in% names(x)) x$external_apps <- ""
  if (!"marker_status" %in% names(x)) x$marker_status <- ""
  
  if (!"n_direct" %in% names(x)) x$n_direct <- 0L
  if (!"n_ld" %in% names(x)) x$n_ld <- 0L
  if (!"n_block" %in% names(x)) x$n_block <- 0L
  if (!"direct_hits" %in% names(x)) x$direct_hits <- ""
  if (!"ld_hits" %in% names(x)) x$ld_hits <- ""
  if (!"block_hits" %in% names(x)) x$block_hits <- ""
  if (!"source_apps" %in% names(x)) x$source_apps <- ""
  
  x <- x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = as.character(block_label)
    )
  
  idx_empty_lab <- is.na(x$block_label) | !nzchar(trimws(x$block_label))
  x$block_label[idx_empty_lab] <- safe_block_label_from_id(
    block_id = x$block_id[idx_empty_lab],
    cluster_id = x$cluster_id[idx_empty_lab]
  )
  
  x <- x %>%
    dplyr::mutate(
      block_label = gsub("_block_", "_B", block_label, ignore.case = TRUE),
      block_label = sub("^block_", "B", block_label, ignore.case = TRUE),
      block_label = sub("^.*_B([0-9]+)$", "B\\1", block_label),
      block_id = paste0(cluster_id, "_", block_label),
      
      start = safe_num0(block_start),
      end = safe_num0(block_end),
      size_kb = dplyr::if_else(
        is.finite(suppressWarnings(as.numeric(block_size_kb))) &
          suppressWarnings(as.numeric(block_size_kb)) > 0,
        suppressWarnings(as.numeric(block_size_kb)),
        pmax(safe_num0(block_size_bp) / 1000, 1)
      ),
      
      support_apps = safe_int0(block_support_apps),
      
      gwas_top_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      gwas_min_p = dplyr::if_else(
        is.finite(gwas_top_logp),
        10^(-gwas_top_logp),
        NA_real_
      ),
      
      block_max_ld = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), ""),
      
      genes_name = dplyr::coalesce(as.character(genes_in_block), ""),
      
      lead_positions = dplyr::coalesce(as.character(lead_positions), ""),
      
      catalog_hits = dplyr::coalesce(as.character(catalog_hits), ""),
      gtex_hits = dplyr::coalesce(as.character(gtex_hits), ""),
      nonsyn_hits = dplyr::coalesce(as.character(nonsyn_hits), ""),
      ewasdis_hits = dplyr::coalesce(as.character(ewasdis_hits), ""),
      ewastum_hits = dplyr::coalesce(as.character(ewastum_hits), ""),
      
      external_apps = dplyr::coalesce(as.character(external_apps), ""),
      marker_status = dplyr::coalesce(as.character(marker_status), ""),
      
      n_direct = safe_int0(n_direct),
      n_ld = safe_int0(n_ld),
      n_block = safe_int0(n_block),
      direct_hits = dplyr::coalesce(as.character(direct_hits), ""),
      ld_hits = dplyr::coalesce(as.character(ld_hits), ""),
      block_hits = dplyr::coalesce(as.character(block_hits), ""),
      source_apps = dplyr::coalesce(as.character(source_apps), "")
    ) %>%
    dplyr::select(
      cluster_id,
      block_id,
      block_label,
      start,
      end,
      size_kb,
      lead_positions,
      support_apps,
      catalog_hits,
      gtex_hits,
      nonsyn_hits,
      ewasdis_hits,
      ewastum_hits,
      external_apps,
      marker_status,
      gwas_hits,
      gwas_top_logp,
      gwas_min_p,
      block_max_ld,
      block_mean_ld,
      ld_proxy_snps,
      genes_name,
      n_direct,
      n_ld,
      n_block,
      direct_hits,
      ld_hits,
      block_hits,
      source_apps
    )
  
  x
}

make_block_extra_table_display <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  x <- df
  
  if (!"cluster_id" %in% names(x)) x$cluster_id <- NA_character_
  if (!"block_id" %in% names(x)) x$block_id <- NA_character_
  if (!"block_label" %in% names(x)) x$block_label <- NA_character_
  if (!"block_start" %in% names(x)) x$block_start <- NA_real_
  if (!"block_end" %in% names(x)) x$block_end <- NA_real_
  if (!"block_size_kb" %in% names(x)) x$block_size_kb <- NA_real_
  if (!"block_size_bp" %in% names(x)) x$block_size_bp <- NA_real_
  
  if (!"lead_positions" %in% names(x)) x$lead_positions <- ""
  if (!"external_apps" %in% names(x)) x$external_apps <- ""
  if (!"marker_status" %in% names(x)) x$marker_status <- ""
  
  if (!"n_direct" %in% names(x)) x$n_direct <- 0L
  if (!"n_ld" %in% names(x)) x$n_ld <- 0L
  if (!"n_block" %in% names(x)) x$n_block <- 0L
  
  if (!"direct_hits" %in% names(x)) x$direct_hits <- ""
  if (!"ld_hits" %in% names(x)) x$ld_hits <- ""
  if (!"block_hits" %in% names(x)) x$block_hits <- ""
  if (!"source_apps" %in% names(x)) x$source_apps <- ""
  
  x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = as.character(block_label),
      start = safe_num0(block_start),
      end = safe_num0(block_end),
      size_kb = dplyr::if_else(
        is.finite(suppressWarnings(as.numeric(block_size_kb))) &
          suppressWarnings(as.numeric(block_size_kb)) > 0,
        suppressWarnings(as.numeric(block_size_kb)),
        pmax(safe_num0(block_size_bp) / 1000, 1)
      ),
      lead_positions = dplyr::coalesce(as.character(lead_positions), ""),
      external_apps = dplyr::coalesce(as.character(external_apps), ""),
      marker_status = dplyr::coalesce(as.character(marker_status), ""),
      n_direct = safe_int0(n_direct),
      n_ld = safe_int0(n_ld),
      n_block = safe_int0(n_block),
      direct_hits = dplyr::coalesce(as.character(direct_hits), ""),
      ld_hits = dplyr::coalesce(as.character(ld_hits), ""),
      block_hits = dplyr::coalesce(as.character(block_hits), ""),
      source_apps = dplyr::coalesce(as.character(source_apps), "")
    ) %>%
    dplyr::select(
      cluster_id,
      block_id,
      block_label,
      start,
      end,
      size_kb,
      lead_positions,
      external_apps,
      marker_status,
      n_direct,
      n_ld,
      n_block,
      direct_hits,
      ld_hits,
      block_hits,
      source_apps
    )
}

# ============================================================
# SHARED block summary builder
# Put this OUTSIDE ld_integrator_module_server(...)
# Same level as make_block_table_display()
# ============================================================
# abans build_block_overlap_summary()
build_block_ld_summary <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble())
  }
  
  x <- df
  
  if (!"cluster_id" %in% names(x)) x$cluster_id <- NA_character_
  if (!"block_id" %in% names(x)) x$block_id <- NA_character_
  if (!"block_label" %in% names(x)) x$block_label <- NA_character_
  
  if (!"block_start" %in% names(x) && "start" %in% names(x)) x$block_start <- x$start
  if (!"block_end" %in% names(x) && "end" %in% names(x)) x$block_end <- x$end
  if (!"block_start" %in% names(x)) x$block_start <- NA_real_
  if (!"block_end" %in% names(x)) x$block_end <- NA_real_
  
  if (!"block_size_bp" %in% names(x) && "size_bp" %in% names(x)) x$block_size_bp <- x$size_bp
  if (!"block_size_bp" %in% names(x)) x$block_size_bp <- NA_real_
  
  if (!"block_size_kb" %in% names(x) && "size_kb" %in% names(x)) x$block_size_kb <- x$size_kb
  if (!"block_size_kb" %in% names(x)) x$block_size_kb <- NA_real_
  
  if (!"block_support_any" %in% names(x) && "support_any" %in% names(x)) x$block_support_any <- x$support_any
  if (!"block_support_any" %in% names(x) && "has_external_support" %in% names(x)) x$block_support_any <- as.integer(x$has_external_support)
  if (!"block_support_any" %in% names(x)) x$block_support_any <- 0L
  
  if (!"block_support_apps" %in% names(x) && "support_apps" %in% names(x)) x$block_support_apps <- x$support_apps
  if (!"block_support_apps" %in% names(x) && "n_block_apps" %in% names(x)) x$block_support_apps <- x$n_block_apps
  if (!"block_support_apps" %in% names(x) && "n_external_apps" %in% names(x)) x$block_support_apps <- x$n_external_apps
  if (!"block_support_apps" %in% names(x)) x$block_support_apps <- 0L
  
  if (!"block_support_hits" %in% names(x) && "support_hits" %in% names(x)) x$block_support_hits <- x$support_hits
  if (!"block_support_hits" %in% names(x) && "n_block_hits" %in% names(x)) x$block_support_hits <- x$n_block_hits
  if (!"block_support_hits" %in% names(x)) x$block_support_hits <- 0L
  
  if (!"block_hit_density" %in% names(x) && "hit_density" %in% names(x)) x$block_hit_density <- x$hit_density
  if (!"block_hit_density" %in% names(x)) x$block_hit_density <- 0
  
  if (!"block_catalog_hits" %in% names(x) && "catalog_hits" %in% names(x)) x$block_catalog_hits <- x$catalog_hits
  if (!"block_catalog_hits" %in% names(x) && "n_catalog_block_hits" %in% names(x)) x$block_catalog_hits <- x$n_catalog_block_hits
  if (!"block_catalog_hits" %in% names(x)) x$block_catalog_hits <- 0L
  
  if (!"block_gtex_hits" %in% names(x) && "gtex_hits" %in% names(x)) x$block_gtex_hits <- x$gtex_hits
  if (!"block_gtex_hits" %in% names(x) && "n_gtex_block_hits" %in% names(x)) x$block_gtex_hits <- x$n_gtex_block_hits
  if (!"block_gtex_hits" %in% names(x)) x$block_gtex_hits <- 0L
  
  if (!"block_nonsyn_hits" %in% names(x) && "nonsyn_hits" %in% names(x)) x$block_nonsyn_hits <- x$nonsyn_hits
  if (!"block_nonsyn_hits" %in% names(x) && "n_nonsyn_block_hits" %in% names(x)) x$block_nonsyn_hits <- x$n_nonsyn_block_hits
  if (!"block_nonsyn_hits" %in% names(x)) x$block_nonsyn_hits <- 0L
  
  if (!"block_ewasdis_hits" %in% names(x) && "ewasdis_hits" %in% names(x)) x$block_ewasdis_hits <- x$ewasdis_hits
  if (!"block_ewasdis_hits" %in% names(x) && "n_ewasdis_block_hits" %in% names(x)) x$block_ewasdis_hits <- x$n_ewasdis_block_hits
  if (!"block_ewasdis_hits" %in% names(x)) x$block_ewasdis_hits <- 0L
  
  if (!"block_ewastum_hits" %in% names(x) && "ewastum_hits" %in% names(x)) x$block_ewastum_hits <- x$ewastum_hits
  if (!"block_ewastum_hits" %in% names(x) && "n_ewastum_block_hits" %in% names(x)) x$block_ewastum_hits <- x$n_ewastum_block_hits
  if (!"block_ewastum_hits" %in% names(x)) x$block_ewastum_hits <- 0L
  
  if (!"n_genes_in_block" %in% names(x) && "n_block_genes" %in% names(x)) x$n_genes_in_block <- x$n_block_genes
  if (!"n_genes_in_block" %in% names(x) && "n_genes" %in% names(x)) x$n_genes_in_block <- x$n_genes
  if (!"n_genes_in_block" %in% names(x)) x$n_genes_in_block <- 0L
  
  if (!"genes_in_block" %in% names(x) && "genes_name" %in% names(x)) x$genes_in_block <- x$genes_name
  if (!"genes_in_block" %in% names(x)) x$genes_in_block <- ""
  
  if (!"n_gwas_sig_hits" %in% names(x)) x$n_gwas_sig_hits <- 0L
  if (!"max_gwas_logp" %in% names(x)) x$max_gwas_logp <- NA_real_
  if (!"mean_gwas_logp" %in% names(x)) x$mean_gwas_logp <- NA_real_
  if (!"gwas_sig_hits" %in% names(x)) x$gwas_sig_hits <- ""
  
  if (!"n_ld_proxy_hits" %in% names(x)) x$n_ld_proxy_hits <- 0L
  if (!"block_max_ld_value" %in% names(x)) x$block_max_ld_value <- NA_real_
  if (!"block_mean_ld_value" %in% names(x)) x$block_mean_ld_value <- NA_real_
  if (!"ld_proxy_snps" %in% names(x)) x$ld_proxy_snps <- ""
  
  x <- x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = as.character(block_label)
    )
  
  idx_empty_lab <- is.na(x$block_label) | !nzchar(trimws(x$block_label))
  x$block_label[idx_empty_lab] <- safe_block_label_from_id(
    block_id = x$block_id[idx_empty_lab],
    cluster_id = x$cluster_id[idx_empty_lab]
  )
  
  x %>%
    dplyr::mutate(
      block_label = gsub("_block_", "_B", block_label, ignore.case = TRUE),
      block_label = sub("^block_", "B", block_label, ignore.case = TRUE),
      block_label = sub("^.*_B([0-9]+)$", "B\\1", block_label),
      block_id = paste0(cluster_id, "_", block_label),
      block_start = safe_num0(block_start),
      block_end = safe_num0(block_end),
      block_size_bp = safe_num0(block_size_bp),
      block_size_kb = dplyr::if_else(
        is.finite(suppressWarnings(as.numeric(block_size_kb))) &
          suppressWarnings(as.numeric(block_size_kb)) > 0,
        suppressWarnings(as.numeric(block_size_kb)),
        pmax(block_size_bp / 1000, 1)
      ),
      block_support_any = safe_int0(block_support_any),
      block_support_apps = safe_int0(block_support_apps),
      block_support_hits = safe_int0(block_support_hits),
      block_hit_density = safe_num0(block_hit_density),
      block_catalog_hits = safe_int0(block_catalog_hits),
      block_gtex_hits = safe_int0(block_gtex_hits),
      block_nonsyn_hits = safe_int0(block_nonsyn_hits),
      block_ewasdis_hits = safe_int0(block_ewasdis_hits),
      block_ewastum_hits = safe_int0(block_ewastum_hits),
      n_gwas_sig_hits = safe_int0(n_gwas_sig_hits),
      max_gwas_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      mean_gwas_logp = suppressWarnings(as.numeric(mean_gwas_logp)),
      gwas_sig_hits = dplyr::coalesce(as.character(gwas_sig_hits), ""),
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), ""),
      n_genes_in_block = safe_int0(n_genes_in_block),
      genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
    ) %>%
    dplyr::select(
      cluster_id,
      block_id,
      block_label,
      block_start,
      block_end,
      block_size_bp,
      block_size_kb,
      block_support_any,
      block_support_apps,
      block_support_hits,
      block_hit_density,
      block_catalog_hits,
      block_gtex_hits,
      block_nonsyn_hits,
      block_ewasdis_hits,
      block_ewastum_hits,
      n_gwas_sig_hits,
      max_gwas_logp,
      mean_gwas_logp,
      gwas_sig_hits,
      n_ld_proxy_hits,
      block_max_ld_value,
      block_mean_ld_value,
      ld_proxy_snps,
      n_genes_in_block,
      genes_in_block
    )
}

# per garantir compatibilitat al Global LD
build_block_overlap_summary <- function(df) {
  build_block_ld_summary(df)
}

read_gwas_significance_bridge <- function(session_dir) {
  f <- file.path(session_dir, "gwas_significance_bridge.rds")
  
  if (!file.exists(f)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      position = integer(),
      rsid = character(),
      p_value = numeric(),
      logp = numeric()
    ))
  }
  
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  
  if (!is.data.frame(x) || !nrow(x)) {
    return(tibble::tibble(
      cluster_id = character(),
      chr = integer(),
      position = integer(),
      rsid = character(),
      p_value = numeric(),
      logp = numeric()
    ))
  }
  
  x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      p_value = suppressWarnings(as.numeric(p_value)),
      logp = suppressWarnings(as.numeric(logp))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(position)
    ) %>%
    dplyr::distinct()
}

summarize_gwas_significance_by_block <- function(block_df, gwas_bridge_df) {
  if (!is.data.frame(block_df) || !nrow(block_df) ||
      !is.data.frame(gwas_bridge_df) || !nrow(gwas_bridge_df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  if (!all(c("cluster_id", "chr", "block_id", "block_start", "block_end") %in% names(block_df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  bb <- block_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      block_id = as.character(block_id),
      block_start = suppressWarnings(as.integer(block_start)),
      block_end = suppressWarnings(as.integer(block_end))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(block_start), is.finite(block_end)
    )
  
  gg <- gwas_bridge_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      position = suppressWarnings(as.integer(position)),
      rsid = as.character(rsid),
      logp = suppressWarnings(as.numeric(logp))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      is.finite(chr), is.finite(position)
    )
  
  ov <- gg %>%
    dplyr::inner_join(
      bb %>% dplyr::select(cluster_id, chr, block_id, block_start, block_end),
      by = c("cluster_id", "chr"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      position >= block_start,
      position <= block_end
    ) %>%
    dplyr::distinct(cluster_id, block_id, rsid, position, .keep_all = TRUE)
  
  if (!nrow(ov)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_gwas_sig_hits = integer(),
      max_gwas_logp = numeric(),
      mean_gwas_logp = numeric(),
      gwas_sig_hits = character()
    ))
  }
  
  collapse_hits <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- x[nzchar(x)]
    paste(sort(x), collapse = "; ")
  }
  
  ov %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      n_gwas_sig_hits = dplyr::n_distinct(rsid[!is.na(rsid) & nzchar(rsid)]),
      max_gwas_logp = if (any(is.finite(logp))) max(logp, na.rm = TRUE) else NA_real_,
      mean_gwas_logp = if (any(is.finite(logp))) mean(logp[is.finite(logp)], na.rm = TRUE) else NA_real_,
      gwas_sig_hits = collapse_hits(rsid),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_gwas_sig_hits = safe_int0(n_gwas_sig_hits),
      max_gwas_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      mean_gwas_logp = suppressWarnings(as.numeric(mean_gwas_logp)),
      gwas_sig_hits = dplyr::coalesce(gwas_sig_hits, "")
    )
}

summarize_ld_by_block <- function(block_df, proxy_tbl) {
  if (!is.data.frame(block_df) || !nrow(block_df) ||
      !is.data.frame(proxy_tbl) || !nrow(proxy_tbl)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  summarize_ld_structure_by_block <- function(block_df, ld_pairs_df) {
    if (!is.data.frame(block_df) || !nrow(block_df) ||
        !is.data.frame(ld_pairs_df) || !nrow(ld_pairs_df)) {
      return(tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_ld_proxy_hits = integer(),
        block_max_ld_value = numeric(),
        block_mean_ld_value = numeric(),
        ld_proxy_snps = character()
      ))
    }
    
    req_block <- c("cluster_id", "block_id", "block_start", "block_end")
    if (!all(req_block %in% names(block_df))) {
      return(tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_ld_proxy_hits = integer(),
        block_max_ld_value = numeric(),
        block_mean_ld_value = numeric(),
        ld_proxy_snps = character()
      ))
    }
    
    req_ld <- c("SNP_A", "SNP_B", "value")
    if (!all(req_ld %in% names(ld_pairs_df))) {
      return(tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_ld_proxy_hits = integer(),
        block_max_ld_value = numeric(),
        block_mean_ld_value = numeric(),
        ld_proxy_snps = character()
      ))
    }
    
    # necessitem posicions. si ja hi són, perfecte; si no, no podem resumir
    has_pos <- all(c("BPA", "BPB") %in% names(ld_pairs_df))
    if (!has_pos) {
      return(tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_ld_proxy_hits = integer(),
        block_max_ld_value = numeric(),
        block_mean_ld_value = numeric(),
        ld_proxy_snps = character()
      ))
    }
    
    bb <- block_df %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        block_start = suppressWarnings(as.integer(block_start)),
        block_end = suppressWarnings(as.integer(block_end))
      ) %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(block_id), nzchar(block_id),
        is.finite(block_start), is.finite(block_end)
      )
    
    ll <- ld_pairs_df %>%
      dplyr::transmute(
        SNP_A = as.character(SNP_A),
        SNP_B = as.character(SNP_B),
        BPA = suppressWarnings(as.integer(BPA)),
        BPB = suppressWarnings(as.integer(BPB)),
        value = suppressWarnings(as.numeric(value))
      ) %>%
      dplyr::filter(
        !is.na(SNP_A), nzchar(SNP_A),
        !is.na(SNP_B), nzchar(SNP_B),
        is.finite(BPA), is.finite(BPB), is.finite(value)
      )
    
    if (!nrow(bb) || !nrow(ll)) {
      return(tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_ld_proxy_hits = integer(),
        block_max_ld_value = numeric(),
        block_mean_ld_value = numeric(),
        ld_proxy_snps = character()
      ))
    }
    
    collapse_hits <- function(x) {
      x <- unique(stats::na.omit(as.character(x)))
      x <- x[nzchar(x)]
      paste(sort(x), collapse = "; ")
    }
    
    out_list <- lapply(seq_len(nrow(bb)), function(i) {
      b0 <- bb[i, , drop = FALSE]
      
      ld_sub <- ll %>%
        dplyr::filter(
          BPA >= b0$block_start[1],
          BPA <= b0$block_end[1],
          BPB >= b0$block_start[1],
          BPB <= b0$block_end[1]
        )
      
      if (!nrow(ld_sub)) {
        return(data.frame(
          cluster_id = b0$cluster_id[1],
          block_id = b0$block_id[1],
          n_ld_proxy_hits = 0L,
          block_max_ld_value = NA_real_,
          block_mean_ld_value = NA_real_,
          ld_proxy_snps = "",
          stringsAsFactors = FALSE
        ))
      }
      
      snps <- sort(unique(c(ld_sub$SNP_A, ld_sub$SNP_B)))
      snps <- snps[!is.na(snps) & nzchar(snps)]
      
      data.frame(
        cluster_id = b0$cluster_id[1],
        block_id = b0$block_id[1],
        n_ld_proxy_hits = length(snps),
        block_max_ld_value = max(ld_sub$value, na.rm = TRUE),
        block_mean_ld_value = mean(ld_sub$value, na.rm = TRUE),
        ld_proxy_snps = collapse_hits(snps),
        stringsAsFactors = FALSE
      )
    })
    
    dplyr::bind_rows(out_list) %>%
      dplyr::mutate(
        n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
        block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
        block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
        ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
      )
  }
  
  req_block_cols <- c("cluster_id", "block_id", "block_start", "block_end")
  if (!all(req_block_cols %in% names(block_df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  req_proxy_cols <- c("cluster_id", "proxy_snp", "proxy_pos", "ld_value")
  if (!all(req_proxy_cols %in% names(proxy_tbl))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  bb <- block_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_start = suppressWarnings(as.integer(block_start)),
      block_end = suppressWarnings(as.integer(block_end))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      !is.na(block_id), nzchar(block_id),
      is.finite(block_start), is.finite(block_end)
    )
  
  pp <- proxy_tbl %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      proxy_snp = as.character(proxy_snp),
      proxy_pos = suppressWarnings(as.integer(proxy_pos)),
      ld_value = suppressWarnings(as.numeric(ld_value))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      !is.na(proxy_snp), nzchar(proxy_snp),
      is.finite(proxy_pos),
      is.finite(ld_value)
    ) %>%
    dplyr::distinct(cluster_id, proxy_snp, proxy_pos, .keep_all = TRUE)
  
  ov <- pp %>%
    dplyr::inner_join(
      bb,
      by = "cluster_id",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      proxy_pos >= block_start,
      proxy_pos <= block_end
    ) %>%
    dplyr::distinct(cluster_id, block_id, proxy_snp, proxy_pos, .keep_all = TRUE)
  
  if (!nrow(ov)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  collapse_hits <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- x[nzchar(x)]
    paste(sort(x), collapse = "; ")
  }
  
  ov %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      n_ld_proxy_hits = dplyr::n_distinct(proxy_snp),
      block_max_ld_value = if (any(is.finite(ld_value))) max(ld_value, na.rm = TRUE) else NA_real_,
      block_mean_ld_value = if (any(is.finite(ld_value))) mean(ld_value[is.finite(ld_value)], na.rm = TRUE) else NA_real_,
      ld_proxy_snps = collapse_hits(proxy_snp),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
    )
}

summarize_ld_structure_by_block <- function(block_df, ld_pairs_df) {
  if (!is.data.frame(block_df) || !nrow(block_df) ||
      !is.data.frame(ld_pairs_df) || !nrow(ld_pairs_df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  req_block <- c("cluster_id", "block_id", "block_start", "block_end")
  if (!all(req_block %in% names(block_df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  req_ld <- c("SNP_A", "SNP_B", "value")
  if (!all(req_ld %in% names(ld_pairs_df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  has_pos <- all(c("BPA", "BPB") %in% names(ld_pairs_df))
  if (!has_pos) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  bb <- block_df %>%
    dplyr::transmute(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_start = suppressWarnings(as.integer(block_start)),
      block_end = suppressWarnings(as.integer(block_end))
    ) %>%
    dplyr::filter(
      !is.na(cluster_id), nzchar(cluster_id),
      !is.na(block_id), nzchar(block_id),
      is.finite(block_start), is.finite(block_end)
    )
  
  ll <- ld_pairs_df %>%
    dplyr::transmute(
      SNP_A = as.character(SNP_A),
      SNP_B = as.character(SNP_B),
      BPA = suppressWarnings(as.integer(BPA)),
      BPB = suppressWarnings(as.integer(BPB)),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(
      !is.na(SNP_A), nzchar(SNP_A),
      !is.na(SNP_B), nzchar(SNP_B),
      is.finite(BPA), is.finite(BPB), is.finite(value)
    )
  
  if (!nrow(bb) || !nrow(ll)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    ))
  }
  
  collapse_hits <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- x[nzchar(x)]
    paste(sort(x), collapse = "; ")
  }
  
  out_list <- lapply(seq_len(nrow(bb)), function(i) {
    b0 <- bb[i, , drop = FALSE]
    
    ld_sub <- ll %>%
      dplyr::filter(
        BPA >= b0$block_start[1],
        BPA <= b0$block_end[1],
        BPB >= b0$block_start[1],
        BPB <= b0$block_end[1]
      )
    
    if (!nrow(ld_sub)) {
      return(data.frame(
        cluster_id = b0$cluster_id[1],
        block_id = b0$block_id[1],
        n_ld_proxy_hits = 0L,
        block_max_ld_value = NA_real_,
        block_mean_ld_value = NA_real_,
        ld_proxy_snps = "",
        stringsAsFactors = FALSE
      ))
    }
    
    snps <- sort(unique(c(ld_sub$SNP_A, ld_sub$SNP_B)))
    snps <- snps[!is.na(snps) & nzchar(snps)]
    
    data.frame(
      cluster_id = b0$cluster_id[1],
      block_id = b0$block_id[1],
      n_ld_proxy_hits = length(snps),
      block_max_ld_value = max(ld_sub$value, na.rm = TRUE),
      block_mean_ld_value = mean(ld_sub$value, na.rm = TRUE),
      ld_proxy_snps = collapse_hits(snps),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out_list) %>%
    dplyr::mutate(
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
    )
}

summarize_block_hit_links_for_block <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_direct = integer(),
      n_ld = integer(),
      n_block = integer(),
      direct_hits = character(),
      ld_hits = character(),
      block_hits = character(),
      source_apps = character()
    ))
  }
  
  if (!all(c("cluster_id", "block_id") %in% names(df))) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_direct = integer(),
      n_ld = integer(),
      n_block = integer(),
      direct_hits = character(),
      ld_hits = character(),
      block_hits = character(),
      source_apps = character()
    ))
  }
  
  x <- df
  
  if (!"source_app" %in% names(x)) x$source_app <- NA_character_
  if (!"classe" %in% names(x)) x$classe <- NA_character_
  if (!"hit_rsid" %in% names(x)) x$hit_rsid <- NA_character_
  if (!"hit_id" %in% names(x)) x$hit_id <- NA_character_
  if (!"relation_to_lead" %in% names(x)) x$relation_to_lead <- NA_character_
  
  collapse_semicolon_unique <- function(v) {
    v <- as.character(v)
    v <- trimws(v)
    v <- v[!is.na(v) & nzchar(v)]
    v <- unique(v)
    paste(sort(v), collapse = "; ")
  }
  
  x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      source_app = tolower(trimws(as.character(source_app))),
      source_app = dplyr::case_when(
        source_app %in% c("catalog") ~ "catalog",
        source_app %in% c("gtex") ~ "gtex",
        source_app %in% c("nonsyn", "dbnsfp", "non_syn", "non-syn") ~ "nonsyn",
        source_app %in% c("ewasdis", "ewas_disease", "disease") ~ "ewasdis",
        source_app %in% c("ewastum", "ewas_tumor", "tumor") ~ "ewastum",
        TRUE ~ source_app
      ),
      classe = normalize_classe(classe),
      relation_to_lead = as.character(relation_to_lead),
      hit_show = dplyr::coalesce(
        dplyr::na_if(as.character(hit_rsid), ""),
        dplyr::na_if(as.character(hit_id), "")
      ),
      is_external_hit = classe %in% c(
        "catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit"
      ),
      source_app_external = dplyr::case_when(
        classe == "catalog_hit" ~ "catalog",
        classe == "gtex_hit" ~ "gtex",
        classe == "nonsyn_hit" ~ "nonsyn",
        classe == "ewasdis_hit" ~ "ewasdis",
        classe == "ewastum_hit" ~ "ewastum",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      n_direct = sum(relation_to_lead == "DIRECT", na.rm = TRUE),
      n_ld     = sum(relation_to_lead == "LD", na.rm = TRUE),
      n_block  = sum(relation_to_lead == "BLOCK", na.rm = TRUE),
      direct_hits = collapse_semicolon_unique(hit_show[relation_to_lead == "DIRECT"]),
      ld_hits     = collapse_semicolon_unique(hit_show[relation_to_lead == "LD"]),
      block_hits  = collapse_semicolon_unique(hit_show[relation_to_lead == "BLOCK"]),
      source_apps = collapse_semicolon_unique(source_app_external[is_external_hit]),
      .groups = "drop"
    )
}


# Canonical block hit summary:
# compute support by normalized source_app directly from block_hits.
# Do not rely on classe-only helpers here, because source_app is the
# stable cross-app identifier used by the integrator.

# funcio priority score by tertiles
classify_priority_tertiles <- function(score) {
  score <- suppressWarnings(as.numeric(score))
  
  q33 <- stats::quantile(score, 0.33, na.rm = TRUE, names = FALSE, type = 7)
  q66 <- stats::quantile(score, 0.66, na.rm = TRUE, names = FALSE, type = 7)
  
  dplyr::case_when(
    is.na(score) ~ NA_character_,
    q33 == q66 & score >= q66 ~ "High",
    q33 == q66 ~ "Low",
    score >= q66 ~ "High",
    score >= q33 ~ "Medium",
    TRUE ~ "Low"
  )
}

block_overlap_summary_df <- function(
    block_ranges,
    block_hits,
    block_genes = NULL,
    gwas_bridge_df = NULL,
    proxy_tbl = NULL
) {
  # --------------------------------------------------
  # 0) base buida
  # --------------------------------------------------
  empty_out <- tibble::tibble(
    cluster_id = character(),
    block_id = character(),
    block_label = character(),
    block_priority_score = numeric(),
    priority_class = character(),
    block_start = numeric(),
    block_end = numeric(),
    block_size_bp = numeric(),
    block_size_kb = numeric(),
    block_support_any = integer(),
    block_support_apps = integer(),
    block_support_hits = integer(),
    block_hit_density = numeric(),
    block_catalog_hits = integer(),
    block_gtex_hits = integer(),
    block_nonsyn_hits = integer(),
    block_ewasdis_hits = integer(),
    block_ewastum_hits = integer(),
    n_gwas_sig_hits = integer(),
    max_gwas_logp = numeric(),
    mean_gwas_logp = numeric(),
    gwas_sig_hits = character(),
    n_genes_in_block = integer(),
    n_ld_proxy_hits = integer(),
    block_max_ld_value = numeric(),
    block_mean_ld_value = numeric(),
    ld_proxy_snps = character(),
    genes_in_block = character()
  )
  
  if (!is.data.frame(block_ranges) || !nrow(block_ranges)) {
    return(empty_out)
  }
  
  # --------------------------------------------------
  # 1) rangs base
  # --------------------------------------------------
  base_df <- block_ranges %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = dplyr::coalesce(as.character(block_label), ""),
      chr = suppressWarnings(as.integer(chr)),
      block_start = suppressWarnings(as.numeric(block_start)),
      block_end = suppressWarnings(as.numeric(block_end)),
      block_size_bp = suppressWarnings(as.numeric(block_size_bp))
    ) %>%
    dplyr::mutate(
      block_label = dplyr::if_else(
        !nzchar(block_label),
        safe_block_label_from_id(block_id, cluster_id),
        block_label
      ),
      block_size_bp = dplyr::if_else(
        is.finite(block_size_bp) & block_size_bp >= 0,
        block_size_bp,
        pmax(block_end - block_start, 0)
      )
    ) %>%
    dplyr::select(cluster_id, chr, block_id, block_label, block_start, block_end, block_size_bp) %>%
    dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
  
  # --------------------------------------------------
  # 2) resum hits únics per block
  # --------------------------------------------------
  hits_sum <- if (is.data.frame(block_hits) && nrow(block_hits)) {
    
    bh <- block_hits %>%
      dplyr::mutate(
        cluster_id = as.character(cluster_id),
        block_id   = as.character(block_id),
        source_app = tolower(trimws(as.character(source_app))),
        source_app = dplyr::case_when(
          source_app %in% c("catalog") ~ "catalog",
          source_app %in% c("gtex") ~ "gtex",
          source_app %in% c("nonsyn", "dbnsfp", "non_syn", "non-syn") ~ "nonsyn",
          source_app %in% c("ewasdis", "ewas_disease", "disease") ~ "ewasdis",
          source_app %in% c("ewastum", "ewas_tumor", "tumor") ~ "ewastum",
          TRUE ~ source_app
        )
      )
    
    pick_existing_col <- function(df, candidates, default = NA_character_) {
      nm <- candidates[candidates %in% names(df)][1]
      if (length(nm) == 0 || is.na(nm)) {
        return(rep(default, nrow(df)))
      }
      as.character(df[[nm]])
    }
    
    # columna de classe real si existeix
    classe_vec <- pick_existing_col(bh, c("classe", "class"), default = NA_character_)
    classe_vec <- tolower(trimws(classe_vec))
    classe_vec[classe_vec == ""] <- NA_character_
    
    # identificador únic del hit, amb fallback robust
    chr_vec      <- pick_existing_col(bh, c("chr", "CHR", "chrom", "chromosome"), default = NA_character_)
    pos_vec      <- pick_existing_col(bh, c("position", "pos", "bp", "BP"), default = NA_character_)
    pos_ini_vec  <- pick_existing_col(bh, c("pos_ini", "start"), default = NA_character_)
    pos_end_vec  <- pick_existing_col(bh, c("pos_end", "end"), default = NA_character_)
    id_hit_vec   <- pick_existing_col(bh, c("id_hit", "hit_id", "variant_id", "id"), default = NA_character_)
    rsid_vec     <- pick_existing_col(bh, c("rsid", "rsID", "snp", "SNP"), default = NA_character_)
    
    loc1 <- ifelse(
      !is.na(chr_vec) & nzchar(chr_vec) & !is.na(pos_vec) & nzchar(pos_vec),
      paste0(chr_vec, ":", pos_vec),
      NA_character_
    )
    
    loc2 <- ifelse(
      !is.na(chr_vec) & nzchar(chr_vec) &
        !is.na(pos_ini_vec) & nzchar(pos_ini_vec) &
        !is.na(pos_end_vec) & nzchar(pos_end_vec),
      paste0(chr_vec, ":", pos_ini_vec, ":", pos_end_vec),
      NA_character_
    )
    
    pick_first_nonempty <- function(...) {
      xs <- list(...)
      out <- rep(NA_character_, length(xs[[1]]))
      for (x in xs) {
        x <- as.character(x)
        x[trimws(x) == ""] <- NA_character_
        ok <- is.na(out) & !is.na(x)
        out[ok] <- trimws(x[ok])
      }
      out
    }
    
    hit_key_vec <- pick_first_nonempty(id_hit_vec, rsid_vec, loc1, loc2)
    
    bh <- bh %>%
      dplyr::mutate(
        classe = classe_vec,
        hit_key = hit_key_vec
      )
    
    # Regla robusta:
    # - si classe està informada, comptem només hits reals d'app
    # - si classe no existeix / és NA, no filtrem per classe
    bh_app <- bh %>%
      dplyr::filter(
        !is.na(cluster_id), nzchar(cluster_id),
        !is.na(block_id), nzchar(block_id),
        !is.na(source_app), nzchar(source_app),
        source_app %in% c("catalog", "gtex", "nonsyn", "ewasdis", "ewastum"),
        !is.na(hit_key), nzchar(hit_key)
      ) %>%
      dplyr::filter(
        (
          !is.na(classe) &
            (
              (source_app == "catalog"  & classe == "catalog_hit")  |
                (source_app == "gtex"     & classe == "gtex_hit")     |
                (source_app == "nonsyn"   & classe == "nonsyn_hit")   |
                (source_app == "ewasdis"  & classe == "ewasdis_hit")  |
                (source_app == "ewastum"  & classe == "ewastum_hit")
            )
        ) |
          is.na(classe)
      )
    
    bh_unique <- bh_app %>%
      dplyr::distinct(cluster_id, block_id, source_app, hit_key, .keep_all = FALSE)
    
    
    bh_unique %>%
      dplyr::group_by(cluster_id, block_id) %>%
      dplyr::summarise(
        n_block_hits = dplyr::n(),
        n_catalog_block_hits = sum(source_app == "catalog", na.rm = TRUE),
        n_gtex_block_hits    = sum(source_app == "gtex", na.rm = TRUE),
        n_nonsyn_block_hits  = sum(source_app == "nonsyn", na.rm = TRUE),
        n_ewasdis_block_hits = sum(source_app == "ewasdis", na.rm = TRUE),
        n_ewastum_block_hits = sum(source_app == "ewastum", na.rm = TRUE),
        n_block_apps = dplyr::n_distinct(source_app),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        n_block_hits = safe_int0(n_block_hits),
        n_catalog_block_hits = safe_int0(n_catalog_block_hits),
        n_gtex_block_hits = safe_int0(n_gtex_block_hits),
        n_nonsyn_block_hits = safe_int0(n_nonsyn_block_hits),
        n_ewasdis_block_hits = safe_int0(n_ewasdis_block_hits),
        n_ewastum_block_hits = safe_int0(n_ewastum_block_hits),
        n_block_apps = safe_int0(n_block_apps)
      )
    
  } else {
    tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_block_hits = integer(),
      n_catalog_block_hits = integer(),
      n_gtex_block_hits = integer(),
      n_nonsyn_block_hits = integer(),
      n_ewasdis_block_hits = integer(),
      n_ewastum_block_hits = integer(),
      n_block_apps = integer()
    )
  }
  # --------------------------------------------------
  # 3) gens per block
  # --------------------------------------------------
  gene_sum <- if (is.data.frame(block_genes) && nrow(block_genes)) {
    bg <- block_genes
    
    if (!"cluster_id" %in% names(bg)) bg$cluster_id <- NA_character_
    if (!"block_id" %in% names(bg)) bg$block_id <- NA_character_
    
    if ("n_genes_in_block" %in% names(bg) && "genes_in_block" %in% names(bg)) {
      bg %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          n_genes_in_block = safe_int0(n_genes_in_block),
          genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
        ) %>%
        dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
      
    } else {
      tmp_gene_sum <- summarize_genes_by_block_simple(bg)
      
      if (!is.data.frame(tmp_gene_sum)) {
        tmp_gene_sum <- tibble::tibble()
      }
      
      if (!"cluster_id" %in% names(tmp_gene_sum)) tmp_gene_sum$cluster_id <- NA_character_
      if (!"block_id" %in% names(tmp_gene_sum)) tmp_gene_sum$block_id <- NA_character_
      if (!"n_genes_in_block" %in% names(tmp_gene_sum)) tmp_gene_sum$n_genes_in_block <- 0L
      if (!"genes_in_block" %in% names(tmp_gene_sum)) tmp_gene_sum$genes_in_block <- ""
      
      cat("[DBG gene_sum tmp_gene_sum] names = ",
          paste(names(tmp_gene_sum), collapse = ", "), "\n", sep = "")
      cat("[DBG gene_sum tmp_gene_sum] nrow = ",
          if (is.data.frame(tmp_gene_sum)) nrow(tmp_gene_sum) else NA, "\n", sep = "")
      if (is.data.frame(tmp_gene_sum) && nrow(tmp_gene_sum)) {
        print(utils::head(tmp_gene_sum, 10))
      }
      
      tmp_gene_sum %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          n_genes_in_block = safe_int0(n_genes_in_block),
          genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
        ) %>%
        dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
    }
  } else {
    tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_genes_in_block = integer(),
      genes_in_block = character()
    )
  }
  # --------------------------------------------------
  # 4) GWAS significació per block
  # --------------------------------------------------
  gwas_sum <- tibble::tibble(
    cluster_id = character(),
    block_id = character(),
    n_gwas_sig_hits = integer(),
    max_gwas_logp = numeric(),
    mean_gwas_logp = numeric(),
    gwas_sig_hits = character()
  )
  
  # --------------------------------------------------
  # 4b) LD resum per block
  # --------------------------------------------------
  ld_sum <- if (is.data.frame(proxy_tbl) && nrow(proxy_tbl)) {
    summarize_ld_by_block(
      block_df = base_df,
      proxy_tbl = proxy_tbl
    ) %>%
      dplyr::transmute(
        cluster_id = as.character(cluster_id),
        block_id = as.character(block_id),
        n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
        block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
        block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
        ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
      )
  } else {
    tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      n_ld_proxy_hits = integer(),
      block_max_ld_value = numeric(),
      block_mean_ld_value = numeric(),
      ld_proxy_snps = character()
    )
  }
  
  # --------------------------------------------------
  # 5) unir-ho tot
  # --------------------------------------------------
  out0 <- base_df %>%
    dplyr::left_join(hits_sum, by = c("cluster_id", "block_id")) %>%
    dplyr::left_join(gene_sum, by = c("cluster_id", "block_id")) %>%
    dplyr::left_join(gwas_sum, by = c("cluster_id", "block_id")) %>%
    dplyr::left_join(ld_sum, by = c("cluster_id", "block_id"))
  
  # crear columnes si falten
  if (!"n_block_hits" %in% names(out0)) out0$n_block_hits <- 0L
  if (!"n_catalog_block_hits" %in% names(out0)) out0$n_catalog_block_hits <- 0L
  if (!"n_gtex_block_hits" %in% names(out0)) out0$n_gtex_block_hits <- 0L
  if (!"n_nonsyn_block_hits" %in% names(out0)) out0$n_nonsyn_block_hits <- 0L
  if (!"n_ewasdis_block_hits" %in% names(out0)) out0$n_ewasdis_block_hits <- 0L
  if (!"n_ewastum_block_hits" %in% names(out0)) out0$n_ewastum_block_hits <- 0L
  if (!"n_block_apps" %in% names(out0)) out0$n_block_apps <- 0L
  
  if (!"n_genes_in_block" %in% names(out0)) out0$n_genes_in_block <- 0L
  if (!"genes_in_block" %in% names(out0)) out0$genes_in_block <- ""
  
  if (!"n_gwas_sig_hits" %in% names(out0)) out0$n_gwas_sig_hits <- 0L
  if (!"max_gwas_logp" %in% names(out0)) out0$max_gwas_logp <- NA_real_
  if (!"mean_gwas_logp" %in% names(out0)) out0$mean_gwas_logp <- NA_real_
  if (!"gwas_sig_hits" %in% names(out0)) out0$gwas_sig_hits <- ""
  
  if (!"n_ld_proxy_hits" %in% names(out0)) out0$n_ld_proxy_hits <- 0L
  if (!"block_max_ld_value" %in% names(out0)) out0$block_max_ld_value <- NA_real_
  if (!"block_mean_ld_value" %in% names(out0)) out0$block_mean_ld_value <- NA_real_
  if (!"ld_proxy_snps" %in% names(out0)) out0$ld_proxy_snps <- ""
  
  out0 <- out0 %>%
    dplyr::mutate(
      n_block_hits = safe_int0(n_block_hits),
      n_catalog_block_hits = safe_int0(n_catalog_block_hits),
      n_gtex_block_hits = safe_int0(n_gtex_block_hits),
      n_nonsyn_block_hits = safe_int0(n_nonsyn_block_hits),
      n_ewasdis_block_hits = safe_int0(n_ewasdis_block_hits),
      n_ewastum_block_hits = safe_int0(n_ewastum_block_hits),
      n_block_apps = safe_int0(n_block_apps),
      n_genes_in_block = safe_int0(n_genes_in_block),
      genes_in_block = dplyr::coalesce(as.character(genes_in_block), ""),
      n_gwas_sig_hits = safe_int0(n_gwas_sig_hits),
      max_gwas_logp = suppressWarnings(as.numeric(max_gwas_logp)),
      mean_gwas_logp = suppressWarnings(as.numeric(mean_gwas_logp)),
      gwas_sig_hits = dplyr::coalesce(as.character(gwas_sig_hits), ""),
      n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
      block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
      block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
      ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
    )
  
  # --------------------------------------------------
  # 6) passar per la pipeline canònica única
  # --------------------------------------------------
  out1 <- build_block_overlap_summary(out0)
  
  if ("block_priority_score" %in% names(out1)) {
    out1 <- out1 %>%
      dplyr::mutate(
        priority_class = classify_priority_tertiles(block_priority_score)
      )
  }
  
  # --------------------------------------------------
  # 7) ordre final
  # --------------------------------------------------
  out1 %>%
    dplyr::arrange(
      cluster_id,
      block_start,
      block_id
    )
}

summarize_block_overlap_detail_for_dt <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(tibble::tibble(
      cluster_id = character(),
      block_id = character(),
      block_label = character(),
      block_start = integer(),
      block_end = integer(),
      block_size_bp = integer(),
      n_lead_snps = integer(),
      lead_snps = character(),
      lead_positions = character(),
      gwas_hits = character(),
      catalog_hits = character(),
      gtex_hits = character(),
      nonsyn_hits = character(),
      ewasdis_hits = character(),
      ewastum_hits = character(),
      n_external_apps = integer(),
      external_apps = character(),
      marker_status = character()
    ))
  }
  
  collapse_chr <- function(x) {
    x <- unique(stats::na.omit(as.character(x)))
    x <- trimws(x)
    x <- x[nzchar(x)]
    paste(x, collapse = "; ")
  }
  
  collapse_int <- function(x) {
    x <- unique(stats::na.omit(suppressWarnings(as.integer(x))))
    x <- x[is.finite(x)]
    paste(sort(x), collapse = "; ")
  }
  
  df %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      block_id = as.character(block_id),
      block_label = as.character(block_label),
      lead_snp = as.character(lead_snp),
      lead_pos = suppressWarnings(as.integer(lead_pos)),
      gwas_hits = dplyr::coalesce(as.character(gwas_hits), ""),
      catalog_hits = dplyr::coalesce(as.character(catalog_hits), ""),
      gtex_hits = dplyr::coalesce(as.character(gtex_hits), ""),
      nonsyn_hits = dplyr::coalesce(as.character(nonsyn_hits), ""),
      ewasdis_hits = dplyr::coalesce(as.character(ewasdis_hits), ""),
      ewastum_hits = dplyr::coalesce(as.character(ewastum_hits), ""),
      n_external_apps = safe_int0(n_external_apps),
      external_apps = dplyr::coalesce(as.character(external_apps), ""),
      marker_status = dplyr::coalesce(as.character(marker_status), "")
    ) %>%
    dplyr::group_by(cluster_id, block_id) %>%
    dplyr::summarise(
      block_label = dplyr::first(block_label),
      block_start = dplyr::first(block_start),
      block_end = dplyr::first(block_end),
      block_size_bp = dplyr::first(block_size_bp),
      n_lead_snps = dplyr::n_distinct(lead_snp[!is.na(lead_snp) & nzchar(lead_snp)]),
      lead_snps = collapse_chr(lead_snp),
      lead_positions = collapse_int(lead_pos),
      gwas_hits = dplyr::first(gwas_hits),
      catalog_hits = dplyr::first(catalog_hits),
      gtex_hits = dplyr::first(gtex_hits),
      nonsyn_hits = dplyr::first(nonsyn_hits),
      ewasdis_hits = dplyr::first(ewasdis_hits),
      ewastum_hits = dplyr::first(ewastum_hits),
      n_external_apps = max(n_external_apps, na.rm = TRUE),
      external_apps = dplyr::first(external_apps),
      marker_status = dplyr::first(marker_status),
      .groups = "drop"
    ) %>%
    dplyr::arrange(block_start, block_id)
}

collapse_semicolon_unique <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- unique(x)
  paste(sort(x), collapse = "; ")
}

collapse_items_html <- function(x, max_visible = 5L) {
  vals <- split_semicolon(x)
  
  if (!length(vals)) return("")
  
  vals_html <- html_escape_basic(vals)
  
  if (length(vals_html) <= max_visible) {
    return(paste(vals_html, collapse = "; "))
  }
  
  first_part <- paste(vals_html[seq_len(max_visible)], collapse = "; ")
  rest_part  <- paste(vals_html[(max_visible + 1):length(vals_html)], collapse = "; ")
  
  paste0(
    "<details>",
    "<summary>", length(vals_html), " elements</summary>",
    "<div style='margin-top:6px;'>",
    first_part,
    if (nzchar(rest_part)) paste0("; ", rest_part) else "",
    "</div>",
    "</details>"
  )
}

format_block_hits_grouped_for_dt <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  
  x <- df
  
  if (!"cluster_id" %in% names(x)) x$cluster_id <- NA_character_
  if (!"chr" %in% names(x)) x$chr <- NA_integer_
  if (!"lead_snp" %in% names(x)) x$lead_snp <- NA_character_
  if (!"lead_pos" %in% names(x)) x$lead_pos <- NA_integer_
  if (!"block_id" %in% names(x)) x$block_id <- NA_character_
  if (!"block_label" %in% names(x)) x$block_label <- NA_character_
  if (!"block_start" %in% names(x)) x$block_start <- NA_integer_
  if (!"block_end" %in% names(x)) x$block_end <- NA_integer_
  if (!"block_size_bp" %in% names(x)) x$block_size_bp <- NA_integer_
  if (!"source_app" %in% names(x)) x$source_app <- NA_character_
  if (!"classe" %in% names(x)) x$classe <- NA_character_
  if (!"hit_id" %in% names(x)) x$hit_id <- NA_character_
  if (!"hit_rsid" %in% names(x)) x$hit_rsid <- NA_character_
  if (!"hit_pos" %in% names(x)) x$hit_pos <- NA_integer_
  if (!"relation_to_lead" %in% names(x)) x$relation_to_lead <- NA_character_
  
  x <- x %>%
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      chr = suppressWarnings(as.integer(chr)),
      lead_snp = as.character(lead_snp),
      lead_pos = suppressWarnings(as.integer(lead_pos)),
      block_id = as.character(block_id),
      block_label = as.character(block_label),
      block_start = suppressWarnings(as.integer(block_start)),
      block_end = suppressWarnings(as.integer(block_end)),
      block_size_bp = suppressWarnings(as.integer(block_size_bp)),
      source_app = as.character(source_app),
      classe = as.character(classe),
      hit_id = as.character(hit_id),
      hit_rsid = as.character(hit_rsid),
      hit_pos = suppressWarnings(as.integer(hit_pos)),
      relation_to_lead = as.character(relation_to_lead),
      hit_show = dplyr::coalesce(
        dplyr::na_if(hit_rsid, ""),
        dplyr::na_if(hit_id, "")
      )
    )
  
  out <- x %>%
    dplyr::group_by(
      cluster_id, chr, lead_snp, lead_pos,
      block_id, block_label, block_start, block_end, block_size_bp
    ) %>%
    dplyr::summarise(
      direct_hits = collapse_semicolon_unique(hit_show[relation_to_lead == "DIRECT"]),
      ld_hits     = collapse_semicolon_unique(hit_show[relation_to_lead == "LD"]),
      block_hits  = collapse_semicolon_unique(hit_show[relation_to_lead == "BLOCK"]),
      
      gwas_hits     = collapse_semicolon_unique(hit_show[classe == "GWAS"]),
      catalog_hits  = collapse_semicolon_unique(hit_show[classe == "catalog_hit"]),
      gtex_hits     = collapse_semicolon_unique(hit_show[classe == "gtex_hit"]),
      nonsyn_hits   = collapse_semicolon_unique(hit_show[classe == "nonsyn_hit"]),
      ewasdis_hits  = collapse_semicolon_unique(hit_show[classe == "ewasdis_hit"]),
      ewastum_hits  = collapse_semicolon_unique(hit_show[classe == "ewastum_hit"]),
      
      source_apps = collapse_semicolon_unique(
        source_app[
          classe %in% c("catalog_hit", "gtex_hit", "nonsyn_hit", "ewasdis_hit", "ewastum_hit")
        ]
      ),
      
      n_direct = sum(relation_to_lead == "DIRECT", na.rm = TRUE),
      n_ld     = sum(relation_to_lead == "LD", na.rm = TRUE),
      n_block  = sum(relation_to_lead == "BLOCK", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      direct_hits   = vapply(direct_hits, collapse_items_html, character(1)),
      ld_hits       = vapply(ld_hits, collapse_items_html, character(1)),
      block_hits    = vapply(block_hits, collapse_items_html, character(1)),
      gwas_hits     = vapply(gwas_hits, collapse_items_html, character(1)),
      catalog_hits  = vapply(catalog_hits, collapse_items_html, character(1)),
      gtex_hits     = vapply(gtex_hits, collapse_items_html, character(1)),
      nonsyn_hits   = vapply(nonsyn_hits, collapse_items_html, character(1)),
      ewasdis_hits  = vapply(ewasdis_hits, collapse_items_html, character(1)),
      ewastum_hits  = vapply(ewastum_hits, collapse_items_html, character(1)),
      source_apps   = vapply(source_apps, collapse_items_html, character(1))
    ) %>%
    dplyr::arrange(lead_pos, block_start, block_id)
  
  out
}

compact_shared_terms <- function(x, n_show = 5) {
  x <- as.character(x %||% "")
  x <- trimws(x)
  
  empty_out <- list(
    terms_clean = "",
    terms_short = "",
    n_terms_clean = 0L
  )
  
  if (!nzchar(x)) return(empty_out)
  
  # --------------------------------------------------
  # 1) split raw terms
  # --------------------------------------------------
  parts <- unlist(strsplit(x, "|", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  
  if (!length(parts)) return(empty_out)
  
  # --------------------------------------------------
  # 2) helper: clean one term
  # --------------------------------------------------
  clean_one <- function(s) {
    s <- as.character(s %||% "")
    s <- trimws(s)
    
    # remove PMID/citation blocks like [12345;67890]
    s <- gsub("\\[[^\\]]*\\]", "", s)
    
    # remove repeated spaces
    s <- gsub("\\s+", " ", s)
    
    # remove common noisy suffixes
    s <- gsub("\\bmeasurement\\b", "", s, ignore.case = TRUE)
    s <- gsub("\\blevels\\b", "", s, ignore.case = TRUE)
    s <- gsub("\\blevel\\b", "", s, ignore.case = TRUE)
    s <- gsub("\\bamount\\b", "", s, ignore.case = TRUE)
    s <- gsub("\\bquantity\\b", "", s, ignore.case = TRUE)
    
    # remove hanging punctuation / separators
    s <- gsub("\\s+,", ",", s)
    s <- gsub(",\\s*$", "", s)
    s <- gsub(";\\s*$", "", s)
    s <- gsub("\\s+\\)", ")", s)
    s <- gsub("\\(\\s+", "(", s)
    
    # remove empty parentheses left after cleanup
    s <- gsub("\\(\\s*\\)", "", s)
    
    # collapse spaces again
    s <- gsub("\\s+", " ", s)
    s <- trimws(s)
    
    s
  }
  
  cleaned <- vapply(parts, clean_one, character(1))
  cleaned <- cleaned[nzchar(cleaned)]
  
  if (!length(cleaned)) return(empty_out)
  
  # --------------------------------------------------
  # 3) canonical key for deduplication
  # --------------------------------------------------
  canonical_one <- function(s) {
    z <- tolower(s)
    z <- trimws(z)
    
    # normalize punctuation spacing
    z <- gsub("\\s+", " ", z)
    z <- gsub("[“”]", "\"", z)
    z <- gsub("[‘’]", "'", z)
    
    # soften some common variants
    z <- gsub("\\balzheimer disease\\b", "alzheimer's disease", z)
    z <- gsub("\\blate onset\\b", "late-onset", z)
    z <- gsub("\\bcolour\\b", "color", z)
    
    z <- trimws(z)
    z
  }
  
  canon <- vapply(cleaned, canonical_one, character(1))
  
  # keep first presentable version for each canonical term
  keep_idx <- !duplicated(canon)
  clean_unique <- cleaned[keep_idx]
  canon_unique <- canon[keep_idx]
  
  # optional stable sort by canonical form
  ord <- order(canon_unique)
  clean_unique <- clean_unique[ord]
  
  n_clean <- length(clean_unique)
  
  # --------------------------------------------------
  # 4) short display
  # --------------------------------------------------
  if (n_clean <= n_show) {
    short_txt <- paste(clean_unique, collapse = " | ")
  } else {
    short_txt <- paste0(
      paste(clean_unique[seq_len(n_show)], collapse = " | "),
      " … (+", n_clean - n_show, " més)"
    )
  }
  
  list(
    terms_clean = paste(clean_unique, collapse = " | "),
    terms_short = short_txt,
    n_terms_clean = as.integer(n_clean)
  )
}

shared_terms_compact_df <- function(df, n_show = 5) {
  if (!is.data.frame(df) || !nrow(df) || !"shared_terms" %in% names(df)) {
    return(df)
  }
  
  cc <- lapply(df$shared_terms, compact_shared_terms)
  
  df$shared_terms_clean <- vapply(cc, `[[`, character(1), "terms_clean")
  df$shared_terms_short <- vapply(cc, `[[`, character(1), "terms_short")
  df$n_terms_clean <- vapply(cc, `[[`, integer(1), "n_terms_clean")
  
  df
}

summarize_shared_terms_for_dt <- function(x, max_terms = 6) {
  x <- as.character(x %||% "")
  x <- trimws(x)
  
  if (!nzchar(x) || identical(x, "NA")) {
    return("")
  }
  
  parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  
  if (!length(parts)) {
    return("")
  }
  
  parts_unique <- unique(parts)
  n_all <- length(parts_unique)
  
  shown <- head(parts_unique, max_terms)
  txt <- paste(shown, collapse = "; ")
  
  if (n_all > max_terms) {
    paste0(
      n_all, " terms: ",
      txt,
      "; ... (+", n_all - max_terms, " més)"
    )
  } else {
    paste0(n_all, " terms: ", txt)
  }
}

# helper link GeneCard
make_genecards_link <- function(x) {
  x2 <- trimws(as.character(x))
  if (!nzchar(x2)) return("")
  
  url <- paste0(
    "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
    utils::URLencode(x2, reserved = TRUE)
  )
  
  paste0("<a href='", url, "' target='_blank'>", html_escape_basic(x2), "</a>")
}

build_gene_link_column_html <- function(genes_str, max_visible = 5L) {
  genes <- split_semicolon(genes_str)
  if (!length(genes)) return("")
  
  links <- vapply(genes, make_genecards_link, character(1))
  
  collapse_link_list_html(
    links = links,
    max_visible = max_visible,
    summary_label = paste0(length(links), " genes")
  )
}

build_snp_link_column_html <- function(snps_str, max_visible = 5L) {
  snps <- split_semicolon(snps_str)
  if (!length(snps)) return("")
  
  links <- vapply(snps, make_dbsnp_link, character(1))
  
  collapse_link_list_html(
    links = links,
    max_visible = max_visible,
    summary_label = paste0(length(links), " SNPs")
  )
}

# ------------------------------------------------------------------
# UI
# ------------------------------------------------------------------
ld_integrator_module_ui <- function(id) {
  ns <- shiny::NS(id)
  ld_def <- gi_ld_defaults()
  
  shiny::tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        
        shiny::tags$details(
          shiny::tags$summary("Advanced LD paths"),
          shiny::wellPanel(
            shiny::textInput(ns("plink_bin"), "PLINK binary path", value = ld_def$plink %||% ""),
            shiny::textInput(ns("bfile_ref_prefix"), "Reference bfile prefix", value = ld_def$bfile %||% ""),
            shiny::textInput(ns("ld_pops_dir"), "POP keep-files directory (.txt)", value = ld_def$popdir %||% ""),
            shiny::tags$div(
              style = "font-size:12px; color:#666; margin-top:6px;",
              "Integrator session folder is controlled by the main Integrator app."
            )
          )
        ),
        
        shiny::uiOutput(ns("ld_pop_ui")),
        shiny::tags$hr(),
        
        shiny::selectInput(ns("ld_chr"), "Chromosome", choices = character(0)),
        shiny::selectInput(ns("ld_cluster_id"), "Integrated cluster", choices = character(0)),
        shiny::tags$hr(),
        
        shiny::radioButtons(
          ns("ld_metric"), "LD metric",
          choices = c("r²" = "R2", "D′" = "Dprime"),
          selected = "R2", inline = TRUE
        ),
        
        shiny::radioButtons(
          ns("ld_x_mode"), "X spacing",
          choices = c("Genomic (BP)" = "bp", "Equal spacing" = "equal"),
          selected = "bp"
        ),
        
        shiny::numericInput(ns("ld_max_snps_interval"), "Max SNPs in interval (plot safety)", value = 400, min = 50, step = 50),
        shiny::numericInput(ns("ld_proxy_r2_min"), "Min r² for priority proxies", value = 0.6, min = 0, max = 1, step = 0.05),
        
        shiny::tags$hr(),
        shiny::actionButton(ns("run_cluster_ld"), "Compute LD & blocks", icon = shiny::icon("play")),
        shiny::tags$hr(),
        
        shiny::verbatimTextOutput(ns("ld_log"))
      ),
      
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          
          shiny::tabPanel(
            "LD plot",
            
            shinycssloaders::withSpinner(
              plotly::plotlyOutput(ns("ld_plot"), height = "800px")
            ),
            
            shiny::tags$hr(),
            
            shiny::h4("LD block summary"),
            # shiny::h5("See 'Prioritized LD blocks' at 'Prioritize evidence'"),
            shiny::tags$p(
              style = "font-size:13px; color:#555; margin-bottom:8px;",
              "One row per LD block, summarizing LD structure, integrated hits, GWAS signal and genes within the block."
            ),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("block_overlap_summary"))
            ),
            
            shiny::tags$hr(),
            
            shiny::h4("Block-level hit details"),
            shiny::tags$p(
              style = "font-size:13px; color:#555; margin-bottom:8px;",
              "Detailed lead/block hit table showing DIRECT, LD and BLOCK relationships for the selected cluster."
            ),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("ld_block_hits_dt"))
            )
          ),
          
          shiny::tabPanel(
            "Seed & proxy summary",
            shiny::h4("GWAS seed summary"),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("ld_seed_dt"))
            ),
            shiny::tags$hr(),
            shiny::h4("Proxy table"),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("ld_proxy_dt"))
            )
          ),
          
          shiny::tabPanel(
            "Debug LD pairs",
            shiny::h4("SNPs in interval"),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("ld_fl_dt"))
            ),
            shiny::tags$hr(),
            shiny::h4("Top LD pairs"),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("ld_pairs_dt"))
            )
          )
        )
      )
    )
  )
}

# ------------------------------------------------------------------
# SERVER
# ------------------------------------------------------------------
ld_integrator_module_server <- function(
    id,
    activate_r = NULL,
    integrator_dir_r = NULL,
    integrated_clusters_r = NULL,
    integrated_candidates_r = NULL,
    gwas_bridge_r = NULL,
    apps = c("catalog", "gtex", "nonsyn", "ewastum", "ewasdis"),
    ...
) {
  moduleServer(id, function(input, output, session) {
    
    clusters_base_r <- shiny::reactive({
      if (!is.null(integrated_clusters_r)) {
        x <- tryCatch(integrated_clusters_r(), error = function(e) NULL)
        if (is.data.frame(x)) return(x)
      }
      NULL
    })
    
    candidates_base_r <- shiny::reactive({
      if (!is.null(integrated_candidates_r)) {
        x <- tryCatch(integrated_candidates_r(), error = function(e) NULL)
        if (is.data.frame(x)) return(x)
      }
      NULL
    })
    
    gwas_bridge_base_r <- shiny::reactive({
      if (!is.null(gwas_bridge_r)) {
        x <- tryCatch(gwas_bridge_r(), error = function(e) NULL)
        if (is.data.frame(x)) return(x)
      }
      tibble::tibble(
        cluster_id = character(),
        chr = integer(),
        position = integer(),
        rsid = character(),
        p_value = numeric(),
        logp = numeric()
      )
    })

    # alias locals de compatibilitat
    integrated_clusters <- clusters_base_r
    integrated_candidates <- candidates_base_r
    gwas_bridge_df <- gwas_bridge_base_r
    
    
    log_buf <- character(0)
    ld_log  <- shiny::reactiveVal("")
    
    append_log <- function(...) {
      txt <- paste(..., collapse = "")
      cat(txt, "\n")
      log_buf <<- c(log_buf, txt)
      ld_log(paste(log_buf, collapse = "\n"))
    }
    
    output$ld_log <- shiny::renderText(ld_log())
    
    append_log("[INIT] ld_integrator_module_server started")
    
    is_active <- shiny::reactive({
      if (is.null(activate_r)) return(TRUE)
      isTRUE(activate_r())
    })
    
    
    active_integrator_dir <- shiny::reactive({
      if (is.null(integrator_dir_r)) {
        return("")
      }
      
      idir <- tryCatch(integrator_dir_r(), error = function(e) "")
      idir <- as.character(idir %||% "")
      
      if (!nzchar(idir) || !dir.exists(idir)) {
        return("")
      }
      
      normalizePath(idir, winslash = "/", mustWork = FALSE)
    })
    
    observeEvent(active_integrator_dir(), {
      idir <- active_integrator_dir()
      if (!nzchar(idir)) return()
      
      mf <- read_manifest_safe_mod(idir)
      if (is.null(mf)) {
        append_log("[MANIFEST] missing/unreadable in ", idir)
      } else {
        append_log("[MANIFEST] session_id=", mf$session_id %||% "")
        append_log("[MANIFEST] gwas_session_file=", mf$gwas_session_file %||% "")
        append_log("[MANIFEST] cluster_method=", mf$cluster_method %||% "")
        append_log("[MANIFEST] threshold_used=", as.character(mf$threshold_used %||% ""))
      }
    }, ignoreInit = FALSE)
    
    workdir <- file.path(tempdir(), paste0("ldint_", session$token))
    dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
    
    
    output$ld_debug_clusters_dt <- DT::renderDT({
      cl <- integrated_clusters()
      if (is.null(cl)) {
        return(DT::datatable(data.frame(Message = "Integrated clusters not available"), options = list(dom = "t"), rownames = FALSE))
      }
      DT::datatable(cl, extensions = "Buttons",
                    options = list(dom = "Bfrtip", 
                                   buttons = list(
                                     list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
                                     list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
                                     list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
                                   ), 
                                   pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    }, server = FALSE)
    
    
    output$ld_debug_candidates_dt <- DT::renderDT({
      ca <- integrated_candidates()
      if (is.null(ca)) {
        return(DT::datatable(data.frame(Message = "Integrated candidates not available"), options = list(dom = "t"), rownames = FALSE))
      }
      DT::datatable(ca, extensions = "Buttons",
                    options = list(dom = "Bfrtip", buttons = list(
                      list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
                      list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
                      list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
                    ),
                    pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    }, server = FALSE)
    
    clusters_available <- shiny::reactive({
      cl <- clusters_base_r()
      ca <- candidates_base_r()
      
      if (!isTRUE(is_active())) return(NULL)
      if (!is.data.frame(cl) || !nrow(cl)) return(NULL)
      if (!is.data.frame(ca) || !nrow(ca)) return(NULL)
      
      cls_sum <- ca %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          has_gwas = any(classe == "GWAS", na.rm = TRUE),
          classes  = paste(sort(unique(classe)), collapse = "; "),
          .groups = "drop"
        )
      
      keep_ids <- cls_sum %>%
        dplyr::filter(has_gwas) %>%
        dplyr::pull(cluster_id) %>%
        unique()
      
      out <- cl %>%
        dplyr::filter(cluster_id %in% keep_ids) %>%
        dplyr::arrange(chr, start, end)
      
      append_log("[AVAILABLE] integrated clusters with GWAS n=", nrow(out))
      out
    })
    
    # chromosome selector
    observe({
      if (!isTRUE(is_active())) return()
      
      cl <- clusters_available()
      
      append_log("[OBS ld_chr] entered")
      append_log("[OBS ld_chr] is.data.frame(cl)=", is.data.frame(cl))
      append_log("[OBS ld_chr] nrow(cl)=", if (is.data.frame(cl)) nrow(cl) else NA)
      
      if (!is.data.frame(cl) || !nrow(cl)) {
        shiny::updateSelectInput(
          session,
          "ld_chr",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      chr_vals <- suppressWarnings(as.integer(cl$chr))
      chr_vals <- chr_vals[is.finite(chr_vals)]
      chr_vals <- sort(unique(chr_vals))
      
      append_log("[OBS ld_chr] chr_vals=", paste(chr_vals, collapse = ", "))
      
      if (!length(chr_vals)) {
        shiny::updateSelectInput(
          session,
          "ld_chr",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      chr_choices <- as.character(chr_vals)
      chr_selected <- chr_choices[1]
      
      append_log("[OBS ld_chr] updating choices=", paste(chr_choices, collapse = ", "))
      append_log("[OBS ld_chr] selected=", chr_selected)
      
      shiny::updateSelectInput(
        session,
        "ld_chr",
        choices = chr_choices,
        selected = chr_selected
      )
    })
    
    # cluster selector
    observe({
      if (!isTRUE(is_active())) return()
      
      cl <- clusters_available()
      
      if (!is.data.frame(cl) || !nrow(cl)) {
        shiny::updateSelectInput(
          session,
          "ld_cluster_id",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      chr_sel <- suppressWarnings(as.integer(input$ld_chr))
      
      chr_available <- suppressWarnings(as.integer(cl$chr))
      chr_available <- chr_available[is.finite(chr_available)]
      chr_available <- sort(unique(chr_available))
      
      if (!length(chr_available)) {
        shiny::updateSelectInput(
          session,
          "ld_cluster_id",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      if (!is.finite(chr_sel) || !chr_sel %in% chr_available) {
        chr_sel <- chr_available[1]
      }
      
      cl_chr <- cl %>%
        dplyr::filter(chr == chr_sel) %>%
        dplyr::arrange(start, end)
      
      if (!nrow(cl_chr)) {
        shiny::updateSelectInput(
          session,
          "ld_cluster_id",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      if (!"cluster_id" %in% names(cl_chr)) cl_chr$cluster_id <- NA_character_
      if (!"cluster_key" %in% names(cl_chr)) cl_chr$cluster_key <- NA_character_
      if (!"source_apps" %in% names(cl_chr)) cl_chr$source_apps <- ""
      
      cl_chr <- cl_chr %>%
        dplyr::mutate(
          cluster_id = as.character(cluster_id),
          cluster_key = as.character(cluster_key),
          source_apps = as.character(source_apps),
          cluster_key = dplyr::if_else(
            !is.na(cluster_key) & nzchar(cluster_key),
            cluster_key,
            paste0(
              cluster_id, "|chr", chr_label_plink(chr), ":",
              start, "-", end
            )
          )
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(cluster_key), nzchar(cluster_key)
        )
      
      if (!nrow(cl_chr)) {
        shiny::updateSelectInput(
          session,
          "ld_cluster_id",
          choices = character(0),
          selected = character(0)
        )
        return()
      }
      
      labs <- paste0(
        cl_chr$cluster_id,
        " (chr", chr_label_plink(cl_chr$chr), ":", cl_chr$start, "-", cl_chr$end, ")",
        " [", cl_chr$source_apps, "]"
      )
      
      vals <- cl_chr$cluster_key
      
      choices <- as.list(vals)
      names(choices) <- labs
      
      shiny::updateSelectInput(
        session,
        "ld_cluster_id",
        choices = choices,
        selected = vals[1]
      )
    })
    
    
    selected_cluster <- shiny::reactive({
      cl <- clusters_available()
      
      if (!is.data.frame(cl) || !nrow(cl)) {
        return(NULL)
      }
      
      if (!"cluster_key" %in% names(cl)) {
        cl <- cl %>%
          dplyr::mutate(
            cluster_key = paste0(
              as.character(cluster_id),
              "|chr", chr_label_plink(chr), ":",
              start, "-", end
            )
          )
      }
      
      sel <- as.character(input$ld_cluster_id %||% "")
      if (!nzchar(sel)) {
        return(NULL)
      }
      
      one <- cl %>%
        dplyr::filter(cluster_key == sel)
      
      if (!nrow(one)) {
        return(NULL)
      }
      
      one[1, , drop = FALSE]
    })
    
    # ---------------------
    selected_cluster_id_r <- shiny::reactive({
      cl <- selected_cluster()
      if (!is.data.frame(cl) || !nrow(cl) || !"cluster_id" %in% names(cl)) return("")
      as.character(cl$cluster_id[[1]] %||% "")
    })
    
    cluster_details_rds_path_r <- shiny::reactive({
      idir <- active_integrator_dir()
      cid <- selected_cluster_id_r()
      if (!nzchar(idir) || !dir.exists(idir) || !nzchar(cid)) return("")
      file.path(idir, paste0("ld_cluster_", cid, "_details.rds"))
    })
    
    cluster_summary_rds_path_r <- shiny::reactive({
      idir <- active_integrator_dir()
      cid <- selected_cluster_id_r()
      if (!nzchar(idir) || !dir.exists(idir) || !nzchar(cid)) return("")
      file.path(idir, paste0("ld_cluster_", cid, "_summary.rds"))
    })
    
    cluster_details_saved_r <- shiny::reactive({
      fp <- cluster_details_rds_path_r()
      if (!nzchar(fp) || !file.exists(fp)) return(NULL)
      obj <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (!is.list(obj)) return(NULL)
      obj
    })
    
    cluster_summary_saved_r <- shiny::reactive({
      fp <- cluster_summary_rds_path_r()
      if (!nzchar(fp) || !file.exists(fp)) return(NULL)
      obj <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (!is.data.frame(obj) || !nrow(obj)) return(NULL)
      obj
    })
    # ----------------------
    
    
    
    selected_candidates <- shiny::reactive({
      cl <- selected_cluster()
      ca <- candidates_base_r()
      
      if (!is.data.frame(cl) || !nrow(cl)) return(NULL)
      if (!is.data.frame(ca) || !nrow(ca)) return(NULL)
      
      chr_sel <- as.integer(cl$chr[1])
      st <- as.integer(cl$start[1])
      en <- as.integer(cl$end[1])
      cid <- as.character(cl$cluster_id[1])
      
      ca_sub <- ca %>%
        dplyr::filter(
          chr == chr_sel,
          position >= st,
          position <= en
        ) %>%
        dplyr::mutate(cluster_id = normalize_cluster_id(cluster_id))
      
      ca_sub2 <- ca_sub %>% dplyr::filter(cluster_id == cid)
      if (nrow(ca_sub2) > 0) ca_sub <- ca_sub2
      
      ca_sub <- ca_sub %>%
        dplyr::mutate(classe = normalize_classe(classe)) %>%
        dplyr::distinct(source_app, cluster_id, id_hit, rsid, chr, position, classe, .keep_all = TRUE) %>%
        dplyr::arrange(position, classe, id_hit)
      
      append_log("[SELECTED_CANDIDATES] cluster_id=", cid, " | n=", nrow(ca_sub))
      ca_sub
    })
    
    selected_gwas_hits <- shiny::reactive({
      ca <- selected_candidates()
      if (!is.data.frame(ca) || !nrow(ca)) return(NULL)
      
      out <- ca %>%
        dplyr::filter(classe == "GWAS") %>%
        dplyr::distinct(cluster_id, rsid, position, .keep_all = TRUE) %>%
        dplyr::arrange(position, rsid)
      
      append_log("[SELECTED_GWAS] n=", nrow(out))
      out
    })
    
    output$ld_pop_ui <- shiny::renderUI({
      dd <- input$ld_pops_dir %||% ""
      d  <- tryCatch(normalizePath(dd, winslash = "/", mustWork = FALSE), error = function(e) dd)
      
      ff <- if (nzchar(d) && dir.exists(d)) list.files(d, pattern="\\.[Tt][Xx][Tt]$", full.names = FALSE) else character(0)
      pops <- sub("\\.[Tt][Xx][Tt]$", "", ff)
      
      if (length(pops) == 0) {
        return(shiny::tags$div(
          style="color:#b00020; font-weight:600;",
          paste0("No .txt keep files found in: ", dd)
        ))
      }
      
      sel <- if (!is.null(input$ld_pop) && input$ld_pop %in% pops) input$ld_pop else {
        if ("EUR" %in% pops) "EUR" else pops[1]
      }
      
      append_log("[LD-POP-UI] pops=", paste(sort(pops), collapse = ", "))
      append_log("[LD-POP-UI] input$ld_pop=", input$ld_pop %||% "NULL")
      append_log("[LD-POP-UI] sel=", sel)
      
      shiny::selectInput(session$ns("ld_pop"), "Population (keep file)", choices = sort(pops), selected = sel)
    })
    
    ld_state <- shiny::reactiveValues(
      fl = NULL,
      ld_pairs = NULL,
      blocks_ij = NULL,
      subset_prefix = NULL,
      subset_bim = NULL,
      tag = NULL,
      chr_sel = NULL,
      st = NULL,
      en = NULL,
      seeds = NULL,
      proxy_table = NULL,
      cand_tracks = NULL,
      block_summary = NULL,
      block_genes = NULL
    )
    
    get_ld_pop_selected <- function() {
      dd <- input$ld_pops_dir %||% ""
      d  <- tryCatch(normalizePath(dd, winslash = "/", mustWork = FALSE), error = function(e) dd)
      
      ff <- if (nzchar(d) && dir.exists(d)) {
        list.files(d, pattern = "\\.[Tt][Xx][Tt]$", full.names = FALSE)
      } else {
        character(0)
      }
      
      pops <- sort(sub("\\.[Tt][Xx][Tt]$", "", ff))
      if (!length(pops)) return("")
      
      cur <- trimws(as.character(input$ld_pop %||% ""))
      if (nzchar(cur) && cur %in% pops) return(cur)
      
      if ("EUR" %in% pops) return("EUR")
      pops[1]
    }
    
    ##########################################################################
    ##################### LD PLOT SELECTION FROM LIVE or RDS #################
    ##########################################################################
    ld_cluster_details_for_display <- reactive({
      obj <- cluster_details_saved_r()
      
      required_fields <- c(
        "cluster_id", "fl_plot", "ld_pairs", "blocks_ij", "cand_tracks",
        "chr", "start", "end", "ld_metric"
      )
      
      if (is.list(obj) && all(required_fields %in% names(obj))) {
        return(obj)
      }
      
      NULL
    })
    
    ld_plot_payload <- reactive({
      saved <- ld_cluster_details_for_display()
      
      if (is.list(saved)) {
        return(list(
          cluster_id = saved$cluster_id %||% NA_character_,
          chr = saved$chr %||% NA_character_,
          start = saved$start %||% NA_real_,
          end = saved$end %||% NA_real_,
          ld_metric = saved$ld_metric %||% "R2",
          fl_plot = saved$fl_plot,
          ld_pairs = saved$ld_pairs,
          blocks_ij = saved$blocks_ij,
          cand_tracks = saved$cand_tracks,
          source = "saved"
        ))
      }
      
      if (!is.null(ld_state$fl) && !is.null(ld_state$ld_pairs)) {
        return(list(
          cluster_id = ld_state$cluster_id %||% NA_character_,
          chr = ld_state$chr_sel %||% NA_character_,
          start = ld_state$st %||% NA_real_,
          end = ld_state$en %||% NA_real_,
          ld_metric = input$ld_metric %||% "R2",
          fl_plot = ld_state$fl,
          ld_pairs = ld_state$ld_pairs,
          blocks_ij = ld_state$blocks_ij,
          cand_tracks = ld_state$cand_tracks,
          source = "live"
        ))
      }
      
      NULL
    })
    
    ld_pairs_for_summary_r <- shiny::reactive({
      saved <- ld_cluster_details_for_display()
      
      if (is.list(saved) &&
          is.data.frame(saved$ld_pairs) &&
          nrow(saved$ld_pairs) > 0) {
        return(saved$ld_pairs)
      }
      
      if (is.data.frame(ld_state$ld_pairs) && nrow(ld_state$ld_pairs) > 0) {
        return(ld_state$ld_pairs)
      }
      
      tibble::tibble()
    })
    
    ld_plot_is_stale <- reactive({
      saved <- ld_cluster_details_for_display()
      if (is.list(saved)) return(FALSE)
      
      cl <- selected_cluster()
      sel_id <- if (is.data.frame(cl) && nrow(cl) == 1 && "cluster_id" %in% names(cl)) {
        as.character(cl$cluster_id[[1]])
      } else {
        NA_character_
      }
      
      calc_id <- ld_state$cluster_id %||% NA_character_
      
      isTRUE(!is.na(sel_id) && !is.na(calc_id) && nzchar(sel_id) && nzchar(calc_id) && sel_id != calc_id)
    })
    
    ##########################################################################
 
    observeEvent(input$run_cluster_ld, {
      if (!isTRUE(is_active())) return()
  
      
      ld_log("")
      append_log("[LD-ALL] Starting...")
      
      tryCatch({
        withProgress(message = "Running LD, blocks and prioritization...", value = 0, {
          
          # --------------------------------------------------
          # 0) reset state
          # --------------------------------------------------
          ld_state$fl <- NULL
          ld_state$ld_pairs <- NULL
          ld_state$blocks_ij <- NULL
          ld_state$subset_prefix <- NULL
          ld_state$subset_bim <- NULL
          ld_state$tag <- NULL
          ld_state$chr_sel <- NULL
          ld_state$st <- NULL
          ld_state$en <- NULL
          ld_state$cand_tracks <- NULL
          ld_state$seeds <- NULL
          ld_state$proxy_table <- NULL
          ld_state$block_summary <- NULL
          ld_state$block_genes <- NULL
          
          # --------------------------------------------------
          # 1) inputs / checks
          # --------------------------------------------------
          incProgress(0.10, detail = "Preparing inputs")
          
          cl_row <- selected_cluster()
          ca_all <- integrated_candidates_r()
          
          shiny::validate(shiny::need(
            is.data.frame(cl_row) && nrow(cl_row) > 0,
            "No selected cluster available."
          ))
          shiny::validate(shiny::need(
            is.data.frame(ca_all) && nrow(ca_all) > 0,
            "Integrated candidates not available."
          ))
          
          plink_bin <- input$plink_bin %||% ""
          bfile_ref <- input$bfile_ref_prefix %||% ""
          pop_sel   <- get_ld_pop_selected()
          keep_dir  <- input$ld_pops_dir %||% ""
          keep_path <- read_keep_file_dir(pop_sel, keep_dir)
          
          maxn <- suppressWarnings(as.integer(input$ld_max_snps_interval %||% 400))
          shiny::validate(shiny::need(is.finite(maxn) && maxn >= 2, "Invalid max SNPs limit."))
          
          r2_min <- suppressWarnings(as.numeric(input$ld_proxy_r2_min %||% 0.6))
          ld_metric_now <- input$ld_metric %||% "R2"
          
          # --------------------------------------------------
          # 2) common canonical bundle
          # --------------------------------------------------
          incProgress(0.45, detail = "Computing canonical LD bundle")
          
          cat("[DBG RUN] before compute_ld_bundle_common\n")
          
          bundle <- compute_ld_bundle_common(
            cluster_row = cl_row,
            candidates_df = ca_all,
            bfile_ref = bfile_ref,
            keep_path = keep_path,
            plink_bin = plink_bin,
            workdir = workdir,
            pop = pop_sel,
            ld_metric = ld_metric_now,
            r2_min = r2_min,
            max_snps_interval = maxn,
            compute_blocks = TRUE,
            append_log = append_log
          )
          
          cat("[DBG RUN] after compute_ld_bundle_common\n")
          
          fl <- bundle$fl
          bim_df <- bundle$bim
          chr_sel <- bundle$chr
          st <- bundle$start
          en <- bundle$end
          cid <- bundle$cluster_id
          subset_prefix <- bundle$subset_prefix
          
          cand_sub <- bundle$candidates
          hits_sub <- bundle$gwas_hits
          seeds_df <- bundle$seeds
          proxy_tbl <- bundle$proxies
          ij <- bundle$blocks_ij
          
          append_log("[LD-ALL] cluster_id=", cid, " | query_hits=", nrow(hits_sub))
          
          shiny::validate(shiny::need(
            is.data.frame(hits_sub) && nrow(hits_sub) > 0,
            "No GWAS hits available in selected integrated cluster."
          ))
          
          # --------------------------------------------------
          # 3) full LD matrix for plot/debug
          # --------------------------------------------------
          incProgress(0.70, detail = "Computing full LD matrix")
          
          ld_prefix <- file.path(workdir, paste0("ldint_", cid, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_ld"))
          r2_args <- c(
            "--bfile", subset_prefix,
            "--ld-window", "999999",
            "--ld-window-kb", "999999",
            "--ld-window-r2", "0"
          )
          
          if (identical(ld_metric_now, "Dprime")) {
            r2_args <- c(r2_args, "--r2", "dprime", "gz")
          } else {
            r2_args <- c(r2_args, "--r2", "gz")
          }
          
          r2_args <- c(r2_args, "--out", ld_prefix)
          r2 <- run_plink(r2_args, ld_prefix, plink_bin)
          if (length(r2$stdout)) append_log(paste(r2$stdout, collapse = "\n"))
          if (is.null(r2$status) || r2$status != 0) stop("PLINK LD failed.")
          
          ld_raw <- read_plink_ld(ld_prefix)
          snpA <- pick_col(ld_raw, c("SNP_A","SNP_A1","SNP1"))
          snpB <- pick_col(ld_raw, c("SNP_B","SNP_B1","SNP2"))
          vcol <- if (identical(ld_metric_now, "Dprime")) {
            pick_col(ld_raw, c("Dprime","D'","DP","DPRIME"))
          } else {
            pick_col(ld_raw, c("R2","r2"))
          }
          
          shiny::validate(shiny::need(
            !is.null(snpA) && !is.null(snpB) && !is.null(vcol),
            "LD file missing required columns."
          ))
          
          fl_map <- fl %>% dplyr::select(SNP, BP, ix)
          
          ld2 <- ld_raw %>%
            dplyr::transmute(
              SNP_A = as.character(.data[[snpA]]),
              SNP_B = as.character(.data[[snpB]]),
              value = suppressWarnings(as.numeric(.data[[vcol]]))
            ) %>%
            dplyr::filter(!is.na(SNP_A), !is.na(SNP_B), is.finite(value)) %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixA = ix, BPA = BP), by = c("SNP_A" = "SNP")) %>%
            dplyr::left_join(fl_map %>% dplyr::rename(ixB = ix, BPB = BP), by = c("SNP_B" = "SNP")) %>%
            dplyr::filter(is.finite(ixA), is.finite(ixB), ixA != ixB) %>%
            dplyr::mutate(i = pmax(ixA, ixB), j = pmin(ixA, ixB)) %>%
            dplyr::distinct(i, j, .keep_all = TRUE) %>%
            dplyr::filter(i > j)
          
          # --------------------------------------------------
          # 4) candidate tracks (same candidates used in canonical bundle)
          # --------------------------------------------------
          incProgress(0.82, detail = "Building candidate tracks")
          
          cand_tracks <- tibble::tibble(
            BP = integer(),
            classe = character(),
            label_id = character(),
            source_app = character()
          )
          
          if (is.data.frame(cand_sub) && nrow(cand_sub) > 0) {
            cand_tracks <- cand_sub %>%
              dplyr::transmute(
                BP = suppressWarnings(as.integer(position)),
                classe = normalize_classe(classe),
                label_id = dplyr::coalesce(
                  as.character(rsid),
                  as.character(id_hit),
                  as.character(cluster_id)
                ),
                source_app = as.character(source_app)
              ) %>%
              dplyr::mutate(
                label_id = trimws(label_id),
                label_id = dplyr::if_else(!nzchar(label_id), paste0("pos_", BP), label_id)
              ) %>%
              dplyr::filter(is.finite(BP), !is.na(classe), nzchar(classe)) %>%
              dplyr::distinct(BP, classe, label_id, source_app, .keep_all = TRUE) %>%
              dplyr::arrange(BP, classe, label_id)
          }
          
          fl2 <- fl %>%
            dplyr::mutate(is_candidate = FALSE, classe = NA_character_) %>%
            dplyr::select(SNP, BP, ix, is_candidate, classe)
          
          if (nrow(cand_tracks) > 0) {
            fl2 <- fl2 %>%
              dplyr::mutate(
                is_candidate = BP %in% cand_tracks$BP
              ) %>%
              dplyr::left_join(
                cand_tracks %>%
                  dplyr::group_by(BP) %>%
                  dplyr::summarise(
                    classe_track = paste(sort(unique(classe)), collapse = "; "),
                    .groups = "drop"
                  ),
                by = "BP"
              ) %>%
              dplyr::mutate(
                classe = dplyr::coalesce(classe_track, classe)
              ) %>%
              dplyr::select(SNP, BP, ix, is_candidate, classe)
          }
          
          
          # --------------------------------------------------
          # 4b) canonical block summary for table + persistence
          # --------------------------------------------------
          session_dir <- active_integrator_dir()
          
          gwas_bridge_save <- gwas_bridge_base_r()
          
          canonical_block_summary <- block_overlap_summary_df(
            block_ranges = bundle$block_ranges,
            block_hits = bundle$block_hits,
            block_genes = bundle$block_genes,
            gwas_bridge_df = gwas_bridge_save,
            proxy_tbl = bundle$proxies
          )
          
          if (!is.data.frame(canonical_block_summary)) {
            canonical_block_summary <- tibble::tibble()
          }
          
          cat("[DBG RUN] canonical_block_summary nrow = ", nrow(canonical_block_summary), "\n", sep = "")
          if (nrow(canonical_block_summary)) {
            cat("[DBG RUN] canonical_block_summary cols = ",
                paste(names(canonical_block_summary), collapse = ", "), "\n", sep = "")
          }
          # --------------------------------------------------
          # 5) save state from canonical bundle
          # --------------------------------------------------
          incProgress(0.95, detail = "Saving state")
          
          ld_state$cand_tracks <- cand_tracks
          ld_state$fl <- fl2
          ld_state$ld_pairs <- ld2
          ld_state$blocks_ij <- ij
          ld_state$subset_prefix <- subset_prefix
          ld_state$subset_bim <- bim_df
          ld_state$tag <- paste0("ldint_", cid)
          ld_state$cluster_id <- cid
          ld_state$chr_sel <- chr_sel
          ld_state$st <- st
          ld_state$en <- en
          ld_state$seeds <- seeds_df
          ld_state$proxy_table <- proxy_tbl
          ld_state$block_summary <- canonical_block_summary
          ld_state$block_genes <- bundle$block_genes
          
          # --------------------------------------------------
          # 6) save cluster-specific RDS
          # --------------------------------------------------
          if (nzchar(session_dir) && dir.exists(session_dir)) {
            
            cid_safe <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(cid))
            
            summary_path_one <- file.path(
              session_dir,
              paste0("ld_cluster_", cid_safe, "_summary.rds")
            )
            
            details_path_one <- file.path(
              session_dir,
              paste0("ld_cluster_", cid_safe, "_details.rds")
            )
            
            summary_one <- data.frame(
              cluster_id = bundle$cluster_id,
              chr = as.character(chr_sel),
              start = as.integer(st),
              end = as.integer(en),
              population = as.character(pop_sel),
              ld_metric = as.character(ld_metric_now),
              n_interval_snps = nrow(bundle$fl),
              n_candidate_snps = nrow(bundle$candidates),
              n_gwas_hits = nrow(bundle$gwas_hits),
              gwas_hits = if (is.data.frame(bundle$gwas_hits) && nrow(bundle$gwas_hits) && "rsid" %in% names(bundle$gwas_hits)) {
                paste(unique(stats::na.omit(bundle$gwas_hits$rsid)), collapse = "; ")
              } else {
                ""
              },
              n_seeds_exact = if (is.data.frame(bundle$seeds) && nrow(bundle$seeds) && "seed_type" %in% names(bundle$seeds)) {
                sum(bundle$seeds$seed_type == "exact", na.rm = TRUE)
              } else {
                0L
              },
              n_seeds_nearest = if (is.data.frame(bundle$seeds) && nrow(bundle$seeds) && "seed_type" %in% names(bundle$seeds)) {
                sum(bundle$seeds$seed_type == "nearest", na.rm = TRUE)
              } else {
                0L
              },
              seed_snps = if (is.data.frame(bundle$seeds) && nrow(bundle$seeds) && "seed_snp" %in% names(bundle$seeds)) {
                paste(unique(stats::na.omit(bundle$seeds$seed_snp)), collapse = "; ")
              } else {
                ""
              },
              n_proxy_snps = if (is.data.frame(bundle$proxies) && nrow(bundle$proxies) && "proxy_snp" %in% names(bundle$proxies)) {
                dplyr::n_distinct(bundle$proxies$proxy_snp)
              } else {
                0L
              },
              proxy_snps = if (is.data.frame(bundle$proxies) && nrow(bundle$proxies) && "proxy_snp" %in% names(bundle$proxies)) {
                paste(unique(stats::na.omit(bundle$proxies$proxy_snp)), collapse = "; ")
              } else {
                ""
              },
              max_ld_value = if (is.data.frame(bundle$proxies) && nrow(bundle$proxies) && "ld_value" %in% names(bundle$proxies)) {
                max(bundle$proxies$ld_value, na.rm = TRUE)
              } else {
                NA_real_
              },
              mean_ld_value = if (is.data.frame(bundle$proxies) && nrow(bundle$proxies) && "ld_value" %in% names(bundle$proxies)) {
                mean(bundle$proxies$ld_value, na.rm = TRUE)
              } else {
                NA_real_
              },
              n_blocks = if (is.data.frame(bundle$block_ranges)) nrow(bundle$block_ranges) else 0L,
              n_block_hits = if (is.data.frame(canonical_block_summary) && nrow(canonical_block_summary) && "block_support_hits" %in% names(canonical_block_summary)) {
                sum(canonical_block_summary$block_support_hits, na.rm = TRUE)
              } else if (is.data.frame(canonical_block_summary) && nrow(canonical_block_summary) && "n_block_hits" %in% names(canonical_block_summary)) {
                sum(canonical_block_summary$n_block_hits, na.rm = TRUE)
              } else {
                0L
              },
              
              n_block_genes = if (is.data.frame(canonical_block_summary) && nrow(canonical_block_summary) && "n_genes_in_block" %in% names(canonical_block_summary)) {
                sum(canonical_block_summary$n_genes_in_block, na.rm = TRUE)
              } else if (is.data.frame(canonical_block_summary) && nrow(canonical_block_summary) && "n_block_genes" %in% names(canonical_block_summary)) {
                sum(canonical_block_summary$n_block_genes, na.rm = TRUE)
              } else {
                0L
              },
              ld_has_signal = as.integer(is.data.frame(bundle$proxies) && nrow(bundle$proxies) > 0),
              ld_computed = TRUE,
              status = "ok",
              stringsAsFactors = FALSE
            )
            
            details_one <- list(
              # identificació bàsica
              cluster_id = bundle$cluster_id,
              chr = as.character(chr_sel),
              start = as.integer(st),
              end = as.integer(en),
              population = as.character(pop_sel),
              ld_metric = as.character(ld_metric_now),
              x_mode = as.character(input$ld_x_mode %||% "bp"),
              
              # objectes canònics
              fl = bundle$fl,
              bim = bundle$bim,
              candidates = bundle$candidates,
              gwas_hits = bundle$gwas_hits,
              seeds = bundle$seeds,
              proxies = bundle$proxies,
              blocks = bundle$blocks,
              blocks_ij = bundle$blocks_ij,
              block_ranges = bundle$block_ranges,
              block_hits = bundle$block_hits,
              block_genes = bundle$block_genes,
              block_summary = NULL,
              
              # objectes específics del plot
              fl_plot = fl2,
              ld_pairs = ld2,
              cand_tracks = cand_tracks,
              
              # metadades de render
              plot_meta = list(
                x_mode = as.character(input$ld_x_mode %||% "bp"),
                ld_metric = as.character(ld_metric_now),
                max_snps_interval = as.integer(maxn),
                ld_proxy_r2_min = as.numeric(r2_min),
                saved_at = as.character(Sys.time())
              ),
              
              # resum final
              summary = summary_one
            )
            
            summary_tbl <- canonical_block_summary
            
            if (!is.data.frame(summary_tbl)) {
              summary_tbl <- data.frame()
            }
            
            cat("[LD-SAVE SUMMARY] class(summary_tbl) = ", paste(class(summary_tbl), collapse = ", "), "\n", sep = "")
            cat("[LD-SAVE SUMMARY] is.data.frame(summary_tbl) = ", is.data.frame(summary_tbl), "\n", sep = "")
            cat("[LD-SAVE SUMMARY] nrow(summary_tbl) = ", if (is.data.frame(summary_tbl)) nrow(summary_tbl) else NA, "\n", sep = "")
            cat("[LD-SAVE SUMMARY] ncol(summary_tbl) = ", if (is.data.frame(summary_tbl)) ncol(summary_tbl) else NA, "\n", sep = "")
            if (is.data.frame(summary_tbl)) {
              cat("[LD-SAVE SUMMARY] names = ", paste(names(summary_tbl), collapse = ", "), "\n", sep = "")
            }
            
            details_one$block_summary <- summary_tbl
            
            saveRDS(summary_tbl, summary_path_one)
            saveRDS(details_one, details_path_one)
            
            append_log("[LD-ALL] saved summary table: ", summary_path_one)
            append_log("[LD-ALL] saved details: ", details_path_one)
          }
          
          incProgress(1, detail = "Done")
          append_log("[LD-ALL] Done")
        })
      }, error = function(e) {
        append_log("[LD-ALL][ERROR] ", conditionMessage(e))
        showNotification(
          paste("LD run error:", conditionMessage(e)),
          type = "error",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)
  
    
    
    output$ld_fl_dt <- DT::renderDT({
      fl  <- ld_state$fl
      blk <- block_ranges_selected()
      
      # Cas 1: encara no s'ha fet RUN LD/blocks
      if (!is.data.frame(blk) || !nrow(blk)) {
        return(DT::datatable(
          data.frame(Message = "Run LD & blocks to populate SNP list."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      # Cas 2: sí que hi ha LD/blocks, però aquest cluster no té SNP list
      if (is.null(fl) || !is.data.frame(fl) || !nrow(fl)) {
        return(DT::datatable(
          data.frame(Message = "No SNP list available for this cluster."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      show <- fl %>%
        dplyr::mutate(
          candidate = dplyr::if_else(
            isTRUE(is_candidate) | (is.logical(is_candidate) & is_candidate),
            "YES",
            ""
          ),
          classe = dplyr::coalesce(as.character(classe), "")
        ) %>%
        dplyr::transmute(ix, SNP, BP, candidate, classe)
      
      DT::datatable(
        show,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          pageLength = 15,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_pairs_dt <- DT::renderDT({
      lp  <- ld_state$ld_pairs
      fl  <- ld_state$fl
      blk <- block_ranges_selected()
      
      # Cas 1: encara no s'ha fet RUN LD/blocks
      if (!is.data.frame(blk) || !nrow(blk)) {
        return(DT::datatable(
          data.frame(Message = "Run LD & blocks to populate pairs."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      # Cas 2: sí que hi ha LD/blocks, però aquest cluster no té parelles LD
      if (is.null(lp) || !is.data.frame(lp) || !nrow(lp) ||
          is.null(fl) || !is.data.frame(fl) || !nrow(fl)) {
        return(DT::datatable(
          data.frame(Message = "No LD pairs found for this cluster."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      DT::datatable(
        lp %>% dplyr::arrange(dplyr::desc(value)),
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          pageLength = 12,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_seed_dt <- DT::renderDT({
      dd  <- ld_state$seeds
      blk <- block_ranges_selected()
      
      # Cas 1: encara no s'ha fet RUN LD/blocks
      if (!is.data.frame(blk) || !nrow(blk)) {
        return(DT::datatable(
          data.frame(Message = "Run LD & blocks to populate seeds."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      # Cas 2: sí que hi ha LD/blocks, però aquest cluster no té seeds
      if (is.null(dd) || !is.data.frame(dd) || !nrow(dd)) {
        return(DT::datatable(
          data.frame(Message = "No seed SNPs found for this cluster."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      DT::datatable(
        dd,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          pageLength = 12,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$ld_proxy_dt <- DT::renderDT({
      dd  <- ld_state$proxy_table
      blk <- block_ranges_selected()
      
      # Cas 1: encara no s'ha fet RUN LD/blocks
      if (!is.data.frame(blk) || !nrow(blk)) {
        return(DT::datatable(
          data.frame(Message = "Run LD & blocks to populate proxy table."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      # Cas 2: sí que hi ha LD/blocks, però aquest cluster no té proxies
      if (is.null(dd) || !is.data.frame(dd) || !nrow(dd)) {
        return(DT::datatable(
          data.frame(Message = "No proxy SNPs found for this cluster."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      DT::datatable(
        dd,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          pageLength = 15,
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    
    make_gene_track_plot <- function(chr_sel, st, en, fl, x_limits2, x_mode, append_log = NULL) {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomicRanges", quietly = TRUE)) return(NULL)
      if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) return(NULL)
      if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(NULL)
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
      if (!requireNamespace("S4Vectors", quietly = TRUE)) return(NULL)
      if (!requireNamespace("IRanges", quietly = TRUE)) return(NULL)
      if (!requireNamespace("dplyr", quietly = TRUE)) return(NULL)
      if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
      
      if (is.function(append_log)) {
        append_log("[GENE] building collapsed gene track | chr=", chr_sel, " start=", st, " end=", en)
      }
      
      chr_ucsc <- paste0("chr", chr_label_plink(chr_sel))
      
      region_gr <- GenomicRanges::GRanges(
        seqnames = chr_ucsc,
        ranges   = IRanges::IRanges(start = st, end = en)
      )
      
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      
      genes_gr <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
      GenomeInfoDb::seqlevelsStyle(genes_gr) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
      
      ov <- GenomicRanges::findOverlaps(genes_gr, region_gr, ignore.strand = TRUE)
      keep_idx <- unique(S4Vectors::queryHits(ov))
      
      if (!length(keep_idx)) {
        if (is.function(append_log)) append_log("[GENE] no genes overlapping region")
        return(NULL)
      }
      
      gsub_gr <- genes_gr[keep_idx]
      
      gene_id <- S4Vectors::mcols(gsub_gr)$gene_id
      if (is.null(gene_id)) gene_id <- names(gsub_gr)
      gene_id <- as.character(gene_id)
      
      gene_df <- data.frame(
        gene_id = gene_id,
        gene_start = GenomicRanges::start(gsub_gr),
        gene_end   = GenomicRanges::end(gsub_gr),
        stringsAsFactors = FALSE
      )
      
      mp <- tryCatch(
        suppressMessages(
          AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys    = unique(gene_df$gene_id),
            keytype = "ENTREZID",
            columns = c("SYMBOL")
          )
        ),
        error = function(e) NULL
      )
      
      if (is.data.frame(mp) && nrow(mp)) {
        mp <- mp %>%
          dplyr::mutate(gene_id = as.character(ENTREZID)) %>%
          dplyr::select(gene_id, SYMBOL) %>%
          dplyr::distinct(gene_id, .keep_all = TRUE)
        gene_df <- gene_df %>% dplyr::left_join(mp, by = "gene_id")
      } else {
        gene_df$SYMBOL <- gene_df$gene_id
      }
      
      gene_df <- gene_df %>%
        dplyr::mutate(label = dplyr::coalesce(SYMBOL, gene_id)) %>%
        dplyr::arrange(gene_start, gene_end) %>%
        dplyr::mutate(track_y = -seq_len(dplyr::n()))
      
      if (is.function(append_log)) {
        append_log("[GENE] genes overlapping region n=", nrow(gene_df))
      }
      
      ex_all <- suppressMessages(GenomicFeatures::exons(txdb, columns = c("gene_id")))
      GenomeInfoDb::seqlevelsStyle(ex_all) <- GenomeInfoDb::seqlevelsStyle(region_gr)[1]
      
      ov_ex <- GenomicRanges::findOverlaps(ex_all, region_gr, ignore.strand = TRUE)
      ex_sub <- ex_all[unique(S4Vectors::queryHits(ov_ex))]
      
      if (!length(ex_sub)) {
        if (is.function(append_log)) append_log("[GENE] no exons overlapping region")
        return(NULL)
      }
      
      ex_gene_id <- S4Vectors::mcols(ex_sub)$gene_id
      if (is.null(ex_gene_id)) {
        if (is.function(append_log)) append_log("[GENE] exon gene_id missing")
        return(NULL)
      }
      
      ex_df <- data.frame(
        gene_id = as.character(ex_gene_id),
        exon_start = GenomicRanges::start(ex_sub),
        exon_end   = GenomicRanges::end(ex_sub),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(gene_id %in% gene_df$gene_id)
      
      if (!nrow(ex_df)) {
        if (is.function(append_log)) append_log("[GENE] exon_df empty after gene filter")
        return(NULL)
      }
      
      collapse_exons_one_gene <- function(df_one) {
        df_one <- df_one[order(df_one$exon_start, df_one$exon_end), , drop = FALSE]
        
        starts <- df_one$exon_start
        ends   <- df_one$exon_end
        
        out_start <- integer(0)
        out_end   <- integer(0)
        
        cur_s <- starts[1]
        cur_e <- ends[1]
        
        if (length(starts) > 1) {
          for (k in 2:length(starts)) {
            if (starts[k] <= (cur_e + 1)) {
              cur_e <- max(cur_e, ends[k])
            } else {
              out_start <- c(out_start, cur_s)
              out_end   <- c(out_end, cur_e)
              cur_s <- starts[k]
              cur_e <- ends[k]
            }
          }
        }
        
        out_start <- c(out_start, cur_s)
        out_end   <- c(out_end, cur_e)
        
        data.frame(
          gene_id = df_one$gene_id[1],
          exon_start = out_start,
          exon_end = out_end,
          stringsAsFactors = FALSE
        )
      }
      
      ex_collapsed <- do.call(
        rbind,
        lapply(split(ex_df, ex_df$gene_id), collapse_exons_one_gene)
      )
      
      ex_collapsed <- ex_collapsed %>%
        dplyr::left_join(gene_df %>% dplyr::select(gene_id, label, track_y), by = "gene_id")
      
      bp_to_x <- function(bp_vec) {
        bp_vec <- suppressWarnings(as.numeric(bp_vec))
        if (!identical(x_mode, "equal")) return(bp_vec)
        
        stats::approx(
          x = fl$BP,
          y = fl$X,
          xout = bp_vec,
          rule = 2,
          ties = "ordered"
        )$y
      }
      
      gene_df <- gene_df %>%
        dplyr::mutate(
          x0 = bp_to_x(gene_start),
          x1 = bp_to_x(gene_end)
        )
      
      ex_collapsed <- ex_collapsed %>%
        dplyr::mutate(
          x0 = bp_to_x(exon_start),
          x1 = bp_to_x(exon_end)
        )
      
      p_gene <- plotly::plot_ly()
      
      for (i in seq_len(nrow(gene_df))) {
        p_gene <- p_gene %>%
          plotly::add_segments(
            x = gene_df$x0[i], xend = gene_df$x1[i],
            y = gene_df$track_y[i], yend = gene_df$track_y[i],
            inherit = FALSE,
            hoverinfo = "text",
            text = paste0(
              "Gene: ", gene_df$label[i],
              "<br>Start: ", gene_df$gene_start[i],
              "<br>End: ", gene_df$gene_end[i]
            ),
            showlegend = FALSE,
            line = list(width = 1)
          ) %>%
          plotly::add_text(
            x = x_limits2[2],
            y = gene_df$track_y[i],
            text = gene_df$label[i],
            textposition = "middle left",
            inherit = FALSE,
            hoverinfo = "skip",
            showlegend = FALSE
          )
      }
      
      exon_half_height <- 0.28
      
      for (i in seq_len(nrow(ex_collapsed))) {
        y0 <- ex_collapsed$track_y[i] - exon_half_height
        y1 <- ex_collapsed$track_y[i] + exon_half_height
        
        p_gene <- p_gene %>%
          plotly::add_trace(
            type = "scatter",
            mode = "lines",
            x = c(ex_collapsed$x0[i], ex_collapsed$x1[i], ex_collapsed$x1[i], ex_collapsed$x0[i], ex_collapsed$x0[i]),
            y = c(y0, y0, y1, y1, y0),
            fill = "toself",
            inherit = FALSE,
            hoverinfo = "text",
            text = paste0(
              "Gene: ", ex_collapsed$label[i],
              "<br>Exon: ", ex_collapsed$exon_start[i], "-", ex_collapsed$exon_end[i]
            ),
            showlegend = FALSE,
            line = list(width = 1)
          )
      }
      
      p_gene %>%
        plotly::layout(
          xaxis = list(
            showticklabels = FALSE,
            title = "",
            range = x_limits2
          ),
          yaxis = list(
            showticklabels = FALSE,
            title = ""
          ),
          margin = list(l = 40, r = 120, t = 10, b = 0)
        )
    }
    
    
    
    
    output$ld_plot <- plotly::renderPlotly({
      
      payload <- ld_plot_payload()
      
      validate(
        need(!is.null(payload), "No LD plot available yet. Run 'Compute LD / blocks'"),
        need(!isTRUE(ld_plot_is_stale()), "You selected a different cluster. Run 'Compute LD / blocks' to refresh the plot for the new cluster.")
      )
      
      fl_plot <- payload$fl_plot
      ld_pairs <- payload$ld_pairs
      blocks_ij <- payload$blocks_ij
      cand_tracks <- payload$cand_tracks
      cluster_lab <- payload$cluster_id %||% "LD cluster"
      ld_metric_plot <- payload$ld_metric %||% "R2"
      
      validate(
        need(is.data.frame(fl_plot) && nrow(fl_plot) > 0, "No SNP coordinates available."),
        need(is.data.frame(ld_pairs) && nrow(ld_pairs) > 0, "No LD pairs available.")
      )
      
      fl0   <- fl_plot
      ld    <- ld_pairs
      cand0 <- cand_tracks
      bd0   <- blocks_ij
      
      shiny::validate(
        shiny::need(is.data.frame(fl0) && nrow(fl0) >= 2, "Run LD plot first."),
        shiny::need(is.data.frame(ld) && nrow(ld) > 0, "No LD pairs available.")
      )
      
      metric <- input$ld_metric %||% "R2"
      x_mode <- input$ld_x_mode %||% "bp"
      
      append_log("[LD-PLOT][RENDER] entered")
      
      fl <- fl0 %>%
        dplyr::mutate(
          SNP = as.character(SNP),
          BP  = suppressWarnings(as.integer(BP)),
          is_candidate = dplyr::coalesce(as.logical(is_candidate), FALSE),
          classe = dplyr::coalesce(as.character(classe), "")
        ) %>%
        dplyr::distinct(SNP, .keep_all = TRUE) %>%
        dplyr::arrange(BP, SNP)
      
      if (identical(x_mode, "equal")) {
        fl <- fl %>% dplyr::mutate(X = dplyr::row_number())
        x_title <- "Equal spacing"
      } else {
        fl <- fl %>% dplyr::mutate(X = BP)
        x_title <- "Genomic (BP)"
      }
      
      xmap   <- stats::setNames(fl$X, fl$SNP)
      bpmap  <- stats::setNames(fl$BP, fl$SNP)
      idxmap <- stats::setNames(seq_len(nrow(fl)), fl$SNP)
      
      x_limits <- range(fl$X, na.rm = TRUE)
      x_span   <- max(1, diff(x_limits))
      
      bp_span <- diff(range(fl$BP, na.rm = TRUE))
      if (!is.finite(bp_span) || bp_span <= 0) bp_span <- 1
      
      tri <- ld %>%
        dplyr::transmute(
          SNP_A = as.character(SNP_A),
          SNP_B = as.character(SNP_B),
          val   = suppressWarnings(as.numeric(value))
        ) %>%
        dplyr::filter(is.finite(val), SNP_A %in% fl$SNP, SNP_B %in% fl$SNP) %>%
        dplyr::mutate(
          X_A = unname(xmap[SNP_A]),
          X_B = unname(xmap[SNP_B]),
          i   = unname(idxmap[SNP_A]),
          j   = unname(idxmap[SNP_B]),
          bpA = unname(bpmap[SNP_A]),
          bpB = unname(bpmap[SNP_B]),
          x_mid = (X_A + X_B) / 2,
          y = if (identical(x_mode, "equal")) {
            -(abs(j - i) / max(1, nrow(fl) - 1))
          } else {
            -(abs(bpB - bpA) / bp_span)
          },
          hover = paste0(
            "SNP_A: ", SNP_A,
            "<br>SNP_B: ", SNP_B,
            "<br>", metric, ": ", signif(val, 3),
            "<br>BP_A: ", bpA,
            "<br>BP_B: ", bpB
          )
        ) %>%
        dplyr::filter(is.finite(x_mid), is.finite(y))
      
      shiny::validate(shiny::need(nrow(tri) > 0, "No LD pairs to plot."))
      
      cand <- if (is.data.frame(cand0) && nrow(cand0)) {
        
        bp_to_x_track <- function(bp_vec, fl_df, x_mode_cur) {
          bp_vec <- suppressWarnings(as.numeric(bp_vec))
          
          if (!identical(x_mode_cur, "equal")) {
            return(bp_vec)
          }
          
          stats::approx(
            x = fl_df$BP,
            y = fl_df$X,
            xout = bp_vec,
            rule = 2,
            ties = "ordered"
          )$y
        }
        
        cand0 %>%
          dplyr::mutate(
            BP = suppressWarnings(as.integer(BP)),
            classe = normalize_classe(classe),
            label_id = as.character(label_id),
            source_app = as.character(source_app),
            X = bp_to_x_track(BP, fl, x_mode)
          ) %>%
          dplyr::filter(is.finite(BP), is.finite(X)) %>%
          dplyr::mutate(
            hover = paste0(
              "ID: ", label_id,
              "<br>BP: ", BP,
              "<br>classe: ", classe
           #   "<br>source_app: ", source_app
            )
          ) %>%
          dplyr::distinct(BP, classe, label_id, source_app, .keep_all = TRUE)
        
      } else {
        tibble::tibble()
      }
      
      append_log("[LD-PLOT][RENDER] cand nrow=", nrow(cand))
      if (nrow(cand) > 0) {
        cls_tab3 <- table(cand$classe, useNA = "ifany")
        append_log(
          "[LD-PLOT][RENDER] classes: ",
          paste(names(cls_tab3), as.integer(cls_tab3), sep = "=", collapse = " | ")
        )
      }
      
      cand_classes <- unique(cand$classe)
      cand_classes <- cand_classes[!is.na(cand_classes) & nzchar(cand_classes)]
      preferred_order <- c("GWAS", "catalog_hit", "gtex_hit", "nonsyn_hit", "ewastum_hit", "ewasdis_hit", "candidate")
      cand_classes <- unique(c(
        preferred_order[preferred_order %in% cand_classes],
        cand_classes[!(cand_classes %in% preferred_order)]
      ))
      if (!length(cand_classes)) cand_classes <- "candidate"
      
      track_base_y <- seq(0.04, by = 0.10, length.out = length(cand_classes))
      names(track_base_y) <- cand_classes
      
      if (nrow(cand) > 0) {
        cand <- cand %>%
          dplyr::mutate(
            track_y0 = unname(track_base_y[classe]),
            track_y1 = track_y0 + 0.06
          )
      }
      
      x_lab <- if (identical(x_mode, "equal")) {
        x_limits[2] + max(2, ceiling(0.08 * nrow(fl)))
      } else {
        x_limits[2] + 0.08 * x_span
      }
      x_limits2 <- c(x_limits[1], x_lab)
      
      g0 <- NULL
      if (is.finite(as.numeric(payload$chr)) &&
          is.finite(as.numeric(payload$start)) &&
          is.finite(as.numeric(payload$end))) {
        
        g0 <- tryCatch(
          make_gene_track_plot(
            chr_sel   = as.integer(payload$chr),
            st        = as.integer(payload$start),
            en        = as.integer(payload$end),
            fl        = fl,
            x_limits2 = x_limits2,
            x_mode    = x_mode,
            append_log = append_log
          ),
          error = function(e) {
            append_log("[GENE][ERROR] ", conditionMessage(e))
            NULL
          }
        )
      }
      
      
      p <- plotly::plot_ly()
      
      # ===============================
      # LD TRIANGLE
      # ===============================
      p <- p %>%
        plotly::add_trace(
          data = tri,
          x = ~x_mid,
          y = ~y,
          type = "scatter",
          mode = "markers",
          text = ~hover,
          hoverinfo = "text",
          marker = list(
            symbol = "square",
            size = 9,
            color = tri$val,
            cmin = 0,
            cmax = 1,
            colorscale = list(
              c(0.0, "white"),
              c(0.2, "orange"),
              c(1.0, "red")
            ),
            colorbar = list(title = metric)
          ),
          showlegend = FALSE
        )
      
      # ===============================
      # LD BLOCKS
      # ===============================
      bd <- NULL
      
      if (is.data.frame(bd0) && nrow(bd0) > 0) {
        
        bd <- bd0 %>%
          dplyr::mutate(
            i = suppressWarnings(as.integer(i)),
            j = suppressWarnings(as.integer(j))
          ) %>%
          dplyr::filter(
            is.finite(i), is.finite(j),
            i >= 1, j >= 1,
            i <= nrow(fl), j <= nrow(fl),
            j > i
          )
        
        if (nrow(bd) > 0) {
          bd <- bd %>%
            dplyr::mutate(
              xL  = fl$X[i],
              xR  = fl$X[j],
              xM  = (xL + xR) / 2,
              bpL = fl$BP[i],
              bpR = fl$BP[j]
            ) %>%
            dplyr::filter(is.finite(xL), is.finite(xR), is.finite(xM), is.finite(bpL), is.finite(bpR)) %>%
            dplyr::mutate(
              block_label = paste0("B", dplyr::row_number())
            )
          
          if (nrow(bd) > 0) {
            y0_blk <- -0.01
            
            if (identical(x_mode, "equal")) {
              denom <- max(1, nrow(fl) - 1)
              bd <- bd %>% dplyr::mutate(yA = -(abs(j - i) / denom))
            } else {
              bd <- bd %>% dplyr::mutate(yA = -(abs(bpR - bpL) / bp_span))
            }
            
            bd <- bd %>% dplyr::mutate(yA = pmin(yA, y0_blk - 0.02))
            
            p <- p %>%
              plotly::add_segments(
                data = bd,
                x = ~xL, xend = ~xR,
                y = y0_blk, yend = y0_blk,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              ) %>%
              plotly::add_segments(
                data = bd,
                x = ~xL, xend = ~xM,
                y = y0_blk, yend = ~yA,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              ) %>%
              plotly::add_segments(
                data = bd,
                x = ~xR, xend = ~xM,
                y = y0_blk, yend = ~yA,
                inherit = FALSE,
                hoverinfo = "skip",
                showlegend = FALSE,
                line = list(width = 2)
              ) %>%
              plotly::add_text(
                data = bd,
                x = ~xM,
                y = y0_blk + 0.03,
                text = ~block_label,
                textposition = "middle center",
                inherit = FALSE,
                hovertext = ~paste0(block_label, "<br>BP: ", bpL, "-", bpR),
                hoverinfo = "text",
                textfont = list(size = 11),
                showlegend = FALSE
              )
          }
        }
      }
      
      # ===============================
      # TRACKS
      # ===============================
      if (nrow(cand) > 0) {
        for (cls in cand_classes) {
          dsub <- cand %>% dplyr::filter(classe == cls)
          if (!nrow(dsub)) next
          
          y0 <- track_base_y[[cls]]
          y1 <- y0 + 0.06
          
          p <- p %>%
            plotly::add_segments(
              x = x_limits[1],
              xend = x_limits[2],
              y = y0,
              yend = y0,
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE,
              line = list(width = 1)
            ) %>%
            plotly::add_segments(
              data = dsub,
              x = ~X,
              xend = ~X,
              y = y0,
              yend = y1,
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE,
              line = list(width = 1)
            ) %>%
            plotly::add_markers(
              data = dsub,
              x = ~X,
              y = ~track_y1,
              text = ~hover,
              hoverinfo = "text",
              name = cls,
              marker = list(size = 7),
              showlegend = FALSE
            ) %>%
            plotly::add_text(
              x = x_lab,
              y = y1,
              text = cls,
              textposition = "middle left",
              inherit = FALSE,
              hoverinfo = "skip",
              showlegend = FALSE
            )
        }
      }
      
      ymax_tracks <- if (length(track_base_y)) max(track_base_y) + 0.12 else 0.25
      
   #  cl <- selected_cluster()
   #  cluster_lab <- if (is.data.frame(cl) && nrow(cl) == 1 && "cluster_id" %in% names(cl)) {
   #    as.character(cl$cluster_id[[1]])
   #  } else {
   #    "selected cluster"
   #  }
      
      p <- p %>%
        plotly::layout(
          title = list(
            text = paste0("LD plot cluster: — ", cluster_lab),
            y = 0.95
          ),
          hovermode = "closest",
          xaxis = list(
            title = x_title,
            automargin = TRUE,
            range = x_limits2
          ),
          yaxis = list(
            title = "",
            range = c(-1.05, max(0.25, ymax_tracks))
          ),
          margin = list(l = 40, r = 120, t = 70, b = 60)
        )
      
      if (!is.null(g0)) {
        plotly::subplot(
          g0, p,
          nrows = 2,
          heights = c(0.20, 0.80),
          shareX = TRUE,
          titleX = TRUE
        ) %>%
          plotly::layout(hovermode = "closest")
      } else {
        p
      }
    })
    
 
    # ============================================================
    # Block prioritization
    # ============================================================
    block_ranges_selected <- shiny::reactive({
      saved <- ld_cluster_details_for_display()
      
      # --------------------------------------------------
      # 1) prioritat: details.rds carregat
      # --------------------------------------------------
      if (is.list(saved) &&
          is.data.frame(saved$block_ranges) &&
          nrow(saved$block_ranges) > 0) {
        
        br <- tibble::as_tibble(saved$block_ranges)
        
        if (!"cluster_id" %in% names(br)) br$cluster_id <- saved$cluster_id %||% NA_character_
        if (!"chr" %in% names(br)) br$chr <- suppressWarnings(as.integer(saved$chr %||% NA))
        if (!"block_id" %in% names(br)) br$block_id <- NA_character_
        if (!"block_label" %in% names(br)) br$block_label <- NA_character_
        if (!"i" %in% names(br)) br$i <- NA_integer_
        if (!"j" %in% names(br)) br$j <- NA_integer_
        if (!"block_start" %in% names(br)) br$block_start <- NA_integer_
        if (!"block_end" %in% names(br)) br$block_end <- NA_integer_
        if (!"block_size_bp" %in% names(br)) br$block_size_bp <- NA_integer_
        
        return(
          br %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              chr = suppressWarnings(as.integer(chr)),
              block_id = as.character(block_id),
              block_label = dplyr::coalesce(as.character(block_label), ""),
              i = suppressWarnings(as.integer(i)),
              j = suppressWarnings(as.integer(j)),
              block_start = suppressWarnings(as.integer(block_start)),
              block_end = suppressWarnings(as.integer(block_end)),
              block_size_bp = suppressWarnings(as.integer(block_size_bp))
            ) %>%
            dplyr::mutate(
              block_label = dplyr::if_else(
                !nzchar(block_label),
                safe_block_label_from_id(block_id, cluster_id),
                block_label
              ),
              block_size_bp = dplyr::if_else(
                is.finite(block_size_bp),
                block_size_bp,
                pmax(block_end - block_start, 0L)
              )
            ) %>%
            dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
        )
      }
      
      # --------------------------------------------------
      # 2) fallback live state
      # --------------------------------------------------
      cl <- selected_cluster()
      fl <- ld_state$fl
      bij <- ld_state$blocks_ij
      
      if (!is.data.frame(cl) || !nrow(cl) ||
          !is.data.frame(fl) || !nrow(fl) ||
          !is.data.frame(bij) || !nrow(bij)) {
        return(tibble::tibble(
          cluster_id = character(),
          chr = integer(),
          block_id = character(),
          block_label = character(),
          i = integer(),
          j = integer(),
          block_start = integer(),
          block_end = integer(),
          block_size_bp = integer()
        ))
      }
      
      cid <- as.character(cl$cluster_id[1])
      chr_sel <- suppressWarnings(as.integer(cl$chr[1]))
      
      bij %>%
        dplyr::transmute(
          i = suppressWarnings(as.integer(i)),
          j = suppressWarnings(as.integer(j))
        ) %>%
        dplyr::filter(
          is.finite(i), is.finite(j),
          i >= 1, j >= 1,
          i <= nrow(fl), j <= nrow(fl),
          j > i
        ) %>%
        dplyr::mutate(
          cluster_id = cid,
          chr = chr_sel,
          block_start = as.integer(fl$BP[i]),
          block_end   = as.integer(fl$BP[j]),
          block_size_bp = as.integer(block_end - block_start),
          block_label = paste0("B", dplyr::row_number()),
          block_id = paste0(cluster_id, "_", block_label)
        ) %>%
        dplyr::select(cluster_id, chr, block_id, block_label, i, j, block_start, block_end, block_size_bp)
    })
    
    ld_block_hits_selected <- shiny::reactive({
      saved <- ld_cluster_details_for_display()
      
      empty_hits <- tibble::tibble(
        cluster_id = character(),
        chr = integer(),
        lead_snp = character(),
        lead_pos = integer(),
        block_id = character(),
        block_label = character(),
        block_start = integer(),
        block_end = integer(),
        block_size_bp = integer(),
        source_app = character(),
        classe = character(),
        hit_id = character(),
        hit_rsid = character(),
        hit_pos = integer(),
        relation_to_lead = character()
      )
      
      # --------------------------------------------------
      # 1) prioritat: block_hits guardat a details.rds
      # --------------------------------------------------
      if (is.list(saved) &&
          is.data.frame(saved$block_hits) &&
          nrow(saved$block_hits) > 0) {
        
        dd <- tibble::as_tibble(saved$block_hits)
        
        if (!"cluster_id" %in% names(dd)) dd$cluster_id <- saved$cluster_id %||% NA_character_
        if (!"chr" %in% names(dd)) dd$chr <- suppressWarnings(as.integer(saved$chr %||% NA))
        if (!"lead_snp" %in% names(dd)) dd$lead_snp <- NA_character_
        if (!"lead_pos" %in% names(dd)) dd$lead_pos <- NA_integer_
        if (!"block_id" %in% names(dd)) dd$block_id <- NA_character_
        if (!"block_label" %in% names(dd)) dd$block_label <- NA_character_
        if (!"block_start" %in% names(dd)) dd$block_start <- NA_integer_
        if (!"block_end" %in% names(dd)) dd$block_end <- NA_integer_
        if (!"block_size_bp" %in% names(dd)) dd$block_size_bp <- NA_integer_
        if (!"source_app" %in% names(dd)) dd$source_app <- NA_character_
        if (!"classe" %in% names(dd)) dd$classe <- NA_character_
        if (!"hit_id" %in% names(dd)) dd$hit_id <- NA_character_
        if (!"hit_rsid" %in% names(dd)) dd$hit_rsid <- NA_character_
        if (!"hit_pos" %in% names(dd)) dd$hit_pos <- NA_integer_
        if (!"relation_to_lead" %in% names(dd)) dd$relation_to_lead <- NA_character_
        
        return(
          dd %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              chr = suppressWarnings(as.integer(chr)),
              lead_snp = as.character(lead_snp),
              lead_pos = suppressWarnings(as.integer(lead_pos)),
              block_id = as.character(block_id),
              block_label = dplyr::coalesce(as.character(block_label), ""),
              block_start = suppressWarnings(as.integer(block_start)),
              block_end = suppressWarnings(as.integer(block_end)),
              block_size_bp = suppressWarnings(as.integer(block_size_bp)),
              source_app = as.character(source_app),
              classe = as.character(classe),
              hit_id = as.character(hit_id),
              hit_rsid = as.character(hit_rsid),
              hit_pos = suppressWarnings(as.integer(hit_pos)),
              relation_to_lead = as.character(relation_to_lead)
            ) %>%
            dplyr::distinct() %>%
            dplyr::arrange(lead_pos, block_start, hit_pos, classe, hit_rsid)
        )
      }
      
      # --------------------------------------------------
      # 2) reconstrucció des de saved details
      # --------------------------------------------------
      if (is.list(saved)) {
        out_saved <- rebuild_block_hits_from_components(
          cl = selected_cluster(),
          ca = if (is.data.frame(saved$candidates)) saved$candidates else NULL,
          hits = if (is.data.frame(saved$gwas_hits)) saved$gwas_hits else NULL,
          blk = if (is.data.frame(saved$block_ranges)) saved$block_ranges else NULL,
          proxy_tbl = if (is.data.frame(saved$proxies)) saved$proxies else NULL
        )
        
        if (is.data.frame(out_saved) && nrow(out_saved)) {
          return(out_saved)
        }
      }
      
      # --------------------------------------------------
      # 3) fallback live state
      # --------------------------------------------------
      out_live <- rebuild_block_hits_from_components(
        cl = selected_cluster(),
        ca = selected_candidates(),
        hits = selected_gwas_hits(),
        blk = block_ranges_selected(),
        proxy_tbl = ld_state$proxy_table
      )
      
      if (!is.data.frame(out_live) || !nrow(out_live)) {
        return(empty_hits)
      }
      
      out_live
    })
  
    observe({
      x <- ld_block_hits_selected()
      req(is.data.frame(x))
      cat("\n[ld_block_hits_selected] column names:\n")
      print(names(x))
      cat("\n")
    })
    
    
    ld_block_overlap_detail_selected <- shiny::reactive({
      dd <- ld_block_hits_selected()
      summarize_hits_by_block(dd)
    })
    
    ld_block_overlap_summary_selected <- shiny::reactive({
      live_block_summary_r()
    })
    
    
    ld_block_summary_merged_r <- shiny::reactive({
      ld_df <- live_block_summary_r()
      
      detail_df <- ld_block_overlap_detail_selected()
      detail_sum <- summarize_block_overlap_detail_for_dt(detail_df)
      
      hit_link_sum <- summarize_block_hit_links_for_block(
        ld_block_hits_selected()
      )
      
      if (!is.data.frame(ld_df) || !nrow(ld_df)) {
        return(tibble::tibble())
      }
      
      out <- ld_df
      
      if (is.data.frame(detail_sum) && nrow(detail_sum)) {
        out <- out %>%
          dplyr::left_join(
            detail_sum %>%
              dplyr::select(
                cluster_id,
                block_id,
                n_lead_snps,
                lead_snps,
                lead_positions,
                gwas_hits,
                catalog_hits,
                gtex_hits,
                nonsyn_hits,
                ewasdis_hits,
                ewastum_hits,
                n_external_apps,
                external_apps,
                marker_status
              ),
            by = c("cluster_id", "block_id")
          )
      }
      
      if (is.data.frame(hit_link_sum) && nrow(hit_link_sum)) {
        out <- out %>%
          dplyr::left_join(
            hit_link_sum,
            by = c("cluster_id", "block_id")
          )
      }
      
      out
    })
    
    # ============================================================
    # HAPLOBLOCK SUPPORT 
    # ============================================================
    
    # for cluster
    block_gene_overlap_summary <- reactive({
      saved <- ld_cluster_details_for_display()
      
      empty_genes <- tibble::tibble(
        cluster_id = character(),
        block_id = character(),
        n_genes_in_block = integer(),
        genes_in_block = character()
      )
      
      # --------------------------------------------------
      # 1) prioritat: details.rds carregat
      # --------------------------------------------------
      if (is.list(saved) &&
          is.data.frame(saved$block_genes) &&
          nrow(saved$block_genes) > 0) {
        
        bg <- tibble::as_tibble(saved$block_genes)
        cid <- as.character(saved$cluster_id %||% NA_character_)
        
        if (!"cluster_id" %in% names(bg)) {
          bg$cluster_id <- cid
        }
        
        if (!"block_id" %in% names(bg)) {
          return(empty_genes)
        }
        
        if (all(c("n_genes_in_block", "genes_in_block") %in% names(bg))) {
          return(
            bg %>%
              dplyr::transmute(
                cluster_id = as.character(cluster_id),
                block_id = as.character(block_id),
                n_genes_in_block = safe_int0(n_genes_in_block),
                genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
              ) %>%
              dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
          )
        }
        
        gene_col <- intersect(c("gene", "gene_symbol", "symbol", "label"), names(bg))
        gene_col <- if (length(gene_col)) gene_col[1] else NA_character_
        
        if (is.na(gene_col)) {
          return(empty_genes)
        }
        
        return(
          bg %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              block_id = as.character(block_id),
              gene = trimws(as.character(.data[[gene_col]]))
            ) %>%
            dplyr::filter(
              !is.na(cluster_id), nzchar(cluster_id),
              !is.na(block_id), nzchar(block_id),
              !is.na(gene), nzchar(gene)
            ) %>%
            dplyr::distinct(cluster_id, block_id, gene) %>%
            dplyr::group_by(cluster_id, block_id) %>%
            dplyr::summarise(
              n_genes_in_block = dplyr::n(),
              genes_in_block = paste(sort(unique(gene)), collapse = "; "),
              .groups = "drop"
            )
        )
      }
      
      # --------------------------------------------------
      # 2) fallback live state
      # --------------------------------------------------
      bg <- ld_state$block_genes
      
      if (!is.data.frame(bg) || !nrow(bg)) {
        return(empty_genes)
      }
      
      cl <- selected_cluster()
      cid <- if (is.data.frame(cl) && nrow(cl)) as.character(cl$cluster_id[1]) else NA_character_
      
      bg <- tibble::as_tibble(bg)
      
      if (!"cluster_id" %in% names(bg)) {
        bg$cluster_id <- cid
      }
      
      if (!"block_id" %in% names(bg)) {
        return(empty_genes)
      }
      
      if (all(c("n_genes_in_block", "genes_in_block") %in% names(bg))) {
        return(
          bg %>%
            dplyr::transmute(
              cluster_id = as.character(cluster_id),
              block_id = as.character(block_id),
              n_genes_in_block = safe_int0(n_genes_in_block),
              genes_in_block = dplyr::coalesce(as.character(genes_in_block), "")
            ) %>%
            dplyr::distinct(cluster_id, block_id, .keep_all = TRUE)
        )
      }
      
      gene_col <- intersect(c("gene", "gene_symbol", "symbol", "label"), names(bg))
      gene_col <- if (length(gene_col)) gene_col[1] else NA_character_
      
      if (is.na(gene_col)) {
        return(empty_genes)
      }
      
      bg %>%
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          block_id = as.character(block_id),
          gene = trimws(as.character(.data[[gene_col]]))
        ) %>%
        dplyr::filter(
          !is.na(cluster_id), nzchar(cluster_id),
          !is.na(block_id), nzchar(block_id),
          !is.na(gene), nzchar(gene)
        ) %>%
        dplyr::distinct(cluster_id, block_id, gene) %>%
        dplyr::group_by(cluster_id, block_id) %>%
        dplyr::summarise(
          n_genes_in_block = dplyr::n(),
          genes_in_block = paste(sort(unique(gene)), collapse = "; "),
          .groups = "drop"
        )
    })
    
    ld_pairs_for_summary_r <- shiny::reactive({
      saved <- ld_cluster_details_for_display()
      
      if (is.list(saved) &&
          is.data.frame(saved$ld_pairs) &&
          nrow(saved$ld_pairs) > 0) {
        return(saved$ld_pairs)
      }
      
      if (is.data.frame(ld_state$ld_pairs) && nrow(ld_state$ld_pairs) > 0) {
        return(ld_state$ld_pairs)
      }
      
      tibble::tibble()
    })
    
      
    live_block_summary_r <- shiny::reactive({
      base_df <- block_overlap_summary_df(
        block_ranges = block_ranges_selected(),
        block_hits = ld_block_hits_selected(),
        block_genes = block_gene_overlap_summary(),
        gwas_bridge_df = gwas_bridge_base_r(),
        proxy_tbl = ld_state$proxy_table
      )
      
      base_df <- build_block_ld_summary(base_df)
      
      ld_struct_df <- summarize_ld_structure_by_block(
        block_df = block_ranges_selected(),
        ld_pairs_df = ld_pairs_for_summary_r()
      )
      
      if (!is.data.frame(base_df) || !nrow(base_df)) {
        return(tibble::tibble())
      }
      
      if (!is.data.frame(ld_struct_df) || !nrow(ld_struct_df)) {
        return(base_df)
      }
      
      base_df %>%
        dplyr::select(
          -dplyr::any_of(c(
            "n_ld_proxy_hits",
            "block_max_ld_value",
            "block_mean_ld_value",
            "ld_proxy_snps"
          ))
        ) %>%
        dplyr::left_join(
          ld_struct_df,
          by = c("cluster_id", "block_id")
        ) %>%
        dplyr::mutate(
          n_ld_proxy_hits = safe_int0(n_ld_proxy_hits),
          block_max_ld_value = suppressWarnings(as.numeric(block_max_ld_value)),
          block_mean_ld_value = suppressWarnings(as.numeric(block_mean_ld_value)),
          ld_proxy_snps = dplyr::coalesce(as.character(ld_proxy_snps), "")
        )
    })  
    
    # pendent de revisar
    ld_block_overlap_summary_for_display <- reactive({
      obj <- cluster_summary_saved_r()
      
      required_cols <- c(
        "cluster_id",
        "block_id",
        "block_label",
        "block_priority_score",
        "priority_class",
        "block_start",
        "block_end",
        "block_size_bp",
        "block_size_kb",
        "block_support_any",
        "block_support_apps",
        "block_support_hits",
        "block_hit_density",
        "block_catalog_hits",
        "block_gtex_hits",
        "block_nonsyn_hits",
        "block_ewasdis_hits",
        "block_ewastum_hits",
        "n_gwas_sig_hits",
        "max_gwas_logp",
        "mean_gwas_logp",
        "gwas_sig_hits",
        "n_genes_in_block",
        "n_ld_proxy_hits",
        "block_max_ld_value",
        "block_mean_ld_value",
        "ld_proxy_snps",
        "genes_in_block"
      )
      
      if (is.data.frame(obj) && nrow(obj) > 0 && length(setdiff(required_cols, names(obj))) == 0) {
        return(obj)
      }
      
      ld_block_overlap_summary_selected()
    })
    
# Table clusters LD
    output$block_overlap_summary <- DT::renderDT({
      df <- ld_block_overlap_summary_for_display()
      
      validate(
        need(is.data.frame(df) && nrow(df) > 0, "No LD block summary available.")
      )
      
      show_df <- make_block_canonical_table_display(df)
      
      if (!"genes_name" %in% names(show_df)) show_df$genes_name <- ""
      if (!"gwas_hits" %in% names(show_df)) show_df$gwas_hits <- ""
      if (!"ld_proxy_snps" %in% names(show_df)) show_df$ld_proxy_snps <- ""
      
      show_df <- show_df %>%
        dplyr::mutate(
          gwas_min_p = dplyr::if_else(
            is.finite(gwas_min_p),
            formatC(gwas_min_p, format = "e", digits = 2),
            ""
          ),
          genes_name = vapply(genes_name, build_gene_link_column_html, character(1)),
          gwas_hits = vapply(
            gwas_hits,
            function(z) collapse_plain_list_html(
              z,
              max_visible = 5L,
              summary_label = paste0(length(split_semicolon(z)), " GWAS hits")
            ),
            character(1)
          ),
          ld_proxy_snps = vapply(ld_proxy_snps, build_snp_link_column_html, character(1))
        )
      
      DT::datatable(
        show_df,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        escape = FALSE,
        extensions = "Buttons",
        callback = DT::JS("
      table.columns.adjust();
      $(window).on('resize', function() {
        table.columns.adjust();
      });
      $(document).on('shown.bs.tab', 'button[data-bs-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function() {
        $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
      });
    "),
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          autoWidth = FALSE,
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          order = list(list(which(names(show_df) == "block_mean_ld") - 1L, "desc"))
        )
      ) %>%
        DT::formatRound(
          c(
            "start", "end",
            "support_apps",
            "catalog_hits", "gtex_hits", "nonsyn_hits", "ewasdis_hits", "ewastum_hits",
            "n_ld_proxy_hits"
          ),
          0
        ) %>%
        DT::formatRound(
          c(
            "size_kb",
            "gwas_top_logp",
            "block_max_ld",
            "block_mean_ld"
          ),
          2
        )
    }, server = FALSE)
    

    output$ld_block_hits_dt <- DT::renderDT({
      dd <- ld_block_hits_selected()
      blk <- block_ranges_selected()
      
      # Cas 1: encara no s'ha fet cap RUN LD/blocks per al cluster
      if (!is.data.frame(blk) || !nrow(blk)) {
        return(DT::datatable(
          data.frame(Message = "Run LD & blocks to populate block-level hits."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      
      # Cas 2: LD/blocks existeix però no hi ha hits per mostrar
      if (!is.data.frame(dd) || !nrow(dd)) {
        return(DT::datatable(
          data.frame(Message = "No block-level hits found for this cluster."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      
      show_dt <- format_block_hits_grouped_for_dt(dd) %>%
        dplyr::select(
          cluster_id,
          lead_snp,
          lead_pos,
          block_id,
          block_label,
          block_start,
          block_end,
          block_size_bp,
          n_direct,
          n_ld,
          n_block,
          direct_hits,
          ld_hits,
          block_hits,
          gwas_hits,
          catalog_hits,
          gtex_hits,
          nonsyn_hits,
          ewasdis_hits,
          ewastum_hits,
          source_apps
        )
      
      # Cas 3: hi ha blocs, però després del format no queda cap fila visible
      if (!is.data.frame(show_dt) || !nrow(show_dt)) {
        return(DT::datatable(
          data.frame(Message = "No block-level hits found for this cluster."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      
      DT::datatable(
        show_dt,
        rownames = FALSE,
        escape = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy",  exportOptions = list(modifier = list(page = "all"))),
            list(extend = "csv",   exportOptions = list(modifier = list(page = "all"))),
            list(extend = "excel", exportOptions = list(modifier = list(page = "all")))
          ),
          pageLength = 15,
          scrollX = TRUE,
          scrollY = "350px",
          autoWidth = FALSE
        ),
        class = "compact stripe hover order-column",
        fillContainer = TRUE
      )
    }, server = FALSE)
    
 #####################   
    invisible(list(
      integrated_clusters = integrated_clusters_r,
      integrated_candidates = integrated_candidates_r,
      selected_cluster = selected_cluster,
      selected_candidates = selected_candidates,
      selected_gwas_hits = selected_gwas_hits,
      seed_table = shiny::reactive(ld_state$seeds),
      proxy_table = shiny::reactive(ld_state$proxy_table),
      block_ranges_selected = block_ranges_selected,
      ld_block_hits_selected = ld_block_hits_selected,
      block_gene_overlap_summary = block_gene_overlap_summary
    ))
  })
 
  
}