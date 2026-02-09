# R/format_dbnsfp.R
# Normaliza/da formato a la salida de dbNSFP (.out) y escribe CSV+RDS.
# - Convierte a numérico un conjunto de columnas por índice.
# - Renombra columnas si coincide el recuento.
# - Deduplica valores separados por ';' en columnas multivalor.
# Devuelve una lista con rutas y un preview (data.frame).

# ---------- helpers ----------
.make_numeric_columns <- function(df, idxs) {
  idxs <- idxs[idxs %in% seq_len(ncol(df))]
  if (!length(idxs)) return(df)
  df[, idxs] <- lapply(df[, idxs, drop = FALSE], function(x) {
    if (is.numeric(x)) return(x)
    suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", as.character(x))))
  })
  df
}

# Deduplica valores separados por ';' dentro de cada celda, manteniendo el orden.
# - trim de espacios
# - quita vacíos y "NA" literales
# - opcionalmente fuerza a minúsculas (tolower_vals = TRUE)
dedup_semicol <- function(x, sep = ";", tolower_vals = FALSE) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  x <- ifelse(is.na(x), "", x)
  vapply(strsplit(x, sep, fixed = TRUE), FUN.VALUE = character(1), function(parts){
    parts <- trimws(parts)
    if (tolower_vals) parts <- tolower(parts)
    parts <- parts[nzchar(parts) & parts != "NA"]
    if (!length(parts)) return("")
    parts <- parts[!duplicated(parts)]  # mantiene el orden
    paste(parts, collapse = sep)
  })
}

# ---------- configuración por defecto ----------
# Índices de columnas a numeric (ajústalo si tu dbNSFP cambia de versión)
.dbnsfp_numeric_idxs <- c(
  1, 2, 8:12, 30, 31, 37, 38, 40, 41, 43, 44, 46, 47,
  49, 50, 53:56, 58, 59, 61:64, 66, 67, 69:71, 73, 74,
  76:79, 83:90, 92, 93, 95, 96, 98, 99, 101, 102, 104, 105,
  107:116, 118, 119, 121:126, 133:139, 142:180,
  182:370, 385, 417, 418, 425:440, 444:446, 448, 449,
  460, 465, 468:517
)

# Nombres de columnas esperados (si coincide ncol, se renombra en bloque)
.dbnsfp_colnames <- c(
  "chr","BP","ref","alt","aaref","aaalt","SNP","hg19_chr","hg19_pos_1_based","hg18_chr","hg18_pos_1_based",
  "aapos","genename","Ensembl_geneid","Ensembl_transcriptid","Ensembl_proteinid","Uniprot_acc1","Uniprot_entry",
  "HGVSc_snpEff","HGVSp_snpEff","HGVSc_VEP","HGVSp_VEP","APPRIS","GENCODE_basic","TSL","VEP_canonical","MANE",
  "cds_strand","refcodon","codonpos","codon_degeneracy","Ancestral_allele","AltaiNeandertal","Denisova",
  "VindijiaNeandertal","ChagyrskayaNeandertal","SIFT_score","SIFT_converted_rankscore","SIFT_pred","SIFT4G_score",
  "SIFT4G_converted_rankscore","SIFT4G_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_rankscore",
  "Polyphen2_HDIV_pred","Polyphen2_HVAR_score","Polyphen2_HVAR_rankscore","Polyphen2_HVAR_pred",
  "MutationTaster_score","MutationTaster_converted_rankscore","MutationTaster_pred","MutationTaster_model",
  "MutationTaster_trees_benign","MutationTaster_trees_deleterious","MutationAssessor_score",
  "MutationAssessor_rankscore","MutationAssessor_pred","PROVEAN_score","PROVEAN_converted_rankscore",
  "PROVEAN_pred","VEST4_score","VEST4_rankscore","MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred",
  "MetaLR_score","MetaLR_rankscore","MetaLR_pred","Reliability_index","MetaRNN_score","MetaRNN_rankscore",
  "MetaRNN_pred","M_CAP_score","M_CAP_rankscore","M_CAP_pred","REVEL_score","REVEL_rankscore","MutPred_score",
  "MutPred_rankscore","MutPred_protID","MutPred_AAchange","MutPred_Top5features","MVP_score","MVP_rankscore",
  "gMVP_score","gMVP_rankscore","MPC_score","MPC_rankscore","PrimateAI_score","PrimateAI_rankscore",
  "PrimateAI_pred","DEOGEN2_score","DEOGEN2_rankscore","DEOGEN2_pred","BayesDel_addAF_score",
  "BayesDel_addAF_rankscore","BayesDel_addAF_pred","BayesDel_noAF_score","BayesDel_noAF_rankscore",
  "BayesDel_noAF_pred","ClinPred_score","ClinPred_rankscore","ClinPred_pred","LIST_S2_score",
  "LIST_S2_rankscore","LIST_S2_pred","VARITY_R_score","VARITY_R_rankscore","VARITY_ER_score",
  "VARITY_ER_rankscore","VARITY_R_LOO_score","VARITY_R_LOO_rankscore","VARITY_ER_LOO_score",
  "VARITY_ER_LOO_rankscore","ESM1b_score","ESM1b_rankscore","ESM1b_pred","AlphaMissense_score",
  "AlphaMissense_rankscore","AlphaMissense_pred","PHACTboost_score","PHACTboost_rankscore","MutFormer_score",
  "MutFormer_rankscore","MutScore_score","MutScore_rankscore","Aloft_Fraction_transcripts_affected",
  "Aloft_prob_Tolerant","Aloft_prob_Recessive","Aloft_prob_Dominant","Aloft_pred","Aloft_Confidence",
  "CADD_raw","CADD_raw_rankscore","CADD_phred","DANN_score","DANN_rankscore","fathmm_XF_coding_score",
  "fathmm_XF_coding_rankscore","fathmm_XF_coding_pred","Eigen_raw_coding","Eigen_raw_coding_rankscore",
  "Eigen_phred_coding","Eigen_PC_raw_coding","Eigen_PC_raw_coding_rankscore","Eigen_PC_phred_coding",
  "GERP++_NR","GERP++_RS","GERP++_RS_rankscore","GERP_91_mammals","GERP_91_mammals_rankscore",
  "phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore","phyloP470way_mammalian",
  "phyloP470way_mammalian_rankscore","phyloP17way_primate","phyloP17way_primate_rankscore",
  "phastCons100way_vertebrate","phastCons100way_vertebrate_rankscore","phastCons470way_mammalian",
  "phastCons470way_mammalian_rankscore","phastCons17way_primate","phastCons17way_primate_rankscore",
  "bStatistic","bStatistic_converted_rankscore","KGp3_AC","KGp3_AF","KGp3_AFR_AC","KGp3_AFR_AF",
  "KGp3_EUR_AC","KGp3_EUR_AF","KGp3_AMR_AC","KGp3_AMR_AF","KGp3_EAS_AC","KGp3_EAS_AF","KGp3_SAS_AC",
  "KGp3_SAS_AF","TOPMed_frz8_AC","TOPMed_frz8_AN","TOPMed_frz8_AF","gnomAD2_1_1_exomes_flag",
  "gnomAD2_1_1_exomes_controls_AC","gnomAD2_1_1_exomes_controls_AN","gnomAD2_1_1_exomes_controls_AF",
  "gnomAD2_1_1_exomes_controls_nhomalt","gnomAD2_1_1_exomes_non_neuro_AC","gnomAD2_1_1_exomes_non_neuro_AN",
  "gnomAD2_1_1_exomes_non_neuro_AF","gnomAD2_1_1_exomes_non_neuro_nhomalt","gnomAD2_1_1_exomes_non_cancer_AC",
  "gnomAD2_1_1_exomes_non_cancer_AN","gnomAD2_1_1_exomes_non_cancer_AF","gnomAD2_1_1_exomes_non_cancer_nhomalt",
  "gnomAD2_1_1_exomes_controls_AFR_AC","gnomAD2_1_1_exomes_controls_AFR_AN","gnomAD2_1_1_exomes_controls_AFR_AF",
  "gnomAD2_1_1_exomes_controls_AFR_nhomalt","gnomAD2_1_1_exomes_controls_AMR_AC",
  "gnomAD2_1_1_exomes_controls_AMR_AN","gnomAD2_1_1_exomes_controls_AMR_AF",
  "gnomAD2_1_1_exomes_controls_AMR_nhomalt","gnomAD2_1_1_exomes_controls_ASJ_AC",
  "gnomAD2_1_1_exomes_controls_ASJ_AN","gnomAD2_1_1_exomes_controls_ASJ_AF",
  "gnomAD2_1_1_exomes_controls_ASJ_nhomalt","gnomAD2_1_1_exomes_controls_EAS_AC",
  "gnomAD2_1_1_exomes_controls_EAS_AN","gnomAD2_1_1_exomes_controls_EAS_AF",
  "gnomAD2_1_1_exomes_controls_EAS_nhomalt","gnomAD2_1_1_exomes_controls_FIN_AC",
  "gnomAD2_1_1_exomes_controls_FIN_AN","gnomAD2_1_1_exomes_controls_FIN_AF",
  "gnomAD2_1_1_exomes_controls_FIN_nhomalt","gnomAD2_1_1_exomes_controls_NFE_AC",
  "gnomAD2_1_1_exomes_controls_NFE_AN","gnomAD2_1_1_exomes_controls_NFE_AF",
  "gnomAD2_1_1_exomes_controls_NFE_nhomalt","gnomAD2_1_1_exomes_controls_SAS_AC",
  "gnomAD2_1_1_exomes_controls_SAS_AN","gnomAD2_1_1_exomes_controls_SAS_AF",
  "gnomAD2_1_1_exomes_controls_SAS_nhomalt","gnomAD2_1_1_exomes_controls_POPMAX_AC",
  "gnomAD2_1_1_exomes_controls_POPMAX_AN","gnomAD2_1_1_exomes_controls_POPMAX_AF",
  "gnomAD2_1_1_exomes_controls_POPMAX_nhomalt","gnomAD2_1_1_exomes_non_neuro_AFR_AC",
  "gnomAD2_1_1_exomes_non_neuro_AFR_AN","gnomAD2_1_1_exomes_non_neuro_AFR_AF",
  "gnomAD2_1_1_exomes_non_neuro_AFR_nhomalt","gnomAD2_1_1_exomes_non_neuro_AMR_AC",
  "gnomAD2_1_1_exomes_non_neuro_AMR_AN","gnomAD2_1_1_exomes_non_neuro_AMR_AF",
  "gnomAD2_1_1_exomes_non_neuro_AMR_nhomalt","gnomAD2_1_1_exomes_non_neuro_ASJ_AC",
  "gnomAD2_1_1_exomes_non_neuro_ASJ_AN","gnomAD2_1_1_exomes_non_neuro_ASJ_AF",
  "gnomAD2_1_1_exomes_non_neuro_ASJ_nhomalt","gnomAD2_1_1_exomes_non_neuro_EAS_AC",
  "gnomAD2_1_1_exomes_non_neuro_EAS_AN","gnomAD2_1_1_exomes_non_neuro_EAS_AF",
  "gnomAD2_1_1_exomes_non_neuro_EAS_nhomalt","gnomAD2_1_1_exomes_non_neuro_FIN_AC",
  "gnomAD2_1_1_exomes_non_neuro_FIN_AN","gnomAD2_1_1_exomes_non_neuro_FIN_AF",
  "gnomAD2_1_1_exomes_non_neuro_FIN_nhomalt","gnomAD2_1_1_exomes_non_neuro_NFE_AC",
  "gnomAD2_1_1_exomes_non_neuro_NFE_AN","gnomAD2_1_1_exomes_non_neuro_NFE_AF",
  "gnomAD2_1_1_exomes_non_neuro_NFE_nhomalt","gnomAD2_1_1_exomes_non_neuro_SAS_AC",
  "gnomAD2_1_1_exomes_non_neuro_SAS_AN","gnomAD2_1_1_exomes_non_neuro_SAS_AF",
  "gnomAD2_1_1_exomes_non_neuro_SAS_nhomalt","gnomAD2_1_1_exomes_non_neuro_POPMAX_AC",
  "gnomAD2_1_1_exomes_non_neuro_POPMAX_AN","gnomAD2_1_1_exomes_non_neuro_POPMAX_AF",
  "gnomAD2_1_1_exomes_non_neuro_POPMAX_nhomalt","gnomAD2_1_1_exomes_non_cancer_AFR_AC",
  "gnomAD2_1_1_exomes_non_cancer_AFR_AN","gnomAD2_1_1_exomes_non_cancer_AFR_AF",
  "gnomAD2_1_1_exomes_non_cancer_AFR_nhomalt","gnomAD2_1_1_exomes_non_cancer_AMR_AC",
  "gnomAD2_1_1_exomes_non_cancer_AMR_AN","gnomAD2_1_1_exomes_non_cancer_AMR_AF",
  "gnomAD2_1_1_exomes_non_cancer_AMR_nhomalt","gnomAD2_1_1_exomes_non_cancer_ASJ_AC",
  "gnomAD2_1_1_exomes_non_cancer_ASJ_AN","gnomAD2_1_1_exomes_non_cancer_ASJ_AF",
  "gnomAD2_1_1_exomes_non_cancer_ASJ_nhomalt","gnomAD2_1_1_exomes_non_cancer_EAS_AC",
  "gnomAD2_1_1_exomes_non_cancer_EAS_AN","gnomAD2_1_1_exomes_non_cancer_EAS_AF",
  "gnomAD2_1_1_exomes_non_cancer_EAS_nhomalt","gnomAD2_1_1_exomes_non_cancer_FIN_AC",
  "gnomAD2_1_1_exomes_non_cancer_FIN_AN","gnomAD2_1_1_exomes_non_cancer_FIN_AF",
  "gnomAD2_1_1_exomes_non_cancer_FIN_nhomalt","gnomAD2_1_1_exomes_non_cancer_NFE_AC",
  "gnomAD2_1_1_exomes_non_cancer_NFE_AN","gnomAD2_1_1_exomes_non_cancer_NFE_AF",
  "gnomAD2_1_1_exomes_non_cancer_NFE_nhomalt","gnomAD2_1_1_exomes_non_cancer_SAS_AC",
  "gnomAD2_1_1_exomes_non_cancer_SAS_AN","gnomAD2_1_1_exomes_non_cancer_SAS_AF",
  "gnomAD2_1_1_exomes_non_cancer_SAS_nhomalt","gnomAD2_1_1_exomes_non_cancer_POPMAX_AC",
  "gnomAD2_1_1_exomes_non_cancer_POPMAX_AN","gnomAD2_1_1_exomes_non_cancer_POPMAX_AF",
  "gnomAD2_1_1_exomes_non_cancer_POPMAX_nhomalt","gnomAD4_1_joint_flag","gnomAD4_1_joint_AC",
  "gnomAD4_1_joint_AN","gnomAD4_1_joint_AF","gnomAD4_1_joint_nhomalt","gnomAD4_1_joint_POPMAX_AC",
  "gnomAD4_1_joint_POPMAX_AN","gnomAD4_1_joint_POPMAX_AF","gnomAD4_1_joint_POPMAX_nhomalt",
  "gnomAD4_1_joint_AFR_AC","gnomAD4_1_joint_AFR_AN","gnomAD4_1_joint_AFR_AF","gnomAD4_1_joint_AFR_nhomalt",
  "gnomAD4_1_joint_AMI_AC","gnomAD4_1_joint_AMI_AN","gnomAD4_1_joint_AMI_AF","gnomAD4_1_joint_AMI_nhomalt",
  "gnomAD4_1_joint_AMR_AC","gnomAD4_1_joint_AMR_AN","gnomAD4_1_joint_AMR_AF","gnomAD4_1_joint_AMR_nhomalt",
  "gnomAD4_1_joint_ASJ_AC","gnomAD4_1_joint_ASJ_AN","gnomAD4_1_joint_ASJ_AF","gnomAD4_1_joint_ASJ_nhomalt",
  "gnomAD4_1_joint_EAS_AC","gnomAD4_1_joint_EAS_AN","gnomAD4_1_joint_EAS_AF","gnomAD4_1_joint_EAS_nhomalt",
  "gnomAD4_1_joint_FIN_AC","gnomAD4_1_joint_FIN_AN","gnomAD4_1_joint_FIN_AF","gnomAD4_1_joint_FIN_nhomalt",
  "gnomAD4_1_joint_MID_AC","gnomAD4_1_joint_MID_AN","gnomAD4_1_joint_MID_AF","gnomAD4_1_joint_MID_nhomalt",
  "gnomAD4_1_joint_NFE_AC","gnomAD4_1_joint_NFE_AN","gnomAD4_1_joint_NFE_AF","gnomAD4_1_joint_NFE_nhomalt",
  "gnomAD4_1_joint_SAS_AC","gnomAD4_1_joint_SAS_AN","gnomAD4_1_joint_SAS_AF","gnomAD4_1_joint_SAS_nhomalt",
  "ALFA_European_AC","ALFA_European_AN","ALFA_European_AF","ALFA_African_Others_AC","ALFA_African_Others_AN",
  "ALFA_African_Others_AF","ALFA_East_Asian_AC","ALFA_East_Asian_AN","ALFA_East_Asian_AF",
  "ALFA_African_American_AC","ALFA_African_American_AN","ALFA_African_American_AF",
  "ALFA_Latin_American_1_AC","ALFA_Latin_American_1_AN","ALFA_Latin_American_1_AF",
  "ALFA_Latin_American_2_AC","ALFA_Latin_American_2_AN","ALFA_Latin_American_2_AF",
  "ALFA_Other_Asian_AC","ALFA_Other_Asian_AN","ALFA_Other_Asian_AF","ALFA_South_Asian_AC",
  "ALFA_South_Asian_AN","ALFA_South_Asian_AF","ALFA_Other_AC","ALFA_Other_AN","ALFA_Other_AF",
  "ALFA_African_AC","ALFA_African_AN","ALFA_African_AF","ALFA_Asian_AC","ALFA_Asian_AN","ALFA_Asian_AF",
  "ALFA_Total_AC","ALFA_Total_AN","ALFA_Total_AF","clinvar_id","clinvar_clnsig","clinvar_trait",
  "clinvar_review","clinvar_hgvs","clinvar_var_source","clinvar_MedGen_id","clinvar_OMIM_id",
  "clinvar_Orphanet_id","Interpro_domain","Gene_old_names","Gene_other_names","Uniprot_acc2",
  "Uniprot_id","Entrez_gene_id","CCDS_id","Refseq_id","ucsc_id","MIM_id","OMIM_id","Gene_full_name",
  "Pathway_Uniprot","Pathway_BioCarta_short","Pathway_BioCarta_full","Pathway_ConsensusPathDB",
  "Pathway_KEGG_id","Pathway_KEGG_full","Function_description","Disease_description","MIM_phenotype_id",
  "MIM_disease","Orphanet_disorder_id","Orphanet_disorder","Orphanet_association_type",
  "Trait_association_GWAS","MGI_mouse_gene","MGI_mouse_phenotype","ZFIN_zebrafish_gene",
  "ZFIN_zebrafish_structure","ZFIN_zebrafish_phenotype_quality","ZFIN_zebrafish_phenotype_tag",
  "HPO_id","HPO_name","GO_biological_process","GO_cellular_component","GO_molecular_function",
  "P_HI","HIPred_score","HIPred","GHIS","ClinGen_Haploinsufficiency_Score","ClinGen_Haploinsufficiency_Description",
  "ClinGen_Haploinsufficiency_PMID","ClinGen_Haploinsufficiency_Disease","P_rec","Known_rec_info",
  "RVIS_EVS","RVIS_percentile_EVS","LoF_FDR_ExAC","RVIS_ExAC","RVIS_percentile_ExAC","ExAC_pLI",
  "ExAC_pRec","ExAC_pNull","ExAC_nonTCGA_pLI","ExAC_nonTCGA_pRec","ExAC_nonTCGA_pNull","ExAC_nonpsych_pLI",
  "ExAC_nonpsych_pRec","ExAC_nonpsych_pNull","gnomAD_pLI","gnomAD_pRec","gnomAD_pNull","ExAC_del_score",
  "ExAC_dup_score","ExAC_cnv_score","ExAC_cnv_flag","GDI","GDI_Phred","Gene_damage_prediction_all_disease_causing_genes",
  "Gene_damage_prediction_all_Mendelian_disease_causing_genes",
  "Gene_damage_prediction_Mendelian_AD_disease_causing_genes",
  "Gene_damage_prediction_Mendelian_AR_disease_causing_genes",
  "Gene_damage_prediction_all_PID_disease_causing_genes",
  "Gene_damage_prediction_PID_AD_disease_causing_genes",
  "Gene_damage_prediction_PID_AR_disease_causing_genes",
  "Gene_damage_prediction_all_cancer_disease_causing_genes",
  "Gene_damage_prediction_cancer_recessive_disease_causing_genes",
  "Gene_damage_prediction_cancer_dominant_disease_causing_genes",
  "LoFtool_score","Essential_gene","Essential_gene_CRISPR",
  "Essential_gene_CRISPR2","Essential_gene_gene_trap","Gene_indispensability_score",
  "Gene_indispensability_pred","Tissue_specificity_Uniprot","HPA_consensus_adipose_tissue",
  "HPA_consensus_adrenal_gland","HPA_consensus_amygdala","HPA_consensus_appendix",
  "HPA_consensus_basal_ganglia","HPA_consensus_bone_marrow","HPA_consensus_breast",
  "HPA_consensus_cerebellum","HPA_consensus_cerebral_cortex","HPA_consensus_cervix",
  "HPA_consensus_choroid_plexus","HPA_consensus_colon","HPA_consensus_duodenum",
  "HPA_consensus_endometrium","HPA_consensus_epididymis","HPA_consensus_esophagus",
  "HPA_consensus_fallopian_tube","HPA_consensus_gallbladder","HPA_consensus_heart_muscle",
  "HPA_consensus_hippocampal_formation","HPA_consensus_hypothalamus","HPA_consensus_kidney",
  "HPA_consensus_liver","HPA_consensus_lung","HPA_consensus_lymph_node","HPA_consensus_midbrain",
  "HPA_consensus_ovary","HPA_consensus_pancreas","HPA_consensus_parathyroid_gland",
  "HPA_consensus_pituitary_gland","HPA_consensus_placenta","HPA_consensus_prostate",
  "HPA_consensus_rectum","HPA_consensus_retina","HPA_consensus_salivary_gland",
  "HPA_consensus_seminal_vesicle","HPA_consensus_skeletal_muscle","HPA_consensus_skin",
  "HPA_consensus_small_intestine","HPA_consensus_smooth_muscle","HPA_consensus_spinal_cord",
  "HPA_consensus_spleen","HPA_consensus_stomach","HPA_consensus_testis","HPA_consensus_thymus",
  "HPA_consensus_thyroid_gland","HPA_consensus_tongue","HPA_consensus_tonsil","HPA_consensus_urinary_bladder",
  "HPA_consensus_vagina","HPA_consensus_highly_expressed"
)

# ---------- función principal ----------
normalize_dbnsfp <- function(
    in_file, out_dir,
    numeric_idxs     = .dbnsfp_numeric_idxs,
    colnames_vec     = .dbnsfp_colnames,
    preview_n        = 100,
    dedup_tolower    = FALSE,
    extra_multi_cols = NULL
) {
  
  stopifnot(file.exists(in_file))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # 1) Lectura robusta:
  #   - Si data.table está disponible: fread(fill=TRUE) (más robusto + rápido)
  #   - Si no: read.delim (fallback)
  # ------------------------------------------------------------
  if (requireNamespace("data.table", quietly = TRUE)) {
    
    df <- data.table::fread(
      in_file,
      sep = "\t",
      header = TRUE,
      quote = "",
      data.table = FALSE,
      colClasses = "character",
      na.strings = c("", "NA", "."),
      fill = TRUE,
      showProgress = FALSE
    )
    
    # fread suele devolver nombres sintácticos si hay duplicados raros; mantenemos tal cual
    # (check.names no aplica aquí)
    
  } else {
    
    df <- utils::read.delim(
      in_file,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = "",
      colClasses = "character",       # <- CLAVE: no adivinar tipos
      na.strings = c("", "NA", ".")   # <- recomendable con dbNSFP
    )
  }
  
  if (!nrow(df)) stop("El fichero dbNSFP está vacío: ", in_file)
  
  # ------------------------------------------------------------
  # 2) Numeric para índices indicados (solo los que EXISTEN)
  #    (evita errores si numeric_idxs contiene índices > ncol(df))
  # ------------------------------------------------------------
  if (!is.null(numeric_idxs) && length(numeric_idxs)) {
    numeric_idxs <- numeric_idxs[is.finite(numeric_idxs)]
    numeric_idxs <- numeric_idxs[numeric_idxs >= 1 & numeric_idxs <= ncol(df)]
  }
  df <- .make_numeric_columns(df, numeric_idxs)
  
  # ------------------------------------------------------------
  # 3) Renombrado si coincide el nº de columnas
  # ------------------------------------------------------------
  if (!is.null(colnames_vec) && length(colnames_vec) == ncol(df)) {
    colnames(df) <- colnames_vec
  }
  
  # ------------------------------------------------------------
  # 4) Deduplicación de columnas multivalor separadas por ';'
  # ------------------------------------------------------------
  priority_multi <- c(
    "genename","Ensembl_transcriptid","Ensembl_proteinid",
    "SIFT_pred","SIFT4G_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred",
    "MutationTaster_pred","MutationAssessor_pred","PROVEAN_pred",
    "ClinVar_sig","KEGG_pathway",
    "GO_term_BP","GO_term_CC","GO_term_MF",
    "cosmic_tumor_site"
  )
  
  if (!is.null(extra_multi_cols)) {
    priority_multi <- unique(c(priority_multi, extra_multi_cols))
  }
  
  priority_multi <- intersect(priority_multi, names(df))
  
  is_charfac <- vapply(df, function(v) is.character(v) || is.factor(v), logical(1))
  cands <- names(df)[is_charfac]
  has_semi <- vapply(df[cands], function(v) any(grepl(";", v, fixed = TRUE), na.rm = TRUE), logical(1))
  auto_multi <- cands[has_semi]
  
  # proteger columnas clave que no queremos tocar
  protect <- intersect(c("SNP","rsid_final","chr","BP","ref","alt"), names(df))
  
  cols_multi <- setdiff(union(priority_multi, auto_multi), protect)
  if (length(cols_multi)) {
    df[cols_multi] <- lapply(df[cols_multi], dedup_semicol, tolower_vals = dedup_tolower)
  }
  
  # ------------------------------------------------------------
  # 5) Guardar
  # ------------------------------------------------------------
  out_csv <- file.path(out_dir, "dbnsfp_normalized.csv")
  out_rds <- file.path(out_dir, "dbnsfp_normalized.rds")
  
  utils::write.csv(df, out_csv, row.names = FALSE)
  saveRDS(df, out_rds)
  
  # ------------------------------------------------------------
  # 6) Preview
  # ------------------------------------------------------------
  prev <- utils::head(df, preview_n)
  
  list(
    csv     = out_csv,
    rds     = out_rds,
    preview = prev,
    n_rows  = nrow(df),
    n_cols  = ncol(df)
  )
}
