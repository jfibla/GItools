get_clinical_metric_sets <- function() {
  list(
    "Default" = c(
      "SIFT_score","SIFT4G_score",
      "Polyphen2_HDIV_score","Polyphen2_HVAR_score",
      "MutationTaster_score","MutationAssessor_score",
      "PROVEAN_score","REVEL_score",
      "GERP_91_mammals","phyloP17way_primate","phastCons17way_primate",
      "SIFT_converted_rankscore","SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_rankscore","Polyphen2_HVAR_rankscore",
      "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
      "REVEL_rankscore","GERP_91_mammals_rankscore",
      "phyloP17way_primate_rankscore","phastCons17way_primate_rankscore"
    ),
   
    "Pathogenicity" = c(
      "SIFT_score","Polyphen2_HDIV_score",
      "Polyphen2_HVAR_score","MutationTaster_score",
      "MutationAssessor_score",
      "PROVEAN_score","VEST4_score","MetaSVM_score",
      "MetaLR_score","MetaRNN_score",
      "M_CAP_score","REVEL_score",
      "MutPred_score","MVP_score","gMVP_score","MPC_score",
      "PrimateAI_score","DEOGEN2_score",
      "BayesDel_addAF_score",
      "BayesDel_noAF_score",
      "ClinPred_score",
      "LIST_S2_score",
      "VARITY_R_score","VARITY_ER_score",
      "ESM1b_score","AlphaMissense_score",
      "PHACTboost_score","MutFormer_score",
      "MutScore_score"
    ),
    
    "Conservation" = c(
      "GERP++_RS","phyloP100way_vertebrate",
      "phyloP470way_mammalian","phyloP17way_primate",
      "phastCons100way_vertebrate","phastCons470way_mammalian",
      "phastCons17way_primate"
    ),
    
    "Functional impact" = c(
      "CADD_raw","CADD_phred",
      "DANN_score","Eigen_raw_coding","Eigen_phred_coding"
    ),
    
    "LoF / Haploinsufficiency" = c(
      "P_HI","GHIS","ClinGen_Haploinsufficiency_Score",
      "ExAC_del_score","ExAC_dup_score","ExAC_cnv_score"
    ),
    
    "Gene damage" = c(
      "GDI","GDI_Phred",
      "LoFtool_score","Gene_indispensability_score"
    )
  )
}
