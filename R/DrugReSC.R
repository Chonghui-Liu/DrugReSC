#' Title
#'
#' @param sc_matrix Single-cell transcriptome expression profile
#' @param drug_signature_list signature list of drug instances
#' @param PACSI_result Results of the PACSI method
#' @param D2C_transformation Methods for transforming gene-by-cell matrices into drug-by-cell matrices
#' @param improtance_method Feature importance calculation method
#' @param ncore The number of parallel threads
#'
#' @return Ranked list of candidate drugs
#' @export
#'

DrugReSC <- function(sc_matrix,drug_signature_list,PACSI_result,D2C_transformation="ssgsea",improtance_method = "rf",ncore=8){

  print("Extraction of up-regulated and down-regulated signature")

  D2C_matrix_final <-list()
  up_signature_list <- list()
  down_signature_list <- list()
  for (i in 1:length(drug_signature_list)) {
    up_signature_list[[i]] <- drug_signature_list[[i]]$up
    down_signature_list[[i]] <- drug_signature_list[[i]]$down
  }
  names(up_signature_list) <- names(drug_signature_list)
  names(down_signature_list) <- names(drug_signature_list)

  print("Construction of the D2C matrix")

  if(D2C_transformation=="ssgsea"){
    D2C_matrix_final <- ssGSEA_transformation(sc_matrix, up_signature_list, down_signature_list,ncore=ncore)
  }else if(D2C_transformation=="AUCell"){
    D2C_matrix_final <- AUCell_transformation(sc_matrix, up_signature_list, down_signature_list)
  }else if(D2C_transformation=="UniPath"){
    D2C_matrix_final <- UniPath_transformation(sc_matrix, up_signature_list, down_signature_list,ncore=ncore)
  }else{
    stop("Incorrect 'D2C_transformation' parameter setting")
  }

  print("Calculation of drug scores")

  cell_labels <- PACSI_result[["result"]]$cell_phenotype_labels
  if(improtance_method == "lgr" ) {
    raw_drug_scores <- cal_lgr_drug_score(D2C_matrix_final,cell_labels)
  } else if( improtance_method == "svm" ) {
    raw_drug_scores <- cal_svm_drug_score(D2C_matrix_final,cell_labels)
  } else if( improtance_method == "aov" ) {
    raw_drug_scores <- cal_aov_drug_score(D2C_matrix_final,cell_labels)
  } else if( improtance_method == "rf" ) {
    raw_drug_scores <- cal_rf_drug_score(D2C_matrix_final,cell_labels,ncore=ncore)
  } else if( improtance_method == "wilcoxon" ) {
    raw_drug_scores <- cal_wilcoxon_drug_score(D2C_matrix_final,cell_labels)
  } else if(improtance_method == "xgboost"){
    raw_drug_scores <- cal_xgboost_drug_score(D2C_matrix_final,cell_labels)
  }else {
    stop("Incorrect 'improtance_method' parameter setting")
  }

  print("normalization and output")

  drug_scores <- min_max_normalize(raw_drug_scores)
  output <- TopDrug(drug_scores)
  return(output)
}
