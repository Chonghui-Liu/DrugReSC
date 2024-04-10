#' Title
#'
#' @param iLINCS_signature_list Results of drug differences analysis extracted from the iLINCS database
#' @param pvalue_cutoff Thresholds for p values
#' @param logFC_cutoff  Thresholds for logFC
#'
#' @return signature list of drug instances
#' @export
#'

get_drug_sig <- function(iLINCS_signature_list, pvalue_cutoff=0.01, logFC_cutoff=0){
  drug_signature_list <- list()
  for(i in 1:length(iLINCS_signature_list)){
    up_signature <- list(iLINCS_signature_list[[i]][iLINCS_signature_list[[i]]$Significance_pvalue < pvalue_cutoff&iLINCS_signature_list[[i]]$Value_LogDiffExp  > logFC_cutoff,]$Name_GeneSymbol)
    down_signature <- list(iLINCS_signature_list[[i]][iLINCS_signature_list[[i]]$Significance_pvalue < pvalue_cutoff&iLINCS_signature_list[[i]]$Value_LogDiffExp  < -logFC_cutoff,]$Name_GeneSymbol)
    drug_signature_list[[i]] <- c(up=up_signature, down=down_signature)
    names(drug_signature_list)[i] <- names(iLINCS_signature_list)[i]
  }
  names(drug_signature_list) <- names(iLINCS_signature_list)
  return(drug_signature_list)
}
