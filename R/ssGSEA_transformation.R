#' Title
#'
#' @param sc_matrix Single-cell transcriptome expression profile
#' @param up_signature_list The list of Gene sets up-regulated by drug instances
#' @param down_signature_list The list of Gene sets down-regulated by drug instances
#' @param ncore The number of parallel threads
#'
#' @return The drug-by-cell matrix
#' @export
#'
#' @importFrom GSVA gsva
#'
ssGSEA_transformation <- function(sc_matrix, up_signature_list, down_signature_list,ncore){
  print("Runing the  ssGSEA")
  D2C_matrix_up <- gsva(as.matrix(sc_matrix), up_signature_list, method='ssgsea',ssgsea.norm = TRUE,kcdf='Gaussian',parallel.sz=ncore,verbose=F)
  D2C_matrix_down <- gsva(as.matrix(sc_matrix), down_signature_list, method='ssgsea',ssgsea.norm = TRUE,kcdf='Gaussian',parallel.sz=ncore,verbose=F)
  if(all(rownames(D2C_matrix_up) == rownames(D2C_matrix_down)) & length(rownames(D2C_matrix_up)) == length(up_signature_list)){
    D2C_matrix_final <- D2C_matrix_up-D2C_matrix_down
  }else{
    D2C_matrix_up_final <- matrix(NA,ncol = ncol(D2C_matrix_up),nrow = length(up_signature_list))
    for (i in 1:length(names(up_signature_list))) {
      if(names(up_signature_list)[i] %in% rownames(D2C_matrix_up)){
        D2C_matrix_up_final[i,] <- D2C_matrix_up[match(names(up_signature_list)[i],rownames(D2C_matrix_up)),]
      }else{
        D2C_matrix_up_final[i,] <- c(rep(0,ncol(D2C_matrix_up)))
      }
    }
    rownames(D2C_matrix_up_final) = names(up_signature_list)
    D2C_matrix_down_final <- matrix(NA,ncol = ncol(D2C_matrix_down),nrow = length(down_signature_list))
    for (i in 1:length(names(down_signature_list))) {
      if(names(down_signature_list)[i] %in% rownames(D2C_matrix_down)){
        D2C_matrix_down_final[i,] <- D2C_matrix_down[match(names(down_signature_list)[i],rownames(D2C_matrix_down)),]
      }else{
        D2C_matrix_down_final[i,] <- c(rep(0,ncol(D2C_matrix_down)))
      }
    }
    rownames(D2C_matrix_down_final) = names(down_signature_list)
    D2C_matrix_final <- D2C_matrix_up_final-D2C_matrix_down_final
  }
  names(D2C_matrix_final) <- names(up_signature_list)
  return(D2C_matrix_final)
}
