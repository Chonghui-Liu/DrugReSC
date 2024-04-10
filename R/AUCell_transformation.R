#' Title
#'
#' @param sc_matrix Single-cell transcriptome expression profile
#' @param up_signature_list The list of Gene sets up-regulated by drug instances
#' @param down_signature_list The list of Gene sets down-regulated by drug instances
#'
#' @return The drug-by-cell matrix
#' @export
#'
#' @importFrom AUCell AUCell_run
#' @importFrom BiocParallel MulticoreParam
#' @importFrom methods as
AUCell_transformation <- function(sc_matrix, up_signature_list, down_signature_list){
  print("Runing the  AUCell")
  D2C_matrix_up <- AUCell_run(as(sc_matrix,"dgCMatrix"), up_signature_list)@assays@data@listData[["AUC"]]
  D2C_matrix_down <- AUCell_run(as(sc_matrix,"dgCMatrix"), down_signature_list)@assays@data@listData[["AUC"]]
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
