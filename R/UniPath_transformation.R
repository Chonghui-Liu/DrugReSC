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
#' @importFrom UniPath binorm
#' @importFrom UniPath combine
#' @importFrom UniPath adjust
#' @importFrom utils data
#'


UniPath_transformation <- function(sc_matrix, up_signature_list, down_signature_list,ncore){
  print("Runing the  UniPath")
  human_null_data <- UniPath::human_null_data
  Pval = binorm(human_null_data)
  Pval1 = binorm(sc_matrix)
  max_length_up <- max(sapply(up_signature_list, length))
  up_signature_list_padded <- lapply(up_signature_list, function(x) c(x, rep(NA, max_length_up - length(x))))
  up_df <- as.data.frame(do.call(rbind, up_signature_list_padded));rownames(up_df) <- names(up_signature_list)
  up_combp_ref = combine(up_df,human_null_data,rownames(human_null_data),Pval,thr=ncore)
  up_combp = combine(up_df,sc_matrix,rownames(sc_matrix),Pval1,thr=ncore)
  D2C_matrix_up <- adjust(up_combp,up_combp_ref)$adjpvaraw

  max_length_down <- max(sapply(down_signature_list, length))
  down_signature_list_padded <- lapply(down_signature_list, function(x) c(x, rep(NA, max_length_down - length(x))))
  down_df <- as.data.frame(do.call(rbind, down_signature_list_padded));rownames(down_df) <- names(down_signature_list)
  down_combp_ref = combine(down_df,human_null_data,rownames(human_null_data),Pval,thr=ncore)
  down_combp = combine(down_df,sc_matrix,rownames(sc_matrix),Pval1,thr=ncore)
  D2C_matrix_down <- adjust(down_combp,down_combp_ref)$adjpvaraw

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
