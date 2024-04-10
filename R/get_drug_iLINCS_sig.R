#' Title
#'
#' @param true_drug_sig_id_list ID of the drug instance in the iLINCS database
#'
#' @return Results of drug differences analysis extracted from the iLINCS database
#' @export
#'
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom utils read.table

get_drug_iLINCS_sig <- function(true_drug_sig_id_list){
  iLINCS_signature_list <- list()
  for (i in 1:length(true_drug_sig_id_list)) {
    ilincs_signatureId <- true_drug_sig_id_list[i]
    if(is.na(ilincs_signatureId)){
      iLINCS_signature_list[[i]] <- NA
      names(iLINCS_signature_list)[i] <- true_drug_sig_id_list[i]
    }else{
      req <- httr::POST("http://www.ilincs.org/api/ilincsR/downloadSignature", body = list(sigID = paste(ilincs_signatureId), display = FALSE), encode = "json")
      ilincs_sessionId<-unlist(httr::content(req))
      fileUrl=paste("http://www.ilincs.org/tmp/",ilincs_sessionId,".xls",sep="")
      signatureData<-read.table(fileUrl,sep="\t",header=T,stringsAsFactors = F)
      iLINCS_signature_list[[i]] <- signatureData
      names(iLINCS_signature_list)[i] <- true_drug_sig_id_list[i]
    }
  }
  names(iLINCS_signature_list) <- true_drug_sig_id_list
  return(iLINCS_signature_list)
}
