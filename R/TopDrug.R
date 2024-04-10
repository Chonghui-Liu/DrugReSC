#' Title
#'
#' @param drug_scores The normalized drug scores
#'
#' @return Ranked list of candidate drugs
#' @export
#'@importFrom httr GET
#'@importFrom httr content
#'@importFrom jsonlite fromJSON
#'

TopDrug <- function(drug_scores){
  predicted_table <- data.frame()
  tmp_imp_score <- drug_scores
  treatment <- NULL
  drug_info <- NULL
  for (j in 1:length(names(tmp_imp_score))) {
    ilincs_signatureId <- names(tmp_imp_score)[j]
    apiUrl <- paste("http://www.ilincs.org/api/SignatureMeta/",ilincs_signatureId,sep="")
    req <- GET(apiUrl)
    ilincsJSON<-content(req,type="text")
    #prettify(ilincsJSON)
    ilincsSigMetaData<-fromJSON(ilincsJSON)
    treatment[j] <- ilincsSigMetaData$lincsSigID
    drug_info[j] <- paste0(ilincsSigMetaData[["compound"]],"_",ilincsSigMetaData[["cellline"]],"_",ilincsSigMetaData[["concentration"]],"_",ilincsSigMetaData[["time"]])
  }
  tmp_rank <- rank(-drug_scores)
  tmp_predicted_table <- predicted_table <- data.frame(rank=tmp_rank,L1000_id=treatment,drug_infor = drug_info,scores=drug_scores)
  predicted_table <- tmp_predicted_table[order(tmp_predicted_table$rank), ]
  return(predicted_table)
}
