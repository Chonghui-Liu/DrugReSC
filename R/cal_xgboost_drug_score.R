#' Title
#'
#' @param D2C_matrix The drug-by-cell matrix generated by DrugReSC
#' @param cell_labels The cell phenotype labels
#'
#' @return The raw drug scores
#' @export
#'@importFrom xgboost xgboost
#'@importFrom xgboost xgb.importance


cal_xgboost_drug_score <- function(D2C_matrix,cell_labels){
  print("Running xgboost")
  zero_rows <- which(apply(D2C_matrix, 1, function(x) all(x == 0)))
  if(length(zero_rows) != 0){
    D2C_matrix <- D2C_matrix[-zero_rows,]
  }
  tmp_matrix <- as.matrix(t(D2C_matrix))
  xgboost_model <- xgboost(data = tmp_matrix, label = cell_labels, nthread = 2, nrounds=500,objective = "binary:logistic", verbose=0)
  sc_predict <- xgb.importance(feature_names = colnames(tmp_matrix), model = xgboost_model)$Gain
  names(sc_predict) <- xgb.importance(feature_names = colnames(tmp_matrix), model = xgboost_model)$Feature
  sorted_sc_predict <- numeric(length(colnames(tmp_matrix)))
  for (i in seq_along(colnames(tmp_matrix))) {
    colname <- colnames(tmp_matrix)[i]
    if (colname %in% names(sc_predict)) {
      sorted_sc_predict[i] <- sc_predict[colname]
    } else {
      sorted_sc_predict[i] <- 0
    }
  }
  return(sorted_sc_predict)
}
