#' Title
#'
#' @param x The raw drug score vector
#'
#' @return The normalized drug score vectors
#' @export
#'

min_max_normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
