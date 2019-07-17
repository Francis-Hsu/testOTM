#' Semi-discrete Optimal Transport Quantile
#'
#' Compute the optimal transport quantiles.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a list containing the quantiles of the data, and the indices of the corresponding cells.
#' @export
otm.quantile = function(object, Q, ...) {
  UseMethod("otm.quantile")
}

#' 2D Semi-discrete Optimal Transport Quantile
#'
#' Compute the 2D optimal transport quantiles.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a list containing the quantiles of the data, and the indices of the corresponding cells.
#' @export
otm.quantile.OTM.2D = function(object, Q) {
  cell.id = locateRVD2D(Q, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, ], Cell.Id = cell.id))
}