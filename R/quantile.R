#' Semi-discrete Optimal Transport Quantile
#'
#' \code{otm.quantile} computes the optimal transport quantiles.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries, as well as the indices of the corresponding cells.
#' @keywords multivariate
#' @export
otm.quantile = function(object, Q, ...) {
  UseMethod("otm.quantile")
}

#' 2D Semi-discrete Optimal Transport Quantile
#'
#' The 2D implementation for \code{otm.quantile}.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries, as well as the indices of the corresponding cells.
#' @keywords internal
#' @export
otm.quantile.otm.2d = function(object, Q, ...) {
  cell.id = locateRVD2D(Q, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, ], Cell.Id = cell.id))
}