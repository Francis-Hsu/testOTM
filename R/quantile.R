#' Semi-discrete Optimal Transport Quantiles
#'
#' \code{tos.quantile} computes the optimal transport quantiles.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries, as well as the indices of the corresponding cells.
#' @keywords multivariate
#' @export
tos.quantile = function(object, query, scale = TRUE, ...) {
  UseMethod("tos.quantile")
}

#' 2D Semi-discrete Optimal Transport Quantiles
#'
#' The 2D implementation for \code{tos.quantile}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries, as well as the indices of the corresponding cells.
#' @keywords internal
#' @export
tos.quantile.tos.2d = function(object, query, scale = TRUE, ...) {
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD2D(query, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, ], Cell.Id = cell.id))
}

#' 3D Semi-discrete Optimal Transport Quantiles
#'
#' The 3D implementation for \code{tos.quantile}.
#' @param object a fitted 3D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries, as well as the indices of the corresponding cells.
#' @keywords internal
#' @export
tos.quantile.tos.3d = function(object, query, scale = TRUE, ...) {
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD3D(query, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, ], Cell.Id = cell.id))
}