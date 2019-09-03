#' Semi-discrete Optimal Transport Potential
#'
#' \code{otm.potential} computes the Alexandrov potential function.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @keywords multivariate
#' @export
otm.potential = function(object, query, scale = TRUE, ...) {
  UseMethod("otm.potential")
}

#' 2D Semi-discrete Optimal Transport Potential
#'
#' The 2D implementation of \code{otm.potential}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @keywords internal
#' @export
otm.potential.otm.2d = function(object, query, scale = TRUE, ...) {
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD2D(query, object$Data, object$Weight)
  potential = rowSums(query * object$Data[cell.id, , drop = F]) + object$Height[cell.id]
  
  return(potential)
}