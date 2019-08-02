#' Semi-discrete Optimal Transport Potential
#'
#' \code{otm.potential} computes the Alexandrov potential function.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @keywords multivariate
#' @export
otm.potential = function(object, Q, ...) {
  UseMethod("otm.potential")
}

#' 2D Semi-discrete Optimal Transport Potential
#'
#' The 2D implementation of \code{otm.potential}.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a matrix containing the potentials of the queries.
#' @keywords internal
#' @export
otm.potential.otm.2d = function(object, Q) {
  cell.id = locateRVD2D(Q, object$Data, object$Weight)
  potential = rowSums(Q * object$Data[cell.id, , drop = F]) + object$Height[cell.id]
  
  return(potential)
}