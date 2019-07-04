#' Semi-discrete Optimal Transport Depth
#' 
#' Compute the MK depth relative to the multivariate uniform distribution U[0, 1]^d.
#' @param object a fitted optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @return a vector containing depth of the data.
#' @export
otm.depth = function(object, X, ...) {
  UseMethod("otm.depth")
}

#' 2D Semi-discrete Optimal Transport Depth
#' 
#' Compute the MK depth relative to the 2D standard uniform distribution.
#' @param object a fitted 2D optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @return a vector containing depth of the data.
#' @export
otm.depth.OTM_2D = function(object, X) {
  ranks = otm.rank.OTM_2D(object, X)
  otm.depth = depth.uniform(ranks)
  
  return(otm.depth)
}

#' Multivariate Uniform Depth
#' 
#' Compute the depth associated with U[0, 1]^d.
#' @param X a \eqn{n} by \eqn{d} numeric matrix.
#' @return a vector containing depth of the data.
#' @export
depth.uniform = function(X) {
  depth = 0.5 - apply(abs(t(t(X) - rep(1, ncol(X)) / 2)), 1, max)
  
  return(depth)
}