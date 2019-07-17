#' Semi-discrete Optimal Transport Depth
#'
#' Compute the optimal transport depth relative to the multivariate uniform distribution U[0, 1]^d.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a vector containing the depths of the data.
#' @export
otm.depth = function(object, Q, ...) {
  UseMethod("otm.depth")
}

#' 2D Semi-discrete Optimal Transport Depth
#'
#' Compute the optimal transport depth relative to the 2D standard uniform distribution.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param use.geo logical indicating if the geometric method should be used to compute the ranks.
#' @return a vector containing the depths of the data.
#' @export
otm.depth.OTM.2D = function(object, Q, use.geo = FALSE) {
  ranks = otm.rank.OTM.2D(object, Q, use.geo)$Rank
  otm.depth = depth.uniform(ranks)
  
  return(otm.depth)
}

#' Multivariate Uniform Depth Function
#' 
#' Compute the depths associated with U[0, 1]^d.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a vector containing the depths of the data.
#' @export
depth.uniform = function(Q) {
  depth = 0.5 - apply(abs(t(t(Q) - rep(0.5, ncol(Q)))), 1, max)
  
  return(depth)
}