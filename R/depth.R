#' Semi-discrete Optimal Transport Depth
#'
#' \code{otm.depth} computes the optimal transport depth relative to the multivariate uniform distribution \eqn{U[0, 1]^d}.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' @param \dots additional arguments.
#' @return a vector containing the depths of the data.
#' @keywords multivariate
#' @export
otm.depth = function(object, query, scale = TRUE, rank.data = "uniform", ...) {
  UseMethod("otm.depth")
}

#' 2D Semi-discrete Optimal Transport Depth
#'
#' The 2D implementation of \code{otm.depth}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' @param \dots additional arguments, currently without effect.
#' @return a vector containing the depths of the data.
#' @keywords internal
otm.depth.otm.2d = function(object, query, scale = TRUE, rank.data = "uniform", ...) {
  ranks = otm.rank.otm.2d(object, query, scale = scale, rank.data = rank.data, rank.algo = "lp")$Rank
  otm.depth = depth.uniform(ranks)
  
  return(otm.depth)
}

#' Multivariate Uniform Depth Function
#' 
#' \code{otm.depth} computes the depths associated with the \eqn{U[0, 1]^d} measure.
#' @param query a numeric matrix where each row represents a query point.
#' @return a vector containing the depths of the data.
#' @keywords internal
depth.uniform = function(query) {
  depth = 0.5 - apply(abs(t(t(query) - rep(0.5, ncol(query)))), 1, max)
  
  return(depth)
}