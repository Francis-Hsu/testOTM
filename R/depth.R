#' Semi-discrete Optimal Transport Depth
#'
#' \code{tos.depth} computes the optimal transport depth relative to the multivariate uniform distribution \eqn{U[0, 1]^d}.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 or 3 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' @param \dots additional arguments.
#' @return a vector containing depths of the queries.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the depth of point (0, 0).
#' tos.depth(XY.OTM, c(0, 0))
#' @keywords multivariate
#' @references Victor Chernozhukov, Alfred Galichon, Marc Hallin, and Marc Henry (2017).
#' \emph{Monge–Kantorovich Depth, Quantiles, Ranks and Signs}.
#' The Annals of Statistics 45 (1): 223–56.
#' \url{https://doi.org/10.1214/16-AOS1450}.
#' @export
tos.depth = function(object, query, scale = TRUE, rank.data = "uniform", ...) {
  UseMethod("tos.depth")
}

#' 2D Semi-discrete Optimal Transport Depth
#'
#' The 2D implementation of \code{tos.depth}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 representing a single point).
#' @param rank.data choose the method for assigning ranks to the data points. 
#' @param \dots additional arguments, currently without effect.
#' @return a vector containing depths of the queries.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the depth of point (0, 0).
#' tos.depth(XY.OTM, c(0, 0))
#' @keywords internal
#' @export
tos.depth.tos.2d = function(object, query, scale = TRUE, rank.data = "uniform", ...) {
  ranks = tos.rank.tos.2d(object, query, scale = scale, rank.data = rank.data, rank.algo = "lp")$Rank
  tos.depth = depth.uniform(ranks)
  
  return(tos.depth)
}

#' 3D Semi-discrete Optimal Transport Depth
#'
#' The 3D implementation of \code{tos.depth}.
#' @param object a fitted 3D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 3 representing a single point).
#' @param rank.data choose the method for assigning ranks to the data points. 
#' @param \dots additional arguments, currently without effect.
#' @return a vector containing the depths of the queries.
#' @examples 
#' # generate some data
#' X = c(2.3, 0.8, 1.2, 0.6, 0.2, 0.5, -0.7, -0.6, -0.2, 3.2)
#' Y = c(1.5, -0.8, 0.6, -2.1, 1.3, 1.3, -1.1, 1.0, -0.2, 1.0)
#' Z = c(0.8, -0.3, -1.3, -1.8, 1.2, -0.7, -0.8, 1.3, -0.9, 2.0)
#' 
#' # compute the optimal transport map from U[0, 1]^3 to the data
#' XYZ.OTM = tos.fit(cbind(X, Y, Z))
#' 
#' # compute the depth of point (0, 0, 0).
#' tos.depth(XYZ.OTM, c(0, 0, 0))
#' @keywords internal
#' @export
tos.depth.tos.3d = function(object, query, scale = TRUE, rank.data = "center", ...) {
  ranks = tos.rank.tos.3d(object, query, scale = scale, rank.data = rank.data)$Rank
  tos.depth = depth.uniform(ranks)
  
  return(tos.depth)
}

#' Multivariate Uniform Depth Function
#' 
#' \code{depth.uniform} computes the depths associated with the \eqn{U[0, 1]^d} measure.
#' @param query a numeric matrix where each row represents a query point.
#' @return a vector containing the depths of the queries.
#' @keywords internal
depth.uniform = function(query) {
  depth = 0.5 - apply(abs(t(t(query) - rep(0.5, NCOL(query)))), 1, max)
  
  return(depth)
}