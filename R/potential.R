#' Semi-discrete Optimal Transport Potential
#'
#' \code{tos.potential} computes the Alexandrov potential function.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the potential of point (0, 0)
#' tos.potential(XY.OTM, c(0, 0))
#' @keywords multivariate
#' @references Xianfeng Gu, Feng Luo, Jian Sun, and Shing-Tung Yau (2015).
#' \emph{Variational Principles for Minkowski Type Problems, Discrete Optimal Transport, and Discrete Monge-Ampere Equations}.
#' Asian Journal of Mathematics 20 (January).
#' \url{https://doi.org/10.4310/AJM.2016.v20.n2.a7}.
#' @export
tos.potential = function(object, query, scale = TRUE, ...) {
  UseMethod("tos.potential")
}

#' 2D Semi-discrete Optimal Transport Potential
#'
#' The 2D implementation of \code{tos.potential}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the potential of point (0, 0)
#' tos.potential(XY.OTM, c(0, 0))
#' @keywords internal
#' @export
tos.potential.tos.2d = function(object, query, scale = TRUE, ...) {
  # input validation
  if (is.vector(query) && length(query) == 2) {
    query = t(matrix(query))
  }
  
  if (!is.matrix(query) || NCOL(query) != 2) {
    stop("Input query must be a matrix with ncol = 2.")
  }
  
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD2D(query, object$Data, object$Weight)
  potential = rowSums(query * object$Data[cell.id, , drop = FALSE]) + object$Height[cell.id]
  
  return(potential)
}

#' 3D Semi-discrete Optimal Transport Potential
#'
#' The 3D implementation of \code{tos.potential}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 3 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a matrix containing the potentials of the queries.
#' @examples 
#' # generate some data
#' X = c(2.3, 0.8, 1.2, 0.6, 0.2, 0.5, -0.7, -0.6, -0.2, 3.2)
#' Y = c(1.5, -0.8, 0.6, -2.1, 1.3, 1.3, -1.1, 1.0, -0.2, 1.0)
#' Z = c(0.8, -0.3, -1.3, -1.8, 1.2, -0.7, -0.8, 1.3, -0.9, 2.0)
#' 
#' # compute the optimal transport map from U[0, 1]^3 to the data
#' XYZ.OTM = tos.fit(cbind(X, Y, Z))
#' 
#' # compute the potential of point (0, 0, 0)
#' tos.potential(XYZ.OTM, c(0, 0, 0))
#' @keywords internal
#' @export
tos.potential.tos.3d = function(object, query, scale = TRUE, ...) {
  # input validation
  if (is.vector(query) && length(query) == 3) {
    query = t(matrix(query))
  }
  
  if (!is.matrix(query) || NCOL(query) != 3) {
    stop("Input query must be a matrix with ncol = 3.")
  }
  
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD3D(query, object$Data, object$Weight)
  potential = rowSums(query * object$Data[cell.id, , drop = FALSE]) + object$Height[cell.id]
  
  return(potential)
}