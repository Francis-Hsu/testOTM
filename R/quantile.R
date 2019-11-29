#' Semi-discrete Optimal Transport Quantiles
#'
#' \code{tos.quantile} computes the optimal transport quantiles.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 or 3 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries and the indices of the corresponding cells.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the quantile of point (0, 0).
#' tos.quantile(XY.OTM, c(0, 0))
#' @keywords multivariate
#' @references Promit Ghosal and Bodhisattva Sen (2019). 
#' \emph{Multivariate Ranks and Quantiles Using Optimal Transportation and Applications to Goodness-of-Fit Testing}.
#' \url{http://arxiv.org/abs/1905.05340}.
#' @export
tos.quantile = function(object, query, scale = TRUE, ...) {
  UseMethod("tos.quantile")
}

#' 2D Semi-discrete Optimal Transport Quantiles
#'
#' The 2D implementation for \code{tos.quantile}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 2 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries and the indices of the corresponding cells.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # compute the quantile of point (0, 0).
#' tos.quantile(XY.OTM, c(0, 0))
#' @keywords internal
#' @export
tos.quantile.tos.2d = function(object, query, scale = TRUE, ...) {
  # input validation
  if (is.vector(query) && length(query) == 2) {
    query = t(matrix(query))
  }
  if (!is.matrix(query) || NCOL(query) != 2) {
    stop("Input query must be a matrix with ncol = 2.")
  }
  
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD2D(query, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, , drop = FALSE], Cell.Id = cell.id))
}

#' 3D Semi-discrete Optimal Transport Quantiles
#'
#' The 3D implementation for \code{tos.quantile}.
#' @param object a fitted 3D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point
#' (or a vector of length 3 representing a single point).
#' @param scale logical indicating if the queries should be scaled.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the quantiles of the queries and the indices of the corresponding cells.
#' @examples 
#' # generate some data
#' X = c(2.3, 0.8, 1.2, 0.6, 0.2, 0.5, -0.7, -0.6, -0.2, 3.2)
#' Y = c(1.5, -0.8, 0.6, -2.1, 1.3, 1.3, -1.1, 1.0, -0.2, 1.0)
#' Z = c(0.8, -0.3, -1.3, -1.8, 1.2, -0.7, -0.8, 1.3, -0.9, 2.0)
#' 
#' # compute the optimal transport map from U[0, 1]^3 to the data
#' XYZ.OTM = tos.fit(cbind(X, Y, Z))
#' 
#' # compute the quantile of point (0, 0, 0).
#' tos.quantile(XYZ.OTM, c(0, 0, 0))
#' @keywords internal
#' @export
tos.quantile.tos.3d = function(object, query, scale = TRUE, ...) {
  # input validation
  if (is.vector(query) && length(query) == 3) {
    query = t(matrix(query))
  }
  if (!is.matrix(query) || NCOL(query) != 3) {
    stop("Input query must be a matrix with ncol = 3.")
  }
  
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  cell.id = locateRVD3D(query, object$Data, object$Weight)
  
  return(list(Quantiles = object$Data[cell.id, , drop = FALSE], Cell.Id = cell.id))
}