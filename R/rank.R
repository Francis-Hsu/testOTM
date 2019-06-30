#' Semi-discrete optimal transport rank
#' 
#' Compute the rank function based on optimal transport
#' Semi-discrete optimal transport rank.
#' @param object a fitted optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @param h step size for the numerical differentiation.
#' @return a matrix containing the rank of the data.
#' @export
otm.rank = function(object, X, h, ...) {
  UseMethod("otm.rank")
}

#' 2D semi-discrete optimal transport rank
#' 
#' Compute the rank function based on optimal transport
#' 2D semi-discrete optimal transport rank.
#' @param object a fitted 2D optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @param h step size for the numerical differentiation.
#' @return a matrix containing the rank of the data.
#' @export
otm.rank.OTM_2D = function(object, X, h = 1e-7) {
  n = nrow(X)
  d = 2
  
  G = matrix(0, n, d)
  for (i in 1:d) {
    Xp = X
    Xn = X
    Xp[, i] = Xp[, i] + h
    Xn[, i] = Xn[, i] - h
    G[, i] = (otm.potential(object, Xp) - otm.potential(object, Xn)) / (2 * h)
  }
  
  return(G)
}