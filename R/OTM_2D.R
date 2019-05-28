#' 2D semi-continuous optimal transport map
#' 
#' Computes the semi-continuous optimal transport map from a uniform measure on [0, 1]^2 to a given data set.
#' @param data Input coordinate matrix, of size \eqn{n} by \eqn{2}.
#' @param epsilon Convergence threshold for optimization.
#' @param maxit Max number of iterations before termination.
#' @param verbose Wether to display messages during optimization.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return A list describing the resulting Laguerre diagram.
#' @keywords Voronoi, Laguerre
#' @export
OTM_2D = function(data, epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  if (!is.matrix(data) || ncol(data) < 2) {
    stop("Data must be a matrix with ncol >= 2.")
  }
  
  if (na.rm) {
    data = data[complete.cases(data), ]
  }
  
  otm_list = powerDiag2D(data, epsilon, maxit, verbose)
  colnames(otm_list$Vert) = c("vert.x", "vert.y", "cell")
  colnames(otm_list$Edge) = c("source.x", "source.y", "target.x", "target.y")
  colnames(otm_list$Centroids) = c("cent.x", "cent.y")
  
  class(otm_list) = "OTM_2D"
  
  return(otm_list)
}