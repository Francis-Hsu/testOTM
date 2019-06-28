#' 2D semi-continuous optimal transport map
#' 
#' Computes the semi-discrete optimal transport map from a uniform measure on [0, 1]^2 to a given data set.
#' @param data input coordinate matrix, of size \eqn{n} by \eqn{2}.
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether NA values should be stripped before the computation proceeds.
#' @return a list describing the resulting Laguerre diagram.
#' @keywords Voronoi, Laguerre
#' @export
OTM_2D = function(data, epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  if (!is.matrix(data) || ncol(data) < 2) {
    stop("Data must be a matrix with ncol >= 2.")
  }
  
  if (na.rm) {
    data = data[complete.cases(data), ]
  }
  
  object = dualGraphs2D(data, epsilon, maxit, verbose)
  colnames(object$Centroid) = c("x", "y")
  colnames(object$Vertex.RDT) = c("cell", "x", "y")
  colnames(object$Vertex.RVD) = c("cell", "x", "y")
  
  object$Vertex.RDT = as.data.frame(object$Vertex.RDT)
  object$Vertex.RVD = as.data.frame(object$Vertex.RVD)
  object$Height = -(rowSums(object$Data^2) + object$Weight) / 2 # use for computing Alexandrov's potential
  
  class(object) = "OTM_2D"
  
  return(object)
}