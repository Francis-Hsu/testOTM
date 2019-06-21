#' 2D semi-continuous optimal transport map
#' 
#' Computes the semi-continuous optimal transport map from a uniform measure on [0, 1]^2 to a given data set.
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
  
  otm_list = dualGraphs2D(data, epsilon, maxit, verbose)
  colnames(otm_list$Centroid) = c("x", "y")
  colnames(otm_list$Vertex.RDT) = c("cell", "x", "y")
  colnames(otm_list$Vertex.RVD) = c("cell", "x", "y")
  
  otm_list$Vertex.RDT = as.data.frame(otm_list$Vertex.RDT)
  otm_list$Vertex.RVD = as.data.frame(otm_list$Vertex.RVD)
  
  class(otm_list) = "OTM_2D"
  
  return(otm_list)
}