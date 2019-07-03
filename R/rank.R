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
  
  otm.ranks = matrix(0, n, d)
  
  # compute the ranks geometrically
  location.id = locateRDT2D(X, as.matrix(object$Vertex.RDT[, 2:3]))
  inside.id = (1:n)[location.id > 0]
  outside.id = (1:n)[location.id <= 0]
  for (i in inside.id) {
    # get cells that are dual to the triangle we found
    cell.id = unlist(subset(object$Vertex.RDT, cell == location.id[i], select = id))
    
    # we round to the 8th digit for the vertices comparison
    rvd.vert.freq = as.data.frame(table(round(subset(object$Vertex.RVD, 
                                                     cell %in% cell.id, select = c("x", "y")), 8)), 
                                  stringsAsFactors = F)
    
    # we search for the common vertex shared by all three cells from the RVD
    # there should always be such a vertex since we located the query point
    # on the dual (RDT) of RVD
    otm.ranks[i, ] = as.numeric(rvd.vert.freq[match(3, rvd.vert.freq$Freq), 1:2])
  }
  otm.ranks[outside.id, ] = object$Centroid[locateRVD2D(X[outside.id, ], object$Data, object$Weight), ]
  
  # compute the ranks numerically
  # for (j in 1:d) {
  #   Xp = X[outside.id, ]
  #   Xn = X[outside.id, ]
  #   Xp[, j] = Xp[, j] + h
  #   Xn[, j] = Xn[, j] - h
  #   otm.ranks[outside.id, j] = (otm.potential(object, Xp) - otm.potential(object, Xn)) / (2 * h)
  # }
  
  return(otm.ranks)
}