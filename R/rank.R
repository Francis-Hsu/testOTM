#' Semi-discrete Optimal Transport Rank
#' 
#' Compute the ranks based on optimal transport.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a matrix containing the ranks of the data.
#' @export
otm.rank = function(object, Q, ...) {
  UseMethod("otm.rank")
}

#' 2D Semi-discrete Optimal Transport Rank
#' 
#' Compute the 2D ranks based on optimal transport.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @return a matrix containing the ranks of the data.
#' @export
otm.rank.OTM_2D = function(object, Q) {
  n = nrow(Q)
  d = 2
  
  # container for the result
  otm.ranks = matrix(0, n, d)
  
  # compute the ranks geometrically
  location.id = locateRDT2D(Q, as.matrix(object$Vertex.RDT[, 2:3]))
  inside.id = (1:n)[location.id > 0]
  outside.id = (1:n)[location.id <= 0]
  for (i in inside.id) {
    # get cells that are dual to the triangle we found
    cell.id = unlist(subset(object$Vertex.RDT, cell == location.id[i], select = id))
    
    # we round to the 8th digit for the vertices comparison
    round.verts = round(subset(object$Vertex.RVD, cell %in% cell.id, select = c("x", "y")), 8)
    rvd.vert.freq = as.matrix(aggregate(row.names(round.verts) ~ ., data = round.verts, length))
    
    # we search for the common vertex shared by all three cells from the RVD
    # there should always be such a vertex since we located the query point
    # on the dual (RDT) of RVD
    otm.ranks[i, ] = rvd.vert.freq[match(3, rvd.vert.freq[, 3]), 1:2]
  }
  
  # numerical mapping
  if (length(outside.id)) {
    acc.verts = c(0, cumsum(as.vector(table(object$Vertex.RVD$cell))))
    otm.ranks[outside.id, ] = dualPotential2D(Q[outside.id, , drop = F], object$Data, 
                                              as.matrix(object$Vertex.RVD[, 2:3]), object$Height, 
                                              acc.verts)$optimal.vertex
  }
  
  return(otm.ranks)
}