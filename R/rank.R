#' Semi-discrete Optimal Transport Rank
#'
#' \code{otm.rank} computes the optimal transport ranks.
#' @param object a fitted optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param use.geo logical indicating if the geometric method should be used to compute the ranks.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the ranks of the data and the corresponding convex conjugate potential values.
#' If \code{use.geo = TRUE} then the convex conjugates will not be computed (\code{NA}s will be returned).
#' @keywords multivariate
#' @export
otm.rank = function(object, Q, use.geo = FALSE, ...) {
  UseMethod("otm.rank")
}

#' 2D Semi-discrete Optimal Transport Rank
#'
#' The 2D implementation of \code{otm.rank}.
#' @param object a fitted 2D optimal transport map object.
#' @param Q a numeric matrix where each row represents a query point.
#' @param use.geo logical indicating if the geometric method should be used to compute the ranks.
#' @return a list containing the ranks of the data and the corresponding convex conjugate potential values.
#' If \code{use.geo == TRUE} then the convex conjugates will not be computed (\code{NA}s will be returned).
#' @keywords internal
#' @export
otm.rank.otm.2d = function(object, Q, use.geo = FALSE) {
  n = nrow(Q)
  d = 2
  
  # containers for the result
  otm.dual.potential = rep(NA_integer_, n)
  otm.ranks = matrix(0, n, d)
  use.num.id = !logical(n)
  
  # compute the ranks geometrically
  if (use.geo) {
    location.id = locateRDT2D(Q, as.matrix(object$Vertex.RDT[, 2:3]))
    inside.id = (1:n)[location.id > 0]
    use.num.id[inside.id] = F
    for (i in inside.id) {
      # get cells that are dual to the triangle we found
      cell.id = unlist(subset(object$Vertex.RDT, cell == location.id[i], select = id))
      
      # we round to the 8th digit for the vertices comparison
      round.verts = round(subset(object$Vertex.RVD, cell %in% cell.id, select = c("x", "y")),
                          8)
      rvd.vert.freq = as.matrix(aggregate(row.names(round.verts) ~ ., data = round.verts, length))
      
      # we search for the common vertex shared by all three cells from the RVD
      # there should always be such a vertex since we located the query point
      # on the dual (RDT) of RVD (unless with pathological data)
      # need to catch exceptions
      otm.ranks[i, ] = rvd.vert.freq[match(3, rvd.vert.freq[, 3]), 1:2]
    }
  }
  
  # numerical mapping
  if (any(use.num.id)) {
    acc.verts = c(0, cumsum(as.vector(table(
      object$Vertex.RVD$cell
    ))))
    dual.potential = dualPotential2D(Q[use.num.id, , drop = F],
                                     object$Data,
                                     as.matrix(object$Vertex.RVD[, 2:3]),
                                     object$Height,
                                     acc.verts)
    otm.dual.potential[use.num.id] = dual.potential$dual.potential
    otm.ranks[use.num.id, ] = dual.potential$optimal.vertex
  }
  
  return(list(Rank = otm.ranks, Dual.Potential = otm.dual.potential))
}