#' Semi-discrete Optimal Transport Rank
#'
#' \code{tos.rank} computes the optimal transport ranks.
#' @param object a fitted optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be \code{max}, \code{min}, \code{center}, or \code{uniform}.
#' @param rank.algo choose the algorithm for computing the ranks. Can be \code{lp} or \code{geom}.
#' @param \dots additional arguments.
#' @details The rank of the a data points used to fit the optimal transport map is not uniquely defined. 
#' In principle it can be any point within the corresponding Laguerre cell. 
#' Set \code{rank.data} to \code{max} to use the vertices with the maximum l2 norm (\code{min} for minimum l2 norm),
#' which is analogous to the usual ECDF rank in 1D. Option \code{center} will assign the centroid of the cell as the rank.
#' Finally, option \code{uniform} will sample a point uniformly from the Laguerre cell as the rank.
#' 
#' The rank for a non-data query point can be computed by solving a series of linear programs, 
#' or by finding the intersection of Laguerre cells that correspond to the simplex containing the query.
#' Set \code{rank.algo} to \code{lp} to use the former and \code{geom} for the latter.
#' Note the geometric method is slow and does not work for points lie outside the RDT, it is also not implemented for 3D, 
#' \code{lp} method will be used in these cases.
#' @return a list containing the ranks of the data and the corresponding convex conjugate potential values.
#' If \code{rank.algo = "geom"} or data points are presented in the queries,
#' then the corresponding conjugate potentials will not be computed (\code{NA}s will be returned).
#' @seealso \code{\link{tos.gof.test}} and \code{\link{tos.dep.test}} for statistical tests using optimal transport rank.
#' @keywords multivariate
#' @export
tos.rank = function(object, query, scale = TRUE, rank.data = "uniform", rank.algo = "lp", ...) {
  UseMethod("tos.rank")
}

#' 2D Semi-discrete Optimal Transport Rank
#'
#' The 2D implementation of \code{tos.rank}.
#' @param object a fitted 2D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be \code{max}, \code{min}, \code{center}, or \code{uniform}.
#' @param rank.algo choose the algorithm for computing the ranks. Can be \code{lp} or \code{geom}.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the ranks of the data and the corresponding convex conjugate potential values.
#' If \code{rank.algo = "geom"} then the conjugate potentials will not be computed (\code{NA}s will be returned).
#' @keywords internal
#' @importFrom stats aggregate na.omit
#' @export
tos.rank.tos.2d = function(object, query, scale = TRUE, rank.data = "uniform", rank.algo = "lp", ...) {
  n = nrow(query)
  d = 2
  
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  algo.id = switch(rank.algo,
                   lp = 0,
                   geom = 1,
                   stop("Unknown rank algorithm!"))
  
  # containers for the result
  tos.dual.potential = rep(NA_integer_, n)
  tos.ranks = matrix(0, n, d)
  lp.id = !logical(n)
  
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  # find data points (if any) in the queries
  data.match = match(split(round(query, 8), 1:n),
                     split(round(object$Data, 8), 1:nrow(object$Data)))
  if (!all(is.na(data.match))) {
    query.data.id = which(!is.na(data.match))
    data.id = data.match[query.data.id]
    
    # assign ranks to data points
    if (rank.id == 0) {
      tos.ranks[query.data.id, ] = object$Centroid[data.id, ]
    } else {
      data.rank = subset(object$Vertex.RVD, object$Vertex.RVD[, 1] %in% data.id)
      data.rank = lapply(split(data.rank[, -1], data.rank[, 1]), matrix, ncol = 2)
      
      if (rank.id == 3) {
        tos.ranks[query.data.id, ] = uniform.rank(data.rank)
      } else {
        tos.ranks[query.data.id, ] = t(sapply(data.rank, choose.vert, type = rank.id))
      }
    }
    
    lp.id[query.data.id] = F
  }
  
  # compute the ranks geometrically
  if (algo.id == 1) {
    location.id = locateRDT2D(query, as.matrix(object$Vertex.RDT[, 3:4]))
    inside.id = (1:n)[location.id > 0]
    lp.id[inside.id] = F
    for (i in inside.id) {
      # get cells that are dual to the triangle we found
      cell.id = unlist(subset(object$Vertex.RDT, object$Vertex.RDT[, 1] == location.id[i], select = 4))
      
      # we round to the 8th digit for the vertices comparison
      round.verts = round(subset(object$Vertex.RVD, object$Vertex.RVD[, 1] %in% cell.id, select = 2:3),
                          8)
      rvd.vert.freq = as.matrix(aggregate(1:nrow(round.verts) ~ ., data = round.verts, length))
      
      # we search for the common vertex shared by all three cells from the RVD
      # there should always be such a vertex since we located the query point
      # on the dual (RDT) of RVD (unless with pathological data)
      # need to catch exceptions
      tos.ranks[i, ] = rvd.vert.freq[match(3, rvd.vert.freq[, 3]), 1:2]
    }
  }
  
  # compute the ranks through LP
  if (any(lp.id)) {
    acc.verts = c(0, cumsum(as.vector(table(object$Vertex.RVD[, 1]))))
    dual.potential = dualPotential2D(query[lp.id, , drop = F],
                                     object$Data,
                                     object$Vertex.RVD[, 2:3],
                                     object$Height,
                                     acc.verts)
    tos.dual.potential[lp.id] = dual.potential$dual.potential
    tos.ranks[lp.id, ] = dual.potential$optimal.vertex
  }
  
  return(list(Rank = tos.ranks, Dual.Potential = tos.dual.potential))
}

#' 3D Semi-discrete Optimal Transport Rank
#'
#' The 3D implementation of \code{tos.rank}.
#' @param object a fitted 3D optimal transport map object.
#' @param query a numeric matrix where each row represents a query point.
#' @param scale logical indicating if the queries should be scaled.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be \code{max}, \code{min}, or \code{center}.
#' @param \dots additional arguments, currently without effect.
#' @return a list containing the ranks of the data and the corresponding convex conjugate potential values.
#' @keywords internal
#' @importFrom stats aggregate na.omit
#' @export
tos.rank.tos.3d = function(object,
                           query,
                           scale = TRUE,
                           rank.data = "center",
                           ...) {
  n = nrow(query)
  d = 3
  
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   stop("Unknown method of rank mapping!"))
  
  # containers for the result
  tos.dual.potential = rep(NA_integer_, n)
  tos.ranks = matrix(0, n, d)
  lp.id = !logical(n)
  
  # scale the queries
  if (scale) {
    query = scale(query, object$Location, object$Scale)
    query = query * diff(range(object$Data)) + min(object$Data)
  }
  
  # find data points (if any) in the queries
  data.match = match(split(round(query, 8), 1:n),
                     split(round(object$Data, 8), 1:nrow(object$Data)))
  if (!all(is.na(data.match))) {
    query.data.id = which(!is.na(data.match))
    data.id = data.match[query.data.id]
    
    # assign ranks to data points
    if (rank.id == 0) {
      tos.ranks[query.data.id, ] = object$Centroid[data.id, ]
    } else {
      data.rank = subset(object$Vertex.RVD, object$Vertex.RVD[, 1] %in% data.id)[, c(1, 4:6)]
      data.rank = lapply(split(data.rank[, -1], data.rank[, 1]), matrix, ncol = 3)
      tos.ranks[query.data.id, ] = t(sapply(data.rank, choose.vert, type = rank.id))
    }
    
    lp.id[query.data.id] = F
  }
  
  # get unique vertices from each cell
  unique.vert = lapply(unique(object$Vertex.RVD[, 1]), function(i) {
    sub.cell = subset(object$Vertex.RVD, object$Vertex.RVD[, 1] == i)
    sub.cell[!duplicated(sub.cell[, 3]),]
  })
  unique.vert = do.call(rbind, unique.vert)
  
  # compute the ranks through LP
  if (any(lp.id)) {
    acc.verts = c(0, cumsum(as.vector(table(unique.vert[, 1]))))
    dual.potential = dualPotential3D(query[lp.id, , drop = F],
                                     object$Data,
                                     unique.vert[, 4:6],
                                     object$Height,
                                     acc.verts)
    tos.dual.potential[lp.id] = dual.potential$dual.potential
    tos.ranks[lp.id,] = dual.potential$optimal.vertex
  }
  
  return(list(Rank = tos.ranks, Dual.Potential = tos.dual.potential))
}