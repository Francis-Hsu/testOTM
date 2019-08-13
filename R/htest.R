#' 2D Goodness-of-fit Test
#' 
#' \code{otm.gof.test} computes the 2D goodness-of-test statistic using ranks defined through the semi-discrete optimal transport map.
#' @param X input data matrix, of size \eqn{n} by \eqn{2}.
#' @param Y input data matrix, of size \eqn{m} by \eqn{2}.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param rank choose the method for assigning ranks to the data points. 
#' Can be "\code{max}", "\code{min}", "\code{center}", or "\code{uniform}".
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return the value of the goodness-of-fit test statistic.
#' @keywords htest, multivariate
#' @importFrom randtoolbox sobol
#' @export
otm.gof.test = function(X, Y, mc = 1000, rank = "center", epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  # validate inputs
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) < 2 || ncol(Y) < 2) {
    stop("Input data must be matrices with more than 2 columns.")
  }
  
  if (ncol(X) !=  ncol(Y)) {
    stop("Mismatch in dimensions of the input data.")
  }
  
  if (na.rm) {
    X = X[complete.cases(X), ]
    Y = Y[complete.cases(Y), ]
  }
  
  rank.id = switch(rank,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  # compute quantiles
  d = ncol(X)
  U = sobol(mc, d)
  XY = rbind(X, Y)
  if (d == 2) {
    gof.list = gof2D(X, Y, U, rank.id == 0, epsilon, maxit, verbose)
  }
  
  # assign ranks
  if (rank.id == 0) {
    gof.rank = gof.list$Elem_XY
  } else if (rank.id == 3) {
    gof.rank = lapply(split(gof.list$Elem_XY[, -1], gof.list$Elem_XY[, 1]), matrix, ncol = 2)
    gof.rank = uniform.rank(gof.rank)
  } else {
    gof.rank = lapply(split(gof.list$Elem_XY[, -1], gof.list$Elem_XY[, 1]), matrix, ncol = 2)
    gof.rank = t(sapply(gof.rank, choose.vert, type = rank.id))
  }
  
  # compute the test statistic
  gof.stat = mean(rowSums((gof.rank[gof.list$U_Map_X, ] - gof.rank[gof.list$U_Map_Y + nrow(X), ])^2))
  
  return(gof.stat)
}

#' 1D Test of Independence
#' 
#' \code{otm.dep.test} computes the 1D mutual independence test statistic using ranks defined through the semi-discrete optimal transport map.
#' @param X input data vector.
#' @param Y input data vector.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param rank choose the method for assigning ranks to the data points. 
#' Can be "\code{max}", "\code{min}", "\code{center}", or "\code{uniform}".
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return the value of the  mutual independence test statistic.
#' @keywords htest, multivariate
#' @importFrom randtoolbox sobol
#' @importFrom stats ecdf
#' @export
otm.dep.test = function(X, Y, mc = 1000, rank = "center", epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  # validate inputs
  if ((!is.vector(X) && ncol(X) != 1) || (!is.vector(Y) && ncol(Y) != 1)) {
    stop("Input data must be vectors or 1D matrices.")
  }
  
  if (length(X) != length(Y)) {
    stop("Mismatch in lengths of the input data.")
  }
  
  XY = cbind(X, Y)
  if (na.rm) {
    X = X[complete.cases(XY)]
    Y = Y[complete.cases(XY)]
    XY = XY[complete.cases(XY),]
  }
  
  rank.id = switch(rank,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  # compute quantiles
  d = ifelse(is.vector(X) && is.vector(Y), 1, max(ncol(X), ncol(Y)))
  U = sobol(mc, d + 1)
  if (d == 1) {
    dep.list = dep1D(XY, U, rank.id == 0, epsilon, maxit, verbose)
  }
  
  # assign ranks
  if (rank.id == 0) {
    dep.rank = dep.list$Elem_XY
  } else if (rank.id == 3) {
    dep.rank = lapply(split(dep.list$Elem_XY[, -1], dep.list$Elem_XY[, 1]), matrix, ncol = 2)
    dep.rank = uniform.rank(dep.rank)
  } else {
    dep.rank = lapply(split(dep.list$Elem_XY[, -1], dep.list$Elem_XY[, 1]), matrix, ncol = 2)
    dep.rank = t(sapply(dep.rank, choose.vert, type = rank.id))
  }
  
  # compute the test statistic
  dep.stat = mean(rowSums((dep.rank[dep.list$U_Map_XY, ] - cbind(ecdf(X)(XY[dep.list$U_Map_XY, 1]), 
                                                                 ecdf(Y)(XY[dep.list$U_Map_XY, 2])))^2))
  
  return(dep.stat)
}

# helper for choosing the vertex with the min/max norm.
choose.vert = function(V, type = 1) {
  if (type == 1) {
    as.numeric(V[which.max(rowSums(V^2)), , drop = F])
  } else if (type == 2) {
    as.numeric(V[which.min(rowSums(V^2)), , drop = F])
  } else {
    stop("Unknown method of rank mapping!")
  }
}

# helper for uniform sampling over a convex polygon
uniform.rank = function(V) {
  n.cell = length(V)
  
  # random numbers for barycentric sampling
  r1 = sqrt(runif(n.cell))
  r2 = runif(n.cell)
  
  # sample 1 point from each cell using fan triangulation implicitly
  unif.rank = matrix(0, n.cell, 2)
  for (i in 1:n.cell) {
    n.vert = nrow(V[[i]])
    
    # area of the fan triangles
    tri.area = rep(0, n.vert - 2)
    A = V[[i]][1, ]
    for (j in 1:(n.vert - 2)) {
      B = V[[i]][j + 1, ]
      C = V[[i]][j + 2, ]
      tri.area[j] = 0.5 * abs(A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]))
    }
    
    # sample a triangle, then sample uniformly from that triangle
    tri.id = sample(1:(n.vert - 2), 1, prob = tri.area)
    B = V[[i]][tri.id + 1, ]
    C = V[[i]][tri.id + 2, ]
    unif.rank[i, 1] = (1 - r1[i]) * A[1] + r1[i] * (1 - r2[i]) * B[1] + r1[i] * r2[i] * C[1]
    unif.rank[i, 2] = (1 - r1[i]) * A[2] + r1[i] * (1 - r2[i]) * B[2] + r1[i] * r2[i] * C[2]
  }
  
  return(unif.rank)
}