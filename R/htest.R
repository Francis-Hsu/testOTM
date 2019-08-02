#' 2D Goodness-of-fit Test
#' 
#' \code{otm.gof.test} computes the 2D goodness-of-test statistic using ranks defined through the semi-discrete optimal transport map.
#' @param X input data matrix, of size \eqn{n} by \eqn{2}.
#' @param Y input data matrix, of size \eqn{m} by \eqn{2}.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param rank choose the method for assigning ranks to the data points. Can be "\code{max}", "\code{min}", or "\code{center}".
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
  } else {
    gof.rank = split(data.frame(gof.list$Elem_XY)[, -1], gof.list$Elem_XY[, 1])
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
#' @param rank choose the method for assigning ranks to the data points. Can be "\code{max}", "\code{min}", or "\code{center}".
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return the value of the  mutual independence test statistic.
#' @keywords htest, multivariate
#' @importFrom randtoolbox sobol
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
  } else {
    dep.rank = split(data.frame(dep.list$Elem_XY)[, -1], dep.list$Elem_XY[, 1])
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