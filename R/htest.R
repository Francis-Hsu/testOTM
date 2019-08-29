#' 2D Goodness-of-fit Test
#' 
#' \code{otm.gof.test} computes the 2D goodness-of-test statistic using ranks defined through the semi-discrete optimal transport map.
#' @param X input data matrix, of size \eqn{n} by \eqn{2}.
#' @param Y input data matrix, of size \eqn{m} by \eqn{2}.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be "\code{max}", "\code{min}", "\code{center}", or "\code{uniform}".
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return the value of the goodness-of-fit test statistic.
#' @keywords htest, multivariate
#' @importFrom randtoolbox sobol
#' @export
otm.gof.test = function(X,
                        Y,
                        mc = 10000,
                        scale = c(0, 1),
                        rank.data = "uniform",
                        epsilon = 1e-6,
                        maxit = 100,
                        verbose = F,
                        na.rm = F) {
  # validate inputs
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) < 2 || ncol(Y) < 2) {
    stop("Input data must be matrices with more than 2 columns.")
  }
  
  if (NCOL(X) !=  NCOL(Y)) {
    stop("Mismatch in dimensions of the input data.")
  }
  
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  if (!is.null(scale)) {
    if (scale[1] < 0 || scale[2] > 1) {
      stop("Scaling range must be within [0, 1].")
    }
    
    # scale X and Y together
    tempXY = rbind(X, Y)
    tempXY = scaling.min.max(tempXY, scale[1], scale[2])
    X = tempXY[1:NROW(X), ]
    Y = tempXY[1:NROW(Y) + NROW(X), ]
  }
  
  if (na.rm) {
    X = X[complete.cases(X), ]
    Y = Y[complete.cases(Y), ]
  }
  XY = rbind(X, Y)
  
  # compute quantiles
  d = ncol(X)
  U = sobol(mc, d)
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
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be "\code{max}", "\code{min}", "\code{center}", or "\code{uniform}".
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return the value of the  mutual independence test statistic.
#' @keywords htest, multivariate
#' @importFrom stats ecdf
#' @export
otm.dep.test = function(X,
                        Y,
                        scale = c(0, 1),
                        rank.data = "uniform",
                        epsilon = 1e-6,
                        maxit = 100,
                        verbose = F,
                        na.rm = F) {
  # validate inputs
  if (NCOL(X) != 1 || NCOL(Y) != 1) {
    stop("Input data must be vectors or 1D matrices.")
  }
  
  if (length(X) != length(Y)) {
    stop("Mismatch in lengths of the input data.")
  }
  
  if (!is.null(scale)) {
    if (scale[1] < 0 || scale[2] > 1) {
      stop("Scaling range must be within [0, 1].")
    }
    
    # scale X and Y together
    tempXY = cbind(X, Y)
    tempXY = scaling.min.max(tempXY, scale[1], scale[2])
    X = tempXY[, NCOL(X)]
    Y = tempXY[, NCOL(Y) + NCOL(X)]
  }
  
  XY = cbind(X, Y)
  if (na.rm) {
    X = X[complete.cases(XY)]
    Y = Y[complete.cases(XY)]
    XY = XY[complete.cases(XY),]
  }
  N = nrow(XY)
  
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  # compute quantiles
  d = ifelse(is.vector(X) && is.vector(Y), 1, max(ncol(X), ncol(Y)))
  if (d == 1) {
    dep.elem = dep1D(XY, rank.id == 0, epsilon, maxit, verbose)
  }
  
  # assign ranks
  if (rank.id == 0) {
    dep.rank.hat = dep.elem
    dep.rank.tilde = cbind((2 * rank(X) - 1) / (2 * N), 
                           (2 * rank(Y) - 1) / (2 * N))
  } else if (rank.id == 3) {
    dep.rank.hat = lapply(split(dep.elem[, -1], dep.elem[, 1]), matrix, ncol = 2)
    dep.rank.hat = uniform.rank(dep.rank.hat)
    dep.rank.tilde = cbind(runif(N, min = rank(X) - 1, max = rank(X)) / N, 
                           runif(N, min = rank(Y) - 1, max = rank(Y)) / N)
  } else {
    dep.rank.hat = lapply(split(dep.elem[, -1], dep.elem[, 1]), matrix, ncol = 2)
    dep.rank.hat = t(sapply(dep.rank.hat, choose.vert, type = rank.id))
    dep.rank.tilde = cbind((rank(X) - (rank.id == 2)) / N,
                           (rank(Y) - (rank.id == 2)) / N)
  }
  
  # compute the test statistic
  dep.stat = mean(rowSums((dep.rank.hat - dep.rank.tilde)^2))
  
  return(dep.stat)
}