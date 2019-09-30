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
#' @details Given samples \eqn{X_1, \dots, X_n} and \eqn{Y_1, \dots, Y_m}, we use the following statistic for goodness-of-fit testing:
#' \deqn{T_{X,Y}=\int_{[0, 1]^d}\|\hat{R}_{X,Y}[\hat{Q}_{X}(u)]-\hat{R}_{X,Y}[\hat{Q}_{Y}(u)]\|^2\,d\mu(u),} 
#' where \eqn{\mu\sim U[0, 1]^d}. Evaluation of this integral is done through quasi-Monte-Carlo using Sobol sequence.
#' @return the value of the goodness-of-fit test statistic.
#' @seealso \code{\link{otm.rank}} for optimal transport rank.
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

#' 1D (Permutation) Test of Independence
#' 
#' \code{otm.dep.test} computes the 1D mutual independence test statistic using ranks defined through 
#' the semi-discrete optimal transport map.
#' @param X input data vector.
#' @param Y input data vector.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param n.perm number of permutations used for computing \eqn{p}-value.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be "\code{max}", "\code{min}", "\code{center}", or "\code{uniform}".
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @details \code{otm.dep.test} tests the null hypothesis that \eqn{X} and \eqn{Y} are independent. 
#' The \eqn{p}-value are computed through permuatation. For very samll samlple (size less than \eqn{8}), 
#' all possible permutations are generated (assuming a suitable \code{n.perm} is provided). 
#' For larger sample size Monte Carlo permutation sampling is used to approximate the \eqn{p}-value.
#' 
#' Given samples \eqn{(X_1, Y_1), \dots, (X_n, Y_n)}, the following statistic is used for test of independence:
#' \deqn{T_n=\sum_{i=1}^n\|\hat{R}(X_i, Y_i)-\tilde{R}(X_i, Y_i)\|^2.} 
#' This statistics will converge to \eqn{0} under the null hypothesis.
#' @return a list, which contains the values of the independence test statistics,
#' together with the \eqn{p}-value computed based on that sample.
#' @seealso \code{\link{otm.rank}} for semi-discrete optimal transport rank.
#' @keywords htest, multivariate
#' @importFrom stats ecdf
#' @export
otm.dep.test = function(X,
                        Y,
                        scale = c(0, 1),
                        n.perm = 1,
                        rank.data = "uniform",
                        epsilon = 1e-6,
                        maxit = 100,
                        verbose = F) {
  # validate inputs
  if (NCOL(X) != 1 || NCOL(Y) != 1) {
    stop("Input data must be vectors or 1D matrices.")
  }
  
  if (NROW(X) == 1 || NROW(Y) == 1) {
    stop("Input data contains only one sample.")
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
  }
  N = nrow(tempXY)
  
  if (n.perm < 1) {
    stop("n.perm must be an integera greater than or equal to 1.")
  }
  
  if (n.perm < 499) {
    warning("n.perm is relatively small. The p-value given may not be reliable.")
  }
  
  # switch method for rank mapping
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  # for small sample size we use all possible permutations
  # otherwise we generate distinct random permutations
  if (2 * N < 8) {
    n.perm = min(n.perm, factorial(N))
    perm.mat = pc.perm(2 * N, n.perm)
  } else {
    perm.mat = rand.perm(2 * N, n.perm)
  }
  
  # permute and rescale the data
  # the first one is always the original
  XY = vector("list", n.perm)
  XY[[1]] = tempXY
  if (n.perm > 1) {
    if (is.null(scale)) {
      scale = range(X, Y)
    }
    
    for (i in 2:n.perm) {
      tempXY = matrix(c(X, Y)[perm.mat[i, ]], N)
      tempXY = scaling.min.max(tempXY, scale[1], scale[2])
      XY[[i]] = tempXY
    }
  }
  
  # compute quantiles
  # d = ifelse(is.vector(X) && is.vector(Y), 1, max(ncol(X), ncol(Y)))
  dep.stat = rep(0, n.perm)
  for (i in 1:n.perm) {
    dep.elem = dep1D(XY[[i]], rank.id == 0, epsilon, maxit, verbose)
    
    # assign ranks
    rpX = rank(XY[[i]][, 1])
    rpY = rank(XY[[i]][, 2])
    if (rank.id == 0) {
      dep.rank.hat = dep.elem
      dep.rank.tilde = cbind((2 * rpX - 1) / (2 * N), 
                             (2 * rpY - 1) / (2 * N))
    } else if (rank.id == 3) {
      dep.rank.hat = lapply(split(dep.elem[, -1], dep.elem[, 1]), matrix, ncol = 2)
      dep.rank.hat = uniform.rank(dep.rank.hat)
      dep.rank.tilde = cbind(runif(N, min = rpX - 1, max = rpX) / N, 
                             runif(N, min = rpY - 1, max = rpY) / N)
    } else {
      dep.rank.hat = lapply(split(dep.elem[, -1], dep.elem[, 1]), matrix, ncol = 2)
      dep.rank.hat = t(sapply(dep.rank.hat, choose.vert, type = rank.id))
      dep.rank.tilde = cbind((rpX - (rank.id == 2)) / N,
                             (rpY - (rank.id == 2)) / N)
    }
    
    # compute the test statistic
    dep.stat[i] = mean(rowSums((dep.rank.hat - dep.rank.tilde)^2))
  }
  
  return(list(perm.stat = dep.stat,
              p.value = (1 + sum(dep.stat[-1] >= dep.stat[1])) / (1 + n.perm)))
}