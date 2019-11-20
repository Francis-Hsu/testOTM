#' 2D/3D Goodness-of-fit Test
#' 
#' \code{tos.gof.test} computes the 2D/3D goodness-of-test statistic using ranks defined through the semi-discrete optimal transport map.
#' @param X input data matrix, of size \eqn{n} by \eqn{2} or \eqn{3}.
#' @param Y input data matrix, of size \eqn{m} by \eqn{2} or \eqn{3}.
#' @param mc number of quasi-Monte Carlo samples used to evaluate the test statistic.
#' @param n.perm number of permutations used for computing \eqn{p}-value.
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be \code{max}, \code{min}, \code{center}, or \code{uniform}. 
#' \code{uniform} is not currently implemented for 3D, and \code{center} will be used if it is chosen.
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @details \code{tos.gof.test} tests whether the samples come from the same population. 
#' The \eqn{p}-value is computed through permutations. For a very small sample (size less than \eqn{8}), 
#' all possible permutations are generated (assuming a suitable \code{n.perm} is provided). 
#' For larger sample size, Monte Carlo permutation sampling is used to approximate the permutations \eqn{p}-value.
#' 
#' Given samples \eqn{X_1, \dots, X_n} and \eqn{Y_1, \dots, Y_m}, we use the following statistic for goodness-of-fit testing:
#' \deqn{T_{X,Y}=\int_{[0, 1]^d}\|\hat{R}_{X,Y}[\hat{Q}_{X}(u)]-\hat{R}_{X,Y}[\hat{Q}_{Y}(u)]\|_2^2\,d\mu(u),} 
#' where \eqn{\mu\sim U[0, 1]^d}. Evaluation of this integral is done through quasi-Monte Carlo with Sobol sequence.
#' @return a list which contains the permutation test statistics and the \eqn{p}-value.
#' @seealso \code{\link{tos.rank}} for optimal transport rank.
#' @keywords htest, multivariate
#' @importFrom randtoolbox sobol
#' @export
tos.gof.test = function(X, Y, mc = 10000, n.perm = 0, scale = c(0, 1), rank.data = "uniform", epsilon = 1e-6, maxit = 100, verbose = F, na.rm = F) {
  # get number of elements
  D = NCOL(X)
  NX = NROW(X)
  NY = NROW(Y)
  N = NX + NY
  
  # validate inputs
  if (D != NCOL(Y)) {
    stop("Mismatch in dimensions of the input data.")
  }
  
  if (!is.matrix(X) || !is.matrix(Y) || D < 2 || D > 3) {
    stop("Input data must be matrices with 2 or 3 columns.")
  }
  
  if (D == 2) {
    if (NX < 3 || NY < 3) {
      stop("Input data should contain at least 3 samples for 2D testing.")
    }
  } else {
    if (NX < 4 || NY < 4) {
      stop("Input data should contain at least 4 samples for 3D testing.")
    }
  }

  if (n.perm < 0 || n.perm %% 1 != 0) {
    stop("n.perm must be an integera greater than or equal to 0.")
  }
  
  if (n.perm != 0 && n.perm < 500) {
    message("Consider setting n.perm to a larger value, otherwise the permutation p-value may not be reliable.")
  }
  
  if (!is.null(scale)) {
    if (scale[1] < 0 || scale[2] > 1) {
      stop("Scaling range must be within [0, 1].")
    }
  }
  
  if (na.rm) {
    X = X[complete.cases(X), ]
    Y = Y[complete.cases(Y), ]
    NX = NROW(X)
    NY = NROW(Y)
  }
  
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  if (D == 3 && rank.id == 3) {
    rank.id = 0
  }
  
  # scale X and Y together
  XY = scaling.min.max(rbind(X, Y), scale[1], scale[2])
  
  # generate quasi-MC sequence
  U = sobol(mc, D)
  
  # for small sample size we use all possible permutations, 
  # otherwise we generate distinct random permutations for moderate sample size
  # for large sample size it's nigh impossible to see collision
  if (n.perm != 0) {
    # check if n.perm too large
    if (n.perm > factorial(N) - 1) {
      n.perm = factorial(N) - 1
      warning("n.perm exceeds maximal number of permutations.")
    }
    
    if (N < 8) {
      perm.mat = pc.perm(N, n.perm + 1)[-1, ] # first row is the original order
    } else if (N <= 100) {
      perm.mat = rand.perm(N, n.perm)
    }
  }
  
  # storage for test statistics
  gof.stats = rep(0, n.perm + 1)
  
  # first test statistic (for the original samples)
  if (D == 2) {
    # compute the first joint RVD
    gof.elem = jointRankHelper2D(XY, rank.id == 0, epsilon, maxit, verbose)
    gof.rank = gof.assign.rank(gof.elem, rank.id, D)
    gof.list = gof2DHelper(XY[1:NX, ], XY[(NX + 1):(NX + NY), ], U, epsilon, maxit, verbose)
  } else {
    gof.elem = jointRankHelper3D(XY, rank.id == 0, epsilon, maxit, verbose)
    if (rank.id == 0) {
      gof.elem = gof.elem[[1]]
    } else {
      gof.elem = do.call(rbind, lapply(seq_along(gof.elem), function(i) {
        cbind(i, gof.elem[[i]][, 3:5])
      }))
    }
    gof.rank = gof.assign.rank(gof.elem, rank.id, D)
    gof.list = gof3DHelper(XY[1:NX, ], XY[(NX + 1):(NX + NY), ], U, epsilon, maxit, verbose)
  }
  gof.stats[1] = mean(rowSums((gof.rank[gof.list$U_Map_X, ] - gof.rank[gof.list$U_Map_Y + NX, ])^2))
  
  # compute the permutation test statistics
  if (n.perm != 0) {
    for (i in 1:n.perm) {
      # get permutation id
      if (N > 100) {
        perm.id = sample(1:N)
      } else {
        perm.id = perm.mat[i, ]
      }
      
      # permute the input
      tempXY = XY[perm.id, ]
      tempX = tempXY[1:NX, ]
      tempY = tempXY[-(1:NX), ]
      
      # compute permutation RVD
      if (rank.id == 0) {
        gof.elem.perm = gof.elem[perm.id, ]
      } else {
        # get the inverse permutation
        perm.inv = perm.id
        perm.inv[perm.inv] = seq_along(perm.inv)
        
        # the index column of the RVD needs to be mapped to the post-permutation order
        # which can be achieved by expanding the inverse permutation indices
        gof.elem.perm = gof.elem
        gof.elem.perm[, 1] = perm.inv[gof.elem[, 1]]
      }
      
      if (D == 2) {
        # assign ranks
        gof.rank = gof.assign.rank(gof.elem.perm, rank.id, D)
        
        # compute the quantiles of quasi-MC
        gof.list = gof2DHelper(tempX, tempY, U, epsilon, maxit, verbose)
      } else {
        gof.rank = gof.assign.rank(gof.elem.perm, rank.id, D)
        gof.list = gof3DHelper(tempX, tempY, U, epsilon, maxit, verbose)
      }
      
      # compute the test statistic
      gof.stats[i + 1] = mean(rowSums((gof.rank[gof.list$U_Map_X, ] - gof.rank[gof.list$U_Map_Y + NX, ])^2))
    }
  }
  
  return(list(perm.stats = gof.stats,
              p.value = ifelse(n.perm == 0, NA, (1 + sum(gof.stats[-1] >= gof.stats[1])) / (1 + n.perm))))
}

#' 1D (Permutation) Test of Independence
#' 
#' \code{tos.dep.test} computes the 1D mutual independence test statistic using ranks defined through 
#' the semi-discrete optimal transport map.
#' @param X input data vector.
#' @param Y input data vector.
#' @param n.perm number of permutations used for computing \eqn{p}-value.
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param rank.data choose the method for assigning ranks to the data points. 
#' Can be \code{max}, \code{min}, \code{center}, or \code{uniform}.
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @details \code{tos.dep.test} tests the null hypothesis that \eqn{X} and \eqn{Y} are independent. 
#' The \eqn{p}-value is computed through permutations. For a very small sample (size less than \eqn{8}), 
#' all possible permutations are generated (assuming a suitable \code{n.perm} is provided). 
#' For larger sample size, Monte Carlo permutation sampling is used to approximate the permutations \eqn{p}-value.
#' 
#' Given samples \eqn{(X_1, Y_1), \dots, (X_n, Y_n)}, the following statistic is used for test of independence:
#' \deqn{T_n=\sum_{i=1}^n\|\hat{R}(X_i, Y_i)-\tilde{R}(X_i, Y_i)\|^2,} 
#' where \eqn{\tilde{R}} is the joining of the usual 1D ranks.
#' @return a list which contains the permutation test statistics and the \eqn{p}-value.
#' @seealso \code{\link{tos.rank}} for semi-discrete optimal transport rank.
#' @keywords htest, multivariate
#' @export
tos.dep.test = function(X, Y, n.perm = 0, scale = c(0, 1), rank.data = "uniform", epsilon = 1e-6, maxit = 100, verbose = F) {
  # get number of elements
  NX = NROW(X)
  NY = NROW(Y)
  
  # validate inputs
  if (NCOL(X) != 1 || NCOL(Y) != 1) {
    stop("Input data must be vectors or 1D matrices.")
  }
  
  if (NX < 3 || NY < 3) {
    stop("Input data should contain more than two samples.")
  }
  
  if (NX != NY) {
    stop("Mismatch in lengths of the input data.")
  }
  N = NX
  
  if (!is.null(scale)) {
    if (scale[1] < 0 || scale[2] > 1) {
      stop("Scaling range must be within [0, 1].")
    }
  }
  
  if (n.perm < 0 || n.perm %% 1 != 0) {
    stop("n.perm must be an integera greater than or equal to 0.")
  }
  
  if (n.perm != 0 && n.perm < 500) {
    message("Consider setting n.perm to a larger value, otherwise the permutation p-value may not be reliable.")
  }
  
  # switch method for rank mapping
  rank.id = switch(rank.data,
                   center = 0,
                   max = 1,
                   min = 2,
                   uniform = 3,
                   stop("Unknown method of rank mapping!"))
  
  # for small sample size we use all possible permutations, 
  # otherwise we generate distinct random permutations for moderate sample size
  # for large sample size it's nigh impossible to see collision
  if (n.perm != 0) {
    # check if n.perm too large
    if (n.perm > factorial(N) - 1) {
      n.perm = factorial(N) - 1
      warning("n.perm exceeds maximal number of permutations.")
    }
    
    if (2 * N < 8) {
      perm.mat = pc.perm(2 * N, n.perm + 1)[-1, ] # first row is the original order
    } else if (N <= 50) {
      perm.mat = rand.perm(2 * N, n.perm)
    }
  }
  
  # scale X and Y
  tempXY = matrix(scaling.min.max(c(X, Y), scale[1], scale[2]), N)
  sX = tempXY[, 1]
  sY = tempXY[, 2]
  
  # storage for test statistics
  dep.stats = rep(0, n.perm + 1)
  
  # first test statistic (for the original samples)
  dep.elem = jointRankHelper2D(tempXY, rank.id == 0, epsilon, maxit, verbose)
  dep.rank = dep.assign.rank(tempXY, dep.elem, rank.id)
  dep.stats[1] = mean(rowSums((dep.rank$h - dep.rank$t)^2))
  
  # permute the data and compute statistics
  if (n.perm != 0) {
    for (i in 1:n.perm) {
      # get permutation id
      if (N > 50) {
        perm.id = sample(1:(2 * N))
      } else {
        perm.id = perm.mat[i, ]
      }
      
      # permute the input
      tempXY = matrix(c(sX, sY)[perm.id], N)
      
      # compute permutation RVD
      dep.elem = jointRankHelper2D(tempXY, rank.id == 0, epsilon, maxit, verbose)
      
      # assign ranks
      dep.rank = dep.assign.rank(tempXY, dep.elem, rank.id)
      
      # compute the test statistic
      dep.stats[i + 1] = mean(rowSums((dep.rank$h - dep.rank$t)^2))
    }
  }
  
  return(list(perm.stats = dep.stats,
              p.value = ifelse(n.perm == 0, NA, (1 + sum(dep.stats[-1] >= dep.stats[1])) / (1 + n.perm))))
}