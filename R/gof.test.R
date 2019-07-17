#' 2D Goodness-of-test Statistic
#' 
#' Computes the 2D goodness-of-test statistic using rank defined through the semi-discrete optimal transport map.
#' @param X sample data matrix, of size \eqn{n} by \eqn{2}.
#' @param Y sample data matrix, of size \eqn{m} by \eqn{2}.
#' @param mc number of quasi-Monte-Carlo samples used to evaluate the test statistic.
#' @param rank choose the method for assigning ranks to the data points. Can be "max", "min", or "center".
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether NA values should be stripped before the computation proceeds.
#' @keywords multivariate, htest
#' @importFrom randtoolbox sobol
#' @export
otm.gof.test = function(X, Y, mc = 1000, rank = "center", epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  # validate inputs
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) < 2 || ncol(Y) < 2) {
    stop("Input data must be matrices with ncol > 2.")
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
    gof_list = GoF2D(X, Y, U, rank.id == 0, epsilon, maxit, verbose)
  }
  
  # assign ranks
  if (rank.id == 0) {
    gof_rank = gof_list$Cent_XY
  } else {
    gof_rank = split(data.frame(gof_list$Vert_XY)[, -1], gof_list$Vert_XY[, 1])
    gof_rank = t(sapply(gof_rank, choose_vert, type = rank.id))
  }
  
  
  # compute the test statistics
  gof_stat = mean(rowSums((gof_rank[gof_list$U_Map_X, ] - gof_rank[gof_list$U_Map_Y + nrow(X), ])^2))
  
  return(gof_stat)
}

# helper for choosing the vertex with the min/max norm.
choose_vert = function(V, type = 1) {
  if (type == 1) {
    as.numeric(V[which.max(rowSums(V^2)), , drop = F])
  } else if (type == 2) {
    as.numeric(V[which.min(rowSums(V^2)), , drop = F])
  } else {
    stop("Unknown method of rank mapping!")
  }
}