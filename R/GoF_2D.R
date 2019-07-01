#' 2D goodness-of-test statistic
#' 
#' Computes the 2D goodness-of-test statistic using rank defined through the semi-discrete optimal transport map.
#' @param X sample data matrix, of size \eqn{n} by \eqn{2}.
#' @param Y sample data matrix, of size \eqn{m} by \eqn{2}.
#' @param mc number of Monte Carlo iteration samples used to evaluate the test statistic.
#' @param type method used to choose the rank over the Voronoi cell boundary.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether NA values should be stripped before the computation proceeds.
#' @keywords multivariate, htest
#' @importFrom randtoolbox sobol
#' @export
GoF_2D = function(X, Y, mc = 1000, type = "max", epsilon = 1e-3, maxit = 100, verbose = F, na.rm = F) {
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) != 2 || ncol(Y) != 2) {
    stop("Input data must be matrices with ncol = 2.")
  }
  
  type_id = switch(type,
                   max = 1,
                   min = 2,
                   stop("Unknown type of rank method!"))
  
  if (na.rm) {
    X = X[complete.cases(X), ]
    Y = Y[complete.cases(Y), ]
  }
  
  U = sobol(mc, 2)
  XY = rbind(X, Y)
  gof_list = GoF2D(X, Y, XY, U, epsilon, maxit, verbose)
  cell_id_x = gof_list$U_Map_X + 1
  cell_id_y = gof_list$U_Map_Y + 1 + nrow(X)
  colnames(gof_list$Vert_XY) = c("vert.x", "vert.y", "cell")
  gof_verts = gof_list$Vert_XY
  gof_rank = t(sapply(split(data.frame(gof_verts)[, -1], gof_verts[, 1]), choose_vert, type = type_id))
  gof_stat = mean(rowSums((gof_rank[cell_id_x, ] - gof_rank[cell_id_y, ])^2))
  
  return(gof_stat)
}

choose_vert = function(V, type = 1) {
  if (type == 1) {
    as.numeric(V[which.max(rowSums(V^2)), , drop = F])
  } else if (type == 2) {
    as.numeric(V[which.min(rowSums(V^2)), , drop = F])
  } else {
    stop("Unknown type of rank method!")
  }
}
