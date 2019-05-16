GoF_2D = function(X, Y, mc, type = "max", epsilon = 1e-3, maxit = 100, na.rm = F) {
  N = nrow(X) + nrow(Y)
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) < 2 || ncol(Y) < 2) {
    stop("Data must be matrices with ncol >= 2.")
  }
  
  if (na.rm) {
    X = X[complete.cases(X), ]
    Y = Y[complete.cases(Y), ]
  }
  
  U = sobol(mc, 2)
  gof_list = GoF2D(X, Y, XY, U)
  cell_id = c(gof_list$U_Map_X, nrow(X) + gof_list$U_Map_Y)
  colnames(gof_list$Vert_XY) = c("vert.x", "vert.y", "cell")
  l2norm = rowSums(gof_list$Vert_XY[, 1:2]^2)
  
  max_norm = rep(0, N)
  for (i in 1:N) {
    vert_ids = which(gof_list$Vert_XY[, 3] == i - 1)
    max_norm[i] = max(l2norm[vert_ids])
  }
  class(otm_list) = "OTM_2D"
  
  return(otm_list)
}