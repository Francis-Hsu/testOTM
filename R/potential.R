#' Semi-discrete Optimal Transport Dual Potential
#' 
#' Compute the Legendre-Fenchel dual of the Alexandrov potential.
#' @param object a fitted optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @return a matrix containing dual potentials of the data.
#' @export
otm.potential = function(object, X, ...) {
  UseMethod("otm.potential")
}

#' 2D Semi-discrete Optimal Transport Dual Potential
#' 
#' Compute the Legendre-Fenchel dual of the Alexandrov potential.
#' @param object a fitted 2D optimal transport map object.
#' @param X a \eqn{n} by \eqn{2} numeric matrix.
#' @return a matrix containing dual potentials of the data.
#' @importFrom lpSolveAPI make.lp set.objfn add.constraint delete.constraint get.objective get.variables solve.lpExtPtr
#' @export
otm.potential.OTM_2D = function(object, X) {
  m = nrow(X)
  n = nrow(object$Data)
  d = 2
  h = object$Height
  
  lp.mdl = make.lp(0, n)
  set.objfn(lp.mdl, -h)
  add.constraint(lp.mdl, rep(1, n), type = "=", rhs = 1)
  
  lp.status = rep(0, m)
  potential = rep(0, m)
  for (i in 1:m) {
    for (j in 1:d) {
      add.constraint(lp.mdl, object$Data[, j], type = "=", rhs = X[i, j])
    }
    
    lp.status[i] = solve(lp.mdl)
    potential[i] = get.objective(lp.mdl)
    
    for (j in 1:d) {
      delete.constraint(lp.mdl, 2)
    }
  }
  
  # compute the infeasible set with another method
  if (any(lp.status != 0)) {
    invalid.id = which(lp.status != 0)
    acc.verts = c(0, cumsum(as.vector(table(object$Vertex.RVD$cell))))
    potential[invalid.id] = dualPotential2D(X[invalid.id, , drop = F], object$Data, 
                                            as.matrix(object$Vertex.RVD[, 2:3]), object$Height, 
                                            acc.verts)$dual.potential
  }
  
  return(potential)
}