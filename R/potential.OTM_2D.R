#' 2D semi-continuous optimal transport dual potential
#' 
#' Compute the Legendre-Fenchel dual of the Alexandrov potential
#' 2D semi-continuous optimal transport dual potential.
#' @param object a fitted 2D optimal transport map object.
#' @param x a \eqn{n} by \eqn{2} numeric matrix.
#' @return a matrix containing dual potentials of the data.
#' @keywords Voronoi, Laguerre
#' @importFrom lpSolveAPI make.lp set.objfn add.constraint delete.constraint get.objective get.variables solve
#' @export
rank.OTM_2D = function(object, X) {
  m = nrow(X)
  n = nrow(object$Data)
  d = 2
  h = object$Height
  
  lp.mdl = make.lp(0, n)
  set.objfn(lp.mdl, -h)
  add.constraint(lp.mdl, rep(1, n), type = "=", rhs = 1)
  
  lp.status = rep(0, m)
  lp.objective = rep(0, m)
  for (i in 1:m) {
    for (j in 1:d) {
      add.constraint(lp.mdl, object$Data[, j], type = "=", rhs = X[i, j])
    }
    
    lp.status[i] = solve(lp.mdl)
    lp.objective[i] = get.objective(lp.mdl)
    
    for (j in 1:d) {
      delete.constraint(lp.mdl, 2)
    }
  }
  delete.lp(lp.mdl)
  
  # compute the infeasible set with other method
  
  return(lp.objective)
}