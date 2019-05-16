#' 2D semi-continuous optimal transport map plot
#' 
#' Plot OTM
#' @param Obj Fitted optimal transport map.
#' @keywords Voronoi, Laguerre
#' @export
plot.OTM_2D = function(Obj, col.data = "cornflowerblue", col.center = "firebrick", col.edge = "black", connect = F, ...) {
  plot(Obj$Data[, 1], Obj$Data[, 2], xlim = c(0, 1), ylim = c(0, 1), col = col.data, pch = 20, xlab = "x", ylab = "y", ...)
  points(Obj$Centroids[, 1], Obj$Centroids[, 2], xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = col.center, ...)
  for (i in 1:nrow(Obj$Edge)) {
    segments(Obj$Edge[i, 1], Obj$Edge[i, 2], Obj$Edge[i, 3], Obj$Edge[i, 4], col = col.edge, ...)
  }
  if (connect) {
    for (i in 1:nrow(Obj$Data)) {
      segments(Obj$Data[i, 1], Obj$Data[i, 2], Obj$Centroids[i, 1], Obj$Centroids[i, 2], lty = 2, lwd = 0.5)
    }
  }
}