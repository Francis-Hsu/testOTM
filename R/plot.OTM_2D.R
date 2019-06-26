#' 2D semi-continuous optimal transport map plot
#' 
#' Plot the restricted Voronoi diagram (RVD) and the restricted Delaunay triangulation (RDT) of the 
#' 2D semi-continuous optimal transport map.
#' @param object a fitted 2D optimal transport map object.
#' @param which specify which graph to plot. Can be "RVD", "RDT", or "Both".
#' @param draw.center logical indicating if the centroids should be plotted.
#' @param draw.map logical indicating if dashed lines should be added to show mapping between the data and the Voronoi cells.
#' @param col.data color of the data points
#' @param col.center color of the Voronoi centroids
#' @param col.edge color of the edges in plotting RVD and RDT
#' @keywords Voronoi, Laguerre
#' @export
plot.OTM_2D = function(object, 
                       which = "Both", 
                       col.data = "cornflowerblue", 
                       col.center = "firebrick", 
                       col.edge = "black", 
                       draw.center = T, 
                       draw.map = F, ...) {
  type = match.arg(which, c("RVD", "RDT", "Both"))
  
  # plot the restricted Voronoi diagram
  if (which == "RVD" || which == "Both") {
    # plot data
    plot(object$Data[, 1], object$Data[, 2], 
         xlim = c(0, 1), ylim = c(0, 1), col = col.data, pch = 20, 
         xlab = expression('u'[1]), ylab = expression('u'[2]), ...)
    
    # plot Laguerre cells
    n.edge = nrow(object$Vertex.RVD)
    rvd.edges = matrix(0, n.edge, 4)
    curr.row = 1
    
    # extract edges
    for (i in 1:object$N.Cells) {
      curr.cell = subset(object$Vertex.RVD, cell == i, select = c(x, y))
      for (j in 1:nrow(curr.cell)) {
        p.id = j
        q.id = 1 + j %% nrow(curr.cell)
        
        # sort by l1 norm, everything positive!
        if (curr.cell[p.id, 1] + curr.cell[p.id, 2] <= curr.cell[q.id, 1] + curr.cell[q.id, 2]) {
          rvd.edges[curr.row, ] = c(curr.cell[p.id, 1], curr.cell[p.id, 2], 
                                    curr.cell[q.id, 1], curr.cell[q.id, 2])
        } else {
          rvd.edges[curr.row, ] = c(curr.cell[q.id, 1], curr.cell[q.id, 2],
                                    curr.cell[p.id, 1], curr.cell[p.id, 2])
        }
        curr.row = curr.row + 1
      }
    }
    
    # remove duplicated edges
    rvd.edges = unique(round(rvd.edges, 8))
    
    for (i in 1:nrow(rvd.edges)) {
      segments(rvd.edges[i, 1], rvd.edges[i, 2], 
               rvd.edges[i, 3], rvd.edges[i, 4], 
               col = col.edge, ...)
    }
    
    # plot centroids
    if (draw.center) {
      points(object$Centroid[, 1], object$Centroid[, 2], 
             xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = col.center, ...)
    }
    
    # plot mappings
    if (draw.map) {
      for (i in 1:nrow(object$Data)) {
        segments(object$Data[i, 1], object$Data[i, 2], 
                 object$Centroid[i, 1], object$Centroid[i, 2], lty = 2, lwd = 0.5)
      }
    }
  }
  
  # plot the restricted Delaunay triangulation
  if (which == "RDT" || which == "Both") {
    # plot data
    plot(object$Data[, 1], object$Data[, 2], 
         xlim = c(0, 1), ylim = c(0, 1), col = col.data, pch = 20, 
         xlab = expression('u'[1]), ylab = expression('u'[2]), ...)
    
    # plot triangles
    n.edge = nrow(object$Vertex.RDT)
    rdt.edges = matrix(0, n.edge, 4)
    curr.row = 1
    
    # extract edges
    for (i in 1:object$N.Triangles) {
      curr.cell = subset(object$Vertex.RDT, cell == i, select = c(x, y))
      for (j in 1:nrow(curr.cell)) {
        p.id = j
        q.id = 1 + j %% nrow(curr.cell)
        
        # sort by l1 norm, everything positive!
        if (curr.cell[p.id, 1] + curr.cell[p.id, 2] <= curr.cell[q.id, 1] + curr.cell[q.id, 2]) {
          rdt.edges[curr.row, ] = c(curr.cell[p.id, 1], curr.cell[p.id, 2], 
                                    curr.cell[q.id, 1], curr.cell[q.id, 2])
        } else {
          rdt.edges[curr.row, ] = c(curr.cell[q.id, 1], curr.cell[q.id, 2],
                                    curr.cell[p.id, 1], curr.cell[p.id, 2])
        }
        curr.row = curr.row + 1
      }
    }
    
    # remove duplicated edges
    rdt.edges = unique(round(rdt.edges, 8))
    
    for (i in 1:nrow(rdt.edges)) {
      segments(rdt.edges[i, 1], rdt.edges[i, 2], 
               rdt.edges[i, 3], rdt.edges[i, 4], 
               col = col.edge, ...)
    }
  }
}