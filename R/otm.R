#' Uniform Semi-discrete Optimal Transport Map
#'
#' \code{otm.fit} computes the semi-discrete optimal transport map from the \eqn{U[0, 1]^d} measure to the input data set.
#' @param data input data matrix, of size \eqn{n} by \eqn{d}.
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating wether to display optimization messages.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return \code{otm.fit} returns an object of class "\code{otm.2d}" or "\code{otm.3d}", depending on the dimesion of input data.
#' An object of class "\code{otm}" is a list describing the resulting optimal transport map.
#' @keywords optimize, graphs
#' @importFrom stats complete.cases
#' @export
otm.fit = function(data,
                   epsilon = 1e-3,
                   maxit = 100,
                   verbose = F,
                   na.rm = F) {
  # input validation
  if (!is.matrix(data) || ncol(data) < 2) {
    stop("Data must be a matrix with ncol >= 2.")
  }
  
  if (na.rm) {
    data = data[complete.cases(data),]
  }
  
  if (min(data) < 0 || max(data) > 1) {
    stop("Data must be scaled to [0, 1] range.")
  }
  
  d = ncol(data)
  if (d == 2) {
    object = dualGraphs2D(data, epsilon, maxit, verbose)
    class(object) = "otm.2d"
    colnames(object$Centroid) = c("x", "y")
    colnames(object$Vertex.RDT) = c("cell", "x", "y", "id")
    colnames(object$Vertex.RVD) = c("cell", "x", "y")
  } else if (d == 3) {
    stop("3-dimension are not supported currently.")
  } else {
    stop("Dimension higher than 3 are not supported.")
  }
  object$Vertex.RDT = as.data.frame(object$Vertex.RDT)
  object$Vertex.RVD = as.data.frame(object$Vertex.RVD)
  
  # for computing Alexandrov's potential
  # noted the signs of weights
  object$Height = -(rowSums(object$Data ^ 2) - object$Weight) / 2
  
  return(object)
}

#' Plot the 2D Semi-discrete Optimal Transport Map
#'
#' Plots the restricted Voronoi diagram (RVD) and the restricted Delaunay triangulation (RDT) of a given
#' 2D semi-discrete optimal transport map.
#' @param x a fitted \code{otm.2d} object.
#' @param which specify which graph(s) to plot. Can be "\code{RVD}", "\code{RDT}", or "\code{Both}".
#' @param col.data color of the data points.
#' @param col.center color of the Voronoi centroids.
#' @param col.edge color of the edges in plotting RVD and RDT.
#' @param draw.center logical indicating if the centroids should be plotted.
#' @param draw.map logical indicating if dashed lines should be added to show mapping between the data and the Voronoi cells.
#' @param \dots other graphical parameters to plot.
#' @keywords hplot
#' @importFrom graphics plot.default segments points
#' @export
plot.otm.2d = function(x,
                       which = "Both",
                       col.data = "cornflowerblue",
                       col.center = "firebrick",
                       col.edge = "black",
                       draw.center = T,
                       draw.map = F,
                       ...) {
  type = match.arg(which, c("RVD", "RDT", "Both"))
  
  # plot the restricted Voronoi diagram
  if (which == "RVD" || which == "Both") {
    # plot data
    plot.default(
      x$Data[, 1],
      x$Data[, 2],
      xlim = c(0, 1),
      ylim = c(0, 1),
      col = col.data,
      pch = 20,
      xlab = expression('u'[1]),
      ylab = expression('u'[2]),
      ...
    )
    
    # plot Laguerre cells
    n.edge = nrow(x$Vertex.RVD)
    rvd.edges = matrix(0, n.edge, 4)
    curr.row = 1
    
    # extract edges
    for (i in 1:x$N.Cells) {
      curr.cell = subset(x$Vertex.RVD, x$Vertex.RVD[, 1] == i, select = 2:3)
      for (j in 1:nrow(curr.cell)) {
        p.id = j
        q.id = 1 + j %% nrow(curr.cell)
        
        # sort by l1 norm, everything positive!
        if (curr.cell[p.id, 1] + curr.cell[p.id, 2] <= curr.cell[q.id, 1] + curr.cell[q.id, 2]) {
          rvd.edges[curr.row,] = c(curr.cell[p.id, 1], curr.cell[p.id, 2],
                                   curr.cell[q.id, 1], curr.cell[q.id, 2])
        } else {
          rvd.edges[curr.row,] = c(curr.cell[q.id, 1], curr.cell[q.id, 2],
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
      points(
        x$Centroid[, 1],
        x$Centroid[, 2],
        xlim = c(0, 1),
        ylim = c(0, 1),
        pch = 20,
        col = col.center,
        ...
      )
    }
    
    # plot mappings
    if (draw.map) {
      for (i in 1:nrow(x$Data)) {
        segments(
          x$Data[i, 1],
          x$Data[i, 2],
          x$Centroid[i, 1],
          x$Centroid[i, 2],
          lty = 2,
          lwd = 0.5
        )
      }
    }
  }
  
  # plot the restricted Delaunay triangulation
  if (which == "RDT" || which == "Both") {
    # plot data
    plot.default(
      x$Data[, 1],
      x$Data[, 2],
      xlim = c(0, 1),
      ylim = c(0, 1),
      col = col.data,
      pch = 20,
      xlab = expression('u'[1]),
      ylab = expression('u'[2]),
      ...
    )
    
    # plot triangles
    n.edge = nrow(x$Vertex.RDT)
    rdt.edges = matrix(0, n.edge, 4)
    curr.row = 1
    
    # extract edges
    for (i in 1:x$N.Triangles) {
      curr.cell = subset(x$Vertex.RDT, x$Vertex.RDT[, 1] == i, select = 2:3)
      for (j in 1:nrow(curr.cell)) {
        p.id = j
        q.id = 1 + j %% nrow(curr.cell)
        
        # sort by l1 norm, everything positive!
        if (curr.cell[p.id, 1] + curr.cell[p.id, 2] <= curr.cell[q.id, 1] + curr.cell[q.id, 2]) {
          rdt.edges[curr.row,] = c(curr.cell[p.id, 1], curr.cell[p.id, 2],
                                   curr.cell[q.id, 1], curr.cell[q.id, 2])
        } else {
          rdt.edges[curr.row,] = c(curr.cell[q.id, 1], curr.cell[q.id, 2],
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