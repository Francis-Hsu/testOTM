#' Fitting Uniform Semi-discrete Optimal Transport Map
#'
#' \code{tos.fit} computes the semi-discrete optimal transport map from the \eqn{U[0, 1]^d} measure to the input data set.
#' @param data a numeric matrix for the input data, of size \eqn{n} by \eqn{2} or \eqn{3}.
#' @param scale a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @param epsilon convergence threshold for optimization.
#' @param maxit max number of iterations before termination.
#' @param verbose logical indicating whether to display optimization messages.
#' @details Input data needs be scaled within the [0, 1] range for computation. 
#' Supply \code{scale} to let \code{tos.fit} handles the scaling internally.
#' Centering and scaling factors will be returned to help transforming the results back to their original range.
#' @return \code{tos.fit} returns an object of class \code{tos.2d} or \code{tos.3d}, depending on the dimension of input data.
#' An object of class \code{tos.2d} or \code{tos.3d} is a list describing the optimal transport map. It contains the following elements:
#' \item{Data}{Scaled data.}
#' \item{Centroid}{Centroid of the Laguerre cells.}
#' \item{Weight}{Weights of the Laguerre cells.}
#' \item{Vertex.RDT}{Vertices of the restricted Delaunay triangulation (RDT), organized by cells.}
#' \item{Vertex.RVD}{Vertices of the restricted Voronoi diagram (RVD), organized by cells and facets.}
#' \item{N.Triangles}{Number of tetrahedra in the RDT.}
#' \item{N.Cells}{Number of Laguerre cells in the RVD.}
#' \item{Height}{Height of the Laguerre cells, used for computing Alexandrov's potential.}
#' \item{Location}{The numeric centering applied to the data.}
#' \item{Scale}{The numeric scalings applied to the data.}
#' @examples
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' @keywords optimize graphs
#' @references Bruno Levy and Erica L. Schwindt (2018).
#' \emph{Notions of Optimal Transport Theory and How to Implement Them on a Computer}.
#' Computers & Graphics 72: 135-48.
#' \url{https://doi.org/10.1016/j.cag.2018.01.009}.
#' @importFrom stats complete.cases
#' @export
tos.fit = function(data,
                   scale = c(0.05, 0.95),
                   na.rm = FALSE,
                   epsilon = 1e-6,
                   maxit = 100,
                   verbose = F) {
  # input validation
  if (!is.matrix(data) || NCOL(data) < 2) {
    stop("Input data must be a matrix with ncol >= 2.")
  }
  
  if (!is.null(scale)) {
    if (scale[1] < 0 || scale[2] > 1) {
      stop("Scaling range must be within [0, 1].")
    }
    
    data = scaling.min.max(data, scale[1], scale[2])
  }
  
  if (min(data, na.rm = T) < 0 || max(data, na.rm = T) > 1) {
    stop("Data must be within the [0, 1] range.")
  }
  
  if (na.rm) {
    data = data[complete.cases(data),]
  }
  
  if (min(data) < 0 || max(data) > 1) {
    stop("Data must be scaled to [0, 1] range.")
  }
  
  d = NCOL(data)
  if (d == 2) {
    if (NROW(data) < 3) {
      stop("Too few data, must have more than 3 points for 2D transport.")
    }
    object = dualGraphs2D(data, epsilon, maxit, verbose)
    class(object) = "tos.2d"
    colnames(object$Centroid) = c("x", "y")
    colnames(object$Vertex.RDT) = c("cell", "id", "x", "y")
    colnames(object$Vertex.RVD) = c("cell", "x", "y")
  } else if (d == 3) {
    if (NROW(data) < 4) {
      stop("Too few data, must have more than 4 points for 3D transport.")
    }
    object = dualGraphs3D(data, epsilon, maxit, verbose)
    class(object) = "tos.3d"
    colnames(object$Centroid) = c("x", "y", "z")
    colnames(object$Vertex.RDT) = c("cell", "facet", "id", "x", "y", "z")
    object$Vertex.RVD = do.call(rbind, lapply(seq_along(object$Vertex.RVD), function(i) {
      cbind(i, object$Vertex.RVD[[i]])
    }))
    colnames(object$Vertex.RVD) = c("cell", "facet", "id", "x", "y", "z")
  } else {
    stop("Dimension higher than 3 are not supported.")
  }
  
  # for computing Alexandrov's potential
  # note the sign of weights
  object$Height = -(rowSums(object$Data^2) - object$Weight) / 2
  
  if (!is.null(scale)) {
    object$Location = attr(data, "scaled:center")
    object$Scale = c(attr(data, "scaled:scale"))
  }
  
  return(object)
}

#' Plotting 2D Semi-discrete Optimal Transport Map
#'
#' \code{plot.tos.2d} plots the restricted Voronoi diagram (RVD) and the restricted Delaunay triangulation (RDT) of a given
#' 2D semi-discrete optimal transport map.
#' @param x a fitted \code{tos.2d} object.
#' @param which specify which graph(s) to plot. Can be \code{None}, \code{RVD}, \code{RDT}, or \code{Both}.
#' @param col.data color of the data points.
#' @param col.center color of the Voronoi centroids.
#' @param col.edge color of the edges in plotting RVD and RDT.
#' @param draw.data logical indicating if the data points should be plotted.
#' @param draw.center logical indicating if the centroids should be plotted.
#' @param draw.map logical indicating if dashed lines should be added to show mapping between the data and the Voronoi cells.
#' @param xlim the x limits of the plot.
#' @param ylim the y limits of the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param pch a vector of plotting characters or symbols. 
#' @param \dots other graphical parameters to plot.
#' @examples 
#' # generate some data
#' X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
#' Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)
#' 
#' # compute the optimal transport map from U[0, 1]^2 to the data
#' XY.OTM = tos.fit(cbind(X, Y))
#' 
#' # plot the RVD
#' plot(XY.OTM, which = "RVD", draw.map = TRUE)
#' @keywords hplot
#' @importFrom graphics plot.default segments points
#' @export
plot.tos.2d = function(x,
                       which = "Both",
                       col.data = "dodgerblue3",
                       col.center = "darkorange1",
                       col.edge = "black",
                       draw.data = TRUE,
                       draw.center = TRUE,
                       draw.map = FALSE,
                       xlim = c(0, 1),
                       ylim = c(0, 1),
                       xlab = expression('x'),
                       ylab = expression('y'),
                       pch = 20,
                       ...) {
  # validate plot type
  type = match.arg(which, c("None", "RVD", "RDT", "Both"))
  
  # plot thins other than RVD and RDT
  if (type != "RDT") {
    # plot data
    plot.default(
      x$Data[, 1],
      x$Data[, 2],
      type = ifelse(draw.data, "p", "n"),
      col = col.data,
      pch = pch,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      ...
    )
    
    # plot centroids
    if (draw.center) {
      points(
        x$Centroid[, 1],
        x$Centroid[, 2],
        col = col.center,
        pch = pch,
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
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
  
  # plot the restricted Voronoi diagram
  if (type == "RVD" || type == "Both") {
    rvd.edges = extract.segment(x$Vertex.RVD, 1, 2:3)
    for (i in 1:nrow(rvd.edges)) {
      segments(rvd.edges[i, 1], rvd.edges[i, 2],
               rvd.edges[i, 3], rvd.edges[i, 4],
               col = col.edge, ...)
    }
  }
  
  # plot the restricted Delaunay triangulation
  if (type == "RDT" || type == "Both") {
    # plot data
    plot.default(
      x$Data[, 1],
      x$Data[, 2],
      type = ifelse(draw.data, "p", "n"),
      col = col.data,
      pch = pch,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      ...
    )
    
    rdt.edges = extract.segment(x$Vertex.RDT, 1, 3:4)
    for (i in 1:nrow(rdt.edges)) {
      segments(rdt.edges[i, 1], rdt.edges[i, 2],
               rdt.edges[i, 3], rdt.edges[i, 4],
               col = col.edge, ...)
    }
  }
}

#' Plotting the 3D Semi-discrete Optimal Transport Map
#'
#' \code{plot.tos.3d} plots the restricted Voronoi diagram (RVD) and the restricted Delaunay triangulation (RDT) of a given
#' 3D semi-discrete optimal transport map.
#' @param x a fitted \code{tos.3d} object.
#' @param inter logical indicating if the plot should be interactive.
#' @param which specify which graph(s) to plot. Can be \code{None}, \code{RVD}, \code{RDT}, or \code{Both}.
#' @param col.data color of the data points.
#' @param col.center color of the Voronoi centroids.
#' @param col.rvd color of the edges in RVD plot.
#' @param col.rdt color of the edges in RDT plot.
#' @param draw.data logical indicating if the data points should be plotted.
#' @param draw.center logical indicating if the centroids should be plotted.
#' @param draw.map logical indicating if lines should be added to show mapping between the data and the Voronoi cells.
#' @param draw.id vector of indices specifying which data points (and the corresponding cells in RVD and RDT) will be plotted. 
#' Set to \code{NULL} to plot all points. If it is not \code{NULL} and RDT is to be drawn, 
#' \code{plot.tos.3d} will plot all tetrahedra that have points in \code{draw.id} as one of its vertices.
#' @param xlim the x limits of the plot.
#' @param ylim the y limits of the plot.
#' @param zlim the z limits of the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param zlab a label for the z axis.
#' @param size size of the plotted points in the interactive plots. 
#' @param theta,phi the angles defining the viewing direction for non-interactive plots. 
#' \code{theta} gives the azimuthal direction and \code{phi} the colatitude.
#' @param pch a vector of plotting characters or symbols.
#' @param \dots other graphical parameters to plot.
#' @examples 
#' # generate some data
#' X = c(2.3, 0.8, 1.2, 0.6, 0.2, 0.5, -0.7, -0.6, -0.2, 3.2)
#' Y = c(1.5, -0.8, 0.6, -2.1, 1.3, 1.3, -1.1, 1.0, -0.2, 1.0)
#' Z = c(0.8, -0.3, -1.3, -1.8, 1.2, -0.7, -0.8, 1.3, -0.9, 2.0)
#' 
#' # compute the optimal transport map from U[0, 1]^3 to the data
#' XYZ.OTM = tos.fit(cbind(X, Y, Z))
#' 
#' # plot the (non-interactive) RVD
#' plot(XYZ.OTM, interactive = FALSE, which = "RVD", draw.map = TRUE)
#' @keywords hplot
#' @importFrom plot3D points3D segments3D
#' @importFrom rgl clear3d plot3d points3d segments3d
#' @export
plot.tos.3d = function(x,
                       inter = interactive(),
                       which = "Both",
                       col.data = "dodgerblue3",
                       col.center = "darkorange1",
                       col.rvd = "black",
                       col.rdt = "firebrick3",
                       draw.data = TRUE,
                       draw.center = TRUE,
                       draw.map = FALSE,
                       draw.id = NULL,
                       xlim = c(0, 1),
                       ylim = c(0, 1),
                       zlim = c(0, 1),
                       xlab = expression('x'),
                       ylab = expression('y'),
                       zlab = expression('z'),
                       size = 10,
                       phi = 45,
                       theta = 45,
                       pch = 20,
                       ...) {
  # validate plot type
  type = match.arg(which, c("None", "RVD", "RDT", "Both"))
  
  # indices of elements to plot
  if (is.null(draw.id)) {
    id.draw = 1:nrow(x$Data)
  } else {
    id.draw = draw.id
  }
  
  if (inter) {
    # clear exisiting plot, if any
    clear3d(type = "all")
    
    # plot data
    plot3d(
      x$Data[id.draw, , drop = FALSE],
      type = ifelse(draw.data, "p", "n"),
      col = col.data,
      pch = pch,
      size = size,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      ...
    )
    
    # plot centroids
    if (draw.center) {
      points3d(x$Centroid[id.draw, , drop = FALSE],
               col = col.center,
               pch = pch,
               size = size,
               ...)
    }
    
    # plot mappings
    if (draw.map) {
      for (i in id.draw) {
        segments3d(rbind(x$Data[i, ], x$Centroid[i, ]), col = "gray50", lwd = 1.0)
      }
    }
    
    # plot RVD
    if (type == "RVD" || type == "Both") {
      rvd.edges = lapply(id.draw, function(i) {
        # extract segments cell by cell
        extract.segment(subset(x$Vertex.RVD, x$Vertex.RVD[, 1] == i), 2, 4:6) 
      })
      rvd.edges = do.call(rbind, rvd.edges)
      rvd.edges = matrix(t(rvd.edges), ncol = 3, byrow = T)
      segments3d(rvd.edges, col = col.rvd, ...)
    }
    
    # plot RDT
    if (type == "RDT" || type == "Both") {
      # plot RVD segments
      if (is.null(draw.id)) {
        # draw all the cells in RDT
        id.draw = unique(x$Vertex.RDT[, 1])
      } else {
        # get the cells that contain vertices in the draw list
        id.draw = unique(subset(x$Vertex.RDT, x$Vertex.RDT[, 6] %in% id.draw)[, 1])
      }
      rdt.edges = lapply(id.draw, function(i) {
        # extract segments cell by cell
        extract.segment(subset(x$Vertex.RDT, x$Vertex.RDT[, 1] == i), 2, 4:6) 
      })
      rdt.edges = do.call(rbind, rdt.edges)
      rdt.edges = matrix(t(rdt.edges), ncol = 3, byrow = T)
      segments3d(rdt.edges, col = col.rdt, ...)
    }
  } else {
    # plot data
    points3D(
      x = x$Data[id.draw, 1, drop = FALSE],
      y = x$Data[id.draw, 2, drop = FALSE],
      z = x$Data[id.draw, 3, drop = FALSE],
      alpha = ifelse(draw.data, 1, 0),
      col = col.data,
      pch = pch,
      phi = phi,
      theta = theta,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      ...
    )
    
    # plot centroids
    if (draw.center) {
      points3D(
        x = x$Centroid[id.draw, 1, drop = FALSE],
        y = x$Centroid[id.draw, 2, drop = FALSE],
        z = x$Centroid[id.draw, 3, drop = FALSE],
        col = col.center,
        pch = pch,
        add = T,
        ...
      )
    }
    
    # plot mappings
    if (draw.map) {
      for (i in id.draw) {
        segments3D(
          x0 = x$Data[i, 1, drop = F],
          y0 = x$Data[i, 2, drop = F],
          z0 = x$Data[i, 3, drop = F],
          x1 = x$Centroid[i, 1, drop = F],
          y1 = x$Centroid[i, 2, drop = F],
          z1 = x$Centroid[i, 3, drop = F],
          lty = 2,
          lwd = 1.0,
          add = T
        )
      }
    }
    
    # plot RVD
    if (type == "RVD" || type == "Both") {
      # plot RVD segments
      rvd.edges = lapply(id.draw, function(i) {
        # extract segments cell by cell
        extract.segment(subset(x$Vertex.RVD, x$Vertex.RVD[, 1] == i), 2, 4:6) 
      })
      rvd.edges = do.call(rbind, rvd.edges)
      segments3D(
        rvd.edges[, 1],
        rvd.edges[, 2],
        rvd.edges[, 3],
        rvd.edges[, 4],
        rvd.edges[, 5],
        rvd.edges[, 6],
        col = col.rvd,
        add = T,
        ...
      )
    }
    
    # plot RDT
    if (type == "RDT" || type == "Both") {
      # plot RVD segments
      if (is.null(draw.id)) {
        # draw all the cells in RDT
        id.draw = unique(x$Vertex.RDT[, 1])
      } else {
        # get the cells that contain vertices in the draw list
        id.draw = unique(subset(x$Vertex.RDT, x$Vertex.RDT[, 6] %in% id.draw)[, 1])
      }
      rdt.edges = lapply(id.draw, function(i) {
        # extract segments cell by cell
        extract.segment(subset(x$Vertex.RDT, x$Vertex.RDT[, 1] == i), 2, 4:6) 
      })
      rdt.edges = do.call(rbind, rdt.edges)
      segments3D(
        rdt.edges[, 1],
        rdt.edges[, 2],
        rdt.edges[, 3],
        rdt.edges[, 4],
        rdt.edges[, 5],
        rdt.edges[, 6],
        col = col.rdt,
        add = T,
        ...
      )
    }
  }
}
  