% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otm.R
\name{tos.fit}
\alias{tos.fit}
\title{Fitting Uniform Semi-discrete Optimal Transport Map}
\usage{
tos.fit(
  data,
  scale = c(0.05, 0.95),
  na.rm = FALSE,
  epsilon = 1e-06,
  maxit = 100,
  verbose = F
)
}
\arguments{
\item{data}{a numeric matrix for the input data, of size \eqn{n} by \eqn{2} or \eqn{3}.}

\item{scale}{a numeric vector indicating the minimum and maximum of the scaled data. Set to \code{NULL} to skip scaling.}

\item{na.rm}{logical indicating whether \code{NA} values should be stripped before the computation proceeds.}

\item{epsilon}{convergence threshold for optimization.}

\item{maxit}{max number of iterations before termination.}

\item{verbose}{logical indicating whether to display optimization messages.}
}
\value{
\code{tos.fit} returns an object of class \code{tos.2d} or \code{tos.3d}, depending on the dimension of input data.
An object of class \code{tos.2d} or \code{tos.3d} is a list describing the optimal transport map. It contains the following elements:
\item{Data}{Scaled data.}
\item{Centroid}{Centroid of the Laguerre cells.}
\item{Weight}{Weights of the Laguerre cells.}
\item{Vertex.RDT}{Vertices of the restricted Delaunay triangulation (RDT), organized by cells.}
\item{Vertex.RVD}{Vertices of the restricted Voronoi diagram (RVD), organized by cells and facets.}
\item{N.Triangles}{Number of tetrahedra in the RDT.}
\item{N.Cells}{Number of Laguerre cells in the RVD.}
\item{Height}{Height of the Laguerre cells, used for computing Alexandrov's potential.}
\item{Location}{The numeric centering applied to the data.}
\item{Scale}{The numeric scalings applied to the data.}
}
\description{
\code{tos.fit} computes the semi-discrete optimal transport map from the \eqn{U[0, 1]^d} measure to the input data set.
}
\details{
Input data needs be scaled within the [0, 1] range for computation. 
Supply \code{scale} to let \code{tos.fit} handles the scaling internally.
Centering and scaling factors will be returned to help transforming the results back to their original range.
}
\examples{
# generate some data
X = c(0.5, 0.8, -0.2, -1.5, 1.4, 0.5, -1.1, -0.1, -1.1, -2.6)
Y = c(1.6, -1.0, -0.1, 0.5, -1.3, 2.9, -0.4, 1.3, -1.8, -2.5)

# compute the optimal transport map from U[0, 1]^2 to the data
XY.OTM = tos.fit(cbind(X, Y))
}
\references{
Bruno Levy and Erica L. Schwindt (2018).
\emph{Notions of Optimal Transport Theory and How to Implement Them on a Computer}.
Computers & Graphics 72: 135-48.
\url{https://doi.org/10.1016/j.cag.2018.01.009}.
}
\keyword{graphs}
\keyword{optimize}
