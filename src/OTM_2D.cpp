#include "utilities.h"
using namespace Rcpp;

void OTM2D(GEO::OptimalTransportMap2d &OTM, const arma::mat &X, double epsilon, int maxit, bool verbose) {
  const int n = X.n_rows;
  
  // create a mesh for data points
  GEO::Mesh dataMesh(2, false);
  dataMesh.vertices.create_vertices(n);
  setMeshPoint(dataMesh, X);
  
  // embed in 3-dimension
  dataMesh.vertices.set_dimension(3);
  
  // setup OTM
  OTM.set_points(n, dataMesh.vertices.point_ptr(0), dataMesh.vertices.dimension());
  
  // optmizer settings
  OTM.set_verbose(verbose);
  OTM.set_epsilon(epsilon);
  OTM.set_regularization(0.0);
  OTM.set_Newton(true);
  
  OTM.optimize(maxit);
}

//' 2D Semi-discrete Optimal Transport Map
//' 
//' @param X input data matrix.
//' @param epsilon maximum error for the optimization algorithm.
//' @param maxit maximum number of solver iterations.
//' @param verbose logical indicating whether to display optimization messages.
//' @return a list describing the resulting optimal transport map.
//' @keywords internal
// [[Rcpp::export]]
List dualGraphs2D(const arma::mat &X, double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(2, false);
  setUniMesh(unifMesh, 2, true);
  
  // compute OTM
  GEO::OptimalTransportMap2d OTM(&unifMesh);
  OTM2D(OTM, X, epsilon, maxit, verbose);
  
  // must save centroids before remesh
  arma::mat Centroid = getCentroids(OTM);
  
  // store weights
  arma::vec W = getWeights(OTM);
  
  // create a squared uniform mesh
  setUniMesh(unifMesh, 2, false);
  
  // this will not return the correct cell order
  // GEO::Mesh otmRVD;
  // OTM.get_RVD(otmRVD);
  
  // get the RDT
  GEO::Mesh otmRDT;
  OTM.RVD()->set_volumetric(false);
  OTM.RVD()->compute_RDT(otmRDT);
  
  // collect objects to return
  List lst;
  lst["Data"] = X;
  lst["Centroid"] = Centroid;
  lst["Weight"] = W;
  lst["Vertex.RDT"] = getVertices2D(otmRDT);
  lst["Vertex.RVD"] = getVerticesGen2D(unifMesh, OTM);
  lst["N.Triangles"] = otmRDT.facets.nb();
  lst["N.Cells"] = OTM.nb_points();
  
  return lst;
}

// [[Rcpp::export]]
List dualPotential2D(const arma::mat &Y, const arma::mat &X, const arma::mat &V, const arma::vec &h, 
                     const arma::uvec accuVerts) {
  const int m = Y.n_rows;
  const int n = X.n_rows;
  
  // seperate the vertices out for each cell
  // then shift each by its Voronoi site
  arma::field<arma::mat> cellVerts(n);
  for (int i = 0; i < n; i++) {
    cellVerts(i) = V.rows(accuVerts(i), accuVerts(i + 1) - 1);
  }
  
  // to store indices of the maximizers
  arma::ivec cellMaxInd(n);
  int psiInd;
  
  // to store the maximum values
  arma::vec cellMax(n, arma::fill::zeros);
  arma::vec psi(m, arma::fill::zeros);
  
  // to recover the optimal vertices
  // which corresponds to ranks
  arma::mat argmaxVert(m, 2, arma::fill::zeros);
  
  // search through the cells to find vertices that maximizes the inner product
  arma::vec tempProd;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      tempProd = cellVerts(j) * trans(Y.row(i) - X.row(j));
      cellMaxInd(j) = arma::index_max(tempProd);
      cellMax(j) = tempProd(cellMaxInd(j));
    }
    psiInd = arma::index_max(cellMax - h);
    psi(i) = cellMax(psiInd) - h(psiInd);
    argmaxVert.row(i) = cellVerts(psiInd).row(cellMaxInd(psiInd));
  }
  
  List lst;
  lst["optimal.vertex"] = argmaxVert;
  lst["dual.potential"] = psi;
  
  return lst;
}

//' Locate within a 2D RVD
//' 
//' Find the RVD cells where a set Q of query points is transported to.
//' @param Q input query matrix.
//' @param X input data used to build the RVD.
//' @param w weights of the RVD cells.
//' @return a vector of indices indicating where the query points are transported to.
//' @keywords internal
// [[Rcpp::export]]
arma::ivec locateRVD2D(const arma::mat &Q, const arma::mat &X, const arma::vec &w) {
  const int m = Q.n_rows;
  const int n = X.n_rows;
  
  // setup weighted vertices for X
  arma::vec wX(3 * n);
  getWeightedVerts(X, w, wX.memptr());
  
  // build a Kd-tree
  GEO::NearestNeighborSearch* treeX = GEO::NearestNeighborSearch::create(3, "BNN");
  treeX->set_points(n, wX.memptr());
  
  double p[3] = {};
  arma::ivec cellID(m, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    p[0] = Q(i, 0);
    p[1] = Q(i, 1);
    
    cellID(i) = treeX->get_nearest_neighbor(p) + 1;
  }
  
  return cellID;
}

//' Locate Points within a 2D RDT
//' 
//' Find the RDT triangles where a set Q of query points belongs to.
//' @param Q input query matrix.
//' @param V matrix of sorted triangle vertices.
//' @return a vector of indices indicating which triangles the query points belong.
//' @keywords internal
// [[Rcpp::export]]
arma::ivec locateRDT2D(const arma::mat &Q, const arma::mat &V) {
  const int m = Q.n_rows; // number of query points
  const int n = V.n_rows / 3; // number of triangles
  
  // vector recording the indices of triangles that contains the queries
  // 0 means q is not in any triangle, negative index means q is on the edge
  arma::ivec triID(m, arma::fill::zeros);
  
  // memoization for barycentric coordinates
  arma::mat memBaryCoord(n, 5);
  for (int j = 0; j < n; j++) {
    memBaryCoord(j, 0) = V(3 * j, 1) - V(3 * j + 2, 1);
    memBaryCoord(j, 1) = V(3 * j + 1, 1) - V(3 * j + 2, 1);
    memBaryCoord(j, 2) = V(3 * j, 0) - V(3 * j + 2, 0);
    memBaryCoord(j, 3) = V(3 * j + 2, 0) - V(3 * j + 1, 0);
    memBaryCoord(j, 4) = memBaryCoord(j, 1) * memBaryCoord(j, 2);
    memBaryCoord(j, 4) += memBaryCoord(j, 3) * memBaryCoord(j, 0);
  }
  
  arma::vec lambda(3, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    // check for interruption every 1000 iterations
    if (i % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    for (int j = 0; j < n; j++) {
      // compute the barycentric coordinates of q_i
      // first coordinate
      lambda(0) = memBaryCoord(j, 1) * (Q(i, 0) - V(3 * j + 2, 0));
      lambda(0) += memBaryCoord(j, 3) * (Q(i, 1) - V(3 * j + 2, 1));
      lambda(0) /= memBaryCoord(j, 4);
      
      // second coordinate
      lambda(1) = -memBaryCoord(j, 0) * (Q(i, 0) - V(3 * j + 2, 0));
      lambda(1) += memBaryCoord(j, 2) * (Q(i, 1) - V(3 * j + 2, 1));
      lambda(1) /= memBaryCoord(j, 4);
      
      // third coordinate
      lambda(2) = 1.0 - lambda(0) - lambda(1);
      
      // check the relative position
      if (all(lambda >= 0) && all(lambda <= 1)) {
        if (any(lambda == 0)) {
          // point q_i is on the edge
          triID(i) = -j - 1;
        } else {
          // point q_i is inside the triangle
          triID(i) = j + 1;
        }
        break; // didn't find anything
      }
    }
  }
  
  return triID;
}

//' 2D Goodness-of-fit Test Helper
//' 
//' Compute the SDOT quantiles of quasi-MC sequence \eqn{U} with respect to \eqn{X} and \eqn{Y}.
//' @param X input data matrix.
//' @param Y input data matrix.
//' @param U sequence used to evaluate the integral.
//' @param epsilon convergence threshold for optimization.
//' @param maxit max number of iterations before termination.
//' @param verbose logical indicating whether to display optimization messages.
//' @return a list containing the quantile indices for \eqn{X} and \eqn{Y}.
//' @keywords internal
// [[Rcpp::export]]
List gof2DHelper(const arma::mat &X, const arma::mat &Y, const arma::mat &U, 
                 double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(2, false);
  setUniMesh(unifMesh, 2, true);
  
  // compute OTMs
  GEO::OptimalTransportMap2d OTMX(&unifMesh);
  GEO::OptimalTransportMap2d OTMY(&unifMesh);
  OTM2D(OTMX, X, epsilon, maxit, verbose);
  OTM2D(OTMY, Y, epsilon, maxit, verbose);
  
  // get weights
  arma::vec wX = getWeights(OTMX);
  arma::vec wY = getWeights(OTMY);
  
  // get the optimal transport mapping
  arma::ivec uX = locateRVD2D(U, X, wX);
  arma::ivec uY = locateRVD2D(U, Y, wY);
  
  // collect objects to return
  List lst;
  lst["U_Map_X"] = uX;
  lst["U_Map_Y"] = uY;
  
  return lst;
}

//' 2D Joint Samples Rank Helper
//' 
//' Helps computing the Empirical Rank for 2D Joint Samples \eqn{(X, Y)}.
//' @param XY 2D input data matrix.
//' @param center logical indicating if the centroids should be computed.
//' @param epsilon convergence threshold for optimization.
//' @param maxit max number of iterations before termination.
//' @param verbose logical indicating whether to display optimization messages.
//' @return a matrix, represents either the centroids/cells of the RVD of joint samples.
//' @keywords internal
// [[Rcpp::export]]
arma::mat jointRankHelper2D(const arma::mat &XY, bool center, double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(2, false);
  setUniMesh(unifMesh, 2, true);
  
  // embed in 3-dimension
  unifMesh.vertices.set_dimension(3);
  
  // compute OTMs
  GEO::OptimalTransportMap2d OTMXY(&unifMesh);
  OTM2D(OTMXY, XY, epsilon, maxit, verbose);
  
  // get elements from corresponding cells to help computing ranks
  arma::mat elemXY;
  if (center) {
    elemXY = getCentroids(OTMXY);
  } else {
    // create a squared uniform mesh
    setUniMesh(unifMesh, 2, false);
    
    // get Voronoi cells
    elemXY = getVerticesGen2D(unifMesh, OTMXY);
  }
  
  return elemXY;
}