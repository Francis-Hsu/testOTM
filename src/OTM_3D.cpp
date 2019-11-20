#include "utilities.h"
using namespace Rcpp;

void OTM3D(GEO::OptimalTransportMap3d &OTM, const arma::mat &X, double epsilon, int maxit, bool verbose) {
  const int n = X.n_rows;
  
  // create a mesh for data points
  GEO::Mesh dataMesh(3, false);
  dataMesh.vertices.create_vertices(n);
  setMeshPoint(dataMesh, X);
  
  // embed in 4-dimension
  dataMesh.vertices.set_dimension(4);
  
  // setup OTM
  OTM.set_points(n, dataMesh.vertices.point_ptr(0), dataMesh.vertices.dimension());
  
  // optmizer settings
  OTM.set_verbose(verbose);
  OTM.set_epsilon(epsilon);
  OTM.set_regularization(0.0);
  OTM.set_Newton(true);
  OTM.optimize(maxit);
}

//' 3D Semi-discrete Optimal Transport Map
//' 
//' @param X input data matrix.
//' @param epsilon maximum error for the optimization algorithm.
//' @param maxit maximum number of solver iterations.
//' @param verbose logical indicating whether to display optimization messages.
//' @return a list describing the resulting optimal transport map.
//' @keywords internal
// [[Rcpp::export]]
List dualGraphs3D(const arma::mat &X, double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 3D uniform measure
  GEO::Mesh unifMesh(3, false);
  setUniMesh(unifMesh, 3, true);
  
  // compute OTM
  GEO::OptimalTransportMap3d OTM(&unifMesh);
  OTM3D(OTM, X, epsilon, maxit, verbose);
  
  // must save centroids before remesh
  arma::mat Centroid = getCentroids(OTM);
  
  // store weights
  arma::vec W = getWeights(OTM);
  
  // get the RDT
  GEO::Mesh otmRDT;
  OTM.RVD()->compute_RDT(otmRDT);
  
  // compute the RVCs and save the vertices to a field
  GEO::Mesh otmRVD;
  arma::field<arma::mat> rvcVerts(OTM.nb_points());
  setUniMesh(unifMesh, 3, false);
  for (unsigned int i = 0; i < OTM.nb_points(); i++) {
    otmRVD.clear();
    OTM.RVD()->compute_RVC(i, unifMesh, otmRVD);
    rvcVerts(i) = getVertices3D(otmRVD, false);
  }
  
  // save the RVD to a mesh
  // set integration_simplices to true so that 
  // the tetrahedra have the data as the first vertex
  // GEO::Mesh otmRVD;
  // GEO::Attribute<GEO::index_t> tet_region(otmRVD.cells.attributes(), "region");
  // OTM.RVD()->compute_RVD(otmRVD, 0, false, true);
  
  // collect objects to return
  List lst;
  lst["Data"] = X;
  lst["Centroid"] = Centroid;
  lst["Weight"] = W;
  lst["Vertex.RDT"] = getVertices3D(otmRDT, true);
  lst["Vertex.RVD"] = rvcVerts;
  lst["N.Triangles"] = otmRDT.cells.nb();
  lst["N.Cells"] = OTM.nb_points();
  
  return lst;
}

// [[Rcpp::export]]
List dualPotential3D(const arma::mat &Y, const arma::mat &X, const arma::mat &V, const arma::vec &h, 
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
  arma::mat argmaxVert(m, 3, arma::fill::zeros);
  
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


//' Locate within a 3D RVD
//' 
//' Find the RVD seeds where a set Q of query points is transported to.
//' @param Q input query matrix.
//' @param X input data used to build the RVD.
//' @param w weights of the RVD cells.
//' @return a vector of indices indicating where the query points are transported to.
//' @keywords internal
// [[Rcpp::export]]
arma::ivec locateRVD3D(const arma::mat &Q, const arma::mat &X, const arma::vec &w) {
  const int m = Q.n_rows;
  const int n = X.n_rows;
  
  // setup weighted vertices for X
  arma::vec wX(4 * n);
  getWeightedVerts(X, w, wX.memptr());
  
  // build a Kd-tree
  GEO::NearestNeighborSearch* treeX = GEO::NearestNeighborSearch::create(4, "BNN");
  treeX->set_points(n, wX.memptr());
  
  double p[4] = {};
  arma::ivec cellID(m, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    p[0] = Q(i, 0);
    p[1] = Q(i, 1);
    p[2] = Q(i, 2);
    
    cellID(i) = treeX->get_nearest_neighbor(p) + 1;
  }
  
  return cellID;
}

//' 3D Goodness-of-fit Test Helper
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
List gof3DHelper(const arma::mat &X, const arma::mat &Y, const arma::mat &U, 
                 double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(3, false);
  setUniMesh(unifMesh, 3, true);
  
  // compute OTMs
  GEO::OptimalTransportMap3d OTMX(&unifMesh);
  GEO::OptimalTransportMap3d OTMY(&unifMesh);
  OTM3D(OTMX, X, epsilon, maxit, verbose);
  OTM3D(OTMY, Y, epsilon, maxit, verbose);
  
  // get weights
  arma::vec wX = getWeights(OTMX);
  arma::vec wY = getWeights(OTMY);
  
  // get the optimal transport mapping
  arma::ivec uX = locateRVD3D(U, X, wX);
  arma::ivec uY = locateRVD3D(U, Y, wY);
  
  // collect objects to return
  List lst;
  lst["U_Map_X"] = uX;
  lst["U_Map_Y"] = uY;
  
  return lst;
}

//' 3D Joint Samples Rank Helper
//' 
//' Helps computing the Empirical Rank for 3D Joint Samples \eqn{(X, Y)}.
//' @param XY 3D input data matrix.
//' @param center logical indicating if the centroids should be computed.
//' @param epsilon convergence threshold for optimization.
//' @param maxit max number of iterations before termination.
//' @param verbose logical indicating whether to display optimization messages.
//' @return a matrix, represents either the centroids or the cells of the RVD of joint samples.
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::mat> jointRankHelper3D(const arma::mat &XY, bool center, double epsilon, int maxit, bool verbose) {
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(3, false);
  setUniMesh(unifMesh, 3, true);
  
  // embed in 4-dimension
  unifMesh.vertices.set_dimension(4);
  
  // compute OTMs
  GEO::OptimalTransportMap3d OTMXY(&unifMesh);
  OTM3D(OTMXY, XY, epsilon, maxit, verbose);
  
  // get elements from corresponding cells to help computing ranks
  arma::field<arma::mat> elemXY;
  if (center) {
    elemXY.set_size(1);
    elemXY(0) = getCentroids(OTMXY);
  } else {
    elemXY.set_size(OTMXY.nb_points());
    // compute the RVCs and save the vertices to a field
    GEO::Mesh otmRVD;
    setUniMesh(unifMesh, 3, false);
    for (unsigned int i = 0; i < OTMXY.nb_points(); i++) {
      otmRVD.clear();
      OTMXY.RVD()->compute_RVC(i, unifMesh, otmRVD);
      elemXY(i) = getVertices3D(otmRVD, false);
    }
  }
  
  return elemXY;
}