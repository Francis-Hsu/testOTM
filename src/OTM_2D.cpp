#include "utilities.h"
using namespace Rcpp;

void OTM2D(GEO::OptimalTransportMap2d &OTM, const arma::mat &X, double epsilon, int maxit, bool verbose) {
  int n = X.n_rows;
  int d = 2;
  
  // create a mesh for data points
  GEO::Mesh dataMesh(d, false);
  dataMesh.vertices.create_vertices(n);
  setMeshPoint(dataMesh, X);
  
  // embed in 3-dimension
  dataMesh.vertices.set_dimension(d + 1);
  
  // setup OTM
  OTM.set_points(n, dataMesh.vertices.point_ptr(0), dataMesh.vertices.dimension());
  
  // optmizer settings
  OTM.set_verbose(verbose);
  OTM.set_epsilon(epsilon);
  OTM.set_regularization(0.0);
  OTM.set_Newton(true);
  
  OTM.optimize(maxit);
}

// [[Rcpp::export]]
List dualGraphs2D(const arma::mat &X, double epsilon, int maxit, bool verbose) {
  int d = 2;
  
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const arma::mat cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.n_rows);
  setMeshPoint(unifMesh, cubeVertices);
  
  // must triangulate for OTM
  unifMesh.facets.create_triangle(0, 1, 2);
  unifMesh.facets.create_triangle(2, 1, 3);
  unifMesh.facets.connect();
  
  // embed in 3-dimension
  unifMesh.vertices.set_dimension(d + 1);
  
  // compute OTM
  GEO::OptimalTransportMap2d OTM(&unifMesh);
  OTM2D(OTM, X, epsilon, maxit, verbose);
  
  // must save centroids before remesh
  arma::mat Centroid = getCentroids(OTM);

  // store weights
  arma::vec W = getWeights(OTM);

  // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.n_rows);
  setMeshPoint(unifMesh, cubeVertices);
  unifMesh.facets.create_quad(0, 2, 3, 1);
  
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
  lst["Vertex.RDT"] = getVertices(otmRDT);
  lst["Vertex.RVD"] = getVertices(unifMesh, OTM);
  lst["N.Triangles"] = otmRDT.facets.nb();
  lst["N.Cells"] = OTM.nb_points();
  
  return lst;
}

// under construction
// [[Rcpp::export]]
void OTMRank2D(const arma::mat &X, arma::vec &weight, double wMax, const arma::mat &Q) {
  int n = X.n_rows;
  int d = 2;
  
  // initialize the Geogram library.
  initializeGeogram();
  
  double wV[(d + 1) * n];
  getWeightedVerts(X, weight, wV);
  
  GEO::RegularWeightedDelaunay2d* otmRWD = dynamic_cast<GEO::RegularWeightedDelaunay2d*>(GEO::Delaunay::create(d + 1, "BPOW2d"));
  //otmRWD->set_keeps_infinite(true);
  otmRWD->set_vertices(n, wV);
  
  GEO::index_t cellID;
  Rcout << "#########" << std::endl;
  double q[2];
  for (unsigned int i = 0; i < Q.n_rows; i++) {
    q[0] = Q(i, 0);
    q[1] = Q(i, 1);
    
    cellID = otmRWD->locate(q);
    Rcout << cellID << std::endl;
    //Rcout << otmRWD->nearest_vertex(q) << std::endl;
    Rcout << X(otmRWD->cell_vertex(cellID, 0), 0) << " " << X(otmRWD->triangle_vertex(cellID, 0), 1) << std::endl;
    Rcout << X(otmRWD->cell_vertex(cellID, 1), 0) << " " << X(otmRWD->triangle_vertex(cellID, 1), 1) << std::endl;
    Rcout << X(otmRWD->cell_vertex(cellID, 2), 0) << " " << X(otmRWD->triangle_vertex(cellID, 2), 1) << std::endl;
  }
  
  Rcout << "#####" << std::endl;
  for (unsigned int i = 0; i < otmRWD->nb_cells(); i++) {
    //Rcout << otmRWD->triangle_adjacent(i, 0) << " " << 
    //  otmRWD->triangle_adjacent(i, 1) << " " << otmRWD->triangle_adjacent(i, 2) << std::endl;
    // Rcout << i << std::endl;
    // for (int j = 0; j <= 2; j++) {
    //   Rcout << X(otmRWD->cell_vertex(i, j), 0) << " "
    //         << X(otmRWD->cell_vertex(i, j), 1) << std::endl;
    // }
  }
}

// [[Rcpp::export]]
List GoF2D(const arma::mat &X, const arma::mat &Y, const arma::mat &XY, const arma::mat &U, 
           double epsilon, int maxit, bool verbose) {
  int d = 2;
  int n = X.n_rows;
  int m = Y.n_rows;
  int k = U.n_rows;

  // initialize the Geogram library.
  initializeGeogram();

  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const arma::mat cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.n_rows);
  setMeshPoint(unifMesh, cubeVertices);
  
  // must triangulate for OTM
  unifMesh.facets.create_triangle(0, 1, 2);
  unifMesh.facets.create_triangle(2, 1, 3);
  unifMesh.facets.connect();
  
  // embed in 3-dimension
  unifMesh.vertices.set_dimension(d + 1);

  // compute OTMs
  GEO::OptimalTransportMap2d OTMX(&unifMesh);
  GEO::OptimalTransportMap2d OTMY(&unifMesh);
  GEO::OptimalTransportMap2d OTMXY(&unifMesh);
  OTM2D(OTMX, X, epsilon, maxit, verbose);
  OTM2D(OTMY, Y, epsilon, maxit, verbose);
  OTM2D(OTMXY, XY, epsilon, maxit, verbose);
  
  // get weights
  arma::vec wX = getWeights(OTMX);
  arma::vec wY = getWeights(OTMY);

  // setup weighted vertices for X and Y
  double wVX[(d + 1) * n];
  double wVY[(d + 1) * m];
  getWeightedVerts(X, wX, wVX);
  getWeightedVerts(Y, wY, wVY);
  
  // build Kd-tree
  GEO::NearestNeighborSearch* treeX = GEO::NearestNeighborSearch::create(d + 1, "BNN");
  GEO::NearestNeighborSearch* treeY = GEO::NearestNeighborSearch::create(d + 1, "BNN");
  treeX->set_points(n, wVX);
  treeY->set_points(m, wVY);
  
  double p[3] = {};
  IntegerVector uX(k), uY(k);
  for (int i = 0; i < k; i++) {
    p[0] = U(i, 0);
    p[1] = U(i, 1);
    
    uX(i) = treeX->get_nearest_neighbor(p);
    uY(i) = treeY->get_nearest_neighbor(p);
  }
  
  // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.n_rows);
  setMeshPoint(unifMesh, cubeVertices);
  unifMesh.facets.create_quad(0, 2, 3, 1);
  
  // collect objects to return
  List lst;
  lst["U_Map_X"] = uX;
  lst["U_Map_Y"] = uY;
  lst["Vert_XY"] = getVertices(unifMesh, OTMXY);
  lst["Cent_XY"] = getCentroids(OTMXY);

  return lst;
}