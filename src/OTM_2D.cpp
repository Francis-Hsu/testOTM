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
arma::vec dualPotential2D(const arma::mat &Y, const arma::mat &X, const arma::mat &V, const arma::vec &h, 
                          const arma::uvec accuVerts) {
  int m = Y.n_rows;
  int n = X.n_rows;
  
  arma::mat currCellVerts;
  arma::vec cellMax(n);
  arma::vec psi(m, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      currCellVerts = V.rows(accuVerts(j), accuVerts(j + 1) - 1);
      currCellVerts.each_row() -= X.row(j);
      cellMax(j) = arma::max(currCellVerts * Y.row(i).t());
    }
    
    psi(i) = arma::max(cellMax - h);
  }
  
  return psi;
}

// find the cells that a set Q of query points was transported to
// [[Rcpp::export]]
arma::ivec locateRVD2D(const arma::mat &Q, const arma::mat &X, const arma::vec &w) {
  int m = Q.n_rows;
  int n = X.n_rows;
  int d = 2;
  
  // setup weighted vertices for X
  double wX[(d + 1) * n];
  getWeightedVerts(X, w, wX);
  
  // build a Kd-tree
  GEO::NearestNeighborSearch* treeX = GEO::NearestNeighborSearch::create(d + 1, "BNN");
  treeX->set_points(n, wX);
  
  double p[d + 1] = {};
  arma::ivec cellID(m, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    p[0] = Q(i, 0);
    p[1] = Q(i, 1);
    
    cellID(i) = treeX->get_nearest_neighbor(p) + 1;
  }
  
  return cellID;
}

// find the triangles that contain a set Q of query points
// [[Rcpp::export]]
arma::ivec locateRDT2D(const arma::mat &Q, const arma::mat &V) {
  int m = Q.n_rows; // number of query points
  int n = V.n_rows / 3; // number of triangles
  
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


// [[Rcpp::export]]
List GoF2D(const arma::mat &X, const arma::mat &Y, const arma::mat &XY, const arma::mat &U, 
           double epsilon, int maxit, bool verbose) {
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
  
  // get the optimal transport mapping
  arma::ivec uX = locateRVD2D(U, X, wX);
  arma::ivec uY = locateRVD2D(U, Y, wY);
  
  // must save centroids before remesh
  arma::mat centroidXY = getCentroids(OTMXY);
  
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
  lst["Cent_XY"] = centroidXY;
  
  return lst;
}