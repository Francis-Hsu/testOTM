#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
 
#include <exploragram/optimal_transport/optimal_transport.h>
#include <exploragram/optimal_transport/optimal_transport_2d.h>
#include <exploragram/optimal_transport/optimal_transport_3d.h>
#include <exploragram/optimal_transport/sampling.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/generic_RVD_cell.h>
 
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/delaunay/delaunay_3d.h>

#include <geogram/points/nn_search.h>
#include <geogram/numerics/optimizer.h>
#include <geogram/numerics/lbfgs_optimizers.h>
 
#include <geogram/basic/command_line.h>
#include <geogram/basic/memory.h>

#ifdef Realloc
#undef Realloc
#endif
 
#ifdef Free
#undef Free
#endif
 
#ifdef ERROR
#undef ERROR
#endif

#include <Rcpp.h>
using namespace Rcpp;

namespace GEO { 
  inline void set_mesh_point(Mesh &M, const NumericMatrix &X, int row) {
   geo_debug_assert(M.vertices.dimension() >= data.ncol());
   double* p = M.vertices.point_ptr(row);
   for (int col = 0; col < X.ncol(); col++) {
     p[col] = X(row, col);
   }
  }
}

// write a wrapper struct?

// R-like modulo operator
inline int divisorModulo(int a, int b) {
  return ((a % b) + b) % b;
}

// compute the vertices of a unit hypercube using binary expansion
// this will generate a matrix with 2^d rows
NumericMatrix CubeVert(int d) {
  if (d > 20) {
    stop("Dimesion is too high!");
  }
  
  int x = 0;
  NumericMatrix rtn(std::pow(2, d), d);
  for (int  i = 0; i < std::pow(2, d); i++) {
    x = i;
    for (int j = 0; j < d; j++) {
      rtn(i, d - j - 1) = x % 2;
      x /= 2;
    }
  }
  
  return rtn;
}

void extractWeights(GEO::OptimalTransportMap &OTM, double &wMax, double* weights) {
  wMax = OTM.weight(0);
  for (unsigned int i = 0; i < OTM.nb_points(); i++) {
    weights[i] = OTM.weight(i);
    if (wMax < weights[i]) {
      wMax = weights[i];
    }
  }
}

NumericMatrix extractVertices(GEO::Mesh &M) {
  GEO::vec3 v;
  unsigned int nFacets = M.facets.nb();
  unsigned int vertID;

  int totalNbVert = 0;
  int* nbVert = new int[nFacets];
  int* accuVert = new int[nFacets];
  for (unsigned int i = 0; i < nFacets; i++) {
    nbVert[i] = M.facets.nb_vertices(i);
    accuVert[i] = totalNbVert;
    totalNbVert += nbVert[i];
  }

  NumericMatrix Vert(totalNbVert, 3);
  for (unsigned int i = 0; i < nFacets; i++) {
    for (int j = 0; j < nbVert[i]; j++) {
      vertID = M.facets.vertex(i, j);
      v = M.vertices.point(vertID);
      Vert(accuVert[i] + j, 0) = v.x;
      Vert(accuVert[i] + j, 1) = v.y;
      Vert(accuVert[i] + j, 2) = i;
    }
  }

  return Vert;
}

NumericMatrix extractEdges(GEO::Mesh &M) {
  GEO::vec3 v;
  unsigned int nVerts = M.vertices.nb();
  unsigned int nFacets = M.facets.nb();

  // extract vertices in order
  NumericMatrix Vert(nVerts, 2);
  for (unsigned int i = 0; i < nVerts; i++) {
    v = M.vertices.point(i);
    Vert(i, 0) = v.x;
    Vert(i, 1) = v.y;
  }

  // find out number of vertices in each cell
  int totalNbVert = 0;
  int* nbVert = new int[nFacets];
  for (unsigned int i = 0; i < nFacets; i++) {
    nbVert[i] = M.facets.nb_vertices(i);
    totalNbVert += nbVert[i];
  }

  // store edges
  int currVertID, currNbVert, currAccuVert = 0;
  NumericMatrix Edge(totalNbVert, 4);
  for (unsigned int i = 0; i < nFacets; i++) {
    currNbVert = nbVert[i];
    for (int j = 0; j < currNbVert; j++) {
      currVertID = M.facets.vertex(i, j);
      Edge(currAccuVert + j, 0) = Vert(currVertID, 0);
      Edge(currAccuVert + j, 1) = Vert(currVertID, 1);
      Edge(currAccuVert + divisorModulo(j - 1, currNbVert), 2) = Vert(currVertID, 0);
      Edge(currAccuVert + divisorModulo(j - 1, currNbVert), 3) = Vert(currVertID, 1);
    }
    currAccuVert += currNbVert;
  }

  return Edge;
}

void OTM2D(GEO::OptimalTransportMap2d &OTM, const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int d = 2;
  int n = X.nrow();

  // create a mesh for data points
  GEO::Mesh dataMesh(d, false);
  dataMesh.vertices.create_vertices(n);
  for (int i = 0; i < n; i++) {
    GEO::set_mesh_point(dataMesh, X, i);
  }

  // embed in 3-dimension
  dataMesh.vertices.set_dimension(d + 1);

  // setup OTM
  OTM.set_points(n, dataMesh.vertices.point_ptr(0), dataMesh.vertices.dimension());

  // optmizer settings
  OTM.set_verbose(verbose);
  OTM.set_epsilon(epsilon);
  OTM.set_regularization(0.0);
  OTM.set_Newton(true); // BFGS not implemented?
  OTM.optimize(maxit);
}

// [[Rcpp::export]]
List powerDiag2D(const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int n = X.nrow();
  int d = 2;
  
  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "exact");
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = CubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    GEO::set_mesh_point(unifMesh, cubeVertices, i);
  }
  
  // must triangulate for OTM
  unifMesh.facets.create_triangle(0, 1, 2);
  unifMesh.facets.create_triangle(2, 1, 3);
  unifMesh.facets.connect();
  
  // embed in 3-dimension
  unifMesh.vertices.set_dimension(d + 1);
  
  // compute OTM
  GEO::OptimalTransportMap2d OTM(&unifMesh);
  OTM2D(OTM, X, epsilon, maxit, verbose);
  
  // store centroids
  double* Centroids = new double[n * d];
  OTM.compute_Laguerre_centroids(Centroids);
  NumericMatrix centroids(n, d);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      centroids(i, j) = Centroids[i * d + j];
    }
  }

  // store weights
  double* Weights = new double[n];
  double wMax = OTM.weight(0);
  extractWeights(OTM, wMax, Weights);

  // store the normalized weights
  NumericVector weight(n, wMax);
  for (int i = 0; i < n; i++) {
    weight(i) -= Weights[i];
  }

  // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    GEO::set_mesh_point(unifMesh, cubeVertices, i);
  }
  unifMesh.facets.create_quad(0, 2, 3, 1);

  // setup weighted vertices
  double* weightedVerts = new double[(d + 1) * n];
  for (int i = 0; i < (d + 1) * n; i++) {
    if (i % (d + 1) == d) {
      weightedVerts[i] = (double) std::sqrt(wMax - Weights[i / (d + 1)]);
    } else {
      weightedVerts[i] = X(i / (d + 1), i % (d + 1));
    }
  }

  // generate a new regular triangulation based on weights computed
  GEO::Delaunay_var transTri = GEO::Delaunay2d::create(d + 1);
  transTri->set_vertices(n, weightedVerts);

  // create a new RVD
  GEO::RestrictedVoronoiDiagram_var RVD = GEO::RestrictedVoronoiDiagram::create(transTri, &unifMesh);
  GEO::Mesh transMap;
  RVD->compute_RVD(transMap);

  // extract all vertices
  NumericMatrix Vert = extractVertices(transMap);

  // store edges
  NumericMatrix Edge = extractEdges(transMap);
  
  List ret;
  ret["Data"] = X;
  ret["MaxWeight"] = wMax;
  ret["Weight"] = weight;
  ret["Vert"] = Vert;
  ret["Edge"] = Edge;
  ret["Centroids"] = centroids;

  return ret;
}

// [[Rcpp::export]]
List GoF2D(const NumericMatrix &X, const NumericMatrix &Y, const NumericMatrix &XY, const NumericMatrix &U) {
  int d = 2;
  int n = X.nrow();
  int m = Y.nrow();
  int k = U.nrow();

  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = CubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    GEO::set_mesh_point(unifMesh, cubeVertices, i);
  }

  // must triangulate for OTM
  unifMesh.facets.create_triangle(0, 1, 2);
  unifMesh.facets.create_triangle(2, 1, 3);
  unifMesh.facets.connect();
  unifMesh.vertices.set_dimension(d + 1);

  // compute OTMs
  GEO::OptimalTransportMap2d OTMX(&unifMesh);
  GEO::OptimalTransportMap2d OTMY(&unifMesh);
  GEO::OptimalTransportMap2d OTMXY(&unifMesh);
  OTM2D(OTMX, X, 1e-3, 250, false);
  OTM2D(OTMY, Y, 1e-3, 250, false);
  OTM2D(OTMXY, XY, 1e-3, 250, false);
  
  // get weights
  double* wX = new double[n];
  double* wY = new double[m];
  double* wXY = new double[n + m];
  double xwMax = OTMX.weight(0);
  double ywMax = OTMY.weight(0);
  double xywMax = OTMXY.weight(0);
  extractWeights(OTMX, xwMax, wX);
  extractWeights(OTMY, ywMax, wY);
  extractWeights(OTMXY, xywMax, wXY);

  // // setup weighted vertices for X and Y
  double* wVX = new double[(d + 1) * n];
  double* wVY = new double[(d + 1) * m];
  double* wVXY = new double[(d + 1) * (n + m)];
  //
  for (int i = 0; i < (d + 1) * n; i++) {
    if (i % (d + 1) == d) {
      wVX[i] = (double) std::sqrt(xwMax - wX[i / (d + 1)]);
    } else {
      wVX[i] = X(i / (d + 1), i % (d + 1));
    }
  }
  //
  for (int i = 0; i < (d + 1) * m; i++) {
    if (i % (d + 1) == d) {
      wVY[i] = (double) std::sqrt(ywMax - wY[i / (d + 1)]);
    } else {
      wVY[i] = Y(i / (d + 1), i % (d + 1));
    }
  }

  for (int i = 0; i < (d + 1) * (n + m); i++) {
    if (i % (d + 1) == d) {
      wVXY[i] = (double) std::sqrt(xywMax - wXY[i / (d + 1)]);
    } else {
      wVXY[i] = XY(i / (d + 1), i % (d + 1));
    }
  }

  // triangulate
  GEO::Delaunay_var transTriX = GEO::Delaunay2d::create(d + 1);
  GEO::Delaunay_var transTriY = GEO::Delaunay2d::create(d + 1);
  GEO::Delaunay_var transTriXY = GEO::Delaunay2d::create(d + 1);
  transTriX->set_vertices(n, wVX);
  transTriY->set_vertices(m, wVY);
  transTriXY->set_vertices(n + m, wVXY);
  //
  double* p = new double[2];
  IntegerVector uX(k), uY(k);
  for (int i = 0; i < k; i++) {
    p[0] = U(i, 0);
    p[1] = U(i, 1);
    uX(i) = transTriX->nearest_vertex(p);
    uY(i) = transTriY->nearest_vertex(p);
  }
  //
  // // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    GEO::set_mesh_point(unifMesh, cubeVertices, i);
  }
  unifMesh.facets.create_quad(0, 2, 3, 1);

  // create new RVDs
  GEO::RestrictedVoronoiDiagram_var RVDXY = GEO::RestrictedVoronoiDiagram::create(transTriXY, &unifMesh);
  GEO::Mesh transMapXY;
  RVDXY->compute_RVD(transMapXY);

  // extract all vertices
  NumericMatrix Vert = extractVertices(transMapXY);

  List rtn;
  rtn["U_Map_X"] = uX;
  rtn["U_Map_Y"] = uY;
  rtn["Vert_XY"] = Vert;

  return rtn;
}


// 3d CCW orientation
// unifMesh.cells.create_tet(1, 6, 7, 5);
// unifMesh.cells.create_tet(1, 4, 7, 6);
// unifMesh.cells.create_tet(1, 2, 4, 6);
// unifMesh.cells.create_tet(6, 4, 8, 7);
// unifMesh.cells.create_tet(1, 3, 4, 7);
// unifMesh.cells.create_hex(1, 2, 4, 3, 7, 8, 6, 5);