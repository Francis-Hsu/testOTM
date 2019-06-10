#include "utilities.h"
using namespace Rcpp;

void OTM2D(GEO::OptimalTransportMap2d &OTM, const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int n = X.nrow();
  int d = 2;

  // create a mesh for data points
  GEO::Mesh dataMesh(d, false);
  dataMesh.vertices.create_vertices(n);
  for (int i = 0; i < n; i++) {
    setMeshPoint(dataMesh, X, i);
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
List dualGraphes2D(const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int n = X.nrow();
  int d = 2;
  
  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "fast");
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    setMeshPoint(unifMesh, cubeVertices, i);
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
  double ct[n * d];
  OTM.compute_Laguerre_centroids(ct);
  NumericMatrix Centroid(n, d);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      Centroid(i, j) = ct[i * d + j];
    }
  }

  // store weights
  double w[n];
  double wMax;
  getWeights(OTM, wMax, w);

  // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    setMeshPoint(unifMesh, cubeVertices, i);
  }
  unifMesh.facets.create_quad(0, 2, 3, 1);

  // setup weighted vertices
  double wV[(d + 1) * n];
  getWeightedVerts(X, wMax, w, n, d, wV);
  
  // generate a new regular triangulation based on weights computed
  GEO::Delaunay_var otmRWD = GEO::Delaunay::create(d + 1);
  otmRWD->set_vertices(n, wV);
  
  // create a new RVD and RDT
  GEO::RestrictedVoronoiDiagram_var RVD = GEO::RestrictedVoronoiDiagram::create(otmRWD, &unifMesh);
  GEO::Mesh otmRDT, otmRVD;
  RVD->compute_RDT(otmRDT);
  RVD->compute_RVD(otmRVD);
  
  // get the weights
  NumericVector Weight(w, w + sizeof(w) / sizeof(*w));
  
  // extract all vertices of RVD and RDT
  NumericMatrix rdtVert = getVertices(otmRDT);
  NumericMatrix rvdVert = getVertices(otmRVD);
  
  List ret;
  ret["Data"] = X;
  ret["Centroid"] = Centroid;
  ret["Weight"] = Weight;
  ret["Vertex.RDT"] = rdtVert;
  ret["Vertex.RVD"] = rvdVert;
  ret["N.Triangles"] = otmRDT.facets.nb();
  ret["N.Cells"] = otmRVD.facets.nb();
  
  return ret;
}

// [[Rcpp::export]]
List GoF2D(const NumericMatrix &X, const NumericMatrix &Y, const NumericMatrix &XY, const NumericMatrix &U) {
  int d = 2;
  int n = X.nrow();
  int m = Y.nrow();
  int k = U.nrow();

  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "exact");

  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    setMeshPoint(unifMesh, cubeVertices, i);
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
  double wX[n];
  double wY[m];
  double wXY[n + m];
  double xwMax, ywMax, xywMax;
  getWeights(OTMX, xwMax, wX);
  getWeights(OTMY, ywMax, wY);
  getWeights(OTMXY, xywMax, wXY);

  // setup weighted vertices for X and Y
  double wVX[(d + 1) * n];
  double wVY[(d + 1) * m];
  double wVXY[(d + 1) * (n + m)];
  getWeightedVerts(X, xwMax, wX, n, d, wVX);
  getWeightedVerts(Y, ywMax, wY, m, d, wVY);
  getWeightedVerts(XY, xywMax, wXY, n + m, d, wVXY);

  // triangulate
  GEO::Delaunay_var transTriX = GEO::Delaunay::create(d + 1);
  GEO::Delaunay_var transTriY = GEO::Delaunay::create(d + 1);
  GEO::Delaunay_var transTriXY = GEO::Delaunay::create(d + 1);
  transTriX->set_vertices(n, wVX);
  transTriY->set_vertices(m, wVY);
  transTriXY->set_vertices(n + m, wVXY);

  double p[2];
  IntegerVector uX(k), uY(k);
  for (int i = 0; i < k; i++) {
    p[0] = U(i, 0);
    p[1] = U(i, 1);
    uX(i) = transTriX->nearest_vertex(p);
    uY(i) = transTriY->nearest_vertex(p);
  }
  
  // create a squared uniform mesh
  unifMesh.clear();
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  for (int i = 0; i < cubeVertices.nrow(); i++) {
    setMeshPoint(unifMesh, cubeVertices, i);
  }
  unifMesh.facets.create_quad(0, 2, 3, 1);

  // create new RVDs
  GEO::RestrictedVoronoiDiagram_var RVDXY = GEO::RestrictedVoronoiDiagram::create(transTriXY, &unifMesh);
  GEO::Mesh transMapXY;
  RVDXY->compute_RVD(transMapXY);

  // extract all vertices
  NumericMatrix Vert = getVertices(transMapXY);
  
  List rtn;
  rtn["U_Map_X"] = uX;
  rtn["U_Map_Y"] = uY;
  rtn["Vert_XY"] = Vert;

  return rtn;
}