#include "utilities.h"
using namespace Rcpp;

void OTM2D(GEO::OptimalTransportMap2d &OTM, const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int n = X.nrow();
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
List dualGraphs2D(const NumericMatrix &X, double epsilon, int maxit, bool verbose) {
  int n = X.nrow();
  int d = 2;
  
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
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
  setMeshPoint(unifMesh, cubeVertices);
  unifMesh.facets.create_quad(0, 2, 3, 1);
  
  // get RVD and RDT
  GEO::Mesh otmRDT, otmRVD;
  OTM.get_RVD(otmRVD); // 2D ordered incorrectly? See get_RVD()
  OTM.RVD()->set_volumetric(false);
  OTM.RVD()->compute_RDT(otmRDT);
  
  // extract vertices from RVD and RDT
  NumericMatrix rdtVert = getVertices(otmRDT);
  NumericMatrix rvdVert = getVertices(otmRVD);

  List ret;
  ret["Data"] = X;
  ret["Centroid"] = Centroid;
  ret["Weight"] = NumericVector(w, w + sizeof(w) / sizeof(*w));
  ret["Vertex.RDT"] = rdtVert;
  ret["Vertex.RVD"] = rvdVert;
  ret["N.Triangles"] = otmRDT.facets.nb();
  ret["N.Cells"] = otmRVD.facets.nb();
  
  return ret;
}

// under construction
// [[Rcpp::export]]
void OTMRank2D(const NumericMatrix &X, NumericVector &weight, double wMax, const NumericMatrix &Q) {
  int n = X.nrow();
  int d = 2;
  
  // initialize the Geogram library.
  initializeGeogram();
  
  double* w = weight.begin();
  double wV[(d + 1) * n];
  getWeightedVerts(X, wMax, w, wV);
  
  GEO::RegularWeightedDelaunay2d* otmRWD = dynamic_cast<GEO::RegularWeightedDelaunay2d*>(GEO::Delaunay::create(d + 1, "BPOW2d"));
  //otmRWD->set_keeps_infinite(true);
  otmRWD->set_vertices(n, wV);
  
  GEO::index_t cellID;
  Rcout << "#########" << std::endl;
  double q[2];
  for (int i = 0; i < Q.nrow(); i++) {
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
List GoF2D(const NumericMatrix &X, const NumericMatrix &Y, const NumericMatrix &XY, const NumericMatrix &U, 
           double epsilon, int maxit, bool verbose) {
  int d = 2;
  int n = X.nrow();
  int m = Y.nrow();
  int k = U.nrow();

  // initialize the Geogram library.
  initializeGeogram();

  // create a mesh for 2D uniform measure
  GEO::Mesh unifMesh(d, false);
  const NumericMatrix cubeVertices = cubeVert(2);
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
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
  double wX[n];
  double wY[m];
  double xwMax, ywMax;
  getWeights(OTMX, xwMax, wX);
  getWeights(OTMY, ywMax, wY);

  // setup weighted vertices for X and Y
  double wVX[(d + 1) * n];
  double wVY[(d + 1) * m];
  getWeightedVerts(X, xwMax, wX, wVX);
  getWeightedVerts(Y, ywMax, wY, wVY);
  
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
  unifMesh.vertices.create_vertices(cubeVertices.nrow());
  setMeshPoint(unifMesh, cubeVertices);
  unifMesh.facets.create_quad(0, 2, 3, 1);
  
  // extract all vertices
  std::vector<double> RVDVerts;
  GEO::Attribute<double> vertex_weight;
  vertex_weight.bind_if_is_defined(
    unifMesh.vertices.attributes(), "weight"
  );
  GEOGen::Polygon P;
  P.initialize_from_mesh_facet(&unifMesh, 0, false, vertex_weight);
  GEOGen::Polygon* cellIntersect;
  GEOGen::RestrictedVoronoiDiagram<3> transMapXY(OTMXY.RVD()->delaunay(), &unifMesh);
  
  // extract vertices
  int totalNbVert = 0;
  int nbVert[n + m];
  int accuVert[n + m];
  for (int i = 0; i < n + m; i++) {
    cellIntersect = transMapXY.intersect_cell_facet(i, P);
    nbVert[i] = cellIntersect->nb_vertices();
    accuVert[i] = totalNbVert;
    totalNbVert += nbVert[i];
    for (unsigned int j = 0; j < cellIntersect->nb_vertices(); j++) {
      RVDVerts.push_back(cellIntersect->vertex(j).point()[0]);
      RVDVerts.push_back(cellIntersect->vertex(j).point()[1]);
    }
  }

  int currID;
  NumericMatrix Vert(totalNbVert, d + 1);
  for (int i = 0; i < n + m; i++) {
    for (int j = 0; j < nbVert[i]; j++) {
      currID = accuVert[i] + j;
      Vert(currID, 0) = i + 1;
      Vert(currID, 1) = RVDVerts[2 * currID];
      Vert(currID, 2) = RVDVerts[2 * currID + 1];
    }
  }
  
  List rtn;
  rtn["U_Map_X"] = uX;
  rtn["U_Map_Y"] = uY;
  rtn["Vert_XY"] = Vert;

  return rtn;
}