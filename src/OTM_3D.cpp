#include "utilities.h"
using namespace Rcpp;

void OTM3D(GEO::OptimalTransportMap3d &OTM, const arma::mat &X, double epsilon, int maxit, bool verbose) {
  const int n = X.n_rows;
  const int d = 3;
  
  // create a mesh for data points
  GEO::Mesh dataMesh(d, false);
  dataMesh.vertices.create_vertices(n);
  setMeshPoint(dataMesh, X);
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
List dualGraphs3D(const arma::mat &X, double epsilon, int maxit, bool verbose) {
  const int d = 3;
  
  // initialize the Geogram library.
  initializeGeogram();
  
  // create a mesh for 3D uniform measure
  GEO::Mesh unifMesh(d, false);
  setUniMesh(unifMesh, d, true);
  
  // compute OTM
  GEO::OptimalTransportMap3d OTM(&unifMesh);
  OTM3D(OTM, X, epsilon, maxit, verbose);
  
  // must save centroids before remesh
  arma::mat Centroid = getCentroids(OTM);
  
  // store weights
  arma::vec W = getWeights(OTM);
  
  //GEO::Mesh otmRVD;
  //OTM.get_RVD(otmRVD);
  
  // get the RDT
  GEO::Mesh otmRDT;
  // OTM.RVD()->set_volumetric(true);
  OTM.RVD()->compute_RDT(otmRDT);
  
  // create a squared uniform mesh
  // setUniMesh(unifMesh, d, false);
  
  // collect objects to return
  List lst;
  lst["Data"] = X;
  lst["Centroid"] = Centroid;
  lst["Weight"] = W;
  lst["Vertex.RDT"] = getVertices3D(otmRDT);
  // lst["Vertex.RVD"] = getVerticesGen3D(unifMesh, OTM);
  // lst["Vertex.RVD"] = getVerticesGen3D(unifMesh, OTM);
  lst["N.Triangles"] = otmRDT.cells.nb();
  // lst["N.Cells"] = OTM.nb_points();
  
  return lst;
}