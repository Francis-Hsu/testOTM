#include <Rcpp.h>
using namespace Rcpp;

#ifdef Realloc
#undef Realloc
#endif

#ifdef Free
#undef Free
#endif
 
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
 
namespace GEO { 
  inline void set_mesh_point(Mesh &M, const NumericMatrix &X, int row) {
   geo_debug_assert(M.vertices.dimension() >= data.ncol());
   double* p = M.vertices.point_ptr(row);
   for (int col = 0; col < X.ncol(); col++) {
     p[col] = X(row, col);
   }
  }
}

// write a wrapper class/struct?

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

GEO::OptimalTransportMap2d OTM2D(const NumericMatrix &X) {
  int d = 2;
  int n = X.nrow();

  // create a mesh for data points
  GEO::Mesh dataMesh(d, false);
  dataMesh.vertices.create_vertices(n);
  for (int i = 0; i < n; i++) {
    GEO::set_mesh_point(dataMesh, X, i);
  }

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
  dataMesh.vertices.set_dimension(d + 1);

  // compute the 2D optimal transportation map
  GEO::OptimalTransportMap2d OTM(&unifMesh);
  OTM.set_points(n, dataMesh.vertices.point_ptr(0), dataMesh.vertices.dimension());

  // optmizer settings
  OTM.set_verbose(false);
  OTM.set_epsilon(1e-6);
  OTM.set_regularization(0.0);
  OTM.set_Newton(true); // BFGS not implemented?
  OTM.optimize(1000);

  return OTM;
}


// List OTMGoF2D(const NumericMatrix &X, const NumericMatrix &Y, const NumericMatrix &U) {
//   int n = X.nrow();
//   int m = Y.nrow();
//   
//   GEO::OptimalTransportMap2d OTMX = OTM2D(X);
//   GEO::OptimalTransportMap2d OTMY = OTM2D(Y);
//   
//   List rtn;
//   
//   return rtn;
// }

// [[Rcpp::export]]
List powDiag(const NumericMatrix &X) {
  int n = X.nrow();
  int d = 2;
  
  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "exact");
  
  GEO::OptimalTransportMap2d OTM = OTM2D(X);
  
  Rcout << OTM.nb_points() << std::endl;
  Rcout << OTM.dimension() << std::endl;
  // 
  double* Centroids = new double[n * d];
  // double* Weights = new double[n];
  // 
  // save centroids
  OTM.compute_Laguerre_centroids(Centroids);
  NumericMatrix centroids(n, d);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      centroids(i, j) = Centroids[i * d + j];
    }
  }
  // 
  // // save weights
  // double wMax = OTM.weight(0);
  // for (int i = 0; i < n; i++) {
  //   Weights[i] = OTM.weight(i);
  //   if (wMax < Weights[i]) {
  //     wMax = Weights[i];
  //   }
  // }
  // 
  // // save the normalized weights
  // NumericVector weight(n, wMax);
  // for (int i = 0; i < n; i++) {
  //   weight(i) -= Weights[i];
  // }
  // 
  // // create a squared uniform mesh
  // GEO::Mesh unifMesh(d, false);
  // const NumericMatrix cubeVertices = CubeVert(d);
  // unifMesh.vertices.create_vertices(cubeVertices.nrow());
  // for (int i = 0; i < cubeVertices.nrow(); i++) {
  //   GEO::set_mesh_point(unifMesh, cubeVertices, i);
  // }
  // unifMesh.facets.create_quad(0, 2, 3, 1);
  // 
  // // setup weighted vertices
  // double* weightedVerts = new double[(d + 1) * n];
  // for (int i = 0; i < (d + 1) * n; i++) {
  //   if (i % (d + 1) == d) {
  //     weightedVerts[i] = (double) std::sqrt(wMax - Weights[i / (d + 1)]);
  //   } else {
  //     weightedVerts[i] = X(i / (d + 1), i % (d + 1));
  //   }
  // }
  // 
  // // generate a new regular triangulation based on weights computed
  // GEO::Delaunay_var transTri = GEO::Delaunay2d::create(d + 1);
  // transTri->set_vertices(n, weightedVerts);
  // 
  // // create a new RVD
  // GEO::RestrictedVoronoiDiagram_var RVD = GEO::RestrictedVoronoiDiagram::create(transTri, &unifMesh);
  // GEO::Mesh transMap;
  // RVD->compute_RVD(transMap);
  // 
  // // extract all vertices
  // GEO::vec3 v;
  // unsigned int nVert = transMap.vertices.nb();
  // NumericMatrix Vert(nVert, 2);
  // for (unsigned int i = 0; i < nVert; i++) {
  //   v = transMap.vertices.point(i);
  //   Vert(i, 0) = v.x;
  //   Vert(i, 1) = v.y;
  // }
  // 
  // // store edges
  // int nFacets = transMap.facets.nb();
  // int totalNbVert = 0;
  // int* nbVert = new int[nFacets];
  // for (int i = 0; i < nFacets; i++) {
  //   nbVert[i] = transMap.facets.nb_vertices(i);
  //   totalNbVert += nbVert[i];
  // }
  // 
  // int currVertID, currNbVert, currAccuVert = 0;
  // NumericMatrix Edge(totalNbVert, 4);
  // for (int i = 0; i < nFacets; i++) {
  //   currNbVert = nbVert[i];
  //   for (int j = 0; j < currNbVert; j++) {
  //     currVertID = transMap.facets.vertex(i, j);
  //     Edge(currAccuVert + j, 0) = Vert(currVertID, 0);
  //     Edge(currAccuVert + j, 1) = Vert(currVertID, 1);
  //     Edge(currAccuVert + divisorModulo(j - 1, currNbVert), 2) = Vert(currVertID, 0);
  //     Edge(currAccuVert + divisorModulo(j - 1, currNbVert), 3) = Vert(currVertID, 1);
  //   }
  //   currAccuVert += currNbVert;
  // }
  
  List ret;
  // ret["MaxWeight"] = wMax;
  // ret["Weight"] = weight;
  // ret["Edge"] = Edge;
  ret["Centroids"] = centroids;

  return ret;
}


// 3d CCW orientation
// unifMesh.cells.create_tet(1, 6, 7, 5);
// unifMesh.cells.create_tet(1, 4, 7, 6);
// unifMesh.cells.create_tet(1, 2, 4, 6);
// unifMesh.cells.create_tet(6, 4, 8, 7);
// unifMesh.cells.create_tet(1, 3, 4, 7);
// unifMesh.cells.create_hex(1, 2, 4, 3, 7, 8, 6, 5);