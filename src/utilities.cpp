#include "utilities.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

void initializeGeogram() {
  // use std::call_once?
  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "fast");
}

void setMeshPoint(GEO::Mesh &M, const arma::mat &X) {
  for (unsigned int i = 0; i < X.n_rows; i++) {
    double* p = M.vertices.point_ptr(i);
    for (unsigned int j = 0; j < X.n_cols; j++) {
      p[j] = X(i, j);
    }
  }
}

void getWeightedVerts(const arma::mat &X, const arma::vec &w, double* wV) {
  const int n = X.n_rows;
  const int d = X.n_cols;
  double wMax = w.max();
  
  for (int i = 0; i < (d + 1) * n; i++) {
    if (i % (d + 1) == d) {
      wV[i] = std::sqrt(wMax - w(i / (d + 1)));
    } else {
      wV[i] = X(i / (d + 1), i % (d + 1));
    }
  }
}

arma::vec getWeights(GEO::OptimalTransportMap &OTM) {
  arma::vec w(OTM.nb_points());
  for (unsigned int i = 0; i < OTM.nb_points(); i++) {
    w(i) = OTM.weight(i);
  }
  
  return w;
}

arma::mat getCentroids(GEO::OptimalTransportMap &OTM) {
  const unsigned int n = OTM.nb_points();
  const unsigned int d = OTM.dimension();
  
  // extract the centroids
  double ct[n * d];
  OTM.compute_Laguerre_centroids(ct);
  
  // save the centroids to a matrix
  arma::mat Centroid(n, d);
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < d; j++) {
      Centroid(i, j) = ct[i * d + j];
    }
  }
  
  return Centroid;
}

arma::mat getVertices(GEO::Mesh &M) {
  const unsigned int nFacets = M.facets.nb();
  GEO::vec3 v;
  
  
  int totalNbVert = 0;
  int nbVert[nFacets];
  int accuVert[nFacets];
  for (unsigned int i = 0; i < nFacets; i++) {
    nbVert[i] = M.facets.nb_vertices(i);
    accuVert[i] = totalNbVert;
    totalNbVert += nbVert[i];
  }
  
  unsigned int id;
  arma::mat Vert(totalNbVert, 4);
  for (unsigned int i = 0; i < nFacets; i++) {
    for (int j = 0; j < nbVert[i]; j++) {
      id = M.facets.vertex(i, j);
      v = M.vertices.point(id);
      Vert(accuVert[i] + j, 0) = i + 1;
      Vert(accuVert[i] + j, 1) = v.x;
      Vert(accuVert[i] + j, 2) = v.y;
      Vert(accuVert[i] + j, 3) = id + 1;
    }
  }
  
  return Vert;
}

arma::mat getVertices(GEO::Mesh &S, GEO::OptimalTransportMap &OTM) {
  const unsigned int n = OTM.nb_points();
  const unsigned int d = OTM.dimension();
  
  // construct a polygon from the source mesh
  GEO::Attribute<double> vertex_weight;
  vertex_weight.bind_if_is_defined(
    S.vertices.attributes(), "weight"
  );
  GEOGen::Polygon P;
  P.initialize_from_mesh_facet(&S, 0, false, vertex_weight);
  
  // some containers
  GEOGen::Polygon* cP; // polygon that stores the intersection
  std::vector<double> RVDVerts; 
  int totalNbVert = 0;
  int nbVert[n];
  int accuVert[n];
  
  // construct a generic RVD
  if (d == 2) {
    GEOGen::RVDHelper<3> transMapGen(OTM.RVD()->delaunay(), &S);
    
    // compute intersections
    for (unsigned int i = 0; i < n; i++) {
      cP = transMapGen.intersect_cell_facet(i, P);
      nbVert[i] = cP->nb_vertices();
      accuVert[i] = totalNbVert;
      totalNbVert += nbVert[i];
      for (unsigned int j = 0; j < cP->nb_vertices(); j++) {
        RVDVerts.push_back(cP->vertex(j).point()[0]); 
        RVDVerts.push_back(cP->vertex(j).point()[1]);
      }
    }
  } else {
    // take care 3d case in future
    // GEOGen::RestrictedVoronoiDiagram<4> transMapGen(OTM.RVD()->delaunay(), &S);
  }
  
  // get the vertices of each cell
  // should be in ccw orientation
  int currID;
  arma::mat Vert(totalNbVert, d + 1);
  for (unsigned int i = 0; i < n; i++) {
    for (int j = 0; j < nbVert[i]; j++) {
      currID = accuVert[i] + j;
      Vert(currID, 0) = i + 1;
      Vert(currID, 1) = RVDVerts[2 * currID];
      Vert(currID, 2) = RVDVerts[2 * currID + 1];
    }
  }
  
  return Vert;
}

// compute the vertices of a unit hypercube using binary expansion
// this will generate a matrix with 2^d rows
arma::mat cubeVert(int d) {
  if (d > 20) {
    stop("Dimesion is too high!");
  }
  
  int x = 0;
  arma::mat verts(std::pow(2, d), d);
  for (unsigned int i = 0; i < verts.n_rows; i++) {
    x = i;
    for (int j = 0; j < d; j++) {
      verts(i, d - j - 1) = x % 2;
      x /= 2;
    }
  }
  
  return verts;
}