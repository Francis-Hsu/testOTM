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

// setup a standard uniform mesh in 2D or 3D
void setUniMesh(GEO::Mesh &M, unsigned int d, bool tri) {
  const arma::mat cubeVertices = cubeVert(d);
  M.clear();
  M.vertices.create_vertices(cubeVertices.n_rows);
  setMeshPoint(M, cubeVertices);
  
  // embed in (d + 1)-dimension
  M.vertices.set_dimension(d + 1);
  
  if (d == 2) {
    if (tri) {
      M.facets.create_triangle(0, 1, 2);
      M.facets.create_triangle(2, 1, 3);
      M.facets.connect();
    } else {
      M.facets.create_quad(0, 2, 3, 1);
    }
  } else if (d == 3) {
    if (tri) {
      // divided the cube into 5 tetrahedra
      M.cells.create_tet(0, 2, 1, 4);
      M.cells.create_tet(3, 7, 1, 2);
      M.cells.create_tet(6, 4, 2, 7);
      M.cells.create_tet(5, 1, 7, 4);
      M.cells.create_tet(1, 2, 4, 7);
      M.cells.connect();
    } else {
      // took from initialize_mesh_with_box()
      M.facets.create_quad(7, 6, 2, 3);
      M.facets.create_quad(1, 3, 2, 0);
      M.facets.create_quad(5, 7, 3, 1);
      M.facets.create_quad(4, 6, 7, 5);
      M.facets.create_quad(4, 5, 1, 0);
      M.facets.create_quad(6, 4, 0, 2);
      M.facets.connect();
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

// extract vertices from a 2D mesh (embedded in 3D)
arma::mat getVertices2D(GEO::Mesh &M) {
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

// extract vertices from a 3D mesh
arma::mat getVertices3D(GEO::Mesh &M, bool volumetric) {
  // count total number of elements to traverse
  unsigned int nElems;
  if (volumetric) {
    // in volumetric the elements are the facets of cells
    nElems = 0;
    for (unsigned int i = 0; i < M.cells.nb(); i++) {
      nElems += M.cells.nb_facets(i);
    }
  } else {
    // in surfacic mode the elements are just facets
    nElems = M.facets.nb();
  }
  
  // count total number of (repeated) vertices to store
  unsigned int accuFacets = 0;
  int totalNbVert = 0;
  int nbVert[nElems];
  int accuVert[nElems];
  if (volumetric) {
    for (unsigned int i = 0; i < M.cells.nb(); i++) {
      for (unsigned int j = 0; j < M.cells.nb_facets(i); j++) {
        nbVert[accuFacets] = M.cells.facet_nb_vertices(i, j);
        accuVert[accuFacets] = totalNbVert;
        totalNbVert += nbVert[accuFacets];
        accuFacets++;
      }
    }
  } else {
    for (unsigned int i = 0; i < nElems; i++) {
      nbVert[i] = M.facets.nb_vertices(i);
      accuVert[i] = totalNbVert;
      totalNbVert += nbVert[i];
    }
  }
  
  
  GEO::vec3 v;
  arma::mat Vert(totalNbVert, 5 + volumetric);
  unsigned int id;
  if (volumetric) {
    // we iterate through facets of each cell and save the vertices in order
    accuFacets = 0;
    for (unsigned int i = 0; i < M.cells.nb(); i++) {
      for (unsigned int j = 0; j < M.cells.nb_facets(i); j++) {
        for (unsigned int k = 0; k < M.cells.facet_nb_vertices(i, j); k++) {
          id = M.cells.facet_vertex(i, j, k);
          v = M.vertices.point(id);
          Vert(accuVert[accuFacets] + k, 0) = i + 1;
          Vert(accuVert[accuFacets] + k, 1) = j + 1;
          Vert(accuVert[accuFacets] + k, 2) = v.x;
          Vert(accuVert[accuFacets] + k, 3) = v.y;
          Vert(accuVert[accuFacets] + k, 4) = v.z;
          Vert(accuVert[accuFacets] + k, 5) = id + 1;
        }
        accuFacets++;
      }
    }
  } else {
    for (unsigned int i = 0; i < nElems; i++) {
      for (int j = 0; j < nbVert[i]; j++) {
        id = M.facets.vertex(i, j);
        v = M.vertices.point(id);
        Vert(accuVert[i] + j, 0) = i + 1;
        Vert(accuVert[i] + j, 1) = v.x;
        Vert(accuVert[i] + j, 2) = v.y;
        Vert(accuVert[i] + j, 3) = v.z;
        Vert(accuVert[i] + j, 4) = id + 1;
      }
    }
  }

  return Vert;
}

// use the generic RVD class to extract vertices from a RVD
arma::mat getVerticesGen2D(GEO::Mesh &S, GEO::OptimalTransportMap &OTM) {
  const unsigned int n = OTM.nb_points();
  const unsigned int d = 2;
  
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
arma::mat cubeVert(unsigned int d) {
  if (d > 20) {
    stop("Input dimesion is too high!");
  }
  
  int x = 0;
  arma::mat verts(std::pow(2, d), d);
  for (unsigned int i = 0; i < verts.n_rows; i++) {
    x = i;
    for (unsigned int j = 0; j < d; j++) {
      verts(i, d - j - 1) = x % 2;
      x /= 2;
    }
  }
  
  return verts;
}