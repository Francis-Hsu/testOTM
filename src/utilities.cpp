#include "utilities.h"
using namespace Rcpp;

void initializeGeogram() {
  // use std::call_once?
  
  // initialize the Geogram library.
  GEO::initialize();
  
  // import command line arguments.
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("algo:predicates", "fast");
}

void setMeshPoint(GEO::Mesh &M, const NumericMatrix &X) {
  for (int i = 0; i < X.nrow(); i++) {
    double* p = M.vertices.point_ptr(i);
    for (int j = 0; j < X.ncol(); j++) {
      p[j] = X(i, j);
    }
  }
}

void getWeights(GEO::OptimalTransportMap &OTM, double &wMax, double* w) {
  wMax = OTM.weight(0);
  for (unsigned int i = 0; i < OTM.nb_points(); i++) {
    w[i] = OTM.weight(i);
    if (wMax < w[i]) {
      wMax = w[i];
    }
  }
}

void getWeightedVerts(const NumericMatrix &X, double &wMax, double* w, double* wV) {
  int n = X.nrow();
  int d = X.ncol();
  
  for (int i = 0; i < (d + 1) * n; i++) {
    if (i % (d + 1) == d) {
      wV[i] = std::sqrt(wMax - w[i / (d + 1)]);
    } else {
      wV[i] = X(i / (d + 1), i % (d + 1));
    }
  }
}

NumericMatrix getVertices(GEO::Mesh &M) {
  GEO::vec3 v;
  unsigned int nFacets = M.facets.nb();

  int totalNbVert = 0;
  int nbVert[nFacets];
  int accuVert[nFacets];
  for (unsigned int i = 0; i < nFacets; i++) {
    nbVert[i] = M.facets.nb_vertices(i);
    accuVert[i] = totalNbVert;
    totalNbVert += nbVert[i];
  }

  NumericMatrix Vert(totalNbVert, 3);
  for (unsigned int i = 0; i < nFacets; i++) {
    for (int j = 0; j < nbVert[i]; j++) {
      v = M.vertices.point(M.facets.vertex(i, j));
      Vert(accuVert[i] + j, 0) = i + 1;
      Vert(accuVert[i] + j, 1) = v.x;
      Vert(accuVert[i] + j, 2) = v.y;
    }
  }
  
  return Vert;
}

Rcpp::NumericMatrix getVertices(GEO::Mesh &S, GEO::OptimalTransportMap &OTM) {
  unsigned int n = OTM.nb_points();
  unsigned int d = OTM.dimension();
  
  // construct polygon from source mesh

  GEO::Attribute<double> vertex_weight;
  vertex_weight.bind_if_is_defined(
    S.vertices.attributes(), "weight"
  );
  GEOGen::Polygon P;
  P.initialize_from_mesh_facet(&S, 0, false, vertex_weight);
  
  // some containers
  GEOGen::Polygon* cP; // polygon storing intersection
  std::vector<double> RVDVerts; 
  int totalNbVert = 0;
  int nbVert[n];
  int accuVert[n];
  
  // construct generic RVD
  if (d == 2) {
    GEOGen::RestrictedVoronoiDiagram<3> transMapGen(OTM.RVD()->delaunay(), &S);
    // compute intersections
    for (int i = 0; i < n; i++) {
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
    GEOGen::RestrictedVoronoiDiagram<4> transMapGen(OTM.RVD()->delaunay(), &S);
  }
  
  int currID;
  NumericMatrix Vert(totalNbVert, d + 1);
  for (int i = 0; i < n; i++) {
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
NumericMatrix cubeVert(int d) {
  if (d > 20) {
    stop("Dimesion is too high!");
  }
  
  int x = 0;
  NumericMatrix verts(std::pow(2, d), d);
  for (int i = 0; i < verts.nrow(); i++) {
    x = i;
    for (int j = 0; j < d; j++) {
      verts(i, d - j - 1) = x % 2;
      x /= 2;
    }
  }
  
  return verts;
}