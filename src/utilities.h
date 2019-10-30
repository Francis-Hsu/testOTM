#ifndef UTILITIES_H
#define UTILITIES_H

#include <RcppArmadillo.h>

#include "geogram.h"

#ifdef Realloc
#undef Realloc
#endif

#ifdef Free
#undef Free
#endif

#ifdef ERROR
#undef ERROR
#endif

void initializeGeogram();
void setMeshPoint(GEO::Mesh &M, const arma::mat &X);
void getWeightedVerts(const arma::mat &X, const arma::vec &w, double* wV);
void setSqUniMesh(GEO::Mesh &M, unsigned int d, bool tri);
arma::vec getWeights(GEO::OptimalTransportMap &OTM);
arma::mat getCentroids(GEO::OptimalTransportMap &OTM);
arma::mat getVertices(GEO::Mesh &M);
arma::mat getVertices(GEO::Mesh &S, GEO::OptimalTransportMap &OTM);
arma::mat cubeVert(unsigned int d);

namespace GEOGen {
  template <GEO::index_t DIM> class RVDHelper: public RestrictedVoronoiDiagram<DIM> {
  public:
    using RestrictedVoronoiDiagram<DIM>::RestrictedVoronoiDiagram;
    using RestrictedVoronoiDiagram<DIM>::intersect_cell_facet;
  };
}

#endif