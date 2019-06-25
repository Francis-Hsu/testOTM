#ifndef UTILITIES_H
#define UTILITIES_H

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

#include <RcppArmadillo.h>

void initializeGeogram();
void setMeshPoint(GEO::Mesh &M, const arma::mat &X);
void getWeightedVerts(const arma::mat &X, const arma::vec &w, double* wV);
arma::vec getWeights(GEO::OptimalTransportMap &OTM);
arma::mat getCentroids(GEO::OptimalTransportMap &OTM);
arma::mat getVertices(GEO::Mesh &M);
arma::mat getVertices(GEO::Mesh &S, GEO::OptimalTransportMap &OTM);
arma::mat cubeVert(int d);

#endif