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

#include <Rcpp.h>

void setMeshPoint(GEO::Mesh &M, const Rcpp::NumericMatrix &X, int row);
void getWeights(GEO::OptimalTransportMap &OTM, double &wMax, double* w);
void getWeightedVerts(const Rcpp::NumericMatrix &X, double &wMax, double* w, int n, int d, double* wV);
Rcpp::NumericMatrix getVertices(GEO::Mesh &M);
Rcpp::NumericMatrix cubeVert(int d);

#endif