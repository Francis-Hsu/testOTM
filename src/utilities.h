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

void initializeGeogram();
void setMeshPoint(GEO::Mesh &M, const Rcpp::NumericMatrix &X);
void getWeights(GEO::OptimalTransportMap &OTM, double &wMax, double* w);
void getWeightedVerts(const Rcpp::NumericMatrix &X, double &wMax, double* w, double* wV);
Rcpp::NumericMatrix getVertices(GEO::Mesh &M);
Rcpp::NumericMatrix getVertices(GEO::Mesh &S, GEO::OptimalTransportMap &OTM);
Rcpp::NumericMatrix cubeVert(int d);

#endif