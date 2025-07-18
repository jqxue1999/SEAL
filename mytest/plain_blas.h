#ifndef PLAIN_BLAS_H
#define PLAIN_BLAS_H

#include <vector>
#include <cstdint>
#include <cblas.h>
#include <omp.h>

using std::vector;

void cvps_plain_blas(const vector<double>& x, double scalar, vector<double>& result);
void cvpv_plain_blas(const vector<double>& x, const vector<double>& y, vector<vector<double>>& result);
void pvcm_plain_blas(const vector<double>& x, const vector<vector<double>>& Y, vector<double>& result);
void pmcm_plain_blas(const vector<vector<double>>& X, const vector<vector<double>>& Y, vector<vector<double>>& result);

#endif // PLAIN_BLAS_H 