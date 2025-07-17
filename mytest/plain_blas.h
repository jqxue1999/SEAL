#ifndef PLAIN_BLAS_H
#define PLAIN_BLAS_H

#include <vector>
#include <cstdint>

using std::vector;
using uint64 = uint64_t;

void cvps_blas(const vector<uint64>& x, uint64 scalar, vector<uint64>& result);
void cvpv_blas(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result);
void pvcm_blas(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result);
void pmcm_blas(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result);

#endif // PLAIN_BLAS_H 