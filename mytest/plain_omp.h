#ifndef PLAIN_OMP_H
#define PLAIN_OMP_H

#include <vector>
#include <cstdint>
#include <omp.h>
#include <algorithm>
#include <random>

using std::vector;
using uint64 = uint64_t;

void cvps_plain_omp(const vector<uint64>& x, uint64 scalar, vector<uint64>& result);
void cvpv_plain_omp(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result);
void pvcm_plain_omp(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result);
void pmcm_plain_omp(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result);

#endif // PLAIN_OMP_H 