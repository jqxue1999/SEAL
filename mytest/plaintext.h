#ifndef PLAINTEXT_H
#define PLAINTEXT_H

#include <vector>
#include <cstdint>

using std::vector;
using uint64 = uint64_t;

// 1. Vector-scalar multiplication
vector<uint64> vec_scalar_mul_blas(const vector<uint64>& x, uint64 scalar);
vector<uint64> vec_scalar_mul_omp(const vector<uint64>& x, uint64 scalar);

// 2. Vector-vector outer product
vector<vector<uint64>> vec_outer_product_blas(const vector<uint64>& x, const vector<uint64>& y);
vector<vector<uint64>> vec_outer_product_omp(const vector<uint64>& x, const vector<uint64>& y);

// 3. Vector-matrix multiplication
vector<uint64> vec_mat_mul_blas(const vector<uint64>& x, const vector<vector<uint64>>& Y);
vector<uint64> vec_mat_mul_omp(const vector<uint64>& x, const vector<vector<uint64>>& Y);

// 4. Matrix-matrix multiplication
vector<vector<uint64>> mat_mat_mul_blas(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y);
vector<vector<uint64>> mat_mat_mul_omp(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y);

#endif // PLAINTEXT_H
