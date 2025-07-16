#include "plaintext.h"
#include <omp.h>
#include <cblas.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

// 1. Vector-scalar multiplication
vector<uint64> vec_scalar_mul_blas(const vector<uint64>& x, uint64 scalar) {
    size_t n = x.size();
    std::vector<double> x_d(n);
    for (size_t i = 0; i < n; ++i) x_d[i] = static_cast<double>(x[i]);
    double s = static_cast<double>(scalar);
    cblas_dscal(n, s, x_d.data(), 1);
    vector<uint64> result(n);
    for (size_t i = 0; i < n; ++i) result[i] = static_cast<uint64>(x_d[i]);
    return result;
}

vector<uint64> vec_scalar_mul_omp(const vector<uint64>& x, uint64 scalar) {
    size_t n = x.size();
    vector<uint64> result(n);
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        result[i] = x[i] * scalar;
    }
    return result;
}

// 2. Vector-vector outer product
vector<vector<uint64>> vec_outer_product_blas(const vector<uint64>& x, const vector<uint64>& y) {
    size_t m = x.size(), n = y.size();
    std::vector<double> x_d(m), y_d(n);
    for (size_t i = 0; i < m; ++i) x_d[i] = static_cast<double>(x[i]);
    for (size_t j = 0; j < n; ++j) y_d[j] = static_cast<double>(y[j]);
    std::vector<double> out_d(m * n);
    cblas_dger(CblasRowMajor, m, n, 1.0, x_d.data(), 1, y_d.data(), 1, out_d.data(), n);
    vector<vector<uint64>> result(m, vector<uint64>(n));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            result[i][j] = static_cast<uint64>(out_d[i * n + j]);
    return result;
}

vector<vector<uint64>> vec_outer_product_omp(const vector<uint64>& x, const vector<uint64>& y) {
    size_t m = x.size(), n = y.size();
    vector<vector<uint64>> result(m, vector<uint64>(n));
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            result[i][j] = x[i] * y[j];
    return result;
}

// 3. Vector-matrix multiplication (x * Y, Y is matrix with cols)
vector<uint64> vec_mat_mul_blas(const vector<uint64>& x, const vector<vector<uint64>>& Y) {
    size_t m = x.size();
    if (m == 0) return {};
    size_t n = Y[0].size();
    std::vector<double> x_d(m);
    std::vector<double> Y_d(m * n);
    for (size_t i = 0; i < m; ++i) x_d[i] = static_cast<double>(x[i]);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            Y_d[i * n + j] = static_cast<double>(Y[i][j]);
    std::vector<double> res_d(n, 0.0);
    cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, Y_d.data(), n, x_d.data(), 1, 0.0, res_d.data(), 1);
    vector<uint64> result(n);
    for (size_t j = 0; j < n; ++j) result[j] = static_cast<uint64>(res_d[j]);
    return result;
}

vector<uint64> vec_mat_mul_omp(const vector<uint64>& x, const vector<vector<uint64>>& Y) {
    size_t m = x.size();
    if (m == 0) return {};
    size_t n = Y[0].size();
    vector<uint64> result(n, 0);
    #pragma omp parallel for
    for (size_t j = 0; j < n; ++j) {
        uint64 sum = 0;
        for (size_t i = 0; i < m; ++i) {
            sum += x[i] * Y[i][j];
        }
        result[j] = sum;
    }
    return result;
}

// 4. Matrix-matrix multiplication (X * Y)
vector<vector<uint64>> mat_mat_mul_blas(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y) {
    size_t m = X.size();
    if (m == 0) return {};
    size_t k = X[0].size();
    size_t n = Y[0].size();
    std::vector<double> X_d(m * k);
    std::vector<double> Y_d(k * n);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < k; ++j)
            X_d[i * k + j] = static_cast<double>(X[i][j]);
    for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < n; ++j)
            Y_d[i * n + j] = static_cast<double>(Y[i][j]);
    std::vector<double> res_d(m * n, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, X_d.data(), k, Y_d.data(), n, 0.0, res_d.data(), n);
    vector<vector<uint64>> result(m, vector<uint64>(n));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            result[i][j] = static_cast<uint64>(res_d[i * n + j]);
    return result;
}

vector<vector<uint64>> mat_mat_mul_omp(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y) {
    size_t m = X.size();
    if (m == 0) return {};
    size_t k = X[0].size();
    size_t n = Y[0].size();
    vector<vector<uint64>> result(m, vector<uint64>(n, 0));
    #pragma omp parallel for
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            uint64 sum = 0;
            for (size_t l = 0; l < k; ++l) {
                sum += X[i][l] * Y[l][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

void flush_cache() {
    static const size_t cache_size = 32 * 1024 * 1024; // 32MB
    static std::vector<char> cache(cache_size, 1);
    volatile char sum = 0;
    for (size_t i = 0; i < cache_size; ++i) {
        sum += cache[i];
    }
}

template<typename F, typename... Args>
double bench(F func, int repeat, Args&&... args) {
    double total = 0.0;
    for (int i = 0; i < repeat; ++i) {
        // cout << "Running iteration " << i + 1 << " of " << repeat << endl;
        flush_cache();
        auto start = high_resolution_clock::now();
        volatile auto res = func(std::forward<Args>(args)...);
        auto end = high_resolution_clock::now();
        total += duration_cast<microseconds>(end - start).count() / 1000.0;
    }
    return total / repeat;
}

vector<uint64> random_vec(size_t n) {
    static mt19937_64 rng(42);
    uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> v(n);
    for (auto& x : v) x = dist(rng);
    return v;
}

vector<vector<uint64>> random_mat(size_t m, size_t n) {
    vector<vector<uint64>> mat(m);
    for (auto& row : mat) row = random_vec(n);
    return mat;
}

int main() {
    vector<size_t> dims = {8192, 16384};
    int repeat = 10;
    uint64 scalar = 123;
    cout << "Benchmarking vector/matrix operations (avg ms over " << repeat << " runs)\n";
    for (auto dim : dims) {
        double t1, t2;
        vector<uint64> x, y;
        vector<vector<uint64>> X, Y;
        cout << "\n==== Dimension: " << dim << " ====" << endl;
        // 1. Vector-scalar
        x = random_vec(dim);
        t1 = bench(vec_scalar_mul_blas, repeat, x, scalar);
        t2 = bench(vec_scalar_mul_omp, repeat, x, scalar);
        cout << "vec_scalar_mul_blas: " << t1 << " ms\n";
        cout << "vec_scalar_mul_omp:  " << t2 << " ms\n";
        // 2. Vector-vector outer product
        y = random_vec(dim);
        t1 = bench(vec_outer_product_blas, repeat, x, y);
        t2 = bench(vec_outer_product_omp, repeat, x, y);
        cout << "vec_outer_product_blas: " << t1 << " ms\n";
        cout << "vec_outer_product_omp:  " << t2 << " ms\n";
        // 3. Vector-matrix
        Y = random_mat(dim, dim);
        t1 = bench(vec_mat_mul_blas, repeat, x, Y);
        t2 = bench(vec_mat_mul_omp, repeat, x, Y);
        cout << "vec_mat_mul_blas: " << t1 << " ms\n";
        cout << "vec_mat_mul_omp:  " << t2 << " ms\n";
        // 4. Matrix-matrix
        X = random_mat(dim, dim);
        Y = random_mat(dim, dim);
        t1 = bench(mat_mat_mul_blas, repeat, X, Y);
        t2 = bench(mat_mat_mul_omp, repeat, X, Y);
        cout << "mat_mat_mul_blas: " << t1 << " ms\n";
        cout << "mat_mat_mul_omp:  " << t2 << " ms\n";
    }
    return 0;
}
