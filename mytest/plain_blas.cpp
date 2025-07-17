#include <benchmark/benchmark.h>
#include "plain_blas.h"
#include <omp.h>
#include <cblas.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

void cvps_blas(const vector<uint64>& x, uint64 scalar, vector<uint64>& result) {
    std::vector<double> x_d(x.size());
    result.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        x_d[i] = static_cast<double>(x[i]);
    double s = static_cast<double>(scalar);
    cblas_dscal(x.size(), s, x_d.data(), 1);
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        result[i] = static_cast<uint64>(x_d[i]);
}
void cvpv_blas(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        cvps_blas(y, x[i], result[i]);
}
void pvcm_blas(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        vector<uint64> partial;
        cvps_blas(Y[i], x[i], partial);
        for (size_t j = 0; j < Y[i].size(); ++j)
            result[j] += partial[j];
    }
}
void pmcm_blas(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result) {
    result.resize(X.size());
    #pragma omp parallel for
    for (size_t i = 0; i < X.size(); ++i)
        pvcm_blas(X[i], Y, result[i]);
}

// Google Benchmark for cvps_blas
static void BM_cvps_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    uint64 scalar = 123;
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        cvps_blas(x, scalar, result);
    }
}
BENCHMARK(BM_cvps_blas)->Arg(4096)->Arg(8192)->Arg(16384)->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv_blas
static void BM_cvpv_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim), y(dim);
    for (auto& v : x) v = dist(rng);
    for (auto& v : y) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        cvpv_blas(x, y, result);
    }
}
BENCHMARK(BM_cvpv_blas)->Arg(4096)->Arg(8192)->Arg(16384)->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm_blas
static void BM_pvcm_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<vector<uint64>> Y(dim, vector<uint64>(dim));
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        pvcm_blas(x, Y, result);
    }
}
BENCHMARK(BM_pvcm_blas)->Arg(4096)->Arg(8192)->Arg(16384)->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm_blas
static void BM_pmcm_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<vector<uint64>> X(dim, vector<uint64>(dim)), Y(dim, vector<uint64>(dim));
    for (auto& row : X) for (auto& v : row) v = dist(rng);
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        pmcm_blas(X, Y, result);
    }
}
BENCHMARK(BM_pmcm_blas)->Arg(4096)->Arg(8192)->Arg(16384)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN(); 