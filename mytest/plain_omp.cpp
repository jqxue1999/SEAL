#include <benchmark/benchmark.h>
#include "plain_omp.h"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <random>

using namespace std;

void cvps_omp(const vector<uint64>& x, uint64 scalar, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        result[i] = x[i] * scalar;
}
void cvpv_omp(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        cvps_omp(y, x[i], result[i]);
}
void pvcm_omp(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        vector<uint64> partial;
        cvps_omp(Y[i], x[i], partial);
        for (size_t j = 0; j < Y[i].size(); ++j)
            result[j] += partial[j];
    }
}
void pmcm_omp(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result) {
    result.resize(X.size());
    #pragma omp parallel for
    for (size_t i = 0; i < X.size(); ++i)
        pvcm_omp(X[i], Y, result[i]);
}

// Google Benchmark for cvps_omp
static void BM_cvps_omp(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    uint64 scalar = 123;
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        cvps_omp(x, scalar, result);
    }
}


// Google Benchmark for cvpv_omp
static void BM_cvpv_omp(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim), y(dim);
    for (auto& v : x) v = dist(rng);
    for (auto& v : y) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        cvpv_omp(x, y, result);
    }
}


// Google Benchmark for pvcm_omp
static void BM_pvcm_omp(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<vector<uint64>> Y(dim, vector<uint64>(dim));
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        pvcm_omp(x, Y, result);
    }
}


// Google Benchmark for pmcm_omp
static void BM_pmcm_omp(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<vector<uint64>> X(dim, vector<uint64>(dim)), Y(dim, vector<uint64>(dim));
    for (auto& row : X) for (auto& v : row) v = dist(rng);
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        pmcm_omp(X, Y, result);
    }
}

// Google Benchmark for cvps_omp
BENCHMARK(BM_cvps_omp)
    ->Arg(1024)
    ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv_omp
BENCHMARK(BM_cvpv_omp)
    ->Arg(1024)
    ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm_omp
BENCHMARK(BM_pvcm_omp)
    ->Arg(1024)
    ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm_omp
BENCHMARK(BM_pmcm_omp)
    ->Arg(1024)
    ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN(); 