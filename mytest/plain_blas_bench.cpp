#include <benchmark/benchmark.h>
#include "plain_blas.h"
#include <random>

// Google Benchmark for cvps_blas
static void BM_cvps_plain_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    uint64 scalar = 123;
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        cvps_plain_blas(x, scalar, result);
    }
}

// Google Benchmark for cvpv_blas
static void BM_cvpv_plain_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim), y(dim);
    for (auto& v : x) v = dist(rng);
    for (auto& v : y) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        cvpv_plain_blas(x, y, result);
    }
}

// Google Benchmark for pvcm_blas
static void BM_pvcm_plain_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<uint64> x(dim);
    for (auto& v : x) v = dist(rng);
    vector<vector<uint64>> Y(dim, vector<uint64>(dim));
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<uint64> result;
    for (auto _ : state) {
        pvcm_plain_blas(x, Y, result);
    }
}

// Google Benchmark for pmcm_blas
static void BM_pmcm_plain_blas(benchmark::State& state) {
    size_t dim = static_cast<size_t>(state.range(0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64> dist(0, 1000);
    vector<vector<uint64>> X(dim, vector<uint64>(dim)), Y(dim, vector<uint64>(dim));
    for (auto& row : X) for (auto& v : row) v = dist(rng);
    for (auto& row : Y) for (auto& v : row) v = dist(rng);
    vector<vector<uint64>> result;
    for (auto _ : state) {
        pmcm_plain_blas(X, Y, result);
    }
}

// Google Benchmark for cvps_blas
BENCHMARK(BM_cvps_plain_blas)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv_blas
BENCHMARK(BM_cvpv_plain_blas)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm_blas
BENCHMARK(BM_pvcm_plain_blas)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm_blas
BENCHMARK(BM_pmcm_plain_blas)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN(); 