#include <benchmark/benchmark.h>
#include "seal/seal.h"
#include "coeff_modulus_config.h"
#include <vector>
#include <random>
#include <cblas.h>
#include <omp.h>
#include "cpmm.h"

using namespace seal;
using namespace std;

static void BM_encrypt_matrix_rowwise(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0)); // 作为矩阵大小n
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    BatchEncoder encoder(context);

    // 随机生成矩阵
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(encoder.slot_count(), 0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64_t> dist(0, parms.plain_modulus().value() - 1);
    for (size_t i = 0; i < poly_modulus_degree; ++i) {
        for (size_t j = 0; j < encoder.slot_count(); ++j) {
            matrix[i][j] = dist(rng);
        }
    }

    vector<Ciphertext> encrypted_rows(poly_modulus_degree);
    for (auto _ : state) {
        #pragma omp parallel for
        for (size_t i = 0; i < poly_modulus_degree; ++i) {
            Plaintext plain;
            encoder.encode(matrix[i], plain);
            encryptor.encrypt(plain, encrypted_rows[i]);
        }
    }
}


static void BM_encrypt_matrix_rowwise_coeff(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0)); // 作为矩阵大小n
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    BatchEncoder encoder(context);

    // 随机生成矩阵
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(encoder.slot_count(), 0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64_t> dist(0, parms.plain_modulus().value() - 1);
    for (size_t i = 0; i < poly_modulus_degree; ++i) {
        for (size_t j = 0; j < encoder.slot_count(); ++j) {
            matrix[i][j] = dist(rng);
        }
    }

    for (auto _ : state) {
        vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
        #pragma omp parallel for
        for (size_t i = 0; i < poly_modulus_degree; i++) {
            Plaintext plaintext;
            encode_vector_to_plaintext(matrix[i], context, plaintext);
            encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
        }
    }
}

// BLAS矩阵减法
static void BM_matrix_sub_blas(benchmark::State& state) {
    size_t n = static_cast<size_t>(state.range(0));
    size_t N = n * n;
    std::vector<double> A(N), B(N), C(N);
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1000.0);
    for (size_t i = 0; i < N; ++i) {
        A[i] = dist(rng);
        B[i] = dist(rng);
    }
    for (auto _ : state) {
        std::copy(A.begin(), A.end(), C.begin());
        // C = A - B
        cblas_daxpy(N, -1.0, B.data(), 1, C.data(), 1);
        benchmark::DoNotOptimize(C);
    }
}

// OpenMP for循环矩阵减法
static void BM_matrix_sub_omp(benchmark::State& state) {
    size_t n = static_cast<size_t>(state.range(0));
    size_t N = n * n;
    std::vector<double> A(N), B(N), C(N);
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1000.0);
    for (size_t i = 0; i < N; ++i) {
        A[i] = dist(rng);
        B[i] = dist(rng);
    }
    for (auto _ : state) {
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            C[i] = A[i] - B[i];
        }
        benchmark::DoNotOptimize(C);
    }
}

BENCHMARK(BM_encrypt_matrix_rowwise)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_encrypt_matrix_rowwise_coeff)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_matrix_sub_blas)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_matrix_sub_omp)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
