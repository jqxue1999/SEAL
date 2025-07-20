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

// BLAS矩阵乘法
static void BM_matrix_mul_blas(benchmark::State& state) {
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
        // C = A * B (C = A x B, row-major)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A.data(), n, B.data(), n, 0.0, C.data(), n);
        benchmark::DoNotOptimize(C);
    }
}

// 密文矩阵加法（逐行）
static void BM_encrypted_matrix_add(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
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
    Evaluator evaluator(context);

    // 随机生成两个明文矩阵
    vector<vector<uint64_t>> matrix1(poly_modulus_degree, vector<uint64_t>(encoder.slot_count(), 0));
    vector<vector<uint64_t>> matrix2(poly_modulus_degree, vector<uint64_t>(encoder.slot_count(), 0));
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64_t> dist(0, parms.plain_modulus().value() - 1);
    for (size_t i = 0; i < poly_modulus_degree; ++i) {
        for (size_t j = 0; j < encoder.slot_count(); ++j) {
            matrix1[i][j] = dist(rng);
            matrix2[i][j] = dist(rng);
        }
    }

    // 按行加密
    vector<Ciphertext> encrypted_matrix1(poly_modulus_degree), encrypted_matrix2(poly_modulus_degree);
    #pragma omp parallel for
    for (size_t i = 0; i < poly_modulus_degree; ++i) {
        Plaintext plain1, plain2;
        encode_vector_to_plaintext(matrix1[i], context, plain1);
        encode_vector_to_plaintext(matrix2[i], context, plain2);
        encryptor.encrypt(plain1, encrypted_matrix1[i]);
        encryptor.encrypt(plain2, encrypted_matrix2[i]);
    }

    for (auto _ : state) {
        std::vector<Ciphertext> encrypted_result(poly_modulus_degree);
        #pragma omp parallel for
        for (size_t i = 0; i < poly_modulus_degree; ++i) {
            encrypted_result[i] = encrypted_matrix1[i];
            evaluator.add_inplace(encrypted_result[i], encrypted_matrix2[i]);
        }
        benchmark::DoNotOptimize(encrypted_result);
    }
}

// 密文矩阵解密+解码测试
static void BM_decrypt_matrix_rowwise_coeff(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
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

    // 统计 matrix 的总大小（以MB为单位）
    size_t total_bytes_matrix = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
        total_bytes_matrix += matrix[i].size() * sizeof(uint64_t);
    }
    double total_MB_matrix = static_cast<double>(total_bytes_matrix) / (1024.0 * 1024.0);
    std::cout << "matrix total size: " << total_MB_matrix << " MB" << std::endl;

    // 先加密
    vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
    #pragma omp parallel for
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        Plaintext plaintext;
        encode_vector_to_plaintext(matrix[i], context, plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
    }

    // 统计 encrypted_matrix_coeff 的总大小（以MB为单位）
    size_t total_bytes_encrypted_matrix_coeff = 0;
    for (size_t i = 0; i < encrypted_matrix_coeff.size(); ++i) {
        total_bytes_encrypted_matrix_coeff += encrypted_matrix_coeff[i].save_size(compr_mode_type::none);
    }
    double total_MB_encrypted_matrix_coeff = static_cast<double>(total_bytes_encrypted_matrix_coeff) / (1024.0 * 1024.0);
    std::cout << "encrypted_matrix_coeff total size: " << total_MB_encrypted_matrix_coeff << " MB" << std::endl;

    for (auto _ : state) {
        vector<vector<uint64_t>> decoded_matrix(poly_modulus_degree);
        #pragma omp parallel for
        for (size_t i = 0; i < poly_modulus_degree; i++) {
            Plaintext pt;
            decryptor.decrypt(encrypted_matrix_coeff[i], pt);
            decode_plaintext_to_vector(pt, context, decoded_matrix[i]);
        }
        benchmark::DoNotOptimize(decoded_matrix);
    }
}

// BENCHMARK(BM_encrypt_matrix_rowwise_coeff)
//     ->Arg(4096)
//     ->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_encrypted_matrix_add)
//     ->Arg(4096)
//     ->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_decrypt_matrix_rowwise_coeff)
//     ->Arg(4096)
//     ->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_matrix_sub_blas)
//     ->Arg(4096)
//     ->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_matrix_sub_omp)
//     ->Arg(4096)
//     ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_matrix_mul_blas)
    ->Arg(4096)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
