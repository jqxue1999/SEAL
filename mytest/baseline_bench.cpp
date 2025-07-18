#include <benchmark/benchmark.h>
#include "baseline.h"
#include <iostream>
#include "coeff_modulus_config.h"

// Google Benchmark for cvps baseline
static void BM_cvps_baseline(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    uint64_t plain_modulus_value = parms.plain_modulus().value();

    vector<uint64_t> input_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        input_vector[i] = abs(rand()) % plain_modulus_value / 2 + 1;
    }
    uint64_t scalar = abs(rand()) % plain_modulus_value / 2 + 1;
    Plaintext plaintext;
    batch_encoder.encode(input_vector, plaintext);
    Ciphertext encrypted_vector_coeff;
    encryptor.encrypt(plaintext, encrypted_vector_coeff);
    Ciphertext encrypted_vector_coeff_result;

    for (auto _ : state) {
        cvps_baseline(context, batch_encoder, evaluator, scalar, encrypted_vector_coeff, encrypted_vector_coeff_result);
    }
}

// Google Benchmark for cvpv baseline
static void BM_cvpv_baseline(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    uint64_t plain_modulus_value = parms.plain_modulus().value();

    vector<uint64_t> input_vector(poly_modulus_degree);
    vector<uint64_t> clear_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        input_vector[i] = abs(rand()) % plain_modulus_value / 2 + 1;
        clear_vector[i] = abs(rand()) % plain_modulus_value / 2 + 1;
    }
    Plaintext plaintext;
    batch_encoder.encode(input_vector, plaintext);
    Ciphertext encrypted_vector_coeff;
    encryptor.encrypt(plaintext, encrypted_vector_coeff);
    vector<Ciphertext> encrypted_vector_coeff_result;

    for (auto _ : state) {
        cvpv_baseline(context, batch_encoder, evaluator, clear_vector, encrypted_vector_coeff, encrypted_vector_coeff_result);
    }
}

// Google Benchmark for pvcm baseline
static void BM_pvcm_baseline(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    uint64_t plain_modulus_value = parms.plain_modulus().value();

    vector<uint64_t> clear_vector(poly_modulus_degree);
    vector<vector<uint64_t>> input_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        clear_vector[i] = abs(rand()) % plain_modulus_value / 2 + 1;
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            input_matrix[i][j] = abs(rand()) % plain_modulus_value / 2 + 1;
        }
    }
    Plaintext plaintext;
    vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        batch_encoder.encode(input_matrix[i], plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
    }
    Ciphertext encrypted_matrix_coeff_result;

    for (auto _ : state) {
        pvcm_baseline(context, batch_encoder, evaluator, clear_vector, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for pmcm baseline
static void BM_pmcm_baseline(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);

    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    uint64_t plain_modulus_value = parms.plain_modulus().value();

    vector<vector<uint64_t>> clear_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    vector<vector<uint64_t>> input_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            clear_matrix[i][j] = abs(rand()) % plain_modulus_value / 2 + 1;
            input_matrix[i][j] = abs(rand()) % plain_modulus_value / 2 + 1;
        }
    }
    Plaintext plaintext;
    vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        batch_encoder.encode(input_matrix[i], plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
    }
    vector<Ciphertext> encrypted_matrix_coeff_result;

    for (auto _ : state) {
        pmcm_baseline(context, batch_encoder, evaluator, clear_matrix, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for cvps baseline
BENCHMARK(BM_cvps_baseline)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv baseline
BENCHMARK(BM_cvpv_baseline)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm baseline
BENCHMARK(BM_pvcm_baseline)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm baseline
BENCHMARK(BM_pmcm_baseline)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();