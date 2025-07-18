#include "coefficients_seal.h"
#include <benchmark/benchmark.h>

// Google Benchmark for cvps
static void BM_cvps_coefficients_seal(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
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

    vector<uint64_t> input_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        input_vector[i] = rand();
    }
    uint64_t scalar = rand();
    Plaintext plaintext;
    encode_vector_to_plaintext(input_vector, context, plaintext);
    Ciphertext encrypted_vector_coeff;
    encryptor.encrypt(plaintext, encrypted_vector_coeff);
    Ciphertext encrypted_vector_coeff_result;

    for (auto _ : state) {
        cvps_coefficients_seal(context, scalar, encrypted_vector_coeff, encrypted_vector_coeff_result);
        benchmark::DoNotOptimize(encrypted_vector_coeff_result);
    }
}

// Google Benchmark for cvpv
static void BM_cvpv_coefficients_seal(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
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

    vector<uint64_t> input_vector(poly_modulus_degree);
    vector<uint64_t> clear_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        input_vector[i] = rand();
        clear_vector[i] = rand();
    }
    Plaintext plaintext;
    encode_vector_to_plaintext(input_vector, context, plaintext);
    Ciphertext encrypted_vector_coeff;
    encryptor.encrypt(plaintext, encrypted_vector_coeff);
    vector<Ciphertext> encrypted_vector_coeff_result;

    for (auto _ : state) {
        cvpv_coefficients_seal(context, clear_vector, encrypted_vector_coeff, encrypted_vector_coeff_result);
        benchmark::DoNotOptimize(encrypted_vector_coeff_result);
    }
}

// Google Benchmark for pvcm
static void BM_pvcm_coefficients_seal(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
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

    vector<uint64_t> clear_vector(poly_modulus_degree);
    vector<vector<uint64_t>> input_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        clear_vector[i] = rand();
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            input_matrix[i][j] = rand();
        }
    }
    Plaintext plaintext;
    vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        encode_vector_to_plaintext(input_matrix[i], context, plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
    }
    Ciphertext encrypted_matrix_coeff_result;

    for (auto _ : state) {
        pvcm_coefficients_seal(context, evaluator, clear_vector, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
        benchmark::DoNotOptimize(encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for pmcm
static void BM_pmcm_coefficients_seal(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
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

    vector<vector<uint64_t>> clear_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    vector<vector<uint64_t>> input_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            clear_matrix[i][j] = rand();
            input_matrix[i][j] = rand();
        }
    }
    Plaintext plaintext;
    vector<Ciphertext> encrypted_matrix_coeff(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        encode_vector_to_plaintext(input_matrix[i], context, plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix_coeff[i]);
    }
    vector<Ciphertext> encrypted_matrix_coeff_result;

    for (auto _ : state) {
        pmcm_coefficients_seal(context, evaluator, clear_matrix, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
        benchmark::DoNotOptimize(encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for cvps
BENCHMARK(BM_cvps_coefficients_seal)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv
BENCHMARK(BM_cvpv_coefficients_seal)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm
BENCHMARK(BM_pvcm_coefficients_seal)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm
BENCHMARK(BM_pmcm_coefficients_seal)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096)
    ->Arg(8192)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();