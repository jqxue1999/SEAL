#include <benchmark/benchmark.h>
#include "digits.h"
#include <iostream>
#include "coeff_modulus_config.h"

static void BM_cvps_digits(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    int num_bits = static_cast<int>(state.range(1));
    
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    uint64_t plain_modulus_value = 1ULL << num_bits;
    
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    vector<uint64_t> input_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        input_vector[i] = rand() % plain_modulus_value;
    }
    uint64_t multiplier = rand() % plain_modulus_value;
    vector<vector<uint64_t>> bit_vectors;
    decompose_to_bit_vectors(input_vector, bit_vectors, num_bits);
    vector<Ciphertext> encrypted_bit_vectors;
    encrypt_bit_vectors(context, encryptor, bit_vectors, encrypted_bit_vectors, num_bits);
    Ciphertext zero_ciphertext;
    initialize_zero_ciphertext(context, encryptor, zero_ciphertext);

    for (auto _ : state) {
        vector<Ciphertext> result_vectors;
        cvps_digits(context, encryptor, evaluator, zero_ciphertext, encrypted_bit_vectors, multiplier, result_vectors, num_bits);
    }
}

static void BM_cvpv_digits(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    int num_bits = static_cast<int>(state.range(1));
    
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    uint64_t plain_modulus_value = 1ULL << num_bits;
    
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    // clear_vector
    vector<uint64_t> clear_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++)
        clear_vector[i] = rand() % plain_modulus_value;
    // plain_vector
    vector<uint64_t> plain_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++)
        plain_vector[i] = rand() % plain_modulus_value;
    // bit_vectors_plaintext
    vector<vector<uint64_t>> bit_vectors_plaintext;
    decompose_to_bit_vectors(plain_vector, bit_vectors_plaintext, num_bits);
    // encrypted_bit_vectors
    vector<Ciphertext> encrypted_bit_vectors;
    encrypt_bit_vectors(context, encryptor, bit_vectors_plaintext, encrypted_bit_vectors, num_bits);

    for (auto _ : state) {
        vector<vector<Ciphertext>> outer_product_results;
        cvpv_digits(context, encryptor, evaluator, decryptor, num_bits, clear_vector, encrypted_bit_vectors, outer_product_results);
    }
}

static void BM_pvcm_digits(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    int num_bits = static_cast<int>(state.range(1));
    
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    uint64_t plain_modulus_value = 1ULL << num_bits;
    
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    // clear_vector
    vector<uint64_t> clear_vector(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++)
        clear_vector[i] = rand() % 2;
    // plain_matrix
    vector<vector<uint64_t>> plain_matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++)
        for (size_t j = 0; j < poly_modulus_degree; j++)
            plain_matrix[i][j] = rand() % 2;
    // encrypted_matrix
    vector<vector<Ciphertext>> encrypted_matrix(poly_modulus_degree);
    #pragma omp parallel for
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        vector<vector<uint64_t>> bit_vectors_plaintext;
        decompose_to_bit_vectors(plain_matrix[i], bit_vectors_plaintext, num_bits);
        encrypt_bit_vectors(context, encryptor, bit_vectors_plaintext, encrypted_matrix[i], num_bits);
    }

    for (auto _ : state) {
        vector<Ciphertext> result;
        pvcm_digits(context, encryptor, evaluator, clear_vector, encrypted_matrix, result, num_bits);
    }
}

static void BM_pmcm_digits(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    int num_bits = static_cast<int>(state.range(1));
    
    std::vector<int> coeff_modulus_params = CoeffModulusConfig::get_coeff_modulus_params(poly_modulus_degree);
    scheme_type scheme = scheme_type::bfv;
    EncryptionParameters parms(scheme);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    uint64_t plain_modulus_value = 1ULL << num_bits;
    
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    // plain_matrix_A
    vector<vector<uint64_t>> plain_matrix_A(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            plain_matrix_A[i][j] = rand() % 2;
        }
    }
    // plain_matrix_B
    vector<vector<uint64_t>> plain_matrix_B(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        for (size_t j = 0; j < poly_modulus_degree; j++) {
            plain_matrix_B[i][j] = rand() % 2;
        }
    }
    // encrypted_matrix_A
    vector<vector<Ciphertext>> encrypted_matrix_A(poly_modulus_degree);
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        vector<vector<uint64_t>> bit_vectors_plaintext;
        decompose_to_bit_vectors(plain_matrix_A[i], bit_vectors_plaintext, num_bits);
        encrypt_bit_vectors(context, encryptor, bit_vectors_plaintext, encrypted_matrix_A[i], num_bits);
    }

    for (auto _ : state) {
        vector<vector<Ciphertext>> result;
        pmcm_digits(context, encryptor, evaluator, encrypted_matrix_A, plain_matrix_B, result, num_bits);
    }
}

BENCHMARK(BM_cvps_digits)
    ->Args({1024, 8})
    ->Args({2048, 8})
    ->Args({4096, 8})
    ->Args({8192, 8})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_cvpv_digits)
    ->Args({1024, 8})
    ->Args({2048, 8})
    ->Args({4096, 8})
    ->Args({8192, 8})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_pvcm_digits)
    ->Args({1024, 8})
    ->Args({2048, 8})
    ->Args({4096, 8})
    ->Args({8192, 8})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_pmcm_digits)
    ->Args({1024, 8})
    ->Args({2048, 8})
    ->Args({4096, 8})
    ->Args({8192, 8})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN(); 