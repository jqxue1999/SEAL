#include "coefficients_blas.h"

void cvps(
    const SEALContext& context,
    const uint64_t scalar,
    const Ciphertext& encrypted_vector_coeff,
    Ciphertext& encrypted_vector_coeff_result
) {        
    scale_vector_blas(context, encrypted_vector_coeff, encrypted_vector_coeff_result, scalar);
}

void cvpv(
    const SEALContext& context,
    const vector<uint64_t>& clear_vector,
    const Ciphertext& encrypted_vector_coeff, 
    vector<Ciphertext>& encrypted_vector_coeff_result
) {
    encrypted_vector_coeff_result.resize(clear_vector.size());
    #pragma omp parallel for
    for (size_t i = 0; i < clear_vector.size(); ++i) {
        cvps(context, clear_vector[i], encrypted_vector_coeff, encrypted_vector_coeff_result[i]);
    }
}

void pvcm(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    Ciphertext& encrypted_matrix_coeff_result
) {
    cvps(context, clear_vector[0], encrypted_matrix_coeff[0], encrypted_matrix_coeff_result);

    #pragma omp parallel for
    for (size_t i = 1; i < clear_vector.size(); i++) {
        Ciphertext partial;
        cvps(context, clear_vector[i], encrypted_matrix_coeff[i], partial);
        evaluator.add_inplace(encrypted_matrix_coeff_result, partial);
    }
}

void pmcm(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<vector<uint64_t>>& clear_matrix,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    vector<Ciphertext>& encrypted_matrix_coeff_result
) {
    encrypted_matrix_coeff_result.resize(clear_matrix.size());
    #pragma omp parallel for
    for (size_t i = 0; i < clear_matrix.size(); i++)
        pvcm(context, evaluator, clear_matrix[i], encrypted_matrix_coeff, encrypted_matrix_coeff_result[i]);
}

// Google Benchmark for cvps
static void BM_cvps(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    json config = read_seal_config("/home/jiaq/Research/SEAL/mytest/seal_config.json");
    vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
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
        cvps(context, scalar, encrypted_vector_coeff, encrypted_vector_coeff_result);
    }
}

// Google Benchmark for cvpv
static void BM_cvpv(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    json config = read_seal_config("/home/jiaq/Research/SEAL/mytest/seal_config.json");
    vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
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
        cvpv(context, clear_vector, encrypted_vector_coeff, encrypted_vector_coeff_result);
    }
}

// Google Benchmark for pvcm
static void BM_pvcm(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    json config = read_seal_config("/home/jiaq/Research/SEAL/mytest/seal_config.json");
    vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
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
        pvcm(context, evaluator, clear_vector, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for pmcm
static void BM_pmcm(benchmark::State& state) {
    size_t poly_modulus_degree = static_cast<size_t>(state.range(0));
    json config = read_seal_config("/home/jiaq/Research/SEAL/mytest/seal_config.json");
    vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
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
        pmcm(context, evaluator, clear_matrix, encrypted_matrix_coeff, encrypted_matrix_coeff_result);
    }
}

// Google Benchmark for cvps
BENCHMARK(BM_cvps)
    ->Arg(1024)
    // ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for cvpv
BENCHMARK(BM_cvpv)
    ->Arg(1024)
    // ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pvcm
BENCHMARK(BM_pvcm)
    ->Arg(1024)
    // ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

// Google Benchmark for pmcm
BENCHMARK(BM_pmcm)
    ->Arg(1024)
    // ->Arg(2048)
    // ->Arg(4096)
    // ->Arg(8192)
    // ->Arg(16384)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();