// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "test.h"

void print_first_k_elements(const vector<uint64_t>& vec, size_t k) {
    cout << "First " << k << " elements: ";
    for (size_t i = 0; i < min(k, vec.size()); i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}

void example_batch_single_multiply()
{
    // =============================== Cryptographic Setup ===============================
    // Set up encryption parameters
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    
    // Set plain_modulus to enable batching
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

    SEALContext context(parms);
    print_parameters(context);
    cout << endl;

    // Create keys
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    // Create encryptor, evaluator, and decryptor
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // Create batch encoder
    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    size_t row_size = slot_count / 2;

    cout << "slot_count: " << slot_count << endl;
    cout << "coeff_count: " << parms.poly_modulus_degree() << endl;

    // =============================== Prepare the clear vector clear_vector_a ===============================
    // Create a vector of values to batch encode
    vector<uint64_t> clear_vector_a(slot_count, 0ULL);

    clear_vector_a[0] = 1ULL;
    clear_vector_a[1] = 2ULL;
    clear_vector_a[2] = 3ULL;
    clear_vector_a[3] = 4ULL;
    clear_vector_a[4] = 5ULL;
    clear_vector_a[5] = 6ULL;
    clear_vector_a[6] = 7ULL;
    clear_vector_a[7] = 8ULL;

    print_first_k_elements(clear_vector_a, 8);

    // =============================== Encrypt the clear_vector_a using CRT encoding ===============================
    // Encode the matrix into a plaintext
    Plaintext plain_vector_a_CRT;
    batch_encoder.encode(clear_vector_a, plain_vector_a_CRT);

    // Encrypt the batched plaintext
    Ciphertext encrypted_vector_a_CRT;
    encryptor.encrypt(plain_vector_a_CRT, encrypted_vector_a_CRT);

    // =============================== Encrypt the clear_vector_a on the coefficient ===============================
    // Create a plaintext with coefficients directly from clear_vector_a
    Plaintext plain_vector_a_coeff;
    plain_vector_a_coeff.resize(poly_modulus_degree);
    plain_vector_a_coeff.set_zero();
    // Use modulo_poly_coeffs to encode the coefficients
    util::modulo_poly_coeffs(clear_vector_a.data(), slot_count, parms.plain_modulus(), plain_vector_a_coeff.data());
    cout << "original plaintext polynomial: " << plain_vector_a_coeff.to_string() << endl;

    // Encrypt the coefficient-encoded plaintext
    Ciphertext encrypted_vector_a_coeff;
    // =============================== Multiply the encrypted_vector_a_CRT with a vector encrypted by CRT encoding ===============================
    encryptor.encrypt(plain_vector_a_coeff, encrypted_vector_a_coeff);
    vector<uint64_t> clear_vector_b(slot_count, 0);
    for (size_t i = 0; i < slot_count; i++) {
        clear_vector_b[i] = rand() % parms.plain_modulus().value();
    }
    Plaintext plain_vector_b_CRT;
    batch_encoder.encode(clear_vector_b, plain_vector_b_CRT);

    // Multiply encrypted_vector_a_CRT with plain_vector_b_CRT and measure latency
    Ciphertext encrypted_result_CRT;
    auto start1 = chrono::high_resolution_clock::now();
    evaluator.multiply_plain(encrypted_vector_a_CRT, plain_vector_b_CRT, encrypted_result_CRT);
    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
    auto time_CRT = duration1.count();


    // =============================== Multiply the encrypted_vector_a_coeff with a scalar on the coefficient ===============================
    uint64_t scalar =  5;    
    MemoryPoolHandle pool = MemoryManager::GetPool();
    auto start2 = chrono::high_resolution_clock::now();
    util::negacyclic_multiply_poly_mono_coeffmod(
        encrypted_vector_a_coeff, encrypted_vector_a_coeff.size(), scalar, 0, parms.coeff_modulus(), encrypted_vector_a_coeff, pool
    );
    
    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(end2 - start2);
    auto time_coeff = duration2.count();

    // Decrypt the result and check if it is correct
    Plaintext plain_result_coeff;
    decryptor.decrypt(encrypted_vector_a_coeff, plain_result_coeff);
    cout << "decrypted plaintext polynomial: " << plain_result_coeff.to_string() << endl;


    cout << "time_CRT: " << time_CRT << " microseconds" << endl;
    cout << "time_coeff: " << time_coeff << " microseconds" << endl;

    // =============================== Multiply the encrypted_vector_a_coeff with a scalar on the coefficient ===============================

    cout << "-------------------------------- Multiply the encrypted_vector_a_coeff with a vector of scalars on the coefficient --------------------------------\n";
    cout << "scalar: " << scalar << endl;
    encryptor.encrypt(plain_vector_a_coeff, encrypted_vector_a_coeff);
    cout << "original plaintext polynomial: " << plain_vector_a_coeff.to_string() << endl;

    auto &context_data = *context.key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t N = poly_modulus_degree * encrypted_vector_a_coeff.coeff_modulus_size();
    std::vector<uint64_t> scale_vector(N, scalar);

    auto total_start = chrono::high_resolution_clock::now();
    chrono::microseconds total_mul_time(0);
    chrono::microseconds total_mod_time(0);

    // Process each polynomial in the ciphertext
    for (int i = 0; i < encrypted_vector_a_coeff.size(); i++) {
        uint64_t* coeffs = encrypted_vector_a_coeff.data(i);
        std::vector<uint64_t> temp(N);
        
        // Time the multiplication operation
        auto mul_start = chrono::high_resolution_clock::now();
        // Perform element-wise multiplication across all RNS layers
        transform(
            coeffs, coeffs + N,
            scale_vector.begin(),
            temp.begin(),
            [](uint64_t a, uint64_t b) { return a * b; }
        );
        auto mul_end = chrono::high_resolution_clock::now();
        total_mul_time += chrono::duration_cast<chrono::microseconds>(mul_end - mul_start);

        auto mod_start = chrono::high_resolution_clock::now();
        // Apply modulo operation for each RNS layer separately
        for (int j = 0; j < encrypted_vector_a_coeff.coeff_modulus_size(); j++) {
            uint64_t modulus = coeff_modulus[j].value();
            uint64_t* layer_start = temp.data() + j * encrypted_vector_a_coeff.poly_modulus_degree();
            transform(
                layer_start, layer_start + encrypted_vector_a_coeff.poly_modulus_degree(),
                layer_start,
                [modulus](uint64_t a) { return a % modulus; }
            );
        }
        auto mod_end = chrono::high_resolution_clock::now();
        total_mod_time += chrono::duration_cast<chrono::microseconds>(mod_end - mod_start);

        // Copy results back to the original coefficients
        copy(temp.begin(), temp.end(), coeffs);
    }

    auto total_end = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::microseconds>(total_end - total_start);

    cout << "Total multiplication time: " << total_mul_time.count() << " microseconds" << endl;
    cout << "Total modulo time: " << total_mod_time.count() << " microseconds" << endl;
    cout << "Total time: " << total_duration.count() << " microseconds" << endl;

    decryptor.decrypt(encrypted_vector_a_coeff, plain_vector_a_coeff);
    cout << "decrypted plaintext polynomial: " << plain_vector_a_coeff.to_string() << endl;
}

int main()
{
    example_batch_single_multiply();
}