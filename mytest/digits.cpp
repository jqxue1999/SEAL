#include "digits.h"
#include <omp.h>


void decompose_to_bit_vectors(
    const vector<uint64_t>& input_vector,
    vector<vector<uint64_t>>& bit_vectors,
    int num_bits)
{    
    size_t vector_size = input_vector.size();
    
    bit_vectors.resize(num_bits);
    for (int bit = 0; bit < num_bits; bit++) {
        bit_vectors[bit].resize(vector_size, 0);
    }
    
    #pragma omp parallel for
    for (size_t i = 0; i < vector_size; i++) {
        uint64_t value = input_vector[i];
        
        for (int bit = 0; bit < num_bits; bit++) {
            uint64_t bit_value = (value >> bit) & 1;
            bit_vectors[bit][i] = bit_value;
        }
    }
}

void compose_from_bit_vectors(
    const vector<vector<uint64_t>>& bit_vectors,
    vector<uint64_t>& output_vector,
    int num_bits,
    uint64_t plain_modulus_value)
{
    if (bit_vectors.empty() || bit_vectors[0].empty()) {
        cerr << "位级向量为空" << endl;
        return;
    }
    
    size_t vector_size = bit_vectors[0].size();
    
    // 检查所有位级向量的大小是否一致
    for (int bit = 0; bit < num_bits; bit++) {
        if (bit_vectors[bit].size() != vector_size) {
            cerr << "位级向量大小不一致" << endl;
            return;
        }
    }
    
    output_vector.resize(vector_size, 0);
    
    // 对每个位置，使用加权求和（支持非标准二进制）
    for (size_t i = 0; i < vector_size; i++) {
        uint64_t result = 0;
        
        // 加权求和：sum(d[i] * 2^i)
        for (int bit = 0; bit < num_bits; bit++) {
            uint64_t bit_value = bit_vectors[bit][i];
            result += bit_value * (1ULL << bit);
        }
        
        output_vector[i] = result % plain_modulus_value;
    }
}

void encrypt_bit_vectors(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& bit_vectors,
    vector<Ciphertext>& encrypted_bit_vectors,
    int num_bits)
{   
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
        
    encrypted_bit_vectors.resize(num_bits);
    
    for (int bit = 0; bit < num_bits; bit++) {
        const vector<uint64_t>& bit_vector = bit_vectors[bit];
        
        Plaintext plaintext;
        encode_vector_to_plaintext(bit_vector, context, plaintext);
        encryptor.encrypt(plaintext, encrypted_bit_vectors[bit]);
    }
}

void decrypt_bit_vectors(
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_bit_vectors,
    vector<vector<uint64_t>>& decrypted_bit_vectors,
    int num_bits)
{
    if (encrypted_bit_vectors.empty()) {
        cerr << "加密的位级向量为空" << endl;
        return;
    }    
    // 初始化解密向量数组
    decrypted_bit_vectors.resize(num_bits);
    
    // 对每个加密的位级向量进行解密
    for (int bit = 0; bit < num_bits; bit++) {        
        // 解密
        Plaintext plaintext;
        decryptor.decrypt(encrypted_bit_vectors[bit], plaintext);
        
        // 解码为向量
        vector<uint64_t> decrypted_vector;
        decode_plaintext_to_vector(plaintext, context, decrypted_vector);
        
        decrypted_bit_vectors[bit] = decrypted_vector;
    }
}

void multiply_by_power_of_2(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<Ciphertext>& encrypted_bit_vectors,
    int power,
    vector<Ciphertext>& result_vectors,
    int num_bits)
{
    if (encrypted_bit_vectors.empty()) {
        cerr << "加密的位级向量为空" << endl;
        return;
    }
    
    if (power < 0 || power >= num_bits) {
        cerr << "幂次方必须在0-" << (num_bits-1) << "范围内" << endl;
        return;
    }
    
    result_vectors.resize(num_bits);
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    vector<uint64_t> zero_vector(poly_modulus_degree, 0);
    
    Plaintext zero_plaintext;
    encode_vector_to_plaintext(zero_vector, context, zero_plaintext);
    
    Ciphertext zero_ciphertext;
    encryptor.encrypt(zero_plaintext, zero_ciphertext);
    
    #pragma omp parallel for
    for (int bit = 0; bit < num_bits; bit++) {
        result_vectors[bit] = zero_ciphertext;
    }
    
    #pragma omp parallel for
    for (int bit = 0; bit < num_bits; bit++) {
        int new_bit = bit + power;
        
        if (new_bit < num_bits)
            result_vectors[new_bit] = encrypted_bit_vectors[bit];
    }
}


vector<int> decompose_to_powers_of_2(uint64_t value) {
    vector<int> powers;
    int exponent = 0;

    while (value > 0) {
        if ((value & 1) == 1) {
            powers.push_back(exponent);
        }
        value >>= 1;
        exponent++;
    }

    return powers;
}

double initialize_zero_ciphertext(
    const SEALContext& context,
    Encryptor& encryptor,
    Ciphertext& zero_ciphertext)
{
    auto start = chrono::high_resolution_clock::now();
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    vector<uint64_t> zero_vector(poly_modulus_degree, 0);
    Plaintext zero_plaintext;
    encode_vector_to_plaintext(zero_vector, context, zero_plaintext);
    encryptor.encrypt(zero_plaintext, zero_ciphertext);
    auto end = chrono::high_resolution_clock::now();

    auto duration = chrono::duration<double>(end - start);

    return duration.count();
}

void cvps_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const Ciphertext& zero_ciphertext,
    const vector<Ciphertext>& encrypted_bit_vectors,
    uint64_t multiplier,
    vector<Ciphertext>& result_vectors,
    int num_bits)
{
    vector<int> powers = decompose_to_powers_of_2(multiplier);
    
    result_vectors.resize(num_bits);
    if (powers.empty()) {
        #pragma omp parallel for
        for (int bit = 0; bit < num_bits; bit++)
            result_vectors[bit] = zero_ciphertext;
    } else {
        int p_min = *min_element(powers.begin(), powers.end());
        #pragma omp parallel for
        for (int target_bit = 0; target_bit < num_bits; target_bit++) {
            int source_bit = target_bit - p_min;
            if (source_bit >= 0 && source_bit < num_bits)
                result_vectors[target_bit] = encrypted_bit_vectors[source_bit];
            else
                result_vectors[target_bit] = zero_ciphertext;

        }
        int number_of_add_inplace = 0;
        #pragma omp parallel for reduction(+:number_of_add_inplace)
        for (int target_bit = 0; target_bit < num_bits; target_bit++) {
            for (int power : powers) {
                if (power == p_min) continue;
                int source_bit2 = target_bit - power;
                if (source_bit2 >= 0 && source_bit2 < num_bits) {
                    evaluator.add_inplace(result_vectors[target_bit], encrypted_bit_vectors[source_bit2]);
                    number_of_add_inplace++;
                }
            }
        }
    }
}

void cvpv_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    Decryptor& decryptor,
    int num_bits,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& bit_vectors_ciphertext,
    vector<vector<Ciphertext>>& outer_product_results)
{
    size_t poly_modulus_degree = clear_vector.size();

    outer_product_results.resize(clear_vector.size());
    Ciphertext zero_ciphertext;

    #pragma omp parallel for
    for (size_t i = 0; i < clear_vector.size(); ++i) {
        uint64_t scalar = clear_vector[i];
        vector<Ciphertext> result_vectors;
        cvps_digits(context, encryptor, evaluator, zero_ciphertext, bit_vectors_ciphertext, scalar, result_vectors, num_bits);
        outer_product_results[i] = result_vectors;
    }
}

double pvcm_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<vector<Ciphertext>>& encrypted_matrix,
    vector<Ciphertext>& result,
    int num_bits)
{
    auto start_time = chrono::high_resolution_clock::now();
    size_t n = clear_vector.size();
    size_t m = encrypted_matrix[0].size();

    Ciphertext zero_ciphertext;
    initialize_zero_ciphertext(context, encryptor, zero_ciphertext);
    result.resize(m);
    for (size_t j = 0; j < m; ++j) {
        result[j] = zero_ciphertext;
    }
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        uint64_t scalar = clear_vector[i];
        const vector<Ciphertext>& row = encrypted_matrix[i];
        vector<Ciphertext> partial;
        cvps_digits(context, encryptor, evaluator, zero_ciphertext, row, scalar, partial, num_bits);
        for (size_t j = 0; j < m; ++j) {
            evaluator.add_inplace(result[j], partial[j]);
        }
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>(end_time - start_time);
    return duration.count();
}

double pmcm_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const vector<vector<Ciphertext>>& encrypted_matrix,
    const vector<vector<uint64_t>>& clear_matrix,
    vector<vector<Ciphertext>>& result,
    int num_bits)
{
    auto start_time = chrono::high_resolution_clock::now();
    size_t n = encrypted_matrix.size();
    result.resize(n);

    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        pvcm_digits(context, encryptor, evaluator, clear_matrix[i], encrypted_matrix, result[i], num_bits);
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>(end_time - start_time);
    return duration.count();
}
