#include "digits.h"
#include <omp.h>


void decompose_to_bit_vectors(
    const vector<uint64_t>& input_vector,
    vector<vector<uint64_t>>& bit_vectors,
    int num_bits)
{
    if (input_vector.empty()) {
        cerr << "输入向量为空" << endl;
        return;
    }
    
    size_t vector_size = input_vector.size();
    
    // 初始化64个位级向量
    bit_vectors.resize(num_bits);
    for (int bit = 0; bit < num_bits; bit++) {
        bit_vectors[bit].resize(vector_size, 0);
    }
    
    // 对每个输入整数，提取其64位（并行化）
    #pragma omp parallel for
    for (size_t i = 0; i < vector_size; i++) {
        uint64_t value = input_vector[i];
        
        // 提取每一位
        for (int bit = 0; bit < num_bits; bit++) {
            // 使用位掩码提取第bit位
            uint64_t bit_value = (value >> bit) & 1;
            bit_vectors[bit][i] = bit_value;
        }
    }
}

void compose_from_bit_vectors(
    const vector<vector<uint64_t>>& bit_vectors,
    vector<uint64_t>& output_vector,
    int num_bits)
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
        
        output_vector[i] = result;
    }
}

void encrypt_bit_vectors(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& bit_vectors,
    vector<Ciphertext>& encrypted_bit_vectors,
    int num_bits)
{
    if (bit_vectors.empty()) {
        cerr << "位级向量为空" << endl;
        return;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
        
    // 初始化加密向量数组 - 64个Ciphertext，每个对应一个位级
    encrypted_bit_vectors.resize(num_bits);
    
    // 对每个位级向量进行加密
    for (int bit = 0; bit < num_bits; bit++) {
        const vector<uint64_t>& bit_vector = bit_vectors[bit];
        
        // 检查向量大小是否匹配多项式模数次数
        if (bit_vector.size() != poly_modulus_degree) {
            cerr << "\n错误: 位级向量 " << bit << " 的大小 (" << bit_vector.size() 
                 << ") 与多项式模数次数 (" << poly_modulus_degree << ") 不匹配" << endl;
            return;
        }
        
        // 编码为明文多项式
        Plaintext plaintext;
        encode_vector_to_plaintext(bit_vector, context, plaintext);
        
        // 加密
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
    
    // 初始化结果向量
    result_vectors.resize(num_bits);
    
    // 通过交换向量位置来实现乘以2的幂次方
    // 例如：乘以2^3意味着所有位向左移动3位
    
    // 创建全零的明文多项式作为占位符
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    // 创建全零向量
    vector<uint64_t> zero_vector(poly_modulus_degree, 0);
    
    // 编码为明文多项式
    Plaintext zero_plaintext;
    encode_vector_to_plaintext(zero_vector, context, zero_plaintext);
    
    // 加密为全零密文
    Ciphertext zero_ciphertext;
    encryptor.encrypt(zero_plaintext, zero_ciphertext);
    
    // 首先，将所有位置初始化为全零密文
    #pragma omp parallel for
    for (int bit = 0; bit < num_bits; bit++) {
        result_vectors[bit] = zero_ciphertext;
    }
    
    // 然后，执行位级左移（并行化）
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

void multiply_by_general_scalar(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const vector<Ciphertext>& encrypted_bit_vectors,
    uint64_t multiplier,
    vector<Ciphertext>& result_vectors,
    int num_bits)
{
    if (encrypted_bit_vectors.empty()) {
        cerr << "加密的位级向量为空" << endl;
        return;
    }

    // 分解乘数为2的幂次方之和
    auto start = chrono::high_resolution_clock::now();
    vector<int> powers = decompose_to_powers_of_2(multiplier);
    cout << "乘数" << multiplier << "分解为: ";
    for (size_t i = 0; i < powers.size(); i++) {
        cout << "2^" << powers[i];
        if (i < powers.size() - 1) cout << " + ";
    }
    cout << endl;
    auto end = chrono::high_resolution_clock::now();
    auto step_1_duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    // 初始化全零密文（只在需要时用）
    start = chrono::high_resolution_clock::now();
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    vector<uint64_t> zero_vector(poly_modulus_degree, 0);
    Plaintext zero_plaintext;
    encode_vector_to_plaintext(zero_vector, context, zero_plaintext);
    Ciphertext zero_ciphertext;
    encryptor.encrypt(zero_plaintext, zero_ciphertext);
    end = chrono::high_resolution_clock::now();
    auto step_2_duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    start = chrono::high_resolution_clock::now();
    result_vectors.resize(num_bits);
    if (powers.empty()) {
        // multiplier为0，全部置零
        for (int bit = 0; bit < num_bits; bit++)
            result_vectors[bit] = zero_ciphertext;
        end = chrono::high_resolution_clock::now();
        auto step_3_duration = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "乘以" << multiplier << "操作完成！耗时: " << step_1_duration.count() << " + " 
        << step_2_duration.count() << " + " << step_3_duration.count() << " = " << 
        step_1_duration.count() + step_2_duration.count() + step_3_duration.count() << " microseconds" << endl;
        cout << "add_inplace 操作次数: 0" << endl;
    } else {
        int p_min = *min_element(powers.begin(), powers.end());
        for (int target_bit = 0; target_bit < num_bits; target_bit++) {
            int source_bit = target_bit - p_min;
            if (source_bit >= 0 && source_bit < num_bits) {
                result_vectors[target_bit] = encrypted_bit_vectors[source_bit];
            } else {
                result_vectors[target_bit] = zero_ciphertext;
            }
        }
        // 对剩下的 powers 做 add_inplace
        int number_of_add_inplace = 0;
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
        end = chrono::high_resolution_clock::now();
        auto step_3_duration = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "乘以" << multiplier << "操作完成！耗时: " << step_1_duration.count() << " + " 
        << step_2_duration.count() << " + " << step_3_duration.count() << " = " << 
        step_1_duration.count() + step_2_duration.count() + step_3_duration.count() << " microseconds" << endl;
        cout << "add_inplace 操作次数: " << number_of_add_inplace << endl;
    }
    
    return;
}

bool verify_general_multiplication(
    const vector<uint64_t>& input_vector,
    uint64_t multiplier,
    const vector<uint64_t>& output_vector)
{
    if (input_vector.size() != output_vector.size()) {
        cerr << "输入和输出向量大小不匹配" << endl;
        return false;
    }
        
    bool all_correct = true;
    
    for (size_t i = 0; i < input_vector.size(); i++) {
        uint64_t expected = input_vector[i] * multiplier;
        uint64_t actual = output_vector[i];
        
        if (expected != actual) {
            cout << "位置 " << i << ": 期望 " << expected 
                 << ", 实际 " << actual << endl;
            all_correct = false;
        }
    }
    
    if (all_correct) {
        cout << "验证成功！所有结果都正确。" << endl;
    } else {
        cout << "验证失败！存在错误结果。" << endl;
    }
    
    return all_correct;
} 