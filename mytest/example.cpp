#include "example.h"
#include <chrono>
#include <cstdlib>
#include <cblas.h>
#include <unistd.h>
#include "cpscale.h"

using namespace seal;
using namespace std;

json read_seal_config(const string& config_file) {
    json config;
    
    // 打印当前工作目录
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        cout << "当前工作目录: " << cwd << endl;
    } else {
        cout << "无法获取当前工作目录" << endl;
    }
    
    // 直接使用固定路径
    string config_path = "../../" + config_file;
    cout << "尝试读取配置文件: " << config_path << endl;
    
    try {
        ifstream file(config_path);
        if (file.is_open()) {
            file >> config;
            file.close();
            cout << "成功读取配置文件: " << config_path << endl;
            return config;
        } else {
            cerr << "错误: 无法打开配置文件 " << config_path << endl;
        }
    } catch (const exception& e) {
        cerr << "错误: 解析配置文件失败: " << e.what() << endl;
    }
    
    return config;
}

size_t get_user_poly_modulus_degree(const json& config) {
    vector<size_t> options = config["poly_modulus_degree_options"].get<vector<size_t>>();
    
    cout << "\n=== 选择多项式模数次数 ===" << endl;
    cout << "可选的次数: ";
    for (size_t i = 0; i < options.size(); i++) {
        cout << options[i];
        if (i < options.size() - 1) cout << ", ";
    }
    cout << endl;
    
    while (true) {
        cout << "请输入多项式模数次数: ";
        size_t user_choice;
        cin >> user_choice;
        
        // 检查输入是否在有效范围内
        for (size_t option : options) {
            if (option == user_choice) {
                cout << "已选择: " << user_choice << endl;
                return user_choice;
            }
        }
        
        cout << "无效选择！请输入以下值之一: ";
        for (size_t i = 0; i < options.size(); i++) {
            cout << options[i];
            if (i < options.size() - 1) cout << ", ";
        }
        cout << endl;
    }
}

vector<int> get_coeff_modulus_params(const json& config, size_t poly_modulus_degree) {
    string degree_str = to_string(poly_modulus_degree);
    if (config["coeff_modulus_configs"].contains(degree_str)) {
        return config["coeff_modulus_configs"][degree_str]["coeff_modulus"].get<vector<int>>();
    } else {
        cerr << "错误: 配置文件中未找到多项式模数次数 " << poly_modulus_degree << " 的参数" << endl;
        return vector<int>();
    }
}

void print_menu() {
    cout << "\n=== SEAL 矩阵乘法测试程序 ===" << endl;
    cout << "请选择要测试的功能：" << endl;
    cout << "1. CPMM (Ciphertext-Plaintext Matrix Multiplication) - 密文-明文矩阵乘法" << endl;
    cout << "2. CCMM (Ciphertext-Ciphertext Matrix Multiplication) - 密文-密文矩阵乘法" << endl;
    cout << "3. CMT (Ciphertext Matrix Transpose) - 密文矩阵转置" << endl;
    cout << "4. BLAS Performance Test - 测试BLAS浮点矩阵乘法性能" << endl;
    cout << "5. General Multiplication & Carry Recovery - 通用乘法和进位恢复" << endl;
    cout << "6. Ciphertext Scale Multiplication - 密文scale乘法" << endl;
    cout << "7. 退出程序" << endl;
    cout << "请输入选择 (1-7): ";
}

int test_cpmm() {
    try {
        // 读取配置文件
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "无法读取配置文件，使用默认参数" << endl;
            return 1;
        }
        
        // 获取用户输入的多项式模数次数
        size_t poly_modulus_degree = get_user_poly_modulus_degree(config);
        
        // 获取对应的系数模数参数
        vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
        if (coeff_modulus_params.empty()) {
            cerr << "无法获取系数模数参数" << endl;
            return 1;
        }
        
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));

        uint64_t plain_modulus_value = 1;
        if (scheme == scheme_type::ckks) {
            for (const auto &mod : parms.coeff_modulus())
                plain_modulus_value *= mod.value();
        } else {
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
            plain_modulus_value = parms.plain_modulus().value();            
        }
        

        SEALContext context(parms);
        print_parameters(context);
        
        // 生成密钥
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        
        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
                
        size_t d = poly_modulus_degree;
        cout << "\n=== 测试参数 ===" << endl;
        cout << "多项式模数次数: " << d << endl;
        cout << "系数模数层数: " << parms.coeff_modulus().size() << endl;
        
        // 创建测试矩阵
        cout << "\n=== 创建测试矩阵 ===" << endl;
        
        // 创建1×d的密文矩阵（每个密文编码一个矩阵行）
        vector<Ciphertext> encrypted_matrix_2;
        
        // 创建d×d的明文矩阵
        vector<vector<uint64_t>> plain_matrix_1(d, vector<uint64_t>(d));
        vector<vector<uint64_t>> plain_matrix_2(d, vector<uint64_t>(d));
        
        // 初始化明文矩阵（简单的测试模式）
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) {
                plain_matrix_1[i][j] = rand() % 10;
                plain_matrix_2[i][j] = rand() % 10;
            }
        }
        
        // 创建并加密矩阵行
        encrypt_matrix(true, context, encryptor, plain_matrix_2, encrypted_matrix_2);

        print_matrix(plain_matrix_1, "明文矩阵 1");
        print_matrix(plain_matrix_2, "明文矩阵 2");
        print_ciphertext_info(encrypted_matrix_2, "密文矩阵 2");
        
        // 执行矩阵乘法
        cout << "\n=== 执行密文-明文矩阵乘法 ===" << endl;
        vector<Ciphertext> result_matrix;
        
        auto start_time = chrono::high_resolution_clock::now();
        
        bool success = ciphertext_plaintext_matrix_multiply(
            context, plain_matrix_1, encrypted_matrix_2, result_matrix);

        bool success_2 = ciphertext_plaintext_matrix_multiply(
            context, plain_matrix_1, encrypted_matrix_2, result_matrix);
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        
        if (!success) {
            cerr << "矩阵乘法失败！" << endl;
            return 1;
        }
        
        cout << "矩阵乘法完成，耗时: " << duration.count() << " ms" << endl;
        print_ciphertext_info(result_matrix, "结果密文矩阵");
        
        // 验证结果
        cout << "\n=== 验证结果 ===" << endl;
        
        // 计算期望结果（明文计算）
        vector<vector<uint64_t>> expected_result;
        matrix_multiply_plain_blas(plain_matrix_1, plain_matrix_2, expected_result);
        
        print_matrix(expected_result, "期望结果矩阵");
        
        // 解密并验证结果
        cout << "\n解密结果验证:" << endl;
        bool all_correct = true;

        vector<vector<uint64_t>> decrypted_result_matrix;
        decrypt_ciphertexts_to_matrix(true, context, decryptor, result_matrix, decrypted_result_matrix);
        print_matrix(decrypted_result_matrix, "解密结果矩阵");
        
        check_matrix_equal(expected_result, decrypted_result_matrix);

        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

int test_ccmm() {
    try {
        // 读取配置文件
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "无法读取配置文件，使用默认参数" << endl;
            return 1;
        }
        
        // 获取用户输入的多项式模数次数
        size_t poly_modulus_degree = get_user_poly_modulus_degree(config);
        
        // 获取对应的系数模数参数
        vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
        if (coeff_modulus_params.empty()) {
            cerr << "无法获取系数模数参数" << endl;
            return 1;
        }
        
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));

        uint64_t plain_modulus_value = 1;
        if (scheme == scheme_type::ckks) {
            for (const auto &mod : parms.coeff_modulus())
                plain_modulus_value *= mod.value();
        } else {
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
            plain_modulus_value = parms.plain_modulus().value();            
        }
        

        SEALContext context(parms);
        print_parameters(context);
        
        // 生成密钥
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        RelinKeys relin_keys;
        keygen.create_relin_keys(relin_keys);
        
        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
        Evaluator evaluator(context);

/*
        Plaintext x_plain("6");
        cout << "x_plain: " << x_plain.to_string() << endl;
        Ciphertext x_encrypted;
        encryptor.encrypt(x_plain, x_encrypted);

        Plaintext plain_2("2");
        Ciphertext x_encrypted_2;
        encryptor.encrypt(plain_2, x_encrypted_2);

        Ciphertext x_encrypted_3;
        evaluator.multiply(x_encrypted, x_encrypted_2, x_encrypted_3);

        Ciphertext x_encrypted_3_my;
        x_encrypted_3_my.resize(context, 3);

        for (size_t i = 0; i < x_encrypted_3.size(); i++) {
            // 获取源密文和目标密文的数据指针
            const uint64_t* src_data = x_encrypted_3.data(i);
            uint64_t* dst_data = x_encrypted_3_my.data(i);
            
            // 获取该多项式的系数数量
            size_t coeff_count = x_encrypted_3.poly_modulus_degree() * x_encrypted_3.coeff_modulus_size();
            
            // 复制所有系数
            for (size_t j = 0; j < coeff_count; j++) {
                dst_data[j] = src_data[j];
            }
        }

        evaluator.relinearize_inplace(x_encrypted_3_my, relin_keys);
        evaluator.relinearize_inplace(x_encrypted_3, relin_keys);

        Plaintext x_plain_my;
        decryptor.decrypt(x_encrypted_3_my, x_plain_my);
        cout << "x_plain_my: " << x_plain_my.to_string() << endl;
        decryptor.decrypt(x_encrypted_3, x_plain);
        cout << "x_plain: " << x_plain.to_string() << endl;
*/

        size_t d = poly_modulus_degree;
        cout << "\n=== 测试参数 ===" << endl;
        cout << "多项式模数次数: " << d << endl;
        cout << "系数模数层数: " << parms.coeff_modulus().size() << endl;

        // 创建测试矩阵
        cout << "\n=== 创建测试矩阵 ===" << endl;
        
        // 创建密文矩阵（每个密文编码一个矩阵行/列）
        vector<Ciphertext> encrypted_matrix_1;
        vector<Ciphertext> encrypted_matrix_2;
        
        // 创建d×d的明文矩阵
        vector<vector<uint64_t>> plain_matrix_1(d, vector<uint64_t>(d));
        vector<vector<uint64_t>> plain_matrix_2(d, vector<uint64_t>(d));
        
        // 初始化明文矩阵（简单的测试模式）
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) {
                plain_matrix_1[i][j] = rand() % 10;
                plain_matrix_2[i][j] = rand() % 10;
            }
        }
        
        // 创建并加密矩阵行和列
        encrypt_matrix(true, context, encryptor, plain_matrix_1, encrypted_matrix_1);
        encrypt_matrix(true, context, encryptor, plain_matrix_2, encrypted_matrix_2);

        print_matrix(plain_matrix_1, "明文矩阵 1");
        print_matrix(plain_matrix_2, "明文矩阵 2");
        print_ciphertext_info(encrypted_matrix_1, "密文矩阵 1");
        print_ciphertext_info(encrypted_matrix_2, "密文矩阵 2");
        
        // 执行矩阵乘法
        cout << "\n=== 执行密文-密文矩阵乘法 ===" << endl;
        vector<Ciphertext> result_matrix;
        
        auto start_time = chrono::high_resolution_clock::now();
        
        bool success = ciphertext_ciphertext_matrix_multiply(
            context, relin_keys, evaluator, encrypted_matrix_1, encrypted_matrix_2, result_matrix);
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        
        if (!success) {
            cerr << "矩阵乘法失败！" << endl;
            return 1;
        }
        
        cout << "矩阵乘法完成，耗时: " << duration.count() << " ms" << endl;
        print_ciphertext_info(result_matrix, "结果密文矩阵");
        
        // 验证结果
        cout << "\n=== 验证结果 ===" << endl;
        
        // 计算期望结果（明文计算）
        vector<vector<uint64_t>> expected_result;
        matrix_multiply_plain_blas(plain_matrix_1, plain_matrix_2, expected_result, plain_modulus_value);
        
        print_matrix(expected_result, "期望结果矩阵");
        
        // 解密并验证结果
        cout << "\n解密结果验证:" << endl;
        bool all_correct = true;

        vector<vector<uint64_t>> decrypted_result_matrix;
        decrypt_ciphertexts_to_matrix(true, context, decryptor, result_matrix, decrypted_result_matrix);
        print_matrix(decrypted_result_matrix, "解密结果矩阵");
        
        check_matrix_equal(expected_result, decrypted_result_matrix);

        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

int test_CMT() {
    try {
        // 读取配置文件
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "无法读取配置文件，使用默认参数" << endl;
            return 1;
        }
        
        // 获取用户输入的多项式模数次数
        size_t poly_modulus_degree = get_user_poly_modulus_degree(config);
        
        // 获取对应的系数模数参数
        vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
        if (coeff_modulus_params.empty()) {
            cerr << "无法获取系数模数参数" << endl;
            return 1;
        }
        
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
        parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

        SEALContext context(parms);
        print_parameters(context);
        
        // 生成密钥
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        
        // 生成Galois密钥（用于automorphism操作）
        GaloisKeys galois_keys;
        
        // 为CM-T算法生成所需的Galois元素：2*t+1，其中t从0到N-1
        vector<uint32_t> galois_elts;
        for (size_t t = 0; t < poly_modulus_degree; t++) {
            galois_elts.push_back(2 * t + 1);
        }
        keygen.create_galois_keys(galois_elts, galois_keys);
        
        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
        Evaluator evaluator(context);

        size_t d = poly_modulus_degree;
        cout << "\n=== 测试参数 ===" << endl;
        cout << "多项式模数次数: " << d << endl;
        cout << "系数模数层数: " << parms.coeff_modulus().size() << endl;
        
        // 创建测试矩阵
        cout << "\n=== 创建测试矩阵 ===" << endl;
        
        // 创建d×d的明文矩阵
        vector<vector<uint64_t>> plain_matrix(d, vector<uint64_t>(d));
        
        // 初始化明文矩阵（简单的测试模式，便于验证转置）
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) {
                // 使用简单的模式：第i行第j列 = i*1000 + j
                plain_matrix[i][j] = i * 1000 + j;
            }
        }
        
        print_matrix(plain_matrix, "原始明文矩阵（前10x10）");
        
        // 按行加密矩阵
        cout << "\n=== 按行加密矩阵 ===" << endl;
        vector<Ciphertext> row_encrypted_matrix;
        encrypt_matrix(true, context, encryptor, plain_matrix, row_encrypted_matrix);
        print_ciphertext_info(row_encrypted_matrix, "按行加密的密文矩阵");
        
        // 执行密文矩阵转置
        cout << "\n=== 执行密文矩阵转置 ===" << endl;
        vector<Ciphertext> col_encrypted_matrix;
        
        auto start_time = chrono::high_resolution_clock::now();
        
        ciphertext_matrix_transpose(context, galois_keys, evaluator, row_encrypted_matrix, col_encrypted_matrix);
        // ciphertext_matrix_transpose_tweak(context, galois_keys, evaluator, row_encrypted_matrix, col_encrypted_matrix);

        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        
        cout << "密文矩阵转置完成，耗时: " << duration.count() << " ms" << endl;
        print_ciphertext_info(col_encrypted_matrix, "按列加密的密文矩阵");
        
        // 解密转置后的矩阵
        cout << "\n=== 解密转置后的矩阵 ===" << endl;
        vector<vector<uint64_t>> decrypted_transposed_matrix;
        decrypt_ciphertexts_to_matrix(false, context, decryptor, 
                                     col_encrypted_matrix, decrypted_transposed_matrix);
        print_matrix(decrypted_transposed_matrix, "解密后的转置矩阵（前10x10）");
        
        // 验证转置结果
        cout << "\n=== 验证转置结果 ===" << endl;
        
        // 计算期望的转置矩阵
        vector<vector<uint64_t>> expected_transposed_matrix(d, vector<uint64_t>(d));
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) {
                expected_transposed_matrix[i][j] = plain_matrix[j][i];
            }
        }
        
        print_matrix(expected_transposed_matrix, "期望的转置矩阵（前10x10）");
        
        // 检查结果是否正确
        bool all_correct = true;
        size_t check_size = min(size_t(10), d); // 只检查前10x10
        
        cout << "\n检查前" << check_size << "x" << check_size << "的结果:" << endl;
        for (size_t i = 0; i < check_size && all_correct; i++) {
            for (size_t j = 0; j < check_size && all_correct; j++) {
                if (decrypted_transposed_matrix[i][j] != expected_transposed_matrix[i][j]) {
                    cout << "错误: 位置[" << i << "][" << j << "] "
                         << "期望=" << expected_transposed_matrix[i][j] 
                         << ", 实际=" << decrypted_transposed_matrix[i][j] << endl;
                    all_correct = false;
                }
            }
        }
        
        if (all_correct) {
            cout << "✓ 密文矩阵转置测试成功！前" << check_size << "x" << check_size << "元素完全正确。" << endl;
        } else {
            cout << "✗ 密文矩阵转置测试失败！" << endl;
        }

        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

int test_general_multiplication(int num_bits) {
    try {
        // 读取配置文件
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "无法读取配置文件，使用默认参数" << endl;
            return 1;
        }
        // 获取用户输入的多项式模数次数
        size_t poly_modulus_degree = get_user_poly_modulus_degree(config);
        // 获取对应的系数模数参数
        vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
        if (coeff_modulus_params.empty()) {
            cerr << "无法获取系数模数参数" << endl;
            return 1;
        }
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
        uint64_t plain_modulus_value = 1;
        if (scheme == scheme_type::ckks) {
            for (const auto &mod : parms.coeff_modulus())
                plain_modulus_value *= mod.value();
        } else {
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
            plain_modulus_value = parms.plain_modulus().value();            
        }
        SEALContext context(parms);
        print_parameters(context);
        // 生成密钥
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
        Evaluator evaluator(context);
        cout << "\n=== 测试参数 ===" << endl;
        cout << "多项式模数次数: " << poly_modulus_degree << endl;
        cout << "系数模数层数: " << parms.coeff_modulus().size() << endl;
        cout << "位数: " << num_bits << endl;
        // 创建测试向量
        cout << "\n=== 创建测试向量 ===" << endl;
        vector<uint64_t> input_vector(poly_modulus_degree);
        for (size_t i = 0; i < poly_modulus_degree; i++) {
            input_vector[i] = rand() % 8;
        }
        cout << "输入向量大小: " << input_vector.size() << endl;
        cout << "前10个元素: ";
        for (size_t i = 0; i < min(size_t(10), input_vector.size()); i++) {
            cout << input_vector[i] << " ";
        }
        cout << endl;
        vector<uint64_t> test_multipliers = {3, 7, 15, 31, 63, 127, 255, 511, 1023};
        for (uint64_t multiplier : test_multipliers) {
            cout << "\n=== 测试乘以" << multiplier << " ===" << endl;
            vector<vector<uint64_t>> bit_vectors;
            decompose_to_bit_vectors(input_vector, bit_vectors, num_bits);
            vector<Ciphertext> encrypted_bit_vectors;
            encrypt_bit_vectors(context, encryptor, bit_vectors, encrypted_bit_vectors, num_bits);
            vector<Ciphertext> result_vectors;
            multiply_by_general_scalar(context, encryptor, evaluator, encrypted_bit_vectors, multiplier, result_vectors, num_bits);
            vector<vector<uint64_t>> decrypted_bit_vectors;
            decrypt_bit_vectors(context, decryptor, result_vectors, decrypted_bit_vectors, num_bits);
            vector<uint64_t> output_vector;
            compose_from_bit_vectors(decrypted_bit_vectors, output_vector, num_bits);
            bool verify_success = verify_general_multiplication(input_vector, multiplier, output_vector);
            if (!verify_success) {
                cerr << "验证失败！" << endl;
                return 1;
            }
            cout << "✓ 乘以" << multiplier << "测试成功！" << endl;
        }
        cout << "\n=== 通用向量乘法测试完成 ===" << endl;
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    return 0;
}

int test_ciphertext_scale_multiplication() {
    try {
        // 读取配置文件
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "无法读取配置文件，使用默认参数" << endl;
            return 1;
        }
        
        // 获取用户输入的多项式模数次数
        size_t poly_modulus_degree = get_user_poly_modulus_degree(config);
        
        // 获取对应的系数模数参数
        vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
        if (coeff_modulus_params.empty()) {
            cerr << "无法获取系数模数参数" << endl;
            return 1;
        }
        
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
        parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

        SEALContext context(parms);
        print_parameters(context);
        
        // 生成密钥
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        
        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
        Evaluator evaluator(context);
        
        cout << "\n=== 测试参数 ===" << endl;
        cout << "多项式模数次数: " << poly_modulus_degree << endl;
        cout << "系数模数层数: " << parms.coeff_modulus().size() << endl;
        
        // 第一步：创建测试向量并加密
        cout << "\n=== 第一步：创建测试向量并加密 ===" << endl;
        vector<uint64_t> input_vector(poly_modulus_degree);
        
         for (size_t i = 0; i < poly_modulus_degree; i++)
            input_vector[i] = rand() % 100;
        
        cout << "输入向量大小: " << input_vector.size() << endl;
        cout << "前10个元素: ";
        for (size_t i = 0; i < min(size_t(10), input_vector.size()); i++) {
            cout << input_vector[i] << " ";
        }
        cout << endl;
        
        // 编码为明文多项式
        Plaintext plaintext;
        encode_vector_to_plaintext(input_vector, context, plaintext);
        
        // 加密
        Ciphertext encrypted;
        encryptor.encrypt(plaintext, encrypted);
        
        cout << "加密完成，密文大小: " << encrypted.size() << endl;
        
        // 第二步：提取系数、乘以scale、重新打包
        cout << "\n=== 第二步：提取系数、乘以scale、重新打包 ===" << endl;
        
        // 获取用户输入的scale值
        uint64_t scale;
        cout << "请输入scale值: ";
        cin >> scale;
        cin.ignore();
        
        cout << "使用scale值: " << scale << endl;
        
        // 统计时间
        auto total_start = chrono::high_resolution_clock::now();
        
        // 选择使用哪种方法
        cout << "请选择scale乘法方法:" << endl;
        cout << "1. BLAS加速版本（提取系数->BLAS乘法->重新打包）" << endl;
        cout << "2. 直接版本（使用util::negacyclic_multiply_poly_mono_coeffmod）" << endl;
        cout << "请输入选择 (1-2): ";
        int version_choice;
        cin >> version_choice;
        cin.ignore();
        
        vector<vector<uint64_t>> coeff_matrix_a, coeff_matrix_b;
        vector<uint64_t> modulus_vector;
        chrono::microseconds extract_duration(0);
        
        if (version_choice == 1) {
            // 2.1 提取系数 - 使用cpscale.cpp中的函数
            auto extract_start = chrono::high_resolution_clock::now();
            
            double extract_time = extract_coefficients_from_single_ciphertext(
                context, encrypted, coeff_matrix_a, coeff_matrix_b, modulus_vector);
            
            auto extract_end = chrono::high_resolution_clock::now();
            extract_duration = chrono::duration_cast<chrono::microseconds>(extract_end - extract_start);
            
            cout << "系数提取完成，耗时: " << extract_duration.count() << " microseconds" << endl;
        }
        
        // 2.2 乘以scale
        auto multiply_start = chrono::high_resolution_clock::now();
        
        double multiply_time = 0;
        Ciphertext scaled_encrypted;
        
        if (version_choice == 1) {
            multiply_time = scale_coefficients_blas(coeff_matrix_a, coeff_matrix_b, modulus_vector, scale);
            
            // 重新打包为密文
            auto repack_start = chrono::high_resolution_clock::now();
            double repack_time = build_single_ciphertext_from_result(
                context, coeff_matrix_a, coeff_matrix_b, scaled_encrypted);
            auto repack_end = chrono::high_resolution_clock::now();
            auto repack_duration = chrono::duration_cast<chrono::microseconds>(repack_end - repack_start);
            cout << "重新打包完成，耗时: " << repack_duration.count() << " microseconds" << endl;
            
        } else if (version_choice == 2) {
            // 直接对密文进行scale，不需要提取和重新打包
            for (size_t i = 0; i < 100; i++) {
                scaled_encrypted = encrypted; // 复制原始密文
                MemoryPoolHandle pool = MemoryManager::GetPool();
                multiply_time += scale_ciphertext_direct(context, scaled_encrypted, scale, pool);
            }
            multiply_time /= 100;
            
        } else {
            throw std::runtime_error("Invalid version choice");
        }
        
        auto multiply_end = chrono::high_resolution_clock::now();
        auto multiply_duration = chrono::duration_cast<chrono::microseconds>(multiply_end - multiply_start);
        
        cout << "系数乘法完成，耗时: " << multiply_duration.count() << " microseconds" << endl;
        
        auto total_end = chrono::high_resolution_clock::now();
        auto total_duration = chrono::duration_cast<chrono::microseconds>(total_end - total_start);
        
        cout << "总耗时: " << total_duration.count() << " microseconds" << endl;
        cout << "详细统计:" << endl;
        if (version_choice == 2) {
            cout << "  直接scale乘法: " << multiply_duration.count() << " microseconds" << endl;
        } else {
            cout << "  系数提取: " << extract_duration.count() << " microseconds" << endl;
            cout << "  系数乘法: " << multiply_duration.count() << " microseconds" << endl;
        }
        
        // 第三步：解密并验证
        cout << "\n=== 第三步：解密并验证 ===" << endl;
        
        // 解密原始密文
        Plaintext decrypted_original;
        decryptor.decrypt(encrypted, decrypted_original);
        vector<uint64_t> original_vector;
        decode_plaintext_to_vector(decrypted_original, context, original_vector);
        
        // 解密scaled密文
        Plaintext decrypted_scaled;
        decryptor.decrypt(scaled_encrypted, decrypted_scaled);
        vector<uint64_t> scaled_vector;
        decode_plaintext_to_vector(decrypted_scaled, context, scaled_vector);
        
        cout << "原始向量前10个元素: ";
        for (size_t i = 0; i < min(size_t(10), original_vector.size()); i++) {
            cout << original_vector[i] << " ";
        }
        cout << endl;
        
        cout << "Scaled向量前10个元素: ";
        for (size_t i = 0; i < min(size_t(10), scaled_vector.size()); i++) {
            cout << scaled_vector[i] << " ";
        }
        cout << endl;
        
        // 验证结果
        bool all_correct = true;
        size_t check_size = min(size_t(10), input_vector.size());
        
        cout << "\n验证前" << check_size << "个元素:" << endl;
        for (size_t i = 0; i < check_size && all_correct; i++) {
            uint64_t expected = (input_vector[i] * scale) % parms.plain_modulus().value();
            uint64_t actual = scaled_vector[i];
            
            if (expected != actual) {
                cout << "错误: 位置[" << i << "] "
                     << "期望=" << expected 
                     << ", 实际=" << actual << endl;
                all_correct = false;
            }
        }
        
        if (all_correct) {
            cout << "✓ 密文scale乘法测试成功！前" << check_size << "个元素完全正确。" << endl;
        } else {
            cout << "✗ 密文scale乘法测试失败！" << endl;
        }
        
        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

void test_function() {
    using namespace std::chrono;
    constexpr size_t N = 4096;
    cout << "\n=== vector<vector<uint64_t>> 到 vector<double> 两种转换方式性能对比 ===" << endl;
    cout << "矩阵大小: " << N << "x" << N << endl;

    // 生成测试数据
    vector<vector<uint64_t>> A(N, vector<uint64_t>(N));
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            A[i][j] = rand() % 100;

    // 方式1：直接for循环转换
    cout << "\n方式1: 直接for循环转换为vector<double>..." << endl;
    auto t1_start = high_resolution_clock::now();
    vector<double> A_1(N * N);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            A_1[i * N + j] = static_cast<double>(A[i][j]);
    auto t1_end = high_resolution_clock::now();
    auto t1 = duration_cast<milliseconds>(t1_end - t1_start).count();
    cout << "耗时: " << t1 << " ms" << endl;

    // 方式2：先展开为一维再用构造函数
    cout << "\n方式2: 先展开为vector<uint64_t>再用构造函数转换为vector<double>..." << endl;
    auto t2_start = high_resolution_clock::now();
    
    // 步骤2a：展开为vector<uint64_t>
    auto t2a_start = high_resolution_clock::now();
    vector<uint64_t> A_flat;
    A_flat.reserve(N * N);
    for (const auto& row : A) A_flat.insert(A_flat.end(), row.begin(), row.end());
    auto t2a_end = high_resolution_clock::now();
    auto t2a = duration_cast<milliseconds>(t2a_end - t2a_start).count();
    cout << "  步骤2a (展开为vector<uint64_t>): " << t2a << " ms" << endl;
    
    // 步骤2b：用构造函数转换为vector<double>
    auto t2b_start = high_resolution_clock::now();
    vector<double> A_2(A_flat.begin(), A_flat.end());
    auto t2b_end = high_resolution_clock::now();
    auto t2b = duration_cast<milliseconds>(t2b_end - t2b_start).count();
    cout << "  步骤2b (构造函数转换为vector<double>): " << t2b << " ms" << endl;
    
    auto t2_end = high_resolution_clock::now();
    auto t2 = duration_cast<milliseconds>(t2_end - t2_start).count();
    cout << "  方式2总耗时: " << t2 << " ms" << endl;

    cout << "\n方式1/方式2的A_1[0] = " << A_1[0] << ", A_2[0] = " << A_2[0] << endl;
    cout << "\n=== 性能对比结束 ===" << endl;
}

int main() {
    string choice;
    
    while (true) {
        print_menu();
        getline(cin, choice);
        
        if (choice == "1") {
            cout << "\n开始测试 CPMM (密文-明文矩阵乘法)..." << endl;
            cout << "==========================================" << endl;
            int result = test_cpmm();
            if (result == 0) {
                cout << "\nCPMM 测试成功完成！" << endl;
            } else {
                cout << "\nCPMM 测试失败！" << endl;
            }
        }
        else if (choice == "2") {
            cout << "\n开始测试 CCMM (密文-密文矩阵乘法)..." << endl;
            cout << "==========================================" << endl;
            int result = test_ccmm();
            if (result == 0) {
                cout << "\nCCMM 测试成功完成！" << endl;
            } else {
                cout << "\nCCMM 测试失败！" << endl;
            }
        }
        else if (choice == "3") {
            cout << "\n开始测试 CMT (密文矩阵转置)..." << endl;
            cout << "==========================================" << endl;
            int result = test_CMT();
            if (result == 0) {
                cout << "\nCMT 测试成功完成！" << endl;
            } else {
                cout << "\nCMT 测试失败！" << endl;
            }
        }
        else if (choice == "4") {
            cout << "\n开始测试 BLAS浮点矩阵乘法性能..." << endl;
            cout << "==========================================" << endl;
            test_function();
        }
        else if (choice == "5") {
            cout << "\n开始自定义位数测试..." << endl;
            cout << "==========================================" << endl;
            
            int num_bits = 64;
            cout << "请输入要测试的bit宽度（如8、16、32、64等）: ";
            cin >> num_bits;
            cin.ignore();
            int result = test_general_multiplication(num_bits);
            if (result == 0)
                cout << "\n通用乘法测试成功完成！（bit宽度=" << num_bits << ")" << endl;
            else
                cout << "\n通用乘法测试失败！（bit宽度=" << num_bits << ")" << endl;
        }
        else if (choice == "6") {
            cout << "\n开始测试密文scale乘法..." << endl;
            cout << "==========================================" << endl;
            int result = test_ciphertext_scale_multiplication();
            if (result == 0)
                cout << "\n密文scale乘法测试成功完成！" << endl;
            else
                cout << "\n密文scale乘法测试失败！" << endl;
        }
        else if (choice == "7") {
            cout << "程序退出。" << endl;
            break;
        }
        else {
            cout << "无效选择，请输入 1、2、3、4、5、6、7。" << endl;
        }
        
        cout << "\n按回车键继续...";
        getline(cin, choice);
    }
    
    return 0;
}