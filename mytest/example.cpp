#include "example.h"
#include <chrono>
#include <cstdlib>
#include <cblas.h>
#include <unistd.h>
#include "cpscale.h"

using namespace seal;
using namespace std;


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



void print_menu() {
    cout << "\n=== SEAL 矩阵乘法测试程序 ===" << endl;
    cout << "请选择要测试的功能：" << endl;
    cout << "1. CPMM (Ciphertext-Plaintext Matrix Multiplication) - 密文-明文矩阵乘法" << endl;
    cout << "2. CCMM (Ciphertext-Ciphertext Matrix Multiplication) - 密文-密文矩阵乘法" << endl;
    cout << "3. CMT (Ciphertext Matrix Transpose) - 密文矩阵转置" << endl;
    cout << "10. 退出程序" << endl;
    cout << "请输入选择 (1-10): ";
}

void test_cpmm(bool verbose) {
    try {
        json config = read_seal_config("/home/jiaq/Research/SEAL/mytest/seal_config.json");
        if (config.empty()) {
            throw std::runtime_error("Error: cannot read config file");
        }
        
        vector<size_t> poly_modulus_degree_options = {4096, 8192, 16384};
        for (size_t poly_modulus_degree : poly_modulus_degree_options) {
            cout << "degree: " << poly_modulus_degree << endl;
            vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
            if (coeff_modulus_params.empty()) {
                throw std::runtime_error("Error: cannot get coeff modulus params");
            }
            
            scheme_type scheme = scheme_type::bfv;
            EncryptionParameters parms(scheme);
            parms.set_poly_modulus_degree(poly_modulus_degree);
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

            SEALContext context(parms);
            if (verbose)
                print_parameters(context);
            
            // 生成密钥
            KeyGenerator keygen(context);
            SecretKey secret_key = keygen.secret_key();
            PublicKey public_key;
            keygen.create_public_key(public_key);
            
            Encryptor encryptor(context, public_key);
            Decryptor decryptor(context, secret_key);
            
            size_t d = poly_modulus_degree;
            vector<double> time_vec{0, 0, 0};
            
            int num_runs = 10;
            for (size_t run = 0; run < num_runs; run++) {
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
                vector<Ciphertext> encrypted_matrix_2;
                encrypt_matrix(true, context, encryptor, plain_matrix_2, encrypted_matrix_2);
                
                // 执行矩阵乘法
                vector<Ciphertext> result_matrix;
                vector<double> time_vec_tmp = ciphertext_plaintext_matrix_multiply(
                    context, plain_matrix_1, encrypted_matrix_2, result_matrix, verbose);
                
                time_vec[0] += time_vec_tmp[0];
                time_vec[1] += time_vec_tmp[1];
                time_vec[2] += time_vec_tmp[2];
                
                // 验证结果
                vector<vector<uint64_t>> expected_result;
                matrix_multiply_plain_blas(plain_matrix_1, plain_matrix_2, expected_result);
                
                vector<vector<uint64_t>> decrypted_result_matrix;
                decrypt_ciphertexts_to_matrix(true, context, decryptor, result_matrix, decrypted_result_matrix);
                
                // 检查结果是否正确
                bool all_correct = true;
                size_t check_size = min(size_t(10), d); // 只检查前10x10
                
                for (size_t i = 0; i < check_size && all_correct; i++) {
                    for (size_t j = 0; j < check_size && all_correct; j++) {
                        if (decrypted_result_matrix[i][j] != expected_result[i][j]) {
                            cout << "Error: position[" << i << "][" << j << "] "
                                 << "expected=" << expected_result[i][j] 
                                 << ", actual=" << decrypted_result_matrix[i][j] << endl;
                            all_correct = false;
                        }
                    }
                }
                
                if (!all_correct) {
                    throw runtime_error("CPMM verification failed!");
                }
            }
            
            cout << "\033[33mCoefficients extraction time: " << time_vec[0] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mMatrix multiplication time: " << time_vec[1] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mCiphertext repackaging time: " << time_vec[2] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[32mTotal time: " << (time_vec[0] + time_vec[1] + time_vec[2]) * 1000 / num_runs << " ms\033[0m" << endl;
        }
    } catch (const exception& e) {
        throw std::runtime_error("Error: " + string(e.what()));
    }
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


int test_ciphertext_scale_multiplication(bool verbose) {
    try {
        json config = read_seal_config();
        if (config.empty()) {
            cerr << "Error: cannot read config file" << endl;
            return 1;
        }
        
        vector<size_t> poly_modulus_degree_options = {16384, 8192, 4096};
        for (size_t poly_modulus_degree : poly_modulus_degree_options) {
            cout << "degree: " << poly_modulus_degree << endl;
            vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
            if (coeff_modulus_params.empty()) {
                cerr << "Error: cannot get coeff modulus params" << endl;
                return 1;
            }
            
            scheme_type scheme = scheme_type::bfv;
            EncryptionParameters parms(scheme);
            parms.set_poly_modulus_degree(poly_modulus_degree);
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

            SEALContext context(parms);
            if (verbose)
                print_parameters(context);
            
            // 生成密钥
            KeyGenerator keygen(context);
            SecretKey secret_key = keygen.secret_key();
            PublicKey public_key;
            keygen.create_public_key(public_key);
            
            Encryptor encryptor(context, public_key);
            Decryptor decryptor(context, secret_key);
            Evaluator evaluator(context);

            vector<uint64_t> input_vector(poly_modulus_degree);

            vector<double> time_vec{0, 0, 0, 0};
            
            int num_runs = 1000;
            for (size_t run = 0; run < num_runs; run++)
            {
                for (size_t i = 0; i < poly_modulus_degree; i++)
                    input_vector[i] = rand() % 1000;

                uint64_t scale = rand() % 1000;
                
                Plaintext plaintext;
                encode_vector_to_plaintext(input_vector, context, plaintext);
                Ciphertext encrypted;
                encryptor.encrypt(plaintext, encrypted);

                Ciphertext scaled_encrypted;
                vector<double> time_vec_tmp = scale_vector_blas(context, encrypted, scaled_encrypted, scale, verbose);
                time_vec[0] += time_vec_tmp[0];
                time_vec[1] += time_vec_tmp[1];
                time_vec[2] += time_vec_tmp[2];
                time_vec[3] += time_vec_tmp[3];
                Plaintext scaled_decrypted;
                decryptor.decrypt(scaled_encrypted, scaled_decrypted);
                vector<uint64_t> scaled_vector;
                decode_plaintext_to_vector(scaled_decrypted, context, scaled_vector);

                for (size_t i = 0; i < input_vector.size(); i++) {
                    uint64_t expected = (input_vector[i] * scale) % parms.plain_modulus().value();
                    uint64_t actual = scaled_vector[i];
                    
                    if (expected != actual) {
                        throw runtime_error("Verification failed!");
                    }
                }
            }
            cout << "\033[33mCoefficients extraction time: " << time_vec[0] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mCoefficients multiplication time: " << time_vec[1] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mCiphertext repackaging time: " << time_vec[2] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[32mTotal time: " << time_vec[3] * 1000 / num_runs << " ms\033[0m" << endl;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

int test_vector_vector_outer_multiplication(bool verbose) {
    try {
        json config = read_seal_config("seal_config.json", verbose);
        if (config.empty()) {
            cerr << "Error: cannot read config file" << endl;
            return 1;
        }
        
        vector<size_t> poly_modulus_degree_options = {4096, 8192, 16384};
        for (size_t poly_modulus_degree : poly_modulus_degree_options) {
            cout << "degree: " << poly_modulus_degree << endl;
            vector<int> coeff_modulus_params = get_coeff_modulus_params(config, poly_modulus_degree);
            if (coeff_modulus_params.empty()) {
                cerr << "Error: cannot get coeff modulus params" << endl;
                return 1;
            }
            
            scheme_type scheme = scheme_type::bfv;
            EncryptionParameters parms(scheme);
            parms.set_poly_modulus_degree(poly_modulus_degree);
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus_params));
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

            SEALContext context(parms);
            if (verbose)
                print_parameters(context);
            
            // 生成密钥
            KeyGenerator keygen(context);
            SecretKey secret_key = keygen.secret_key();
            PublicKey public_key;
            keygen.create_public_key(public_key);
            
            Encryptor encryptor(context, public_key);
            Decryptor decryptor(context, secret_key);
            Evaluator evaluator(context);

            vector<uint64_t> input_vector(poly_modulus_degree);
            vector<double> time_vec{0, 0, 0, 0};
            
            int num_runs = 1;
            int plain_vector_size = 4096;
            for (size_t run = 0; run < num_runs; run++) {
                for (size_t i = 0; i < poly_modulus_degree; i++) {
                    input_vector[i] = rand() % 1000;
                }
                
                vector<uint64_t> plain_vector(plain_vector_size);
                for (size_t i = 0; i < plain_vector_size; i++) {
                    plain_vector[i] = rand() % 1000;
                }
                
                Plaintext plaintext;
                encode_vector_to_plaintext(input_vector, context, plaintext);
                Ciphertext encrypted;
                encryptor.encrypt(plaintext, encrypted);
                
                vector<Ciphertext> result_vector;
                vector<double> time_vec_tmp = vector_vector_outer_multiply_blas(
                    context, plain_vector, encrypted, result_vector, verbose);
                
                time_vec[0] += time_vec_tmp[0];
                time_vec[1] += time_vec_tmp[1];
                time_vec[2] += time_vec_tmp[2];
                time_vec[3] += time_vec_tmp[3];
                
                vector<vector<uint64_t>> expected_outer_product(plain_vector_size, vector<uint64_t>(poly_modulus_degree));
                for (size_t i = 0; i < plain_vector_size; i++) {
                    for (size_t j = 0; j < poly_modulus_degree; j++) {
                        expected_outer_product[i][j] = (input_vector[j] * plain_vector[i]) % parms.plain_modulus().value();
                    }
                }
                
                for (size_t i = 0; i < plain_vector_size; i++) {
                    Plaintext decrypted_plaintext;
                    decryptor.decrypt(result_vector[i], decrypted_plaintext);
                    vector<uint64_t> decrypted_vector;
                    decode_plaintext_to_vector(decrypted_plaintext, context, decrypted_vector);
                    
                    for (size_t j = 0; j < poly_modulus_degree; j++) {
                        uint64_t expected = expected_outer_product[i][j];
                        uint64_t actual = decrypted_vector[j];
                        
                        if (expected != actual) {
                            throw runtime_error("Vector-vector outer multiplication verification failed!");
                        }
                    }
                }
            }
            
            cout << "\033[33mCoefficients extraction time: " << time_vec[0] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mCoefficients multiplication time: " << time_vec[1] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[33mCiphertext repackaging time: " << time_vec[2] * 1000 / num_runs << " ms\033[0m" << endl;
            cout << "\033[32mTotal time: " << time_vec[3] * 1000 / num_runs << " ms\033[0m" << endl;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

int main() {
    string choice;
    
    while (true) {
        print_menu();
        getline(cin, choice);
        
        if (choice == "1") {
            cout << "\n开始测试 CPMM (密文-明文矩阵乘法)..." << endl;
            cout << "==========================================" << endl;
            test_cpmm(false);
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
            cout << "\n开始测试密文scale乘法..." << endl;
            cout << "==========================================" << endl;
            test_ciphertext_scale_multiplication(false);
        }
        else if (choice == "5") {
            cout << "\n开始测试明文vector外乘密文vector..." << endl;
            cout << "==========================================" << endl;
            test_vector_vector_outer_multiplication(false);
        }
        else if (choice == "6") {     
            cout << "\ntest_digits_cvps..." << endl;
            cout << "==========================================" << endl;       
            int num_bits;
            bool check_correctness;
            cout << "Input num_bits: ";
            cin >> num_bits;
            cout << "Input check_correctness: ";
            cin >> check_correctness;
            cin.ignore();            
            test_digits_cvps(num_bits, check_correctness);
        }
        else if (choice == "7") {
            cout << "\ntest_digits_cvpv..." << endl;
            cout << "==========================================" << endl;
            
            int num_bits;
            cout << "Input num_bits: ";
            cin >> num_bits;

            bool check_correctness;
            cout << "Input check_correctness: ";
            cin >> check_correctness;
            cin.ignore();
            
            int result = test_digits_cvpv(num_bits, check_correctness);
            if (result == 0)
                cout << "\n通用乘法测试成功完成！（bit宽度=" << num_bits << ")" << endl;
            else
                cout << "\n通用乘法测试失败！（bit宽度=" << num_bits << ")" << endl;
        }
        else if (choice == "8") {
            cout << "\ntest_digits_pvcm..." << endl;
            cout << "==========================================" << endl;
            int num_bits;
            cout << "Input num_bits: ";
            cin >> num_bits;
            bool check_correctness;
            cout << "Input check_correctness: ";
            cin >> check_correctness;
            cin.ignore();
            int result = test_digits_pvcm(num_bits, check_correctness);
            if (result == 0)
                cout << "\n明文向量与密文矩阵乘法测试成功完成！（bit宽度=" << num_bits << ")" << endl;
            else
                cout << "\n明文向量与密文矩阵乘法测试失败！（bit宽度=" << num_bits << ")" << endl;
        }
        else if (choice == "9") {
            cout << "\ntest_digits_pmcm..." << endl;
            cout << "==========================================" << endl;
            int num_bits;
            cout << "Input num_bits: ";
            cin >> num_bits;
            bool check_correctness;
            cout << "Input check_correctness: ";
            cin >> check_correctness;
            cin.ignore();
            int result = test_digits_pmcm(num_bits, check_correctness);
            if (result == 0)
                cout << "\n密文矩阵与明文矩阵乘法测试成功完成！（bit宽度=" << num_bits << ")" << endl;
            else
                cout << "\n密文矩阵与明文矩阵乘法测试失败！（bit宽度=" << num_bits << ")" << endl;
        }
        else if (choice == "10") {
            cout << "程序退出。" << endl;
            break;
        }
        
        else {
            cout << "无效选择，请输入 1、2、3、4、5、6、7、8、9、10。" << endl;
        }
        
        cout << "\n按回车键继续...";
        getline(cin, choice);
    }
    
    return 0;
}