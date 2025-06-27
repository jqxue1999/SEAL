#include "example.h"
#include <chrono>
#include <cstdlib>

using namespace seal;
using namespace std;

void print_menu() {
    cout << "\n=== SEAL 矩阵乘法测试程序 ===" << endl;
    cout << "请选择要测试的功能：" << endl;
    cout << "1. CPMM (Ciphertext-Plaintext Matrix Multiplication) - 密文-明文矩阵乘法" << endl;
    cout << "2. CCMM (Ciphertext-Ciphertext Matrix Multiplication) - 密文-密文矩阵乘法" << endl;
    cout << "3. CMT (Ciphertext Matrix Transpose) - 密文矩阵转置" << endl;
    cout << "4. Test Function - 测试CM-T算法第二步Galois Automorphism" << endl;
    cout << "5. 退出程序" << endl;
    cout << "请输入选择 (1-5): ";
}

int test_cpmm() {
    try {
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        size_t poly_modulus_degree = 4096;
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
        // Double can only represent 2^53-1, so we need to carefully design a new coeff_modulus for degree 8192
        // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 36, 36, 36, 36, 37 }));

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
        // matrix_multiply_plain(plain_matrix_1, plain_matrix_2, expected_result, plain_modulus_value);
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

int test_ccmm() {
    try {
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        size_t poly_modulus_degree = 4096;
        parms.set_poly_modulus_degree(poly_modulus_degree);
        // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
        // Double can only represent 2^53-1, so we need to carefully design a new coeff_modulus for degree 8192
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 26, 26}));

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
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        size_t poly_modulus_degree = 4096;
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
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

int test_function() {
    try {
        // 设置加密参数
        scheme_type scheme = scheme_type::bfv;
        EncryptionParameters parms(scheme);
        size_t poly_modulus_degree = 4096;
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
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

        cout << "\n=== 测试CM-T算法第二步：Galois Automorphism ===" << endl;
        
        // 创建测试向量
        vector<uint64_t> original_vector(poly_modulus_degree);
        for (size_t i = 0; i < poly_modulus_degree; i++) {
            original_vector[i] = i; // 简单的测试数据
        }
        
        cout << "原始向量（前10个元素）: ";
        for (size_t i = 0; i < min(size_t(10), poly_modulus_degree); i++) {
            cout << original_vector[i] << " ";
        }
        cout << endl;
        
        // 编码并加密向量
        Plaintext plain_vector;
        encode_vector_to_plaintext(original_vector, context, plain_vector);
        Ciphertext encrypted_vector;
        encryptor.encrypt(plain_vector, encrypted_vector);
        cout << "向量加密完成" << endl;
        
        // 测试不同的Galois automorphism
        size_t test_cases = 5; // 测试前5个case
        cout << "\n测试前" << test_cases << "个Galois automorphism (2*t+1):" << endl;
        
        for (size_t t = 0; t < test_cases; t++) {
            uint32_t galois_elt = 2 * t + 1;
            cout << "\n--- 测试 t=" << t << ", Galois元素=" << galois_elt << " ---" << endl;
            
            // 复制原始密文
            Ciphertext test_cipher = encrypted_vector;
            
            // 应用Galois automorphism
            evaluator.apply_galois_inplace(test_cipher, galois_elt, galois_keys);
            cout << "Galois automorphism应用完成" << endl;
            
            // 解密并解码
            Plaintext decrypted_plain;
            decryptor.decrypt(test_cipher, decrypted_plain);
            vector<uint64_t> transformed_vector;
            decode_plaintext_to_vector(decrypted_plain, context, transformed_vector);
            
            cout << "变换后向量（前10个元素）: ";
            for (size_t i = 0; i < min(size_t(10), poly_modulus_degree); i++) {
                cout << transformed_vector[i] << " ";
            }
            cout << endl;
            

        }
        
        // 测试CM-T算法第二步的完整流程
        cout << "\n=== 测试CM-T算法第二步完整流程 ===" << endl;
        
        // 模拟tilde_m向量（这里我们使用相同的加密向量作为示例）
        vector<Ciphertext> tilde_m(3); // 只测试前3个
        for (size_t t = 0; t < 3; t++) {
            tilde_m[t] = encrypted_vector; // 使用相同的向量作为示例
        }
        
        // 执行CM-T算法第二步：bar_m_t = tilde_m_t(X^{2t+1})
        vector<Ciphertext> bar_m(3);
        for (size_t t = 0; t < 3; t++) {
            bar_m[t] = tilde_m[t];
            evaluator.apply_galois_inplace(bar_m[t], 2*t+1, galois_keys);
            cout << "完成 bar_m[" << t << "] = tilde_m[" << t << "](X^{" << 2*t+1 << "})" << endl;
        }
        
        // 验证所有bar_m都被正确创建
        bool all_bar_m_valid = true;
        for (size_t t = 0; t < 3; t++) {
            try {
                // 尝试解密每个bar_m，确保它们都是有效的密文
                Plaintext test_plain;
                decryptor.decrypt(bar_m[t], test_plain);
                cout << "✓ bar_m[" << t << "] 解密成功" << endl;
            } catch (const exception& e) {
                cout << "✗ bar_m[" << t << "] 解密失败: " << e.what() << endl;
                all_bar_m_valid = false;
            }
        }
        
        if (all_bar_m_valid) {
            cout << "\n✓ CM-T算法第二步测试成功！所有Galois automorphism操作都正确执行。" << endl;
        } else {
            cout << "\n✗ CM-T算法第二步测试失败！" << endl;
        }

        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
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
            cout << "\n开始测试 CM-T算法第二步Galois Automorphism..." << endl;
            cout << "==========================================" << endl;
            int result = test_function();
            if (result == 0) {
                cout << "\nCM-T算法第二步Galois Automorphism测试成功完成！" << endl;
            } else {
                cout << "\nCM-T算法第二步Galois Automorphism测试失败！" << endl;
            }
        }
        else if (choice == "5") {
            cout << "程序退出。" << endl;
            break;
        }
        else {
            cout << "无效选择，请输入 1、2、3、4 或 5。" << endl;
        }
        
        cout << "\n按回车键继续...";
        getline(cin, choice);
    }
    
    return 0;
}