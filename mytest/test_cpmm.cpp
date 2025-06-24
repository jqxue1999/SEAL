#include "cpmm.h"
#include "utils.h"


using namespace seal;
using namespace std;



int main() {
    try {
        // 设置加密参数
        EncryptionParameters parms(scheme_type::bfv);
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
        encrypt_matrix_rows(context, encryptor, plain_matrix_2, encrypted_matrix_2);

        print_matrix(plain_matrix_1, "明文矩阵 1");
        print_matrix(plain_matrix_2, "明文矩阵 2");
        print_ciphertext_info(encrypted_matrix_2, "密文矩阵 2");
        
        vector<vector<uint64_t>> plain_matrix_decrypted_2;
        decrypt_ciphertexts_to_matrix(context, decryptor, encrypted_matrix_2, plain_matrix_decrypted_2);
        print_matrix(plain_matrix_decrypted_2, "解密后的明文矩阵 2");
        check_matrix_equal(plain_matrix_2, plain_matrix_decrypted_2);
        
        
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
        
        // 获取明文模数
        uint64_t plain_modulus_value = parms.plain_modulus().value();
        cout << "明文模数: " << plain_modulus_value << endl;
        
        // 计算期望结果（明文计算）
        vector<vector<uint64_t>> expected_result;
        matrix_multiply_plain(plain_matrix_1, plain_matrix_2, expected_result, plain_modulus_value);
        
        print_matrix(expected_result, "期望结果矩阵");
        
        // 解密并验证结果
        cout << "\n解密结果验证:" << endl;
        bool all_correct = true;

        vector<vector<uint64_t>> decrypted_result_matrix;
        decrypt_ciphertexts_to_matrix(context, decryptor, result_matrix, decrypted_result_matrix);
        print_matrix(decrypted_result_matrix, "解密结果矩阵");
        
        check_matrix_equal(expected_result, decrypted_result_matrix);

        cout << "\n=== 测试完成 ===" << endl;
        
    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }
    
    return 0;
} 