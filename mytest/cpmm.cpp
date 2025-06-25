#include "cpmm.h"
#include <cblas.h>
#include <vector>

using namespace seal;
using namespace std;

bool ciphertext_plaintext_matrix_multiply(
    const SEALContext& context,
        const vector<vector<uint64_t>>& plain_matrix,
    const vector<Ciphertext>& encrypted_matrix,
    vector<Ciphertext>& result_matrix)
{
    try {
        cout << "\n=== 开始密文-明文矩阵乘法 ===" << endl;
        
        // 获取参数
        auto context_data = context.get_context_data(encrypted_matrix[0].parms_id());
        if (!context_data) {
            cerr << "Invalid context data" << endl;
            return false;
        }
        
        size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
        size_t coeff_modulus_size = encrypted_matrix[0].coeff_modulus_size();
        
        // 检查矩阵维度
        if (encrypted_matrix.empty() || plain_matrix.empty() || plain_matrix[0].empty()) {
            cerr << "Empty matrix" << endl;
            return false;
        }
        
        size_t encrypted_cols = encrypted_matrix.size(); // 应该是 d
        size_t plain_rows = plain_matrix.size();
        size_t plain_cols = plain_matrix[0].size();
        
        // 检查矩阵大小
        if (encrypted_cols != poly_modulus_degree) {
            cerr << "Encrypted matrix size must be " << poly_modulus_degree 
                 << ", but got " << encrypted_cols << endl;
            return false;
        }
        
        if (plain_rows != poly_modulus_degree || plain_cols != poly_modulus_degree) {
            cerr << "Plain matrix size must be " << poly_modulus_degree << "x" << poly_modulus_degree 
                 << ", but got " << plain_rows << "x" << plain_cols << endl;
            return false;
        }
        
        cout << "矩阵大小验证通过:" << endl;
        cout << "  密文矩阵: 1x" << encrypted_cols << endl;
        cout << "  明文矩阵: " << plain_rows << "x" << plain_cols << endl;
        cout << "  预期结果: 1x" << plain_cols << endl;
        cout << "  多项式模数次数: " << poly_modulus_degree << endl;
        cout << "  系数模数层数: " << coeff_modulus_size << endl;
        
        // 初始化结果矩阵
        result_matrix.resize(plain_cols);
        
        cout << "\n步骤1/3: 提取系数..." << endl;
        // 步骤0: 从ciphertext矩阵每一行的ciphertext提取系数，组成矩阵
        vector<vector<vector<uint64_t>>> coeff_matrix_a, coeff_matrix_b;
        vector<uint64_t> modulus_vector;
        extract_coefficients_from_ciphertext_vector(context, encrypted_matrix, coeff_matrix_a, coeff_matrix_b, modulus_vector);
        
        cout << "\n步骤2/3: 执行矩阵乘法..." << endl;
        // 步骤1: 对每个RNS片，计算plaintext matrix * a和plaintext matrix * b
        vector<vector<vector<uint64_t>>> result_a_matrix, result_b_matrix;
        // matrix_multiply_for_all_rns_layers(plain_matrix, coeff_matrix_a, coeff_matrix_b, result_a_matrix, result_b_matrix, modulus_vector);
        matrix_multiply_for_all_rns_layers_blas(plain_matrix, coeff_matrix_a, coeff_matrix_b, result_a_matrix, result_b_matrix, modulus_vector);

        cout << "\n步骤3/3: 构建密文..." << endl;
        // 步骤2: 构建d个ciphertext，按行构建
        build_ciphertexts_from_result_matrices(context, result_a_matrix, result_b_matrix, result_matrix);
        
        cout << "\n=== 密文-明文矩阵乘法完成 ===" << endl;
        
        return true;
    }
    catch (const exception& e) {
        cerr << "Error in ciphertext_plaintext_matrix_multiply: " << e.what() << endl;
        return false;
    }
}

void extract_coefficients_from_ciphertext_vector(
    const SEALContext& context,
    const vector<Ciphertext>& encrypted_vector,
    vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector)
{
    if (encrypted_vector.empty()) {
        cerr << "Empty encrypted vector" << endl;
        return;
    }
    
    auto context_data = context.get_context_data(encrypted_vector[0].parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = encrypted_vector[0].coeff_modulus_size();
    size_t num_ciphertexts = encrypted_vector.size();
    
    // 检查向量大小
    if (num_ciphertexts != poly_modulus_degree) {
        cerr << "Encrypted vector size must be " << poly_modulus_degree 
             << ", but got " << num_ciphertexts << endl;
        return;
    }
    
    cout << "开始提取系数..." << endl;
    cout << "密文数量: " << num_ciphertexts << ", RNS层数: " << coeff_modulus_size << endl;
    
    // 获取modulus vector
    auto coeff_modulus = context_data->parms().coeff_modulus();
    modulus_vector.resize(coeff_modulus_size);
    for (size_t i = 0; i < coeff_modulus_size; i++) {
        modulus_vector[i] = coeff_modulus[i].value();
    }
    // 打印modulus_vector
    cout << "modulus_vector: ";
    for (size_t i = 0; i < coeff_modulus_size; i++) {
        cout << modulus_vector[i] << " ";
    }
    cout << endl;
    
    // 初始化系数矩阵 - 第一维是RNS层，第二维是密文索引，第三维是该密文在该RNS层上的系数
    coeff_matrix_a.resize(coeff_modulus_size);
    coeff_matrix_b.resize(coeff_modulus_size);
    
    for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
        coeff_matrix_a[rns_layer].resize(num_ciphertexts);
        coeff_matrix_b[rns_layer].resize(num_ciphertexts);
        
        for (size_t cipher_idx = 0; cipher_idx < num_ciphertexts; cipher_idx++) {
            coeff_matrix_a[rns_layer][cipher_idx].resize(poly_modulus_degree);
            coeff_matrix_b[rns_layer][cipher_idx].resize(poly_modulus_degree);
            
            const Ciphertext& cipher = encrypted_vector[cipher_idx];
            
            // 提取a多项式的系数（第一个多项式）
            if (cipher.size() > 0) {
                const uint64_t* coeffs_a = cipher.data(0);
                // 提取第rns_layer层的系数
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeff_matrix_a[rns_layer][cipher_idx][coeff_idx] = 
                        coeffs_a[coeff_idx + rns_layer * poly_modulus_degree];
                }
            }
            
            // 提取b多项式的系数（第二个多项式）
            if (cipher.size() > 1) {
                const uint64_t* coeffs_b = cipher.data(1);
                // 提取第rns_layer层的系数
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeff_matrix_b[rns_layer][cipher_idx][coeff_idx] = 
                        coeffs_b[coeff_idx + rns_layer * poly_modulus_degree];
                }
            }
        }
        
        // 显示进度
        cout << "\r提取系数进度: RNS层 " << (rns_layer + 1) << "/" << coeff_modulus_size 
             << " [" << (rns_layer + 1) * 100 / coeff_modulus_size << "%]" << flush;
    }
    
    cout << "\n系数提取完成！" << endl;
}

void matrix_multiply_for_all_rns_layers(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<vector<vector<uint64_t>>>& result_a_matrix,
    vector<vector<vector<uint64_t>>>& result_b_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (plain_matrix.empty() || coeff_matrix_a.empty() || coeff_matrix_b.empty()) {
        cerr << "Empty matrices in multiplication" << endl;
        return;
    }
    
    size_t plain_rows = plain_matrix.size();
    size_t plain_cols = plain_matrix[0].size();
    size_t rns_layers = coeff_matrix_a.size();
    
    if (coeff_matrix_b.size() != rns_layers) {
        cerr << "RNS layer count mismatch" << endl;
        return;
    }
    
    if (modulus_vector.size() != rns_layers) {
        cerr << "Modulus vector size mismatch" << endl;
        return;
    }
    
    cout << "开始矩阵乘法计算..." << endl;
    cout << "RNS层数: " << rns_layers << ", 矩阵大小: " << plain_rows << "x" << plain_cols << endl;
    
    // 初始化结果矩阵 - 每个RNS层一个d×d矩阵
    result_a_matrix.resize(rns_layers);
    result_b_matrix.resize(rns_layers);
    
    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        // 显示RNS层进度
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;
        
        result_a_matrix[rns_layer].resize(plain_rows);
        result_b_matrix[rns_layer].resize(plain_rows);
        
        for (size_t i = 0; i < plain_rows; i++) {
            result_a_matrix[rns_layer][i].resize(plain_cols, 0);
            result_b_matrix[rns_layer][i].resize(plain_cols, 0);
        }
        
        uint64_t modulus = modulus_vector[rns_layer];
        
        // 对当前RNS层执行矩阵乘法：plain_matrix * coeff_matrix_a[rns_layer] 和 plain_matrix * coeff_matrix_b[rns_layer]
        for (size_t i = 0; i < plain_rows; i++) {
            // 显示当前处理的行
            cout << "\rRNS层 " << (rns_layer + 1) << "/" << rns_layers 
                 << " - 处理第 " << (i + 1) << "/" << plain_rows << " 行" << flush;
            
            for (size_t j = 0; j < plain_cols; j++) {
                uint64_t sum_a = 0;
                uint64_t sum_b = 0;
                
                for (size_t k = 0; k < plain_cols; k++) {
                    // 对明文矩阵的值也需要进行模运算
                    uint64_t plain_val = plain_matrix[i][k];
                    uint64_t coeff_a_val = coeff_matrix_a[rns_layer][k][j];
                    uint64_t coeff_b_val = coeff_matrix_b[rns_layer][k][j];
                    
                    // 计算乘积并累加，使用模运算避免溢出
                    uint64_t product_a = plain_val * coeff_a_val;
                    uint64_t product_b = plain_val * coeff_b_val;
                    
                    sum_a = (sum_a + product_a) % modulus;
                    sum_b = (sum_b + product_b) % modulus;
                }
                
                result_a_matrix[rns_layer][i][j] = sum_a;
                result_b_matrix[rns_layer][i][j] = sum_b;
            }
        }
    }
    
    cout << "\n矩阵乘法计算完成！" << endl;
}

void matrix_multiply_for_all_rns_layers_blas(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<vector<vector<uint64_t>>>& result_a_matrix,
    vector<vector<vector<uint64_t>>>& result_b_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (plain_matrix.empty() || coeff_matrix_a.empty() || coeff_matrix_b.empty()) {
        cerr << "Empty matrices in multiplication" << endl;
        return;
    }

    size_t plain_rows = plain_matrix.size();
    size_t plain_cols = plain_matrix[0].size();
    size_t rns_layers = coeff_matrix_a.size();

    if (coeff_matrix_b.size() != rns_layers || modulus_vector.size() != rns_layers) {
        cerr << "RNS layer count or modulus vector size mismatch" << endl;
        return;
    }

    cout << "开始使用BLAS进行矩阵乘法计算..." << endl;
    cout << "RNS层数: " << rns_layers << ", 矩阵大小: " << plain_rows << "x" << plain_cols << endl;

    result_a_matrix.resize(rns_layers);
    result_b_matrix.resize(rns_layers);

    std::vector<double> plain_matrix_double(plain_rows * plain_cols);
    for(size_t i = 0; i < plain_rows; ++i) {
        for(size_t j = 0; j < plain_cols; ++j) {
            plain_matrix_double[i * plain_cols + j] = static_cast<double>(plain_matrix[i][j]);
        }
    }

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        result_a_matrix[rns_layer].resize(plain_rows, vector<uint64_t>(plain_cols));
        result_b_matrix[rns_layer].resize(plain_rows, vector<uint64_t>(plain_cols));

        size_t matrix_a_rows = coeff_matrix_a[rns_layer].size();
        size_t matrix_a_cols = coeff_matrix_a[rns_layer][0].size();
        std::vector<double> coeff_a_double(matrix_a_rows * matrix_a_cols);
        for(size_t i = 0; i < matrix_a_rows; ++i) {
            for (size_t j = 0; j < matrix_a_cols; ++j) {
                coeff_a_double[i * matrix_a_cols + j] = static_cast<double>(coeff_matrix_a[rns_layer][i][j]);
            }
        }
        
        size_t matrix_b_rows = coeff_matrix_b[rns_layer].size();
        size_t matrix_b_cols = coeff_matrix_b[rns_layer][0].size();
        std::vector<double> coeff_b_double(matrix_b_rows * matrix_b_cols);
        for(size_t i = 0; i < matrix_b_rows; ++i) {
            for (size_t j = 0; j < matrix_b_cols; ++j) {
                coeff_b_double[i * matrix_b_cols + j] = static_cast<double>(coeff_matrix_b[rns_layer][i][j]);
            }
        }

        std::vector<double> result_a_double(plain_rows * matrix_a_cols);
        std::vector<double> result_b_double(plain_rows * matrix_b_cols);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    plain_rows, matrix_a_cols, plain_cols, 1.0,
                    plain_matrix_double.data(), plain_cols,
                    coeff_a_double.data(), matrix_a_cols, 0.0,
                    result_a_double.data(), matrix_a_cols);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    plain_rows, matrix_b_cols, plain_cols, 1.0,
                    plain_matrix_double.data(), plain_cols,
                    coeff_b_double.data(), matrix_b_cols, 0.0,
                    result_b_double.data(), matrix_b_cols);
        
        uint64_t modulus = modulus_vector[rns_layer];
        for(size_t i = 0; i < plain_rows; ++i) {
            for(size_t j = 0; j < matrix_a_cols; ++j) {
                result_a_matrix[rns_layer][i][j] = static_cast<uint64_t>(fmod(result_a_double[i * matrix_a_cols + j], modulus));
            }
        }
        for(size_t i = 0; i < plain_rows; ++i) {
            for(size_t j = 0; j < matrix_b_cols; ++j) {
                result_b_matrix[rns_layer][i][j] = static_cast<uint64_t>(fmod(result_b_double[i * matrix_b_cols + j], modulus));
            }
        }
    }

    cout << "\nBLAS矩阵乘法计算完成！" << endl;
}

void build_ciphertexts_from_result_matrices(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_a_matrix,
    const vector<vector<vector<uint64_t>>>& result_b_matrix,
    vector<Ciphertext>& result_matrix)
{
    if (result_a_matrix.empty() || result_b_matrix.empty()) {
        cerr << "Empty result matrices" << endl;
        return;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = result_a_matrix.size();
    size_t num_ciphertexts = result_a_matrix[0].size();
    
    // 检查维度
    if (result_b_matrix.size() != coeff_modulus_size) {
        cerr << "RNS layer count mismatch in result matrices" << endl;
        return;
    }
    
    if (result_matrix.size() != num_ciphertexts) {
        cerr << "Result matrix size mismatch" << endl;
        return;
    }
    
    cout << "开始构建密文..." << endl;
    cout << "密文数量: " << num_ciphertexts << ", RNS层数: " << coeff_modulus_size << endl;
    
    // 对每个ciphertext位置构建密文
    for (size_t cipher_idx = 0; cipher_idx < num_ciphertexts; cipher_idx++) {
        // 显示进度
        cout << "\r构建密文进度: " << (cipher_idx + 1) << "/" << num_ciphertexts 
             << " [" << (cipher_idx + 1) * 100 / num_ciphertexts << "%]" << flush;
        
        // 创建新的密文
        result_matrix[cipher_idx].resize(context, 2); // 2个多项式：a和b
        
        // 设置a多项式的系数
        if (result_matrix[cipher_idx].size() > 0) {
            uint64_t* coeffs_a = result_matrix[cipher_idx].data(0);
            for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeffs_a[coeff_idx + rns_layer * poly_modulus_degree] = 
                        result_a_matrix[rns_layer][cipher_idx][coeff_idx];
                }
            }
        }
        
        // 设置b多项式的系数
        if (result_matrix[cipher_idx].size() > 1) {
            uint64_t* coeffs_b = result_matrix[cipher_idx].data(1);
            for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeffs_b[coeff_idx + rns_layer * poly_modulus_degree] = 
                        result_b_matrix[rns_layer][cipher_idx][coeff_idx];
                }
            }
        }
    }
    
    cout << "\n密文构建完成！" << endl;
}

// 将矩阵行编码到明文多项式系数上的函数
void encode_matrix_row_to_plaintext(
    const vector<uint64_t>& matrix_row,
    const SEALContext& context,
    Plaintext& plaintext)
{
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    // 检查行大小
    if (matrix_row.size() != poly_modulus_degree) {
        cerr << "Matrix row size must be " << poly_modulus_degree 
             << ", but got " << matrix_row.size() << endl;
        return;
    }
    
    // 创建明文多项式
    if (context_data->parms().scheme() == scheme_type::ckks)
        plaintext.resize(poly_modulus_degree * 2);
    else
        plaintext.resize(poly_modulus_degree);
    
    plaintext.set_zero();
    
    // 使用 modulo_poly_coeffs 将矩阵行编码到多项式系数上
    util::modulo_poly_coeffs(
        matrix_row.data(), 
        poly_modulus_degree, 
        (context_data->parms().scheme() == scheme_type::ckks 
            ? context_data->parms().coeff_modulus()[0] 
            : context_data->parms().plain_modulus()),
        plaintext.data()
    );
    if (context_data->parms().scheme() == scheme_type::ckks) {
        util::ntt_negacyclic_harvey(plaintext.data(), context_data->small_ntt_tables()[0]);
        plaintext.parms_id() = context.key_parms_id();
        // assert(plaintext.parms_id() == context.key_parms_id());
    }
}

// 从明文多项式系数解码回矩阵行
void decode_plaintext_to_matrix_row(
    const Plaintext& plaintext,
    const SEALContext& context,
    vector<uint64_t>& matrix_row)
{
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    matrix_row.resize(poly_modulus_degree);
    
    // 直接复制系数
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        matrix_row[i] = plaintext.data()[i];
    }
}

void encrypt_matrix_rows(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& plain_matrix,
    vector<Ciphertext>& encrypted_matrix)
{
    size_t d = plain_matrix.size();
    if (d == 0) {
        cerr << "明文矩阵为空" << endl;
        return;
    }
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    // 检查每行长度
    for (size_t i = 0; i < d; i++) {
        if (plain_matrix[i].size() != poly_modulus_degree) {
            cerr << "明文矩阵第" << i << "行长度错误，应为" << poly_modulus_degree << endl;
            return;
        }
    }
    encrypted_matrix.resize(d);
    for (size_t i = 0; i < d; i++) {
        Plaintext plaintext;
        encode_matrix_row_to_plaintext(plain_matrix[i], context, plaintext);
        encryptor.encrypt(plaintext, encrypted_matrix[i]);
    }
}

void decrypt_ciphertexts_to_matrix(
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_matrix,
    vector<vector<uint64_t>>& plain_matrix)
{
    size_t d = encrypted_matrix.size();
    plain_matrix.resize(d);
    for (size_t i = 0; i < d; i++) {
        Plaintext pt;
        decryptor.decrypt(encrypted_matrix[i], pt);
        decode_plaintext_to_matrix_row(pt, context, plain_matrix[i]);
    }
}

void matrix_multiply_plain(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus)
{
    size_t n = A.size();
    size_t m = A[0].size();
    size_t p = B[0].size();
    // 检查维度
    if (B.size() != m) {
        cerr << "矩阵乘法维度不匹配: A is " << n << "x" << m << ", B is " << B.size() << "x" << p << endl;
        C.clear();
        return;
    }
    C.assign(n, vector<uint64_t>(p, 0));
    for (size_t i = 0; i < n; i++) {
        // 显示进度
        cout << "\r矩阵乘法进度: 处理第 " << (i + 1) << "/" << n << " 行 [" 
             << (i + 1) * 100 / n << "%]" << flush;
        
        for (size_t j = 0; j < p; j++) {
            uint64_t sum = 0;
            for (size_t k = 0; k < m; k++) {
                uint64_t a = A[i][k];
                uint64_t b = B[k][j];
                uint64_t prod = (a * b);
                sum = (sum + prod) % modulus;
            }
            C[i][j] = sum;
        }
    }
    cout << endl; // 换行
}

void matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus)
{
    size_t m = A.size();
    if (m == 0) return;
    size_t k = A[0].size();
    if (B.size() != k) {
        cerr << "Matrix dimensions do not match for multiplication." << endl;
        return;
    }
    size_t n = B[0].size();

    vector<double> A_double(m * k);
    vector<double> B_double(k * n);
    vector<double> C_double(m * n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < k; ++j) {
            A_double[i * k + j] = static_cast<double>(A[i][j]);
        }
    }

    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < n; ++j) {
            B_double[i * n + j] = static_cast<double>(B[i][j]);
        }
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0,
                A_double.data(), k,
                B_double.data(), n, 0.0,
                C_double.data(), n);
    
    C.assign(m, vector<uint64_t>(n));
    double d_modulus = static_cast<double>(modulus);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = static_cast<uint64_t>(fmod(C_double[i * n + j], d_modulus));
        }
    }
}
