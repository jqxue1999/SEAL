#include "cpmm.h"

// AVX辅助函数：将uint64_t数组转换为double数组
#ifdef __AVX2__
void convert_uint64_to_double_avx(const uint64_t* src, double* dst, size_t size) {
    size_t vec_size = size / 4 * 4; // 处理4个元素一组
    
    // 并行处理AVX向量化部分
    #pragma omp parallel for
    for (size_t i = 0; i < vec_size; i += 4) {
        // 加载4个uint64_t
        __m256i uint64_vec = _mm256_loadu_si256((__m256i*)(src + i));
        
        // 转换为double (使用AVX2兼容的方法)
        __m128i low = _mm256_extracti128_si256(uint64_vec, 0);
        __m128i high = _mm256_extracti128_si256(uint64_vec, 1);
        
        // 使用标量转换，因为_mm_cvtepi64_pd需要AVX-512
        // 必须使用编译时常量作为选择器
        int64_t low_val_0 = _mm_extract_epi64(low, 0);
        int64_t low_val_1 = _mm_extract_epi64(low, 1);
        dst[i + 0] = static_cast<double>(low_val_0);
        dst[i + 1] = static_cast<double>(low_val_1);
        
        int64_t high_val_0 = _mm_extract_epi64(high, 0);
        int64_t high_val_1 = _mm_extract_epi64(high, 1);
        dst[i + 2] = static_cast<double>(high_val_0);
        dst[i + 3] = static_cast<double>(high_val_1);
    }
    
    // 并行处理剩余元素
    #pragma omp parallel for
    for (size_t i = vec_size; i < size; ++i) {
        dst[i] = static_cast<double>(src[i]);
    }
}

void convert_double_to_uint64_avx(const double* src, uint64_t* dst, size_t size) {
    size_t vec_size = size / 4 * 4;
    
    // 并行处理AVX向量化部分
    #pragma omp parallel for
    for (size_t i = 0; i < vec_size; i += 4) {
        // 加载4个double
        __m256d double_vec = _mm256_loadu_pd(src + i);
        
        // 使用标量转换，因为_mm256_cvtpd_epi64需要AVX-512
        // 提取double值并转换
        double vals[4];
        _mm256_storeu_pd(vals, double_vec);
        
        for (int j = 0; j < 4; ++j) {
            dst[i + j] = static_cast<uint64_t>(vals[j]);
        }
    }
    
    // 并行处理剩余元素
    #pragma omp parallel for
    for (size_t i = vec_size; i < size; ++i) {
        dst[i] = static_cast<uint64_t>(src[i]);
    }
}
#endif

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
        
        cout << "  多项式模数次数: " << poly_modulus_degree << endl;
        cout << "  系数模数层数: " << coeff_modulus_size << endl;

        // 初始化结果矩阵
        result_matrix.resize(encrypted_matrix.size());
        
        // 计时变量
        double step1_time = 0;
        double step2_time = 0;
        double step3_time = 0;
        
        cout << "\n步骤1/3: 提取系数..." << endl;
        
        // 步骤1: 从ciphertext矩阵每一行的ciphertext提取系数，组成矩阵
        vector<vector<vector<uint64_t>>> coeff_matrix_a, coeff_matrix_b;
        vector<uint64_t> modulus_vector;
        step1_time += extract_coefficients_from_ciphertext_vector(context, encrypted_matrix, coeff_matrix_a, coeff_matrix_b, modulus_vector);
                
        cout << "\n步骤2/3: 执行矩阵乘法..." << endl;
        
        // 步骤2: 对每个RNS片，计算plaintext matrix * a和plaintext matrix * b
        vector<vector<vector<uint64_t>>> result_a_matrix, result_b_matrix;
        
        cout << "处理a矩阵..." << endl;
        step2_time += Normal_RNS_multiply(plain_matrix, coeff_matrix_a, result_a_matrix, modulus_vector);
        cout << "处理b矩阵..." << endl;
        step2_time += Normal_RNS_multiply(plain_matrix, coeff_matrix_b, result_b_matrix, modulus_vector);

        cout << "\n步骤3/3: 构建密文..." << endl;        
        // 步骤3: 构建d个ciphertext，按行构建
        step3_time += build_ciphertexts_from_result_matrices_2(context, result_a_matrix, result_b_matrix, result_matrix);
        
        cout << "\n[CPMM] ========== 密文-明文矩阵乘法时间统计 ==========" << endl;
        cout << "Step 1 (提取系数): " << step1_time << " seconds" << endl;
        cout << "Step 2 (矩阵乘法): " << step2_time << " seconds" << endl;
        cout << "Step 3 (构建密文): " << step3_time << " seconds" << endl;
        cout << "----------------------------------------" << endl;
        cout << "总时间:              " << step1_time + step2_time + step3_time << " seconds" << endl;
        cout << "=========================================" << endl;
        
        cout << "\n=== 密文-明文矩阵乘法完成 ===" << endl;
        
        return true;
    }
    catch (const exception& e) {
        cerr << "Error in ciphertext_plaintext_matrix_multiply: " << e.what() << endl;
        return false;
    }
}

double extract_coefficients_from_ciphertext_vector(
    const SEALContext& context,
    const vector<Ciphertext>& encrypted_vector,
    vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector)
{
    auto start_time = chrono::high_resolution_clock::now();
    if (encrypted_vector.empty()) {
        cerr << "Empty encrypted vector" << endl;
        return 0;
    }
    
    auto context_data = context.get_context_data(encrypted_vector[0].parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = encrypted_vector[0].coeff_modulus_size();
    size_t num_ciphertexts = encrypted_vector.size();
    
    // 检查向量大小
    if (num_ciphertexts != poly_modulus_degree) {
        cerr << "Encrypted vector size must be " << poly_modulus_degree 
             << ", but got " << num_ciphertexts << endl;
        return 0;
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

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}

double Normal_RNS_multiply(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (plain_matrix.empty() || coeff_matrix.empty()) {
        cerr << "Empty matrices in multiplication" << endl;
        return 0;
    }

    size_t plain_rows = plain_matrix.size();
    size_t plain_cols = plain_matrix[0].size();
    size_t rns_layers = coeff_matrix.size();

    if (modulus_vector.size() != rns_layers) {
        cerr << "RNS layer count or modulus vector size mismatch" << endl;
        return 0;
    }

    cout << "RNS层数: " << rns_layers << ", 矩阵大小: " << plain_rows << "x" << plain_cols << endl;

    double total_time = 0;
    result_matrix.resize(rns_layers);

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        size_t matrix_rows = coeff_matrix[rns_layer].size();
        size_t matrix_cols = coeff_matrix[rns_layer][0].size();
        
        total_time += matrix_multiply_plain_blas(
            plain_matrix, 
            coeff_matrix[rns_layer], 
            result_matrix[rns_layer], 
            modulus_vector[rns_layer]
        );
    }

    return total_time;
}

double build_ciphertexts_from_result_matrices_2(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_a_matrix,
    const vector<vector<vector<uint64_t>>>& result_b_matrix,
    vector<Ciphertext>& result_matrix)
{
    auto start_time = chrono::high_resolution_clock::now();
    if (result_a_matrix.empty() || result_b_matrix.empty()) {
        cerr << "Empty result matrices" << endl;
        return 0;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = result_a_matrix.size();
    size_t num_ciphertexts = result_a_matrix[0].size();
    
    // 检查维度
    if (result_b_matrix.size() != coeff_modulus_size) {
        cerr << "RNS layer count mismatch in result matrices" << endl;
        return 0;
    }
    
    if (result_matrix.size() != num_ciphertexts) {
        cerr << "Result matrix size mismatch" << endl;
        return 0;
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

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}

// 将向量编码到明文多项式系数上的函数
void encode_vector_to_plaintext(
    const vector<uint64_t>& vector_data,
    const SEALContext& context,
    Plaintext& plaintext)
{
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    // 检查向量大小
    if (vector_data.size() != poly_modulus_degree) {
        cerr << "Vector size must be " << poly_modulus_degree 
             << ", but got " << vector_data.size() << endl;
        return;
    }
    
    // 创建明文多项式
    if (context_data->parms().scheme() == scheme_type::ckks)
        plaintext.resize(poly_modulus_degree * 2);
    else
        plaintext.resize(poly_modulus_degree);
    
    plaintext.set_zero();
    
    // 使用 modulo_poly_coeffs 将向量编码到多项式系数上
    util::modulo_poly_coeffs(
        vector_data.data(), 
        poly_modulus_degree, 
        (context_data->parms().scheme() == scheme_type::ckks 
            ? context_data->parms().coeff_modulus()[0] 
            : context_data->parms().plain_modulus()),
        plaintext.data()
    );
    if (context_data->parms().scheme() == scheme_type::ckks) {
        util::ntt_negacyclic_harvey(plaintext.data(), context_data->small_ntt_tables()[0]);
        plaintext.parms_id() = context.key_parms_id();
    }
}

// 从明文多项式系数解码回向量
void decode_plaintext_to_vector(
    const Plaintext& plaintext,
    const SEALContext& context,
    vector<uint64_t>& vector_data)
{
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    vector_data.resize(poly_modulus_degree);
    
    // 直接复制系数
    for (size_t i = 0; i < poly_modulus_degree; i++) {
        vector_data[i] = plaintext.data()[i];
    }
}

void encrypt_matrix(
    const bool is_row,
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& plain_matrix,
    vector<Ciphertext>& encrypted_matrix)
{
    if (plain_matrix.empty()) {
        cerr << "明文矩阵为空" << endl;
        return;
    }
    
    size_t rows = plain_matrix.size();
    size_t cols = plain_matrix[0].size();
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    if (is_row) {
        // 按行编码
        // 检查每行长度
        for (size_t i = 0; i < rows; i++) {
            if (plain_matrix[i].size() != poly_modulus_degree) {
                cerr << "明文矩阵第" << i << "行长度错误，应为" << poly_modulus_degree << endl;
                return;
            }
        }
        encrypted_matrix.resize(rows);
        for (size_t i = 0; i < rows; i++) {
            Plaintext plaintext;
            encode_vector_to_plaintext(plain_matrix[i], context, plaintext);
            encryptor.encrypt(plaintext, encrypted_matrix[i]);
        }
    } else {
        // 按列编码
        // 检查矩阵维度
        if (rows != poly_modulus_degree) {
            cerr << "明文矩阵行数错误，应为" << poly_modulus_degree << "，实际为" << rows << endl;
            return;
        }
        
        // 检查每行长度
        for (size_t i = 0; i < rows; i++) {
            if (plain_matrix[i].size() != cols) {
                cerr << "明文矩阵第" << i << "行长度不一致" << endl;
                return;
            }
        }
        
        encrypted_matrix.resize(cols);
        
        // 对每一列进行编码和加密
        for (size_t j = 0; j < cols; j++) {
            // 提取第j列
            vector<uint64_t> matrix_col(rows);
            for (size_t i = 0; i < rows; i++) {
                matrix_col[i] = plain_matrix[i][j];
            }
            
            // 编码并加密
            Plaintext plaintext;
            encode_vector_to_plaintext(matrix_col, context, plaintext);
            encryptor.encrypt(plaintext, encrypted_matrix[j]);
        }
    }
}

void decrypt_ciphertexts_to_matrix(
    const bool is_row,
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_matrix,
    vector<vector<uint64_t>>& plain_matrix)
{
    if (encrypted_matrix.empty()) {
        cerr << "密文矩阵为空" << endl;
        return;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    
    if (is_row) {
        // 按行解码
        size_t d = encrypted_matrix.size();
        plain_matrix.resize(d);
        for (size_t i = 0; i < d; i++) {
            Plaintext pt;
            decryptor.decrypt(encrypted_matrix[i], pt);
            decode_plaintext_to_vector(pt, context, plain_matrix[i]);
        }
    } else {
        // 按列解码
        size_t cols = encrypted_matrix.size();
        
        plain_matrix.resize(poly_modulus_degree);
        for (size_t i = 0; i < poly_modulus_degree; i++) {
            plain_matrix[i].resize(cols);
        }
        
        // 对每个密文进行解密并提取对应的列
        for (size_t j = 0; j < cols; j++) {
            Plaintext pt;
            decryptor.decrypt(encrypted_matrix[j], pt);
            
            vector<uint64_t> matrix_col;
            decode_plaintext_to_vector(pt, context, matrix_col);
            
            // 将列数据填充到矩阵中
            for (size_t i = 0; i < poly_modulus_degree; i++) {
                plain_matrix[i][j] = matrix_col[i];
            }
        }
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

void transpose_matrix_blas(
    const vector<vector<vector<uint64_t>>>& input_matrix,
    vector<vector<vector<uint64_t>>>& output_matrix)
{
    if (input_matrix.empty()) {
        cerr << "Input matrix is empty" << endl;
        return;
    }
    
    size_t rns_layers = input_matrix.size();
    size_t rows = input_matrix[0].size();
    size_t cols = input_matrix[0][0].size();
    
    // 检查所有RNS层的矩阵大小是否一致
    for (size_t rns = 0; rns < rns_layers; rns++) {
        if (input_matrix[rns].size() != rows || input_matrix[rns][0].size() != cols) {
            cerr << "Matrix dimensions are not consistent across RNS layers" << endl;
            return;
        }
    }
    
    output_matrix.resize(rns_layers);
    
    cout << "开始使用BLAS进行矩阵转置..." << endl;
    cout << "RNS层数: " << rns_layers << ", 矩阵大小: " << rows << "x" << cols << endl;
    
    for (size_t rns = 0; rns < rns_layers; rns++) {
        cout << "\r处理RNS层: " << (rns + 1) << "/" << rns_layers 
             << " [" << (rns + 1) * 100 / rns_layers << "%]" << flush;
        
        // 初始化输出矩阵
        output_matrix[rns].resize(cols, vector<uint64_t>(rows));
        
        // 转换为double数组
        vector<double> input_double(rows * cols);
        vector<double> output_double(cols * rows);
        
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                input_double[i * cols + j] = static_cast<double>(input_matrix[rns][i][j]);
            }
        }
        
        // 使用BLAS的domatcopy进行转置
        // CblasRowMajor: 行主序
        // CblasTrans: 转置操作
        // rows, cols: 源矩阵的行数和列数
        // 1.0: alpha值（缩放因子）
        // input_double.data(), cols: 源矩阵数据，源矩阵的列数（步长）
        // output_double.data(), rows: 目标矩阵数据，目标矩阵的列数（步长）
        cblas_domatcopy(CblasRowMajor, CblasTrans, rows, cols, 1.0,
                        input_double.data(), cols, output_double.data(), rows);
        
        // 转换回uint64_t
        for (size_t i = 0; i < cols; i++) {
            for (size_t j = 0; j < rows; j++) {
                output_matrix[rns][i][j] = static_cast<uint64_t>(output_double[i * rows + j]);
            }
        }
    }
    
    cout << "\nBLAS矩阵转置完成！" << endl;
}

double matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus)
{
    chrono::duration<double> convert_time(0);
    chrono::duration<double> blas_time(0);
    chrono::duration<double> mod_time(0);
    
    size_t m = A.size();
    if (m == 0) return 0;
    size_t k = A[0].size();
    if (B.size() != k) {
        cerr << "Matrix dimensions do not match for multiplication." << endl;
        return 0;
    }
    size_t n = B[0].size();

    // 转换为 double 数组
    auto convert_start = chrono::high_resolution_clock::now();
    vector<double> A_double(m * k);
    vector<double> B_double(k * n);
    vector<double> C_double(m * n);

    #pragma omp parallel
    {
        #ifdef __AVX2__
        // 使用AVX优化的转换
        #pragma omp for
        for (size_t i = 0; i < m; ++i) {
            // 将A的第i行转换为double数组
            convert_uint64_to_double_avx(A[i].data(), A_double.data() + i * k, k);
        }
        
        #pragma omp for
        for (size_t i = 0; i < k; ++i) {
            // 将B的第i行转换为double数组
            convert_uint64_to_double_avx(B[i].data(), B_double.data() + i * n, n);
        }
        #else
        cout << "使用原始转换方法" << endl;
        // 原始转换方法（作为fallback）
        #pragma omp for collapse(2)
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < k; ++j)
                A_double[i * k + j] = static_cast<double>(A[i][j]);
        
        #pragma omp for collapse(2)
        for (size_t i = 0; i < k; ++i)
            for (size_t j = 0; j < n; ++j)
                B_double[i * n + j] = static_cast<double>(B[i][j]);
        #endif
    }

    auto convert_end = chrono::high_resolution_clock::now();
    convert_time = chrono::duration<double>(convert_end - convert_start);

    // 执行 BLAS 矩阵乘法
    auto blas_start = chrono::high_resolution_clock::now();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0,
                A_double.data(), k,
                B_double.data(), n, 0.0,
                C_double.data(), n);
    auto blas_end = chrono::high_resolution_clock::now();
    blas_time = chrono::duration<double>(blas_end - blas_start);

    // 先转换为uint64_t，统计到convert_time
    auto convert_start2 = chrono::high_resolution_clock::now();
    vector<vector<uint64_t>> C_int(m, vector<uint64_t>(n));
    
    #ifdef __AVX2__
    // 使用AVX优化的转换
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < m; ++i) {
            convert_double_to_uint64_avx(C_double.data() + i * n, C_int[i].data(), n);
        }
    }
    #else
    // 原始转换方法（作为fallback）
    cout << "使用原始转换方法" << endl;

    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C_int[i][j] = static_cast<uint64_t>(C_double[i * n + j]);
        }
    }

    #endif
    auto convert_end2 = chrono::high_resolution_clock::now();
    convert_time += chrono::duration<double>(convert_end2 - convert_start2);

    // 应用模运算并存储结果，统计到mod_time
    auto mod_start = chrono::high_resolution_clock::now();
    C.assign(m, vector<uint64_t>(n));
    Modulus mod_obj(modulus);
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j)
                C[i][j] = util::barrett_reduce_64(C_int[i][j], mod_obj);
        }
    }
    auto mod_end = chrono::high_resolution_clock::now();
    mod_time = chrono::duration<double>(mod_end - mod_start);

    cout << "convert_time: " << convert_time.count() << " | blas_time: " << blas_time.count() << " | mod_time: " << mod_time.count() << endl;

    return convert_time.count() + blas_time.count() + mod_time.count();
}


double matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C)
{
    chrono::duration<double> convert_time(0);
    chrono::duration<double> blas_time(0);
    
    size_t m = A.size();
    if (m == 0) return 0;
    size_t k = A[0].size();
    if (B.size() != k) {
        cerr << "Matrix dimensions do not match for multiplication." << endl;
        return 0;
    }
    size_t n = B[0].size();

    // 转换为 double 数组
    auto convert_start = chrono::high_resolution_clock::now();
    vector<double> A_double(m * k);
    vector<double> B_double(k * n);
    vector<double> C_double(m * n);

    #pragma omp parallel
    {
        #ifdef __AVX2__
        // 使用AVX优化的转换
        #pragma omp for
        for (size_t i = 0; i < m; ++i) {
            // 将A的第i行转换为double数组
            convert_uint64_to_double_avx(A[i].data(), A_double.data() + i * k, k);
        }
        
        #pragma omp for
        for (size_t i = 0; i < k; ++i) {
            // 将B的第i行转换为double数组
            convert_uint64_to_double_avx(B[i].data(), B_double.data() + i * n, n);
        }
        #else
        cout << "使用原始转换方法" << endl;
        // 原始转换方法（作为fallback）
        #pragma omp for collapse(2)
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < k; ++j)
                A_double[i * k + j] = static_cast<double>(A[i][j]);
        
        #pragma omp for collapse(2)
        for (size_t i = 0; i < k; ++i)
            for (size_t j = 0; j < n; ++j)
                B_double[i * n + j] = static_cast<double>(B[i][j]);
        #endif
    }

    auto convert_end = chrono::high_resolution_clock::now();
    convert_time = chrono::duration<double>(convert_end - convert_start);

    // 执行 BLAS 矩阵乘法
    auto blas_start = chrono::high_resolution_clock::now();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0,
                A_double.data(), k,
                B_double.data(), n, 0.0,
                C_double.data(), n);
    auto blas_end = chrono::high_resolution_clock::now();
    blas_time = chrono::duration<double>(blas_end - blas_start);

    // 先转换为uint64_t，统计到convert_time
    auto convert_start2 = chrono::high_resolution_clock::now();
    C.assign(m, vector<uint64_t>(n));
    
    #ifdef __AVX2__
    // 使用AVX优化的转换
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < m; ++i) {
            convert_double_to_uint64_avx(C_double.data() + i * n, C[i].data(), n);
        }
    }
    #else
    // 原始转换方法（作为fallback）
    cout << "使用原始转换方法" << endl;

    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = static_cast<uint64_t>(C_double[i * n + j]);
        }
    }

    #endif
    auto convert_end2 = chrono::high_resolution_clock::now();
    convert_time += chrono::duration<double>(convert_end2 - convert_start2);

    cout << "convert_time: " << convert_time.count() << " | blas_time: " << blas_time.count() << endl;

    return convert_time.count() + blas_time.count();
}