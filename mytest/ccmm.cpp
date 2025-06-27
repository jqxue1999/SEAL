#include "ccmm.h"
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include "seal/util/uintarithsmallmod.h"
#include <chrono>

void RNS_RNS_multiply(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (coeff_matrix_1.empty() || coeff_matrix_2.empty()) {
        cerr << "Empty matrices in multiplication" << endl;
        return;
    }

    size_t rns_layers = coeff_matrix_1.size();
    
    // 检查RNS层数是否一致
    if (coeff_matrix_2.size() != rns_layers) {
        cerr << "RNS layer count mismatch: matrix1 has " << rns_layers 
             << " layers, matrix2 has " << coeff_matrix_2.size() << " layers" << endl;
        return;
    }
    
    if (modulus_vector.size() != rns_layers) {
        cerr << "Modulus vector size mismatch: expected " << rns_layers 
             << ", got " << modulus_vector.size() << endl;
        return;
    }

    cout << "开始使用BLAS进行RNS-RNS矩阵乘法计算..." << endl;
    cout << "RNS层数: " << rns_layers << endl;

    result_matrix.resize(rns_layers);

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        size_t matrix1_rows = coeff_matrix_1[rns_layer].size();
        size_t matrix1_cols = coeff_matrix_1[rns_layer][0].size();
        size_t matrix2_rows = coeff_matrix_2[rns_layer].size();
        size_t matrix2_cols = coeff_matrix_2[rns_layer][0].size();
        
        // 检查矩阵维度是否匹配
        if (matrix1_cols != matrix2_rows) {
            cerr << "\nMatrix dimensions do not match for multiplication at RNS layer " << rns_layer 
                 << ": matrix1 is " << matrix1_rows << "x" << matrix1_cols 
                 << ", matrix2 is " << matrix2_rows << "x" << matrix2_cols << endl;
            return;
        }

        result_matrix[rns_layer].resize(matrix1_rows, vector<uint64_t>(matrix2_cols));

        // 转换为double数组
        std::vector<double> matrix1_double(matrix1_rows * matrix1_cols);
        std::vector<double> matrix2_double(matrix2_rows * matrix2_cols);
        
        for(size_t i = 0; i < matrix1_rows; ++i) {
            for (size_t j = 0; j < matrix1_cols; ++j) {
                matrix1_double[i * matrix1_cols + j] = static_cast<double>(coeff_matrix_1[rns_layer][i][j]);
            }
        }
        
        for(size_t i = 0; i < matrix2_rows; ++i) {
            for (size_t j = 0; j < matrix2_cols; ++j) {
                matrix2_double[i * matrix2_cols + j] = static_cast<double>(coeff_matrix_2[rns_layer][i][j]);
            }
        }

        std::vector<double> result_double(matrix1_rows * matrix2_cols);

        // 使用BLAS进行矩阵乘法: result = matrix1 * matrix2
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    matrix1_rows, matrix2_cols, matrix1_cols, 1.0,
                    matrix1_double.data(), matrix1_cols,
                    matrix2_double.data(), matrix2_cols, 0.0,
                    result_double.data(), matrix2_cols);
        
        // 转换回uint64_t并应用模运算
        uint64_t modulus = modulus_vector[rns_layer];
        for(size_t i = 0; i < matrix1_rows; ++i) {
            for(size_t j = 0; j < matrix2_cols; ++j) {
                result_matrix[rns_layer][i][j] = static_cast<uint64_t>(fmod(result_double[i * matrix2_cols + j], modulus));
            }
        }
    }

    cout << "\nBLAS RNS-RNS矩阵乘法计算完成！" << endl;
}

bool ciphertext_ciphertext_matrix_multiply(
    const SEALContext& context,
    const RelinKeys& relin_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& encrypted_matrix_1,
    const vector<Ciphertext>& encrypted_matrix_2,
    vector<Ciphertext>& result_matrix)
{
    try {
        cout << "\n=== 开始密文-密文矩阵乘法 ===" << endl;
        
        // 获取参数
        auto context_data = context.get_context_data(encrypted_matrix_1[0].parms_id());
        if (!context_data) {
            cerr << "Invalid context data" << endl;
            return false;
        }
        
        size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
        size_t coeff_modulus_size = encrypted_matrix_1[0].coeff_modulus_size();
        
        size_t encrypted_rows = encrypted_matrix_1.size();
        size_t encrypted_cols = encrypted_matrix_2.size();
        
        // 检查矩阵大小
        if (encrypted_rows != poly_modulus_degree || encrypted_cols != poly_modulus_degree) {
            cerr << "Encrypted matrix size must be " << poly_modulus_degree 
                 << ", but got " << encrypted_rows << endl;
            return false;
        }
        
        cout << "多项式模数次数: " << poly_modulus_degree << endl;
        cout << "系数模数层数: " << coeff_modulus_size << endl;
        
        // 初始化结果矩阵
        result_matrix.resize(encrypted_rows);
        
        cout << "\n步骤1/4: 提取系数..." << endl;
        // 步骤0: 从ciphertext矩阵每一行的ciphertext提取系数，组成矩阵
        vector<vector<vector<uint64_t>>> coeff_matrix_1_a, coeff_matrix_1_b, coeff_matrix_2_a, coeff_matrix_2_b;
        vector<uint64_t> modulus_vector;
        extract_coefficients_from_ciphertext_vector(context, encrypted_matrix_1, coeff_matrix_1_a, coeff_matrix_1_b, modulus_vector);
        extract_coefficients_from_ciphertext_vector(context, encrypted_matrix_2, coeff_matrix_2_a, coeff_matrix_2_b, modulus_vector);
        
        cout << "\n步骤2/4: 转置系数矩阵..." << endl;
        // 步骤1: 转置系数矩阵
        vector<vector<vector<uint64_t>>> coeff_matrix_1_a_T, coeff_matrix_1_b_T;
        transpose_matrix_blas(coeff_matrix_1_a, coeff_matrix_1_a_T);
        transpose_matrix_blas(coeff_matrix_1_b, coeff_matrix_1_b_T);

        // vector<vector<vector<uint64_t>>> coeff_matrix_1_a_T_flint, coeff_matrix_1_b_T_flint;
        // transpose_matrix_flint(coeff_matrix_1_a, coeff_matrix_1_a_T_flint);
        // transpose_matrix_flint(coeff_matrix_1_b, coeff_matrix_1_b_T_flint);

        cout << "\n步骤3/4: 执行矩阵乘法..." << endl;
        vector<vector<vector<uint64_t>>> result_c00_matrix, result_c01_matrix, result_c10_matrix, result_c11_matrix;
        RNS_RNS_multiply(coeff_matrix_1_a_T, coeff_matrix_2_a, result_c00_matrix, modulus_vector);
        RNS_RNS_multiply(coeff_matrix_1_a_T, coeff_matrix_2_b, result_c01_matrix, modulus_vector);
        RNS_RNS_multiply(coeff_matrix_1_b_T, coeff_matrix_2_a, result_c10_matrix, modulus_vector);
        RNS_RNS_multiply(coeff_matrix_1_b_T, coeff_matrix_2_b, result_c11_matrix, modulus_vector);

        // 可选：使用FLINT版本进行大整数矩阵运算
        // RNS_RNS_multiply_flint(coeff_matrix_1_a_T, coeff_matrix_2_a, result_c00_matrix, modulus_vector);
        // RNS_RNS_multiply_flint(coeff_matrix_1_a_T, coeff_matrix_2_b, result_c01_matrix, modulus_vector);
        // RNS_RNS_multiply_flint(coeff_matrix_1_b_T, coeff_matrix_2_a, result_c10_matrix, modulus_vector);
        // RNS_RNS_multiply_flint(coeff_matrix_1_b_T, coeff_matrix_2_b, result_c11_matrix, modulus_vector);

        vector<vector<vector<uint64_t>>> result_d0_matrix, result_d1_matrix, result_d2_matrix, result_d3_matrix;
        transpose_matrix_blas(result_c00_matrix, result_d0_matrix);
        // transpose_matrix_flint(result_c00_matrix, result_d0_matrix);
        result_d1_matrix.resize(result_c00_matrix.size(), vector<vector<uint64_t>>(result_c00_matrix[0].size(), vector<uint64_t>(result_c00_matrix[0][0].size(), 0)));
        transpose_matrix_blas(result_c01_matrix, result_d2_matrix);
        // transpose_matrix_flint(result_c01_matrix, result_d2_matrix);
        result_d3_matrix.resize(result_c01_matrix.size(), vector<vector<uint64_t>>(result_c01_matrix[0].size(), vector<uint64_t>(result_c01_matrix[0][0].size(), 0)));

        // C0 = D3 + C11
        vector<vector<vector<uint64_t>>> result_c0_matrix;
        matrix_add_plain_blas(result_d3_matrix, result_c11_matrix, result_c0_matrix, modulus_vector);
        // matrix_add_plain_flint(result_d3_matrix, result_c11_matrix, result_c0_matrix, modulus_vector);

        // C1 = D1 + D2 + C10
        vector<vector<vector<uint64_t>>> result_c1_matrix;
        vector<vector<vector<uint64_t>>> temp_matrix;
        matrix_add_plain_blas(result_d1_matrix, result_d2_matrix, temp_matrix, modulus_vector);
        // matrix_add_plain_flint(result_d1_matrix, result_d2_matrix, temp_matrix, modulus_vector);
        matrix_add_plain_blas(temp_matrix, result_c10_matrix, result_c1_matrix, modulus_vector);
        // matrix_add_plain_flint(temp_matrix, result_c10_matrix, result_c1_matrix, modulus_vector);
        
        // C2 = D0
        vector<vector<vector<uint64_t>>> result_c2_matrix = result_d0_matrix;
        
        cout << "\n步骤4/4: 构建密文..." << endl;
        build_ciphertexts_from_result_matrices_3(context, result_c0_matrix, result_c1_matrix, result_c2_matrix, result_matrix);
        for (size_t i = 0; i < result_matrix.size(); i++) {
            evaluator.relinearize_inplace(result_matrix[i], relin_keys);
        }
        cout << "\n=== 密文-密文矩阵乘法完成 ===" << endl;
        return true;
    }
    catch (const exception& e) {
        cerr << "Error in ciphertext_plaintext_matrix_multiply: " << e.what() << endl;
        return false;
    }
}

void matrix_add_plain_blas(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (coeff_matrix_1.empty() || coeff_matrix_2.empty()) {
        cerr << "Empty matrices in addition" << endl;
        return;
    }

    size_t rns_layers = coeff_matrix_1.size();
    
    // 检查RNS层数是否一致
    if (coeff_matrix_2.size() != rns_layers) {
        cerr << "RNS layer count mismatch: matrix1 has " << rns_layers 
             << " layers, matrix2 has " << coeff_matrix_2.size() << " layers" << endl;
        return;
    }
    
    if (modulus_vector.size() != rns_layers) {
        cerr << "Modulus vector size mismatch: expected " << rns_layers 
             << ", got " << modulus_vector.size() << endl;
        return;
    }

    cout << "开始使用BLAS进行RNS-RNS矩阵加法计算..." << endl;
    cout << "RNS层数: " << rns_layers << endl;

    result_matrix.resize(rns_layers);

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        size_t matrix1_rows = coeff_matrix_1[rns_layer].size();
        size_t matrix1_cols = coeff_matrix_1[rns_layer][0].size();
        size_t matrix2_rows = coeff_matrix_2[rns_layer].size();
        size_t matrix2_cols = coeff_matrix_2[rns_layer][0].size();
        
        // 检查矩阵维度是否匹配
        if (matrix1_rows != matrix2_rows || matrix1_cols != matrix2_cols) {
            cerr << "\nMatrix dimensions do not match for addition at RNS layer " << rns_layer 
                 << ": matrix1 is " << matrix1_rows << "x" << matrix1_cols 
                 << ", matrix2 is " << matrix2_rows << "x" << matrix2_cols << endl;
            return;
        }

        result_matrix[rns_layer].resize(matrix1_rows, vector<uint64_t>(matrix1_cols));

        // 转换为double数组
        std::vector<double> matrix1_double(matrix1_rows * matrix1_cols);
        std::vector<double> matrix2_double(matrix2_rows * matrix2_cols);
        std::vector<double> result_double(matrix1_rows * matrix1_cols);
        
        for(size_t i = 0; i < matrix1_rows; ++i) {
            for (size_t j = 0; j < matrix1_cols; ++j) {
                matrix1_double[i * matrix1_cols + j] = static_cast<double>(coeff_matrix_1[rns_layer][i][j]);
            }
        }
        
        for(size_t i = 0; i < matrix2_rows; ++i) {
            for (size_t j = 0; j < matrix2_cols; ++j) {
                matrix2_double[i * matrix2_cols + j] = static_cast<double>(coeff_matrix_2[rns_layer][i][j]);
            }
        }

        // 使用BLAS进行矩阵加法: result = matrix1 + matrix2
        // 先复制matrix1到结果
        cblas_dcopy(matrix1_rows * matrix1_cols, matrix1_double.data(), 1, result_double.data(), 1);
        // 然后加上matrix2
        cblas_daxpy(matrix2_rows * matrix2_cols, 1.0, matrix2_double.data(), 1, result_double.data(), 1);
        
        // 转换回uint64_t并应用模运算
        uint64_t modulus = modulus_vector[rns_layer];
        for(size_t i = 0; i < matrix1_rows; ++i) {
            for(size_t j = 0; j < matrix1_cols; ++j) {
                result_matrix[rns_layer][i][j] = static_cast<uint64_t>(fmod(result_double[i * matrix1_cols + j], modulus));
            }
        }
    }

    cout << "\nBLAS RNS-RNS矩阵加法计算完成！" << endl;
}

void build_ciphertexts_from_result_matrices_3(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_c0_matrix,
    const vector<vector<vector<uint64_t>>>& result_c1_matrix,
    const vector<vector<vector<uint64_t>>>& result_c2_matrix,
    vector<Ciphertext>& result_matrix)
{
    if (result_c0_matrix.empty() || result_c1_matrix.empty() || result_c2_matrix.empty()) {
        cerr << "Empty result matrices" << endl;
        return;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = result_c0_matrix.size();
    size_t num_ciphertexts = result_c0_matrix[0].size();
    
    // 检查维度
    if (result_c1_matrix.size() != coeff_modulus_size || result_c2_matrix.size() != coeff_modulus_size) {
        cerr << "RNS layer count mismatch in result matrices" << endl;
        return;
    }
    
    if (result_matrix.size() != num_ciphertexts) {
        cerr << "Result matrix size mismatch" << endl;
        return;
    }
    
    cout << "开始构建密文（3个多项式）..." << endl;
    cout << "密文数量: " << num_ciphertexts << ", RNS层数: " << coeff_modulus_size << endl;
    
    // 对每个ciphertext位置构建密文
    for (size_t cipher_idx = 0; cipher_idx < num_ciphertexts; cipher_idx++) {
        // 显示进度
        cout << "\r构建密文进度: " << (cipher_idx + 1) << "/" << num_ciphertexts 
             << " [" << (cipher_idx + 1) * 100 / num_ciphertexts << "%]" << flush;
        
        // 创建新的密文（3个多项式：c0, c1, c2）
        result_matrix[cipher_idx].resize(context, 3);
        
        // 设置c0多项式的系数
        if (result_matrix[cipher_idx].size() > 0) {
            uint64_t* coeffs_c0 = result_matrix[cipher_idx].data(0);
            for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeffs_c0[coeff_idx + rns_layer * poly_modulus_degree] = 
                        result_c0_matrix[rns_layer][cipher_idx][coeff_idx];
                }
            }
        }
        
        // 设置c1多项式的系数
        if (result_matrix[cipher_idx].size() > 1) {
            uint64_t* coeffs_c1 = result_matrix[cipher_idx].data(1);
            for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeffs_c1[coeff_idx + rns_layer * poly_modulus_degree] = 
                        result_c1_matrix[rns_layer][cipher_idx][coeff_idx];
                }
            }
        }
        
        // 设置c2多项式的系数
        if (result_matrix[cipher_idx].size() > 2) {
            uint64_t* coeffs_c2 = result_matrix[cipher_idx].data(2);
            for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
                for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                    coeffs_c2[coeff_idx + rns_layer * poly_modulus_degree] = 
                        result_c2_matrix[rns_layer][cipher_idx][coeff_idx];
                }
            }
        }
    }
    
    cout << "\n密文构建完成（3个多项式）！" << endl;
}

void RNS_RNS_multiply_flint(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (coeff_matrix_1.empty() || coeff_matrix_2.empty()) {
        cerr << "Empty matrices in multiplication" << endl;
        return;
    }

    size_t rns_layers = coeff_matrix_1.size();
    
    // 检查RNS层数是否一致
    if (coeff_matrix_2.size() != rns_layers) {
        cerr << "RNS layer count mismatch: matrix1 has " << rns_layers 
             << " layers, matrix2 has " << coeff_matrix_2.size() << " layers" << endl;
        return;
    }
    
    if (modulus_vector.size() != rns_layers) {
        cerr << "Modulus vector size mismatch: expected " << rns_layers 
             << ", got " << modulus_vector.size() << endl;
        return;
    }

    cout << "开始使用FLINT进行RNS-RNS矩阵乘法计算..." << endl;
    cout << "RNS层数: " << rns_layers << endl;

    result_matrix.resize(rns_layers);

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        size_t matrix1_rows = coeff_matrix_1[rns_layer].size();
        size_t matrix1_cols = coeff_matrix_1[rns_layer][0].size();
        size_t matrix2_rows = coeff_matrix_2[rns_layer].size();
        size_t matrix2_cols = coeff_matrix_2[rns_layer][0].size();
        
        // 检查矩阵维度是否匹配
        if (matrix1_cols != matrix2_rows) {
            cerr << "\nMatrix dimensions do not match for multiplication at RNS layer " << rns_layer 
                 << ": matrix1 is " << matrix1_rows << "x" << matrix1_cols 
                 << ", matrix2 is " << matrix2_rows << "x" << matrix2_cols << endl;
            return;
        }

        result_matrix[rns_layer].resize(matrix1_rows, vector<uint64_t>(matrix2_cols));

        // 使用FLINT进行矩阵乘法
        fmpz_t modulus;
        fmpz_init_set_ui(modulus, modulus_vector[rns_layer]);
        
        fmpz_mod_mat_t mat1, mat2, result;
        fmpz_mod_mat_init(mat1, matrix1_rows, matrix1_cols, modulus);
        fmpz_mod_mat_init(mat2, matrix2_rows, matrix2_cols, modulus);
        fmpz_mod_mat_init(result, matrix1_rows, matrix2_cols, modulus);
        
        // 填充矩阵1
        for (size_t i = 0; i < matrix1_rows; i++) {
            for (size_t j = 0; j < matrix1_cols; j++) {
                fmpz_t val;
                fmpz_init_set_ui(val, coeff_matrix_1[rns_layer][i][j]);
                fmpz_mod_mat_set_entry(mat1, i, j, val);
                fmpz_clear(val);
            }
        }
        
        // 填充矩阵2
        for (size_t i = 0; i < matrix2_rows; i++) {
            for (size_t j = 0; j < matrix2_cols; j++) {
                fmpz_t val;
                fmpz_init_set_ui(val, coeff_matrix_2[rns_layer][i][j]);
                fmpz_mod_mat_set_entry(mat2, i, j, val);
                fmpz_clear(val);
            }
        }
        
        // 执行矩阵乘法
        fmpz_mod_mat_mul(result, mat1, mat2);
        
        // 提取结果
        for (size_t i = 0; i < matrix1_rows; i++) {
            for (size_t j = 0; j < matrix2_cols; j++) {
                fmpz_t val;
                fmpz_init(val);
                fmpz_mod_mat_get_entry(val, result, i, j);
                result_matrix[rns_layer][i][j] = fmpz_get_ui(val);
                fmpz_clear(val);
            }
        }
        
        // 清理FLINT对象
        fmpz_mod_mat_clear(mat1);
        fmpz_mod_mat_clear(mat2);
        fmpz_mod_mat_clear(result);
        fmpz_clear(modulus);
    }

    cout << "\nFLINT RNS-RNS矩阵乘法计算完成！" << endl;
}

void matrix_add_plain_flint(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector)
{
    if (coeff_matrix_1.empty() || coeff_matrix_2.empty()) {
        cerr << "Empty matrices in addition" << endl;
        return;
    }

    size_t rns_layers = coeff_matrix_1.size();
    
    // 检查RNS层数是否一致
    if (coeff_matrix_2.size() != rns_layers) {
        cerr << "RNS layer count mismatch: matrix1 has " << rns_layers 
             << " layers, matrix2 has " << coeff_matrix_2.size() << " layers" << endl;
        return;
    }
    
    if (modulus_vector.size() != rns_layers) {
        cerr << "Modulus vector size mismatch: expected " << rns_layers 
             << ", got " << modulus_vector.size() << endl;
        return;
    }

    cout << "开始使用FLINT进行RNS-RNS矩阵加法计算..." << endl;
    cout << "RNS层数: " << rns_layers << endl;

    result_matrix.resize(rns_layers);

    for (size_t rns_layer = 0; rns_layer < rns_layers; rns_layer++) {
        cout << "\r处理RNS层: " << (rns_layer + 1) << "/" << rns_layers 
             << " [" << (rns_layer + 1) * 100 / rns_layers << "%]" << flush;

        size_t matrix1_rows = coeff_matrix_1[rns_layer].size();
        size_t matrix1_cols = coeff_matrix_1[rns_layer][0].size();
        size_t matrix2_rows = coeff_matrix_2[rns_layer].size();
        size_t matrix2_cols = coeff_matrix_2[rns_layer][0].size();
        
        // 检查矩阵维度是否匹配
        if (matrix1_rows != matrix2_rows || matrix1_cols != matrix2_cols) {
            cerr << "\nMatrix dimensions do not match for addition at RNS layer " << rns_layer 
                 << ": matrix1 is " << matrix1_rows << "x" << matrix1_cols 
                 << ", matrix2 is " << matrix2_rows << "x" << matrix2_cols << endl;
            return;
        }

        result_matrix[rns_layer].resize(matrix1_rows, vector<uint64_t>(matrix1_cols));

        // 使用FLINT进行矩阵加法
        fmpz_t modulus;
        fmpz_init_set_ui(modulus, modulus_vector[rns_layer]);
        
        fmpz_mod_mat_t mat1, mat2, result;
        fmpz_mod_mat_init(mat1, matrix1_rows, matrix1_cols, modulus);
        fmpz_mod_mat_init(mat2, matrix2_rows, matrix2_cols, modulus);
        fmpz_mod_mat_init(result, matrix1_rows, matrix1_cols, modulus);
        
        // 填充矩阵1
        for (size_t i = 0; i < matrix1_rows; i++) {
            for (size_t j = 0; j < matrix1_cols; j++) {
                fmpz_t val;
                fmpz_init_set_ui(val, coeff_matrix_1[rns_layer][i][j]);
                fmpz_mod_mat_set_entry(mat1, i, j, val);
                fmpz_clear(val);
            }
        }
        
        // 填充矩阵2
        for (size_t i = 0; i < matrix2_rows; i++) {
            for (size_t j = 0; j < matrix2_cols; j++) {
                fmpz_t val;
                fmpz_init_set_ui(val, coeff_matrix_2[rns_layer][i][j]);
                fmpz_mod_mat_set_entry(mat2, i, j, val);
                fmpz_clear(val);
            }
        }
        
        // 执行矩阵加法
        fmpz_mod_mat_add(result, mat1, mat2);
        
        // 提取结果
        for (size_t i = 0; i < matrix1_rows; i++) {
            for (size_t j = 0; j < matrix1_cols; j++) {
                fmpz_t val;
                fmpz_init(val);
                fmpz_mod_mat_get_entry(val, result, i, j);
                result_matrix[rns_layer][i][j] = fmpz_get_ui(val);
                fmpz_clear(val);
            }
        }
        
        // 清理FLINT对象
        fmpz_mod_mat_clear(mat1);
        fmpz_mod_mat_clear(mat2);
        fmpz_mod_mat_clear(result);
        fmpz_clear(modulus);
    }

    cout << "\nFLINT RNS-RNS矩阵加法计算完成！" << endl;
}

void ciphertext_matrix_transpose(
    const SEALContext& context,
    const GaloisKeys& galois_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& row_encrypted_matrix,
    vector<Ciphertext>& col_encrypted_matrix)
{
    // 获取参数
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t N = context_data->parms().poly_modulus_degree();
    size_t num_ciphertexts = row_encrypted_matrix.size();
    if (N != num_ciphertexts) {
        cerr << "CM-T: 输入密文数量与poly_modulus_degree不符" << endl;
        return;
    }
    cout << "[CM-T] N = " << N << endl;
    uint64_t N_inv_mod_q = modinv(N, context_data->parms().plain_modulus().value());

    // 计时变量
    auto start_time = chrono::high_resolution_clock::now();
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> inv_time(0);
    chrono::duration<double> shift_time(0);
    chrono::duration<double> add_time(0);
    chrono::duration<double> galois_time(0);

    // Step 1: 计算 tilde_m_t = sum_i m_i(X) * X^{i*(2t+1)^{-1}}
    vector<Ciphertext> tilde_m(N);
    for (size_t t = 0; t < N; t++) {
        cout << "\r[CM-T] 处理tilde_m_t: " << (t + 1) << "/" << N << " [" << (t + 1) * 100 / N << "%]" << flush;
        uint64_t two_t_plus_one = (2 * t + 1) % (2 * N);

        start_time = chrono::high_resolution_clock::now();
        uint64_t inv_two_t_plus_one = modinv(two_t_plus_one, 2 * N);
        end_time = chrono::high_resolution_clock::now();
        inv_time += chrono::duration<double>(end_time - start_time);

        bool first = true;
        for (size_t i = 0; i < N; i++) {
            // 计算 i * (2t+1)^{-1} mod 2N
            uint64_t shift = (i * inv_two_t_plus_one) % (2 * N);

            start_time = chrono::high_resolution_clock::now();
            // 使用 negacyclic_shift_poly_coeffmod 直接处理整个密文
            Ciphertext shifted_ct = row_encrypted_matrix[i];
            util::negacyclic_shift_poly_coeffmod(
                shifted_ct, shifted_ct.size(), shift, 
                context_data->parms().coeff_modulus(), shifted_ct
            );
            end_time = chrono::high_resolution_clock::now();
            shift_time += chrono::duration<double>(end_time - start_time);
            
            start_time = chrono::high_resolution_clock::now();
            if (first) {
                tilde_m[t] = shifted_ct;
                first = false;
            } else {
                evaluator.add_inplace(tilde_m[t], shifted_ct);
            }
            end_time = chrono::high_resolution_clock::now();
            add_time += chrono::duration<double>(end_time - start_time);
        }
    }
    cout << endl;
    // Step 2: 计算 bar_m_t = tilde_m_t(X^{2t+1})
    vector<Ciphertext> bar_m(N);
    for (size_t t = 0; t < N; t++) {
        cout << "\r[CM-T] 处理bar_m_t: " << (t + 1) << "/" << N << " [" << (t + 1) * 100 / N << "%]" << flush;
        bar_m[t] = tilde_m[t];
        start_time = chrono::high_resolution_clock::now();
        evaluator.apply_galois_inplace(bar_m[t], (2 * t + 1) % (2 * N), galois_keys);
        end_time = chrono::high_resolution_clock::now();
        galois_time += chrono::duration<double>(end_time - start_time);
    }
    cout << endl;
    // Step 3: 计算 m'_j = N^{-1} * sum_t bar_m_t(X) * X^{-j*(2t+1)}
    col_encrypted_matrix.resize(N);
    
    for (size_t j = 0; j < N; j++) {
        cout << "\r[CM-T] 处理m'_j: " << (j + 1) << "/" << N << " [" << (j + 1) * 100 / N << "%]" << flush;
        bool first = true;
        for (size_t t = 0; t < N; t++) {
            uint64_t two_t_plus_one = (2 * t + 1) % (2 * N);
            // 计算 -j*(2t+1) mod N
            uint64_t shift = (2 * N - (j * two_t_plus_one) % (2 * N)) % (2 * N);
            Ciphertext shifted_ct = bar_m[t];
            start_time = chrono::high_resolution_clock::now();
            util::negacyclic_shift_poly_coeffmod(
                shifted_ct, shifted_ct.size(), shift, 
                context_data->parms().coeff_modulus(), shifted_ct
            );
            end_time = chrono::high_resolution_clock::now();
            shift_time += chrono::duration<double>(end_time - start_time);
            
            start_time = chrono::high_resolution_clock::now();
            if (first) {
                col_encrypted_matrix[j] = shifted_ct;
                first = false;
            } else {
                evaluator.add_inplace(col_encrypted_matrix[j], shifted_ct);
            }
            end_time = chrono::high_resolution_clock::now();
            add_time += chrono::duration<double>(end_time - start_time);
        }
        util::multiply_poly_scalar_coeffmod(
            col_encrypted_matrix[j], col_encrypted_matrix[j].size(), N_inv_mod_q,
            context_data->parms().coeff_modulus(), col_encrypted_matrix[j]
        );
    }
    
    cout << "\n[CM-T] 时间统计:" << endl;
    cout << "  模逆元计算时间: " << inv_time.count() << " 秒" << endl;
    cout << "  多项式移位时间: " << shift_time.count() << " 秒" << endl;
    cout << "  密文加法时间: " << add_time.count() << " 秒" << endl;
    cout << "  Galois自同构时间: " << galois_time.count() << " 秒" << endl;
    cout << "  总时间: " << (inv_time + shift_time + add_time + galois_time).count() << " 秒" << endl;
    
    cout << "[CM-T] 密文矩阵转置完成！" << endl;
}

void ciphertext_matrix_transpose_tweak(
    const SEALContext& context,
    const GaloisKeys& galois_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& row_encrypted_matrix,
    vector<Ciphertext>& col_encrypted_matrix)
{
    size_t N = row_encrypted_matrix.size();
    if ((N & (N-1)) != 0) {
        cerr << "TWEAK算法要求N为2的幂" << endl;
        return;
    }
    // Step 1: aux = TWEAK(N, {X^i * ct_i})
    cout << "[TWEAK] 开始计算TWEAK(N, {X^i * ct_i})" << endl;
    vector<Ciphertext> x_ct(N);
    auto context_data = context.get_context_data(context.first_parms_id());
    for (size_t i = 0; i < N; ++i) {
        cout << "\r[TWEAK] 处理X^i * ct_i: " << (i + 1) << "/" << N << " [" << (i + 1) * 100 / N << "%]" << flush;
        x_ct[i] = row_encrypted_matrix[i];
        util::negacyclic_shift_poly_coeffmod(
            x_ct[i], x_ct[i].size(), i % N,
            context_data->parms().coeff_modulus(), x_ct[i]
        );
    }
    vector<Ciphertext> aux;
    ciphertext_tweak(context, x_ct, evaluator, aux);

    // Step 2: aux'_j = (N^{-1} mod Q) * aux_{((2j+1)^{-1} mod 2N-1)/2}
    //          aux'_j = Auto(aux'_j, 2j+1)
    cout << "\n[TWEAK] 开始计算aux'_j = (N^{-1} mod Q) * aux_{((2j+1)^{-1} mod 2N-1)/2}" << endl;
    vector<Ciphertext> aux_prime(N);
    uint64_t N_inv_mod_q = 0;
    // 取第一个modulus
    auto modulus = context_data->parms().coeff_modulus()[0];
    for (uint64_t k = 1; k < modulus.value(); ++k) {
        if ((k * N) % modulus.value() == 1) {
            N_inv_mod_q = k;
            break;
        }
    }
    for (size_t j = 0; j < N; ++j) {
        cout << "\r[TWEAK] 处理aux'_j: " << (j + 1) << "/" << N << " [" << (j + 1) * 100 / N << "%]" << flush;
        // 计算(2j+1)^{-1} mod 2N
        uint64_t two_j_plus_1 = (2 * j + 1) % (N);
        uint64_t inv_two_j_plus_1 = 1;
        for (uint64_t k = 1; k < N; ++k) {
            if ((k * two_j_plus_1) % (N) == 1) {
                inv_two_j_plus_1 = k;
                break;
            }
        }
        size_t idx = ((inv_two_j_plus_1 % (2*N)) / 2) % N;
        aux_prime[j] = aux[idx];
        // 乘以N^{-1} mod Q
        util::multiply_poly_scalar_coeffmod(
            aux_prime[j], aux_prime[j].size(), N_inv_mod_q,
            context_data->parms().coeff_modulus(), aux_prime[j]
        );
        // Auto(aux'_j, 2j+1)
        evaluator.apply_galois_inplace(aux_prime[j], 2*j+1, galois_keys);
    }
    // Step 3: ct'' = TWEAK(N, aux')
    cout << "\n[TWEAK] 开始计算ct'' = TWEAK(N, aux')" << endl;
    vector<Ciphertext> ct2;
    ciphertext_tweak(context, aux_prime, evaluator, ct2);
    // Step 4: ct'_j = ct2_{j} - X^{-j} * ct2_{N-j mod N}
    col_encrypted_matrix.resize(N);
    for (size_t j = 0; j < N; ++j) {
        cout << "\r[TWEAK] 处理ct'_j: " << (j + 1) << "/" << N << " [" << (j + 1) * 100 / N << "%]" << flush;
        col_encrypted_matrix[j] = ct2[j];
        Ciphertext neg_j_ct = ct2[(N-j)%N];
        util::negacyclic_shift_poly_coeffmod(
            neg_j_ct, neg_j_ct.size(), (N-j)%N,
            context_data->parms().coeff_modulus(), neg_j_ct
        );
        evaluator.sub_inplace(col_encrypted_matrix[j], neg_j_ct);
    }
}

void ciphertext_tweak(const SEALContext& context, const vector<Ciphertext>& ct, Evaluator& evaluator, vector<Ciphertext>& ct_out) {
    size_t n = ct.size();
    if (n == 1) {
        ct_out = ct;
        return;
    }
    // 递归分治
    size_t half = n / 2;
    vector<Ciphertext> ct_even(half), ct_odd(half);
    for (size_t j = 0; j < half; ++j) {
        ct_even[j] = ct[j];
        ct_odd[j] = ct[j + half];
    }
    vector<Ciphertext> aux_even, aux_odd;
    ciphertext_tweak(context, ct_even, evaluator, aux_even);
    ciphertext_tweak(context, ct_odd, evaluator, aux_odd);
    
    // 获取参数
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t N = context_data->parms().poly_modulus_degree();
    
    ct_out.resize(n);
    for (size_t j = 0; j < half; ++j) {
        // X^{N/2 * j} * aux_odd[j]
        Ciphertext shifted;
        shifted = aux_odd[j];
        util::negacyclic_shift_poly_coeffmod(
            shifted, shifted.size(), (N/2) * j % N,
            context_data->parms().coeff_modulus(), shifted
        );
        ct_out[j] = aux_even[j];
        evaluator.add_inplace(ct_out[j], shifted);
        
        // X^{N/2 * j} * aux_odd[j]
        Ciphertext shifted2 = aux_odd[j];
        util::negacyclic_shift_poly_coeffmod(
            shifted2, shifted2.size(), (N/2) * j % N,
            context_data->parms().coeff_modulus(), shifted2
        );
        ct_out[j + half] = aux_even[j];
        evaluator.sub_inplace(ct_out[j + half], shifted2);
    }
}

uint64_t modinv(uint64_t a, uint64_t modulus) {
    uint64_t result = 0;
    if (!seal::util::try_invert_uint_mod(a, modulus, result)) {
        throw std::invalid_argument("No modular inverse exists for the given input.");
    }
    return result;
}