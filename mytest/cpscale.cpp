#include "cpscale.h"
#include "seal/seal.h"
#include "seal/util/uintarith.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cblas.h>
#include <omp.h>

using namespace seal;
using namespace std;
using namespace seal::util;

double extract_coefficients_from_single_ciphertext(
    const SEALContext& context,
    const Ciphertext& encrypted,
    vector<vector<uint64_t>>& coeff_matrix_a,
    vector<vector<uint64_t>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector)
{
    auto start_time = chrono::high_resolution_clock::now();
    
    auto context_data = context.get_context_data(encrypted.parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = encrypted.coeff_modulus_size();
    
    cout << "开始提取系数..." << endl;
    cout << "RNS层数: " << coeff_modulus_size << endl;
    
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
    
    // 初始化系数矩阵 - 第一维是RNS层，第二维是该密文在该RNS层上的系数
    coeff_matrix_a.resize(coeff_modulus_size);
    coeff_matrix_b.resize(coeff_modulus_size);
    
    for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
        coeff_matrix_a[rns_layer].resize(poly_modulus_degree);
        coeff_matrix_b[rns_layer].resize(poly_modulus_degree);

        // 提取a多项式的系数（第一个多项式）
        if (encrypted.size() > 0) {
            const uint64_t* coeffs_a = encrypted.data(0);
            // 提取第rns_layer层的系数
            for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
                coeff_matrix_a[rns_layer][coeff_idx] = coeffs_a[coeff_idx + rns_layer * poly_modulus_degree];
        }
        
        // 提取b多项式的系数（第二个多项式）
        if (encrypted.size() > 1) {
            const uint64_t* coeffs_b = encrypted.data(1);
            // 提取第rns_layer层的系数
            for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
                coeff_matrix_b[rns_layer][coeff_idx] = coeffs_b[coeff_idx + rns_layer * poly_modulus_degree];
        }
        
        // 显示进度
        cout << "\r提取系数进度: RNS层 " << (rns_layer + 1) << "/" << coeff_modulus_size 
             << " [" << (rns_layer + 1) * 100 / coeff_modulus_size << "%]" << flush;
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}


double build_single_ciphertext_from_result(
    const SEALContext& context,
    const vector<vector<uint64_t>>& result_a_matrix,
    const vector<vector<uint64_t>>& result_b_matrix,
    Ciphertext& result_encrypted)
{
    auto start_time = chrono::high_resolution_clock::now();
    
    if (result_a_matrix.empty() || result_b_matrix.empty()) {
        cerr << "Empty result matrices" << endl;
        return 0;
    }
    
    auto context_data = context.get_context_data(context.first_parms_id());
    size_t poly_modulus_degree = context_data->parms().poly_modulus_degree();
    size_t coeff_modulus_size = result_a_matrix.size();
    
    // 检查维度
    if (result_b_matrix.size() != coeff_modulus_size) {
        cerr << "RNS layer count mismatch in result matrices" << endl;
        return 0;
    }
    
    cout << "开始构建单个密文..." << endl;
    cout << "RNS层数: " << coeff_modulus_size << endl;
    
    // 创建新的密文
    result_encrypted.resize(context, 2); // 2个多项式：a和b
    
    // 设置a多项式的系数
    if (result_encrypted.size() > 0) {
        uint64_t* coeffs_a = result_encrypted.data(0);
        for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
            for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                coeffs_a[coeff_idx + rns_layer * poly_modulus_degree] = 
                    result_a_matrix[rns_layer][coeff_idx];
            }
        }
    }
    
    // 设置b多项式的系数
    if (result_encrypted.size() > 1) {
        uint64_t* coeffs_b = result_encrypted.data(1);
        for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
            for (size_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++) {
                coeffs_b[coeff_idx + rns_layer * poly_modulus_degree] = 
                    result_b_matrix[rns_layer][coeff_idx];
            }
        }
    }
    
    cout << "单个密文构建完成！" << endl;

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}

double scale_coefficients_blas(
    vector<vector<uint64_t>>& coeff_matrix_a,
    vector<vector<uint64_t>>& coeff_matrix_b,
    const vector<uint64_t>& modulus_vector,
    uint64_t scale)
{
    auto start_time = chrono::high_resolution_clock::now();
    
    if (coeff_matrix_a.empty() || coeff_matrix_b.empty()) {
        cerr << "Empty coefficient matrices" << endl;
        return 0;
    }
    
    size_t coeff_modulus_size = coeff_matrix_a.size();
    size_t coeff_count = coeff_matrix_a[0].size();
    
    cout << "开始对系数进行scale乘法（BLAS版本），scale=" << scale << endl;
    
    // 统计数据类型转换时间
    chrono::duration<double> convert_time(0);
    chrono::duration<double> blas_time(0);
    chrono::duration<double> mod_time(0);
    
    // 对每个RNS层的系数进行scale乘法
    #pragma omp parallel for
    for (size_t rns_layer = 0; rns_layer < coeff_modulus_size; rns_layer++) {
        uint64_t modulus = modulus_vector[rns_layer];
        Modulus mod_obj(modulus);
        
        // 使用BLAS加速：将uint64_t转换为double，进行BLAS乘法，再转换回uint64_t
        vector<double> coeff_a_double(coeff_count);
        vector<double> coeff_b_double(coeff_count);
        
        // 转换为double数组（使用OpenMP并行化）
        auto convert_start = chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < coeff_count; i++) {
            coeff_a_double[i] = static_cast<double>(coeff_matrix_a[rns_layer][i]);
            coeff_b_double[i] = static_cast<double>(coeff_matrix_b[rns_layer][i]);
        }
        auto convert_end = chrono::high_resolution_clock::now();
        convert_time += chrono::duration_cast<chrono::duration<double>>(convert_end - convert_start);
        
        // 使用BLAS的dscal进行向量缩放（相当于乘以scale）
        auto blas_start = chrono::high_resolution_clock::now();
        cblas_dscal(coeff_count, static_cast<double>(scale), coeff_a_double.data(), 1);
        cblas_dscal(coeff_count, static_cast<double>(scale), coeff_b_double.data(), 1);
        auto blas_end = chrono::high_resolution_clock::now();
        blas_time += chrono::duration_cast<chrono::duration<double>>(blas_end - blas_start);
        
        // 转换回uint64_t并应用模运算（使用OpenMP并行化）
        auto mod_start = chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < coeff_count; i++) {
            coeff_matrix_a[rns_layer][i] = util::barrett_reduce_64(
                static_cast<uint64_t>(coeff_a_double[i]), mod_obj);
            coeff_matrix_b[rns_layer][i] = util::barrett_reduce_64(
                static_cast<uint64_t>(coeff_b_double[i]), mod_obj);
        }
        auto mod_end = chrono::high_resolution_clock::now();
        mod_time += chrono::duration_cast<chrono::duration<double>>(mod_end - mod_start);
        
        // 显示进度
        cout << "\r系数乘法进度: RNS层 " << (rns_layer + 1) << "/" << coeff_modulus_size 
             << " [" << (rns_layer + 1) * 100 / coeff_modulus_size << "%]" << flush;
    }
    
    cout << "\n系数乘法完成！" << endl;
    cout << "BLAS版本时间统计:" << endl;
    cout << "  数据类型转换时间: " << convert_time.count() << " seconds" << endl;
    cout << "  BLAS运算时间: " << blas_time.count() << " seconds" << endl;
    cout << "  模运算时间: " << mod_time.count() << " seconds" << endl;

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}

double scale_ciphertext_direct(
    const SEALContext& context,
    Ciphertext& encrypted,
    uint64_t scale,
    MemoryPoolHandle pool)
{
    auto start_time = chrono::high_resolution_clock::now();
    
    cout << "开始直接对密文进行scale乘法，scale=" << scale << endl;
    
    // 获取加密参数
    auto context_data = context.get_context_data(encrypted.parms_id());
    const auto& parms = context_data->parms();
    
    // 使用util::negacyclic_multiply_poly_mono_coeffmod直接对密文进行scale
    util::negacyclic_multiply_poly_mono_coeffmod(
        encrypted, 
        encrypted.size(), 
        scale, 
        0, 
        parms.coeff_modulus(), 
        encrypted, 
        pool);
    
    cout << "直接scale乘法完成！" << endl;

    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);

    return total_time.count();
}