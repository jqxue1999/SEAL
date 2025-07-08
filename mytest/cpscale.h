#ifndef CPSCALE_H
#define CPSCALE_H

#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include <vector>
#include <chrono>

using namespace seal;
using namespace std;

/**
 * 从单个密文中提取系数矩阵
 * 
 * @param context SEAL上下文
 * @param encrypted 单个密文
 * @param coeff_matrix_a 提取的a多项式系数矩阵 [RNS层][多项式系数]
 * @param coeff_matrix_b 提取的b多项式系数矩阵 [RNS层][多项式系数]
 * @param modulus_vector 每个RNS层的模数
 * @return 提取时间（秒）
 */
double extract_coefficients_from_single_ciphertext(
    const SEALContext& context,
    const Ciphertext& encrypted,
    vector<vector<uint64_t>>& coeff_matrix_a,
    vector<vector<uint64_t>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector);

/**
 * 从结果矩阵构建单个密文
 * 
 * @param context SEAL上下文
 * @param result_a_matrix 结果a矩阵 [RNS层][多项式系数]
 * @param result_b_matrix 结果b矩阵 [RNS层][多项式系数]
 * @param result_encrypted 输出的单个密文
 * @return 构建时间（秒）
 */
double build_single_ciphertext_from_result(
    const SEALContext& context,
    const vector<vector<uint64_t>>& result_a_matrix,
    const vector<vector<uint64_t>>& result_b_matrix,
    Ciphertext& result_encrypted);
    
/**
 * 对系数矩阵进行scale乘法（BLAS加速版本）
 * 
 * @param coeff_matrix_a a多项式系数矩阵 [RNS层][多项式系数]
 * @param coeff_matrix_b b多项式系数矩阵 [RNS层][多项式系数]
 * @param modulus_vector 每个RNS层的模数
 * @param scale 缩放因子
 * @return 乘法时间（秒）
 */
double scale_coefficients_blas(
    vector<vector<uint64_t>>& coeff_matrix_a,
    vector<vector<uint64_t>>& coeff_matrix_b,
    const vector<uint64_t>& modulus_vector,
    uint64_t scale);

/**
 * 直接对密文进行scale乘法（使用util::negacyclic_multiply_poly_mono_coeffmod）
 * 
 * @param context SEAL上下文
 * @param encrypted 输入密文
 * @param scale 缩放因子
 * @param pool 内存池
 * @return 乘法时间（秒）
 */
double scale_ciphertext_direct(
    const SEALContext& context,
    Ciphertext& encrypted,
    uint64_t scale,
    MemoryPoolHandle pool);

#endif // CPSCALE_H
