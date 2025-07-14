#ifndef CPMM_H
#define CPMM_H

#include "seal/seal.h"
#include "seal/modulus.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/common.h"
#include "seal/util/polycore.h"
#include "seal/util/scalingvariant.h"
#include "seal/util/uintarith.h"
#include "seal/util/uintcore.h"
#include "seal/util/uintarithsmallmod.h"
#include "utils.h"
#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cblas.h>
#include <cassert>

// AVX相关头文件
#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <omp.h>

using namespace seal;
using namespace std;

/**
 * 执行密文矩阵与明文矩阵的乘法
 * 
 * @param context SEAL上下文
 * @param plain_matrix 明文矩阵
 * @param encrypted_matrix 密文矩阵
 * @param result_matrix 结果密文矩阵
 * @param verbose 是否显示详细输出
 * @return 时间向量 [extract_time, multiply_time, repack_time]
 */
vector<double> ciphertext_plaintext_matrix_multiply(
    const SEALContext& context,
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<Ciphertext>& encrypted_matrix,
    vector<Ciphertext>& result_matrix,
    bool verbose = false);

/**
 * 从密文向量中提取系数矩阵
 * 
 * @param context SEAL上下文
 * @param encrypted_vector 1×d的密文向量
 * @param coeff_matrix_a 提取的a多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param coeff_matrix_b 提取的b多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param modulus_vector 每个RNS层的模数
 * @param verbose 是否显示详细输出
 * @return 提取时间
 */
double extract_coefficients_from_ciphertext_vector(
    const SEALContext& context,
    const vector<Ciphertext>& encrypted_vector,
    vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector,
    bool verbose = false);

/**
 * 从结果矩阵构建密文向量
 * 
 * @param context SEAL上下文
 * @param result_a_matrix 结果a矩阵 [RNS层][行][列]
 * @param result_b_matrix 结果b矩阵 [RNS层][行][列]
 * @param result_matrix 输出的密文向量
 * @param verbose 是否显示详细输出
 * @return 构建时间
 */
double build_ciphertexts_from_result_matrices_2(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_a_matrix,
    const vector<vector<vector<uint64_t>>>& result_b_matrix,
    vector<Ciphertext>& result_matrix,
    bool verbose = false);

/**
 * 将向量编码到明文多项式系数上
 * 
 * @param vector_data 向量数据
 * @param context SEAL上下文
 * @param plaintext 输出的明文多项式
 */
void encode_vector_to_plaintext(
    const vector<uint64_t>& vector_data,
    const SEALContext& context,
    Plaintext& plaintext);

/**
 * 从明文多项式系数解码回向量
 * 
 * @param plaintext 明文多项式
 * @param context SEAL上下文
 * @param vector_data 输出的向量数据
 */
void decode_plaintext_to_vector(
    const Plaintext& plaintext,
    const SEALContext& context,
    vector<uint64_t>& vector_data);

/**
 * 将明文矩阵按行或列编码并加密为密文向量
 * 
 * @param is_row 如果为true则按行编码，否则按列编码
 * @param context SEAL上下文
 * @param encryptor 加密器
 * @param plain_matrix 明文矩阵
 * @param encrypted_matrix 输出的密文向量
 */
void encrypt_matrix(
    const bool is_row,
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& plain_matrix,
    vector<Ciphertext>& encrypted_matrix);

/**
 * 批量解密密文向量为明文矩阵（按行或列）
 *
 * @param is_row 如果为true则按行解码，否则按列解码
 * @param context SEAL上下文
 * @param decryptor 解密器
 * @param encrypted_matrix 密文向量
 * @param plain_matrix 输出的明文矩阵
 */
void decrypt_ciphertexts_to_matrix(
    const bool is_row,
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_matrix,
    vector<vector<uint64_t>>& plain_matrix);

/**
 * 计算明文矩阵乘法 C = A * B (带模运算)
 * @param A 左矩阵
 * @param B 右矩阵
 * @param C 结果矩阵
 * @param modulus 模数
 * @param verbose 是否显示详细输出
 * @return 计算时间
 */
double matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus,
    bool verbose = false);

/**
 * 计算明文矩阵乘法 C = A * B (不带模运算)
 * @param A 左矩阵
 * @param B 右矩阵
 * @param C 结果矩阵
 * @param verbose 是否显示详细输出
 * @return 计算时间
 */
double matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    bool verbose = false);

/**
 * 对单个系数矩阵进行RNS层矩阵乘法
 * 
 * @param plain_matrix 明文矩阵
 * @param coeff_matrix 系数矩阵 [RNS层][行][列]
 * @param result_matrix 结果矩阵 [RNS层][行][列]
 * @param modulus_vector 每个RNS层的模数
 * @param verbose 是否显示详细输出
 * @return 计算时间
 */
double Normal_RNS_multiply(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector,
    bool verbose = false);

#endif // CPMM_H
