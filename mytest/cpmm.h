#ifndef CPMM_H
#define CPMM_H

#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/common.h"
#include "seal/util/polycore.h"
#include "seal/util/scalingvariant.h"
#include "seal/util/uintarith.h"
#include "seal/util/uintcore.h"
#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cblas.h>
#include <cassert>

using namespace seal;
using namespace std;

/**
 * 执行密文矩阵与明文矩阵的乘法
 * 
 * @param context SEAL上下文
 * @param encrypted_matrix 1×d的密文向量，每个密文编码一个矩阵行
 * @param plain_matrix d×d的明文矩阵
 * @param result_matrix 1×d的结果密文向量
 * @return 成功返回true，失败返回false
 */
bool ciphertext_plaintext_matrix_multiply(
    const SEALContext& context,
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<Ciphertext>& encrypted_matrix,
    vector<Ciphertext>& result_matrix);

/**
 * 从密文向量中提取系数矩阵
 * 
 * @param context SEAL上下文
 * @param encrypted_vector 1×d的密文向量
 * @param coeff_matrix_a 提取的a多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param coeff_matrix_b 提取的b多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param modulus_vector 每个RNS层的模数
 */
void extract_coefficients_from_ciphertext_vector(
    const SEALContext& context,
    const vector<Ciphertext>& encrypted_vector,
    vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<uint64_t>& modulus_vector);

/**
 * 对所有RNS层执行矩阵乘法
 * 
 * @param plain_matrix d×d明文矩阵
 * @param coeff_matrix_a a多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param coeff_matrix_b b多项式系数矩阵 [RNS层][密文索引][多项式系数]
 * @param result_a_matrix 结果a矩阵 [RNS层][行][列]
 * @param result_b_matrix 结果b矩阵 [RNS层][行][列]
 * @param modulus_vector 每个RNS层的模数
 */
void matrix_multiply_for_all_rns_layers(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<vector<vector<uint64_t>>>& result_a_matrix,
    vector<vector<vector<uint64_t>>>& result_b_matrix,
    const vector<uint64_t>& modulus_vector);

void matrix_multiply_for_all_rns_layers_blas(
    const vector<vector<uint64_t>>& plain_matrix,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_a,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_b,
    vector<vector<vector<uint64_t>>>& result_a_matrix,
    vector<vector<vector<uint64_t>>>& result_b_matrix,
    const vector<uint64_t>& modulus_vector);

/**
 * 从结果矩阵构建密文向量
 * 
 * @param context SEAL上下文
 * @param result_a_matrix 结果a矩阵 [RNS层][行][列]
 * @param result_b_matrix 结果b矩阵 [RNS层][行][列]
 * @param result_matrix 输出的密文向量
 */
void build_ciphertexts_from_result_matrices(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_a_matrix,
    const vector<vector<vector<uint64_t>>>& result_b_matrix,
    vector<Ciphertext>& result_matrix);

/**
 * 将矩阵行编码到明文多项式系数上
 * 
 * @param matrix_row 矩阵行
 * @param context SEAL上下文
 * @param plaintext 输出的明文多项式
 */
void encode_matrix_row_to_plaintext(
    const vector<uint64_t>& matrix_row,
    const SEALContext& context,
    Plaintext& plaintext);

/**
 * 从明文多项式系数解码回矩阵行
 * 
 * @param plaintext 明文多项式
 * @param context SEAL上下文
 * @param matrix_row 输出的矩阵行
 */
void decode_plaintext_to_matrix_row(
    const Plaintext& plaintext,
    const SEALContext& context,
    vector<uint64_t>& matrix_row);

/**
 * 将明文矩阵按行编码并加密为密文向量
 * 
 * @param context SEAL上下文
 * @param encryptor 加密器
 * @param plain_matrix 明文矩阵
 * @param encrypted_matrix 输出的密文向量
 */
void encrypt_matrix_rows(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& plain_matrix,
    vector<Ciphertext>& encrypted_matrix);

/**
 * 批量解密密文向量为明文矩阵
 *
 * @param context SEAL上下文
 * @param decryptor 解密器
 * @param encrypted_matrix 密文向量
 * @param plain_matrix 输出的明文矩阵
 */
void decrypt_ciphertexts_to_matrix(
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_matrix,
    vector<vector<uint64_t>>& plain_matrix);

/**
 * 计算明文矩阵乘法 C = A * B
 * @param A 左矩阵
 * @param B 右矩阵
 * @param C 结果矩阵
 * @param modulus 若非0则对结果取模
 */
void matrix_multiply_plain(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus = 0);

void matrix_multiply_plain_blas(
    const vector<vector<uint64_t>>& A,
    const vector<vector<uint64_t>>& B,
    vector<vector<uint64_t>>& C,
    uint64_t modulus);

#endif // CPMM_H
