#ifndef CCMM_H
#define CCMM_H

#include "cpmm.h"


using namespace seal;
using namespace std; 

/**
 * 执行密文矩阵与密文矩阵的乘法
 * 
 * @param context SEAL上下文
 * @param relin_keys 重线性化密钥
 * @param evaluator 评估器
 * @param encrypted_matrix_1_row 按行编码的密文矩阵
 * @param encrypted_matrix_2_col 按列编码的密文矩阵
 * @param result_matrix 结果密文矩阵
 * @return 成功返回true，失败返回false
 */
bool ciphertext_ciphertext_matrix_multiply(
    const SEALContext& context,
    const RelinKeys& relin_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& encrypted_matrix_1_row,
    const vector<Ciphertext>& encrypted_matrix_2_col,
    vector<Ciphertext>& result_matrix);

/**
 * 对两个RNS系数矩阵进行乘法运算
 * 
 * @param coeff_matrix_1 第一个系数矩阵 [RNS层][行][列]
 * @param coeff_matrix_2 第二个系数矩阵 [RNS层][行][列]
 * @param result_matrix 结果矩阵 [RNS层][行][列]
 * @param modulus_vector 每个RNS层的模数
 */
void RNS_RNS_multiply(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector);

/**
 * 使用BLAS对两个RNS系数矩阵进行加法运算
 * 
 * @param coeff_matrix_1 第一个系数矩阵 [RNS层][行][列]
 * @param coeff_matrix_2 第二个系数矩阵 [RNS层][行][列]
 * @param result_matrix 结果矩阵 [RNS层][行][列]
 * @param modulus_vector 每个RNS层的模数
 */
void matrix_add_plain_blas(
    const vector<vector<vector<uint64_t>>>& coeff_matrix_1,
    const vector<vector<vector<uint64_t>>>& coeff_matrix_2,
    vector<vector<vector<uint64_t>>>& result_matrix,
    const vector<uint64_t>& modulus_vector);

/**
 * 从三个结果矩阵构建密文向量（3个多项式）
 * 
 * @param context SEAL上下文
 * @param result_c0_matrix 结果c0矩阵 [RNS层][行][列]
 * @param result_c1_matrix 结果c1矩阵 [RNS层][行][列]
 * @param result_c2_matrix 结果c2矩阵 [RNS层][行][列]
 * @param result_matrix 输出的密文向量
 */
void build_ciphertexts_from_result_matrices_3(
    const SEALContext& context,
    const vector<vector<vector<uint64_t>>>& result_c0_matrix,
    const vector<vector<vector<uint64_t>>>& result_c1_matrix,
    const vector<vector<vector<uint64_t>>>& result_c2_matrix,
    vector<Ciphertext>& result_matrix);

/**
 * 密文矩阵转置（CM-T算法）：将按行加密的密文向量转为按列加密的密文向量
 * @param context SEAL上下文
 * @param galois_keys Galois密钥
 * @param evaluator 评估器
 * @param row_encrypted_matrix 按行加密的密文向量
 * @param col_encrypted_matrix 输出：按列加密的密文向量
 */
void ciphertext_matrix_transpose(
    const SEALContext& context,
    const GaloisKeys& galois_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& row_encrypted_matrix,
    vector<Ciphertext>& col_encrypted_matrix);

/**
 * TWEAK算法：高效计算 { sum_i X^{2ij} * ct_i }_{j∈[N] }
 * @param context SEAL上下文
 * @param ct 输入密文向量
 * @param evaluator 评估器
 * @param ct_out 输出密文向量
 */
void ciphertext_tweak(const SEALContext& context, const vector<Ciphertext>& ct, Evaluator& evaluator, vector<Ciphertext>& ct_out);

/**
 * 基于TWEAK算法的密文矩阵转置（Algorithm 2）
 * @param context SEAL上下文
 * @param galois_keys Galois密钥
 * @param evaluator 评估器
 * @param row_encrypted_matrix 按行加密的密文向量
 * @param col_encrypted_matrix 输出：按列加密的密文向量
 */
void ciphertext_matrix_transpose_tweak(
    const SEALContext& context,
    const GaloisKeys& galois_keys,
    Evaluator& evaluator,
    const vector<Ciphertext>& row_encrypted_matrix,
    vector<Ciphertext>& col_encrypted_matrix);

/**
 * 计算模逆元：a^(-1) mod modulus
 * @param a 输入值
 * @param modulus 模数
 * @return a在模modulus下的逆元
 * @throws std::invalid_argument 如果逆元不存在
 */
uint64_t modinv(uint64_t a, uint64_t modulus);

#endif // CCMM_H