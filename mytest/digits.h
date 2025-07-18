#ifndef DIGITS_H
#define DIGITS_H

#include "seal/seal.h"
#include "seal/modulus.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/common.h"
#include "seal/util/polycore.h"
#include "seal/util/scalingvariant.h"
#include "seal/util/uintarith.h"
#include "seal/util/uintcore.h"
#include "seal/util/uintarithsmallmod.h"
#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "cpmm.h"
#include <chrono>
#include <cmath>

using namespace seal;
using namespace std;

/**
 * 将整数向量分解为位级向量
 * 每个位级向量包含所有输入整数的对应位
 * 
 * @param input_vector 输入的整数向量
 * @param bit_vectors 输出的位级向量数组 [64个向量]
 */
void decompose_to_bit_vectors(
    const vector<uint64_t>& input_vector,
    vector<vector<uint64_t>>& bit_vectors,
    int num_bits = 64);

/**
 * 将位级向量重新组合为整数向量
 * 
 * @param bit_vectors 位级向量数组 [64个向量]
 * @param plain_modulus_value 明文模数
 * @param output_vector 输出的整数向量
 */
void compose_from_bit_vectors(
    const vector<vector<uint64_t>>& bit_vectors,
    vector<uint64_t>& output_vector,
    int num_bits = 64,
    uint64_t plain_modulus_value = 0);

/**
 * 加密位级向量
 * 
 * @param context SEAL上下文
 * @param encryptor 加密器
 * @param bit_vectors 位级向量数组
 * @param encrypted_bit_vectors 加密后的位级向量数组（64个Ciphertext）
 */
void encrypt_bit_vectors(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<vector<uint64_t>>& bit_vectors,
    vector<Ciphertext>& encrypted_bit_vectors,
    int num_bits = 64);

/**
 * 解密位级向量
 * 
 * @param context SEAL上下文
 * @param decryptor 解密器
 * @param encrypted_bit_vectors 加密的位级向量数组（64个Ciphertext）
 * @param decrypted_bit_vectors 解密后的位级向量数组
 */
void decrypt_bit_vectors(
    const SEALContext& context,
    Decryptor& decryptor,
    const vector<Ciphertext>& encrypted_bit_vectors,
    vector<vector<uint64_t>>& decrypted_bit_vectors,
    int num_bits = 64);

/**
 * 对加密的位级向量执行乘以2的幂次方操作
 * 通过交换向量位置来实现
 * 
 * @param context SEAL上下文
 * @param encryptor 加密器（用于创建全零密文）
 * @param encrypted_bit_vectors 加密的位级向量数组（64个Ciphertext）
 * @param power 2的幂次方 (0-63)
 * @param result_vectors 结果向量数组（64个Ciphertext）
 */
void multiply_by_power_of_2(
    const SEALContext& context,
    Encryptor& encryptor,
    const vector<Ciphertext>& encrypted_bit_vectors,
    int power,
    vector<Ciphertext>& result_vectors,
    int num_bits = 64);

/**
 * 将整数分解为2的幂次方之和
 * 
 * @param value 要分解的整数
 * @return 包含幂次方的向量
 */
vector<int> decompose_to_powers_of_2(uint64_t value);

/**
 * 初始化全零密文
 * 
 * @param context SEAL上下文
 * @param encryptor 加密器
 * @param zero_ciphertext 全零密文
 */
double initialize_zero_ciphertext(
    const SEALContext& context,
    Encryptor& encryptor,
    Ciphertext& zero_ciphertext);

void cvps_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const Ciphertext& zero_ciphertext,
    const vector<Ciphertext>& encrypted_bit_vectors,
    uint64_t multiplier,
    vector<Ciphertext>& result_vectors,
    int num_bits);

/**
 * 验证通用向量乘法的正确性
 * 
 * @param input_vector 输入向量
 * @param multiplier 乘数
 * @param output_vector 输出向量
 * @return 验证成功返回true，失败返回false
 */
bool verify_general_multiplication(
    const vector<uint64_t>& input_vector,
    uint64_t multiplier,
    const vector<uint64_t>& output_vector);

void cvpv_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    Decryptor& decryptor,
    int num_bits,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& bit_vectors_ciphertext,
    vector<vector<Ciphertext>>& outer_product_results);

// 计算明文向量和密文矩阵的乘积，结果为密文bit向量
double pvcm_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<vector<Ciphertext>>& encrypted_matrix,
    vector<Ciphertext>& result,
    int num_bits);

// 计算密文矩阵（按行密文bit向量）与明文矩阵的乘积，结果为密文bit向量矩阵
double pmcm_digits(
    const SEALContext& context,
    Encryptor& encryptor,
    Evaluator& evaluator,
    const vector<vector<Ciphertext>>& encrypted_matrix,
    const vector<vector<uint64_t>>& clear_matrix,
    vector<vector<Ciphertext>>& result,
    int num_bits);

#endif // DIGITS_H 