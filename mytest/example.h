#ifndef EXAMPLE_H
#define EXAMPLE_H

#include "cpmm.h"
#include "ccmm.h"
#include "digits.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

/**
 * 读取SEAL配置文件
 * @param config_file 配置文件路径
 * @param verbose 是否显示详细输出
 * @return 配置JSON对象，失败时返回空对象
 */
json read_seal_config(const string& config_file = "seal_config.json", bool verbose = false);

/**
 * 获取用户输入的多项式模数次数
 * @param config 配置对象
 * @return 用户选择的多项式模数次数
 */
size_t get_user_poly_modulus_degree(const json& config);

/**
 * 根据多项式模数次数获取对应的系数模数参数
 * @param config 配置对象
 * @param poly_modulus_degree 多项式模数次数
 * @return 系数模数参数向量
 */
vector<int> get_coeff_modulus_params(const json& config, size_t poly_modulus_degree);

/**
 * 测试CPMM（Ciphertext-Plaintext Matrix Multiplication）功能
 * @param verbose 是否显示详细输出
 * @return 成功返回0，失败返回1
 */
void test_cpmm(bool verbose = false);

/**
 * 测试CCMM（Ciphertext-Ciphertext Matrix Multiplication）功能
 * @return 成功返回0，失败返回1
 */
int test_ccmm();

/**
 * 测试CMT（Ciphertext Matrix Transpose）功能
 * @return 成功返回0，失败返回1
 */
int test_CMT();


/**
 * 测试通用向量乘法功能
 * @return 成功返回0，失败返回1
 */
int test_general_multiplication(int num_bits = 64);

/**
 * 测试密文scale乘法功能
 * @param verbose 是否显示详细输出
 * @return 成功返回0，失败返回1
 */
int test_ciphertext_scale_multiplication(bool verbose = false);

/**
 * 测试明文vector外乘密文vector功能
 * @param verbose 是否显示详细输出
 * @return 成功返回0，失败返回1
 */
int test_vector_vector_outer_multiplication(bool verbose = false);

/**
 * 打印菜单选项
 */
void print_menu();

#endif // EXAMPLE_H 