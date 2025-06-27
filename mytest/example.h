#ifndef EXAMPLE_H
#define EXAMPLE_H

#include "cpmm.h"
#include "ccmm.h"
#include <iostream>
#include <string>

using namespace std;

/**
 * 测试CPMM（Ciphertext-Plaintext Matrix Multiplication）功能
 * @return 成功返回0，失败返回1
 */
int test_cpmm();

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
 * 测试negacyclic_multiply_poly_mono_coeffmod的移位功能
 * @return 成功返回0，失败返回1
 */
int test_function();

/**
 * 打印菜单选项
 */
void print_menu();

#endif // EXAMPLE_H 