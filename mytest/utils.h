#ifndef UTILS_H
#define UTILS_H

#include "seal/seal.h"
#include <vector>
#include <string>
#include <iostream>
#include <nlohmann/json.hpp>

using namespace seal;
using namespace std;
using json = nlohmann::json;

void print_parameters(const seal::SEALContext &context);

void print_matrix(const vector<vector<uint64_t>>& matrix, const string& name);

void print_ciphertext_info(const vector<Ciphertext>& ciphertexts, const string& name);

void check_matrix_equal(const vector<vector<uint64_t>>& matrix1, const vector<vector<uint64_t>>& matrix2);

// 新增的配置读取函数
json read_seal_config(const string& config_file = "seal_config.json", bool verbose = false);
vector<int> get_coeff_modulus_params(const json& config, size_t poly_modulus_degree);

#endif // UTILS_H