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

json read_seal_config(const string& config_file = "seal_config.json", bool verbose = false);

size_t get_user_poly_modulus_degree(const json& config);

vector<int> get_coeff_modulus_params(const json& config, size_t poly_modulus_degree);

void test_cpmm(bool verbose = false);

int test_ccmm();

int test_CMT();

int test_digits_cvps(int num_bits = 64);

int test_digits_cvpv(int num_bits = 64, bool check_correctness = false);

int test_digits_pvcm(int num_bits = 64, bool check_correctness = false);

int test_ciphertext_scale_multiplication(bool verbose = false);

int test_vector_vector_outer_multiplication(bool verbose = false);

int test_encrypted_matrix_times_clear_matrix(int num_bits = 64, bool check_correctness = false);

void print_menu();

#endif // EXAMPLE_H 