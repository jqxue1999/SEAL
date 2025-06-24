#ifndef UTILS_H
#define UTILS_H

#include "seal/seal.h"
#include <vector>
#include <string>
#include <iostream>

using namespace seal;
using namespace std;

void print_parameters(const seal::SEALContext &context);

void print_matrix(const vector<vector<uint64_t>>& matrix, const string& name);

void print_ciphertext_info(const vector<Ciphertext>& ciphertexts, const string& name);

void check_matrix_equal(const vector<vector<uint64_t>>& matrix1, const vector<vector<uint64_t>>& matrix2);

#endif // UTILS_H