#ifndef BASELINE_H
#define BASELINE_H

#include "seal/seal.h"
#include <vector>

using namespace seal;
using namespace std;

void cvps_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const uint64_t scalar,
    const Ciphertext& encrypted_vector_coeff,
    Ciphertext& encrypted_vector_coeff_result
);

void cvpv_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const Ciphertext& encrypted_vector_coeff, 
    vector<Ciphertext>& encrypted_vector_coeff_result
);

void pvcm_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    Ciphertext& encrypted_matrix_coeff_result
);

void pmcm_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<vector<uint64_t>>& clear_matrix,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    vector<Ciphertext>& encrypted_matrix_coeff_result
);

#endif // BASELINE_H
