#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include "cpmm.h"
#include "coeff_modulus_config.h"
#include <vector>
#include <iostream>

using namespace seal;
using namespace std;

void cvps_coefficients_seal(const SEALContext& context,
    const uint64_t scalar,
    const Ciphertext& encrypted_vector_coeff,
    Ciphertext& encrypted_vector_coeff_resul
);

void cvpv_coefficients_seal(
    const SEALContext& context,
    const vector<uint64_t>& clear_vector,
    const Ciphertext& encrypted_vector_coeff,
    vector<Ciphertext>& encrypted_vector_coeff_resul
);

void pvcm_coefficients_seal(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    Ciphertext& encrypted_matrix_coeff_result
);

void pmcm_coefficients_seal(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<vector<uint64_t>>& clear_matrix,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    vector<Ciphertext>& encrypted_matrix_coeff_result
);

#endif // COEFFICIENTS_H
