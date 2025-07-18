#include "coefficients_seal.h"

void cvps_coefficients_seal(
    const SEALContext& context,
    const uint64_t scalar,
    const Ciphertext& encrypted_vector_coeff,
    Ciphertext& encrypted_vector_coeff_result
) {    
    encrypted_vector_coeff_result = encrypted_vector_coeff;
    EncryptionParameters parms = context.first_context_data()->parms();
    
    MemoryPoolHandle pool = MemoryManager::GetPool();
    util::negacyclic_multiply_poly_mono_coeffmod(
        encrypted_vector_coeff, encrypted_vector_coeff.size(), scalar, 0, parms.coeff_modulus(), encrypted_vector_coeff_result, pool
    );
}

void cvpv_coefficients_seal(
    const SEALContext& context,
    const vector<uint64_t>& clear_vector,
    const Ciphertext& encrypted_vector_coeff, 
    vector<Ciphertext>& encrypted_vector_coeff_result
) {
    encrypted_vector_coeff_result.resize(clear_vector.size());
    #pragma omp parallel for
    for (size_t i = 0; i < clear_vector.size(); ++i) {
        cvps_coefficients_seal(context, clear_vector[i], encrypted_vector_coeff, encrypted_vector_coeff_result[i]);
    }
}

void pvcm_coefficients_seal(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    Ciphertext& encrypted_matrix_coeff_result
) {
    cvps_coefficients_seal(context, clear_vector[0], encrypted_matrix_coeff[0], encrypted_matrix_coeff_result);

    #pragma omp parallel for
    for (size_t i = 1; i < clear_vector.size(); i++) {
        Ciphertext partial;
        cvps_coefficients_seal(context, clear_vector[i], encrypted_matrix_coeff[i], partial);
        evaluator.add_inplace(encrypted_matrix_coeff_result, partial);
    }
}

void pmcm_coefficients_seal(
    const SEALContext& context,
    const Evaluator& evaluator,
    const vector<vector<uint64_t>>& clear_matrix,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    vector<Ciphertext>& encrypted_matrix_coeff_result
) {
    encrypted_matrix_coeff_result.resize(clear_matrix.size());
    #pragma omp parallel for
    for (size_t i = 0; i < clear_matrix.size(); i++)
        pvcm_coefficients_seal(context, evaluator, clear_matrix[i], encrypted_matrix_coeff, encrypted_matrix_coeff_result[i]);
}

