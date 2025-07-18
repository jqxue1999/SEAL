#include "baseline.h"

void cvps_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const uint64_t scalar,
    const Ciphertext& encrypted_vector_coeff,
    Ciphertext& encrypted_vector_coeff_result
) {
    vector<uint64_t> scalar_vector(context.first_context_data()->parms().poly_modulus_degree(), scalar);
    Plaintext plain_scalar;
    batch_encoder.encode(scalar_vector, plain_scalar);

    evaluator.multiply_plain(encrypted_vector_coeff, plain_scalar, encrypted_vector_coeff_result);
}

void cvpv_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const Ciphertext& encrypted_vector_coeff, 
    vector<Ciphertext>& encrypted_vector_coeff_result
) {
    encrypted_vector_coeff_result.resize(clear_vector.size());
    // #pragma omp parallel for
    for (size_t i = 0; i < clear_vector.size(); ++i) {
        cvps_baseline(context, batch_encoder, evaluator, clear_vector[i], encrypted_vector_coeff, encrypted_vector_coeff_result[i]);
    }
}

void pvcm_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<uint64_t>& clear_vector,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    Ciphertext& encrypted_matrix_coeff_result
) {
    cvps_baseline(context, batch_encoder, evaluator, clear_vector[0], encrypted_matrix_coeff[0], encrypted_matrix_coeff_result);

    // #pragma omp parallel for
    for (size_t i = 1; i < clear_vector.size(); i++) {
        Ciphertext partial;
        cvps_baseline(context, batch_encoder, evaluator, clear_vector[i], encrypted_matrix_coeff[i], partial);
        evaluator.add_inplace(encrypted_matrix_coeff_result, partial);
    }
}

void pmcm_baseline(
    const SEALContext& context,
    const BatchEncoder& batch_encoder,
    const Evaluator& evaluator,
    const vector<vector<uint64_t>>& clear_matrix,
    const vector<Ciphertext>& encrypted_matrix_coeff,
    vector<Ciphertext>& encrypted_matrix_coeff_result
) {
    encrypted_matrix_coeff_result.resize(clear_matrix.size());
    // #pragma omp parallel for
    for (size_t i = 0; i < clear_matrix.size(); i++)
        pvcm_baseline(context, batch_encoder, evaluator, clear_matrix[i], encrypted_matrix_coeff, encrypted_matrix_coeff_result[i]);
}