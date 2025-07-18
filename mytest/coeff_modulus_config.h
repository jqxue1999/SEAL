#ifndef COEFF_MODULUS_CONFIG_H
#define COEFF_MODULUS_CONFIG_H

#include <vector>
#include <cstddef>

class CoeffModulusConfig {
public:
    static std::vector<int> get_coeff_modulus_params(size_t poly_modulus_degree);
};

#endif // COEFF_MODULUS_CONFIG_H 