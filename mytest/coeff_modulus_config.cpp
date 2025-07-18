#include "coeff_modulus_config.h"
#include <stdexcept>

std::vector<int> CoeffModulusConfig::get_coeff_modulus_params(size_t poly_modulus_degree) {
    switch (poly_modulus_degree) {
        case 1024:
            return {27}; // TODO: 替换为实际参数
        case 2048:
            return {27, 27}; // TODO: 替换为实际参数
        case 4096:
            return {36, 36, 37}; // TODO: 替换为实际参数
        case 8192:
            return {36, 36, 36, 36, 37}; // TODO: 替换为实际参数
        case 16384:
            return {36, 36, 36, 36, 37}; // TODO: 替换为实际参数
        default:
            throw std::invalid_argument("Unsupported poly_modulus_degree");
    }
} 