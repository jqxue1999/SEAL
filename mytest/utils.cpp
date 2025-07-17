#include "utils.h"
#include <algorithm>
#include <unistd.h>  // for getcwd
#include <fstream>

using namespace seal;
using namespace std;

void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::bgv:
        scheme_name = "BGV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

void print_matrix(const vector<vector<uint64_t>>& matrix, const string& name) {
    cout << name << ":" << endl;
    for (size_t i = 0; i < min(size_t(5), matrix.size()); i++) {
        cout << "  Row " << i << ": ";
        for (size_t j = 0; j < min(size_t(5), matrix[i].size()); j++) {
            cout << matrix[i][j] << " ";
        }
        if (matrix[i].size() > 5) cout << "...";
        cout << endl;
    }
    if (matrix.size() > 5) cout << "  ..." << endl;
    cout << "  Size: " << matrix.size() << "x" << (matrix.empty() ? 0 : matrix[0].size()) << endl;
}

void print_ciphertext_info(const vector<Ciphertext>& ciphertexts, const string& name) {
    cout << name << ":" << endl;
    cout << "  Count: " << ciphertexts.size() << endl;
    if (!ciphertexts.empty()) {
        cout << "  Size: " << ciphertexts[0].size() << " polynomials" << endl;
        cout << "  Coeff modulus size: " << ciphertexts[0].coeff_modulus_size() << endl;
        cout << "  Poly modulus degree: " << ciphertexts[0].poly_modulus_degree() << endl;
    }
}

void check_matrix_equal(const vector<vector<uint64_t>>& matrix1, const vector<vector<uint64_t>>& matrix2) {
    if (matrix1.size() != matrix2.size()) {
        cout << "矩阵大小不一致" << endl;
        return;
    }
    for (size_t i = 0; i < matrix1.size(); i++) {
        for (size_t j = 0; j < matrix1[i].size(); j++) {
            if (matrix1[i][j] != matrix2[i][j]) {
                cout << "矩阵不一致" << endl;
                return;
            }
        }
    }
    cout << "矩阵一致" << endl;
}

json read_seal_config(const string& config_file, bool verbose) {
    json config;
    
    // 打印当前工作目录
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        if (verbose)
            cout << "Current working directory: " << cwd << endl;
    } else {
            cout << "Cannot get current working directory" << endl;
    }
    
    // 直接使用固定路径
    string config_path = config_file;
    if (verbose)
        cout << "Try to read config file: " << config_path << endl;
    
    try {
        ifstream file(config_path);
        if (file.is_open()) {
            file >> config;
            file.close();
            if (verbose)
                cout << "Successfully read config file: " << config_path << endl;
            return config;
        } else {
            cerr << "Error: cannot open config file " << config_path << endl;
        }
    } catch (const exception& e) {
        cerr << "Error: cannot parse config file: " << e.what() << endl;
    }
    
    return config;
}

vector<int> get_coeff_modulus_params(const json& config, size_t poly_modulus_degree) {
    string degree_str = to_string(poly_modulus_degree);
    if (config["coeff_modulus_configs"].contains(degree_str)) {
        return config["coeff_modulus_configs"][degree_str]["coeff_modulus"].get<vector<int>>();
    } else {
        cerr << "错误: 配置文件中未找到多项式模数次数 " << poly_modulus_degree << " 的参数" << endl;
        return vector<int>();
    }
}
