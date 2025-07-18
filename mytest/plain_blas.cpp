#include "plain_blas.h"

using namespace std;

void cvps_plain_blas(const vector<uint64>& x, uint64 scalar, vector<uint64>& result) {
    std::vector<double> x_d(x.size());
    result.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        x_d[i] = static_cast<double>(x[i]);
    double s = static_cast<double>(scalar);
    cblas_dscal(x.size(), s, x_d.data(), 1);
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        result[i] = static_cast<uint64>(x_d[i]);
}
void cvpv_plain_blas(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        cvps_plain_blas(y, x[i], result[i]);
}
void pvcm_plain_blas(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        vector<uint64> partial;
        cvps_plain_blas(Y[i], x[i], partial);
        for (size_t j = 0; j < Y[i].size(); ++j)
            result[j] += partial[j];
    }
}
void pmcm_plain_blas(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result) {
    result.resize(X.size());
    #pragma omp parallel for
    for (size_t i = 0; i < X.size(); ++i)
        pvcm_plain_blas(X[i], Y, result[i]);
}