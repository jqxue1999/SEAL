#include "plain_blas.h"

using namespace std;

void cvps_plain_blas(const vector<double>& x, double scalar, vector<double>& result) {
    result.resize(x.size());
    cblas_dscal(x.size(), scalar, const_cast<double*>(x.data()), 1);
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        result[i] = x[i];
}
void cvpv_plain_blas(const vector<double>& x, const vector<double>& y, vector<vector<double>>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        cvps_plain_blas(y, x[i], result[i]);
}
void pvcm_plain_blas(const vector<double>& x, const vector<vector<double>>& Y, vector<double>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        vector<double> partial;
        cvps_plain_blas(Y[i], x[i], partial);
        for (size_t j = 0; j < Y[i].size(); ++j)
            result[j] += partial[j];
    }
}
void pmcm_plain_blas(const vector<vector<double>>& X, const vector<vector<double>>& Y, vector<vector<double>>& result) {
    result.resize(X.size());
    #pragma omp parallel for
    for (size_t i = 0; i < X.size(); ++i)
        pvcm_plain_blas(X[i], Y, result[i]);
}