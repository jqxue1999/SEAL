#include "plain_omp.h"


using namespace std;

void cvps_plain_omp(const vector<uint64>& x, uint64 scalar, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        result[i] = x[i] * scalar;
}
void cvpv_plain_omp(const vector<uint64>& x, const vector<uint64>& y, vector<vector<uint64>>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i)
        cvps_plain_omp(y, x[i], result[i]);
}
void pvcm_plain_omp(const vector<uint64>& x, const vector<vector<uint64>>& Y, vector<uint64>& result) {
    result.resize(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        vector<uint64> partial;
        cvps_plain_omp(Y[i], x[i], partial);
        for (size_t j = 0; j < Y[i].size(); ++j)
            result[j] += partial[j];
    }
}
void pmcm_plain_omp(const vector<vector<uint64>>& X, const vector<vector<uint64>>& Y, vector<vector<uint64>>& result) {
    result.resize(X.size());
    #pragma omp parallel for
    for (size_t i = 0; i < X.size(); ++i)
        pvcm_plain_omp(X[i], Y, result[i]);
}