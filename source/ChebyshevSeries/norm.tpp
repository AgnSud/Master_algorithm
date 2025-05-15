#pragma once

#include "Norm.hpp"


template <typename T, int DIM>
//T Norm<T, DIM>::computeNorm(const ChebyshevSeries<T, DIM>& a) const {
template<class V>
T Norm<T, DIM>::computeNorm(const V& a) const {
    T sum = capd::abs<T>(a[0]); // Zaczynamy od |a_0|
    for (int k = 1; k < a.dimension(); ++k) {
        sum += 2 * capd::abs<T>(a[k]) * std::pow(nu, k); // 2 * |a_k| * nu^k
    }
    return sum;
}


// TODO: Czy to jest dobrze? Czy dobrze zrozumiałam wzór?
template <typename T, int DIM>
template<class V>
T Norm<T, DIM>::computeDualNorm(const V& a) const {
    T maxTerm = capd::abs(a[0]);
//    cout << "term for " << a[0] << "= " << maxTerm << endl;
    for (int k = 1; k < a.dimension(); ++k) {
        T term = capd::abs(a[k]) / (2 * std::pow(nu, k));
//        cout << "term for " << a[k] << "= " << term << endl;
        maxTerm = capd::max(maxTerm, term); // zastepuje sup -> max, poniewaz i tak numerycznie jest max?
    }
    return maxTerm;
}

template <typename T, int DIM>
template<class V>
T Norm<T, DIM>::computeNorm_n(const capd::vectalg::Vector<V, DIM>& vec) const {
    T result = 0;

    // Maksymalizujemy po każdej składowej ChebyshevSeries
    for (int j = 0; j < n; ++j) {
        T sum = computeNorm(vec[j]);
        result = capd::max<T>(result, sum); // Maksymalizujemy dla każdego szeregu
    }
    return result;
}

template <typename T, int DIM>
T Norm<T, DIM>::computeOperatorNorm_Pi0(const MatrixType& C) {
    T result_norm = C[0][0];
    for(int j_tilde = 0; j_tilde < n; j_tilde++){
        result_norm += mu_j(C, j_tilde);
    }
    return result_norm;
}

template <typename T, int DIM>
T Norm<T, DIM>::computeOperatorNorm_Pi1j(const MatrixType& C, int j) {
    T result_norm = eta_j(C, j);
    for(int j_tilde = 0; j_tilde < n; j_tilde++){
        result_norm += ksi_j_j(C, j, j_tilde);
    }
    return result_norm;
}


template <typename T, int DIM>
T Norm<T, DIM>::computeOperatorNorm(const MatrixType& C) {
    auto Pi0_C = computeOperatorNorm_Pi0(C);
    cout << "||Pi0_C|| = " << Pi0_C << endl;

    for(int j = 0; j < n; j++){
        auto Pi1j_C = computeOperatorNorm_Pi1j(C, j);
        cout << "||Pi1" << j << "_C|| = " << Pi1j_C << endl;
    }
    return C[0][0];
}


template <typename T, int DIM>
T Norm<T, DIM>::eta_j(const MatrixType &C, int j) {
    VectorType C_a_0(N);
    for(int k = 0; k < N; k++){
        C_a_0[k] = C[1 + n * k + j][0];
    }
    T result = computeNorm(C_a_0);
//    cout << "(C_a^0)_j. = " << C_a_0 << endl;
    return result;
}

template <typename T, int DIM>
T Norm<T, DIM>::mu_j(const MatrixType &C, int j_tilde) {
    VectorType C_0_a(N);
    for (int k_tilde = 0; k_tilde < N; k_tilde++) {
        C_0_a[k_tilde] = C[0][1 + n * k_tilde + j_tilde];
    }
    T result = computeDualNorm(C_0_a);
//    cout << "(C_0^a)^j. = " << C_0_a << endl;
    return result;
}

template <typename T, int DIM>
T Norm<T, DIM>::ksi_tilde_j_jk(const MatrixType &C, int j, int j_tilde, int k_tilde) {
    VectorType C_a_a(N);
    for (int k = 0; k < N; k++) {
        C_a_a[k] = C[1 + n * k + j][1 + n * k_tilde + j_tilde];
    }
    T result = computeNorm(C_a_a);
//    cout << "(C_a^a)_" << j << ".^" << j_tilde << k_tilde << " RESULT = " << result << endl;
    return result;
}

template <typename T, int DIM>
T Norm<T, DIM>::ksi_j_j(const MatrixType &C, int j, int j_tilde) {
    VectorType ksi_tilde(N);
    for (int k_tilde = 0; k_tilde < N; k_tilde++) {
        ksi_tilde[k_tilde] = ksi_tilde_j_jk(C, j, j_tilde, k_tilde);
    }
    T result = computeDualNorm(ksi_tilde);
//    cout << "ksi_tilde_" << j << "^" << j_tilde << ". RESULT = " << result << endl;
    return result;
}
