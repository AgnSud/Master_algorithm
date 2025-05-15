#pragma once

#include "RadiiPolynomials.hpp"

template<typename T>
RadiiPolynomials<T>::RadiiPolynomials(int N_, int n_, double nu_, const ChebyshevOperatorFinite<T>& finiteOp_)
        : N(N_), n(n_), nu(nu_), finiteOp(finiteOp_) {}

template <typename T>
T RadiiPolynomials<T>::Pi0(const VectorType& x) {
    return x[0];
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::Pi1(const VectorType& x) {
    VectorType result(x.dimension() - 1);
    for (int i = 1; i < x.dimension(); ++i) {
        result[i - 1] = x[i];
    }
    return result;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::PiN(const VectorType& a, int N_, int n_) {
    VectorType result(N_ * n_);
    for (int i = 0; i < N_ * n_; ++i) {
        result[i] = a[i];
    }
    return result;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::PiN_x(const VectorType& x, int N_, int n_) {
    VectorType result(N_ * n_ + 1);
    result[0] = x[0];
    VectorType a(N_ * n_);
    for (int i = 0; i < a.dimension(); ++i) {
        a[i] = x[i + 1];
    }

    VectorType projected = PiN(a, N_, n_);
    for (int i = 0; i < projected.dimension(); ++i) {
        result[i + 1] = projected[i];
    }
    return result;
}

template<typename T>
template<class V>
V RadiiPolynomials<T>::Pi1_j(const V &x, int j, int N_, int n_) {
    return ChebyshevOperatorFinite<T>::getCoeffVectorI_thSquareParan(x, j, N_, n_);
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::get_kth(const VectorType& a, int k){
    VectorType a_k(n);
    for (int i = 0; i < n; ++i){
        a_k[i] = a[k * n + i];
    }
    return a_k;
}



template <typename T>
T RadiiPolynomials<T>::computeY0() {

    VectorType x_approx = finiteOp.convertToXVector();          // x^*
//    cout << "x_approx = " << x_approx << endl;
    VectorType F_x_approx = finiteOp.getF_x_approx();                        // F_N(x^*)
//    cout << "F_x_approx = " << F_x_approx << endl;
    MatrixType A_N = finiteOp.getInverseDerivativeFinite();   // A_N
//    cout << "A_N = " << A_N << endl;
    VectorType AF = A_N * F_x_approx;                         // A_N F_N(x^*)
//    cout << "A_N * F_x_approx = " << AF << endl;

    return Pi0(AF); // Pi_0 — tylko pierwsza współrzędna
}

template <typename T>
T RadiiPolynomials<T>::computeY1j(int j, int N_g) {

    // Parametry
    T omega = finiteOp.getOmega();
    auto a_series = finiteOp.getASeries();
    auto c_series = finiteOp.getCSeries();

    T sum = 0;

    for (int k = N; k <= N_g * (N - 1) + 1; ++k) {
        T diff = c_series[j][k - 1];
        if (k + 1 < c_series[j].dimension())
            diff -= c_series[j][k + 1];
        diff = capd::abs<T>(diff);
//        cout << "k=" << k << ", [c_{k-1} - c_{k+1}]_" << j << " = " << diff << endl;
        T weighted = diff * std::pow(nu, k) / k;
//        cout << "weighted by nu^k/k =" << std::pow(nu, k) / k << " results: " << weighted << endl;
        sum += weighted;
    }

    // część druga: || Pi_{1,j} A_N F_N(x^*) ||_nu
    VectorType x_approx = finiteOp.convertToXVector();
    VectorType F_x_approx = finiteOp.getF_x_approx();
    auto A_N = finiteOp.getInverseDerivativeFinite();
    VectorType AF = Pi1_j(A_N * F_x_approx, j, N, n);

    // Pi_1,j — współczynniki j-tego ciągu a_k (czyli co n-ta współrzędna od j)
    Norm<T> weighted_norm(nu, N, n);
    auto w_norm_AF = weighted_norm.computeNorm(AF);
//    cout << "w_norm_AF = " << w_norm_AF << endl;
//    cout << "sum = " << sum << endl;
    T result = sum / finiteOp.getOmega() + w_norm_AF;

    return result;
}

template <typename T>
void RadiiPolynomials<T>::testOperatorNorm() {
    Norm<T> weighted_norm(nu, N, n);
    auto tmp = weighted_norm.computeOperatorNorm(finiteOp.getInverseDerivativeFinite());
}
