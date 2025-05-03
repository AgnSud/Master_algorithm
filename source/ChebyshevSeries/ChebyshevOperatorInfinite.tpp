#pragma once

#include "ChebyshevOperatorInfinite.hpp"

template <typename T>
ChebyshevOperatorInfinite<T>::ChebyshevOperatorInfinite(int N, int n, const ChebyshevOperatorFinite<T>& finiteOp)
        : N(N), n(n), finiteOp(finiteOp) {}

//ok - zwraca nieskonczone
template <typename T>
T ChebyshevOperatorInfinite<T>::Pi0(const VectorType& x) const {
    return x[0];
}

//ok - zwraca nieskonczone
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::Pi1(const VectorType& x) const {
    VectorType result(x.dimension() - 1);
    for (int i = 1; i < x.dimension(); ++i) {
        result[i - 1] = x[i];
    }
    return result;
}

//wynik zwracany to wektor 1 x nN, przy czym reszta domyślnie jest 0, ale tu utozsamiamy ze skonczonym wektorem
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::PiN(const VectorType& a) const {
    VectorType result(N*n);
    for (int i = 0; i < N*n; ++i) {
        if (i < N * n) {
            result[i] = a[i];
        }
    }
    return result;
}

//ok
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::PiN_x(const VectorType& x) const {
    VectorType result(N*n+1);
    result[0] = x[0];
    VectorType a(N*n);
    for (int i = 0; i < a.dimension(); ++i) {
        a[i] = x[i + 1];
    }

    VectorType projected = PiN(a);
    for (int i = 0; i < projected.dimension(); ++i) {
        result[i + 1] = projected[i];
    }
    return result;
}

template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::get_kth(const VectorType& a, int k){
    VectorType a_k(n);
    for (int i = 0; i < n; i++){
        a_k[i] = a[k * n + i];
    }
    return a_k;
}

template <typename T>
T ChebyshevOperatorInfinite<T>::Pi0_HatA(const VectorType& x) {
    MatrixType derivative_finite = finiteOp.getDerivativeFinite();
    VectorType multiply_derivative_finite_Pi_N = derivative_finite * PiN_x(x);
    return Pi0(multiply_derivative_finite_Pi_N);
}

//wyciagnac WEKTOR odpowiadajacy pozycji k - WEKTOR jest dlugosci n!!!
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::Pi1_HatA_k(const VectorType &x, int k) {
    VectorType result(x.dimension());

    VectorType multiply_derivative_finite_Pi_N = finiteOp.getDerivativeFinite() * PiN_x(x);

    if (k < N){
//        cout << multiply_derivative_finite_Pi_N << endl;
        auto tmp = Pi1(multiply_derivative_finite_Pi_N);
        result = get_kth(tmp, k);
//        cout << result << endl;
    }
    else{
        auto tmp = get_kth(Pi1(x), k);
        result = finiteOp.getOmega() * k * get_kth(Pi1(x), k);
    }

    return result;
}