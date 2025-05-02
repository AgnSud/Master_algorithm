#pragma once

#include "ChebyshevOperatorInfinite.hpp"

template <typename T>
ChebyshevOperatorInfinite<T>::ChebyshevOperatorInfinite(int N, const ChebyshevOperatorFinite<T>& finiteOp)
        : N(N), finiteOp(finiteOp) {}

//ok
template <typename T>
T ChebyshevOperatorInfinite<T>::Pi0(const VectorType& x) const {
    return x[0];
}

//ok - ale Pi1(omega, a) = Pi1(x) to to samo co PiN(x) - to tylko moze zrezygnowac z Pi1? - zostaiwam Pi1
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
    VectorType result(a.dimension());
    for (int i = 0; i < a.dimension(); ++i) {
        if (i < N) {
            result[i] = a[i];
        } else {
            result[i] = 0;
        }
    }
    return result;
}

//ok
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::PiN_x(const VectorType& x) const {
    VectorType result(x.dimension());
    result[0] = x[0];
    VectorType a(x.dimension() - 1);
    for (int i = 0; i < x.dimension() -1; ++i) {
        a[i] = x[i + 1];
    }

    VectorType projected = PiN(a);
    for (int i = 0; i < projected.dimension(); ++i) {
        result[i + 1] = projected[i];
    }
    return result;
}

template <typename T>
T ChebyshevOperatorInfinite<T>::Pi0_HatA(const VectorType& x) const {
    MatrixType derivative_finite = finiteOp.getDerivativeFinite();
    VectorType multiply_derivative_finite_Pi_N = derivative_finite * PiN_x(x);
    return Pi0(multiply_derivative_finite_Pi_N);
}


//template <typename T>
//typename ChebyshevOperatorInfinite<T>::VectorType
//ChebyshevOperatorInfinite<T>::applyHatA(const VectorType& x) const {
//    VectorType result(x.dimension());
//
//    // Apply Π₀Â(x) = Π₀ DF_N(x*) Π_N(x)
//    T omega_proj = Pi0(x);
//    VectorType a_proj = Pi1(x);
//    VectorType PiN_x = PiN(a_proj);
//
//    VectorType df_x = finiteOp.getDerivativeFinite() * PiN_x;
//    result[0] = df_x[0]; // Only Π₀ component
//
//    // Apply Π₁Â(x)
//    VectorType Pi1_x = Pi1(x);
//    for (int k = 0; k < Pi1_x.dimension(); ++k) {
//        if (k < N * n) {
//            result[k + 1] = df_x[k + 1]; // project DF_N(x*) Π_N(x)
//        } else {
//            result[k + 1] = omega_proj * (k + 1 - N * n) * Pi1_x[k]; // ω * k * (Π₁(x))_k
//        }
//    }
//
//    return result;
//}