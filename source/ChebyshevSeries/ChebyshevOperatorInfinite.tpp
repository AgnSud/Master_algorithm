#pragma once

#include "ChebyshevOperatorInfinite.hpp"

template <typename T>
ChebyshevOperatorInfinite<T>::ChebyshevOperatorInfinite(int N, int n, const ChebyshevOperatorFinite<T>& finiteOp)
        : N(N), n(n), finiteOp(finiteOp), pi0_hatA(0), pi1_hatA(N*n) {}


template<typename T>
void ChebyshevOperatorInfinite<T>::set_pi0_hatA(const T& pi0_hatA_input) {
    this->pi0_hatA = pi0_hatA_input;
}
template <typename T>
T ChebyshevOperatorInfinite<T>::get_pi0_hatA() const {
    return this->pi0_hatA;
}

template<typename T>
void ChebyshevOperatorInfinite<T>::set_pi1_hatA(const VectorType& pi1_hatA_input) {
    this->pi1_hatA = pi1_hatA_input;
}
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::get_pi1_hatA() const {
    return this->pi1_hatA;
}

//ok - zwraca nieskonczone
template <typename T>
T ChebyshevOperatorInfinite<T>::Pi0(const VectorType& x) {
    return x[0];
}

//ok - zwraca nieskonczone
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::Pi1(const VectorType& x) {
    VectorType result(x.dimension() - 1);
    for (int i = 1; i < x.dimension(); ++i) {
        result[i - 1] = x[i];
    }
    return result;
}

template<typename T>
template<class V>
V ChebyshevOperatorInfinite<T>::Pi1_j(const V &x, int j, int N_, int n_) {
    return ChebyshevOperatorFinite<T>::getCoeffVectorI_thSquareParan(x, j, N_, n_);
}

//wynik zwracany to wektor 1 x nN, przy czym reszta domyślnie jest 0, ale tu utozsamiamy ze skonczonym wektorem
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::PiN(const VectorType& a, int N_, int n_) {
    VectorType result(N_*n_);
    for (int i = 0; i < N_*n_; ++i) {
        if (i < N_ * n_) {
            result[i] = a[i];
        }
    }
    return result;
}

//ok
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::PiN_x(const VectorType& x, int N_, int n_) {
    VectorType result(N_*n_+1);
    result[0] = x[0];
    VectorType a(N_*n_);
    for (int i = 0; i < a.dimension(); ++i) {
        a[i] = x[i + 1];
    }

    VectorType projected = PiN(a, N_, n_);
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
    VectorType multiply_derivative_finite_Pi_N = derivative_finite * PiN_x(x, this->N, this->n);
    return Pi0(multiply_derivative_finite_Pi_N);
}

//wyciagnac WEKTOR odpowiadajacy pozycji k - WEKTOR jest dlugosci n!!!
template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::Pi1_HatA_k(const VectorType &x, int k) {
    VectorType result(x.dimension());

    VectorType multiply_derivative_finite_Pi_N = finiteOp.getDerivativeFinite() * PiN_x(x, this->N, this->n);

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

template <typename T>
typename ChebyshevOperatorInfinite<T>::VectorType
ChebyshevOperatorInfinite<T>::Pi1_HatA(const VectorType& x) {
    VectorType a = Pi1(x);
    int K = a.dimension() / n;

    VectorType result(a.dimension());

    for (int k = 0; k < K; ++k) {
        VectorType kth = Pi1_HatA_k(x, k);
        for (int i = 0; i < n; ++i) {
            result[k * n + i] = kth[i];
        }
    }

    return result;
}

template <typename T>
std::pair<T, typename ChebyshevOperatorInfinite<T>::VectorType>
ChebyshevOperatorInfinite<T>::HatA(const VectorType& x) {
    this->set_pi0_hatA(Pi0_HatA(x));
    this->set_pi1_hatA(Pi1_HatA(x));
    return std::make_pair(pi0_hatA, pi1_hatA);
}