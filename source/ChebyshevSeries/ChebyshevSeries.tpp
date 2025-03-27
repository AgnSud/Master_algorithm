#pragma once

#include "ChebyshevSeries.hpp"

template <typename T, int DIM>
T ChebyshevSeries<T, DIM>::evaluateFirstKind(int k, T x) {
    if (k == 0) return 1;
    if (k == 1) return x;
    return 2 * x * evaluateFirstKind(k - 1, x) - evaluateFirstKind(k - 2, x);
}

template<typename T, int DIM>
typename ChebyshevSeries<T, DIM>::CVector
ChebyshevSeries<T, DIM>::getCoefficients() const {
    return *this;
}

template<typename T, int DIM>
int ChebyshevSeries<T, DIM>::getN() const {
    return this->N;
}

template <typename T, int DIM>
T ChebyshevSeries<T, DIM>::operator()(T x) const {
    T sum = (*this)[0];  // a_0
    for (int k = 1; k < this->N; ++k) {
        sum += 2 * (*this)[k] * evaluateFirstKind(k, x);  // 2 * a_k * T_k(x)
    }
    return sum;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::operator+(const ChebyshevSeries<T, DIM>& other) const {
    ChebyshevSeries<T, DIM> result(std::max(this->N, other.N));
    for (int i = 0; i <= result.N; ++i) {
        if (i <= this->N) result[i] += (*this)[i];
        if (i <= other.N) result[i] += other[i];
    }
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::convolve(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b) {
    CVector a_coefs = a.getCoefficients();
    CVector b_coefs = b.getCoefficients();
    int size = a_coefs.dimension() + b_coefs.dimension() - 1;
    ChebyshevSeries<T, DIM> result(size);
    int n = std::max(a_coefs.dimension() - 1, b_coefs.dimension() - 1);
    int negate = -1 * n;
    for (int k = 0; k < size; ++k) {
        result[k] = 0;
        for (int k1 = negate; k1 <= n; k1++) {
//            std::cout << "k1=" << k1 << "\n";
            int k2 = k - k1;
            if (abs(k1) < a_coefs.dimension() && abs(k2) < b_coefs.dimension()) {
                result[k] += a_coefs[abs(k1)] * b_coefs[abs(k2)];
            }
        }
    }
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::operator*(const ChebyshevSeries<T, DIM>& other) const {
    ChebyshevSeries<T, DIM> result(this->N + other.N - 1);
    result = convolve(*this, other);
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::power(int n) const {
//    std::cout << "N: " << this->N << '\n';
    ChebyshevSeries<T, DIM> result(this->N);
//    std::cout << result;
    for (int i = 0; i < this->N; ++i) {
//        std::cout << "Iteracja " << i << ": ";
        result[i] = std::pow((*this)[i], n);
    }
    return result;
}


template <typename T, int DIM>
void ChebyshevSeries<T, DIM>::prettyPrint() const {
    std::cout << "Chebyshev Series Expansion:\n";
    std::cout << "{";
    for (int i = 0; i < this->N-1; ++i) {
        std::cout << "a_" << i << ": " << (*this)[i] << ", ";
    }
    std::cout << "a_" << this->N - 1 << ": " << (*this)[this->N-1] << "}";
}