#pragma once

#include "ChebyshevSeries.hpp"

template <typename T, int DIM>
ChebyshevSeries<T, DIM>::ChebyshevSeries(std::initializer_list<T> list)
        : N(static_cast<int>(list.size())), CVector(list.size()) {

    int i = 0;
    for (auto val : list) {
        (*this)[i++] = val;
    }
}

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
    for (int i = 0; i < result.N; ++i) {
        if (i <= this->N) result[i] += (*this)[i];
        if (i <= other.N) result[i] += other[i];
    }
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::operator-(const ChebyshevSeries<T, DIM>& other) const {
    ChebyshevSeries<T, DIM> result(std::max(this->N, other.N));
    for (int i = 0; i < result.N; ++i) {
        if (i <= this->N) result[i] -= (*this)[i];
        if (i <= other.N) result[i] -= other[i];
    }
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::convolve(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b) {
//    CVector a_coefs = a.getCoefficients();
//    CVector b_coefs = b.getCoefficients();
    int size = a.getN() + b.getN() - 1;
    ChebyshevSeries<T, DIM> result(size);
    int n = std::max(a.getN() - 1, b.getN() - 1);
    int negate = -1 * n;
    for (int k = 0; k < size; ++k) {
        result[k] = 0;
        for (int k1 = negate; k1 <= n; k1++) {
//            std::cout << "k1=" << k1 << "\n";
            int k2 = k - k1;
            if (abs(k1) < a.getN() && abs(k2) < b.getN()) {
                result[k] += a[abs(k1)] * b[abs(k2)];
            }
        }
    }
    return result;
}

//template <typename T, int DIM>
//T ChebyshevSeries<T, DIM>::dot(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b) {
//    T result = 0;
//    int N = std::min(a.getN(), b.getN());
//
//    for (int k = 0; k < N; ++k) {
//        result += a[k] * b[k];
//    }
//
//    return result;
//}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::operator*(const ChebyshevSeries<T, DIM>& other) const {
    ChebyshevSeries<T, DIM> result(this->N + other.N - 1);
    result = convolve(*this, other);
    return result;
}

template <typename T, int DIM>
ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::power(int n) const {
    ChebyshevSeries<T, DIM> tmp(this->N);
    tmp[0] = 1;
    for (int i = 0; i < n; ++i) {
        tmp = tmp * (*this);
//        result[i] = std::pow((*this)[i], n);
    }

    ChebyshevSeries<T, DIM> result(this->N);
    for (int i = 0; i < this->N; i++){
        result[i] = tmp[i];
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