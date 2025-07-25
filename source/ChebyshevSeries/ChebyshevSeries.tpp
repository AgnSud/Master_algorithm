#pragma once

#include "ChebyshevSeries.hpp"

using namespace capd;
using namespace std;

template <typename T, int DIM>
ChebyshevSeries<T, DIM>::ChebyshevSeries(initializer_list<T> list)
        : N(static_cast<int>(list.size())), VectorType(list.size()) {

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
typename ChebyshevSeries<T, DIM>::VectorType
ChebyshevSeries<T, DIM>::getCoefficients() const {
    return *this;
}

template<typename T, int DIM>
void ChebyshevSeries<T, DIM>::setCoefficients(const VectorType& x) {
    if (x.dimension() != this->N) {
        cerr << "Error: The size of the input vector does not match the number of coefficients.\n";
        return;
    }
    for (int i = 0; i < this->N; ++i) {
        (*this)[i] = x[i];
    }
}

template<typename T, int DIM>
int ChebyshevSeries<T, DIM>::getN() const {
    return this->N;
}

template <typename T, int DIM>
T ChebyshevSeries<T, DIM>::operator()(T t) const {
    T sum = (*this)[0];  // a_0
    for (int k = 1; k < this->N; ++k) {
        T T_k = evaluateFirstKind(k, t);
        auto a_k = (*this)[k];
        sum += 2 * a_k * T_k;  // 2 * a_k * T_k(x)
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

//template <typename T, int DIM>
//ChebyshevSeries<T, DIM> ChebyshevSeries<T, DIM>::convolve(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b) {
////    CVector a_coefs = a.getCoefficients();
////    CVector b_coefs = b.getCoefficients();
//    int size = a.getN() + b.getN() - 1;
//    ChebyshevSeries<T, DIM> result(size);
//    int n = std::max(a.getN() - 1, b.getN() - 1);
//    int negate = -1 * n;
//    for (int k = 0; k < size; ++k) {
//        result[k] = 0;
//        for (int k1 = negate; k1 <= n; k1++) {
////            cout << "k1=" << k1 << "\n";
//            int k2 = k - k1;
//            if (abs(k1) < a.getN() && abs(k2) < b.getN()) {
//                result[k] += a[abs(k1)] * b[abs(k2)];
//            }
//        }
//    }
//    return result;
//}

template <typename T, int DIM>
template<class V>
V ChebyshevSeries<T, DIM>::convolve(const V& a, const V& b) {

    int size = a.dimension() + b.dimension() - 1;
    V result(size);
    int negate = -1 * std::max(a.dimension() -1, b.dimension() - 1);
    int upper=std::max(a.dimension() -1, b.dimension() - 1);
    for (int k = 0; k < size; ++k) {
        result[k] = 0;
        for (int k1 = negate; k1 <= upper; k1++) {
            int k2 = k - k1;
            if (abs(k1) < a.dimension() && abs(k2) < b.dimension()) {
                result[k] += a[abs(k1)] * b[abs(k2)];
            }
        }
    }
    return result;
}

template <typename T, int DIM>
template<class V>
V ChebyshevSeries<T, DIM>::operator*(const V& other) const {
    V result(this->N + other.N - 1);
    result = convolve(*this, other);
    return result;
}

//template <typename T, int DIM>
//template<class V>
//V ChebyshevSeries<T, DIM>::power(int n) const {
//    V tmp(this->N);
//    tmp[0] = 1;
//    for (int i = 0; i < n; ++i) {
//        tmp = tmp * (*this);
////        result[i] = std::pow((*this)[i], n);
//    }
//
//    return tmp;
//}


template <typename T, int DIM>
void ChebyshevSeries<T, DIM>::prettyPrint() const {
    cout << "Chebyshev Series Expansion:\n";
    cout << "{";
    for (int i = 0; i < this->N-1; ++i) {
        cout << "a_" << i << ": " << (*this)[i] << ", ";
    }
    cout << "a_" << this->N - 1 << ": " << (*this)[this->N-1] << "}";
}