#pragma once

#include "norm.hpp"


template <typename T, int DIM>
//T norm<T, DIM>::computeNorm(const ChebyshevSeries<T, DIM>& a) const {
template<class V>
T norm<T, DIM>::computeNorm(const V& a) const {
    T sum = std::abs(a[0]); // Zaczynamy od |a_0|
    for (int k = 1; k < a.dimension(); ++k) {
        sum += 2 * std::abs(a[k]) * std::pow(nu, k); // 2 * |a_k| * nu^k
    }
    return sum;
}

template <typename T, int DIM>
template<class V>
T norm<T, DIM>::computeNorm(const capd::vectalg::Vector<V, DIM>& vec, int n) const {
    T result = 0;

    // Maksymalizujemy po każdej składowej ChebyshevSeries
    for (int j = 0; j < n; ++j) {
        T sum = computeNorm(vec[j]);
        result = std::max(result, sum); // Maksymalizujemy dla każdego szeregu
    }
    return result;
}

template <typename T, int DIM>
template<class V>
T norm<T, DIM>::computeOperatorNorm(const V& a) const {
    //TBD
}
