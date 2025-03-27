#pragma once

#include "ChebyshevOperatorFinite.hpp"

template<typename T>
ChebyshevOperatorFinite<T>::ChebyshevOperatorFinite(
        const int N, const int n,
        const capd::vectalg::Vector<T, 0>& y0_init,
        const capd::vectalg::Matrix<T, 0, 0>& g_init)
        : N(N), n(n), y0(y0_init), g(g_init)
{}

//template <typename T>
//std::vector<std::vector<int>> ChebyshevOperatorFinite<T>::computeMultiIndexMatrix(int maxDegree) const {
////    std::vector<std::vector<int>> multiIndices = generateMultiIndices(n, maxDegree);
////    return multiIndices;
//}

template <typename T>
capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> ChebyshevOperatorFinite<T>::computeC(
        const std::vector<std::vector<int>>& multiIndices,
        const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series){

//    using CVec = ChebyshevSeries<T, DIMENSION>;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> c_result(this->n);

    // kazdy element c_result (czyli wektor) to będzie ChebyshevSeries i ich będzie tyle, jaki wymiar ma a (małe n),
    // ale kazdy z nich bedzie wymiaru DUŻE N
    for (int k = 0; k < this->n; ++k)
        c_result[k] = ChebyshevSeries<T, DIMENSION>(this->N);


    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alhpa - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        ChebyshevSeries<T, DIMENSION> conv(this->N);  // stały 1
        conv[0] = 1.0;
        std::cout << conv;

//        std::cout << "numAlphas: " << numAlphas << '\n';
        for (int alpha_idx = 0; alpha_idx < this->n; ++alpha_idx) {
            int exponent = multiIndices[alpha][alpha_idx];
            std::cout << "Exponent: "<< exponent << ", alpha: " << alpha << ", alpha_idx: " << alpha_idx << '\n';
            std::cout << "a_series[alpha]: " << a_series[alpha_idx];
            ChebyshevSeries<T, DIMENSION> aj_pow = a_series[alpha_idx].power(exponent);
            conv = ChebyshevSeries<T>::convolve(conv, aj_pow);
            std::cout << "Iteracja " << alpha_idx << ": ";
//            std::cout << "a_series[alpha] ^ exponent: " << aj_pow;
            std::cout << conv << '\n';
        }
        std::cout << '\n';

        //--------------------------------------------------------------------
//        std::cout << conv << '\n';
//        auto term = conv.getCoefficients();

//        // Teraz dodaj g_alpha * term do c
//        for (int k = 0; k < std::min(N, term.dimension()); ++k) {
//            for (int d = 0; d < DIMENSION; ++d) {
//                c_result[k][d] += g[d][alpha_idx] * term[k];
//            }
//        }
    }

    return c_result;
}
