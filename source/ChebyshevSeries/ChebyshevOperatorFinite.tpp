#pragma once

#include "ChebyshevOperatorFinite.hpp"

template<typename T>
ChebyshevOperatorFinite<T>::ChebyshevOperatorFinite(
        const int N, const int n,
        const capd::vectalg::Vector<T, 0>& y0_init,
        const capd::vectalg::Matrix<T, 0, 0>& g_init
) : N(N), n(n), y0(y0_init), g(g_init), a_series(n), c_series(n) {}

template<typename T>
void ChebyshevOperatorFinite<T>::setASeries(const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_input) {
    this->a_series = a_input;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setCSeries(const std::vector<std::vector<int>>& multiIndices) {
    this->c_series = computeC(multiIndices);
}


template <typename T>
capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION>
ChebyshevOperatorFinite<T>::computeC(const std::vector<std::vector<int>>& multiIndices) {

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
//            std::cout << "Iteracja " << alpha_idx << ": ";
//            std::cout << "a_series[alpha] ^ exponent: " << aj_pow;
//            std::cout << conv << '\n';
        }
//        std::cout << '\n';

        //--------------------------------------------------------------------
//        std::cout << conv << '\n';
//        auto term = conv.getCoefficients();

        // Teraz dodaj g_alpha * term do c
        // Czy tutaj max() czy min()? Czy c będzie mogło miec wiecej wspolczynnikow od a?
        // Znaczy jako wynik splotu sam w sobie to tak, pytanie czy tego chcemy, czy c tez obcinamy do N?
        for (int k = 0; k < this->n; ++k) {
            for (int d = 0; d < std::min(this->N, conv.getN()); ++d) {
                c_result[k][d] += g[alpha][d] * conv[k];
            }
        }
        std::cout << "c_result: "<< c_result << '\n';
    }

    return c_result;
}


template<typename T>
T ChebyshevOperatorFinite<T>::computeF0(
        const ChebyshevSeries<T, DIMENSION>& v,
        const ChebyshevSeries<T, DIMENSION>& u) {

    T result = ChebyshevSeries<T, DIMENSION>::dot(v, u) - ChebyshevSeries<T, DIMENSION>::dot(v, this->a_series[0]);

    for (int k = 1; k < N; ++k) {
        result -= 2 * ChebyshevSeries<T>::dot(v, this->a_series[k]);
    }

    return result;
}

template<typename T>
capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> ChebyshevOperatorFinite<T>::computeF1(
        T omega,
        const std::vector<std::vector<int>>& multiIndices) {

    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> result(n);
    setCSeries(multiIndices);

    // 1. Obliczenie f_1[0]
    result[0] = y0[0] - this->a_series[0];
    for (int l = 1; l < this->N; l++){
        if (l % 2 == 0){
            result[0] = result[0] - 2.0 * this->a_series[l];
        }
        else{
            result[0] = result[0] + 2.0 * this->a_series[l];
        }
    }

    //COS TU NIZEJ NIE DZIALA
    // 2. Obliczanie f_1 dla k > 0
//    for (int k = 1; k < this->N - 1; ++k) {
////        T c_k = 0.0;  // Zmienna pomocnicza do obliczenia c_k
////        auto c = computeC(multiIndices, this->a_series);  // Wykorzystanie computeC() do obliczenia c_k
//        result[k] = omega * k * this->a_series[k] - 0.5 * (this->c_series[k-1] - this->c_series[k+1]);
//    }

    return result;
}

