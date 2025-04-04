#pragma once

#include "ChebyshevOperatorFinite.hpp"
#include "capd/fadbad/fadbad.h"
#include <type_traits>

template<typename T>
ChebyshevOperatorFinite<T>::ChebyshevOperatorFinite(
        const int N, const int n,
        const capd::vectalg::Vector<T, 0>& y0_init,
        const capd::vectalg::Matrix<T, 0, 0>& g_init,
        const ChebyshevSeries<T, 0>& v,
        const ChebyshevSeries<T, 0>& u,
        const std::vector<std::vector<int>>& multiIndices
) : N(N), n(n), omega(0), y0(y0_init), g(g_init), a_series(n), c_series(n), v(v), u(u), multiIndices(multiIndices) {}

template<typename T>
void ChebyshevOperatorFinite<T>::setASeries(const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_input) {
    this->a_series = a_input;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setOmega(const T& omega_su) {
    this->omega = omega_su;
}


template <typename T>
template <class V>
capd::vectalg::Vector<typename V::VectorType, 0>
ChebyshevOperatorFinite<T>::computeC(const V& x) {
    auto [omega_input, a_series_input] = convertXVectorToOmegaAndASeries(x);

    //konwersja a_series_input na ChebyshevSeries dla funkcji power()
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> converted_a_series_input(this->n);
    for (int i = 0; i < this->n; i++){
        converted_a_series_input[i] = ChebyshevSeries<T, DIMENSION>(this->N);
        for (int k = 0; k < this->N; k++){
            if constexpr (std::is_same<decltype(a_series_input[i][k]), double&>::value){
                converted_a_series_input[i][k] = a_series_input[i][k];
            }
            else{
                converted_a_series_input[i][k] = a_series_input[i][k].x();
            }
        }
    }
    capd::vectalg::Vector<typename V::VectorType, 0> c_series_result(this->n);
    for (int i = 0; i < this->n; ++i)
        c_series_result[i] = typename V::VectorType(this->N + 1);

    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alpha - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        ChebyshevSeries<T, DIMENSION> conv(this->N + 1);  // stały 1
        conv[0] = 1.0;

        for (int alpha_idx = 0; alpha_idx < this->n; ++alpha_idx) {
            int exponent = multiIndices[alpha][alpha_idx];
            ChebyshevSeries<T, DIMENSION> aj_pow = converted_a_series_input[alpha_idx].power(exponent);
            conv = ChebyshevSeries<T>::convolve(conv, aj_pow);
        }

        // Teraz dodaj g_alpha * term do c
        // Potrzebujemy wielkosc c_{k+1}, więc obetniemy c do N+1
        for (int k = 0; k < this->n; ++k) {
            for (int d = 0; d <= this->N; ++d) {
                c_series_result[k][d] += this->g[alpha][k] * conv[d];
            }
        }
    }
    if constexpr (std::is_same<decltype(c_series_result[0][0]), double&>::value)
        std::cout << "c_series_result: " << c_series_result << '\n';

    return c_series_result;
}

template<typename T>
capd::vectalg::Vector<T, 0> ChebyshevOperatorFinite<T>::convertToXVector() {
    capd::vectalg::Vector<T, 0> x(this->n * this->N + 1);

    x[0] = omega;

    int index = 1;
    for (int k = 0; k < this->N; ++k) {
        for (int i = 0; i < this->n; ++i) {
            x[index] = a_series[i][k];
            index++;
        }
    }

    return x;
}

template<typename T>
template <class V>
std::tuple<typename V::ScalarType, capd::vectalg::Vector<typename V::VectorType, 0>>
ChebyshevOperatorFinite<T>::convertXVectorToOmegaAndASeries(const V& x) {
    // Wyciągamy omega z pierwszego elementu wektora x
    typename V::ScalarType omega_result = x[0];

    // Tworzymy a_series z pozostałych elementów wektora x
    capd::vectalg::Vector<typename V::VectorType, 0> a_series_result(this->n);
    for (int i = 0; i < this->n; i++)
        a_series_result[i] = typename V::VectorType(this->N);

    int index = 1;
    for (int k = 0; k < this->N; ++k) {
        for (int i = 0; i < this->n; i++) {
            a_series_result[i][k] = x[index];
            index++;
        }
    }
//    if constexpr (std::is_same<decltype(a_series_result[0][0]), double&>::value) {
//        std::cout << "x: " << x << '\n';
//        std::cout << "a_series_result: " << a_series_result << '\n';
//    }

    return {omega_result, a_series_result};
}


template<typename T>
capd::vectalg::Vector<T, DIMENSION> ChebyshevOperatorFinite<T>::findFiniteSolution(
        T omega_start,
        const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series_start,
        int max_iterations, T tolerance) {
    // Ustaw a i omega
    setOmega(omega_start);
    setASeries(a_series_start);
    int iteration = 0;
    T norm = 1.0;

    // Utwórz wynik jako Vector<T, DIMENSION>
    capd::vectalg::Vector<T, DIMENSION> F(this->n * this->N + 1);
    capd::vectalg::Matrix<T, DIMENSION, DIMENSION> jacobian(this->n * this->N + 1, this->n * this->N + 1);

    //x_k_1 = x_k - F_x_k/DF_x_k
    //co oznacza x_{k+1} = x_k - F(x_k) / DF(x_k)
    while (iteration < max_iterations && norm > tolerance) {
        capd::vectalg::Vector<T, 0> x = convertToXVector();
        std::cout << "x = (omega, a) =" << x << '\n';

        auto F_x_k = (*this)(x);
        std::cout << "F(x_k) = " << F_x_k << '\n';
        computeDerivative(*this, x, jacobian);
        std::cout << "Jacobian: " << jacobian << "\n";

        iteration++;
    }
    return F;
}

template<typename T>
template<class V>
V ChebyshevOperatorFinite<T>::operator() (const V& x) {
    V result(x.dimension());

    result[0] = computeF0(x);
    auto f1_result = computeF1(x);
    int index = 1;
    for (int k = 0; k < this->N; ++k) {
        for (int i = 0; i < this->n; ++i) {
            result[index] = f1_result[i][k];
            index++;
        }
    }

    if constexpr (std::is_same<decltype(f1_result[0][0]), double&>::value) {
        std::cout << "f1_result: " << f1_result << '\n';
        std::cout << "result: " << result << '\n';
    }

    return result;
}


template<typename T>
template<class V>
typename V::ScalarType ChebyshevOperatorFinite<T>::computeF0(const V& x){
    auto [omega_input, a_series_input] = convertXVectorToOmegaAndASeries(x);
    typename V::ScalarType result = 0;
    for (int i = 0; i < this->n; i++){
        // <v, u>
        result += v[i] * u[i];

        // <v, a_0>
        result -= v[i] * a_series_input[i][0];

        // sum_{k=1}^{N-1} 2<v, a_k>
        for (int k = 1; k < this->N; ++k) {
            result -= v[i] * 2 * a_series_input[i][k];
        }
    }
    return result;
}


template<typename T>
template<class V>
capd::vectalg::Vector<typename V::VectorType, 0> ChebyshevOperatorFinite<T>::computeF1(const V& x) {

    capd::vectalg::Vector<typename V::VectorType, 0> result(this->n);
    auto [omega_input, a_series_input] = convertXVectorToOmegaAndASeries(x);
    auto c_series_input = computeC(x);

    // 1. Obliczenie f_1[0]
    for (int i = 0; i < this->n; i++){
        result[i] = typename V::VectorType(this->N);
        result[i][0] = this->y0[i] - this->a_series[i][0];
        for (int l = 1; l < this->N; l++){
            if (l % 2 == 0){
                result[i][0] = result[i][0] - 2.0 * a_series_input[i][l];
            }
            else{
                result[i][0] = result[i][0] + 2.0 * a_series_input[i][l];
            }
        }
    }

//     2. Obliczanie f_1 dla k > 0
    for (int i = 0; i < this->n; i++){
        for (int k = 1; k < this->N; ++k) {
            result[i][k] = omega_input * k * a_series_input[i][k] - 0.5 * (c_series_input[i][k-1] - c_series_input[i][k+1]);
        }
    }
    return result;
}


