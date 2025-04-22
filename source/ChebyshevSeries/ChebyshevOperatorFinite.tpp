#pragma once

#include "ChebyshevOperatorFinite.hpp"
#include "capd/fadbad/fadbad.h"
#include <type_traits>



template<typename T>
ChebyshevOperatorFinite<T>::ChebyshevOperatorFinite(
        const int N, const int n,
        const VectorType& u0_init,
        const MatrixType& g_init,
        const ChebyshevSeries<T, 0>& v,
        const ChebyshevSeries<T, 0>& w,
        const std::vector<std::vector<int>>& multiIndices
) : N(N), n(n), omega(0), u0(u0_init), g(g_init), a_series(n), c_series(n), v(v), w(w), multiIndices(multiIndices) {}

template<typename T>
void ChebyshevOperatorFinite<T>::setASeries(const VectorOfChebyshevsType& a_input) {
    this->a_series = a_input;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setOmega(const T& omega_su) {
    this->omega = omega_su;
}

template <typename T>
template <class V>
V ChebyshevOperatorFinite<T>::multiply(const V& a, const V& b) {
    V result(a.dimension() + b.dimension() - 1);
    result = ChebyshevSeries<T>::convolve(a, b);
    return result;
}

template <typename T>
template <class V>
typename ChebyshevOperatorFinite<T>::template Types<V>::VectorOfVectorsVType ChebyshevOperatorFinite<T>::compute_c(const V& x) {
    typename Types<V>::VectorOfVectorsVType c_series_result(this->n);
    for (int i = 0; i < this->n; ++i)
        c_series_result[i] = typename V::VectorType(this->N + 1);

    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alpha - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        V conv(this->N + 1);  // stały 1
        conv[0] = 1.0;

        for (int alpha_idx = 0; alpha_idx < this->n; ++alpha_idx) {
            int exponent = multiIndices[alpha][alpha_idx];
            //compute power:
            V a_square_paran_alpha_idx = getCoeffVectorI_thSquareParan(x, alpha_idx);
            V aj_pow(this->N);

            aj_pow[0] = 1;
            for (int i = 0; i < exponent; ++i) {
//                aj_pow = aj_pow * a_square_paran_alpha_idx;
                aj_pow = multiply(aj_pow, a_square_paran_alpha_idx);
            }
            conv = ChebyshevSeries<T>::convolve(conv, aj_pow);
//            if constexpr (std::is_same<decltype(conv[0]), double&>::value){
//                std::cout << "aj_pow: " << aj_pow << '\n';
//                std::cout << "conv: " << conv << '\n';
//            }
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
typename ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::convertToXVector() {
    VectorType x(this->n * this->N + 1);

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

//zwraca [a_k]_i
template<typename T>
template<class V>
inline typename V::ScalarType ChebyshevOperatorFinite<T>::getCoeff(const V& x, int i, int k, bool is_omega) const {
    if (is_omega)
        return x[0];
    return x[1 + k * this->n + i];
}


//zwraca [a]_i
template<typename T>
template<class V>
inline typename V::VectorType ChebyshevOperatorFinite<T>::getCoeffVectorI_thSquareParan(const V& x, int i) const {
    V result(this->N);
    for (int k = 0; k < this->N; ++k) {
        result[k] = x[1 + k * this->n + i];
    }
    return result;
}


template<typename T>
typename ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::findFiniteSolution(
        T omega_start,
        const VectorOfChebyshevsType& a_series_start,
        int max_iterations, T tolerance) {
    // Ustaw a i omega
    setOmega(omega_start);
    setASeries(a_series_start);
    int iteration = 0;
    T norm = 1.0;

    // Utwórz wynik jako Vector<T, DIMENSION>
    VectorType F(this->n * this->N + 1);
    MatrixType jacobian(this->n * this->N + 1, this->n * this->N + 1);

    //x_k_1 = x_k - F_x_k/DF_x_k
    //co oznacza x_{k+1} = x_k - F(x_k) / DF(x_k)
    while (iteration < max_iterations && norm > tolerance) {
        VectorType x = convertToXVector();
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

    result[0] = compute_f_0(x);
    auto f1_result = compute_f_1(x);
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
typename V::ScalarType ChebyshevOperatorFinite<T>::compute_f_0(const V& x){
    typename V::ScalarType result = 0;
    for (int i = 0; i < this->n; i++){
        // <v, w>
        result += v[i] * w[i];

        // <v, a_0>
        result -= v[i] * this->getCoeff(x, i, 0);

        // sum_{k=1}^{N-1} 2<v, a_k>
        for (int k = 1; k < this->N; ++k) {
            result -= 2 * v[i] * this->getCoeff(x, i, k);
        }
    }
    return result;
}


template<typename T>
template<class V>
typename ChebyshevOperatorFinite<T>::template Types<V>::VectorOfVectorsVType ChebyshevOperatorFinite<T>::compute_f_1(const V& x) {

    typename Types<V>::VectorOfVectorsVType result(this->n);
    auto c_series_input = compute_c(x);


    // 1. Obliczenie f_1[0]
    for (int i = 0; i < this->n; i++){
        result[i] = typename V::VectorType(this->N);
        result[i][0] = this->u0[i] - this->getCoeff(x, i, 0);
        for (int l = 1; l < this->N; l++){
            if (l % 2 == 0){
                result[i][0] = result[i][0] - 2.0 * this->getCoeff(x, i, l);
            }
            else{
                result[i][0] = result[i][0] + 2.0 * this->getCoeff(x, i, l);
            }
        }
    }

//     2. Obliczanie f_1 dla k > 0
    for (int i = 0; i < this->n; i++){
        for (int k = 1; k < this->N; ++k) {
            result[i][k] = this->getCoeff(x, 0, 0, true) * k * this->getCoeff(x, i, k) - 0.5 * (c_series_input[i][k-1] - c_series_input[i][k+1]);
        }
    }
    return result;
}


