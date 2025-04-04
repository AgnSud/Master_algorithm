#pragma once

#include "ChebyshevOperatorFinite.hpp"
#include "capd/fadbad/fadiff.h"

/**
 * TODO: Teraz zrobie podejscie dla dwoch zmiennych - dla van der Pol i zobacze czy wychodzi, czyli przyjmuje n=2
 * Nie moge zrobic konwersji     capd::vectalg::Vector<T, 0> x_helper(x.dimension()) = x.val...
 * w funkcji convertXVectorToOmegaAndASeries(), bo wtedy nie bedzie sie wyliczała pochodna?
 * dlatego w tej funkcji i innych zmienie a_series, na a_series_first_variable i a_series_second_variable
 * aby nie robic z tego Vector<ChebyshevSeries>, jak to było w poprzedniej wersji.
 * Zmienie tez computeC aby przyjmował (const V& x) i w funkcji juz odpakował. Funkcja computeC będzie również zwracała
 * parę (c_series_first_variable, c_series_second_variable)
 */

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

//template<typename T>
//void ChebyshevOperatorFinite<T>::setCSeries() {
//    this->c_series = computeC();
//}

template<typename T>
void ChebyshevOperatorFinite<T>::setOmega(const T& omega_su) {
    this->omega = omega_su;
}


template <typename T>
template <class V>
std::tuple<typename V::VectorType, typename V::VectorType>
ChebyshevOperatorFinite<T>::computeC(const V& x) {
    auto [omega_input, a_series_input_first_variable, a_series_input_second_variable] = convertXVectorToOmegaAndASeries(x);

    //konwersja a_series na ChebyshevSeries:
    ChebyshevSeries<T, DIMENSION> converted_a_series_input_first_variable(this->N);
    ChebyshevSeries<T, DIMENSION> converted_a_series_input_second_variable(this->N);
    for (int i = 0; i < this->N; i++){
        converted_a_series_input_first_variable[i] = a_series_input_first_variable[i];
        converted_a_series_input_second_variable[i] = a_series_input_second_variable[i];
    }

    // Tworzymy wektory dla c_series
    typename V::VectorType c_series_result_first_variable(this->N + 1);
    typename V::VectorType c_series_result_second_variable(this->N +1);


    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alpha - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        auto exponent = multiIndices[alpha];
        ChebyshevSeries<T, DIMENSION> aj_pow_first_variable = converted_a_series_input_first_variable.power(exponent[0]);
        ChebyshevSeries<T, DIMENSION> aj_pow_second_variable = converted_a_series_input_second_variable.power(exponent[1]);
        ChebyshevSeries<T, DIMENSION> conv = ChebyshevSeries<T>::convolve(aj_pow_first_variable, aj_pow_second_variable);

        // Teraz dodaj g_alpha * term do c
        // Potrzebujemy wielkosc c_{k+1}, więc obetniemy c do N+1
        for (int i = 0; i < this->N + 1; i++){
            c_series_result_first_variable[i] += this->g[alpha][0] * conv[i];
            c_series_result_second_variable[i] += this->g[alpha][1] * conv[i];
        }
    }
    std::cout << "c_series_result_first_variable: " << c_series_result_first_variable << '\n';
    std::cout << "c_series_result_second_variable: " << c_series_result_second_variable << '\n';

    return {c_series_result_first_variable, c_series_result_second_variable};
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
std::tuple<typename V::ScalarType, typename V::VectorType, typename V::VectorType>
ChebyshevOperatorFinite<T>::convertXVectorToOmegaAndASeries(const V& x) {
    // Wyciągamy omega z pierwszego elementu wektora x
    typename V::ScalarType omega_result = x[0];

    // Tworzymy a_series z pozostałych elementów wektora x
    typename V::VectorType a_series_result_first_variable(this->N);
    typename V::VectorType a_series_result_second_variable(this->N);

    int index = 1;
    for (int k = 0; k < this->N; ++k) {
        a_series_result_first_variable[k] = x[index];
        a_series_result_second_variable[k] = x[index+1];
        index += 2;
    }
    return {omega_result, a_series_result_first_variable, a_series_result_second_variable};
}


template<typename T>
capd::vectalg::Vector<T, DIMENSION> ChebyshevOperatorFinite<T>::findFiniteSolution(
        T omega_start,
        const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series_start,
        int max_iterations, T tolerance) {
    // Ustaw a i c
    setOmega(omega_start);
    setASeries(a_series_start);
//    setCSeries();

    int iteration = 0;
    T norm = 1.0;

    // Utwórz wynik jako Vector<T, DIMENSION>
    capd::vectalg::Vector<T, DIMENSION> F(this->n * this->N + 1);
    capd::vectalg::Matrix<T, DIMENSION, DIMENSION> jacobian(this->n * this->N + 1, this->n * this->N + 1);

    //x_k_1 = x_k - F_x_k/DF_x_k
    //co oznacza x_{k+1} = x_k - F(x_k) / DF(x_k)
    while (iteration < max_iterations && norm > tolerance) {
        capd::vectalg::Vector<T, 0> x = convertToXVector();

        std::cout << "Wyliczenia operatora Czebyszewa F:" << '\n';
        std::cout << "a_series: " << this->a_series << '\n';
//        std::cout << "c_series: " << this->c_series << '\n';
        std::cout << "omega: " << this->omega << '\n';
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
capd::vectalg::Vector<T, 0> ChebyshevOperatorFinite<T>::computeF(const capd::vectalg::Vector<T, 0>& x) {
    // Oblicz f_0
    auto f0 = computeF0(x);
    capd::vectalg::Vector<T, 0> result;
    // Oblicz f_1
//    auto [f1_first_variable, f1_second_variable] = computeF1();
//
//    capd::vectalg::Vector<T,0> result = convertToXVector();

//    int index = 1;
//    for (int k = 0; k < this->N; ++k) {
//        result[index]
//        for (int i = 0; i < this->n; ++i) {
//            result[index] = f1[i][k];  // f1 w postaci F
//            index++;
//        }
//    }

    return result;
}

template<typename T>
template<class V>
V ChebyshevOperatorFinite<T>::operator() (const V& x) {
    V result(x.dimension());

    result[0] = computeF0(x);
    auto [f1_first_variable, f1_second_variable] = computeF1(x);
//    std::cout << "f1_first_variable: " << f1_first_variable << '\n';
//    std::cout << "f1_second_variable: " << f1_second_variable << '\n';

    return result;
}



template<typename T>
template<class V>
typename V::ScalarType ChebyshevOperatorFinite<T>::computeF0(const V& x){

    auto [omega_input, a_series_input_first_variable, a_series_input_second_variable] = convertXVectorToOmegaAndASeries(x);
//    <v, u>
    typename V::ScalarType result = v[0] * u[0] + v[1] * u[1];

//    <v, a_0>
    result -= v[0] * a_series_input_first_variable[0] + v[1] * a_series_input_second_variable[0];

//    sum_{k=1}^{N-1} 2<v, a_k>
    for (int k = 1; k < this->N; ++k) {
        result -= 2 * v[0] * a_series_input_first_variable[k] + 2 * v[1] * a_series_input_second_variable[k];
    }

    return result;
}


template<typename T>
template<class V>
std::tuple<typename V::VectorType, typename V::VectorType> ChebyshevOperatorFinite<T>::computeF1(const V& x) {

    typename V::VectorType result_first_variable(this->N);
    typename V::VectorType result_second_variable(this->N);
    auto [omega_input, a_series_input_first_variable, a_series_input_second_variable] = convertXVectorToOmegaAndASeries(x);
    auto [c_series_result_first_variable, c_series_result_second_variable] = computeC(x);


    // 1. Obliczenie f_1[0]
    result_first_variable[0] = this->y0[0] - a_series_input_first_variable[0];
    result_second_variable[0] = this->y0[1] - a_series_input_second_variable[0];

    for (int l = 1; l < this->N; l++){
        if (l % 2 == 0){
            result_first_variable[0] -= 2.0 * a_series_input_first_variable[l];
            result_second_variable[0] -= 2.0 * a_series_input_second_variable[l];
        }
        else{
            result_first_variable[0] += 2.0 * a_series_input_first_variable[l];
            result_second_variable[0] += 2.0 * a_series_input_second_variable[l];
        }
    }

//     2. Obliczanie f_1 dla k > 0
    for (int k = 1; k < this->N; ++k) {
        result_first_variable[k] = omega * k * a_series_input_first_variable[k] -
                0.5 * (c_series_result_first_variable[k-1] - c_series_result_first_variable[k+1]);

        result_second_variable[k] = omega * k * a_series_input_second_variable[k] -
                                   0.5 * (c_series_result_second_variable[k-1] - c_series_result_second_variable[k+1]);
    }


    return {result_first_variable, result_second_variable};
}


