#pragma once

#include "ChebyshevOperatorFinite.hpp"
#include "capd/fadbad/fadbad.h"
#include <type_traits>
#include <iomanip>

using namespace capd;
using namespace std;

template<typename T>
ChebyshevOperatorFinite<T>::ChebyshevOperatorFinite(
        const int N, const int n,
        const VectorType& u0_init,
        const MatrixType& g_init,
        const VectorType& v,
        const VectorType& w,
        const vector<vector<int>>& multiIndices
) : N(N), n(n), omega(0), u0(u0_init), g(g_init), a_series(n), c_series(n), v(v), w(w), multiIndices(multiIndices), F_x_approx(N*n+1), x_approx(N*n+1) {}

template<typename T>
void ChebyshevOperatorFinite<T>::setASeries(const VectorOfChebyshevsType& a_input) {
    this->a_series = a_input;
}
template <typename T>
ChebyshevOperatorFinite<T>::VectorOfChebyshevsType ChebyshevOperatorFinite<T>::getASeries() const {
    return this->a_series;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setACoeff(T coeff, int j, int k) {
    this->a_series[j][k] = coeff;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setOmega(const T& omega_su) {
    this->omega = omega_su;
}
template <typename T>
T ChebyshevOperatorFinite<T>::getOmega() const {
    return this->omega;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setF_x_approx(const VectorType& F_x_approx_) {
    this->F_x_approx = F_x_approx_;
}
template <typename T>
ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::getF_x_approx() const {
    return this->F_x_approx;
}

template<typename T>
void ChebyshevOperatorFinite<T>::setX_approx(const VectorType& x_approx) {
    this->x_approx = x_approx;
}
template <typename T>
ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::getX_approx() const {
    return this->x_approx;
}

template <typename T>
ChebyshevOperatorFinite<T>::MatrixType ChebyshevOperatorFinite<T>::getG() const {
    return this->g;
}

template <typename T>
vector<vector<int>> ChebyshevOperatorFinite<T>::getMultiIndices() const {
    return this->multiIndices;
}


template <typename T>
void ChebyshevOperatorFinite<T>::setCSeries(const VectorOfChebyshevsType& c_input) {
    this->c_series = c_input;
}

template <typename T>
ChebyshevOperatorFinite<T>::VectorOfChebyshevsType ChebyshevOperatorFinite<T>::getCSeries() const {
    return this->c_series;
}

template <typename T>
void ChebyshevOperatorFinite<T>::setInverseDerivativeFinite(const MatrixType& inverse_derivative_finite) {
    this->inverse_derivative_finite = inverse_derivative_finite;
}

template <typename T>
ChebyshevOperatorFinite<T>::MatrixType ChebyshevOperatorFinite<T>::getInverseDerivativeFinite() const {
    return this->inverse_derivative_finite;
}

template <typename T>
void ChebyshevOperatorFinite<T>::setDerivativeFinite(const MatrixType& derivative) {
    this->derivative_finite = derivative;
}

template <typename T>
ChebyshevOperatorFinite<T>::MatrixType ChebyshevOperatorFinite<T>::getDerivativeFinite() const {
    return this->derivative_finite;
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
V ChebyshevOperatorFinite<T>::compute_c(const V& x) {
//    V c_series_result_flatten(this->n * (this->N+1) + 1);
    // TODO: zmiana wynikowego rozmiaru c, bez obcinania po N+1, czy mozna to jakos inaczej zainicjowac, nie na etapie kompliacji?
    V c_series_result_flatten(n * (2 * N -1) + 1);
    c_series_result_flatten[0] = 0;          //ustawienie domyślne dla przesunięc

    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alpha - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        V conv(1);  // stały 1
        conv[0] = 1.0;

//        cout << "multiIndices[alpha] = " << multiIndices[alpha] << endl;
        for (int alpha_idx = 0; alpha_idx < this->n; ++alpha_idx) {
            int exponent = multiIndices[alpha][alpha_idx];
            //compute power:
            V a_square_paran_alpha_idx = getCoeffVectorI_thSquareParan(x, alpha_idx, this->N, this->n);
            V aj_pow(1);

            aj_pow[0] = 1;
            for (int i = 0; i < exponent; ++i) {
//                aj_pow = aj_pow * a_square_paran_alpha_idx;
                aj_pow = multiply(aj_pow, a_square_paran_alpha_idx);
            }
            conv = ChebyshevSeries<T>::convolve(conv, aj_pow);
//            if constexpr (is_same<decltype(conv[0]), double&>::value){
//                cout << "aj_pow: " << aj_pow << '\n';
//                cout << "conv: " << conv << '\n';
//            }
        }

        // Teraz dodaj g_alpha * term do c
        // Potrzebujemy wielkosc c_{k+1}, więc obetniemy c do N+1
        // TODO: nie moge obcinac, bo poxniej w wyliczeniach Y jest to istotne, wiec obliczam do teoretycznego wymiaru c
        int index = 1;
        for (int i = 0; i < this->n; ++i) {
            for (int k = 0; k < conv.dimension(); ++k) {
                c_series_result_flatten[1 + k * this->n + i] += this->g[alpha][i] * conv[k];
            }
        }
    }


//    if constexpr (is_same<decltype(c_series_result_flatten[0]), double&>::value)
//        cout << "c_series_result_flatten: " << c_series_result_flatten << '\n';

    return c_series_result_flatten;
}

template<typename T>
typename ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::convertToXVector() {
    VectorType x(this->n * this->N + 1);
    x[0] = this->omega;

    int index = 1;
    for (int k = 0; k < this->N; ++k) {
        for (int i = 0; i < this->n; ++i) {
            x[index] = this->a_series[i][k];
            index++;
        }
    }
    return x;
}

template<typename T>
template<class V>
inline typename V::ScalarType ChebyshevOperatorFinite<T>::getCoeff(const V &x, int i, int k, int n, bool is_omega) {
    if (is_omega)
        return x[0];
    return x[1 + k * n + i];
}

template<typename T>
template<class V>
inline V ChebyshevOperatorFinite<T>::getCoeffVectorI_thSquareParan(const V& x, int i, int size, int n) {
    V result(size);
    for (int k = 0; k < size; ++k) {
        result[k] = x[1 + k * n + i];
    }
    return result;
}

//size się różni, bo czasem wywoluje na x, a czasem na c
template<typename T>
template<class V>
typename ChebyshevOperatorFinite<T>::VectorOfChebyshevsType ChebyshevOperatorFinite<T>::convertToSeriesFromXForm(const V& x, int size){
    VectorOfChebyshevsType series_form(this->n);
    for (int i = 0; i < this->n; i++){
        ChebyshevSeries<T, DIMENSION> tmp(size);
        tmp.setCoefficients(getCoeffVectorI_thSquareParan(x, i, size, this->n));
//        if (size > this->N){
//            cout << "c_flatten[i] = " << tmp << endl;
//        }
        series_form[i] = tmp;
    }
    return series_form;
}


template<typename T>
std::pair<T, typename ChebyshevOperatorFinite<T>::VectorOfChebyshevsType> ChebyshevOperatorFinite<T>::findFiniteSolution(
        T omega_start,
        const VectorOfChebyshevsType& a_series_start,
        int max_iterations, T tolerance) {
    setOmega(omega_start);
    setASeries(a_series_start);
    int iteration = 0;

    VectorType F_x_k(this->n * this->N + 1);
    MatrixType jacobian(this->n * this->N + 1, this->n * this->N + 1);
    auto x = convertToXVector();

    NormType myNorm;
    T norm_tolerance = myNorm(x);

    while (iteration < max_iterations && norm_tolerance > tolerance) {
        F_x_k = (*this)(x);
        computeDerivative(*this, x, jacobian);
        x = x - matrixAlgorithms::gauss(jacobian, F_x_k);
        norm_tolerance = myNorm(F_x_k);
        iteration++;
    }
    VectorOfChebyshevsType a_series_final = convertToSeriesFromXForm(x, this->N);
    VectorOfChebyshevsType c_series_final = convertToSeriesFromXForm(compute_c(x), 2 * this->N - 1);
    T omega_final = getCoeff(x, 0, 0, this->n, true);

    this->setASeries(a_series_final);
    this->setCSeries(c_series_final);
    this->setOmega(omega_final);
    this->setF_x_approx((*this)(x));
    this->setX_approx(x);
    cout << "Liczba iteracji Newtona dla Czebyszewa: " << iteration << endl;
    cout << "Czas przejścia: " << 1./omega_final << endl;
    //obliczenie odwrotnosci jacobianu
//    cout << "tutaj" << endl;
    computeDerivativeInverse(x);

    return std::make_pair(omega_final, a_series_final);
}

template <typename T>
typename ChebyshevOperatorFinite<T>::VectorType ChebyshevOperatorFinite<T>::NewtonLikeOperatorTx_x(const VectorType& x) {
    /// poniżej test wykorzystujacy wczesniejsze T(x) - x = x - AF(x) - x = -AF(x) - czyli roznica ktora tak naprawde sprawdzamy
    auto F_x = (*this)(x);
    VectorType mult_A_F_x = getInverseDerivativeFinite() * F_x;
    VectorType result = mult_A_F_x;
    return result;
}

template <typename T>
typename ChebyshevOperatorFinite<T>::VectorType
ChebyshevOperatorFinite<T>::applyDT(const VectorType& x, const VectorType& x2) {
    using MatrixType = typename ChebyshevOperatorFinite<T>::MatrixType;

    // DT(x) = I - A * DF(x) - np x = x^* + x1
    MatrixType DF(this->n * this->N + 1, this->n * this->N + 1);;
    computeDerivative(*this, x, DF);

    MatrixType A = getInverseDerivativeFinite();

    // DT = I - A * DF
    MatrixType I(this->n * this->N + 1, this->n * this->N + 1);
    I.setToIdentity();
    MatrixType DT = I - A * DF;

    return DT * x2;
}

template <typename T>
typename ChebyshevOperatorFinite<T>::VectorType
ChebyshevOperatorFinite<T>::applyDT_decomposed(const VectorType& x, const VectorType& x2) {
    // DF(𝑥̂ + x₁)
    MatrixType DF(this->n * this->N + 1, this->n * this->N + 1);
    computeDerivative(*this, x, DF);  // DF(x̂ + x₁)

    // A = [DF(x̂)]^{-1}, hatA = DF(x̂)
    MatrixType A = getInverseDerivativeFinite();
    MatrixType A_hat = getDerivativeFinite();

    // część 1: (I - A A_hat) x2
    MatrixType I(this->n * this->N + 1, this->n * this->N + 1);
    I.setToIdentity();
    VectorType first_term = (I - A * A_hat) * x2;

    // część 2: A (DF(x) - A_hat) x2
    VectorType second_term = A * ((DF - A_hat) * x2);

    // końcowy wynik
    return first_term - second_term;
}



template <typename T>
template<class V>
void ChebyshevOperatorFinite<T>::computeDerivativeInverse(const V& x) {
    // Obliczamy macierz Jacobiego w punkcie x
    MatrixType jacobian(this->n * this->N + 1, this->n * this->N + 1);
    computeDerivative(*this, x, jacobian);
    MatrixType jacobianInversed = matrixAlgorithms::gaussInverseMatrix(jacobian);
    this->setDerivativeFinite(jacobian);
    this->setInverseDerivativeFinite(jacobianInversed);
}

template<typename T>
template<class V>
V ChebyshevOperatorFinite<T>::operator() (const V& x) {
    V result(x.dimension());

    result[0] = compute_f_0(x);
    auto f1_result = compute_f_1(x);
    for (int i = 0; i < f1_result.dimension(); i++){
        result[i + 1] = f1_result[i];
    }

//    if constexpr (is_same<decltype(f1_result[0]), double&>::value) {
//        cout << "f1_result: " << f1_result << '\n';
//        cout << "result: " << result << '\n';
//    }

    return result;
}

template<typename T>
template<class V>
typename V::ScalarType ChebyshevOperatorFinite<T>::compute_f_0(const V& x){
    typename V::ScalarType result = 0;
    for (int i = 0; i < this->n; i++){
        // <v, w>
        result += v[i] * w[i];
//        if constexpr (std::is_same_v<typename V::ScalarType, double> ||
//                      std::is_same_v<typename V::ScalarType, capd::intervals::Interval<double>>){
//            std::cout << "f_0 result: " << result << ", v: " << v[i] << ", w: " << w[i] << std::endl;
//        }

        // <v, a_0>
        result -= v[i] * this->getCoeff(x, i, 0, this->n);
//        if constexpr (std::is_same_v<typename V::ScalarType, double> ||
//                      std::is_same_v<typename V::ScalarType, capd::intervals::Interval<double>>){
//            std::cout << "f_0 result: " << result << ", v: " << v[i] << ", a_0: " << this->getCoeff(x, i, 0, this->n) << std::endl;
//        }

        // sum_{k=1}^{N-1} 2<v, a_k>
        for (int k = 1; k < this->N; ++k) {
            result -= 2 * v[i] * this->getCoeff(x, i, k, this->n);
//            if constexpr (std::is_same_v<typename V::ScalarType, double> ||
//                          std::is_same_v<typename V::ScalarType, capd::intervals::Interval<double>>){
//                std::cout << "f_0 result: " << result << ", v: " << v[i] << ", a_k: " << this->getCoeff(x, i, k, this->n) << std::endl;
//            }
        }
    }
//    if constexpr (std::is_same_v<typename V::ScalarType, double> ||
//                  std::is_same_v<typename V::ScalarType, capd::intervals::Interval<double>>){
//        std::cout << "f_0 result: " << result << std::endl;
//    }
    return result;
}


template<typename T>
template<class V>
V ChebyshevOperatorFinite<T>::compute_f_1(const V& x) {
    V result(this->n * this->N);
    auto c = compute_c(x);
//    if constexpr (is_same<decltype(result[0]), double&>::value) {
//        cout << "x: " << x << '\n';
//        cout << "c: " << c << '\n';
//    }


    // 1. Obliczenie f_1[0]
    for (int i = 0; i < this->n; i++){
        result[i] = this->u0[i] - this->getCoeff(x, i, 0, this->n);
        for (int k = 1; k < this->N; k++){
            if (k % 2 == 0){
                result[i] = result[i] - 2.0 * this->getCoeff(x, i, k, this->n);
            }
            else{
                result[i] = result[i] + 2.0 * this->getCoeff(x, i, k, this->n);
            }
        }
    }
//    if constexpr (is_same<decltype(result[0]), double&>::value) {
//        cout << "result: " << result << '\n';
//    }

//     2. Obliczanie f_1 dla k > 0
    for (int i = 0; i < this->n; i++){
        for (int k = 1; k < this->N; ++k) {
            result[k * this->n + i] = this->getCoeff(x, 0, 0, this->n, true) * k * this->getCoeff(x, i, k, this->n)
                    - 0.25 * (this->getCoeff(c, i, k - 1, this->n) -this->getCoeff(c, i, k + 1, this->n));
        }
    }
    // TODO: jak jest 0.25 to jest po prostu szybciej - różnica w czasie, zapytać o interpretację tego czasu
    return result;
}


