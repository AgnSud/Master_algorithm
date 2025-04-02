#pragma once

//#include "ChebyshevOperatorFinite.hpp"

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
void ChebyshevOperatorFinite<T>::setCSeries() {
    this->c_series = computeC(multiIndices);
}

template<typename T>
void ChebyshevOperatorFinite<T>::setOmega(const T& omega_su) {
    this->omega = omega_su;
}


template <typename T>
capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION>
ChebyshevOperatorFinite<T>::computeC(const std::vector<std::vector<int>>& multiIndices) {
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> c_result(this->n);

    // kazdy element c_result (czyli wektor) to będzie ChebyshevSeries i ich będzie tyle, jaki wymiar ma a (małe n),
    // ale kazdy z nich bedzie wymiaru 2 * N - 1 (bo są wynikiem convolve współczynników a, które są wielkości N)
    for (int k = 0; k < this->n; ++k)
        c_result[k] = ChebyshevSeries<T, DIMENSION>(this->N + 1);


    int numAlphas = multiIndices.size();  // A
    // Oblicz [a]_1^{alpha_1} * [a]_2^{alpha_2} * ... * [a]_n^{alpha_n}
    // alhpa - numer wielowskaznika
    // alpha_idx - elementy wielowskaznika
    for (int alpha = 0; alpha < numAlphas; ++alpha) {
        ChebyshevSeries<T, DIMENSION> conv(this->N + 1);  // stały 1
        conv[0] = 1.0;


        for (int alpha_idx = 0; alpha_idx < this->n; ++alpha_idx) {
            int exponent = multiIndices[alpha][alpha_idx];
            ChebyshevSeries<T, DIMENSION> aj_pow = this->a_series[alpha_idx].power(exponent);
            conv = ChebyshevSeries<T>::convolve(conv, aj_pow);
//            std::cout << "Iteracja " << alpha_idx << ": ";
//            std::cout << "a_series[" << alpha_idx << "] ^" << exponent << ": " << aj_pow;
//            std::cout << conv;
        }
        //--------------------------------------------------------------------

        // Teraz dodaj g_alpha * term do c
        // Czy tutaj max() czy min()? Czy c będzie mogło miec wiecej wspolczynnikow od a?
        // Znaczy jako wynik splotu sam w sobie to tak, pytanie czy tego chcemy, czy c tez obcinamy do N?
        // Potrzebujemy wielkosc c_{k+1}, więc obetniemy c do N+1
        // Robie tak, bo mam mnozenie g[alpha] (vector) * conv (ChebyshevSeries)
        for (int k = 0; k < this->n; ++k) {
            for (int d = 0; d <= this->N; ++d) {
//                std::cout << "g[alpha] = " << this->g[alpha] << "\n";
//                std::cout << "conv[d] =" << conv[d] << "\n";
//                std::cout << "k=" << k << ", d=" << d << "\n";
//                std::cout << "g[alpha][k] * conv[d] = " << this->g[alpha][k] * conv[d] << "\n\n";
                c_result[k][d] += this->g[alpha][k] * conv[d];
            }
        }
    }

    return c_result;
}



template<typename T>
capd::vectalg::Vector<T, DIMENSION> ChebyshevOperatorFinite<T>::findFiniteSolution(
        T omega_start,
        const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series_start,
        int max_iterations, T tolerance) {

    // Ustaw a i c
    setOmega(omega_start);
    setASeries(a_series_start);
    setCSeries();

    int iteration = 0;
    T norm = 1.0;

    // Utwórz wynik jako Vector<T, DIMENSION>
    capd::vectalg::Vector<T, DIMENSION> F(this->n * this->N + 1);

    while (iteration < max_iterations && norm > tolerance) {
        auto [f0, f1] = this->computeF();
        std::cout << "Wyliczenia operatora Czebyszewa F:" << '\n';
        std::cout << "a_series: " << this->a_series << '\n';
        std::cout << "c_series: " << this->c_series << '\n';
        std::cout << "omega: " << this->omega << '\n';
        std::cout << "f0: " << f0 << '\n';
        std::cout << "f1: " << f1 << '\n';

        // Zamiana na jeden wektor F
        F[0] = f0;
        int index = 1;
        for (int k = 0; k < this->N; ++k) {
            for (int i = 0; i < this->n; ++i) {
                F[index] = f1[i][k];  // f1 jest wektorem szeregów Czebyszewa
                index++;
            }
        }

        //Zamiana na jeden punkt x = (\omega, a)
        capd::vectalg::Vector<T, DIMENSION> x(this->n * this->N + 1);
        x[0] = omega;

        index = 1;
        for (int k = 0; k < this->N; ++k) {
            for (int i = 0; i < this->n; ++i) {
                x[index] = a_series[i][k];  // Wypełniamy wartościami z a_series
                index++;
            }
        }
        std::cout << "x=(omega, a): " << x << "\n";
        std::cout << "F(x): " << F << "\n";

        setOmega(f0);
        setASeries(f1);
        setCSeries();

//        std::cout << "Po zaktualizowaniu:" << '\n';
//        std::cout << "a_series: " << this->a_series << '\n';
//        std::cout << "c_series: " << this->c_series << '\n';
//        std::cout << "omega: " << this->omega << '\n';
        iteration++;
    }

    return F;
}


template<typename T>
std::pair<T, capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>>
ChebyshevOperatorFinite<T>::computeF() {


    // Oblicz f_0
    T f0 = this->computeF0();

    // Oblicz f_1
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> f1 = this->computeF1();

    return std::make_pair(f0, f1);
}



template<typename T>
T ChebyshevOperatorFinite<T>::computeF0() {

    //skoro i tak potrzebuje petli, aby wyciagnac a_k, to nie bede robic osobnej funkcji dot, tylko tutaj wykonam iloczyn skalarny
    T result = 0;
    for (int i = 0; i < this->n; i++){
        // <v, u>
//        std::cout << "i=" << i << '\n';
//        std::cout << "v[i] * u[i] = " << v[i] * u[i] << '\n';
        result += v[i] * u[i];

        // <v, a_0>
//        std::cout << "v[i] * a_series[i][0] = " << v[i] * this->a_series[i][0] << '\n';
        result -= v[i] * this->a_series[i][0];

        // sum_{k=1}^{N-1} 2<v, a_k>
        for (int k = 1; k < this->N; ++k) {
//            std::cout << "v[i] * a_series[i][k] = " << v[i] * this->a_series[i][k] << '\n';
            result -= v[i] * 2 * this->a_series[i][k];
        }
    }
    return result;
}

template<typename T>
capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> ChebyshevOperatorFinite<T>::computeF1() {

    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> a_series_next(this->n);

    // 1. Obliczenie f_1[0]
    for (int i = 0; i < this->n; i++){
        a_series_next[i] = ChebyshevSeries<T, DIMENSION>(this->N);
        a_series_next[i][0] = this->y0[i] - this->a_series[i][0];
//        std::cout << "i=" << i << '\n';
//        std::cout << "a_series_next[i][0] = " << a_series_next[i][0] << '\n';
        for (int l = 1; l < this->N; l++){
            if (l % 2 == 0){
                a_series_next[i][0] = a_series_next[i][0] - 2.0 * this->a_series[i][l];
            }
            else{
                a_series_next[i][0] = a_series_next[i][0] + 2.0 * this->a_series[i][l];
            }
        }
    }

//     2. Obliczanie f_1 dla k > 0
    for (int i = 0; i < this->n; i++){
        for (int k = 1; k < this->N; ++k) {
            a_series_next[i][k] = omega * k * this->a_series[i][k] - 0.5 * (this->c_series[i][k-1] - this->c_series[i][k+1]);
        }
    }
    return a_series_next;
}