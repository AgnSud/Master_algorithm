#pragma once

#include <capd/vectalg/lib.h>
#include <capd/fadbad/differentiate.h>
#include <cmath>
#include <iostream>
#include "ChebyshevSeries.hpp"
#include "Norm.hpp"


#ifndef TRUNCATION
#define TRUNCATION 10
#endif

using namespace capd;
using namespace std;

template <typename T>
class ChebyshevOperatorFinite {
public:
    ChebyshevOperatorFinite() : N(1), n(1), omega(1), u0(1), g(1, 1), v(1), w(1), multiIndices(1), a_series(1), c_series(1) {}
    typedef  vectalg::Matrix<T, DIMENSION, DIMENSION> MatrixType;
//    typedef  vectalg::Matrix<double, DIMENSION, DIMENSION> DMatrixType;
    typedef vectalg::Vector<T, DIMENSION> VectorType;
//    typedef vectalg::Vector<double, DIMENSION> DVectorType;
    typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;
    template <class V>
    struct Types {
        typedef vectalg::Vector<V, DIMENSION> VectorOfVType;
    };
    typedef vectalg::SumNorm<VectorType, MatrixType> NormType;

    /**
     * parametry konstruktora.
     * Beda one inicjowane oraz aktualizowane o kolejne iteracje w metodzie w metodzie findFiniteSolutionByNewton
     * metoda findFiniteSolution bedzie najbardziej zewnatrzna do wyznaczenia przyblizonego 0 oraz odwrotnosci pochodnej
     * DODATKOWE POLA, ktore jeszcze przewiduje to macierz pochodnej, macierz odwrotnosci pochodnej i przyblizone 0,
     * ktore beda mialy swoje gettery
     *
     * @param N liczba elementów w każdym szereg Czebyszewa
     * @param n wymiar przestrzeni - liczba zmiennych
     * @param u0_init wektor początkowy
     * @param g_init funkcja w postaci macierzowej
     * @param v fixed vector v (rozmiaru n)
     * @param w fixed vector w (rozmiaru n)
     * @param multiIndices lista wielowskaźników
     */
    ChebyshevOperatorFinite(
            int N, int n,
            const VectorType& u0_init,
            const MatrixType& g_init,
            const ChebyshevSeries<T, 0>& v,
            const ChebyshevSeries<T, 0>& w,
            const vector<vector<int>>& multiIndices
    );

    void setASeries(const VectorOfChebyshevsType& a_input);
    VectorOfChebyshevsType getASeries() const;
    void setOmega(const T& omega);
    T getOmega() const;
    void setF_x_approx(const VectorType& F_x_approx);
    VectorType getF_x_approx() const;
    void setX_approx(const VectorType& x_approx);
    VectorType getX_approx() const;

    MatrixType getG() const;
    vector<vector<int>> getMultiIndices() const;

    void setCSeries(const VectorOfChebyshevsType& c_input);
    VectorOfChebyshevsType getCSeries() const;
    void setInverseDerivativeFinite(const MatrixType& derivative);
    MatrixType getInverseDerivativeFinite() const;
    void setDerivativeFinite(const MatrixType& derivative);
    MatrixType getDerivativeFinite() const;

    /**
     * dziala jak operator Czebyszewa F z pracy - ma przyjmować (omega, a_series) i zwracać (omega, a_series) -
     * wtedy jest zgodnośc z praca
     * drugie podejście, to żeby zwracła void i były gettery do a_series i omega
     *
     * @param omega_start Początkowa wartość omega.
     * @param a_series_start Początkowy wektor szeregów Czebyszewa a_series.
     * @param max_iterations Maksymalna liczba iteracji.
     * @param tolerance Tolerancja obliczeniowa.
     * @return Wektor wynikowy.
     */
    std::pair<T, typename ChebyshevOperatorFinite<T>::VectorOfChebyshevsType> findFiniteSolution(
            T omega_start,
            const VectorOfChebyshevsType& a_series_start,
            int max_iterations = 100,
            T tolerance = 1e-13);

    VectorType NewtonLikeOperatorTx_x(const VectorType& x);
    VectorType applyDT(const VectorType& x, const VectorType& direction);
    VectorType applyDT_decomposed(const VectorType& x, const VectorType& x2);

    // wylicza operator Czebyszewa skonczony F_N
    //ten operator() przyjmuje x i zwraca x (wewnetrzna funkcja)
    template<class V>
    V operator() (const V& x);

    template <class V>
    V compute_c(const V& x);

    /*
     * Konwersja jest następująca: omega typu T, a_series = [ [a]_1, [a]_2, ..., [a]_n ]
     * do tego [a]_k = [ [a_0]_k, [a_1]_k, [a_2]_k, ..., [a_{N-1}]_k ]
     * zamienia się na x = [omega, [a_0]_1, [a_0]_2,  ..., [a_0]_n, [a_1]_1, [a_1]_2, ...]
     */
    VectorType convertToXVector();
    template<class V>
    VectorOfChebyshevsType convertToSeriesFromXForm(const V& x, int size);

    //zwraca [a_k]_i lub [c_k]_i (ustawione c[0] = 0 dla przesunięcia)
    template<class V>
    static inline typename V::ScalarType getCoeff(const V &x, int i, int k, int n, bool is_omega = false);

    //zwraca [a]_i lub [c]_i (ustawione c[0] = 0 dla przesunięcia)
    template<class V>
    static inline V getCoeffVectorI_thSquareParan(const V &x, int i, int size, int n);

    template<class V>
    void computeDerivativeInverse(const V& x);

    friend std::ostream& operator<<(std::ostream& os, const ChebyshevOperatorFinite<T>& op) {
        os << "\n========= ChebyshevOperatorFinite =========\n";
        os << "Omega: " << op.omega << "\n";

        os << "\nA Series:\n";
        os << op.a_series << "\n";

        os << "\nC Series:\n";
        os << op.c_series << "\n";

        os << "\nF_x_approx:\n";
        os << op.F_x_approx << "\n";

        os << "\nDerivative Finite:\n";
        os << op.derivative_finite << "\n";

        os << "\nInverse Derivative Finite:\n";
        os << op.inverse_derivative_finite << "\n";

        os << "===========================================\n";

        return os;    }

private:
    // TODO: Czy wyciagnac tez do zmiennej F_x_approx?
    int N;
    int n;
    VectorType u0;
    MatrixType g;
    ChebyshevSeries<T, 0> v;
    ChebyshevSeries<T, 0> w;
    vector<vector<int>> multiIndices;
    MatrixType derivative_finite;
    MatrixType inverse_derivative_finite;
    VectorType F_x_approx;
    // TODO: dodałam jako dodatkowe pole, bo w sumie uzywam naprzemiennie? Tylko co jest bardziej kosztowne, trzymanie obu form rozwiazania czy konwersja
    VectorType x_approx;

    T omega;
    // prev:    przemyslec czy nie zmienic na po prostu Vector Vectorow, bo  juz nie wiem kiedy uzywam ChebyshevSeries faktycznie, a najwazniejsze to jest do mnozenia w zasadzie
    // TODO: a_series MUSI być DVectorOfChebyshevsType, bo to później jest wykorzystywane w sprawdzeniu wartosci
    VectorOfChebyshevsType a_series;
    VectorOfChebyshevsType c_series;



    template<class V>
    typename V::ScalarType compute_f_0(const V& x);

    template<class V>
    V compute_f_1(const V& x);


    template<class V>
    V multiply(const V& a, const V& b);


};

#include "ChebyshevOperatorFinite.tpp"
