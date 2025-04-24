#pragma once

#include <capd/vectalg/lib.h>
#include <capd/fadbad/differentiate.h>
#include <cmath>
#include "ChebyshevSeries.hpp"
#include "norm.hpp"


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
    typedef vectalg::Vector<T, DIMENSION> VectorType;
    typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;
    template <class V>
    struct Types {
        typedef vectalg::Vector<V, DIMENSION> VectorOfVType;
    };

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
    void setOmega(const T& omega);
    void setCSeries(const VectorOfChebyshevsType& c_input);
    VectorOfChebyshevsType getCSeries() const;
    void setDerivativeFinite(const MatrixType& derivative);
    MatrixType getDerivativeFinite() const;

    /**
     * Oblicza przybliżone rozwiązanie w metodzie Newtona.
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
            T tolerance = 1e-8);

    // wylicza operator Czebyszewa skonczony F_N
    template<class V>
    V operator() (const V& x);

private:
    int N;
    int n;
    VectorType u0;
    MatrixType g;
    ChebyshevSeries<T, 0> v;
    ChebyshevSeries<T, 0> w;
    vector<vector<int>> multiIndices;
    MatrixType derivative_finite;

    T omega;
    VectorOfChebyshevsType a_series;
    VectorOfChebyshevsType c_series;

    template <class V>
    V compute_c(const V& x);

    template<class V>
    typename V::ScalarType compute_f_0(const V& x);

    template<class V>
    V compute_f_1(const V& x);

    template<class V>
    VectorOfChebyshevsType convertToSeriesFromXForm(const V& x, int size);
    VectorType convertToXVector();

    template<class V>
    inline typename V::ScalarType getCoeff(const V &x, int i, int k, bool is_omega=false) const;

    template<class V>
    inline V getCoeffVectorI_thSquareParan(const V &x, int i, int size) const;

    template<class V>
    V multiply(const V& a, const V& b);

    template<class V>
    void computeDerivativeInverse(const V& x);
};

#include "ChebyshevOperatorFinite.tpp"
