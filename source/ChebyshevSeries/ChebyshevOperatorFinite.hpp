#pragma once

#include <capd/vectalg/lib.h>
#include <capd/fadbad/differentiate.h>
#include <cmath>
#include "ChebyshevSeries.hpp"

#ifndef TRUNCATION
#define TRUNCATION 10
#endif

template <typename T>
class ChebyshevOperatorFinite {
public:
    ChebyshevOperatorFinite() : N(1), n(1), omega(1), y0(1), g(1, 1), v(1), u(1), multiIndices(1), a_series(1), c_series(1) {}
    typedef  capd::vectalg::Matrix<T, DIMENSION, DIMENSION> MatrixType;
    typedef capd::vectalg::Vector<T, DIMENSION> VectorType;
//    TODO

    /**
     * parametry konstruktora.
     * Beda one inicjowane oraz aktualizowane o kolejne iteracje w metodzie w metodzie findFiniteSolutionByNewton
     * metoda findFiniteSolution bedzie najbardziej zewnatrzna do wyznaczenia przyblizonego 0 oraz odwrotnosci pochodnej
     * DODATKOWE POLA, ktore jeszcze przewiduje to macierz pochodnej, macierz odwrotnosci pochodnej i przyblizone 0,
     * ktore beda mialy swoje gettery
     *
     * @param N liczba elementów w każdym szereg Czebyszewa
     * @param n wymiar przestrzeni - liczba zmiennych
     * @param y0_init wektor początkowy
     * @param g_init funkcja w postaci macierzowej
     * @param v fixed vector v (rozmiaru n)
     * @param u fixed vector u (rozmiaru n)
     * @param multiIndices lista wielowskaźników
     */
    ChebyshevOperatorFinite(
            int N, int n,
            const capd::vectalg::Vector<T, 0>& y0_init,
            const capd::vectalg::Matrix<T, 0, 0>& g_init,
            const ChebyshevSeries<T, 0>& v,
            const ChebyshevSeries<T, 0>& u,
            const std::vector<std::vector<int>>& multiIndices
    );

    void setASeries(const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_input);
    void setOmega(const T& omega);

    /**
     * Oblicza przybliżone rozwiązanie w metodzie Newtona.
     *
     * @param omega_start Początkowa wartość omega.
     * @param a_series_start Początkowy wektor szeregów Czebyszewa a_series.
     * @param max_iterations Maksymalna liczba iteracji.
     * @param tolerance Tolerancja obliczeniowa.
     * @return Wektor wynikowy.
     */
    capd::vectalg::Vector<T, DIMENSION> findFiniteSolution(
            T omega_start,
            const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series_start,
            int max_iterations = 100,
            T tolerance = 1e-8);

    template<class V>
    V operator() (const V& x);

private:
    int N;
    int n;
    capd::vectalg::Vector<T, 0> y0;
    capd::vectalg::Matrix<T, 0, 0> g;
    ChebyshevSeries<T, 0> v;
    ChebyshevSeries<T, 0> u;
    std::vector<std::vector<int>> multiIndices;

    T omega;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> a_series;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> c_series;

    template <class V>
    capd::vectalg::Vector<typename V::VectorType, 0> computeC(const V& x);
    template<class V>
    typename V::ScalarType computeF0(const V& x);
    template<class V>
    capd::vectalg::Vector<typename V::VectorType, 0> computeF1(const V& x);

    capd::vectalg::Vector<T, 0> convertToXVector();
    template <class V>
    std::tuple<typename V::ScalarType, capd::vectalg::Vector<typename V::VectorType, 0>> convertXVectorToOmegaAndASeries(const V& x);

    template<class V>
    inline typename V::ScalarType getCoeff(const V &x, int i, int k, bool is_omega=false) const;

    template<class V>
    inline typename V::VectorType getCoeffVectorI_thSquareParan(const V &x, int i) const;

    template<class V>
    V convolve(const V& a, const V& b);

    template<class V>
    V multiply(const V& a, const V& b);
};

//#include "ChebyshevOperatorFinite.old.old.old.tpp"
#include "ChebyshevOperatorFinite.tpp"
//#include "ChebyshevOperatorFiniteCopFor2Variables.tpp"
