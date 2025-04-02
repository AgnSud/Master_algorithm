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
    ChebyshevOperatorFinite() : N(1), n(1), y0(1), g(0, 2), a_series(1), c_series(1) {}

    /**
     * parametry konstruktora. Pola które nie sa inicjowane w konstruktorze to a_series, c_series oraz omega.
     * Beda one inicjowane oraz aktualizowane o kolejne iteracje w metodzie w metodzie findFiniteSolutionByNewton
     * metoda findFiniteSolution bedzie najbardziej zewnatrzna do wyznaczenia przyblizonego 0 oraz odwrotnosci pochodnej
     * DODATKOWE POLA, ktore jeszcze przewiduje to macierz pochodnej, macierz odwrotnosci pochodnej i przyblizone 0,
     * ktore beda mialy swoje gettery
     *
     * @param N
     * @param n
     * @param y0_init
     * @param g_init
     * @param v
     * @param u
     * @param multiIndices
     */
    ChebyshevOperatorFinite(
            int N, int n,
            const capd::vectalg::Vector<T, 0>& y0_init,
            const capd::vectalg::Matrix<T, 0, 0>& g_init,
            const ChebyshevSeries<T, 0>& v,
            const ChebyshevSeries<T, 0>& u,
            const std::vector<std::vector<int>>& multiIndices
    );

    /**
     * Ustawia wektor szeregów Czebyszewa a_series, który reprezentuje przybliżenie funkcji w przestrzeni współrzędnych.
     *
     * @param a_input Wektor szeregów Czebyszewa.
     */
    void setASeries(const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_input);

    /**
     * Ustawia wektor c_series, który jest wynikiem obliczenia szeregu c(a) dla wielowskaźników i splotu a
     * Funkcja ta oblicza wyniki operacji na szeregach Czebyszewa [a]_k oraz multiindeksach reprezentujacych funkcje g.
     *
     * @param multiIndices Wektor wielowskaźników, na podstawie których obliczane są składniki c_series.
     */
    void setCSeries();

    void setOmega(const T& omega);


    /**
     * Oblicza c_series na podstawie wielowskaźników i szeregów Czebyszewa a_series.
     * Dla każdego wielowskaźnika, obliczany jest wynik splotu odpowiednich szeregów Czebyszewa, a następnie mnożony przez odpowiednią wartość z funkcji g.
     *
     * @param multiIndices Wektor wielowskaźników, który jest używany do obliczenia odpowiednich składników dla c_series.
     * @return Wektor szeregów Czebyszewa c_result, który zawiera wyniki obliczeń.
     */
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> computeC(
            const std::vector<std::vector<int>>& multiIndices
    );

    capd::vectalg::Vector<T, DIMENSION> findFiniteSolution(
            T omega_start,
            const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series_start,
            int max_iterations = 100,
            T tolerance = 1e-8);

    /**
     * Oblicza wartość funkcji F. Funkcja ta wywołuje computeF0 i computeF1, aby obliczyć wartości f0 oraz f1.
     *
     * @return Para, zawierająca wartość f0 oraz wektor f1.
     */
    capd::vectalg::Vector<fadbad::F<T>,0> computeF(const capd::vectalg::Vector<fadbad::F<T>,0>& x);


    /**
     * Oblicza wartość funkcji f0, która jest pierwszym składnikiem funkcji F. Jest to obliczenie iloczynu skalarnego między wektorami v, u oraz szeregami Czebyszewa a_series.
     *
     * @return Wartość f0, która jest wynikiem obliczeń, jednowymiarowa
     */
    T computeF0();

    /**
     * Oblicza wartość funkcji f1, która jest drugim składnikiem funkcji F według wzoru z omega, a_series, c_series
     *
     * @param omega Wartość parametru omega, który jest używany w obliczeniu funkcji f1.
     * @return Wektor szeregów Czebyszewa f1, który zawiera wyniki obliczeń dla funkcji f1.
     */
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> computeF1();

private:
    int N; //truncation - tyle elementow bedzie mial kazdy ChebyshevSeries
    int n; //wymiar - liczba zmiennych
    capd::vectalg::Vector<T, 0> y0;
    capd::vectalg::Matrix<T, 0, 0> g;
    ChebyshevSeries<T, 0> v;
    ChebyshevSeries<T, 0> u;
    std::vector<std::vector<int>> multiIndices;

    T omega;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> a_series;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> c_series;
};

#include "ChebyshevOperatorFinite.tpp"
