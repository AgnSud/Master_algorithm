#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"

#ifndef TRUNCATION
#define TRUNCATION 10
#endif

template <typename T>
class ChebyshevOperatorFinite {
public:
    ChebyshevOperatorFinite() : N(1), n(1), y0(1), g(0, 2), a_series(1), c_series(1) {}

    ChebyshevOperatorFinite(
            int N, int n,
            const capd::vectalg::Vector<T, 0>& y0_init,
            const capd::vectalg::Matrix<T, 0, 0>& g_init
    );

    void setASeries(const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_input);
    void setCSeries(const std::vector<std::vector<int>>& multiIndices);

    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> computeC(
            const std::vector<std::vector<int>>& multiIndices
    );

    T computeF0(
            const ChebyshevSeries<T, DIMENSION>& v,
            const ChebyshevSeries<T, DIMENSION>& u
    );

    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> computeF1(T omega, const std::vector<std::vector<int>>& multiIndices);

private:
    int N;
    int n;
    capd::vectalg::Vector<T, 0> y0;
    capd::vectalg::Matrix<T, 0, 0> g;

    // Nowe pola wewnętrzne:
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> a_series;
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> c_series;
};

#include "ChebyshevOperatorFinite.tpp"
