#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"

//ewentualnie jakby sie zmienilo na parametr
#ifndef TRUNCATION
#define TRUNCATION 10
#endif  //TRUNCATION


/*
 * Matrix będzie reprezentować funkcje g w postaci wielowskaznikow jako macierz, np
 * Dla Van der Polla g = g(y) = (y2, -\mu (y1^2 - 1)y2 + y1)
 * zostanie zapisana jako macierz:
 * [[0, 0], [1,\mu], [0,0], [0,0], [0,1], [0,0], [0,0], [0,0], [0, -\mu], [0,0]]
 * gdzie wielowskazniki sa w porzadku leksykograficznym czyli powyzsza macierz odpowiada:
 * [(0,0), (0,1), (0,2), (0,3), (1,0), (1,1), (1,2), (2,0), (2,1), (3,0)]
 * */


template <typename T>
class ChebyshevOperatorFinite{
public:
    ChebyshevOperatorFinite() : N(1), n(1), y0(1), g(0, 2) {}

    ChebyshevOperatorFinite(
            int N, int n,
            const capd::vectalg::Vector<T, 0>& y0_init,
            const capd::vectalg::Matrix<T, 0, 0>& g_init);

    //domyslnie innego typu
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> computeC(
            const std::vector<std::vector<int>>& multiIndices,
            const capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0>& a_series);


private:
    int N; //N - przyblizenie skonczone
    int n; //wymiar y = wymiar a = wymiar c np dla Van der Poll n = 2
    capd::vectalg::Vector<T, 0> y0;
    capd::vectalg::Matrix<T, 0, 0> g;  //macierz wielowskaznikow, wymiaru n x liczba kombinacji ktore sum=n
};


#include "ChebyshevOperatorFinite.tpp"