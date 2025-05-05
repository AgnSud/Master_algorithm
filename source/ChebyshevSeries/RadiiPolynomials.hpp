#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "Norm.hpp"
#include "ChebyshevSeries.hpp"
#include "ChebyshevOperatorFinite.hpp"

template <typename T>
class RadiiPolynomials {
public:
    typedef typename ChebyshevOperatorFinite<T>::VectorType VectorType;
    typedef typename ChebyshevOperatorFinite<T>::MatrixType MatrixType;
    typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;

    RadiiPolynomials(int N_, int n_, double nu_, const ChebyshevOperatorFinite<T>& finiteOp_);

    T computeY0();
    //N_g przekazany jako argument
    T computeY1j(int j, int N_g);


private:
    // TODO: wgaa nie będzie przedziałem dla interval, prawda?
    int N;
    int n;
    double nu;
    ChebyshevOperatorFinite<T> finiteOp;

};



#include "RadiiPolynomials.tpp"