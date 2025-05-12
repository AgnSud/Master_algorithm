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

    static T Pi0(const VectorType& x);
    static VectorType Pi1(const VectorType& x);
    static VectorType PiN(const VectorType& a, int N_, int n_);
    static VectorType PiN_x(const VectorType& x, int N_, int n_);
    template<class V>
    static V Pi1_j(const V& x, int j, int N_, int n_);
    VectorType get_kth(const VectorType& a, int k);

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