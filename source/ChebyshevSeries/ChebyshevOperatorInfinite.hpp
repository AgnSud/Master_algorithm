#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"
#include "ChebyshevOperatorFinite.hpp"

template <typename T>
class ChebyshevOperatorInfinite {
public:
    typedef typename ChebyshevOperatorFinite<T>::VectorType VectorType;
    typedef typename ChebyshevOperatorFinite<T>::MatrixType MatrixType;
    typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;

    ChebyshevOperatorInfinite(int N, int n, const ChebyshevOperatorFinite<T>& finiteOp);

    T Pi0(const VectorType& x) const;
    VectorType Pi1(const VectorType& x) const;
    VectorType PiN(const VectorType& a) const;
    VectorType PiN_x(const VectorType& x) const;

    T Pi0_HatA(const VectorType& x);

    // Operator \hat{A}(x)
    VectorType Pi1_HatA_k(const VectorType &x, int k);

    //Ten operator wyciaga k-ty element z Pi1  (czyli wektor {a_{k*n}, ..., a_{(k+1)*n}}
    VectorType get_kth(const VectorType& a, int k);
private:
    int N;
    int n;
    ChebyshevOperatorFinite<T> finiteOp;
};

#include "ChebyshevOperatorInfinite.tpp"