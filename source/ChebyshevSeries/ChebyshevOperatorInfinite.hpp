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

    ChebyshevOperatorInfinite(int N, const ChebyshevOperatorFinite<T>& finiteOp);

    T Pi0(const VectorType& x) const;
    VectorType Pi1(const VectorType& x) const;
    VectorType PiN(const VectorType& a) const;
    VectorType PiN_x(const VectorType& x) const;

    T Pi0_HatA(const VectorType& x) const;

    // Operator \hat{A}(x)
    VectorType applyHatA(const VectorType& x) const;
private:
    int N;
    ChebyshevOperatorFinite<T> finiteOp;
};

#include "ChebyshevOperatorInfinite.tpp"