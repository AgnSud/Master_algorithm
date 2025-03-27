#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"

template <typename T, int DIM = DIMENSION>
class norm {
public:
    explicit norm(T nu) : nu(nu) {}

    // Norma w przestrzeni l1_nu dla n=1
    T computeNorm(const ChebyshevSeries<T, DIM>& a) const;

    // Norma w przestrzeni l1_nu dla n>1
    T computeNorm(const capd::vectalg::Vector<ChebyshevSeries<T, DIM>, DIM>& vec, int n) const;

    //Norma operatorowa
    T computeOperatorNorm(const ChebyshevSeries<T, DIM>& a) const;


private:
    T nu;
};

#include "norm.tpp"
