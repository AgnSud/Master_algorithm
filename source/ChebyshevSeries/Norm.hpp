#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"

template <typename T, int DIM = DIMENSION>
class Norm {
public:
    //TODO: dopisac norme operatorowa
    explicit Norm(T nu) : nu(nu) {}

    // Norma w przestrzeni l1_nu dla n=1
//    T computeNorm(const ChebyshevSeries<T, DIM>& a) const;
    template<class V>
    T computeNorm(const V& a) const;


    // Norma w przestrzeni l1_nu dla n>1
    template<class V>
    T computeNorm(const capd::vectalg::Vector<V, DIM>& vec, int n) const;

    //Norma operatorowa
    template<class V>
    T computeOperatorNorm(const V& a) const;


private:
    T nu;
};

#include "norm.tpp"
