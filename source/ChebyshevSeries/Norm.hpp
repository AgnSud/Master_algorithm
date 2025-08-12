#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "ChebyshevSeries.hpp"
#include "ChebyshevOperatorFinite.hpp"

using namespace std;
using namespace capd;

template <typename T, int DIM = DIMENSION>
class Norm {
public:
    typedef typename ChebyshevOperatorFinite<T>::MatrixType MatrixType;
    typedef typename ChebyshevOperatorFinite<T>::VectorType VectorType;

    Norm(double nu, int N, int n) : nu(nu), N(N), n(n) {}

    template<class V>
    T computeNorm(const V& a) const;

    template<class V>
    T computeDualNorm(const V& a) const;

    template<class V>
    T computeNorm_n(const capd::vectalg::Vector<V, DIM>& vec) const;

    T computeOperatorNorm_Pi0(const MatrixType& C);
    T computeOperatorNorm_Pi1j(const MatrixType& C, int j);

    T computeOperatorNorm(const MatrixType& C);

private:
    double nu;
    int n;
    int N;

    T eta_j(const MatrixType& C, int j);
    T mu_j(const MatrixType& C, int j_tilde);
    T ksi_tilde_j_jk(const MatrixType& C, int j, int j_tilde, int k_tilde);
    T ksi_j_j(const MatrixType& C, int j, int j_tilde);
};

#include "norm.tpp"
