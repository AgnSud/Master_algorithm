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

    // Konstruktor z inicjalizacją nu oraz n i N
    Norm(double nu, int N, int n) : nu(nu), N(N), n(n) {}

    // Norma w przestrzeni l1_nu dla n=1
    template<class V>
    T computeNorm(const V& a) const;

    // dualna norma jest tylko dla n=1
    template<class V>
    T computeDualNorm(const V& a) const;

    // Norma w przestrzeni l1_nu dla n>1
    template<class V>
    T computeNorm_n(const capd::vectalg::Vector<V, DIM>& vec) const;

    T computeOperatorNorm_Pi0(const MatrixType& C);
    T computeOperatorNorm_Pi1j(const MatrixType& C, int j);


    // Norma operatorowa
    T computeOperatorNorm(const MatrixType& C);

private:
    double nu;
    int n;
    int N;

    // eta_j
    T eta_j(const MatrixType& C, int j);

    // mu^j
    T mu_j(const MatrixType& C, int j_tilde);

    // ksi_tilde_j^{jk}
    T ksi_tilde_j_jk(const MatrixType& C, int j, int j_tilde, int k_tilde);

    // ksi_j^j
    T ksi_j_j(const MatrixType& C, int j, int j_tilde);
};

#include "norm.tpp"
