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
    typedef typename ChebyshevOperatorFinite<double>::VectorType DVectorType;
    typedef typename ChebyshevOperatorFinite<double>::MatrixType DMatrixType;
    typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;

    RadiiPolynomials(int N_, int n_, double nu_, const ChebyshevOperatorFinite<T>& finiteOp_);

    T Pi0(const VectorType& x);
    VectorType Pi1(const VectorType& x);
    VectorType PiN(const VectorType& a, int N_, int n_);
    VectorType PiN_x(const VectorType& x, int N_, int n_);
    template<class V>
    V Pi1_j(const V& x, int j, int N_, int n_);
    VectorType get_kth(const VectorType& a, int k);

    /// obliczy wszystkie h = [h0, h1x, h1y, hjz]
    VectorType compute_h();

    /// gamma bedzie double, jesli g bedzie double - jesli g bedzie interval to gamma bedzie interval
    /// ale g chyba powinno zostac double, tak?
    T compute_gamma();
    T compute_GammaMinus_a(const VectorType& a);
    T compute_GammaPlus_a(const VectorType& a);

    /// B_k jest interval, bo operatorNormPsi_ak jest interval oraz B_k jest Vectorem
    VectorType computeB_k(int k);
    /// to jest rowniez w takim formacie jak x
    VectorType compute_Z1_tilde();
    VectorType compute_Z1();

    T operatorNormPsi_ak(VectorType& a, int k);

    VectorType compute_d1();
    VectorType compute_d2();


    T computeY0();
    //N_g przekazany jako argument
    T computeY1j(int j, int N_g);

    T compute_Z1j(T r, int j);
    T compute_Z0(T r);

    void testOperatorNorm();
    double fast_pow(double number, int exp);

private:
    // TODO: waga nie będzie przedziałem dla interval, prawda?
    int N;
    int n;
    double nu;
    VectorType h;
    ChebyshevOperatorFinite<T> finiteOp;

    VectorType g_unit_vector(int j);
    VectorType g_ls(int l, int s);

    /// jesli g bedzie interval to to bedzie tez interval
    template <class V>
    double vector_sum(const V & v);
    template <class V>
    V vector_abs(const V & v);

};

#include "RadiiPolynomials.tpp"