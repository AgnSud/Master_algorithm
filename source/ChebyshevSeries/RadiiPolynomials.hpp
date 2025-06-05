#pragma once

#include <capd/vectalg/lib.h>
#include <cmath>
#include "Norm.hpp"
#include "ChebyshevSeries.hpp"
#include "ChebyshevOperatorFinite.hpp"

template <typename T>
class RadiiPolynomials {
public:
    typedef  ChebyshevOperatorFinite<T>::VectorType VectorType;
    typedef  ChebyshevOperatorFinite<T>::MatrixType MatrixType;
    typedef  vectalg::SumNorm<VectorType, MatrixType> NormType;
    typedef  ChebyshevOperatorFinite<double>::VectorType DVectorType;
    typedef  ChebyshevOperatorFinite<double>::MatrixType DMatrixType;
    typedef  vectalg::SumNorm<DVectorType, DMatrixType> DNormType;

    RadiiPolynomials(int N_, int n_, double nu_, const ChebyshevOperatorFinite<T>& finiteOp_);

    DVectorType getYBounds();
    DVectorType getZBounds();
//    VectorType getYBounds();
//    VectorType getZBounds();

    T Pi0(const VectorType& x);
    VectorType Pi1(const VectorType& x);
    VectorType PiN(const VectorType& a, int N_, int n_);
    VectorType PiN_x(const VectorType& x, int N_, int n_);
    template<class V>
    V Pi1_j(const V& x, int j, int N_, int n_);

    void compute_YBounds(int N_g);
    T computeY0();
    //N_g przekazany jako argument
    T computeY1j(int j, int N_g);

//    Czy V moze byc inne dla kazdej funkcji?
//    VectorType compute_ZBounds(VectorType& r);
//    template <class V>
//    V compute_Z0(V r);
//    template <class V>
//    V compute_Z1j(V r, int j);

    void compute_ZBounds(DVectorType& r);
    T compute_Z0(double r);
    T compute_Z1j(double r, int j);


    /// obliczy wszystkie h = [h0, h1x, h1y, hjz]
    VectorType compute_h();

    /// B_k jest interval, bo operatorNormPsi_ak jest interval oraz B_k jest Vectorem
    VectorType compute_Z1();
    VectorType compute_Z1_tilde();
    VectorType computeB_k(int k);
    T operatorNormPsi_ak(VectorType& a, int k);

    /// gamma bedzie double, jesli g bedzie double - jesli g bedzie interval to gamma bedzie interval
    /// ale g chyba powinno zostac double, tak?
    T compute_gamma();

    VectorType compute_d1();
    VectorType compute_d2();
    T compute_GammaMinus_a(const VectorType& a);
    T compute_GammaPlus_a(const VectorType& a);


    T findRIntervalForRadiiPolynomials();
    T findRootNewtonForRadiiPolynomials(T r_start);
    template <class V>
    V operator()(V& r);

    void testOperatorNorm();



private:
    // TODO: waga nie będzie przedziałem dla interval, prawda?
    int N;
    int n;
    double nu;
    ChebyshevOperatorFinite<T> finiteOp;
    DVectorType Y_bounds;  // length n + 1, [Y_0, Y_11, Y_12, Y_13] - wektor juz samych bound <double>
    DVectorType Z_bounds;  // length n + 1, [Z_0, Z_11, Z_12, Z_13] - wektor już samych bound <double>

    VectorType g_unit_vector(int j);
    VectorType g_ls(int l, int s);

    template <class V>
    V vector_abs(const V & v);

};

#include "RadiiPolynomials.tpp"