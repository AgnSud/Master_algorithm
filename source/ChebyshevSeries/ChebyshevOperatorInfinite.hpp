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

    //TODO: W tej klasie wyliczane i przechowywane mają być przybliżenie pochodnej i odwrotnosc, ale narazie jest tylko przyblizenie pochodnej, po wyprostowaniu dopisac odwrotność
    ChebyshevOperatorInfinite(int N, int n, const ChebyshevOperatorFinite<T>& finiteOp);

    void set_pi0_hatA(const T& pi0_hatA_input);
    T get_pi0_hatA() const;

    void set_pi1_hatA(const VectorType& pi1_hatA_input);
    VectorType get_pi1_hatA() const;

    static T Pi0(const VectorType& x);
    static VectorType Pi1(const VectorType& x);
    static VectorType PiN(const VectorType& a, int N_, int n_);
    static VectorType PiN_x(const VectorType& x, int N_, int n_);
    template<class V>
    static V Pi1_j(const V& x, int j, int N_, int n_);

    T Pi0_HatA(const VectorType& x);

    // Operator \hat{A}(x)
    VectorType Pi1_HatA_k(const VectorType &x, int k);
    VectorType Pi1_HatA(const VectorType& x);
    std::pair<T, VectorType> HatA(const VectorType& x);

    //Ten operator wyciaga k-ty element z Pi1  (czyli wektor {a_{k*n}, ..., a_{(k+1)*n}}
    VectorType get_kth(const VectorType& a, int k);
private:
    int N;
    int n;
    ChebyshevOperatorFinite<T> finiteOp;

    T pi0_hatA;
    VectorType pi1_hatA;
};

#include "ChebyshevOperatorInfinite.tpp"