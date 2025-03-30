#pragma once

#include <capd/vectalg/lib.h>
#ifndef DIMENSION
#define DIMENSION 0
#endif  //DIMENSION

/*
 * ChebyshevSeries będzie klasą reprezentujaca szereg Czebyszewa na pojedynczej wspolrzednej
 * np a_k \in R^2, a_k = (b_k, c_k), rozwiniecie y = (y_1, y_2) = a_0 + \sum a_k T_k
 * Wtedy mamy dwa ChebyshevSeries:
 * ChebyshevSeries B: y_1 = b_0 + \sum b_k T_k
 * ChebyshevSeries C: y_2 = c_0 + \sum c_k T_k
 * Zatem samo y = DVector<ChebyshevSeries, DIM>
 * */

template <typename T, int DIM = DIMENSION>
class ChebyshevSeries : private capd::vectalg::Vector<T, DIM>{
public:
    //docelowo Vector<ChebyshevSeries, DIM>
    using CVector = capd::vectalg::Vector<T, DIM>;
    using capd::vectalg::Vector<T, DIM>::operator=;  // Używamy operatora przypisania z klasy bazowej
    using capd::vectalg::Vector<T, DIM>::operator[];

    ChebyshevSeries() : N(1), CVector(1) {}

    explicit ChebyshevSeries(int N) : N(N), CVector(N) {}

    // Zwraca wartosc wielomianu T_k(x)
    static T evaluateFirstKind(int k, T x);

    CVector getCoefficients() const;
    int getN() const;
    void prettyPrint() const;


    // wylicza wartosc w punkcie (do zmiany, narazie doslowny wzor zadany)
    T operator()(T x) const;

    ChebyshevSeries<T, DIM> operator+(const ChebyshevSeries<T, DIM>& other) const;
    ChebyshevSeries<T, DIM> operator-(const ChebyshevSeries<T, DIM>& other) const;
    ChebyshevSeries<T, DIM> operator*(const ChebyshevSeries<T, DIM>& other) const;
    ChebyshevSeries<T, DIM> power(int n) const;



    // Splot
    static ChebyshevSeries<T, DIM> convolve(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b);
    static T dot(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b);

    friend std::ostream& operator<<(std::ostream& os, const ChebyshevSeries<T, DIM>& a){
        os << "{";
        for (int i = 0; i < a.N-1; ++i) {
            os << a[i] << ", ";
        }
        os << a[a.N-1] << "}\n";
        return os;
    }

    friend ChebyshevSeries<T, DIM> operator*(T scalar, const ChebyshevSeries<T, DIM>& a){
        ChebyshevSeries<T, DIM> result(a.N);
        for (int k = 0; k < a.N; ++k) {
            result[k] = a[k] * scalar;
        }
        return result;
    }


    friend ChebyshevSeries<T, DIM> operator*(const ChebyshevSeries<T, DIM>& a, T scalar){
        return scalar * a;
    }

    friend ChebyshevSeries<T, DIM> operator-(T scalar, const ChebyshevSeries<T, DIM>& a) {
        ChebyshevSeries<T, DIM> result(a.getN());
        for (int i = 0; i < a.getN(); ++i) {
            result[i] = scalar - a[i];
        }
        return result;
    }

/*
* TODO:
 * efektywniejsze wyliczanie wartosci w punkcie
 * */

private:
    int N;  // Stopień wielomianu + 1 (liczba wspolczynnikow)
};

#include "ChebyshevSeries.tpp"