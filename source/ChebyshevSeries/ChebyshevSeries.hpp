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
class ChebyshevSeries : public capd::vectalg::Vector<T, DIM>{
public:
    //docelowo Vector<ChebyshevSeries, DIM>
    using VectorType = capd::vectalg::Vector<T, DIM>;
    using capd::vectalg::Vector<T, DIM>::operator=;  // Używamy operatora przypisania z klasy bazowej
    using capd::vectalg::Vector<T, DIM>::operator[];

    ChebyshevSeries() : N(1), VectorType(1) {}

    ChebyshevSeries(int N) : N(N), VectorType(N) {}
    ChebyshevSeries(std::initializer_list<T> list);

    // Zwraca wartosc wielomianu T_k(x)
    static T evaluateFirstKind(int k, T x);

    VectorType getCoefficients() const;
    int getN() const;
    void prettyPrint() const;


    // wylicza wartosc w punkcie (do zmiany, narazie doslowny wzor zadany)
    T operator()(T x) const;

    ChebyshevSeries<T, DIM> operator+(const ChebyshevSeries<T, DIM>& other) const;
    ChebyshevSeries<T, DIM> operator-(const ChebyshevSeries<T, DIM>& other) const;
    template<class V>
    V operator*(const V& other) const;

//    template<class V>
//    V power(int n) const;


    // Splot
    template<class V>
    static V convolve(const V& a, const V& b);

    //funkcja dot - ALE nie uzywam,  bo przeciez iloczyn skalarny jest na WEKTORACH z R^n nie na ChebyshevSeries
    //static T dot(const ChebyshevSeries<T, DIM>& a, const ChebyshevSeries<T, DIM>& b);

    friend std::ostream& operator<<(std::ostream& os, const ChebyshevSeries<T, DIM>& a){
        os << "{";
        for (int i = 0; i < a.N-1; ++i) {
            os << a[i] << ", ";
        }
        os << a[a.N-1] << "}";
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