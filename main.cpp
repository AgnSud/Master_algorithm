#include <iostream>
#include "capd/capdlib.h"
#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"

using namespace capd;

//funkcje pomocnicze
// Rekurencyjna funkcja generująca wszystkie multiindeksy a = (a_1, ..., a_order) o sumie <= maxDegree
void generateMultiIndicesRecursive(
        int order, int maxDegree, int pos, int remainingSum,
        std::vector<int>& current, std::vector<std::vector<int>>& result){
    if (pos == order) {
        if (remainingSum == 0) {
            result.push_back(current);
        }
        return;
    }

    for (int i = 0; i <= remainingSum; ++i) {
        current[pos] = i;
        generateMultiIndicesRecursive(order, maxDegree, pos + 1, remainingSum - i, current, result);
    }
}

std::vector<std::vector<int>> generateMultiIndices(int order, int maxDegree) {
    std::vector<std::vector<int>> result;
    std::vector<int> current(order, 0);
    for (int d = 0; d <= maxDegree; ++d)
        generateMultiIndicesRecursive(order, d, 0, d, current, result);
    return result;
}

int main() {
//    std::cout << "========================= TEST ChebyshevSeries =========================\n";
//
//    ChebyshevSeries<double> poly1(3);  // Wielomian o stopniu 3
//    ChebyshevSeries<double> poly2(4);  // Wielomian o stopniu 4
//
//    DVector coefficients1{1.0, 1, 1};
//    poly1 = coefficients1;
//
//    DVector coefficients2{2, 2, 2, 2};
//    poly2 = coefficients2;
//
//    std::cout << "\nWielomian 1:\n";
//    std::cout << poly1;
//    poly1.prettyPrint();
//
//    std::cout << "\nWielomian 2:\n";
//    std::cout << poly2;
//
//
//    ChebyshevSeries<double> sumPoly = poly1 + poly2;
//    std::cout << "\nWynik dodawania (poly1 + poly2):\n";
//    std::cout << sumPoly;
//
//    ChebyshevSeries<double> productPoly = poly1 * poly2;
//    std::cout << "\nWynik mnożenia (poly1 * poly2):\n";
//    std::cout << productPoly;
//
//    ChebyshevSeries<double> productPolyScalarLeft = 4.0 * poly2;
//    std::cout << "\nWynik mnożenia przez skalar 4.0 * poly2):\n";
//    std::cout << productPolyScalarLeft;
//
//    ChebyshevSeries<double> productPolyScalarRIght = poly2 * 2.0;
//    std::cout << "\nWynik mnożenia przez skalar poly2 * 2.0):\n";
//    std::cout << productPolyScalarRIght;
//
//    double x = 0.5;
//    std::cout << "\nTest wartości w punkcie x = " << x << ":\n";
//    std::cout << "poly1(x) = " << poly1(x) << "\n";
//    std::cout << "poly2(x) = " << poly2(x) << "\n";
//    std::cout << "sumPoly(x) = " << sumPoly(x) << "\n";
//    std::cout << "productPoly(x) = " << productPoly(x) << "\n";
//
//    if (typeid(poly1) == typeid(ChebyshevSeries<double>)) {
//        std::cout << "Obiekt poly1 jest typu ChebyshevSeries\n";
//    }
//
//    if (typeid(sumPoly) == typeid(ChebyshevSeries<double>)) {
//        std::cout << "Obiekt sumPoly jest typu ChebyshevSeries\n";
//    }
//
//    std::cout << "\n======================================================================\n\n";
//    std::cout << "========================= TEST Norm =========================\n";
//    vectalg::Vector<ChebyshevSeries<double>, 0> vec(2);
//    //pytanie -> czemu dopiero po stworzeniu konstruktora pustego (bez żadnych argumentów) w ChebyshevSeries
//    //powyzsza linijka nie zwróciła błędu, a wcześniej zwracała
//
//    vec[0] = poly1;
//    vec[1] = poly2;
//
//    std::cout << vec << "\n";
//
//    norm<double> myNorm(2.0);  // Definiujemy normę z wagą nu = 2.0
//
//    // Testujemy computeNorm dla pojedynczego szeregu
//    double norm1 = myNorm.computeNorm(poly1);
//    std::cout << "Norma pojedynczego szeregu: " << norm1 << std::endl;
//
//    // Testujemy computeNorm dla wektora szeregów
//    double norm2 = myNorm.computeNorm(vec, 2);
//    std::cout << "Norma wektora szeregów: " << norm2 << std::endl;
//    std::cout << "\n======================================================================\n\n";

    using T = double;
    int N = 5;
    int n = 2;

    capd::vectalg::Vector<T, 0> y0(n);
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, 0> a_series(n);

    for (int i = 0; i < n; ++i) {
        a_series[i] = ChebyshevSeries<T, DIMENSION>(N);  // N współczynników
        for (int k = 0; k < N; ++k) {
            a_series[i][k] = 2.0;
        }
    }
//    std::cout << a_series[0].power(2);


    std::vector<std::vector<int>> multiIndices = generateMultiIndices(n, 3);
    std::cout << "Multi-indices (up to degree 3):\n";
    std::cout << multiIndices << "\n";
    std::cout << multiIndices[6][1] << "\n\n\n";


    double mu = 0.5;  // lub dowolna inna wartość
    int A = multiIndices.size();
    capd::vectalg::Matrix<T, 0, 0> g(A, n);

    for (int j = 0; j < A; ++j) {
        const auto& alpha = multiIndices[j];
        if (alpha[0] == 0 && alpha[1] == 1) {
            g[j][0] = 1.0;         // y2
            g[j][1] = mu;          // +μ y2
        }
        else if (alpha[0] == 1 && alpha[1] == 0) {
            g[j][1] = 1.0;         // +y1
        }
        else if (alpha[0] == 2 && alpha[1] == 1) {
            g[j][1] = -mu;         // -μ y1^2 y2
        }
        else {
            g[j][0] = 0.0;
            g[j][1] = 0.0;
        }
    }

    std::cout << g << "\n\n";

    ChebyshevOperatorFinite<double> test1(N, n, y0,g);
    test1.computeC(multiIndices, a_series);


    return 0;
}
