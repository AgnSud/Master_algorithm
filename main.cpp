#include <iostream>
#include "capd/capdlib.h"
#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"

using namespace capd;
using T = double;


// ---------- Pomocnicze ----------
void generateMultiIndicesRecursive(int order, int maxDegree, int pos, int remainingSum,
                                   std::vector<int>& current, std::vector<std::vector<int>>& result) {
    if (pos == order) {
        if (remainingSum == 0)
            result.push_back(current);
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

capd::vectalg::Matrix<T, DIMENSION, DIMENSION> defineFunctionG(std::vector<std::vector<int>> multiIndices, int n){
    int A = multiIndices.size();
    double rho = 28;
    double beta = 8/3.;
    double sigma = 10;
    capd::vectalg::Matrix<T, DIMENSION, DIMENSION> g(A, n);

    //definiuje poprzez wielowskazniki
    //Lorenz
    for (int j = 0; j < A; ++j) {
        const auto& alpha = multiIndices[j];
        if (alpha[0] == 0 && alpha[1] == 0 && alpha[2] == 1) {
            g[j][0] = 0;
            g[j][1] = 0;
            g[j][2] = -beta;
        }
        else if (alpha[0] == 0 && alpha[2] == 0 && alpha[1] == 1) {
            g[j][0] = sigma;
            g[j][1] = -1;
            g[j][2] = 0;
        }
        else if (alpha[1] == 0 && alpha[2] == 0 && alpha[0] == 1) {
            g[j][0] = -sigma;
            g[j][1] = rho;
            g[j][2] = 0;
        }
        else if (alpha[0] == 1 && alpha[2] == 1 && alpha[1] == 0) {
            g[j][0] = 0;
            g[j][1] = -1;
            g[j][2] = 0;
        }
        else if (alpha[0] == 1 && alpha[1] == 1 && alpha[2] == 0) {
            g[j][0] = 0;
            g[j][1] = 0;
            g[j][2] = 1;
        }
        else{
            g[j][0] = g[j][1] = g[j][2] = 0;
        }
    }


    //To ponizej do Van der Poll
//    for (int j = 0; j < A; ++j) {
//        const auto& alpha = multiIndices[j];
//        if (alpha[0] == 0 && alpha[1] == 1) {
//            g[j][0] = 1.0;
//            g[j][1] = mu;
//        }
//        else if (alpha[0] == 1 && alpha[1] == 0) {
//            g[j][1] = 1.0;
//        }
//        else if (alpha[0] == 2 && alpha[1] == 1) {
//            g[j][1] = -mu;
//        }
//    }

    // g(y) = (y1 + y2, y2)
//    for (int i = 0; i < A; ++i) {
//        const auto& alpha = multiIndices[i];
//        if (alpha[0] == 1 && alpha[1] == 0)
//            g[i][0] = 1.0;
//        if (alpha[0] == 0 && alpha[1] == 1) {
//            g[i][0] = 1.0;
//            g[i][1] = 1.0;
//        }
//    }

    return g;
}



// ---------- TEST 1: ChebyshevSeries ----------
void testChebyshevSeries() {
    std::cout << "\n========== TEST: ChebyshevSeries ==========\n";

    ChebyshevSeries<double> poly1(3);  // a_0 + a_1 T_1 + a_2 T_2
    poly1[0] = 1.0;
    poly1[1] = 1.0;
    poly1[2] = 1.0;

    ChebyshevSeries<double> poly2(4);
    poly2[0] = 2.0;
    poly2[1] = 2.0;
    poly2[2] = 2.0;
    poly2[3] = 2.0;

    std::cout << "poly1: "; poly1.prettyPrint(); std::cout << "\n";
    std::cout << "poly2: "; poly2.prettyPrint(); std::cout << "\n";

    std::cout << "poly1 + poly2: " << poly1 + poly2;
    std::cout << "poly1 * poly2: " << poly1 * poly2;
    std::cout << "4.0 * poly2: " << 4.0 * poly2;
    std::cout << "poly2 * 2.0: " << poly2 * 2.0;
//    std::cout << "poly1 ^ 2: " << poly1.power(2);

    double x = 0.5;
    std::cout << "\nEvaluations at x = " << x << ":\n";
    std::cout << "poly1(x) = " << poly1(x) << "\n";
    std::cout << "poly2(x) = " << poly2(x) << "\n";
    std::cout << "sum(x)   = " << (poly1 + poly2)(x) << "\n";
    std::cout << "prod(x)  = " << (poly1 * poly2)(x) << "\n";
}

// ---------- TEST 2: Norm ----------
void testNorms() {
    std::cout << "\n========== TEST: Norm ==========\n";

    ChebyshevSeries<double> poly1(3);
    poly1[0] = 1.0; poly1[1] = 1.0; poly1[2] = 1.0;

    ChebyshevSeries<double> poly2(4);
    poly2[0] = 2.0; poly2[1] = 2.0; poly2[2] = 2.0; poly2[3] = 2.0;

    vectalg::Vector<ChebyshevSeries<double>, 0> vec(2);
    vec[0] = poly1;
    vec[1] = poly2;

    norm<double> myNorm(2.0);

    std::cout << "Norm of poly1: " << myNorm.computeNorm(poly1) << "\n";
    std::cout << "Norm of vector: " << myNorm.computeNorm(vec, 2) << "\n";
}

// ---------- TEST 3: ChebyshevOperatorFinite (Van der Pol) ----------
// Przyjeta jest konwencja, ktora tez jest zgodna z tym co jest w pracy, ze N oznacza
// liczbe wspolczynnikow niezerowych, startujac od a_0, czyli od k=N mamy a_k=0
// ALE w tym celu, do wyliczenia operatora Czebyszewa F_N potrzebujemy wyznaczyc c_{k+1}
// co za tym idzie, mamy dwa rozne N -> ustalajac liczbe wspolczynnikow niezerowych przyblizen jako N
// czyli a_0, a_1, a_2, ..., a_N oraz c_0, c_1, c_2, ..., c_N
// bedziemy wyznaczac F_{N-1}
void testChebyshevOperatorFinite() {
    std::cout << "\n========== TEST: ChebyshevOperatorFinite Van der Pol ==========\n";

    constexpr int N = 2;
    constexpr int n = 3;
    constexpr int N_g = 2;

    //Ponizej startowe parametry
    constexpr T omega_start = 1; //omega jest czescia rownania
    capd::vectalg::Vector<T, 0> u0{5.0, 5.0, 5.0};
    capd::vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> a_series_start(n);

    //kazdy jest rozmiaru N
    a_series_start[0] = ChebyshevSeries<T, DIMENSION>{5,0};
    a_series_start[1] = ChebyshevSeries<T, DIMENSION>{5,0};
    a_series_start[2] = ChebyshevSeries<T, DIMENSION>{5,0};

    std::cout << "Ustawienia startowe (omega_0, a_series_0, u0):" << '\n';
    std::cout << "omega_0: " << omega_start << '\n';
    std::cout << "a_series_start: " << a_series_start << '\n';
    std::cout << "u0: " << u0 << '\n';
    //-------------------------------------------------------------------

    //Definicja v i u - takiego rozmiaru jak n
    ChebyshevSeries<T, DIMENSION> v{1.0, 0, 0};
    ChebyshevSeries<T, DIMENSION> w{6.0, 0, 0};
    //-------------------------------------------------------------------

    //Wielowskazniki
    auto multiIndices = generateMultiIndices(n, N_g);
    capd::vectalg::Matrix<T, DIMENSION, DIMENSION> g = defineFunctionG(multiIndices, n);
    std::cout << "Funkcja g: " << g << '\n';

    std::cout << "Wielowskaźniki:\n";
    for (const auto& mi : multiIndices) {
        std::cout << "(";
        for (int i = 0; i < mi.size(); ++i) {
            std::cout << mi[i] << (i < mi.size() - 1 ? "," : "");
        }
        std::cout << ") ";
    }
    std::cout << "\n-----------------------------------------------------------------\n\n";


    ChebyshevOperatorFinite<T> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 1;
    auto F = op.findFiniteSolution(omega_start, a_series_start, max_iterations);
}

// ---------- MAIN ----------
int main() {
//    testChebyshevSeries();
//    testNorms();
    testChebyshevOperatorFinite();
    return 0;
}
