#include <iostream>
#include "capd/capdlib.h"
#include "capd/map/Map.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/dynset/C1DoubletonSet.hpp"

#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorInfinite.hpp"

using namespace capd;
using T = double;
using namespace std;


// ---------- Pomocnicze ----------
void generateMultiIndicesRecursive(int order, int maxDegree, int pos, int remainingSum,
                                   vector<int>& current, vector<vector<int>>& result) {
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

vector<vector<int>> generateMultiIndices(int order, int maxDegree) {
    vector<vector<int>> result;
    vector<int> current(order, 0);
    for (int d = 0; d <= maxDegree; ++d)
        generateMultiIndicesRecursive(order, d, 0, d, current, result);
    return result;
}

vectalg::Matrix<T, DIMENSION, DIMENSION> defineFunctionG(vector<vector<int>> multiIndices, int n){
    int A = multiIndices.size();
    double rho = 28;
    double beta = 8/3.;
    double sigma = 10;
    vectalg::Matrix<T, DIMENSION, DIMENSION> g(A, n);

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
//    double mu = 0.5;
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

void checkSolution(const vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION>& a_series_approx,
                   T argument){
    vectalg::Vector<T, DIMENSION> value(a_series_approx.dimension());
//    argument = argument * 2 - 1.;
    for(int i = 0; i < a_series_approx.dimension(); i++){
        value[i] = a_series_approx[i](argument * 2 - 1);
    }
    cout << "u(" << argument << ")= " << value << '\n';
}



// ---------- TEST 1: ChebyshevSeries ----------
void testChebyshevSeries() {
    cout << "\n========== TEST: ChebyshevSeries ==========\n";

    ChebyshevSeries<double> poly1(2);  // a_0 + a_1 T_1 + a_2 T_2
    ChebyshevSeries<double> poly2(2);
    vectalg::Vector<double, 0> tmp = {5.481443923147439, 0.2407219615737195};
    poly1.setCoefficients(tmp);
    tmp = {12.339608554714449, 3.6698042773572244};
    poly2.setCoefficients(tmp);
//    poly2[0] = 2.0;
//    poly2[1] = 2.0;
//    poly2[2] = 2.0;
//    poly2[3] = 2.0;

    cout << "poly1: "; poly1.prettyPrint(); cout << "\n";
    cout << "poly2: "; poly2.prettyPrint(); cout << "\n";
    cout << "convolve: " << ChebyshevSeries<T>::convolve(poly1, poly2) << endl;

    cout << "poly1 + poly2: " << poly1 + poly2;
    cout << "poly1 * poly2: " << poly1 * poly2;
    cout << "4.0 * poly2: " << 4.0 * poly2;
    cout << "poly2 * 2.0: " << poly2 * 2.0;

    double x = 0.5;
    cout << "\nEvaluations at x = " << x << ":\n";
    cout << "poly1(x) = " << poly1(x) << "\n";
    cout << "poly2(x) = " << poly2(x) << "\n";
    cout << "sum(x)   = " << (poly1 + poly2)(x) << "\n";
    cout << "prod(x)  = " << (poly1 * poly2)(x) << "\n";
}

// ---------- TEST 2: Norm ----------
void testNorms() {
    cout << "\n========== TEST: Norm ==========\n";

    ChebyshevSeries<double> poly1(3);
    poly1[0] = 1.0; poly1[1] = 1.0; poly1[2] = 1.0;

    ChebyshevSeries<double> poly2(4);
    poly2[0] = 2.0; poly2[1] = 2.0; poly2[2] = 2.0; poly2[3] = 2.0;

    vectalg::Vector<ChebyshevSeries<double>, 0> vec(2);
    vec[0] = poly1;
    vec[1] = poly2;

    norm<double> myNorm(2.0);

    cout << "Norm of poly1: " << myNorm.computeNorm(poly1) << "\n";
    cout << "Norm of vector: " << myNorm.computeNorm(vec, 2) << "\n";
}

// ---------- TEST 3: ChebyshevOperatorFinite (Van der Pol) ----------
// Przyjeta jest konwencja, ktora tez jest zgodna z tym co jest w pracy, ze N oznacza
// liczbe wspolczynnikow niezerowych, startujac od a_0, czyli od k=N mamy a_k=0
// ALE w tym celu, do wyliczenia operatora Czebyszewa F_N potrzebujemy wyznaczyc c_{k+1}
// co za tym idzie, mamy dwa rozne N -> ustalajac liczbe wspolczynnikow niezerowych przyblizen jako N
// czyli a_0, a_1, a_2, ..., a_N oraz c_0, c_1, c_2, ..., c_N
// bedziemy wyznaczac F_{N-1}
ChebyshevOperatorFinite<T> testChebyshevOperatorFinite(int N, int n) {
    cout << "\n========== TEST: ChebyshevOperatorFinite Van der Pol ==========\n";

    constexpr int N_g = 2;

    //Ponizej startowe parametry
    constexpr T omega_start = 0.27; //omega jest czescia rownania
    vectalg::Vector<T, 0> u0{5., 5., 23.};
    vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> a_series_start(n);

    //kazdy jest rozmiaru N
    a_series_start[0] = ChebyshevSeries<T, DIMENSION>(N);
    a_series_start[1] = ChebyshevSeries<T, DIMENSION>(N);
    a_series_start[2] = ChebyshevSeries<T, DIMENSION>(N);
    a_series_start[0][0] = 5.0;
    a_series_start[0][1] = 1e-8;
    a_series_start[1][0] = 5.0;
    a_series_start[2][0] = 23.0;


    //Definicja v i u - takiego rozmiaru jak n
    ChebyshevSeries<T, DIMENSION> v{0, 0, 1.};
    ChebyshevSeries<T, DIMENSION> w{0, 0, 27.};
    //-------------------------------------------------------------------

    cout << "Ustawienia startowe:" << '\n';
    cout << "omega_0: " << omega_start << '\n';
    cout << "a_series_start: " << a_series_start << '\n';
    cout << "u_0: " << u0 << '\n';
    cout << "v: " << v << '\n';
    cout << "w: " << w << '\n';
    cout << "n:" << n << '\n';
    cout << "N: " << N << '\n';
    cout << "N_0: " << N_g << '\n';
    //-------------------------------------------------------------------


    //Wielowskazniki
    auto multiIndices = generateMultiIndices(n, N_g);
    vectalg::Matrix<T, DIMENSION, DIMENSION> g = defineFunctionG(multiIndices, n);
    cout << "Funkcja g: " << g << '\n';

    cout << "Wielowskaźniki:\n";
    for (const auto& mi : multiIndices) {
        cout << "(";
        for (int i = 0; i < mi.size(); ++i) {
            cout << mi[i] << (i < mi.size() - 1 ? "," : "");
        }
        cout << ") ";
    }
    cout << "\n-----------------------------------------------------------------\n\n";


    ChebyshevOperatorFinite<T> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 100;
    auto solution_approx = op.findFiniteSolution(omega_start, a_series_start, max_iterations);
    auto omega_approx = solution_approx.first;
    auto a_series_approx = solution_approx.second;
    cout << "-------------------------------------------------------" << '\n';
    cout << "omega_approx: " << omega_approx << '\n';
    cout << "1/omega: " << 1/omega_approx << '\n';
    cout << "a_series_approx: " << a_series_approx << '\n';
    cout << "c_series_approx: " << op.getCSeries() << '\n';
    cout << "jacobianDerivative: " << op.getDerivativeFinite() << '\n';
    cout << "-------------------------------------------------------" << '\n';
    cout << "jacobianInverseDerivative: " << op.getInverseDerivativeFinite() << '\n';
    cout << "-------------------------------------------------------" << '\n';

//    T t = 0;
//    while (t < 1){
//        checkSolution(a_series_approx, t);
//        t += 0.1;
//    }
    return op;
}

void testChebyshevOperatorInfinite(int N, int n, ChebyshevOperatorFinite<T>& finiteOp) {
    cout << "\n========== TEST: ChebyshevOperatorInfinite (Sztuczne dane) ==========\n";

    int totalLength = 2 * N * n + 1;

    using VectorType = vectalg::Vector<T, DIMENSION>;
    using MatrixType = vectalg::Matrix<T, DIMENSION, DIMENSION>;

    ChebyshevOperatorInfinite<T> op(N, n, finiteOp);

    // Przykładowy wektor x = [omega, a0, a1, ..., a_{2*N*n+1}]
//    vectalg::Vector<T, DIMENSION> x(totalLength);
//    for (int i = 0; i < totalLength; ++i)
//        x[i] = static_cast<T>(i + 1);  // x = [1, 2, 3, ..., totalLength]

    auto x = finiteOp.convertToXVector();
    cout << x.dimension() << endl;
    VectorType x_expanded(100 * N * n + 1);
    for (int i = 0; i < x.dimension(); i++){
        x_expanded[i] = x[i];
    }
    cout << "x = " << x_expanded << "\n";

    // Test Pi0
    T pi0 = op.Pi0(x_expanded);
    cout << "Pi0(x_expanded) = " << pi0 << "\n";

    // Test Pi1
    auto pi1 = op.Pi1(x_expanded);
    cout << "Pi1(x_expanded) = " << pi1 << "\n";

    // Test PiN (obetnij po N = 3 => 6 elementów z pi1, reszta 0)
    auto piN = op.PiN(pi1);
    cout << "PiN(a) = " << piN << "\n";

    // Test PiN_x
    auto piN_x = op.PiN_x(x_expanded);
    cout << "PiN_x(x_expanded) = " << piN_x << "\n";

    // Pi0_HatA test
    T Pi0_hatA = op.Pi0_HatA(x_expanded);
    cout << "Pi0_HatA(x_expanded) = " << Pi0_hatA << "\n";

    cout << "Pi1_HatA(x, 0) = " << op.Pi1_HatA_k(x_expanded, 0) << "\n";
    cout << "Pi1_HatA(x, 1) = " << op.Pi1_HatA_k(x_expanded, 1) << "\n";
    cout << "Pi1_HatA(x, 2) = " << op.Pi1_HatA_k(x_expanded, 2) << "\n";


    cout << "========== KONIEC TESTU ==========\n";
}


// ---------- MAIN ----------
int main() {
    cout.precision(10);

    constexpr int N = 2;
    constexpr int n = 3;

//    testChebyshevSeries();
//    testNorms();
    ChebyshevOperatorFinite<T> finiteOp = testChebyshevOperatorFinite(N, n);
    testChebyshevOperatorInfinite(N, n, finiteOp);
    cout << "##########################################################################################\n";

//    int order = 20;
//    {
//        DMap vectorField("par:q;var:x,y,z;fun:10*(y-x),x*(28-z)-y,x*y-8*z/3;");
//        DOdeSolver solver(vectorField,order);
//        DCoordinateSection section(3,2,27.);
//        DPoincareMap pm(solver,section);
//        DTimeMap tm(solver);
//
//        DVector u{5,5,23};
//        DVector u1 = u;
//        double rt=0;
//        cout << pm(u,rt) << endl;
//        cout << rt << endl;
//
//        // TODO: Czy tutaj solution daje rozwiazanie po czasie zadanym? Chyba inaczej zrozumialam poprzednio,
//        // ze wynik powyzszego, powinien byc taki jak ponizszego po czasie 1?
//        // EDIT: NIE, czas zadany w ponizszym to czas przykladowy, aby byly te same wyniki potrzeba
//        // dac wyjscie czasu powyzszego
//        tm.stopAfterStep(true);
//        int counter = 0;
//
//        DTimeMap::SolutionCurve solution(0.);
//        do {
//            tm(rt,u1,solution);
//            counter++;
//        }while(!tm.completed());
//
//        cout << "counter= " << counter << endl;
//        T t = 0;
//        rt = 0.27713911300544136;
//        T del = rt/10.;
//        while (t < rt){
//            cout << "sol(" << t << ") = " << solution(t) << endl;
//            t += del;
//        }
//
//    }
    return 0;
}
