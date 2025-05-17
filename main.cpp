#include <iostream>
#include <cmath>
#include "capd/capdlib.h"
#include "capd/map/Map.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/dynset/C1DoubletonSet.hpp"

#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/Norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"
#include "source/ChebyshevSeries/RadiiPolynomials.hpp"


using namespace capd;
using T = double;
using namespace std;

typedef vectalg::Vector<T, DIMENSION> VectorType;
typedef vectalg::Matrix<T, DIMENSION, DIMENSION> MatrixType;
typedef vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> VectorOfChebyshevsType;
typedef vectalg::Vector<Interval, DIMENSION> IVectorType;
typedef vectalg::Matrix<Interval, DIMENSION, DIMENSION> IMatrixType;
typedef vectalg::Vector<ChebyshevSeries<Interval, DIMENSION>, DIMENSION> IVectorOfChebyshevsType;
typedef vectalg::SumNorm<VectorType, MatrixType> NormType;


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

vectalg::Vector<T, DIMENSION> checkSolution(const vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION>& a_series_approx,
                   T argument){
    vectalg::Vector<T, DIMENSION> value(a_series_approx.dimension());
//    argument = argument * 2 - 1.;
    for(int i = 0; i < a_series_approx.dimension(); i++){
        value[i] = a_series_approx[i](argument * 2 - 1);
    }
    return value;
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

    ChebyshevSeries<double> poly2(3);
    poly2[0] = 2.0; poly2[1] = 2.0; poly2[2] = 2.0;

    vectalg::Vector<ChebyshevSeries<double>, DIMENSION> vec(2);
    vec[0] = poly1;
    vec[1] = poly2;

    Norm<double> myNorm(2.0, 3, 2);

    cout << "Norm of poly1: " << myNorm.computeNorm(poly1) << "\n";
//    cout << "Norm of vector: " << myNorm.computeNorm(vec) << "\n";
}

void printPreparation(int N, int n, int N_g,
                      vectalg::Vector<T, 0>& u0,
                      ChebyshevSeries<T>& v,
                      ChebyshevSeries<T>& w,
                      vectalg::Matrix<T, DIMENSION, DIMENSION>& g,
                      vector<vector<int>>& multiIndices){
    cout << "=============PRZYGOTOWANIE=============" << endl;
    cout << "Znalezienie rozwiązania dla układu Lorenza przy użyciu wielomianów Czebyszewa." << endl;
    cout << "Układ Lorenza jest układem " << n << " zmiennych, stopnia " << N_g << endl;
    cout << "Zadany problem brzegowy ma postać:" << endl;
    cout << "\t du/dt = 1/omega (10(y-x), 28x - xz - y, xy - 8/3z) \n\t u(0) = "<< u0 << "\n\t <" << v << ", " << w << " - u(1)> = 0" << endl;
    cout << "Reprezentowany jest za pomocą wielowskaźników -> g(u) = ";
    for (int i = 0; i < multiIndices.size(); i++){
        cout << g[i] << "*u^" << multiIndices[i] << (i < multiIndices.size() - 1 ? " + " : "");
    }
    cout << endl;
    cout << "Liczba wyliczanych współczynników szeregu Czebyszewa ustawiona została na N=" << N << endl << endl;
}




// Przyjeta jest konwencja, ktora tez jest zgodna z tym co jest w pracy, ze N oznacza
// liczbe wspolczynnikow niezerowych, startujac od a_0, czyli od k=N mamy a_k=0
// ALE w tym celu, do wyliczenia operatora Czebyszewa F_N potrzebujemy wyznaczyc c_{k+1}
// co za tym idzie, mamy dwa rozne N -> ustalajac liczbe wspolczynnikow niezerowych przyblizen jako N
// czyli a_0, a_1, a_2, ..., a_N oraz c_0, c_1, c_2, ..., c_N
// bedziemy wyznaczac F_{N-1}
ChebyshevOperatorFinite<T> prepareChebyshevOperatorAndFindFiniteSolution(int N, int n,
                                                                         VectorType u0,
                                                                         ChebyshevSeries<T>& v,
                                                                         ChebyshevSeries<T>& w,
                                                                         MatrixType& g,
                                                                         vector<vector<int>>& multiIndices) {
    cout << "============= WYLICZANIE PRZYBLIŻONEGO ROZWIĄZANIA =============" << endl;

    //Ponizej startowe parametry
    constexpr T omega_start = 1.; //omega jest czescia rownania
    vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> a_series_start(n);
    for (int i = 0; i < n; i++){
        a_series_start[i] = ChebyshevSeries<T, DIMENSION>(N);
        a_series_start[i][0] = u0[i];
    }
    a_series_start[0][1] = 1e-8; //zadanie, aby macierz pochodnej byla odwracalna (pierwsza kolumna byla niezerowa)


    ChebyshevOperatorFinite<T> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 100;
    auto solution_approx = op.findFiniteSolution(omega_start, a_series_start, max_iterations);

    cout << "Found approximate finite solution (omega, a) =" << endl;
    cout << "omega = " << solution_approx.first << endl;
    cout << "a_series_approx = " << solution_approx.second << endl;
    cout << "with value F(omega, a) = " << op.getF_x_approx() << endl << endl;
    return op;
}


// USUWAM CAŁA KLASĘ ChebyshevOperatorInfinite
void testChebyshevOperatorInfinite(int N, int n, ChebyshevOperatorFinite<T>& finiteOp) {
    cout << "\n========== TEST: ChebyshevOperatorInfinite ==========\n";
    cout << "========== KONIEC TESTU ==========\n";
}

ChebyshevOperatorFinite<Interval> convertToInterval(int N, int n, const VectorType& u0,
                                                    ChebyshevSeries<T>& v,
                                                    ChebyshevSeries<T>& w,
                                                    MatrixType& g,
                                                    vector<vector<int>>& multiIndices,
                                                    const ChebyshevOperatorFinite<double>& finiteOp){
    // zamiana na arytmetyke przedziałowa
    ChebyshevOperatorFinite<Interval> IFiniteOp(N, n,
                                                capd::vectalg::convertObject<IVectorType, VectorType>(u0),
                                                capd::vectalg::convertObject<IMatrixType, MatrixType>(g),
                                                capd::vectalg::convertObject<ChebyshevSeries<Interval>, ChebyshevSeries<T>>(v),
                                                capd::vectalg::convertObject<ChebyshevSeries<Interval>, ChebyshevSeries<T>>(w),
                                                multiIndices);

    IFiniteOp.setOmega( Interval(finiteOp.getOmega()));
    IFiniteOp.setASeries( capd::vectalg::convertObject<IVectorOfChebyshevsType, VectorOfChebyshevsType>(finiteOp.getASeries()) );
    IFiniteOp.setX_approx( capd::vectalg::convertObject<IVectorType, VectorType>(finiteOp.getX_approx()) );

    /// do powyzszego momentu uzywam convertObject, a reszta obliczyc już dla ŚCISŁEJ werjsji przybliżonego punktu?
    IVectorOfChebyshevsType IC_series = IFiniteOp.convertToSeriesFromXForm(IFiniteOp.compute_c(IFiniteOp.getX_approx()), 2*N-1);
    IFiniteOp.setCSeries(IC_series);

    IVectorType IF_x_approx = IFiniteOp(IFiniteOp.getX_approx());
    IFiniteOp.setF_x_approx(IF_x_approx);

    IFiniteOp.computeDerivativeInverse(IFiniteOp.getX_approx());
//    IFiniteOp.setCSeries( capd::vectalg::convertObject<IVectorOfChebyshevsType, VectorOfChebyshevsType>(finiteOp.getCSeries()) );
//    IFiniteOp.setF_x_approx( capd::vectalg::convertObject<IVectorType, VectorType>(finiteOp.getF_x_approx()) );
//    IFiniteOp.setDerivativeFinite( capd::vectalg::convertObject<IMatrixType, MatrixType>(finiteOp.getDerivativeFinite()) );
//    IFiniteOp.setInverseDerivativeFinite( capd::vectalg::convertObject<IMatrixType, MatrixType>(finiteOp.getInverseDerivativeFinite()) );

    return IFiniteOp;
}


void testRadiiPolynomials(int N, int n, int N_g, double nu, ChebyshevOperatorFinite<Interval>& IFiniteOp) {
    cout << "\n========== TEST: ChebyshevOperatorInfinite ==========\n";

    RadiiPolynomials<Interval> radii_pol(N, n, nu, IFiniteOp);

    //test na odwrocenie macierzy odwrotnej - odwrócona, powinna być zbliżona do IFiniteOp.derivative_finite
    IMatrixType interval_inverse_of_inverse_of_derivative_test = matrixAlgorithms::gaussInverseMatrix(IFiniteOp.getInverseDerivativeFinite());
//    cout << interval_inverse_of_inverse_of_derivative_test << endl;

    auto diff_T_x =  IFiniteOp.NewtonLikeOperatorTx_x(IFiniteOp.getX_approx());
    cout << "T(x*) - x* = " << diff_T_x << endl;
    auto Pi_0_diff_T_x = radii_pol.Pi0(diff_T_x);
    auto abs_Pi_0_diff_T_x = capd::abs(Pi_0_diff_T_x);

    auto Y_0 = radii_pol.computeY0();
    cout << "|Pi_0 (T(x*) - x*)| = " << abs_Pi_0_diff_T_x << ", Y_0 = " << Y_0 << endl;
    if (abs_Pi_0_diff_T_x == Y_0) //czemu dziala ==, ale nie dziala <= ???
        cout << "Y_0 bound OK" << endl;
    else
        cout << "Y_0 bound NOT OK" << endl;

    /// zawsze 0 będzie na dole przedziału i trochę nie wiem jak to poprawić (bo  wspolczynniki sa przedzialami [-,+])
    for(int i = 0; i < n; i++){
        Norm<Interval> weighted_norm(nu, N, n);
        auto Y_1_i = radii_pol.computeY1j(i, N_g);
        auto Pi_1i_diff_T_x = radii_pol.Pi1_j(diff_T_x, i, N, n);
        auto norm_Pi_1i_diff_T_x = weighted_norm.computeNorm(Pi_1i_diff_T_x);
        cout << "||Pi_1" << i << "(T(x*) - x*)|| = " << norm_Pi_1i_diff_T_x << ", Y_1_" << i << "= " << Y_1_i << endl;
    }

    auto h = radii_pol.compute_h();
    cout << "h=" << h << endl;
    cout << "========== KONIEC TESTU ==========\n";
}

// ---------- MAIN ----------
int main() {
    cout.precision(23);

    constexpr int N = 25;
    constexpr int n = 3;
    constexpr int N_g = 2;
    double nu = 1.01;
    vectalg::Vector<T, 0> u0{5., 5., 23.};
    ChebyshevSeries<T, DIMENSION> v{0, 0, 1.};
    ChebyshevSeries<T, DIMENSION> w{0, 0, 27.};
    auto multiIndices = generateMultiIndices(n, N_g);
    vectalg::Matrix<T, DIMENSION, DIMENSION> g = defineFunctionG(multiIndices, n);

    printPreparation(N, n, N_g, u0, v, w, g, multiIndices);

//    testChebyshevSeries();
//    testNorms();
    ChebyshevOperatorFinite<T> finiteOp = prepareChebyshevOperatorAndFindFiniteSolution(N, n, u0, v, w, g, multiIndices);
    ChebyshevOperatorFinite<Interval> IFiniteOp = convertToInterval(N, n, u0, v, w, g, multiIndices, finiteOp);
    cout << IFiniteOp << endl;
//    Interval omegaI = capd::vectalg::convertObject<Interval, double>(finiteOp.getOmega()); - czemu to nie działa?
    testRadiiPolynomials(N, n, N_g, nu, IFiniteOp);
    cout << "##########################################################################################\n";

    {
        DMap vectorField("par:q;var:x,y,z;fun:10*(y-x),x*(28-z)-y,x*y-8*z/3;");
        DOdeSolver solver(vectorField,N);
        DCoordinateSection section(3,2,27.);
        DPoincareMap pm(solver,section);
        DTimeMap tm(solver);

        DVector u{5,5,23};
        DVector u1 = u;
        double rt=0;
        cout << pm(u,rt) << endl;
        cout << rt << endl;

        // ze wynik powyzszego, powinien byc taki jak ponizszego po czasie 1?
        // EDIT: NIE, czas zadany w ponizszym to czas przykladowy, aby byly te same wyniki potrzeba
        // dac wyjscie czasu powyzszego
        tm.stopAfterStep(true);
        int counter = 0;

        DTimeMap::SolutionCurve solution(0.);
        do {
            tm(rt,u1,solution);
            counter++;
        }while(!tm.completed());

        cout << "counter= " << counter << endl;
//        T t = 0;
//        rt = 3.50971460024166;
//        T del = rt/10.;
//        while (t < rt){
//            cout << "sol(" << t << ") = " << solution(t) << endl;
//            t += del;
//        }

        auto omega_approx = finiteOp.getOmega();
        auto a_series_approx = finiteOp.getASeries();
        T t_chebyshev = 0;
        T t_taylor = 0;
        T del = 1 / (omega_approx * 10);
        cout << "rt= " << rt << endl;
        while (t_chebyshev < 1){
            auto diff = checkSolution(a_series_approx, t_chebyshev) - solution(t_taylor);
            cout << "u("<< t_chebyshev << ") - sol(" << t_taylor << ") = " << diff << endl;
            t_chebyshev += 0.1;
            t_taylor += del;
            t_taylor = capd::min(t_taylor, rt);
        }
    }
    return 0;
}
