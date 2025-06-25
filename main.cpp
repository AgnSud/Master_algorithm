#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include "capd/capdlib.h"
#include "capd/map/Map.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"

#define LOGGER(A) cout << (#A) << " = " << (A) << endl

#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/Norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"
#include "source/ChebyshevSeries/RadiiPolynomials.hpp"


using namespace capd;
//using T = double;
using namespace std;

typedef vectalg::Vector<double, DIMENSION> DVectorType;
typedef vectalg::Matrix<double, DIMENSION, DIMENSION> DMatrixType;
typedef vectalg::Vector<ChebyshevSeries<double, DIMENSION>, DIMENSION> DVectorOfChebyshevsType;
typedef ChebyshevSeries<double, DIMENSION> DChebyshevsVectorType;
typedef vectalg::SumNorm<DVectorType, DMatrixType> DSumNormType;


typedef vectalg::Vector<Interval, DIMENSION> IVectorType;
typedef vectalg::Matrix<Interval, DIMENSION, DIMENSION> IMatrixType;
typedef vectalg::Vector<ChebyshevSeries<Interval, DIMENSION>, DIMENSION> IVectorOfChebyshevsType;
typedef ChebyshevSeries<Interval, DIMENSION> IChebyshevsVectorType;
typedef vectalg::SumNorm<IVectorType, IMatrixType> ISumNormType;



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

DMatrixType defineFunctionG(vector<vector<int>> multiIndices, int n){
    int A = multiIndices.size();
    double rho = 28;
    double beta = 8/3.;
    double sigma = 10;
    DMatrixType g(A, n);

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

IVectorType checkSolution(const IVectorOfChebyshevsType& a_series_approx,
                   double argument){
    IVectorType value(a_series_approx.dimension());
//    argument = argument * 2 - 1.;
    for(int i = 0; i < a_series_approx.dimension(); i++){
        value[i] = a_series_approx[i](argument * 2 - 1);
    }
    return value;
}

IVectorType randomVectorInBallWeightedNorm(int dim, Interval r, const Norm<Interval>& weighted_norm) {
    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> dist(0.0, 1.0);
    std::uniform_real_distribution<> u_dist(0.0, 1.0);  // dla losowego skalowania <1

    DVectorType vec(dim);
    for (int i = 0; i < dim; ++i)
        vec[i] = dist(gen);

    // Przekształcenie na przedziałowy wektor (z dokładnymi punktami)
    IVectorType ivec = capd::vectalg::convertObject<IVectorType, DVectorType>(vec);

    // Oblicz normę ważoną
    Interval norm_val = weighted_norm.computeNorm(ivec);

    // Skaluj tak, aby ||x||_ν ≤ r
    double u = std::pow(u_dist(gen), 1.0 / dim);
    Interval scale = Interval(r * u) / norm_val;

    // Skalowanie wektora przedziałowego
    for (int i = 0; i < dim; ++i)
        ivec[i] = ivec[i] * scale;

    return ivec;
}
void checkRadiiPolynomialsCoeffs(int n, double nu, RadiiPolynomials<Interval> radii_pol){
    string name = "radii_polynomials_coeffs_nu_" + std::to_string(nu) + ".csv";
    std::ofstream coeff_out(name);
    coeff_out << std::setprecision(17) << std::scientific;
    coeff_out << "name,nu,A_coeff,B_coeff,C_coeff\n";
    auto Y_bounds = radii_pol.getYBounds();

    {
        auto [A, B] = radii_pol.compute_Z0_terms();
        B -= 1;
        auto C = Y_bounds[0];
        coeff_out << "p_0," << nu << "," << A.mid().leftBound() << "," << B.mid().leftBound() << "," << C.mid().leftBound() << "\n";
    }

    for (int j = 0; j < n; ++j) {
        auto [A, B] = radii_pol.compute_Z1j_terms(j);
        B -= 1;
        auto C = Y_bounds[j + 1];
        coeff_out << "p_1" << j << "," << nu << "," << A.mid().leftBound() << "," << B.mid().leftBound() << "," << C.mid().leftBound() << "\n";
    }

    coeff_out.close();
}





// ---------- TEST 1: ChebyshevSeries ----------
void testChebyshevSeries() {
    cout << "\n========== TEST: ChebyshevSeries ==========\n";

    DChebyshevsVectorType poly1(2);  // a_0 + a_1 T_1 + a_2 T_2
    DChebyshevsVectorType poly2(2);
    DVectorType tmp = {5.481443923147439, 0.2407219615737195};
    poly1.setCoefficients(tmp);
    tmp = {12.339608554714449, 3.6698042773572244};
    poly2.setCoefficients(tmp);

    cout << "poly1: "; poly1.prettyPrint(); cout << "\n";
    cout << "poly2: "; poly2.prettyPrint(); cout << "\n";
    cout << "convolve: " << ChebyshevSeries<double>::convolve(poly1, poly2) << endl;

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

    DChebyshevsVectorType poly1(3);
    poly1[0] = 1.0; poly1[1] = 1.0; poly1[2] = 1.0;

    DChebyshevsVectorType poly2(3);
    poly2[0] = 2.0; poly2[1] = 2.0; poly2[2] = 2.0;

    vectalg::Vector<DChebyshevsVectorType, DIMENSION> vec(2);
    vec[0] = poly1;
    vec[1] = poly2;

    Norm<double> myNorm(2.0, 3, 2);

    cout << "Norm of poly1: " << myNorm.computeNorm(poly2) << "\n";
//    cout << "Norm of vector: " << myNorm.computeNorm(vec) << "\n";
}

void printPreparation(int N, int n, int N_g,
                      DVectorType & u0,
                      DVectorType & v,
                      DVectorType & w,
                      DMatrixType& g,
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

// TODO: r zadawany dla Z też ma być double? Bo to Z i Y to finalnie prawe konce przedziałow
void verify_bounds(int N, int n, double nu, Interval r, RadiiPolynomials<Interval> radii_pol,
                   ChebyshevOperatorFinite<Interval>& IFiniteOp) {
    IMatrixType A_N = IFiniteOp.getInverseDerivativeFinite();
    IVectorType x_approx = IFiniteOp.getX_approx();
    IVectorType F_x_approx = IFiniteOp.getF_x_approx();
    IVectorType AF = A_N * F_x_approx;
    int l = 1;
    IVectorType Pi1j_Z1j_norm_supremum(n);
    Interval Pi0_DT_supremum(0.);

//    IVectorType x1(n*N + 1), x2(n*N + 1);
//    x1[1] = r; x2[1] = r;
//    IVectorType x = x_approx + x1;
//    IVectorType DT_result = IFiniteOp.applyDT(x, x2);
//    std::cout << "size of DT(x + x1) * x2 = " << DT_result.dimension() << "\n DT(x + x1) * x2 = " << DT_result << std::endl;

    Norm<Interval> weighted_norm(nu, N, n);

    std::cout << "============ Verifying bounds ============\n";
    auto Y_bounds = radii_pol.getYBounds();
    auto Z_bounds = radii_pol.getZBounds();

    // TODO: Check: ||Pi_{1,j}(T(x) - x)|| ≤ Y1_j -> Y1_j jest prawym koncem przedzialu i porównujemy z prawym końcem przedziału po lewej, tak?
    for (int j = 0; j < n; ++j) {
        auto Pi1j_Y1j = radii_pol.Pi1_j(AF, j, N, n);
        auto Pi1j_Y1j_norm = weighted_norm.computeNorm(Pi1j_Y1j);
        std::cout << "Y_bound for j=" << j << ": " << Pi1j_Y1j_norm << " ≤ " << Y_bounds[j+1].rightBound()
                  << " => " << (Pi1j_Y1j_norm.rightBound() <= Y_bounds[j+1].rightBound() ? "OK" : "FAIL") << "\n";
    }

    // Check: |Pi_0(T(x) - x)| ≤ Y0
    Interval Pi0_val = capd::abs(radii_pol.Pi0(AF));
    std::cout << "Y0: " << Pi0_val << " ≤ " << Y_bounds[0].rightBound()
              << " => " << (Pi0_val.rightBound() <= Y_bounds[0].rightBound() ? "OK" : "FAIL") << "\n";


    for (int sample = 0; sample < l; ++sample) {
        IVectorType x1 = randomVectorInBallWeightedNorm(n * N + 1, r, weighted_norm);
        IVectorType x2 = randomVectorInBallWeightedNorm(n * N + 1, r, weighted_norm);
        IVectorType x = x_approx + x1;
        IVectorType DT_result = IFiniteOp.applyDT_decomposed(x, x2);

        for (int j = 0; j < n; ++j) {
            auto Pi1j_Z1j = radii_pol.Pi1_j(DT_result, j, N, n);
            auto Pi1j_Z1j_norm = weighted_norm.computeNorm(Pi1j_Z1j);
            if (Pi1j_Z1j_norm.rightBound() >= Pi1j_Z1j_norm_supremum[j].rightBound()){
                Pi1j_Z1j_norm_supremum[j] = Pi1j_Z1j_norm;
            }
        }

        Interval Pi0_DT = capd::abs(radii_pol.Pi0(DT_result));
        if (Pi0_DT.rightBound() >= Pi0_DT_supremum.rightBound()) {
            Pi0_DT_supremum = Pi0_DT;
        }
    }

    for (int j = 0; j < n; ++j) {
        std::cout << "sup_j ||Pi_{1," << j << "} D(T(x + x1)) x2|| = " << Pi1j_Z1j_norm_supremum[j]
                  << " ≤ " << Z_bounds[j+1] << " => " << (Pi1j_Z1j_norm_supremum[j].rightBound() <= Z_bounds[j+1].rightBound() ? "OK" : "FAIL") << "\n";
    }
    std::cout << "sup ||Pi_0 D(T(x + x1)) x2|| = " << Pi0_DT_supremum
              << " ≤ " << Z_bounds[0] << " => " << (Pi0_DT_supremum.rightBound() <= Z_bounds[0].rightBound() ? "OK" : "FAIL") << "\n";

}

// Przyjeta jest konwencja, ktora tez jest zgodna z tym co jest w pracy, ze N oznacza
// liczbe wspolczynnikow niezerowych, startujac od a_0, czyli od k=N mamy a_k=0
// ALE w tym celu, do wyliczenia operatora Czebyszewa F_N potrzebujemy wyznaczyc c_{k+1}
// co za tym idzie, mamy dwa rozne N -> ustalajac liczbe wspolczynnikow niezerowych przyblizen jako N
// czyli a_0, a_1, a_2, ..., a_N oraz c_0, c_1, c_2, ..., c_N
// bedziemy wyznaczac F_{N-1}
ChebyshevOperatorFinite<double> prepareChebyshevOperatorAndFindFiniteSolution(int N, int n,
                                                                              const IVectorType& _u0,
                                                                              DVectorType & v,
                                                                              DVectorType & w,
                                                                              DMatrixType& g,
                                                                              vector<vector<int>>& multiIndices) {
    cout << "Wyliczam przybliżone rozwiązanie..." << endl;
    auto u0 = capd::vectalg::convertObject<DVectorType>(_u0);
    LOGGER(u0);
    //Ponizej startowe parametry
    constexpr double omega_start = 1.; //omega jest czescia rownania
    DVectorOfChebyshevsType a_series_start(n);
    for (int i = 0; i < n; i++){
        a_series_start[i] = DChebyshevsVectorType(N);
        a_series_start[i][0] = u0[i];
    }
    a_series_start[0][1] = 1e-13; //zadane, aby macierz pochodnej byla odwracalna (pierwsza kolumna byla niezerowa)
    a_series_start[n-1][0] = w[n-1];
//    LOGGER(a_series_start);


    ChebyshevOperatorFinite<double> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 100;
    auto solution_approx = op.findFiniteSolution(omega_start, a_series_start, max_iterations);

//    cout << "Found approximate finite solution (omega, a) =" << endl;
//    cout << "omega = " << solution_approx.first << endl;
//    cout << "a_series_approx = " << solution_approx.second << endl;
//    cout << "with value F(omega, a) = " << op.getF_x_approx() << endl << endl;
    return op;
}


ChebyshevOperatorFinite<Interval> convertToInterval(int N, int n, const IVectorType& u0,
                                                    DVectorType & v,
                                                    DVectorType & w,
                                                    DMatrixType& g,
                                                    vector<vector<int>>& multiIndices,
                                                    const ChebyshevOperatorFinite<double>& finiteOp){
    /// to co zostało przekonwertowane na Interval to:
    /// * u0 (punkt startowy)
    /// * v i w (sekcja Poincare)
    /// * omega, a, x (dosłownie, używajac convertObject z a na [a,a]
    /// * g

    /// to co zostało wyliczone już w interval to:
    /// * c_series
    /// * F(x*)
    /// * derivative i inverse_derivative

    /// to co zostało zostawione bez zmiany na interval to:
    /// * N i n
    /// * multiIndeces
    // zamiana na arytmetyke przedziałowa
    cout << "Konwersja na arytmetykę przedziałową..." << endl;
    ChebyshevOperatorFinite<Interval> IFiniteOp(N, n, u0,
                                                capd::vectalg::convertObject<IMatrixType, DMatrixType>(g),
                                                capd::vectalg::convertObject<IVectorType, DVectorType>(v),
                                                capd::vectalg::convertObject<IVectorType, DVectorType>(w),
                                                multiIndices);

    IFiniteOp.setOmega( Interval(finiteOp.getOmega()));
    IFiniteOp.setASeries( capd::vectalg::convertObject<IVectorOfChebyshevsType, DVectorOfChebyshevsType>(finiteOp.getASeries()) );
    auto converted_X_approx = capd::vectalg::convertObject<IVectorType, DVectorType>(finiteOp.getX_approx());
//    LOGGER(IFiniteOp.getASeries());
//    LOGGER(converted_X_approx);
    for (int j = 0; j < n; j++){
        Interval tmp, delta;
        u0[j].split(tmp, delta);
        tmp = IFiniteOp.getASeries()[j][0] + delta;
        IFiniteOp.setACoeff(tmp, j, 0);
        converted_X_approx[1 + j] = tmp;
//        LOGGER(tmp);
    }
//    LOGGER(converted_X_approx);
//    LOGGER(IFiniteOp.getASeries());
    IFiniteOp.setX_approx(converted_X_approx);

    IVectorOfChebyshevsType IC_series = IFiniteOp.convertToSeriesFromXForm(IFiniteOp.compute_c(IFiniteOp.getX_approx()), 2*N-1);
    IFiniteOp.setCSeries(IC_series);
    IVectorType IF_x_approx = IFiniteOp(IFiniteOp.getX_approx());
    IFiniteOp.setF_x_approx(IF_x_approx);
    IFiniteOp.computeDerivativeInverse(IFiniteOp.getX_approx());

    return IFiniteOp;
}


double testRadiiPolynomials(int N, int n, int N_g, double nu, ChebyshevOperatorFinite<Interval>& IFiniteOp) {
    cout << "Wyliczam Y i Z bounds i znajduję promień r..." << endl;

    RadiiPolynomials<Interval> radii_pol(N, n, nu, IFiniteOp);
    Norm<Interval> weighted_norm(nu, N, n);

    //test na odwrocenie macierzy odwrotnej - odwrócona, powinna być zbliżona do IFiniteOp.derivative_finite
    IMatrixType interval_inverse_of_inverse_of_derivative_test = matrixAlgorithms::gaussInverseMatrix(IFiniteOp.getInverseDerivativeFinite());
//    cout << interval_inverse_of_inverse_of_derivative_test << endl;


//    auto gamma = radii_pol.compute_gamma();
//    LOGGER(gamma);
//
//    IMatrixType A_N = IFiniteOp.getInverseDerivativeFinite();
//    auto Pi_0_AN_op_norm = weighted_norm.computeOperatorNorm_Pi0(A_N);
//    LOGGER(Pi_0_AN_op_norm);
//
//    for (int j = 0; j < n; j++){
//        auto AN_op_norm = weighted_norm.computeOperatorNorm_Pi1j(A_N, j);
//        cout << "For j=" << j << " ";
//        LOGGER(AN_op_norm);
//    }
//
//    auto d1 = radii_pol.compute_d1();
//    LOGGER(d1);
//
//    auto h = radii_pol.compute_h();
//    LOGGER(h);
//
//    auto Z1 = radii_pol.compute_Z1();
//    LOGGER(Z1);
//    auto Pi_0_Z1 = radii_pol.Pi0(Z1);
//    LOGGER(Pi_0_Z1);
//
//    for (int j = 0; j < n; j++){
//        auto Pi_1_j_Z1 = radii_pol.Pi1_j(Z1, j, N, n);
//        auto norm_Pi_1_j_Z1 = weighted_norm.computeNorm(Pi_1_j_Z1);
//        cout << "For j=" << j << " ";
//        LOGGER(norm_Pi_1_j_Z1);
//    }
//
//    auto d2 = radii_pol.compute_d2();
//    LOGGER(d2);
//
//    cout << "================ Y - bounds ================" << endl;
//    radii_pol.compute_YBounds(N_g);
//    auto Y_bounds = radii_pol.getYBounds();
//    LOGGER(Y_bounds);
//
//    cout << "================ Z - bounds ================" << endl;


    double r = radii_pol.findRForRadiiPolynomials();
    checkRadiiPolynomialsCoeffs(n, nu, radii_pol);
//    Interval r_scalar = Interval(9.9197764446251516e-13);
//    auto polynomials = radii_pol(r_scalar);
//    LOGGER(polynomials);

//    cout << "========== KONIEC TESTU ==========\n";
    return r;
}

void compareWithTaylorAndSaveResults(int N, double z, int step, const string& points, const IVectorType& u0,
                                     const IVectorOfChebyshevsType& a_series_approx,
                                     double omega_approx_mid){
    omega_approx_mid = std::abs(omega_approx_mid);
    DMap vectorField("par:q;var:x,y,z;fun:10*(y-x),x*(28-z)-y,x*y-8*z/3;");
    DOdeSolver solver(vectorField, N);
    DCoordinateSection section(3, 2, z);
    DPoincareMap pm(solver, section);
    DTimeMap tm(solver);

    double rt = 0;
    auto u0_mid = capd::vectalg::convertObject<DVectorType>(u0);
    cout << "Result for Taylor: " << pm(u0_mid, rt) << endl;
    cout << "Time for Taylor: " << rt << endl;

    tm.stopAfterStep(true);
    int counter = 0;
    DTimeMap::SolutionCurve solution(0.);
    do {
        tm(rt,u0_mid,solution);
        counter++;
    }while(!tm.completed());

    double t_chebyshev = 0;
    double t_taylor = 0;
    double del = 1.0 / (omega_approx_mid * 100);
    string name = "trajectoryN_" + std::to_string(N) + "_step_" + std::to_string(step) + points + "_start_point1" + ".csv";

    std::ofstream fout(name);
    fout << std::setprecision(17) << std::scientific;
    fout << "t_chebyshev,t_taylor,x_cheb,y_cheb,z_cheb,x_tay,y_tay,z_tay,diff_norm\n";

    DSumNormType sumNorm;
    IVectorType cheb;
    while (true) {
        cheb = checkSolution(a_series_approx, t_chebyshev);
        auto cheb_d = capd::vectalg::convertObject<DVectorType>(cheb);
        auto tay = solution(t_taylor);
        double diff = sumNorm(cheb_d - tay);

        fout << t_chebyshev << "," << t_taylor << ","
             << cheb_d[0] << "," << cheb_d[1] << "," << cheb_d[2] << ","
             << tay[0] << "," << tay[1] << "," << tay[2] << "," << diff << "\n";

        if (t_chebyshev == 1.0 && t_taylor == rt) break;
        t_chebyshev = capd::min(t_chebyshev + 0.01, 1.0);
        t_taylor = capd::min(t_taylor + del, rt);
    }

    std::cout << "Zapisano plik: " << std::filesystem::absolute(name) << std::endl;
}

// ---------- MAIN ----------
int main() {
    cout.precision(17);

    int N = 28;
    double nu = 1.1;
    constexpr int n = 3;
    constexpr int N_g = 2;
    auto multiIndices = generateMultiIndices(n, N_g);
    DMatrixType g = defineFunctionG(multiIndices, n);
    int steps = 5;

    DVectorType v{0, 0, 1.};

    std::vector<DVectorType> ws = {
            {0, 0, 25.},
            {0, 0, 27.},
            {0, 0, 29.},
            {0, 0, 31.},
            {0, 0, 33.}
    };
    IVectorType u0{5.,5.,23.};
    auto u0_mid = capd::vectalg::convertObject<DVector>(u0);
//    printPreparation(N, n, N_g, u0, v, ws[0], g, multiIndices);

    for (int step = 0; step < steps; step++){
        cout << "======================== Iteracja " << step << "========================" << endl;
        auto w = ws[step];
        ChebyshevOperatorFinite<double> finiteOp = prepareChebyshevOperatorAndFindFiniteSolution(N, n, u0, v, w, g, multiIndices);
        LOGGER(u0);
        ChebyshevOperatorFinite<Interval> IFiniteOp = convertToInterval(N, n, u0, v, w, g, multiIndices, finiteOp);
        auto r = testRadiiPolynomials(N, n, N_g, nu, IFiniteOp);
        LOGGER(r);

        auto omega_approx = IFiniteOp.getOmega();
        auto a_series_approx = IFiniteOp.getASeries();
        Interval r_omega;
        double omega_approx_mid;
        omega_approx.split(omega_approx_mid, r_omega);

        auto cheb = checkSolution(a_series_approx, 1.0);
        cout << "Result at t = 1.0 for Chebyshev: " << cheb << endl;
        cout << "Scale time for Chebyshev: " << 1./omega_approx << endl;

        string points;
        for (int j = 0; j < steps; j++){
            points += "_" + std::to_string(int(ws[j][2]));
        }
        compareWithTaylorAndSaveResults(N, w[n-1], step, points, u0, a_series_approx, omega_approx_mid);

        u0 = checkSolution(a_series_approx, 1.0);
        LOGGER(u0);
        for (int j = 0; j < n-1; j++){
            u0[j] = Interval(u0[j].leftBound() - r, u0[j].rightBound() + r);
        }
        u0[n-1] = Interval(w[n-1]);
        u0_mid = capd::vectalg::convertObject<DVector>(u0);
        cout << "Punkt na sekcji Poincare o z=" << w[2] << ", punkt znaleziony to u0= " << u0 << endl;
        cout << "Środek u0= " << u0_mid << endl;
    }

    return 0;
}

// start_point2 = {}
// start_point1 = {5,5,23}