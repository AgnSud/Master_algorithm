#include <iostream>
#include <filesystem>
#include <random>
#include "capd/capdlib.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"

#define LOGGER(A) cout << (#A) << " = " << (A) << endl

#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/Norm.hpp"
#include "source/ChebyshevSeries/ChebyshevOperatorFinite.hpp"
#include "source/ChebyshevSeries/RadiiPolynomials.hpp"

#include "helpers.hpp"
#include "testing.hpp"


using namespace capd;
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
    ChebyshevOperatorFinite<Interval> IFiniteOp(N, n, u0,
                                                capd::vectalg::convertObject<IMatrixType, DMatrixType>(g),
                                                capd::vectalg::convertObject<IVectorType, DVectorType>(v),
                                                capd::vectalg::convertObject<IVectorType, DVectorType>(w),
                                                multiIndices);

    IFiniteOp.setOmega( Interval(finiteOp.getOmega()));
    IFiniteOp.setASeries( capd::vectalg::convertObject<IVectorOfChebyshevsType, DVectorOfChebyshevsType>(finiteOp.getASeries()) );
    auto converted_X_approx = capd::vectalg::convertObject<IVectorType, DVectorType>(finiteOp.getX_approx());

    for (int j = 0; j < n; j++){
        Interval tmp, delta;
        u0[j].split(tmp, delta);
        tmp = IFiniteOp.getASeries()[j][0] + delta;
        IFiniteOp.setACoeff(tmp, j, 0);
        converted_X_approx[1 + j] = tmp;
    }

    IFiniteOp.setX_approx(converted_X_approx);

    //test na odwrocenie macierzy odwrotnej - odwrócona, powinna być zbliżona do IFiniteOp.derivative_finite
    try{
        auto I_inverse_of_derivative_test = capd::vectalg::convertObject<IMatrixType>(finiteOp.getInverseDerivativeFinite());
        IMatrixType interval_inverse_of_inverse_of_derivative_test = matrixAlgorithms::gaussInverseMatrix(I_inverse_of_derivative_test);
        cout << "Test na odwrócenie macierzy odwrotnej... - ZDANY" << endl;
    }catch(const std::runtime_error& e){
        std::cerr << "Błąd: " << e.what() << "\n";
        cout << "Test na odwrócenie macierzy odwrotnej... - NIEZDANY" << endl;
    }

    IVectorOfChebyshevsType IC_series = IFiniteOp.convertToSeriesFromXForm(IFiniteOp.compute_c(IFiniteOp.getX_approx()), 2*N-1);
    IFiniteOp.setCSeries(IC_series);
    IVectorType IF_x_approx = IFiniteOp(IFiniteOp.getX_approx());
    IFiniteOp.setF_x_approx(IF_x_approx);
    IFiniteOp.computeDerivativeInverse(IFiniteOp.getX_approx());

    return IFiniteOp;
}


double computeBoundsAndFindRadius(int N, int n, int N_g, double nu, ChebyshevOperatorFinite<Interval>& IFiniteOp) {

    RadiiPolynomials<Interval> radii_pol(N, n, nu, IFiniteOp);
    Norm<Interval> weighted_norm(nu, N, n);
    double r = radii_pol.findRForRadiiPolynomials();
    checkRadiiPolynomialsCoeffs(n, nu, radii_pol);

    return r;
}

void compareWithTaylorAndSaveResults(int N, double z, int step, int steps, double length_of_time_period, double a, double b,
                                     const IVectorType& u0,
                                     const IVectorOfChebyshevsType& a_series_approx,
                                     double omega_approx_mid,
                                     const DTimeMap::SolutionCurve& solution){
    double t_chebyshev = 0;
    double t_taylor = a;
    double del = 1.0 / (omega_approx_mid * 100);
    string name = "trajectoryN_" + std::to_string(N) + "_step_" + std::to_string(step)
            + "_NRsteps_" + std::to_string(steps) + "_delta_" + std::to_string(length_of_time_period)
            + "_start_point1" + ".csv";

    std::ofstream fout(name);
    fout << std::setprecision(17) << std::scientific;
    fout << "t_chebyshev,t_taylor,x_cheb,y_cheb,z_cheb,x_tay,y_tay,z_tay,diff_norm\n";

    DSumNormType sumNorm;
    IVectorType cheb;
    //C0HORect2Set
    while (true) {
        cheb = checkSolution(a_series_approx, t_chebyshev);
        auto cheb_d = capd::vectalg::convertObject<DVectorType>(cheb);
        auto tay = solution(t_taylor);
        double diff = sumNorm(cheb_d - tay);

        fout << t_chebyshev << "," << t_taylor << ","
             << cheb_d[0] << "," << cheb_d[1] << "," << cheb_d[2] << ","
             << tay[0] << "," << tay[1] << "," << tay[2] << "," << diff << "\n";

        if (t_chebyshev == 1.0 && t_taylor == b) break;
        t_chebyshev = capd::min(t_chebyshev + 0.01, 1.0);
        t_taylor = capd::min(t_taylor + del, b);
    }

    std::cout << "Zapisano plik z porównaniem Czebyszewa i Taylora dla kroku " << step << ": "  << name << std::endl;
}


ChebyshevOperatorFinite<double> findStartApproximation_PrepareChebyshevOperator_FindFiniteSolution(int N, int n, double rt,
                                                double a, double b,
                                                const IVectorType& _u0,
                                                DVectorType & v,
                                                DVectorType & w,
                                                DMatrixType& g,
                                                vector<vector<int>>& multiIndices,
                                                const DTimeMap::SolutionCurve& solution){
    double length_of_time_period = b - a;
    double del = rt/(N*20);

    std::vector<DVectorType> f_vals(N);
    std::vector<double> f_vals_x(N), f_vals_y(N), f_vals_z(N);

    LOGGER(a);
    LOGGER(b);
    for (int j = 0; j < N; ++j) {
        double theta = M_PI * (j + 0.5) / N;
        double xj = 0.5 * (a + b) + 0.5 * (b - a) * std::cos(theta);
        f_vals[j] = solution(xj);
        f_vals_x[j] = f_vals[j][0];
        f_vals_y[j] = f_vals[j][1];
        f_vals_z[j] = f_vals[j][2];
    }

    DChebyshevsVectorType coeffs_x = czebyszewCoefficientsFromAbValues(f_vals_x);
    DChebyshevsVectorType coeffs_y = czebyszewCoefficientsFromAbValues(f_vals_y);
    DChebyshevsVectorType coeffs_z = czebyszewCoefficientsFromAbValues(f_vals_z);
    double x = a;
    for (int i = 0; i < N; i++) {
        double t = (2 * x - (a + b)) / (b - a);
        double approx_x = coeffs_x(t);
        double approx_y = coeffs_y(t);
        double approx_z = coeffs_z(t);
        std::vector<double> approx_interpolation = {approx_x, approx_y, approx_z};
//        std::cout << "t=" << t << ", approx_inteprolation=" << approx_interpolation << ", solution" << solution(x) << "\n";
        x = capd::min(x + del, rt);
    }

    DVectorOfChebyshevsType a_series_start(n);
    a_series_start[0] = coeffs_x;
    a_series_start[1] = coeffs_y;
    a_series_start[2] = coeffs_z;

    // TODO: to zostawić do podrozdziału o wpływie a_series_start (wtedy dać tutaj rózne wartości i wytłumaczyć dlaczego postanowiliśmy poczatkowo wyliczać
//    DVectorOfChebyshevsType a_series_start(n);
//    for (int i = 0; i < n; i++){
//        a_series_start[i] = DChebyshevsVectorType(N);
//        a_series_start[i][0] = _u0[i];
//    }
//    a_series_start[0][1] = 1e-8; //zadanie, aby macierz pochodnej byla odwracalna (pierwsza kolumna byla niezerowa)
//    a_series_start[0][1] = 1e-13; //zadane, aby macierz pochodnej byla odwracalna (pierwsza kolumna byla niezerowa)
//    a_series_start[n-1][0] = w[n-1];

    auto u0 = capd::vectalg::convertObject<DVectorType>(_u0);
    constexpr double omega_start = 1.;

    w[n-1] = solution(b)[n-1]; // zadanie sekcji (po wyliczeniu po czasie)
    ChebyshevOperatorFinite<double> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 100;
    auto solution_approx = op.findFiniteSolution(omega_start, a_series_start, max_iterations);

//    cout << "Found approximate finite solution (omega, a) =" << endl;
//    cout << "omega = " << solution_approx.first << endl;
    cout << "a_series_approx = " << solution_approx.second << endl;
    IVectorType check_solution_for_interpolated = checkSolution(capd::vectalg::convertObject<IVectorOfChebyshevsType>(a_series_start), 1.);  // zakłada, że checkSolution sam robi 2*t - 1
    LOGGER(check_solution_for_interpolated);
    IVectorType check_solution_for_approx = checkSolution(capd::vectalg::convertObject<IVectorOfChebyshevsType>(solution_approx.second), 1.);  // zakłada, że checkSolution sam robi 2*t - 1
    LOGGER(check_solution_for_approx);
    return op;
}

// ---------- MAIN ----------
int main() {
    cout.precision(17);
    //parametry niezmienne
    constexpr int n = 3;
    constexpr int N_g = 2;
    auto multiIndices = generateMultiIndices(n, N_g);
    DMatrixType g = defineFunctionG(multiIndices, n);
    DVectorType v{0, 0, 1.};
    DVectorType w{0, 0, 0};

    //parametry zmienne
    int N = 28;
    double nu = 1.1;
    IVectorType u0{5.,5.,23.};
    double rt = 5.;
    double length_of_time_period = 0.25;
    int steps = rt/length_of_time_period;
    steps = 7; //TEMPORARY

    //wartosci wyliczane
    std::vector<IVectorType> poincare_section_target_points(steps);      //punkty po kazdym kroku
    IVectorType list_of_omegas(steps);                      //czasy po kazdym kroku
    DVectorType list_of_rs(steps);                          //promienie po kazdym kroku

    double a, b;
    auto u0_mid = capd::vectalg::convertObject<DVector>(u0);
    auto solution = findStartTaylorApproximation(N, rt, u0);

    for (int step = 0; step < steps; step++){
        cout << "======================== Iteracja " << step << "========================" << endl;
        a = step * length_of_time_period;
        b = (step + 1) * length_of_time_period;
        cout << "Interpoluję szereg Czebyszewa... \nWyznaczam operator Czebyszewa... \nWyliczam przybliżone rozwiązanie F(x)=0..." << endl;
        ChebyshevOperatorFinite<double> finiteOp =
                findStartApproximation_PrepareChebyshevOperator_FindFiniteSolution(N, n, rt, a, b, u0, v, w, g, multiIndices, solution);
        auto omega_approx = finiteOp.getOmega();
        auto a_series_approx = finiteOp.getASeries();
        cout << "\tPrzybliżony punkt po czasie " << 1./omega_approx << ": "
             << checkSolution(capd::vectalg::convertObject<IVectorOfChebyshevsType>(a_series_approx), 1.0) << endl;


        cout << "Konwersja na arytmetykę przedziałową..." << endl;
        ChebyshevOperatorFinite<Interval> IFiniteOp = convertToInterval(N, n, u0, v, w, g, multiIndices, finiteOp);

        cout << "Wyliczam Y i Z bounds i znajduję promień r..." << endl;
        auto r = computeBoundsAndFindRadius(N, n, N_g, nu, IFiniteOp);
        list_of_rs[step] = r;

        auto omega_interval = IFiniteOp.getOmega();
        auto a_series_interval = IFiniteOp.getASeries();
        Interval r_omega;
        double omega_interval_mid;
        omega_interval.split(omega_interval_mid, r_omega);
        
        u0 = checkSolution(a_series_interval, 1.0);
        for (int j = 0; j < n-1; j++){
            u0[j] = Interval(u0[j].leftBound() - r, u0[j].rightBound() + r);
        }
        u0[n-1] = Interval(w[n-1]);
        u0_mid = capd::vectalg::convertObject<DVector>(u0);

        compareWithTaylorAndSaveResults(N, w[n-1], step, steps, length_of_time_period, a, b, u0, a_series_interval, omega_interval_mid, solution);
        poincare_section_target_points[step] = u0;
        list_of_omegas[step] = 1./omega_approx;
        cout << "\tRozwiązanie dokładne znajduje się w kuli o promieniu: " << r << endl;
        cout << "\tPrzedział zawierajacy rozwiązanie dokładnie: " << u0 << endl;
        cout << "\tŚrodek przedziału: " << u0_mid << endl;



        std::string filename_intervals = "interval_and_midpoint_resultsN_" + std::to_string(N) +
                                         "_NRsteps_" + std::to_string(steps) + "_delta_" + std::to_string(length_of_time_period) + "_start_point1.csv";

        std::ofstream intervalOut(filename_intervals);
        if (intervalOut.is_open()) {
            intervalOut << std::setprecision(17) << std::scientific;
            intervalOut << "step,"
                        << "x_interval_left,x_interval_right,"
                        << "y_interval_left,y_interval_right,"
                        << "z_interval_left,z_interval_right,"
                        << "x_mid,y_mid,z_mid\n";

            for (int i = 0; i < steps; ++i) {
                intervalOut << i << ","
                            << u0[0].leftBound() << "," << u0[0].rightBound() << ","
                            << u0[1].leftBound() << "," << u0[1].rightBound() << ","
                            << u0[2].leftBound() << "," << u0[2].rightBound() << ","
                            << u0_mid[0] << "," << u0_mid[1] << "," << u0_mid[2] << "\n";
            }

            intervalOut.close();
            std::cout << "Zapisano dane przedziałów i środków do pliku: " << filename_intervals << std::endl;
        } else {
            std::cerr << "Nie udało się otworzyć pliku do zapisu przedziałów." << std::endl;
        }
    }

// =====================================================================================================================
    string filename = "radii_polynomial_resultsN_" + std::to_string(N) + "_NRsteps_" + std::to_string(steps) + "_delta_" + std::to_string(length_of_time_period) + "_start_point1" + ".csv";;
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        outFile << std::setprecision(17) << std::scientific;
        outFile << "r\n";
        for (int i = 0; i < steps; ++i) {
            outFile << list_of_rs[i] << "\n";
//            if (i != steps - 1) {
//                outFile << "";
//            }
        }
        outFile << std::endl;

        outFile.close();
        std::cout << "Zapisano dane do pliku: " << filename << std::endl;
    } else {
        std::cerr << "Nie udało się otworzyć pliku do zapisu." << std::endl;
    }
    return 0;
}
