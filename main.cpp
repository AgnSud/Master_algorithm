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

int attempt_nr = 27;


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

void compareWithTaylorAndSaveResults(int N, double z, double nu, int step, double base_time_step, double length_of_time_period, double a, double b,
                                     const IVectorType& u0,
                                     const IVectorOfChebyshevsType& a_series_approx,
                                     double omega_approx_mid,
                                     const DTimeMap::SolutionCurve& solution){
    double t_chebyshev = 0;
    double t_taylor = a;
    int number_of_points = 30;
    double del_tay = 1.0 / (omega_approx_mid * number_of_points);
    double del_cheb = 1.0 / number_of_points;

    string name = "trajectoryN_" + std::to_string(N) + "_step_" + std::to_string(step)
            + "_base_time_step_" + std::to_string(base_time_step) +
            + "_nu_" + std::to_string(nu) + "_attempt_nr_" + std::to_string(attempt_nr) + ".csv";

    std::ofstream fout(name);
    fout << std::setprecision(17) << std::scientific;
    fout << "t_chebyshev,t_taylor,x_cheb,y_cheb,z_cheb,x_tay,y_tay,z_tay,diff_norm\n";

    DSumNormType sumNorm;
    while (true) {
        auto cheb_d_X_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[0]);
        ChebyshevSeries<double> cheb_d_X(N);
        cheb_d_X.setCoefficients(cheb_d_X_tmp);

        auto cheb_d_Y_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[1]);
        ChebyshevSeries<double> cheb_d_Y(N);
        cheb_d_Y.setCoefficients(cheb_d_Y_tmp);

        auto cheb_d_Z_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[2]);
        ChebyshevSeries<double> cheb_d_Z(N);
        cheb_d_Z.setCoefficients(cheb_d_Z_tmp);

        double time_cheb_scaled = 2 * t_chebyshev - 1;
        auto cheb_x = cheb_d_X(time_cheb_scaled);
        auto cheb_y = cheb_d_Y(time_cheb_scaled);
        auto cheb_z = cheb_d_Z(time_cheb_scaled);
        DVectorType cheb_d = {cheb_x, cheb_y, cheb_z};
        auto tay = solution(t_taylor);
        double diff = sumNorm(cheb_d - tay);

        fout << t_chebyshev << "," << t_taylor << ","
             << cheb_d[0] << "," << cheb_d[1] << "," << cheb_d[2] << ","
             << tay[0] << "," << tay[1] << "," << tay[2] << "," << diff << "\n";

        if (t_chebyshev == 1.0 && t_taylor == b) break;
        t_chebyshev = capd::min(t_chebyshev + del_cheb, 1.0);
        t_taylor = capd::min(t_taylor + del_tay, b);
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
    double del = (b-a)/(N);

    std::vector<DVectorType> f_vals(N);
    std::vector<double> f_vals_x(N), f_vals_y(N), f_vals_z(N);

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
    cout << "===== Interpolacja dla przedziału [" << a << "," << b << "] =======" << endl;
    for (int i = 0; i < N; i++) {
        double t = (2 * x - (a + b)) / (b - a);
        double approx_x = coeffs_x(t);
        double approx_y = coeffs_y(t);
        double approx_z = coeffs_z(t);
        std::vector<double> approx_interpolation = {approx_x, approx_y, approx_z};
        x = capd::min(x + del, rt);
    }

    DVectorOfChebyshevsType a_series_start(n);
    a_series_start[0] = coeffs_x;
    a_series_start[1] = coeffs_y;
    a_series_start[2] = coeffs_z;

    auto u0 = capd::vectalg::convertObject<DVectorType>(_u0);
    constexpr double omega_start = 1.;

    w[n-1] = solution(b)[n-1]; // zadanie sekcji (po wyliczeniu po czasie)
    ChebyshevOperatorFinite<double> op(N, n, u0, g, v, w, multiIndices);
    int max_iterations = 100;
    auto solution_approx = op.findFiniteSolution(omega_start, a_series_start, max_iterations);

    cout << "a_series_approx = " << solution_approx.second << endl;
    IVectorType check_solution_for_interpolated = checkSolution(capd::vectalg::convertObject<IVectorOfChebyshevsType>(a_series_start), 1.);  // zakłada, że checkSolution sam robi 2*t - 1
    IVectorType check_solution_for_approx = checkSolution(capd::vectalg::convertObject<IVectorOfChebyshevsType>(solution_approx.second), 1.);  // zakłada, że checkSolution sam robi 2*t - 1
    return op;
}

int main() {
    cout.precision(17);
    constexpr int n = 3;
    constexpr int N_g = 2;
    auto multiIndices = generateMultiIndices(n, N_g);
    DMatrixType g = defineFunctionG(multiIndices, n);
    LOGGER(g);
    DVectorType v{0, 0, 1.};
    DVectorType w{0, 0, 0};

    int N = 80;
    double nu = 1.075;
    double rt = 5.0;

    IVectorType u0{5., 5., 23.};

    DVectorType list_of_rs(1000);
    std::vector<IVectorType> list_of_u0(1000);
    std::vector<DVectorType> list_of_u0_mid(1000);
    std::vector<double> list_of_omegas(1000);
    std::vector<double> list_of_used_dts(1000);

    double base_time_step = 0.55;
    std::vector<double> fallback_steps{0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05};

    auto solution = findStartTaylorApproximation(N, rt, u0);

    double current_time = 0.0;
    int global_step = 0;
    while (current_time < rt) {
        bool success = false;
        double dt = base_time_step;
        for (int try_idx = -1; try_idx < (int)fallback_steps.size(); ++try_idx) {
            if (try_idx >= 0)
                dt = fallback_steps[try_idx];
            double a = current_time;
            double b = std::min(current_time + dt, rt);
            try {
                cout << "==== Global step " << global_step << " (" << a << " -> " << b << ", dt=" << dt << ") ====" << endl;
                auto u0_mid = capd::vectalg::convertObject<DVector>(u0);

                ChebyshevOperatorFinite<double> finiteOp =
                        findStartApproximation_PrepareChebyshevOperator_FindFiniteSolution(N, n, rt, a, b, u0, v, w, g, multiIndices, solution);

                auto omega_approx = finiteOp.getOmega();
                auto a_series_approx = finiteOp.getASeries();

                auto cheb_d_X_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[0]);
                ChebyshevSeries<double> cheb_d_X(N);
                cheb_d_X.setCoefficients(cheb_d_X_tmp);

                auto cheb_d_Y_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[1]);
                ChebyshevSeries<double> cheb_d_Y(N);
                cheb_d_Y.setCoefficients(cheb_d_Y_tmp);

                auto cheb_d_Z_tmp = capd::vectalg::convertObject<DVectorType>(a_series_approx[2]);
                ChebyshevSeries<double> cheb_d_Z(N);
                cheb_d_Z.setCoefficients(cheb_d_Z_tmp);

                auto cheb_x = cheb_d_X(1.0);
                auto cheb_y = cheb_d_Y(1.0);
                auto cheb_z = cheb_d_Z(1.0);
                DVectorType cheb_d = {cheb_x, cheb_y, cheb_z};

                cout << "\tPrzybliżony punkt po czasie " << 1. / omega_approx << ": "
                     << cheb_d << endl;

                ChebyshevOperatorFinite<Interval> IFiniteOp = convertToInterval(N, n, u0, v, w, g, multiIndices, finiteOp);
                auto r = computeBoundsAndFindRadius(N, n, N_g, nu, IFiniteOp);

                auto omega_interval = IFiniteOp.getOmega();
                auto a_series_interval = IFiniteOp.getASeries();
                Interval r_omega;
                double omega_interval_mid;
                omega_interval.split(omega_interval_mid, r_omega);

                u0 = checkSolution(a_series_interval, 1.0);
                for (int j = 0; j < n - 1; j++) {
                    u0[j] = Interval(u0[j].leftBound() - r, u0[j].rightBound() + r);
                }
                u0[n - 1] = Interval(w[n - 1]);
                u0_mid = capd::vectalg::convertObject<DVector>(u0);
                cout << "\tRozwiązanie dokładne znajduje się w kuli o promieniu: " << r << endl;
                cout << "\tPrzedział zawierający rozwiązanie dokładnie: " << u0 << endl;

                compareWithTaylorAndSaveResults(N, w[n - 1], nu, global_step, base_time_step, dt, a, b, u0, a_series_interval, omega_interval_mid, solution);

                list_of_rs[global_step] = r;
                list_of_u0[global_step] = u0;
                list_of_u0_mid[global_step] = u0_mid;
                list_of_omegas[global_step] = 1. / omega_approx;
                list_of_used_dts[global_step] = dt;


                success = true;
                break;
            } catch (const std::exception& e) {
                std::cerr << "Błąd przy próbie z krokiem dt = " << dt << ": " << e.what() << std::endl;
            }
        }

        if (!success) {
            std::cout << "Nie udało się wykonać kroku z żadnym dopuszczalnym dt. Kończę program." << std::endl;
            break;
        }
        current_time += dt;
        global_step++;
    }


    // ======= ZAPISUJEMY WYNIKI =======

    std::string filename_intervals = "u0_interval_and_midpoint_N" + std::to_string(N) +
            "_base_time_step" + std::to_string(base_time_step) + "_global_step_" + std::to_string(global_step) +
            "_nu_" + std::to_string(nu) + "_attempt_nr_" + std::to_string(attempt_nr) +  ".csv";
    std::ofstream intervalOut(filename_intervals);
    if (intervalOut.is_open()) {
        intervalOut << std::setprecision(17) << std::scientific;
        intervalOut << "step,x_interval_left,x_interval_right,y_interval_left,y_interval_right,z_interval_left,z_interval_right,x_mid,y_mid,z_mid,dt\n";
        for (int i = 0; i < global_step; ++i) {
            intervalOut << i;
            if (std::isnan(list_of_rs[i])) {
                intervalOut << ",NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN\n";
                break;
            } else {
                intervalOut << "," << list_of_u0[i][0].leftBound() << "," << list_of_u0[i][0].rightBound()
                            << "," << list_of_u0[i][1].leftBound() << "," << list_of_u0[i][1].rightBound()
                            << "," << list_of_u0[i][2].leftBound() << "," << list_of_u0[i][2].rightBound()
                            << "," << list_of_u0_mid[i][0]
                            << "," << list_of_u0_mid[i][1]
                            << "," << list_of_u0_mid[i][2]
                            << "," << list_of_used_dts[i]
                            << "\n";
            }
        }
        intervalOut.close();
        std::cout << "Zapisano dane przedziałów i środków do pliku: " << filename_intervals << std::endl;
    } else {
        std::cerr << "Nie udało się otworzyć pliku do zapisu przedziałów." << std::endl;
    }

    std::string filename_radii = "radii_polynomial_N" + std::to_string(N) +
                                 "_base_time_step" + std::to_string(base_time_step) + "_global_step_" + std::to_string(global_step) +
                                 "_nu_" + std::to_string(nu) + "_attempt_nr_" + std::to_string(attempt_nr) + ".csv";
    std::ofstream outFile(filename_radii);
    if (outFile.is_open()) {
        outFile << std::setprecision(17) << std::scientific;
        outFile << "r\n";
        for (const auto& r : list_of_rs) {
            outFile << r << "\n";
        }
        outFile.close();
        std::cout << "Zapisano dane promieni do pliku: " << filename_radii << std::endl;
    } else {
        std::cerr << "Nie udało się otworzyć pliku do zapisu promieni." << std::endl;
    }

    return 0;
}

