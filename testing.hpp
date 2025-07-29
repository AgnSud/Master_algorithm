#ifndef TESTING_HPP
#define TESTING_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <capd/capdlib.h>
#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/Norm.hpp"
#include "source/ChebyshevSeries/RadiiPolynomials.hpp"

using namespace capd;
using namespace std;

typedef vectalg::Vector<double, DIMENSION> DVectorType;
typedef vectalg::Vector<Interval, DIMENSION> IVectorType;
typedef vectalg::Vector<ChebyshevSeries<double, DIMENSION>, DIMENSION> DVectorOfChebyshevsType;
typedef ChebyshevSeries<double, DIMENSION> DChebyshevsVectorType;
typedef vectalg::Vector<ChebyshevSeries<Interval, DIMENSION>, DIMENSION> IVectorOfChebyshevsType;
typedef vectalg::Matrix<Interval, DIMENSION, DIMENSION> IMatrixType;



IVectorType randomVectorInBallWeightedNorm(int dim, Interval r, const Norm<Interval>& weighted_norm) {
    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> dist(0.0, 1.0);
    std::uniform_real_distribution<> u_dist(0.0, 1.0);  // dla losowego skalowania <1

    DVectorType vec(dim);
    for (int i = 0; i < dim; ++i)
        vec[i] = dist(gen);

    IVectorType ivec = capd::vectalg::convertObject<IVectorType, DVectorType>(vec);
    Interval norm_val = weighted_norm.computeNorm(ivec);
    double u = std::pow(u_dist(gen), 1.0 / dim);
    Interval scale = Interval(r * u) / norm_val;

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
}

void verify_bounds(int N, int n, double nu, Interval r, RadiiPolynomials<Interval> radii_pol,
                   ChebyshevOperatorFinite<Interval>& IFiniteOp) {
    IMatrixType A_N = IFiniteOp.getInverseDerivativeFinite();
    IVectorType x_approx = IFiniteOp.getX_approx();
    IVectorType F_x_approx = IFiniteOp.getF_x_approx();
    IVectorType AF = A_N * F_x_approx;
    int l = 1;
    IVectorType Pi1j_Z1j_norm_supremum(n);
    Interval Pi0_DT_supremum(0.);

    Norm<Interval> weighted_norm(nu, N, n);

    std::cout << "============ Verifying bounds ============\n";
    auto Y_bounds = radii_pol.getYBounds();
    auto Z_bounds = radii_pol.getZBounds();

    for (int j = 0; j < n; ++j) {
        auto Pi1j_Y1j = radii_pol.Pi1_j(AF, j, N, n);
        auto Pi1j_Y1j_norm = weighted_norm.computeNorm(Pi1j_Y1j);
        std::cout << "Y_bound for j=" << j << ": " << Pi1j_Y1j_norm << " ≤ " << Y_bounds[j+1].rightBound()
                  << " => " << (Pi1j_Y1j_norm.rightBound() <= Y_bounds[j+1].rightBound() ? "OK" : "FAIL") << "\n";
    }

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

#endif //TESTING_HPP
