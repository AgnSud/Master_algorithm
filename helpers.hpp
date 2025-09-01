#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <capd/capdlib.h>
#include "source/ChebyshevSeries/ChebyshevSeries.hpp"
#include "source/ChebyshevSeries/Norm.hpp"

using namespace capd;
using namespace std;

typedef vectalg::Vector<double, DIMENSION> DVectorType;
typedef vectalg::Matrix<double, DIMENSION, DIMENSION> DMatrixType;
typedef vectalg::Vector<ChebyshevSeries<double, DIMENSION>, DIMENSION> DVectorOfChebyshevsType;
typedef ChebyshevSeries<double, DIMENSION> DChebyshevsVectorType;

typedef vectalg::Vector<long double, DIMENSION> LDVectorType;
typedef vectalg::Matrix<long double, DIMENSION, DIMENSION> LDMatrixType;
typedef vectalg::Vector<ChebyshevSeries<long double, DIMENSION>, DIMENSION> LDVectorOfChebyshevsType;
typedef ChebyshevSeries<long double, DIMENSION> LDChebyshevsVectorType;

typedef vectalg::Vector<Interval, DIMENSION> IVectorType;
typedef vectalg::Matrix<Interval, DIMENSION, DIMENSION> IMatrixType;
typedef vectalg::Vector<ChebyshevSeries<Interval, DIMENSION>, DIMENSION> IVectorOfChebyshevsType;
typedef ChebyshevSeries<Interval, DIMENSION> IChebyshevsVectorType;

typedef DTimeMap::SolutionCurve DTimeMapSolution;


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

LDMatrixType defineFunctionG(vector<vector<int>> multiIndices, int n){
    int A = multiIndices.size();
    long double rho = 28;
    long double beta = 8/3.;
    long double sigma = 10;
    LDMatrixType g(A, n);

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

    return g;
}

IVectorType checkSolution(const IVectorOfChebyshevsType& a_series_approx,
                          long double argument){
    IVectorType value(a_series_approx.dimension());
    for(int i = 0; i < a_series_approx.dimension(); i++){
        value[i] = a_series_approx[i](argument * 2 - 1);
    }
    return value;
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

DChebyshevsVectorType czebyszewCoefficientsFromAbValues(const std::vector<long double>& f_vals) {
    int N = f_vals.size();
    LDChebyshevsVectorType a_k(N);

    std::vector<long double> theta(N);
    for (int j = 0; j < N; ++j) {
        theta[j] = M_PI * (j + 0.5) / N;
    }

    for (int k = 0; k < N; ++k) {
        long double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            sum += f_vals[j] * std::cos(k * theta[j]);
        }
        a_k[k] = sum / N;
    }

    return a_k;
}

LDTimeMap::SolutionCurve findStartTaylorApproximation(int N, long double rt, const IVectorType& _u0){
    LDMap vectorField("par:q;var:x,y,z;fun:10*(y-x),x*(28-z)-y,x*y-8*z/3;");
    LDOdeSolver solver(vectorField, N);
    LDTimeMap tm(solver);

    auto u0_mid = capd::vectalg::convertObject<LDVectorType>(_u0);
    cout << "Starting computing to time " << rt << endl;
    tm.stopAfterStep(true);
    int counter = 0;
    LDTimeMap::SolutionCurve solution(0.);
    do {
        tm(rt, u0_mid, solution);
        counter++;
        if (counter % 50 == 0)
            cout << "check, counter=" << counter << endl;
    } while (!tm.completed());
    int nr_of_points = 250;
    long double acc_t = 0;

    // === ZAPIS DO PLIKU ===
    std::ofstream outFile("solution_output.csv");
    if (outFile.is_open()) {
        outFile << "x_tay,y_tay,z_tay\n";
        for (int i = 0; i < nr_of_points; i++) {
            for (int j = 0; j < 3; ++j) {
                outFile << solution(acc_t)[j];
                if (j < 2)
                    outFile << ",";
            }
            acc_t = std::min(rt, acc_t + (rt/nr_of_points));
            outFile << "\n";
        }
        outFile.close();
        std::cout << "Solution saved to solution_output.txt" << std::endl;
    } else {
        std::cerr << "Error opening file for writing." << std::endl;
    }

    return solution;
}

#endif // HELPERS_HPP
