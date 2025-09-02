#pragma once

#include <numeric>
#include "RadiiPolynomials.hpp"

template<typename T>
RadiiPolynomials<T>::RadiiPolynomials(int N_, int n_, long double nu_, const ChebyshevOperatorFinite<T>& finiteOp_)
        : N(N_), n(n_), nu(nu_), finiteOp(finiteOp_), Y_bounds(n_+1), Z_bounds(n_+1) {}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
void RadiiPolynomials<T>::testOperatorNorm() {
    Norm<T> weighted_norm(nu, N, n);
    auto tmp = weighted_norm.computeOperatorNorm(finiteOp.getInverseDerivativeFinite());
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::g_unit_vector(int j) {
    std::vector<int> e_j(n, 0);
    e_j[j] = 1;
    auto m_idx = finiteOp.getMultiIndices();
    for (int i = 0; i < m_idx.size(); i++){
        if (e_j == m_idx[i]){
            return finiteOp.getG()[i];
        }
    }
    throw std::runtime_error("g_unit_vector: unit vector not found in multiIndices");
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::g_ls(int l, int s) {
    std::vector<int> e_ls(n, 0);
    e_ls[l] += 1;
    e_ls[s] += 1;
    auto m_idx = finiteOp.getMultiIndices();
    for (int i = 0; i < m_idx.size(); i++){
        if (e_ls == m_idx[i]){
            return finiteOp.getG()[i];
        }
    }
    throw std::runtime_error("g_ls: unit vector not found in multiIndices");
}

template <typename T>
template <class V>
V RadiiPolynomials<T>::vector_abs(const V& v) {
    V result(v.dimension());
    for (int j = 0; j < v.dimension(); j++)
        result[j] = capd::abs(v[j]);
    return result;
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::getYBounds() {
    return Y_bounds;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::getZBounds() {
    return Z_bounds;
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
T RadiiPolynomials<T>::Pi0(const VectorType& x) {
    return x[0];
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::Pi1(const VectorType& x) {
    VectorType result(x.dimension() - 1);
    for (int i = 1; i < x.dimension(); ++i) {
        result[i - 1] = x[i];
    }
    return result;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::PiN(const VectorType& a, int N_, int n_) {
    VectorType result(N_ * n_);
    for (int i = 0; i < N_ * n_; ++i) {
        result[i] = a[i];
    }
    return result;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::PiN_x(const VectorType& x, int N_, int n_) {
    VectorType result(N_ * n_ + 1);
    result[0] = x[0];
    VectorType a(N_ * n_);
    for (int i = 0; i < a.dimension(); ++i) {
        a[i] = x[i + 1];
    }

    VectorType projected = PiN(a, N_, n_);
    for (int i = 0; i < projected.dimension(); ++i) {
        result[i + 1] = projected[i];
    }
    return result;
}

template<typename T>
template<class V>
V RadiiPolynomials<T>::Pi1_j(const V &x, int j, int N_, int n_) {
    return ChebyshevOperatorFinite<T>::getCoeffVectorI_thSquareParan(x, j, N_, n_);
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
void RadiiPolynomials<T>::compute_YBounds(int N_g) {
    Y_bounds[0] = computeY0().right();
    for (int j = 0; j < n; ++j) {
        Y_bounds[j+1] = computeY1j(j, N_g).right();
    }
}

template <typename T>
T RadiiPolynomials<T>::computeY0() {
    VectorType x_approx = finiteOp.getX_approx();               // x^*
    VectorType F_x_approx = finiteOp.getF_x_approx();           // F_N(x^*)
    MatrixType A_N = finiteOp.getInverseDerivativeFinite();     // A_N
    VectorType AF = A_N * F_x_approx;                           // A_N F_N(x^*)

    return capd::abs(Pi0(AF)); // Pi_0 — tylko pierwsza współrzędna
}

template <typename T>
T RadiiPolynomials<T>::computeY1j(int j, int N_g) {
    T omega = finiteOp.getOmega();
    auto c_series = finiteOp.getCSeries();
    const int k_start = N;
    const int k_end = N_g * (N - 1) + 1;

    T horner_sum = capd::abs<T>(c_series[j][k_end - 1]) / k_end;
    for (int k = k_end - 1; k >= k_start; --k) {
        T diff = c_series[j][k - 1];
        if (k + 1 < c_series[j].dimension())
            diff -= c_series[j][k + 1];
        diff = capd::abs<T>(diff);

        horner_sum = horner_sum * nu + diff / k;
    }
    horner_sum = horner_sum * std::pow(nu, k_start);

    VectorType x_approx = finiteOp.getX_approx();
    VectorType F_x_approx = finiteOp.getF_x_approx();
    auto A_N = finiteOp.getInverseDerivativeFinite();
    VectorType AF = Pi1_j(A_N * F_x_approx, j, N, n);

    Norm<T> weighted_norm(nu, N, n);
    auto w_norm_AF = weighted_norm.computeNorm(AF);

    return horner_sum / (2 * omega) + w_norm_AF;
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
void RadiiPolynomials<T>::compute_ZBounds(T r) {
    auto [Z0_r2, Z0_r1] = compute_Z0_terms();
    Z_bounds[0] = Z0_r2 * r * r + Z0_r1 * r;

    for (int j = 0; j < n; ++j) {
        auto [Z1j_r2, Z1j_r1] = compute_Z1j_terms(j);
        Z_bounds[j + 1] = Z1j_r2 * r * r + Z1j_r1 * r;
    }
}

template <typename T>
std::pair<T, T> RadiiPolynomials<T>::compute_Z0_terms() {
    T gamma = compute_gamma();
    VectorType h = compute_h();

    MatrixType A_N = finiteOp.getInverseDerivativeFinite();
    Norm<T> weighted_norm(nu, N, n);
    T op_norm = weighted_norm.computeOperatorNorm_Pi0(A_N);

    VectorType Z1 = compute_Z1();
    T Pi0_Z1 = Pi0(Z1);

    T r2_coeff = gamma * op_norm;
    T r_coeff = h[0] + Pi0_Z1;

    return std::make_pair(r2_coeff, r_coeff);
}

template <typename T>
std::pair<T, T> RadiiPolynomials<T>::compute_Z1j_terms(int j) {
    T gamma = compute_gamma();
    T omega = finiteOp.getOmega();

    VectorType d1_vec = compute_d1();
    VectorType d2_vec = compute_d2();
    VectorType Z1 = compute_Z1();
    VectorType h = compute_h();

    MatrixType A_N = finiteOp.getInverseDerivativeFinite();
    Norm<T> weighted_norm(nu, N, n);
    T AN_op_norm = weighted_norm.computeOperatorNorm_Pi1j(A_N, j);

    T Z1_norm = weighted_norm.computeNorm(Pi1_j(Z1, j, N, n));

    T r2_coeff = gamma * AN_op_norm + d2_vec[j] / (omega * 4.0) + 2 / omega;
    T r_coeff = h[j + 1] + Z1_norm + (d1_vec[j] / (4.0 * omega));

    return std::make_pair(r2_coeff, r_coeff);
}

///////////////////////////////////////////////////////////////////////////////

/// obliczy wszystkie h = [h0, h1x, h1y, hjz]
template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_h() {
    VectorType result_h(1 + n);
    MatrixType identity_matrix(1 + n * N, 1 + n * N);
    identity_matrix.setToIdentity();

    MatrixType DF_N = finiteOp.getDerivativeFinite();
    MatrixType A_N = finiteOp.getInverseDerivativeFinite();
    MatrixType diff = identity_matrix - A_N * DF_N;
    Norm<T> weighted_norm(nu, N, n);

    result_h[0] = weighted_norm.computeOperatorNorm_Pi0(diff);
    for (int j = 0; j < n; j++){
        result_h[j+1] = weighted_norm.computeOperatorNorm_Pi1j(diff, j);
    }

    return result_h;
}

///////////////////////////////////////////////////////////////////////////////

// Final step: Z1 = |A_N| * Z1_tilde
template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_Z1() {
    auto Z1_tilde = compute_Z1_tilde();
    auto AN = finiteOp.getInverseDerivativeFinite();
    for (int i = 0; i < AN[0].dimension(); i++){
        for (int j = 0; j < AN[0].dimension(); j++) {
            AN[i][j] = capd::abs(AN[i][j]);
        }
    }
    auto Z1 = AN * Z1_tilde;
    return Z1;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_Z1_tilde() {
    int totalDim = 1 + n * N;
    VectorType Z1_tilde(totalDim);
    Z1_tilde[0] = T(0);
    for (int j = 0; j < n; j++){
        Z1_tilde[j + 1] = 1 / std::pow(nu, N);
    }

    int index = 1 + n;
    for (int k = 1; k < N; ++k) {
        auto B_k = computeB_k(k);
        for (int j = 0; j < n; ++j) {
            Z1_tilde[index] = B_k[j] * 0.25;
            index++;
        }
    }
    return Z1_tilde;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::computeB_k(int k){
    VectorType first_sum(n);
    auto a_series = finiteOp.getASeries();
    if (k == N - 1){
        for (int j = 0; j < n; j++){
            first_sum += vector_abs(g_unit_vector(j));
        }
        first_sum = first_sum / (std::pow(nu, N) * 2);
    }
    VectorType result = first_sum;
    for (int l = 0; l < n; l++){
        for (int s = l; s < n; s++){
            auto abs_g_ls = vector_abs(g_ls(l, s));
            if (is_nonzero(abs_g_ls)){
                auto operatorNormPsi_al_k = operatorNormPsi_ak(a_series[l], k);
                auto operatorNormPsi_as_k = operatorNormPsi_ak(a_series[s], k);
                result = result + abs_g_ls * (operatorNormPsi_al_k + operatorNormPsi_as_k);
            }
        }
    }
    return result;
}

template <typename T>
T RadiiPolynomials<T>::operatorNormPsi_ak(VectorType& a, int k) {
    T max_result = capd::max(capd::abs(a[N-1]) / std::pow(nu, k+N),
                             capd::abs(a[capd::abs(N-2)]) / std::pow(nu, k+N-1));

    for (int l = N; l <= k+N-2; l++){
        auto diff = capd::abs(a[capd::abs(k-l-1)] - a[capd::abs(k-l+1)]);
        diff = diff / std::pow(nu, l);
        max_result = capd::max(max_result, diff);
    }
    return max_result / 2;
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
T RadiiPolynomials<T>::compute_gamma() {
    VectorType sum_g_for_alpha_2(n);
    auto multiIndices = finiteOp.getMultiIndices();
    auto g = finiteOp.getG();

    for (int alpha = 0; alpha < multiIndices.size(); alpha++){
        auto sum_alpha = std::reduce(multiIndices[alpha].begin(), multiIndices[alpha].end());
        if (sum_alpha == 2){
            VectorType abs_g_alpha(n);
            for (int j = 0; j < n; j++)
                abs_g_alpha[j] = capd::abs(g[alpha][j]);
            sum_g_for_alpha_2 += abs_g_alpha;
        }
    }
    sum_g_for_alpha_2 = ((2 * nu * nu + 1) * sum_g_for_alpha_2) / (2 * nu);
    for (int j = 0; j < n; j++)
        sum_g_for_alpha_2[j] += 2*(N-1);

    return *std::max_element(sum_g_for_alpha_2.begin(), sum_g_for_alpha_2.end());
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_d1() {
    auto multiIndices = finiteOp.getMultiIndices();
    auto a_series = finiteOp.getASeries();

    // (nu^2 + 1) / (nu * N)
    T coeff = (nu * nu + 1) / (nu * N);

    // Suma po |g_{e_j}|
    VectorType sum_g_ej(n);
    for (int j = 0; j < n; ++j) {
        sum_g_ej += vector_abs(g_unit_vector(j));
    }

    VectorType part1 = coeff * sum_g_ej;

    // Suma po g_{ls}
    VectorType part2(n);
    for (int l = 0; l < n; ++l) {
        for (int s = l; s < n; ++s) {
            auto g_ls_abs = vector_abs(g_ls(l, s));

            T Gamma_l_minus = compute_GammaMinus_a(a_series[l]);
            T Gamma_s_minus = compute_GammaMinus_a(a_series[s]);
            T Gamma_l_plus  = compute_GammaPlus_a(a_series[l]);
            T Gamma_s_plus  = compute_GammaPlus_a(a_series[s]);

            T inner = nu * (Gamma_l_minus + Gamma_s_minus)
                      + (1.0 / nu) * (Gamma_l_plus + Gamma_s_plus);
            part2 += g_ls_abs * inner;
        }
    }

    return part1 + part2;
}

template <typename T>
RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_d2() {
    auto g = finiteOp.getG();
    auto multiIndices = finiteOp.getMultiIndices();

    T coeff = 2 * (nu * nu + 1) / (nu * N);
    VectorType sum(n);

    for (int alpha = 0; alpha < multiIndices.size(); ++alpha) {
        int deg = std::accumulate(multiIndices[alpha].begin(), multiIndices[alpha].end(), 0);
        if (deg == 2) {
            VectorType abs_g_alpha(n);
            for (int j = 0; j < n; j++)
                abs_g_alpha[j] = capd::abs(g[alpha][j]);
            sum += abs_g_alpha;
        }
    }

    return coeff * sum;
}

template <typename T>
T RadiiPolynomials<T>::compute_GammaMinus_a(const VectorType& a) {
    T max_val = 2 * capd::abs(a[N - 1]) * std::pow(nu, N - 1) / N;

    for (int l = 1; l <= 2 * (N - 1); ++l) {
        T sum = 0;
        for (int k = N - 1 - l; k <= N - 1; ++k) {
            T term = capd::abs( a[capd::abs(k)] ) * std::pow(nu, k) / (k + l - 1);
            sum += term;
        }
        max_val = capd::max(max_val, sum);
    }

    return max_val;
}

template <typename T>
T RadiiPolynomials<T>::compute_GammaPlus_a(const VectorType& a) {
    T max_val = 0;

    for (int l = 2; l <= 2 * N; ++l) {
        T sum = 0;
        for (int k = N + 1 - l; k <= N - 1; ++k) {
            T term = capd::abs( a[capd::abs(k)] ) * std::pow(nu, k) / (k + l - 1);
            sum += term;
        }
        max_val = capd::max(max_val, sum);
    }

    return max_val;
}

///////////////////////////////////////////////////////////////////////////////

template <typename T>
long double RadiiPolynomials<T>::findRForRadiiPolynomials(){
    VectorType intervals(1 + n);
    intervals[0] = findRIntervalForRadiiPolynomials_0();
    for (int j = 0; j < n; j++){
        intervals[1 + j] = findRIntervalForRadiiPolynomials_1j(j);
    }
    long double r = intervals[0].leftBound();
    for (int j = 0; j < n; j++){
        auto left_r = intervals[j + 1].leftBound();
        if (left_r > r)
            r = left_r;
    }
    if (r < 0){
        throw std::runtime_error("RadiiPolynomials::findRForRadiiPolynomials: "
                                     "promien jest ujemny.");
    }
    return r;
}

template <typename T>
T RadiiPolynomials<T>::findRIntervalForRadiiPolynomials_0(){
    compute_YBounds(2);
    auto [A_coeff, B_coeff] = compute_Z0_terms();
    auto C_coeff = Y_bounds[0];
//        zamiana A_coeff, B_coeff, C_coeff na prawe końce przedziałów tych Z_bounds i Y_bounds
    A_coeff = A_coeff.right();
    B_coeff = B_coeff.right();
    C_coeff = C_coeff.right();

    B_coeff -= 1;
    auto delta = B_coeff * B_coeff - (4 * A_coeff * C_coeff);
    if (delta < 0)
        throw std::logic_error("Delta < 0 -> wielomian radii nie ma pierwiastkow");

    auto x1 = (-B_coeff - sqrt(delta)) / (2 * A_coeff);
    auto x2 = (-B_coeff + sqrt(delta)) / (2 * A_coeff);
    return Interval(x1.rightBound(), x2.leftBound());
}

template <typename T>
T RadiiPolynomials<T>::findRIntervalForRadiiPolynomials_1j(int j){
    compute_YBounds(2);
    auto [A_coeff, B_coeff] = compute_Z1j_terms(j);
    auto C_coeff = Y_bounds[j + 1];
    //    zamiana A_coeff, B_coeff, C_coeff na prawe końce przedziałów tych Z_bounds i Y_bounds
    A_coeff = A_coeff.right();
    B_coeff = B_coeff.right();
    C_coeff = C_coeff.right();

    B_coeff -= 1;
    auto delta = B_coeff * B_coeff - 4 * A_coeff * C_coeff;
    if (delta < 0)
        throw std::logic_error("Delta < 0 -> wielomian radii nie ma pierwiastkow");

    auto x1 = (-B_coeff - sqrt(delta)) / (2 * A_coeff);
    auto x2 = (-B_coeff + sqrt(delta)) / (2 * A_coeff);
    return Interval(x1.rightBound(), x2.leftBound());
}


template <typename T>
RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::operator()(T r) {
    compute_ZBounds(r);
    VectorType result(n + 1); // [p_0, p_{1,1}, ..., p_{1,n}]

    result[0] = Z_bounds[0] + Y_bounds[0] - r;
    for (int j = 0; j < n; ++j) {
        result[j + 1] = Z_bounds[j+1] + Y_bounds[j + 1] - r;
    }
    return result;
}

