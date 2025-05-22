#pragma once

#include <numeric>
#include "RadiiPolynomials.hpp"

template<typename T>
RadiiPolynomials<T>::RadiiPolynomials(int N_, int n_, double nu_, const ChebyshevOperatorFinite<T>& finiteOp_)
        : N(N_), n(n_), nu(nu_), h(1+n), finiteOp(finiteOp_) {}

template <typename T>
typename RadiiPolynomials<T>::DVectorType RadiiPolynomials<T>::g_unit_vector(int j) {
    std::vector<int> e_j(n, 0);
    e_j[j] = 1.0;
    auto m_idx = finiteOp.getMultiIndices();
    for (int i = 0; i < m_idx.size(); i++){
        if (e_j == m_idx[i]){
           return finiteOp.getG()[i];
        }
    }
    throw std::runtime_error("g_unit_vector: unit vector not found in multiIndices");
}

template <typename T>
typename RadiiPolynomials<T>::DVectorType RadiiPolynomials<T>::g_ls(int l, int s) {
//    DVectorType e_l(n, 0.0), e_s(n, 0.0);
//    e_l[l] = 1.0;
//    e_s[s] = 1.0;
    return g_unit_vector(l) + g_unit_vector(s);
}

template <typename T>
template <class V>
double RadiiPolynomials<T>::vector_sum(const V& v) {
//    typename V::ScalarType sum = 0;
    double sum = 0;
    /// czy tutaj niżej ma być abs???? odnosze sie do wzoru |alpha| = alpha_1 + ... + alpha_n
    for (const auto& val : v) sum += val;
//    cout << "sum = " << sum << endl;
    return capd::abs(sum);
}

template <typename T>
template <class V>
V RadiiPolynomials<T>::vector_abs(const V& v) {
    V result(v.dimension());
    for (int j = 0; j < v.dimension(); j++)
        result[j] = capd::abs(v[j]);
    return result;
}

template <typename T>
T RadiiPolynomials<T>::operatorNormPsi_ak(VectorType& a, int k) {
    T max_result = capd::max(capd::abs(a[N-1]) / std::pow(nu, k+N),
                             capd::abs(a[capd::abs(N-2)]) / std::pow(nu, k+N-1));
//    cout << "a_{N-2}/nu = " << capd::abs(a[N-1]) / std::pow(nu, k+N) << endl;
//    cout << "a_{N-1}/nu = " << capd::abs(a[capd::abs(N-2)]) / std::pow(nu, k+N-1)  << endl;
    for (int l = N; l <= k+N-2; l++){
        auto diff = capd::abs(a[capd::abs(k-l-1)] - a[capd::abs(k-l+1)]);
//        cout << "diff for l = " << l << " = " << diff << "BEFORE DIV" << endl;
        diff = diff / std::pow(nu, l);
//        cout << "diff for l = " << l << " = " << diff << "AFTER DIV" << endl;
        max_result = capd::max(max_result, diff);
    }
    return max_result / 2;
}



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

template <typename T>
typename RadiiPolynomials<T>::VectorType
RadiiPolynomials<T>::get_kth(const VectorType& a, int k){
    VectorType a_k(n);
    for (int i = 0; i < n; ++i){
        a_k[i] = a[k * n + i];
    }
    return a_k;
}

/// obliczy wszystkie h = [h0, h1x, h1y, hjz]
template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_h() {
    VectorType result_h(1 + n);
    MatrixType identity_matrix(1 + n * N, 1 + n * N);  // Utworzenie macierzy identycznościowej o odpowiednich rozmiarach
    identity_matrix.setToIdentity();  // Inicjalizacja macierzy jako macierz identycznościowa

    MatrixType DF_N = finiteOp.getDerivativeFinite();
    MatrixType A_N = finiteOp.getInverseDerivativeFinite();
//    DMatrixType identity_matrix = capd::vectalg::Matrix<T, 1+n*N, 1+n*N>::Identity(1+n*N);
    MatrixType diff = identity_matrix - A_N * DF_N;
    Norm<T> weighted_norm(nu, N, n);

    result_h[0] = weighted_norm.computeOperatorNorm_Pi0(diff);
    for (int j = 0; j < n; j++){
        result_h[j+1] = weighted_norm.computeOperatorNorm_Pi1j(diff, j);
    }

    this->h = result_h;
    return result_h;
}

/// zapis na gamme rozumiem rozlacznie? Jako że zbiór i maksymalny element
/// (np \Gamma_{1,2} + 2(N-1) + c \sum_{alpha = 2} |g_alpha|
/// element (2(N-1) jest jednowymiarowy, suma będzie n-wymiarowa)
template <typename T>
double RadiiPolynomials<T>::compute_gamma() {
    DVectorType sum_g_for_alpha_2(n);
    auto multiIndices = finiteOp.getMultiIndices();
    auto g = finiteOp.getG();

    for (int alpha = 0; alpha < multiIndices.size(); alpha++){
        auto sum_alpha = std::reduce(multiIndices[alpha].begin(), multiIndices[alpha].end());
        if (sum_alpha == 2){
//            cout << "|alpha| == 2 : " << multiIndices[alpha] << endl;
//            cout << "g_alpha : " << g[alpha] << endl;
            DVectorType abs_g_alpha(n);
            for (int j = 0; j < n; j++)
                abs_g_alpha[j] = capd::abs(g[alpha][j]);
//            cout << "abs_g_alpha : " << abs_g_alpha << endl;
            sum_g_for_alpha_2 += abs_g_alpha;
        }
    }
    sum_g_for_alpha_2 = ((2 * nu * nu + 1) * sum_g_for_alpha_2) / (2 * nu);
    for (int j = 0; j < n; j++)
        sum_g_for_alpha_2[j] += 2*(N-1);

    cout << "sum_g_for_alpha_2 = " << sum_g_for_alpha_2 << endl;
    return *std::max_element(sum_g_for_alpha_2.begin(), sum_g_for_alpha_2.end());
}


template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::computeB_k(int k){
    DVectorType first_sum(n);
    vectalg::Vector<ChebyshevSeries<T, DIMENSION>, DIMENSION> a_series = finiteOp.getASeries();
    if (k == N-1){
        for (int j = 0; j < n; j++){
//            first_sum = vector_sum(g_unit_vector(j));
            first_sum = vector_abs(g_unit_vector(j));
//            for (int j = 0; j < n; j++)
//                first_sum[j] = capd:i \:abs(g[alpha][j]);
            first_sum = first_sum / (std::pow(nu, N) * 2);
        }
    }

//    T result = T(first_sum);
    VectorType result = vectalg::convertObject<VectorType, DVectorType>(first_sum);
    cout << "first_sum = " << first_sum << endl;
    for (int l = 0; l < n; l++){
        auto operatorNormPsi_al_k = operatorNormPsi_ak(a_series[l], k);
        for (int s = l; s < n; s++){
            auto operatorNormPsi_as_k = operatorNormPsi_ak(a_series[s], k);
//            double sum_g_ls = vector_sum( g_ls(l, s) );
            auto abs_g_ls = vector_abs(g_ls(l, s));
            auto I_abs_g_ls = vectalg::convertObject<VectorType, DVectorType>(abs_g_ls);
//            cout << "||Psi_al_k|| + ||Psi_as_k|| =" << (operatorNormPsi_al_k + operatorNormPsi_as_k) << endl;
//            cout << "I_abs_g_ls, for l, s =  " << l << " " << s << " = "  << I_abs_g_ls << endl;
            result += I_abs_g_ls * (operatorNormPsi_al_k + operatorNormPsi_as_k);
        }
    }
    cout << "B_k for k = " << k << " = " << result << endl;
    cout << "===========================" << endl;
    return result;
}

template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_Z1_tilde() {
    int totalDim = 1 + n * N;
    VectorType Z1_tilde(totalDim);
    Z1_tilde[0] = T(0);

    int index = 1 + n; //od tego momentu zaczynamy
    for (int k = 1; k < N; ++k) {
        auto B_k = computeB_k(k);
        cout << "B_k = " << B_k << endl;
        for (int j = 0; j < n; ++j) {
            Z1_tilde[index++] = B_k[j] / 4;
        }
    }
    cout << "Z1_tilde = " << Z1_tilde << endl;
    return Z1_tilde;
}

// Final step: Z1 = |A_N| * Z1_tilde
template <typename T>
typename RadiiPolynomials<T>::VectorType RadiiPolynomials<T>::compute_Z1() {
    auto Z1_tilde = compute_Z1_tilde();
    auto AN = finiteOp.getInverseDerivativeFinite();  // assuming A_N is defined
    for (int i = 0; i < AN[0].dimension(); i++){
        for (int j = 0; j < AN[0].dimension(); j++) {
            AN[i][j] = capd::abs(AN[i][j]);
        }
    }
    return AN * Z1_tilde;
}


template <typename T>
T RadiiPolynomials<T>::computeY0() {

    VectorType x_approx = finiteOp.getX_approx();          // x^*
//    cout << "x_approx = " << x_approx << endl;
    VectorType F_x_approx = finiteOp.getF_x_approx();                        // F_N(x^*)
//    cout << "F_x_approx = " << F_x_approx << endl;
    MatrixType A_N = finiteOp.getInverseDerivativeFinite();   // A_N
//    cout << "A_N = " << A_N << endl;
    VectorType AF = A_N * F_x_approx;                         // A_N F_N(x^*)
//    cout << "A_N * F_x_approx = " << AF << endl;

    return capd::abs(Pi0(AF)); // Pi_0 — tylko pierwsza współrzędna
}

template <typename T>
T RadiiPolynomials<T>::computeY1j(int j, int N_g) {

    // Parametry
    T omega = finiteOp.getOmega();
    auto a_series = finiteOp.getASeries();
    auto c_series = finiteOp.getCSeries();

    T sum = 0;

    for (int k = N; k <= N_g * (N - 1) + 1; ++k) {
        T diff = c_series[j][k - 1];
        if (k + 1 < c_series[j].dimension())
            diff -= c_series[j][k + 1];
        diff = capd::abs<T>(diff);
//        cout << "k=" << k << ", [c_{k-1} - c_{k+1}]_" << j << " = " << diff << endl;
        T weighted = diff * std::pow(nu, k) / k;
//        cout << "weighted by nu^k/k =" << std::pow(nu, k) / k << " results: " << weighted << endl;
        sum += weighted;
    }

    // część druga: || Pi_{1,j} A_N F_N(x^*) ||_nu
    VectorType x_approx = finiteOp.getX_approx();
    VectorType F_x_approx = finiteOp.getF_x_approx();
    auto A_N = finiteOp.getInverseDerivativeFinite();
    VectorType AF = Pi1_j(A_N * F_x_approx, j, N, n);

    // Pi_1,j — współczynniki j-tego ciągu a_k (czyli co n-ta współrzędna od j)
    Norm<T> weighted_norm(nu, N, n);
    auto w_norm_AF = weighted_norm.computeNorm(AF);
//    cout << "w_norm_AF = " << w_norm_AF << endl;
//    cout << "sum = " << sum << endl;
    T result = sum / 2 * finiteOp.getOmega() + w_norm_AF;
    /// tutaj czesc sum/2*omega sprawia, ze przedzial nie jest od 0 (jest jak najmniejszy ale dodatni)

    return result;
}

template <typename T>
void RadiiPolynomials<T>::testOperatorNorm() {
    Norm<T> weighted_norm(nu, N, n);
    auto tmp = weighted_norm.computeOperatorNorm(finiteOp.getInverseDerivativeFinite());
}


