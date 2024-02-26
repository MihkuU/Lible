#include <lible/mcmurchie_davidson.hpp>
#include <lible/ints_util.hpp>
#include <lible/cart_exps.hpp>

namespace LIMD = lible::ints::MD;

using std::array, std::vector;

void LIMD::calcHCoeffs(const int &l, const vector<double> &exps,
                       vector<arma::dmat> &h_coeffs_out)
{
    // Calculates Hermite expansion coefficients for a diagonal shell pair located on the same
    // atom.

    int dim = dimCartesians(l);
    size_t k = exps.size();

    h_coeffs_out.resize(k * k, arma::zeros(dim, dim));

    lible::vec3d E_x(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d E_y(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d E_z(l + 1, l + 1, 2 * l + 1, 0);

    array<double, 3> A{0, 0, 0}; // dummy center

    auto &cart_exps_a = cart_exps[l];

    size_t iab = 0;
    for (size_t ia = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++)
        {
            double a = exps[ia];
            double b = exps[ib];

            std::array<double, 3> Kab{1, 1, 1};

            coeffs(a, b, l, l, A, A, Kab, E_x, E_y, E_z);

            for (auto &[i, j, k, mu] : cart_exps_a)
                for (auto &[i_, j_, k_, nu] : cart_exps_a)
                {
                    double E = E_x(i, i_, 0) * E_y(j, j_, 0) * E_z(k, k_, 0);
                    h_coeffs_out[iab](mu, nu) = E;
                }
            iab++;
        }
}

void LIMD::calcHCoeffs(const int &la, const int &lb, const ShellPairData &shell_pair_data,
                       vector<vector<arma::dmat>> &h_coeffs_out)
{
    // Calculates Hermite expansion coefficients for all shell-pairs corresponding to the
    // shell-pair data.

    h_coeffs_out.resize(shell_pair_data.n_pairs);

    int dim_a = dimCartesians(la);
    int dim_b = dimCartesians(lb);

    lible::vec3d E_x(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d E_y(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d E_z(la + 1, lb + 1, la + lb + 1, 0);

    auto &cart_exps_a = cart_exps[la];
    auto &cart_exps_b = cart_exps[lb];

    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        auto &exps = shell_pair_data.exps[ipair];

        array<double, 3> A = shell_pair_data.coords[ipair].first;
        array<double, 3> B = shell_pair_data.coords[ipair].second;

        size_t ka = exps.first.size();
        size_t kb = exps.second.size();

        vector<arma::dmat> h_coeffs_pair(ka * kb, arma::zeros(dim_a, dim_b));

        size_t iab = 0;
        for (size_t ia = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++)
            {
                double a = exps.first[ia];
                double b = exps.second[ib];

                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};                                   

                coeffs(a, b, la, lb, A, B, Kab, E_x, E_y, E_z);

                for (auto &[i, j, k, mu] : cart_exps_a)
                    for (auto &[i_, j_, k_, nu] : cart_exps_b)
                    {
                        double E = E_x(i, i_, 0) * E_y(j, j_, 0) * E_z(k, k_, 0);
                        h_coeffs_pair[iab](mu, nu) = E;
                    }
                iab++;
            }
        h_coeffs_out[ipair] = h_coeffs_pair;
    }
}

void LIMD::coeffs(const double &a, const double &b, const double &PA, const double &PB,
                  const double one_o_2p, const int &la, const int &lb, vec3d &E)
{
    for (int i = 1; i <= la; i++)
    {
        E(i, 0, 0) = PA * E(i - 1, 0, 0) + E(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            E(i, 0, t) = one_o_2p * E(i - 1, 0, t - 1) +
                         PA * E(i - 1, 0, t) +
                         (t + 1) * E(i - 1, 0, t + 1);

        E(i, 0, i) = one_o_2p * E(i - 1, 0, i - 1) + PA * E(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            E(i, j, 0) = PB * E(i, j - 1, 0) + E(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                E(i, j, t) = one_o_2p * E(i, j - 1, t - 1) +
                             PB * E(i, j - 1, t) +
                             (t + 1) * E(i, j - 1, t + 1);

            E(i, j, i + j) = one_o_2p * E(i, j - 1, i + j - 1) + PB * E(i, j - 1, i + j);
        }
}

void LIMD::coeffs(const double &a, const double &b, const int &la, const int &lb,
                  const std::array<double, 3> &A, const std::array<double, 3> &B,
                  const std::array<double, 3> &Kab,
                  vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    E_x.set(0);
    E_y.set(0);
    E_z.set(0);

    E_x(0, 0, 0) = Kab[0];
    E_y(0, 0, 0) = Kab[1];
    E_z(0, 0, 0) = Kab[2];

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs(a, b, PA[0], PB[0], one_o_2p, la, lb, E_x);
    coeffs(a, b, PA[1], PB[1], one_o_2p, la, lb, E_y);
    coeffs(a, b, PA[2], PB[2], one_o_2p, la, lb, E_z);
}