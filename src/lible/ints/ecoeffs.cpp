#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <array>

namespace LI = lible::ints;

using std::array, std::vector;

void LI::coeffs(const double a, const double b, const double PA, const double PB,
                const double one_o_2p, const int la, const int lb, vec3d &E)
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

void LI::coeffs(const double a, const double b, const int la, const int lb,
                const array<double, 3> &A, const array<double, 3> &B,
                const array<double, 3> &Kab, vec3d &Ex, vec3d &Ey, vec3d &Ez)
{
    Ex.set(0);
    Ey.set(0);
    Ez.set(0);

    Ex(0, 0, 0) = Kab[0];
    Ey(0, 0, 0) = Kab[1];
    Ez(0, 0, 0) = Kab[2];

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs(a, b, PA[0], PB[0], one_o_2p, la, lb, Ex);
    coeffs(a, b, PA[1], PB[1], one_o_2p, la, lb, Ey);
    coeffs(a, b, PA[2], PB[2], one_o_2p, la, lb, Ez);
}

void LI::calcECoeffs(const int l, const vector<double> &exps,
                     vector<arma::dmat> &ecoeffs_out)
{
    // Calculates Hermite expansion coefficients for a diagonal shell pair located on the same
    // atom.

    int dim = dimCartesians(l);
    size_t k = exps.size();

    ecoeffs_out.resize(k * k, arma::zeros(dim, dim));

    lible::vec3d Ex(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d Ey(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d Ez(l + 1, l + 1, 2 * l + 1, 0);

    array<double, 3> A{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            std::array<double, 3> Kab{1, 1, 1};

            coeffs(a, b, l, l, A, A, Kab, Ex, Ey, Ez);

            for (const auto [i, j, k, mu] : cart_exps_a)
                for (const auto [i_, j_, k_, nu] : cart_exps_a)
                    ecoeffs_out[iab](mu, nu) = Ex(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
        }
}

void LI::calcECoeffs(const int la, const int lb, const ShellPairData_new &sp_data,
                     vector<vector<vec3d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec3d> ecoeffs_pair(dim_a * dim_b, vec3d(3, la + 1, lb + 1, 0));
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                coeffs(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](0, i, j) = Ex(i, j, 0);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](1, i, j) = Ey(i, j, 0);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](2, i, j) = Ez(i, j, 0);
            }

        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LI::calcECoeffs(const int la, const int lb, const ShellPairData_new &sp_data,
                     vector<vector<vec4d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec4d> ecoeffs_pair(dim_a * dim_b, vec4d(3, la + 1, lb + 1, la + lb + 1, 0));
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                coeffs(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](0, i, j, t) = Ex(i, j, t);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](1, i, j, t) = Ey(i, j, t);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](2, i, j, t) = Ez(i, j, t);
            }

        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LI::calcECoeffsSpherical(const int la, const int lb, const ShellPairData_new &sp_data,
                              vector<double> &ecoeffs_out, vector<double> &ecoeffs_tsp_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    arma::dmat sph_trafo_a = returnSphericalTrafo(la);
    arma::dmat sph_trafo_b = returnSphericalTrafo(lb);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int dim_a_cart = dimCartesians(la);
    int dim_b_cart = dimCartesians(lb);
    int dim_a_sph = dimSphericals(la);
    int dim_b_sph = dimSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = dimHermiteGaussians(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vec3i tuv_poss = returnTUVPoss(lab);
    vector<IdxsTUV> idxs_tuv = returnIdxsTUV(lab);

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];

        vec3d ecoeffs_ppair_cc(dim_a_cart, dim_b_cart, dim_tuv, 0);
        vec3d ecoeffs_ppair_sc(dim_a_sph, dim_b_cart, dim_tuv, 0);
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                coeffs(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                ecoeffs_ppair_cc.set(0);
                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);
                                    ecoeffs_ppair_cc(mu, nu, tuv) =
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }

                ecoeffs_ppair_sc.set(0);
                for (int mu = 0; mu < dim_a_sph; mu++)
                    for (int mu_ = 0; mu_ < dim_a_cart; mu_++)
                        for (int nu = 0; nu < dim_b_cart; nu++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_ppair_sc(mu, nu, tuv) +=
                                    ecoeffs_ppair_cc(mu_, nu, tuv) * sph_trafo_a(mu, mu_);

                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];
                int pos = offset_ecoeffs + iab * n_ecoeffs;
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int nu_ = 0; nu_ < dim_b_cart; nu_++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                            {
                                double ecoeff = da * db * ecoeffs_ppair_sc(mu, nu_, tuv) *
                                                sph_trafo_b(nu, nu_);

                                int idx = pos + munu * dim_tuv + tuv;
                                ecoeffs_out[idx] += ecoeff;

                                int idx_tsp = pos + tuv * dim_ab + munu;
                                ecoeffs_tsp_out[idx_tsp] += ecoeff;
                            }
            }
    }
}