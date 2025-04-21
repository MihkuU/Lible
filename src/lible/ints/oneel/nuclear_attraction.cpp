#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using std::array, std::vector;

namespace lible::ints
{
    // Forward declaration
    void externalChargesKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                               const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                               const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                               const vector<array<double, 4>> &charges, const BoysGrid &boys_grid,
                               double *ints_contracted);

    // Forward declaration
    void externalChargesDerivKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                                    const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                                    const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                                    const vector<array<double, 4>> charges, const BoysGrid &boys_grid,
                                    double *ints_contracted);

    // Forward declaration
    void externalChargesOperatorDerivKernel(const int la, const int lb, const int cdepth_a,
                                            const int cdepth_b, const double *cexps_a,
                                            const double *cexps_b, const double *ccoeffs_a,
                                            const double *ccoeffs_b, const double *xyz_a,
                                            const double *xyz_b, 
                                            const vector<array<double, 4>> charges,
                                            const BoysGrid &boys_grid,
                                            double *intderivs_contracted);
}

void LI::externalChargesKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                               const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                               const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                               const vector<array<double, 4>> &charges, const BoysGrid &boys_grid,
                               double *ints_contracted)
{
    int lab = la + lb;
    int n_a_cart = numCartesians(la);
    int n_b_cart = numCartesians(lb);
    int n_ab_cart = n_a_cart * n_b_cart;

    std::fill(ints_contracted, ints_contracted + n_ab_cart, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = 2 * (M_PI / p) * dadb;

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                   (a * xyz_a[1] + b * xyz_b[1]) / p,
                                   (a * xyz_a[2] + b * xyz_b[2]) / p};

            vec3d rints_sum(lab + 1, 0);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab, x, boys_grid);

                vec3d rints = calcRInts(lab, p, &xyz_pc[0], &fnx[0]);

                for (int t = 0; t <= lab; t++)
                    for (int u = 0; u <= lab; u++)
                        for (int v = 0; v <= lab; v++)
                            rints_sum(t, u, v) += charge * rints(t, u, v); 
            }

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                                ints_contracted[mu * n_b_cart + nu] += (-1) * fac * // -1 = charge of electron
                                                                       Ex(i, i_, t) *
                                                                       Ey(j, j_, u) *
                                                                       Ez(k, k_, v) * 
                                                                       rints_sum(t, u, v);
        }
}

void LI::externalChargesDerivKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                                    const double *cexps_a, const double *cexps_b,
                                    const double *ccoeffs_a, const double *ccoeffs_b,
                                    const double *xyz_a, const double *xyz_b,
                                    const vector<array<double, 4>> charges, const BoysGrid &boys_grid,
                                    double *intderivs_contracted)
{
    int lab = la + lb;
    int n_a_cart = numCartesians(la);
    int n_b_cart = numCartesians(lb);
    int n_ab_cart = n_a_cart * n_b_cart;

    std::fill(intderivs_contracted, intderivs_contracted + 6 * n_ab_cart, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = 2 * (M_PI / p) * dadb;

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                           {Ex, Ey, Ez});

            array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                   (a * xyz_a[1] + b * xyz_b[1]) / p,
                                   (a * xyz_a[2] + b * xyz_b[2]) / p};

            vec3d rints_sum(lab + 2, 0);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

                vec3d rints = calcRInts(lab + 1, p, &xyz_pc[0], &fnx[0]);

                for (int t = 0; t <= lab + 1; t++)
                    for (int u = 0; u <= lab + 1; u++)
                        for (int v = 0; v <= lab + 1; v++)
                            rints_sum(t, u, v) += charge * rints(t, u, v);
            }

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                double dpx = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t + 1, u, v);
                                double dpy = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u + 1, v);
                                double dpz = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v + 1);

                                double drx = fac * E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v);
                                double dry = fac * Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v);
                                double drz = fac * Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v) * rints_sum(t, u, v);

                                int munu = mu * n_b_cart + nu;

                                // d/dA
                                intderivs_contracted[0 * n_ab_cart + munu] += -1 * ((a / p) * dpx + drx); // -1 = charge of electron
                                intderivs_contracted[1 * n_ab_cart + munu] += -1 * ((a / p) * dpy + dry);
                                intderivs_contracted[2 * n_ab_cart + munu] += -1 * ((a / p) * dpz + drz);

                                // d/dB
                                intderivs_contracted[3 * n_ab_cart + munu] += -1 * ((b / p) * dpx - drx);
                                intderivs_contracted[4 * n_ab_cart + munu] += -1 * ((b / p) * dpy - dry);
                                intderivs_contracted[5 * n_ab_cart + munu] += -1 * ((b / p) * dpz - drz);
                            }
        }
}

void LI::externalChargesOperatorDerivKernel(const int la, const int lb, const int cdepth_a,
                                            const int cdepth_b, const double *cexps_a,
                                            const double *cexps_b, const double *ccoeffs_a,
                                            const double *ccoeffs_b, const double *xyz_a,
                                            const double *xyz_b,
                                            const vector<array<double, 4>> charges,
                                            const BoysGrid &boys_grid,
                                            double *intderivs_contracted)
{
    int lab = la + lb;
    int n_a_cart = numCartesians(la);
    int n_b_cart = numCartesians(lb);
    int n_ab_cart = n_a_cart * n_b_cart;

    std::fill(intderivs_contracted, intderivs_contracted + 6 * n_ab_cart, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = 2 * (M_PI / p) * dadb;

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                           {Ex, Ey, Ez});

            array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                   (a * xyz_a[1] + b * xyz_b[1]) / p,
                                   (a * xyz_a[2] + b * xyz_b[2]) / p};

            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

                vec3d rints = calcRInts(lab + 1, p, &xyz_pc[0], &fnx[0]);

                for (const auto &[i, j, k, mu] : cart_exps_a)
                    for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    double Exyz = Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);

                                    int munu = mu * n_b_cart + nu;
                                    int idx1 = icharge * n_ab_cart + munu + 0;
                                    int idx2 = icharge * n_ab_cart + munu + 1;
                                    int idx3 = icharge * n_ab_cart + munu + 2;

                                    intderivs_contracted[idx1] += charge * fac * Exyz * rints(t + 1, u, v);
                                    intderivs_contracted[idx2] += charge * fac * Exyz * rints(t, u + 1, v);
                                    intderivs_contracted[idx3] += charge * fac * Exyz * rints(t, u, v + 1);
                                }
            }


        }
}

template <>
void LIO::kernel<LIO::Option::nuclear_attraction>(const int la, const int lb,
                                                  const ShellPairData &sp_data,
                                                  vec2d &ints_out)
{
    vector<array<double, 4>> charges(sp_data.n_atoms);
    for (int iatom = 0; iatom < sp_data.n_atoms; iatom++)
    {
        const auto &coords = sp_data.atomic_coords;
        array<double, 3> xyz_c{coords[3 * iatom],
                               coords[3 * iatom + 1],
                               coords[3 * iatom + 2]};

        double Z = sp_data.atomic_nrs[iatom];
        charges[iatom] = {xyz_c[0], xyz_c[1], xyz_c[2], Z};
    }

    int lab = la + lb;
    BoysGrid boys_grid(lab);

    int dim_a_cart = numCartesians(sp_data.la);
    int dim_b_cart = numCartesians(sp_data.lb);

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    arma::dmat ints_contracted(dim_b_cart, dim_a_cart);
    arma::dmat ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        externalChargesKernel(la, lb, cdepth_a, cdepth_b, &sp_data.exps[cofs_a], &sp_data.exps[cofs_b],
                              &sp_data.coeffs[cofs_a], &sp_data.coeffs[cofs_b],
                              &sp_data.coords[6 * ipair + 0], &sp_data.coords[6 * ipair + 3],
                              charges, boys_grid, ints_contracted.memptr());

        ints_sph = sph_trafo_bra * ints_contracted.t() * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph, ints_out);
    }
}