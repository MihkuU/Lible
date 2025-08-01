#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/oneel/oneel_detail.hpp>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using std::array, std::vector;

namespace lible::ints
{
    vec2d externalCharges(const vector<array<double, 4>> &charges,
                          const Structure &structure);

    vec2d externalChargesErf(const vector<array<double, 4>> &charges,
                             const std::vector<double> &omegas,
                             const Structure &structure);

    vec2d externalChargesKernel(const int ipair, const vector<array<double, 4>> &charges,
                                const BoysGrid &boys_grid, const ShellPairData &sp_data);

    vec2d externalChargesErfKernel(const int ipair,
                                   const vector<array<double, 4>> &charges,
                                   const vector<double> &omegas,
                                   const BoysGrid &boys_grid,
                                   const ShellPairData &sp_data);

    array<vec2d, 3> spinOrbitCoupling1El(const Structure &structure);

    array<vec2d, 3> spinOrbitCoupling1ElKernel(const int ipair,
                                               const vector<array<double, 4>> &charges,
                                               const BoysGrid &boys_grid,
                                               const ShellPairData &sp_data);

    array<vec2d, 3> spinOrbitCoupling1ElKernelDebug(const int ipair,
                                                    const vector<array<double, 4>> &charges,
                                                    const BoysGrid &boys_grid,
                                                    const ShellPairData &sp_data);

    array<vec2d, 6> externalChargesD1Kernel(const int ipair,
                                            const vector<array<double, 4>> &charges,
                                            const BoysGrid &boys_grid,
                                            const ShellPairData &sp_data);

    vector<array<vec2d, 3>>
    externalChargesOperatorD1Kernel(const int ipair, const vector<array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);
}

lible::vec2d LI::externalChargesKernel(const int ipair, const vector<array<double, 4>> &charges,
                                       const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    const int cofs_a = sp_data.coffsets[2 * ipair + 0];
    const int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int n_cart_a = numCartesians(la);
    int n_cart_b = numCartesians(lb);

    vec2d ints_cart(Fill(0), n_cart_a, n_cart_b);
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

            vec3d rints_sum(Fill(0), lab + 1);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab, x, boys_grid);

                vec3d rints = calcRInts3D(lab, p, &xyz_pc[0], &fnx[0]);

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
                                ints_cart(mu, nu) += (-1) * fac * // -1 = charge of electron
                                                     Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                                     rints_sum(t, u, v);
        }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms[ofs_norm_a + mu];
            double norm_b = sp_data.norms[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

lible::vec2d LI::externalChargesErfKernel(const int ipair,
                                          const vector<array<double, 4>> &charges,
                                          const vector<double> &omegas,
                                          const BoysGrid &boys_grid,
                                          const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;
    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int n_cart_a = numCartesians(la);
    int n_cart_b = numCartesians(lb);

    vec2d ints_cart(Fill(0), n_cart_a, n_cart_b);
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

            vec3d rints_sum(Fill(0), lab + 1);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;

                double omega_factor = omegas[icharge] / std::pow(omegas[icharge] * omegas[icharge] + p, 0.5);

                double omega_squared = omegas[icharge] * omegas[icharge];

                // second part of (52) in https://doi.org/10.1039/B605188J
                double x = p * xyz_pc_dot * omega_squared / (omega_squared + p);

                vector<double> fnx = calcBoysF(lab, x, boys_grid);

                vec3d rints = calcRInts3DErf(lab, p, omegas[icharge], &xyz_pc[0], &fnx[0]);

                for (int t = 0; t <= lab; t++)
                    for (int u = 0; u <= lab; u++)
                        for (int v = 0; v <= lab; v++)
                            rints_sum(t, u, v) += charge * omega_factor // erf-related factors
                                                  * rints(t, u, v);
            }

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                                ints_cart(mu, nu) += (-1) * fac * // -1 = charge of electron
                                                     Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                                     rints_sum(t, u, v);
        }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms[ofs_norm_a + mu];
            double norm_b = sp_data.norms[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

array<lible::vec2d, 6>
LI::externalChargesD1Kernel(const int ipair, const vector<array<double, 4>> &charges,
                            const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;
    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int n_cart_a = numCartesians(la);
    int n_cart_b = numCartesians(lb);

    array<vec2d, 6> ints_cart;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), n_cart_a, n_cart_b);

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

            vec3d rints_sum(Fill(0), lab + 2);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

                vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

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

                                // d/dA
                                ints_cart[0](mu, nu) += -1 * ((a / p) * dpx + drx); // -1 = charge of electron
                                ints_cart[1](mu, nu) += -1 * ((a / p) * dpy + dry);
                                ints_cart[2](mu, nu) += -1 * ((a / p) * dpz + drz);

                                // d/dB
                                ints_cart[3](mu, nu) += -1 * ((b / p) * dpx - drx);
                                ints_cart[4](mu, nu) += -1 * ((b / p) * dpy - dry);
                                ints_cart[5](mu, nu) += -1 * ((b / p) * dpz - drz);
                            }
        }

    array<vec2d, 6> ints_sph;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];
                ints_sph[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

vector<array<lible::vec2d, 3>>
LI::externalChargesOperatorD1Kernel(const int ipair, const vector<array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    const int cofs_a = sp_data.coffsets[2 * ipair + 0];
    const int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    const int n_cart_a = numCartesians(la);
    const int n_cart_b = numCartesians(lb);

    int n_charges = charges.size();
    vector<array<vec2d, 3>> ints_cart(n_charges);
    for (int icharge = 0; icharge < n_charges; icharge++)
        for (int icoord = 0; icoord < 3; icoord++)
            ints_cart[icharge][icoord] = vec2d(Fill(0), n_cart_a, n_cart_b);

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

            for (int icharge = 0; icharge < n_charges; icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

                vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

                for (const auto &[i, j, k, mu] : cart_exps_a)
                    for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    double Exyz = Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);

                                    ints_cart[icharge][0](mu, nu) += charge * fac * Exyz * rints(t + 1, u, v);
                                    ints_cart[icharge][1](mu, nu) += charge * fac * Exyz * rints(t, u + 1, v);
                                    ints_cart[icharge][2](mu, nu) += charge * fac * Exyz * rints(t, u, v + 1);
                                }
            }
        }

    vector<array<vec2d, 3>> ints_sph(n_charges);
    for (int icharge = 0; icharge < n_charges; icharge++)
    {
        ints_sph[icharge][0] = trafo2Spherical(la, lb, ints_cart[icharge][0]);
        ints_sph[icharge][1] = trafo2Spherical(la, lb, ints_cart[icharge][1]);
        ints_sph[icharge][2] = trafo2Spherical(la, lb, ints_cart[icharge][2]);
    }

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (int icharge = 0; icharge < n_charges; icharge++)
        for (size_t mu = 0; mu < ints_sph[icharge][0].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icharge][0].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];

                ints_sph[icharge][0](mu, nu) *= norm_a * norm_b;
                ints_sph[icharge][1](mu, nu) *= norm_a * norm_b;
                ints_sph[icharge][2](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

lible::vec2d LI::externalCharges(const vector<array<double, 4>> &charges,
                                 const Structure &structure)
{
    int l_max = structure.getMaxL();
    int dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data = shellPairDataSymm(la, lb, structure);

            int lab = la + lb;
            BoysGrid boys_grid(lab);

            for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
            {
                vec2d ints_ipair = externalChargesKernel(ipair, charges, boys_grid, sp_data);

                int n_sph_a = 2 * la + 1;
                int n_sph_b = 2 * lb + 1;
                int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
                int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                    {
                        int mu = ofs_a + a;
                        int nu = ofs_b + b;
                        ints(mu, nu) = ints_ipair(a, b);
                        ints(nu, mu) = ints_ipair(a, b);
                    }
            }
        }

    return ints;
}

lible::vec2d LI::externalChargesErf(const vector<array<double, 4>> &charges,
                                    const std::vector<double> &omegas,
                                    const Structure &structure)
{
    int l_max = structure.getMaxL();
    int dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data = shellPairDataSymm(la, lb, structure);

            int lab = la + lb;
            BoysGrid boys_grid(lab);

            for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
            {
                vec2d ints_ipair = externalChargesErfKernel(ipair, charges,
                                                            omegas, boys_grid,
                                                            sp_data);

                int n_sph_a = 2 * la + 1;
                int n_sph_b = 2 * lb + 1;
                int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
                int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                    {
                        int mu = ofs_a + a;
                        int nu = ofs_b + b;
                        ints(mu, nu) = ints_ipair(a, b);
                        ints(nu, mu) = ints_ipair(a, b);
                    }
            }
        }

    return ints;
}

array<lible::vec2d, 3> LI::spinOrbitCoupling1El(const Structure &structure)
{
    int l_max = structure.getMaxL();
    int dim_ao = structure.getDimAO();

    vector<array<double, 4>> charges = structure.getZs();

    array<vec2d, 3> ints;
    for (int i = 0; i < 3; i++)
        ints[i] = vec2d(Fill(0), dim_ao, dim_ao);

    for (int la = 0; la <= l_max; la++)
        for (int lb = 0; lb <= l_max; lb++)
        {
            ShellPairData sp_data = shellPairDataNoSymm(la, lb, structure);

            int lab = la + lb;
            BoysGrid boys_grid(lab + 1);

            for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
            {
                array<vec2d, 3> ints_batch = spinOrbitCoupling1ElKernel(ipair, charges, boys_grid,
                                                                        sp_data);

                int n_sph_a = 2 * la + 1;
                int n_sph_b = 2 * lb + 1;
                int ofs_a = sp_data.offsets_sph[2 * ipair];
                int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                    {
                        int mu = ofs_a + a;
                        int nu = ofs_b + b;
                        ints[0](mu, nu) = ints_batch[0](a, b);
                        ints[1](mu, nu) = ints_batch[1](a, b);
                        ints[2](mu, nu) = ints_batch[2](a, b);
                    }
            }
        }

    return ints;
}

array<lible::vec2d, 3> LI::spinOrbitCoupling1ElKernel(const int ipair,
                                                      const vector<array<double, 4>> &charges,
                                                      const BoysGrid &boys_grid,
                                                      const ShellPairData &sp_data)
{
    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    const int cofs_a = sp_data.coffsets[2 * ipair + 0];
    const int cofs_b = sp_data.coffsets[2 * ipair + 1];
    const int n_cart_a = numCartesians(la);
    const int n_cart_b = numCartesians(lb);

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    array<vec2d, 3> ints_cart;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), n_cart_a, n_cart_b);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                           {Ex, Ey, Ez});

            // R integrals
            double p = a + b;
            array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                   (a * xyz_a[1] + b * xyz_b[1]) / p,
                                   (a * xyz_a[2] + b * xyz_b[2]) / p};

            vec3d rints_sum(Fill(0), lab + 2);
            for (size_t icharge = 0; icharge < charges.size(); icharge++)
            {
                auto [xc, yc, zc, charge] = charges[icharge];

                array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                double x = p * xyz_pc_dot;

                vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

                vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

                for (int t = 0; t <= lab + 1; t++)
                    for (int u = 0; u <= lab + 1; u++)
                        for (int v = 0; v <= lab + 1; v++)
                            rints_sum(t, u, v) += charge * rints(t, u, v);
            }

            // (\nabla_A x \nabla_B) integrals
            double dadb = da * db;
            double fac = (2 * M_PI / p) * dadb;
            double a_p = a / p;
            double b_p = b / p;
            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                double r100 = rints_sum(t + 1, u, v);
                                double r010 = rints_sum(t, u + 1, v);
                                double r001 = rints_sum(t, u, v + 1);
                                double e100 = E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                double e010 = Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                                double e001 = Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);

                                double pr_yz = r010 * e001;
                                double pr_zy = r001 * e010;
                                double pr_zx = r001 * e100;
                                double pr_xz = r100 * e001;
                                double pr_xy = r100 * e010;
                                double pr_yx = r010 * e100;

                                double ab_yz = -a_p * pr_yz + b_p * pr_zy;
                                double ab_zy = -a_p * pr_zy + b_p * pr_yz;
                                double ab_zx = -a_p * pr_zx + b_p * pr_xz;
                                double ab_xz = -a_p * pr_xz + b_p * pr_zx;
                                double ab_xy = -a_p * pr_xy + b_p * pr_yx;
                                double ab_yx = -a_p * pr_yx + b_p * pr_xy;

                                // -1 = charge of electron
                                ints_cart[0](mu, nu) += -1 * fac * (ab_yz - ab_zy);
                                ints_cart[1](mu, nu) += -1 * fac * (ab_zx - ab_xz);
                                ints_cart[2](mu, nu) += -1 * fac * (ab_xy - ab_yx);
                            }
        }

    array<vec2d, 3> ints_sph;
    for (int icoord = 0; icoord < 3; icoord++)
        ints_sph[icoord] = trafo2Spherical(la, lb, ints_cart[icoord]);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (int icoord = 0; icoord < 3; icoord++)
        for (size_t mu = 0; mu < ints_sph[icoord].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icoord].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];
                ints_sph[icoord](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

array<lible::vec2d, 3> LI::spinOrbitCoupling1ElKernelDebug(const int ipair,
                                                           const vector<array<double, 4>> &charges,
                                                           const BoysGrid &boys_grid,
                                                           const ShellPairData &sp_data)
{
    const double delta = 1e-5; // Finite differences with a fixed value.

    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    const int cofs_a = sp_data.coffsets[2 * ipair + 0];
    const int cofs_b = sp_data.coffsets[2 * ipair + 1];
    const int n_cart_a = numCartesians(la);
    const int n_cart_b = numCartesians(lb);
    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    auto kernel_fun = [&](const array<double, 6> &coords) -> vec2d
    {
        array<double, 3> xyz_a{coords[0], coords[1], coords[2]};
        array<double, 3> xyz_b{coords[4], coords[5], coords[6]};

        vec2d ints_cart(Fill(0), n_cart_a, n_cart_b);
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

                auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, &xyz_a[0], &xyz_b[0]);

                array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                       (a * xyz_a[1] + b * xyz_b[1]) / p,
                                       (a * xyz_a[2] + b * xyz_b[2]) / p};

                vec3d rints_sum(Fill(0), lab + 1);
                for (size_t icharge = 0; icharge < charges.size(); icharge++)
                {
                    auto [xc, yc, zc, charge] = charges[icharge];

                    array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

                    double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                    double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                    double x = p * xyz_pc_dot;

                    vector<double> fnx = calcBoysF(lab, x, boys_grid);

                    vec3d rints = calcRInts3D(lab, p, &xyz_pc[0], &fnx[0]);

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
                                    ints_cart(mu, nu) += (-1) * fac * // -1 = charge of electron
                                                         Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                                         rints_sum(t, u, v);
            }

        vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

        int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
        int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
        for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];
                ints_sph(mu, nu) *= norm_a * norm_b;
            }

        return ints_sph;
    };

    array<vec2d, 3> ints_sph;
    for (int icart = 0; icart < 3; icart++)
        ints_sph[icart] = vec2d(Fill(0), n_sph_a, n_sph_b);

    array<double, 6> coords{xyz_a[0], xyz_a[1], xyz_a[2],
                            xyz_b[0], xyz_b[1], xyz_b[2]};

    array<std::pair<int, int>, 6> idxs{{
        {1, 5},
        {2, 4},
        {2, 3},
        {0, 5},
        {1, 3},
        {0, 4}}};

    array<vec2d, 6> terms;
    for (int iterm = 0; iterm < 6; iterm++)
    {
        auto [di, dj] = idxs[iterm];

        array<double, 6> coords_upup = coords;
        coords_upup[di] += delta;
        coords_upup[dj] += delta;

        array<double, 6> coords_dnup = coords;
        coords_dnup[di] -= delta;
        coords_dnup[dj] += delta;

        array<double, 6> coords_updn = coords;
        coords_updn[di] += delta;
        coords_updn[dj] -= delta;

        array<double, 6> coords_dndn = coords;
        coords_dndn[di] -= delta;
        coords_dndn[dj] -= delta;

        vec2d batch_upup = kernel_fun(coords_upup);
        vec2d batch_dnup = kernel_fun(coords_dnup);
        vec2d batch_updn = kernel_fun(coords_updn);
        vec2d batch_dndn = kernel_fun(coords_dndn);

        double denom = 4 * delta * delta;

        terms[iterm] = (1.0 / denom) * (batch_upup - batch_dnup -
                                        batch_updn + batch_dndn);
    }

    ints_sph[0] = terms[0] - terms[1];
    // ints_sph[1] = terms[2] - terms[3];
    // ints_sph[2] = terms[4] - terms[5];    

    return ints_sph;
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

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        vec2d ints_ipair = externalChargesKernel(ipair, charges, boys_grid, sp_data);

        int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
        int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
        for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
            {
                ints_out(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                ints_out(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
            }
    }
}

template <>
void LIO::kernel<LIO::Option::nuclear_attraction_erf>(const int la,
                                                      const int lb,
                                                      const ShellPairData &sp_data,
                                                      const std::vector<double> &omegas,
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

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        vec2d ints_ipair = externalChargesErfKernel(ipair, charges, omegas, boys_grid, sp_data);

        int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
        int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
        for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
            {
                ints_out(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                ints_out(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
            }
    }
}