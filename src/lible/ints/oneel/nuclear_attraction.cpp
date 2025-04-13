#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>

namespace LIO = lible::ints::one;

using std::array, std::vector;

template <>
void LIO::kernel<LIO::Option::nuclear_attraction>(const int la, const int lb,
                                                  const ShellPairData &sp_data,
                                                  vec2d &ints_out)
{
    int lab = la + lb;
    // vec3d rints(lab + 1, 0);
    // vec3d rints_sum(lab + 1, 0);
    // vec4d rints_tmp(lab + 1, 0);
    // vector<double> fnx(lab + 1, 0);

    vector<vector<vec4d>> ecoeffs = ecoeffsSPData_Eijt(la, lb, sp_data);

    BoysGrid boys_grid(lab);

    int dim_a_cart = numCartesians(sp_data.la);
    int dim_b_cart = numCartesians(sp_data.lb);

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    arma::dmat ints_contracted(dim_a_cart, dim_b_cart);
    arma::dmat ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

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

        const vector<vec4d> &Exyz = ecoeffs[ipair];
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double p = a + b;

                array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                       (a * xyz_a[1] + b * xyz_b[1]) / p,
                                       (a * xyz_a[2] + b * xyz_b[2]) / p};

                // rints_sum.set(0);
                vec3d rints_sum(lab + 1, 0);
                for (int iatom = 0; iatom < sp_data.n_atoms; iatom++)
                {
                    const auto &coords = sp_data.atomic_coords;
                    arma::vec::fixed<3> xyz_c{coords[3 * iatom],
                                              coords[3 * iatom + 1],
                                              coords[3 * iatom + 2]};

                    array<double, 3> xyz_pc{xyz_p[0] - xyz_c[0], xyz_p[1] - xyz_c[1],
                                            xyz_p[2] - xyz_c[2]};

                    double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                    double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
                    double x = p * xyz_pc_dot;

                    vector<double> fnx = calcBoysF(lab, x, boys_grid);

                    vec3d rints = calcRInts(lab, p, &xyz_pc[0], &fnx[0]);

                    // calcRInts_(la, lb, p, xyz_pc, fnx, rints_tmp, rints);

                    int Z = sp_data.atomic_nrs[iatom];
                    for (int t = 0; t <= lab; t++)
                        for (int u = 0; u <= lab; u++)
                            for (int v = 0; v <= lab; v++)
                                rints_sum(t, u, v) += Z * rints(t, u, v);
                }

                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];
                double dadb = da * db;
                double fac = 2 * M_PI / p * dadb;

                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                    ints_contracted(mu, nu) += -fac *
                                                               Exyz[iab](0, i, i_, t) *
                                                               Exyz[iab](1, j, j_, u) *
                                                               Exyz[iab](2, k, k_, v) *
                                                               rints_sum(t, u, v);
            }

        ints_sph = sph_trafo_bra * ints_contracted * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph, ints_out);
    }
}