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
                                                  const ShellPairData &shell_pair_data,
                                                  vec2d &ints_out)
{
    int lab = la + lb;
    vec3d rints(lab + 1, 0);
    vec3d rints_sum(lab + 1, 0);
    vec4d rints_tmp(lab + 1, 0);
    vector<double> fnx(lab + 1, 0);

    vector<vector<vec4d>> ecoeffs;
    calcECoeffs(la, lb, shell_pair_data, ecoeffs);

    BoysF boys_f(lab);

    int dim_a_cart = dimCartesians(shell_pair_data.la);
    int dim_b_cart = dimCartesians(shell_pair_data.lb);

    auto cart_exps_a = cart_exps[la];
    auto cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(shell_pair_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(shell_pair_data.lb).t();

    arma::dmat ints_contracted(dim_a_cart, dim_b_cart);
    arma::dmat ints_sph;
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
        const auto &[xyz_a, xyz_b] = shell_pair_data.coords[ipair];

        arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
        arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        const vector<vec4d> &Exyz = ecoeffs[ipair];

        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double p = a + b;

                arma::vec::fixed<3> P = (a * A + b * B) / p;

                rints_sum.set(0);
                for (size_t iatom = 0; iatom < shell_pair_data.structure->getNAtoms(); iatom++)
                {
                    int Z = shell_pair_data.structure->getZ(iatom);
                    array<double, 3> xyz_c = shell_pair_data.structure->getCoordsAtom(iatom);
                    arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};

                    arma::vec::fixed<3> RPC = P - C;
                    double RPC2 = arma::dot(RPC, RPC);
                    double x = p * RPC2;

                    boys_f.calcFnx(lab, x, fnx);

                    calcRInts(la, lb, p, RPC, fnx, rints_tmp, rints);

                    for (int t = 0; t <= lab; t++)
                        for (int u = 0; u <= lab; u++)
                            for (int v = 0; v <= lab; v++)
                                rints_sum(t, u, v) += Z * rints(t, u, v);
                }

                double dadb = ccoeffs_a[ia] * ccoeffs_b[ib];
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

        transferIntegrals(ipair, shell_pair_data, ints_sph, ints_out);
    }
}

template <>
void LIO::kernel_new<LIO::Option::nuclear_attraction>(const int la, const int lb,
                                                      const ShellPairData_new &sp_data,
                                                      vec2d &ints_out)
{
    int lab = la + lb;
    vec3d rints(lab + 1, 0);
    vec3d rints_sum(lab + 1, 0);
    vec4d rints_tmp(lab + 1, 0);
    vector<double> fnx(lab + 1, 0);

    vector<vector<vec4d>> ecoeffs;
    calcECoeffs(la, lb, sp_data, ecoeffs);

    BoysF boys_f(lab);

    int dim_a_cart = dimCartesians(sp_data.la);
    int dim_b_cart = dimCartesians(sp_data.lb);

    auto cart_exps_a = cart_exps[la];
    auto cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    arma::dmat ints_contracted(dim_a_cart, dim_b_cart);
    arma::dmat ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        arma::vec::fixed<3> xyz_a{sp_data.coords[6 * ipair + 0],
                                  sp_data.coords[6 * ipair + 1],
                                  sp_data.coords[6 * ipair + 2]};

        arma::vec::fixed<3> xyz_b{sp_data.coords[6 * ipair + 3],
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

                arma::vec::fixed<3> xyz_p = (a * xyz_a + b * xyz_b) / p;

                rints_sum.set(0);
                for (int iatom = 0; iatom < sp_data.n_atoms; iatom++)
                {
                    const auto &coords = sp_data.atomic_coords;
                    arma::vec::fixed<3> xyz_c{coords[3 * iatom],
                                              coords[3 * iatom + 1],
                                              coords[3 * iatom + 2]};

                    arma::vec::fixed<3> xyz_pc = xyz_p - xyz_c;
                    double xyz_pc_dot = arma::dot(xyz_pc, xyz_pc);
                    double x = p * xyz_pc_dot;

                    boys_f.calcFnx(lab, x, fnx);

                    calcRInts(la, lb, p, xyz_pc, fnx, rints_tmp, rints);

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

        transferIntegrals(ipair, sp_data, ints_sph, ints_out);
    }
}