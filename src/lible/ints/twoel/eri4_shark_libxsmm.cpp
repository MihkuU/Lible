#include <lible/ints/twoel/eri4_shark_libxsmm.hpp>

#include <libxsmm_source.h>

namespace LIT = lible::ints::two;

using std::vector;

// void LIT::kernelERI4SharkFlat_libxsmm(const int lalb, const int lcld,
//                                       const int lab, const int lcd,
//                                       const int dim_ab_sph, const int dim_cd_sph,
//                                       const vector<double> &ecoeffs_lalb,
//                                       const vector<double> &ecoeffs_lcld,
//                                       const vector<double> &ecoeffs_lcld_tsp,
//                                       const vector<IdxsTUV> &idxs_tuv_ab,
//                                       const vector<IdxsTUV> &idxs_tuv_cd,
//                                       const ShellPairData &shell_pair_data_ab,
//                                       const ShellPairData &shell_pair_data_cd,
//                                       const BoysF &boys_f, vec4d &eri4)
// {
//     size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
//     size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

//     int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
//     int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

//     libxsmm_mmfunction<double> kernel_R_X_E(LIBXSMM_GEMM_FLAG_NONE, dim_tuv_ab,
//                                             dim_cd_sph, dim_tuv_cd, 1.0, 1.0);
//     libxsmm_mmfunction<double> kernel_E_X_T(LIBXSMM_GEMM_FLAG_NONE, dim_ab_sph,
//                                             dim_cd_sph, dim_tuv_ab, 1.0, 1.0);

//     int labcd = lab + lcd;

//     vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
//     vector<double> fnx(labcd + 1, 0);
//     vec4d rints_tmp(labcd + 1, 0);

//     if (lalb == lcld)
//         for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
//             for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
//             {
//                 vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
//             }
//     else
//         for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
//             for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
//             {
//                 vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
//             }

//     // int labcd = lab + lcd;

//     // const auto &[exps_a, exps_b] = shell_pair_data_ab.exps[ipair_ab];
//     // const auto &[exps_c, exps_d] = shell_pair_data_cd.exps[ipair_cd];
//     // const auto &[xyz_a, xyz_b] = shell_pair_data_ab.coords[ipair_ab];
//     // const auto &[xyz_c, xyz_d] = shell_pair_data_cd.coords[ipair_cd];

//     // const auto &offsets_ecoeffs_ab = shell_pair_data_ab.offsets_ecoeffs;
//     // const auto &offsets_ecoeffs_cd = shell_pair_data_cd.offsets_ecoeffs;

//     // size_t offset_prim_ab = shell_pair_data_ab.offsets_prims[ipair_ab];
//     // size_t offset_prim_cd = shell_pair_data_cd.offsets_prims[ipair_cd];

//     // size_t ka = exps_a.size();
//     // size_t kb = exps_b.size();
//     // size_t kc = exps_c.size();
//     // size_t kd = exps_d.size();

//     // arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
//     // arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};
//     // arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};
//     // arma::vec::fixed<3> D{xyz_d[0], xyz_d[1], xyz_d[2]};

//     // int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
//     // int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

//     // int la = shell_pair_data_ab.la, lb = shell_pair_data_ab.lb;
//     // int lc = shell_pair_data_cd.la, ld = shell_pair_data_cd.lb;
//     // int dim_sph_a = dimSphericals(la);
//     // int dim_sph_b = dimSphericals(lb);
//     // int dim_sph_c = dimSphericals(lc);
//     // int dim_sph_d = dimSphericals(ld);
//     // int dim_sph_ab = dim_sph_a * dim_sph_b;
//     // int dim_sph_cd = dim_sph_c * dim_sph_d;
//     // int dim_E_x_R = dim_tuv_ab * dim_sph_cd;
//     // vector<double> X(ka * kb * dim_E_x_R, 0);

//     // size_t iab = 0;
//     // for (size_t ia = 0; ia < ka; ia++)
//     //     for (size_t ib = 0; ib < kb; ib++)
//     //     {
//     //         size_t pos_X = iab * dim_E_x_R;

//     //         size_t icd = 0;
//     //         for (size_t ic = 0; ic < kc; ic++)
//     //             for (size_t id = 0; id < kd; id++)
//     //             {
//     //                 double a = exps_a[ia];
//     //                 double b = exps_b[ib];
//     //                 double c = exps_c[ic];
//     //                 double d = exps_d[id];

//     //                 double p = a + b;
//     //                 double q = c + d;
//     //                 double alpha = p * q / (p + q);

//     //                 arma::vec::fixed<3> P = (a * A + b * B) / p;
//     //                 arma::vec::fixed<3> Q = (c * C + d * D) / q;
//     //                 arma::vec::fixed<3> RPQ = P - Q;

//     //                 double x = alpha * arma::dot(RPQ, RPQ);

//     //                 boys_f.calcFnx(labcd, x, fnx);

//     //                 double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
//     //                 MD::calcRIntsCM(lab, lcd, alpha, fac, RPQ, fnx, idxs_tuv_ab, idxs_tuv_cd,
//     //                                 rints_tmp, rints);

//     //                 size_t iprim_cd = offset_prim_cd + icd;
//     //                 size_t pos_ecoeffs_cd = offsets_ecoeffs_cd[iprim_cd];

//     //                 kernel_R_X_E(&rints[0], &ecoeffs_lcld_tsp[pos_ecoeffs_cd], &X[pos_X]);
//     //                 icd++;
//     //             }
//     //         iab++;
//     //     }

//     // for (size_t ia = 0, iab = 0; ia < ka; ia++)
//     //     for (size_t ib = 0; ib < kb; ib++, iab++)
//     //     {
//     //         size_t iprim_ab = offset_prim_ab + iab;
//     //         size_t pos_ecoeffs_ab = offsets_ecoeffs_ab[iprim_ab];
//     //         size_t pos_X = iab * dim_E_x_R;

//     //         kernel_E_X_T(&ecoeffs_lalb[pos_ecoeffs_ab], &X[pos_X], &eri4_shells_sph[0]);
//     //     }
// }