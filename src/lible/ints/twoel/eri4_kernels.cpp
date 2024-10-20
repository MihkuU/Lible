#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/rints_meta.hpp>
#include <lible/ints/utils.hpp>

#include <format>

#include <fmt/core.h>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using LIT::kernel_eri4_t;

using std::array, std::string, std::vector;

template <int la, int lb, int lc, int ld>
void LIT::eri4Kernel(const int cdepth_a, const int cdepth_b,
                     const int cdepth_c, const int cdepth_d,
                     const double *exps_a, const double *exps_b,
                     const double *exps_c, const double *exps_d,
                     const double *coords_a, const double *coords_b,
                     const double *coords_c, const double *coords_d,
                     const double *ecoeffs_ab,
                     const double *ecoeffs_cd_tsp,
                     double *eri4_batch)
{
    constexpr int lab = la + lb;
    constexpr int lcd = lc + ld;
    constexpr int labcd = lab + lcd;

    constexpr int n_sph_a = numSphericalsC(la);
    constexpr int n_sph_b = numSphericalsC(lb);
    constexpr int n_sph_c = numSphericalsC(lc);
    constexpr int n_sph_d = numSphericalsC(ld);
    constexpr int n_hermite_ab = numHermitesC(lab);
    constexpr int n_hermite_cd = numHermitesC(lcd);
    constexpr int n_sph_ab = n_sph_a * n_sph_b;
    constexpr int n_sph_cd = n_sph_c * n_sph_d;
    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;

    std::array<double, labcd + 1> fnx;
    BoysF2<labcd> boys_f;

    constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);
    array<double, n_hermites_abcd> rints;

    constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;
    vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double p = a + b;
                    double q = c + d;

                    array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                           (a * coords_a[1] + b * coords_b[1]) / p,
                                           (a * coords_a[2] + b * coords_b[2]) / p};

                    array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,
                                           (c * coords_c[1] + d * coords_d[1]) / q,
                                           (c * coords_c[2] + d * coords_d[2]) / q};

                    array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                            xyz_p[2] - xyz_q[2]};

                    double alpha = p * q / (p + q);
                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                    double x = alpha * (dx * dx + dy * dy + dz * dz);
                    boys_f.calcFnx(x, &fnx[0]);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    calcRInts<lab, lcd>(fac, &fnx[0], p, &xyz_pq[0], &rints[0]);

                    int pos_ecoeffs_cd = icd * n_ecoeffs_cd;

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab,
                                n_sph_cd, n_hermite_cd, 1.0, &rints[0], n_hermite_cd,
                                &ecoeffs_cd_tsp[pos_ecoeffs_cd], n_sph_cd, 1.0,
                                &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd);
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
            int pos_ecoeffs_ab = iab * n_ecoeffs_ab;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd,
                        n_hermite_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], n_hermite_ab,
                        &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd, 1.0,
                        &eri4_batch[0], n_sph_cd);
        }
}

// Table of available kernels
namespace lible::ints::two
{
    const std::array<kernel_eri4_t, 120> eri4_kernels{
        eri4Kernel<0, 0, 0, 0>,
        eri4Kernel<1, 0, 0, 0>,
        eri4Kernel<1, 0, 1, 0>,
        eri4Kernel<1, 1, 0, 0>,
        eri4Kernel<1, 1, 1, 0>,
        eri4Kernel<1, 1, 1, 1>,
        eri4Kernel<2, 0, 0, 0>,
        eri4Kernel<2, 0, 1, 0>,
        eri4Kernel<2, 0, 1, 1>,
        eri4Kernel<2, 0, 2, 0>,
        eri4Kernel<2, 1, 0, 0>,
        eri4Kernel<2, 1, 1, 0>,
        eri4Kernel<2, 1, 1, 1>,
        eri4Kernel<2, 1, 2, 0>,
        eri4Kernel<2, 1, 2, 1>,
        eri4Kernel<2, 2, 0, 0>,
        eri4Kernel<2, 2, 1, 0>,
        eri4Kernel<2, 2, 1, 1>,
        eri4Kernel<2, 2, 2, 0>,
        eri4Kernel<2, 2, 2, 1>,
        eri4Kernel<2, 2, 2, 2>,
        eri4Kernel<3, 0, 0, 0>,
        eri4Kernel<3, 0, 1, 0>,
        eri4Kernel<3, 0, 1, 1>,
        eri4Kernel<3, 0, 2, 0>,
        eri4Kernel<3, 0, 2, 1>,
        eri4Kernel<3, 0, 2, 2>,
        eri4Kernel<3, 0, 3, 0>,
        eri4Kernel<3, 1, 0, 0>,
        eri4Kernel<3, 1, 1, 0>,
        eri4Kernel<3, 1, 1, 1>,
        eri4Kernel<3, 1, 2, 0>,
        eri4Kernel<3, 1, 2, 1>,
        eri4Kernel<3, 1, 2, 2>,
        eri4Kernel<3, 1, 3, 0>,
        eri4Kernel<3, 1, 3, 1>,
        eri4Kernel<3, 2, 0, 0>,
        eri4Kernel<3, 2, 1, 0>,
        eri4Kernel<3, 2, 1, 1>,
        eri4Kernel<3, 2, 2, 0>,
        eri4Kernel<3, 2, 2, 1>,
        eri4Kernel<3, 2, 2, 2>,
        eri4Kernel<3, 2, 3, 0>,
        eri4Kernel<3, 2, 3, 1>,
        eri4Kernel<3, 2, 3, 2>,
        eri4Kernel<3, 3, 0, 0>,
        eri4Kernel<3, 3, 1, 0>,
        eri4Kernel<3, 3, 1, 1>,
        eri4Kernel<3, 3, 2, 0>,
        eri4Kernel<3, 3, 2, 1>,
        eri4Kernel<3, 3, 2, 2>,
        eri4Kernel<3, 3, 3, 0>,
        eri4Kernel<3, 3, 3, 1>,
        eri4Kernel<3, 3, 3, 2>,
        eri4Kernel<3, 3, 3, 3>,
        eri4Kernel<4, 0, 0, 0>,
        eri4Kernel<4, 0, 1, 0>,
        eri4Kernel<4, 0, 1, 1>,
        eri4Kernel<4, 0, 2, 0>,
        eri4Kernel<4, 0, 2, 1>,
        eri4Kernel<4, 0, 2, 2>,
        eri4Kernel<4, 0, 3, 0>,
        eri4Kernel<4, 0, 3, 1>,
        eri4Kernel<4, 0, 3, 2>,
        eri4Kernel<4, 0, 3, 3>,
        eri4Kernel<4, 0, 4, 0>,
        eri4Kernel<4, 1, 0, 0>,
        eri4Kernel<4, 1, 1, 0>,
        eri4Kernel<4, 1, 1, 1>,
        eri4Kernel<4, 1, 2, 0>,
        eri4Kernel<4, 1, 2, 1>,
        eri4Kernel<4, 1, 2, 2>,
        eri4Kernel<4, 1, 3, 0>,
        eri4Kernel<4, 1, 3, 1>,
        eri4Kernel<4, 1, 3, 2>,
        eri4Kernel<4, 1, 3, 3>,
        eri4Kernel<4, 1, 4, 0>,
        eri4Kernel<4, 1, 4, 1>,
        eri4Kernel<4, 2, 0, 0>,
        eri4Kernel<4, 2, 1, 0>,
        eri4Kernel<4, 2, 1, 1>,
        eri4Kernel<4, 2, 2, 0>,
        eri4Kernel<4, 2, 2, 1>,
        eri4Kernel<4, 2, 2, 2>,
        eri4Kernel<4, 2, 3, 0>,
        eri4Kernel<4, 2, 3, 1>,
        eri4Kernel<4, 2, 3, 2>,
        eri4Kernel<4, 2, 3, 3>,
        eri4Kernel<4, 2, 4, 0>,
        eri4Kernel<4, 2, 4, 1>,
        eri4Kernel<4, 2, 4, 2>,
        eri4Kernel<4, 3, 0, 0>,
        eri4Kernel<4, 3, 1, 0>,
        eri4Kernel<4, 3, 1, 1>,
        eri4Kernel<4, 3, 2, 0>,
        eri4Kernel<4, 3, 2, 1>,
        eri4Kernel<4, 3, 2, 2>,
        eri4Kernel<4, 3, 3, 0>,
        eri4Kernel<4, 3, 3, 1>,
        eri4Kernel<4, 3, 3, 2>,
        eri4Kernel<4, 3, 3, 3>,
        eri4Kernel<4, 3, 4, 0>,
        eri4Kernel<4, 3, 4, 1>,
        eri4Kernel<4, 3, 4, 2>,
        eri4Kernel<4, 3, 4, 3>,
        eri4Kernel<4, 4, 0, 0>,
        eri4Kernel<4, 4, 1, 0>,
        eri4Kernel<4, 4, 1, 1>,
        eri4Kernel<4, 4, 2, 0>,
        eri4Kernel<4, 4, 2, 1>,
        eri4Kernel<4, 4, 2, 2>,
        eri4Kernel<4, 4, 3, 0>,
        eri4Kernel<4, 4, 3, 1>,
        eri4Kernel<4, 4, 3, 2>,
        eri4Kernel<4, 4, 3, 3>,
        eri4Kernel<4, 4, 4, 0>,
        eri4Kernel<4, 4, 4, 1>,
        eri4Kernel<4, 4, 4, 2>,
        eri4Kernel<4, 4, 4, 3>,
        eri4Kernel<4, 4, 4, 4>
    };
}

kernel_eri4_t LIT::deployERI4Kernel(const int la, const int lb, const int lc, const int ld)
{
    int idx_ab = la * (la + 1) / 2 + lb;
    int idx_cd = lc * (lc + 1) / 2 + ld;

    if (idx_ab < idx_cd)
        throw std::runtime_error(std::format("(la, lb) >= (lc, ld) condition must be satisfied, given was: {} vs {}!\n",
                                             idx_ab, idx_cd));

    int labcd = la + lb + lc + ld;
    int l_max = _eri_kernel_max_l_;
    if (labcd > _eri_kernel_max_l_)
        throw std::runtime_error(std::format("lab + lcd = {} is larger than the allowed max: {}!\n",
                                             labcd, l_max));

    int idx_abcd = idx_ab * (idx_ab + 1) / 2 + idx_cd;

    return eri4_kernels.at(idx_abcd);
}