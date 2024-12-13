#include <lible/ints/twoel/eri_kernels.hpp>

template<> void lible::ints::two::eri4Kernel<1, 1, 0, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 1, lb = 1, lc = 0, ld = 0;
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

    std::fill(eri4_batch, eri4_batch + n_sph_ab * n_sph_cd, 0);

    std::array<double, labcd + 1> fnx;
    BoysF2<labcd> boys_f;

    constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);
    std::array<double, n_hermites_abcd> rints;

    constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;
    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

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
                    std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                                (a * coords_a[1] + b * coords_b[1]) / p,
                                                (a * coords_a[2] + b * coords_b[2]) / p};

                    std::array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,
                                                (c * coords_c[1] + d * coords_d[1]) / q,
                                                (c * coords_c[2] + d * coords_d[2]) / q};

                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                 xyz_p[2] - xyz_q[2]};

                    double alpha = p * q / (p + q);
                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                    double x = alpha * (dx * dx + dy * dy + dz * dz);
                    boys_f.calcFnx(x, &fnx[0]);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    calcRInts<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                    const double* p_ecoeffs_cd_tsp = &ecoeffs_cd_tsp[icd * n_ecoeffs_cd];
                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_cd_tsp[0];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch[1] += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[1];
            eri4_batch[1] += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[3];
            eri4_batch[1] += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[6];
            eri4_batch[2] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch[2] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri4_batch[2] += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[8];
            eri4_batch[3] += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[1];
            eri4_batch[3] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri4_batch[3] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[3];
            eri4_batch[3] += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[6];
            eri4_batch[4] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri4_batch[4] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch[4] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri4_batch[5] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[1];
            eri4_batch[5] += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[0];
            eri4_batch[5] += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[5];
            eri4_batch[5] += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[2];
            eri4_batch[6] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri4_batch[6] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri4_batch[6] += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[2];
            eri4_batch[6] += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[8];
            eri4_batch[7] += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri4_batch[7] += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri4_batch[7] += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri4_batch[7] += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri4_batch[8] += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[7];
            eri4_batch[8] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri4_batch[8] += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[2];
        }
}

template<> void lible::ints::two::eri4Kernel<2, 0, 0, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 2, lb = 0, lc = 0, ld = 0;
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

    std::fill(eri4_batch, eri4_batch + n_sph_ab * n_sph_cd, 0);

    std::array<double, labcd + 1> fnx;
    BoysF2<labcd> boys_f;

    constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);
    std::array<double, n_hermites_abcd> rints;

    constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;
    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

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
                    std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                                (a * coords_a[1] + b * coords_b[1]) / p,
                                                (a * coords_a[2] + b * coords_b[2]) / p};

                    std::array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,
                                                (c * coords_c[1] + d * coords_d[1]) / q,
                                                (c * coords_c[2] + d * coords_d[2]) / q};

                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                 xyz_p[2] - xyz_q[2]};

                    double alpha = p * q / (p + q);
                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                    double x = alpha * (dx * dx + dy * dy + dz * dz);
                    boys_f.calcFnx(x, &fnx[0]);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    calcRInts<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                    const double* p_ecoeffs_cd_tsp = &ecoeffs_cd_tsp[icd * n_ecoeffs_cd];
                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_cd_tsp[0];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch[1] += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[1];
            eri4_batch[1] += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[3];
            eri4_batch[1] += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[6];
            eri4_batch[2] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch[2] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri4_batch[2] += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[8];
            eri4_batch[3] += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[2];
            eri4_batch[3] += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[7];
            eri4_batch[3] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri4_batch[3] += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[1];
            eri4_batch[3] += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[4];
            eri4_batch[4] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri4_batch[4] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch[4] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri4_batch[4] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
        }
}

template<> void lible::ints::two::eri3Kernel<1, 1, 0>(const int cdepth_a, const int cdepth_b, const int cdepth_c,
                                                      const double* exps_a, const double* exps_b,
                                                      const double* exps_c, const double* coords_a,
                                                      const double* coords_b, const double* coords_c,
                                                      const double* ecoeffs_ab, const double* ecoeffs_c,
                                                      double* eri3_batch)
{
    constexpr int la = 1, lb = 1, lc = 0;
    constexpr int lab = la + lb;
    constexpr int labc = lab + lc;
    constexpr int n_sph_a = numSphericalsC(la);
    constexpr int n_sph_b = numSphericalsC(lb);
    constexpr int n_sph_c = numSphericalsC(lc);
    constexpr int n_hermite_ab = numHermitesC(lab);
    constexpr int n_hermite_c = numHermitesC(lc);
    constexpr int n_sph_ab = n_sph_a * n_sph_b;
    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    constexpr int n_ecoeffs_c = n_sph_c * n_hermite_c;

    std::fill(eri3_batch, eri3_batch + n_sph_ab * n_sph_c, 0);

    std::array<double, labc + 1> fnx;
    BoysF2<labc> boys_f;

    constexpr int n_hermites_abc = numHermitesC(lab) * numHermitesC(lc);
    std::array<double, n_hermites_abc> rints;

    constexpr int n_rints_x_ecoeffs = n_sph_c * n_hermite_ab;
    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double c = exps_c[ic];

                double p = a + b;
                double alpha = p * c / (p + c);

                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                            (a * coords_a[2] + b * coords_b[2]) / p};

                std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                             xyz_p[1] - coords_c[1],
                                             xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                boys_f.calcFnx(x, &fnx[0]);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                calcRInts<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                    const double* p_ecoeffs_c = &ecoeffs_c[ic * n_ecoeffs_c];
                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_c[0];
            }
}

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch[1] += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[1];
            eri3_batch[1] += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri3_batch[1] += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[3];
            eri3_batch[1] += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[6];
            eri3_batch[2] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch[2] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri3_batch[2] += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[8];
            eri3_batch[3] += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[1];
            eri3_batch[3] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri3_batch[3] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[3];
            eri3_batch[3] += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[6];
            eri3_batch[4] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri3_batch[4] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch[4] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri3_batch[5] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[1];
            eri3_batch[5] += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[0];
            eri3_batch[5] += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[5];
            eri3_batch[5] += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[2];
            eri3_batch[6] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri3_batch[6] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri3_batch[6] += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[2];
            eri3_batch[6] += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[8];
            eri3_batch[7] += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri3_batch[7] += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri3_batch[7] += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri3_batch[7] += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri3_batch[8] += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[7];
            eri3_batch[8] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri3_batch[8] += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[2];
        }
}

template<> void lible::ints::two::eri3Kernel<2, 0, 0>(const int cdepth_a, const int cdepth_b, const int cdepth_c,
                                                      const double* exps_a, const double* exps_b,
                                                      const double* exps_c, const double* coords_a,
                                                      const double* coords_b, const double* coords_c,
                                                      const double* ecoeffs_ab, const double* ecoeffs_c,
                                                      double* eri3_batch)
{
    constexpr int la = 2, lb = 0, lc = 0;
    constexpr int lab = la + lb;
    constexpr int labc = lab + lc;
    constexpr int n_sph_a = numSphericalsC(la);
    constexpr int n_sph_b = numSphericalsC(lb);
    constexpr int n_sph_c = numSphericalsC(lc);
    constexpr int n_hermite_ab = numHermitesC(lab);
    constexpr int n_hermite_c = numHermitesC(lc);
    constexpr int n_sph_ab = n_sph_a * n_sph_b;
    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    constexpr int n_ecoeffs_c = n_sph_c * n_hermite_c;

    std::fill(eri3_batch, eri3_batch + n_sph_ab * n_sph_c, 0);

    std::array<double, labc + 1> fnx;
    BoysF2<labc> boys_f;

    constexpr int n_hermites_abc = numHermitesC(lab) * numHermitesC(lc);
    std::array<double, n_hermites_abc> rints;

    constexpr int n_rints_x_ecoeffs = n_sph_c * n_hermite_ab;
    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double c = exps_c[ic];

                double p = a + b;
                double alpha = p * c / (p + c);

                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                            (a * coords_a[2] + b * coords_b[2]) / p};

                std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                             xyz_p[1] - coords_c[1],
                                             xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                boys_f.calcFnx(x, &fnx[0]);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                calcRInts<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                    const double* p_ecoeffs_c = &ecoeffs_c[ic * n_ecoeffs_c];
                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_c[0];
            }
}

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch[1] += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[1];
            eri3_batch[1] += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri3_batch[1] += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[3];
            eri3_batch[1] += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[6];
            eri3_batch[2] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch[2] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri3_batch[2] += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[8];
            eri3_batch[3] += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[2];
            eri3_batch[3] += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[7];
            eri3_batch[3] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri3_batch[3] += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[1];
            eri3_batch[3] += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[4];
            eri3_batch[4] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri3_batch[4] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch[4] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri3_batch[4] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
        }
}

template<> void lible::ints::two::eri2Kernel<2, 0>(const int cdepth_a, const int cdepth_b,
                                                   const double* exps_a, const double* exps_b,
                                                   const double* coords_a, const double* coords_b,
                                                   const double* ecoeffs_a, const double* ecoeffs_b_tsp,
                                                   double* eri2_batch)
{
    constexpr int la = 2, lb = 0;
    constexpr int lab = la + lb;
    constexpr int n_sph_a = numSphericalsC(la);
    constexpr int n_sph_b = numSphericalsC(lb);
    constexpr int n_hermite_a = numHermitesC(la);
    constexpr int n_hermite_b = numHermitesC(lb);
    constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
    constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

    std::fill(eri2_batch, eri2_batch + n_sph_a * n_sph_b, 0);

    std::array<double, lab + 1> fnx;
    BoysF2<lab> boys_f;

    constexpr int n_hermites_ab = numHermitesC(la) * numHermitesC(lb);
    std::array<double, n_hermites_ab> rints;

    constexpr int n_rints_x_ecoeffs = n_hermite_a * n_sph_b;
    std::vector<double> rints_x_ecoeffs(cdepth_a * n_rints_x_ecoeffs, 0);

    std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
    coords_a[1] - coords_b[1],
    coords_a[2] - coords_b[2]};

    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;	

    for (int ia = 0; ia < cdepth_a; ia++)
    {
        int pos_rints_x_ecoeffs = ia * n_rints_x_ecoeffs;
        for (int ib = 0; ib < cdepth_b; ib++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            double alpha = a * b / (a + b);
            double x = alpha * xyz_ab_dot;
            boys_f.calcFnx(x, &fnx[0]);

            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
            calcRInts<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

                    const double* p_ecoeffs_b = &ecoeffs_b_tsp[ib * n_ecoeffs_b];
                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_b[0];
        }
    }

    for (int ia = 0; ia < cdepth_a; ia++)
    {
            const double* p_ecoeffs_a = &ecoeffs_a[ia * n_ecoeffs_a];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];

            eri2_batch[0] += p_ecoeffs_a[2] * p_rints_x_ecoeffs[2];
            eri2_batch[0] += p_ecoeffs_a[7] * p_rints_x_ecoeffs[7];
            eri2_batch[0] += p_ecoeffs_a[0] * p_rints_x_ecoeffs[0];
            eri2_batch[0] += p_ecoeffs_a[1] * p_rints_x_ecoeffs[1];
            eri2_batch[0] += p_ecoeffs_a[9] * p_rints_x_ecoeffs[9];
            eri2_batch[0] += p_ecoeffs_a[4] * p_rints_x_ecoeffs[4];
            eri2_batch[0] += p_ecoeffs_a[3] * p_rints_x_ecoeffs[3];
            eri2_batch[1] += p_ecoeffs_a[11] * p_rints_x_ecoeffs[1];
            eri2_batch[1] += p_ecoeffs_a[10] * p_rints_x_ecoeffs[0];
            eri2_batch[1] += p_ecoeffs_a[13] * p_rints_x_ecoeffs[3];
            eri2_batch[1] += p_ecoeffs_a[16] * p_rints_x_ecoeffs[6];
            eri2_batch[2] += p_ecoeffs_a[20] * p_rints_x_ecoeffs[0];
            eri2_batch[2] += p_ecoeffs_a[23] * p_rints_x_ecoeffs[3];
            eri2_batch[2] += p_ecoeffs_a[22] * p_rints_x_ecoeffs[2];
            eri2_batch[2] += p_ecoeffs_a[28] * p_rints_x_ecoeffs[8];
            eri2_batch[3] += p_ecoeffs_a[32] * p_rints_x_ecoeffs[2];
            eri2_batch[3] += p_ecoeffs_a[37] * p_rints_x_ecoeffs[7];
            eri2_batch[3] += p_ecoeffs_a[30] * p_rints_x_ecoeffs[0];
            eri2_batch[3] += p_ecoeffs_a[31] * p_rints_x_ecoeffs[1];
            eri2_batch[3] += p_ecoeffs_a[34] * p_rints_x_ecoeffs[4];
            eri2_batch[4] += p_ecoeffs_a[41] * p_rints_x_ecoeffs[1];
            eri2_batch[4] += p_ecoeffs_a[40] * p_rints_x_ecoeffs[0];
            eri2_batch[4] += p_ecoeffs_a[45] * p_rints_x_ecoeffs[5];
            eri2_batch[4] += p_ecoeffs_a[42] * p_rints_x_ecoeffs[2];
    }
}
