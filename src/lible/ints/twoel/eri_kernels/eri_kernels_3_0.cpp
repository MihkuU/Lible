#include <lible/ints/twoel/eri_kernels.hpp>

template<> void lible::ints::two::eri4Kernel<2, 1, 0, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 2, lb = 1, lc = 0, ld = 0;
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
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[12] += rints[12] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[13] += rints[13] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[16] += rints[16] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[18] += rints[18] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[19] += rints[19] * p_ecoeffs_cd_tsp[0];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri4_batch[0] += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri4_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch[0] += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri4_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch[0] += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri4_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch[0] += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri4_batch[1] += p_ecoeffs_ab[26] * p_rints_x_ecoeffs[6];
            eri4_batch[1] += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[5];
            eri4_batch[1] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri4_batch[1] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[13];
            eri4_batch[1] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[10];
            eri4_batch[1] += p_ecoeffs_ab[27] * p_rints_x_ecoeffs[7];
            eri4_batch[1] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[21] * p_rints_x_ecoeffs[1];
            eri4_batch[1] += p_ecoeffs_ab[29] * p_rints_x_ecoeffs[9];
            eri4_batch[1] += p_ecoeffs_ab[24] * p_rints_x_ecoeffs[4];
            eri4_batch[1] += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[15];
            eri4_batch[1] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri4_batch[2] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
            eri4_batch[2] += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[7];
            eri4_batch[2] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[11];
            eri4_batch[2] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch[2] += p_ecoeffs_ab[56] * p_rints_x_ecoeffs[16];
            eri4_batch[2] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri4_batch[2] += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[9];
            eri4_batch[2] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri4_batch[2] += p_ecoeffs_ab[58] * p_rints_x_ecoeffs[18];
            eri4_batch[2] += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[8];
            eri4_batch[3] += p_ecoeffs_ab[66] * p_rints_x_ecoeffs[6];
            eri4_batch[3] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri4_batch[3] += p_ecoeffs_ab[61] * p_rints_x_ecoeffs[1];
            eri4_batch[3] += p_ecoeffs_ab[69] * p_rints_x_ecoeffs[9];
            eri4_batch[3] += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[15];
            eri4_batch[3] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri4_batch[4] += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[6];
            eri4_batch[4] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri4_batch[4] += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[1];
            eri4_batch[4] += p_ecoeffs_ab[92] * p_rints_x_ecoeffs[12];
            eri4_batch[4] += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[4];
            eri4_batch[4] += p_ecoeffs_ab[83] * p_rints_x_ecoeffs[3];
            eri4_batch[5] += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[6];
            eri4_batch[5] += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[5];
            eri4_batch[5] += p_ecoeffs_ab[102] * p_rints_x_ecoeffs[2];
            eri4_batch[5] += p_ecoeffs_ab[100] * p_rints_x_ecoeffs[0];
            eri4_batch[5] += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[1];
            eri4_batch[5] += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[3];
            eri4_batch[5] += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[14];
            eri4_batch[5] += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[8];
            eri4_batch[6] += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[2];
            eri4_batch[6] += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[0];
            eri4_batch[6] += p_ecoeffs_ab[129] * p_rints_x_ecoeffs[9];
            eri4_batch[6] += p_ecoeffs_ab[138] * p_rints_x_ecoeffs[18];
            eri4_batch[6] += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[3];
            eri4_batch[6] += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[8];
            eri4_batch[7] += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri4_batch[7] += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri4_batch[7] += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri4_batch[7] += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri4_batch[7] += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri4_batch[7] += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri4_batch[7] += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri4_batch[7] += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri4_batch[8] += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[17];
            eri4_batch[8] += p_ecoeffs_ab[162] * p_rints_x_ecoeffs[2];
            eri4_batch[8] += p_ecoeffs_ab[167] * p_rints_x_ecoeffs[7];
            eri4_batch[8] += p_ecoeffs_ab[160] * p_rints_x_ecoeffs[0];
            eri4_batch[8] += p_ecoeffs_ab[163] * p_rints_x_ecoeffs[3];
            eri4_batch[8] += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[8];
            eri4_batch[9] += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[17];
            eri4_batch[9] += p_ecoeffs_ab[186] * p_rints_x_ecoeffs[6];
            eri4_batch[9] += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[2];
            eri4_batch[9] += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[7];
            eri4_batch[9] += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[0];
            eri4_batch[9] += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[1];
            eri4_batch[9] += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[12];
            eri4_batch[9] += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[4];
            eri4_batch[9] += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[3];
            eri4_batch[9] += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[8];
            eri4_batch[10] += p_ecoeffs_ab[205] * p_rints_x_ecoeffs[5];
            eri4_batch[10] += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[2];
            eri4_batch[10] += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[13];
            eri4_batch[10] += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[10];
            eri4_batch[10] += p_ecoeffs_ab[207] * p_rints_x_ecoeffs[7];
            eri4_batch[10] += p_ecoeffs_ab[200] * p_rints_x_ecoeffs[0];
            eri4_batch[10] += p_ecoeffs_ab[201] * p_rints_x_ecoeffs[1];
            eri4_batch[10] += p_ecoeffs_ab[204] * p_rints_x_ecoeffs[4];
            eri4_batch[11] += p_ecoeffs_ab[225] * p_rints_x_ecoeffs[5];
            eri4_batch[11] += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[2];
            eri4_batch[11] += p_ecoeffs_ab[231] * p_rints_x_ecoeffs[11];
            eri4_batch[11] += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[7];
            eri4_batch[11] += p_ecoeffs_ab[220] * p_rints_x_ecoeffs[0];
            eri4_batch[11] += p_ecoeffs_ab[236] * p_rints_x_ecoeffs[16];
            eri4_batch[11] += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[1];
            eri4_batch[11] += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[4];
            eri4_batch[12] += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[6];
            eri4_batch[12] += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[5];
            eri4_batch[12] += p_ecoeffs_ab[242] * p_rints_x_ecoeffs[2];
            eri4_batch[12] += p_ecoeffs_ab[240] * p_rints_x_ecoeffs[0];
            eri4_batch[12] += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[1];
            eri4_batch[12] += p_ecoeffs_ab[243] * p_rints_x_ecoeffs[3];
            eri4_batch[12] += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[14];
            eri4_batch[12] += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[8];
            eri4_batch[13] += p_ecoeffs_ab[265] * p_rints_x_ecoeffs[5];
            eri4_batch[13] += p_ecoeffs_ab[262] * p_rints_x_ecoeffs[2];
            eri4_batch[13] += p_ecoeffs_ab[271] * p_rints_x_ecoeffs[11];
            eri4_batch[13] += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[0];
            eri4_batch[13] += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[1];
            eri4_batch[13] += p_ecoeffs_ab[264] * p_rints_x_ecoeffs[4];
            eri4_batch[14] += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri4_batch[14] += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri4_batch[14] += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri4_batch[14] += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri4_batch[14] += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri4_batch[14] += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
        }
}

template<> void lible::ints::two::eri4Kernel<3, 0, 0, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 3, lb = 0, lc = 0, ld = 0;
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
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[12] += rints[12] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[13] += rints[13] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[16] += rints[16] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[18] += rints[18] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[19] += rints[19] * p_ecoeffs_cd_tsp[0];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri4_batch[0] += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri4_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch[0] += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri4_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch[0] += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri4_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch[0] += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri4_batch[1] += p_ecoeffs_ab[26] * p_rints_x_ecoeffs[6];
            eri4_batch[1] += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[5];
            eri4_batch[1] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri4_batch[1] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[13];
            eri4_batch[1] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[10];
            eri4_batch[1] += p_ecoeffs_ab[27] * p_rints_x_ecoeffs[7];
            eri4_batch[1] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[21] * p_rints_x_ecoeffs[1];
            eri4_batch[1] += p_ecoeffs_ab[29] * p_rints_x_ecoeffs[9];
            eri4_batch[1] += p_ecoeffs_ab[24] * p_rints_x_ecoeffs[4];
            eri4_batch[1] += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[15];
            eri4_batch[1] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri4_batch[2] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
            eri4_batch[2] += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[7];
            eri4_batch[2] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[11];
            eri4_batch[2] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch[2] += p_ecoeffs_ab[56] * p_rints_x_ecoeffs[16];
            eri4_batch[2] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri4_batch[2] += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[9];
            eri4_batch[2] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri4_batch[2] += p_ecoeffs_ab[58] * p_rints_x_ecoeffs[18];
            eri4_batch[2] += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[3];
            eri4_batch[2] += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[8];
            eri4_batch[3] += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[17];
            eri4_batch[3] += p_ecoeffs_ab[66] * p_rints_x_ecoeffs[6];
            eri4_batch[3] += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[2];
            eri4_batch[3] += p_ecoeffs_ab[67] * p_rints_x_ecoeffs[7];
            eri4_batch[3] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri4_batch[3] += p_ecoeffs_ab[61] * p_rints_x_ecoeffs[1];
            eri4_batch[3] += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[12];
            eri4_batch[3] += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[4];
            eri4_batch[3] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri4_batch[3] += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[8];
            eri4_batch[4] += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[6];
            eri4_batch[4] += p_ecoeffs_ab[85] * p_rints_x_ecoeffs[5];
            eri4_batch[4] += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[2];
            eri4_batch[4] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri4_batch[4] += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[1];
            eri4_batch[4] += p_ecoeffs_ab[83] * p_rints_x_ecoeffs[3];
            eri4_batch[4] += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[14];
            eri4_batch[4] += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[8];
            eri4_batch[5] += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[5];
            eri4_batch[5] += p_ecoeffs_ab[102] * p_rints_x_ecoeffs[2];
            eri4_batch[5] += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[13];
            eri4_batch[5] += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[10];
            eri4_batch[5] += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[7];
            eri4_batch[5] += p_ecoeffs_ab[100] * p_rints_x_ecoeffs[0];
            eri4_batch[5] += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[1];
            eri4_batch[5] += p_ecoeffs_ab[104] * p_rints_x_ecoeffs[4];
            eri4_batch[6] += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[5];
            eri4_batch[6] += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[2];
            eri4_batch[6] += p_ecoeffs_ab[131] * p_rints_x_ecoeffs[11];
            eri4_batch[6] += p_ecoeffs_ab[127] * p_rints_x_ecoeffs[7];
            eri4_batch[6] += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[0];
            eri4_batch[6] += p_ecoeffs_ab[136] * p_rints_x_ecoeffs[16];
            eri4_batch[6] += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[1];
            eri4_batch[6] += p_ecoeffs_ab[124] * p_rints_x_ecoeffs[4];
        }
}

template<> void lible::ints::two::eri3Kernel<2, 1, 0>(const int cdepth_a, const int cdepth_b, const int cdepth_c,
                                                      const double* exps_a, const double* exps_b,
                                                      const double* exps_c, const double* coords_a,
                                                      const double* coords_b, const double* coords_c,
                                                      const double* ecoeffs_ab, const double* ecoeffs_c,
                                                      double* eri3_batch)
{
    constexpr int la = 2, lb = 1, lc = 0;
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
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[12] += rints[12] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[13] += rints[13] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[16] += rints[16] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[18] += rints[18] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[19] += rints[19] * p_ecoeffs_c[0];
            }
}

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch[0] += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri3_batch[0] += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri3_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch[0] += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri3_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch[0] += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri3_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch[0] += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri3_batch[1] += p_ecoeffs_ab[26] * p_rints_x_ecoeffs[6];
            eri3_batch[1] += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[5];
            eri3_batch[1] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri3_batch[1] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[13];
            eri3_batch[1] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[10];
            eri3_batch[1] += p_ecoeffs_ab[27] * p_rints_x_ecoeffs[7];
            eri3_batch[1] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch[1] += p_ecoeffs_ab[21] * p_rints_x_ecoeffs[1];
            eri3_batch[1] += p_ecoeffs_ab[29] * p_rints_x_ecoeffs[9];
            eri3_batch[1] += p_ecoeffs_ab[24] * p_rints_x_ecoeffs[4];
            eri3_batch[1] += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[15];
            eri3_batch[1] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri3_batch[2] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
            eri3_batch[2] += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[7];
            eri3_batch[2] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[11];
            eri3_batch[2] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch[2] += p_ecoeffs_ab[56] * p_rints_x_ecoeffs[16];
            eri3_batch[2] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri3_batch[2] += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[9];
            eri3_batch[2] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri3_batch[2] += p_ecoeffs_ab[58] * p_rints_x_ecoeffs[18];
            eri3_batch[2] += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[8];
            eri3_batch[3] += p_ecoeffs_ab[66] * p_rints_x_ecoeffs[6];
            eri3_batch[3] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri3_batch[3] += p_ecoeffs_ab[61] * p_rints_x_ecoeffs[1];
            eri3_batch[3] += p_ecoeffs_ab[69] * p_rints_x_ecoeffs[9];
            eri3_batch[3] += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[15];
            eri3_batch[3] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri3_batch[4] += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[6];
            eri3_batch[4] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri3_batch[4] += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[1];
            eri3_batch[4] += p_ecoeffs_ab[92] * p_rints_x_ecoeffs[12];
            eri3_batch[4] += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[4];
            eri3_batch[4] += p_ecoeffs_ab[83] * p_rints_x_ecoeffs[3];
            eri3_batch[5] += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[6];
            eri3_batch[5] += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[5];
            eri3_batch[5] += p_ecoeffs_ab[102] * p_rints_x_ecoeffs[2];
            eri3_batch[5] += p_ecoeffs_ab[100] * p_rints_x_ecoeffs[0];
            eri3_batch[5] += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[1];
            eri3_batch[5] += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[3];
            eri3_batch[5] += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[14];
            eri3_batch[5] += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[8];
            eri3_batch[6] += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[2];
            eri3_batch[6] += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[0];
            eri3_batch[6] += p_ecoeffs_ab[129] * p_rints_x_ecoeffs[9];
            eri3_batch[6] += p_ecoeffs_ab[138] * p_rints_x_ecoeffs[18];
            eri3_batch[6] += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[3];
            eri3_batch[6] += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[8];
            eri3_batch[7] += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri3_batch[7] += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri3_batch[7] += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri3_batch[7] += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri3_batch[7] += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri3_batch[7] += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri3_batch[7] += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri3_batch[7] += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri3_batch[8] += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[17];
            eri3_batch[8] += p_ecoeffs_ab[162] * p_rints_x_ecoeffs[2];
            eri3_batch[8] += p_ecoeffs_ab[167] * p_rints_x_ecoeffs[7];
            eri3_batch[8] += p_ecoeffs_ab[160] * p_rints_x_ecoeffs[0];
            eri3_batch[8] += p_ecoeffs_ab[163] * p_rints_x_ecoeffs[3];
            eri3_batch[8] += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[8];
            eri3_batch[9] += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[17];
            eri3_batch[9] += p_ecoeffs_ab[186] * p_rints_x_ecoeffs[6];
            eri3_batch[9] += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[2];
            eri3_batch[9] += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[7];
            eri3_batch[9] += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[0];
            eri3_batch[9] += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[1];
            eri3_batch[9] += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[12];
            eri3_batch[9] += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[4];
            eri3_batch[9] += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[3];
            eri3_batch[9] += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[8];
            eri3_batch[10] += p_ecoeffs_ab[205] * p_rints_x_ecoeffs[5];
            eri3_batch[10] += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[2];
            eri3_batch[10] += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[13];
            eri3_batch[10] += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[10];
            eri3_batch[10] += p_ecoeffs_ab[207] * p_rints_x_ecoeffs[7];
            eri3_batch[10] += p_ecoeffs_ab[200] * p_rints_x_ecoeffs[0];
            eri3_batch[10] += p_ecoeffs_ab[201] * p_rints_x_ecoeffs[1];
            eri3_batch[10] += p_ecoeffs_ab[204] * p_rints_x_ecoeffs[4];
            eri3_batch[11] += p_ecoeffs_ab[225] * p_rints_x_ecoeffs[5];
            eri3_batch[11] += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[2];
            eri3_batch[11] += p_ecoeffs_ab[231] * p_rints_x_ecoeffs[11];
            eri3_batch[11] += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[7];
            eri3_batch[11] += p_ecoeffs_ab[220] * p_rints_x_ecoeffs[0];
            eri3_batch[11] += p_ecoeffs_ab[236] * p_rints_x_ecoeffs[16];
            eri3_batch[11] += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[1];
            eri3_batch[11] += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[4];
            eri3_batch[12] += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[6];
            eri3_batch[12] += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[5];
            eri3_batch[12] += p_ecoeffs_ab[242] * p_rints_x_ecoeffs[2];
            eri3_batch[12] += p_ecoeffs_ab[240] * p_rints_x_ecoeffs[0];
            eri3_batch[12] += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[1];
            eri3_batch[12] += p_ecoeffs_ab[243] * p_rints_x_ecoeffs[3];
            eri3_batch[12] += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[14];
            eri3_batch[12] += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[8];
            eri3_batch[13] += p_ecoeffs_ab[265] * p_rints_x_ecoeffs[5];
            eri3_batch[13] += p_ecoeffs_ab[262] * p_rints_x_ecoeffs[2];
            eri3_batch[13] += p_ecoeffs_ab[271] * p_rints_x_ecoeffs[11];
            eri3_batch[13] += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[0];
            eri3_batch[13] += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[1];
            eri3_batch[13] += p_ecoeffs_ab[264] * p_rints_x_ecoeffs[4];
            eri3_batch[14] += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri3_batch[14] += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri3_batch[14] += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri3_batch[14] += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri3_batch[14] += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri3_batch[14] += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
        }
}

template<> void lible::ints::two::eri3Kernel<3, 0, 0>(const int cdepth_a, const int cdepth_b, const int cdepth_c,
                                                      const double* exps_a, const double* exps_b,
                                                      const double* exps_c, const double* coords_a,
                                                      const double* coords_b, const double* coords_c,
                                                      const double* ecoeffs_ab, const double* ecoeffs_c,
                                                      double* eri3_batch)
{
    constexpr int la = 3, lb = 0, lc = 0;
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
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[12] += rints[12] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[13] += rints[13] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[16] += rints[16] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[18] += rints[18] * p_ecoeffs_c[0];
                    p_rints_x_ecoeffs[19] += rints[19] * p_ecoeffs_c[0];
            }
}

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch[0] += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri3_batch[0] += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri3_batch[0] += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch[0] += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri3_batch[0] += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch[0] += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch[0] += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri3_batch[0] += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch[0] += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch[0] += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch[0] += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri3_batch[1] += p_ecoeffs_ab[26] * p_rints_x_ecoeffs[6];
            eri3_batch[1] += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[5];
            eri3_batch[1] += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[2];
            eri3_batch[1] += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[13];
            eri3_batch[1] += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[10];
            eri3_batch[1] += p_ecoeffs_ab[27] * p_rints_x_ecoeffs[7];
            eri3_batch[1] += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch[1] += p_ecoeffs_ab[21] * p_rints_x_ecoeffs[1];
            eri3_batch[1] += p_ecoeffs_ab[29] * p_rints_x_ecoeffs[9];
            eri3_batch[1] += p_ecoeffs_ab[24] * p_rints_x_ecoeffs[4];
            eri3_batch[1] += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[15];
            eri3_batch[1] += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[5];
            eri3_batch[2] += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[2];
            eri3_batch[2] += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[7];
            eri3_batch[2] += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[11];
            eri3_batch[2] += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch[2] += p_ecoeffs_ab[56] * p_rints_x_ecoeffs[16];
            eri3_batch[2] += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[1];
            eri3_batch[2] += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[9];
            eri3_batch[2] += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[4];
            eri3_batch[2] += p_ecoeffs_ab[58] * p_rints_x_ecoeffs[18];
            eri3_batch[2] += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[3];
            eri3_batch[2] += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[8];
            eri3_batch[3] += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[17];
            eri3_batch[3] += p_ecoeffs_ab[66] * p_rints_x_ecoeffs[6];
            eri3_batch[3] += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[2];
            eri3_batch[3] += p_ecoeffs_ab[67] * p_rints_x_ecoeffs[7];
            eri3_batch[3] += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri3_batch[3] += p_ecoeffs_ab[61] * p_rints_x_ecoeffs[1];
            eri3_batch[3] += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[12];
            eri3_batch[3] += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[4];
            eri3_batch[3] += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[3];
            eri3_batch[3] += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[8];
            eri3_batch[4] += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[6];
            eri3_batch[4] += p_ecoeffs_ab[85] * p_rints_x_ecoeffs[5];
            eri3_batch[4] += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[2];
            eri3_batch[4] += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri3_batch[4] += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[1];
            eri3_batch[4] += p_ecoeffs_ab[83] * p_rints_x_ecoeffs[3];
            eri3_batch[4] += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[14];
            eri3_batch[4] += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[8];
            eri3_batch[5] += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[5];
            eri3_batch[5] += p_ecoeffs_ab[102] * p_rints_x_ecoeffs[2];
            eri3_batch[5] += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[13];
            eri3_batch[5] += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[10];
            eri3_batch[5] += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[7];
            eri3_batch[5] += p_ecoeffs_ab[100] * p_rints_x_ecoeffs[0];
            eri3_batch[5] += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[1];
            eri3_batch[5] += p_ecoeffs_ab[104] * p_rints_x_ecoeffs[4];
            eri3_batch[6] += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[5];
            eri3_batch[6] += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[2];
            eri3_batch[6] += p_ecoeffs_ab[131] * p_rints_x_ecoeffs[11];
            eri3_batch[6] += p_ecoeffs_ab[127] * p_rints_x_ecoeffs[7];
            eri3_batch[6] += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[0];
            eri3_batch[6] += p_ecoeffs_ab[136] * p_rints_x_ecoeffs[16];
            eri3_batch[6] += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[1];
            eri3_batch[6] += p_ecoeffs_ab[124] * p_rints_x_ecoeffs[4];
        }
}

template<> void lible::ints::two::eri2Kernel<3, 0>(const int cdepth_a, const int cdepth_b,
                                                   const double* exps_a, const double* exps_b,
                                                   const double* coords_a, const double* coords_b,
                                                   const double* ecoeffs_a, const double* ecoeffs_b_tsp,
                                                   double* eri2_batch)
{
    constexpr int la = 3, lb = 0;
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
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[12] += rints[12] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[13] += rints[13] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[16] += rints[16] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[18] += rints[18] * p_ecoeffs_b[0];
                    p_rints_x_ecoeffs[19] += rints[19] * p_ecoeffs_b[0];
        }
    }

    for (int ia = 0; ia < cdepth_a; ia++)
    {
            const double* p_ecoeffs_a = &ecoeffs_a[ia * n_ecoeffs_a];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];

            eri2_batch[0] += p_ecoeffs_a[17] * p_rints_x_ecoeffs[17];
            eri2_batch[0] += p_ecoeffs_a[6] * p_rints_x_ecoeffs[6];
            eri2_batch[0] += p_ecoeffs_a[2] * p_rints_x_ecoeffs[2];
            eri2_batch[0] += p_ecoeffs_a[19] * p_rints_x_ecoeffs[19];
            eri2_batch[0] += p_ecoeffs_a[7] * p_rints_x_ecoeffs[7];
            eri2_batch[0] += p_ecoeffs_a[0] * p_rints_x_ecoeffs[0];
            eri2_batch[0] += p_ecoeffs_a[1] * p_rints_x_ecoeffs[1];
            eri2_batch[0] += p_ecoeffs_a[12] * p_rints_x_ecoeffs[12];
            eri2_batch[0] += p_ecoeffs_a[9] * p_rints_x_ecoeffs[9];
            eri2_batch[0] += p_ecoeffs_a[4] * p_rints_x_ecoeffs[4];
            eri2_batch[0] += p_ecoeffs_a[3] * p_rints_x_ecoeffs[3];
            eri2_batch[0] += p_ecoeffs_a[8] * p_rints_x_ecoeffs[8];
            eri2_batch[1] += p_ecoeffs_a[26] * p_rints_x_ecoeffs[6];
            eri2_batch[1] += p_ecoeffs_a[25] * p_rints_x_ecoeffs[5];
            eri2_batch[1] += p_ecoeffs_a[22] * p_rints_x_ecoeffs[2];
            eri2_batch[1] += p_ecoeffs_a[33] * p_rints_x_ecoeffs[13];
            eri2_batch[1] += p_ecoeffs_a[30] * p_rints_x_ecoeffs[10];
            eri2_batch[1] += p_ecoeffs_a[27] * p_rints_x_ecoeffs[7];
            eri2_batch[1] += p_ecoeffs_a[20] * p_rints_x_ecoeffs[0];
            eri2_batch[1] += p_ecoeffs_a[21] * p_rints_x_ecoeffs[1];
            eri2_batch[1] += p_ecoeffs_a[29] * p_rints_x_ecoeffs[9];
            eri2_batch[1] += p_ecoeffs_a[24] * p_rints_x_ecoeffs[4];
            eri2_batch[1] += p_ecoeffs_a[35] * p_rints_x_ecoeffs[15];
            eri2_batch[1] += p_ecoeffs_a[23] * p_rints_x_ecoeffs[3];
            eri2_batch[2] += p_ecoeffs_a[45] * p_rints_x_ecoeffs[5];
            eri2_batch[2] += p_ecoeffs_a[42] * p_rints_x_ecoeffs[2];
            eri2_batch[2] += p_ecoeffs_a[47] * p_rints_x_ecoeffs[7];
            eri2_batch[2] += p_ecoeffs_a[51] * p_rints_x_ecoeffs[11];
            eri2_batch[2] += p_ecoeffs_a[40] * p_rints_x_ecoeffs[0];
            eri2_batch[2] += p_ecoeffs_a[56] * p_rints_x_ecoeffs[16];
            eri2_batch[2] += p_ecoeffs_a[41] * p_rints_x_ecoeffs[1];
            eri2_batch[2] += p_ecoeffs_a[49] * p_rints_x_ecoeffs[9];
            eri2_batch[2] += p_ecoeffs_a[44] * p_rints_x_ecoeffs[4];
            eri2_batch[2] += p_ecoeffs_a[58] * p_rints_x_ecoeffs[18];
            eri2_batch[2] += p_ecoeffs_a[43] * p_rints_x_ecoeffs[3];
            eri2_batch[2] += p_ecoeffs_a[48] * p_rints_x_ecoeffs[8];
            eri2_batch[3] += p_ecoeffs_a[77] * p_rints_x_ecoeffs[17];
            eri2_batch[3] += p_ecoeffs_a[66] * p_rints_x_ecoeffs[6];
            eri2_batch[3] += p_ecoeffs_a[62] * p_rints_x_ecoeffs[2];
            eri2_batch[3] += p_ecoeffs_a[67] * p_rints_x_ecoeffs[7];
            eri2_batch[3] += p_ecoeffs_a[60] * p_rints_x_ecoeffs[0];
            eri2_batch[3] += p_ecoeffs_a[61] * p_rints_x_ecoeffs[1];
            eri2_batch[3] += p_ecoeffs_a[72] * p_rints_x_ecoeffs[12];
            eri2_batch[3] += p_ecoeffs_a[64] * p_rints_x_ecoeffs[4];
            eri2_batch[3] += p_ecoeffs_a[63] * p_rints_x_ecoeffs[3];
            eri2_batch[3] += p_ecoeffs_a[68] * p_rints_x_ecoeffs[8];
            eri2_batch[4] += p_ecoeffs_a[86] * p_rints_x_ecoeffs[6];
            eri2_batch[4] += p_ecoeffs_a[85] * p_rints_x_ecoeffs[5];
            eri2_batch[4] += p_ecoeffs_a[82] * p_rints_x_ecoeffs[2];
            eri2_batch[4] += p_ecoeffs_a[80] * p_rints_x_ecoeffs[0];
            eri2_batch[4] += p_ecoeffs_a[81] * p_rints_x_ecoeffs[1];
            eri2_batch[4] += p_ecoeffs_a[83] * p_rints_x_ecoeffs[3];
            eri2_batch[4] += p_ecoeffs_a[94] * p_rints_x_ecoeffs[14];
            eri2_batch[4] += p_ecoeffs_a[88] * p_rints_x_ecoeffs[8];
            eri2_batch[5] += p_ecoeffs_a[105] * p_rints_x_ecoeffs[5];
            eri2_batch[5] += p_ecoeffs_a[102] * p_rints_x_ecoeffs[2];
            eri2_batch[5] += p_ecoeffs_a[113] * p_rints_x_ecoeffs[13];
            eri2_batch[5] += p_ecoeffs_a[110] * p_rints_x_ecoeffs[10];
            eri2_batch[5] += p_ecoeffs_a[107] * p_rints_x_ecoeffs[7];
            eri2_batch[5] += p_ecoeffs_a[100] * p_rints_x_ecoeffs[0];
            eri2_batch[5] += p_ecoeffs_a[101] * p_rints_x_ecoeffs[1];
            eri2_batch[5] += p_ecoeffs_a[104] * p_rints_x_ecoeffs[4];
            eri2_batch[6] += p_ecoeffs_a[125] * p_rints_x_ecoeffs[5];
            eri2_batch[6] += p_ecoeffs_a[122] * p_rints_x_ecoeffs[2];
            eri2_batch[6] += p_ecoeffs_a[131] * p_rints_x_ecoeffs[11];
            eri2_batch[6] += p_ecoeffs_a[127] * p_rints_x_ecoeffs[7];
            eri2_batch[6] += p_ecoeffs_a[120] * p_rints_x_ecoeffs[0];
            eri2_batch[6] += p_ecoeffs_a[136] * p_rints_x_ecoeffs[16];
            eri2_batch[6] += p_ecoeffs_a[121] * p_rints_x_ecoeffs[1];
            eri2_batch[6] += p_ecoeffs_a[124] * p_rints_x_ecoeffs[4];
    }
}
