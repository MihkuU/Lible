#include <lible/ints/twoel/eri_kernels.hpp>

template<> void lible::ints::two::eri4Kernel<0, 0, 2, 1>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 0, lb = 0, lc = 2, ld = 1;
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

                    p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_cd_tsp[255];
                    p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_cd_tsp[90];
                    p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_cd_tsp[30];
                    p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_cd_tsp[285];
                    p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_cd_tsp[105];
                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_cd_tsp[15];
                    p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_cd_tsp[180];
                    p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_cd_tsp[135];
                    p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_cd_tsp[60];
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[45];
                    p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_cd_tsp[120];
                    p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_cd_tsp[91];
                    p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_cd_tsp[76];
                    p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_cd_tsp[31];
                    p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_cd_tsp[196];
                    p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_cd_tsp[151];
                    p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_cd_tsp[106];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[16];
                    p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_cd_tsp[136];
                    p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_cd_tsp[61];
                    p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_cd_tsp[226];
                    p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_cd_tsp[46];
                    p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_cd_tsp[77];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[32];
                    p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_cd_tsp[107];
                    p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_cd_tsp[167];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_cd_tsp[242];
                    p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_cd_tsp[17];
                    p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_cd_tsp[137];
                    p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_cd_tsp[62];
                    p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_cd_tsp[272];
                    p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_cd_tsp[47];
                    p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_cd_tsp[122];
                    p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_cd_tsp[93];
                    p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_cd_tsp[3];
                    p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_cd_tsp[18];
                    p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_cd_tsp[138];
                    p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_cd_tsp[228];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[48];
                    p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_cd_tsp[94];
                    p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_cd_tsp[19];
                    p_rints_x_ecoeffs[4] += rints[12] * p_ecoeffs_cd_tsp[184];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[64];
                    p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_cd_tsp[49];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[95];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[80];
                    p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_cd_tsp[35];
                    p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_cd_tsp[5];
                    p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_cd_tsp[20];
                    p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_cd_tsp[50];
                    p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_cd_tsp[215];
                    p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_cd_tsp[125];
                    p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_cd_tsp[36];
                    p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_cd_tsp[6];
                    p_rints_x_ecoeffs[6] += rints[9] * p_ecoeffs_cd_tsp[141];
                    p_rints_x_ecoeffs[6] += rints[18] * p_ecoeffs_cd_tsp[276];
                    p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_cd_tsp[51];
                    p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_cd_tsp[126];
                    p_rints_x_ecoeffs[7] += rints[6] * p_ecoeffs_cd_tsp[97];
                    p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_cd_tsp[82];
                    p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_cd_tsp[37];
                    p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_cd_tsp[7];
                    p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_cd_tsp[22];
                    p_rints_x_ecoeffs[7] += rints[3] * p_ecoeffs_cd_tsp[52];
                    p_rints_x_ecoeffs[7] += rints[14] * p_ecoeffs_cd_tsp[217];
                    p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_cd_tsp[127];
                    p_rints_x_ecoeffs[8] += rints[17] * p_ecoeffs_cd_tsp[263];
                    p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_cd_tsp[38];
                    p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_cd_tsp[113];
                    p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[8] += rints[3] * p_ecoeffs_cd_tsp[53];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[128];
                    p_rints_x_ecoeffs[9] += rints[17] * p_ecoeffs_cd_tsp[264];
                    p_rints_x_ecoeffs[9] += rints[6] * p_ecoeffs_cd_tsp[99];
                    p_rints_x_ecoeffs[9] += rints[2] * p_ecoeffs_cd_tsp[39];
                    p_rints_x_ecoeffs[9] += rints[7] * p_ecoeffs_cd_tsp[114];
                    p_rints_x_ecoeffs[9] += rints[0] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[9] += rints[1] * p_ecoeffs_cd_tsp[24];
                    p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_cd_tsp[189];
                    p_rints_x_ecoeffs[9] += rints[4] * p_ecoeffs_cd_tsp[69];
                    p_rints_x_ecoeffs[9] += rints[3] * p_ecoeffs_cd_tsp[54];
                    p_rints_x_ecoeffs[9] += rints[8] * p_ecoeffs_cd_tsp[129];
                    p_rints_x_ecoeffs[10] += rints[5] * p_ecoeffs_cd_tsp[85];
                    p_rints_x_ecoeffs[10] += rints[2] * p_ecoeffs_cd_tsp[40];
                    p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_cd_tsp[205];
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_cd_tsp[160];
                    p_rints_x_ecoeffs[10] += rints[7] * p_ecoeffs_cd_tsp[115];
                    p_rints_x_ecoeffs[10] += rints[0] * p_ecoeffs_cd_tsp[10];
                    p_rints_x_ecoeffs[10] += rints[1] * p_ecoeffs_cd_tsp[25];
                    p_rints_x_ecoeffs[10] += rints[4] * p_ecoeffs_cd_tsp[70];
                    p_rints_x_ecoeffs[11] += rints[5] * p_ecoeffs_cd_tsp[86];
                    p_rints_x_ecoeffs[11] += rints[2] * p_ecoeffs_cd_tsp[41];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_cd_tsp[176];
                    p_rints_x_ecoeffs[11] += rints[7] * p_ecoeffs_cd_tsp[116];
                    p_rints_x_ecoeffs[11] += rints[0] * p_ecoeffs_cd_tsp[11];
                    p_rints_x_ecoeffs[11] += rints[16] * p_ecoeffs_cd_tsp[251];
                    p_rints_x_ecoeffs[11] += rints[1] * p_ecoeffs_cd_tsp[26];
                    p_rints_x_ecoeffs[11] += rints[4] * p_ecoeffs_cd_tsp[71];
                    p_rints_x_ecoeffs[12] += rints[6] * p_ecoeffs_cd_tsp[102];
                    p_rints_x_ecoeffs[12] += rints[5] * p_ecoeffs_cd_tsp[87];
                    p_rints_x_ecoeffs[12] += rints[2] * p_ecoeffs_cd_tsp[42];
                    p_rints_x_ecoeffs[12] += rints[0] * p_ecoeffs_cd_tsp[12];
                    p_rints_x_ecoeffs[12] += rints[1] * p_ecoeffs_cd_tsp[27];
                    p_rints_x_ecoeffs[12] += rints[3] * p_ecoeffs_cd_tsp[57];
                    p_rints_x_ecoeffs[12] += rints[14] * p_ecoeffs_cd_tsp[222];
                    p_rints_x_ecoeffs[12] += rints[8] * p_ecoeffs_cd_tsp[132];
                    p_rints_x_ecoeffs[13] += rints[5] * p_ecoeffs_cd_tsp[88];
                    p_rints_x_ecoeffs[13] += rints[2] * p_ecoeffs_cd_tsp[43];
                    p_rints_x_ecoeffs[13] += rints[11] * p_ecoeffs_cd_tsp[178];
                    p_rints_x_ecoeffs[13] += rints[0] * p_ecoeffs_cd_tsp[13];
                    p_rints_x_ecoeffs[13] += rints[1] * p_ecoeffs_cd_tsp[28];
                    p_rints_x_ecoeffs[13] += rints[4] * p_ecoeffs_cd_tsp[73];
                    p_rints_x_ecoeffs[14] += rints[5] * p_ecoeffs_cd_tsp[89];
                    p_rints_x_ecoeffs[14] += rints[2] * p_ecoeffs_cd_tsp[44];
                    p_rints_x_ecoeffs[14] += rints[13] * p_ecoeffs_cd_tsp[209];
                    p_rints_x_ecoeffs[14] += rints[7] * p_ecoeffs_cd_tsp[119];
                    p_rints_x_ecoeffs[14] += rints[0] * p_ecoeffs_cd_tsp[14];
                    p_rints_x_ecoeffs[14] += rints[1] * p_ecoeffs_cd_tsp[29];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch[2] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch[3] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri4_batch[4] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri4_batch[5] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri4_batch[6] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
            eri4_batch[7] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[7];
            eri4_batch[8] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[8];
            eri4_batch[9] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[9];
            eri4_batch[10] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[10];
            eri4_batch[11] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[11];
            eri4_batch[12] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[12];
            eri4_batch[13] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[13];
            eri4_batch[14] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[14];
        }
}

template<> void lible::ints::two::eri4Kernel<0, 0, 3, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 0, lb = 0, lc = 3, ld = 0;
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

                    p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_cd_tsp[119];
                    p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_cd_tsp[42];
                    p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_cd_tsp[14];
                    p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_cd_tsp[133];
                    p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_cd_tsp[49];
                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_cd_tsp[7];
                    p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_cd_tsp[84];
                    p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_cd_tsp[63];
                    p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_cd_tsp[28];
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[21];
                    p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_cd_tsp[56];
                    p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_cd_tsp[43];
                    p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_cd_tsp[36];
                    p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_cd_tsp[15];
                    p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_cd_tsp[92];
                    p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_cd_tsp[71];
                    p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_cd_tsp[50];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_cd_tsp[64];
                    p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_cd_tsp[29];
                    p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_cd_tsp[106];
                    p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_cd_tsp[22];
                    p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_cd_tsp[37];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[16];
                    p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_cd_tsp[51];
                    p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_cd_tsp[79];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_cd_tsp[114];
                    p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_cd_tsp[65];
                    p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_cd_tsp[30];
                    p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_cd_tsp[128];
                    p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_cd_tsp[23];
                    p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_cd_tsp[58];
                    p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_cd_tsp[122];
                    p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_cd_tsp[45];
                    p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_cd_tsp[17];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[52];
                    p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_cd_tsp[3];
                    p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_cd_tsp[10];
                    p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_cd_tsp[87];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[31];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[24];
                    p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_cd_tsp[59];
                    p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_cd_tsp[46];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[39];
                    p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_cd_tsp[18];
                    p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_cd_tsp[11];
                    p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_cd_tsp[25];
                    p_rints_x_ecoeffs[4] += rints[14] * p_ecoeffs_cd_tsp[102];
                    p_rints_x_ecoeffs[4] += rints[8] * p_ecoeffs_cd_tsp[60];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[40];
                    p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_cd_tsp[19];
                    p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_cd_tsp[96];
                    p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_cd_tsp[75];
                    p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_cd_tsp[54];
                    p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_cd_tsp[5];
                    p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_cd_tsp[12];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[33];
                    p_rints_x_ecoeffs[6] += rints[5] * p_ecoeffs_cd_tsp[41];
                    p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_cd_tsp[20];
                    p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_cd_tsp[83];
                    p_rints_x_ecoeffs[6] += rints[7] * p_ecoeffs_cd_tsp[55];
                    p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_cd_tsp[6];
                    p_rints_x_ecoeffs[6] += rints[16] * p_ecoeffs_cd_tsp[118];
                    p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_cd_tsp[13];
                    p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_cd_tsp[34];
                }
        }

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch[0] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch[1] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch[2] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch[3] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri4_batch[4] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri4_batch[5] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri4_batch[6] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
        }
}

template void lible::ints::two::eri3Kernel<0, 0, 3>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<0, 3>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);

