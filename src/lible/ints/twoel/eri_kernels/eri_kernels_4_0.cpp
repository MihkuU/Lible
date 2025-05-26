#include <lible/ints/twoel/eri_kernels.hpp>

template<> lible::vec4d
lible::ints::two::eri4Kernel<2, 2, 0, 0>(const int ipair_ab, const int ipair_cd,
                                         const std::vector<double> &ecoeffs_ab,
                                         const std::vector<double> &ecoeffs_cd_tsp,
                                         const ShellPairData &sp_data_ab,
                                         const ShellPairData &sp_data_cd)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sp_data_cd.coords[6 * ipair_cd];
    const double *coords_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]];

    constexpr int la = 2, lb = 2, lc = 0, ld = 0;
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

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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
                    p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_cd_tsp[0];
                }
        }

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[5] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[20];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[30];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[23];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[130] * p_rints_x_ecoeffs[25];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[116] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[20];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[135] * p_rints_x_ecoeffs[30];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[23];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[137] * p_rints_x_ecoeffs[32];
            eri4_batch(0, 3, 0, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[156] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[166] * p_rints_x_ecoeffs[26];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[158] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[28];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[161] * p_rints_x_ecoeffs[21];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 4, 0, 0) += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[17];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[27];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[194] * p_rints_x_ecoeffs[19];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[204] * p_rints_x_ecoeffs[29];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[22];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[190] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[235] * p_rints_x_ecoeffs[25];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[219] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[225] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[251] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[273] * p_rints_x_ecoeffs[28];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[263] * p_rints_x_ecoeffs[18];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[259] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[253] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[307] * p_rints_x_ecoeffs[27];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[297] * p_rints_x_ecoeffs[17];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[286] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[290] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[292] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[302] * p_rints_x_ecoeffs[22];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[283] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[294] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 3, 0, 0) += p_ecoeffs_ab[288] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[339] * p_rints_x_ecoeffs[24];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[321] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[320] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[317] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[326] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[315] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[316] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[327] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[319] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[318] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[329] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 4, 0, 0) += p_ecoeffs_ab[323] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[352] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[361] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[383] * p_rints_x_ecoeffs[33];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[366] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[353] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[367] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[356] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[355] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[362] * p_rints_x_ecoeffs[12];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[368] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[369] * p_rints_x_ecoeffs[19];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[357] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[350] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[351] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[354] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[358] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[374] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[381] * p_rints_x_ecoeffs[31];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[359] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[364] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[391] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[390] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[387] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[385] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[413] * p_rints_x_ecoeffs[28];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[386] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[394] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[403] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[400] * p_rints_x_ecoeffs[15];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[388] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[399] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[393] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[437] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[422] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[427] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[420] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[429] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[438] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[452] * p_rints_x_ecoeffs[32];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[423] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[428] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[479] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[461] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[472] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[460] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[486] * p_rints_x_ecoeffs[31];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[457] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[466] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[462] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[455] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[471] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[456] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[467] * p_rints_x_ecoeffs[12];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[459] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[458] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[469] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 3, 0, 0) += p_ecoeffs_ab[463] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[517] * p_rints_x_ecoeffs[27];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[507] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[496] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[495] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[492] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[503] * p_rints_x_ecoeffs[13];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[497] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[490] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[491] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[493] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[504] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 4, 0, 0) += p_ecoeffs_ab[498] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[550] * p_rints_x_ecoeffs[25];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[527] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[536] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[541] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[528] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[542] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[531] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[530] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[535] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[545] * p_rints_x_ecoeffs[20];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[537] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[543] * p_rints_x_ecoeffs[18];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[532] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[525] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[555] * p_rints_x_ecoeffs[30];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[526] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[529] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[548] * p_rints_x_ecoeffs[23];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[533] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[538] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[534] * p_rints_x_ecoeffs[9];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[557] * p_rints_x_ecoeffs[32];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[540] * p_rints_x_ecoeffs[15];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[587] * p_rints_x_ecoeffs[27];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[577] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[566] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[565] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[562] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[573] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[567] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[570] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[560] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[561] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[572] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[564] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[582] * p_rints_x_ecoeffs[22];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[563] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[574] * p_rints_x_ecoeffs[14];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[568] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[619] * p_rints_x_ecoeffs[24];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[601] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[612] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[600] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[626] * p_rints_x_ecoeffs[31];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[597] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[606] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[602] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[595] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[611] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[596] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[607] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[599] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[598] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[609] * p_rints_x_ecoeffs[14];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[603] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[635] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[632] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[643] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[637] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[641] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[630] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[646] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[640] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[660] * p_rints_x_ecoeffs[30];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[650] * p_rints_x_ecoeffs[20];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[631] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[634] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 3, 0, 0) += p_ecoeffs_ab[653] * p_rints_x_ecoeffs[23];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[670] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[667] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[678] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[672] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[676] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[665] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[681] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[675] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[686] * p_rints_x_ecoeffs[21];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[666] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[669] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 4, 0, 0) += p_ecoeffs_ab[691] * p_rints_x_ecoeffs[26];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[702] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[711] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[716] * p_rints_x_ecoeffs[16];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[726] * p_rints_x_ecoeffs[26];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[703] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[706] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[705] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[710] * p_rints_x_ecoeffs[10];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[718] * p_rints_x_ecoeffs[18];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[707] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[700] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[728] * p_rints_x_ecoeffs[28];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[701] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[704] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[708] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[713] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[721] * p_rints_x_ecoeffs[21];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[709] * p_rints_x_ecoeffs[9];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[715] * p_rints_x_ecoeffs[15];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[714] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[759] * p_rints_x_ecoeffs[24];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[741] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[740] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[737] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[746] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[735] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[736] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[747] * p_rints_x_ecoeffs[12];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[739] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[738] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[749] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[743] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[797] * p_rints_x_ecoeffs[27];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[787] * p_rints_x_ecoeffs[17];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[776] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[775] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[772] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[783] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[777] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[770] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[771] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[773] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[784] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[778] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[810] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[807] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[818] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[812] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[816] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[805] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[821] * p_rints_x_ecoeffs[16];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[815] * p_rints_x_ecoeffs[10];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[826] * p_rints_x_ecoeffs[21];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[806] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[809] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 3, 0, 0) += p_ecoeffs_ab[831] * p_rints_x_ecoeffs[26];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[845] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[842] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[853] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[847] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[851] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[840] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[841] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[844] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 4, 0, 0) += p_ecoeffs_ab[863] * p_rints_x_ecoeffs[23];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
                for (int ta = 0; ta < n_sph_d; ta++)
                {
                    double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                    double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                    double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                    double norm_d = sp_data_cd.norms[ofs_norm_d + ta];
                    eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;
                }

    return eri4_batch;
}

template<> lible::vec4d
lible::ints::two::eri4Kernel<3, 1, 0, 0>(const int ipair_ab, const int ipair_cd,
                                         const std::vector<double> &ecoeffs_ab,
                                         const std::vector<double> &ecoeffs_cd_tsp,
                                         const ShellPairData &sp_data_ab,
                                         const ShellPairData &sp_data_cd)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sp_data_cd.coords[6 * ipair_cd];
    const double *coords_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]];

    constexpr int la = 3, lb = 1, lc = 0, ld = 0;
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

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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
                    p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_cd_tsp[0];
                }
        }

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[132] * p_rints_x_ecoeffs[27];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[124] * p_rints_x_ecoeffs[19];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[134] * p_rints_x_ecoeffs[29];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[127] * p_rints_x_ecoeffs[22];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[119] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[165] * p_rints_x_ecoeffs[25];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[160] * p_rints_x_ecoeffs[20];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[152] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[163] * p_rints_x_ecoeffs[23];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[186] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[191] * p_rints_x_ecoeffs[16];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[201] * p_rints_x_ecoeffs[26];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[193] * p_rints_x_ecoeffs[18];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[203] * p_rints_x_ecoeffs[28];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[196] * p_rints_x_ecoeffs[21];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[190] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[212] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[243] * p_rints_x_ecoeffs[33];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[226] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[215] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[228] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[229] * p_rints_x_ecoeffs[19];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[217] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[218] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[234] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[31];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[219] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[256] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[271] * p_rints_x_ecoeffs[26];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[251] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[255] * p_rints_x_ecoeffs[10];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[263] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[252] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[273] * p_rints_x_ecoeffs[28];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[249] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[253] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[258] * p_rints_x_ecoeffs[13];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[266] * p_rints_x_ecoeffs[21];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[15];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[259] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[297] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[303] * p_rints_x_ecoeffs[23];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[291] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[296] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[310] * p_rints_x_ecoeffs[30];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[289] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[298] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[312] * p_rints_x_ecoeffs[32];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[283] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[288] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[332] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[321] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[340] * p_rints_x_ecoeffs[25];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[317] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[322] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[315] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[316] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[327] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[324] * p_rints_x_ecoeffs[9];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[319] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[333] * p_rints_x_ecoeffs[18];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[347] * p_rints_x_ecoeffs[32];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[330] * p_rints_x_ecoeffs[15];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[318] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[323] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[377] * p_rints_x_ecoeffs[27];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[367] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[356] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[355] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[352] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[363] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[357] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[360] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[350] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[351] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[362] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[354] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[372] * p_rints_x_ecoeffs[22];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[353] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[364] * p_rints_x_ecoeffs[14];
            eri4_batch(3, 1, 0, 0) += p_ecoeffs_ab[358] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[409] * p_rints_x_ecoeffs[24];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[391] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[402] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[390] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[416] * p_rints_x_ecoeffs[31];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[387] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[396] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[392] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[385] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[401] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[386] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[397] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[389] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[388] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[399] * p_rints_x_ecoeffs[14];
            eri4_batch(3, 2, 0, 0) += p_ecoeffs_ab[393] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[426] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[425] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[422] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[420] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[448] * p_rints_x_ecoeffs[28];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[421] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[429] * p_rints_x_ecoeffs[9];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[438] * p_rints_x_ecoeffs[18];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[435] * p_rints_x_ecoeffs[15];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[423] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[434] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[428] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[479] * p_rints_x_ecoeffs[24];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[461] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[460] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[457] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[466] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[455] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[456] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[467] * p_rints_x_ecoeffs[12];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[459] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[458] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[469] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 1, 0, 0) += p_ecoeffs_ab[463] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[517] * p_rints_x_ecoeffs[27];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[507] * p_rints_x_ecoeffs[17];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[496] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[495] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[492] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[503] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[497] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[490] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[491] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[493] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[504] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 2, 0, 0) += p_ecoeffs_ab[498] * p_rints_x_ecoeffs[8];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[552] * p_rints_x_ecoeffs[27];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[542] * p_rints_x_ecoeffs[17];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[531] * p_rints_x_ecoeffs[6];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[530] * p_rints_x_ecoeffs[5];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[527] * p_rints_x_ecoeffs[2];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[538] * p_rints_x_ecoeffs[13];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[532] * p_rints_x_ecoeffs[7];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[535] * p_rints_x_ecoeffs[10];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[525] * p_rints_x_ecoeffs[0];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[526] * p_rints_x_ecoeffs[1];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[537] * p_rints_x_ecoeffs[12];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[529] * p_rints_x_ecoeffs[4];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[547] * p_rints_x_ecoeffs[22];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[528] * p_rints_x_ecoeffs[3];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[539] * p_rints_x_ecoeffs[14];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[533] * p_rints_x_ecoeffs[8];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[565] * p_rints_x_ecoeffs[5];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[562] * p_rints_x_ecoeffs[2];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[573] * p_rints_x_ecoeffs[13];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[567] * p_rints_x_ecoeffs[7];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[571] * p_rints_x_ecoeffs[11];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[560] * p_rints_x_ecoeffs[0];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[570] * p_rints_x_ecoeffs[10];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[580] * p_rints_x_ecoeffs[20];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[561] * p_rints_x_ecoeffs[1];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[564] * p_rints_x_ecoeffs[4];
            eri4_batch(5, 1, 0, 0) += p_ecoeffs_ab[583] * p_rints_x_ecoeffs[23];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[600] * p_rints_x_ecoeffs[5];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[597] * p_rints_x_ecoeffs[2];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[608] * p_rints_x_ecoeffs[13];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[602] * p_rints_x_ecoeffs[7];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[606] * p_rints_x_ecoeffs[11];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[595] * p_rints_x_ecoeffs[0];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[611] * p_rints_x_ecoeffs[16];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[605] * p_rints_x_ecoeffs[10];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[616] * p_rints_x_ecoeffs[21];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[596] * p_rints_x_ecoeffs[1];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[599] * p_rints_x_ecoeffs[4];
            eri4_batch(5, 2, 0, 0) += p_ecoeffs_ab[621] * p_rints_x_ecoeffs[26];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[654] * p_rints_x_ecoeffs[24];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[636] * p_rints_x_ecoeffs[6];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[647] * p_rints_x_ecoeffs[17];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[635] * p_rints_x_ecoeffs[5];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[661] * p_rints_x_ecoeffs[31];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[632] * p_rints_x_ecoeffs[2];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[641] * p_rints_x_ecoeffs[11];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[637] * p_rints_x_ecoeffs[7];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[630] * p_rints_x_ecoeffs[0];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[646] * p_rints_x_ecoeffs[16];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[631] * p_rints_x_ecoeffs[1];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[642] * p_rints_x_ecoeffs[12];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[634] * p_rints_x_ecoeffs[4];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[633] * p_rints_x_ecoeffs[3];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[644] * p_rints_x_ecoeffs[14];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[638] * p_rints_x_ecoeffs[8];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[670] * p_rints_x_ecoeffs[5];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[667] * p_rints_x_ecoeffs[2];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[678] * p_rints_x_ecoeffs[13];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[672] * p_rints_x_ecoeffs[7];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[676] * p_rints_x_ecoeffs[11];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[665] * p_rints_x_ecoeffs[0];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[681] * p_rints_x_ecoeffs[16];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[675] * p_rints_x_ecoeffs[10];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[686] * p_rints_x_ecoeffs[21];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[666] * p_rints_x_ecoeffs[1];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[669] * p_rints_x_ecoeffs[4];
            eri4_batch(6, 1, 0, 0) += p_ecoeffs_ab[691] * p_rints_x_ecoeffs[26];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[705] * p_rints_x_ecoeffs[5];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[702] * p_rints_x_ecoeffs[2];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[713] * p_rints_x_ecoeffs[13];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[707] * p_rints_x_ecoeffs[7];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[711] * p_rints_x_ecoeffs[11];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[700] * p_rints_x_ecoeffs[0];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[716] * p_rints_x_ecoeffs[16];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[730] * p_rints_x_ecoeffs[30];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[701] * p_rints_x_ecoeffs[1];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[704] * p_rints_x_ecoeffs[4];
            eri4_batch(6, 2, 0, 0) += p_ecoeffs_ab[723] * p_rints_x_ecoeffs[23];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
                for (int ta = 0; ta < n_sph_d; ta++)
                {
                    double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                    double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                    double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                    double norm_d = sp_data_cd.norms[ofs_norm_d + ta];
                    eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;
                }

    return eri4_batch;
}

template<> lible::vec4d
lible::ints::two::eri4Kernel<4, 0, 0, 0>(const int ipair_ab, const int ipair_cd,
                                         const std::vector<double> &ecoeffs_ab,
                                         const std::vector<double> &ecoeffs_cd_tsp,
                                         const ShellPairData &sp_data_ab,
                                         const ShellPairData &sp_data_cd)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sp_data_cd.coords[6 * ipair_cd];
    const double *coords_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]];

    constexpr int la = 4, lb = 0, lc = 0, ld = 0;
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

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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
                    p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_cd_tsp[0];
                }
        }

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[5] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[20];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[30];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[23];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[130] * p_rints_x_ecoeffs[25];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[116] * p_rints_x_ecoeffs[11];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[16];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[20];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[18];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[135] * p_rints_x_ecoeffs[30];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[23];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[137] * p_rints_x_ecoeffs[32];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[156] * p_rints_x_ecoeffs[16];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[166] * p_rints_x_ecoeffs[26];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[158] * p_rints_x_ecoeffs[18];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[28];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[161] * p_rints_x_ecoeffs[21];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[27];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[17];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[12];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[22];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri4_batch(5, 0, 0, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[234] * p_rints_x_ecoeffs[24];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[17];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[215] * p_rints_x_ecoeffs[5];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[31];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[212] * p_rints_x_ecoeffs[2];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[11];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[217] * p_rints_x_ecoeffs[7];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[226] * p_rints_x_ecoeffs[16];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[14];
            eri4_batch(6, 0, 0, 0) += p_ecoeffs_ab[218] * p_rints_x_ecoeffs[8];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[258] * p_rints_x_ecoeffs[13];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[252] * p_rints_x_ecoeffs[7];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[256] * p_rints_x_ecoeffs[11];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[16];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[255] * p_rints_x_ecoeffs[10];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[275] * p_rints_x_ecoeffs[30];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[265] * p_rints_x_ecoeffs[20];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[249] * p_rints_x_ecoeffs[4];
            eri4_batch(7, 0, 0, 0) += p_ecoeffs_ab[268] * p_rints_x_ecoeffs[23];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[291] * p_rints_x_ecoeffs[11];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[296] * p_rints_x_ecoeffs[16];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[290] * p_rints_x_ecoeffs[10];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[301] * p_rints_x_ecoeffs[21];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri4_batch(8, 0, 0, 0) += p_ecoeffs_ab[306] * p_rints_x_ecoeffs[26];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
                for (int ta = 0; ta < n_sph_d; ta++)
                {
                    double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                    double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                    double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                    double norm_d = sp_data_cd.norms[ofs_norm_d + ta];
                    eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;
                }

    return eri4_batch;
}

template<> lible::vec3d
lible::ints::two::eri3Kernel<2, 2, 0>(const int ipair_ab, const int ishell_c,
                                      const std::vector<double> &ecoeffs_ab,
                                      const std::vector<double> &ecoeffs_c,
                                      const ShellPairData &sp_data_ab,
                                      const ShellData &sh_data_c)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_c = &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]];

    constexpr int la = 2, lb = 2, lc = 0;
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

                const double* p_ecoeffs_c = &pecoeffs_c[ic * n_ecoeffs_c];
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
                p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_c[0];
            }
}

    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[5] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[20];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[30];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[23];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[130] * p_rints_x_ecoeffs[25];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[116] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[20];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[135] * p_rints_x_ecoeffs[30];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[23];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[137] * p_rints_x_ecoeffs[32];
            eri3_batch(0, 3, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[156] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[166] * p_rints_x_ecoeffs[26];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[158] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[28];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[161] * p_rints_x_ecoeffs[21];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 4, 0) += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[17];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[27];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[194] * p_rints_x_ecoeffs[19];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[204] * p_rints_x_ecoeffs[29];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[22];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[190] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[235] * p_rints_x_ecoeffs[25];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[219] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[225] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[251] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[273] * p_rints_x_ecoeffs[28];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[263] * p_rints_x_ecoeffs[18];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[259] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[253] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[307] * p_rints_x_ecoeffs[27];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[297] * p_rints_x_ecoeffs[17];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[286] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[290] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[292] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[302] * p_rints_x_ecoeffs[22];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[283] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[294] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 3, 0) += p_ecoeffs_ab[288] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[339] * p_rints_x_ecoeffs[24];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[321] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[320] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[317] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[326] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[315] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[316] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[327] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[319] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[318] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[329] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 4, 0) += p_ecoeffs_ab[323] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[352] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[361] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[383] * p_rints_x_ecoeffs[33];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[366] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[353] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[367] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[356] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[355] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[362] * p_rints_x_ecoeffs[12];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[368] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[369] * p_rints_x_ecoeffs[19];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[357] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[350] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[351] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[354] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[358] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[374] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[381] * p_rints_x_ecoeffs[31];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[359] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[364] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[391] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[390] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[387] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[385] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[413] * p_rints_x_ecoeffs[28];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[386] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[394] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[403] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[400] * p_rints_x_ecoeffs[15];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[388] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[399] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[393] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[437] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[422] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[427] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[420] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[429] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[438] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[452] * p_rints_x_ecoeffs[32];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[423] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[428] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[479] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[461] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[472] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[460] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[486] * p_rints_x_ecoeffs[31];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[457] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[466] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[462] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[455] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[471] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[456] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[467] * p_rints_x_ecoeffs[12];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[459] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[458] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[469] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 3, 0) += p_ecoeffs_ab[463] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[517] * p_rints_x_ecoeffs[27];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[507] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[496] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[495] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[492] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[503] * p_rints_x_ecoeffs[13];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[497] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[490] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[491] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[493] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[504] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 4, 0) += p_ecoeffs_ab[498] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[550] * p_rints_x_ecoeffs[25];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[527] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[536] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[541] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[528] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[542] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[531] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[530] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[535] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[545] * p_rints_x_ecoeffs[20];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[537] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[543] * p_rints_x_ecoeffs[18];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[532] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[525] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[555] * p_rints_x_ecoeffs[30];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[526] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[529] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[548] * p_rints_x_ecoeffs[23];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[533] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[538] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[534] * p_rints_x_ecoeffs[9];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[557] * p_rints_x_ecoeffs[32];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[540] * p_rints_x_ecoeffs[15];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[587] * p_rints_x_ecoeffs[27];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[577] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[566] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[565] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[562] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[573] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[567] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[570] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[560] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[561] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[572] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[564] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[582] * p_rints_x_ecoeffs[22];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[563] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[574] * p_rints_x_ecoeffs[14];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[568] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[619] * p_rints_x_ecoeffs[24];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[601] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[612] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[600] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[626] * p_rints_x_ecoeffs[31];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[597] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[606] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[602] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[595] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[611] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[596] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[607] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[599] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[598] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[609] * p_rints_x_ecoeffs[14];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[603] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[635] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[632] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[643] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[637] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[641] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[630] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[646] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[640] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[660] * p_rints_x_ecoeffs[30];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[650] * p_rints_x_ecoeffs[20];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[631] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[634] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 3, 0) += p_ecoeffs_ab[653] * p_rints_x_ecoeffs[23];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[670] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[667] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[678] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[672] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[676] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[665] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[681] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[675] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[686] * p_rints_x_ecoeffs[21];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[666] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[669] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 4, 0) += p_ecoeffs_ab[691] * p_rints_x_ecoeffs[26];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[702] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[711] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[716] * p_rints_x_ecoeffs[16];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[726] * p_rints_x_ecoeffs[26];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[703] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[706] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[705] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[710] * p_rints_x_ecoeffs[10];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[718] * p_rints_x_ecoeffs[18];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[707] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[700] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[728] * p_rints_x_ecoeffs[28];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[701] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[704] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[708] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[713] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[721] * p_rints_x_ecoeffs[21];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[709] * p_rints_x_ecoeffs[9];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[715] * p_rints_x_ecoeffs[15];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[714] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[759] * p_rints_x_ecoeffs[24];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[741] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[740] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[737] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[746] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[735] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[736] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[747] * p_rints_x_ecoeffs[12];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[739] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[738] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[749] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[743] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[797] * p_rints_x_ecoeffs[27];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[787] * p_rints_x_ecoeffs[17];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[776] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[775] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[772] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[783] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[777] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[770] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[771] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[773] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[784] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[778] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[810] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[807] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[818] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[812] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[816] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[805] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[821] * p_rints_x_ecoeffs[16];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[815] * p_rints_x_ecoeffs[10];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[826] * p_rints_x_ecoeffs[21];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[806] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[809] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 3, 0) += p_ecoeffs_ab[831] * p_rints_x_ecoeffs[26];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[845] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[842] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[853] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[847] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[851] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[840] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[841] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[844] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 4, 0) += p_ecoeffs_ab[863] * p_rints_x_ecoeffs[23];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
            {
                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                double norm_c = sh_data_c.norms[ofs_norm_c + ka];
                eri3_batch(mu, nu, ka) *= norm_a * norm_b * norm_c;
            }

    return eri3_batch;
}

template<> lible::vec3d
lible::ints::two::eri3Kernel<3, 1, 0>(const int ipair_ab, const int ishell_c,
                                      const std::vector<double> &ecoeffs_ab,
                                      const std::vector<double> &ecoeffs_c,
                                      const ShellPairData &sp_data_ab,
                                      const ShellData &sh_data_c)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_c = &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]];

    constexpr int la = 3, lb = 1, lc = 0;
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

                const double* p_ecoeffs_c = &pecoeffs_c[ic * n_ecoeffs_c];
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
                p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_c[0];
            }
}

    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[132] * p_rints_x_ecoeffs[27];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[124] * p_rints_x_ecoeffs[19];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[134] * p_rints_x_ecoeffs[29];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[127] * p_rints_x_ecoeffs[22];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[119] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[165] * p_rints_x_ecoeffs[25];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[160] * p_rints_x_ecoeffs[20];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[152] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[163] * p_rints_x_ecoeffs[23];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[186] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[191] * p_rints_x_ecoeffs[16];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[201] * p_rints_x_ecoeffs[26];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[193] * p_rints_x_ecoeffs[18];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[203] * p_rints_x_ecoeffs[28];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[196] * p_rints_x_ecoeffs[21];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[184] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[190] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[212] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[243] * p_rints_x_ecoeffs[33];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[226] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[215] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[228] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[229] * p_rints_x_ecoeffs[19];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[217] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[218] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[234] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[31];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[219] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[256] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[271] * p_rints_x_ecoeffs[26];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[248] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[251] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[255] * p_rints_x_ecoeffs[10];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[263] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[252] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[273] * p_rints_x_ecoeffs[28];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[249] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[253] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[258] * p_rints_x_ecoeffs[13];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[266] * p_rints_x_ecoeffs[21];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[254] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[260] * p_rints_x_ecoeffs[15];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[259] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[297] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[303] * p_rints_x_ecoeffs[23];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[291] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[296] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[310] * p_rints_x_ecoeffs[30];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[289] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[298] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[312] * p_rints_x_ecoeffs[32];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[283] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[288] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[332] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[321] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[340] * p_rints_x_ecoeffs[25];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[317] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[322] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[315] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[316] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[327] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[324] * p_rints_x_ecoeffs[9];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[319] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[333] * p_rints_x_ecoeffs[18];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[347] * p_rints_x_ecoeffs[32];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[330] * p_rints_x_ecoeffs[15];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[318] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[323] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[377] * p_rints_x_ecoeffs[27];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[367] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[356] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[355] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[352] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[363] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[357] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[360] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[350] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[351] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[362] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[354] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[372] * p_rints_x_ecoeffs[22];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[353] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[364] * p_rints_x_ecoeffs[14];
            eri3_batch(3, 1, 0) += p_ecoeffs_ab[358] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[409] * p_rints_x_ecoeffs[24];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[391] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[402] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[390] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[416] * p_rints_x_ecoeffs[31];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[387] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[396] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[392] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[385] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[401] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[386] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[397] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[389] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[388] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[399] * p_rints_x_ecoeffs[14];
            eri3_batch(3, 2, 0) += p_ecoeffs_ab[393] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[426] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[425] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[422] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[420] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[448] * p_rints_x_ecoeffs[28];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[421] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[429] * p_rints_x_ecoeffs[9];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[438] * p_rints_x_ecoeffs[18];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[435] * p_rints_x_ecoeffs[15];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[423] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[434] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[428] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[479] * p_rints_x_ecoeffs[24];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[461] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[460] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[457] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[466] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[455] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[456] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[467] * p_rints_x_ecoeffs[12];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[459] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[458] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[469] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 1, 0) += p_ecoeffs_ab[463] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[517] * p_rints_x_ecoeffs[27];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[507] * p_rints_x_ecoeffs[17];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[496] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[495] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[492] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[503] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[497] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[490] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[491] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[493] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[504] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 2, 0) += p_ecoeffs_ab[498] * p_rints_x_ecoeffs[8];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[552] * p_rints_x_ecoeffs[27];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[542] * p_rints_x_ecoeffs[17];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[531] * p_rints_x_ecoeffs[6];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[530] * p_rints_x_ecoeffs[5];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[527] * p_rints_x_ecoeffs[2];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[538] * p_rints_x_ecoeffs[13];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[532] * p_rints_x_ecoeffs[7];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[535] * p_rints_x_ecoeffs[10];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[525] * p_rints_x_ecoeffs[0];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[526] * p_rints_x_ecoeffs[1];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[537] * p_rints_x_ecoeffs[12];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[529] * p_rints_x_ecoeffs[4];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[547] * p_rints_x_ecoeffs[22];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[528] * p_rints_x_ecoeffs[3];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[539] * p_rints_x_ecoeffs[14];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[533] * p_rints_x_ecoeffs[8];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[565] * p_rints_x_ecoeffs[5];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[562] * p_rints_x_ecoeffs[2];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[573] * p_rints_x_ecoeffs[13];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[567] * p_rints_x_ecoeffs[7];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[571] * p_rints_x_ecoeffs[11];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[560] * p_rints_x_ecoeffs[0];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[570] * p_rints_x_ecoeffs[10];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[580] * p_rints_x_ecoeffs[20];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[561] * p_rints_x_ecoeffs[1];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[564] * p_rints_x_ecoeffs[4];
            eri3_batch(5, 1, 0) += p_ecoeffs_ab[583] * p_rints_x_ecoeffs[23];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[600] * p_rints_x_ecoeffs[5];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[597] * p_rints_x_ecoeffs[2];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[608] * p_rints_x_ecoeffs[13];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[602] * p_rints_x_ecoeffs[7];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[606] * p_rints_x_ecoeffs[11];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[595] * p_rints_x_ecoeffs[0];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[611] * p_rints_x_ecoeffs[16];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[605] * p_rints_x_ecoeffs[10];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[616] * p_rints_x_ecoeffs[21];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[596] * p_rints_x_ecoeffs[1];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[599] * p_rints_x_ecoeffs[4];
            eri3_batch(5, 2, 0) += p_ecoeffs_ab[621] * p_rints_x_ecoeffs[26];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[654] * p_rints_x_ecoeffs[24];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[636] * p_rints_x_ecoeffs[6];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[647] * p_rints_x_ecoeffs[17];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[635] * p_rints_x_ecoeffs[5];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[661] * p_rints_x_ecoeffs[31];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[632] * p_rints_x_ecoeffs[2];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[641] * p_rints_x_ecoeffs[11];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[637] * p_rints_x_ecoeffs[7];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[630] * p_rints_x_ecoeffs[0];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[646] * p_rints_x_ecoeffs[16];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[631] * p_rints_x_ecoeffs[1];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[642] * p_rints_x_ecoeffs[12];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[634] * p_rints_x_ecoeffs[4];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[633] * p_rints_x_ecoeffs[3];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[644] * p_rints_x_ecoeffs[14];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[638] * p_rints_x_ecoeffs[8];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[670] * p_rints_x_ecoeffs[5];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[667] * p_rints_x_ecoeffs[2];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[678] * p_rints_x_ecoeffs[13];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[672] * p_rints_x_ecoeffs[7];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[676] * p_rints_x_ecoeffs[11];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[665] * p_rints_x_ecoeffs[0];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[681] * p_rints_x_ecoeffs[16];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[675] * p_rints_x_ecoeffs[10];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[686] * p_rints_x_ecoeffs[21];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[666] * p_rints_x_ecoeffs[1];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[669] * p_rints_x_ecoeffs[4];
            eri3_batch(6, 1, 0) += p_ecoeffs_ab[691] * p_rints_x_ecoeffs[26];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[705] * p_rints_x_ecoeffs[5];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[702] * p_rints_x_ecoeffs[2];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[713] * p_rints_x_ecoeffs[13];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[707] * p_rints_x_ecoeffs[7];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[711] * p_rints_x_ecoeffs[11];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[700] * p_rints_x_ecoeffs[0];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[716] * p_rints_x_ecoeffs[16];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[730] * p_rints_x_ecoeffs[30];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[701] * p_rints_x_ecoeffs[1];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[704] * p_rints_x_ecoeffs[4];
            eri3_batch(6, 2, 0) += p_ecoeffs_ab[723] * p_rints_x_ecoeffs[23];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
            {
                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                double norm_c = sh_data_c.norms[ofs_norm_c + ka];
                eri3_batch(mu, nu, ka) *= norm_a * norm_b * norm_c;
            }

    return eri3_batch;
}

template<> lible::vec3d
lible::ints::two::eri3Kernel<4, 0, 0>(const int ipair_ab, const int ishell_c,
                                      const std::vector<double> &ecoeffs_ab,
                                      const std::vector<double> &ecoeffs_c,
                                      const ShellPairData &sp_data_ab,
                                      const ShellData &sh_data_c)
{
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
    const double *pecoeffs_c = &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]];

    constexpr int la = 4, lb = 0, lc = 0;
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

                const double* p_ecoeffs_c = &pecoeffs_c[ic * n_ecoeffs_c];
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
                p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_c[0];
            }
}

    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[25] * p_rints_x_ecoeffs[25];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[16];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[34];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[17] * p_rints_x_ecoeffs[17];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[6] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[5] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[20];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[12] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[18] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[19] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[30];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[23];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[8] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[32];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[15] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[38] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[17];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[47] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[27];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[54] * p_rints_x_ecoeffs[19];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[35] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[64] * p_rints_x_ecoeffs[29];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[39] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[43] * p_rints_x_ecoeffs[8];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[48] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[57] * p_rints_x_ecoeffs[22];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[49] * p_rints_x_ecoeffs[14];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[81] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[103] * p_rints_x_ecoeffs[33];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[86] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[73] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[76] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[12];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[88] * p_rints_x_ecoeffs[18];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[89] * p_rints_x_ecoeffs[19];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[77] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[74] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[78] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[94] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[101] * p_rints_x_ecoeffs[31];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[79] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[84] * p_rints_x_ecoeffs[14];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[130] * p_rints_x_ecoeffs[25];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[107] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[116] * p_rints_x_ecoeffs[11];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[121] * p_rints_x_ecoeffs[16];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[108] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[122] * p_rints_x_ecoeffs[17];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[111] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[110] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[115] * p_rints_x_ecoeffs[10];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[125] * p_rints_x_ecoeffs[20];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[117] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[123] * p_rints_x_ecoeffs[18];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[112] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[105] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[135] * p_rints_x_ecoeffs[30];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[106] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[109] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[128] * p_rints_x_ecoeffs[23];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[113] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[118] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[114] * p_rints_x_ecoeffs[9];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[137] * p_rints_x_ecoeffs[32];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[120] * p_rints_x_ecoeffs[15];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[142] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[151] * p_rints_x_ecoeffs[11];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[156] * p_rints_x_ecoeffs[16];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[166] * p_rints_x_ecoeffs[26];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[143] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[146] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[145] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[150] * p_rints_x_ecoeffs[10];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[158] * p_rints_x_ecoeffs[18];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[147] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[140] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[168] * p_rints_x_ecoeffs[28];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[141] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[144] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[148] * p_rints_x_ecoeffs[8];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[153] * p_rints_x_ecoeffs[13];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[161] * p_rints_x_ecoeffs[21];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[149] * p_rints_x_ecoeffs[9];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[155] * p_rints_x_ecoeffs[15];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[154] * p_rints_x_ecoeffs[14];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[202] * p_rints_x_ecoeffs[27];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[192] * p_rints_x_ecoeffs[17];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[181] * p_rints_x_ecoeffs[6];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[180] * p_rints_x_ecoeffs[5];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[177] * p_rints_x_ecoeffs[2];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[188] * p_rints_x_ecoeffs[13];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[182] * p_rints_x_ecoeffs[7];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[185] * p_rints_x_ecoeffs[10];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[175] * p_rints_x_ecoeffs[0];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[176] * p_rints_x_ecoeffs[1];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[187] * p_rints_x_ecoeffs[12];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[179] * p_rints_x_ecoeffs[4];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[197] * p_rints_x_ecoeffs[22];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[178] * p_rints_x_ecoeffs[3];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[189] * p_rints_x_ecoeffs[14];
            eri3_batch(5, 0, 0) += p_ecoeffs_ab[183] * p_rints_x_ecoeffs[8];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[234] * p_rints_x_ecoeffs[24];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[216] * p_rints_x_ecoeffs[6];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[227] * p_rints_x_ecoeffs[17];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[215] * p_rints_x_ecoeffs[5];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[241] * p_rints_x_ecoeffs[31];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[212] * p_rints_x_ecoeffs[2];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[221] * p_rints_x_ecoeffs[11];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[217] * p_rints_x_ecoeffs[7];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[210] * p_rints_x_ecoeffs[0];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[226] * p_rints_x_ecoeffs[16];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[211] * p_rints_x_ecoeffs[1];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[222] * p_rints_x_ecoeffs[12];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[214] * p_rints_x_ecoeffs[4];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[213] * p_rints_x_ecoeffs[3];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[224] * p_rints_x_ecoeffs[14];
            eri3_batch(6, 0, 0) += p_ecoeffs_ab[218] * p_rints_x_ecoeffs[8];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[250] * p_rints_x_ecoeffs[5];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[247] * p_rints_x_ecoeffs[2];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[258] * p_rints_x_ecoeffs[13];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[252] * p_rints_x_ecoeffs[7];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[256] * p_rints_x_ecoeffs[11];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[245] * p_rints_x_ecoeffs[0];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[261] * p_rints_x_ecoeffs[16];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[255] * p_rints_x_ecoeffs[10];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[275] * p_rints_x_ecoeffs[30];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[265] * p_rints_x_ecoeffs[20];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[246] * p_rints_x_ecoeffs[1];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[249] * p_rints_x_ecoeffs[4];
            eri3_batch(7, 0, 0) += p_ecoeffs_ab[268] * p_rints_x_ecoeffs[23];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[285] * p_rints_x_ecoeffs[5];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[282] * p_rints_x_ecoeffs[2];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[293] * p_rints_x_ecoeffs[13];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[287] * p_rints_x_ecoeffs[7];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[291] * p_rints_x_ecoeffs[11];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[280] * p_rints_x_ecoeffs[0];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[296] * p_rints_x_ecoeffs[16];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[290] * p_rints_x_ecoeffs[10];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[301] * p_rints_x_ecoeffs[21];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[281] * p_rints_x_ecoeffs[1];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[284] * p_rints_x_ecoeffs[4];
            eri3_batch(8, 0, 0) += p_ecoeffs_ab[306] * p_rints_x_ecoeffs[26];
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
            for (int ka = 0; ka < n_sph_c; ka++)
            {
                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                double norm_c = sh_data_c.norms[ofs_norm_c + ka];
                eri3_batch(mu, nu, ka) *= norm_a * norm_b * norm_c;
            }

    return eri3_batch;
}

template<> lible::vec2d
lible::ints::two::eri2Kernel<4, 0>(const int ishell_a, const int ishell_b,
                                   const std::vector<double> &ecoeffs_a,
                                   const std::vector<double> &ecoeffs_b_tsp,
                                   const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    const int cdepth_a = sh_data_a.cdepths[ishell_a];
    const int cdepth_b = sh_data_b.cdepths[ishell_b];
    const int cofs_a = sh_data_a.coffsets[ishell_a];
    const int cofs_b = sh_data_b.coffsets[ishell_b];

    const double *exps_a = &sh_data_a.exps[cofs_a];
    const double *exps_b = &sh_data_b.exps[cofs_b];
    const double *coords_a = &sh_data_a.coords[3 * ishell_a];
    const double *coords_b = &sh_data_b.coords[3 * ishell_b];
    const double *pecoeffs_a = &ecoeffs_a[sh_data_a.offsets_ecoeffs[ishell_a]];
    const double *pecoeffs_b_tsp = &ecoeffs_b_tsp[sh_data_b.offsets_ecoeffs[ishell_b]];

    constexpr int la = 4, lb = 0;
    constexpr int lab = la + lb;
    constexpr int n_sph_a = numSphericalsC(la);
    constexpr int n_sph_b = numSphericalsC(lb);
    constexpr int n_hermite_a = numHermitesC(la);
    constexpr int n_hermite_b = numHermitesC(lb);
    constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
    constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

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

            const double* p_ecoeffs_b = &pecoeffs_b_tsp[ib * n_ecoeffs_b];
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
            p_rints_x_ecoeffs[20] += rints[20] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[21] += rints[21] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[22] += rints[22] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[23] += rints[23] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[24] += rints[24] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[25] += rints[25] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[26] += rints[26] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[27] += rints[27] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[28] += rints[28] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[29] += rints[29] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[30] += rints[30] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[31] += rints[31] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[32] += rints[32] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[33] += rints[33] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[34] += rints[34] * p_ecoeffs_b[0];
        }
    }

    vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);
    for (int ia = 0; ia < cdepth_a; ia++)
    {
        const double* p_ecoeffs_a = &pecoeffs_a[ia * n_ecoeffs_a];
        const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];

        eri2_batch(0, 0) += p_ecoeffs_a[25] * p_rints_x_ecoeffs[25];
        eri2_batch(0, 0) += p_ecoeffs_a[2] * p_rints_x_ecoeffs[2];
        eri2_batch(0, 0) += p_ecoeffs_a[11] * p_rints_x_ecoeffs[11];
        eri2_batch(0, 0) += p_ecoeffs_a[16] * p_rints_x_ecoeffs[16];
        eri2_batch(0, 0) += p_ecoeffs_a[3] * p_rints_x_ecoeffs[3];
        eri2_batch(0, 0) += p_ecoeffs_a[34] * p_rints_x_ecoeffs[34];
        eri2_batch(0, 0) += p_ecoeffs_a[17] * p_rints_x_ecoeffs[17];
        eri2_batch(0, 0) += p_ecoeffs_a[6] * p_rints_x_ecoeffs[6];
        eri2_batch(0, 0) += p_ecoeffs_a[5] * p_rints_x_ecoeffs[5];
        eri2_batch(0, 0) += p_ecoeffs_a[10] * p_rints_x_ecoeffs[10];
        eri2_batch(0, 0) += p_ecoeffs_a[20] * p_rints_x_ecoeffs[20];
        eri2_batch(0, 0) += p_ecoeffs_a[12] * p_rints_x_ecoeffs[12];
        eri2_batch(0, 0) += p_ecoeffs_a[18] * p_rints_x_ecoeffs[18];
        eri2_batch(0, 0) += p_ecoeffs_a[19] * p_rints_x_ecoeffs[19];
        eri2_batch(0, 0) += p_ecoeffs_a[7] * p_rints_x_ecoeffs[7];
        eri2_batch(0, 0) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[0];
        eri2_batch(0, 0) += p_ecoeffs_a[30] * p_rints_x_ecoeffs[30];
        eri2_batch(0, 0) += p_ecoeffs_a[1] * p_rints_x_ecoeffs[1];
        eri2_batch(0, 0) += p_ecoeffs_a[4] * p_rints_x_ecoeffs[4];
        eri2_batch(0, 0) += p_ecoeffs_a[23] * p_rints_x_ecoeffs[23];
        eri2_batch(0, 0) += p_ecoeffs_a[8] * p_rints_x_ecoeffs[8];
        eri2_batch(0, 0) += p_ecoeffs_a[13] * p_rints_x_ecoeffs[13];
        eri2_batch(0, 0) += p_ecoeffs_a[9] * p_rints_x_ecoeffs[9];
        eri2_batch(0, 0) += p_ecoeffs_a[32] * p_rints_x_ecoeffs[32];
        eri2_batch(0, 0) += p_ecoeffs_a[15] * p_rints_x_ecoeffs[15];
        eri2_batch(1, 0) += p_ecoeffs_a[37] * p_rints_x_ecoeffs[2];
        eri2_batch(1, 0) += p_ecoeffs_a[38] * p_rints_x_ecoeffs[3];
        eri2_batch(1, 0) += p_ecoeffs_a[52] * p_rints_x_ecoeffs[17];
        eri2_batch(1, 0) += p_ecoeffs_a[41] * p_rints_x_ecoeffs[6];
        eri2_batch(1, 0) += p_ecoeffs_a[40] * p_rints_x_ecoeffs[5];
        eri2_batch(1, 0) += p_ecoeffs_a[45] * p_rints_x_ecoeffs[10];
        eri2_batch(1, 0) += p_ecoeffs_a[47] * p_rints_x_ecoeffs[12];
        eri2_batch(1, 0) += p_ecoeffs_a[62] * p_rints_x_ecoeffs[27];
        eri2_batch(1, 0) += p_ecoeffs_a[54] * p_rints_x_ecoeffs[19];
        eri2_batch(1, 0) += p_ecoeffs_a[42] * p_rints_x_ecoeffs[7];
        eri2_batch(1, 0) += p_ecoeffs_a[35] * p_rints_x_ecoeffs[0];
        eri2_batch(1, 0) += p_ecoeffs_a[36] * p_rints_x_ecoeffs[1];
        eri2_batch(1, 0) += p_ecoeffs_a[64] * p_rints_x_ecoeffs[29];
        eri2_batch(1, 0) += p_ecoeffs_a[39] * p_rints_x_ecoeffs[4];
        eri2_batch(1, 0) += p_ecoeffs_a[43] * p_rints_x_ecoeffs[8];
        eri2_batch(1, 0) += p_ecoeffs_a[48] * p_rints_x_ecoeffs[13];
        eri2_batch(1, 0) += p_ecoeffs_a[44] * p_rints_x_ecoeffs[9];
        eri2_batch(1, 0) += p_ecoeffs_a[57] * p_rints_x_ecoeffs[22];
        eri2_batch(1, 0) += p_ecoeffs_a[50] * p_rints_x_ecoeffs[15];
        eri2_batch(1, 0) += p_ecoeffs_a[49] * p_rints_x_ecoeffs[14];
        eri2_batch(2, 0) += p_ecoeffs_a[72] * p_rints_x_ecoeffs[2];
        eri2_batch(2, 0) += p_ecoeffs_a[81] * p_rints_x_ecoeffs[11];
        eri2_batch(2, 0) += p_ecoeffs_a[103] * p_rints_x_ecoeffs[33];
        eri2_batch(2, 0) += p_ecoeffs_a[86] * p_rints_x_ecoeffs[16];
        eri2_batch(2, 0) += p_ecoeffs_a[73] * p_rints_x_ecoeffs[3];
        eri2_batch(2, 0) += p_ecoeffs_a[87] * p_rints_x_ecoeffs[17];
        eri2_batch(2, 0) += p_ecoeffs_a[76] * p_rints_x_ecoeffs[6];
        eri2_batch(2, 0) += p_ecoeffs_a[75] * p_rints_x_ecoeffs[5];
        eri2_batch(2, 0) += p_ecoeffs_a[82] * p_rints_x_ecoeffs[12];
        eri2_batch(2, 0) += p_ecoeffs_a[88] * p_rints_x_ecoeffs[18];
        eri2_batch(2, 0) += p_ecoeffs_a[89] * p_rints_x_ecoeffs[19];
        eri2_batch(2, 0) += p_ecoeffs_a[77] * p_rints_x_ecoeffs[7];
        eri2_batch(2, 0) += p_ecoeffs_a[70] * p_rints_x_ecoeffs[0];
        eri2_batch(2, 0) += p_ecoeffs_a[71] * p_rints_x_ecoeffs[1];
        eri2_batch(2, 0) += p_ecoeffs_a[74] * p_rints_x_ecoeffs[4];
        eri2_batch(2, 0) += p_ecoeffs_a[78] * p_rints_x_ecoeffs[8];
        eri2_batch(2, 0) += p_ecoeffs_a[94] * p_rints_x_ecoeffs[24];
        eri2_batch(2, 0) += p_ecoeffs_a[101] * p_rints_x_ecoeffs[31];
        eri2_batch(2, 0) += p_ecoeffs_a[79] * p_rints_x_ecoeffs[9];
        eri2_batch(2, 0) += p_ecoeffs_a[84] * p_rints_x_ecoeffs[14];
        eri2_batch(3, 0) += p_ecoeffs_a[130] * p_rints_x_ecoeffs[25];
        eri2_batch(3, 0) += p_ecoeffs_a[107] * p_rints_x_ecoeffs[2];
        eri2_batch(3, 0) += p_ecoeffs_a[116] * p_rints_x_ecoeffs[11];
        eri2_batch(3, 0) += p_ecoeffs_a[121] * p_rints_x_ecoeffs[16];
        eri2_batch(3, 0) += p_ecoeffs_a[108] * p_rints_x_ecoeffs[3];
        eri2_batch(3, 0) += p_ecoeffs_a[122] * p_rints_x_ecoeffs[17];
        eri2_batch(3, 0) += p_ecoeffs_a[111] * p_rints_x_ecoeffs[6];
        eri2_batch(3, 0) += p_ecoeffs_a[110] * p_rints_x_ecoeffs[5];
        eri2_batch(3, 0) += p_ecoeffs_a[115] * p_rints_x_ecoeffs[10];
        eri2_batch(3, 0) += p_ecoeffs_a[125] * p_rints_x_ecoeffs[20];
        eri2_batch(3, 0) += p_ecoeffs_a[117] * p_rints_x_ecoeffs[12];
        eri2_batch(3, 0) += p_ecoeffs_a[123] * p_rints_x_ecoeffs[18];
        eri2_batch(3, 0) += p_ecoeffs_a[112] * p_rints_x_ecoeffs[7];
        eri2_batch(3, 0) += p_ecoeffs_a[105] * p_rints_x_ecoeffs[0];
        eri2_batch(3, 0) += p_ecoeffs_a[135] * p_rints_x_ecoeffs[30];
        eri2_batch(3, 0) += p_ecoeffs_a[106] * p_rints_x_ecoeffs[1];
        eri2_batch(3, 0) += p_ecoeffs_a[109] * p_rints_x_ecoeffs[4];
        eri2_batch(3, 0) += p_ecoeffs_a[128] * p_rints_x_ecoeffs[23];
        eri2_batch(3, 0) += p_ecoeffs_a[113] * p_rints_x_ecoeffs[8];
        eri2_batch(3, 0) += p_ecoeffs_a[118] * p_rints_x_ecoeffs[13];
        eri2_batch(3, 0) += p_ecoeffs_a[114] * p_rints_x_ecoeffs[9];
        eri2_batch(3, 0) += p_ecoeffs_a[137] * p_rints_x_ecoeffs[32];
        eri2_batch(3, 0) += p_ecoeffs_a[120] * p_rints_x_ecoeffs[15];
        eri2_batch(4, 0) += p_ecoeffs_a[142] * p_rints_x_ecoeffs[2];
        eri2_batch(4, 0) += p_ecoeffs_a[151] * p_rints_x_ecoeffs[11];
        eri2_batch(4, 0) += p_ecoeffs_a[156] * p_rints_x_ecoeffs[16];
        eri2_batch(4, 0) += p_ecoeffs_a[166] * p_rints_x_ecoeffs[26];
        eri2_batch(4, 0) += p_ecoeffs_a[143] * p_rints_x_ecoeffs[3];
        eri2_batch(4, 0) += p_ecoeffs_a[146] * p_rints_x_ecoeffs[6];
        eri2_batch(4, 0) += p_ecoeffs_a[145] * p_rints_x_ecoeffs[5];
        eri2_batch(4, 0) += p_ecoeffs_a[150] * p_rints_x_ecoeffs[10];
        eri2_batch(4, 0) += p_ecoeffs_a[158] * p_rints_x_ecoeffs[18];
        eri2_batch(4, 0) += p_ecoeffs_a[147] * p_rints_x_ecoeffs[7];
        eri2_batch(4, 0) += p_ecoeffs_a[140] * p_rints_x_ecoeffs[0];
        eri2_batch(4, 0) += p_ecoeffs_a[168] * p_rints_x_ecoeffs[28];
        eri2_batch(4, 0) += p_ecoeffs_a[141] * p_rints_x_ecoeffs[1];
        eri2_batch(4, 0) += p_ecoeffs_a[144] * p_rints_x_ecoeffs[4];
        eri2_batch(4, 0) += p_ecoeffs_a[148] * p_rints_x_ecoeffs[8];
        eri2_batch(4, 0) += p_ecoeffs_a[153] * p_rints_x_ecoeffs[13];
        eri2_batch(4, 0) += p_ecoeffs_a[161] * p_rints_x_ecoeffs[21];
        eri2_batch(4, 0) += p_ecoeffs_a[149] * p_rints_x_ecoeffs[9];
        eri2_batch(4, 0) += p_ecoeffs_a[155] * p_rints_x_ecoeffs[15];
        eri2_batch(4, 0) += p_ecoeffs_a[154] * p_rints_x_ecoeffs[14];
        eri2_batch(5, 0) += p_ecoeffs_a[202] * p_rints_x_ecoeffs[27];
        eri2_batch(5, 0) += p_ecoeffs_a[192] * p_rints_x_ecoeffs[17];
        eri2_batch(5, 0) += p_ecoeffs_a[181] * p_rints_x_ecoeffs[6];
        eri2_batch(5, 0) += p_ecoeffs_a[180] * p_rints_x_ecoeffs[5];
        eri2_batch(5, 0) += p_ecoeffs_a[177] * p_rints_x_ecoeffs[2];
        eri2_batch(5, 0) += p_ecoeffs_a[188] * p_rints_x_ecoeffs[13];
        eri2_batch(5, 0) += p_ecoeffs_a[182] * p_rints_x_ecoeffs[7];
        eri2_batch(5, 0) += p_ecoeffs_a[185] * p_rints_x_ecoeffs[10];
        eri2_batch(5, 0) += p_ecoeffs_a[175] * p_rints_x_ecoeffs[0];
        eri2_batch(5, 0) += p_ecoeffs_a[176] * p_rints_x_ecoeffs[1];
        eri2_batch(5, 0) += p_ecoeffs_a[187] * p_rints_x_ecoeffs[12];
        eri2_batch(5, 0) += p_ecoeffs_a[179] * p_rints_x_ecoeffs[4];
        eri2_batch(5, 0) += p_ecoeffs_a[197] * p_rints_x_ecoeffs[22];
        eri2_batch(5, 0) += p_ecoeffs_a[178] * p_rints_x_ecoeffs[3];
        eri2_batch(5, 0) += p_ecoeffs_a[189] * p_rints_x_ecoeffs[14];
        eri2_batch(5, 0) += p_ecoeffs_a[183] * p_rints_x_ecoeffs[8];
        eri2_batch(6, 0) += p_ecoeffs_a[234] * p_rints_x_ecoeffs[24];
        eri2_batch(6, 0) += p_ecoeffs_a[216] * p_rints_x_ecoeffs[6];
        eri2_batch(6, 0) += p_ecoeffs_a[227] * p_rints_x_ecoeffs[17];
        eri2_batch(6, 0) += p_ecoeffs_a[215] * p_rints_x_ecoeffs[5];
        eri2_batch(6, 0) += p_ecoeffs_a[241] * p_rints_x_ecoeffs[31];
        eri2_batch(6, 0) += p_ecoeffs_a[212] * p_rints_x_ecoeffs[2];
        eri2_batch(6, 0) += p_ecoeffs_a[221] * p_rints_x_ecoeffs[11];
        eri2_batch(6, 0) += p_ecoeffs_a[217] * p_rints_x_ecoeffs[7];
        eri2_batch(6, 0) += p_ecoeffs_a[210] * p_rints_x_ecoeffs[0];
        eri2_batch(6, 0) += p_ecoeffs_a[226] * p_rints_x_ecoeffs[16];
        eri2_batch(6, 0) += p_ecoeffs_a[211] * p_rints_x_ecoeffs[1];
        eri2_batch(6, 0) += p_ecoeffs_a[222] * p_rints_x_ecoeffs[12];
        eri2_batch(6, 0) += p_ecoeffs_a[214] * p_rints_x_ecoeffs[4];
        eri2_batch(6, 0) += p_ecoeffs_a[213] * p_rints_x_ecoeffs[3];
        eri2_batch(6, 0) += p_ecoeffs_a[224] * p_rints_x_ecoeffs[14];
        eri2_batch(6, 0) += p_ecoeffs_a[218] * p_rints_x_ecoeffs[8];
        eri2_batch(7, 0) += p_ecoeffs_a[250] * p_rints_x_ecoeffs[5];
        eri2_batch(7, 0) += p_ecoeffs_a[247] * p_rints_x_ecoeffs[2];
        eri2_batch(7, 0) += p_ecoeffs_a[258] * p_rints_x_ecoeffs[13];
        eri2_batch(7, 0) += p_ecoeffs_a[252] * p_rints_x_ecoeffs[7];
        eri2_batch(7, 0) += p_ecoeffs_a[256] * p_rints_x_ecoeffs[11];
        eri2_batch(7, 0) += p_ecoeffs_a[245] * p_rints_x_ecoeffs[0];
        eri2_batch(7, 0) += p_ecoeffs_a[261] * p_rints_x_ecoeffs[16];
        eri2_batch(7, 0) += p_ecoeffs_a[255] * p_rints_x_ecoeffs[10];
        eri2_batch(7, 0) += p_ecoeffs_a[275] * p_rints_x_ecoeffs[30];
        eri2_batch(7, 0) += p_ecoeffs_a[265] * p_rints_x_ecoeffs[20];
        eri2_batch(7, 0) += p_ecoeffs_a[246] * p_rints_x_ecoeffs[1];
        eri2_batch(7, 0) += p_ecoeffs_a[249] * p_rints_x_ecoeffs[4];
        eri2_batch(7, 0) += p_ecoeffs_a[268] * p_rints_x_ecoeffs[23];
        eri2_batch(8, 0) += p_ecoeffs_a[285] * p_rints_x_ecoeffs[5];
        eri2_batch(8, 0) += p_ecoeffs_a[282] * p_rints_x_ecoeffs[2];
        eri2_batch(8, 0) += p_ecoeffs_a[293] * p_rints_x_ecoeffs[13];
        eri2_batch(8, 0) += p_ecoeffs_a[287] * p_rints_x_ecoeffs[7];
        eri2_batch(8, 0) += p_ecoeffs_a[291] * p_rints_x_ecoeffs[11];
        eri2_batch(8, 0) += p_ecoeffs_a[280] * p_rints_x_ecoeffs[0];
        eri2_batch(8, 0) += p_ecoeffs_a[296] * p_rints_x_ecoeffs[16];
        eri2_batch(8, 0) += p_ecoeffs_a[290] * p_rints_x_ecoeffs[10];
        eri2_batch(8, 0) += p_ecoeffs_a[301] * p_rints_x_ecoeffs[21];
        eri2_batch(8, 0) += p_ecoeffs_a[281] * p_rints_x_ecoeffs[1];
        eri2_batch(8, 0) += p_ecoeffs_a[284] * p_rints_x_ecoeffs[4];
        eri2_batch(8, 0) += p_ecoeffs_a[306] * p_rints_x_ecoeffs[26];
    }

    int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];
    int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];
    for (int mu = 0; mu < n_sph_a; mu++)
        for (int nu = 0; nu < n_sph_b; nu++)
        {
            double norm_a = sh_data_a.norms[ofs_norm_a + mu];
            double norm_b = sh_data_b.norms[ofs_norm_b + nu];
            eri2_batch(mu, nu) *= norm_a * norm_b;
        }

    return eri2_batch;
}
