#include <lible/ints/twoel/eri_kernels.hpp>

template<> lible::vec4d
lible::ints::two::eri4Kernel<1, 1, 1, 0>(const int ipair_ab, const int ipair_cd,
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

    constexpr int la = 1, lb = 1, lc = 1, ld = 0;
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
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[12] += rints[16] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[12] += rints[19] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[14] += rints[16] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[14] += rints[18] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[16] += rints[21] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[17] += rints[20] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[17] += rints[22] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[18] += rints[27] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[19] += rints[25] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[19] += rints[24] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[20] += rints[24] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[21] += rints[28] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[21] += rints[31] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[22] += rints[29] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[22] += rints[28] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[23] += rints[28] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[23] += rints[30] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[24] += rints[32] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[24] += rints[35] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[25] += rints[33] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[25] += rints[32] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[26] += rints[32] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[26] += rints[34] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[27] += rints[36] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[27] += rints[39] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[28] += rints[37] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[28] += rints[36] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[29] += rints[36] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[29] += rints[38] * p_ecoeffs_cd_tsp[8];
                }
        }

    vec4d eri4_batch(n_sph_a, n_sph_b, n_sph_c, n_sph_d, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[27];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[28];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[29];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 1, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 1, 1, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 1, 1, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 1, 1, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 1, 1, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 1, 2, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 1, 2, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 1, 2, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 1, 2, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[20];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 2, 0, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[24];
            eri4_batch(0, 2, 1, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 2, 1, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 2, 1, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 2, 1, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[25];
            eri4_batch(0, 2, 2, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 2, 2, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 2, 2, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 2, 2, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[26];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[18];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[19];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[20];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 1, 0, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[12];
            eri4_batch(1, 1, 1, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 1, 1, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 1, 1, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[13];
            eri4_batch(1, 1, 2, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 1, 2, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 1, 2, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[14];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[15];
            eri4_batch(1, 2, 0, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[6];
            eri4_batch(1, 2, 1, 0) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 2, 1, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 2, 1, 0) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[16];
            eri4_batch(1, 2, 1, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[7];
            eri4_batch(1, 2, 2, 0) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 2, 2, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 2, 2, 0) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[17];
            eri4_batch(1, 2, 2, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[10];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[25];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[26];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[3];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[15];
            eri4_batch(2, 1, 0, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 1, 1, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[4];
            eri4_batch(2, 1, 1, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 1, 1, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[16];
            eri4_batch(2, 1, 1, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 1, 2, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[5];
            eri4_batch(2, 1, 2, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 1, 2, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[17];
            eri4_batch(2, 1, 2, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[21];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 2, 0, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 2, 1, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[22];
            eri4_batch(2, 2, 1, 0) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 2, 1, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 2, 2, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[23];
            eri4_batch(2, 2, 2, 0) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 2, 2, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[8];
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
lible::ints::two::eri4Kernel<2, 0, 1, 0>(const int ipair_ab, const int ipair_cd,
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

    constexpr int la = 2, lb = 0, lc = 1, ld = 0;
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
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[12] += rints[16] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[12] += rints[19] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[14] += rints[16] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[14] += rints[18] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[16] += rints[21] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[17] += rints[20] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[17] += rints[22] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[18] += rints[27] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[19] += rints[25] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[19] += rints[24] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[20] += rints[24] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[21] += rints[28] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[21] += rints[31] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[22] += rints[29] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[22] += rints[28] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[23] += rints[28] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[23] += rints[30] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[24] += rints[32] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[24] += rints[35] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[25] += rints[33] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[25] += rints[32] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[26] += rints[32] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[26] += rints[34] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[27] += rints[36] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[27] += rints[39] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[28] += rints[37] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[28] += rints[36] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[29] += rints[36] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[29] += rints[38] * p_ecoeffs_cd_tsp[8];
                }
        }

    vec4d eri4_batch(n_sph_a, n_sph_b, n_sph_c, n_sph_d, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[21];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[27];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[22];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[28];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[23];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[29];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[3];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[9];
            eri4_batch(1, 0, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[18];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[4];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[1];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[10];
            eri4_batch(1, 0, 1, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[19];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[5];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[2];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[11];
            eri4_batch(1, 0, 2, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[20];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[9];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[6];
            eri4_batch(2, 0, 0, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[24];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[1];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[10];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[7];
            eri4_batch(2, 0, 1, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[25];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[2];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[11];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[8];
            eri4_batch(2, 0, 2, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[26];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[6];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[21];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[3];
            eri4_batch(3, 0, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[12];
            eri4_batch(3, 0, 1, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[7];
            eri4_batch(3, 0, 1, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[22];
            eri4_batch(3, 0, 1, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[1];
            eri4_batch(3, 0, 1, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[4];
            eri4_batch(3, 0, 1, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[13];
            eri4_batch(3, 0, 2, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[8];
            eri4_batch(3, 0, 2, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[23];
            eri4_batch(3, 0, 2, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[2];
            eri4_batch(3, 0, 2, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[5];
            eri4_batch(3, 0, 2, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[14];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[3];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[15];
            eri4_batch(4, 0, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[6];
            eri4_batch(4, 0, 1, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[4];
            eri4_batch(4, 0, 1, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[1];
            eri4_batch(4, 0, 1, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[16];
            eri4_batch(4, 0, 1, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri4_batch(4, 0, 2, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[5];
            eri4_batch(4, 0, 2, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[2];
            eri4_batch(4, 0, 2, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[17];
            eri4_batch(4, 0, 2, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[8];
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
lible::ints::two::eri3Kernel<1, 1, 1>(const int ipair_ab, const int ishell_c,
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

    constexpr int la = 1, lb = 1, lc = 1;
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
                p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[12] += rints[16] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[12] += rints[19] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[14] += rints[16] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[14] += rints[18] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[16] += rints[21] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[17] += rints[20] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[17] += rints[22] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[18] += rints[27] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[19] += rints[25] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[19] += rints[24] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[20] += rints[24] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[21] += rints[28] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[21] += rints[31] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[22] += rints[29] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[22] += rints[28] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[23] += rints[28] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[23] += rints[30] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[24] += rints[32] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[24] += rints[35] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[25] += rints[33] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[25] += rints[32] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[26] += rints[32] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[26] += rints[34] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[27] += rints[36] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[27] += rints[39] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[28] += rints[37] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[28] += rints[36] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[29] += rints[36] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[29] += rints[38] * p_ecoeffs_c[10];
            }
}

    vec3d eri3_batch(n_sph_a, n_sph_b, n_sph_c, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[27];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[28];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[29];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 1, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[18];
            eri3_batch(0, 1, 1) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 1, 1) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 1, 1) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 1, 1) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[19];
            eri3_batch(0, 1, 2) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 1, 2) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 1, 2) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 1, 2) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[20];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 2, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[24];
            eri3_batch(0, 2, 1) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 2, 1) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 2, 1) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 2, 1) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[25];
            eri3_batch(0, 2, 2) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 2, 2) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[11];
            eri3_batch(0, 2, 2) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 2, 2) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[26];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[18];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[19];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[33] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[36] * p_rints_x_ecoeffs[20];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 1, 0) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[12];
            eri3_batch(1, 1, 1) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 1, 1) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 1, 1) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[13];
            eri3_batch(1, 1, 2) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 1, 2) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 1, 2) += p_ecoeffs_ab[44] * p_rints_x_ecoeffs[14];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[15];
            eri3_batch(1, 2, 0) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[6];
            eri3_batch(1, 2, 1) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 2, 1) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 2, 1) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[16];
            eri3_batch(1, 2, 1) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[7];
            eri3_batch(1, 2, 2) += p_ecoeffs_ab[51] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 2, 2) += p_ecoeffs_ab[50] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 2, 2) += p_ecoeffs_ab[55] * p_rints_x_ecoeffs[17];
            eri3_batch(1, 2, 2) += p_ecoeffs_ab[52] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[10];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[25];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[60] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[63] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[62] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[68] * p_rints_x_ecoeffs[26];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[3];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[15];
            eri3_batch(2, 1, 0) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 1, 1) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[4];
            eri3_batch(2, 1, 1) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 1, 1) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[16];
            eri3_batch(2, 1, 1) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 1, 2) += p_ecoeffs_ab[71] * p_rints_x_ecoeffs[5];
            eri3_batch(2, 1, 2) += p_ecoeffs_ab[70] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 1, 2) += p_ecoeffs_ab[75] * p_rints_x_ecoeffs[17];
            eri3_batch(2, 1, 2) += p_ecoeffs_ab[72] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[21];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 2, 0) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 2, 1) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[22];
            eri3_batch(2, 2, 1) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 2, 1) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 2, 2) += p_ecoeffs_ab[87] * p_rints_x_ecoeffs[23];
            eri3_batch(2, 2, 2) += p_ecoeffs_ab[80] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 2, 2) += p_ecoeffs_ab[82] * p_rints_x_ecoeffs[8];
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
lible::ints::two::eri3Kernel<2, 0, 1>(const int ipair_ab, const int ishell_c,
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

    constexpr int la = 2, lb = 0, lc = 1;
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
                p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[12] += rints[16] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[12] += rints[19] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[14] += rints[16] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[14] += rints[18] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[16] += rints[21] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[17] += rints[20] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[17] += rints[22] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[18] += rints[27] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[19] += rints[25] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[19] += rints[24] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[20] += rints[24] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[21] += rints[28] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[21] += rints[31] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[22] += rints[29] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[22] += rints[28] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[23] += rints[28] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[23] += rints[30] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[24] += rints[32] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[24] += rints[35] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[25] += rints[33] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[25] += rints[32] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[26] += rints[32] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[26] += rints[34] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[27] += rints[36] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[27] += rints[39] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[28] += rints[37] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[28] += rints[36] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[29] += rints[36] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[29] += rints[38] * p_ecoeffs_c[10];
            }
}

    vec3d eri3_batch(n_sph_a, n_sph_b, n_sph_c, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[21];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[27];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[12];
            eri3_batch(0, 0, 0) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[9];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[22];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[28];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[13];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[10];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[2] * p_rints_x_ecoeffs[8];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[7] * p_rints_x_ecoeffs[23];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[1] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[9] * p_rints_x_ecoeffs[29];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[4] * p_rints_x_ecoeffs[14];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[3] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[3];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[0];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[9];
            eri3_batch(1, 0, 0) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[18];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[4];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[1];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[10];
            eri3_batch(1, 0, 1) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[19];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[11] * p_rints_x_ecoeffs[5];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[10] * p_rints_x_ecoeffs[2];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[13] * p_rints_x_ecoeffs[11];
            eri3_batch(1, 0, 2) += p_ecoeffs_ab[16] * p_rints_x_ecoeffs[20];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[0];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[9];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[6];
            eri3_batch(2, 0, 0) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[24];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[1];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[10];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[7];
            eri3_batch(2, 0, 1) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[25];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[20] * p_rints_x_ecoeffs[2];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[23] * p_rints_x_ecoeffs[11];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[22] * p_rints_x_ecoeffs[8];
            eri3_batch(2, 0, 2) += p_ecoeffs_ab[28] * p_rints_x_ecoeffs[26];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[6];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[21];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[0];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[3];
            eri3_batch(3, 0, 0) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[12];
            eri3_batch(3, 0, 1) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[7];
            eri3_batch(3, 0, 1) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[22];
            eri3_batch(3, 0, 1) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[1];
            eri3_batch(3, 0, 1) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[4];
            eri3_batch(3, 0, 1) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[13];
            eri3_batch(3, 0, 2) += p_ecoeffs_ab[32] * p_rints_x_ecoeffs[8];
            eri3_batch(3, 0, 2) += p_ecoeffs_ab[37] * p_rints_x_ecoeffs[23];
            eri3_batch(3, 0, 2) += p_ecoeffs_ab[30] * p_rints_x_ecoeffs[2];
            eri3_batch(3, 0, 2) += p_ecoeffs_ab[31] * p_rints_x_ecoeffs[5];
            eri3_batch(3, 0, 2) += p_ecoeffs_ab[34] * p_rints_x_ecoeffs[14];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[3];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[0];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[15];
            eri3_batch(4, 0, 0) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[6];
            eri3_batch(4, 0, 1) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[4];
            eri3_batch(4, 0, 1) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[1];
            eri3_batch(4, 0, 1) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[16];
            eri3_batch(4, 0, 1) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[7];
            eri3_batch(4, 0, 2) += p_ecoeffs_ab[41] * p_rints_x_ecoeffs[5];
            eri3_batch(4, 0, 2) += p_ecoeffs_ab[40] * p_rints_x_ecoeffs[2];
            eri3_batch(4, 0, 2) += p_ecoeffs_ab[45] * p_rints_x_ecoeffs[17];
            eri3_batch(4, 0, 2) += p_ecoeffs_ab[42] * p_rints_x_ecoeffs[8];
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
lible::ints::two::eri2Kernel<2, 1>(const int ishell_a, const int ishell_b,
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

    constexpr int la = 2, lb = 1;
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
            p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[12] += rints[16] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[12] += rints[19] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[14] += rints[16] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[14] += rints[18] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[16] += rints[21] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[17] += rints[20] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[17] += rints[22] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[18] += rints[27] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[19] += rints[25] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[19] += rints[24] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[20] += rints[24] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[21] += rints[28] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[21] += rints[31] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[22] += rints[29] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[22] += rints[28] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[23] += rints[28] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[23] += rints[30] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[24] += rints[32] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[24] += rints[35] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[25] += rints[33] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[25] += rints[32] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[26] += rints[32] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[26] += rints[34] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[27] += rints[36] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[27] += rints[39] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[28] += rints[37] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[28] += rints[36] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[29] += rints[36] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[29] += rints[38] * p_ecoeffs_b[8];
        }
    }

    vec2d eri2_batch(n_sph_a, n_sph_b, 0);
    for (int ia = 0; ia < cdepth_a; ia++)
    {
        const double* p_ecoeffs_a = &pecoeffs_a[ia * n_ecoeffs_a];
        const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];

        eri2_batch(0, 0) += p_ecoeffs_a[2] * p_rints_x_ecoeffs[6];
        eri2_batch(0, 0) += p_ecoeffs_a[7] * p_rints_x_ecoeffs[21];
        eri2_batch(0, 0) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[0];
        eri2_batch(0, 0) += p_ecoeffs_a[1] * p_rints_x_ecoeffs[3];
        eri2_batch(0, 0) += p_ecoeffs_a[9] * p_rints_x_ecoeffs[27];
        eri2_batch(0, 0) += p_ecoeffs_a[4] * p_rints_x_ecoeffs[12];
        eri2_batch(0, 0) += p_ecoeffs_a[3] * p_rints_x_ecoeffs[9];
        eri2_batch(0, 1) += p_ecoeffs_a[2] * p_rints_x_ecoeffs[7];
        eri2_batch(0, 1) += p_ecoeffs_a[7] * p_rints_x_ecoeffs[22];
        eri2_batch(0, 1) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[1];
        eri2_batch(0, 1) += p_ecoeffs_a[1] * p_rints_x_ecoeffs[4];
        eri2_batch(0, 1) += p_ecoeffs_a[9] * p_rints_x_ecoeffs[28];
        eri2_batch(0, 1) += p_ecoeffs_a[4] * p_rints_x_ecoeffs[13];
        eri2_batch(0, 1) += p_ecoeffs_a[3] * p_rints_x_ecoeffs[10];
        eri2_batch(0, 2) += p_ecoeffs_a[2] * p_rints_x_ecoeffs[8];
        eri2_batch(0, 2) += p_ecoeffs_a[7] * p_rints_x_ecoeffs[23];
        eri2_batch(0, 2) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[2];
        eri2_batch(0, 2) += p_ecoeffs_a[1] * p_rints_x_ecoeffs[5];
        eri2_batch(0, 2) += p_ecoeffs_a[9] * p_rints_x_ecoeffs[29];
        eri2_batch(0, 2) += p_ecoeffs_a[4] * p_rints_x_ecoeffs[14];
        eri2_batch(0, 2) += p_ecoeffs_a[3] * p_rints_x_ecoeffs[11];
        eri2_batch(1, 0) += p_ecoeffs_a[11] * p_rints_x_ecoeffs[3];
        eri2_batch(1, 0) += p_ecoeffs_a[10] * p_rints_x_ecoeffs[0];
        eri2_batch(1, 0) += p_ecoeffs_a[13] * p_rints_x_ecoeffs[9];
        eri2_batch(1, 0) += p_ecoeffs_a[16] * p_rints_x_ecoeffs[18];
        eri2_batch(1, 1) += p_ecoeffs_a[11] * p_rints_x_ecoeffs[4];
        eri2_batch(1, 1) += p_ecoeffs_a[10] * p_rints_x_ecoeffs[1];
        eri2_batch(1, 1) += p_ecoeffs_a[13] * p_rints_x_ecoeffs[10];
        eri2_batch(1, 1) += p_ecoeffs_a[16] * p_rints_x_ecoeffs[19];
        eri2_batch(1, 2) += p_ecoeffs_a[11] * p_rints_x_ecoeffs[5];
        eri2_batch(1, 2) += p_ecoeffs_a[10] * p_rints_x_ecoeffs[2];
        eri2_batch(1, 2) += p_ecoeffs_a[13] * p_rints_x_ecoeffs[11];
        eri2_batch(1, 2) += p_ecoeffs_a[16] * p_rints_x_ecoeffs[20];
        eri2_batch(2, 0) += p_ecoeffs_a[20] * p_rints_x_ecoeffs[0];
        eri2_batch(2, 0) += p_ecoeffs_a[23] * p_rints_x_ecoeffs[9];
        eri2_batch(2, 0) += p_ecoeffs_a[22] * p_rints_x_ecoeffs[6];
        eri2_batch(2, 0) += p_ecoeffs_a[28] * p_rints_x_ecoeffs[24];
        eri2_batch(2, 1) += p_ecoeffs_a[20] * p_rints_x_ecoeffs[1];
        eri2_batch(2, 1) += p_ecoeffs_a[23] * p_rints_x_ecoeffs[10];
        eri2_batch(2, 1) += p_ecoeffs_a[22] * p_rints_x_ecoeffs[7];
        eri2_batch(2, 1) += p_ecoeffs_a[28] * p_rints_x_ecoeffs[25];
        eri2_batch(2, 2) += p_ecoeffs_a[20] * p_rints_x_ecoeffs[2];
        eri2_batch(2, 2) += p_ecoeffs_a[23] * p_rints_x_ecoeffs[11];
        eri2_batch(2, 2) += p_ecoeffs_a[22] * p_rints_x_ecoeffs[8];
        eri2_batch(2, 2) += p_ecoeffs_a[28] * p_rints_x_ecoeffs[26];
        eri2_batch(3, 0) += p_ecoeffs_a[32] * p_rints_x_ecoeffs[6];
        eri2_batch(3, 0) += p_ecoeffs_a[37] * p_rints_x_ecoeffs[21];
        eri2_batch(3, 0) += p_ecoeffs_a[30] * p_rints_x_ecoeffs[0];
        eri2_batch(3, 0) += p_ecoeffs_a[31] * p_rints_x_ecoeffs[3];
        eri2_batch(3, 0) += p_ecoeffs_a[34] * p_rints_x_ecoeffs[12];
        eri2_batch(3, 1) += p_ecoeffs_a[32] * p_rints_x_ecoeffs[7];
        eri2_batch(3, 1) += p_ecoeffs_a[37] * p_rints_x_ecoeffs[22];
        eri2_batch(3, 1) += p_ecoeffs_a[30] * p_rints_x_ecoeffs[1];
        eri2_batch(3, 1) += p_ecoeffs_a[31] * p_rints_x_ecoeffs[4];
        eri2_batch(3, 1) += p_ecoeffs_a[34] * p_rints_x_ecoeffs[13];
        eri2_batch(3, 2) += p_ecoeffs_a[32] * p_rints_x_ecoeffs[8];
        eri2_batch(3, 2) += p_ecoeffs_a[37] * p_rints_x_ecoeffs[23];
        eri2_batch(3, 2) += p_ecoeffs_a[30] * p_rints_x_ecoeffs[2];
        eri2_batch(3, 2) += p_ecoeffs_a[31] * p_rints_x_ecoeffs[5];
        eri2_batch(3, 2) += p_ecoeffs_a[34] * p_rints_x_ecoeffs[14];
        eri2_batch(4, 0) += p_ecoeffs_a[41] * p_rints_x_ecoeffs[3];
        eri2_batch(4, 0) += p_ecoeffs_a[40] * p_rints_x_ecoeffs[0];
        eri2_batch(4, 0) += p_ecoeffs_a[45] * p_rints_x_ecoeffs[15];
        eri2_batch(4, 0) += p_ecoeffs_a[42] * p_rints_x_ecoeffs[6];
        eri2_batch(4, 1) += p_ecoeffs_a[41] * p_rints_x_ecoeffs[4];
        eri2_batch(4, 1) += p_ecoeffs_a[40] * p_rints_x_ecoeffs[1];
        eri2_batch(4, 1) += p_ecoeffs_a[45] * p_rints_x_ecoeffs[16];
        eri2_batch(4, 1) += p_ecoeffs_a[42] * p_rints_x_ecoeffs[7];
        eri2_batch(4, 2) += p_ecoeffs_a[41] * p_rints_x_ecoeffs[5];
        eri2_batch(4, 2) += p_ecoeffs_a[40] * p_rints_x_ecoeffs[2];
        eri2_batch(4, 2) += p_ecoeffs_a[45] * p_rints_x_ecoeffs[17];
        eri2_batch(4, 2) += p_ecoeffs_a[42] * p_rints_x_ecoeffs[8];
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
