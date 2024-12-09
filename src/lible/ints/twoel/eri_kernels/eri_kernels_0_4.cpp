#include <lible/ints/twoel/eri_kernels.hpp>

template<> void lible::ints::two::eri4Kernel<0, 0, 2, 2>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 0, lb = 0, lc = 2, ld = 2;
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

                    p_rints_x_ecoeffs[0] += rints[25] * p_ecoeffs_cd_tsp[625];
                    p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_cd_tsp[50];
                    p_rints_x_ecoeffs[0] += rints[11] * p_ecoeffs_cd_tsp[275];
                    p_rints_x_ecoeffs[0] += rints[16] * p_ecoeffs_cd_tsp[400];
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[75];
                    p_rints_x_ecoeffs[0] += rints[34] * p_ecoeffs_cd_tsp[850];
                    p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_cd_tsp[425];
                    p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_cd_tsp[150];
                    p_rints_x_ecoeffs[0] += rints[5] * p_ecoeffs_cd_tsp[125];
                    p_rints_x_ecoeffs[0] += rints[10] * p_ecoeffs_cd_tsp[250];
                    p_rints_x_ecoeffs[0] += rints[20] * p_ecoeffs_cd_tsp[500];
                    p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_cd_tsp[300];
                    p_rints_x_ecoeffs[0] += rints[18] * p_ecoeffs_cd_tsp[450];
                    p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_cd_tsp[475];
                    p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_cd_tsp[175];
                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[0] += rints[30] * p_ecoeffs_cd_tsp[750];
                    p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_cd_tsp[25];
                    p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_cd_tsp[100];
                    p_rints_x_ecoeffs[0] += rints[23] * p_ecoeffs_cd_tsp[575];
                    p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_cd_tsp[200];
                    p_rints_x_ecoeffs[0] += rints[13] * p_ecoeffs_cd_tsp[325];
                    p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_cd_tsp[225];
                    p_rints_x_ecoeffs[0] += rints[32] * p_ecoeffs_cd_tsp[800];
                    p_rints_x_ecoeffs[0] += rints[15] * p_ecoeffs_cd_tsp[375];
                    p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_cd_tsp[51];
                    p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_cd_tsp[76];
                    p_rints_x_ecoeffs[1] += rints[17] * p_ecoeffs_cd_tsp[426];
                    p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_cd_tsp[151];
                    p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_cd_tsp[126];
                    p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_cd_tsp[251];
                    p_rints_x_ecoeffs[1] += rints[12] * p_ecoeffs_cd_tsp[301];
                    p_rints_x_ecoeffs[1] += rints[27] * p_ecoeffs_cd_tsp[676];
                    p_rints_x_ecoeffs[1] += rints[19] * p_ecoeffs_cd_tsp[476];
                    p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_cd_tsp[176];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[26];
                    p_rints_x_ecoeffs[1] += rints[29] * p_ecoeffs_cd_tsp[726];
                    p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_cd_tsp[101];
                    p_rints_x_ecoeffs[1] += rints[8] * p_ecoeffs_cd_tsp[201];
                    p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_cd_tsp[326];
                    p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_cd_tsp[226];
                    p_rints_x_ecoeffs[1] += rints[22] * p_ecoeffs_cd_tsp[551];
                    p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_cd_tsp[376];
                    p_rints_x_ecoeffs[1] += rints[14] * p_ecoeffs_cd_tsp[351];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[52];
                    p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_cd_tsp[277];
                    p_rints_x_ecoeffs[2] += rints[33] * p_ecoeffs_cd_tsp[827];
                    p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_cd_tsp[402];
                    p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_cd_tsp[77];
                    p_rints_x_ecoeffs[2] += rints[17] * p_ecoeffs_cd_tsp[427];
                    p_rints_x_ecoeffs[2] += rints[6] * p_ecoeffs_cd_tsp[152];
                    p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_cd_tsp[127];
                    p_rints_x_ecoeffs[2] += rints[12] * p_ecoeffs_cd_tsp[302];
                    p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_cd_tsp[452];
                    p_rints_x_ecoeffs[2] += rints[19] * p_ecoeffs_cd_tsp[477];
                    p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_cd_tsp[177];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_cd_tsp[27];
                    p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_cd_tsp[102];
                    p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_cd_tsp[202];
                    p_rints_x_ecoeffs[2] += rints[24] * p_ecoeffs_cd_tsp[602];
                    p_rints_x_ecoeffs[2] += rints[31] * p_ecoeffs_cd_tsp[777];
                    p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_cd_tsp[227];
                    p_rints_x_ecoeffs[2] += rints[14] * p_ecoeffs_cd_tsp[352];
                    p_rints_x_ecoeffs[3] += rints[25] * p_ecoeffs_cd_tsp[628];
                    p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_cd_tsp[53];
                    p_rints_x_ecoeffs[3] += rints[11] * p_ecoeffs_cd_tsp[278];
                    p_rints_x_ecoeffs[3] += rints[16] * p_ecoeffs_cd_tsp[403];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[78];
                    p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_cd_tsp[428];
                    p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_cd_tsp[153];
                    p_rints_x_ecoeffs[3] += rints[5] * p_ecoeffs_cd_tsp[128];
                    p_rints_x_ecoeffs[3] += rints[10] * p_ecoeffs_cd_tsp[253];
                    p_rints_x_ecoeffs[3] += rints[20] * p_ecoeffs_cd_tsp[503];
                    p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_cd_tsp[303];
                    p_rints_x_ecoeffs[3] += rints[18] * p_ecoeffs_cd_tsp[453];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[178];
                    p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_cd_tsp[3];
                    p_rints_x_ecoeffs[3] += rints[30] * p_ecoeffs_cd_tsp[753];
                    p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_cd_tsp[28];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[103];
                    p_rints_x_ecoeffs[3] += rints[23] * p_ecoeffs_cd_tsp[578];
                    p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_cd_tsp[203];
                    p_rints_x_ecoeffs[3] += rints[13] * p_ecoeffs_cd_tsp[328];
                    p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_cd_tsp[228];
                    p_rints_x_ecoeffs[3] += rints[32] * p_ecoeffs_cd_tsp[803];
                    p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_cd_tsp[378];
                    p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_cd_tsp[54];
                    p_rints_x_ecoeffs[4] += rints[11] * p_ecoeffs_cd_tsp[279];
                    p_rints_x_ecoeffs[4] += rints[16] * p_ecoeffs_cd_tsp[404];
                    p_rints_x_ecoeffs[4] += rints[26] * p_ecoeffs_cd_tsp[654];
                    p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_cd_tsp[79];
                    p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_cd_tsp[154];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[129];
                    p_rints_x_ecoeffs[4] += rints[10] * p_ecoeffs_cd_tsp[254];
                    p_rints_x_ecoeffs[4] += rints[18] * p_ecoeffs_cd_tsp[454];
                    p_rints_x_ecoeffs[4] += rints[7] * p_ecoeffs_cd_tsp[179];
                    p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[28] * p_ecoeffs_cd_tsp[704];
                    p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_cd_tsp[29];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[104];
                    p_rints_x_ecoeffs[4] += rints[8] * p_ecoeffs_cd_tsp[204];
                    p_rints_x_ecoeffs[4] += rints[13] * p_ecoeffs_cd_tsp[329];
                    p_rints_x_ecoeffs[4] += rints[21] * p_ecoeffs_cd_tsp[529];
                    p_rints_x_ecoeffs[4] += rints[9] * p_ecoeffs_cd_tsp[229];
                    p_rints_x_ecoeffs[4] += rints[15] * p_ecoeffs_cd_tsp[379];
                    p_rints_x_ecoeffs[4] += rints[14] * p_ecoeffs_cd_tsp[354];
                    p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_cd_tsp[55];
                    p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_cd_tsp[80];
                    p_rints_x_ecoeffs[5] += rints[17] * p_ecoeffs_cd_tsp[430];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[155];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[130];
                    p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_cd_tsp[255];
                    p_rints_x_ecoeffs[5] += rints[12] * p_ecoeffs_cd_tsp[305];
                    p_rints_x_ecoeffs[5] += rints[27] * p_ecoeffs_cd_tsp[680];
                    p_rints_x_ecoeffs[5] += rints[19] * p_ecoeffs_cd_tsp[480];
                    p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_cd_tsp[180];
                    p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_cd_tsp[5];
                    p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_cd_tsp[30];
                    p_rints_x_ecoeffs[5] += rints[29] * p_ecoeffs_cd_tsp[730];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[105];
                    p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_cd_tsp[205];
                    p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_cd_tsp[330];
                    p_rints_x_ecoeffs[5] += rints[9] * p_ecoeffs_cd_tsp[230];
                    p_rints_x_ecoeffs[5] += rints[22] * p_ecoeffs_cd_tsp[555];
                    p_rints_x_ecoeffs[5] += rints[15] * p_ecoeffs_cd_tsp[380];
                    p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_cd_tsp[355];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_cd_tsp[156];
                    p_rints_x_ecoeffs[6] += rints[25] * p_ecoeffs_cd_tsp[631];
                    p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_cd_tsp[6];
                    p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_cd_tsp[31];
                    p_rints_x_ecoeffs[6] += rints[12] * p_ecoeffs_cd_tsp[306];
                    p_rints_x_ecoeffs[6] += rints[9] * p_ecoeffs_cd_tsp[231];
                    p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_cd_tsp[106];
                    p_rints_x_ecoeffs[6] += rints[15] * p_ecoeffs_cd_tsp[381];
                    p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_cd_tsp[81];
                    p_rints_x_ecoeffs[7] += rints[6] * p_ecoeffs_cd_tsp[157];
                    p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_cd_tsp[132];
                    p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_cd_tsp[57];
                    p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_cd_tsp[7];
                    p_rints_x_ecoeffs[7] += rints[28] * p_ecoeffs_cd_tsp[707];
                    p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_cd_tsp[32];
                    p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_cd_tsp[232];
                    p_rints_x_ecoeffs[7] += rints[18] * p_ecoeffs_cd_tsp[457];
                    p_rints_x_ecoeffs[7] += rints[15] * p_ecoeffs_cd_tsp[382];
                    p_rints_x_ecoeffs[7] += rints[3] * p_ecoeffs_cd_tsp[82];
                    p_rints_x_ecoeffs[7] += rints[14] * p_ecoeffs_cd_tsp[357];
                    p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_cd_tsp[207];
                    p_rints_x_ecoeffs[8] += rints[27] * p_ecoeffs_cd_tsp[683];
                    p_rints_x_ecoeffs[8] += rints[17] * p_ecoeffs_cd_tsp[433];
                    p_rints_x_ecoeffs[8] += rints[6] * p_ecoeffs_cd_tsp[158];
                    p_rints_x_ecoeffs[8] += rints[5] * p_ecoeffs_cd_tsp[133];
                    p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_cd_tsp[58];
                    p_rints_x_ecoeffs[8] += rints[13] * p_ecoeffs_cd_tsp[333];
                    p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_cd_tsp[183];
                    p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_cd_tsp[258];
                    p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[8] += rints[1] * p_ecoeffs_cd_tsp[33];
                    p_rints_x_ecoeffs[8] += rints[12] * p_ecoeffs_cd_tsp[308];
                    p_rints_x_ecoeffs[8] += rints[4] * p_ecoeffs_cd_tsp[108];
                    p_rints_x_ecoeffs[8] += rints[22] * p_ecoeffs_cd_tsp[558];
                    p_rints_x_ecoeffs[8] += rints[3] * p_ecoeffs_cd_tsp[83];
                    p_rints_x_ecoeffs[8] += rints[14] * p_ecoeffs_cd_tsp[358];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[208];
                    p_rints_x_ecoeffs[9] += rints[24] * p_ecoeffs_cd_tsp[609];
                    p_rints_x_ecoeffs[9] += rints[6] * p_ecoeffs_cd_tsp[159];
                    p_rints_x_ecoeffs[9] += rints[5] * p_ecoeffs_cd_tsp[134];
                    p_rints_x_ecoeffs[9] += rints[2] * p_ecoeffs_cd_tsp[59];
                    p_rints_x_ecoeffs[9] += rints[11] * p_ecoeffs_cd_tsp[284];
                    p_rints_x_ecoeffs[9] += rints[0] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[9] += rints[1] * p_ecoeffs_cd_tsp[34];
                    p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_cd_tsp[309];
                    p_rints_x_ecoeffs[9] += rints[4] * p_ecoeffs_cd_tsp[109];
                    p_rints_x_ecoeffs[9] += rints[3] * p_ecoeffs_cd_tsp[84];
                    p_rints_x_ecoeffs[9] += rints[14] * p_ecoeffs_cd_tsp[359];
                    p_rints_x_ecoeffs[9] += rints[8] * p_ecoeffs_cd_tsp[209];
                    p_rints_x_ecoeffs[10] += rints[2] * p_ecoeffs_cd_tsp[60];
                    p_rints_x_ecoeffs[10] += rints[11] * p_ecoeffs_cd_tsp[285];
                    p_rints_x_ecoeffs[10] += rints[33] * p_ecoeffs_cd_tsp[835];
                    p_rints_x_ecoeffs[10] += rints[16] * p_ecoeffs_cd_tsp[410];
                    p_rints_x_ecoeffs[10] += rints[3] * p_ecoeffs_cd_tsp[85];
                    p_rints_x_ecoeffs[10] += rints[17] * p_ecoeffs_cd_tsp[435];
                    p_rints_x_ecoeffs[10] += rints[6] * p_ecoeffs_cd_tsp[160];
                    p_rints_x_ecoeffs[10] += rints[5] * p_ecoeffs_cd_tsp[135];
                    p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_cd_tsp[310];
                    p_rints_x_ecoeffs[10] += rints[18] * p_ecoeffs_cd_tsp[460];
                    p_rints_x_ecoeffs[10] += rints[19] * p_ecoeffs_cd_tsp[485];
                    p_rints_x_ecoeffs[10] += rints[7] * p_ecoeffs_cd_tsp[185];
                    p_rints_x_ecoeffs[10] += rints[0] * p_ecoeffs_cd_tsp[10];
                    p_rints_x_ecoeffs[10] += rints[1] * p_ecoeffs_cd_tsp[35];
                    p_rints_x_ecoeffs[10] += rints[4] * p_ecoeffs_cd_tsp[110];
                    p_rints_x_ecoeffs[10] += rints[8] * p_ecoeffs_cd_tsp[210];
                    p_rints_x_ecoeffs[10] += rints[24] * p_ecoeffs_cd_tsp[610];
                    p_rints_x_ecoeffs[10] += rints[31] * p_ecoeffs_cd_tsp[785];
                    p_rints_x_ecoeffs[10] += rints[9] * p_ecoeffs_cd_tsp[235];
                    p_rints_x_ecoeffs[10] += rints[14] * p_ecoeffs_cd_tsp[360];
                    p_rints_x_ecoeffs[11] += rints[6] * p_ecoeffs_cd_tsp[161];
                    p_rints_x_ecoeffs[11] += rints[5] * p_ecoeffs_cd_tsp[136];
                    p_rints_x_ecoeffs[11] += rints[2] * p_ecoeffs_cd_tsp[61];
                    p_rints_x_ecoeffs[11] += rints[0] * p_ecoeffs_cd_tsp[11];
                    p_rints_x_ecoeffs[11] += rints[28] * p_ecoeffs_cd_tsp[711];
                    p_rints_x_ecoeffs[11] += rints[1] * p_ecoeffs_cd_tsp[36];
                    p_rints_x_ecoeffs[11] += rints[9] * p_ecoeffs_cd_tsp[236];
                    p_rints_x_ecoeffs[11] += rints[18] * p_ecoeffs_cd_tsp[461];
                    p_rints_x_ecoeffs[11] += rints[15] * p_ecoeffs_cd_tsp[386];
                    p_rints_x_ecoeffs[11] += rints[3] * p_ecoeffs_cd_tsp[86];
                    p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_cd_tsp[361];
                    p_rints_x_ecoeffs[11] += rints[8] * p_ecoeffs_cd_tsp[211];
                    p_rints_x_ecoeffs[12] += rints[17] * p_ecoeffs_cd_tsp[437];
                    p_rints_x_ecoeffs[12] += rints[2] * p_ecoeffs_cd_tsp[62];
                    p_rints_x_ecoeffs[12] += rints[7] * p_ecoeffs_cd_tsp[187];
                    p_rints_x_ecoeffs[12] += rints[0] * p_ecoeffs_cd_tsp[12];
                    p_rints_x_ecoeffs[12] += rints[9] * p_ecoeffs_cd_tsp[237];
                    p_rints_x_ecoeffs[12] += rints[18] * p_ecoeffs_cd_tsp[462];
                    p_rints_x_ecoeffs[12] += rints[32] * p_ecoeffs_cd_tsp[812];
                    p_rints_x_ecoeffs[12] += rints[3] * p_ecoeffs_cd_tsp[87];
                    p_rints_x_ecoeffs[12] += rints[8] * p_ecoeffs_cd_tsp[212];
                    p_rints_x_ecoeffs[13] += rints[24] * p_ecoeffs_cd_tsp[613];
                    p_rints_x_ecoeffs[13] += rints[6] * p_ecoeffs_cd_tsp[163];
                    p_rints_x_ecoeffs[13] += rints[17] * p_ecoeffs_cd_tsp[438];
                    p_rints_x_ecoeffs[13] += rints[5] * p_ecoeffs_cd_tsp[138];
                    p_rints_x_ecoeffs[13] += rints[31] * p_ecoeffs_cd_tsp[788];
                    p_rints_x_ecoeffs[13] += rints[2] * p_ecoeffs_cd_tsp[63];
                    p_rints_x_ecoeffs[13] += rints[11] * p_ecoeffs_cd_tsp[288];
                    p_rints_x_ecoeffs[13] += rints[7] * p_ecoeffs_cd_tsp[188];
                    p_rints_x_ecoeffs[13] += rints[0] * p_ecoeffs_cd_tsp[13];
                    p_rints_x_ecoeffs[13] += rints[16] * p_ecoeffs_cd_tsp[413];
                    p_rints_x_ecoeffs[13] += rints[1] * p_ecoeffs_cd_tsp[38];
                    p_rints_x_ecoeffs[13] += rints[12] * p_ecoeffs_cd_tsp[313];
                    p_rints_x_ecoeffs[13] += rints[4] * p_ecoeffs_cd_tsp[113];
                    p_rints_x_ecoeffs[13] += rints[3] * p_ecoeffs_cd_tsp[88];
                    p_rints_x_ecoeffs[13] += rints[14] * p_ecoeffs_cd_tsp[363];
                    p_rints_x_ecoeffs[13] += rints[8] * p_ecoeffs_cd_tsp[213];
                    p_rints_x_ecoeffs[14] += rints[27] * p_ecoeffs_cd_tsp[689];
                    p_rints_x_ecoeffs[14] += rints[17] * p_ecoeffs_cd_tsp[439];
                    p_rints_x_ecoeffs[14] += rints[6] * p_ecoeffs_cd_tsp[164];
                    p_rints_x_ecoeffs[14] += rints[5] * p_ecoeffs_cd_tsp[139];
                    p_rints_x_ecoeffs[14] += rints[2] * p_ecoeffs_cd_tsp[64];
                    p_rints_x_ecoeffs[14] += rints[13] * p_ecoeffs_cd_tsp[339];
                    p_rints_x_ecoeffs[14] += rints[7] * p_ecoeffs_cd_tsp[189];
                    p_rints_x_ecoeffs[14] += rints[0] * p_ecoeffs_cd_tsp[14];
                    p_rints_x_ecoeffs[14] += rints[1] * p_ecoeffs_cd_tsp[39];
                    p_rints_x_ecoeffs[14] += rints[3] * p_ecoeffs_cd_tsp[89];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_cd_tsp[364];
                    p_rints_x_ecoeffs[14] += rints[8] * p_ecoeffs_cd_tsp[214];
                    p_rints_x_ecoeffs[15] += rints[25] * p_ecoeffs_cd_tsp[640];
                    p_rints_x_ecoeffs[15] += rints[2] * p_ecoeffs_cd_tsp[65];
                    p_rints_x_ecoeffs[15] += rints[11] * p_ecoeffs_cd_tsp[290];
                    p_rints_x_ecoeffs[15] += rints[16] * p_ecoeffs_cd_tsp[415];
                    p_rints_x_ecoeffs[15] += rints[3] * p_ecoeffs_cd_tsp[90];
                    p_rints_x_ecoeffs[15] += rints[17] * p_ecoeffs_cd_tsp[440];
                    p_rints_x_ecoeffs[15] += rints[6] * p_ecoeffs_cd_tsp[165];
                    p_rints_x_ecoeffs[15] += rints[5] * p_ecoeffs_cd_tsp[140];
                    p_rints_x_ecoeffs[15] += rints[10] * p_ecoeffs_cd_tsp[265];
                    p_rints_x_ecoeffs[15] += rints[20] * p_ecoeffs_cd_tsp[515];
                    p_rints_x_ecoeffs[15] += rints[12] * p_ecoeffs_cd_tsp[315];
                    p_rints_x_ecoeffs[15] += rints[18] * p_ecoeffs_cd_tsp[465];
                    p_rints_x_ecoeffs[15] += rints[7] * p_ecoeffs_cd_tsp[190];
                    p_rints_x_ecoeffs[15] += rints[0] * p_ecoeffs_cd_tsp[15];
                    p_rints_x_ecoeffs[15] += rints[30] * p_ecoeffs_cd_tsp[765];
                    p_rints_x_ecoeffs[15] += rints[1] * p_ecoeffs_cd_tsp[40];
                    p_rints_x_ecoeffs[15] += rints[4] * p_ecoeffs_cd_tsp[115];
                    p_rints_x_ecoeffs[15] += rints[23] * p_ecoeffs_cd_tsp[590];
                    p_rints_x_ecoeffs[15] += rints[8] * p_ecoeffs_cd_tsp[215];
                    p_rints_x_ecoeffs[15] += rints[13] * p_ecoeffs_cd_tsp[340];
                    p_rints_x_ecoeffs[15] += rints[9] * p_ecoeffs_cd_tsp[240];
                    p_rints_x_ecoeffs[15] += rints[32] * p_ecoeffs_cd_tsp[815];
                    p_rints_x_ecoeffs[15] += rints[15] * p_ecoeffs_cd_tsp[390];
                    p_rints_x_ecoeffs[16] += rints[27] * p_ecoeffs_cd_tsp[691];
                    p_rints_x_ecoeffs[16] += rints[17] * p_ecoeffs_cd_tsp[441];
                    p_rints_x_ecoeffs[16] += rints[6] * p_ecoeffs_cd_tsp[166];
                    p_rints_x_ecoeffs[16] += rints[5] * p_ecoeffs_cd_tsp[141];
                    p_rints_x_ecoeffs[16] += rints[2] * p_ecoeffs_cd_tsp[66];
                    p_rints_x_ecoeffs[16] += rints[13] * p_ecoeffs_cd_tsp[341];
                    p_rints_x_ecoeffs[16] += rints[7] * p_ecoeffs_cd_tsp[191];
                    p_rints_x_ecoeffs[16] += rints[10] * p_ecoeffs_cd_tsp[266];
                    p_rints_x_ecoeffs[16] += rints[0] * p_ecoeffs_cd_tsp[16];
                    p_rints_x_ecoeffs[16] += rints[1] * p_ecoeffs_cd_tsp[41];
                    p_rints_x_ecoeffs[16] += rints[12] * p_ecoeffs_cd_tsp[316];
                    p_rints_x_ecoeffs[16] += rints[4] * p_ecoeffs_cd_tsp[116];
                    p_rints_x_ecoeffs[16] += rints[22] * p_ecoeffs_cd_tsp[566];
                    p_rints_x_ecoeffs[16] += rints[3] * p_ecoeffs_cd_tsp[91];
                    p_rints_x_ecoeffs[16] += rints[14] * p_ecoeffs_cd_tsp[366];
                    p_rints_x_ecoeffs[16] += rints[8] * p_ecoeffs_cd_tsp[216];
                    p_rints_x_ecoeffs[17] += rints[24] * p_ecoeffs_cd_tsp[617];
                    p_rints_x_ecoeffs[17] += rints[6] * p_ecoeffs_cd_tsp[167];
                    p_rints_x_ecoeffs[17] += rints[17] * p_ecoeffs_cd_tsp[442];
                    p_rints_x_ecoeffs[17] += rints[5] * p_ecoeffs_cd_tsp[142];
                    p_rints_x_ecoeffs[17] += rints[31] * p_ecoeffs_cd_tsp[792];
                    p_rints_x_ecoeffs[17] += rints[2] * p_ecoeffs_cd_tsp[67];
                    p_rints_x_ecoeffs[17] += rints[11] * p_ecoeffs_cd_tsp[292];
                    p_rints_x_ecoeffs[17] += rints[7] * p_ecoeffs_cd_tsp[192];
                    p_rints_x_ecoeffs[17] += rints[0] * p_ecoeffs_cd_tsp[17];
                    p_rints_x_ecoeffs[17] += rints[16] * p_ecoeffs_cd_tsp[417];
                    p_rints_x_ecoeffs[17] += rints[1] * p_ecoeffs_cd_tsp[42];
                    p_rints_x_ecoeffs[17] += rints[12] * p_ecoeffs_cd_tsp[317];
                    p_rints_x_ecoeffs[17] += rints[4] * p_ecoeffs_cd_tsp[117];
                    p_rints_x_ecoeffs[17] += rints[3] * p_ecoeffs_cd_tsp[92];
                    p_rints_x_ecoeffs[17] += rints[14] * p_ecoeffs_cd_tsp[367];
                    p_rints_x_ecoeffs[17] += rints[8] * p_ecoeffs_cd_tsp[217];
                    p_rints_x_ecoeffs[18] += rints[5] * p_ecoeffs_cd_tsp[143];
                    p_rints_x_ecoeffs[18] += rints[2] * p_ecoeffs_cd_tsp[68];
                    p_rints_x_ecoeffs[18] += rints[13] * p_ecoeffs_cd_tsp[343];
                    p_rints_x_ecoeffs[18] += rints[7] * p_ecoeffs_cd_tsp[193];
                    p_rints_x_ecoeffs[18] += rints[11] * p_ecoeffs_cd_tsp[293];
                    p_rints_x_ecoeffs[18] += rints[0] * p_ecoeffs_cd_tsp[18];
                    p_rints_x_ecoeffs[18] += rints[16] * p_ecoeffs_cd_tsp[418];
                    p_rints_x_ecoeffs[18] += rints[10] * p_ecoeffs_cd_tsp[268];
                    p_rints_x_ecoeffs[18] += rints[30] * p_ecoeffs_cd_tsp[768];
                    p_rints_x_ecoeffs[18] += rints[20] * p_ecoeffs_cd_tsp[518];
                    p_rints_x_ecoeffs[18] += rints[1] * p_ecoeffs_cd_tsp[43];
                    p_rints_x_ecoeffs[18] += rints[4] * p_ecoeffs_cd_tsp[118];
                    p_rints_x_ecoeffs[18] += rints[23] * p_ecoeffs_cd_tsp[593];
                    p_rints_x_ecoeffs[19] += rints[5] * p_ecoeffs_cd_tsp[144];
                    p_rints_x_ecoeffs[19] += rints[2] * p_ecoeffs_cd_tsp[69];
                    p_rints_x_ecoeffs[19] += rints[13] * p_ecoeffs_cd_tsp[344];
                    p_rints_x_ecoeffs[19] += rints[7] * p_ecoeffs_cd_tsp[194];
                    p_rints_x_ecoeffs[19] += rints[11] * p_ecoeffs_cd_tsp[294];
                    p_rints_x_ecoeffs[19] += rints[0] * p_ecoeffs_cd_tsp[19];
                    p_rints_x_ecoeffs[19] += rints[16] * p_ecoeffs_cd_tsp[419];
                    p_rints_x_ecoeffs[19] += rints[10] * p_ecoeffs_cd_tsp[269];
                    p_rints_x_ecoeffs[19] += rints[21] * p_ecoeffs_cd_tsp[544];
                    p_rints_x_ecoeffs[19] += rints[1] * p_ecoeffs_cd_tsp[44];
                    p_rints_x_ecoeffs[19] += rints[4] * p_ecoeffs_cd_tsp[119];
                    p_rints_x_ecoeffs[19] += rints[26] * p_ecoeffs_cd_tsp[669];
                    p_rints_x_ecoeffs[20] += rints[2] * p_ecoeffs_cd_tsp[70];
                    p_rints_x_ecoeffs[20] += rints[11] * p_ecoeffs_cd_tsp[295];
                    p_rints_x_ecoeffs[20] += rints[16] * p_ecoeffs_cd_tsp[420];
                    p_rints_x_ecoeffs[20] += rints[26] * p_ecoeffs_cd_tsp[670];
                    p_rints_x_ecoeffs[20] += rints[3] * p_ecoeffs_cd_tsp[95];
                    p_rints_x_ecoeffs[20] += rints[6] * p_ecoeffs_cd_tsp[170];
                    p_rints_x_ecoeffs[20] += rints[5] * p_ecoeffs_cd_tsp[145];
                    p_rints_x_ecoeffs[20] += rints[10] * p_ecoeffs_cd_tsp[270];
                    p_rints_x_ecoeffs[20] += rints[18] * p_ecoeffs_cd_tsp[470];
                    p_rints_x_ecoeffs[20] += rints[7] * p_ecoeffs_cd_tsp[195];
                    p_rints_x_ecoeffs[20] += rints[0] * p_ecoeffs_cd_tsp[20];
                    p_rints_x_ecoeffs[20] += rints[28] * p_ecoeffs_cd_tsp[720];
                    p_rints_x_ecoeffs[20] += rints[1] * p_ecoeffs_cd_tsp[45];
                    p_rints_x_ecoeffs[20] += rints[4] * p_ecoeffs_cd_tsp[120];
                    p_rints_x_ecoeffs[20] += rints[8] * p_ecoeffs_cd_tsp[220];
                    p_rints_x_ecoeffs[20] += rints[13] * p_ecoeffs_cd_tsp[345];
                    p_rints_x_ecoeffs[20] += rints[21] * p_ecoeffs_cd_tsp[545];
                    p_rints_x_ecoeffs[20] += rints[9] * p_ecoeffs_cd_tsp[245];
                    p_rints_x_ecoeffs[20] += rints[15] * p_ecoeffs_cd_tsp[395];
                    p_rints_x_ecoeffs[20] += rints[14] * p_ecoeffs_cd_tsp[370];
                    p_rints_x_ecoeffs[21] += rints[24] * p_ecoeffs_cd_tsp[621];
                    p_rints_x_ecoeffs[21] += rints[6] * p_ecoeffs_cd_tsp[171];
                    p_rints_x_ecoeffs[21] += rints[5] * p_ecoeffs_cd_tsp[146];
                    p_rints_x_ecoeffs[21] += rints[2] * p_ecoeffs_cd_tsp[71];
                    p_rints_x_ecoeffs[21] += rints[11] * p_ecoeffs_cd_tsp[296];
                    p_rints_x_ecoeffs[21] += rints[0] * p_ecoeffs_cd_tsp[21];
                    p_rints_x_ecoeffs[21] += rints[1] * p_ecoeffs_cd_tsp[46];
                    p_rints_x_ecoeffs[21] += rints[12] * p_ecoeffs_cd_tsp[321];
                    p_rints_x_ecoeffs[21] += rints[4] * p_ecoeffs_cd_tsp[121];
                    p_rints_x_ecoeffs[21] += rints[3] * p_ecoeffs_cd_tsp[96];
                    p_rints_x_ecoeffs[21] += rints[14] * p_ecoeffs_cd_tsp[371];
                    p_rints_x_ecoeffs[21] += rints[8] * p_ecoeffs_cd_tsp[221];
                    p_rints_x_ecoeffs[22] += rints[27] * p_ecoeffs_cd_tsp[697];
                    p_rints_x_ecoeffs[22] += rints[17] * p_ecoeffs_cd_tsp[447];
                    p_rints_x_ecoeffs[22] += rints[6] * p_ecoeffs_cd_tsp[172];
                    p_rints_x_ecoeffs[22] += rints[5] * p_ecoeffs_cd_tsp[147];
                    p_rints_x_ecoeffs[22] += rints[2] * p_ecoeffs_cd_tsp[72];
                    p_rints_x_ecoeffs[22] += rints[13] * p_ecoeffs_cd_tsp[347];
                    p_rints_x_ecoeffs[22] += rints[7] * p_ecoeffs_cd_tsp[197];
                    p_rints_x_ecoeffs[22] += rints[0] * p_ecoeffs_cd_tsp[22];
                    p_rints_x_ecoeffs[22] += rints[1] * p_ecoeffs_cd_tsp[47];
                    p_rints_x_ecoeffs[22] += rints[3] * p_ecoeffs_cd_tsp[97];
                    p_rints_x_ecoeffs[22] += rints[14] * p_ecoeffs_cd_tsp[372];
                    p_rints_x_ecoeffs[22] += rints[8] * p_ecoeffs_cd_tsp[222];
                    p_rints_x_ecoeffs[23] += rints[5] * p_ecoeffs_cd_tsp[148];
                    p_rints_x_ecoeffs[23] += rints[2] * p_ecoeffs_cd_tsp[73];
                    p_rints_x_ecoeffs[23] += rints[13] * p_ecoeffs_cd_tsp[348];
                    p_rints_x_ecoeffs[23] += rints[7] * p_ecoeffs_cd_tsp[198];
                    p_rints_x_ecoeffs[23] += rints[11] * p_ecoeffs_cd_tsp[298];
                    p_rints_x_ecoeffs[23] += rints[0] * p_ecoeffs_cd_tsp[23];
                    p_rints_x_ecoeffs[23] += rints[16] * p_ecoeffs_cd_tsp[423];
                    p_rints_x_ecoeffs[23] += rints[10] * p_ecoeffs_cd_tsp[273];
                    p_rints_x_ecoeffs[23] += rints[21] * p_ecoeffs_cd_tsp[548];
                    p_rints_x_ecoeffs[23] += rints[1] * p_ecoeffs_cd_tsp[48];
                    p_rints_x_ecoeffs[23] += rints[4] * p_ecoeffs_cd_tsp[123];
                    p_rints_x_ecoeffs[23] += rints[26] * p_ecoeffs_cd_tsp[673];
                    p_rints_x_ecoeffs[24] += rints[5] * p_ecoeffs_cd_tsp[149];
                    p_rints_x_ecoeffs[24] += rints[2] * p_ecoeffs_cd_tsp[74];
                    p_rints_x_ecoeffs[24] += rints[13] * p_ecoeffs_cd_tsp[349];
                    p_rints_x_ecoeffs[24] += rints[7] * p_ecoeffs_cd_tsp[199];
                    p_rints_x_ecoeffs[24] += rints[11] * p_ecoeffs_cd_tsp[299];
                    p_rints_x_ecoeffs[24] += rints[0] * p_ecoeffs_cd_tsp[24];
                    p_rints_x_ecoeffs[24] += rints[1] * p_ecoeffs_cd_tsp[49];
                    p_rints_x_ecoeffs[24] += rints[4] * p_ecoeffs_cd_tsp[124];
                    p_rints_x_ecoeffs[24] += rints[23] * p_ecoeffs_cd_tsp[599];
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
            eri4_batch[15] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[15];
            eri4_batch[16] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[16];
            eri4_batch[17] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[17];
            eri4_batch[18] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[18];
            eri4_batch[19] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[19];
            eri4_batch[20] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[20];
            eri4_batch[21] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[21];
            eri4_batch[22] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[22];
            eri4_batch[23] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[23];
            eri4_batch[24] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[24];
        }
}

template<> void lible::ints::two::eri4Kernel<0, 0, 3, 1>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 0, lb = 0, lc = 3, ld = 1;
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

                    p_rints_x_ecoeffs[0] += rints[34] * p_ecoeffs_cd_tsp[714];
                    p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_cd_tsp[357];
                    p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_cd_tsp[126];
                    p_rints_x_ecoeffs[0] += rints[25] * p_ecoeffs_cd_tsp[525];
                    p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_cd_tsp[42];
                    p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_cd_tsp[399];
                    p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_cd_tsp[147];
                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_cd_tsp[21];
                    p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_cd_tsp[252];
                    p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_cd_tsp[189];
                    p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_cd_tsp[84];
                    p_rints_x_ecoeffs[0] += rints[18] * p_ecoeffs_cd_tsp[378];
                    p_rints_x_ecoeffs[0] += rints[32] * p_ecoeffs_cd_tsp[672];
                    p_rints_x_ecoeffs[0] += rints[15] * p_ecoeffs_cd_tsp[315];
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[63];
                    p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_cd_tsp[168];
                    p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_cd_tsp[43];
                    p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_cd_tsp[64];
                    p_rints_x_ecoeffs[1] += rints[17] * p_ecoeffs_cd_tsp[358];
                    p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_cd_tsp[127];
                    p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_cd_tsp[106];
                    p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_cd_tsp[211];
                    p_rints_x_ecoeffs[1] += rints[12] * p_ecoeffs_cd_tsp[253];
                    p_rints_x_ecoeffs[1] += rints[27] * p_ecoeffs_cd_tsp[568];
                    p_rints_x_ecoeffs[1] += rints[19] * p_ecoeffs_cd_tsp[400];
                    p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_cd_tsp[148];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[22];
                    p_rints_x_ecoeffs[1] += rints[29] * p_ecoeffs_cd_tsp[610];
                    p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_cd_tsp[85];
                    p_rints_x_ecoeffs[1] += rints[8] * p_ecoeffs_cd_tsp[169];
                    p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_cd_tsp[274];
                    p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_cd_tsp[190];
                    p_rints_x_ecoeffs[1] += rints[22] * p_ecoeffs_cd_tsp[463];
                    p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_cd_tsp[316];
                    p_rints_x_ecoeffs[1] += rints[14] * p_ecoeffs_cd_tsp[295];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[44];
                    p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_cd_tsp[233];
                    p_rints_x_ecoeffs[2] += rints[33] * p_ecoeffs_cd_tsp[695];
                    p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_cd_tsp[338];
                    p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_cd_tsp[65];
                    p_rints_x_ecoeffs[2] += rints[17] * p_ecoeffs_cd_tsp[359];
                    p_rints_x_ecoeffs[2] += rints[6] * p_ecoeffs_cd_tsp[128];
                    p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_cd_tsp[107];
                    p_rints_x_ecoeffs[2] += rints[12] * p_ecoeffs_cd_tsp[254];
                    p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_cd_tsp[380];
                    p_rints_x_ecoeffs[2] += rints[19] * p_ecoeffs_cd_tsp[401];
                    p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_cd_tsp[149];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_cd_tsp[23];
                    p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_cd_tsp[86];
                    p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_cd_tsp[170];
                    p_rints_x_ecoeffs[2] += rints[24] * p_ecoeffs_cd_tsp[506];
                    p_rints_x_ecoeffs[2] += rints[31] * p_ecoeffs_cd_tsp[653];
                    p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_cd_tsp[191];
                    p_rints_x_ecoeffs[2] += rints[14] * p_ecoeffs_cd_tsp[296];
                    p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_cd_tsp[45];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[66];
                    p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_cd_tsp[360];
                    p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_cd_tsp[129];
                    p_rints_x_ecoeffs[3] += rints[5] * p_ecoeffs_cd_tsp[108];
                    p_rints_x_ecoeffs[3] += rints[10] * p_ecoeffs_cd_tsp[213];
                    p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_cd_tsp[255];
                    p_rints_x_ecoeffs[3] += rints[27] * p_ecoeffs_cd_tsp[570];
                    p_rints_x_ecoeffs[3] += rints[19] * p_ecoeffs_cd_tsp[402];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[150];
                    p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_cd_tsp[3];
                    p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_cd_tsp[24];
                    p_rints_x_ecoeffs[3] += rints[29] * p_ecoeffs_cd_tsp[612];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[87];
                    p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_cd_tsp[171];
                    p_rints_x_ecoeffs[3] += rints[13] * p_ecoeffs_cd_tsp[276];
                    p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_cd_tsp[192];
                    p_rints_x_ecoeffs[3] += rints[22] * p_ecoeffs_cd_tsp[465];
                    p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_cd_tsp[318];
                    p_rints_x_ecoeffs[3] += rints[14] * p_ecoeffs_cd_tsp[297];
                    p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_cd_tsp[130];
                    p_rints_x_ecoeffs[4] += rints[25] * p_ecoeffs_cd_tsp[529];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[109];
                    p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_cd_tsp[46];
                    p_rints_x_ecoeffs[4] += rints[13] * p_ecoeffs_cd_tsp[277];
                    p_rints_x_ecoeffs[4] += rints[7] * p_ecoeffs_cd_tsp[151];
                    p_rints_x_ecoeffs[4] += rints[11] * p_ecoeffs_cd_tsp[235];
                    p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[10] * p_ecoeffs_cd_tsp[214];
                    p_rints_x_ecoeffs[4] += rints[20] * p_ecoeffs_cd_tsp[424];
                    p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_cd_tsp[25];
                    p_rints_x_ecoeffs[4] += rints[12] * p_ecoeffs_cd_tsp[256];
                    p_rints_x_ecoeffs[4] += rints[9] * p_ecoeffs_cd_tsp[193];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[88];
                    p_rints_x_ecoeffs[4] += rints[15] * p_ecoeffs_cd_tsp[319];
                    p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_cd_tsp[67];
                    p_rints_x_ecoeffs[4] += rints[23] * p_ecoeffs_cd_tsp[487];
                    p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_cd_tsp[47];
                    p_rints_x_ecoeffs[5] += rints[11] * p_ecoeffs_cd_tsp[236];
                    p_rints_x_ecoeffs[5] += rints[16] * p_ecoeffs_cd_tsp[341];
                    p_rints_x_ecoeffs[5] += rints[26] * p_ecoeffs_cd_tsp[551];
                    p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_cd_tsp[68];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[131];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[110];
                    p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_cd_tsp[215];
                    p_rints_x_ecoeffs[5] += rints[18] * p_ecoeffs_cd_tsp[383];
                    p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_cd_tsp[152];
                    p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_cd_tsp[5];
                    p_rints_x_ecoeffs[5] += rints[28] * p_ecoeffs_cd_tsp[593];
                    p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_cd_tsp[26];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[89];
                    p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_cd_tsp[173];
                    p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_cd_tsp[278];
                    p_rints_x_ecoeffs[5] += rints[21] * p_ecoeffs_cd_tsp[446];
                    p_rints_x_ecoeffs[5] += rints[9] * p_ecoeffs_cd_tsp[194];
                    p_rints_x_ecoeffs[5] += rints[15] * p_ecoeffs_cd_tsp[320];
                    p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_cd_tsp[299];
                    p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_cd_tsp[48];
                    p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_cd_tsp[237];
                    p_rints_x_ecoeffs[6] += rints[33] * p_ecoeffs_cd_tsp[699];
                    p_rints_x_ecoeffs[6] += rints[16] * p_ecoeffs_cd_tsp[342];
                    p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_cd_tsp[69];
                    p_rints_x_ecoeffs[6] += rints[17] * p_ecoeffs_cd_tsp[363];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_cd_tsp[132];
                    p_rints_x_ecoeffs[6] += rints[5] * p_ecoeffs_cd_tsp[111];
                    p_rints_x_ecoeffs[6] += rints[12] * p_ecoeffs_cd_tsp[258];
                    p_rints_x_ecoeffs[6] += rints[18] * p_ecoeffs_cd_tsp[384];
                    p_rints_x_ecoeffs[6] += rints[19] * p_ecoeffs_cd_tsp[405];
                    p_rints_x_ecoeffs[6] += rints[7] * p_ecoeffs_cd_tsp[153];
                    p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_cd_tsp[6];
                    p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_cd_tsp[27];
                    p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_cd_tsp[90];
                    p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_cd_tsp[174];
                    p_rints_x_ecoeffs[6] += rints[24] * p_ecoeffs_cd_tsp[510];
                    p_rints_x_ecoeffs[6] += rints[31] * p_ecoeffs_cd_tsp[657];
                    p_rints_x_ecoeffs[6] += rints[9] * p_ecoeffs_cd_tsp[195];
                    p_rints_x_ecoeffs[6] += rints[14] * p_ecoeffs_cd_tsp[300];
                    p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_cd_tsp[49];
                    p_rints_x_ecoeffs[7] += rints[11] * p_ecoeffs_cd_tsp[238];
                    p_rints_x_ecoeffs[7] += rints[16] * p_ecoeffs_cd_tsp[343];
                    p_rints_x_ecoeffs[7] += rints[26] * p_ecoeffs_cd_tsp[553];
                    p_rints_x_ecoeffs[7] += rints[3] * p_ecoeffs_cd_tsp[70];
                    p_rints_x_ecoeffs[7] += rints[6] * p_ecoeffs_cd_tsp[133];
                    p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_cd_tsp[112];
                    p_rints_x_ecoeffs[7] += rints[10] * p_ecoeffs_cd_tsp[217];
                    p_rints_x_ecoeffs[7] += rints[18] * p_ecoeffs_cd_tsp[385];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_cd_tsp[154];
                    p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_cd_tsp[7];
                    p_rints_x_ecoeffs[7] += rints[28] * p_ecoeffs_cd_tsp[595];
                    p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_cd_tsp[28];
                    p_rints_x_ecoeffs[7] += rints[4] * p_ecoeffs_cd_tsp[91];
                    p_rints_x_ecoeffs[7] += rints[8] * p_ecoeffs_cd_tsp[175];
                    p_rints_x_ecoeffs[7] += rints[13] * p_ecoeffs_cd_tsp[280];
                    p_rints_x_ecoeffs[7] += rints[21] * p_ecoeffs_cd_tsp[448];
                    p_rints_x_ecoeffs[7] += rints[9] * p_ecoeffs_cd_tsp[196];
                    p_rints_x_ecoeffs[7] += rints[15] * p_ecoeffs_cd_tsp[322];
                    p_rints_x_ecoeffs[7] += rints[14] * p_ecoeffs_cd_tsp[301];
                    p_rints_x_ecoeffs[8] += rints[17] * p_ecoeffs_cd_tsp[365];
                    p_rints_x_ecoeffs[8] += rints[23] * p_ecoeffs_cd_tsp[491];
                    p_rints_x_ecoeffs[8] += rints[5] * p_ecoeffs_cd_tsp[113];
                    p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_cd_tsp[50];
                    p_rints_x_ecoeffs[8] += rints[13] * p_ecoeffs_cd_tsp[281];
                    p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_cd_tsp[155];
                    p_rints_x_ecoeffs[8] += rints[11] * p_ecoeffs_cd_tsp[239];
                    p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[8] += rints[16] * p_ecoeffs_cd_tsp[344];
                    p_rints_x_ecoeffs[8] += rints[30] * p_ecoeffs_cd_tsp[638];
                    p_rints_x_ecoeffs[8] += rints[1] * p_ecoeffs_cd_tsp[29];
                    p_rints_x_ecoeffs[8] += rints[9] * p_ecoeffs_cd_tsp[197];
                    p_rints_x_ecoeffs[8] += rints[4] * p_ecoeffs_cd_tsp[92];
                    p_rints_x_ecoeffs[8] += rints[18] * p_ecoeffs_cd_tsp[386];
                    p_rints_x_ecoeffs[8] += rints[32] * p_ecoeffs_cd_tsp[680];
                    p_rints_x_ecoeffs[8] += rints[3] * p_ecoeffs_cd_tsp[71];
                    p_rints_x_ecoeffs[8] += rints[8] * p_ecoeffs_cd_tsp[176];
                    p_rints_x_ecoeffs[9] += rints[17] * p_ecoeffs_cd_tsp[366];
                    p_rints_x_ecoeffs[9] += rints[6] * p_ecoeffs_cd_tsp[135];
                    p_rints_x_ecoeffs[9] += rints[25] * p_ecoeffs_cd_tsp[534];
                    p_rints_x_ecoeffs[9] += rints[2] * p_ecoeffs_cd_tsp[51];
                    p_rints_x_ecoeffs[9] += rints[7] * p_ecoeffs_cd_tsp[156];
                    p_rints_x_ecoeffs[9] += rints[0] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[9] += rints[1] * p_ecoeffs_cd_tsp[30];
                    p_rints_x_ecoeffs[9] += rints[12] * p_ecoeffs_cd_tsp[261];
                    p_rints_x_ecoeffs[9] += rints[9] * p_ecoeffs_cd_tsp[198];
                    p_rints_x_ecoeffs[9] += rints[4] * p_ecoeffs_cd_tsp[93];
                    p_rints_x_ecoeffs[9] += rints[18] * p_ecoeffs_cd_tsp[387];
                    p_rints_x_ecoeffs[9] += rints[32] * p_ecoeffs_cd_tsp[681];
                    p_rints_x_ecoeffs[9] += rints[15] * p_ecoeffs_cd_tsp[324];
                    p_rints_x_ecoeffs[9] += rints[3] * p_ecoeffs_cd_tsp[72];
                    p_rints_x_ecoeffs[9] += rints[8] * p_ecoeffs_cd_tsp[177];
                    p_rints_x_ecoeffs[10] += rints[27] * p_ecoeffs_cd_tsp[577];
                    p_rints_x_ecoeffs[10] += rints[17] * p_ecoeffs_cd_tsp[367];
                    p_rints_x_ecoeffs[10] += rints[6] * p_ecoeffs_cd_tsp[136];
                    p_rints_x_ecoeffs[10] += rints[5] * p_ecoeffs_cd_tsp[115];
                    p_rints_x_ecoeffs[10] += rints[2] * p_ecoeffs_cd_tsp[52];
                    p_rints_x_ecoeffs[10] += rints[13] * p_ecoeffs_cd_tsp[283];
                    p_rints_x_ecoeffs[10] += rints[7] * p_ecoeffs_cd_tsp[157];
                    p_rints_x_ecoeffs[10] += rints[10] * p_ecoeffs_cd_tsp[220];
                    p_rints_x_ecoeffs[10] += rints[0] * p_ecoeffs_cd_tsp[10];
                    p_rints_x_ecoeffs[10] += rints[1] * p_ecoeffs_cd_tsp[31];
                    p_rints_x_ecoeffs[10] += rints[12] * p_ecoeffs_cd_tsp[262];
                    p_rints_x_ecoeffs[10] += rints[4] * p_ecoeffs_cd_tsp[94];
                    p_rints_x_ecoeffs[10] += rints[22] * p_ecoeffs_cd_tsp[472];
                    p_rints_x_ecoeffs[10] += rints[3] * p_ecoeffs_cd_tsp[73];
                    p_rints_x_ecoeffs[10] += rints[14] * p_ecoeffs_cd_tsp[304];
                    p_rints_x_ecoeffs[10] += rints[8] * p_ecoeffs_cd_tsp[178];
                    p_rints_x_ecoeffs[11] += rints[24] * p_ecoeffs_cd_tsp[515];
                    p_rints_x_ecoeffs[11] += rints[6] * p_ecoeffs_cd_tsp[137];
                    p_rints_x_ecoeffs[11] += rints[17] * p_ecoeffs_cd_tsp[368];
                    p_rints_x_ecoeffs[11] += rints[5] * p_ecoeffs_cd_tsp[116];
                    p_rints_x_ecoeffs[11] += rints[31] * p_ecoeffs_cd_tsp[662];
                    p_rints_x_ecoeffs[11] += rints[2] * p_ecoeffs_cd_tsp[53];
                    p_rints_x_ecoeffs[11] += rints[11] * p_ecoeffs_cd_tsp[242];
                    p_rints_x_ecoeffs[11] += rints[7] * p_ecoeffs_cd_tsp[158];
                    p_rints_x_ecoeffs[11] += rints[0] * p_ecoeffs_cd_tsp[11];
                    p_rints_x_ecoeffs[11] += rints[16] * p_ecoeffs_cd_tsp[347];
                    p_rints_x_ecoeffs[11] += rints[1] * p_ecoeffs_cd_tsp[32];
                    p_rints_x_ecoeffs[11] += rints[12] * p_ecoeffs_cd_tsp[263];
                    p_rints_x_ecoeffs[11] += rints[4] * p_ecoeffs_cd_tsp[95];
                    p_rints_x_ecoeffs[11] += rints[3] * p_ecoeffs_cd_tsp[74];
                    p_rints_x_ecoeffs[11] += rints[14] * p_ecoeffs_cd_tsp[305];
                    p_rints_x_ecoeffs[11] += rints[8] * p_ecoeffs_cd_tsp[179];
                    p_rints_x_ecoeffs[12] += rints[6] * p_ecoeffs_cd_tsp[138];
                    p_rints_x_ecoeffs[12] += rints[5] * p_ecoeffs_cd_tsp[117];
                    p_rints_x_ecoeffs[12] += rints[2] * p_ecoeffs_cd_tsp[54];
                    p_rints_x_ecoeffs[12] += rints[0] * p_ecoeffs_cd_tsp[12];
                    p_rints_x_ecoeffs[12] += rints[28] * p_ecoeffs_cd_tsp[600];
                    p_rints_x_ecoeffs[12] += rints[1] * p_ecoeffs_cd_tsp[33];
                    p_rints_x_ecoeffs[12] += rints[9] * p_ecoeffs_cd_tsp[201];
                    p_rints_x_ecoeffs[12] += rints[18] * p_ecoeffs_cd_tsp[390];
                    p_rints_x_ecoeffs[12] += rints[15] * p_ecoeffs_cd_tsp[327];
                    p_rints_x_ecoeffs[12] += rints[3] * p_ecoeffs_cd_tsp[75];
                    p_rints_x_ecoeffs[12] += rints[14] * p_ecoeffs_cd_tsp[306];
                    p_rints_x_ecoeffs[12] += rints[8] * p_ecoeffs_cd_tsp[180];
                    p_rints_x_ecoeffs[13] += rints[24] * p_ecoeffs_cd_tsp[517];
                    p_rints_x_ecoeffs[13] += rints[6] * p_ecoeffs_cd_tsp[139];
                    p_rints_x_ecoeffs[13] += rints[5] * p_ecoeffs_cd_tsp[118];
                    p_rints_x_ecoeffs[13] += rints[2] * p_ecoeffs_cd_tsp[55];
                    p_rints_x_ecoeffs[13] += rints[11] * p_ecoeffs_cd_tsp[244];
                    p_rints_x_ecoeffs[13] += rints[0] * p_ecoeffs_cd_tsp[13];
                    p_rints_x_ecoeffs[13] += rints[1] * p_ecoeffs_cd_tsp[34];
                    p_rints_x_ecoeffs[13] += rints[12] * p_ecoeffs_cd_tsp[265];
                    p_rints_x_ecoeffs[13] += rints[4] * p_ecoeffs_cd_tsp[97];
                    p_rints_x_ecoeffs[13] += rints[3] * p_ecoeffs_cd_tsp[76];
                    p_rints_x_ecoeffs[13] += rints[14] * p_ecoeffs_cd_tsp[307];
                    p_rints_x_ecoeffs[13] += rints[8] * p_ecoeffs_cd_tsp[181];
                    p_rints_x_ecoeffs[14] += rints[27] * p_ecoeffs_cd_tsp[581];
                    p_rints_x_ecoeffs[14] += rints[17] * p_ecoeffs_cd_tsp[371];
                    p_rints_x_ecoeffs[14] += rints[6] * p_ecoeffs_cd_tsp[140];
                    p_rints_x_ecoeffs[14] += rints[5] * p_ecoeffs_cd_tsp[119];
                    p_rints_x_ecoeffs[14] += rints[2] * p_ecoeffs_cd_tsp[56];
                    p_rints_x_ecoeffs[14] += rints[13] * p_ecoeffs_cd_tsp[287];
                    p_rints_x_ecoeffs[14] += rints[7] * p_ecoeffs_cd_tsp[161];
                    p_rints_x_ecoeffs[14] += rints[0] * p_ecoeffs_cd_tsp[14];
                    p_rints_x_ecoeffs[14] += rints[1] * p_ecoeffs_cd_tsp[35];
                    p_rints_x_ecoeffs[14] += rints[3] * p_ecoeffs_cd_tsp[77];
                    p_rints_x_ecoeffs[14] += rints[14] * p_ecoeffs_cd_tsp[308];
                    p_rints_x_ecoeffs[14] += rints[8] * p_ecoeffs_cd_tsp[182];
                    p_rints_x_ecoeffs[15] += rints[27] * p_ecoeffs_cd_tsp[582];
                    p_rints_x_ecoeffs[15] += rints[17] * p_ecoeffs_cd_tsp[372];
                    p_rints_x_ecoeffs[15] += rints[6] * p_ecoeffs_cd_tsp[141];
                    p_rints_x_ecoeffs[15] += rints[5] * p_ecoeffs_cd_tsp[120];
                    p_rints_x_ecoeffs[15] += rints[2] * p_ecoeffs_cd_tsp[57];
                    p_rints_x_ecoeffs[15] += rints[13] * p_ecoeffs_cd_tsp[288];
                    p_rints_x_ecoeffs[15] += rints[7] * p_ecoeffs_cd_tsp[162];
                    p_rints_x_ecoeffs[15] += rints[10] * p_ecoeffs_cd_tsp[225];
                    p_rints_x_ecoeffs[15] += rints[0] * p_ecoeffs_cd_tsp[15];
                    p_rints_x_ecoeffs[15] += rints[1] * p_ecoeffs_cd_tsp[36];
                    p_rints_x_ecoeffs[15] += rints[12] * p_ecoeffs_cd_tsp[267];
                    p_rints_x_ecoeffs[15] += rints[4] * p_ecoeffs_cd_tsp[99];
                    p_rints_x_ecoeffs[15] += rints[22] * p_ecoeffs_cd_tsp[477];
                    p_rints_x_ecoeffs[15] += rints[3] * p_ecoeffs_cd_tsp[78];
                    p_rints_x_ecoeffs[15] += rints[14] * p_ecoeffs_cd_tsp[309];
                    p_rints_x_ecoeffs[15] += rints[8] * p_ecoeffs_cd_tsp[183];
                    p_rints_x_ecoeffs[16] += rints[5] * p_ecoeffs_cd_tsp[121];
                    p_rints_x_ecoeffs[16] += rints[2] * p_ecoeffs_cd_tsp[58];
                    p_rints_x_ecoeffs[16] += rints[13] * p_ecoeffs_cd_tsp[289];
                    p_rints_x_ecoeffs[16] += rints[7] * p_ecoeffs_cd_tsp[163];
                    p_rints_x_ecoeffs[16] += rints[11] * p_ecoeffs_cd_tsp[247];
                    p_rints_x_ecoeffs[16] += rints[0] * p_ecoeffs_cd_tsp[16];
                    p_rints_x_ecoeffs[16] += rints[10] * p_ecoeffs_cd_tsp[226];
                    p_rints_x_ecoeffs[16] += rints[20] * p_ecoeffs_cd_tsp[436];
                    p_rints_x_ecoeffs[16] += rints[1] * p_ecoeffs_cd_tsp[37];
                    p_rints_x_ecoeffs[16] += rints[4] * p_ecoeffs_cd_tsp[100];
                    p_rints_x_ecoeffs[16] += rints[23] * p_ecoeffs_cd_tsp[499];
                    p_rints_x_ecoeffs[17] += rints[5] * p_ecoeffs_cd_tsp[122];
                    p_rints_x_ecoeffs[17] += rints[2] * p_ecoeffs_cd_tsp[59];
                    p_rints_x_ecoeffs[17] += rints[13] * p_ecoeffs_cd_tsp[290];
                    p_rints_x_ecoeffs[17] += rints[7] * p_ecoeffs_cd_tsp[164];
                    p_rints_x_ecoeffs[17] += rints[11] * p_ecoeffs_cd_tsp[248];
                    p_rints_x_ecoeffs[17] += rints[0] * p_ecoeffs_cd_tsp[17];
                    p_rints_x_ecoeffs[17] += rints[16] * p_ecoeffs_cd_tsp[353];
                    p_rints_x_ecoeffs[17] += rints[10] * p_ecoeffs_cd_tsp[227];
                    p_rints_x_ecoeffs[17] += rints[21] * p_ecoeffs_cd_tsp[458];
                    p_rints_x_ecoeffs[17] += rints[1] * p_ecoeffs_cd_tsp[38];
                    p_rints_x_ecoeffs[17] += rints[4] * p_ecoeffs_cd_tsp[101];
                    p_rints_x_ecoeffs[17] += rints[26] * p_ecoeffs_cd_tsp[563];
                    p_rints_x_ecoeffs[18] += rints[24] * p_ecoeffs_cd_tsp[522];
                    p_rints_x_ecoeffs[18] += rints[6] * p_ecoeffs_cd_tsp[144];
                    p_rints_x_ecoeffs[18] += rints[17] * p_ecoeffs_cd_tsp[375];
                    p_rints_x_ecoeffs[18] += rints[5] * p_ecoeffs_cd_tsp[123];
                    p_rints_x_ecoeffs[18] += rints[31] * p_ecoeffs_cd_tsp[669];
                    p_rints_x_ecoeffs[18] += rints[2] * p_ecoeffs_cd_tsp[60];
                    p_rints_x_ecoeffs[18] += rints[11] * p_ecoeffs_cd_tsp[249];
                    p_rints_x_ecoeffs[18] += rints[7] * p_ecoeffs_cd_tsp[165];
                    p_rints_x_ecoeffs[18] += rints[0] * p_ecoeffs_cd_tsp[18];
                    p_rints_x_ecoeffs[18] += rints[16] * p_ecoeffs_cd_tsp[354];
                    p_rints_x_ecoeffs[18] += rints[1] * p_ecoeffs_cd_tsp[39];
                    p_rints_x_ecoeffs[18] += rints[12] * p_ecoeffs_cd_tsp[270];
                    p_rints_x_ecoeffs[18] += rints[4] * p_ecoeffs_cd_tsp[102];
                    p_rints_x_ecoeffs[18] += rints[3] * p_ecoeffs_cd_tsp[81];
                    p_rints_x_ecoeffs[18] += rints[14] * p_ecoeffs_cd_tsp[312];
                    p_rints_x_ecoeffs[18] += rints[8] * p_ecoeffs_cd_tsp[186];
                    p_rints_x_ecoeffs[19] += rints[5] * p_ecoeffs_cd_tsp[124];
                    p_rints_x_ecoeffs[19] += rints[2] * p_ecoeffs_cd_tsp[61];
                    p_rints_x_ecoeffs[19] += rints[13] * p_ecoeffs_cd_tsp[292];
                    p_rints_x_ecoeffs[19] += rints[7] * p_ecoeffs_cd_tsp[166];
                    p_rints_x_ecoeffs[19] += rints[11] * p_ecoeffs_cd_tsp[250];
                    p_rints_x_ecoeffs[19] += rints[0] * p_ecoeffs_cd_tsp[19];
                    p_rints_x_ecoeffs[19] += rints[16] * p_ecoeffs_cd_tsp[355];
                    p_rints_x_ecoeffs[19] += rints[10] * p_ecoeffs_cd_tsp[229];
                    p_rints_x_ecoeffs[19] += rints[21] * p_ecoeffs_cd_tsp[460];
                    p_rints_x_ecoeffs[19] += rints[1] * p_ecoeffs_cd_tsp[40];
                    p_rints_x_ecoeffs[19] += rints[4] * p_ecoeffs_cd_tsp[103];
                    p_rints_x_ecoeffs[19] += rints[26] * p_ecoeffs_cd_tsp[565];
                    p_rints_x_ecoeffs[20] += rints[5] * p_ecoeffs_cd_tsp[125];
                    p_rints_x_ecoeffs[20] += rints[2] * p_ecoeffs_cd_tsp[62];
                    p_rints_x_ecoeffs[20] += rints[13] * p_ecoeffs_cd_tsp[293];
                    p_rints_x_ecoeffs[20] += rints[7] * p_ecoeffs_cd_tsp[167];
                    p_rints_x_ecoeffs[20] += rints[11] * p_ecoeffs_cd_tsp[251];
                    p_rints_x_ecoeffs[20] += rints[0] * p_ecoeffs_cd_tsp[20];
                    p_rints_x_ecoeffs[20] += rints[16] * p_ecoeffs_cd_tsp[356];
                    p_rints_x_ecoeffs[20] += rints[30] * p_ecoeffs_cd_tsp[650];
                    p_rints_x_ecoeffs[20] += rints[1] * p_ecoeffs_cd_tsp[41];
                    p_rints_x_ecoeffs[20] += rints[4] * p_ecoeffs_cd_tsp[104];
                    p_rints_x_ecoeffs[20] += rints[23] * p_ecoeffs_cd_tsp[503];
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
            eri4_batch[15] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[15];
            eri4_batch[16] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[16];
            eri4_batch[17] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[17];
            eri4_batch[18] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[18];
            eri4_batch[19] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[19];
            eri4_batch[20] += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[20];
        }
}

template<> void lible::ints::two::eri4Kernel<0, 0, 4, 0>(const int cdepth_a, const int cdepth_b,
                                                         const int cdepth_c, const int cdepth_d,
                                                         const double *exps_a, const double *exps_b,
                                                         const double *exps_c, const double *exps_d,
                                                         const double *coords_a, const double *coords_b,
                                                         const double *coords_c, const double *coords_d,
                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,
                                                         double *eri4_batch)
{
    constexpr int la = 0, lb = 0, lc = 4, ld = 0;
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

                    p_rints_x_ecoeffs[0] += rints[25] * p_ecoeffs_cd_tsp[225];
                    p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_cd_tsp[18];
                    p_rints_x_ecoeffs[0] += rints[11] * p_ecoeffs_cd_tsp[99];
                    p_rints_x_ecoeffs[0] += rints[16] * p_ecoeffs_cd_tsp[144];
                    p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_cd_tsp[27];
                    p_rints_x_ecoeffs[0] += rints[34] * p_ecoeffs_cd_tsp[306];
                    p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_cd_tsp[153];
                    p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_cd_tsp[54];
                    p_rints_x_ecoeffs[0] += rints[5] * p_ecoeffs_cd_tsp[45];
                    p_rints_x_ecoeffs[0] += rints[10] * p_ecoeffs_cd_tsp[90];
                    p_rints_x_ecoeffs[0] += rints[20] * p_ecoeffs_cd_tsp[180];
                    p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_cd_tsp[108];
                    p_rints_x_ecoeffs[0] += rints[18] * p_ecoeffs_cd_tsp[162];
                    p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_cd_tsp[171];
                    p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_cd_tsp[63];
                    p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_cd_tsp[0];
                    p_rints_x_ecoeffs[0] += rints[30] * p_ecoeffs_cd_tsp[270];
                    p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_cd_tsp[9];
                    p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_cd_tsp[36];
                    p_rints_x_ecoeffs[0] += rints[23] * p_ecoeffs_cd_tsp[207];
                    p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_cd_tsp[72];
                    p_rints_x_ecoeffs[0] += rints[13] * p_ecoeffs_cd_tsp[117];
                    p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_cd_tsp[81];
                    p_rints_x_ecoeffs[0] += rints[32] * p_ecoeffs_cd_tsp[288];
                    p_rints_x_ecoeffs[0] += rints[15] * p_ecoeffs_cd_tsp[135];
                    p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_cd_tsp[19];
                    p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_cd_tsp[28];
                    p_rints_x_ecoeffs[1] += rints[17] * p_ecoeffs_cd_tsp[154];
                    p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_cd_tsp[55];
                    p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_cd_tsp[46];
                    p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_cd_tsp[91];
                    p_rints_x_ecoeffs[1] += rints[12] * p_ecoeffs_cd_tsp[109];
                    p_rints_x_ecoeffs[1] += rints[27] * p_ecoeffs_cd_tsp[244];
                    p_rints_x_ecoeffs[1] += rints[19] * p_ecoeffs_cd_tsp[172];
                    p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_cd_tsp[64];
                    p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_cd_tsp[1];
                    p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_cd_tsp[10];
                    p_rints_x_ecoeffs[1] += rints[29] * p_ecoeffs_cd_tsp[262];
                    p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_cd_tsp[37];
                    p_rints_x_ecoeffs[1] += rints[8] * p_ecoeffs_cd_tsp[73];
                    p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_cd_tsp[118];
                    p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_cd_tsp[82];
                    p_rints_x_ecoeffs[1] += rints[22] * p_ecoeffs_cd_tsp[199];
                    p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_cd_tsp[136];
                    p_rints_x_ecoeffs[1] += rints[14] * p_ecoeffs_cd_tsp[127];
                    p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_cd_tsp[20];
                    p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_cd_tsp[101];
                    p_rints_x_ecoeffs[2] += rints[33] * p_ecoeffs_cd_tsp[299];
                    p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_cd_tsp[146];
                    p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_cd_tsp[29];
                    p_rints_x_ecoeffs[2] += rints[17] * p_ecoeffs_cd_tsp[155];
                    p_rints_x_ecoeffs[2] += rints[6] * p_ecoeffs_cd_tsp[56];
                    p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_cd_tsp[47];
                    p_rints_x_ecoeffs[2] += rints[12] * p_ecoeffs_cd_tsp[110];
                    p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_cd_tsp[164];
                    p_rints_x_ecoeffs[2] += rints[19] * p_ecoeffs_cd_tsp[173];
                    p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_cd_tsp[65];
                    p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_cd_tsp[2];
                    p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_cd_tsp[11];
                    p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_cd_tsp[38];
                    p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_cd_tsp[74];
                    p_rints_x_ecoeffs[2] += rints[24] * p_ecoeffs_cd_tsp[218];
                    p_rints_x_ecoeffs[2] += rints[31] * p_ecoeffs_cd_tsp[281];
                    p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_cd_tsp[83];
                    p_rints_x_ecoeffs[2] += rints[14] * p_ecoeffs_cd_tsp[128];
                    p_rints_x_ecoeffs[3] += rints[25] * p_ecoeffs_cd_tsp[228];
                    p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_cd_tsp[21];
                    p_rints_x_ecoeffs[3] += rints[11] * p_ecoeffs_cd_tsp[102];
                    p_rints_x_ecoeffs[3] += rints[16] * p_ecoeffs_cd_tsp[147];
                    p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_cd_tsp[30];
                    p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_cd_tsp[156];
                    p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_cd_tsp[57];
                    p_rints_x_ecoeffs[3] += rints[5] * p_ecoeffs_cd_tsp[48];
                    p_rints_x_ecoeffs[3] += rints[10] * p_ecoeffs_cd_tsp[93];
                    p_rints_x_ecoeffs[3] += rints[20] * p_ecoeffs_cd_tsp[183];
                    p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_cd_tsp[111];
                    p_rints_x_ecoeffs[3] += rints[18] * p_ecoeffs_cd_tsp[165];
                    p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_cd_tsp[66];
                    p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_cd_tsp[3];
                    p_rints_x_ecoeffs[3] += rints[30] * p_ecoeffs_cd_tsp[273];
                    p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_cd_tsp[12];
                    p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_cd_tsp[39];
                    p_rints_x_ecoeffs[3] += rints[23] * p_ecoeffs_cd_tsp[210];
                    p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_cd_tsp[75];
                    p_rints_x_ecoeffs[3] += rints[13] * p_ecoeffs_cd_tsp[120];
                    p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_cd_tsp[84];
                    p_rints_x_ecoeffs[3] += rints[32] * p_ecoeffs_cd_tsp[291];
                    p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_cd_tsp[138];
                    p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_cd_tsp[22];
                    p_rints_x_ecoeffs[4] += rints[11] * p_ecoeffs_cd_tsp[103];
                    p_rints_x_ecoeffs[4] += rints[16] * p_ecoeffs_cd_tsp[148];
                    p_rints_x_ecoeffs[4] += rints[26] * p_ecoeffs_cd_tsp[238];
                    p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_cd_tsp[31];
                    p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_cd_tsp[58];
                    p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_cd_tsp[49];
                    p_rints_x_ecoeffs[4] += rints[10] * p_ecoeffs_cd_tsp[94];
                    p_rints_x_ecoeffs[4] += rints[18] * p_ecoeffs_cd_tsp[166];
                    p_rints_x_ecoeffs[4] += rints[7] * p_ecoeffs_cd_tsp[67];
                    p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_cd_tsp[4];
                    p_rints_x_ecoeffs[4] += rints[28] * p_ecoeffs_cd_tsp[256];
                    p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_cd_tsp[13];
                    p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_cd_tsp[40];
                    p_rints_x_ecoeffs[4] += rints[8] * p_ecoeffs_cd_tsp[76];
                    p_rints_x_ecoeffs[4] += rints[13] * p_ecoeffs_cd_tsp[121];
                    p_rints_x_ecoeffs[4] += rints[21] * p_ecoeffs_cd_tsp[193];
                    p_rints_x_ecoeffs[4] += rints[9] * p_ecoeffs_cd_tsp[85];
                    p_rints_x_ecoeffs[4] += rints[15] * p_ecoeffs_cd_tsp[139];
                    p_rints_x_ecoeffs[4] += rints[14] * p_ecoeffs_cd_tsp[130];
                    p_rints_x_ecoeffs[5] += rints[27] * p_ecoeffs_cd_tsp[248];
                    p_rints_x_ecoeffs[5] += rints[17] * p_ecoeffs_cd_tsp[158];
                    p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_cd_tsp[59];
                    p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_cd_tsp[50];
                    p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_cd_tsp[23];
                    p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_cd_tsp[122];
                    p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_cd_tsp[68];
                    p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_cd_tsp[95];
                    p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_cd_tsp[5];
                    p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_cd_tsp[14];
                    p_rints_x_ecoeffs[5] += rints[12] * p_ecoeffs_cd_tsp[113];
                    p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_cd_tsp[41];
                    p_rints_x_ecoeffs[5] += rints[22] * p_ecoeffs_cd_tsp[203];
                    p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_cd_tsp[32];
                    p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_cd_tsp[131];
                    p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_cd_tsp[77];
                    p_rints_x_ecoeffs[6] += rints[24] * p_ecoeffs_cd_tsp[222];
                    p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_cd_tsp[60];
                    p_rints_x_ecoeffs[6] += rints[17] * p_ecoeffs_cd_tsp[159];
                    p_rints_x_ecoeffs[6] += rints[5] * p_ecoeffs_cd_tsp[51];
                    p_rints_x_ecoeffs[6] += rints[31] * p_ecoeffs_cd_tsp[285];
                    p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_cd_tsp[24];
                    p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_cd_tsp[105];
                    p_rints_x_ecoeffs[6] += rints[7] * p_ecoeffs_cd_tsp[69];
                    p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_cd_tsp[6];
                    p_rints_x_ecoeffs[6] += rints[16] * p_ecoeffs_cd_tsp[150];
                    p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_cd_tsp[15];
                    p_rints_x_ecoeffs[6] += rints[12] * p_ecoeffs_cd_tsp[114];
                    p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_cd_tsp[42];
                    p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_cd_tsp[33];
                    p_rints_x_ecoeffs[6] += rints[14] * p_ecoeffs_cd_tsp[132];
                    p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_cd_tsp[78];
                    p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_cd_tsp[52];
                    p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_cd_tsp[25];
                    p_rints_x_ecoeffs[7] += rints[13] * p_ecoeffs_cd_tsp[124];
                    p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_cd_tsp[70];
                    p_rints_x_ecoeffs[7] += rints[11] * p_ecoeffs_cd_tsp[106];
                    p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_cd_tsp[7];
                    p_rints_x_ecoeffs[7] += rints[16] * p_ecoeffs_cd_tsp[151];
                    p_rints_x_ecoeffs[7] += rints[10] * p_ecoeffs_cd_tsp[97];
                    p_rints_x_ecoeffs[7] += rints[30] * p_ecoeffs_cd_tsp[277];
                    p_rints_x_ecoeffs[7] += rints[20] * p_ecoeffs_cd_tsp[187];
                    p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_cd_tsp[16];
                    p_rints_x_ecoeffs[7] += rints[4] * p_ecoeffs_cd_tsp[43];
                    p_rints_x_ecoeffs[7] += rints[23] * p_ecoeffs_cd_tsp[214];
                    p_rints_x_ecoeffs[8] += rints[5] * p_ecoeffs_cd_tsp[53];
                    p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_cd_tsp[26];
                    p_rints_x_ecoeffs[8] += rints[13] * p_ecoeffs_cd_tsp[125];
                    p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_cd_tsp[71];
                    p_rints_x_ecoeffs[8] += rints[11] * p_ecoeffs_cd_tsp[107];
                    p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_cd_tsp[8];
                    p_rints_x_ecoeffs[8] += rints[16] * p_ecoeffs_cd_tsp[152];
                    p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_cd_tsp[98];
                    p_rints_x_ecoeffs[8] += rints[21] * p_ecoeffs_cd_tsp[197];
                    p_rints_x_ecoeffs[8] += rints[1] * p_ecoeffs_cd_tsp[17];
                    p_rints_x_ecoeffs[8] += rints[4] * p_ecoeffs_cd_tsp[44];
                    p_rints_x_ecoeffs[8] += rints[26] * p_ecoeffs_cd_tsp[242];
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
        }
}

template void lible::ints::two::eri3Kernel<0, 0, 4>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<0, 4>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);

