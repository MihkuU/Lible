#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<0, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template<> lible::vec4d
lible::ints::two::eri4Kernel<0, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
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
                    calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 0, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 0, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 0, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 1, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 1, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 1, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 0, 1, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 2, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 0, 2, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 2, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 0, 2, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 0, 3, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 0, 3, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 0, 3, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 0, 3, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 0, 3, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 0, 4, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[20];
            eri4_batch(0, 0, 4, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[21];
            eri4_batch(0, 0, 4, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[22];
            eri4_batch(0, 0, 4, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[23];
            eri4_batch(0, 0, 4, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[24];
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

template lible::vec4d lible::ints::eri4KernelFun<0, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 0, 1, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template<> lible::vec4d
lible::ints::two::eri4Kernel<0, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
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
                    calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 0, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 0, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 1, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 1, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 2, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 2, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[8];
            eri4_batch(0, 0, 3, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[9];
            eri4_batch(0, 0, 3, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[10];
            eri4_batch(0, 0, 3, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[11];
            eri4_batch(0, 0, 4, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[12];
            eri4_batch(0, 0, 4, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[13];
            eri4_batch(0, 0, 4, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[14];
            eri4_batch(0, 0, 5, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[15];
            eri4_batch(0, 0, 5, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[16];
            eri4_batch(0, 0, 5, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[17];
            eri4_batch(0, 0, 6, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[18];
            eri4_batch(0, 0, 6, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[19];
            eri4_batch(0, 0, 6, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[20];
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

template lible::vec4d lible::ints::eri4KernelFun<0, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 0, 0, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template<> lible::vec4d
lible::ints::two::eri4Kernel<0, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
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
                    calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                    const double* p_ecoeffs_cd_tsp = &pecoeffs_cd_tsp[icd * n_ecoeffs_cd];
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

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri4_batch(0, 0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri4_batch(0, 0, 1, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri4_batch(0, 0, 2, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri4_batch(0, 0, 3, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri4_batch(0, 0, 4, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri4_batch(0, 0, 5, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri4_batch(0, 0, 6, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
            eri4_batch(0, 0, 7, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[7];
            eri4_batch(0, 0, 8, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[8];
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

template lible::vec3d lible::ints::eri3KernelFun<0, 0, 4>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<0, 0, 4>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template<> lible::vec3d
lible::ints::two::eri3Kernel<0, 0, 4>(const int ipair_ab, const int ishell_c,
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

    constexpr int la = 0, lb = 0, lc = 4;
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
                calcRInts_ERI<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                const double* p_ecoeffs_c = &pecoeffs_c[ic * n_ecoeffs_c];
                double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

                p_rints_x_ecoeffs[0] += rints[25] * p_ecoeffs_c[25];
                p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_c[2];
                p_rints_x_ecoeffs[0] += rints[11] * p_ecoeffs_c[11];
                p_rints_x_ecoeffs[0] += rints[16] * p_ecoeffs_c[16];
                p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_c[3];
                p_rints_x_ecoeffs[0] += rints[34] * p_ecoeffs_c[34];
                p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_c[17];
                p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_c[6];
                p_rints_x_ecoeffs[0] += rints[5] * p_ecoeffs_c[5];
                p_rints_x_ecoeffs[0] += rints[10] * p_ecoeffs_c[10];
                p_rints_x_ecoeffs[0] += rints[20] * p_ecoeffs_c[20];
                p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_c[12];
                p_rints_x_ecoeffs[0] += rints[18] * p_ecoeffs_c[18];
                p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_c[19];
                p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_c[7];
                p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_c[0];
                p_rints_x_ecoeffs[0] += rints[30] * p_ecoeffs_c[30];
                p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_c[1];
                p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_c[4];
                p_rints_x_ecoeffs[0] += rints[23] * p_ecoeffs_c[23];
                p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_c[8];
                p_rints_x_ecoeffs[0] += rints[13] * p_ecoeffs_c[13];
                p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_c[9];
                p_rints_x_ecoeffs[0] += rints[32] * p_ecoeffs_c[32];
                p_rints_x_ecoeffs[0] += rints[15] * p_ecoeffs_c[15];
                p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_c[37];
                p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_c[38];
                p_rints_x_ecoeffs[1] += rints[17] * p_ecoeffs_c[52];
                p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_c[41];
                p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_c[40];
                p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_c[45];
                p_rints_x_ecoeffs[1] += rints[12] * p_ecoeffs_c[47];
                p_rints_x_ecoeffs[1] += rints[27] * p_ecoeffs_c[62];
                p_rints_x_ecoeffs[1] += rints[19] * p_ecoeffs_c[54];
                p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_c[42];
                p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_c[35];
                p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_c[36];
                p_rints_x_ecoeffs[1] += rints[29] * p_ecoeffs_c[64];
                p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_c[39];
                p_rints_x_ecoeffs[1] += rints[8] * p_ecoeffs_c[43];
                p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_c[48];
                p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_c[44];
                p_rints_x_ecoeffs[1] += rints[22] * p_ecoeffs_c[57];
                p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_c[50];
                p_rints_x_ecoeffs[1] += rints[14] * p_ecoeffs_c[49];
                p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_c[72];
                p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_c[81];
                p_rints_x_ecoeffs[2] += rints[33] * p_ecoeffs_c[103];
                p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_c[86];
                p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_c[73];
                p_rints_x_ecoeffs[2] += rints[17] * p_ecoeffs_c[87];
                p_rints_x_ecoeffs[2] += rints[6] * p_ecoeffs_c[76];
                p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_c[75];
                p_rints_x_ecoeffs[2] += rints[12] * p_ecoeffs_c[82];
                p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_c[88];
                p_rints_x_ecoeffs[2] += rints[19] * p_ecoeffs_c[89];
                p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_c[77];
                p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_c[70];
                p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_c[71];
                p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_c[74];
                p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_c[78];
                p_rints_x_ecoeffs[2] += rints[24] * p_ecoeffs_c[94];
                p_rints_x_ecoeffs[2] += rints[31] * p_ecoeffs_c[101];
                p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_c[79];
                p_rints_x_ecoeffs[2] += rints[14] * p_ecoeffs_c[84];
                p_rints_x_ecoeffs[3] += rints[25] * p_ecoeffs_c[130];
                p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_c[107];
                p_rints_x_ecoeffs[3] += rints[11] * p_ecoeffs_c[116];
                p_rints_x_ecoeffs[3] += rints[16] * p_ecoeffs_c[121];
                p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_c[108];
                p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_c[122];
                p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_c[111];
                p_rints_x_ecoeffs[3] += rints[5] * p_ecoeffs_c[110];
                p_rints_x_ecoeffs[3] += rints[10] * p_ecoeffs_c[115];
                p_rints_x_ecoeffs[3] += rints[20] * p_ecoeffs_c[125];
                p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_c[117];
                p_rints_x_ecoeffs[3] += rints[18] * p_ecoeffs_c[123];
                p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_c[112];
                p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_c[105];
                p_rints_x_ecoeffs[3] += rints[30] * p_ecoeffs_c[135];
                p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_c[106];
                p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_c[109];
                p_rints_x_ecoeffs[3] += rints[23] * p_ecoeffs_c[128];
                p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_c[113];
                p_rints_x_ecoeffs[3] += rints[13] * p_ecoeffs_c[118];
                p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_c[114];
                p_rints_x_ecoeffs[3] += rints[32] * p_ecoeffs_c[137];
                p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_c[120];
                p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_c[142];
                p_rints_x_ecoeffs[4] += rints[11] * p_ecoeffs_c[151];
                p_rints_x_ecoeffs[4] += rints[16] * p_ecoeffs_c[156];
                p_rints_x_ecoeffs[4] += rints[26] * p_ecoeffs_c[166];
                p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_c[143];
                p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_c[146];
                p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_c[145];
                p_rints_x_ecoeffs[4] += rints[10] * p_ecoeffs_c[150];
                p_rints_x_ecoeffs[4] += rints[18] * p_ecoeffs_c[158];
                p_rints_x_ecoeffs[4] += rints[7] * p_ecoeffs_c[147];
                p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_c[140];
                p_rints_x_ecoeffs[4] += rints[28] * p_ecoeffs_c[168];
                p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_c[141];
                p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_c[144];
                p_rints_x_ecoeffs[4] += rints[8] * p_ecoeffs_c[148];
                p_rints_x_ecoeffs[4] += rints[13] * p_ecoeffs_c[153];
                p_rints_x_ecoeffs[4] += rints[21] * p_ecoeffs_c[161];
                p_rints_x_ecoeffs[4] += rints[9] * p_ecoeffs_c[149];
                p_rints_x_ecoeffs[4] += rints[15] * p_ecoeffs_c[155];
                p_rints_x_ecoeffs[4] += rints[14] * p_ecoeffs_c[154];
                p_rints_x_ecoeffs[5] += rints[27] * p_ecoeffs_c[202];
                p_rints_x_ecoeffs[5] += rints[17] * p_ecoeffs_c[192];
                p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_c[181];
                p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_c[180];
                p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_c[177];
                p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_c[188];
                p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_c[182];
                p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_c[185];
                p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_c[175];
                p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_c[176];
                p_rints_x_ecoeffs[5] += rints[12] * p_ecoeffs_c[187];
                p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_c[179];
                p_rints_x_ecoeffs[5] += rints[22] * p_ecoeffs_c[197];
                p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_c[178];
                p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_c[189];
                p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_c[183];
                p_rints_x_ecoeffs[6] += rints[24] * p_ecoeffs_c[234];
                p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_c[216];
                p_rints_x_ecoeffs[6] += rints[17] * p_ecoeffs_c[227];
                p_rints_x_ecoeffs[6] += rints[5] * p_ecoeffs_c[215];
                p_rints_x_ecoeffs[6] += rints[31] * p_ecoeffs_c[241];
                p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_c[212];
                p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_c[221];
                p_rints_x_ecoeffs[6] += rints[7] * p_ecoeffs_c[217];
                p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_c[210];
                p_rints_x_ecoeffs[6] += rints[16] * p_ecoeffs_c[226];
                p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_c[211];
                p_rints_x_ecoeffs[6] += rints[12] * p_ecoeffs_c[222];
                p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_c[214];
                p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_c[213];
                p_rints_x_ecoeffs[6] += rints[14] * p_ecoeffs_c[224];
                p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_c[218];
                p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_c[250];
                p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_c[247];
                p_rints_x_ecoeffs[7] += rints[13] * p_ecoeffs_c[258];
                p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_c[252];
                p_rints_x_ecoeffs[7] += rints[11] * p_ecoeffs_c[256];
                p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_c[245];
                p_rints_x_ecoeffs[7] += rints[16] * p_ecoeffs_c[261];
                p_rints_x_ecoeffs[7] += rints[10] * p_ecoeffs_c[255];
                p_rints_x_ecoeffs[7] += rints[30] * p_ecoeffs_c[275];
                p_rints_x_ecoeffs[7] += rints[20] * p_ecoeffs_c[265];
                p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_c[246];
                p_rints_x_ecoeffs[7] += rints[4] * p_ecoeffs_c[249];
                p_rints_x_ecoeffs[7] += rints[23] * p_ecoeffs_c[268];
                p_rints_x_ecoeffs[8] += rints[5] * p_ecoeffs_c[285];
                p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_c[282];
                p_rints_x_ecoeffs[8] += rints[13] * p_ecoeffs_c[293];
                p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_c[287];
                p_rints_x_ecoeffs[8] += rints[11] * p_ecoeffs_c[291];
                p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_c[280];
                p_rints_x_ecoeffs[8] += rints[16] * p_ecoeffs_c[296];
                p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_c[290];
                p_rints_x_ecoeffs[8] += rints[21] * p_ecoeffs_c[301];
                p_rints_x_ecoeffs[8] += rints[1] * p_ecoeffs_c[281];
                p_rints_x_ecoeffs[8] += rints[4] * p_ecoeffs_c[284];
                p_rints_x_ecoeffs[8] += rints[26] * p_ecoeffs_c[306];
            }
}

    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            const double* p_ecoeffs_ab = &pecoeffs_ab[iab * n_ecoeffs_ab];
            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];

            eri3_batch(0, 0, 0) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[0];
            eri3_batch(0, 0, 1) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[1];
            eri3_batch(0, 0, 2) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[2];
            eri3_batch(0, 0, 3) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[3];
            eri3_batch(0, 0, 4) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[4];
            eri3_batch(0, 0, 5) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[5];
            eri3_batch(0, 0, 6) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[6];
            eri3_batch(0, 0, 7) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[7];
            eri3_batch(0, 0, 8) += p_ecoeffs_ab[0] * p_rints_x_ecoeffs[8];
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

template lible::vec2d lible::ints::eri2KernelFun<0, 4>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<0, 4>(const int ishell_a, const int ishell_b,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b,
                                                                        const ERI2D1Kernel *eri2d1_kernel);

template<> lible::vec2d
lible::ints::two::eri2Kernel<0, 4>(const int ishell_a, const int ishell_b,
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

    constexpr int la = 0, lb = 4;
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
            calcRInts_ERI<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

            const double* p_ecoeffs_b = &pecoeffs_b_tsp[ib * n_ecoeffs_b];
            double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];

            p_rints_x_ecoeffs[0] += rints[25] * p_ecoeffs_b[225];
            p_rints_x_ecoeffs[0] += rints[2] * p_ecoeffs_b[18];
            p_rints_x_ecoeffs[0] += rints[11] * p_ecoeffs_b[99];
            p_rints_x_ecoeffs[0] += rints[16] * p_ecoeffs_b[144];
            p_rints_x_ecoeffs[0] += rints[3] * p_ecoeffs_b[27];
            p_rints_x_ecoeffs[0] += rints[34] * p_ecoeffs_b[306];
            p_rints_x_ecoeffs[0] += rints[17] * p_ecoeffs_b[153];
            p_rints_x_ecoeffs[0] += rints[6] * p_ecoeffs_b[54];
            p_rints_x_ecoeffs[0] += rints[5] * p_ecoeffs_b[45];
            p_rints_x_ecoeffs[0] += rints[10] * p_ecoeffs_b[90];
            p_rints_x_ecoeffs[0] += rints[20] * p_ecoeffs_b[180];
            p_rints_x_ecoeffs[0] += rints[12] * p_ecoeffs_b[108];
            p_rints_x_ecoeffs[0] += rints[18] * p_ecoeffs_b[162];
            p_rints_x_ecoeffs[0] += rints[19] * p_ecoeffs_b[171];
            p_rints_x_ecoeffs[0] += rints[7] * p_ecoeffs_b[63];
            p_rints_x_ecoeffs[0] += rints[0] * p_ecoeffs_b[0];
            p_rints_x_ecoeffs[0] += rints[30] * p_ecoeffs_b[270];
            p_rints_x_ecoeffs[0] += rints[1] * p_ecoeffs_b[9];
            p_rints_x_ecoeffs[0] += rints[4] * p_ecoeffs_b[36];
            p_rints_x_ecoeffs[0] += rints[23] * p_ecoeffs_b[207];
            p_rints_x_ecoeffs[0] += rints[8] * p_ecoeffs_b[72];
            p_rints_x_ecoeffs[0] += rints[13] * p_ecoeffs_b[117];
            p_rints_x_ecoeffs[0] += rints[9] * p_ecoeffs_b[81];
            p_rints_x_ecoeffs[0] += rints[32] * p_ecoeffs_b[288];
            p_rints_x_ecoeffs[0] += rints[15] * p_ecoeffs_b[135];
            p_rints_x_ecoeffs[1] += rints[2] * p_ecoeffs_b[19];
            p_rints_x_ecoeffs[1] += rints[3] * p_ecoeffs_b[28];
            p_rints_x_ecoeffs[1] += rints[17] * p_ecoeffs_b[154];
            p_rints_x_ecoeffs[1] += rints[6] * p_ecoeffs_b[55];
            p_rints_x_ecoeffs[1] += rints[5] * p_ecoeffs_b[46];
            p_rints_x_ecoeffs[1] += rints[10] * p_ecoeffs_b[91];
            p_rints_x_ecoeffs[1] += rints[12] * p_ecoeffs_b[109];
            p_rints_x_ecoeffs[1] += rints[27] * p_ecoeffs_b[244];
            p_rints_x_ecoeffs[1] += rints[19] * p_ecoeffs_b[172];
            p_rints_x_ecoeffs[1] += rints[7] * p_ecoeffs_b[64];
            p_rints_x_ecoeffs[1] += rints[0] * p_ecoeffs_b[1];
            p_rints_x_ecoeffs[1] += rints[1] * p_ecoeffs_b[10];
            p_rints_x_ecoeffs[1] += rints[29] * p_ecoeffs_b[262];
            p_rints_x_ecoeffs[1] += rints[4] * p_ecoeffs_b[37];
            p_rints_x_ecoeffs[1] += rints[8] * p_ecoeffs_b[73];
            p_rints_x_ecoeffs[1] += rints[13] * p_ecoeffs_b[118];
            p_rints_x_ecoeffs[1] += rints[9] * p_ecoeffs_b[82];
            p_rints_x_ecoeffs[1] += rints[22] * p_ecoeffs_b[199];
            p_rints_x_ecoeffs[1] += rints[15] * p_ecoeffs_b[136];
            p_rints_x_ecoeffs[1] += rints[14] * p_ecoeffs_b[127];
            p_rints_x_ecoeffs[2] += rints[2] * p_ecoeffs_b[20];
            p_rints_x_ecoeffs[2] += rints[11] * p_ecoeffs_b[101];
            p_rints_x_ecoeffs[2] += rints[33] * p_ecoeffs_b[299];
            p_rints_x_ecoeffs[2] += rints[16] * p_ecoeffs_b[146];
            p_rints_x_ecoeffs[2] += rints[3] * p_ecoeffs_b[29];
            p_rints_x_ecoeffs[2] += rints[17] * p_ecoeffs_b[155];
            p_rints_x_ecoeffs[2] += rints[6] * p_ecoeffs_b[56];
            p_rints_x_ecoeffs[2] += rints[5] * p_ecoeffs_b[47];
            p_rints_x_ecoeffs[2] += rints[12] * p_ecoeffs_b[110];
            p_rints_x_ecoeffs[2] += rints[18] * p_ecoeffs_b[164];
            p_rints_x_ecoeffs[2] += rints[19] * p_ecoeffs_b[173];
            p_rints_x_ecoeffs[2] += rints[7] * p_ecoeffs_b[65];
            p_rints_x_ecoeffs[2] += rints[0] * p_ecoeffs_b[2];
            p_rints_x_ecoeffs[2] += rints[1] * p_ecoeffs_b[11];
            p_rints_x_ecoeffs[2] += rints[4] * p_ecoeffs_b[38];
            p_rints_x_ecoeffs[2] += rints[8] * p_ecoeffs_b[74];
            p_rints_x_ecoeffs[2] += rints[24] * p_ecoeffs_b[218];
            p_rints_x_ecoeffs[2] += rints[31] * p_ecoeffs_b[281];
            p_rints_x_ecoeffs[2] += rints[9] * p_ecoeffs_b[83];
            p_rints_x_ecoeffs[2] += rints[14] * p_ecoeffs_b[128];
            p_rints_x_ecoeffs[3] += rints[25] * p_ecoeffs_b[228];
            p_rints_x_ecoeffs[3] += rints[2] * p_ecoeffs_b[21];
            p_rints_x_ecoeffs[3] += rints[11] * p_ecoeffs_b[102];
            p_rints_x_ecoeffs[3] += rints[16] * p_ecoeffs_b[147];
            p_rints_x_ecoeffs[3] += rints[3] * p_ecoeffs_b[30];
            p_rints_x_ecoeffs[3] += rints[17] * p_ecoeffs_b[156];
            p_rints_x_ecoeffs[3] += rints[6] * p_ecoeffs_b[57];
            p_rints_x_ecoeffs[3] += rints[5] * p_ecoeffs_b[48];
            p_rints_x_ecoeffs[3] += rints[10] * p_ecoeffs_b[93];
            p_rints_x_ecoeffs[3] += rints[20] * p_ecoeffs_b[183];
            p_rints_x_ecoeffs[3] += rints[12] * p_ecoeffs_b[111];
            p_rints_x_ecoeffs[3] += rints[18] * p_ecoeffs_b[165];
            p_rints_x_ecoeffs[3] += rints[7] * p_ecoeffs_b[66];
            p_rints_x_ecoeffs[3] += rints[0] * p_ecoeffs_b[3];
            p_rints_x_ecoeffs[3] += rints[30] * p_ecoeffs_b[273];
            p_rints_x_ecoeffs[3] += rints[1] * p_ecoeffs_b[12];
            p_rints_x_ecoeffs[3] += rints[4] * p_ecoeffs_b[39];
            p_rints_x_ecoeffs[3] += rints[23] * p_ecoeffs_b[210];
            p_rints_x_ecoeffs[3] += rints[8] * p_ecoeffs_b[75];
            p_rints_x_ecoeffs[3] += rints[13] * p_ecoeffs_b[120];
            p_rints_x_ecoeffs[3] += rints[9] * p_ecoeffs_b[84];
            p_rints_x_ecoeffs[3] += rints[32] * p_ecoeffs_b[291];
            p_rints_x_ecoeffs[3] += rints[15] * p_ecoeffs_b[138];
            p_rints_x_ecoeffs[4] += rints[2] * p_ecoeffs_b[22];
            p_rints_x_ecoeffs[4] += rints[11] * p_ecoeffs_b[103];
            p_rints_x_ecoeffs[4] += rints[16] * p_ecoeffs_b[148];
            p_rints_x_ecoeffs[4] += rints[26] * p_ecoeffs_b[238];
            p_rints_x_ecoeffs[4] += rints[3] * p_ecoeffs_b[31];
            p_rints_x_ecoeffs[4] += rints[6] * p_ecoeffs_b[58];
            p_rints_x_ecoeffs[4] += rints[5] * p_ecoeffs_b[49];
            p_rints_x_ecoeffs[4] += rints[10] * p_ecoeffs_b[94];
            p_rints_x_ecoeffs[4] += rints[18] * p_ecoeffs_b[166];
            p_rints_x_ecoeffs[4] += rints[7] * p_ecoeffs_b[67];
            p_rints_x_ecoeffs[4] += rints[0] * p_ecoeffs_b[4];
            p_rints_x_ecoeffs[4] += rints[28] * p_ecoeffs_b[256];
            p_rints_x_ecoeffs[4] += rints[1] * p_ecoeffs_b[13];
            p_rints_x_ecoeffs[4] += rints[4] * p_ecoeffs_b[40];
            p_rints_x_ecoeffs[4] += rints[8] * p_ecoeffs_b[76];
            p_rints_x_ecoeffs[4] += rints[13] * p_ecoeffs_b[121];
            p_rints_x_ecoeffs[4] += rints[21] * p_ecoeffs_b[193];
            p_rints_x_ecoeffs[4] += rints[9] * p_ecoeffs_b[85];
            p_rints_x_ecoeffs[4] += rints[15] * p_ecoeffs_b[139];
            p_rints_x_ecoeffs[4] += rints[14] * p_ecoeffs_b[130];
            p_rints_x_ecoeffs[5] += rints[27] * p_ecoeffs_b[248];
            p_rints_x_ecoeffs[5] += rints[17] * p_ecoeffs_b[158];
            p_rints_x_ecoeffs[5] += rints[6] * p_ecoeffs_b[59];
            p_rints_x_ecoeffs[5] += rints[5] * p_ecoeffs_b[50];
            p_rints_x_ecoeffs[5] += rints[2] * p_ecoeffs_b[23];
            p_rints_x_ecoeffs[5] += rints[13] * p_ecoeffs_b[122];
            p_rints_x_ecoeffs[5] += rints[7] * p_ecoeffs_b[68];
            p_rints_x_ecoeffs[5] += rints[10] * p_ecoeffs_b[95];
            p_rints_x_ecoeffs[5] += rints[0] * p_ecoeffs_b[5];
            p_rints_x_ecoeffs[5] += rints[1] * p_ecoeffs_b[14];
            p_rints_x_ecoeffs[5] += rints[12] * p_ecoeffs_b[113];
            p_rints_x_ecoeffs[5] += rints[4] * p_ecoeffs_b[41];
            p_rints_x_ecoeffs[5] += rints[22] * p_ecoeffs_b[203];
            p_rints_x_ecoeffs[5] += rints[3] * p_ecoeffs_b[32];
            p_rints_x_ecoeffs[5] += rints[14] * p_ecoeffs_b[131];
            p_rints_x_ecoeffs[5] += rints[8] * p_ecoeffs_b[77];
            p_rints_x_ecoeffs[6] += rints[24] * p_ecoeffs_b[222];
            p_rints_x_ecoeffs[6] += rints[6] * p_ecoeffs_b[60];
            p_rints_x_ecoeffs[6] += rints[17] * p_ecoeffs_b[159];
            p_rints_x_ecoeffs[6] += rints[5] * p_ecoeffs_b[51];
            p_rints_x_ecoeffs[6] += rints[31] * p_ecoeffs_b[285];
            p_rints_x_ecoeffs[6] += rints[2] * p_ecoeffs_b[24];
            p_rints_x_ecoeffs[6] += rints[11] * p_ecoeffs_b[105];
            p_rints_x_ecoeffs[6] += rints[7] * p_ecoeffs_b[69];
            p_rints_x_ecoeffs[6] += rints[0] * p_ecoeffs_b[6];
            p_rints_x_ecoeffs[6] += rints[16] * p_ecoeffs_b[150];
            p_rints_x_ecoeffs[6] += rints[1] * p_ecoeffs_b[15];
            p_rints_x_ecoeffs[6] += rints[12] * p_ecoeffs_b[114];
            p_rints_x_ecoeffs[6] += rints[4] * p_ecoeffs_b[42];
            p_rints_x_ecoeffs[6] += rints[3] * p_ecoeffs_b[33];
            p_rints_x_ecoeffs[6] += rints[14] * p_ecoeffs_b[132];
            p_rints_x_ecoeffs[6] += rints[8] * p_ecoeffs_b[78];
            p_rints_x_ecoeffs[7] += rints[5] * p_ecoeffs_b[52];
            p_rints_x_ecoeffs[7] += rints[2] * p_ecoeffs_b[25];
            p_rints_x_ecoeffs[7] += rints[13] * p_ecoeffs_b[124];
            p_rints_x_ecoeffs[7] += rints[7] * p_ecoeffs_b[70];
            p_rints_x_ecoeffs[7] += rints[11] * p_ecoeffs_b[106];
            p_rints_x_ecoeffs[7] += rints[0] * p_ecoeffs_b[7];
            p_rints_x_ecoeffs[7] += rints[16] * p_ecoeffs_b[151];
            p_rints_x_ecoeffs[7] += rints[10] * p_ecoeffs_b[97];
            p_rints_x_ecoeffs[7] += rints[30] * p_ecoeffs_b[277];
            p_rints_x_ecoeffs[7] += rints[20] * p_ecoeffs_b[187];
            p_rints_x_ecoeffs[7] += rints[1] * p_ecoeffs_b[16];
            p_rints_x_ecoeffs[7] += rints[4] * p_ecoeffs_b[43];
            p_rints_x_ecoeffs[7] += rints[23] * p_ecoeffs_b[214];
            p_rints_x_ecoeffs[8] += rints[5] * p_ecoeffs_b[53];
            p_rints_x_ecoeffs[8] += rints[2] * p_ecoeffs_b[26];
            p_rints_x_ecoeffs[8] += rints[13] * p_ecoeffs_b[125];
            p_rints_x_ecoeffs[8] += rints[7] * p_ecoeffs_b[71];
            p_rints_x_ecoeffs[8] += rints[11] * p_ecoeffs_b[107];
            p_rints_x_ecoeffs[8] += rints[0] * p_ecoeffs_b[8];
            p_rints_x_ecoeffs[8] += rints[16] * p_ecoeffs_b[152];
            p_rints_x_ecoeffs[8] += rints[10] * p_ecoeffs_b[98];
            p_rints_x_ecoeffs[8] += rints[21] * p_ecoeffs_b[197];
            p_rints_x_ecoeffs[8] += rints[1] * p_ecoeffs_b[17];
            p_rints_x_ecoeffs[8] += rints[4] * p_ecoeffs_b[44];
            p_rints_x_ecoeffs[8] += rints[26] * p_ecoeffs_b[242];
        }
    }

    vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);
    for (int ia = 0; ia < cdepth_a; ia++)
    {
        const double* p_ecoeffs_a = &pecoeffs_a[ia * n_ecoeffs_a];
        const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];

        eri2_batch(0, 0) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[0];
        eri2_batch(0, 1) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[1];
        eri2_batch(0, 2) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[2];
        eri2_batch(0, 3) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[3];
        eri2_batch(0, 4) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[4];
        eri2_batch(0, 5) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[5];
        eri2_batch(0, 6) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[6];
        eri2_batch(0, 7) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[7];
        eri2_batch(0, 8) += p_ecoeffs_a[0] * p_rints_x_ecoeffs[8];
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

template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<0, 4>(const int ishell_a, const int ishell_b,
                                                           const std::vector<double> &ecoeffs_a,
                                                           const std::vector<double> &ecoeffs_b_tsp,
                                                           const ShellData &sh_data_a,
                                                           const ShellData &sh_data_b);

