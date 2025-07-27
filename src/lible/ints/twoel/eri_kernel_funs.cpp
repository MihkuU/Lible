#include <lible/ints/rints.hpp>
#include <lible/ints/twoel/eri_kernel_funs.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

lible::vec4d LI::eri4KernelFun(const int ipair_ab, const int ipair_cd,
                               const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                               const ERI4Kernel *eri4_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sp_data_cd.la;
    const int ld = sp_data_cd.lb;

    const int lab = la + lb;
    const int lcd = lc + ld;
    const int labcd = lab + lcd;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_sph_d = numSphericals(ld);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_cd = n_sph_c * n_sph_d;
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_cd = numHermites(lcd);
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
    const int ofs_E_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_E_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *ecoeffs_ab = &eri4_kernel->ecoeffs_bra[ofs_E_ab];
    const double *ecoeffs_cd = &eri4_kernel->ecoeffs_ket[ofs_E_cd];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lcd);

    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            vector<double> R_x_E(n_hermite_ab * n_sph_cd, 0);
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double a = sp_data_ab.exps[cofs_a + ia];
                    double b = sp_data_ab.exps[cofs_b + ib];
                    double c = sp_data_cd.exps[cofs_c + ic];
                    double d = sp_data_cd.exps[cofs_d + id];

                    double p = a + b;
                    double q = c + d;
                    double alpha = p * q / (p + q);

                    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                                           (a * xyz_a[2] + b * xyz_b[2]) / p};

                    array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                           (c * xyz_c[1] + d * xyz_d[1]) / q,
                                           (c * xyz_c[2] + d * xyz_d[2]) / q};

                    array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                            xyz_p[1] - xyz_q[1],
                                            xyz_p[2] - xyz_q[2]};

                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                    double x = alpha * (dx * dx + dy * dy + dz * dz);
                    vector<double> fnx = calcBoysF(labcd, x, eri4_kernel->boys_grid);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    vector<double> rints = calcRIntsMatrix(labcd, fac, alpha, &xyz_pq[0], &fnx[0],
                                                           hermite_idxs_bra, hermite_idxs_ket);

                    int ofs_e1_cd = icd * n_ecoeffs_cd;
                    shark_mm_ket(n_hermite_ab, n_sph_cd, n_hermite_cd, &rints[0], &ecoeffs_cd[ofs_e1_cd], &R_x_E[0]);
                }
            int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
            shark_mm_bra(n_sph_ab, n_sph_cd, n_hermite_ab, &ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri4_batch[0]);
        }

    return eri4_batch;
}

lible::vec3d LI::eri3KernelFun(const int ipair_ab, const int ishell_c,
                               const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                               const ERI3Kernel *eri3_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sh_data_c.l;
    const int lab = la + lb;
    const int labc = lab + lc;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_c = numHermites(lc);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_c = n_sph_c * n_hermite_c;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];
    const int ofs_E_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_E_c = sh_data_c.offsets_ecoeffs[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *ecoeffs_ab = &eri3_kernel->ecoeffs_bra[ofs_E_ab];
    const double *ecoeffs_c = &eri3_kernel->ecoeffs_ket[ofs_E_c];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lc);

    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            vector<double> R_x_E(n_hermite_ab * n_sph_c, 0);
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double c = exps_c[ic];

                double p = a + b;
                double alpha = p * c / (p + c);

                array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                            (a * coords_a[2] + b * coords_b[2]) / p};

                array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                             xyz_p[1] - coords_c[1],
                                             xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                vector<double> fnx = calcBoysF(labc, x, eri3_kernel->boys_grid);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                vector<double> rints = calcRIntsMatrix(labc, fac, alpha, &xyz_pc[0], &fnx[0],
                                                       hermite_idxs_bra, hermite_idxs_ket);

                int ofs_e0_c = ic * n_ecoeffs_c;
                shark_mm_ket(n_hermite_ab, n_sph_c, n_hermite_c, &rints[0], &ecoeffs_c[ofs_e0_c], &R_x_E[0]);
            }
            int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
            shark_mm_bra(n_sph_ab, n_sph_c, n_hermite_ab, &ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri3_batch[0]);
        }

    return eri3_batch;
}

lible::vec2d LI::eri2KernelFun(const int ishell_a, const int ishell_b,
                               const ShellData &sh_data_a, const ShellData &sh_data_b,
                               const ERI2Kernel *eri2_kernel)
{
    const int la = sh_data_a.l;
    const int lb = sh_data_b.l;
    const int lab = la + lb;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_hermite_a = numHermites(la);
    const int n_hermite_b = numHermites(lb);
    const int n_ecoeffs_a = n_sph_a * n_hermite_a;
    const int n_ecoeffs_b = n_sph_b * n_hermite_b;

    // Read-in data
    const int cdepth_a = sh_data_a.cdepths[ishell_a];
    const int cdepth_b = sh_data_b.cdepths[ishell_b];
    const int cofs_a = sh_data_a.coffsets[ishell_a];
    const int cofs_b = sh_data_b.coffsets[ishell_b];
    const int ofs_E_a = sh_data_a.offsets_ecoeffs[ishell_a];
    const int ofs_E_b = sh_data_b.offsets_ecoeffs[ishell_b];

    const double *exps_a = &sh_data_a.exps[cofs_a];
    const double *exps_b = &sh_data_b.exps[cofs_b];
    const double *coords_a = &sh_data_a.coords[3 * ishell_a];
    const double *coords_b = &sh_data_b.coords[3 * ishell_b];
    const double *ecoeffs_a = &eri2_kernel->ecoeffs_bra[ofs_E_a];
    const double *ecoeffs_b = &eri2_kernel->ecoeffs_ket[ofs_E_b];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(la);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lb);

    array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                            coords_a[1] - coords_b[1],
                            coords_a[2] - coords_b[2]};
    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

    vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);
    for (int ia = 0; ia < cdepth_a; ia++)
    {
        vector<double> R_x_E(n_hermite_a * n_sph_b, 0);
        for (int ib = 0; ib < cdepth_b; ib++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            double alpha = a * b / (a + b);
            double x = alpha * xyz_ab_dot;
            vector<double> fnx = calcBoysF(lab, x, eri2_kernel->boys_grid);

            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
            vector<double> rints = calcRIntsMatrix(lab, fac, alpha, &xyz_ab[0], &fnx[0],
                                                   hermite_idxs_bra, hermite_idxs_ket);

            int ofs_ecoeffs_b = ib * n_ecoeffs_b;
            shark_mm_ket(n_hermite_a, n_sph_b, n_hermite_b, &rints[0], &ecoeffs_b[ofs_ecoeffs_b], &R_x_E[0]);
        }
        int ofs_ecoeffs_a = ia * n_ecoeffs_a;
        shark_mm_bra(n_sph_a, n_sph_b, n_hermite_a, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[0], &eri2_batch[0]);
    }

    return eri2_batch;
}

array<lible::vec2d, 6> LI::eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                           const ShellData &sh_data_a,
                                           const ShellData &sh_data_b,
                                           const ERI2D1Kernel *eri2d1_kernel)
{
    const int la = sh_data_a.l;
    const int lb = sh_data_b.l;
    const int lab = la + lb;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_hermite_a = numHermites(la);
    const int n_hermite_b = numHermites(lb);
    const int n_ecoeffs_a = n_sph_a * n_hermite_a;
    const int n_ecoeffs_b = n_sph_b * n_hermite_b;
    const int n_rints = n_hermite_a * n_hermite_b;
    const int n_R_x_E = n_hermite_a * n_sph_b;

    // Read-in data
    const int cdepth_a = sh_data_a.cdepths[ishell_a];
    const int cdepth_b = sh_data_b.cdepths[ishell_b];
    const int cofs_a = sh_data_a.coffsets[ishell_a];
    const int cofs_b = sh_data_b.coffsets[ishell_b];
    const int ofs_E_a = sh_data_a.offsets_ecoeffs[ishell_a];
    const int ofs_E_b = sh_data_b.offsets_ecoeffs[ishell_b];

    const double *exps_a = &sh_data_a.exps[cofs_a];
    const double *exps_b = &sh_data_b.exps[cofs_b];
    const double *coords_a = &sh_data_a.coords[3 * ishell_a];
    const double *coords_b = &sh_data_b.coords[3 * ishell_b];
    const double *ecoeffs_a = &eri2d1_kernel->ecoeffs_bra[ofs_E_a];
    const double *ecoeffs_b = &eri2d1_kernel->ecoeffs_ket[ofs_E_b];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(la);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lb);

    array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                            coords_a[1] - coords_b[1],
                            coords_a[2] - coords_b[2]};
    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

    array<vec2d, 6> eri2_batch;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        eri2_batch[ideriv] = vec2d(Fill(0), n_sph_a, n_sph_b);

    for (int ia = 0; ia < cdepth_a; ia++)
    {        
        vector<double> R_x_E(3 * n_R_x_E, 0);
        for (int ib = 0; ib < cdepth_b; ib++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            double alpha = a * b / (a + b);
            double x = alpha * xyz_ab_dot;
            vector<double> fnx = calcBoysF(lab + 1, x, eri2d1_kernel->boys_grid);

            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
            vector<double> rints = calcRInts_ERI2D1(lab, alpha, fac, &fnx[0], &xyz_ab[0],
                                                    hermite_idxs_bra, hermite_idxs_ket);            

            int m = n_hermite_a, n = n_sph_b, k = n_hermite_b;
            int ofs_e0_b = ib * n_ecoeffs_b;

            shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[0 * n_R_x_E]);
            shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[1 * n_R_x_E]);
            shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[2 * n_R_x_E]);
        }

        // A
        int m = n_sph_a, n = n_sph_b, k = n_hermite_a;
        int ofs_e0_a = ia * n_ecoeffs_a;

        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[1 * n_R_x_E], &eri2_batch[1][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[2 * n_R_x_E], &eri2_batch[2][0]);
    }

    // B
    for (int ideriv = 3; ideriv < 6; ideriv++)
        eri2_batch[ideriv] = -1 * eri2_batch[ideriv - 3];

    return eri2_batch;
}

lible::arr2d<lible::vec2d, 6, 6> LI::eri2d2KernelFun(const int ishell_a, const int ishell_b,
                                                        const ShellData &sh_data_a,
                                                        const ShellData &sh_data_b,
                                                        const ERI2D2Kernel *eri2d2_kernel)
{
    const int la = sh_data_a.l;
    const int lb = sh_data_b.l;
    const int lab = la + lb;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_hermite_a = numHermites(la);
    const int n_hermite_b = numHermites(lb);
    const int n_ecoeffs_a = n_sph_a * n_hermite_a;
    const int n_ecoeffs_b = n_sph_b * n_hermite_b;
    const int n_rints = n_hermite_a * n_hermite_b;
    const int n_R_x_E = n_hermite_a * n_sph_b;

    // Read-in data
    const int cdepth_a = sh_data_a.cdepths[ishell_a];
    const int cdepth_b = sh_data_b.cdepths[ishell_b];
    const int cofs_a = sh_data_a.coffsets[ishell_a];
    const int cofs_b = sh_data_b.coffsets[ishell_b];
    const int ofs_E_a = sh_data_a.offsets_ecoeffs[ishell_a];
    const int ofs_E_b = sh_data_b.offsets_ecoeffs[ishell_b];

    const double *exps_a = &sh_data_a.exps[cofs_a];
    const double *exps_b = &sh_data_b.exps[cofs_b];
    const double *coords_a = &sh_data_a.coords[3 * ishell_a];
    const double *coords_b = &sh_data_b.coords[3 * ishell_b];
    const double *ecoeffs_a = &eri2d2_kernel->ecoeffs_bra[ofs_E_a];
    const double *ecoeffs_b = &eri2d2_kernel->ecoeffs_ket[ofs_E_b];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(la);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lb);

    array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                            coords_a[1] - coords_b[1],
                            coords_a[2] - coords_b[2]};
    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

    arr2d<vec2d, 6, 6> eri2_batch;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            eri2_batch[i][j] = vec2d(Fill(0), n_sph_a, n_sph_b);

    for (int ia = 0; ia < cdepth_a; ia++)
    {
        vector<double> R_x_E(6 * n_R_x_E, 0);
        for (int ib = 0; ib < cdepth_b; ib++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            double alpha = a * b / (a + b);
            double x = alpha * xyz_ab_dot;
            vector<double> fnx = calcBoysF(lab + 2, x, eri2d2_kernel->boys_grid);

            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
            vector<double> rints = calcRInts_ERI2D2(lab, alpha, fac, &fnx[0], &xyz_ab[0],
                                                    hermite_idxs_bra, hermite_idxs_ket);

            int m = n_hermite_a, n = n_sph_b, k = n_hermite_b;
            int ofs_e0_b = ib * n_ecoeffs_b;

            shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[0 * n_R_x_E]); // 200
            shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[1 * n_R_x_E]); // 110
            shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[2 * n_R_x_E]); // 101
            shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[3 * n_R_x_E]); // 020
            shark_mm_ket(m, n, k, &rints[4 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[4 * n_R_x_E]); // 011
            shark_mm_ket(m, n, k, &rints[5 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[5 * n_R_x_E]); // 002
        }

        int m = n_sph_a, n = n_sph_b, k = n_hermite_a;
        int ofs_e0_a = ia * n_ecoeffs_a;

        // AA upper triangle
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[1 * n_R_x_E], &eri2_batch[0][1][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[2 * n_R_x_E], &eri2_batch[0][2][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[3 * n_R_x_E], &eri2_batch[1][1][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[4 * n_R_x_E], &eri2_batch[1][2][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_e0_a], &R_x_E[5 * n_R_x_E], &eri2_batch[2][2][0]);
    }

    // Complete AA
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            eri2_batch[j][i] = eri2_batch[i][j];

    // AB und BA
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            eri2_batch[j][3 + i] = -1 * eri2_batch[i][j];
            eri2_batch[3 + i][j] = -1 * eri2_batch[i][j];
        }

    // BB
    for (int i = 3; i < 6; i++)
        for (int j = 3; j < 6; j++)
            eri2_batch[i][j] = -1 * eri2_batch[i - 3][j];

    return eri2_batch;
}

array<lible::vec3d, 9> LI::eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                           const ShellPairData &sp_data_ab,
                                           const ShellData &sh_data_c,
                                           const ERI3D1Kernel *eri3d1_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sh_data_c.l;
    const int lab = la + lb;
    const int labc = lab + lc;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_c = numHermites(lc);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_abc = n_sph_a * n_sph_b * n_sph_c;
    const int n_rints = n_hermite_ab * n_hermite_c;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_c = n_sph_c * n_hermite_c;
    const int n_R_x_E = n_hermite_ab * n_sph_c;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];
    const int ofs_E0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_E0_c = sh_data_c.offsets_ecoeffs[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *ecoeffs0_ab = &eri3d1_kernel->ecoeffs0_bra[ofs_E0_ab];
    const double *ecoeffs1_ab = &eri3d1_kernel->ecoeffs1_bra[ofs_E1_ab];
    const double *ecoeffs0_c = &eri3d1_kernel->ecoeffs0_ket[ofs_E0_c];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lc);

    array<vec3d, 9> eri3_batch;
    for (int ideriv = 0; ideriv < 9; ideriv++)
        eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(4 * n_R_x_E, 0);
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double c = exps_c[ic];

                double alpha = p * c / (p + c);

                array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                       (a * coords_a[1] + b * coords_b[1]) / p,
                                       (a * coords_a[2] + b * coords_b[2]) / p};

                array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                        xyz_p[1] - coords_c[1],
                                        xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                vector<double> fnx = calcBoysF(labc + 1, x, eri3d1_kernel->boys_grid);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                vector<double> rints = calcRInts_ERI3D1(labc, alpha, fac, &fnx[0], &xyz_pc[0],
                                                        hermite_idxs_bra, hermite_idxs_ket);

                int m = n_hermite_ab, n = n_sph_c, k = n_hermite_c;
                int ofs_e0_c = ic * n_ecoeffs_c;

                shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[0 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[1 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[2 * n_R_x_E]);                

                shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[3 * n_R_x_E]);
            }

            int m = n_sph_ab, n = n_sph_c, k = n_hermite_ab;
            int ofs_e0_ab = iab * n_ecoeffs_ab;
            int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;

            // P & R
            vector<double> P(3 * n_sph_abc, 0);
            vector<double> R(3 * n_sph_abc, 0);

            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abc]);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abc]);

            // A
            cblas_daxpy(n_sph_abc, (a / p), &P[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, (a / p), &P[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, (a / p), &P[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

            cblas_daxpy(n_sph_abc, 1.0, &R[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &R[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &R[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

            // C
            cblas_daxpy(n_sph_abc, -1.0, &P[0 * n_sph_abc], 1, &eri3_batch[6][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &P[1 * n_sph_abc], 1, &eri3_batch[7][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &P[2 * n_sph_abc], 1, &eri3_batch[8][0], 1);
        }

    // B
    for (int ideriv = 3; ideriv < 6; ideriv++)
        eri3_batch[ideriv] = -1 * (eri3_batch[ideriv - 3] + eri3_batch[ideriv + 3]);

    return eri3_batch;
}

lible::arr2d<lible::vec4d, 12, 12> LI::eri4d2KernelFun(const int ipair_ab, const int ipair_cd,
                                                       const ShellPairData &sp_data_ab,
                                                       const ShellPairData &sp_data_cd,
                                                       const ERI4D2Kernel *eri4d2_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sp_data_cd.la;
    const int ld = sp_data_cd.lb;
    const int lab = la + lb;
    const int lcd = lc + ld;
    const int labcd = lab + lcd;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_sph_d = numSphericals(ld);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_cd = n_sph_c * n_sph_d;
    const int n_sph_abcd = n_sph_ab * n_sph_cd;
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_cd = numHermites(lcd);
    const int n_rints = n_hermite_ab * n_hermite_cd;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
    const int n_R_x_E = n_hermite_ab * n_sph_cd;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    const int ofs_e0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_e1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_e2_ab = sp_data_ab.offsets_ecoeffs_deriv2[ipair_ab];
    const int ofs_e0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];
    const int ofs_e1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];
    const int ofs_e2_cd = sp_data_cd.offsets_ecoeffs_deriv2[ipair_cd];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *ecoeffs0_ab = &eri4d2_kernel->ecoeffs0_bra[ofs_e0_ab];
    const double *ecoeffs1_ab = &eri4d2_kernel->ecoeffs1_bra[ofs_e1_ab];
    const double *ecoeffs2_ab = &eri4d2_kernel->ecoeffs1_bra[ofs_e2_ab];
    const double *ecoeffs0_cd = &eri4d2_kernel->ecoeffs0_ket[ofs_e0_cd];
    const double *ecoeffs1_cd = &eri4d2_kernel->ecoeffs1_ket[ofs_e1_cd];
    const double *ecoeffs2_cd = &eri4d2_kernel->ecoeffs1_ket[ofs_e2_cd];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lcd);

    arr2d<vec4d, 12, 12> eri4_batch;
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++)
            eri4_batch[i][j] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(10 * n_R_x_E, 0);
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double q = c + d;
                    double alpha = p * q / (p + q);

                    std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                (a * xyz_a[2] + b * xyz_b[2]) / p};

                    std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                (c * xyz_c[2] + d * xyz_d[2]) / q};

                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                 xyz_p[1] - xyz_q[1],
                                                 xyz_p[2] - xyz_q[2]};

                    double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    vector<double> fnx = calcBoysF(labcd + 2, x, eri4d2_kernel->boys_grid);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    vector<double> rints = calcRInts_ERI3D2(labcd, alpha, fac, &fnx[0], &xyz_pq[0],
                                                            hermite_idxs_bra, hermite_idxs_ket);
                }
        }

    return eri4_batch;    
}

lible::arr2d<lible::vec3d, 9, 9> LI::eri3d2KernelFun(const int ipair_ab, const int ishell_c,
                                                     const ShellPairData &sp_data_ab,
                                                     const ShellData &sh_data_c,
                                                     const ERI3D2Kernel *eri3d2_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sh_data_c.l;
    const int lab = la + lb;
    const int labc = lab + lc;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_c = numHermites(lc);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_abc = n_sph_a * n_sph_b * n_sph_c;
    const int n_rints = n_hermite_ab * n_hermite_c;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_c = n_sph_c * n_hermite_c;
    const int n_R_x_E = n_hermite_ab * n_sph_c;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];
    const int ofs_E0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_E2_ab = sp_data_ab.offsets_ecoeffs_deriv2[ipair_ab];
    const int ofs_E0_c = sh_data_c.offsets_ecoeffs[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *ecoeffs0_ab = &eri3d2_kernel->ecoeffs0_bra[ofs_E0_ab];
    const double *ecoeffs1_ab = &eri3d2_kernel->ecoeffs1_bra[ofs_E1_ab];
    const double *ecoeffs2_ab = &eri3d2_kernel->ecoeffs2_bra[ofs_E2_ab];
    const double *ecoeffs0_c = &eri3d2_kernel->ecoeffs0_ket[ofs_E0_c];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lc);

    arr2d<vec3d, 9, 9> eri3_batch;
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            eri3_batch[i][j] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            vector<double> R_x_E(10 * n_R_x_E, 0);
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double c = exps_c[ic];

                double alpha = p * c / (p + c);

                array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                       (a * coords_a[1] + b * coords_b[1]) / p,
                                       (a * coords_a[2] + b * coords_b[2]) / p};

                array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                        xyz_p[1] - coords_c[1],
                                        xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                vector<double> fnx = calcBoysF(labc + 2, x, eri3d2_kernel->boys_grid);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                vector<double> rints = calcRInts_ERI3D2(labc, alpha, fac, &fnx[0], &xyz_pc[0],
                                                        hermite_idxs_bra, hermite_idxs_ket);

                int m = n_hermite_ab, n = n_sph_c, k = n_hermite_c;
                int ofs_e0_c = ic * n_ecoeffs_c;

                vector<double> R_x_E(10 * n_R_x_E, 0);
                for (int i = 0; i < 10; i++)
                    shark_mm_ket(m, n, k, &rints[i * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[i * n_R_x_E]);
            }
        }

    return eri3_batch;
}

array<lible::vec4d, 12> LI::eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                            const ShellPairData &sp_data_ab,
                                            const ShellPairData &sp_data_cd,
                                            const ERI4D1Kernel *eri4d1_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sp_data_cd.la;
    const int ld = sp_data_cd.lb;
    const int lab = la + lb;
    const int lcd = lc + ld;
    const int labcd = lab + lcd;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_sph_d = numSphericals(ld);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_cd = n_sph_c * n_sph_d;
    const int n_sph_abcd = n_sph_ab * n_sph_cd;
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_cd = numHermites(lcd);
    const int n_rints = n_hermite_ab * n_hermite_cd;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
    const int n_R_x_E = n_hermite_ab * n_sph_cd;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    const int ofs_E0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_E0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];
    const int ofs_E1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *ecoeffs0_ab = &eri4d1_kernel->ecoeffs0_bra[ofs_E0_ab];
    const double *ecoeffs1_ab = &eri4d1_kernel->ecoeffs1_bra[ofs_E1_ab];
    const double *ecoeffs0_cd = &eri4d1_kernel->ecoeffs0_ket[ofs_E0_cd];
    const double *ecoeffs1_cd = &eri4d1_kernel->ecoeffs1_ket[ofs_E1_cd];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lcd);

    array<vec4d, 12> eri4_batch;
    for (int ideriv = 0; ideriv < 12; ideriv++)
        eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(7 * n_R_x_E, 0);
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double q = c + d;
                    double alpha = p * q / (p + q);

                    std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                (a * xyz_a[2] + b * xyz_b[2]) / p};

                    std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                (c * xyz_c[2] + d * xyz_d[2]) / q};

                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                 xyz_p[1] - xyz_q[1],
                                                 xyz_p[2] - xyz_q[2]};

                    double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    vector<double> fnx = calcBoysF(labcd + 1, x, eri4d1_kernel->boys_grid);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    vector<double> rints = calcRInts_ERI3D1(labcd, alpha, fac, &fnx[0], &xyz_pq[0],
                                                            hermite_idxs_bra, hermite_idxs_ket);

                    vector<double> I(3 * n_R_x_E, 0);

                    int ofs_e0_cd = icd * n_ecoeffs_cd;
                    int ofs_e1_cd = 3 * icd * n_ecoeffs_cd;
                    int m = n_hermite_ab, n = n_sph_cd, k = n_hermite_cd;

                    shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &I[0 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &I[1 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &I[2 * n_R_x_E]);

                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[3 * n_R_x_E]);

                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_e1_cd + 0 * n_ecoeffs_cd], &R_x_E[4 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_e1_cd + 1 * n_ecoeffs_cd], &R_x_E[5 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_e1_cd + 2 * n_ecoeffs_cd], &R_x_E[6 * n_R_x_E]);                    

                    cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
                    cblas_daxpy(3 * n_R_x_E, -(c / q), &I[0], 1, &R_x_E[4 * n_R_x_E], 1);
                }

            int ofs_e0_ab = iab * n_ecoeffs_ab;
            int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;

            int m = n_sph_ab, n = n_sph_cd, k = n_hermite_ab;

            // P & R
            vector<double> P(3 * n_sph_abcd, 0);
            vector<double> R(3 * n_sph_abcd, 0);

            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abcd]);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abcd]);

            // A
            cblas_daxpy(n_sph_abcd, (a / p), &P[0 * n_sph_abcd], 1, &eri4_batch[0][0], 1);
            cblas_daxpy(n_sph_abcd, (a / p), &P[1 * n_sph_abcd], 1, &eri4_batch[1][0], 1);
            cblas_daxpy(n_sph_abcd, (a / p), &P[2 * n_sph_abcd], 1, &eri4_batch[2][0], 1);

            cblas_daxpy(n_sph_abcd, 1.0, &R[0 * n_sph_abcd], 1, &eri4_batch[0][0], 1);
            cblas_daxpy(n_sph_abcd, 1.0, &R[1 * n_sph_abcd], 1, &eri4_batch[1][0], 1);
            cblas_daxpy(n_sph_abcd, 1.0, &R[2 * n_sph_abcd], 1, &eri4_batch[2][0], 1);

            // B
            cblas_daxpy(n_sph_abcd, (b / p), &P[0 * n_sph_abcd], 1, &eri4_batch[3][0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &P[1 * n_sph_abcd], 1, &eri4_batch[4][0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &P[2 * n_sph_abcd], 1, &eri4_batch[5][0], 1);

            cblas_daxpy(n_sph_abcd, -1.0, &R[0 * n_sph_abcd], 1, &eri4_batch[3][0], 1);
            cblas_daxpy(n_sph_abcd, -1.0, &R[1 * n_sph_abcd], 1, &eri4_batch[4][0], 1);
            cblas_daxpy(n_sph_abcd, -1.0, &R[2 * n_sph_abcd], 1, &eri4_batch[5][0], 1);

            // C         
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[4 * n_R_x_E], &eri4_batch[6][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[5 * n_R_x_E], &eri4_batch[7][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_e0_ab], &R_x_E[6 * n_R_x_E], &eri4_batch[8][0]);
        }

    // D
    for (int ideriv = 9; ideriv < 12; ideriv++)
        eri4_batch[ideriv] = -1 * (eri4_batch[ideriv - 9] + eri4_batch[ideriv - 6] + eri4_batch[ideriv - 3]);

    return eri4_batch;
}

array<lible::vec4d, 3> LI::eri4socKernelFun(const int ipair_ab, const int ipair_cd,
                                            const ShellPairData &sp_data_ab,
                                            const ShellPairData &sp_data_cd,
                                            const ERI4SOCKernel *eri4soc_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sp_data_cd.la;
    const int ld = sp_data_cd.lb;
    const int lab = la + lb;
    const int lcd = lc + ld;
    const int labcd = lab + lcd;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_sph_d = numSphericals(ld);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_cd = n_sph_c * n_sph_d;
    const int n_sph_abcd = n_sph_ab * n_sph_cd;
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_cd = numHermites(lcd);
    const int n_rints = n_hermite_ab * n_hermite_cd;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
    const int n_R_x_E = n_hermite_ab * n_sph_cd;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_E0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *ecoeffs1_ab = &eri4soc_kernel->ecoeffs1_bra[ofs_E1_ab];
    const double *ecoeffs0_cd = &eri4soc_kernel->ecoeffs0_ket[ofs_E0_cd];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lcd);

    array<vec4d, 3> eri4_batch;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            vector<double> R_x_E(3 * n_R_x_E, 0);
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double q = c + d;
                    double alpha = p * q / (p + q);

                    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                (a * xyz_a[2] + b * xyz_b[2]) / p};

                    array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                (c * xyz_c[2] + d * xyz_d[2]) / q};

                    array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                            xyz_p[1] - xyz_q[1],
                                            xyz_p[2] - xyz_q[2]};

                    double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    vector<double> fnx = calcBoysF(labcd + 1, x, eri4soc_kernel->boys_grid);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    vector<double> rints = calcRInts_ERISOC(labcd, fac, alpha, &xyz_pq[0], &fnx[0],
                                                            hermite_idxs_bra, hermite_idxs_ket);

                    int m = n_hermite_ab, n = n_sph_cd, k = n_hermite_cd;
                    int ofs_e0_cd = icd * n_ecoeffs_cd;

                    shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[0 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[1 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[2 * n_R_x_E]);
                }

            int m = n_sph_ab, n = n_sph_cd, k = n_hermite_ab;
            int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;
            int ofs_e1_0 = ofs_e1_ab + 0 * n_ecoeffs_ab;
            int ofs_e1_1 = ofs_e1_ab + 1 * n_ecoeffs_ab;
            int ofs_e1_2 = ofs_e1_ab + 2 * n_ecoeffs_ab;

            // P & R
            vector<double> PR_yz(n_sph_abcd, 0);
            vector<double> PR_zy(n_sph_abcd, 0);
            vector<double> PR_zx(n_sph_abcd, 0);
            vector<double> PR_xz(n_sph_abcd, 0);
            vector<double> PR_xy(n_sph_abcd, 0);
            vector<double> PR_yx(n_sph_abcd, 0);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_2], &R_x_E[1 * n_R_x_E], &PR_yz[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_1], &R_x_E[2 * n_R_x_E], &PR_zy[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_0], &R_x_E[2 * n_R_x_E], &PR_zx[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_2], &R_x_E[0 * n_R_x_E], &PR_xz[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_1], &R_x_E[0 * n_R_x_E], &PR_xy[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_0], &R_x_E[1 * n_R_x_E], &PR_yx[0]);

            // A & B
            vector<double> AB_yz(n_sph_abcd, 0);
            vector<double> AB_zy(n_sph_abcd, 0);
            vector<double> AB_zx(n_sph_abcd, 0);
            vector<double> AB_xz(n_sph_abcd, 0);
            vector<double> AB_xy(n_sph_abcd, 0);
            vector<double> AB_yx(n_sph_abcd, 0);

            cblas_daxpy(n_sph_abcd, -(a / p), &PR_yz[0], 1, &AB_yz[0], 1);
            cblas_daxpy(n_sph_abcd, -(a / p), &PR_zy[0], 1, &AB_zy[0], 1);
            cblas_daxpy(n_sph_abcd, -(a / p), &PR_zx[0], 1, &AB_zx[0], 1);
            cblas_daxpy(n_sph_abcd, -(a / p), &PR_xz[0], 1, &AB_xz[0], 1);
            cblas_daxpy(n_sph_abcd, -(a / p), &PR_xy[0], 1, &AB_xy[0], 1);
            cblas_daxpy(n_sph_abcd, -(a / p), &PR_yx[0], 1, &AB_yx[0], 1);

            cblas_daxpy(n_sph_abcd, (b / p), &PR_yz[0], 1, &AB_yz[0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &PR_zy[0], 1, &AB_zy[0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &PR_zx[0], 1, &AB_zx[0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &PR_xz[0], 1, &AB_xz[0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &PR_xy[0], 1, &AB_xy[0], 1);
            cblas_daxpy(n_sph_abcd, (b / p), &PR_yx[0], 1, &AB_yx[0], 1);

            cblas_daxpy(n_sph_abcd, 1.0, &AB_yz[0], 1, &eri4_batch[0][0], 1);
            cblas_daxpy(n_sph_abcd, 1.0, &AB_zx[0], 1, &eri4_batch[1][0], 1);
            cblas_daxpy(n_sph_abcd, 1.0, &AB_xy[0], 1, &eri4_batch[2][0], 1);

            cblas_daxpy(n_sph_abcd, -1.0, &AB_zy[0], 1, &eri4_batch[0][0], 1);
            cblas_daxpy(n_sph_abcd, -1.0, &AB_xz[0], 1, &eri4_batch[1][0], 1);
            cblas_daxpy(n_sph_abcd, -1.0, &AB_yx[0], 1, &eri4_batch[2][0], 1);
        }

    return eri4_batch;
}

array<lible::vec3d, 3> LI::eri3socKernelFun(const int ipair_ab, const int ishell_c,
                                            const ShellPairData &sp_data_ab,
                                            const ShellData &sh_data_c,
                                            const ERI3SOCKernel *eri3soc_kernel)
{
    const int la = sp_data_ab.la;
    const int lb = sp_data_ab.lb;
    const int lc = sh_data_c.l;
    const int lab = la + lb;
    const int labc = lab + lc;

    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_c = numSphericals(lc);
    const int n_hermite_ab = numHermites(lab);
    const int n_hermite_c = numHermites(lc);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_sph_abc = n_sph_a * n_sph_b * n_sph_c;
    const int n_rints = n_hermite_ab * n_hermite_c;
    const int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_c = n_sph_c * n_hermite_c;
    const int n_R_x_E = n_hermite_ab * n_sph_c;

    // Read-in data
    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    const int cdepth_c = sh_data_c.cdepths[ishell_c];
    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    const int cofs_c = sh_data_c.coffsets[ishell_c];
    const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    const int ofs_E0_c = sh_data_c.offsets_ecoeffs[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *ecoeffs1_ab = &eri3soc_kernel->ecoeffs1_bra[ofs_E1_ab];
    const double *ecoeffs0_c = &eri3soc_kernel->ecoeffs0_ket[ofs_E0_c];

    // SHARK integrals
    vector<array<int, 3>> hermite_idxs_bra = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> hermite_idxs_ket = getHermiteGaussianIdxs(lc);

    array<vec3d, 3> eri3_batch;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(3 * n_R_x_E, 0);
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double c = exps_c[ic];

                double alpha = p * c / (p + c);

                array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                       (a * coords_a[1] + b * coords_b[1]) / p,
                                       (a * coords_a[2] + b * coords_b[2]) / p};

                array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                        xyz_p[1] - coords_c[1],
                                        xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                vector<double> fnx = calcBoysF(labc + 1, x, eri3soc_kernel->boys_grid);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                vector<double> rints = calcRInts_ERISOC(labc, fac, alpha, &xyz_pc[0], &fnx[0], 
                                                        hermite_idxs_bra, hermite_idxs_ket);

                int m = n_hermite_ab, n = n_sph_c, k = n_hermite_c;
                int ofs_e0_c = ic * n_ecoeffs_c;

                shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[0 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[1 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[2 * n_R_x_E]);
            }

            int m = n_sph_ab, n = n_sph_c, k = n_hermite_ab;
            int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;
            int ofs_e1_0 = ofs_e1_ab + 0 * n_ecoeffs_ab;
            int ofs_e1_1 = ofs_e1_ab + 1 * n_ecoeffs_ab;
            int ofs_e1_2 = ofs_e1_ab + 2 * n_ecoeffs_ab;            

            // P & R
            vector<double> PR_yz(n_sph_abc, 0);
            vector<double> PR_zy(n_sph_abc, 0);
            vector<double> PR_zx(n_sph_abc, 0);
            vector<double> PR_xz(n_sph_abc, 0);
            vector<double> PR_xy(n_sph_abc, 0);
            vector<double> PR_yx(n_sph_abc, 0);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_2], &R_x_E[1 * n_R_x_E], &PR_yz[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_1], &R_x_E[2 * n_R_x_E], &PR_zy[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_0], &R_x_E[2 * n_R_x_E], &PR_zx[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_2], &R_x_E[0 * n_R_x_E], &PR_xz[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_1], &R_x_E[0 * n_R_x_E], &PR_xy[0]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_e1_0], &R_x_E[1 * n_R_x_E], &PR_yz[0]);

            // A & B
            vector<double> AB_yz(n_sph_abc, 0);
            vector<double> AB_zy(n_sph_abc, 0);
            vector<double> AB_zx(n_sph_abc, 0);
            vector<double> AB_xz(n_sph_abc, 0);
            vector<double> AB_xy(n_sph_abc, 0);
            vector<double> AB_yx(n_sph_abc, 0);

            cblas_daxpy(n_sph_abc, -(a / p), &PR_yz[0], 1, &AB_yz[0], 1);
            cblas_daxpy(n_sph_abc, -(a / p), &PR_zy[0], 1, &AB_zy[0], 1);
            cblas_daxpy(n_sph_abc, -(a / p), &PR_zx[0], 1, &AB_zx[0], 1);
            cblas_daxpy(n_sph_abc, -(a / p), &PR_xz[0], 1, &AB_xz[0], 1);
            cblas_daxpy(n_sph_abc, -(a / p), &PR_xy[0], 1, &AB_xy[0], 1);
            cblas_daxpy(n_sph_abc, -(a / p), &PR_yx[0], 1, &AB_yx[0], 1);

            cblas_daxpy(n_sph_abc, (b / p), &PR_yz[0], 1, &AB_yz[0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &PR_zy[0], 1, &AB_zy[0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &PR_zx[0], 1, &AB_zx[0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &PR_xz[0], 1, &AB_xz[0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &PR_xy[0], 1, &AB_xy[0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &PR_yx[0], 1, &AB_yx[0], 1);

            cblas_daxpy(n_sph_abc, 1.0, &AB_yz[0], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &AB_zx[0], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &AB_xy[0], 1, &eri3_batch[2][0], 1);

            cblas_daxpy(n_sph_abc, -1.0, &AB_zy[0], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &AB_xz[0], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &AB_yx[0], 1, &eri3_batch[2][0], 1);
        }

    return eri3_batch;
}