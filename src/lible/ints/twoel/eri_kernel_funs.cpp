#include <lible/ints/rints.hpp>
#include <lible/ints/twoel/eri_kernel_funs.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

lible::vec4d LI::eri4KernelFun(const int ipair_ab, const int ipair_cd,
                               const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                               const ERI4Kernel *eri4_kernel)
{
    // Angmom-related data
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

                    std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                (a * xyz_a[2] + b * xyz_b[2]) / p};

                    std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                (c * xyz_c[2] + d * xyz_d[2]) / q};

                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                 xyz_p[1] - xyz_q[1],
                                                 xyz_p[2] - xyz_q[2]};

                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                    double x = alpha * (dx * dx + dy * dy + dz * dz);
                    vector<double> fnx = calcBoysF(labcd, x, eri4_kernel->boys_grid);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                    vector<double> rints = calcRIntsMatrix(labcd, fac, alpha, &xyz_pq[0], &fnx[0],
                                                           hermite_idxs_bra, hermite_idxs_ket);

                    int ofs_ecoeffs_cd = icd * n_ecoeffs_cd;
                    shark_mm_ket(n_hermite_ab, n_sph_cd, n_hermite_cd, &rints[0], &ecoeffs_cd[ofs_ecoeffs_cd], &R_x_E[0]);
                }
            int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
            shark_mm_bra(n_sph_ab, n_sph_cd, n_hermite_ab, &ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri4_batch[0]);
        }

    // Norms
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

                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                            (a * coords_a[2] + b * coords_b[2]) / p};

                std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                             xyz_p[1] - coords_c[1],
                                             xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                vector<double> fnx = calcBoysF(labc, x, eri3_kernel->boys_grid);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                vector<double> rints = calcRIntsMatrix(labc, fac, alpha, &xyz_pc[0], &fnx[0],
                                                       hermite_idxs_bra, hermite_idxs_ket);

                int ofs_ecoeffs_c = ic * n_ecoeffs_c;
                shark_mm_ket(n_hermite_ab, n_sph_c, n_hermite_c, &rints[0], &ecoeffs_c[ofs_ecoeffs_c], &R_x_E[0]);
            }
            int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
            shark_mm_bra(n_sph_ab, n_sph_c, n_hermite_ab, &ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri3_batch[0]);
        }

    // Norms
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

    // Norms
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

std::array<lible::vec2d, 6> LI::eri2d1KernelFun(const int ishell_a, const int ishell_b,
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
        vector<double> R_x_E(6 * n_R_x_E, 0);
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

            vector<double> I(3 * n_R_x_E, 0);

            int m = n_hermite_a, n = n_sph_b, k = n_hermite_b;
            int ofs_ecoeffs_b = ib * n_ecoeffs_b;
            
            shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[0 * n_R_x_E]);
            shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[1 * n_R_x_E]);
            shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[2 * n_R_x_E]);

            cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
            cblas_daxpy(3 * n_R_x_E, -1.0, &I[0], 1, &R_x_E[3 * n_R_x_E], 1);
        }

        int m = n_sph_a, n = n_sph_b, k = n_hermite_a;
        int ofs_ecoeffs_a = ia * n_ecoeffs_a;

        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[1 * n_R_x_E], &eri2_batch[1][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[2 * n_R_x_E], &eri2_batch[2][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[3 * n_R_x_E], &eri2_batch[3][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[4 * n_R_x_E], &eri2_batch[4][0]);
        shark_mm_bra(m, n, k, &ecoeffs_a[ofs_ecoeffs_a], &R_x_E[5 * n_R_x_E], &eri2_batch[5][0]);
    }

    // Norms
    int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];
    int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (int mu = 0; mu < n_sph_a; mu++)
            for (int nu = 0; nu < n_sph_b; nu++)
            {
                double norm_a = sh_data_a.norms[ofs_norm_a + mu];
                double norm_b = sh_data_b.norms[ofs_norm_b + nu];

                eri2_batch[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return eri2_batch;
}

std::array<lible::vec3d, 9> LI::eri3d1KernelFun(const int ipair_ab, const int ishell_c,
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

    std::array<vec3d, 9> eri3_batch;
    for (int ideriv = 0; ideriv < 9; ideriv++)
        eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

    vector<double> rints(4 * n_rints, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(7 * n_R_x_E, 0);
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

                vector<double> I(3 * n_R_x_E, 0);

                int m = n_hermite_ab, n = n_sph_c, k = n_hermite_c;
                int ofs_ecoeffs_c = ic * n_ecoeffs_c;

                shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[0 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[1 * n_R_x_E]);
                shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[2 * n_R_x_E]);

                shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[3 * n_R_x_E]);

                cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
                cblas_daxpy(3 * n_R_x_E, -1.0, &I[0], 1, &R_x_E[4 * n_R_x_E], 1);                
            }

            int m = n_sph_ab, n = n_sph_c, k = n_hermite_ab;
            int ofs_ecoeffs0_ab = iab * n_ecoeffs_ab;
            int ofs_ecoeffs1_ab = 3 * iab * n_ecoeffs_ab;

            // P & R
            std::vector<double> P(3 * n_sph_abc, 0);
            std::vector<double> R(3 * n_sph_abc, 0);

            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abc]);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abc]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abc]);

            // A
            cblas_daxpy(n_sph_abc, (a / p), &P[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, (a / p), &P[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, (a / p), &P[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

            cblas_daxpy(n_sph_abc, 1.0, &R[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &R[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
            cblas_daxpy(n_sph_abc, 1.0, &R[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

            // B
            cblas_daxpy(n_sph_abc, (b / p), &P[0 * n_sph_abc], 1, &eri3_batch[3][0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &P[1 * n_sph_abc], 1, &eri3_batch[4][0], 1);
            cblas_daxpy(n_sph_abc, (b / p), &P[2 * n_sph_abc], 1, &eri3_batch[5][0], 1);

            cblas_daxpy(n_sph_abc, -1.0, &R[0 * n_sph_abc], 1, &eri3_batch[3][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &R[1 * n_sph_abc], 1, &eri3_batch[4][0], 1);
            cblas_daxpy(n_sph_abc, -1.0, &R[2 * n_sph_abc], 1, &eri3_batch[5][0], 1);

            // C
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[4 * n_R_x_E], &eri3_batch[6][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[5 * n_R_x_E], &eri3_batch[7][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[6 * n_R_x_E], &eri3_batch[8][0]);
        }

    // Norms
    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
    for (int ideriv = 0; ideriv < 9; ideriv++)
        for (int a = 0; a < n_sph_a; a++)
            for (int b = 0; b < n_sph_b; b++)
                for (int c = 0; c < n_sph_c; c++)
                {
                    double norm_a = sp_data_ab.norms[ofs_norm_a + a];
                    double norm_b = sp_data_ab.norms[ofs_norm_b + b];
                    double norm_c = sh_data_c.norms[ofs_norm_c + c];

                    eri3_batch[ideriv](a, b, c) *= norm_a * norm_b * norm_c;
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

    std::array<vec4d, 12> eri4_batch;
    for (int ideriv = 0; ideriv < 12; ideriv++)
        eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);    
  
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            std::vector<double> R_x_E(13 * n_R_x_E, 0);
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

                    std::vector<double> I1(3 * n_R_x_E, 0);
                    std::vector<double> I2(3 * n_R_x_E, 0);

                    int ofs_ecoeffs0_cd = icd * n_ecoeffs_cd;
                    int ofs_ecoeffs1_cd = 3 * icd * n_ecoeffs_cd;
                    int m = n_hermite_ab, n = n_sph_cd, k = n_hermite_cd;

                    shark_mm_ket(m, n, k, &rints[0 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I1[0 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[1 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I1[1 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[2 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I1[2 * n_R_x_E]);

                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &R_x_E[3 * n_R_x_E]);

                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 0 * n_ecoeffs_cd], &I2[0 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 1 * n_ecoeffs_cd], &I2[1 * n_R_x_E]);
                    shark_mm_ket(m, n, k, &rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 2 * n_ecoeffs_cd], &I2[2 * n_R_x_E]);

                    cblas_daxpy(3 * n_R_x_E, 1.0, &I1[0], 1, &R_x_E[0 * n_R_x_E], 1);
                    cblas_daxpy(3 * n_R_x_E, -1.0, &I1[0], 1, &R_x_E[4 * n_R_x_E], 1);

                    cblas_daxpy(3 * n_R_x_E, -(c / q), &I1[0], 1, &R_x_E[7 * n_R_x_E], 1);
                    cblas_daxpy(3 * n_R_x_E, -(d / q), &I1[0], 1, &R_x_E[10 * n_R_x_E], 1);

                    cblas_daxpy(3 * n_R_x_E, 1.0, &I2[0], 1, &R_x_E[7 * n_R_x_E], 1);
                    cblas_daxpy(3 * n_R_x_E, -1.0, &I2[0], 1, &R_x_E[10 * n_R_x_E], 1);
                }

            int ofs_ecoeffs0_ab = iab * n_ecoeffs_ab;
            int ofs_ecoeffs1_ab = 3 * iab * n_ecoeffs_ab;

            int m = n_sph_ab, n = n_sph_cd, k = n_hermite_ab;

            // P & R
            std::vector<double> P(3 * n_sph_abcd, 0);
            std::vector<double> R(3 * n_sph_abcd, 0);

            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abcd]);

            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abcd]);
            shark_mm_bra(m, n, k, &ecoeffs1_ab[ofs_ecoeffs1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abcd]);

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

            // C & D
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[7 * n_R_x_E], &eri4_batch[6][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[8 * n_R_x_E], &eri4_batch[7][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[9 * n_R_x_E], &eri4_batch[8][0]);

            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[10 * n_R_x_E], &eri4_batch[9][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[11 * n_R_x_E], &eri4_batch[10][0]);
            shark_mm_bra(m, n, k, &ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[12 * n_R_x_E], &eri4_batch[11][0]);
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd + 0];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    for (int ideriv = 0; ideriv < 12; ideriv++)
        for (int a = 0; a < n_sph_a; a++)
            for (int b = 0; b < n_sph_b; b++)
                for (int c = 0; c < n_sph_c; c++)
                    for (int d = 0; d < n_sph_d; d++)
                    {
                        double norm_a = sp_data_ab.norms[ofs_norm_a + a];
                        double norm_b = sp_data_ab.norms[ofs_norm_b + b];
                        double norm_c = sp_data_cd.norms[ofs_norm_c + c];
                        double norm_d = sp_data_cd.norms[ofs_norm_d + d];                        

                        eri4_batch[ideriv](a, b, c, d) *= norm_a * norm_b * norm_c * norm_d;                                                        
                    }

    return eri4_batch;
}