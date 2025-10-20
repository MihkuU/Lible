#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/utils.hpp>

#include <array>
#include <tuple>

namespace lints = lible::ints;

using std::array, std::pair, std::tuple, std::vector;

lible::vec2d lints::ecoeffsRecurrence1(const double one_o_2a, const int l)
{
    vec2d ecoeffs(Fill(0), l + 1, l + 1);

    ecoeffs(0, 0) = 1;

    for (int i = 1; i <= l; i++)
    {
        if (i % 2 == 0)
            ecoeffs(i, 0) = ecoeffs(i - 1, 1);

        for (int t = 1; t < i; t++)
            if ((t + i) % 2 == 0)
                ecoeffs(i, t) = one_o_2a * ecoeffs(i - 1, t - 1) +
                                (t + 1) * ecoeffs(i - 1, t + 1);

        ecoeffs(i, i) = one_o_2a * ecoeffs(i - 1, i - 1);
    }

    return ecoeffs;
}

lible::vec3d lints::ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                       const double PA, const double PB, const double Kab)
{
    const double p = a + b;
    const double one_o_2p = 1.0 / (2 * p);

    vec3d ecoeffs(Fill(0), la + 1, lb + 1, la + lb + 1);
    ecoeffs(0, 0, 0) = Kab;
    for (int i = 1; i <= la; i++)
    {
        ecoeffs(i, 0, 0) = PA * ecoeffs(i - 1, 0, 0) + ecoeffs(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            ecoeffs(i, 0, t) = one_o_2p * ecoeffs(i - 1, 0, t - 1) +
                               PA * ecoeffs(i - 1, 0, t) +
                               (t + 1) * ecoeffs(i - 1, 0, t + 1);

        ecoeffs(i, 0, i) = one_o_2p * ecoeffs(i - 1, 0, i - 1) + PA * ecoeffs(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            ecoeffs(i, j, 0) = PB * ecoeffs(i, j - 1, 0) + ecoeffs(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                ecoeffs(i, j, t) = one_o_2p * ecoeffs(i, j - 1, t - 1) +
                                   PB * ecoeffs(i, j - 1, t) +
                                   (t + 1) * ecoeffs(i, j - 1, t + 1);

            ecoeffs(i, j, i + j) = one_o_2p * ecoeffs(i, j - 1, i + j - 1) +
                                   PB * ecoeffs(i, j - 1, i + j);
        }

    return ecoeffs;
}

lible::vec3d lints::ecoeffsRecurrence2_n1(const double a, const double b, const int la, const int lb,
                                          const double A, const double B, const vec3d &ecoeffs0)
{
    vec3d ecoeffs1(Fill(0), la + 1, lb + 1, la + lb + 1);

    const double p = a + b;
    const double one_o_2p = 1.0 / (2 * p);
    const double R = A - B;

    ecoeffs1(0, 0, 0) = -2 * (a * b) * R * ecoeffs0(0, 0, 0) / p;

    for (int i = 1; i <= la; i++)
    {
        ecoeffs1(i, 0, 0) = -(b / p) * (R * ecoeffs1(i - 1, 0, 0) + ecoeffs0(i - 1, 0, 0)) +
                            ecoeffs1(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            ecoeffs1(i, 0, t) = one_o_2p * ecoeffs1(i - 1, 0, t - 1) -
                                (b / p) * (R * ecoeffs1(i - 1, 0, t) + ecoeffs0(i - 1, 0, t)) +
                                (t + 1) * ecoeffs1(i - 1, 0, t + 1);

        ecoeffs1(i, 0, i) = one_o_2p * ecoeffs1(i - 1, 0, i - 1) -
                            (b / p) * (R * ecoeffs1(i - 1, 0, i) + ecoeffs0(i - 1, 0, i));
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            ecoeffs1(i, j, 0) = (a / p) * (R * ecoeffs1(i, j - 1, 0) + ecoeffs0(i, j - 1, 0)) +
                                ecoeffs1(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                ecoeffs1(i, j, t) = one_o_2p * ecoeffs1(i, j - 1, t - 1) +
                                    (a / p) * (R * ecoeffs1(i, j - 1, t) + ecoeffs0(i, j - 1, t)) +
                                    (t + 1) * ecoeffs1(i, j - 1, t + 1);

            ecoeffs1(i, j, i + j) = one_o_2p * ecoeffs1(i, j - 1, i + j - 1) +
                                    (a / p) * (R * ecoeffs1(i, j - 1, i + j) + ecoeffs0(i, j - 1, i + j));
        }

    return ecoeffs1;
}

lible::vec3d lints::ecoeffsRecurrence2_n2(const double a, const double b, const int la,
                                          const int lb, const double A, const double B,
                                          const vec3d &ecoeffs0, const vec3d &ecoeffs1)
{
    vec3d ecoeffs2(Fill(0), la + 1, lb + 1, la + lb + 1);

    const double p = a + b;
    const double one_o_2p = 1.0 / (2 * p);
    const double R = A - B;

    ecoeffs2(0, 0, 0) = (-2 * a * b / p) * (R * ecoeffs1(0, 0, 0) + ecoeffs0(0, 0, 0));

    for (int i = 1; i <= la; i++)
    {
        ecoeffs2(i, 0, 0) = -(b / p) * (R * ecoeffs2(i - 1, 0, 0) + 2 * ecoeffs1(i - 1, 0, 0)) +
                            ecoeffs2(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            ecoeffs2(i, 0, t) = one_o_2p * ecoeffs2(i - 1, 0, t - 1) -
                                (b / p) * (R * ecoeffs2(i - 1, 0, t) + 2 * ecoeffs1(i - 1, 0, t)) +
                                (t + 1) * ecoeffs2(i - 1, 0, t + 1);

        ecoeffs2(i, 0, i) = one_o_2p * ecoeffs2(i - 1, 0, i - 1) -
                            (b / p) * (R * ecoeffs2(i - 1, 0, i) + 2 * ecoeffs1(i - 1, 0, i));
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            ecoeffs2(i, j, 0) = (a / p) * (R * ecoeffs2(i, j - 1, 0) + 2 * ecoeffs1(i, j - 1, 0)) +
                                ecoeffs2(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                ecoeffs2(i, j, t) = one_o_2p * ecoeffs2(i, j - 1, t - 1) +
                                    (a / p) * (R * ecoeffs2(i, j - 1, t) + 2 * ecoeffs1(i, j - 1, t)) +
                                    (t + 1) * ecoeffs2(i, j - 1, t + 1);

            ecoeffs2(i, j, i + j) = one_o_2p * ecoeffs2(i, j - 1, i + j - 1) +
                                    (a / p) * (R * ecoeffs2(i, j - 1, i + j) + 2 * ecoeffs1(i, j - 1, i + j));
        }

    return ecoeffs2;
}

array<lible::vec3d, 3> lints::ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                   const int lb, const double *xyz_a,
                                                   const double *xyz_b)
{
    const double p = a + b;
    const double mu = a * b / p;

    array<double, 3> xyz_p{};
    array<double, 3> Kab{};
    for (int i = 0; i < 3; i++)
    {
        xyz_p[i] = (a * xyz_a[i] + b * xyz_b[i]) / p;
        Kab[i] = std::exp(-mu * std::pow(xyz_a[i] - xyz_b[i], 2));
    }

    array<double, 3> xyz_pa{};
    array<double, 3> xyz_pb{};
    for (int i = 0; i < 3; i++)
    {
        xyz_pa[i] = xyz_p[i] - xyz_a[i];
        xyz_pb[i] = xyz_p[i] - xyz_b[i];
    }

    vec3d ecoeffs_x = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[0], xyz_pb[0], Kab[0]);
    vec3d ecoeffs_y = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[1], xyz_pb[1], Kab[1]);
    vec3d ecoeffs_z = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[2], xyz_pb[2], Kab[2]);

    return {ecoeffs_x, ecoeffs_y, ecoeffs_z};
}

array<lible::vec3d, 3> lints::ecoeffsPrimitivePair_n1(const double a, const double b, const int la,
                                                      const int lb, const double *xyz_a,
                                                      const double *xyz_b,
                                                      const array<vec3d, 3> &ecoeffs0)
{
    vec3d ecoeffs1_x = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[0], xyz_b[0], ecoeffs0[0]);
    vec3d ecoeffs1_y = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[1], xyz_b[1], ecoeffs0[1]);
    vec3d ecoeffs1_z = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[2], xyz_b[2], ecoeffs0[2]);

    return {ecoeffs1_x, ecoeffs1_y, ecoeffs1_z};
}

array<lible::vec3d, 3> lints::ecoeffsPrimitivePair_n2(const double a, const double b, const int la,
                                                      const int lb, const double *xyz_a,
                                                      const double *xyz_b,
                                                      const array<vec3d, 3> &ecoeffs0,
                                                      const array<vec3d, 3> &ecoeffs1)
{
    vec3d ecoeffs2_x = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[0], xyz_b[0], ecoeffs0[0], ecoeffs1[1]);
    vec3d ecoeffs2_y = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[1], xyz_b[1], ecoeffs0[1], ecoeffs1[1]);
    vec3d ecoeffs2_z = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[2], xyz_b[2], ecoeffs0[2], ecoeffs1[2]);

    return {ecoeffs2_x, ecoeffs2_y, ecoeffs2_z};
}

array<lible::vec2d, 3> lints::ecoeffsPrimitive(const double a, const int l)
{
    double one_o_2a = 1.0 / (2 * a);
    vec2d ecoeffs_x = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_y = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_z = ecoeffsRecurrence1(one_o_2a, l);

    return {ecoeffs_x, ecoeffs_y, ecoeffs_z};
}

vector<vector<double>> lints::ecoeffsShell(const int l, const vector<double> &exps)
{
    const int dim = numCartesians(l);
    const size_t cdepth = exps.size();

    constexpr array<double, 3> xyz_a{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    vector<vector<double>> ecoeffs_out(cdepth * cdepth, vector<double>(dim * dim, 0));
    for (size_t ia = 0, iab = 0; ia < cdepth; ia++)
        for (size_t ib = 0; ib < cdepth; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] =
                    ecoeffsPrimitivePair(a, b, l, l, &xyz_a[0], &xyz_a[0]);

            for (const auto &[i, j, k, mu]: cart_exps_a)
                for (const auto &[i_, j_, k_, nu]: cart_exps_a)
                {
                    int munu = mu * dim + nu;
                    ecoeffs_out[iab][munu] = ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) *
                                             ecoeffs_z(k, k_, 0);
                }
        }

    return ecoeffs_out;
}

vector<double> lints::ecoeffsSHARK(const ShellData &sh_data, const bool transpose)
{
    const int l = sh_data.l;
    const int n_sph = numSphericals(l);
    const int n_hermite = numHermites(l);
    const int n_ecoeffs_pp = n_sph * n_hermite;
    const int n_ecoeffs = sh_data.n_primitives * n_ecoeffs_pp;

    const auto &cart_exps = cartExps(l);
    const vec3i tuv_poss = getHermiteGaussianPositions(l);
    const vector<tuple<int, int, double>> sph_trafo = sphericalTrafo(l);

    vector<double> ecoeffs(n_ecoeffs, 0);
#pragma omp parallel for
    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int cdepth = sh_data.cdepths[ishell];
        int cofs = sh_data.coffsets[ishell];

        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];

        int ofs_norm = sh_data.offsets_norms[ishell];
        const double *norms = &sh_data.norms[ofs_norm];

        for (int ia = 0; ia < cdepth; ia++)
        {
            double a = sh_data.exps[cofs + ia];
            double d = sh_data.coeffs[cofs + ia];

            auto [Ex, Ey, Ez] = ecoeffsPrimitive(a, l);

            int ofs = offset_ecoeffs + ia * n_sph * n_hermite;
            for (const auto &[a, a_, val]: sph_trafo)
            {
                auto [i, j, k] = cart_exps[a_];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            double ecoeff = norms[a] * d * Ex(i, t) * Ey(j, u) * Ez(k, v) * val;

                            int tuv = tuv_poss(t, u, v);

                            int idx;
                            if (transpose)
                                idx = ofs + tuv * n_sph + a;
                            else
                                idx = ofs + a * n_hermite + tuv;

                            ecoeffs[idx] += ecoeff;
                        }
            }
        }
    }

    return ecoeffs;
}

vector<double> lints::ecoeffsSHARK(const ShellPairData &sp_data, const bool transpose)
{
    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int n_cart_b = numCartesians(lb);
    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_hermite = numHermites(lab);
    const int n_ecoeffs = n_sph_ab * n_hermite;

    const vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    const vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);
    const vec3i tuv_poss = getHermiteGaussianPositions(lab);

    vector<double> ecoeffs(n_ecoeffs * sp_data.n_prim_pairs, 0);
#pragma omp parallel for
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        const int cofs_a = sp_data.coffsets[2 * ipair + 0];
        const int cofs_b = sp_data.coffsets[2 * ipair + 1];

        const double *xyz_a = &sp_data.coords[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords[6 * ipair + 3];

        int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
        int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
        const double *norms_a = &sp_data.norms[ofs_norm_a];
        const double *norms_b = &sp_data.norms[ofs_norm_b];

        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                double a = sp_data.exps[cofs_a + ia];
                double b = sp_data.exps[cofs_b + ib];

                auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

                vec3d ecoeffs_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite);
                for (auto &[a, a_, val]: sph_trafo_a)
                    for (size_t b_ = 0; b_ < cart_exps_b.size(); b_++)
                    {
                        auto [i, j, k] = cart_exps_a[a_];
                        auto [i_, j_, k_] = cart_exps_b[b_];
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);

                                    ecoeffs_ppair_sc(a, b_, tuv) +=
                                            val * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }
                    }

                double da = sp_data.coeffs[cofs_a + ia];
                double db = sp_data.coeffs[cofs_b + ib];
                double dadb = da * db;

                int ofs = offset_ecoeffs + iab * n_ecoeffs;
                for (auto &[b, b_, val]: sph_trafo_b)
                    for (int a = 0; a < n_sph_a; a++)
                        for (int tuv = 0; tuv < n_hermite; tuv++)
                        {
                            int ab = a * n_sph_b + b;

                            int idx;
                            if (transpose)
                                idx = ofs + tuv * n_sph_ab + ab;
                            else
                                idx = ofs + ab * n_hermite + tuv;

                            ecoeffs[idx] += norms_a[a] * norms_b[b] * dadb * val *
                                    ecoeffs_ppair_sc(a, b_, tuv);
                        }
            }
    }

    return ecoeffs;
}

vector<double> lints::ecoeffsD1SHARK(const ShellPairData &sp_data, const bool transpose)
{
    const int la = sp_data.la;
    const int lb = sp_data.lb;
    const int lab = la + lb;
    const int n_cart_b = numCartesians(lb);
    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_sph_ab = n_sph_a * n_sph_b;
    const int n_hermite_ab = numHermites(lab);
    const int n_ecoeffs = n_sph_ab * n_hermite_ab;
    const int n_ecoeffs_prims = 3 * n_ecoeffs * sp_data.n_prim_pairs;

    const vec3i tuv_poss = getHermiteGaussianPositions(lab);
    const vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    const vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    vector<double> ecoeffs_100_010_001(n_ecoeffs_prims, 0);
#pragma omp parallel for
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        const int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        const int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        const int cofs_a = sp_data.coffsets[2 * ipair + 0];
        const int cofs_b = sp_data.coffsets[2 * ipair + 1];

        const double *xyz_a = &sp_data.coords[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords[6 * ipair + 3];

        const int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
        const int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
        const double *norms_a = &sp_data.norms[ofs_norm_a];
        const double *norms_b = &sp_data.norms[ofs_norm_b];

        int offset_ecoeffs = sp_data.offsets_ecoeffs_deriv1[ipair];

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                double a = sp_data.exps[cofs_a + ia];
                double b = sp_data.exps[cofs_b + ib];

                auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

                auto [E1x, E1y, E1z] =
                        ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b, {Ex, Ey, Ez});

                vec3d ecoeffs100_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
                vec3d ecoeffs010_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
                vec3d ecoeffs001_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
                for (auto &[a, a_, val]: sph_trafo_a)
                    for (size_t b_ = 0; b_ < cart_exps_b.size(); b_++)
                    {
                        auto [i, j, k] = cart_exps_a[a_];
                        auto [i_, j_, k_] = cart_exps_b[b_];
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);

                                    ecoeffs100_ppair_sc(a, b_, tuv) +=
                                            val * E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                    ecoeffs010_ppair_sc(a, b_, tuv) +=
                                            val * Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                                    ecoeffs001_ppair_sc(a, b_, tuv) +=
                                            val * Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);
                                }
                    }

                double da = sp_data.coeffs[cofs_a + ia];
                double db = sp_data.coeffs[cofs_b + ib];
                double dadb = da * db;

                int ofs_100 = offset_ecoeffs + (3 * iab + 0) * n_ecoeffs;
                int ofs_010 = offset_ecoeffs + (3 * iab + 1) * n_ecoeffs;
                int ofs_001 = offset_ecoeffs + (3 * iab + 2) * n_ecoeffs;
                for (auto &[b, b_, val]: sph_trafo_b)
                    for (int a = 0; a < n_sph_a; a++)
                        for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                        {
                            int ab = a * n_sph_b + b;

                            int idx_100, idx_010, idx_001;
                            if (transpose)
                            {
                                idx_100 = ofs_100 + tuv * n_sph_ab + ab;
                                idx_010 = ofs_010 + tuv * n_sph_ab + ab;
                                idx_001 = ofs_001 + tuv * n_sph_ab + ab;
                            }
                            else
                            {
                                idx_100 = ofs_100 + ab * n_hermite_ab + tuv;
                                idx_010 = ofs_010 + ab * n_hermite_ab + tuv;
                                idx_001 = ofs_001 + ab * n_hermite_ab + tuv;
                            }

                            double NaNb = norms_a[a] * norms_b[b];

                            ecoeffs_100_010_001[idx_100] += NaNb * dadb * val * ecoeffs100_ppair_sc(a, b_, tuv);
                            ecoeffs_100_010_001[idx_010] += NaNb * dadb * val * ecoeffs010_ppair_sc(a, b_, tuv);
                            ecoeffs_100_010_001[idx_001] += NaNb * dadb * val * ecoeffs001_ppair_sc(a, b_, tuv);
                        }
            }
    }

    return ecoeffs_100_010_001;
}