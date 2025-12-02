#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/utils.hpp>

#include <array>
#include <tuple>

namespace lints = lible::ints;

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

std::array<lible::vec3d, 3> lints::ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                        const int lb, const double *xyz_a,
                                                        const double *xyz_b)
{
    const double p = a + b;
    const double mu = a * b / p;

    std::array<double, 3> xyz_p{};
    std::array<double, 3> Kab{};
    for (int i = 0; i < 3; i++)
    {
        xyz_p[i] = (a * xyz_a[i] + b * xyz_b[i]) / p;
        Kab[i] = std::exp(-mu * std::pow(xyz_a[i] - xyz_b[i], 2));
    }

    std::array<double, 3> xyz_pa{};
    std::array<double, 3> xyz_pb{};
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

std::array<lible::vec3d, 3> lints::ecoeffsPrimitivePair_n1(const double a, const double b,
                                                           const int la, const int lb,
                                                           const double *xyz_a, const double *xyz_b,
                                                           const std::array<vec3d, 3> &ecoeffs0)
{
    vec3d ecoeffs1_x = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[0], xyz_b[0], ecoeffs0[0]);
    vec3d ecoeffs1_y = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[1], xyz_b[1], ecoeffs0[1]);
    vec3d ecoeffs1_z = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[2], xyz_b[2], ecoeffs0[2]);

    return {ecoeffs1_x, ecoeffs1_y, ecoeffs1_z};
}

std::array<lible::vec3d, 3> lints::ecoeffsPrimitivePair_n2(const double a, const double b,
                                                           const int la, const int lb,
                                                           const double *xyz_a, const double *xyz_b,
                                                           const std::array<vec3d, 3> &ecoeffs0,
                                                           const std::array<vec3d, 3> &ecoeffs1)
{
    vec3d ecoeffs2_x = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[0], xyz_b[0], ecoeffs0[0], ecoeffs1[1]);
    vec3d ecoeffs2_y = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[1], xyz_b[1], ecoeffs0[1], ecoeffs1[1]);
    vec3d ecoeffs2_z = ecoeffsRecurrence2_n2(a, b, la, lb, xyz_a[2], xyz_b[2], ecoeffs0[2], ecoeffs1[2]);

    return {ecoeffs2_x, ecoeffs2_y, ecoeffs2_z};
}

std::array<lible::vec2d, 3> lints::ecoeffsPrimitive(const double a, const int l)
{
    double one_o_2a = 1.0 / (2 * a);
    vec2d ecoeffs_x = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_y = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_z = ecoeffsRecurrence1(one_o_2a, l);

    return {ecoeffs_x, ecoeffs_y, ecoeffs_z};
}

std::vector<std::vector<double>> lints::ecoeffsShell(const int l, const std::vector<double> &exps)
{
    const int dim = numCartesians(l);
    const size_t cdepth = exps.size();

    constexpr std::array<double, 3> xyz_a{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    std::vector<std::vector<double>> ecoeffs_out(cdepth * cdepth, std::vector<double>(dim * dim, 0));
    for (size_t ia = 0, iab = 0; ia < cdepth; ia++)
        for (size_t ib = 0; ib < cdepth; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] =
                    ecoeffsPrimitivePair(a, b, l, l, &xyz_a[0], &xyz_a[0]);

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_a)
                {
                    int munu = mu * dim + nu;
                    ecoeffs_out[iab][munu] = ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) *
                                             ecoeffs_z(k, k_, 0);
                }
        }

    return ecoeffs_out;
}

std::vector<double> lints::ecoeffsSHARK(const ShellData &sh_data, const bool transpose)
{
    if (omp_in_parallel() == true)
        throw std::runtime_error("ecoeffsSHARK(): cannot be called inside a parallel region");

    int l = sh_data.l_;
    int n_sph = numSphericals(l);
    int n_hermite = numHermites(l);
    int n_ecoeffs_pp = n_sph * n_hermite;
    size_t n_ecoeffs = sh_data.n_primitives_ * n_ecoeffs_pp;

    const auto &ijk = cartExps(l);
    vec3i tuv_poss = getHermiteGaussianPositions(l);
    std::vector<std::tuple<int, int, double>> sph_trafo = sphericalTrafo(l);

    std::vector<double> ecoeffs(n_ecoeffs, 0);
#pragma omp parallel for
    for (size_t ishell = 0; ishell < sh_data.n_shells_; ishell++)
    {
        size_t cdepth = sh_data.cdepths_[ishell];
        size_t cofs = sh_data.coffsets_[ishell];

        size_t offset_ecoeffs = sh_data.offsets_ecoeffs_[ishell];

        size_t ofs_norm = sh_data.offsets_norms_[ishell];
        const double *norms = &sh_data.norms_[ofs_norm];

        for (size_t ia = 0; ia < cdepth; ia++)
        {
            double a = sh_data.exps_[cofs + ia];
            double d = sh_data.coeffs_[cofs + ia];

            auto [Ex, Ey, Ez] = ecoeffsPrimitive(a, l);

            size_t ofs = offset_ecoeffs + ia * n_sph * n_hermite;
            for (const auto &[mu, mu_, val] : sph_trafo)
            {
                auto [i, j, k] = ijk[mu_];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            double ecoeff = norms[mu] * d * Ex(i, t) * Ey(j, u) * Ez(k, v) * val;

                            int tuv = tuv_poss(t, u, v);

                            size_t idx;
                            if (transpose)
                                idx = ofs + tuv * n_sph + mu;
                            else
                                idx = ofs + mu * n_hermite + tuv;

                            ecoeffs[idx] += ecoeff;
                        }
            }
        }
    }

    return ecoeffs;
}

std::vector<double> lints::ecoeffsSHARK(const ShellPairData &sp_data, const bool transpose)
{
    if (omp_in_parallel() == true)
        throw std::runtime_error("ecoeffsSHARK(): cannot be called inside a parallel region");

    int la = sp_data.la_;
    int lb = sp_data.lb_;
    int lab = la + lb;
    int n_cart_b = numCartesians(lb);
    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_hermite = numHermites(lab);
    int n_ecoeffs = n_sph_ab * n_hermite;

    std::vector<std::tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    std::vector<std::tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);
    vec3i tuv_poss = getHermiteGaussianPositions(lab);

    std::vector<double> ecoeffs(n_ecoeffs * sp_data.n_ppairs_, 0);
#pragma omp parallel for
    for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
    {
        const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

        size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
        size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
        const double *norms_a = &sp_data.norms_[ofs_norm_a];
        const double *norms_b = &sp_data.norms_[ofs_norm_b];

        size_t ofs_prim = sp_data.offsets_primitives_[ipair];
        const double *exps = &sp_data.exps_[ofs_prim];
        const double *coeffs = &sp_data.coeffs_[ofs_prim];

        size_t offset_ecoeffs = sp_data.offsets_ecoeffs_[ipair];

        for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
        {
            double a = exps[iab * 2];
            double b = exps[iab * 2 + 1];

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            // first trafo
            vec3d ecoeffs_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite);
            for (auto &[mu, mu_, val] : sph_trafo_a)
                for (size_t nu_ = 0; nu_ < cart_exps_b.size(); nu_++)
                {
                    auto [i, j, k] = cart_exps_a[mu_];
                    auto [i_, j_, k_] = cart_exps_b[nu_];
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                int tuv = tuv_poss(t, u, v);

                                ecoeffs_ppair_sc(mu, nu_, tuv) +=
                                        val * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                            }
                }

            // second trafo
            double da = coeffs[iab * 2];
            double db = coeffs[iab * 2 + 1];
            double dadb = da * db;

            size_t ofs = offset_ecoeffs + iab * n_ecoeffs;
            for (auto &[nu, nu_, val] : sph_trafo_b)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int tuv = 0; tuv < n_hermite; tuv++)
                    {
                        int munu = mu * n_sph_b + nu;

                        size_t idx;
                        if (transpose)
                            idx = ofs + tuv * n_sph_ab + munu;
                        else
                            idx = ofs + munu * n_hermite + tuv;

                        ecoeffs[idx] += norms_a[mu] * norms_b[nu] * dadb * val *
                                ecoeffs_ppair_sc(mu, nu_, tuv);
                    }
        }
    }

    return ecoeffs;
}

std::vector<double> lints::ecoeffsD1SHARK(const ShellPairData &sp_data, const bool transpose)
{
    if (omp_in_parallel() == true)
        throw std::runtime_error("ecoeffsD1SHARK(): cannot be called inside a parallel region");

    int la = sp_data.la_;
    int lb = sp_data.lb_;
    int lab = la + lb;
    int n_cart_b = numCartesians(lb);
    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_hermite_ab = numHermites(lab);
    int n_ecoeffs = n_sph_ab * n_hermite_ab;
    size_t n_ecoeffs_prims = 3 * n_ecoeffs * sp_data.n_ppairs_;

    vec3i tuv_poss = getHermiteGaussianPositions(lab);
    std::vector<std::tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    std::vector<std::tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    std::vector<double> ecoeffs_100_010_001(n_ecoeffs_prims, 0);
#pragma omp parallel for
    for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
    {
        const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

        size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
        size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
        const double *norms_a = &sp_data.norms_[ofs_norm_a];
        const double *norms_b = &sp_data.norms_[ofs_norm_b];

        size_t ofs_prim = sp_data.offsets_primitives_[ipair];
        const double *exps = &sp_data.exps_[ofs_prim];
        const double *coeffs = &sp_data.coeffs_[ofs_prim];

        size_t offset_ecoeffs = sp_data.offsets_ecoeffs_deriv1_[ipair];

        for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
        {
            double a = exps[iab * 2];
            double b = exps[iab * 2 + 1];

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] =
                    ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b, {Ex, Ey, Ez});

            // first trafo
            vec3d ecoeffs100_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
            vec3d ecoeffs010_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
            vec3d ecoeffs001_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
            for (auto &[mu, mu_, val] : sph_trafo_a)
                for (size_t nu_ = 0; nu_ < cart_exps_b.size(); nu_++)
                {
                    auto [i, j, k] = cart_exps_a[mu_];
                    auto [i_, j_, k_] = cart_exps_b[nu_];
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                int tuv = tuv_poss(t, u, v);

                                ecoeffs100_ppair_sc(mu, nu_, tuv) +=
                                        val * E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                ecoeffs010_ppair_sc(mu, nu_, tuv) +=
                                        val * Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                                ecoeffs001_ppair_sc(mu, nu_, tuv) +=
                                        val * Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);
                            }
                }

            // second trafo
            double da = coeffs[iab * 2];
            double db = coeffs[iab * 2 + 1];
            double dadb = da * db;

            size_t ofs_100 = offset_ecoeffs + (3 * iab + 0) * n_ecoeffs;
            size_t ofs_010 = offset_ecoeffs + (3 * iab + 1) * n_ecoeffs;
            size_t ofs_001 = offset_ecoeffs + (3 * iab + 2) * n_ecoeffs;
            for (auto &[nu, nu_, val] : sph_trafo_b)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                    {
                        int munu = mu * n_sph_b + nu;

                        size_t idx_100, idx_010, idx_001;
                        if (transpose)
                        {
                            idx_100 = ofs_100 + tuv * n_sph_ab + munu;
                            idx_010 = ofs_010 + tuv * n_sph_ab + munu;
                            idx_001 = ofs_001 + tuv * n_sph_ab + munu;
                        }
                        else
                        {
                            idx_100 = ofs_100 + munu * n_hermite_ab + tuv;
                            idx_010 = ofs_010 + munu * n_hermite_ab + tuv;
                            idx_001 = ofs_001 + munu * n_hermite_ab + tuv;
                        }

                        double NaNb = norms_a[mu] * norms_b[nu];

                        ecoeffs_100_010_001[idx_100] += NaNb * dadb * val * ecoeffs100_ppair_sc(mu, nu_, tuv);
                        ecoeffs_100_010_001[idx_010] += NaNb * dadb * val * ecoeffs010_ppair_sc(mu, nu_, tuv);
                        ecoeffs_100_010_001[idx_001] += NaNb * dadb * val * ecoeffs001_ppair_sc(mu, nu_, tuv);
                    }
        }
    }

    return ecoeffs_100_010_001;
}
