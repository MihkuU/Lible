#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <array>

namespace LI = lible::ints;

using std::array, std::vector;

void LI::ecoeffsRecurrence2(const double a, const double b, const double PA, const double PB,
                            const double one_o_2p, const int la, const int lb, vec3d &ecoeffs)
{
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
}

void LI::ecoeffsRecurrence1(const double one_o_2a, const int l, vec2d &ecoeffs)
{
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
}

void LI::ecoeffsPrimitive(const double a, const int l, vec2d &ecoeffs_x, vec2d &ecoeffs_y,
                          vec2d &ecoeffs_z)
{
    ecoeffs_x.set(0);
    ecoeffs_y.set(0);
    ecoeffs_z.set(0);

    ecoeffs_x(0, 0) = 1;
    ecoeffs_y(0, 0) = 1;
    ecoeffs_z(0, 0) = 1;

    double one_o_2a = 1.0 / (2 * a);

    ecoeffsRecurrence1(one_o_2a, l, ecoeffs_x);
    ecoeffsRecurrence1(one_o_2a, l, ecoeffs_y);
    ecoeffsRecurrence1(one_o_2a, l, ecoeffs_z);
}

void LI::ecoeffsPrimitivePair(const double a, const double b, const int la, const int lb,
                              const array<double, 3> &xyz_a, const array<double, 3> &xyz_b,
                              const array<double, 3> &Kab, vec3d &ecoeffs_x, vec3d &ecoeffs_y,
                              vec3d &ecoeffs_z)
{
    ecoeffs_x.set(0);
    ecoeffs_y.set(0);
    ecoeffs_z.set(0);

    ecoeffs_x(0, 0, 0) = Kab[0];
    ecoeffs_y(0, 0, 0) = Kab[1];
    ecoeffs_z(0, 0, 0) = Kab[2];

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                           (a * xyz_a[2] + b * xyz_b[2]) / p};

    array<double, 3> xyz_pa{xyz_p[0] - xyz_a[0], xyz_p[1] - xyz_a[1], xyz_p[2] - xyz_a[2]};
    array<double, 3> xyz_pb{xyz_p[0] - xyz_b[0], xyz_p[1] - xyz_b[1], xyz_p[2] - xyz_b[2]};

    ecoeffsRecurrence2(a, b, xyz_pa[0], xyz_pb[0], one_o_2p, la, lb, ecoeffs_x);
    ecoeffsRecurrence2(a, b, xyz_pa[1], xyz_pb[1], one_o_2p, la, lb, ecoeffs_y);
    ecoeffsRecurrence2(a, b, xyz_pa[2], xyz_pb[2], one_o_2p, la, lb, ecoeffs_z);
}

void LI::ecoeffsShell(const int l, const vector<double> &exps, vector<arma::dmat> &ecoeffs_out)
{
    int dim = numCartesians(l);
    size_t k = exps.size();

    ecoeffs_out.resize(k * k, arma::zeros(dim, dim));

    lible::vec3d ecoeffs_x(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d ecoeffs_y(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d ecoeffs_z(l + 1, l + 1, 2 * l + 1, 0);

    array<double, 3> xyz_a{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            std::array<double, 3> Kab{1, 1, 1};

            ecoeffsPrimitivePair(a, b, l, l, xyz_a, xyz_a, Kab, ecoeffs_x, ecoeffs_y, ecoeffs_z);

            for (const auto [i, j, k, mu] : cart_exps_a)
                for (const auto [i_, j_, k_, nu] : cart_exps_a)
                    ecoeffs_out[iab](mu, nu) = ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) *
                                               ecoeffs_z(k, k_, 0);
        }
}

void LI::ecoeffsShellPairs3D(const int la, const int lb, const ShellPairData &sp_data,
                             vector<vector<vec3d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec3d> ecoeffs_pair(dim_a * dim_b, vec3d(3, la + 1, lb + 1, 0));
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](0, i, j) = Ex(i, j, 0);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](1, i, j) = Ey(i, j, 0);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        ecoeffs_pair[iab](2, i, j) = Ez(i, j, 0);
            }

        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LI::ecoeffsShellPairs4D(const int la, const int lb, const ShellPairData &sp_data,
                             vector<vector<vec4d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec4d> ecoeffs_pair(dim_a * dim_b, vec4d(3, la + 1, lb + 1, la + lb + 1, 0));
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](0, i, j, t) = Ex(i, j, t);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](1, i, j, t) = Ey(i, j, t);

                for (int i = 0; i <= la; i++)
                    for (int j = 0; j <= lb; j++)
                        for (int t = 0; t <= i + j; t++)
                            ecoeffs_pair[iab](2, i, j, t) = Ez(i, j, t);
            }

        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LI::ecoeffsShellsSpherical(const int l, const ShellData &sh_data, vector<double> &ecoeffs_out)
{
    lible::vec2d ecoeffs_x(l + 1, l + 1, 0);
    lible::vec2d ecoeffs_y(l + 1, l + 1, 0);
    lible::vec2d ecoeffs_z(l + 1, l + 1, 0);

    arma::dmat sph_trafo = returnSphericalTrafo(l);

    const auto &cart_exps = cartExps(l);

    int dim_cart = numCartesians(l);
    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    int n_ecoeffs = dim_sph * dim_tuv;

    vec3i tuv_poss = returnHermiteGaussianPositions(l);

    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int dim = sh_data.cdepths[ishell];
        int pos = sh_data.coffsets[ishell];

        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];

        vec2d ecoeffs_c(dim_cart, dim_tuv, 0);
        for (int i = 0; i < dim; i++)
        {
            double a = sh_data.exps[pos + i];
            ecoeffsPrimitive(a, l, ecoeffs_x, ecoeffs_y, ecoeffs_z);

            for (size_t mu = 0; mu < cart_exps.size(); mu++)
            {
                auto [i, j, k] = cart_exps[mu];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            int tuv = tuv_poss(t, u, v); // TODO: do this more effectively by looping over (t, u, v)-triplets directly
                            ecoeffs_c(mu, tuv) = ecoeffs_x(i, t) * ecoeffs_y(j, u) *
                                                 ecoeffs_z(k, v);
                        }
            }

            double d = sh_data.coeffs[pos + i];
            int pos_ecoeffs = offset_ecoeffs + i * n_ecoeffs;
            for (int mu = 0; mu < dim_sph; mu++)
                for (int mu_ = 0; mu_ < dim_cart; mu_++)
                    for (int tuv = 0; tuv < dim_tuv; tuv++)
                    {
                        double ecoeff = d * ecoeffs_c(mu_, tuv) * sph_trafo(mu, mu_);

                        int idx = pos_ecoeffs + mu * dim_tuv + tuv;
                        ecoeffs_out[idx] += ecoeff;
                    }
        }
    }
}

void LI::ecoeffsShellsSpherical(const int l, const ShellData &sh_data, vector<double> &ecoeffs_out,
                                vector<double> &ecoeffs_tsp_out)
{
    std::fill(ecoeffs_out.begin(), ecoeffs_out.end(), 0);
    std::fill(ecoeffs_tsp_out.begin(), ecoeffs_tsp_out.end(), 0);

    lible::vec2d ecoeffs_x(l + 1, l + 1, 0);
    lible::vec2d ecoeffs_y(l + 1, l + 1, 0);
    lible::vec2d ecoeffs_z(l + 1, l + 1, 0);

    arma::dmat sph_trafo = returnSphericalTrafo(l);

    const auto &cart_exps = cartExps(l);

    int dim_cart = numCartesians(l);
    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    int n_ecoeffs = dim_sph * dim_tuv;

    vec3i tuv_poss = returnHermiteGaussianPositions(l);

    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int dim = sh_data.cdepths[ishell];
        int pos = sh_data.coffsets[ishell];

        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];

        vec2d ecoeffs_c(dim_cart, dim_tuv, 0);
        for (int i = 0; i < dim; i++)
        {
            double a = sh_data.exps[pos + i];
            ecoeffsPrimitive(a, l, ecoeffs_x, ecoeffs_y, ecoeffs_z);

            for (size_t mu = 0; mu < cart_exps.size(); mu++)
            {
                auto [i, j, k] = cart_exps[mu];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            int tuv = tuv_poss(t, u, v); // TODO: do this more effectively by looping over (t, u, v)-triplets directly
                            ecoeffs_c(mu, tuv) = ecoeffs_x(i, t) * ecoeffs_y(j, u) *
                                                 ecoeffs_z(k, v);
                        }
            }

            double d = sh_data.coeffs[pos + i];
            int pos_ecoeffs = offset_ecoeffs + i * n_ecoeffs;
            for (int mu = 0; mu < dim_sph; mu++)
                for (int mu_ = 0; mu_ < dim_cart; mu_++)
                    for (int tuv = 0; tuv < dim_tuv; tuv++)
                    {
                        double ecoeff = d * ecoeffs_c(mu_, tuv) * sph_trafo(mu, mu_);

                        int idx = pos_ecoeffs + mu * dim_tuv + tuv;
                        ecoeffs_out[idx] += ecoeff;

                        int idx_tsp = pos_ecoeffs + tuv * dim_sph + mu;
                        ecoeffs_tsp_out[idx_tsp] += ecoeff;
                    }
        }
    }
}

void LI::ecoeffsSPsSpherical(const int la, const int lb, const ShellPairData &sp_data,
                             vector<double> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    arma::dmat sph_trafo_a = returnSphericalTrafo(la);
    arma::dmat sph_trafo_b = returnSphericalTrafo(lb);

    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    int dim_a_cart = numCartesians(la);
    int dim_b_cart = numCartesians(lb);
    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vec3i tuv_poss = returnHermiteGaussianPositions(lab);

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];

        vec3d ecoeffs_ppair_cc(dim_a_cart, dim_b_cart, dim_tuv, 0);
        vec3d ecoeffs_ppair_sc(dim_a_sph, dim_b_cart, dim_tuv, 0);
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                ecoeffs_ppair_cc.set(0);
                for (size_t mu = 0; mu < cart_exps_a.size(); mu++)
                {
                    auto [i, j, k] = cart_exps_a[mu];
                    for (size_t nu = 0; nu < cart_exps_b.size(); nu++)
                    {
                        auto [i_, j_, k_] = cart_exps_b[nu];
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);
                                    ecoeffs_ppair_cc(mu, nu, tuv) =
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }
                    }
                }

                ecoeffs_ppair_sc.set(0);
                for (int mu = 0; mu < dim_a_sph; mu++)
                    for (int mu_ = 0; mu_ < dim_a_cart; mu_++)
                        for (int nu = 0; nu < dim_b_cart; nu++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_ppair_sc(mu, nu, tuv) +=
                                    ecoeffs_ppair_cc(mu_, nu, tuv) * sph_trafo_a(mu, mu_);

                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];
                int pos = offset_ecoeffs + iab * n_ecoeffs;
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int nu_ = 0; nu_ < dim_b_cart; nu_++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                            {
                                double ecoeff = da * db * ecoeffs_ppair_sc(mu, nu_, tuv) *
                                                sph_trafo_b(nu, nu_);

                                int idx = pos + munu * dim_tuv + tuv;
                                ecoeffs_out[idx] += ecoeff;
                            }
            }
    }
}

void LI::ecoeffsSPsSpherical(const int la, const int lb, const ShellPairData &sp_data,
                             vector<double> &ecoeffs_out, vector<double> &ecoeffs_tsp_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    arma::dmat sph_trafo_a = returnSphericalTrafo(la);
    arma::dmat sph_trafo_b = returnSphericalTrafo(lb);

    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    int dim_a_cart = numCartesians(la);
    int dim_b_cart = numCartesians(lb);
    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vec3i tuv_poss = returnHermiteGaussianPositions(lab);

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];

        vec3d ecoeffs_ppair_cc(dim_a_cart, dim_b_cart, dim_tuv, 0);
        vec3d ecoeffs_ppair_sc(dim_a_sph, dim_b_cart, dim_tuv, 0);
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                          std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                          std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

                ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b, Kab, Ex, Ey, Ez);

                ecoeffs_ppair_cc.set(0);
                for (size_t mu = 0; mu < cart_exps_a.size(); mu++)
                {
                    auto [i, j, k] = cart_exps_a[mu];
                    for (size_t nu = 0; nu < cart_exps_b.size(); nu++)
                    {
                        auto [i_, j_, k_] = cart_exps_b[nu];
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);
                                    ecoeffs_ppair_cc(mu, nu, tuv) =
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }
                    }
                }

                ecoeffs_ppair_sc.set(0);
                for (int mu = 0; mu < dim_a_sph; mu++)
                    for (int mu_ = 0; mu_ < dim_a_cart; mu_++)
                        for (int nu = 0; nu < dim_b_cart; nu++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_ppair_sc(mu, nu, tuv) +=
                                    ecoeffs_ppair_cc(mu_, nu, tuv) * sph_trafo_a(mu, mu_);

                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];
                int pos = offset_ecoeffs + iab * n_ecoeffs;
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int nu_ = 0; nu_ < dim_b_cart; nu_++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                            {
                                double ecoeff = da * db * ecoeffs_ppair_sc(mu, nu_, tuv) *
                                                sph_trafo_b(nu, nu_);

                                int idx = pos + munu * dim_tuv + tuv;
                                ecoeffs_out[idx] += ecoeff;

                                int idx_tsp = pos + tuv * dim_ab + munu;
                                ecoeffs_tsp_out[idx_tsp] += ecoeff;
                            }
            }
    }
}