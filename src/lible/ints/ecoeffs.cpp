#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <array>
#include <stdexcept>

namespace LI = lible::ints;

using std::array, std::pair, std::vector;

lible::vec3d LI::ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                    const double PA, const double PB, const double Kab)
{
    vec3d ecoeffs(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs(0, 0, 0) = Kab;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

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

lible::vec3d LI::ecoeffsRecurrence2_n1(const double a, const double b, const int la, const int lb,
                                       const double A, const double B, const vec3d &ecoeffs)
{
    vec3d ecoeffs1(la + 1, lb + 1, la + lb + 1, 0);

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);
    double R = A - B;

    ecoeffs1(0, 0, 0) = -2 * (a * b) * R * ecoeffs(0, 0, 0) / p;

    for (int i = 1; i <= la; i++)
    {
        ecoeffs1(i, 0, 0) = -(b / p) * (R * ecoeffs1(i - 1, 0, 0) + ecoeffs(i - 1, 0, 0)) +
                            ecoeffs1(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            ecoeffs1(i, 0, t) = one_o_2p * ecoeffs1(i - 1, 0, t - 1) -
                                (b / p) * (R * ecoeffs1(i - 1, 0, t) + ecoeffs(i - 1, 0, t)) +
                                (t + 1) * ecoeffs1(i - 1, 0, t + 1);

        ecoeffs1(i, 0, i) = one_o_2p * ecoeffs1(i - 1, 0, i - 1) -
                            (b / p) * (R * ecoeffs1(i - 1, 0, i) + ecoeffs(i - 1, 0, i));
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            ecoeffs1(i, j, 0) = -(b / p) * (R * ecoeffs1(i, j - 1, 0) + ecoeffs(i, j - 1, 0)) +
                                ecoeffs1(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                ecoeffs1(i, j, t) = one_o_2p * ecoeffs1(i, j - 1, t - 1) +
                                    (b / p) * R * (ecoeffs1(i, j - 1, t) + ecoeffs(i, j - 1, t)) +
                                    (t + 1) * ecoeffs1(i, j - 1, t + 1);

            ecoeffs1(i, j, i + j) = one_o_2p * ecoeffs1(i, j - 1, i + j - 1) +
                                    (b / p) * R * (ecoeffs1(i, j - 1, i + j) + ecoeffs(i, j - 1, i + j));
        }

    return ecoeffs1;
}

lible::vec2d LI::ecoeffsRecurrence1(const double one_o_2a, const int l)
{
    vec2d ecoeffs(l + 1, l + 1, 0);

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

array<lible::vec3d, 3> LI::ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                const int lb, const double *xyz_a,
                                                const double *xyz_b, const double *Kab)
{
    double p = a + b;
    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                           (a * xyz_a[2] + b * xyz_b[2]) / p};

    array<double, 3> xyz_pa{xyz_p[0] - xyz_a[0], xyz_p[1] - xyz_a[1], xyz_p[2] - xyz_a[2]};
    array<double, 3> xyz_pb{xyz_p[0] - xyz_b[0], xyz_p[1] - xyz_b[1], xyz_p[2] - xyz_b[2]};

    vec3d ecoeffs_x = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[0], xyz_pb[0], Kab[0]);
    vec3d ecoeffs_y = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[1], xyz_pb[1], Kab[1]);
    vec3d ecoeffs_z = ecoeffsRecurrence2(a, b, la, lb, xyz_pa[2], xyz_pb[2], Kab[2]);

    return {ecoeffs_x, ecoeffs_y, ecoeffs_z};
}

array<lible::vec2d, 3> LI::ecoeffsPrimitive(const double a, const int l)
{
    double one_o_2a = 1.0 / (2 * a);

    vec2d ecoeffs_x = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_y = ecoeffsRecurrence1(one_o_2a, l);
    vec2d ecoeffs_z = ecoeffsRecurrence1(one_o_2a, l);

    return {ecoeffs_x, ecoeffs_y, ecoeffs_z};
}

vector<lible::vec3d>
LI::ecoeffsShellPair_Eij0(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                          const double *exps_a, const double *exps_b, const double *xyz_a,
                          const double *xyz_b)
{
    vector<vec3d> ecoeffs(cdepth_a * cdepth_b, vec3d(3, la + 1, lb + 1, 0));
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double mu = a * b / (a + b);

            array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                 std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                 std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb,
                                                                          xyz_a, xyz_b, &Kab[0]);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    ecoeffs[iab](0, i, j) = ecoeffs_x(i, j, 0);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    ecoeffs[iab](1, i, j) = ecoeffs_y(i, j, 0);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    ecoeffs[iab](2, i, j) = ecoeffs_z(i, j, 0);
        }

    return ecoeffs;
}

vector<lible::vec4d>
LI::ecoeffsShellPair_Eijt(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                          const double *exps_a, const double *exps_b, const double *xyz_a,
                          const double *xyz_b)
{
    vector<vec4d> ecoeffs(cdepth_a * cdepth_b, vec4d(3, la + 1, lb + 1, la + lb + 1, 0));
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double mu = a * b / (a + b);

            array<double, 3> Kab{std::exp(-mu * std::pow(xyz_a[0] - xyz_b[0], 2)),
                                 std::exp(-mu * std::pow(xyz_a[1] - xyz_b[1], 2)),
                                 std::exp(-mu * std::pow(xyz_a[2] - xyz_b[2], 2))};

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb,
                                                                          xyz_a, xyz_b, &Kab[0]);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs[iab](0, i, j, t) = ecoeffs_x(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs[iab](1, i, j, t) = ecoeffs_y(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs[iab](2, i, j, t) = ecoeffs_z(i, j, t);
        }

    return ecoeffs;
}

vector<vector<lible::vec3d>>
LI::ecoeffsSPData_Eij0(const int la, const int lb, const ShellPairData &sp_data)
{
    vector<vector<vec3d>> ecoeffs(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec3d> ecoeffs_ipair = ecoeffsShellPair_Eij0(la, lb, cdepth_a, cdepth_b,
                                                            &sp_data.exps[cofs_a],
                                                            &sp_data.exps[cofs_b],
                                                            &sp_data.coords[6 * ipair],
                                                            &sp_data.coords[6 * ipair + 3]);

        ecoeffs[ipair] = ecoeffs_ipair;
    }

    return ecoeffs;
}

vector<vector<lible::vec4d>>
LI::ecoeffsSPData_Eijt(const int la, const int lb, const ShellPairData &sp_data)
{
    vector<vector<vec4d>> ecoeffs(sp_data.n_pairs);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        vector<vec4d> ecoeffs_ipair = ecoeffsShellPair_Eijt(la, lb, cdepth_a, cdepth_b,
                                                            &sp_data.exps[cofs_a],
                                                            &sp_data.exps[cofs_b],
                                                            &sp_data.coords[6 * ipair],
                                                            &sp_data.coords[6 * ipair + 3]);

        ecoeffs[ipair] = ecoeffs_ipair;
    }

    return ecoeffs;
}

vector<vector<double>> LI::ecoeffsShell(const int l, const vector<double> &exps)
{
    int dim = numCartesians(l);
    size_t k = exps.size();

    vector<vector<double>> ecoeffs_out(k * k, vector<double>(dim * dim, 0));

    vec3d ecoeffs_x(l + 1, l + 1, 2 * l + 1, 0);
    vec3d ecoeffs_y(l + 1, l + 1, 2 * l + 1, 0);
    vec3d ecoeffs_z(l + 1, l + 1, 2 * l + 1, 0);

    array<double, 3> xyz_a{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            std::array<double, 3> Kab{1, 1, 1};            

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, l, l, &xyz_a[0],
                                                                          &xyz_a[0], &Kab[0]);

            for (const auto [i, j, k, mu] : cart_exps_a)
                for (const auto [i_, j_, k_, nu] : cart_exps_a)
                {
                    int munu = mu * dim + nu;
                    ecoeffs_out[iab][munu] = ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) *
                                             ecoeffs_z(k, k_, 0);
                }
        }

    return ecoeffs_out;
}

vector<double>
LI::ecoeffsSphericalSPData_Bra(const int la, const int lb, const ShellPairData &sp_data)
{
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

    int n_ecoeffs_sph = numSphericals(la) * numSphericals(lb) * numHermites(lab) *
                        sp_data.n_prim_pairs;

    vector<double> ecoeffs(n_ecoeffs_sph, 0);    
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

                auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] =
                    ecoeffsPrimitivePair(a, b, la, lb, &xyz_a[0], &xyz_b[0], &Kab[0]);

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
                                        ecoeffs_x(i, i_, t) * ecoeffs_y(j, j_, u) * ecoeffs_z(k, k_, v);
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
                                ecoeffs[idx] += ecoeff;
                            }
            }
    }

    return ecoeffs;
}

pair<vector<double>, vector<double>>
LI::ecoeffsSphericalSPData_BraKet(const int la, const int lb, const ShellPairData &sp_data)
{
    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    int n_ecoeffs_sph = numSphericals(la) * numSphericals(lb) * numHermites(lab) *
                        sp_data.n_prim_pairs;

    vector<double> ecoeffs = ecoeffsSphericalSPData_Bra(la, lb, sp_data);

    vector<double> ecoeffs_tsp(n_ecoeffs_sph, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];
        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos = offset_ecoeffs + iab * n_ecoeffs;
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int tuv = 0; tuv < dim_tuv; tuv++)
                        {
                            int idx = pos + munu * dim_tuv + tuv;
                            double ecoeff = ecoeffs[idx];

                            int idx_tsp = pos + tuv * dim_ab + munu;
                            ecoeffs_tsp[idx_tsp] = ecoeff;
                        }
            }
    }

    return {ecoeffs, ecoeffs_tsp};
}

vector<double>
LI::ecoeffsSphericalShellData_Bra(const int l, const ShellData &sh_data)
{
    arma::dmat sph_trafo = returnSphericalTrafo(l);

    const auto &cart_exps = cartExps(l);

    int dim_cart = numCartesians(l);
    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    vec3i tuv_poss = returnHermiteGaussianPositions(l);
    vector<array<int, 3>> tuv_idxs = returnHermiteGaussianIdxs(l);

    int n_ecoeffs = numSphericals(l) * numHermites(l) * sh_data.n_primitives;

    vector<double> ecoeffs(n_ecoeffs, 0);
    vector<double> ecoeffs_tsp(n_ecoeffs, 0);
    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int dim = sh_data.cdepths[ishell];
        int pos = sh_data.coffsets[ishell];

        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];

        vec2d ecoeffs_c(dim_cart, dim_tuv, 0);
        for (int i = 0; i < dim; i++)
        {
            double a = sh_data.exps[pos + i];
            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitive(a, l);

            for (size_t mu = 0; mu < cart_exps.size(); mu++)
            {
                auto [i, j, k] = cart_exps[mu];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            int tuv = tuv_poss(t, u, v);
                            ecoeffs_c(mu, tuv) = ecoeffs_x(i, t) * ecoeffs_y(j, u) *
                                                 ecoeffs_z(k, v);
                        }
            }

            double d = sh_data.coeffs[pos + i];
            int pos_ecoeffs = offset_ecoeffs + i * dim_sph * dim_tuv;
            for (int mu = 0; mu < dim_sph; mu++)
                for (int mu_ = 0; mu_ < dim_cart; mu_++)
                    for (int tuv = 0; tuv < dim_tuv; tuv++)
                    {
                        double ecoeff = d * ecoeffs_c(mu_, tuv) * sph_trafo(mu, mu_);

                        int idx = pos_ecoeffs + mu * dim_tuv + tuv;
                        ecoeffs[idx] += ecoeff;

                        int idx_tsp = pos_ecoeffs + tuv * dim_sph + mu;
                        ecoeffs_tsp[idx_tsp] += ecoeff;
                    }
        }
    }

    return ecoeffs;
}

pair<vector<double>, vector<double>>
LI::ecoeffsSphericalShellData_BraKet(const int l, const ShellData &sh_data)
{
    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    int n_ecoeffs = numSphericals(l) * numHermites(l) * sh_data.n_primitives;

    vector<double> ecoeffs = ecoeffsSphericalShellData_Bra(l, sh_data);

    vector<double> ecoeffs_tsp(n_ecoeffs, 0);
    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int dim = sh_data.cdepths[ishell];
        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];
        for (int i = 0; i < dim; i++)
        {
            int pos_ecoeffs = offset_ecoeffs + i * dim_sph * dim_tuv;
            for (int mu = 0; mu < dim_sph; mu++)
                for (int tuv = 0; tuv < dim_tuv; tuv++)
                {
                    int idx = pos_ecoeffs + mu * dim_tuv + tuv;
                    double ecoeff = ecoeffs[idx];

                    int idx_tsp = pos_ecoeffs + tuv * dim_sph + mu;
                    ecoeffs_tsp[idx_tsp] = ecoeff;
                }
        }
    }

    return {ecoeffs, ecoeffs_tsp};
}

vector<vector<double>>
LI::ecoeffsSphericalSPDatas_Bra(const vector<pair<int, int>> &l_pairs,
                                const vector<ShellPairData> &sp_datas)
{
    if (l_pairs.size() != sp_datas.size())
        throw std::runtime_error("The sizes of sp_datas and l_pairs don't match!\n");

    vector<vector<double>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        auto ecoeffs_ipair = ecoeffsSphericalSPData_Bra(la, lb, sp_datas[ipair]);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }

    return ecoeffs;
}

pair<vector<vector<double>>, vector<vector<double>>>
LI::ecoeffsSphericalSPDatas_BraKet(const vector<pair<int, int>> &l_pairs,
                                   const vector<ShellPairData> &sp_datas)
{
    if (l_pairs.size() != sp_datas.size())
        throw std::runtime_error("The sizes of sp_datas and l_pairs don't match!\n");

    vector<vector<double>> ecoeffs(l_pairs.size());
    vector<vector<double>> ecoeffs_tsp(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        auto [ecoeffs_ipair, ecoeffs_tsp_ipair] = ecoeffsSphericalSPData_BraKet(la, lb, sp_datas[ipair]);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_tsp[ipair] = std::move(ecoeffs_tsp_ipair);
    }

    return {ecoeffs, ecoeffs_tsp};
}

vector<vector<double>>
LI::ecoeffsSphericalShellDatas_Bra(const int l_max_aux, const vector<ShellData> &sh_datas)
{
    if (size_t(l_max_aux + 1) != (size_t)sh_datas.size())
        throw std::runtime_error("The size of sh_datas doesn't match (l_max_aux + 1)\n");

    vector<vector<double>> ecoeffs(sh_datas.size());
    for (int l = 0; l <= l_max_aux; l++)
    {
        auto ecoeffs_l = ecoeffsSphericalShellData_Bra(l, sh_datas[l]);

        ecoeffs[l] = std::move(ecoeffs_l);
    }

    return ecoeffs;
}

pair<vector<vector<double>>, vector<vector<double>>>
LI::ecoeffsSphericalShellDatas_BraKet(const int l_max_aux, const vector<ShellData> &sh_datas)
{
    if (size_t(l_max_aux + 1) != (size_t)sh_datas.size())
        throw std::runtime_error("The size of sh_datas doesn't match (l_max_aux + 1)\n");

    vector<vector<double>> ecoeffs(sh_datas.size());
    vector<vector<double>> ecoeffs_tsp(sh_datas.size());
    for (int l = 0; l <= l_max_aux; l++)
    {
        auto [ecoeffs_l, ecoeffs_tsp_l] = ecoeffsSphericalShellData_BraKet(l, sh_datas[l]);

        ecoeffs[l] = std::move(ecoeffs_l);
        ecoeffs_tsp[l] = std::move(ecoeffs_tsp_l);
    }

    return {ecoeffs, ecoeffs_tsp};
}