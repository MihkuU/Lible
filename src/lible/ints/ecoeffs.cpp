#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <array>
#include <stdexcept>
#include <tuple>

namespace LI = lible::ints;

using std::array, std::pair, std::tuple, std::vector;

lible::vec3d LI::ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                    const double PA, const double PB, const double Kab)
{
    vec3d ecoeffs(Fill(0), la + 1, lb + 1, la + lb + 1);

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
    vec3d ecoeffs1(Fill(0), la + 1, lb + 1, la + lb + 1);

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
            ecoeffs1(i, j, 0) = (a / p) * (R * ecoeffs1(i, j - 1, 0) + ecoeffs(i, j - 1, 0)) +
                                ecoeffs1(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                ecoeffs1(i, j, t) = one_o_2p * ecoeffs1(i, j - 1, t - 1) +
                                    (a / p) * (R * ecoeffs1(i, j - 1, t) + ecoeffs(i, j - 1, t)) +
                                    (t + 1) * ecoeffs1(i, j - 1, t + 1);

            ecoeffs1(i, j, i + j) = one_o_2p * ecoeffs1(i, j - 1, i + j - 1) +
                                    (a / p) * (R * ecoeffs1(i, j - 1, i + j) + ecoeffs(i, j - 1, i + j));
        }

    return ecoeffs1;
}

lible::vec2d LI::ecoeffsRecurrence1(const double one_o_2a, const int l)
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

array<lible::vec3d, 3> LI::ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                const int lb, const double *xyz_a,
                                                const double *xyz_b)
{
    double p = a + b;
    double mu = a * b / p;

    array<double, 3> xyz_p;
    array<double, 3> Kab;
    for (int i = 0; i < 3; i++)
    {
        xyz_p[i] = (a * xyz_a[i] + b * xyz_b[i]) / p;
        Kab[i] = std::exp(-mu * std::pow(xyz_a[i] - xyz_b[i], 2));
    }

    array<double, 3> xyz_pa;
    array<double, 3> xyz_pb;
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

array<lible::vec3d, 3> LI::ecoeffsPrimitivePair_n1(const double a, const double b, const int la,
                                                   const int lb, const double *xyz_a,
                                                   const double *xyz_b,
                                                   const array<lible::vec3d, 3> &ecoeffs)
{
    vec3d ecoeffs1_x = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[0], xyz_b[0], ecoeffs[0]);
    vec3d ecoeffs1_y = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[1], xyz_b[1], ecoeffs[1]);
    vec3d ecoeffs1_z = ecoeffsRecurrence2_n1(a, b, la, lb, xyz_a[2], xyz_b[2], ecoeffs[2]);

    return {ecoeffs1_x, ecoeffs1_y, ecoeffs1_z};
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
    vector<vec3d> ecoeffs(cdepth_a * cdepth_b, vec3d(Fill(0), 3, la + 1, lb + 1));
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb,
                                                                          xyz_a, xyz_b);

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
    vector<vec4d> ecoeffs(cdepth_a * cdepth_b);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb,
                                                                          xyz_a, xyz_b);

            vec4d ecoeffs_xyz(3, la + 1, lb + 1, la + lb + 1);
            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(0, i, j, t) = ecoeffs_x(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(1, i, j, t) = ecoeffs_y(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(2, i, j, t) = ecoeffs_z(i, j, t);

            ecoeffs[iab] = ecoeffs_xyz;
        }

    return ecoeffs;
}

vector<lible::vec4d>
LI::ecoeffsShellPair_Eijt_Debug(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                                const double *exps_a, const double *exps_b, const double *xyz_a,
                                const double *xyz_b)
{
    vector<vec4d> ecoeffs(cdepth_a * cdepth_b);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb,
                                                                          xyz_a, xyz_b);

            vec4d ecoeffs_xyz(3, la + 1, lb + 1, la + lb + 1);
            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(0, i, j, t) = ecoeffs_x(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(1, i, j, t) = ecoeffs_y(i, j, t);

            for (int i = 0; i <= la; i++)
                for (int j = 0; j <= lb; j++)
                    for (int t = 0; t <= i + j; t++)
                        ecoeffs_xyz(2, i, j, t) = ecoeffs_z(i, j, t);

            ecoeffs[iab] = ecoeffs_xyz;
        }

    return ecoeffs;
}

vector<vector<lible::vec3d>>
LI::ecoeffsSPData_Eij0(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.la;

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
LI::ecoeffsSPData_Eijt(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.la;

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

    vec3d ecoeffs_x(Fill(0), l + 1, l + 1, 2 * l + 1);
    vec3d ecoeffs_y(Fill(0), l + 1, l + 1, 2 * l + 1);
    vec3d ecoeffs_z(Fill(0), l + 1, l + 1, 2 * l + 1);

    array<double, 3> xyz_a{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, l, l, &xyz_a[0],
                                                                          &xyz_a[0]);

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

vector<double> LI::ecoeffsSHARKShellPair(const int ipair, const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;

    int n_cart_b = numCartesians(lb);
    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_hermite_ab = numHermites(lab);

    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cdepth_ab = cdepth_a * cdepth_b;
    int n_cols = cdepth_ab * n_hermite_ab;

    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    vector<array<int, 3>> cart_exps_a = cartExps(la);
    vector<array<int, 3>> cart_exps_b = cartExps(lb);
    vec3i hermite_poss = getHermiteGaussianPositions(lab);

    int n_ecoeffs = n_sph_ab * (cdepth_ab * n_hermite_ab);
    vector<double> ecoeffs(n_ecoeffs, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = sp_data.exps[cofs_a + ia];
            double b = sp_data.exps[cofs_b + ib];

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            vec3d ecoeffs_ppair_sc(Fill(0), n_sph_a, n_cart_b, n_hermite_ab);
            for (auto &[mu, mu_, val] : sph_trafo_a)
                for (size_t nu_ = 0; nu_ < cart_exps_b.size(); nu_++)
                {
                    auto [i, j, k] = cart_exps_a[mu_];
                    auto [i_, j_, k_] = cart_exps_b[nu_];
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                int tuv = hermite_poss(t, u, v);

                                ecoeffs_ppair_sc(mu, nu_, tuv) +=
                                    val * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                            }
                }

            double da = sp_data.coeffs[cofs_a + ia];
            double db = sp_data.coeffs[cofs_b + ib];
            double dadb = da * db;

            for (auto &[nu, nu_, val] : sph_trafo_b)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                    {
                        int irow = mu * n_sph_b + nu;
                        int icol = iab * n_hermite_ab + tuv;
                        int idx = irow * n_cols + icol;
                        ecoeffs[idx] += dadb * val * ecoeffs_ppair_sc(mu, nu_, tuv);
                    }
        }

    return ecoeffs;
}

vector<double> LI::ecoeffsSHARKShell(const int ishell, const ShellData &sh_data)
{
    int l = sh_data.l;

    vector<tuple<int, int, double>> sph_trafo = sphericalTrafo(l);

    const auto &cart_exps = cartExps(l);    

    vec3i tuv_poss = getHermiteGaussianPositions(l);
    vector<array<int, 3>> tuv_idxs = getHermiteGaussianIdxs(l);

    int n_sph = numSphericals(l);
    int n_hermite = numHermites(l);
    int cdepth = sh_data.cdepths[ishell];
    int cofs = sh_data.coffsets[ishell];
    int n_cols = cdepth * n_hermite;

    int n_ecoeffs = n_sph * (cdepth * n_hermite);
    vector<double> ecoeffs(n_ecoeffs, 0);
    for (int ia = 0; ia < cdepth; ia++)
    {        
        double a = sh_data.exps[cofs + ia];
        double d = sh_data.coeffs[cofs + ia];

        auto [Ex, Ey, Ez] = ecoeffsPrimitive(a, l);

        for (auto &[mu, mu_, val] : sph_trafo)
        {
            auto [i, j, k] = cart_exps[mu_];
            for (int t = 0; t <= i; t++)
                for (int u = 0; u <= j; u++)
                    for (int v = 0; v <= k; v++)
                    {
                        double ecoeff = d * Ex(i, t) * Ey(j, u) * Ez(k, v) * val;

                        int tuv = tuv_poss(t, u, v);

                        int irow = mu;
                        int icol = ia * n_hermite + tuv;
                        int idx = irow * n_cols + icol;
                        ecoeffs[idx] += ecoeff;
                    }
        }
    }

    return ecoeffs;
}

vector<double> LI::ecoeffsD1SHARKShellPair(const int ipair, const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;

    int n_cart_b = numCartesians(lb);
    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_hermite_ab = numHermites(lab);

    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cdepth_ab = cdepth_a * cdepth_b;
    int n_cols = cdepth_ab * n_hermite_ab;

    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);
    vector<array<int, 3>> cart_exps_a = cartExps(la);
    vector<array<int, 3>> cart_exps_b = cartExps(lb);
    vec3i hermite_poss = getHermiteGaussianPositions(lab);

    int n_ecoeffs0 = n_sph_ab * (cdepth_ab * n_hermite_ab);
    int n_ecoeffsd1 = 3 * n_ecoeffs0;
    int ofs0 = 0 * n_ecoeffs0;
    int ofs1 = 1 * n_ecoeffs0;
    int ofs2 = 2 * n_ecoeffs0;    
    vector<double> ecoeffsd1(n_ecoeffsd1, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = sp_data.exps[cofs_a + ia];
            double b = sp_data.exps[cofs_b + ib];

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b, {Ex, Ey, Ez});

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
                                int tuv = hermite_poss(t, u, v);

                                ecoeffs100_ppair_sc(mu, nu_, tuv) +=
                                    val * E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                ecoeffs010_ppair_sc(mu, nu_, tuv) +=
                                    val * Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                                ecoeffs001_ppair_sc(mu, nu_, tuv) +=
                                    val * Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);
                            }
                }

            double da = sp_data.coeffs[cofs_a + ia];
            double db = sp_data.coeffs[cofs_b + ib];
            double dadb = da * db;

            for (auto &[nu, nu_, val] : sph_trafo_b)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                    {
                        int irow = mu * n_sph_b + nu;
                        int icol = iab * n_hermite_ab + tuv;
                        int idx = irow * n_cols + icol;
                        ecoeffsd1[ofs0 + idx] += val * dadb * ecoeffs100_ppair_sc(mu, nu_, tuv);
                        ecoeffsd1[ofs1 + idx] += val * dadb * ecoeffs010_ppair_sc(mu, nu_, tuv);
                        ecoeffsd1[ofs2 + idx] += val * dadb * ecoeffs001_ppair_sc(mu, nu_, tuv);
                    }            
        }

    return ecoeffsd1;
}

vector<double> LI::ecoeffsSHARK(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;
    int n_sph_ab = numSphericals(la) * numSphericals(lb);
    int n_hermite_ab = numHermites(lab);
    int n_ecoeffs_pp = n_sph_ab * n_hermite_ab;
    int n_ecoeffs = sp_data.n_prim_pairs * n_ecoeffs_pp;

    vector<double> ecoeffs(n_ecoeffs, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        vector<double> ecoeffs_ipair = ecoeffsSHARKShellPair(ipair, sp_data);

        int ofs = sp_data.offsets_ecoeffs[ipair];
        for (size_t i = 0; i < ecoeffs_ipair.size(); i++)
            ecoeffs[ofs + i] = ecoeffs_ipair[i];
    }

    return ecoeffs;
}

vector<double> LI::ecoeffsD1SHARK(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int lab = la + lb;
    int n_sph_ab = numSphericals(la) * numSphericals(lb);
    int n_hermite_ab = numHermites(lab);
    int n_ecoeffsd1_pp = 3 * n_sph_ab * n_hermite_ab;
    int n_ecoeffsd1 = sp_data.n_prim_pairs * n_ecoeffsd1_pp;

    vector<double> ecoeffsd1(n_ecoeffsd1, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        vector<double> ecoeffsd1_ipair = ecoeffsD1SHARKShellPair(ipair, sp_data);

        int ofs = sp_data.offsets_ecoeffs_deriv1[ipair];
        for (size_t i = 0; i < ecoeffsd1_ipair.size(); i++)
            ecoeffsd1[ofs + i] = ecoeffsd1_ipair[i];
    }

    return ecoeffsd1;
}

vector<double> LI::ecoeffsSHARK(const ShellData &sh_data)
{
    int l = sh_data.l;
    int n_sph = numSphericals(l);
    int n_hermite = numHermites(l);
    int n_ecoeffs_pp = n_sph * n_hermite;
    int n_ecoeffs = sh_data.n_primitives * n_ecoeffs_pp;

    vector<double> ecoeffs(n_ecoeffs);
    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        vector<double> ecoeffs_ishell = ecoeffsSHARKShell(ishell, sh_data);

        int ofs = sh_data.offsets_ecoeffs[ishell];
        for (size_t i = 0; i < ecoeffs_ishell.size(); i++)
            ecoeffs[ofs + i] = ecoeffs_ishell[i];
    }

    return ecoeffs;
}

vector<double>
LI::ecoeffsSphericalSPData_Bra(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;

    vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);

    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    int dim_b_cart = numCartesians(lb);
    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vec3i tuv_poss = getHermiteGaussianPositions(lab);

    int n_ecoeffs_sph = numSphericals(la) * numSphericals(lb) * numHermites(lab) *
                        sp_data.n_prim_pairs;

    vector<double> ecoeffs(n_ecoeffs_sph, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        const double *xyz_a = &sp_data.coords[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords[6 * ipair + 3];

        int offset_ecoeffs = sp_data.offsets_ecoeffs[ipair];

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                double a = sp_data.exps[cofs_a + ia];
                double b = sp_data.exps[cofs_b + ib];

                auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

                vec3d ecoeffs_ppair_sc(Fill(0), dim_a_sph, dim_b_cart, dim_tuv);
                for (auto &[a, a_, val] : sph_trafo_a)
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
                for (auto &[b, b_, val] : sph_trafo_b)
                    for (int a = 0; a < dim_a_sph; a++)
                        for (int tuv = 0; tuv < dim_tuv; tuv++)
                        {
                            int ab = a * dim_b_sph + b;
                            int idx = ofs + ab * dim_tuv + tuv;
                            ecoeffs[idx] += dadb * val * ecoeffs_ppair_sc(a, b_, tuv);
                        }
            }
    }

    return ecoeffs;
}

vector<double>
LI::ecoeffsSphericalSPData_Bra_Deriv1(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;

    vector<tuple<int, int, double>> sph_trafo_a = sphericalTrafo(la);
    vector<tuple<int, int, double>> sph_trafo_b = sphericalTrafo(lb);

    const auto &cart_exps_a = cartExps(la);
    const auto &cart_exps_b = cartExps(lb);

    int dim_b_cart = numCartesians(lb);
    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vec3i tuv_poss = getHermiteGaussianPositions(lab);

    int n_ecoeffs_prims = 3 * n_ecoeffs * sp_data.n_prim_pairs;

    vector<double> ecoeffs_100_010_001(n_ecoeffs_prims, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        const double *xyz_a = &sp_data.coords[6 * ipair + 0];
        const double *xyz_b = &sp_data.coords[6 * ipair + 3];

        int offset_ecoeffs = sp_data.offsets_ecoeffs_deriv1[ipair];

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                double a = sp_data.exps[cofs_a + ia];
                double b = sp_data.exps[cofs_b + ib];

                auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

                auto [E1x, E1y, E1z] =
                    ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b, {Ex, Ey, Ez});

                vec3d ecoeffs100_ppair_sc(Fill(0), dim_a_sph, dim_b_cart, dim_tuv);
                vec3d ecoeffs010_ppair_sc(Fill(0), dim_a_sph, dim_b_cart, dim_tuv);
                vec3d ecoeffs001_ppair_sc(Fill(0), dim_a_sph, dim_b_cart, dim_tuv);
                for (auto &[a, a_, val] : sph_trafo_a)
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
                for (auto &[b, b_, val] : sph_trafo_b)
                    for (int a = 0; a < dim_a_sph; a++)
                        for (int tuv = 0; tuv < dim_tuv; tuv++)
                        {
                            int ab = a * dim_b_sph + b;
                            int idx_100 = ofs_100 + ab * dim_tuv + tuv;
                            int idx_010 = ofs_010 + ab * dim_tuv + tuv;
                            int idx_001 = ofs_001 + ab * dim_tuv + tuv;
                            ecoeffs_100_010_001[idx_100] += val * dadb * ecoeffs100_ppair_sc(a, b_, tuv);
                            ecoeffs_100_010_001[idx_010] += val * dadb * ecoeffs010_ppair_sc(a, b_, tuv);
                            ecoeffs_100_010_001[idx_001] += val * dadb * ecoeffs001_ppair_sc(a, b_, tuv);
                        }
            }
    }

    return ecoeffs_100_010_001;
}

pair<vector<double>, vector<double>>
LI::ecoeffsSphericalSPData_BraKet(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;

    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    int n_ecoeffs_sph = numSphericals(la) * numSphericals(lb) * numHermites(lab) *
                        sp_data.n_prim_pairs;

    vector<double> ecoeffs = ecoeffsSphericalSPData_Bra(sp_data);

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

pair<vector<double>, vector<double>>
LI::ecoeffsSphericalSPData_BraKet_Deriv1(const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;

    int dim_a_sph = numSphericals(la);
    int dim_b_sph = numSphericals(lb);
    int dim_ab = dim_a_sph * dim_b_sph;

    int lab = la + lb;
    int dim_tuv = numHermites(lab);
    int n_ecoeffs = dim_ab * dim_tuv;

    vector<double> ecoeffs = ecoeffsSphericalSPData_Bra_Deriv1(sp_data);

    int n_ecoeffs_sph_x1 = dim_ab * dim_tuv * sp_data.n_prim_pairs;
    vector<double> ecoeffs_tsp(3 * n_ecoeffs_sph_x1, 0);
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int offset_ecoeffs = sp_data.offsets_ecoeffs_deriv1[ipair];
        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                int ofs_100 = offset_ecoeffs + (3 * iab + 0) * n_ecoeffs;
                int ofs_010 = offset_ecoeffs + (3 * iab + 1) * n_ecoeffs;
                int ofs_001 = offset_ecoeffs + (3 * iab + 2) * n_ecoeffs;
                for (int a = 0; a < dim_a_sph; a++)
                    for (int b = 0; b < dim_b_sph; b++)
                        for (int tuv = 0; tuv < dim_tuv; tuv++)
                        {
                            int ab = a * dim_b_sph + b;
                            int idx_100 = ofs_100 + ab * dim_tuv + tuv;
                            int idx_010 = ofs_010 + ab * dim_tuv + tuv;
                            int idx_001 = ofs_001 + ab * dim_tuv + tuv;
                            double ecoeff_100 = ecoeffs[idx_100];
                            double ecoeff_010 = ecoeffs[idx_010];
                            double ecoeff_001 = ecoeffs[idx_001];

                            int idx_100_tsp = ofs_100 + tuv * dim_ab + ab;
                            int idx_010_tsp = ofs_010 + tuv * dim_ab + ab;
                            int idx_001_tsp = ofs_001 + tuv * dim_ab + ab;
                            ecoeffs_tsp[idx_100_tsp] = ecoeff_100;
                            ecoeffs_tsp[idx_010_tsp] = ecoeff_010;
                            ecoeffs_tsp[idx_001_tsp] = ecoeff_001;
                        }
            }
    }

    return {ecoeffs, ecoeffs_tsp};
}

vector<double>
LI::ecoeffsSphericalShellData_Bra(const ShellData &sh_data)
{
    int l = sh_data.l;

    vector<tuple<int, int, double>> sph_trafo = sphericalTrafo(l);

    const auto &cart_exps = cartExps(l);

    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    vec3i tuv_poss = getHermiteGaussianPositions(l);
    vector<array<int, 3>> tuv_idxs = getHermiteGaussianIdxs(l);

    int n_ecoeffs = numSphericals(l) * numHermites(l) * sh_data.n_primitives;

    vector<double> ecoeffs(n_ecoeffs, 0);
    for (int ishell = 0; ishell < sh_data.n_shells; ishell++)
    {
        int cdepth = sh_data.cdepths[ishell];
        int cofs = sh_data.coffsets[ishell];

        int offset_ecoeffs = sh_data.offsets_ecoeffs[ishell];

        for (int i = 0; i < cdepth; i++)
        {
            double a = sh_data.exps[cofs + i];

            auto [Ex, Ey, Ez] = ecoeffsPrimitive(a, l);

            double d = sh_data.coeffs[cofs + i];
            int ofs = offset_ecoeffs + i * dim_sph * dim_tuv;
            for (auto &[a, a_, val] : sph_trafo)
            {
                auto [i, j, k] = cart_exps[a_];
                for (int t = 0; t <= i; t++)
                    for (int u = 0; u <= j; u++)
                        for (int v = 0; v <= k; v++)
                        {
                            double ecoeff = d * Ex(i, t) * Ey(j, u) * Ez(k, v) * val;

                            int tuv = tuv_poss(t, u, v);

                            int idx = ofs + a * dim_tuv + tuv;
                            ecoeffs[idx] += ecoeff;
                        }
            }
        }
    }

    return ecoeffs;
}

pair<vector<double>, vector<double>>
LI::ecoeffsSphericalShellData_BraKet(const ShellData &sh_data)
{
    int l = sh_data.l;

    int dim_sph = numSphericals(l);
    int dim_tuv = numHermites(l);

    int n_ecoeffs = numSphericals(l) * numHermites(l) * sh_data.n_primitives;

    vector<double> ecoeffs = ecoeffsSphericalShellData_Bra(sh_data);

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
LI::ecoeffsSphericalSPDatas_Bra(const vector<ShellPairData> &sp_datas)
{
    vector<vector<double>> ecoeffs(sp_datas.size());
    for (size_t idata = 0; idata < sp_datas.size(); idata++)
        ecoeffs[idata] = ecoeffsSphericalSPData_Bra(sp_datas[idata]);

    return ecoeffs;
}

vector<vector<double>>
LI::ecoeffsSphericalSPDatas_Bra_Deriv1(const vector<ShellPairData> &sp_datas)
{
    vector<vector<double>> ecoeffs(sp_datas.size());
    for (size_t idata = 0; idata < sp_datas.size(); idata++)
        ecoeffs[idata] = ecoeffsSphericalSPData_Bra_Deriv1(sp_datas[idata]);

    return ecoeffs;
}

pair<vector<vector<double>>, vector<vector<double>>>
LI::ecoeffsSphericalSPDatas_BraKet(const vector<ShellPairData> &sp_datas)
{
    vector<vector<double>> ecoeffs(sp_datas.size());
    vector<vector<double>> ecoeffs_tsp(sp_datas.size());
    for (size_t idata = 0; idata < sp_datas.size(); idata++)
    {
        auto [ecoeffs_idata, ecoeffs_tsp_idata] = ecoeffsSphericalSPData_BraKet(sp_datas[idata]);

        ecoeffs[idata] = std::move(ecoeffs_idata);
        ecoeffs_tsp[idata] = std::move(ecoeffs_tsp_idata);
    }

    return {ecoeffs, ecoeffs_tsp};
}

pair<vector<vector<double>>, vector<vector<double>>>
LI::ecoeffsSphericalSPDatas_BraKet_Deriv1(const vector<ShellPairData> &sp_datas)
{
    vector<vector<double>> ecoeffs(sp_datas.size());
    vector<vector<double>> ecoeffs_tsp(sp_datas.size());
    for (size_t idata = 0; idata < sp_datas.size(); idata++)
    {
        auto [ecoeffs_idata, ecoeffs_tsp_idata] =
            ecoeffsSphericalSPData_BraKet_Deriv1(sp_datas[idata]);

        ecoeffs[idata] = std::move(ecoeffs_idata);
        ecoeffs_tsp[idata] = std::move(ecoeffs_tsp_idata);
    }

    return {ecoeffs, ecoeffs_tsp};
}

vector<vector<double>>
LI::ecoeffsSphericalShellDatas_Bra(const vector<ShellData> &sh_datas)
{
    vector<vector<double>> ecoeffs(sh_datas.size());
    for (size_t l = 0; l < sh_datas.size(); l++)
        ecoeffs[l] = ecoeffsSphericalShellData_Bra(sh_datas[l]);

    return ecoeffs;
}

pair<vector<vector<double>>, vector<vector<double>>>
LI::ecoeffsSphericalShellDatas_BraKet(const vector<ShellData> &sh_datas)
{
    vector<vector<double>> ecoeffs(sh_datas.size());
    vector<vector<double>> ecoeffs_tsp(sh_datas.size());
    for (size_t l = 0; l < sh_datas.size(); l++)
    {
        auto [ecoeffs_l, ecoeffs_tsp_l] = ecoeffsSphericalShellData_BraKet(sh_datas[l]);

        ecoeffs[l] = std::move(ecoeffs_l);
        ecoeffs_tsp[l] = std::move(ecoeffs_tsp_l);
    }

    return {ecoeffs, ecoeffs_tsp};
}