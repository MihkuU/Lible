#include <lible/ints/mcmurchie_davidson.hpp>
#include <lible/ints/ints_util.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/spherical_trafo.hpp>

namespace LIMD = lible::ints::MD;

using std::array, std::vector;

lible::vec3i LIMD::returnTUVPoss(const int l)
{
    vec3i tuv_poss(l + 1, l + 1, l + 1, -1);
    for (int t = 0, tuv = 0; t <= l; t++)
        for (int u = 0; u <= l - t; u++)
            for (int v = 0; v <= l - t - u; v++, tuv++)                            
                tuv_poss(t, u, v) = tuv;            

    return tuv_poss;
}

vector<LIMD::IdxsTUV> LIMD::returnIdxsTUV(const int l)
{
    vector<IdxsTUV> idxs_tuv((l + 1) * (l + 2) * (l + 3) / 6);
    for (int t = 0, tuv = 0; t <= l; t++)
        for (int u = 0; u <= l - t; u++)
            for (int v = 0; v <= l - t - u; v++, tuv++)
                idxs_tuv[tuv] = {t, u, v};

    return idxs_tuv;
}

void LIMD::calcECoeffs(const int l, const vector<double> &exps,
                       vector<arma::dmat> &ecoeffs_out)
{
    // Calculates Hermite expansion coefficients for a diagonal shell pair located on the same
    // atom.

    int dim = dimCartesians(l);
    size_t k = exps.size();

    ecoeffs_out.resize(k * k, arma::zeros(dim, dim));

    lible::vec3d Ex(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d Ey(l + 1, l + 1, 2 * l + 1, 0);
    lible::vec3d Ez(l + 1, l + 1, 2 * l + 1, 0);

    array<double, 3> A{0, 0, 0}; // dummy center

    const auto &cart_exps_a = cart_exps[l];

    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];

            std::array<double, 3> Kab{1, 1, 1};

            coeffs(a, b, l, l, A, A, Kab, Ex, Ey, Ez);

            for (const auto [i, j, k, mu] : cart_exps_a)
                for (const auto [i_, j_, k_, nu] : cart_exps_a)
                    ecoeffs_out[iab](mu, nu) = Ex(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
        }
}

void LIMD::calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                       vector<vector<arma::dmat>> &ecoeffs_out)
{
    // Calculates Hermite expansion coefficients for all shell-pairs corresponding to the
    // shell-pair data.

    ecoeffs_out.resize(shell_pair_data.n_pairs);

    int dim_a = dimCartesians(la);
    int dim_b = dimCartesians(lb);

    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &exps = shell_pair_data.exps[ipair];

        array<double, 3> A = shell_pair_data.coords[ipair].first;
        array<double, 3> B = shell_pair_data.coords[ipair].second;

        size_t ka = exps.first.size();
        size_t kb = exps.second.size();

        vector<arma::dmat> ecoeffs_pair(ka * kb, arma::zeros(dim_a, dim_b));

        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps.first[ia];
                double b = exps.second[ib];

                double mu = a * b / (a + b);

                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};

                coeffs(a, b, la, lb, A, B, Kab, Ex, Ey, Ez);

                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        ecoeffs_pair[iab](mu, nu) = Ex(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
            }

        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LIMD::calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                       vector<vector<vec3d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(shell_pair_data.n_pairs);
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];

        const auto &[A, B] = shell_pair_data.coords[ipair];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        vector<vec3d> ecoeffs_pair(ka * kb, vec3d(3, la + 1, lb + 1, 0));
        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double mu = a * b / (a + b);
                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};

                coeffs(a, b, la, lb, A, B, Kab, Ex, Ey, Ez);

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

void LIMD::calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                       vector<vector<vec4d>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    ecoeffs_out.resize(shell_pair_data.n_pairs);
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];

        const auto &[A, B] = shell_pair_data.coords[ipair];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        vector<vec4d> ecoeffs_pair(ka * kb, vec4d(3, la + 1, lb + 1, la + lb + 1, 0));
        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double mu = a * b / (a + b);
                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};

                coeffs(a, b, la, lb, A, B, Kab, Ex, Ey, Ez);

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

void LIMD::calcECoeffsSpherical(const int la, const int lb,
                                const ShellPairData &shell_pair_data,
                                vector<vector<arma::dmat>> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    arma::dmat sph_trafo_a = returnSphericalTrafo(la);
    arma::dmat sph_trafo_b = returnSphericalTrafo(lb);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int dim_a_cart = dimCartesians(la);
    int dim_b_cart = dimCartesians(lb);
    int dim_a_sph = dimSphericals(la);
    int dim_b_sph = dimSphericals(lb);

    int lab = la + lb;
    int dim_tuv = (lab + 1) * (lab + 2) * (lab + 3) / 6;

    vec3i tuv_poss = returnTUVPoss(lab);
    vector<IdxsTUV> idxs_tuv = returnIdxsTUV(lab);

    ecoeffs_out.resize(shell_pair_data.n_pairs);
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];

        const auto &[A, B] = shell_pair_data.coords[ipair];

        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();        

        vec3d ecoeffs_ppair_cc(dim_a_cart, dim_b_cart, dim_tuv, 0);
        vec3d ecoeffs_ppair_sc(dim_a_sph, dim_b_cart, dim_tuv, 0);
        vector<arma::dmat> ecoeffs_pair(ka * kb, arma::zeros(dim_a_sph * dim_b_sph, dim_tuv));
        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double mu = a * b / (a + b);
                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};

                coeffs(a, b, la, lb, A, B, Kab, Ex, Ey, Ez);                

                ecoeffs_ppair_cc.set(0);                
                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);                              
                                    ecoeffs_ppair_cc(mu, nu, tuv) =
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }

                ecoeffs_ppair_sc.set(0);
                for (int mu = 0; mu < dim_a_sph; mu++)
                    for (int mu_ = 0; mu_ < dim_a_cart; mu_++)
                        for (int nu = 0; nu < dim_b_cart; nu++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_ppair_sc(mu, nu, tuv) +=
                                    ecoeffs_ppair_cc(mu_, nu, tuv) * sph_trafo_a(mu, mu_);

                double da = ccoeffs_a[ia], db = ccoeffs_b[ib];
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int nu_ = 0; nu_ < dim_b_cart; nu_++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_pair[iab](munu, tuv) += da * db *
                                                                ecoeffs_ppair_sc(mu, nu_, tuv) * 
                                                                sph_trafo_b(nu, nu_);
            }
        ecoeffs_out[ipair] = ecoeffs_pair;
    }
}

void LIMD::calcECoeffsSpherical(const int la, const int lb,
                                const ShellPairData &shell_pair_data,
                                vector<double> &ecoeffs_out)
{
    lible::vec3d Ex(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ey(la + 1, lb + 1, la + lb + 1, 0);
    lible::vec3d Ez(la + 1, lb + 1, la + lb + 1, 0);

    arma::dmat sph_trafo_a = returnSphericalTrafo(la);
    arma::dmat sph_trafo_b = returnSphericalTrafo(lb);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int dim_a_cart = dimCartesians(la);
    int dim_b_cart = dimCartesians(lb);
    int dim_a_sph = dimSphericals(la);
    int dim_b_sph = dimSphericals(lb);

    int lab = la + lb;
    int dim_tuv = (lab + 1) * (lab + 2) * (lab + 3) / 6;

    vec3i tuv_poss = returnTUVPoss(lab);
    vector<IdxsTUV> idxs_tuv = returnIdxsTUV(lab);

    size_t iprim = 0;
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];

        const auto &[A, B] = shell_pair_data.coords[ipair];

        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        vec3d ecoeffs_ppair_cc(dim_a_cart, dim_b_cart, dim_tuv, 0);
        vec3d ecoeffs_ppair_sc(dim_a_sph, dim_b_cart, dim_tuv, 0);
        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double mu = a * b / (a + b);
                std::array<double, 3> Kab{std::exp(-mu * std::pow(A[0] - B[0], 2)),
                                          std::exp(-mu * std::pow(A[1] - B[1], 2)),
                                          std::exp(-mu * std::pow(A[2] - B[2], 2))};

                coeffs(a, b, la, lb, A, B, Kab, Ex, Ey, Ez);

                ecoeffs_ppair_cc.set(0);
                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        for (int t = 0; t <= i + i_; t++)
                            for (int u = 0; u <= j + j_; u++)
                                for (int v = 0; v <= k + k_; v++)
                                {
                                    int tuv = tuv_poss(t, u, v);
                                    ecoeffs_ppair_cc(mu, nu, tuv) =
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                }

                ecoeffs_ppair_sc.set(0);
                for (int mu = 0; mu < dim_a_sph; mu++)
                    for (int mu_ = 0; mu_ < dim_a_cart; mu_++)
                        for (int nu = 0; nu < dim_b_cart; nu++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                                ecoeffs_ppair_sc(mu, nu, tuv) +=
                                    ecoeffs_ppair_cc(mu_, nu, tuv) * sph_trafo_a(mu, mu_);

                double da = ccoeffs_a[ia], db = ccoeffs_b[ib];
                size_t offset_ecoeffs = shell_pair_data.offsets_ecoeffs[iprim];
                for (int mu = 0, munu = 0; mu < dim_a_sph; mu++)
                    for (int nu = 0; nu < dim_b_sph; nu++, munu++)
                        for (int nu_ = 0; nu_ < dim_b_cart; nu_++)
                            for (int tuv = 0; tuv < dim_tuv; tuv++)
                            {
                                int idx = offset_ecoeffs + munu * dim_tuv + tuv;
                                ecoeffs_out[idx] += da * db *
                                                   ecoeffs_ppair_sc(mu, nu_, tuv) *
                                                   sph_trafo_b(nu, nu_);
                            }
                iprim++;                                
            }
    }
}

void LIMD::calcRInts(const int la, const int lb, const double p, const arma::vec::fixed<3> &RPC,
                     const vector<double> &fnx, vec4d &rints_tmp, vec3d &rints_out)
{
    rints_tmp.set(0);
    rints_out.set(0);

    rints_tmp(0, 0, 0, 0) = fnx[0];

    int lab = la + lb;
    double x = -2 * p;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_tmp(n, 0, 0, 0) = fnx[n] * y;
        y *= x;
    }

    // This clever loop is taken from HUMMR:
    for (int n = lab - 1; n >= 0; n--)
    {
        int n_ = lab - n;
        for (int t = 0; t <= n_; t++)
            for (int u = 0; u <= n_ - t; u++)
                for (int v = 0; v <= n_ - t - u; v++)
                {
                    if (t > 0)
                    {
                        rints_tmp(n, t, u, v) = RPC(0) * rints_tmp(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_tmp(n, t, u, v) += (t - 1) * rints_tmp(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(1) * rints_tmp(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_tmp(n, t, u, v) += (u - 1) * rints_tmp(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(2) * rints_tmp(n + 1, t, u, v - 1);
                            if (v > 1)
                                rints_tmp(n, t, u, v) += (v - 1) * rints_tmp(n + 1, t, u, v - 2);
                        }
                    }
                }
    }

    for (int t = 0; t <= lab; t++)
        for (int u = 0; u <= lab - t; u++)
            for (int v = 0; v <= lab - t - u; v++)
                rints_out(t, u, v) = rints_tmp(0, t, u, v);
}

void LIMD::calcRInts(const int la, const int lb, const double p,
                     const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                     const std::vector<IdxsTUV> &tuv_idxs_a,
                     const std::vector<IdxsTUV> &tuv_idxs_b,
                     vec4d &rints_tmp, arma::dmat &rints_out)
{
    rints_tmp.set(0);
    rints_out.zeros();
    
    rints_tmp(0, 0, 0, 0) = fnx[0];

    int lab = la + lb;
    double x = -2 * p;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_tmp(n, 0, 0, 0) = fnx[n] * y;
        y *= x;
    }

    // This clever loop is taken from HUMMR:
    for (int n = lab - 1; n >= 0; n--)
    {
        int n_ = lab - n;
        for (int t = 0; t <= n_; t++)
            for (int u = 0; u <= n_ - t; u++)
                for (int v = 0; v <= n_ - t - u; v++)
                {
                    if (t > 0)
                    {
                        rints_tmp(n, t, u, v) = RPC(0) * rints_tmp(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_tmp(n, t, u, v) += (t - 1) * rints_tmp(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(1) * rints_tmp(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_tmp(n, t, u, v) += (u - 1) * rints_tmp(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(2) * rints_tmp(n + 1, t, u, v - 1);
                            if (v > 1)
                                rints_tmp(n, t, u, v) += (v - 1) * rints_tmp(n + 1, t, u, v - 2);
                        }
                    }
                }
    }    

    for (size_t i = 0; i < tuv_idxs_a.size(); i++)
    {
        auto [t, u, v] = tuv_idxs_a[i];        
        for (size_t j = 0; j < tuv_idxs_b.size(); j++)
        {
            auto [t_, u_, v_] = tuv_idxs_b[j];
            
            double sign = 1.0; 
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;        

            rints_out(i, j) = sign * rints_tmp(0, t + t_, u + u_, v + v_);
        }
    }
}

void LIMD::calcRInts(const int la, const int lb, const double p,
                     const arma::vec::fixed<3> &RPC, const vector<double> &fnx,
                     const vector<IdxsTUV> &tuv_idxs_a,
                     const vector<IdxsTUV> &tuv_idxs_b,
                     vec4d &rints_tmp, vector<double> &rints_out)
{
    rints_tmp.set(0);    
    std::fill(rints_out.begin(), rints_out.end(), 0);

    rints_tmp(0, 0, 0, 0) = fnx[0];

    int lab = la + lb;
    double x = -2 * p;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_tmp(n, 0, 0, 0) = fnx[n] * y;
        y *= x;
    }

    // This clever loop is taken from HUMMR:
    for (int n = lab - 1; n >= 0; n--)
    {
        int n_ = lab - n;
        for (int t = 0; t <= n_; t++)
            for (int u = 0; u <= n_ - t; u++)
                for (int v = 0; v <= n_ - t - u; v++)
                {
                    if (t > 0)
                    {
                        rints_tmp(n, t, u, v) = RPC(0) * rints_tmp(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_tmp(n, t, u, v) += (t - 1) * rints_tmp(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(1) * rints_tmp(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_tmp(n, t, u, v) += (u - 1) * rints_tmp(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_tmp(n, t, u, v) = RPC(2) * rints_tmp(n + 1, t, u, v - 1);
                            if (v > 1)
                                rints_tmp(n, t, u, v) += (v - 1) * rints_tmp(n + 1, t, u, v - 2);
                        }
                    }
                }
    }

    int dim_tuv_b = (lb + 1) * (lb + 2) * (lb + 3) / 6;
    for (size_t i = 0; i < tuv_idxs_a.size(); i++)
    {
        auto [t, u, v] = tuv_idxs_a[i];
        for (size_t j = 0; j < tuv_idxs_b.size(); j++)
        {
            auto [t_, u_, v_] = tuv_idxs_b[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            rints_out[i * dim_tuv_b + j] = sign * rints_tmp(0, t + t_, u + u_, v + v_);
        }
    }
}

void LIMD::coeffs(const double a, const double b, const double PA, const double PB,
                  const double one_o_2p, const int la, const int lb, vec3d &E)
{
    for (int i = 1; i <= la; i++)
    {
        E(i, 0, 0) = PA * E(i - 1, 0, 0) + E(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            E(i, 0, t) = one_o_2p * E(i - 1, 0, t - 1) +
                         PA * E(i - 1, 0, t) +
                         (t + 1) * E(i - 1, 0, t + 1);

        E(i, 0, i) = one_o_2p * E(i - 1, 0, i - 1) + PA * E(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            E(i, j, 0) = PB * E(i, j - 1, 0) + E(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                E(i, j, t) = one_o_2p * E(i, j - 1, t - 1) +
                             PB * E(i, j - 1, t) +
                             (t + 1) * E(i, j - 1, t + 1);

            E(i, j, i + j) = one_o_2p * E(i, j - 1, i + j - 1) + PB * E(i, j - 1, i + j);
        }
}

void LIMD::coeffs(const double a, const double b, const int la, const int lb,
                  const array<double, 3> &A, const array<double, 3> &B,
                  const array<double, 3> &Kab, vec3d &Ex, vec3d &Ey, vec3d &Ez)
{
    Ex.set(0);
    Ey.set(0);
    Ez.set(0);

    Ex(0, 0, 0) = Kab[0];
    Ey(0, 0, 0) = Kab[1];
    Ez(0, 0, 0) = Kab[2];

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs(a, b, PA[0], PB[0], one_o_2p, la, lb, Ex);
    coeffs(a, b, PA[1], PB[1], one_o_2p, la, lb, Ey);
    coeffs(a, b, PA[2], PB[2], one_o_2p, la, lb, Ez);
}

// void LIMD::coeffs(const double one_o_2p, const int la, const int lb, 
//                   const arma::vec::fixed<3> &A, const arma::vec::fixed<3> &B, 
//                   const array<double, 3> &Kab, vec4d &E)
// {
//     E.set(0);

//     E(0, 0, 0, 0) = Kab[0];
//     E(1, 0, 0, 0) = Kab[1];
//     E(2, 0, 0, 0) = Kab[2];

//     for (int i = 1; i <= la; i++)
//     {
//         E(0, i, 0, 0) = A[0] * E(0, i - 1, 0, 0) + E(0, i - 1, 0, 1);
//         E(1, i, 0, 0) = A[1] * E(1, i - 1, 0, 0) + E(1, i - 1, 0, 1);
//         E(2, i, 0, 0) = A[2] * E(2, i - 1, 0, 0) + E(2, i - 1, 0, 1);

//         for (int t = 1; t < i; t++)
//         {
//             E(0, i, 0, t) = one_o_2p * E(0, i - 1, 0, t - 1) + A[0] * E(0, i - 1, 0, t) +
//                             (t + 1) * E(0, i - 1, 0, t + 1);
//             E(1, i, 0, t) = one_o_2p * E(1, i - 1, 0, t - 1) + A[1] * E(1, i - 1, 0, t) +
//                             (t + 1) * E(1, i - 1, 0, t + 1);
//             E(2, i, 0, t) = one_o_2p * E(2, i - 1, 0, t - 1) + A[2] * E(2, i - 1, 0, t) +
//                             (t + 1) * E(2, i - 1, 0, t + 1);
//         }

//         E(0, i, 0, i) = one_o_2p * E(0, i - 1, 0, i - 1) + A[0] * E(0, i - 1, 0, i);
//         E(1, i, 0, i) = one_o_2p * E(1, i - 1, 0, i - 1) + A[1] * E(1, i - 1, 0, i);
//         E(2, i, 0, i) = one_o_2p * E(2, i - 1, 0, i - 1) + A[2] * E(2, i - 1, 0, i);
//     }

//     for (int j = 1; j <= lb; j++)
//         for (int i = 0; i <= la; i++)
//         {
//             E(0, i, j, 0) = B[0] * E(0, i, j - 1, 0) + E(0, i, j - 1, 1);
//             E(1, i, j, 0) = B[1] * E(1, i, j - 1, 0) + E(1, i, j - 1, 1);
//             E(2, i, j, 0) = B[2] * E(2, i, j - 1, 0) + E(2, i, j - 1, 1);

//             for (int t = 1; t < i + j; t++)
//             {
//                 E(0, i, j, t) = one_o_2p * E(0, i, j - 1, t - 1) + B[0] * E(0, i, j - 1, t) +
//                                 (t + 1) * E(0, i, j - 1, t + 1);
//                 E(1, i, j, t) = one_o_2p * E(1, i, j - 1, t - 1) + B[1] * E(1, i, j - 1, t) +
//                                 (t + 1) * E(1, i, j - 1, t + 1);
//                 E(2, i, j, t) = one_o_2p * E(2, i, j - 1, t - 1) + B[2] * E(2, i, j - 1, t) +
//                                 (t + 1) * E(2, i, j - 1, t + 1);
//             }

//             E(0, i, j, i + j) = one_o_2p * E(0, i, j - 1, i + j - 1) + B[0] * E(0, i, j - 1, i + j);
//             E(1, i, j, i + j) = one_o_2p * E(1, i, j - 1, i + j - 1) + B[1] * E(1, i, j - 1, i + j);
//             E(2, i, j, i + j) = one_o_2p * E(2, i, j - 1, i + j - 1) + B[2] * E(2, i, j - 1, i + j);
//         }
// }