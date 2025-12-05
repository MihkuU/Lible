#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

namespace lible::ints
{
    /// Calculates the Hermite integrals as a matrix, R(tuv, t'u'v') = R(t + t', u + u', v + v').
    std::vector<double> calcRIntsMatrix(int l, double fac, double p, const double *xyz_pq,
                                        const double *fnx,
                                        const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                        const std::vector<std::array<int, 3>> &hermite_idxs_b);

    /// Calculates the Hermite integrals as a 3D array, R(t, u, v). Based on eqs. (9.9.18)-(9.9.20)
    /// from DOI:10.1002/9781119019572.
    vec3d calcRInts3D(int l, double p, const double *xyz_ab, const double *fnx);

    /// Calculates the Hermite integrals with an attenuation parameter omega.
    vec3d calcRInts3DErf(int l, double p, double omega, const double *xyz_ab, const double *fnx);

    /// Calculates the Hermite integrals for ERI2 first derivative. Returns the matrix for
    /// each incremented Cartesian direction, (tt' + 1, uu', vv'), (tt', uu' + 1, vv') and
    /// (tt', uu', vv' + 1).
    std::vector<double> calcRInts_ERI2D1(int l, double alpha, double fac, const double *fnx,
                                         const double *xyz_ab,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_b);

    /// Calculates the Hermite integrals for ERI2 second derivative. Returns the matrix for each
    /// incremented Cartesian direction, (tt' + 2, uu', vv'), (tt' + 1, uu' + 1, vv'),
    /// (tt' + 1, uu', vv' + 1), (tt', uu' + 2, vv'), (tt', uu' + 1, vv' + 1) and
    /// (tt', uu' + 2, vv' + 2).
    std::vector<double> calcRInts_ERI2D2(int l, double alpha, double fac, const double *fnx,
                                         const double *xyz_ab,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_b);

    /// Calculates the Hermite integrals for ERI3 first derivative. Returns the matrix for
    /// the Cartesian directions (tt' + 1, uu', vv'), (tt', uu' + 1, vv') (tt', uu', vv' + 1) and
    /// (tt', uu', vv').
    std::vector<double> calcRInts_ERI3D1(int l, double alpha, double fac, const double *fnx,
                                         const double *xyz_pc,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_ket);

    /// Calculates Hermite integrals for ERI4 and ERI3 SOC integrals. Returns the matrix for
    /// each incremented Cartesian direction, (tt' + 1, uu', vv'), (tt', uu' + 1, vv') and
    /// (tt', uu', vv' + 1).
    std::vector<double> calcRInts_ERISOC(int l, double fac, double alpha, const double *xyz_pq,
                                         const double *fnx,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                         const std::vector<std::array<int, 3>> &hermite_idxs_ket);

    /// Helper function to calculate the index of a Hermite Gaussian in the list
    /// {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), ..., (t_max, u_max, v_max)}.
    constexpr int indexHermite(const int t, const int u, const int v)
    {
        int tuv = t + u + v;

        int idx_cart = indexCart(t, u, v);
        int offset = numHermites(tuv - 1);

        return offset + idx_cart;
    }

    /// Helper function to calculate the Hermite Gaussian index.
    constexpr int indexRRollout(const int lab, const int t, const int u, const int v)
    {
        int tuv = t + u + v;
        if (tuv == 0)
            return 0;

        int idx_cart = indexCart(t, u, v);
        int offset = lab + numHermites(tuv - 1);

        return offset + idx_cart;
    }

    /// Compile time helper function to calculate the list of Hermite index triplets,
    /// {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), ..., (t_max, u_max, v_max)}.
    template <int l>
    consteval std::array<std::array<int, 3>, numHermitesC(l)> generateHermiteIdxs()
    {
        constexpr int n_hermites = numHermitesC(l);

        std::array<std::array<int, 3>, n_hermites> hermite_idxs;
        for (int n = 0, ijk = 0; n <= l; n++)
            for (int i = n; i >= 0; i--)
                for (int j = n - i; j >= 0; j--, ijk++)
                {
                    int k = n - i - j;
                    hermite_idxs[ijk] = {i, j, k};
                }

        return hermite_idxs;
    }

    /// Compile time function to calculate the index for (t, u, v) in the list
    /// {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), ..., (t_max, u_max, v_max)}.
    inline int indexRR(const int offset, const int n, const int t, const int u, const int v)
    {
        int tuv = t + u + v;
        if (tuv == 0)
            return n;

        return offset + indexCart(t, u, v);
    }

    /// Templated function to calculate the R-function as R(t + t', u + u', v + v').
    template <int l>
    void calcRInts(const double alpha, const double *fnx, const double *xyz_pq, double *rints_out)
    {
        rints_out[0] = fnx[0];

        double x = -2 * alpha;
        double y = x;
        for (int n = 1; n <= l; n++)
        {
            rints_out[n] = fnx[n] * y;
            y *= x;
        }

        // Calc R-ints.
        for (int n = l - 1; n >= 0; n--)
            for (int m = l - n; m >= 1; m--)
            {
                int offset_lhs = l + numHermites(m - 1);
                int offset_rhs1 = l + numHermites(m - 2);
                int offset_rhs2 = l + numHermites(m - 3);
                for (int t = m; t >= 0; t--)
                    for (int u = m - t; u >= 0; u--)
                    {
                        int v = m - t - u;
                        int idx_lhs = offset_lhs + indexCart(t, u, v);

                        if (t > 0)
                        {
                            int idx_rhs1 = indexRR(offset_rhs1, n + 1, t - 1, u, v);
                            rints_out[idx_lhs] = xyz_pq[0] * rints_out[idx_rhs1];
                            if (t > 1)
                            {
                                int idx_rhs2 = indexRR(offset_rhs2, n + 1, t - 2, u, v);
                                rints_out[idx_lhs] += (t - 1) * rints_out[idx_rhs2];
                            }
                        }
                        else
                        {
                            if (u > 0)
                            {
                                int idx_rhs1 = indexRR(offset_rhs1, n + 1, t, u - 1, v);
                                rints_out[idx_lhs] = xyz_pq[1] * rints_out[idx_rhs1];
                                if (u > 1)
                                {
                                    int idx_rhs2 = indexRR(offset_rhs2, n + 1, t, u - 2, v);
                                    rints_out[idx_lhs] += (u - 1) * rints_out[idx_rhs2];
                                }
                            }
                            else if (v > 0)
                            {
                                int idx_rhs1 = indexRR(offset_rhs1, n + 1, t, u, v - 1);
                                rints_out[idx_lhs] = xyz_pq[2] * rints_out[idx_rhs1];
                                if (v > 1)
                                {
                                    int idx_rhs2 = indexRR(offset_rhs2, n + 1, t, u, v - 2);
                                    rints_out[idx_lhs] += (v - 1) * rints_out[idx_rhs2];
                                }
                            }
                        }
                    }
            }
    }

    /// Templated function to calculate the R-integrals as a matrix,
    /// R(tuv, t'u'v') = R(t + t', u + u', v + v').
    template <int la, int lb>
    void calcRInts_ERI(const double alpha, const double fac, const double *fnx,
                       const double *xyz_pq, double *rints_out)
    {
        constexpr int lab = la + lb;

        constexpr int buff_size = numHermitesC(lab) + lab;
        std::array<double, buff_size> rints_buff{};
        calcRInts<lab>(alpha, fnx, xyz_pq, &rints_buff[0]);

        constexpr int n_hermites_a = numHermitesC(la);
        constexpr int n_hermites_b = numHermitesC(lb);
        constexpr std::array<std::array<int, 3>, n_hermites_a> idxs_a = generateHermiteIdxs<la>();
        constexpr std::array<std::array<int, 3>, n_hermites_b> idxs_b = generateHermiteIdxs<lb>();
        for (int j = 0; j < n_hermites_b; j++)
        {
            auto &[t_, u_, v_] = idxs_b[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            for (int i = 0; i < n_hermites_a; i++)
            {
                auto &[t, u, v] = idxs_a[i];

                int tt_ = t + t_;
                int uu_ = u + u_;
                int vv_ = v + v_;

                int idx_lhs = i * n_hermites_b + j;
                int idx_rhs = indexRRollout(lab, tt_, uu_, vv_);

                rints_out[idx_lhs] = sign * fac * rints_buff[idx_rhs];
            }
        }
    }

    /// Templated function to calculate the R-integrals for ERI2 derivatives.
    /// Returned as a matrix for the list: (tt' + 1, uu', vv'), (tt', uu' + 1, vv') and
    /// (tt', uu', vv' + 1).
    template <int la, int lb>
    void calcRInts_ERI2D1(const double alpha, const double fac, const double *fnx,
                          const double *xyz_ab, double *rints)
    {
        constexpr int lab = la + lb;

        constexpr int buff_size = numHermitesC(lab + 1) + (lab + 1);
        std::array<double, buff_size> rints_buff{};
        calcRInts<lab + 1>(alpha, fnx, xyz_ab, &rints_buff[0]);

        constexpr int n_hermites_a = numHermitesC(la);
        constexpr int n_hermites_b = numHermitesC(lb);
        constexpr int n_rints = n_hermites_a * n_hermites_b;
        constexpr int ofs0 = n_rints * 0;
        constexpr int ofs1 = n_rints * 1;
        constexpr int ofs2 = n_rints * 2;

        constexpr std::array<std::array<int, 3>, n_hermites_a> idxs_a = generateHermiteIdxs<la>();
        constexpr std::array<std::array<int, 3>, n_hermites_b> idxs_b = generateHermiteIdxs<lb>();
        for (int j = 0; j < n_hermites_b; j++)
        {
            auto &[t_, u_, v_] = idxs_b[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            for (int i = 0; i < n_hermites_a; i++)
            {
                auto &[t, u, v] = idxs_a[i];

                int tt_ = t + t_;
                int uu_ = u + u_;
                int vv_ = v + v_;

                int idx_lhs = i * n_hermites_b + j;
                int idx_rhs0 = indexRRollout(lab + 1, tt_ + 1, uu_, vv_); // TODO: precalc offset
                int idx_rhs1 = indexRRollout(lab + 1, tt_, uu_ + 1, vv_); // TODO: precalc offset
                int idx_rhs2 = indexRRollout(lab + 1, tt_, uu_, vv_ + 1); // TODO: precalc offset

                // d/dA
                rints[ofs0 + idx_lhs] = sign * fac * rints_buff[idx_rhs0];
                rints[ofs1 + idx_lhs] = sign * fac * rints_buff[idx_rhs1];
                rints[ofs2 + idx_lhs] = sign * fac * rints_buff[idx_rhs2];
            }
        }
    }

    /// Templated function to calculate the R-integrals for ERI3 derivatives.
    /// Returned as a matrix for the list: (tt' + 1, uu', vv'), (tt', uu' + 1, vv'),
    /// (tt', uu', vv' + 1) and (tt', uu', vv').
    template <int lab, int lc>
    void calcRInts_ERI3D1(const double alpha, const double fac, const double *fnx,
                          const double *xyz_pc, double *rints)
    {
        constexpr int labc = lab + lc;

        constexpr int buff_size = numHermitesC(labc + 1) + (labc + 1);
        std::array<double, buff_size> rints_buff{};
        calcRInts<labc + 1>(alpha, fnx, xyz_pc, &rints_buff[0]);

        constexpr int n_hermites_ab = numHermitesC(lab);
        constexpr int n_hermites_c = numHermitesC(lc);
        constexpr int n_rints = n_hermites_ab * n_hermites_c;
        constexpr int ofs0 = n_rints * 0;
        constexpr int ofs1 = n_rints * 1;
        constexpr int ofs2 = n_rints * 2;
        constexpr int ofs3 = n_rints * 3;

        constexpr std::array<std::array<int, 3>, n_hermites_ab> idxs_ab = generateHermiteIdxs<lab>();
        constexpr std::array<std::array<int, 3>, n_hermites_c> idxs_c = generateHermiteIdxs<lc>();

        for (int j = 0; j < n_hermites_c; j++)
        {
            auto &[t_, u_, v_] = idxs_c[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            for (int i = 0; i < n_hermites_ab; i++)
            {
                auto &[t, u, v] = idxs_ab[i];

                int tt_ = t + t_;
                int uu_ = u + u_;
                int vv_ = v + v_;

                int idx_lhs = i * n_hermites_c + j;
                int idx_rhs0 = indexRRollout(labc + 1, tt_ + 1, uu_, vv_);
                int idx_rhs1 = indexRRollout(labc + 1, tt_, uu_ + 1, vv_);
                int idx_rhs2 = indexRRollout(labc + 1, tt_, uu_, vv_ + 1);
                int idx_rhs3 = indexRRollout(labc + 1, tt_, uu_, vv_);

                // d/dP
                rints[ofs0 + idx_lhs] = sign * fac * rints_buff[idx_rhs0];
                rints[ofs1 + idx_lhs] = sign * fac * rints_buff[idx_rhs1];
                rints[ofs2 + idx_lhs] = sign * fac * rints_buff[idx_rhs2];

                // d/dR
                rints[ofs3 + idx_lhs] = sign * fac * rints_buff[idx_rhs3];
            }
        }
    }

    /// Templated function to calculate the R-integrals for ERI2 derivatives.
    /// Returned as a matrix for the list: (tt' + 1, uu', vv'), (tt', uu' + 1, vv') and
    /// (tt', uu', vv' + 1).
    template <int la, int lb>
    void calcRInts_ERI2D2(const double alpha, const double fac, const double *fnx,
                          const double *xyz_ab, double *rints)
    {
        constexpr int lab = la + lb;

        constexpr int buff_size = numHermitesC(lab + 2) + (lab + 2);
        std::array<double, buff_size> rints_buff{};
        calcRInts<lab + 2>(alpha, fnx, xyz_ab, &rints_buff[0]);

        constexpr int n_hermites_a = numHermitesC(la);
        constexpr int n_hermites_b = numHermitesC(lb);
        constexpr int n_rints = n_hermites_a * n_hermites_b;
        constexpr int ofs0 = n_rints * 0;
        constexpr int ofs1 = n_rints * 1;
        constexpr int ofs2 = n_rints * 2;
        constexpr int ofs3 = n_rints * 3;
        constexpr int ofs4 = n_rints * 4;
        constexpr int ofs5 = n_rints * 5;

        constexpr std::array<std::array<int, 3>, n_hermites_a> idxs_a = generateHermiteIdxs<la>();
        constexpr std::array<std::array<int, 3>, n_hermites_b> idxs_b = generateHermiteIdxs<lb>();
        for (int j = 0; j < n_hermites_b; j++)
        {
            auto &[t_, u_, v_] = idxs_b[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            for (int i = 0; i < n_hermites_a; i++)
            {
                auto &[t, u, v] = idxs_a[i];

                int tt_ = t + t_;
                int uu_ = u + u_;
                int vv_ = v + v_;

                int idx_lhs = i * n_hermites_b + j;
                int idx_rhs0 = indexRRollout(lab + 2, tt_ + 2, uu_, vv_);
                int idx_rhs1 = indexRRollout(lab + 2, tt_ + 1, uu_ + 1, vv_);
                int idx_rhs2 = indexRRollout(lab + 2, tt_ + 1, uu_, vv_ + 1);
                int idx_rhs3 = indexRRollout(lab + 2, tt_, uu_ + 2, vv_);
                int idx_rhs4 = indexRRollout(lab + 2, tt_, uu_ + 1, vv_ + 1);
                int idx_rhs5 = indexRRollout(lab + 2, tt_, uu_, vv_ + 2);

                rints[ofs0 + idx_lhs] = sign * fac * rints_buff[idx_rhs0];
                rints[ofs1 + idx_lhs] = sign * fac * rints_buff[idx_rhs1];
                rints[ofs2 + idx_lhs] = sign * fac * rints_buff[idx_rhs2];
                rints[ofs3 + idx_lhs] = sign * fac * rints_buff[idx_rhs3];
                rints[ofs4 + idx_lhs] = sign * fac * rints_buff[idx_rhs4];
                rints[ofs5 + idx_lhs] = sign * fac * rints_buff[idx_rhs5];
            }
        }
    }

    /// Calculates Hermite integrals for ERI4 and ERI3 SOC integrals. Returns the matrix for
    /// each incremented Cartesian direction, (tt' + 1, uu', vv'), (tt', uu' + 1, vv') and
    /// (tt', uu', vv' + 1).
    template <int lbra, int lket>
    void calcRInts_ERISOC(const double alpha, const double fac, const double *fnx,
                          const double *xyz_pq, double *rints)
    {
        constexpr int l = lbra + lket;

        constexpr int buff_size = numHermitesC(l + 1) + (l + 1);
        std::array<double, buff_size> rints_buff{};
        calcRInts<l + 1>(alpha, fnx, xyz_pq, &rints_buff[0]);

        constexpr int n_hermites_bra = numHermitesC(lbra);
        constexpr int n_hermites_ket = numHermitesC(lket);
        constexpr int n_rints = n_hermites_bra * n_hermites_ket;
        constexpr int ofs0 = n_rints * 0;
        constexpr int ofs1 = n_rints * 1;
        constexpr int ofs2 = n_rints * 2;

        constexpr std::array<std::array<int, 3>, n_hermites_bra> idxs_bra = generateHermiteIdxs<lbra>();
        constexpr std::array<std::array<int, 3>, n_hermites_ket> idxs_ket = generateHermiteIdxs<lket>();

        for (int j = 0; j < n_hermites_ket; j++)
        {
            auto &[t_, u_, v_] = idxs_ket[j];

            double sign = 1.0;
            if ((t_ + u_ + v_) % 2 != 0)
                sign = -1.0;

            for (int i = 0; i < n_hermites_bra; i++)
            {
                auto &[t, u, v] = idxs_bra[i];

                int tt_ = t + t_;
                int uu_ = u + u_;
                int vv_ = v + v_;

                int idx_lhs = i * n_hermites_ket + j;
                int idx_rhs0 = indexRRollout(l + 1, tt_ + 1, uu_, vv_);
                int idx_rhs1 = indexRRollout(l + 1, tt_, uu_ + 1, vv_);
                int idx_rhs2 = indexRRollout(l + 1, tt_, uu_, vv_ + 1);

                // d/dP
                rints[ofs0 + idx_lhs] = sign * fac * rints_buff[idx_rhs0];
                rints[ofs1 + idx_lhs] = sign * fac * rints_buff[idx_rhs1];
                rints[ofs2 + idx_lhs] = sign * fac * rints_buff[idx_rhs2];
            }
        }
    }
}
