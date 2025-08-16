#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         * Calculates the Hermite Coulomb integrals as a flattened matrix,
         *   R(t + t', u + u', v + v') -> R(tuv, t'u'v').
         * This is used in the SHARK method.
         */
        std::vector<double> calcRIntsMatrix(const int l, const double fac, const double p,
                                            const double *xyz_pq, const double *fnx,
                                            const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                            const std::vector<std::array<int, 3>> &hermite_idxs_b);

        /** Calculates the Hermite Coulomb integrals as a 3D array R(t, u, v). */
        vec3d calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx);

        vec3d calcRInts3DErf(const int l, const double p, const double omega, const double *xyz_ab, const double *fnx);

        std::vector<double> calcRInts_ERI2D1(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_ab,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_b);

        std::vector<double> calcRInts_ERI2D2(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_ab,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_b);

        std::vector<double> calcRInts_ERI3D1(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_pc,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_ket);

        std::vector<double> calcRInts_ERI3D2(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_pc,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_ket);

        std::vector<double> calcRInts_ERISOC(const int l, const double fac, const double alpha,
                                             const double *xyz_pq, const double *fnx,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_ket);

        /** */
        constexpr int indexHermite(int t, int u, int v)
        {
            int tuv = t + u + v;

            int idx_cart = indexCart(t, u, v);
            int offset = numHermites(tuv - 1);

            return offset + idx_cart;
        }

        constexpr int indexRRollout(int lab, int t, int u, int v)
        {
            int tuv = t + u + v;
            if (tuv == 0)
                return 0;

            int idx_cart = indexCart(t, u, v);
            int offset = lab + numHermites(tuv - 1);

            return offset + idx_cart;
        }

        /** */
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

        /** Helper function for calculating the index */
        inline int indexRR(int offset, int n, int t, int u, int v)
        {
            int tuv = t + u + v;
            if (tuv == 0)
                return n;
            else
                return offset + indexCart(t, u, v);
        }

        /** */
        template <int la, int lb>
        void calcRInts(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
                       double *rints_out)
        {
            using std::array;

            // Calc R ints (0, 0, 0, 0) to (n, 0, 0, 0).
            constexpr int lab = la + lb;
            constexpr int buff_size = numHermitesC(lab) + lab;
            array<double, buff_size> rints_buff{};

            rints_buff[0] = fnx[0];

            double x = -2 * alpha;
            double y = x;
            for (int n = 1; n <= lab; n++)
            {
                rints_buff[n] = fnx[n] * y;
                y *= x;
            }

            // Calc R-ints.
            for (int n = lab - 1; n >= 0; n--)
            {
                int n_ = lab - n;
                for (int m = n_; m >= 1; m--)
                {
                    int offset_lhs = lab + numHermites(m - 1);
                    int offset_rhs1 = lab + numHermites(m - 2);
                    int offset_rhs2 = lab + numHermites(m - 3);
                    for (int t = m; t >= 0; t--)
                        for (int u = m - t; u >= 0; u--)
                        {
                            int v = m - t - u;
                            int idx_lhs = offset_lhs + indexCart(t, u, v);

                            if (t > 0)
                            {
                                int idx_rhs1 = indexRR(offset_rhs1, n + 1, t - 1, u, v);
                                rints_buff[idx_lhs] = xyz_ab[0] * rints_buff[idx_rhs1];
                                if (t > 1)
                                {
                                    int idx_rhs2 = indexRR(offset_rhs2, n + 1, t - 2, u, v);
                                    rints_buff[idx_lhs] += (t - 1) * rints_buff[idx_rhs2];
                                }
                            }
                            else
                            {
                                if (u > 0)
                                {
                                    int idx_rhs1 = indexRR(offset_rhs1, n + 1, t, u - 1, v);
                                    rints_buff[idx_lhs] = xyz_ab[1] * rints_buff[idx_rhs1];
                                    if (u > 1)
                                    {
                                        int idx_rhs2 = indexRR(offset_rhs2, n + 1, t, u - 2, v);
                                        rints_buff[idx_lhs] += (u - 1) * rints_buff[idx_rhs2];
                                    }
                                }
                                else if (v > 0)
                                {
                                    int idx_rhs1 = indexRR(offset_rhs1, n + 1, t, u, v - 1);
                                    rints_buff[idx_lhs] = xyz_ab[2] * rints_buff[idx_rhs1];
                                    if (v > 1)
                                    {
                                        int idx_rhs2 = indexRR(offset_rhs2, n + 1, t, u, v - 2);
                                        rints_buff[idx_lhs] += (v - 1) * rints_buff[idx_rhs2];
                                    }
                                }
                            }
                        }
                }
            }

            // Roll out R-ints.
            constexpr int n_hermites_a = numHermitesC(la);
            constexpr int n_hermites_b = numHermitesC(lb);
            constexpr array<array<int, 3>, n_hermites_a> idxs_a = generateHermiteIdxs<la>();
            constexpr array<array<int, 3>, n_hermites_b> idxs_b = generateHermiteIdxs<lb>();
            for (int j = 0; j < n_hermites_b; j++)
            {
                auto [t_, u_, v_] = idxs_b[j];

                double sign = 1.0;
                if ((t_ + u_ + v_) % 2 != 0)
                    sign = -1.0;

                for (int i = 0; i < n_hermites_a; i++)
                {
                    auto [t, u, v] = idxs_a[i];

                    int tt_ = t + t_;
                    int uu_ = u + u_;
                    int vv_ = v + v_;

                    int idx_lhs = i * n_hermites_b + j;
                    int idx_rhs = indexRRollout(lab, tt_, uu_, vv_);

                    rints_out[idx_lhs] = sign * fac * rints_buff[idx_rhs];
                }
            }
        }

        template <int l>
        void calcRInts(const double alpha, const double fac, const double *fnx, const double *xyz_pq,
                       double *rints_out)
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

        template <int la, int lb>
        void calcRInts_ERI(const double alpha, const double fac, const double *fnx,
                           const double *xyz_pq, double *rints_out)
        {
            constexpr int lab = la + lb;

            constexpr int buff_size = numHermitesC(lab) + lab;
            std::array<double, buff_size> rints_buff{};
            calcRInts<lab>(alpha, fac, fnx, xyz_pq, &rints_buff[0]);

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

        template <int la, int lb>
        void calcRInts_ERI2D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_ab, double *rints)
        {
            constexpr int lab = la + lb;

            constexpr int buff_size = numHermitesC(lab + 1) + (lab + 1);
            std::array<double, buff_size> rints_buff{};
            calcRInts<lab + 1>(alpha, fac, fnx, xyz_ab, &rints_buff[0]);

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

        template <int lab, int lc>
        void calcRInts_ERI3D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pc, double *rints)
        {
            constexpr int labc = lab + lc;

            constexpr int buff_size = numHermitesC(labc + 1) + (labc + 1);
            std::array<double, buff_size> rints_buff{};
            calcRInts<labc + 1>(alpha, fac, fnx, xyz_pc, &rints_buff[0]);

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

        template <int la, int lb>
        void calcRInts_ERI2D2(const double alpha, const double fac, const double *fnx,
                              const double *xyz_ab, double *rints)
        {
            constexpr int lab = la + lb;

            constexpr int buff_size = numHermitesC(lab + 2) + (lab + 2);
            std::array<double, buff_size> rints_buff{};
            calcRInts<lab + 2>(alpha, fac, fnx, xyz_ab, &rints_buff[0]);

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

        template <int lab, int lc>
        void calcRInts_ERI3D2(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pc, double *rints)
        {
            constexpr int labc = lab + lc;

            constexpr int buff_size = numHermitesC(labc + 2) + (labc + 2);
            std::array<double, buff_size> rints_buff{};
            calcRInts<labc + 2>(alpha, fac, fnx, xyz_pc, &rints_buff[0]);

            constexpr int n_hermites_ab = numHermitesC(lab);
            constexpr int n_hermites_c = numHermitesC(lc);
            constexpr int n_rints = n_hermites_ab * n_hermites_c;
            constexpr int ofs0 = n_rints * 0;
            constexpr int ofs1 = n_rints * 1;
            constexpr int ofs2 = n_rints * 2;
            constexpr int ofs3 = n_rints * 3;
            constexpr int ofs4 = n_rints * 4;
            constexpr int ofs5 = n_rints * 5;
            constexpr int ofs6 = n_rints * 6;
            constexpr int ofs7 = n_rints * 7;
            constexpr int ofs8 = n_rints * 8;
            constexpr int ofs9 = n_rints * 9;

            constexpr std::array<std::array<int, 3>, n_hermites_ab> idxs_ab = generateHermiteIdxs<lab>();
            constexpr std::array<std::array<int, 3>, n_hermites_c> idxs_c = generateHermiteIdxs<lc>();

            for (int j = 0; j < n_hermites_ab; j++)
            {
                auto [t_, u_, v_] = idxs_ab[j];

                double sign = 1.0;
                if ((t_ + u_ + v_) % 2 != 0)
                    sign = -1.0;

                for (int i = 0; i < n_hermites_c; i++)
                {
                    auto [t, u, v] = idxs_c[i];

                    // int idx = i * n_hermite_c + j;
                    int tt_ = t + t_;
                    int uu_ = u + u_;
                    int vv_ = v + v_;

                    int idx_lhs = i * n_hermites_c + j;

                    rints[ofs0 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_, vv_);

                    rints[ofs1 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_ + 1, uu_, vv_);
                    rints[ofs2 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_ + 1, vv_);
                    rints[ofs3 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_, vv_ + 1);

                    rints[ofs4 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_ + 2, uu_, vv_);
                    rints[ofs5 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_ + 1, uu_ + 1, vv_);
                    rints[ofs6 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_ + 1, uu_, vv_ + 1);
                    rints[ofs7 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_ + 2, vv_);
                    rints[ofs8 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_ + 1, vv_ + 1);
                    rints[ofs9 + idx_lhs] = sign * fac * indexRRollout(labc + 2, tt_, uu_, vv_ + 2);
                }
            }
        }

        template <int lbra, int lket>
        void calcRInts_ERISOC(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pq, double *rints)
        {
            constexpr int l = lbra + lket;

            constexpr int buff_size = numHermitesC(l + 1) + (l + 1);
            std::array<double, buff_size> rints_buff{};
            calcRInts<l + 1>(alpha, fac, fnx, xyz_pq, &rints_buff[0]);

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
}