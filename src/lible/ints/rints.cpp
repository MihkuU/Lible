#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

void LI::calcRInts(const int la, const int lb, const double p, const array<double, 3> &xyz_ab,
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
                        rints_tmp(n, t, u, v) = xyz_ab[0] * rints_tmp(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_tmp(n, t, u, v) += (t - 1) * rints_tmp(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_tmp(n, t, u, v) = xyz_ab[1] * rints_tmp(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_tmp(n, t, u, v) += (u - 1) * rints_tmp(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_tmp(n, t, u, v) = xyz_ab[2] * rints_tmp(n + 1, t, u, v - 1);
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

void LI::calcRInts(const int la, const int lb, const double fac, const double p,
                   const array<double, 3> &xyz_pq, const vector<double> &fnx,
                   const vector<array<int, 3>> &tuv_idxs_a,
                   const vector<array<int, 3>> &tuv_idxs_b,
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
                        rints_tmp(n, t, u, v) = xyz_pq[0] * rints_tmp(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_tmp(n, t, u, v) += (t - 1) * rints_tmp(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_tmp(n, t, u, v) = xyz_pq[1] * rints_tmp(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_tmp(n, t, u, v) += (u - 1) * rints_tmp(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_tmp(n, t, u, v) = xyz_pq[2] * rints_tmp(n + 1, t, u, v - 1);
                            if (v > 1)
                                rints_tmp(n, t, u, v) += (v - 1) * rints_tmp(n + 1, t, u, v - 2);
                        }
                    }
                }
    }

    int dim_tuv_b = dimHermiteGaussians(lb);
    for (size_t j = 0; j < tuv_idxs_b.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (size_t i = 0; i < tuv_idxs_a.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_a[i];

            rints_out[i * dim_tuv_b + j] = sign * fac * rints_tmp(0, t + t_, u + u_, v + v_);
        }
    }
}