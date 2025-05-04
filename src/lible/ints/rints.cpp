#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

lible::vec3d LI::calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx)
{
    vec4d rints_buff(l + 1, 0);

    rints_buff(0, 0, 0, 0) = fnx[0];

    double x = -2 * p;
    double y = x;
    for (int n = 1; n <= l; n++)
    {
        rints_buff(n, 0, 0, 0) = fnx[n] * y;
        y *= x;
    }

    // This clever loop is taken from HUMMR:
    for (int n = l - 1; n >= 0; n--)
    {
        int n_ = l - n;
        for (int t = 0; t <= n_; t++)
            for (int u = 0; u <= n_ - t; u++)
                for (int v = 0; v <= n_ - t - u; v++)
                {
                    if (t > 0)
                    {
                        rints_buff(n, t, u, v) = xyz_ab[0] * rints_buff(n + 1, t - 1, u, v);
                        if (t > 1)
                            rints_buff(n, t, u, v) += (t - 1) * rints_buff(n + 1, t - 2, u, v);
                    }
                    else
                    {
                        if (u > 0)
                        {
                            rints_buff(n, t, u, v) = xyz_ab[1] * rints_buff(n + 1, t, u - 1, v);
                            if (u > 1)
                                rints_buff(n, t, u, v) += (u - 1) * rints_buff(n + 1, t, u - 2, v);
                        }
                        else if (v > 0)
                        {
                            rints_buff(n, t, u, v) = xyz_ab[2] * rints_buff(n + 1, t, u, v - 1);
                            if (v > 1)
                                rints_buff(n, t, u, v) += (v - 1) * rints_buff(n + 1, t, u, v - 2);
                        }
                    }
                }
    }

    vec3d rints(l + 1, 0);
    for (int t = 0; t <= l; t++)
        for (int u = 0; u <= l - t; u++)
            for (int v = 0; v <= l - t - u; v++)
                rints(t, u, v) = rints_buff(0, t, u, v);

    return rints;
}

vector<double> LI::calcRIntsMatrix(const int l, const double fac, const double p,
                                   const double *xyz_pq, const double *fnx,
                                   const vector<array<int, 3>> &tuv_idxs_a,
                                   const vector<array<int, 3>> &tuv_idxs_b)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pq, fnx);

    int dim_tuv_a = tuv_idxs_a.size();
    int dim_tuv_b = tuv_idxs_b.size();

    vector<double> rints_out(dim_tuv_a * dim_tuv_b, 0);
    for (size_t j = 0; j < tuv_idxs_b.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (size_t i = 0; i < tuv_idxs_a.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_a[i];

            rints_out[i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }

    return rints_out;
}

vector<double> LI::calcRInts_ERI4_Deriv1(const int l, const double fac, const double p,
                                         const double *xyz_pq, const double *fnx,
                                         const vector<array<int, 3>> &tuv_idxs_a,
                                         const vector<array<int, 3>> &tuv_idxs_b)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pq, fnx);

    int dim_tuv_a = tuv_idxs_a.size();
    int dim_tuv_b = tuv_idxs_b.size();
    int dim_tuv_ab = dim_tuv_a * dim_tuv_b;

    int ofs0 = dim_tuv_ab * 0;
    int ofs1 = dim_tuv_ab * 1;
    int ofs2 = dim_tuv_ab * 2;
    int ofs3 = dim_tuv_ab * 3;
    int ofs4 = dim_tuv_ab * 4;
    int ofs5 = dim_tuv_ab * 5;
    int ofs6 = dim_tuv_ab * 6;
    int ofs9 = dim_tuv_ab * 9;

    vector<double> rints_out(12 * dim_tuv_ab, 0);
    for (size_t j = 0; j < tuv_idxs_b.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        double sign_ = sign * -1.0;

        for (size_t i = 0; i < tuv_idxs_a.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_a[i];

            // d/dP^(ab)
            rints_out[ofs0 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR^(ab)
            rints_out[ofs3 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
            rints_out[ofs4 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
            rints_out[ofs5 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dP^(cd)
            rints_out[ofs6 + (i + 0) * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs6 + (i + 1) * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs6 + (i + 2) * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR^(cd)
            rints_out[ofs9 + (i + 0) * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
            rints_out[ofs9 + (i + 1) * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
            rints_out[ofs9 + (i + 2) * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }

    return rints_out;
}

vector<double> LI::calcRInts_ERI3_deriv1(const int l, const double fac, const double p,
                                         const double *xyz_pc, const double *fnx,
                                         const vector<array<int, 3>> &tuv_idxs_ab,
                                         const vector<array<int, 3>> &tuv_idxs_c)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pc, fnx);

    int dim_tuv_ab = tuv_idxs_ab.size();
    int dim_tuv_c = tuv_idxs_c.size();
    int dim_tuv_abc = dim_tuv_ab * dim_tuv_c;

    int ofs0 = dim_tuv_abc * 0;
    int ofs1 = dim_tuv_abc * 1;
    int ofs2 = dim_tuv_abc * 2;
    int ofs3 = dim_tuv_abc * 3;
    int ofs4 = dim_tuv_abc * 4;
    int ofs5 = dim_tuv_abc * 5;
    int ofs6 = dim_tuv_abc * 5;

    vector<double> rints_out(7 * dim_tuv_abc, 0);
    for (size_t j = 0; j < tuv_idxs_c.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_c[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        double sign_ = sign * -1.0;

        for (size_t i = 0; i < tuv_idxs_ab.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_ab[i];

            // d/dA
            rints_out[ofs0 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dB
            rints_out[ofs3 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dC
            rints_out[ofs4 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs5 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs6 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }

    return rints_out;    
}

vector<double> LI::calcRInts_ERI2_deriv1(const int l, const double fac, const double p,
                                         const double *xyz_ab, const double *fnx,
                                         const vector<array<int, 3>> &tuv_idxs_a,
                                         const vector<array<int, 3>> &tuv_idxs_b)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_ab, fnx);

    int dim_tuv_a = tuv_idxs_a.size();
    int dim_tuv_b = tuv_idxs_b.size();
    int dim_tuv_ab = dim_tuv_a * dim_tuv_b;

    int ofs0 = dim_tuv_ab * 0;
    int ofs1 = dim_tuv_ab * 1;
    int ofs2 = dim_tuv_ab * 2;
    int ofs3 = dim_tuv_ab * 3;
    int ofs4 = dim_tuv_ab * 4;
    int ofs5 = dim_tuv_ab * 5;

    vector<double> rints_out(6 * dim_tuv_ab, 0);
    for (size_t j = 0; j < tuv_idxs_b.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        double sign_ = sign * -1.0;

        for (size_t i = 0; i < tuv_idxs_a.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_a[i];

            // d/dA
            rints_out[ofs0 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + i * dim_tuv_b + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dB
            rints_out[ofs3 + i * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs4 + i * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs5 + i * dim_tuv_b + j] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }

    return rints_out;
}