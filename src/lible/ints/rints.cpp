#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

lible::vec3d LI::calcRInts3DErf(const int l, const double p, const double omega,
                                const double *xyz_ab, const double *fnx)
{
    vec4d rints_buff(Fill(0), l + 1);

    rints_buff(0, 0, 0, 0) = fnx[0];

    double x = -2 * p * omega * omega / (p + omega * omega);
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

    vec3d rints(Fill(0), l + 1);
    for (int t = 0; t <= l; t++)
        for (int u = 0; u <= l - t; u++)
            for (int v = 0; v <= l - t - u; v++)
                rints(t, u, v) = rints_buff(0, t, u, v);

    return rints;
}
lible::vec3d LI::calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx)
{
    vec4d rints_buff(Fill(0), l + 1);

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

    vec3d rints(Fill(0), l + 1);
    for (int t = 0; t <= l; t++)
        for (int u = 0; u <= l - t; u++)
            for (int v = 0; v <= l - t - u; v++)
                rints(t, u, v) = rints_buff(0, t, u, v);

    return rints;
}

vector<double> LI::calcRIntsMatrix(const int l, const double fac, const double alpha,
                                   const double *xyz_pq, const double *fnx,
                                   const vector<array<int, 3>> &hermite_idxs_bra,
                                   const vector<array<int, 3>> &hermite_idxs_ket)
{
    vec3d rints_3d = calcRInts3D(l, alpha, xyz_pq, fnx);

    int n_hermite_bra = hermite_idxs_bra.size();
    int n_hermite_ket = hermite_idxs_ket.size();

    vector<double> rints_out(n_hermite_bra * n_hermite_ket, 0);
    for (size_t j = 0; j < hermite_idxs_ket.size(); j++)
    {
        auto [t_, u_, v_] = hermite_idxs_ket[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (size_t i = 0; i < hermite_idxs_bra.size(); i++)
        {
            auto [t, u, v] = hermite_idxs_bra[i];

            rints_out[i * n_hermite_ket + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }

    return rints_out;
}

vector<double> LI::calcRInts_ERI2D1(const int l, const double alpha, const double fac,
                                    const double *fnx, const double *xyz_ab,
                                    const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                    const std::vector<std::array<int, 3>> &hermite_idxs_b)
{
    vec3d rints_3d = calcRInts3D(l + 1, alpha, xyz_ab, fnx);

    const int n_hermite_a = hermite_idxs_a.size();
    const int n_hermite_b = hermite_idxs_b.size();
    const int n_rints = n_hermite_a * n_hermite_b;
    const int ofs0 = n_rints * 0;
    const int ofs1 = n_rints * 1;
    const int ofs2 = n_rints * 2;  

    vector<double> rints(3 * n_rints);
    for (int j = 0; j < n_hermite_b; j++)
    {
        auto [t_, u_, v_] = hermite_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (int i = 0; i < n_hermite_a; i++)
        {
            auto [t, u, v] = hermite_idxs_a[i];
            
            int idx = i * n_hermite_b + j;

            // d/dA
            rints[ofs0 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs1 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs2 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }

    return rints;
}

vector<double> LI::calcRInts_ERI3D1(const int l, const double alpha, const double fac,
                                    const double *fnx, const double *xyz_pc,
                                    const vector<array<int, 3>> &hermite_idxs_bra,
                                    const vector<array<int, 3>> &hermite_idxs_ket)
{
    vec3d rints_3d = calcRInts3D(l + 1, alpha, xyz_pc, fnx);

    const int n_hermite_ab = hermite_idxs_bra.size();
    const int n_hermite_c = hermite_idxs_ket.size();
    const int n_rints = n_hermite_ab * n_hermite_c;
    const int ofs0 = n_rints * 0;
    const int ofs1 = n_rints * 1;
    const int ofs2 = n_rints * 2;
    const int ofs3 = n_rints * 3;

    vector<double> rints(4 * n_rints);
    for (int j = 0; j < n_hermite_c; j++)
    {
        auto &[t_, u_, v_] = hermite_idxs_ket[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;        

        for (int i = 0; i < n_hermite_ab; i++)
        {
            auto &[t, u, v] = hermite_idxs_bra[i];

            int idx = i * n_hermite_c + j;

            // d/dP
            rints[ofs0 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs1 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs2 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR
            rints[ofs3 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }

    return rints;
}

vector<double> LI::calcRInts_ERI2D2(const int l, const double alpha, const double fac,
                                    const double *fnx, const double *xyz_ab,
                                    const vector<array<int, 3>> &hermite_idxs_a,
                                    const vector<array<int, 3>> &hermite_idxs_b)
{
    vec3d rints_3d = calcRInts3D(l + 2, alpha, xyz_ab, fnx);

    const int n_hermite_a = hermite_idxs_a.size();
    const int n_hermite_b = hermite_idxs_b.size();
    const int n_rints = n_hermite_a * n_hermite_b;
    const int ofs0 = n_rints * 0;
    const int ofs1 = n_rints * 1;
    const int ofs2 = n_rints * 2;
    const int ofs3 = n_rints * 3;
    const int ofs4 = n_rints * 4;
    const int ofs5 = n_rints * 5;

    vector<double> rints(6 * n_rints);
    for (int j = 0; j < n_hermite_b; j++)
    {
        auto [t_, u_, v_] = hermite_idxs_b[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (int i = 0; i < n_hermite_a; i++)
        {
            auto [t, u, v] = hermite_idxs_a[i];
            
            int idx = i * n_hermite_b + j;

            rints[ofs0 + idx] = sign * fac * rints_3d(t + t_ + 2, u + u_, v + v_);
            rints[ofs1 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_ + 1, v + v_);
            rints[ofs2 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_ + 1);
            rints[ofs3 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 2, v + v_);
            rints[ofs4 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_ + 1);
            rints[ofs5 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 2);
        }
    }

    return rints;
}

vector<double> LI::calcRInts_ERI3D2(const int l, const double alpha, const double fac,
                                    const double *fnx, const double *xyz_pc,
                                    const vector<array<int, 3>> &hermite_idxs_bra,
                                    const vector<array<int, 3>> &hermite_idxs_ket)
{
    vec3d rints_3d = calcRInts3D(l + 2, alpha, xyz_pc, fnx);

    const int n_hermite_bra = hermite_idxs_bra.size();
    const int n_hermite_ket = hermite_idxs_ket.size();
    const int n_rints = n_hermite_bra * n_hermite_ket;
    const int ofs0 = n_rints * 0;
    const int ofs1 = n_rints * 1;
    const int ofs2 = n_rints * 2;
    const int ofs3 = n_rints * 3;
    const int ofs4 = n_rints * 4;
    const int ofs5 = n_rints * 5;
    const int ofs6 = n_rints * 6;
    const int ofs7 = n_rints * 7;
    const int ofs8 = n_rints * 8;
    const int ofs9 = n_rints * 9;

    vector<double> rints(10 * n_rints);
    for (int j = 0; j < n_hermite_bra; j++)
    {
        auto [t_, u_, v_] = hermite_idxs_bra[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (int i = 0; i < n_hermite_ket; i++)
        {
            auto [t, u, v] = hermite_idxs_ket[i];

            int idx = i * n_hermite_ket + j;

            rints[ofs0 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
            
            rints[ofs1 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs2 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs3 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            rints[ofs4 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_ + 1);
            rints[ofs5 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 2);
            rints[ofs6 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_ + 1);
            rints[ofs7 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 2, v + v_);
            rints[ofs8 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_ + 1);
            rints[ofs9 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 2);
        }
    }

    return rints;    
}