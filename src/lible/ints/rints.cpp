#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

using std::array, std::vector;

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

void LI::calcRIntsMatrixTest(const int l, const int n_cols, const int ofs_row, const int ofs_col,
                             const double fac, const double alpha, const double *xyz_pq,
                             const double *fnx, const vector<array<int, 3>> &hermite_idxs_bra,
                             const vector<array<int, 3>> &hermite_idxs_ket, double *rints)
{
    vec3d rints_3d = calcRInts3D(l, alpha, xyz_pq, fnx);

    int n_hermite_a = hermite_idxs_bra.size();
    int n_hermite_b = hermite_idxs_ket.size();

    for (int j = 0; j < n_hermite_b; j++)
    {
        auto [t_, u_, v_] = hermite_idxs_ket[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        for (int i = 0; i < n_hermite_a; i++)
        {
            auto [t, u, v] = hermite_idxs_bra[i];
            
            int irow = ofs_row + i;
            int icol = ofs_col + j;
            int idx = irow * n_cols + icol;
            rints[idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }
}

void LI::calcRInts_ERI2D1(const int l, const int n_rints, const int n_cols,
                          const int ofs_row, const int ofs_col, const double alpha,
                          const double fac, const double *fnx, const double *xyz_ab, 
                          const vector<array<int, 3>> &hermite_idxs_bra,
                          const vector<array<int, 3>> &hermite_idxs_ket, 
                          double *rints_out)
{
    vec3d rints_3d = calcRInts3D(l + 1, alpha, xyz_ab, fnx);

    int n_hermite_a = hermite_idxs_bra.size();
    int n_hermite_b = hermite_idxs_ket.size();
    int ofs0 = n_rints * 0;
    int ofs1 = n_rints * 1;
    int ofs2 = n_rints * 2;
    int ofs3 = n_rints * 3;
    int ofs4 = n_rints * 4;
    int ofs5 = n_rints * 5;

    for (int j = 0; j < n_hermite_b; j++)
    {
        auto [t_, u_, v_] = hermite_idxs_ket[j];

        double sign_A = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign_A = -1.0;

        double sign_B = sign_A * -1.0;

        for (int i = 0; i < n_hermite_a; i++)
        {
            auto [t, u, v] = hermite_idxs_bra[i];
            
            int irow = ofs_row + i;
            int icol = ofs_col + j;
            int idx = irow * n_cols + icol;

            // d/dA
            rints_out[ofs0 + idx] = sign_A * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + idx] = sign_A * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + idx] = sign_A * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dB
            rints_out[ofs3 + idx] = sign_B * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs4 + idx] = sign_B * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs5 + idx] = sign_B * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }
}

void LI::calcRInts_ERI3D1(const int l, const int n_rints, const int n_cols,
                          const int ofs_row, const int ofs_col, const double alpha, 
                          const double fac, const double *fnx, const double *xyz_pc,
                          const vector<array<int, 3>> &hermite_idxs_bra,
                          const vector<array<int, 3>> &hermite_idxs_ket,
                          double *rints_out)
{
    vec3d rints_3d = calcRInts3D(l + 1, alpha, xyz_pc, fnx);

    int n_hermite_ab = hermite_idxs_bra.size();
    int n_hermite_c = hermite_idxs_ket.size();
    int ofs0 = n_rints * 0;
    int ofs1 = n_rints * 1;
    int ofs2 = n_rints * 2;
    int ofs3 = n_rints * 3;
    int ofs4 = n_rints * 4;
    int ofs5 = n_rints * 5;    
    int ofs6 = n_rints * 6;

    for (int j = 0; j < n_hermite_c; j++)
    {
        auto &[t_, u_, v_] = hermite_idxs_ket[j];

        double sign_AB = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign_AB = -1.0;

        double sign_C = sign_AB * -1.0;

        for (int i = 0; i < n_hermite_ab; i++)
        {
            auto &[t, u, v] = hermite_idxs_bra[i];

            int irow = ofs_row + i;
            int icol = ofs_col + j;
            int idx = irow * n_cols + icol;

            // d/dP
            rints_out[ofs0 + idx] = sign_AB * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + idx] = sign_AB * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + idx] = sign_AB * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR
            rints_out[ofs3 + idx] = sign_AB * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dC
            rints_out[ofs4 + idx] = sign_C * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs5 + idx] = sign_C * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs6 + idx] = sign_C * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }
}

vector<double> LI::calcRInts_ERI4_Deriv1(const int l, const double fac, const double p,
                                         const double *xyz_pq, const double *fnx,
                                         const vector<array<int, 3>> &tuv_idxs_ab,
                                         const vector<array<int, 3>> &tuv_idxs_cd)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pq, fnx);

    int dim_tuv_ab = tuv_idxs_ab.size();
    int dim_tuv_cd = tuv_idxs_cd.size();
    int dim_tuv_abcd = dim_tuv_ab * dim_tuv_cd;

    int ofs0 = dim_tuv_abcd * 0;
    int ofs1 = dim_tuv_abcd * 1;
    int ofs2 = dim_tuv_abcd * 2;
    int ofs3 = dim_tuv_abcd * 3;
    int ofs4 = dim_tuv_abcd * 4;
    int ofs5 = dim_tuv_abcd * 5;
    int ofs6 = dim_tuv_abcd * 6;

    vector<double> rints_out(7 * dim_tuv_abcd, 0);
    for (size_t j = 0; j < tuv_idxs_cd.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_cd[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        double sign_ = sign * -1.0;

        for (size_t i = 0; i < tuv_idxs_ab.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_ab[i];

            // d/dP^(ab)
            rints_out[ofs0 + i * dim_tuv_cd + j] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + i * dim_tuv_cd + j] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + i * dim_tuv_cd + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR^(ab) and d/dR^(cd)
            rints_out[ofs3 + i * dim_tuv_cd + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dP^(cd)
            rints_out[ofs4 + i * dim_tuv_cd + j] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs5 + i * dim_tuv_cd + j] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs6 + i * dim_tuv_cd + j] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
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
    int ofs6 = dim_tuv_abc * 6;

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

            // d/dP
            rints_out[ofs0 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs1 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs2 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR
            rints_out[ofs3 + i * dim_tuv_c + j] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dC
            rints_out[ofs4 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints_out[ofs5 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints_out[ofs6 + i * dim_tuv_c + j] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }

    return rints_out;
}

void LI::calcRInts_ERI3_deriv1_test(const int l, const double fac, const double p,
                                    const double *xyz_pc, const double *fnx, const int n_rints,
                                    const int ofs_row, const int ofs_col, const int n_cols,
                                    const vector<array<int, 3>> &tuv_idxs_ab,
                                    const vector<array<int, 3>> &tuv_idxs_c,
                                    double *rints)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pc, fnx);

    int ofs0 = n_rints * 0;
    int ofs1 = n_rints * 1;
    int ofs2 = n_rints * 2;
    int ofs3 = n_rints * 3;
    int ofs4 = n_rints * 4;
    int ofs5 = n_rints * 5;
    int ofs6 = n_rints * 6;

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

            int irow = ofs_row + i;
            int icol = ofs_col + j;
            int idx = irow * n_cols + icol;

            // d/dP
            rints[ofs0 + idx] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs1 + idx] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs2 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR
            rints[ofs3 + idx] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            // d/dC
            rints[ofs4 + idx] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs5 + idx] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs6 + idx] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);
        }
    }
}

void LI::calcRInts_ERI4_deriv1_test(const int l, const double fac, const double p,
                                    const double *xyz_pc, const double *fnx, const int n_rints,
                                    const int ofs_row, const int ofs_col, const int n_cols,
                                    const int n_rows, const vector<array<int, 3>> &tuv_idxs_ab,
                                    const vector<array<int, 3>> &tuv_idxs_cd,
                                    double *rints)
{
    vec3d rints_3d = calcRInts3D(l + 1, p, xyz_pc, fnx);

    int ofs0 = n_rints * 0;
    int ofs1 = n_rints * 1;
    int ofs2 = n_rints * 2;
    int ofs3 = n_rints * 3;
    int ofs4 = n_rints * 4;
    int ofs5 = n_rints * 5;
    int ofs6 = n_rints * 6;
    int ofs7 = n_rints * 7;

    for (size_t j = 0; j < tuv_idxs_cd.size(); j++)
    {
        auto [t_, u_, v_] = tuv_idxs_cd[j];

        double sign = 1.0;
        if ((t_ + u_ + v_) % 2 != 0)
            sign = -1.0;

        double sign_ = sign * -1.0;

        for (size_t i = 0; i < tuv_idxs_ab.size(); i++)
        {
            auto [t, u, v] = tuv_idxs_ab[i];

            int irow = ofs_row + i;
            int icol = ofs_col + j;
            int idx_lhs = irow * n_cols + icol;

            // d/dP
            rints[ofs0 + idx_lhs] = sign * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs1 + idx_lhs] = sign * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs2 + idx_lhs] = sign * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dR
            rints[ofs3 + idx_lhs] = sign * fac * rints_3d(t + t_, u + u_, v + v_);

            int idx_lhs_tsp = icol * n_rows + irow;

            // d/dQ
            rints[ofs4 + idx_lhs_tsp] = sign_ * fac * rints_3d(t + t_ + 1, u + u_, v + v_);
            rints[ofs5 + idx_lhs_tsp] = sign_ * fac * rints_3d(t + t_, u + u_ + 1, v + v_);
            rints[ofs6 + idx_lhs_tsp] = sign_ * fac * rints_3d(t + t_, u + u_, v + v_ + 1);

            // d/dS
            rints[ofs7 + idx_lhs_tsp] = sign * fac * rints_3d(t + t_, u + u_, v + v_);
        }
    }
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