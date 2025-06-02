#include <lible/ints/rints_meta.hpp>

template<>
void lible::ints::calcRInts_ERI_new<1, 1>(const double alpha, const double fac, const double *fnx,
                                          const double *xyz_ab, const int n_cols, const int ofs_row, 
                                          const int ofs_col, double *rints_out)
{
    constexpr int lab = 2;
    constexpr int buff_size = 12;
    std::array<double, buff_size> rints_buff{};

    rints_buff[0] = fnx[0];
    double x = -2 * alpha;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_buff[n] = fnx[n] * y;
        y *= x;
    }

    // R-ints recursion
    rints_buff[3] = xyz_ab[0] * rints_buff[2];
    rints_buff[4] = xyz_ab[1] * rints_buff[2];
    rints_buff[5] = xyz_ab[2] * rints_buff[2];
    rints_buff[6] = xyz_ab[0] * rints_buff[3] + 1 * rints_buff[1];
    rints_buff[7] = xyz_ab[0] * rints_buff[4];
    rints_buff[8] = xyz_ab[0] * rints_buff[5];
    rints_buff[9] = xyz_ab[1] * rints_buff[4] + 1 * rints_buff[1];
    rints_buff[10] = xyz_ab[1] * rints_buff[5];
    rints_buff[11] = xyz_ab[2] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[3] = xyz_ab[0] * rints_buff[1];
    rints_buff[4] = xyz_ab[1] * rints_buff[1];
    rints_buff[5] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[(ofs_row + 0) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[0];
    rints_out[(ofs_row + 1) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[3];
    rints_out[(ofs_row + 2) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[4];
    rints_out[(ofs_row + 3) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[5];
    rints_out[(ofs_row + 0) * n_cols + (ofs_col + 1)] = -1.0 * fac * rints_buff[3];
    rints_out[(ofs_row + 1) * n_cols + (ofs_col + 1)] = -1.0 * fac * rints_buff[6];
    rints_out[(ofs_row + 2) * n_cols + (ofs_col + 1)] = -1.0 * fac * rints_buff[7];
    rints_out[(ofs_row + 3) * n_cols + (ofs_col + 1)] = -1.0 * fac * rints_buff[8];
    rints_out[(ofs_row + 0) * n_cols + (ofs_col + 2)] = -1.0 * fac * rints_buff[4];
    rints_out[(ofs_row + 1) * n_cols + (ofs_col + 2)] = -1.0 * fac * rints_buff[7];
    rints_out[(ofs_row + 2) * n_cols + (ofs_col + 2)] = -1.0 * fac * rints_buff[9];
    rints_out[(ofs_row + 3) * n_cols + (ofs_col + 2)] = -1.0 * fac * rints_buff[10];
    rints_out[(ofs_row + 0) * n_cols + (ofs_col + 3)] = -1.0 * fac * rints_buff[5];
    rints_out[(ofs_row + 1) * n_cols + (ofs_col + 3)] = -1.0 * fac * rints_buff[8];
    rints_out[(ofs_row + 2) * n_cols + (ofs_col + 3)] = -1.0 * fac * rints_buff[10];
    rints_out[(ofs_row + 3) * n_cols + (ofs_col + 3)] = -1.0 * fac * rints_buff[11];
}

namespace lible::ints
{
template <int la, int lb>
void calcRInts_ERI(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
                   double *rints_out);

template<>
void calcRInts_ERI<1, 1>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 2;
    constexpr int buff_size = 12;
    std::array<double, buff_size> rints_buff{};

    rints_buff[0] = fnx[0];
    double x = -2 * alpha;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_buff[n] = fnx[n] * y;
        y *= x;
    }

    // R-ints recursion
    rints_buff[3] = xyz_ab[0] * rints_buff[2];
    rints_buff[4] = xyz_ab[1] * rints_buff[2];
    rints_buff[5] = xyz_ab[2] * rints_buff[2];
    rints_buff[6] = xyz_ab[0] * rints_buff[3] + 1 * rints_buff[1];
    rints_buff[7] = xyz_ab[0] * rints_buff[4];
    rints_buff[8] = xyz_ab[0] * rints_buff[5];
    rints_buff[9] = xyz_ab[1] * rints_buff[4] + 1 * rints_buff[1];
    rints_buff[10] = xyz_ab[1] * rints_buff[5];
    rints_buff[11] = xyz_ab[2] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[3] = xyz_ab[0] * rints_buff[1];
    rints_buff[4] = xyz_ab[1] * rints_buff[1];
    rints_buff[5] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[4] = 1.0 * fac * rints_buff[3];
    rints_out[8] = 1.0 * fac * rints_buff[4];
    rints_out[12] = 1.0 * fac * rints_buff[5];
    rints_out[1] = -1.0 * fac * rints_buff[3];
    rints_out[5] = -1.0 * fac * rints_buff[6];
    rints_out[9] = -1.0 * fac * rints_buff[7];
    rints_out[13] = -1.0 * fac * rints_buff[8];
    rints_out[2] = -1.0 * fac * rints_buff[4];
    rints_out[6] = -1.0 * fac * rints_buff[7];
    rints_out[10] = -1.0 * fac * rints_buff[9];
    rints_out[14] = -1.0 * fac * rints_buff[10];
    rints_out[3] = -1.0 * fac * rints_buff[5];
    rints_out[7] = -1.0 * fac * rints_buff[8];
    rints_out[11] = -1.0 * fac * rints_buff[10];
    rints_out[15] = -1.0 * fac * rints_buff[11];
}
}

template void lible::ints::calcRInts_ERI2_deriv1<1, 1>(const double, const double, const double*, const double*, double*);
template void lible::ints::calcRInts_ERI2D1<1, 1>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int n_cols,
                                                  const int ofs_row, const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI3D1<1, 1>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int n_cols,
                                                  const int ofs_row, const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI4D1<1, 1>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int ofs_row,
                                                  const int ofs_col, const int n_cols, const int n_rows,
                                                  double *rints_out);

