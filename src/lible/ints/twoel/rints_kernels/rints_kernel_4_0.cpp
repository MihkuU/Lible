#include <lible/ints/rints_meta.hpp>

template<>
void lible::ints::calcRInts_ERI_new<4, 0>(const double alpha, const double fac, const double *fnx,
                                          const double *xyz_ab, const int n_cols, const int ofs_row, 
                                          const int ofs_col, double *rints_out)
{
    constexpr int lab = 4;
    constexpr int buff_size = 39;
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
    rints_buff[5] = xyz_ab[0] * rints_buff[4];
    rints_buff[6] = xyz_ab[1] * rints_buff[4];
    rints_buff[7] = xyz_ab[2] * rints_buff[4];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[3];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[3];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[3];
    rints_buff[5] = xyz_ab[0] * rints_buff[3];
    rints_buff[6] = xyz_ab[1] * rints_buff[3];
    rints_buff[7] = xyz_ab[2] * rints_buff[3];
    rints_buff[14] = xyz_ab[0] * rints_buff[8] + 2 * rints_buff[5];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 1 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[1] * rints_buff[11] + 2 * rints_buff[6];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 1 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13];
    rints_buff[23] = xyz_ab[2] * rints_buff[13] + 2 * rints_buff[7];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[2];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[2];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[2];
    rints_buff[5] = xyz_ab[0] * rints_buff[2];
    rints_buff[6] = xyz_ab[1] * rints_buff[2];
    rints_buff[7] = xyz_ab[2] * rints_buff[2];
    rints_buff[24] = xyz_ab[0] * rints_buff[14] + 3 * rints_buff[8];
    rints_buff[25] = xyz_ab[0] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 2 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 1 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 1 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20];
    rints_buff[31] = xyz_ab[0] * rints_buff[21];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[1] * rints_buff[20] + 3 * rints_buff[11];
    rints_buff[35] = xyz_ab[1] * rints_buff[21] + 2 * rints_buff[12];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 1 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23];
    rints_buff[38] = xyz_ab[2] * rints_buff[23] + 3 * rints_buff[13];
    rints_buff[14] = xyz_ab[0] * rints_buff[8] + 2 * rints_buff[5];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 1 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[1] * rints_buff[11] + 2 * rints_buff[6];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 1 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13];
    rints_buff[23] = xyz_ab[2] * rints_buff[13] + 2 * rints_buff[7];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[1];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[1];
    rints_buff[5] = xyz_ab[0] * rints_buff[1];
    rints_buff[6] = xyz_ab[1] * rints_buff[1];
    rints_buff[7] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[(ofs_row + 0) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[0];
    rints_out[(ofs_row + 1) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[5];
    rints_out[(ofs_row + 2) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[6];
    rints_out[(ofs_row + 3) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[7];
    rints_out[(ofs_row + 4) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[8];
    rints_out[(ofs_row + 5) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[9];
    rints_out[(ofs_row + 6) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[10];
    rints_out[(ofs_row + 7) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[11];
    rints_out[(ofs_row + 8) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[12];
    rints_out[(ofs_row + 9) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[13];
    rints_out[(ofs_row + 10) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[14];
    rints_out[(ofs_row + 11) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[15];
    rints_out[(ofs_row + 12) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[16];
    rints_out[(ofs_row + 13) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[17];
    rints_out[(ofs_row + 14) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[18];
    rints_out[(ofs_row + 15) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[19];
    rints_out[(ofs_row + 16) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[20];
    rints_out[(ofs_row + 17) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[21];
    rints_out[(ofs_row + 18) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[22];
    rints_out[(ofs_row + 19) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[23];
    rints_out[(ofs_row + 20) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[24];
    rints_out[(ofs_row + 21) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[25];
    rints_out[(ofs_row + 22) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[26];
    rints_out[(ofs_row + 23) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[27];
    rints_out[(ofs_row + 24) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[28];
    rints_out[(ofs_row + 25) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[29];
    rints_out[(ofs_row + 26) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[30];
    rints_out[(ofs_row + 27) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[31];
    rints_out[(ofs_row + 28) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[32];
    rints_out[(ofs_row + 29) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[33];
    rints_out[(ofs_row + 30) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[34];
    rints_out[(ofs_row + 31) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[35];
    rints_out[(ofs_row + 32) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[36];
    rints_out[(ofs_row + 33) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[37];
    rints_out[(ofs_row + 34) * n_cols + (ofs_col + 0)] = 1.0 * fac * rints_buff[38];
}

namespace lible::ints
{
template <int la, int lb>
void calcRInts_ERI(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
                   double *rints_out);

template<>
void calcRInts_ERI<4, 0>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 4;
    constexpr int buff_size = 39;
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
    rints_buff[5] = xyz_ab[0] * rints_buff[4];
    rints_buff[6] = xyz_ab[1] * rints_buff[4];
    rints_buff[7] = xyz_ab[2] * rints_buff[4];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[3];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[3];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[3];
    rints_buff[5] = xyz_ab[0] * rints_buff[3];
    rints_buff[6] = xyz_ab[1] * rints_buff[3];
    rints_buff[7] = xyz_ab[2] * rints_buff[3];
    rints_buff[14] = xyz_ab[0] * rints_buff[8] + 2 * rints_buff[5];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 1 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[1] * rints_buff[11] + 2 * rints_buff[6];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 1 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13];
    rints_buff[23] = xyz_ab[2] * rints_buff[13] + 2 * rints_buff[7];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[2];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[2];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[2];
    rints_buff[5] = xyz_ab[0] * rints_buff[2];
    rints_buff[6] = xyz_ab[1] * rints_buff[2];
    rints_buff[7] = xyz_ab[2] * rints_buff[2];
    rints_buff[24] = xyz_ab[0] * rints_buff[14] + 3 * rints_buff[8];
    rints_buff[25] = xyz_ab[0] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 2 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 1 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 1 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20];
    rints_buff[31] = xyz_ab[0] * rints_buff[21];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[1] * rints_buff[20] + 3 * rints_buff[11];
    rints_buff[35] = xyz_ab[1] * rints_buff[21] + 2 * rints_buff[12];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 1 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23];
    rints_buff[38] = xyz_ab[2] * rints_buff[23] + 3 * rints_buff[13];
    rints_buff[14] = xyz_ab[0] * rints_buff[8] + 2 * rints_buff[5];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 1 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[1] * rints_buff[11] + 2 * rints_buff[6];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 1 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13];
    rints_buff[23] = xyz_ab[2] * rints_buff[13] + 2 * rints_buff[7];
    rints_buff[8] = xyz_ab[0] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[1] * rints_buff[6] + 1 * rints_buff[1];
    rints_buff[12] = xyz_ab[1] * rints_buff[7];
    rints_buff[13] = xyz_ab[2] * rints_buff[7] + 1 * rints_buff[1];
    rints_buff[5] = xyz_ab[0] * rints_buff[1];
    rints_buff[6] = xyz_ab[1] * rints_buff[1];
    rints_buff[7] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[1] = 1.0 * fac * rints_buff[5];
    rints_out[2] = 1.0 * fac * rints_buff[6];
    rints_out[3] = 1.0 * fac * rints_buff[7];
    rints_out[4] = 1.0 * fac * rints_buff[8];
    rints_out[5] = 1.0 * fac * rints_buff[9];
    rints_out[6] = 1.0 * fac * rints_buff[10];
    rints_out[7] = 1.0 * fac * rints_buff[11];
    rints_out[8] = 1.0 * fac * rints_buff[12];
    rints_out[9] = 1.0 * fac * rints_buff[13];
    rints_out[10] = 1.0 * fac * rints_buff[14];
    rints_out[11] = 1.0 * fac * rints_buff[15];
    rints_out[12] = 1.0 * fac * rints_buff[16];
    rints_out[13] = 1.0 * fac * rints_buff[17];
    rints_out[14] = 1.0 * fac * rints_buff[18];
    rints_out[15] = 1.0 * fac * rints_buff[19];
    rints_out[16] = 1.0 * fac * rints_buff[20];
    rints_out[17] = 1.0 * fac * rints_buff[21];
    rints_out[18] = 1.0 * fac * rints_buff[22];
    rints_out[19] = 1.0 * fac * rints_buff[23];
    rints_out[20] = 1.0 * fac * rints_buff[24];
    rints_out[21] = 1.0 * fac * rints_buff[25];
    rints_out[22] = 1.0 * fac * rints_buff[26];
    rints_out[23] = 1.0 * fac * rints_buff[27];
    rints_out[24] = 1.0 * fac * rints_buff[28];
    rints_out[25] = 1.0 * fac * rints_buff[29];
    rints_out[26] = 1.0 * fac * rints_buff[30];
    rints_out[27] = 1.0 * fac * rints_buff[31];
    rints_out[28] = 1.0 * fac * rints_buff[32];
    rints_out[29] = 1.0 * fac * rints_buff[33];
    rints_out[30] = 1.0 * fac * rints_buff[34];
    rints_out[31] = 1.0 * fac * rints_buff[35];
    rints_out[32] = 1.0 * fac * rints_buff[36];
    rints_out[33] = 1.0 * fac * rints_buff[37];
    rints_out[34] = 1.0 * fac * rints_buff[38];
}
}

template void lible::ints::calcRInts_ERI2D1<4, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI2D2<4, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI3D1<4, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

