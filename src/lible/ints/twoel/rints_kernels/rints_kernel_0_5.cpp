#include <lible/ints/rints_meta.hpp>

namespace lible::ints
{
template<>
void calcRInts_ERI<0, 5>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 5;
    constexpr int buff_size = 61;
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
    rints_buff[6] = xyz_ab[0] * rints_buff[5];
    rints_buff[7] = xyz_ab[1] * rints_buff[5];
    rints_buff[8] = xyz_ab[2] * rints_buff[5];
    rints_buff[9] = xyz_ab[0] * rints_buff[6] + 1 * rints_buff[4];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[1] * rints_buff[7] + 1 * rints_buff[4];
    rints_buff[13] = xyz_ab[1] * rints_buff[8];
    rints_buff[14] = xyz_ab[2] * rints_buff[8] + 1 * rints_buff[4];
    rints_buff[6] = xyz_ab[0] * rints_buff[4];
    rints_buff[7] = xyz_ab[1] * rints_buff[4];
    rints_buff[8] = xyz_ab[2] * rints_buff[4];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 2 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 2 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 1 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14];
    rints_buff[24] = xyz_ab[2] * rints_buff[14] + 2 * rints_buff[8];
    rints_buff[9] = xyz_ab[0] * rints_buff[6] + 1 * rints_buff[3];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[1] * rints_buff[7] + 1 * rints_buff[3];
    rints_buff[13] = xyz_ab[1] * rints_buff[8];
    rints_buff[14] = xyz_ab[2] * rints_buff[8] + 1 * rints_buff[3];
    rints_buff[6] = xyz_ab[0] * rints_buff[3];
    rints_buff[7] = xyz_ab[1] * rints_buff[3];
    rints_buff[8] = xyz_ab[2] * rints_buff[3];
    rints_buff[25] = xyz_ab[0] * rints_buff[15] + 3 * rints_buff[9];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 2 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 2 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 1 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20] + 1 * rints_buff[14];
    rints_buff[31] = xyz_ab[0] * rints_buff[21];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[0] * rints_buff[24];
    rints_buff[35] = xyz_ab[1] * rints_buff[21] + 3 * rints_buff[12];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 2 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23] + 1 * rints_buff[14];
    rints_buff[38] = xyz_ab[1] * rints_buff[24];
    rints_buff[39] = xyz_ab[2] * rints_buff[24] + 3 * rints_buff[14];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 2 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 2 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 1 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14];
    rints_buff[24] = xyz_ab[2] * rints_buff[14] + 2 * rints_buff[8];
    rints_buff[9] = xyz_ab[0] * rints_buff[6] + 1 * rints_buff[2];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[1] * rints_buff[7] + 1 * rints_buff[2];
    rints_buff[13] = xyz_ab[1] * rints_buff[8];
    rints_buff[14] = xyz_ab[2] * rints_buff[8] + 1 * rints_buff[2];
    rints_buff[6] = xyz_ab[0] * rints_buff[2];
    rints_buff[7] = xyz_ab[1] * rints_buff[2];
    rints_buff[8] = xyz_ab[2] * rints_buff[2];
    rints_buff[40] = xyz_ab[0] * rints_buff[25] + 4 * rints_buff[15];
    rints_buff[41] = xyz_ab[0] * rints_buff[26] + 3 * rints_buff[16];
    rints_buff[42] = xyz_ab[0] * rints_buff[27] + 3 * rints_buff[17];
    rints_buff[43] = xyz_ab[0] * rints_buff[28] + 2 * rints_buff[18];
    rints_buff[44] = xyz_ab[0] * rints_buff[29] + 2 * rints_buff[19];
    rints_buff[45] = xyz_ab[0] * rints_buff[30] + 2 * rints_buff[20];
    rints_buff[46] = xyz_ab[0] * rints_buff[31] + 1 * rints_buff[21];
    rints_buff[47] = xyz_ab[0] * rints_buff[32] + 1 * rints_buff[22];
    rints_buff[48] = xyz_ab[0] * rints_buff[33] + 1 * rints_buff[23];
    rints_buff[49] = xyz_ab[0] * rints_buff[34] + 1 * rints_buff[24];
    rints_buff[50] = xyz_ab[0] * rints_buff[35];
    rints_buff[51] = xyz_ab[0] * rints_buff[36];
    rints_buff[52] = xyz_ab[0] * rints_buff[37];
    rints_buff[53] = xyz_ab[0] * rints_buff[38];
    rints_buff[54] = xyz_ab[0] * rints_buff[39];
    rints_buff[55] = xyz_ab[1] * rints_buff[35] + 4 * rints_buff[21];
    rints_buff[56] = xyz_ab[1] * rints_buff[36] + 3 * rints_buff[22];
    rints_buff[57] = xyz_ab[1] * rints_buff[37] + 2 * rints_buff[23];
    rints_buff[58] = xyz_ab[1] * rints_buff[38] + 1 * rints_buff[24];
    rints_buff[59] = xyz_ab[1] * rints_buff[39];
    rints_buff[60] = xyz_ab[2] * rints_buff[39] + 4 * rints_buff[24];
    rints_buff[25] = xyz_ab[0] * rints_buff[15] + 3 * rints_buff[9];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 2 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 2 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 1 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20] + 1 * rints_buff[14];
    rints_buff[31] = xyz_ab[0] * rints_buff[21];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[0] * rints_buff[24];
    rints_buff[35] = xyz_ab[1] * rints_buff[21] + 3 * rints_buff[12];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 2 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23] + 1 * rints_buff[14];
    rints_buff[38] = xyz_ab[1] * rints_buff[24];
    rints_buff[39] = xyz_ab[2] * rints_buff[24] + 3 * rints_buff[14];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 2 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 1 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[1] * rints_buff[12] + 2 * rints_buff[7];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 1 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14];
    rints_buff[24] = xyz_ab[2] * rints_buff[14] + 2 * rints_buff[8];
    rints_buff[9] = xyz_ab[0] * rints_buff[6] + 1 * rints_buff[1];
    rints_buff[10] = xyz_ab[0] * rints_buff[7];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[1] * rints_buff[7] + 1 * rints_buff[1];
    rints_buff[13] = xyz_ab[1] * rints_buff[8];
    rints_buff[14] = xyz_ab[2] * rints_buff[8] + 1 * rints_buff[1];
    rints_buff[6] = xyz_ab[0] * rints_buff[1];
    rints_buff[7] = xyz_ab[1] * rints_buff[1];
    rints_buff[8] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[1] = -1.0 * fac * rints_buff[6];
    rints_out[2] = -1.0 * fac * rints_buff[7];
    rints_out[3] = -1.0 * fac * rints_buff[8];
    rints_out[4] = 1.0 * fac * rints_buff[9];
    rints_out[5] = 1.0 * fac * rints_buff[10];
    rints_out[6] = 1.0 * fac * rints_buff[11];
    rints_out[7] = 1.0 * fac * rints_buff[12];
    rints_out[8] = 1.0 * fac * rints_buff[13];
    rints_out[9] = 1.0 * fac * rints_buff[14];
    rints_out[10] = -1.0 * fac * rints_buff[15];
    rints_out[11] = -1.0 * fac * rints_buff[16];
    rints_out[12] = -1.0 * fac * rints_buff[17];
    rints_out[13] = -1.0 * fac * rints_buff[18];
    rints_out[14] = -1.0 * fac * rints_buff[19];
    rints_out[15] = -1.0 * fac * rints_buff[20];
    rints_out[16] = -1.0 * fac * rints_buff[21];
    rints_out[17] = -1.0 * fac * rints_buff[22];
    rints_out[18] = -1.0 * fac * rints_buff[23];
    rints_out[19] = -1.0 * fac * rints_buff[24];
    rints_out[20] = 1.0 * fac * rints_buff[25];
    rints_out[21] = 1.0 * fac * rints_buff[26];
    rints_out[22] = 1.0 * fac * rints_buff[27];
    rints_out[23] = 1.0 * fac * rints_buff[28];
    rints_out[24] = 1.0 * fac * rints_buff[29];
    rints_out[25] = 1.0 * fac * rints_buff[30];
    rints_out[26] = 1.0 * fac * rints_buff[31];
    rints_out[27] = 1.0 * fac * rints_buff[32];
    rints_out[28] = 1.0 * fac * rints_buff[33];
    rints_out[29] = 1.0 * fac * rints_buff[34];
    rints_out[30] = 1.0 * fac * rints_buff[35];
    rints_out[31] = 1.0 * fac * rints_buff[36];
    rints_out[32] = 1.0 * fac * rints_buff[37];
    rints_out[33] = 1.0 * fac * rints_buff[38];
    rints_out[34] = 1.0 * fac * rints_buff[39];
    rints_out[35] = -1.0 * fac * rints_buff[40];
    rints_out[36] = -1.0 * fac * rints_buff[41];
    rints_out[37] = -1.0 * fac * rints_buff[42];
    rints_out[38] = -1.0 * fac * rints_buff[43];
    rints_out[39] = -1.0 * fac * rints_buff[44];
    rints_out[40] = -1.0 * fac * rints_buff[45];
    rints_out[41] = -1.0 * fac * rints_buff[46];
    rints_out[42] = -1.0 * fac * rints_buff[47];
    rints_out[43] = -1.0 * fac * rints_buff[48];
    rints_out[44] = -1.0 * fac * rints_buff[49];
    rints_out[45] = -1.0 * fac * rints_buff[50];
    rints_out[46] = -1.0 * fac * rints_buff[51];
    rints_out[47] = -1.0 * fac * rints_buff[52];
    rints_out[48] = -1.0 * fac * rints_buff[53];
    rints_out[49] = -1.0 * fac * rints_buff[54];
    rints_out[50] = -1.0 * fac * rints_buff[55];
    rints_out[51] = -1.0 * fac * rints_buff[56];
    rints_out[52] = -1.0 * fac * rints_buff[57];
    rints_out[53] = -1.0 * fac * rints_buff[58];
    rints_out[54] = -1.0 * fac * rints_buff[59];
    rints_out[55] = -1.0 * fac * rints_buff[60];
}
}

template void lible::ints::calcRInts_ERI2D1<0, 5>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI2D2<0, 5>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI3D1<0, 5>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

template void lible::ints::calcRInts_ERISOC<0, 5>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

