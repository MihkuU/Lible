#include <lible/ints/utils.hpp>

namespace lible::ints
{
template <int la, int lb>
void calcRInts(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
               double *rints_out);

template<>
void calcRInts<5, 1>(const double alpha, const double fac, const double *fnx,
                     const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 6;
    constexpr int buff_size = 90;
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
    rints_buff[7] = xyz_ab[0] * rints_buff[6];
    rints_buff[8] = xyz_ab[1] * rints_buff[6];
    rints_buff[9] = xyz_ab[2] * rints_buff[6];
    rints_buff[10] = xyz_ab[0] * rints_buff[7] + 1 * rints_buff[5];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[0] * rints_buff[9];
    rints_buff[13] = xyz_ab[1] * rints_buff[8] + 1 * rints_buff[5];
    rints_buff[14] = xyz_ab[1] * rints_buff[9];
    rints_buff[15] = xyz_ab[2] * rints_buff[9] + 1 * rints_buff[5];
    rints_buff[7] = xyz_ab[0] * rints_buff[5];
    rints_buff[8] = xyz_ab[1] * rints_buff[5];
    rints_buff[9] = xyz_ab[2] * rints_buff[5];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 2 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12] + 1 * rints_buff[9];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[0] * rints_buff[15];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 2 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14] + 1 * rints_buff[9];
    rints_buff[24] = xyz_ab[1] * rints_buff[15];
    rints_buff[25] = xyz_ab[2] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[10] = xyz_ab[0] * rints_buff[7] + 1 * rints_buff[4];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[0] * rints_buff[9];
    rints_buff[13] = xyz_ab[1] * rints_buff[8] + 1 * rints_buff[4];
    rints_buff[14] = xyz_ab[1] * rints_buff[9];
    rints_buff[15] = xyz_ab[2] * rints_buff[9] + 1 * rints_buff[4];
    rints_buff[7] = xyz_ab[0] * rints_buff[4];
    rints_buff[8] = xyz_ab[1] * rints_buff[4];
    rints_buff[9] = xyz_ab[2] * rints_buff[4];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 3 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 2 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 2 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20] + 1 * rints_buff[14];
    rints_buff[31] = xyz_ab[0] * rints_buff[21] + 1 * rints_buff[15];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[0] * rints_buff[24];
    rints_buff[35] = xyz_ab[0] * rints_buff[25];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 3 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23] + 2 * rints_buff[14];
    rints_buff[38] = xyz_ab[1] * rints_buff[24] + 1 * rints_buff[15];
    rints_buff[39] = xyz_ab[1] * rints_buff[25];
    rints_buff[40] = xyz_ab[2] * rints_buff[25] + 3 * rints_buff[15];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 2 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12] + 1 * rints_buff[9];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[0] * rints_buff[15];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 2 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14] + 1 * rints_buff[9];
    rints_buff[24] = xyz_ab[1] * rints_buff[15];
    rints_buff[25] = xyz_ab[2] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[10] = xyz_ab[0] * rints_buff[7] + 1 * rints_buff[3];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[0] * rints_buff[9];
    rints_buff[13] = xyz_ab[1] * rints_buff[8] + 1 * rints_buff[3];
    rints_buff[14] = xyz_ab[1] * rints_buff[9];
    rints_buff[15] = xyz_ab[2] * rints_buff[9] + 1 * rints_buff[3];
    rints_buff[7] = xyz_ab[0] * rints_buff[3];
    rints_buff[8] = xyz_ab[1] * rints_buff[3];
    rints_buff[9] = xyz_ab[2] * rints_buff[3];
    rints_buff[41] = xyz_ab[0] * rints_buff[26] + 4 * rints_buff[16];
    rints_buff[42] = xyz_ab[0] * rints_buff[27] + 3 * rints_buff[17];
    rints_buff[43] = xyz_ab[0] * rints_buff[28] + 3 * rints_buff[18];
    rints_buff[44] = xyz_ab[0] * rints_buff[29] + 2 * rints_buff[19];
    rints_buff[45] = xyz_ab[0] * rints_buff[30] + 2 * rints_buff[20];
    rints_buff[46] = xyz_ab[0] * rints_buff[31] + 2 * rints_buff[21];
    rints_buff[47] = xyz_ab[0] * rints_buff[32] + 1 * rints_buff[22];
    rints_buff[48] = xyz_ab[0] * rints_buff[33] + 1 * rints_buff[23];
    rints_buff[49] = xyz_ab[0] * rints_buff[34] + 1 * rints_buff[24];
    rints_buff[50] = xyz_ab[0] * rints_buff[35] + 1 * rints_buff[25];
    rints_buff[51] = xyz_ab[0] * rints_buff[36];
    rints_buff[52] = xyz_ab[0] * rints_buff[37];
    rints_buff[53] = xyz_ab[0] * rints_buff[38];
    rints_buff[54] = xyz_ab[0] * rints_buff[39];
    rints_buff[55] = xyz_ab[0] * rints_buff[40];
    rints_buff[56] = xyz_ab[1] * rints_buff[36] + 4 * rints_buff[22];
    rints_buff[57] = xyz_ab[1] * rints_buff[37] + 3 * rints_buff[23];
    rints_buff[58] = xyz_ab[1] * rints_buff[38] + 2 * rints_buff[24];
    rints_buff[59] = xyz_ab[1] * rints_buff[39] + 1 * rints_buff[25];
    rints_buff[60] = xyz_ab[1] * rints_buff[40];
    rints_buff[61] = xyz_ab[2] * rints_buff[40] + 4 * rints_buff[25];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 3 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 2 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 2 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20] + 1 * rints_buff[14];
    rints_buff[31] = xyz_ab[0] * rints_buff[21] + 1 * rints_buff[15];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[0] * rints_buff[24];
    rints_buff[35] = xyz_ab[0] * rints_buff[25];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 3 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23] + 2 * rints_buff[14];
    rints_buff[38] = xyz_ab[1] * rints_buff[24] + 1 * rints_buff[15];
    rints_buff[39] = xyz_ab[1] * rints_buff[25];
    rints_buff[40] = xyz_ab[2] * rints_buff[25] + 3 * rints_buff[15];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 2 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12] + 1 * rints_buff[9];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[0] * rints_buff[15];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 2 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14] + 1 * rints_buff[9];
    rints_buff[24] = xyz_ab[1] * rints_buff[15];
    rints_buff[25] = xyz_ab[2] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[10] = xyz_ab[0] * rints_buff[7] + 1 * rints_buff[2];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[0] * rints_buff[9];
    rints_buff[13] = xyz_ab[1] * rints_buff[8] + 1 * rints_buff[2];
    rints_buff[14] = xyz_ab[1] * rints_buff[9];
    rints_buff[15] = xyz_ab[2] * rints_buff[9] + 1 * rints_buff[2];
    rints_buff[7] = xyz_ab[0] * rints_buff[2];
    rints_buff[8] = xyz_ab[1] * rints_buff[2];
    rints_buff[9] = xyz_ab[2] * rints_buff[2];
    rints_buff[62] = xyz_ab[0] * rints_buff[41] + 5 * rints_buff[26];
    rints_buff[63] = xyz_ab[0] * rints_buff[42] + 4 * rints_buff[27];
    rints_buff[64] = xyz_ab[0] * rints_buff[43] + 4 * rints_buff[28];
    rints_buff[65] = xyz_ab[0] * rints_buff[44] + 3 * rints_buff[29];
    rints_buff[66] = xyz_ab[0] * rints_buff[45] + 3 * rints_buff[30];
    rints_buff[67] = xyz_ab[0] * rints_buff[46] + 3 * rints_buff[31];
    rints_buff[68] = xyz_ab[0] * rints_buff[47] + 2 * rints_buff[32];
    rints_buff[69] = xyz_ab[0] * rints_buff[48] + 2 * rints_buff[33];
    rints_buff[70] = xyz_ab[0] * rints_buff[49] + 2 * rints_buff[34];
    rints_buff[71] = xyz_ab[0] * rints_buff[50] + 2 * rints_buff[35];
    rints_buff[72] = xyz_ab[0] * rints_buff[51] + 1 * rints_buff[36];
    rints_buff[73] = xyz_ab[0] * rints_buff[52] + 1 * rints_buff[37];
    rints_buff[74] = xyz_ab[0] * rints_buff[53] + 1 * rints_buff[38];
    rints_buff[75] = xyz_ab[0] * rints_buff[54] + 1 * rints_buff[39];
    rints_buff[76] = xyz_ab[0] * rints_buff[55] + 1 * rints_buff[40];
    rints_buff[77] = xyz_ab[0] * rints_buff[56];
    rints_buff[78] = xyz_ab[0] * rints_buff[57];
    rints_buff[79] = xyz_ab[0] * rints_buff[58];
    rints_buff[80] = xyz_ab[0] * rints_buff[59];
    rints_buff[81] = xyz_ab[0] * rints_buff[60];
    rints_buff[82] = xyz_ab[0] * rints_buff[61];
    rints_buff[83] = xyz_ab[1] * rints_buff[56] + 5 * rints_buff[36];
    rints_buff[84] = xyz_ab[1] * rints_buff[57] + 4 * rints_buff[37];
    rints_buff[85] = xyz_ab[1] * rints_buff[58] + 3 * rints_buff[38];
    rints_buff[86] = xyz_ab[1] * rints_buff[59] + 2 * rints_buff[39];
    rints_buff[87] = xyz_ab[1] * rints_buff[60] + 1 * rints_buff[40];
    rints_buff[88] = xyz_ab[1] * rints_buff[61];
    rints_buff[89] = xyz_ab[2] * rints_buff[61] + 5 * rints_buff[40];
    rints_buff[41] = xyz_ab[0] * rints_buff[26] + 4 * rints_buff[16];
    rints_buff[42] = xyz_ab[0] * rints_buff[27] + 3 * rints_buff[17];
    rints_buff[43] = xyz_ab[0] * rints_buff[28] + 3 * rints_buff[18];
    rints_buff[44] = xyz_ab[0] * rints_buff[29] + 2 * rints_buff[19];
    rints_buff[45] = xyz_ab[0] * rints_buff[30] + 2 * rints_buff[20];
    rints_buff[46] = xyz_ab[0] * rints_buff[31] + 2 * rints_buff[21];
    rints_buff[47] = xyz_ab[0] * rints_buff[32] + 1 * rints_buff[22];
    rints_buff[48] = xyz_ab[0] * rints_buff[33] + 1 * rints_buff[23];
    rints_buff[49] = xyz_ab[0] * rints_buff[34] + 1 * rints_buff[24];
    rints_buff[50] = xyz_ab[0] * rints_buff[35] + 1 * rints_buff[25];
    rints_buff[51] = xyz_ab[0] * rints_buff[36];
    rints_buff[52] = xyz_ab[0] * rints_buff[37];
    rints_buff[53] = xyz_ab[0] * rints_buff[38];
    rints_buff[54] = xyz_ab[0] * rints_buff[39];
    rints_buff[55] = xyz_ab[0] * rints_buff[40];
    rints_buff[56] = xyz_ab[1] * rints_buff[36] + 4 * rints_buff[22];
    rints_buff[57] = xyz_ab[1] * rints_buff[37] + 3 * rints_buff[23];
    rints_buff[58] = xyz_ab[1] * rints_buff[38] + 2 * rints_buff[24];
    rints_buff[59] = xyz_ab[1] * rints_buff[39] + 1 * rints_buff[25];
    rints_buff[60] = xyz_ab[1] * rints_buff[40];
    rints_buff[61] = xyz_ab[2] * rints_buff[40] + 4 * rints_buff[25];
    rints_buff[26] = xyz_ab[0] * rints_buff[16] + 3 * rints_buff[10];
    rints_buff[27] = xyz_ab[0] * rints_buff[17] + 2 * rints_buff[11];
    rints_buff[28] = xyz_ab[0] * rints_buff[18] + 2 * rints_buff[12];
    rints_buff[29] = xyz_ab[0] * rints_buff[19] + 1 * rints_buff[13];
    rints_buff[30] = xyz_ab[0] * rints_buff[20] + 1 * rints_buff[14];
    rints_buff[31] = xyz_ab[0] * rints_buff[21] + 1 * rints_buff[15];
    rints_buff[32] = xyz_ab[0] * rints_buff[22];
    rints_buff[33] = xyz_ab[0] * rints_buff[23];
    rints_buff[34] = xyz_ab[0] * rints_buff[24];
    rints_buff[35] = xyz_ab[0] * rints_buff[25];
    rints_buff[36] = xyz_ab[1] * rints_buff[22] + 3 * rints_buff[13];
    rints_buff[37] = xyz_ab[1] * rints_buff[23] + 2 * rints_buff[14];
    rints_buff[38] = xyz_ab[1] * rints_buff[24] + 1 * rints_buff[15];
    rints_buff[39] = xyz_ab[1] * rints_buff[25];
    rints_buff[40] = xyz_ab[2] * rints_buff[25] + 3 * rints_buff[15];
    rints_buff[16] = xyz_ab[0] * rints_buff[10] + 2 * rints_buff[7];
    rints_buff[17] = xyz_ab[0] * rints_buff[11] + 1 * rints_buff[8];
    rints_buff[18] = xyz_ab[0] * rints_buff[12] + 1 * rints_buff[9];
    rints_buff[19] = xyz_ab[0] * rints_buff[13];
    rints_buff[20] = xyz_ab[0] * rints_buff[14];
    rints_buff[21] = xyz_ab[0] * rints_buff[15];
    rints_buff[22] = xyz_ab[1] * rints_buff[13] + 2 * rints_buff[8];
    rints_buff[23] = xyz_ab[1] * rints_buff[14] + 1 * rints_buff[9];
    rints_buff[24] = xyz_ab[1] * rints_buff[15];
    rints_buff[25] = xyz_ab[2] * rints_buff[15] + 2 * rints_buff[9];
    rints_buff[10] = xyz_ab[0] * rints_buff[7] + 1 * rints_buff[1];
    rints_buff[11] = xyz_ab[0] * rints_buff[8];
    rints_buff[12] = xyz_ab[0] * rints_buff[9];
    rints_buff[13] = xyz_ab[1] * rints_buff[8] + 1 * rints_buff[1];
    rints_buff[14] = xyz_ab[1] * rints_buff[9];
    rints_buff[15] = xyz_ab[2] * rints_buff[9] + 1 * rints_buff[1];
    rints_buff[7] = xyz_ab[0] * rints_buff[1];
    rints_buff[8] = xyz_ab[1] * rints_buff[1];
    rints_buff[9] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[4] = 1.0 * fac * rints_buff[7];
    rints_out[8] = 1.0 * fac * rints_buff[8];
    rints_out[12] = 1.0 * fac * rints_buff[9];
    rints_out[16] = 1.0 * fac * rints_buff[10];
    rints_out[20] = 1.0 * fac * rints_buff[11];
    rints_out[24] = 1.0 * fac * rints_buff[12];
    rints_out[28] = 1.0 * fac * rints_buff[13];
    rints_out[32] = 1.0 * fac * rints_buff[14];
    rints_out[36] = 1.0 * fac * rints_buff[15];
    rints_out[40] = 1.0 * fac * rints_buff[16];
    rints_out[44] = 1.0 * fac * rints_buff[17];
    rints_out[48] = 1.0 * fac * rints_buff[18];
    rints_out[52] = 1.0 * fac * rints_buff[19];
    rints_out[56] = 1.0 * fac * rints_buff[20];
    rints_out[60] = 1.0 * fac * rints_buff[21];
    rints_out[64] = 1.0 * fac * rints_buff[22];
    rints_out[68] = 1.0 * fac * rints_buff[23];
    rints_out[72] = 1.0 * fac * rints_buff[24];
    rints_out[76] = 1.0 * fac * rints_buff[25];
    rints_out[80] = 1.0 * fac * rints_buff[26];
    rints_out[84] = 1.0 * fac * rints_buff[27];
    rints_out[88] = 1.0 * fac * rints_buff[28];
    rints_out[92] = 1.0 * fac * rints_buff[29];
    rints_out[96] = 1.0 * fac * rints_buff[30];
    rints_out[100] = 1.0 * fac * rints_buff[31];
    rints_out[104] = 1.0 * fac * rints_buff[32];
    rints_out[108] = 1.0 * fac * rints_buff[33];
    rints_out[112] = 1.0 * fac * rints_buff[34];
    rints_out[116] = 1.0 * fac * rints_buff[35];
    rints_out[120] = 1.0 * fac * rints_buff[36];
    rints_out[124] = 1.0 * fac * rints_buff[37];
    rints_out[128] = 1.0 * fac * rints_buff[38];
    rints_out[132] = 1.0 * fac * rints_buff[39];
    rints_out[136] = 1.0 * fac * rints_buff[40];
    rints_out[140] = 1.0 * fac * rints_buff[41];
    rints_out[144] = 1.0 * fac * rints_buff[42];
    rints_out[148] = 1.0 * fac * rints_buff[43];
    rints_out[152] = 1.0 * fac * rints_buff[44];
    rints_out[156] = 1.0 * fac * rints_buff[45];
    rints_out[160] = 1.0 * fac * rints_buff[46];
    rints_out[164] = 1.0 * fac * rints_buff[47];
    rints_out[168] = 1.0 * fac * rints_buff[48];
    rints_out[172] = 1.0 * fac * rints_buff[49];
    rints_out[176] = 1.0 * fac * rints_buff[50];
    rints_out[180] = 1.0 * fac * rints_buff[51];
    rints_out[184] = 1.0 * fac * rints_buff[52];
    rints_out[188] = 1.0 * fac * rints_buff[53];
    rints_out[192] = 1.0 * fac * rints_buff[54];
    rints_out[196] = 1.0 * fac * rints_buff[55];
    rints_out[200] = 1.0 * fac * rints_buff[56];
    rints_out[204] = 1.0 * fac * rints_buff[57];
    rints_out[208] = 1.0 * fac * rints_buff[58];
    rints_out[212] = 1.0 * fac * rints_buff[59];
    rints_out[216] = 1.0 * fac * rints_buff[60];
    rints_out[220] = 1.0 * fac * rints_buff[61];
    rints_out[1] = -1.0 * fac * rints_buff[7];
    rints_out[5] = -1.0 * fac * rints_buff[10];
    rints_out[9] = -1.0 * fac * rints_buff[11];
    rints_out[13] = -1.0 * fac * rints_buff[12];
    rints_out[17] = -1.0 * fac * rints_buff[16];
    rints_out[21] = -1.0 * fac * rints_buff[17];
    rints_out[25] = -1.0 * fac * rints_buff[18];
    rints_out[29] = -1.0 * fac * rints_buff[19];
    rints_out[33] = -1.0 * fac * rints_buff[20];
    rints_out[37] = -1.0 * fac * rints_buff[21];
    rints_out[41] = -1.0 * fac * rints_buff[26];
    rints_out[45] = -1.0 * fac * rints_buff[27];
    rints_out[49] = -1.0 * fac * rints_buff[28];
    rints_out[53] = -1.0 * fac * rints_buff[29];
    rints_out[57] = -1.0 * fac * rints_buff[30];
    rints_out[61] = -1.0 * fac * rints_buff[31];
    rints_out[65] = -1.0 * fac * rints_buff[32];
    rints_out[69] = -1.0 * fac * rints_buff[33];
    rints_out[73] = -1.0 * fac * rints_buff[34];
    rints_out[77] = -1.0 * fac * rints_buff[35];
    rints_out[81] = -1.0 * fac * rints_buff[41];
    rints_out[85] = -1.0 * fac * rints_buff[42];
    rints_out[89] = -1.0 * fac * rints_buff[43];
    rints_out[93] = -1.0 * fac * rints_buff[44];
    rints_out[97] = -1.0 * fac * rints_buff[45];
    rints_out[101] = -1.0 * fac * rints_buff[46];
    rints_out[105] = -1.0 * fac * rints_buff[47];
    rints_out[109] = -1.0 * fac * rints_buff[48];
    rints_out[113] = -1.0 * fac * rints_buff[49];
    rints_out[117] = -1.0 * fac * rints_buff[50];
    rints_out[121] = -1.0 * fac * rints_buff[51];
    rints_out[125] = -1.0 * fac * rints_buff[52];
    rints_out[129] = -1.0 * fac * rints_buff[53];
    rints_out[133] = -1.0 * fac * rints_buff[54];
    rints_out[137] = -1.0 * fac * rints_buff[55];
    rints_out[141] = -1.0 * fac * rints_buff[62];
    rints_out[145] = -1.0 * fac * rints_buff[63];
    rints_out[149] = -1.0 * fac * rints_buff[64];
    rints_out[153] = -1.0 * fac * rints_buff[65];
    rints_out[157] = -1.0 * fac * rints_buff[66];
    rints_out[161] = -1.0 * fac * rints_buff[67];
    rints_out[165] = -1.0 * fac * rints_buff[68];
    rints_out[169] = -1.0 * fac * rints_buff[69];
    rints_out[173] = -1.0 * fac * rints_buff[70];
    rints_out[177] = -1.0 * fac * rints_buff[71];
    rints_out[181] = -1.0 * fac * rints_buff[72];
    rints_out[185] = -1.0 * fac * rints_buff[73];
    rints_out[189] = -1.0 * fac * rints_buff[74];
    rints_out[193] = -1.0 * fac * rints_buff[75];
    rints_out[197] = -1.0 * fac * rints_buff[76];
    rints_out[201] = -1.0 * fac * rints_buff[77];
    rints_out[205] = -1.0 * fac * rints_buff[78];
    rints_out[209] = -1.0 * fac * rints_buff[79];
    rints_out[213] = -1.0 * fac * rints_buff[80];
    rints_out[217] = -1.0 * fac * rints_buff[81];
    rints_out[221] = -1.0 * fac * rints_buff[82];
    rints_out[2] = -1.0 * fac * rints_buff[8];
    rints_out[6] = -1.0 * fac * rints_buff[11];
    rints_out[10] = -1.0 * fac * rints_buff[13];
    rints_out[14] = -1.0 * fac * rints_buff[14];
    rints_out[18] = -1.0 * fac * rints_buff[17];
    rints_out[22] = -1.0 * fac * rints_buff[19];
    rints_out[26] = -1.0 * fac * rints_buff[20];
    rints_out[30] = -1.0 * fac * rints_buff[22];
    rints_out[34] = -1.0 * fac * rints_buff[23];
    rints_out[38] = -1.0 * fac * rints_buff[24];
    rints_out[42] = -1.0 * fac * rints_buff[27];
    rints_out[46] = -1.0 * fac * rints_buff[29];
    rints_out[50] = -1.0 * fac * rints_buff[30];
    rints_out[54] = -1.0 * fac * rints_buff[32];
    rints_out[58] = -1.0 * fac * rints_buff[33];
    rints_out[62] = -1.0 * fac * rints_buff[34];
    rints_out[66] = -1.0 * fac * rints_buff[36];
    rints_out[70] = -1.0 * fac * rints_buff[37];
    rints_out[74] = -1.0 * fac * rints_buff[38];
    rints_out[78] = -1.0 * fac * rints_buff[39];
    rints_out[82] = -1.0 * fac * rints_buff[42];
    rints_out[86] = -1.0 * fac * rints_buff[44];
    rints_out[90] = -1.0 * fac * rints_buff[45];
    rints_out[94] = -1.0 * fac * rints_buff[47];
    rints_out[98] = -1.0 * fac * rints_buff[48];
    rints_out[102] = -1.0 * fac * rints_buff[49];
    rints_out[106] = -1.0 * fac * rints_buff[51];
    rints_out[110] = -1.0 * fac * rints_buff[52];
    rints_out[114] = -1.0 * fac * rints_buff[53];
    rints_out[118] = -1.0 * fac * rints_buff[54];
    rints_out[122] = -1.0 * fac * rints_buff[56];
    rints_out[126] = -1.0 * fac * rints_buff[57];
    rints_out[130] = -1.0 * fac * rints_buff[58];
    rints_out[134] = -1.0 * fac * rints_buff[59];
    rints_out[138] = -1.0 * fac * rints_buff[60];
    rints_out[142] = -1.0 * fac * rints_buff[63];
    rints_out[146] = -1.0 * fac * rints_buff[65];
    rints_out[150] = -1.0 * fac * rints_buff[66];
    rints_out[154] = -1.0 * fac * rints_buff[68];
    rints_out[158] = -1.0 * fac * rints_buff[69];
    rints_out[162] = -1.0 * fac * rints_buff[70];
    rints_out[166] = -1.0 * fac * rints_buff[72];
    rints_out[170] = -1.0 * fac * rints_buff[73];
    rints_out[174] = -1.0 * fac * rints_buff[74];
    rints_out[178] = -1.0 * fac * rints_buff[75];
    rints_out[182] = -1.0 * fac * rints_buff[77];
    rints_out[186] = -1.0 * fac * rints_buff[78];
    rints_out[190] = -1.0 * fac * rints_buff[79];
    rints_out[194] = -1.0 * fac * rints_buff[80];
    rints_out[198] = -1.0 * fac * rints_buff[81];
    rints_out[202] = -1.0 * fac * rints_buff[83];
    rints_out[206] = -1.0 * fac * rints_buff[84];
    rints_out[210] = -1.0 * fac * rints_buff[85];
    rints_out[214] = -1.0 * fac * rints_buff[86];
    rints_out[218] = -1.0 * fac * rints_buff[87];
    rints_out[222] = -1.0 * fac * rints_buff[88];
    rints_out[3] = -1.0 * fac * rints_buff[9];
    rints_out[7] = -1.0 * fac * rints_buff[12];
    rints_out[11] = -1.0 * fac * rints_buff[14];
    rints_out[15] = -1.0 * fac * rints_buff[15];
    rints_out[19] = -1.0 * fac * rints_buff[18];
    rints_out[23] = -1.0 * fac * rints_buff[20];
    rints_out[27] = -1.0 * fac * rints_buff[21];
    rints_out[31] = -1.0 * fac * rints_buff[23];
    rints_out[35] = -1.0 * fac * rints_buff[24];
    rints_out[39] = -1.0 * fac * rints_buff[25];
    rints_out[43] = -1.0 * fac * rints_buff[28];
    rints_out[47] = -1.0 * fac * rints_buff[30];
    rints_out[51] = -1.0 * fac * rints_buff[31];
    rints_out[55] = -1.0 * fac * rints_buff[33];
    rints_out[59] = -1.0 * fac * rints_buff[34];
    rints_out[63] = -1.0 * fac * rints_buff[35];
    rints_out[67] = -1.0 * fac * rints_buff[37];
    rints_out[71] = -1.0 * fac * rints_buff[38];
    rints_out[75] = -1.0 * fac * rints_buff[39];
    rints_out[79] = -1.0 * fac * rints_buff[40];
    rints_out[83] = -1.0 * fac * rints_buff[43];
    rints_out[87] = -1.0 * fac * rints_buff[45];
    rints_out[91] = -1.0 * fac * rints_buff[46];
    rints_out[95] = -1.0 * fac * rints_buff[48];
    rints_out[99] = -1.0 * fac * rints_buff[49];
    rints_out[103] = -1.0 * fac * rints_buff[50];
    rints_out[107] = -1.0 * fac * rints_buff[52];
    rints_out[111] = -1.0 * fac * rints_buff[53];
    rints_out[115] = -1.0 * fac * rints_buff[54];
    rints_out[119] = -1.0 * fac * rints_buff[55];
    rints_out[123] = -1.0 * fac * rints_buff[57];
    rints_out[127] = -1.0 * fac * rints_buff[58];
    rints_out[131] = -1.0 * fac * rints_buff[59];
    rints_out[135] = -1.0 * fac * rints_buff[60];
    rints_out[139] = -1.0 * fac * rints_buff[61];
    rints_out[143] = -1.0 * fac * rints_buff[64];
    rints_out[147] = -1.0 * fac * rints_buff[66];
    rints_out[151] = -1.0 * fac * rints_buff[67];
    rints_out[155] = -1.0 * fac * rints_buff[69];
    rints_out[159] = -1.0 * fac * rints_buff[70];
    rints_out[163] = -1.0 * fac * rints_buff[71];
    rints_out[167] = -1.0 * fac * rints_buff[73];
    rints_out[171] = -1.0 * fac * rints_buff[74];
    rints_out[175] = -1.0 * fac * rints_buff[75];
    rints_out[179] = -1.0 * fac * rints_buff[76];
    rints_out[183] = -1.0 * fac * rints_buff[78];
    rints_out[187] = -1.0 * fac * rints_buff[79];
    rints_out[191] = -1.0 * fac * rints_buff[80];
    rints_out[195] = -1.0 * fac * rints_buff[81];
    rints_out[199] = -1.0 * fac * rints_buff[82];
    rints_out[203] = -1.0 * fac * rints_buff[84];
    rints_out[207] = -1.0 * fac * rints_buff[85];
    rints_out[211] = -1.0 * fac * rints_buff[86];
    rints_out[215] = -1.0 * fac * rints_buff[87];
    rints_out[219] = -1.0 * fac * rints_buff[88];
    rints_out[223] = -1.0 * fac * rints_buff[89];
}
}