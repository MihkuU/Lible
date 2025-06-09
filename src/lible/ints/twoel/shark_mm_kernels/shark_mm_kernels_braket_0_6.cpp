#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra2<0, 0, 0, 6>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 1, 5>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
    E_x_R_x_ET[13] += E[0] * R_x_ET[13];
    E_x_R_x_ET[14] += E[0] * R_x_ET[14];
    E_x_R_x_ET[15] += E[0] * R_x_ET[15];
    E_x_R_x_ET[16] += E[0] * R_x_ET[16];
    E_x_R_x_ET[17] += E[0] * R_x_ET[17];
    E_x_R_x_ET[18] += E[0] * R_x_ET[18];
    E_x_R_x_ET[19] += E[0] * R_x_ET[19];
    E_x_R_x_ET[20] += E[0] * R_x_ET[20];
    E_x_R_x_ET[21] += E[0] * R_x_ET[21];
    E_x_R_x_ET[22] += E[0] * R_x_ET[22];
    E_x_R_x_ET[23] += E[0] * R_x_ET[23];
    E_x_R_x_ET[24] += E[0] * R_x_ET[24];
    E_x_R_x_ET[25] += E[0] * R_x_ET[25];
    E_x_R_x_ET[26] += E[0] * R_x_ET[26];
    E_x_R_x_ET[27] += E[0] * R_x_ET[27];
    E_x_R_x_ET[28] += E[0] * R_x_ET[28];
    E_x_R_x_ET[29] += E[0] * R_x_ET[29];
    E_x_R_x_ET[30] += E[0] * R_x_ET[30];
    E_x_R_x_ET[31] += E[0] * R_x_ET[31];
    E_x_R_x_ET[32] += E[0] * R_x_ET[32];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 2, 4>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
    E_x_R_x_ET[13] += E[0] * R_x_ET[13];
    E_x_R_x_ET[14] += E[0] * R_x_ET[14];
    E_x_R_x_ET[15] += E[0] * R_x_ET[15];
    E_x_R_x_ET[16] += E[0] * R_x_ET[16];
    E_x_R_x_ET[17] += E[0] * R_x_ET[17];
    E_x_R_x_ET[18] += E[0] * R_x_ET[18];
    E_x_R_x_ET[19] += E[0] * R_x_ET[19];
    E_x_R_x_ET[20] += E[0] * R_x_ET[20];
    E_x_R_x_ET[21] += E[0] * R_x_ET[21];
    E_x_R_x_ET[22] += E[0] * R_x_ET[22];
    E_x_R_x_ET[23] += E[0] * R_x_ET[23];
    E_x_R_x_ET[24] += E[0] * R_x_ET[24];
    E_x_R_x_ET[25] += E[0] * R_x_ET[25];
    E_x_R_x_ET[26] += E[0] * R_x_ET[26];
    E_x_R_x_ET[27] += E[0] * R_x_ET[27];
    E_x_R_x_ET[28] += E[0] * R_x_ET[28];
    E_x_R_x_ET[29] += E[0] * R_x_ET[29];
    E_x_R_x_ET[30] += E[0] * R_x_ET[30];
    E_x_R_x_ET[31] += E[0] * R_x_ET[31];
    E_x_R_x_ET[32] += E[0] * R_x_ET[32];
    E_x_R_x_ET[33] += E[0] * R_x_ET[33];
    E_x_R_x_ET[34] += E[0] * R_x_ET[34];
    E_x_R_x_ET[35] += E[0] * R_x_ET[35];
    E_x_R_x_ET[36] += E[0] * R_x_ET[36];
    E_x_R_x_ET[37] += E[0] * R_x_ET[37];
    E_x_R_x_ET[38] += E[0] * R_x_ET[38];
    E_x_R_x_ET[39] += E[0] * R_x_ET[39];
    E_x_R_x_ET[40] += E[0] * R_x_ET[40];
    E_x_R_x_ET[41] += E[0] * R_x_ET[41];
    E_x_R_x_ET[42] += E[0] * R_x_ET[42];
    E_x_R_x_ET[43] += E[0] * R_x_ET[43];
    E_x_R_x_ET[44] += E[0] * R_x_ET[44];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 3, 3>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
    E_x_R_x_ET[13] += E[0] * R_x_ET[13];
    E_x_R_x_ET[14] += E[0] * R_x_ET[14];
    E_x_R_x_ET[15] += E[0] * R_x_ET[15];
    E_x_R_x_ET[16] += E[0] * R_x_ET[16];
    E_x_R_x_ET[17] += E[0] * R_x_ET[17];
    E_x_R_x_ET[18] += E[0] * R_x_ET[18];
    E_x_R_x_ET[19] += E[0] * R_x_ET[19];
    E_x_R_x_ET[20] += E[0] * R_x_ET[20];
    E_x_R_x_ET[21] += E[0] * R_x_ET[21];
    E_x_R_x_ET[22] += E[0] * R_x_ET[22];
    E_x_R_x_ET[23] += E[0] * R_x_ET[23];
    E_x_R_x_ET[24] += E[0] * R_x_ET[24];
    E_x_R_x_ET[25] += E[0] * R_x_ET[25];
    E_x_R_x_ET[26] += E[0] * R_x_ET[26];
    E_x_R_x_ET[27] += E[0] * R_x_ET[27];
    E_x_R_x_ET[28] += E[0] * R_x_ET[28];
    E_x_R_x_ET[29] += E[0] * R_x_ET[29];
    E_x_R_x_ET[30] += E[0] * R_x_ET[30];
    E_x_R_x_ET[31] += E[0] * R_x_ET[31];
    E_x_R_x_ET[32] += E[0] * R_x_ET[32];
    E_x_R_x_ET[33] += E[0] * R_x_ET[33];
    E_x_R_x_ET[34] += E[0] * R_x_ET[34];
    E_x_R_x_ET[35] += E[0] * R_x_ET[35];
    E_x_R_x_ET[36] += E[0] * R_x_ET[36];
    E_x_R_x_ET[37] += E[0] * R_x_ET[37];
    E_x_R_x_ET[38] += E[0] * R_x_ET[38];
    E_x_R_x_ET[39] += E[0] * R_x_ET[39];
    E_x_R_x_ET[40] += E[0] * R_x_ET[40];
    E_x_R_x_ET[41] += E[0] * R_x_ET[41];
    E_x_R_x_ET[42] += E[0] * R_x_ET[42];
    E_x_R_x_ET[43] += E[0] * R_x_ET[43];
    E_x_R_x_ET[44] += E[0] * R_x_ET[44];
    E_x_R_x_ET[45] += E[0] * R_x_ET[45];
    E_x_R_x_ET[46] += E[0] * R_x_ET[46];
    E_x_R_x_ET[47] += E[0] * R_x_ET[47];
    E_x_R_x_ET[48] += E[0] * R_x_ET[48];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 4, 2>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
    E_x_R_x_ET[13] += E[0] * R_x_ET[13];
    E_x_R_x_ET[14] += E[0] * R_x_ET[14];
    E_x_R_x_ET[15] += E[0] * R_x_ET[15];
    E_x_R_x_ET[16] += E[0] * R_x_ET[16];
    E_x_R_x_ET[17] += E[0] * R_x_ET[17];
    E_x_R_x_ET[18] += E[0] * R_x_ET[18];
    E_x_R_x_ET[19] += E[0] * R_x_ET[19];
    E_x_R_x_ET[20] += E[0] * R_x_ET[20];
    E_x_R_x_ET[21] += E[0] * R_x_ET[21];
    E_x_R_x_ET[22] += E[0] * R_x_ET[22];
    E_x_R_x_ET[23] += E[0] * R_x_ET[23];
    E_x_R_x_ET[24] += E[0] * R_x_ET[24];
    E_x_R_x_ET[25] += E[0] * R_x_ET[25];
    E_x_R_x_ET[26] += E[0] * R_x_ET[26];
    E_x_R_x_ET[27] += E[0] * R_x_ET[27];
    E_x_R_x_ET[28] += E[0] * R_x_ET[28];
    E_x_R_x_ET[29] += E[0] * R_x_ET[29];
    E_x_R_x_ET[30] += E[0] * R_x_ET[30];
    E_x_R_x_ET[31] += E[0] * R_x_ET[31];
    E_x_R_x_ET[32] += E[0] * R_x_ET[32];
    E_x_R_x_ET[33] += E[0] * R_x_ET[33];
    E_x_R_x_ET[34] += E[0] * R_x_ET[34];
    E_x_R_x_ET[35] += E[0] * R_x_ET[35];
    E_x_R_x_ET[36] += E[0] * R_x_ET[36];
    E_x_R_x_ET[37] += E[0] * R_x_ET[37];
    E_x_R_x_ET[38] += E[0] * R_x_ET[38];
    E_x_R_x_ET[39] += E[0] * R_x_ET[39];
    E_x_R_x_ET[40] += E[0] * R_x_ET[40];
    E_x_R_x_ET[41] += E[0] * R_x_ET[41];
    E_x_R_x_ET[42] += E[0] * R_x_ET[42];
    E_x_R_x_ET[43] += E[0] * R_x_ET[43];
    E_x_R_x_ET[44] += E[0] * R_x_ET[44];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 5, 1>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
    E_x_R_x_ET[13] += E[0] * R_x_ET[13];
    E_x_R_x_ET[14] += E[0] * R_x_ET[14];
    E_x_R_x_ET[15] += E[0] * R_x_ET[15];
    E_x_R_x_ET[16] += E[0] * R_x_ET[16];
    E_x_R_x_ET[17] += E[0] * R_x_ET[17];
    E_x_R_x_ET[18] += E[0] * R_x_ET[18];
    E_x_R_x_ET[19] += E[0] * R_x_ET[19];
    E_x_R_x_ET[20] += E[0] * R_x_ET[20];
    E_x_R_x_ET[21] += E[0] * R_x_ET[21];
    E_x_R_x_ET[22] += E[0] * R_x_ET[22];
    E_x_R_x_ET[23] += E[0] * R_x_ET[23];
    E_x_R_x_ET[24] += E[0] * R_x_ET[24];
    E_x_R_x_ET[25] += E[0] * R_x_ET[25];
    E_x_R_x_ET[26] += E[0] * R_x_ET[26];
    E_x_R_x_ET[27] += E[0] * R_x_ET[27];
    E_x_R_x_ET[28] += E[0] * R_x_ET[28];
    E_x_R_x_ET[29] += E[0] * R_x_ET[29];
    E_x_R_x_ET[30] += E[0] * R_x_ET[30];
    E_x_R_x_ET[31] += E[0] * R_x_ET[31];
    E_x_R_x_ET[32] += E[0] * R_x_ET[32];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 6, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[7] += E[0] * R_x_ET[7];
    E_x_R_x_ET[8] += E[0] * R_x_ET[8];
    E_x_R_x_ET[9] += E[0] * R_x_ET[9];
    E_x_R_x_ET[10] += E[0] * R_x_ET[10];
    E_x_R_x_ET[11] += E[0] * R_x_ET[11];
    E_x_R_x_ET[12] += E[0] * R_x_ET[12];
}

