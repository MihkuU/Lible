#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra1<1, 3>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[21];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[22];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[23];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[3] += E[3] * R_x_ET[24];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[4] += E[3] * R_x_ET[25];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[5] += E[3] * R_x_ET[26];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[6] += E[3] * R_x_ET[27];
    E_x_R_x_ET[7] += E[5] * R_x_ET[7];
    E_x_R_x_ET[7] += E[4] * R_x_ET[0];
    E_x_R_x_ET[8] += E[5] * R_x_ET[8];
    E_x_R_x_ET[8] += E[4] * R_x_ET[1];
    E_x_R_x_ET[9] += E[5] * R_x_ET[9];
    E_x_R_x_ET[9] += E[4] * R_x_ET[2];
    E_x_R_x_ET[10] += E[5] * R_x_ET[10];
    E_x_R_x_ET[10] += E[4] * R_x_ET[3];
    E_x_R_x_ET[11] += E[5] * R_x_ET[11];
    E_x_R_x_ET[11] += E[4] * R_x_ET[4];
    E_x_R_x_ET[12] += E[5] * R_x_ET[12];
    E_x_R_x_ET[12] += E[4] * R_x_ET[5];
    E_x_R_x_ET[13] += E[5] * R_x_ET[13];
    E_x_R_x_ET[13] += E[4] * R_x_ET[6];
    E_x_R_x_ET[14] += E[8] * R_x_ET[0];
    E_x_R_x_ET[14] += E[10] * R_x_ET[14];
    E_x_R_x_ET[15] += E[8] * R_x_ET[1];
    E_x_R_x_ET[15] += E[10] * R_x_ET[15];
    E_x_R_x_ET[16] += E[8] * R_x_ET[2];
    E_x_R_x_ET[16] += E[10] * R_x_ET[16];
    E_x_R_x_ET[17] += E[8] * R_x_ET[3];
    E_x_R_x_ET[17] += E[10] * R_x_ET[17];
    E_x_R_x_ET[18] += E[8] * R_x_ET[4];
    E_x_R_x_ET[18] += E[10] * R_x_ET[18];
    E_x_R_x_ET[19] += E[8] * R_x_ET[5];
    E_x_R_x_ET[19] += E[10] * R_x_ET[19];
    E_x_R_x_ET[20] += E[8] * R_x_ET[6];
    E_x_R_x_ET[20] += E[10] * R_x_ET[20];
}

template<> void lible::ints::shark_mm_bra2<0, 1, 3>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[21];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[22];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[23];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[3] += E[3] * R_x_ET[24];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[4] += E[3] * R_x_ET[25];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[5] += E[3] * R_x_ET[26];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[6] += E[3] * R_x_ET[27];
    E_x_R_x_ET[7] += E[5] * R_x_ET[7];
    E_x_R_x_ET[7] += E[4] * R_x_ET[0];
    E_x_R_x_ET[8] += E[5] * R_x_ET[8];
    E_x_R_x_ET[8] += E[4] * R_x_ET[1];
    E_x_R_x_ET[9] += E[5] * R_x_ET[9];
    E_x_R_x_ET[9] += E[4] * R_x_ET[2];
    E_x_R_x_ET[10] += E[5] * R_x_ET[10];
    E_x_R_x_ET[10] += E[4] * R_x_ET[3];
    E_x_R_x_ET[11] += E[5] * R_x_ET[11];
    E_x_R_x_ET[11] += E[4] * R_x_ET[4];
    E_x_R_x_ET[12] += E[5] * R_x_ET[12];
    E_x_R_x_ET[12] += E[4] * R_x_ET[5];
    E_x_R_x_ET[13] += E[5] * R_x_ET[13];
    E_x_R_x_ET[13] += E[4] * R_x_ET[6];
    E_x_R_x_ET[14] += E[8] * R_x_ET[0];
    E_x_R_x_ET[14] += E[10] * R_x_ET[14];
    E_x_R_x_ET[15] += E[8] * R_x_ET[1];
    E_x_R_x_ET[15] += E[10] * R_x_ET[15];
    E_x_R_x_ET[16] += E[8] * R_x_ET[2];
    E_x_R_x_ET[16] += E[10] * R_x_ET[16];
    E_x_R_x_ET[17] += E[8] * R_x_ET[3];
    E_x_R_x_ET[17] += E[10] * R_x_ET[17];
    E_x_R_x_ET[18] += E[8] * R_x_ET[4];
    E_x_R_x_ET[18] += E[10] * R_x_ET[18];
    E_x_R_x_ET[19] += E[8] * R_x_ET[5];
    E_x_R_x_ET[19] += E[10] * R_x_ET[19];
    E_x_R_x_ET[20] += E[8] * R_x_ET[6];
    E_x_R_x_ET[20] += E[10] * R_x_ET[20];
}

template<> void lible::ints::shark_mm_bra2<1, 0, 3>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[21];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[22];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[23];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[3] += E[3] * R_x_ET[24];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
    E_x_R_x_ET[4] += E[3] * R_x_ET[25];
    E_x_R_x_ET[5] += E[0] * R_x_ET[5];
    E_x_R_x_ET[5] += E[3] * R_x_ET[26];
    E_x_R_x_ET[6] += E[0] * R_x_ET[6];
    E_x_R_x_ET[6] += E[3] * R_x_ET[27];
    E_x_R_x_ET[7] += E[5] * R_x_ET[7];
    E_x_R_x_ET[7] += E[4] * R_x_ET[0];
    E_x_R_x_ET[8] += E[5] * R_x_ET[8];
    E_x_R_x_ET[8] += E[4] * R_x_ET[1];
    E_x_R_x_ET[9] += E[5] * R_x_ET[9];
    E_x_R_x_ET[9] += E[4] * R_x_ET[2];
    E_x_R_x_ET[10] += E[5] * R_x_ET[10];
    E_x_R_x_ET[10] += E[4] * R_x_ET[3];
    E_x_R_x_ET[11] += E[5] * R_x_ET[11];
    E_x_R_x_ET[11] += E[4] * R_x_ET[4];
    E_x_R_x_ET[12] += E[5] * R_x_ET[12];
    E_x_R_x_ET[12] += E[4] * R_x_ET[5];
    E_x_R_x_ET[13] += E[5] * R_x_ET[13];
    E_x_R_x_ET[13] += E[4] * R_x_ET[6];
    E_x_R_x_ET[14] += E[8] * R_x_ET[0];
    E_x_R_x_ET[14] += E[10] * R_x_ET[14];
    E_x_R_x_ET[15] += E[8] * R_x_ET[1];
    E_x_R_x_ET[15] += E[10] * R_x_ET[15];
    E_x_R_x_ET[16] += E[8] * R_x_ET[2];
    E_x_R_x_ET[16] += E[10] * R_x_ET[16];
    E_x_R_x_ET[17] += E[8] * R_x_ET[3];
    E_x_R_x_ET[17] += E[10] * R_x_ET[17];
    E_x_R_x_ET[18] += E[8] * R_x_ET[4];
    E_x_R_x_ET[18] += E[10] * R_x_ET[18];
    E_x_R_x_ET[19] += E[8] * R_x_ET[5];
    E_x_R_x_ET[19] += E[10] * R_x_ET[19];
    E_x_R_x_ET[20] += E[8] * R_x_ET[6];
    E_x_R_x_ET[20] += E[10] * R_x_ET[20];
}

