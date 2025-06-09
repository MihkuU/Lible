#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_ket1<0, 2>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[2] * ET[10];
    R_x_ET[0] += R[7] * ET[35];
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[1] * ET[5];
    R_x_ET[0] += R[9] * ET[45];
    R_x_ET[0] += R[4] * ET[20];
    R_x_ET[0] += R[3] * ET[15];
    R_x_ET[1] += R[1] * ET[6];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[1] += R[3] * ET[16];
    R_x_ET[1] += R[6] * ET[31];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[3] * ET[17];
    R_x_ET[2] += R[2] * ET[12];
    R_x_ET[2] += R[8] * ET[42];
    R_x_ET[3] += R[2] * ET[13];
    R_x_ET[3] += R[7] * ET[38];
    R_x_ET[3] += R[0] * ET[3];
    R_x_ET[3] += R[1] * ET[8];
    R_x_ET[3] += R[4] * ET[23];
    R_x_ET[4] += R[1] * ET[9];
    R_x_ET[4] += R[0] * ET[4];
    R_x_ET[4] += R[5] * ET[29];
    R_x_ET[4] += R[2] * ET[14];
}

template<> void lible::ints::shark_mm_ket2<0, 0, 2>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[2] * ET[10];
    R_x_ET[0] += R[7] * ET[35];
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[1] * ET[5];
    R_x_ET[0] += R[9] * ET[45];
    R_x_ET[0] += R[4] * ET[20];
    R_x_ET[0] += R[3] * ET[15];
    R_x_ET[1] += R[1] * ET[6];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[1] += R[3] * ET[16];
    R_x_ET[1] += R[6] * ET[31];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[3] * ET[17];
    R_x_ET[2] += R[2] * ET[12];
    R_x_ET[2] += R[8] * ET[42];
    R_x_ET[3] += R[2] * ET[13];
    R_x_ET[3] += R[7] * ET[38];
    R_x_ET[3] += R[0] * ET[3];
    R_x_ET[3] += R[1] * ET[8];
    R_x_ET[3] += R[4] * ET[23];
    R_x_ET[4] += R[1] * ET[9];
    R_x_ET[4] += R[0] * ET[4];
    R_x_ET[4] += R[5] * ET[29];
    R_x_ET[4] += R[2] * ET[14];
}

template<> void lible::ints::shark_mm_ket2<0, 1, 1>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[3] * ET[27];
    R_x_ET[0] += R[9] * ET[81];
    R_x_ET[1] += R[1] * ET[10];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[1] += R[3] * ET[28];
    R_x_ET[1] += R[6] * ET[55];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[3] * ET[29];
    R_x_ET[2] += R[2] * ET[20];
    R_x_ET[2] += R[8] * ET[74];
    R_x_ET[3] += R[1] * ET[12];
    R_x_ET[3] += R[0] * ET[3];
    R_x_ET[3] += R[3] * ET[30];
    R_x_ET[3] += R[6] * ET[57];
    R_x_ET[4] += R[1] * ET[13];
    R_x_ET[4] += R[0] * ET[4];
    R_x_ET[4] += R[4] * ET[40];
    R_x_ET[5] += R[1] * ET[14];
    R_x_ET[5] += R[0] * ET[5];
    R_x_ET[5] += R[5] * ET[50];
    R_x_ET[5] += R[2] * ET[23];
    R_x_ET[6] += R[0] * ET[6];
    R_x_ET[6] += R[3] * ET[33];
    R_x_ET[6] += R[2] * ET[24];
    R_x_ET[6] += R[8] * ET[78];
    R_x_ET[7] += R[1] * ET[16];
    R_x_ET[7] += R[0] * ET[7];
    R_x_ET[7] += R[5] * ET[52];
    R_x_ET[7] += R[2] * ET[25];
    R_x_ET[8] += R[7] * ET[71];
    R_x_ET[8] += R[0] * ET[8];
    R_x_ET[8] += R[2] * ET[26];
}

template<> void lible::ints::shark_mm_ket2<0, 2, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[2] * ET[10];
    R_x_ET[0] += R[7] * ET[35];
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[1] * ET[5];
    R_x_ET[0] += R[9] * ET[45];
    R_x_ET[0] += R[4] * ET[20];
    R_x_ET[0] += R[3] * ET[15];
    R_x_ET[1] += R[1] * ET[6];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[1] += R[3] * ET[16];
    R_x_ET[1] += R[6] * ET[31];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[3] * ET[17];
    R_x_ET[2] += R[2] * ET[12];
    R_x_ET[2] += R[8] * ET[42];
    R_x_ET[3] += R[2] * ET[13];
    R_x_ET[3] += R[7] * ET[38];
    R_x_ET[3] += R[0] * ET[3];
    R_x_ET[3] += R[1] * ET[8];
    R_x_ET[3] += R[4] * ET[23];
    R_x_ET[4] += R[1] * ET[9];
    R_x_ET[4] += R[0] * ET[4];
    R_x_ET[4] += R[5] * ET[29];
    R_x_ET[4] += R[2] * ET[14];
}

