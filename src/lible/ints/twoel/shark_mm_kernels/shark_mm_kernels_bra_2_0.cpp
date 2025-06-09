#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra1<2, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[2] * R_x_ET[2];
    E_x_R_x_ET[0] += E[7] * R_x_ET[7];
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[1] * R_x_ET[1];
    E_x_R_x_ET[0] += E[9] * R_x_ET[9];
    E_x_R_x_ET[0] += E[4] * R_x_ET[4];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[11] * R_x_ET[1];
    E_x_R_x_ET[1] += E[10] * R_x_ET[0];
    E_x_R_x_ET[1] += E[13] * R_x_ET[3];
    E_x_R_x_ET[1] += E[16] * R_x_ET[6];
    E_x_R_x_ET[2] += E[20] * R_x_ET[0];
    E_x_R_x_ET[2] += E[23] * R_x_ET[3];
    E_x_R_x_ET[2] += E[22] * R_x_ET[2];
    E_x_R_x_ET[2] += E[28] * R_x_ET[8];
    E_x_R_x_ET[3] += E[32] * R_x_ET[2];
    E_x_R_x_ET[3] += E[37] * R_x_ET[7];
    E_x_R_x_ET[3] += E[30] * R_x_ET[0];
    E_x_R_x_ET[3] += E[31] * R_x_ET[1];
    E_x_R_x_ET[3] += E[34] * R_x_ET[4];
    E_x_R_x_ET[4] += E[41] * R_x_ET[1];
    E_x_R_x_ET[4] += E[40] * R_x_ET[0];
    E_x_R_x_ET[4] += E[45] * R_x_ET[5];
    E_x_R_x_ET[4] += E[42] * R_x_ET[2];
}

template<> void lible::ints::shark_mm_bra2<0, 2, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[2] * R_x_ET[2];
    E_x_R_x_ET[0] += E[7] * R_x_ET[7];
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[1] * R_x_ET[1];
    E_x_R_x_ET[0] += E[9] * R_x_ET[9];
    E_x_R_x_ET[0] += E[4] * R_x_ET[4];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[11] * R_x_ET[1];
    E_x_R_x_ET[1] += E[10] * R_x_ET[0];
    E_x_R_x_ET[1] += E[13] * R_x_ET[3];
    E_x_R_x_ET[1] += E[16] * R_x_ET[6];
    E_x_R_x_ET[2] += E[20] * R_x_ET[0];
    E_x_R_x_ET[2] += E[23] * R_x_ET[3];
    E_x_R_x_ET[2] += E[22] * R_x_ET[2];
    E_x_R_x_ET[2] += E[28] * R_x_ET[8];
    E_x_R_x_ET[3] += E[32] * R_x_ET[2];
    E_x_R_x_ET[3] += E[37] * R_x_ET[7];
    E_x_R_x_ET[3] += E[30] * R_x_ET[0];
    E_x_R_x_ET[3] += E[31] * R_x_ET[1];
    E_x_R_x_ET[3] += E[34] * R_x_ET[4];
    E_x_R_x_ET[4] += E[41] * R_x_ET[1];
    E_x_R_x_ET[4] += E[40] * R_x_ET[0];
    E_x_R_x_ET[4] += E[45] * R_x_ET[5];
    E_x_R_x_ET[4] += E[42] * R_x_ET[2];
}

template<> void lible::ints::shark_mm_bra2<1, 1, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[0] += E[9] * R_x_ET[9];
    E_x_R_x_ET[1] += E[11] * R_x_ET[1];
    E_x_R_x_ET[1] += E[10] * R_x_ET[0];
    E_x_R_x_ET[1] += E[13] * R_x_ET[3];
    E_x_R_x_ET[1] += E[16] * R_x_ET[6];
    E_x_R_x_ET[2] += E[20] * R_x_ET[0];
    E_x_R_x_ET[2] += E[23] * R_x_ET[3];
    E_x_R_x_ET[2] += E[22] * R_x_ET[2];
    E_x_R_x_ET[2] += E[28] * R_x_ET[8];
    E_x_R_x_ET[3] += E[31] * R_x_ET[1];
    E_x_R_x_ET[3] += E[30] * R_x_ET[0];
    E_x_R_x_ET[3] += E[33] * R_x_ET[3];
    E_x_R_x_ET[3] += E[36] * R_x_ET[6];
    E_x_R_x_ET[4] += E[41] * R_x_ET[1];
    E_x_R_x_ET[4] += E[40] * R_x_ET[0];
    E_x_R_x_ET[4] += E[44] * R_x_ET[4];
    E_x_R_x_ET[5] += E[51] * R_x_ET[1];
    E_x_R_x_ET[5] += E[50] * R_x_ET[0];
    E_x_R_x_ET[5] += E[55] * R_x_ET[5];
    E_x_R_x_ET[5] += E[52] * R_x_ET[2];
    E_x_R_x_ET[6] += E[60] * R_x_ET[0];
    E_x_R_x_ET[6] += E[63] * R_x_ET[3];
    E_x_R_x_ET[6] += E[62] * R_x_ET[2];
    E_x_R_x_ET[6] += E[68] * R_x_ET[8];
    E_x_R_x_ET[7] += E[71] * R_x_ET[1];
    E_x_R_x_ET[7] += E[70] * R_x_ET[0];
    E_x_R_x_ET[7] += E[75] * R_x_ET[5];
    E_x_R_x_ET[7] += E[72] * R_x_ET[2];
    E_x_R_x_ET[8] += E[87] * R_x_ET[7];
    E_x_R_x_ET[8] += E[80] * R_x_ET[0];
    E_x_R_x_ET[8] += E[82] * R_x_ET[2];
}

template<> void lible::ints::shark_mm_bra2<2, 0, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[2] * R_x_ET[2];
    E_x_R_x_ET[0] += E[7] * R_x_ET[7];
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[1] * R_x_ET[1];
    E_x_R_x_ET[0] += E[9] * R_x_ET[9];
    E_x_R_x_ET[0] += E[4] * R_x_ET[4];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[11] * R_x_ET[1];
    E_x_R_x_ET[1] += E[10] * R_x_ET[0];
    E_x_R_x_ET[1] += E[13] * R_x_ET[3];
    E_x_R_x_ET[1] += E[16] * R_x_ET[6];
    E_x_R_x_ET[2] += E[20] * R_x_ET[0];
    E_x_R_x_ET[2] += E[23] * R_x_ET[3];
    E_x_R_x_ET[2] += E[22] * R_x_ET[2];
    E_x_R_x_ET[2] += E[28] * R_x_ET[8];
    E_x_R_x_ET[3] += E[32] * R_x_ET[2];
    E_x_R_x_ET[3] += E[37] * R_x_ET[7];
    E_x_R_x_ET[3] += E[30] * R_x_ET[0];
    E_x_R_x_ET[3] += E[31] * R_x_ET[1];
    E_x_R_x_ET[3] += E[34] * R_x_ET[4];
    E_x_R_x_ET[4] += E[41] * R_x_ET[1];
    E_x_R_x_ET[4] += E[40] * R_x_ET[0];
    E_x_R_x_ET[4] += E[45] * R_x_ET[5];
    E_x_R_x_ET[4] += E[42] * R_x_ET[2];
}

