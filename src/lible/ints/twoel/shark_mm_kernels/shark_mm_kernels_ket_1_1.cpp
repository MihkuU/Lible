#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_ket1<1, 1>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[3] * ET[9];
    R_x_ET[1] += R[1] * ET[4];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[2] * ET[8];
    R_x_ET[3] += R[4] * ET[0];
    R_x_ET[3] += R[7] * ET[9];
    R_x_ET[4] += R[5] * ET[4];
    R_x_ET[4] += R[4] * ET[1];
    R_x_ET[5] += R[4] * ET[2];
    R_x_ET[5] += R[6] * ET[8];
    R_x_ET[6] += R[8] * ET[0];
    R_x_ET[6] += R[11] * ET[9];
    R_x_ET[7] += R[9] * ET[4];
    R_x_ET[7] += R[8] * ET[1];
    R_x_ET[8] += R[8] * ET[2];
    R_x_ET[8] += R[10] * ET[8];
    R_x_ET[9] += R[12] * ET[0];
    R_x_ET[9] += R[15] * ET[9];
    R_x_ET[10] += R[13] * ET[4];
    R_x_ET[10] += R[12] * ET[1];
    R_x_ET[11] += R[12] * ET[2];
    R_x_ET[11] += R[14] * ET[8];
}

template<> void lible::ints::shark_mm_ket2<1, 0, 1>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[3] * ET[9];
    R_x_ET[1] += R[1] * ET[4];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[2] * ET[8];
    R_x_ET[3] += R[4] * ET[0];
    R_x_ET[3] += R[7] * ET[9];
    R_x_ET[4] += R[5] * ET[4];
    R_x_ET[4] += R[4] * ET[1];
    R_x_ET[5] += R[4] * ET[2];
    R_x_ET[5] += R[6] * ET[8];
    R_x_ET[6] += R[8] * ET[0];
    R_x_ET[6] += R[11] * ET[9];
    R_x_ET[7] += R[9] * ET[4];
    R_x_ET[7] += R[8] * ET[1];
    R_x_ET[8] += R[8] * ET[2];
    R_x_ET[8] += R[10] * ET[8];
    R_x_ET[9] += R[12] * ET[0];
    R_x_ET[9] += R[15] * ET[9];
    R_x_ET[10] += R[13] * ET[4];
    R_x_ET[10] += R[12] * ET[1];
    R_x_ET[11] += R[12] * ET[2];
    R_x_ET[11] += R[14] * ET[8];
}

template<> void lible::ints::shark_mm_ket2<1, 1, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[0] += R[3] * ET[9];
    R_x_ET[1] += R[1] * ET[4];
    R_x_ET[1] += R[0] * ET[1];
    R_x_ET[2] += R[0] * ET[2];
    R_x_ET[2] += R[2] * ET[8];
    R_x_ET[3] += R[4] * ET[0];
    R_x_ET[3] += R[7] * ET[9];
    R_x_ET[4] += R[5] * ET[4];
    R_x_ET[4] += R[4] * ET[1];
    R_x_ET[5] += R[4] * ET[2];
    R_x_ET[5] += R[6] * ET[8];
    R_x_ET[6] += R[8] * ET[0];
    R_x_ET[6] += R[11] * ET[9];
    R_x_ET[7] += R[9] * ET[4];
    R_x_ET[7] += R[8] * ET[1];
    R_x_ET[8] += R[8] * ET[2];
    R_x_ET[8] += R[10] * ET[8];
    R_x_ET[9] += R[12] * ET[0];
    R_x_ET[9] += R[15] * ET[9];
    R_x_ET[10] += R[13] * ET[4];
    R_x_ET[10] += R[12] * ET[1];
    R_x_ET[11] += R[12] * ET[2];
    R_x_ET[11] += R[14] * ET[8];
}

