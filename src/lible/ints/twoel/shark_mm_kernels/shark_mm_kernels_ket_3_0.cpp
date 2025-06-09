#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_ket1<3, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[1] += R[1] * ET[0];
    R_x_ET[2] += R[2] * ET[0];
    R_x_ET[3] += R[3] * ET[0];
    R_x_ET[4] += R[4] * ET[0];
    R_x_ET[5] += R[5] * ET[0];
    R_x_ET[6] += R[6] * ET[0];
    R_x_ET[7] += R[7] * ET[0];
    R_x_ET[8] += R[8] * ET[0];
    R_x_ET[9] += R[9] * ET[0];
    R_x_ET[10] += R[10] * ET[0];
    R_x_ET[11] += R[11] * ET[0];
    R_x_ET[12] += R[12] * ET[0];
    R_x_ET[13] += R[13] * ET[0];
    R_x_ET[14] += R[14] * ET[0];
    R_x_ET[15] += R[15] * ET[0];
    R_x_ET[16] += R[16] * ET[0];
    R_x_ET[17] += R[17] * ET[0];
    R_x_ET[18] += R[18] * ET[0];
    R_x_ET[19] += R[19] * ET[0];
}

template<> void lible::ints::shark_mm_ket2<3, 0, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[1] += R[1] * ET[0];
    R_x_ET[2] += R[2] * ET[0];
    R_x_ET[3] += R[3] * ET[0];
    R_x_ET[4] += R[4] * ET[0];
    R_x_ET[5] += R[5] * ET[0];
    R_x_ET[6] += R[6] * ET[0];
    R_x_ET[7] += R[7] * ET[0];
    R_x_ET[8] += R[8] * ET[0];
    R_x_ET[9] += R[9] * ET[0];
    R_x_ET[10] += R[10] * ET[0];
    R_x_ET[11] += R[11] * ET[0];
    R_x_ET[12] += R[12] * ET[0];
    R_x_ET[13] += R[13] * ET[0];
    R_x_ET[14] += R[14] * ET[0];
    R_x_ET[15] += R[15] * ET[0];
    R_x_ET[16] += R[16] * ET[0];
    R_x_ET[17] += R[17] * ET[0];
    R_x_ET[18] += R[18] * ET[0];
    R_x_ET[19] += R[19] * ET[0];
}

