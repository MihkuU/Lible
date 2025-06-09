#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra1<1, 1>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[9];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[10];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[11];
    E_x_R_x_ET[3] += E[5] * R_x_ET[3];
    E_x_R_x_ET[3] += E[4] * R_x_ET[0];
    E_x_R_x_ET[4] += E[5] * R_x_ET[4];
    E_x_R_x_ET[4] += E[4] * R_x_ET[1];
    E_x_R_x_ET[5] += E[5] * R_x_ET[5];
    E_x_R_x_ET[5] += E[4] * R_x_ET[2];
    E_x_R_x_ET[6] += E[8] * R_x_ET[0];
    E_x_R_x_ET[6] += E[10] * R_x_ET[6];
    E_x_R_x_ET[7] += E[8] * R_x_ET[1];
    E_x_R_x_ET[7] += E[10] * R_x_ET[7];
    E_x_R_x_ET[8] += E[8] * R_x_ET[2];
    E_x_R_x_ET[8] += E[10] * R_x_ET[8];
}

template<> void lible::ints::shark_mm_bra2<0, 1, 1>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[9];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[10];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[11];
    E_x_R_x_ET[3] += E[5] * R_x_ET[3];
    E_x_R_x_ET[3] += E[4] * R_x_ET[0];
    E_x_R_x_ET[4] += E[5] * R_x_ET[4];
    E_x_R_x_ET[4] += E[4] * R_x_ET[1];
    E_x_R_x_ET[5] += E[5] * R_x_ET[5];
    E_x_R_x_ET[5] += E[4] * R_x_ET[2];
    E_x_R_x_ET[6] += E[8] * R_x_ET[0];
    E_x_R_x_ET[6] += E[10] * R_x_ET[6];
    E_x_R_x_ET[7] += E[8] * R_x_ET[1];
    E_x_R_x_ET[7] += E[10] * R_x_ET[7];
    E_x_R_x_ET[8] += E[8] * R_x_ET[2];
    E_x_R_x_ET[8] += E[10] * R_x_ET[8];
}

template<> void lible::ints::shark_mm_bra2<1, 0, 1>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[9];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[1] += E[3] * R_x_ET[10];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[2] += E[3] * R_x_ET[11];
    E_x_R_x_ET[3] += E[5] * R_x_ET[3];
    E_x_R_x_ET[3] += E[4] * R_x_ET[0];
    E_x_R_x_ET[4] += E[5] * R_x_ET[4];
    E_x_R_x_ET[4] += E[4] * R_x_ET[1];
    E_x_R_x_ET[5] += E[5] * R_x_ET[5];
    E_x_R_x_ET[5] += E[4] * R_x_ET[2];
    E_x_R_x_ET[6] += E[8] * R_x_ET[0];
    E_x_R_x_ET[6] += E[10] * R_x_ET[6];
    E_x_R_x_ET[7] += E[8] * R_x_ET[1];
    E_x_R_x_ET[7] += E[10] * R_x_ET[7];
    E_x_R_x_ET[8] += E[8] * R_x_ET[2];
    E_x_R_x_ET[8] += E[10] * R_x_ET[8];
}

