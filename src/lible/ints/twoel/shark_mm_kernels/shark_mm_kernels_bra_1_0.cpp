#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra1<1, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[5] * R_x_ET[1];
    E_x_R_x_ET[1] += E[4] * R_x_ET[0];
    E_x_R_x_ET[2] += E[8] * R_x_ET[0];
    E_x_R_x_ET[2] += E[10] * R_x_ET[2];
}

template<> void lible::ints::shark_mm_bra2<0, 1, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[5] * R_x_ET[1];
    E_x_R_x_ET[1] += E[4] * R_x_ET[0];
    E_x_R_x_ET[2] += E[8] * R_x_ET[0];
    E_x_R_x_ET[2] += E[10] * R_x_ET[2];
}

template<> void lible::ints::shark_mm_bra2<1, 0, 0>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[0] += E[3] * R_x_ET[3];
    E_x_R_x_ET[1] += E[5] * R_x_ET[1];
    E_x_R_x_ET[1] += E[4] * R_x_ET[0];
    E_x_R_x_ET[2] += E[8] * R_x_ET[0];
    E_x_R_x_ET[2] += E[10] * R_x_ET[2];
}

