#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_bra1<0, 2>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
}

template<> void lible::ints::shark_mm_bra2<0, 0, 2>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)
{
    E_x_R_x_ET[0] += E[0] * R_x_ET[0];
    E_x_R_x_ET[1] += E[0] * R_x_ET[1];
    E_x_R_x_ET[2] += E[0] * R_x_ET[2];
    E_x_R_x_ET[3] += E[0] * R_x_ET[3];
    E_x_R_x_ET[4] += E[0] * R_x_ET[4];
}

