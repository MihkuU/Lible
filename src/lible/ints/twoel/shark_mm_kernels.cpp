#include <lible/ints/twoel/shark_mm_kernels.hpp>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LI = lible::ints;

void LI::shark_mm_bra(int m, int n, int k, const double *E, const double *R_x_ET,
                      double *E_x_R_x_ET)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, E, k, R_x_ET, n, 1.0,
                E_x_R_x_ET, n);
}

void LI::shark_mm_ket(int m, int n, int k, const double *R, const double *E,
                      double *R_x_ET)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, R, k, E, n, 1.0,
                R_x_ET, n);
}