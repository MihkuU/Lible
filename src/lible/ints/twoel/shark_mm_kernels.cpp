#include <lible/ints/twoel/shark_mm_kernels.hpp>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace lints = lible::ints;

void lints::shark_mm_bra(const int m, const int n, const int k, const double *E,
                         const double *R_x_ET, double *E_x_R_x_ET)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, E, k, R_x_ET, n, 1.0,
                E_x_R_x_ET, n);
}

void lints::shark_mm_ket(const int m, const int n, const int k, const double *R, const double *ET,
                         double *R_x_ET)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, R, k, ET, n, 1.0,
                R_x_ET, n);
}
