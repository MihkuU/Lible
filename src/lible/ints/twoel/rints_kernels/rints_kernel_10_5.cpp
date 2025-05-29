#include <lible/ints/rints_meta.hpp>

template void lible::ints::calcRInts_ERI_new<10, 5>(const double alpha, const double fac, const double *fnx,
                                                   const double *xyz_pq,  const int n_cols, const int ofs_row, 
                                                   const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI<10, 5>(const double, const double, const double*, const double*, double*);

template void lible::ints::calcRInts_ERI2_deriv1<10, 5>(const double, const double, const double*, const double*, double*);
