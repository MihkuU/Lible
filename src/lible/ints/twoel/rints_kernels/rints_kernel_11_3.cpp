#include <lible/ints/rints_meta.hpp>

template void lible::ints::calcRInts_ERI_new<11, 3>(const double alpha, const double fac, const double *fnx,
                                                   const double *xyz_pq, const int n_cols, const int ofs_row, 
                                                   const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI2D1<11, 3>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int n_cols,
                                                  const int ofs_row, const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI3D1<11, 3>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int n_cols,
                                                  const int ofs_row, const int ofs_col, double *rints_out);

template void lible::ints::calcRInts_ERI4D1<11, 3>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pq, const int n_rints, const int ofs_row,
                                                  const int ofs_col, const int n_cols, const int n_rows,
                                                  double *rints_out);

template void lible::ints::calcRInts_ERI<11, 3>(const double, const double, const double*, const double*, double*);

template void lible::ints::calcRInts_ERI2_deriv1<11, 3>(const double, const double, const double*, const double*, double*);
