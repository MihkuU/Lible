// #include "kernels_1el_ints.h"
#include "ints.h"
#include "spherical_trafo.h"

using namespace lible;
using namespace lible::ints;

using std::array;
using std::size_t;
using std::vector;

Vec2D<double> calcOverlap(const Structure &Structure)
{

}

// template <>
// void Kernels1El::oneElIntKernel<ints::OVERLAP>(const int &la,
//                                                const int &lb,
//                                                const size_t &n_ao,
//                                                const vector<Shells::ShellPair> &shell_pairs,
//                                                vector<double> &one_el_ints)
// {
//     array<double, 3> overlap_00;
//     array<double, 3> xyz_p;
//     array<double, 3> xyz_pa;
//     array<double, 3> xyz_pb;

//     int dim_cart_ints = (la + 1) * (lb + 1);
//     vector<double> overlap_x(dim_cart_ints, 0);
//     vector<double> overlap_y(dim_cart_ints, 0);
//     vector<double> overlap_z(dim_cart_ints, 0);

//     size_t dim_cart_a = Shells::calcShellDimCartesian(la);
//     size_t dim_cart_b = Shells::calcShellDimCartesian(lb);
//     size_t dim_sph_a = Shells::calcShellDimSpherical(la);
//     size_t dim_sph_b = Shells::calcShellDimSpherical(lb);
//     vector<double> one_el_ints_cart(dim_cart_a * dim_cart_b, 0);
//     vector<double> one_el_ints_sph_cart(dim_sph_a * dim_cart_b, 0);
//     vector<double> one_el_ints_sph(dim_sph_a * dim_sph_b, 0);

//     vector<trafo_coeff_tuple> spherical_trafo_first = SphericalTrafo::returnSphericalTrafo(la);
//     vector<trafo_coeff_tuple> spherical_trafo_second = SphericalTrafo::returnSphericalTrafo(lb);
//     for (size_t ishell_pair = 0; ishell_pair < shell_pairs.size(); ishell_pair++)
//     {
//         auto shell_pair = shell_pairs[ishell_pair];
//         array<double, 3> xyz_a = shell_pair.first.xyz_coordinates;
//         array<double, 3> xyz_b = shell_pair.second.xyz_coordinates;

//         std::fill(one_el_ints_cart.begin(), one_el_ints_cart.end(), 0);
//         std::fill(one_el_ints_sph_cart.begin(), one_el_ints_sph_cart.end(), 0);
//         std::fill(one_el_ints_sph.begin(), one_el_ints_sph.end(), 0);

//         for (size_t a = 0; a < shell_pair.first.contraction_coeffs.size(); a++)
//         {
//             double exp_a = shell_pair.first.contraction_exps[a];
//             for (size_t b = 0; b < shell_pair.second.contraction_coeffs.size(); b++)
//             {
//                 double exp_b = shell_pair.second.contraction_exps[b];

//                 double p = exp_a + exp_b;
//                 double one_over_2p = 1.0 / (2 * p);
//                 double mu = exp_a * exp_b / p;

//                 double prefac = sqrt(M_PI / p);
//                 overlap_00[0] = prefac * exp(-mu * pow(xyz_a[0] - xyz_b[0], 2));
//                 overlap_00[1] = prefac * exp(-mu * pow(xyz_a[1] - xyz_b[1], 2));
//                 overlap_00[2] = prefac * exp(-mu * pow(xyz_a[2] - xyz_b[2], 2));

//                 xyz_p[0] = (exp_a * xyz_a[0] + exp_b * xyz_b[0]) / p;
//                 xyz_p[1] = (exp_a * xyz_a[1] + exp_b * xyz_b[1]) / p;
//                 xyz_p[2] = (exp_a * xyz_a[2] + exp_b * xyz_b[2]) / p;
//                 xyz_pa[0] = xyz_p[0] - xyz_a[0];
//                 xyz_pa[1] = xyz_p[1] - xyz_a[1];
//                 xyz_pa[2] = xyz_p[2] - xyz_a[2];
//                 xyz_pb[0] = xyz_p[0] - xyz_b[0];
//                 xyz_pb[1] = xyz_p[1] - xyz_b[1];
//                 xyz_pb[2] = xyz_p[2] - xyz_b[2];

//                 ObaraSaika::overlapX(la, lb, one_over_2p, xyz_pa[0], xyz_pb[0], overlap_00[0], overlap_x);
//                 ObaraSaika::overlapX(la, lb, one_over_2p, xyz_pa[1], xyz_pb[1], overlap_00[1], overlap_y);
//                 ObaraSaika::overlapX(la, lb, one_over_2p, xyz_pa[2], xyz_pb[2], overlap_00[2], overlap_z);

//                 double dadb = shell_pair.first.contraction_coeffs[a] * shell_pair.second.contraction_coeffs[b];
//                 for (size_t ikm = 0, ikmjln = 0; ikm < dim_cart_a; ikm++)
//                 {
//                     auto cart_exps_a = shell_pair.first.cartesian_exps[ikm];
//                     int i = cart_exps_a[0];
//                     int k = cart_exps_a[1];
//                     int m = cart_exps_a[2];
//                     for (size_t jln = 0; jln < dim_cart_b; jln++, ikmjln++)
//                     {
//                         auto cart_exps_b = shell_pair.second.cartesian_exps[jln];
//                         int j = cart_exps_b[0];
//                         int l = cart_exps_b[1];
//                         int n = cart_exps_b[2];

//                         int ij = i * (lb + 1) + j;
//                         int kl = k * (lb + 1) + l;
//                         int mn = m * (lb + 1) + n;
//                         one_el_ints_cart[ikmjln] += dadb * overlap_x[ij] * overlap_y[kl] * overlap_z[mn];
//                     }
//                 }
//             }
//         }

//         SphericalTrafo::transformCartesianIntsToSpherical(shell_pair,
//                                                           one_el_ints_cart,
//                                                           spherical_trafo_first,
//                                                           spherical_trafo_second,
//                                                           one_el_ints_sph_cart,
//                                                           one_el_ints_sph);

//         SphericalTrafo::transferSphericalInts(n_ao,
//                                               shell_pair,
//                                               one_el_ints_sph,
//                                               one_el_ints);
//     }
// }

// void Kernels1El::ObaraSaika::overlapX(const int &angmom_a,
//                                       const int &angmom_b,
//                                       const double &one_over_2p,
//                                       const double &x_pa,
//                                       const double &x_pb,
//                                       const double &overlap_00,
//                                       vector<double> &overlap_x)
// {
//     /*
//      * Implements eq. (9.3.9) from DOI:10.1002/9781119019572.
//      * Calculates the overlap integral along one chosen cartesian direction.
//      */

//     overlap_x[0] = overlap_00;

//     int dim_b = angmom_b + 1;
//     /* Si0 */
//     if (angmom_a > 0)
//     {
//         overlap_x[dim_b] = x_pa * overlap_00;

//         for (int i = 2; i <= angmom_a; i++)
//             overlap_x[i * dim_b] = x_pa * overlap_x[(i - 1) * dim_b] + one_over_2p * (i - 1) * overlap_x[(i - 2) * dim_b];
//     }

//     if (angmom_b > 0)
//     {
//         /* S0j */
//         overlap_x[1] = x_pb * overlap_00;

//         for (int j = 2; j <= angmom_b; j++)
//             overlap_x[j] = x_pb * overlap_x[j - 1] + one_over_2p * (j - 1) * overlap_x[j - 2];

//         /* Si1 */
//         for (int i = 1; i <= angmom_a; i++)
//             overlap_x[i * dim_b + 1] = x_pb * overlap_x[i * dim_b] + one_over_2p * i * overlap_x[(i - 1) * dim_b];

//         /* Sij */
//         if (angmom_b > 1)
//         {
//             for (int i = 1; i <= angmom_a; i++)
//                 for (int j = 2; j <= angmom_b; j++)
//                     overlap_x[i * dim_b + j] = x_pb * overlap_x[i * dim_b + j - 1] +
//                                                one_over_2p * (i * overlap_x[(i - 1) * dim_b + j - 1] + (j - 1) * overlap_x[i * dim_b + j - 2]);
//         }
//     }
// }