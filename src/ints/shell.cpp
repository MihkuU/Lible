#include <lible/shell.h>
#include <lible/spherical_trafo.h>

using namespace lible;
using std::array;
using std::size_t;
using std::vector;

vector<double> shells::calcShellNormalization(const int &angmom,
                                              const vector<double> &contraction_coeffs,
                                              const vector<double> &contraction_exps,
                                              const vector<array<int, 3>> &cartesian_exps)
{
    int dim_cart_ints = (angmom + 1) * (angmom + 1);
    vector<double> overlap_x(dim_cart_ints, 0);
    vector<double> overlap_y(dim_cart_ints, 0);
    vector<double> overlap_z(dim_cart_ints, 0);

    size_t dim_cart = calcShellDimCartesian(angmom);
    size_t dim_sph = calcShellDimSpherical(angmom);

    vector<double> one_el_ints_cart(dim_cart * dim_cart, 0);
    size_t contraction_dim = contraction_exps.size();
    for (size_t a = 0; a < contraction_dim; a++)
    {
        double exp_a = contraction_exps[a];
        for (size_t b = 0; b < contraction_dim; b++)
        {
            double exp_b = contraction_exps[b];

            double p = exp_a + exp_b;
            double one_over_2p = 1.0 / (2 * p);

            double overlap_00 = std::sqrt(M_PI / p);
            // Kernels1El::ObaraSaika::overlapX(angmom, angmom, one_over_2p, 0, 0, overlap_00, overlap_x);
            overlap_y = overlap_x;
            overlap_z = overlap_x;

            double dadb = contraction_coeffs[a] * contraction_coeffs[b];
            for (size_t ikm = 0, ikmjln = 0; ikm < dim_cart; ikm++)
            {
                auto cart_exps_a = cartesian_exps[ikm];
                int i = cart_exps_a[0];
                int k = cart_exps_a[1];
                int m = cart_exps_a[2];
                for (size_t jln = 0; jln < dim_cart; jln++, ikmjln++)
                {
                    auto cart_exps_b = cartesian_exps[jln];
                    int j = cart_exps_b[0];
                    int l = cart_exps_b[1];
                    int n = cart_exps_b[2];

                    int ij = i * (angmom + 1) + j;
                    int kl = k * (angmom + 1) + l;
                    int mn = m * (angmom + 1) + n;
                    one_el_ints_cart[ikmjln] += dadb * overlap_x[ij] * overlap_y[kl] * overlap_z[mn];
                }
            }
        }
    }

    vector<trafo_coeff_tuple> spherical_trafo_first = SphericalTrafo::returnSphericalTrafo(angmom);
    vector<trafo_coeff_tuple> spherical_trafo_second = SphericalTrafo::returnSphericalTrafo(angmom);
    Shell shell(angmom, dim_cart, dim_sph);
    ShellPair shell_pair(shell, shell);
    vector<double> one_el_ints_sph_cart(dim_sph * dim_cart, 0);
    vector<double> one_el_ints_sph(dim_sph * dim_sph, 0);
    
    SphericalTrafo::transformCartesianIntsToSpherical(shell_pair,
                                                      one_el_ints_cart,
                                                      spherical_trafo_first,
                                                      spherical_trafo_second,
                                                      one_el_ints_sph_cart,
                                                      one_el_ints_sph);

    vector<double> normalization(dim_sph);
    for (size_t i = 0; i < dim_sph; i++)
        normalization[i] = 1.0 / std::sqrt(one_el_ints_sph[i * dim_sph + i]);

    return normalization;
}