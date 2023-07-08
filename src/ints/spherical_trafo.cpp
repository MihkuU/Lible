#include <cassert>

#include "defs_ints.h"
#include "spherical_trafo.h"

using namespace Lible;

template <>
void SphericalTrafo::transformAlongIdx<SphericalTrafo::Idx::FIRST>(const Shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_in,
                                                                   const std::vector<trafo_coeff_tuple> &spherical_trafo, std::vector<double> &one_el_ints_out)
{
    std::size_t dim_cart_b = shell_pair.second.dim_cartesian;
    for (auto &item : spherical_trafo)
        for (std::size_t j = 0; j < dim_cart_b; j++)
            one_el_ints_out[std::get<0>(item) * dim_cart_b + j] += std::get<2>(item) * one_el_ints_in[std::get<1>(item) * dim_cart_b + j];
}

template <>
void SphericalTrafo::transformAlongIdx<SphericalTrafo::Idx::SECOND>(const Shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_in,
                                                                    const std::vector<trafo_coeff_tuple> &spherical_trafo, std::vector<double> &one_el_ints_out)
{
    std::size_t dim_cart_b = shell_pair.second.dim_cartesian;
    std::size_t dim_sph_a = shell_pair.first.dim_spherical;
    std::size_t dim_sph_b = shell_pair.second.dim_cartesian;
    for (auto &item : spherical_trafo)
        for (std::size_t i = 0; i < dim_sph_a; i++)
            one_el_ints_out[i * dim_sph_b + std::get<0>(item)] += std::get<2>(item) * one_el_ints_in[i * dim_cart_b + std::get<1>(item)];
}

void SphericalTrafo::transformCartesianIntsToSpherical(const Shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_cart,
                                                       const std::vector<trafo_coeff_tuple> &spherical_trafo_first,
                                                       const std::vector<trafo_coeff_tuple> &spherical_trafo_second,
                                                       std::vector<double> &one_el_ints_sph_cart, std::vector<double> &one_el_ints_sph)
{
    assert((shell_pair.first.angular_momentum <= IntsDefs::max_angular_momentum));
    assert((shell_pair.second.angular_momentum <= IntsDefs::max_angular_momentum));

    transformAlongIdx<SphericalTrafo::Idx::FIRST>(shell_pair, one_el_ints_cart, spherical_trafo_first, one_el_ints_sph_cart);
    transformAlongIdx<SphericalTrafo::Idx::SECOND>(shell_pair, one_el_ints_sph_cart, spherical_trafo_second, one_el_ints_sph);
}