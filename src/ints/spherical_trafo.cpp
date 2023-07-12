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
    std::size_t dim_sph_b = shell_pair.second.dim_spherical;
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

void SphericalTrafo::transferSphericalInts(const std::size_t &n_ao, const Shells::ShellPair &shell_pair, const std::vector<double> &ints_in,
                                           std::vector<double> &ints_out)
{
    assert((shell_pair.first.angular_momentum >= shell_pair.second.angular_momentum));

    std::size_t pos_a = shell_pair.first.pos;
    std::size_t pos_b = shell_pair.second.pos;

    for (std::size_t i = 0; i < shell_pair.first.dim_spherical; i++)
    {
        for (std::size_t j = 0; j < shell_pair.second.dim_spherical; j++)
        {
            double normalized_int = shell_pair.first.normalization[i] * shell_pair.second.normalization[j] *
                                    ints_in[i * shell_pair.second.dim_spherical + j];
            ints_out[(pos_a + i) * n_ao + (pos_b + j)] = normalized_int;
            ints_out[(pos_b + j) * n_ao + (pos_a + i)] = normalized_int;
        }
    }
}
