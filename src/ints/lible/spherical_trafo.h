#pragma once

#include <tuple>

#include <lible/shell.h>

namespace lible
{
    typedef std::tuple<std::size_t, std::size_t, double> trafo_coeff_tuple; // TODO is this place appropriate??

    namespace SphericalTrafo
    {
        /*
         * TODO: explain here the conventions regarding ordering of spherical gaussian functions and
         * cartesian gaussian functions.
         */
        std::vector<trafo_coeff_tuple> returnSphericalTrafo(const int &angmom);

        void transformCartesianIntsToSpherical(const shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_cart,
                                               const std::vector<trafo_coeff_tuple> &spherical_trafo_first, 
                                               const std::vector<trafo_coeff_tuple> &spherical_trafo_second,
                                               std::vector<double> &one_el_ints_sph_cart, std::vector<double> &one_el_ints_sph);

        void transferSphericalInts(const std::size_t &n_ao, const shells::ShellPair &shell_pair, const std::vector<double> &ints_in,
                                   std::vector<double> &ints_out);

        enum class Idx
        {
            FIRST,
            SECOND,
            THIRD,
            FOURTH
        };

        template <Idx idx>
        void transformAlongIdx(const shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_in,
                               const std::vector<trafo_coeff_tuple> &spherical_trafo, std::vector<double> &one_el_ints_out);

    }
}