#pragma once

#include <tuple>

#include "shell.h"

namespace Lible
{
    typedef std::tuple<std::size_t, std::size_t, double> trafo_coeff_tuple; // TODO is this place appropriate??

    namespace SphericalTrafo
    {
        /*
         * TODO: explain here the conventions regarding ordering of spherical gaussian functions and
         * cartesian gaussian functions.
         */
        std::vector<trafo_coeff_tuple> returnSphericalTrafo(const int &angmom);

        void transformCartesianIntsToSpherical(const Shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_cart,
                                               const std::vector<trafo_coeff_tuple> &spherical_trafo_first, 
                                               const std::vector<trafo_coeff_tuple> &spherical_trafo_second,
                                               std::vector<double> &one_el_ints_sph_cart, std::vector<double> &one_el_ints_sph);

        enum class Idx
        {
            FIRST,
            SECOND,
            THIRD,
            FOURTH
        };

        template <Idx idx>
        void transformAlongIdx(const Shells::ShellPair &shell_pair, const std::vector<double> &one_el_ints_in,
                               const std::vector<trafo_coeff_tuple> &spherical_trafo, std::vector<double> &one_el_ints_out);

    }
}