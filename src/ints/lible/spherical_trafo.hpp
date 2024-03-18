#pragma once

#include <lible/types.hpp>
#include <lible/shell.hpp>
#include <lible/shell_pair_data.hpp>

#include <armadillo>
#include <tuple>

namespace lible
{

    typedef std::tuple<std::size_t, std::size_t, double> trafo_coeff_tuple; // TODO is this place appropriate??

    namespace ints
    {
        /*
         * TODO: explain here the conventions regarding ordering of spherical gaussian functions and
         * cartesian gaussian functions.
         */

        arma::dmat returnSphericalTrafo(const int &angmom);

        void transferIntegrals(const size_t &ipair,
                               const ShellPairData &shell_pair_data,
                               const arma::dmat &ints_sph,
                               vec2d &ints_out);
    }
}