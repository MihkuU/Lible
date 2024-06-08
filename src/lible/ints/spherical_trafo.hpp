#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell.hpp>
#include <lible/ints/shell_pair_data.hpp>

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

        arma::dmat returnSphericalTrafo(const int angmom);

        void sphericalTrafo(const arma::dmat &trafo_a, const arma::dmat &trafo_b,
                            const arma::dmat &trafo_c, const arma::dmat &trafo_d,
                            const vec4d &eri4_shells_cart, vec4d &eri4_shells_sph);

        void transferIntegrals(const size_t ipair,
                               const ShellPairData &shell_pair_data,
                               const arma::dmat &ints_sph, vec2d &ints);

        void transferIntegrals(const size_t ipair_ab, const size_t ipair_cd,
                               const ShellPairData &shell_pair_data_ab,
                               const ShellPairData &shell_pair_data_cd,
                               const vec4d &eri4_shells_sph, vec4d &eri4);

        void transferIntegrals(const size_t ipair_ab, const size_t ipair_cd,
                               const ShellPairData &shell_pair_data_ab,
                               const ShellPairData &shell_pair_data_cd,
                               const arma::dmat &eri4_shells_sph, vec4d &eri4);
    }
}