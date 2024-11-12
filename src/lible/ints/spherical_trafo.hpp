#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <armadillo>
#include <tuple>

namespace lible
{
    // typedef std::tuple<std::size_t, std::size_t, double> trafo_coeff_tuple; // TODO is this place appropriate??

    namespace ints
    {
        /*
         * TODO: explain here the conventions regarding ordering of spherical gaussian functions and
         * cartesian gaussian functions.
         */

        /**
         *
         */
        std::vector<std::tuple<int, int, double>> sphericalTrafo(const int l);

        /**
         *
         */
        arma::dmat returnSphericalTrafo(const int l);

        // TODO: rename the eri2_shells to eri2_shellbatch or something alike?

        /**
         *
         */
        void transferInts1El(const int ipair, const ShellPairData &sp_data,
                             const arma::dmat &ints_sph, vec2d &ints);

        /**
         *
         */
        void transferIntsERI2Diag(const int ishell, const ShellData &sh_data,
                                  const std::vector<double> &eri2_shells_sph,
                                  std::vector<double> &eri2_diagonal);

        /**
         *
         */
        void transferIntsERI4Diag(const int ipair_ab, const ShellPairData &sp_data_ab,
                                  const std::vector<double> &eri4_shells_sph,
                                  vec2d &eri4_diagonal);

        /**
         *
         */
        void transferIntsERI2(const int ishell_a, const int ishell_b,
                              const ShellData &sh_data_a, const ShellData &sh_data_b,
                              const std::vector<double> &eri2_batch, vec2d &eri2);

        /**
         *
         */
        void transferIntsERI3(const int ipair_ab, const int ishell_c,
                              const ShellData &sh_data_c, const ShellPairData &sp_data_ab,
                              const std::vector<double> &eri4_batch, vec3d &eri3);

        /**
         *
         */
        void transferIntsERI4(const int ipair_ab, const int ipair_cd,
                              const ShellPairData &sp_data_ab,
                              const ShellPairData &sp_data_cd,
                              const std::vector<double> &eri4_batch, vec4d &eri4);
    }
}