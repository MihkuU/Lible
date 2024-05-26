#pragma once

#include <lible/types.hpp>
#include <lible/structure.hpp>

#include <cassert>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         * Module for calculating molecular integrals.
         */
        struct test
        {
            /**
             * Tseburek
             */

            void testFun();
        };

        /**
         *
         */
        vec2d overlap(const Structure &structure);

        /**
         *
         */
        vec2d kineticEnergy(const Structure &structure);

        /**
         *
         */
        vec2d nuclearAttraction(const Structure &structure);

        /**
         *
         */
        vec2d dipoleMoment(const Structure &structure);

        /**
         *
         */
        vec2d eri2();

        /**
         *
         */
        vec3d eri3();

        /**
         *
         */
        vec4d eri4(const Structure &structure);

        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        /**
         *
         */
        std::set<std::string> returnAvailableBasisSets();

        /**
         *
         */
        std::map<int, std::vector<shell_exps_coeffs_t>>
        returnBasisForAtom(const int atomic_nr, const std::string &basis_set);

        /**
         *
         */
        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        returnBasisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);

        namespace gpu
        {
            /**
             *
             */
            vec2d overlap();
        }
    }
}
