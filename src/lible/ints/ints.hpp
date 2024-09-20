#pragma once

#include <lible/types.hpp>
#include <lible/ints/structure.hpp>

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
        std::vector<double> eri2Diagonal(const Structure &structure);
        
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
        vec2d eri4Diagonal(const Structure &structure);

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

        /**
         *
         */
        void eri4Benchmark(const Structure &structure);

        /**
         *
         */
        std::set<std::string> availableBasisSets();

        /**
         *
         */
        std::set<std::string> availableBasisSetsAux();

        /** */
        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;
        
        /**
         *
         */
        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtom(const int atomic_nr, const std::string &basis_set);

        /**
         *
         */
        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtomAux(const int atomic_nr, const std::string &aux_basis_set);        

        /**
         *
         */
        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);

        /**
         *
         */
        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtomsAux(const std::set<int> &atomic_nrs, const std::string &aux_basis_set);

        /**
         *
         */
        std::vector<std::array<int, 3>> cartExps(const int l);

        /**
         *
         */
        std::vector<std::tuple<int, int, double>> sphericalTrafo(const int l);

#ifdef _LIBLE_USE_HIP_
        namespace gpu
        {
            /**
             *
             */
            vec2d overlap0(const Structure &structure);

            /**
             *
             */
            vec2d overlap(const Structure &structure);
        }
#endif
    }
}
