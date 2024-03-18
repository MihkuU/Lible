#pragma once

#include <lible/types.h>
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
        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        vec2d overlap(const Structure &structure);

        vec2d kineticEnergy(const Structure &structure);

        vec2d nuclearAttraction(const Structure &structure);

        vec2d dipoleMoment(const Structure &structure);
        
        std::map<int, std::vector<shell_exps_coeffs_t>>
        returnBasisForAtom(const int atomic_nr, const std::string &basis_set);

        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        returnBasisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);
    }

    // class ints
    // {
    // public:
    // ints(const std::string &basis_set, const std::vector<double> &coordinates_ang,
    //      const std::vector<std::string> &elements);

    // enum Option1El
    // {
    //     KINETIC,
    //     OVERLAP,
    //     POTENTIAL
    // };

    // enum Option2El
    // {
    //     COULOMB,
    //     EXCHANGE
    // };

    // template <Option1El option>
    // std::vector<double> calcOneElInts();

    // template <Option1El option>
    // std::vector<double> calcOneElIntsGPU();

    // template <Option2El option>
    // std::vector<double> calcTwoElInts(const std::vector<double> &density); // TODO: rethink the name of this function

    // struct Structure;
    // struct StructureGPU;
    // std::unique_ptr<Structure> structure;
    // std::unique_ptr<StructureGPU> structure_gpu;
    // };
}