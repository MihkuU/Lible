#pragma once

#include <lible/types.h>
#include <lible/structure.h>

#include <cassert>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace lible
{
    namespace ints
    {
        vec2d overlap(const Structure &structure);

        vec2d kineticEnergy(const Structure &structure);

        vec2d coulombAttraction(const Structure &structure);

        vec2d dipoleMoment(const Structure &structure);

        namespace kernels
        {
            void overlap();

            void kineticEnergy();

            void coulombAttraction();

            void dipoleMoment();
        }
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