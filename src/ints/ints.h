#pragma once

#include <cassert>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <lible/types.h>
#include "structure.h"

namespace lible
{
    namespace ints
    {
        Vec2D<double> overlap(const Structure &structure);

        Vec2D<double> kineticEnergy(const Structure &structure);

        Vec2D<double> coulombAttraction(const Structure &structure);

        Vec2D<double> dipoleMoment(const Structure &structure);

        namespace Kernels
        {
            Vec2D<double> overlap();
            
            Vec2D<double> kineticEnergy();

            Vec2D<double> coulombAttraction();

            Vec2D<double> dipoleMoment();
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