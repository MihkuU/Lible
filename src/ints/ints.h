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
        vec2d overlap(const Structure &structure);
        
        vec2d kineticEnergy(const Structure &structure);

        vec2d coulombAttraction(const Structure &structure);

        vec2d dipoleMoment(const Structure &structure);        

        namespace Kernels
        {
            vec2d overlap();
            
            vec2d kineticEnergy();

            vec2d coulombAttraction();

            vec2d dipoleMoment();
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