#pragma once

#include <cassert>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace Lible
{
    /*
     * A class that implements functions for calculating molecular integrals.
     *
     */
    class Ints
    {
    public:
        enum Option1El
        {
            KINETIC,
            OVERLAP,
            POTENTIAL
        };

        enum Option2El
        {
            COULOMB,
            EXCHANGE
        };

        Ints(const std::string &basis_set, const std::vector<double> &coordinates, const std::vector<std::string> &elements);

        template <Option1El option>
        std::vector<double> calcOneElInts();

        template <Option2El option>
        std::vector<double> calcTwoElInts(const std::vector<double> &density); // TODO: rethink the name of this function

    private:
        struct Structure;
        std::unique_ptr<Structure> structure;
    };
}