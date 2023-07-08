#pragma once

#include <armadillo>
#include <cassert>
#include <map>
#include <vector>

#include "geomopt.h"

namespace Lible
{

    class GeomOpt::Geometry
    {
    public:
        Geometry(const std::vector<double> &coords, const std::vector<std::string> &atoms);        

    private:
        size_t n_atoms;
        std::vector<std::string> atoms;
        std::vector<arma::dvec> cartesian_coords;

        /*
         * Implementation of redundant internal coordinates as defined in:
         * https://doi.org/10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0.
         */
        struct RedundantInternalCoordinates
        {
            std::vector<double> bonds;
            std::vector<double> angles;
            std::vector<double> dihedrals;
        };
        RedundantInternalCoordinates red_int_coords;

        void constructRedIntCoords();
    };

}