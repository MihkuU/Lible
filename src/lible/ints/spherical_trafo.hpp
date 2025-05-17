#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <tuple>

namespace lible
{
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

        /** */
        vec2d trafo2Spherical(const int la, const int lb, const vec2d &ints_cart);
    }
}