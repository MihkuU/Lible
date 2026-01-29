#pragma once

#include <lible/types.hpp>

#include <tuple>

namespace lible::ints
{
    /// Returns the Cartesian to spherical transformation as {{mu, mu_, val}}, where mu and mu_
    /// refer to the spherical and Cartesian Gaussians, respectively.
    std::vector<std::tuple<int, int, double>> sphericalTrafo(int l);


}
