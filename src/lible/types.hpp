#pragma once

#include <lible/vectormd.hpp>

#include <array>
#include <vector>

namespace lible
{
    using vec2d = VectorMD<double, 2>;
    using vec3d = VectorMD<double, 3>;
    using vec4d = VectorMD<double, 4>;

    using vec3i = VectorMD<int, 3>;

    template <typename T, size_t d1, size_t d2>
    using arr2d = std::array<std::array<T, d2>, d1>;

    template <typename T, size_t d1, size_t d2, size_t d3>
    using arr3d = std::array<std::array<std::array<T, d3>, d2>, d1>;

    template <typename T, size_t d1, size_t d2, size_t d3, size_t d4>
    using arr4d = std::array<std::array<std::array<std::array<T, d4>, d3>, d2>, d1>;

    /// Type alias for a nested (2-dimensional) vector.
    template <typename T>
    using vecvec = std::vector<std::vector<T>>;
}
