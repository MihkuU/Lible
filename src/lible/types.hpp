#pragma once

#include <lible/vectormd.hpp>

#include <array>
#include <vector>

namespace lible
{
    // STL type aliases.

    /// Type alias for a 2D static array of type `T`.
    template <typename T, size_t d1, size_t d2>
    using arr2d = std::array<std::array<T, d2>, d1>;
    /// Type alias for a 3D static array of type `T`.
    template <typename T, size_t d1, size_t d2, size_t d3>
    using arr3d = std::array<std::array<std::array<T, d3>, d2>, d1>;
    /// Type alias for 4D static array of type `T`.
    template <typename T, size_t d1, size_t d2, size_t d3, size_t d4>
    using arr4d = std::array<std::array<std::array<std::array<T, d4>, d3>, d2>, d1>;

    /// Type alias for a nested (2-dimensional) vector.
    template <typename T>
    using vecvec = std::vector<std::vector<T>>;

    // Internal type aliases.

    /// Type alias for a 2D array of doubles.
    using vec2d = VectorMD<double, 2>;
    /// Type alias for a 3D array of doubles.
    using vec3d = VectorMD<double, 3>;
    /// Type alias for a 4D array of doubles.
    using vec4d = VectorMD<double, 4>;
    /// Type alias for a 3D array of integers.
    using vec3i = VectorMD<int, 3>;
}
