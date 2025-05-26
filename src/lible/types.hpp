#pragma once

#include <lible/vec2d.hpp>
#include <lible/vec3d.hpp>
#include <lible/vec4d.hpp>
#include <lible/vectormd.hpp>

#include <array>
#include <vector>

namespace lible
{
    using vec2d = Vec2D<double>;
    using vec3d = Vec3D<double>;
    using vec4d = Vec4D<double>;

    using vec3i = Vec3D<int>;

    // using vec2d = VectorMD<double, 2>;
    // using vec3d = VectorMD<double, 3>;
    // using vec4d = VectorMD<double, 4>;
}