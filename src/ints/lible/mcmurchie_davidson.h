#pragma once

#include <lible/types.h>

namespace lible
{
    namespace ints
    {
        namespace MD
        {
            void coeffs(const double &a, const double &b, const double &PA, const double &PB,
                        const double one_o_2p, const int &la, const int &lb, vec3d &E);

            void coeffs(const double &a, const double &b, const int &la, const int &lb,
                        const std::array<double, 3> &A, const std::array<double, 3> &B,
                        vec3d &E_x, vec3d &E_y, vec3d &E_z);

            void coeffsSpherical();

            void test();
        }
    }
}