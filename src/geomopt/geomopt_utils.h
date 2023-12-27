#pragma once

#include "geomopt.h"
#include "geometry.h"
#include <vector>

namespace Lible
{
    namespace GeomOpt
    {
        double calcGradNorm(const std::vector<double> &grad);

        void optimizePrintPreamble(const Geometry &geometry);
        void optimizePrintEpilogue(const bool &converged);

        template <Option option>
        void optimizePrintIter(const std::size_t &iter);

    }
}