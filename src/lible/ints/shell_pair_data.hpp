#pragma once

#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        struct ShellPairData
        {
            ShellPairData();

            ShellPairData(const int la, const int lb, const Structure &structure);

            const Structure *structure;

            int la;
            int lb;
            size_t n_pairs;

            std::vector<std::pair<size_t, size_t>> offsets;
            std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> coords;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> ccoeffs;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> exps;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> norms;
        };
    }
}