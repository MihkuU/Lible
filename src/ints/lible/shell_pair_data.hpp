#pragma once

#include <lible/structure.hpp>

namespace lible
{
    namespace ints
    {
        struct ShellPairData
        {
            ShellPairData(const int &la, const int &lb, const Structure &structure);
            
            int la;
            int lb;
            size_t n_pairs;

            // std::vector<std::pair<size_t, size_t>> dims_cart;
            // std::vector<std::pair<size_t, size_t>> dims_sph;
            std::vector<std::pair<size_t, size_t>> offsets;
            std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> coords;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> coeffs;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> exps;            
            std::vector<std::pair<std::vector<double>, std::vector<double>>> norms;
        };
    }
}