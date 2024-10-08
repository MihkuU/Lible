#pragma once

#include <functional>
#include <utility>
#include <vector>

namespace lible
{
    namespace davidson
    {        
        std::pair<std::vector<double>, std::vector<std::vector<double>>>
        diagonalize(const size_t &n_roots,
                    const std::function<std::vector<double>()> &calcDiag,
                    const std::function<std::vector<std::vector<double>>(const std::vector<double> &diag)> &calcGuess,
                    const std::function<std::vector<double>(const std::vector<double> &trial)> &calcSigma);

        std::pair<std::vector<double>, std::vector<std::vector<double>>>
        diagonalize(const size_t &n_roots,
                    const std::function<std::vector<double>()> &calcDiag,
                    const std::function<std::vector<std::vector<double>>(const std::vector<double> &diag)> &calcGuess,
                    const std::function<std::vector<double>(const std::vector<double> &trial)> &calcSigma,
                    const std::function<std::vector<double>()> &preConditioner);        

        struct Aux
        {
            static inline size_t n_iters = 0;
        };
    }
}