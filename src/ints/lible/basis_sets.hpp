#pragma once

#include <string>
#include <unordered_map>

namespace lible
{
    namespace ints
    {
        enum class BasisFamily
        {
            ahlrichs
        };

        std::string returnBasisFamily(const std::string &basis_set);

        std::string returnBasisPath(const std::string &basis_set);
    }
}