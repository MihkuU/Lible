#include <lible/basis_sets.hpp>
#include <lible/ints.hpp>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include <fmt/core.h>

namespace LI = lible::ints;

using std::set, std::string, std::unordered_map;

static unordered_map<string, LI::BasisFamily>
    basis_families{{"def2-svp", LI::BasisFamily::ahlrichs},
                   {"def2-tzvp", LI::BasisFamily::ahlrichs},
                   {"ahlrichs-pvdz", LI::BasisFamily::ahlrichs},
                   {"ahlrichs-tzv", LI::BasisFamily::ahlrichs},
                   {"ahlrichs-vdz", LI::BasisFamily::ahlrichs},
                   {"def2-qzvpd", LI::BasisFamily::ahlrichs},
                   {"def2-qzvppd", LI::BasisFamily::ahlrichs},
                   {"def2-qzvpp", LI::BasisFamily::ahlrichs},
                   {"def2-qzvp", LI::BasisFamily::ahlrichs},
                   {"def2-svpd", LI::BasisFamily::ahlrichs},
                   {"def2-sv(p)", LI::BasisFamily::ahlrichs},
                   {"def2-svp", LI::BasisFamily::ahlrichs},
                   {"def2-tzvpd", LI::BasisFamily::ahlrichs},
                   {"def2-tzvppd", LI::BasisFamily::ahlrichs},
                   {"def2-tzvpp", LI::BasisFamily::ahlrichs},
                   {"def2-tzvp", LI::BasisFamily::ahlrichs}};

set<string> LI::returnAvailableBasisSets()
{
    set<string> basis_sets;
    for (auto &item : basis_families)
        basis_sets.insert(item.first);

    return basis_sets;
}


string LI::returnBasisFamily(const string &basis_set)
{
    string basis_set_lc = basis_set;
    std::transform(basis_set_lc.begin(), basis_set_lc.end(), basis_set_lc.begin(),
                   [](auto c)
                   { return std::tolower(c); });

    if (basis_families.find(basis_set_lc) == basis_families.end())
        throw std::runtime_error(fmt::format("The requested basis set {} could not be found!",
                                             basis_set));
    else
    {
        switch (basis_families.at(basis_set_lc))
        {
        case BasisFamily::ahlrichs:
            return "ahlrichs";
        default:
            throw std::runtime_error("");
        }
    }
}