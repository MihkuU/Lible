#include <lible/ints/basis_sets.hpp>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include <fmt/core.h>

namespace LI = lible::ints;

using std::string;

const static std::unordered_map<string, LI::AuxBasisFamily>
    aux_basis_families{{"def2-qzvp-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-qzvpp-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-qzvppd-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-sv(p)-jkfit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-sv(p)-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-svp-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-svpd-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-tzvp-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-tzvpd-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-tzvpp-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-tzvppd-rifit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-universal-jfit", LI::AuxBasisFamily::ahlrichs_fit},
                       {"def2-universal-jkfit", LI::AuxBasisFamily::ahlrichs_fit}};

const static std::unordered_map<string, LI::BasisFamily>
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
                   {"def2-tzvp", LI::BasisFamily::ahlrichs},
                   {"ano-dk3", LI::BasisFamily::ano},
                   {"ano-pv5z", LI::BasisFamily::ano},
                   {"ano-pvdz", LI::BasisFamily::ano},
                   {"ano-pvqz", LI::BasisFamily::ano},
                   {"ano-pvtz", LI::BasisFamily::ano},
                   {"ano-r", LI::BasisFamily::ano},
                   {"ano-r0", LI::BasisFamily::ano},
                   {"ano-r1", LI::BasisFamily::ano},
                   {"ano-r2", LI::BasisFamily::ano},
                   {"ano-r3", LI::BasisFamily::ano},
                   {"no-rcc-mb", LI::BasisFamily::ano},
                   {"ano-rcc-vdz", LI::BasisFamily::ano},
                   {"ano-rcc-vdzp", LI::BasisFamily::ano},
                   {"ano-rcc-vqzp", LI::BasisFamily::ano},
                   {"ano-rcc-vtz", LI::BasisFamily::ano},
                   {"ano-rcc-vtzp", LI::BasisFamily::ano},
                   {"ano-rcc", LI::BasisFamily::ano},
                   {"aug-ano-pv5z", LI::BasisFamily::ano},
                   {"aug-ano-pvdz", LI::BasisFamily::ano},
                   {"aug-ano-pvqz", LI::BasisFamily::ano},
                   {"aug-ano-pvtz", LI::BasisFamily::ano},
                   {"nasa-ames-ano", LI::BasisFamily::ano},
                   {"nasa-ames-ano2", LI::BasisFamily::ano},
                   {"roos-augmented-double-zeta-ano", LI::BasisFamily::ano},
                   {"roos-augmented-triple-zeta-ano", LI::BasisFamily::ano},
                   {"saug-ano-pv5z", LI::BasisFamily::ano},
                   {"saug-ano-pvdz", LI::BasisFamily::ano},
                   {"saug-ano-pvqz", LI::BasisFamily::ano},
                   {"saug-ano-pvtz", LI::BasisFamily::ano},
                   {"sto-2g", LI::BasisFamily::sto},
                   {"sto-3g", LI::BasisFamily::sto},
                   {"sto-3g*", LI::BasisFamily::sto},
                   {"sto-4g", LI::BasisFamily::sto},
                   {"sto-5g", LI::BasisFamily::sto},
                   {"sto-6g", LI::BasisFamily::sto}};

std::set<string> LI::availableBasisSets()
{
    std::set<string> basis_sets;
    for (auto &item : basis_families)
        basis_sets.insert(item.first);

    return basis_sets;
}

std::set<string> LI::availableBasisSetsAux()
{
    std::set<string> aux_basis_sets;
    for (auto &item : aux_basis_families)
        aux_basis_sets.insert(item.first);

    return aux_basis_sets;
}

string LI::returnAuxBasisFamilyString(const string &aux_basis_set)
{
    string basis_set_lc = aux_basis_set;
    std::transform(basis_set_lc.begin(), basis_set_lc.end(), basis_set_lc.begin(),
                   [](auto c)
                   { return std::tolower(c); });

    if (aux_basis_families.find(basis_set_lc) == aux_basis_families.end())
        throw std::runtime_error(fmt::format("The requested auxiliary basis set {} could not be found!",
                                             aux_basis_set));
    else
    {
        switch (aux_basis_families.at(basis_set_lc))
        {
        case AuxBasisFamily::ahlrichs_fit:
            return "ahlrichs";
        default:
            throw std::runtime_error("");
        }
    }
}

string LI::returnBasisFamilyString(const string &basis_set)
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
        case BasisFamily::sto:
            return "sto";
        default:
            throw std::runtime_error("");
        }
    }
}

LI::AuxBasisFamily LI::returnAuxBasisFamily(const std::string &aux_basis_set)
{
    return aux_basis_families.at(aux_basis_set);
}

LI::BasisFamily LI::returnBasisFamily(const std::string &basis_set)
{
    return basis_families.at(basis_set);
}