#include <lible/ints/basis_sets.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/ints.hpp>

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;
namespace lints = lible::ints;

namespace lible::ints
{
    /// Backend function for reading and returning a basis set from the given path.
    BasisAtom basisForAtomImpl(int atomic_nr, const std::string &basis_set,
                               const std::string &basis_path);
}

/// Mapping between auxiliary basis set names and their families.
const static std::unordered_map<std::string, lints::AuxBasisFamily>
aux_basis_families{
    {"def2-qzvp-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-qzvpp-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-qzvppd-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-sv(p)-jkfit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-sv(p)-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-svp-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-svpd-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-tzvp-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-tzvpd-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-tzvpp-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-tzvppd-rifit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-universal-jfit", lints::AuxBasisFamily::ahlrichs_fit},
    {"def2-universal-jkfit", lints::AuxBasisFamily::ahlrichs_fit}
};

/// Mapping between main basis set names and their families.
const static std::unordered_map<std::string, lints::BasisFamily>
basis_families{
    {"def2-svp", lints::BasisFamily::ahlrichs},
    {"def2-tzvp", lints::BasisFamily::ahlrichs},
    {"ahlrichs-pvdz", lints::BasisFamily::ahlrichs},
    {"ahlrichs-tzv", lints::BasisFamily::ahlrichs},
    {"ahlrichs-vdz", lints::BasisFamily::ahlrichs},
    {"ahlrichs-vtz", lints::BasisFamily::ahlrichs},
    {"def2-qzvpd", lints::BasisFamily::ahlrichs},
    {"def2-qzvppd", lints::BasisFamily::ahlrichs},
    {"def2-qzvpp", lints::BasisFamily::ahlrichs},
    {"def2-qzvp", lints::BasisFamily::ahlrichs},
    {"def2-svpd", lints::BasisFamily::ahlrichs},
    {"def2-sv(p)", lints::BasisFamily::ahlrichs},
    {"def2-svp", lints::BasisFamily::ahlrichs},
    {"def2-tzvpd", lints::BasisFamily::ahlrichs},
    {"def2-tzvppd", lints::BasisFamily::ahlrichs},
    {"def2-tzvpp", lints::BasisFamily::ahlrichs},
    {"def2-tzvp", lints::BasisFamily::ahlrichs},
    {"ano-dk3", lints::BasisFamily::ano},
    {"ano-pv5z", lints::BasisFamily::ano},
    {"ano-pvdz", lints::BasisFamily::ano},
    {"ano-pvqz", lints::BasisFamily::ano},
    {"ano-pvtz", lints::BasisFamily::ano},
    {"ano-r", lints::BasisFamily::ano},
    {"ano-r0", lints::BasisFamily::ano},
    {"ano-r1", lints::BasisFamily::ano},
    {"ano-r2", lints::BasisFamily::ano},
    {"ano-r3", lints::BasisFamily::ano},
    {"ano-rcc-mb", lints::BasisFamily::ano},
    {"ano-rcc-vdz", lints::BasisFamily::ano},
    {"ano-rcc-vdzp", lints::BasisFamily::ano},
    {"ano-rcc-vqzp", lints::BasisFamily::ano},
    {"ano-rcc-vtz", lints::BasisFamily::ano},
    {"ano-rcc-vtzp", lints::BasisFamily::ano},
    {"ano-rcc", lints::BasisFamily::ano},
    {"aug-ano-pv5z", lints::BasisFamily::ano},
    {"aug-ano-pvdz", lints::BasisFamily::ano},
    {"aug-ano-pvqz", lints::BasisFamily::ano},
    {"aug-ano-pvtz", lints::BasisFamily::ano},
    {"nasa-ames-ano", lints::BasisFamily::ano},
    {"nasa-ames-ano2", lints::BasisFamily::ano},
    {"roos-augmented-double-zeta-ano", lints::BasisFamily::ano},
    {"roos-augmented-triple-zeta-ano", lints::BasisFamily::ano},
    {"saug-ano-pv5z", lints::BasisFamily::ano},
    {"saug-ano-pvdz", lints::BasisFamily::ano},
    {"saug-ano-pvqz", lints::BasisFamily::ano},
    {"saug-ano-pvtz", lints::BasisFamily::ano},
    {"aug-cc-pcv5z", lints::BasisFamily::dunning},
    {"aug-cc-pv7z", lints::BasisFamily::dunning},
    {"aug-cc-pwcv5z", lints::BasisFamily::dunning},
    {"aug-seg-cc-pvqz-pp", lints::BasisFamily::dunning},
    {"cc-pcvdz", lints::BasisFamily::dunning},
    {"cc-pv9z", lints::BasisFamily::dunning},
    {"cc-pv(t+d)z", lints::BasisFamily::dunning},
    {"d-aug-cc-pv5z", lints::BasisFamily::dunning},
    {"seg-cc-pv5z-pp", lints::BasisFamily::dunning},
    {"seg-cc-pwcvtz-pp", lints::BasisFamily::dunning},
    {"aug-cc-pcvdz", lints::BasisFamily::dunning},
    {"aug-cc-pv(d+d)z", lints::BasisFamily::dunning},
    {"aug-cc-pwcvdz", lints::BasisFamily::dunning},
    {"aug-seg-cc-pvtz-pp", lints::BasisFamily::dunning},
    {"cc-pcvqz", lints::BasisFamily::dunning},
    {"cc-pv(d+d)z", lints::BasisFamily::dunning},
    {"cc-pvtz(seg-opt)", lints::BasisFamily::dunning},
    {"d-aug-cc-pv6z", lints::BasisFamily::dunning},
    {"seg-cc-pvdz-pp", lints::BasisFamily::dunning},
    {"aug-cc-pcvqz", lints::BasisFamily::dunning},
    {"aug-cc-pvdz", lints::BasisFamily::dunning},
    {"aug-cc-pwcvqz", lints::BasisFamily::dunning},
    {"aug-seg-cc-pwcv5z-pp", lints::BasisFamily::dunning},
    {"cc-pcvtz", lints::BasisFamily::dunning},
    {"cc-pvdz(seg-opt)", lints::BasisFamily::dunning},
    {"cc-pvtz", lints::BasisFamily::dunning},
    {"d-aug-cc-pvdz", lints::BasisFamily::dunning},
    {"seg-cc-pvqz-pp", lints::BasisFamily::dunning},
    {"aug-cc-pcvtz", lints::BasisFamily::dunning},
    {"aug-cc-pv(q+d)z", lints::BasisFamily::dunning},
    {"aug-cc-pwcvtz", lints::BasisFamily::dunning},
    {"aug-seg-cc-pwcvdz-pp", lints::BasisFamily::dunning},
    {"cc-pv(5+d)z", lints::BasisFamily::dunning},
    {"cc-pvdz", lints::BasisFamily::dunning},
    {"cc-pwcv5z", lints::BasisFamily::dunning},
    {"d-aug-cc-pvqz", lints::BasisFamily::dunning},
    {"seg-cc-pvtz-pp", lints::BasisFamily::dunning},
    {"aug-cc-pv(5+d)z", lints::BasisFamily::dunning},
    {"aug-cc-pvqz", lints::BasisFamily::dunning},
    {"aug-pv7z", lints::BasisFamily::dunning},
    {"aug-seg-cc-pwcvqz-pp", lints::BasisFamily::dunning},
    {"cc-pv5z", lints::BasisFamily::dunning},
    {"cc-pv(q+d)z", lints::BasisFamily::dunning},
    {"cc-pwcvdz", lints::BasisFamily::dunning},
    {"d-aug-cc-pvtz", lints::BasisFamily::dunning},
    {"seg-cc-pwcv5z-pp", lints::BasisFamily::dunning},
    {"aug-cc-pv5z", lints::BasisFamily::dunning},
    {"aug-cc-pv(t+d)z", lints::BasisFamily::dunning},
    {"aug-seg-cc-pv5z-pp", lints::BasisFamily::dunning},
    {"aug-seg-cc-pwcvtz-pp", lints::BasisFamily::dunning},
    {"cc-pv6z", lints::BasisFamily::dunning},
    {"cc-pvqz(seg-opt)", lints::BasisFamily::dunning},
    {"cc-pwcvqz", lints::BasisFamily::dunning},
    {"pv6z", lints::BasisFamily::dunning},
    {"seg-cc-pwcvdz-pp", lints::BasisFamily::dunning},
    {"aug-cc-pv6z", lints::BasisFamily::dunning},
    {"aug-cc-pvtz", lints::BasisFamily::dunning},
    {"aug-seg-cc-pvdz-pp", lints::BasisFamily::dunning},
    {"cc-pcv5z", lints::BasisFamily::dunning},
    {"cc-pv8z", lints::BasisFamily::dunning},
    {"cc-pvqz", lints::BasisFamily::dunning},
    {"cc-pwcvtz", lints::BasisFamily::dunning},
    {"pv7z", lints::BasisFamily::dunning},
    {"seg-cc-pwcvqz-pp", lints::BasisFamily::dunning},
    {"3-21g", lints::BasisFamily::pople},
    {"6-21g", lints::BasisFamily::pople},
    {"6-311+g(2d,p)", lints::BasisFamily::pople},
    {"6-311+g**", lints::BasisFamily::pople},
    {"6-311g**", lints::BasisFamily::pople},
    {"6-311++g**", lints::BasisFamily::pople},
    {"6-31g(2df,p)", lints::BasisFamily::pople},
    {"6-31+g**", lints::BasisFamily::pople},
    {"6-31g**", lints::BasisFamily::pople},
    {"6-31++g**", lints::BasisFamily::pople},
    {"4-31g", lints::BasisFamily::pople},
    {"6-311++g(2d,2p)", lints::BasisFamily::pople},
    {"6-311++g(3df,3pd)", lints::BasisFamily::pople},
    {"6-311+g*", lints::BasisFamily::pople},
    {"6-311g*", lints::BasisFamily::pople},
    {"6-311++g*", lints::BasisFamily::pople},
    {"6-31g(3df,3pd)", lints::BasisFamily::pople},
    {"6-31+g*", lints::BasisFamily::pople},
    {"6-31g*", lints::BasisFamily::pople},
    {"6-31++g*", lints::BasisFamily::pople},
    {"5-21g", lints::BasisFamily::pople},
    {"6-311g(2df,2pd)", lints::BasisFamily::pople},
    {"6-311g(d,p)", lints::BasisFamily::pople},
    {"6-311+g", lints::BasisFamily::pople},
    {"6-311g", lints::BasisFamily::pople},
    {"6-311++g", lints::BasisFamily::pople},
    {"6-31g(d,p)", lints::BasisFamily::pople},
    {"6-31+g", lints::BasisFamily::pople},
    {"6-31g", lints::BasisFamily::pople},
    {"6-31++g", lints::BasisFamily::pople},
    {"sto-2g", lints::BasisFamily::sto},
    {"sto-3g", lints::BasisFamily::sto},
    {"sto-3g*", lints::BasisFamily::sto},
    {"sto-4g", lints::BasisFamily::sto},
    {"sto-5g", lints::BasisFamily::sto},
    {"sto-6g", lints::BasisFamily::sto}
};

lints::BasisAtom lints::basisForAtomImpl(const int atomic_nr, const std::string &basis_set,
                                         const std::string &basis_path)
{
    std::ifstream basis_file(basis_path, std::ios::in);

    std::vector<std::string> lines;
    {
        bool basis_found{false}, read_lines{false};
        std::string line;
        while (std::getline(basis_file, line))
        {
            if (line == std::format("end={}", atomic_nr))
                break;

            if (read_lines == true)
                lines.push_back(line);

            if (line == std::format("element={}", atomic_nr))
            {
                read_lines = true;
                basis_found = true;
            }
        }

        if (basis_found == false)
        {
            std::string msg = std::format("Basis set {} not found for element {}!",
                                          basis_set, atomic_symbols.at(atomic_nr));
            throw std::runtime_error(msg);
        }
    }

    int l, counter = 0;
    BasisShell basis_shell;
    basis_shells_t basis_shells;
    for (const std::string &line : lines)
    {
        if (line.find("l=") != std::string::npos)
        {
            if (counter != 0)
                basis_shells.push_back(basis_shell);

            l = line[2] - '0';
            basis_shell.l_ = l;
            basis_shell.exps_.clear();
            basis_shell.coeffs_.clear();
            counter = 0;
        }
        else
        {
            double exp, coeff;
            std::stringstream line_ss;
            line_ss << line;

            line_ss >> exp;
            line_ss >> coeff;
            basis_shell.exps_.push_back(exp);
            basis_shell.coeffs_.push_back(coeff);
            counter++;
        }
    }
    basis_shells.push_back(basis_shell);

    return {atomic_nr, basis_shells};
}

lints::BasisAtom lints::basisForAtom(const int atomic_nr, const std::string &basis_set)
{
    std::string basis_path = basisPath(basis_set);
    return basisForAtomImpl(atomic_nr, basis_set, basis_path);
}

lints::BasisAtom lints::basisForAtomAux(const int atomic_nr, const std::string &basis_set)
{
    std::string aux_basis_path = auxBasisPath(basis_set);
    return basisForAtomImpl(atomic_nr, basis_set, aux_basis_path);
}

lints::basis_atoms_t lints::basisForAtoms(const std::vector<int> &atomic_nrs,
                                          const std::string &basis_set)
{
    basis_atoms_t basis_atoms;
    for (const int atomic_nr : atomic_nrs)
    {
        BasisAtom basis_atom = basisForAtom(atomic_nr, basis_set);
        basis_atoms.push_back(basis_atom);
    }

    return basis_atoms;
}

lints::basis_atoms_t lints::basisForAtomsAux(const std::vector<int> &atomic_nrs,
                                             const std::string &aux_basis_set)
{
    basis_atoms_t aux_basis_atoms;
    for (const int atomic_nr : atomic_nrs)
    {
        BasisAtom basis_atom = basisForAtomAux(atomic_nr, aux_basis_set);
        aux_basis_atoms.push_back(basis_atom);
    }

    return aux_basis_atoms;
}

std::string lints::auxBasisPath(const std::string &aux_basis_set)
{
    std::string bs = aux_basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(), [](unsigned char c)
    {
        return std::tolower(c);
    });

    std::string basis_family_str = auxBasisFamilyString(aux_basis_set);
    std::string basis_prefix = BasisPaths::getAuxBasisSetsPath() + "/" + basis_family_str;

    std::string basis_path;
    for (const auto &entry : fs::directory_iterator(basis_prefix))
    {
        basis_path = entry.path();
        std::string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find('.'));
        if (basis_name == bs)
            return basis_path;
    }

    std::string message = std::format("The requested basis set {} could not be found!",
                                      aux_basis_set);
    throw std::runtime_error(message);
}

std::string lints::basisPath(const std::string &basis_set)
{
    std::string bs = basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(), [](unsigned char c)
    {
        return std::tolower(c);
    });

    std::string basis_family_str = basisFamilyString(basis_set);
    std::string basis_prefix = BasisPaths::getMainBasisSetsPath() + "/" + basis_family_str;

    std::string basis_path;
    for (const auto &entry : fs::directory_iterator(basis_prefix))
    {
        basis_path = entry.path();
        std::string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find('.'));
        if (basis_name == bs)
            return basis_path;
    }

    std::string message = std::format("The requested basis set {} could not be found!", basis_set);
    throw std::runtime_error(message);
}

std::set<std::string> lints::availableBasisSets()
{
    std::set<std::string> basis_sets;
    for (const auto &item : basis_families)
        basis_sets.insert(item.first);

    return basis_sets;
}

std::set<std::string> lints::availableBasisSetsAux()
{
    std::set<std::string> aux_basis_sets;
    for (const auto &item : aux_basis_families)
        aux_basis_sets.insert(item.first);

    return aux_basis_sets;
}

std::string lints::auxBasisFamilyString(const std::string &aux_basis_set)
{
    std::string basis_set_lc = aux_basis_set;
    std::transform(basis_set_lc.begin(), basis_set_lc.end(), basis_set_lc.begin(),
                   [](auto c)
                   {
                       return std::tolower(c);
                   });

    if (aux_basis_families.find(basis_set_lc) == aux_basis_families.end())
        throw std::runtime_error(std::format("The requested auxiliary basis set {} could not be found!",
                                             aux_basis_set));

    switch (aux_basis_families.at(basis_set_lc))
    {
        case AuxBasisFamily::ahlrichs_fit :
            return "ahlrichs";
        default :
            throw std::runtime_error("");
    }
}

std::string lints::basisFamilyString(const std::string &basis_set)
{
    std::string basis_set_lc = basis_set;
    std::transform(basis_set_lc.begin(), basis_set_lc.end(), basis_set_lc.begin(),
                   [](auto c)
                   {
                       return std::tolower(c);
                   });

    if (basis_families.find(basis_set_lc) == basis_families.end())
        throw std::runtime_error(std::format("The requested basis set {} could not be found!",
                                             basis_set));

    switch (basis_families.at(basis_set_lc))
    {
        case BasisFamily::ahlrichs :
            return "ahlrichs";
        case BasisFamily::ano :
            return "ano";
        case BasisFamily::dunning :
            return "dunning";
        case BasisFamily::pople :
            return "pople";
        case BasisFamily::sto :
            return "sto";
        default :
            throw std::runtime_error("");
    }
}

lints::AuxBasisFamily lints::auxBasisFamily(const std::string &aux_basis_set)
{
    return aux_basis_families.at(aux_basis_set);
}

lints::BasisFamily lints::basisFamily(const std::string &basis_set)
{
    return basis_families.at(basis_set);
}
