#include <lible/ints/basis_sets.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/utils.hpp>

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace FS = std::filesystem;
namespace LI = lible::ints;

using std::pair, std::string, std::vector;

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#define path_to_aux_basis_sets AUX_BASIS_DIR
#endif

namespace lible::ints
{
    static std::map<int, vector<LI::shell_exps_coeffs_t>>
    basisForAtomImpl(const int atomic_nr, const string &basis_set, const string &basis_path)
    {
        std::ifstream basis_file(basis_path, std::ios::in);

        bool basis_found{false}, read_lines{false};

        string line;
        vector<string> lines;
        while (std::getline(basis_file, line))
        {
            if (line == std::format("end={}", atomic_nr))
                break;

            if (read_lines)
                lines.push_back(line);

            if (line == std::format("element={}", atomic_nr))
            {
                read_lines = true;
                basis_found = true;
            }
        }

        if (!basis_found)
        {
            string msg = std::format("Basis set {} not found for element {}!",
                                     basis_set, atomic_symbols.at(atomic_nr));
            throw std::runtime_error(msg);
        }

        int l, counter = 0;
        std::map<int, vector<shell_exps_coeffs_t>> basis_atom;
        shell_exps_coeffs_t exps_coeffs;
        for (string &line : lines)
        {
            if (line.find("l=") != string::npos)
            {
                if (counter != 0)
                    basis_atom[l].push_back(exps_coeffs);

                l = line[2] - '0';
                exps_coeffs.first.clear();
                exps_coeffs.second.clear();
                counter = 0;
            }
            else
            {
                double exp, coeff;
                std::stringstream line_ss;
                line_ss << line;

                line_ss >> exp;
                line_ss >> coeff;
                exps_coeffs.first.push_back(exp);
                exps_coeffs.second.push_back(coeff);
                counter++;
            }
        }
        basis_atom[l].push_back(exps_coeffs);

        return basis_atom;
    }
}

std::map<int, vector<LI::shell_exps_coeffs_t>>
LI::basisForAtom(const int atomic_nr, const string &basis_set)
{
    string basis_path = returnBasisPath(basis_set);
    return basisForAtomImpl(atomic_nr, basis_set, basis_path);
}

std::map<int, vector<LI::shell_exps_coeffs_t>>
LI::basisForAtomAux(const int atomic_nr, const string &basis_set)
{
    string aux_basis_path = returnAuxBasisPath(basis_set);
    return basisForAtomImpl(atomic_nr, basis_set, aux_basis_path);
}

std::map<int, std::map<int, vector<LI::shell_exps_coeffs_t>>>
LI::basisForAtoms(const std::set<int> &atomic_nrs, const string &basis_set)
{
    std::map<int, std::map<int, vector<shell_exps_coeffs_t>>> basis_atoms;
    for (const int atomic_nr : atomic_nrs)
    {
        auto basis_atom = basisForAtom(atomic_nr, basis_set);
        basis_atoms[atomic_nr] = basis_atom;
    }

    return basis_atoms;
}

std::map<int, std::map<int, vector<LI::shell_exps_coeffs_t>>>
LI::basisForAtomsAux(const std::set<int> &atomic_nrs, const string &aux_basis_set)
{
    std::map<int, std::map<int, vector<shell_exps_coeffs_t>>> aux_basis_atoms;
    for (const int atomic_nr : atomic_nrs)
    {
        auto basis_atom = basisForAtomAux(atomic_nr, aux_basis_set);
        aux_basis_atoms[atomic_nr] = basis_atom;
    }

    return aux_basis_atoms;
}

string LI::returnAuxBasisPath(const string &aux_basis_set)
{
    string bs = aux_basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(), [](unsigned char c)
                   { return std::tolower(c); });

    string basis_family_str = returnAuxBasisFamilyString(aux_basis_set);
    string basis_prefix = string(path_to_aux_basis_sets) + "/" + basis_family_str;

    string basis_path;
    for (const auto &entry : FS::directory_iterator(basis_prefix))
    {
        basis_path = entry.path();
        string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find("."));
        if (basis_name == bs)
            return basis_path;
    }

    string message = std::format("The requested basis set {} could not be found!", aux_basis_set);
    throw std::runtime_error(message);
}

string LI::returnBasisPath(const string &basis_set)
{
    string bs = basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(), [](unsigned char c)
                   { return std::tolower(c); });

    string basis_family_str = returnBasisFamilyString(basis_set);
    string basis_prefix = string(path_to_basis_sets) + "/" + basis_family_str;

    string basis_path;
    for (const auto &entry : FS::directory_iterator(basis_prefix))
    {
        basis_path = entry.path();
        string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find("."));
        if (basis_name == bs)
            return basis_path;
    }

    string message = std::format("The requested basis set {} could not be found!", basis_set);
    throw std::runtime_error(message);
}