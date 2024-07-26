#include <lible/ints/basis_sets.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/util.hpp>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <fmt/core.h>

namespace FS = std::filesystem;
namespace LI = lible::ints;

using std::map, std::pair, std::set, std::string, std::vector;

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#endif

map<int, vector<LI::shell_exps_coeffs_t>>
LI::returnBasisForAtom(const int atomic_nr, const string &basis_set)
{
    string basis_path = returnBasisPath(basis_set);

    std::ifstream basis_file(basis_path, std::ios::in);

    bool basis_found = false, read_lines = false;
            
    string line;
    vector<string> lines;
    while (std::getline(basis_file, line))
    {
        if (line == fmt::format("end={}", atomic_nr))
            break;

        if (read_lines)
            lines.push_back(line);

        if (line == fmt::format("element={}", atomic_nr))
        {
            read_lines = true;
            basis_found = true;
        }
    }

    if (!basis_found)
    {
        string msg = fmt::format("Basis set {} not found for element {}!",
                                 basis_set, atomic_symbols.at(atomic_nr));
        throw std::runtime_error(msg);
    }

    int l, counter = 0;
    map<int, vector<shell_exps_coeffs_t>> basis_atom;
    shell_exps_coeffs_t shell;
    for (string &line : lines)
    {
        if (line.find("l=") != string::npos)
        {
            if (counter != 0)
                basis_atom[l].push_back(shell);

            l = line[2] - '0';
            shell.first.clear();
            shell.second.clear();
            counter = 0;
        }
        else
        {
            double exp, coeff;
            std::stringstream line_ss;
            line_ss << line;

            line_ss >> exp;
            line_ss >> coeff;
            shell.first.push_back(exp);
            shell.second.push_back(coeff);
            counter++;
        }
    }
    basis_atom[l].push_back(shell);

    return basis_atom;
}

map<int, map<int, vector<LI::shell_exps_coeffs_t>>>
LI::returnBasisForAtoms(const set<int> &atomic_nrs, const string &basis_set)
{
    map<int, map<int, vector<shell_exps_coeffs_t>>> basis_atoms;
    for (const int atomic_nr : atomic_nrs)
    {
        auto basis_atom = returnBasisForAtom(atomic_nr, basis_set);
        basis_atoms[atomic_nr] = basis_atom;
    }

    return basis_atoms;
}

string LI::returnBasisPath(const string &basis_set)
{
    string bs = basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    bool basis_found = false;
    string basis_path;

    string basis_family = returnBasisFamily(basis_set);
    string basis_prefix = string(path_to_basis_sets) + "/" + basis_family;

    for (const auto &entry : FS::directory_iterator(basis_prefix))
    {
        basis_path = entry.path();
        string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find("."));
        if (basis_name == bs)
            return basis_path;
    }

    string message = fmt::format("The requested basis set {} could not be found!",
                                 basis_set);

    throw std::runtime_error(message);
}