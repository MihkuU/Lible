#include <lible/ints.hpp>
#include <lible/ints_util.hpp>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <fmt/core.h>
#include <nlohmann/json.hpp>

namespace FS = std::filesystem;
namespace LI = lible::ints;

using std::map, std::pair, std::set, std::string, std::vector;

using json = nlohmann::json;

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#endif

map<int, vector<LI::shell_exps_coeffs_t>>
LI::returnBasisForAtom(const int atomic_nr, const string &basis_set)
{
    /*
     *
     *
     */

    string basis_path = returnBasisPath(basis_set);
    json basis_json = json::parse(std::ifstream{basis_path});

    string atomic_nr_str = std::to_string(atomic_nr);
    auto shells_json = basis_json.at("elements").at(atomic_nr_str).at("electron_shells");

    map<int, vector<shell_exps_coeffs_t>> basis_atom;
    for (auto &shell : shells_json)
    {
        int angmom = shell.at("angular_momentum")[0];

        assert((shell.at("angular_momentum").size() == 1));

        size_t n_exps = shell.at("exponents").size();
        size_t n_coeff_sets = shell.at("coefficients").size();

        auto exps_json = shell.at("exponents");
        auto coeffs_json = shell.at("coefficients");

        vector<shell_exps_coeffs_t> exp_coeff_pairs(n_coeff_sets);
        for (size_t icoeffs = 0; icoeffs < n_coeff_sets; icoeffs++)
        {
            vector<double> exps;
            vector<double> coeffs;
            for (size_t iexp = 0; iexp < n_exps; iexp++)
                if (coeffs_json[icoeffs][iexp] != 0)
                {
                    exps.push_back(exps_json[iexp]);
                    coeffs.push_back(coeffs_json[icoeffs][iexp]);
                }

            exp_coeff_pairs[icoeffs] = std::make_pair(exps, coeffs);
        }

        basis_atom[angmom].insert(basis_atom[angmom].end(), exp_coeff_pairs.begin(),
                                  exp_coeff_pairs.end());
    }

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
    for (const auto &entry : FS::directory_iterator(path_to_basis_sets))
    {
        basis_path = entry.path();
        string basis_name = entry.path().filename();
        basis_name = basis_name.substr(0, basis_name.find("."));
        if (basis_name == bs)        
            return basis_path;        
    }

    string message = fmt::format("The requested basis set {} could not be found!", basis_set);
    throw std::runtime_error(message);
}