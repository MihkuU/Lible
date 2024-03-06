#include <lible/structure.hpp>
#include <lible/ints.hpp>
#include <lible/ints_defs.hpp>
#include <lible/ints_util.hpp>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include <fmt/core.h>

namespace LI = lible::ints;

using std::array, std::map, std::set, std::string, std::vector;

LI::Structure::Structure(const string &basis_set,
                         const vector<int> &atomic_nrs,
                         const vector<double> &coordinates_angstroem)
    : basis_set(basis_set), atomic_nrs(atomic_nrs), coordinates(coordinates_angstroem)
{
    n_atoms = atomic_nrs.size();

    elements.resize(n_atoms);
    for (size_t iatom = 0; iatom < atomic_nrs.size(); iatom)
        elements[iatom] = atomic_symbols.at(atomic_nrs[iatom]);

    for (size_t i = 0; i < coordinates.size(); i++)
        coordinates[i] *= ang_to_bohr;

    readBasis(basis_set, shells);
}

vector<LI::Shell> LI::Structure::parseBasisJSONFile(const string &basis_path)
{
    // json basis_set_json = json::parse(std::ifstream{basis_path});

    // max_angular_momentum = 0;

    // size_t pos = 0;
    // vector<Shell> shells;
    // for (size_t iatom = 0; iatom < n_atoms; iatom++)
    // {
    //     string element = elements[iatom];
    //     string atomic_number = std::to_string(atomic_numbers.at(element));
    //     auto shells_json = basis_set_json.at("elements").at(atomic_number).at("electron_shells");        

    //     // std::cout << shells_json << std::endl;
    //     // printf("\nshells_json.size() = %d\n", shells_json.size());

    //     for (auto &shell_json : shells_json)
    //     {
    //         // TODO: sometimes there are multiple sets of coefficients for one set of exponents.
    //         // Implement a more flexible scheme that can deal with that, for now we only get the
    //         // first ones -- this can lead to big problems!
    //         // TODO: handle pople basis sets.

    //         int angular_momentum = shell_json.at("angular_momentum")[0]; // TODO: above^

    //         size_t n_exps = shell_json.at("exponents").size();
    //         size_t n_coeff_sets = shell_json.at("coefficients").size();
    //         printf("n_coeff_sets = %d\n", n_coeff_sets);

    //         for (auto &coeffs : shell_json.at("coefficients"))
    //             printf("coeffs.size() = %d\n", coeffs.size());

            

    //         assert((angular_momentum >= 0));

    //         if (angular_momentum > max_angular_momentum)
    //             max_angular_momentum = angular_momentum;

    //         int atomic_number = atomic_numbers.at(element);
    //         size_t dim_cartesian = dimCartesians(angular_momentum);
    //         size_t dim_spherical = dimSphericals(angular_momentum);
    //         array<double, 3> xyz_coordinates{coordinates[3 * iatom],
    //                                          coordinates[3 * iatom + 1],
    //                                          coordinates[3 * iatom + 2]};

    //         vector<double> coeffs;
    //         for (const string &coeff : shell_json.at("coefficients")[0]) // TODO: above^
    //             coeffs.push_back(std::stod(coeff));

    //         vector<double> exps;
    //         for (const string &exp : shell_json.at("exponents"))
    //             exps.push_back(std::stod(exp));

    //         vector<double> coeffs_raw = coeffs;
    //         for (size_t i = 0; i < exps.size(); i++)
    //             coeffs[i] *= calcPurePrimitiveNorm(angular_momentum,
    //                                                   exps[i]);

    //         vector<array<int, 3>> cartesian_exps = returnCartesianExps(angular_momentum);

    //         vector<double> norms = calcShellNormalization(angular_momentum,
    //                                                       coeffs,
    //                                                       exps);

    //         Shell shell(angular_momentum, atomic_number,
    //                     dim_cartesian, dim_spherical,
    //                     pos, xyz_coordinates,
    //                     coeffs, coeffs_raw,
    //                     exps, norms, 
    //                     cartesian_exps);

    //         shells.push_back(shell);
    //         pos += dim_spherical;
    //     }
    // }

    // n_atomic_orbitals = pos;

    // return shells;
}

void LI::Structure::constructShells(map<int, vector<Shell>> &shells)
{
    set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());

    auto basis_atoms = returnBasisForAtoms(atomic_nrs_set, basis_set);

    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        int atomic_nr = atomic_nrs[iatom];

        auto &basis_atom = basis_atoms.at(atomic_nr);


    }
}

void LI::Structure::readBasis(const string &basis_set, map<int, vector<Shell>> &shells)
{
    string basis_path = returnBasisPath(basis_set);
    constructShells(shells);
    // vector<Shell> shells_vec = parseBasisJSONFile(basis_path);

    // for (const auto &shell : shells_vec)
    //     shells[shell.angular_momentum].push_back(shell);

    // size_t dim_ao = 0;
    // for (auto &shell : shells_vec)
    //     dim_ao += shell.dim_spherical;
    // printf("dim_ao = %d\n", dim_ao); // TMP

    // for (int la = max_angular_momentum; la >= 0; la--)
    //     for (int lb = la; lb >= 0; lb--)
    //         if (la == lb)
    //         {
    //             for (size_t ishell = 0; ishell < shells.at(la).size(); ishell++)
    //                 for (size_t jshell = ishell; jshell < shells.at(lb).size(); jshell++)
    //                     shell_pairs[std::make_pair(la, lb)].push_back(
    //                         shells::ShellPair(shells.at(la)[ishell], shells.at(lb)[jshell]));
    //         }
    //         else
    //         {
    //             for (const auto &shell_a : shells.at(la))
    //                 for (const auto &shell_b : shells.at(lb))
    //                     shell_pairs[std::make_pair(la, lb)].push_back(shells::ShellPair(shell_a, shell_b));
    //         }
}