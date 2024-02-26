#include <lible/ints_defs.hpp>
#include <lible/ints_util.hpp>
#include <lible/structure.hpp>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <fmt/core.h>
#include <nlohmann/json.hpp>

namespace FS = std::filesystem;
namespace LI = lible::ints;

using std::array, std::string, std::vector;

using json = nlohmann::json;

LI::Structure::Structure(const string &basis_set,
                         const vector<double> &coordinates_angstroem,
                         const vector<string> &elements)
    : basis_set(basis_set), coordinates(coordinates_angstroem), elements(elements)
{
    n_atoms = elements.size();

    for (size_t i = 0; i < coordinates.size(); i++)
        coordinates[i] *= ang_to_bohr;

    readBasis(basis_set);
}

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#endif
string LI::Structure::returnBasisPath(const string &basis_set)
{
    string bs = basis_set;
    transform(bs.begin(), bs.end(), bs.begin(),
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
        {
            return basis_path;
        }
    }

    // TODO: do some kind of exception handling here, e.g. with assert
    // assert(not(fmt::format("Basis set {} in basis-library.", bs), basis_found));

    return basis_path;
}

vector<LI::Shell> LI::Structure::parseBasisJSONFile(const string &basis_path)
{
    json basis_set_json = json::parse(std::ifstream{basis_path});

    max_angular_momentum = 0;

    size_t pos = 0;
    vector<Shell> shells;
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        string element = elements[iatom];
        string atomic_number = std::to_string(atomic_numbers.at(element));
        auto shells_json = basis_set_json.at("elements").at(atomic_number).at("electron_shells");
        for (auto &shell_json : shells_json)
        {
            // TODO: sometimes there are multiple sets of coefficients for one set of exponents.
            // Implement a more flexible scheme that can deal with that, for now we only get the
            // first ones -- this can lead to big problems!

            int angular_momentum = shell_json.at("angular_momentum")[0];
            assert((angular_momentum >= 0));

            if (angular_momentum > max_angular_momentum)
                max_angular_momentum = angular_momentum;

            int atomic_number = atomic_numbers.at(element);
            size_t dim_cartesian = dimCartesians(angular_momentum);
            size_t dim_spherical = dimSphericals(angular_momentum);
            array<double, 3> xyz_coordinates{coordinates[3 * iatom],
                                             coordinates[3 * iatom + 1],
                                             coordinates[3 * iatom + 2]};

            vector<double> coeffs;
            for (const string &coeff : shell_json.at("coefficients")[0])
                coeffs.push_back(std::stod(coeff));

            vector<double> exps;
            for (const string &exp : shell_json.at("exponents"))
                exps.push_back(std::stod(exp));

            vector<double> coeffs_raw = coeffs;
            for (size_t i = 0; i < exps.size(); i++)
                coeffs[i] *= calcPurePrimitiveNorm(angular_momentum,
                                                      exps[i]);

            vector<array<int, 3>> cartesian_exps = returnCartesianExps(angular_momentum);

            vector<double> norms = calcShellNormalization(angular_momentum,
                                                          coeffs,
                                                          exps);

            Shell shell(angular_momentum, atomic_number,
                        dim_cartesian, dim_spherical,
                        pos, xyz_coordinates,
                        coeffs, coeffs_raw,
                        exps, norms, 
                        cartesian_exps);

            shells.push_back(shell);
            pos += dim_spherical;
        }
    }

    n_atomic_orbitals = pos;

    return shells;
}

void LI::Structure::readBasis(const string &basis_set)
{
    string basis_path = returnBasisPath(basis_set);
    vector<Shell> shells_vec = parseBasisJSONFile(basis_path);

    for (const auto &shell : shells_vec)
        shells[shell.angular_momentum].push_back(shell);

    size_t dim_ao = 0;
    for (auto &shell : shells_vec)
        dim_ao += shell.dim_spherical;
    printf("dim_ao = %d\n", dim_ao); // TMP

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