#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <fmt/core.h>
#include <nlohmann/json.hpp>
#include "defs_ints.h"
#include "structure.h"

namespace fs = std::filesystem;

using namespace lible;
using namespace lible::ints;

using std::array;
using std::string;
using std::vector;

using json = nlohmann::json;

Structure::Structure(const string &basis_set,
                     const vector<double> &coordinates_angstroem,
                     const vector<string> &elements)
    : basis_set(basis_set), coordinates(coordinates_angstroem), elements(elements)
{
    n_atoms = elements.size();

    for (size_t i = 0; i < coordinates.size(); i++)
        coordinates[i] *= IntsDefs::ang_to_bohr;

    readBasis(basis_set);
}

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#endif
string Structure::returnBasisPath(const string &basis_set)
{
    string bs = basis_set;
    transform(bs.begin(), bs.end(), bs.begin(),
              [](unsigned char c)
              { return std::tolower(c); });

    bool basis_found = false;
    string basis_path;
    for (const auto &entry : fs::directory_iterator(path_to_basis_sets))
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

vector<Shells::Shell> Structure::parseBasisJSONFile(const string &basis_path)
{
    json basis_set_json = json::parse(std::ifstream{basis_path});

    max_angular_momentum = 0;

    size_t pos = 0;
    vector<Shells::Shell> shells;
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        string element = elements[iatom];
        string atomic_number = std::to_string(IntsDefs::atomic_numbers.at(element));
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

            int atomic_number = IntsDefs::atomic_numbers.at(element);
            size_t dim_cartesian = Shells::calcShellDimCartesian(angular_momentum);
            size_t dim_spherical = Shells::calcShellDimSpherical(angular_momentum);
            array<double, 3> xyz_coordinates{coordinates[3 * iatom], coordinates[3 * iatom + 1], coordinates[3 * iatom + 2]};

            vector<double> contraction_coeffs;
            for (const string &coeff : shell_json.at("coefficients")[0])
                contraction_coeffs.push_back(std::stod(coeff));

            vector<double> contraction_exps;
            for (const string &exp : shell_json.at("exponents"))
                contraction_exps.push_back(std::stod(exp));

            for (size_t i = 0; i < contraction_exps.size(); i++)
                contraction_coeffs[i] *= Shells::calcPureGaussianPrimitiveNorm(angular_momentum, contraction_exps[i]);

            vector<array<int, 3>> cartesian_exps = Shells::calcShellCartesianExps(angular_momentum);

            vector<double> normalization = Shells::calcShellNormalization(angular_momentum, contraction_coeffs, contraction_exps,
                                                                          cartesian_exps);

            Shells::Shell shell(angular_momentum, atomic_number, dim_cartesian, dim_spherical, pos, xyz_coordinates,
                                contraction_coeffs, contraction_exps, normalization, cartesian_exps);
            shells.push_back(shell);
            pos += dim_spherical;
        }
    }

    n_atomic_orbitals = pos;

    return shells;
}

void ints::Structure::readBasis(const string &basis_set)
{
    string basis_path = returnBasisPath(basis_set);
    vector<Shells::Shell> shells_vec = parseBasisJSONFile(basis_path);

    for (const auto &shell : shells_vec)
        shells[shell.angular_momentum].push_back(shell);

    for (int la = max_angular_momentum; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
            if (la == lb)
            {
                for (size_t ishell = 0; ishell < shells.at(la).size(); ishell++)
                    for (size_t jshell = ishell; jshell < shells.at(lb).size(); jshell++)
                        shell_pairs[std::make_pair(la, lb)].push_back(
                            Shells::ShellPair(shells.at(la)[ishell], shells.at(lb)[jshell]));
            }
            else
            {
                for (const auto &shell_a : shells.at(la))
                    for (const auto &shell_b : shells.at(lb))
                        shell_pairs[std::make_pair(la, lb)].push_back(Shells::ShellPair(shell_a, shell_b));
            }
}