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

using namespace Lible;
using json = nlohmann::json;

namespace fs = std::filesystem;

Ints::Structure::Structure(const std::string &basis_set, const std::vector<double> &coordinates, const std::vector<std::string> &elements)
    : basis_set(basis_set), coordinates(coordinates), elements(elements)
{
    n_atoms = elements.size();
    readBasis(basis_set);
}

#ifdef BASIS_DIR
#define path_to_basis_sets BASIS_DIR
#endif
std::string Ints::Structure::returnBasisPath(const std::string &basis_set)
{
    std::string bs = basis_set;
    std::transform(bs.begin(), bs.end(), bs.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    bool basis_found = false;
    std::string basis_path;
    for (const auto &entry : fs::directory_iterator(path_to_basis_sets))
    {
        basis_path = entry.path();
        std::string basis_name = entry.path().filename();
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

std::vector<Shells::Shell> Ints::Structure::parseBasisJSONFile(const std::string &basis_path)
{
    json basis_set_json = json::parse(std::ifstream{basis_path});

    std::size_t pos = 0;
    std::vector<Shells::Shell> shells;
    for (std::size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        std::string element = elements[iatom];
        std::cout << element << std::endl;
        std::string atomic_number = std::to_string(IntsDefs::atomic_numbers.at(element));
        auto shells_json = basis_set_json.at("elements").at(atomic_number).at("electron_shells");
        for (auto &shell_json : shells_json)
        {
            // TODO: sometimes there are multiple sets of coefficients for one set of exponents.
            // Implement a more flexible scheme that can deal with that, for now we only get the
            // first ones -- this can lead to big problems!

            int angular_momentum = shell_json.at("angular_momentum")[0];
            int atomic_number = IntsDefs::atomic_numbers.at(element);
            std::size_t dim_cartesian = Shells::calcShellDimCartesian(angular_momentum);
            std::size_t dim_spherical = Shells::calcShellDimSpherical(angular_momentum);
            std::array<double, 3> xyz_coordinates{coordinates[3 * iatom], coordinates[3 * iatom + 1], coordinates[3 * iatom + 2]};

            std::vector<double> coefficients;
            for (const std::string &coeff : shell_json.at("coefficients")[0])
                coefficients.push_back(std::stod(coeff));

            std::vector<double> exponents;
            for (const std::string &exp : shell_json.at("exponents"))
                exponents.push_back(std::stod(exp));

            for (std::size_t i = 0; i < exponents.size(); i++)
                coefficients[i] *= Shells::calcPureGaussianPrimitiveNorm(angular_momentum, exponents[i]);

            std::vector<std::array<int, 3>> cartesian_exps = Shells::calcShellCartesianExps(angular_momentum);

            std::vector<double> normalization; // TODO
            Shells::Shell shell(angular_momentum, atomic_number, dim_cartesian, dim_spherical, pos, xyz_coordinates,
                                coefficients, exponents, normalization, cartesian_exps);
            shells.push_back(shell);
            pos += dim_spherical;
        }
    }

    return shells;
}

void Ints::Structure::readBasis(const std::string &basis_set)
{
    std::string basis_path = returnBasisPath(basis_set);
    std::vector<Shells::Shell> shells = parseBasisJSONFile(basis_path);
}