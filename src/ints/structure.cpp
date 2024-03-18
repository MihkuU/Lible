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
        coordinates[i] *= _ang_to_bohr_;

    readBasis(basis_set, this->max_angular_momentum, this->n_atomic_orbitals, this->shells);
}

LI::Structure::Structure(const string &basis_set,
                         const vector<string> &elements,
                         const vector<double> &coordinates_angstroem)
    : basis_set(basis_set), elements(elements), coordinates(coordinates_angstroem)
{
    n_atoms = elements.size();

    atomic_nrs.resize(n_atoms);
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        string atomic_symbol = elements[iatom];
        std::transform(atomic_symbol.begin(), atomic_symbol.end(), atomic_symbol.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });
                       
        atomic_symbol[0] = std::toupper(atomic_symbol[0]);
        this->elements[iatom] = atomic_symbol;

        atomic_nrs[iatom] = atomic_numbers.at(atomic_symbol);
    }

    for (size_t i = 0; i < coordinates.size(); i++)
        coordinates[i] *= _ang_to_bohr_;

    readBasis(basis_set, this->max_angular_momentum, this->n_atomic_orbitals, this->shells);
}

void LI::Structure::constructShells(int &max_angular_momentum,
                                    size_t &n_atomic_orbitals,
                                    map<int, vector<Shell>> &shells)
{
    set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());

    auto basis_atoms = returnBasisForAtoms(atomic_nrs_set, basis_set);

    int max_angmom = 0, pos = 0;
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        int atomic_nr = atomic_nrs[iatom];

        auto basis_atom = basis_atoms.at(atomic_nr);

        array<double, 3> xyz_coordinates{coordinates[3 * iatom],
                                         coordinates[3 * iatom + 1],
                                         coordinates[3 * iatom + 2]};

        for (const auto &[angmom, shell_exps_coeffs] : basis_atom)
        {
            if (angmom > max_angmom)
                max_angmom = angmom;

            assert((max_angmom <= _max_angular_momentum_));

            size_t dim_cart = dimCartesians(angmom);
            size_t dim_sphe = dimSphericals(angmom);

            for (const auto &[exps, coeffs_raw] : shell_exps_coeffs)
            {
                vector<double> coeffs = coeffs_raw;
                for (size_t i = 0; i < exps.size(); i++)
                    coeffs[i] *= calcPurePrimitiveNorm(angmom, exps[i]);

                vector<double> norms = calcShellNormalization(angmom, coeffs, exps);

                Shell shell(angmom, atomic_nr,
                            dim_cart, dim_sphe,
                            pos, xyz_coordinates,
                            coeffs, coeffs_raw,
                            exps, norms);

                shells[angmom].push_back(shell);

                pos += dim_sphe;
            }
        }
    }

    n_atomic_orbitals = pos;
    max_angular_momentum = max_angmom;
}

void LI::Structure::readBasis(const string &basis_set,
                              int &max_angular_momentum,
                              size_t &n_atomic_orbitals,
                              map<int, vector<Shell>> &shells)
{
    string basis_path = returnBasisPath(basis_set);

    constructShells(max_angular_momentum, n_atomic_orbitals, shells);
}