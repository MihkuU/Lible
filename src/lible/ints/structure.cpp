#include <lible/ints/structure.hpp>
#include <lible/ints/basis_sets.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/utils.hpp>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

namespace LI = lible::ints;

using std::array, std::string, std::vector;

LI::Structure::Structure()
{
}

// TODO: do error handling where atomic_nrs and coords_angstrom match

LI::Structure::Structure(const string &basis_set, const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : basis_set(basis_set), coords(coords_angstrom), atomic_nrs(atomic_nrs)
{
    n_atoms = atomic_nrs.size();

    elements.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
        elements[iatom] = atomic_symbols.at(atomic_nrs[iatom]);

    for (size_t i = 0; i < coords.size(); i++)
        coords[i] *= _ang_to_bohr_;

    coords_xyz.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
        coords_xyz[iatom] = {coords[3 * iatom], coords[3 * iatom + 1], coords[3 * iatom + 2]};

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms = basisForAtoms(atomic_nrs_set, basis_set);

    constructShells(basis_atoms, max_l, dim_ao, dim_ao_cart, shells);
}

LI::Structure::Structure(const string &basis_set, const vector<string> &elements,
                         const vector<double> &coords_angstrom)
    : basis_set(basis_set), coords(coords_angstrom), elements(elements)
{
    n_atoms = elements.size();

    atomic_nrs.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
    {
        string atomic_symbol = elements[iatom];
        std::transform(atomic_symbol.begin(), atomic_symbol.end(), atomic_symbol.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });

        atomic_symbol[0] = std::toupper(atomic_symbol[0]);
        this->elements[iatom] = atomic_symbol;

        atomic_nrs[iatom] = atomic_numbers.at(atomic_symbol);
    }

    for (size_t i = 0; i < coords.size(); i++)
        coords[i] *= _ang_to_bohr_;

    coords_xyz.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
        coords_xyz[iatom] = {coords[3 * iatom], coords[3 * iatom + 1], coords[3 * iatom + 2]};

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms = basisForAtoms(atomic_nrs_set, basis_set);

    constructShells(basis_atoms, max_l, dim_ao, dim_ao_cart, shells);
}

LI::Structure::Structure(const string &basis_set, const string &basis_set_aux,
                         const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = basis_set_aux;
    use_ri = true;

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_set, basis_set_aux);

    constructShells(basis_atoms_aux, max_l_aux, dim_ao_aux, dim_ao_cart_aux, shells_aux);
}

LI::Structure::Structure(const string &basis_set, const string &basis_set_aux,
                         const vector<string> &elements,
                         const vector<double> &coords_angstrom)
    : Structure(basis_set, elements, coords_angstrom)
{
    this->basis_set_aux = basis_set_aux;
    use_ri = true;

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_set, basis_set_aux);

    constructShells(basis_atoms_aux, max_l_aux, dim_ao_aux, dim_ao_cart_aux, shells_aux);
}

LI::Structure::Structure(const basis_atoms_t &basis_set_custom, const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : basis_set("custom"), coords(coords_angstrom), atomic_nrs(atomic_nrs)
{
    n_atoms = atomic_nrs.size();

    elements.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
        elements[iatom] = atomic_symbols.at(atomic_nrs[iatom]);

    for (size_t i = 0; i < coords.size(); i++)
        coords[i] *= _ang_to_bohr_;

    coords_xyz.resize(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
        coords_xyz[iatom] = {coords[3 * iatom], coords[3 * iatom + 1], coords[3 * iatom + 2]};

    constructShells(basis_set_custom, max_l, dim_ao, dim_ao_cart, shells);
}

LI::Structure::Structure(const string &basis_set, const basis_atoms_t &basis_set_custom_aux,
                         const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = "custom";
    use_ri = true;

    constructShells(basis_set_custom_aux, max_l_aux, dim_ao_aux, dim_ao_cart_aux, shells_aux);
}

LI::Structure::Structure(const basis_atoms_t &basis_set_custom, const string &basis_set_aux,
                         const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : Structure(basis_set_custom, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = basis_set_aux;
    use_ri = true;

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_set, basis_set_aux);

    constructShells(basis_atoms_aux, max_l_aux, dim_ao_aux, dim_ao_cart_aux, shells_aux);
}

LI::Structure::Structure(const basis_atoms_t &basis_set_custom,
                         const basis_atoms_t &basis_set_custom_aux,
                         const vector<int> &atomic_nrs,
                         const vector<double> &coords_angstrom)
    : Structure(basis_set_custom, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = "custom";
    use_ri = true;

    constructShells(basis_set_custom_aux, max_l_aux, dim_ao_aux, dim_ao_cart_aux, shells_aux);
}

bool LI::Structure::getUseRI() const
{
    return use_ri;
}

int LI::Structure::getMaxL() const
{
    return max_l;
}

int LI::Structure::getMaxLAux() const
{
    return max_l_aux;
}

int LI::Structure::getZ(const int iatom) const
{
    return atomic_nrs[iatom];
}

int LI::Structure::getDimAO() const
{
    return dim_ao;
}

int LI::Structure::getDimAOAux() const
{
    return dim_ao_aux;
}

int LI::Structure::getDimAOCart() const
{
    return dim_ao_cart;
}

int LI::Structure::getDimAOCartAux() const
{
    return dim_ao_cart_aux;
}

int LI::Structure::getNAtoms() const
{
    return n_atoms;
}

array<double, 3> LI::Structure::getCoordsAtom(const int iatom) const
{
    return coords_xyz[iatom];
}

vector<array<double, 4>> LI::Structure::getZs() const
{
    int n_atoms = getNAtoms();

    vector<array<double, 4>> atomic_charges(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
    {
        auto [x, y, z] = coords_xyz[iatom];
        double Z = atomic_nrs[iatom];
        atomic_charges[iatom] = {x, y, z, Z};
    }

    return atomic_charges;
}

vector<LI::Shell> LI::Structure::getShellsL(const int l) const
{
    return shells.at(l);
}

vector<LI::Shell> LI::Structure::getShellsLAux(const int l) const
{
    return shells_aux.at(l);
}

std::map<int, vector<LI::Shell>> LI::Structure::getShells() const
{
    return shells;
}

std::map<int, vector<LI::Shell>> LI::Structure::getShellsAux() const
{
    return shells_aux;
}

void LI::Structure::constructShells(const basis_atoms_t &basis_atoms, int &max_l, int &dim_ao,
                                    int &dim_ao_cart, std::map<int, vector<Shell>> &shells)
{
    int max_angmom = 0, pos = 0, pos_cart = 0;
    for (int iatom = 0; iatom < n_atoms; iatom++)
    {
        int atomic_nr = atomic_nrs[iatom];

        auto basis_atom = basis_atoms.at(atomic_nr);

        array<double, 3> coords_xyz{coords[3 * iatom], coords[3 * iatom + 1],
                                    coords[3 * iatom + 2]};
        for (const auto &[angmom, shell_exps_coeffs] : basis_atom)
        {
            if (angmom > max_angmom)
                max_angmom = angmom;

            if (max_angmom > _max_angular_momentum_)
            {
                std::string msg = std::format("Maximum angular momentum {} is larger than allowed {}!",
                                              max_angmom, _max_angular_momentum_);
                throw std::runtime_error(msg);
            }

            int dim_cart = numCartesians(angmom);
            int dim_sphe = numSphericals(angmom);

            for (const auto &[exps, coeffs_raw] : shell_exps_coeffs)
            {
                vector<double> coeffs = coeffs_raw;
                for (size_t i = 0; i < exps.size(); i++)
                    coeffs[i] *= calcPurePrimitiveNorm(angmom, exps[i]);

                vector<double> norms = calcShellNorms(angmom, coeffs, exps);

                Shell shell(angmom, atomic_nr, iatom, dim_cart, dim_sphe, pos, pos_cart,
                            coords_xyz, coeffs, coeffs_raw, exps, norms);

                shells[angmom].push_back(shell);

                pos += dim_sphe;
                pos_cart += dim_cart;
            }
        }
    }

    dim_ao = pos;
    dim_ao_cart = pos_cart;
    max_l = max_angmom;
}