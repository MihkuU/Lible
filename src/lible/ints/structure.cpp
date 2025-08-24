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

namespace lints = lible::ints;

using std::array, std::string, std::vector;

lints::Structure::Structure()
{
}

lints::Structure::Structure(const string &basis_set, const vector<int> &atomic_nrs,
                            const vector<array<double, 3>> &coords_angstrom)
    : basis_set(basis_set), coords(coords_angstrom), atomic_nrs(atomic_nrs)
{
    if (atomic_nrs.size() != coords.size())
        throw std::runtime_error("Structure::Structure(): number of atomic numbers does not match "
                                 "the number of (x, y, z)-oordinates");

    n_atoms = atomic_nrs.size();

    for (int iatom = 0; iatom < n_atoms; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords[iatom][icart] *= _ang_to_bohr_;

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms = basisForAtoms(atomic_nrs_set, basis_set);

    shells = constructShells(basis_atoms, atomic_nrs, coords);

    for (const Shell &shell : shells)
    {
        shells_map[shell.l].push_back(shell);

        if (shell.l > max_l)
            max_l = shell.l;

        dim_ao += shell.dim_sph;
        dim_ao_cart += shell.dim_cart;
    }
}

lints::Structure::Structure(const string &basis_set, const string &basis_set_aux,
                            const vector<int> &atomic_nrs,
                            const vector<array<double, 3>> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = basis_set_aux;
    use_ri = true;

    std::set<int> atomic_nrs_set(atomic_nrs.begin(), atomic_nrs.end());
    basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_set, basis_set_aux);

    shells_aux = constructShells(basis_atoms_aux, atomic_nrs, coords);

    for (const Shell &shell : shells_aux)
    {
        shells_map_aux[shell.l].push_back(shell);

        if (shell.l > max_l_aux)
            max_l_aux = shell.l;

        dim_ao_aux += shell.dim_sph;
        dim_ao_cart_aux += shell.dim_cart;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set, const vector<int> &atomic_nrs,
                            const vector<array<double, 3>> &coords_angstrom)
    : coords(coords_angstrom), atomic_nrs(atomic_nrs)
{
    if (atomic_nrs.size() != coords.size())
        throw std::runtime_error("Structure::Structure(): number of atomic numbers does not match "
                                 "the number of (x, y, z)-oordinates");

    this->basis_set = "custom";
    n_atoms = atomic_nrs.size();

    for (int iatom = 0; iatom < n_atoms; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords[iatom][icart] *= _ang_to_bohr_;

    shells = constructShells(basis_set, atomic_nrs, coords);

    for (const Shell &shell : shells)
    {
        shells_map[shell.l].push_back(shell);

        if (shell.l > max_l)
            max_l = shell.l;

        dim_ao += shell.dim_sph;
        dim_ao_cart += shell.dim_cart;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set,
                            const basis_atoms_t &basis_set_aux,
                            const vector<int> &atomic_nrs,
                            const vector<array<double, 3>> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    this->basis_set_aux = "custom";
    use_ri = true;

    shells_aux = constructShells(basis_set_aux, atomic_nrs, coords);

    for (const Shell &shell : shells_aux)
    {
        shells_map_aux[shell.l].push_back(shell);

        if (shell.l > max_l_aux)
            max_l_aux = shell.l;

        dim_ao_aux += shell.dim_sph;
        dim_ao_cart_aux += shell.dim_cart;
    }
}

bool lints::Structure::getUseRI() const
{
    return use_ri;
}

int lints::Structure::getMaxL() const
{
    return max_l;
}

int lints::Structure::getMaxLAux() const
{
    return max_l_aux;
}

int lints::Structure::getZ(const int iatom) const
{
    return atomic_nrs[iatom];
}

int lints::Structure::getDimAO() const
{
    return dim_ao;
}

int lints::Structure::getDimAOAux() const
{
    return dim_ao_aux;
}

int lints::Structure::getDimAOCart() const
{
    return dim_ao_cart;
}

int lints::Structure::getDimAOCartAux() const
{
    return dim_ao_cart_aux;
}

int lints::Structure::getNAtoms() const
{
    return n_atoms;
}

array<double, 3> lints::Structure::getCoordsAtom(const int iatom) const
{
    return coords[iatom];
}

vector<array<double, 4>> lints::Structure::getZs() const
{
    int n_atoms = getNAtoms();

    vector<array<double, 4>> atomic_charges(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
    {
        auto [x, y, z] = coords[iatom];
        double Z = atomic_nrs[iatom];
        atomic_charges[iatom] = {x, y, z, Z};
    }

    return atomic_charges;
}

vector<lints::Shell> lints::Structure::getShells() const
{
    return shells;
}

vector<lints::Shell> lints::Structure::getShellsAux() const
{
    return shells_aux;
}

vector<lints::Shell> lints::Structure::getShellsL(const int l) const
{
    return shells_map.at(l);
}

vector<lints::Shell> lints::Structure::getShellsLAux(const int l) const
{
    return shells_map_aux.at(l);
}

std::map<int, vector<lints::Shell>> lints::Structure::getShellsMap() const
{
    return shells_map;
}

std::map<int, vector<lints::Shell>> lints::Structure::getShellsMapAux() const
{
    return shells_map_aux;
}