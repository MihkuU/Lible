#include <lible/ints/structure.hpp>
#include <lible/ints/basis_sets.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/utils.hpp>

#include <filesystem>
#include <map>

namespace lints = lible::ints;

lints::Structure::Structure(const std::string &basis_set, const std::vector<int> &atomic_nrs,
                            const std::vector<std::array<double, 3>> &coords_angstrom)
    : coords_(coords_angstrom), atomic_nrs_(atomic_nrs), basis_set_(basis_set)
{
    if (atomic_nrs.size() != coords_.size())
        throw std::runtime_error("Structure::Structure(): number of atomic numbers does not match "
            "the number of (x, y, z)-coordinates");

    n_atoms_ = atomic_nrs.size();

    for (size_t iatom = 0; iatom < n_atoms_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_[iatom][icart] *= _ang_to_bohr_;

    basis_atoms_t basis_atoms = basisForAtoms(atomic_nrs_, basis_set);

    shells_ = constructShells(basis_atoms, coords_);

    for (const Shell &shell : shells_)
    {
        shells_map_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_)
            max_l_ = shell.l_;

        dim_ao_ += shell.dim_sph_;
        dim_ao_cart_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const std::string &basis_set, const std::string &basis_set_aux,
                            const std::vector<int> &atomic_nrs,
                            const std::vector<std::array<double, 3>> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    basis_set_aux_ = basis_set_aux;
    use_ri_ = true;

    const basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_, basis_set_aux);

    shells_aux_ = constructShells(basis_atoms_aux, coords_);

    for (const Shell &shell : shells_aux_)
    {
        shells_map_aux_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_aux_)
            max_l_aux_ = shell.l_;

        dim_ao_aux_ += shell.dim_sph_;
        dim_ao_cart_aux_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set, const std::vector<int> &atomic_nrs,
                            const std::vector<std::array<double, 3>> &coords_angstrom)
    : coords_(coords_angstrom), atomic_nrs_(atomic_nrs)
{
    if (atomic_nrs.size() != coords_.size())
        throw std::runtime_error("Structure::Structure(): number of atomic numbers does not match "
            "the number of (x, y, z)-coordinates");

    basis_set_ = "custom";
    n_atoms_ = atomic_nrs.size();

    for (size_t iatom = 0; iatom < n_atoms_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_[iatom][icart] *= _ang_to_bohr_;

    shells_ = constructShells(basis_set, coords_);

    for (const Shell &shell : shells_)
    {
        shells_map_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_)
            max_l_ = shell.l_;

        dim_ao_ += shell.dim_sph_;
        dim_ao_cart_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set,
                            const basis_atoms_t &basis_set_aux,
                            const std::vector<int> &atomic_nrs,
                            const std::vector<std::array<double, 3>> &coords_angstrom)
    : Structure(basis_set, atomic_nrs, coords_angstrom)
{
    basis_set_aux_ = "custom";
    use_ri_ = true;

    shells_aux_ = constructShells(basis_set_aux, coords_);

    for (const Shell &shell : shells_aux_)
    {
        shells_map_aux_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_aux_)
            max_l_aux_ = shell.l_;

        dim_ao_aux_ += shell.dim_sph_;
        dim_ao_cart_aux_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const std::string &basis_set, const std::string &basis_set_ghost,
                            const std::vector<int> &atomic_nrs, const std::vector<int> &atomic_nrs_ghost,
                            const std::vector<std::array<double, 3>> &coords_angstrom,
                            const std::vector<std::array<double, 3>> &coords_angstrom_ghost)
    : coords_(coords_angstrom), coords_ghost_(coords_angstrom_ghost), atomic_nrs_(atomic_nrs),
      atomic_nrs_ghost_(atomic_nrs_ghost), basis_set_(basis_set), basis_set_ghost_(basis_set_ghost)
{
    if (atomic_nrs_ghost_.size() != coords_ghost_.size() or atomic_nrs_.size() != coords_.size())
        throw std::runtime_error("Structure::Structure(): number of ghost atomic or atomic numbers does not match "
            "the number of corresponding (x, y, z)-coordinates");

    n_atoms_ = atomic_nrs_.size();
    n_atoms_ghost_ = atomic_nrs_ghost_.size();

    for (size_t iatom = 0; iatom < n_atoms_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_[iatom][icart] *= _ang_to_bohr_;
    for (size_t iatom = 0; iatom < n_atoms_ghost_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_ghost_[iatom][icart] *= _ang_to_bohr_;

    const basis_atoms_t basis_atoms = basisForAtoms(atomic_nrs_, basis_set_);
    const basis_atoms_t basis_atoms_ghost = basisForAtoms(atomic_nrs_ghost_, basis_set_ghost_);

    shells_ = constructShellsGhost(basis_atoms, basis_atoms_ghost, coords_, coords_ghost_);

    for (const Shell &shell : shells_)
    {
        shells_map_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_)
            max_l_ = shell.l_;

        dim_ao_ += shell.dim_sph_;
        dim_ao_cart_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const std::string &basis_set, const std::string &basis_set_ghost,
                            const std::string &basis_set_aux, const std::string &basis_set_aux_ghost,
                            const std::vector<int> &atomic_nrs, const std::vector<int> &atomic_nrs_ghost,
                            const std::vector<std::array<double, 3>> &coords_angstrom,
                            const std::vector<std::array<double, 3>> &coords_angstrom_ghost)
    : Structure(basis_set, basis_set_ghost, atomic_nrs, atomic_nrs_ghost, coords_angstrom, coords_angstrom_ghost)
{
    basis_set_aux_ = basis_set_aux;
    basis_set_aux_ghost_ = basis_set_aux_ghost;
    use_ri_ = true;

    const basis_atoms_t basis_atoms_aux = basisForAtomsAux(atomic_nrs_, basis_set_aux_);
    const basis_atoms_t basis_atoms_aux_ghost = basisForAtomsAux(atomic_nrs_ghost_, basis_set_aux_ghost_);

    shells_aux_ = constructShellsGhost(basis_atoms_aux, basis_atoms_aux_ghost, coords_, coords_ghost_);

    for (const Shell &shell : shells_aux_)
    {
        shells_map_aux_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_aux_)
            max_l_aux_ = shell.l_;

        dim_ao_aux_ += shell.dim_sph_;
        dim_ao_cart_aux_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set, const basis_atoms_t &basis_set_ghost,
                            const std::vector<int> &atomic_nrs, const std::vector<int> &atomic_nrs_ghost,
                            const std::vector<std::array<double, 3>> &coords_angstrom,
                            const std::vector<std::array<double, 3>> &coords_angstrom_ghost)
    : coords_(coords_angstrom), coords_ghost_(coords_angstrom_ghost), atomic_nrs_(atomic_nrs), atomic_nrs_ghost_(atomic_nrs_ghost)
{
    if (atomic_nrs_.size() != coords_.size() or atomic_nrs_ghost_.size() != coords_ghost_.size())
        throw std::runtime_error("Structure::Structure(): number of atomic numbers does not match "
            "the number of (x, y, z)-coordinates");

    basis_set_ = "custom";
    n_atoms_ = atomic_nrs_.size();
    n_atoms_ghost_ = atomic_nrs_ghost_.size();

    for (size_t iatom = 0; iatom < n_atoms_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_[iatom][icart] *= _ang_to_bohr_;
    for (size_t iatom = 0; iatom < n_atoms_ghost_; iatom++)
        for (int icart = 0; icart < 3; icart++)
            coords_ghost_[iatom][icart] *= _ang_to_bohr_;

    shells_ = constructShellsGhost(basis_set, basis_set_ghost, coords_, coords_ghost_);

    for (const Shell &shell : shells_)
    {
        shells_map_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_)
            max_l_ = shell.l_;

        dim_ao_ += shell.dim_sph_;
        dim_ao_cart_ += shell.dim_cart_;
    }
}

lints::Structure::Structure(const basis_atoms_t &basis_set, const basis_atoms_t &basis_set_ghost,
                            const basis_atoms_t &basis_set_aux, const basis_atoms_t &basis_set_ghost_aux,
                            const std::vector<int> &atomic_nrs, const std::vector<int> &atomic_nrs_ghost,
                            const std::vector<std::array<double, 3>> &coords_angstrom,
                            const std::vector<std::array<double, 3>> &coords_angstrom_ghost)
    : Structure(basis_set, basis_set_ghost, atomic_nrs, atomic_nrs_ghost, coords_angstrom, coords_angstrom_ghost)
{
    basis_set_aux_ = "custom";
    use_ri_ = true;

    shells_aux_ = constructShellsGhost(basis_set_aux, basis_set_ghost_aux, coords_, coords_ghost_);

    for (const Shell &shell : shells_aux_)
    {
        shells_map_aux_[shell.l_].push_back(shell);

        if (shell.l_ > max_l_aux_)
            max_l_aux_ = shell.l_;

        dim_ao_aux_ += shell.dim_sph_;
        dim_ao_cart_aux_ += shell.dim_cart_;
    }
}

bool lints::Structure::getUseRI() const
{
    return use_ri_;
}

int lints::Structure::getMaxL() const
{
    return max_l_;
}

int lints::Structure::getMaxLAux() const
{
    return max_l_aux_;
}

int lints::Structure::getZ(const size_t iatom) const
{
    return atomic_nrs_[iatom];
}

size_t lints::Structure::getDimAO() const
{
    return dim_ao_;
}

size_t lints::Structure::getDimAOAux() const
{
    return dim_ao_aux_;
}

size_t lints::Structure::getDimAOCart() const
{
    return dim_ao_cart_;
}

size_t lints::Structure::getDimAOCartAux() const
{
    return dim_ao_cart_aux_;
}

size_t lints::Structure::getNAtoms() const
{
    return n_atoms_;
}

std::array<double, 3> lints::Structure::getCoordsAtom(const size_t iatom) const
{
    return coords_[iatom];
}

std::vector<std::array<double, 4>> lints::Structure::getZs() const
{
    size_t n_atoms = getNAtoms();

    std::vector<std::array<double, 4>> atomic_charges(n_atoms);
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        auto [x, y, z] = coords_[iatom];
        double Z = atomic_nrs_[iatom];
        atomic_charges[iatom] = {x, y, z, Z};
    }

    return atomic_charges;
}

std::vector<lints::Shell> lints::Structure::getShells() const
{
    return shells_;
}

std::vector<lints::Shell> lints::Structure::getShellsAux() const
{
    return shells_aux_;
}

std::vector<lints::Shell> lints::Structure::getShellsL(const int l) const
{
    return shells_map_.at(l);
}

std::vector<lints::Shell> lints::Structure::getShellsLAux(const int l) const
{
    return shells_map_aux_.at(l);
}

std::map<int, std::vector<lints::Shell>> lints::Structure::getShellsMap() const
{
    return shells_map_;
}

std::map<int, std::vector<lints::Shell>> lints::Structure::getShellsMapAux() const
{
    return shells_map_aux_;
}
