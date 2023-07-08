#include "defs_geomopt.h"
#include "geometry.h"

using namespace Lible;

GeomOpt::Geometry::Geometry(const std::vector<double> &coords, const std::vector<std::string> &atoms) : atoms(atoms)
{
    assert((cartesian_coords.size() % 3 == 0));
    assert((coords.size() % 3 == atoms.size()));

    n_atoms = cartesian_coords.size() / 3;
    cartesian_coords.resize(n_atoms);
    for (size_t i = 0; i < n_atoms; i++)
        cartesian_coords[i] = arma::dvec({coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]});

    constructRedIntCoords();
}

void GeomOpt::Geometry::constructRedIntCoords()
{
    std::vector<double> covalent_radii(n_atoms);
    for (size_t i = 0; i < n_atoms; i++)    
        covalent_radii[i] = GeomOptDefs::covalent_radii.at(atoms[i]);

    std::vector<std::vector<std::size_t>> atom_bonding_partners;
}