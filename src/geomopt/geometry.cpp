#include <cmath>
#include <limits>
#include <numbers>

#include "defs_geomopt.h"
#include "geometry.h"

using namespace Lible;
using std::pair;
using std::set;
using std::size_t;
using std::vector;

GeomOpt::Geometry::Geometry(const vector<double> &coords_cart_, const vector<std::string> &atoms_) : coords_cart(coords_cart_), atoms(atoms_)
{
    assert((coords_cart.size() % 3 == 0));
    assert((coords_cart.size() / 3 == atoms.size()));

    n_atoms = atoms.size() / 3;
    atom_coords_cart.resize(n_atoms);
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
        atom_coords_cart[iatom] = arma::dvec({coords_cart[3 * iatom], coords_cart[3 * iatom + 1], coords_cart[3 * iatom + 2]});

    red_int_coords = constructRedIntCoords();
}

double GeomOpt::Geometry::calcDistance(const size_t &iatom, const size_t &jatom)
{
    return arma::norm(atom_coords_cart[iatom] - atom_coords_cart[jatom]);
}

double GeomOpt::Geometry::calcAngle(const size_t &iatom, const size_t &jatom, const size_t &katom)
{
    arma::dvec ij_vec = atom_coords_cart[iatom] - atom_coords_cart[jatom];
    arma::dvec kj_vec = atom_coords_cart[katom] - atom_coords_cart[jatom];
    return acos(arma::dot(ij_vec, kj_vec) / (arma::norm(ij_vec) * arma::norm(kj_vec)));
}

double GeomOpt::Geometry::calcDihedral(const size_t &iatom, const size_t &jatom, const size_t &katom, const size_t &latom)
{
    return 0; // TODO
}

size_t GeomOpt::Geometry::findClosestAtom(const size_t &iatom)
{
    double min_distance = std::numeric_limits<double>::max();
    size_t closest_atom;
    for (size_t jatom = iatom + 1; jatom < n_atoms; jatom++)
    {
        double distance = calcDistance(iatom, jatom);
        if (distance < min_distance)
        {
            min_distance = distance;
            closest_atom = jatom;
        }
    }
    return closest_atom;
}

vector<set<size_t>> GeomOpt::Geometry::findAtomBondingPartners()
{
    vector<double> covalent_radii(n_atoms);
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
        covalent_radii[iatom] = GeomOptDefs::covalent_radii.at(atoms[iatom]);

    vector<set<size_t>> atom_bonding_partners(n_atoms);
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        set<size_t> bonding_partners;
        for (size_t jatom = iatom + 1; jatom < n_atoms; jatom++)
        {
            double sum_cov_radii = covalent_radii[iatom] + covalent_radii[jatom];
            double bonding_distance = GeomOptDefs::bonding_factor * sum_cov_radii;
            double distance = calcDistance(iatom, jatom);
            if (distance < bonding_distance)
                bonding_partners.insert(jatom);
        }

        if (bonding_partners.size() == 0)
            bonding_partners.insert(findClosestAtom(iatom));

        atom_bonding_partners[iatom] = bonding_partners;
    }
    return atom_bonding_partners;
}

GeomOpt::Geometry::RedundantInternalCoordinates GeomOpt::Geometry::constructRedIntCoords()
{
    vector<set<size_t>> atom_bonding_partners = findAtomBondingPartners();

    vector<pair<doublet, double>> bonds;
    vector<pair<triplet, double>> angles;
    vector<pair<quartet, double>> dihedrals;
    for (size_t iatom = 0; iatom < n_atoms; iatom++)
    {
        set<size_t> bonding_partners_iatom = atom_bonding_partners[iatom];
        for (const size_t &jatom : bonding_partners_iatom)
        {
            bonds.emplace_back(std::make_pair(std::make_tuple(iatom, jatom), calcDistance(iatom, jatom)));

            set<size_t> bonding_partners_jatom = atom_bonding_partners[jatom];
            for (const size_t &katom : bonding_partners_jatom)
            {
                angles.emplace_back(std::make_pair(std::make_tuple(iatom, jatom, katom),
                                                   calcAngle(iatom, jatom, katom)));

                set<size_t> bonding_partners_katom = atom_bonding_partners[katom];
                for (const size_t &latom : bonding_partners_katom)
                    dihedrals.emplace_back(std::make_pair(std::make_tuple(iatom, jatom, katom, latom),
                                                          calcDihedral(iatom, jatom, katom, latom)));
            }
        }
    }

    RedundantInternalCoordinates red_int_coords;
    red_int_coords.bonds = bonds;
    red_int_coords.angles = angles;
    red_int_coords.dihedrals = dihedrals;

    return red_int_coords;
}

void GeomOpt::Geometry::constructBMatrix(arma::dmat &bmatrix)
{

}

void GeomOpt::Geometry::constructBMatrix(arma::sp_dmat &bmatrix)
{

}