#include <lible/geomopt/geomopt.hpp>

#include <armadillo>

namespace lgopt = lible::geomopt;

namespace lible::geomopt
{
    /// Type alias for a nested (2-dimensional) vector.
    template<typename T>
    using vec2d = std::vector<std::vector<T>>;

    /// Type alias for a list of atomic coordinates.
    using xyz_coords_arma_t = std::vector<arma::vec3>;

    /// Calculates the bond length for two atoms.
    double bondLength(const arma::vec3 &xyz_m, const arma::vec3 &xyz_n);

    /// Calculates a bond angle between three atoms.
    double bondAngle(const arma::vec3 &xyz_m, const arma::vec3 &xyz_o, const arma::vec3 &xyz_n);

    /// Calculates a dihedral angle formed by four atoms.
    double dihedralAngle(const arma::vec3 &xyz_m, const arma::vec3 &xyz_o, const arma::vec3 &xyz_p,
                         const arma::vec3 &xyz_n);

    /// Finds the bonding partner atoms for each atom.
    vec2d<size_t> bondingPartners(const std::vector<int> &atomic_nrs,
                                  const xyz_coords_arma_t &xyz_coords);

    /// Constructs and returns a list of bond length coordinates.
    std::vector<BondLength> bondLengths(const vec2d<size_t> &bonding_partners,
                                        const xyz_coords_arma_t &xyz_coords);

    /// Constructs and returns a list of bond angle coordinates.
    std::vector<BondAngle> bondAngles(const vec2d<size_t> &bonding_partners,
                                      const xyz_coords_arma_t &xyz_coords);

    /// Constructs and returns a list of dihedral angle coordinates.
    std::vector<DihedralAngle> dihedralAngles(const vec2d<size_t> &bonding_partners,
                                              const xyz_coords_arma_t &xyz_coords);

    /// Converts the coordinate list to arma.
    xyz_coords_arma_t conv2arma(const xyz_coords_t &xyz_coords);

    /// Constructs the Wilson's B matrix using formulas from https://doi.org/10.1063/1.1515483.
    arma::dmat builB(const xyz_coords_arma_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Constructs the K matrix defined in eq. (8) from https://doi.org/10.1063/1.1515483.
    arma::dmat buildK(const arma::dmat &b_matrix, const arma::dvec &grad_redint,
                      const xyz_coords_arma_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Covalent radii of all elements up to an atomic number of 96 in Bohr. Note, for Mn, Fe and Co
    /// the values for high spin are given, and for C the value for sp3 is given. All data are
    /// extracted from https://doi.org/10.1039/B801115J.
    const inline std::map<int, double> covalent_radii_table{
        {1, 0.59},
        {2, 0.53},
        {3, 2.42},
        {4, 1.81},
        {5, 1.59},
        {6, 1.44},
        {7, 1.34},
        {8, 1.25},
        {9, 1.08},
        {10, 1.10},
        {11, 3.14},
        {12, 2.66},
        {13, 2.29},
        {14, 2.10},
        {15, 2.02},
        {16, 1.98},
        {17, 1.93},
        {18, 2.00},
        {19, 3.84},
        {20, 3.33},
        {21, 3.21},
        {22, 3.02},
        {23, 2.89},
        {24, 2.63},
        {25, 3.04},
        {26, 2.87},
        {27, 2.83},
        {28, 2.34},
        {29, 2.49},
        {30, 2.31},
        {31, 2.31},
        {32, 2.27},
        {33, 2.25},
        {34, 2.27},
        {35, 2.27},
        {36, 2.19},
        {37, 4.16},
        {38, 3.68},
        {39, 3.59},
        {40, 3.31},
        {41, 3.10},
        {42, 2.91},
        {43, 2.78},
        {44, 2.76},
        {45, 2.68},
        {46, 2.63},
        {47, 2.74},
        {48, 2.72},
        {49, 2.68},
        {50, 2.63},
        {51, 2.63},
        {52, 2.61},
        {53, 2.63},
        {54, 2.65},
        {55, 4.61},
        {56, 4.06},
        {57, 3.91},
        {58, 3.86},
        {59, 3.84},
        {60, 3.80},
        {61, 3.76},
        {62, 3.74},
        {63, 3.74},
        {64, 3.70},
        {65, 3.67},
        {66, 3.63},
        {67, 3.63},
        {68, 3.57},
        {69, 3.59},
        {70, 3.53},
        {71, 3.53},
        {72, 3.31},
        {73, 3.21},
        {74, 3.06},
        {75, 2.85},
        {76, 2.72},
        {77, 2.66},
        {78, 2.57},
        {79, 2.57},
        {80, 2.49},
        {81, 2.74},
        {82, 2.76},
        {83, 2.80},
        {84, 2.65},
        {85, 2.83},
        {86, 2.83},
        {87, 4.91},
        {88, 4.18},
        {89, 4.06},
        {90, 3.89},
        {91, 3.78},
        {92, 3.70},
        {93, 3.59},
        {94, 3.53},
        {95, 3.40},
        {96, 3.19}
    };

    /// Taken from https://doi.org/10.1063/1.1515483.
    const static double bonding_factor = 1.3;
}

lgopt::RedIntCoords lgopt::redIntCoords(const std::vector<int> &atomic_nrs,
                                        const xyz_coords_t &xyz_coords)
{
    xyz_coords_arma_t xyz_coords_arma = conv2arma(xyz_coords);

    vec2d<size_t> bonding_partners = bondingPartners(atomic_nrs, xyz_coords_arma);

    return {
        bondLengths(bonding_partners, xyz_coords_arma),
        bondAngles(bonding_partners, xyz_coords_arma),
        dihedralAngles(bonding_partners, xyz_coords_arma)
    };
}

std::vector<lgopt::BondLength> lgopt::bondLengths(const vec2d<size_t> &bonding_partners,
                                                  const xyz_coords_arma_t &xyz_coords)
{
    std::vector<BondLength> bond_lengths;
    for (size_t m = 0; m < bonding_partners.size(); m++)
        for (size_t n : bonding_partners.at(m))
            if (m < n)
            {
                double bond_length = bondLength(xyz_coords[m], xyz_coords[n]);
                bond_lengths.push_back({m, n, bond_length});
            }

    return bond_lengths;
}

std::vector<lgopt::BondAngle> lgopt::bondAngles(const vec2d<size_t> &bonding_partners,
                                                const xyz_coords_arma_t &xyz_coords)
{
    std::vector<BondAngle> bond_angles;
    for (size_t o = 0; o < bonding_partners.size(); o++)
        for (size_t m : bonding_partners.at(o))
            for (size_t n : bonding_partners.at(o))
                if (m < n)
                {
                    double bond_angle = bondAngle(xyz_coords[m], xyz_coords[o], xyz_coords[n]);
                    bond_angles.push_back({o, m, n, bond_angle});
                }

    return bond_angles;
}

std::vector<lgopt::DihedralAngle> lgopt::dihedralAngles(const vec2d<size_t> &bonding_partners,
                                                        const xyz_coords_arma_t &xyz_coords)
{
    std::vector<std::pair<size_t, size_t>> op_bonds;
    for (size_t o = 0; o < bonding_partners.size(); o++)
    {
        const std::vector<size_t> &bonding_partners_o = bonding_partners.at(o);
        if (bonding_partners_o.size() < 2)
            continue;

        for (size_t p : bonding_partners_o)
        {
            if (bonding_partners.at(p).size() < 2)
                continue;

            if (o < p)
                op_bonds.emplace_back(o, p);
        }
    }

    std::vector<DihedralAngle> dihedral_angles;
    for (const auto &[o, p] : op_bonds)
    {
        const std::vector<size_t> &bonding_partners_o = bonding_partners.at(o);
        const std::vector<size_t> &bonding_partners_p = bonding_partners.at(p);

        for (size_t m : bonding_partners_o)
            for (size_t n : bonding_partners_p)
                if (m != p && n != o)
                {
                    double dihedral_angle = dihedralAngle(xyz_coords[m], xyz_coords[o],
                                                          xyz_coords[p], xyz_coords[n]);
                    dihedral_angles.push_back({o, m, p, n, dihedral_angle});
                }
    }

    return dihedral_angles;
}

size_t lgopt::numRedIntCoords(const RedIntCoords &red_int_coords)
{
    const auto &[bonds, angles, dihedrals] = red_int_coords;
    return bonds.size() + angles.size() + dihedrals.size();
}

lgopt::vec2d<size_t> lgopt::bondingPartners(const std::vector<int> &atomic_nrs,
                                            const xyz_coords_arma_t &xyz_coords)
{
    size_t n_atoms = atomic_nrs.size();

    vec2d<size_t> bonding_partners(n_atoms);
    for (size_t m = 0; m < n_atoms; m++)
    {
        std::vector<size_t> bonding_partners_m;
        for (size_t n = 0; n < n_atoms; n++)
        {
            if (m == n)
                continue;

            int z_m = atomic_nrs[m];
            int z_n = atomic_nrs[n];

            double r_m = covalent_radii_table.at(z_m);
            double r_n = covalent_radii_table.at(z_n);
            double sum_r_mn = r_m + r_n;

            double distance = bondLength(xyz_coords[m], xyz_coords[n]);
            double bonding_distance = bonding_factor * sum_r_mn;

            if (distance < bonding_distance)
                bonding_partners_m.push_back(n);
        }
        bonding_partners[m] = bonding_partners_m;
        // TODO: if we do not have a bonding partner use the least distant atom
    }

    return bonding_partners;
}

double lgopt::bondLength(const arma::vec3 &xyz_m, const arma::vec3 &xyz_n)
{
    return arma::norm(xyz_m - xyz_n);
}

double lgopt::bondAngle(const arma::vec3 &xyz_m, const arma::vec3 &xyz_o, const arma::vec3 &xyz_n)
{
    arma::dvec3 u = xyz_m - xyz_o;
    arma::dvec3 v = xyz_n - xyz_o;

    return std::acos(arma::dot(u, v) / (arma::norm(u) * arma::norm(v)));
}

double lgopt::dihedralAngle(const arma::dvec3 &xyz_m, const arma::dvec3 &xyz_o,
                            const arma::dvec3 &xyz_p, const arma::dvec3 &xyz_n)
{
    arma::dvec3 u = xyz_m - xyz_o;
    arma::dvec3 v = xyz_n - xyz_p;
    arma::dvec3 w = xyz_p - xyz_o;

    u /= arma::norm(u);
    v /= arma::norm(v);
    w /= arma::norm(w);

    double sin_u = std::sqrt(1 - std::pow(arma::dot(u, w), 2));
    double sin_v = std::sqrt(1 - std::pow(arma::dot(v, w), 2));

    return std::acos(arma::dot(arma::cross(u, w), arma::cross(v, w)) / (sin_u * sin_v));
}

lgopt::xyz_coords_arma_t lgopt::conv2arma(const xyz_coords_t &xyz_coords)
{
    xyz_coords_arma_t xyz_coords_arma(xyz_coords.size());
    for (size_t i = 0; i < xyz_coords.size(); i++)
        xyz_coords_arma[i] = {xyz_coords[i][0], xyz_coords[i][1], xyz_coords[i][2]};

    return xyz_coords_arma;
}
