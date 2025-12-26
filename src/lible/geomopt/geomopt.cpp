#include <lible/geomopt/geomopt.hpp>
#include <lible/geomopt/utils.hpp>

#include <numbers>

#include <armadillo>

namespace lgopt = lible::geomopt;

namespace lible::geomopt
{
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
    vec2d<size_t> bonding_partners = bondingPartners(atomic_nrs, xyz_coords);

    return {
        bondLengths(bonding_partners, xyz_coords),
        bondAngles(bonding_partners, xyz_coords),
        dihedralAngles(bonding_partners, xyz_coords)
    };
}

std::vector<lgopt::BondLength> lgopt::bondLengths(const vec2d<size_t> &bonding_partners,
                                                  const xyz_coords_t &xyz_coords)
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
                                                const xyz_coords_t &xyz_coords)
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
                                                        const xyz_coords_t &xyz_coords)
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
                                            const xyz_coords_t &xyz_coords)
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

double lgopt::bondLength(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_n)
{
    return norm(xyz_m - xyz_n);
}

double lgopt::bondAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                        const std::array<double, 3> &xyz_n)
{
    std::array<double, 3> u = xyz_m - xyz_o;
    std::array<double, 3> v = xyz_n - xyz_o;

    return std::acos(dot(u, v) / (norm(u) * norm(v)));
}

double lgopt::dihedralAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                            const std::array<double, 3> &xyz_p, const std::array<double, 3> &xyz_n)
{
    std::array<double, 3> u = xyz_m - xyz_o;
    std::array<double, 3> v = xyz_n - xyz_p;
    std::array<double, 3> w = xyz_p - xyz_o;

    u = u / norm(u);
    v = v / norm(v);
    w = w / norm(w);

    double sin_u = std::sqrt(1 - std::pow(dot(u, w), 2));
    double sin_v = std::sqrt(1 - std::pow(dot(v, w), 2));

    return std::acos(dot(cross(u, w), cross(v, w)) / (sin_u * sin_v));
}

lgopt::vec2d<double> lgopt::builBMatrix(const xyz_coords_t &xyz_coords,
                                        const RedIntCoords &red_int_coords,
                                        const double tol)
{
    vec2d<double> b_matrix_bonds = buildBMatrixBondLengths(xyz_coords, red_int_coords);

    vec2d<double> b_matrix_angles = buildBMatrixBondAngles(xyz_coords, red_int_coords, tol);

    vec2d<double> b_matrix_dihedrals = buildBMatrixDihedralAngles(xyz_coords, red_int_coords, tol);

    vec2d<double> b_matrix = b_matrix_bonds;
    b_matrix.insert(b_matrix.end(), b_matrix_angles.begin(), b_matrix_angles.end());
    b_matrix.insert(b_matrix.end(), b_matrix_dihedrals.begin(), b_matrix_dihedrals.end());

    return b_matrix;
}

lgopt::vec2d<double> lgopt::buildBMatrixBondLengths(const xyz_coords_t &xyz_coords,
                                                    const RedIntCoords &red_int_coords)
{
    size_t n_bond_lenghts = red_int_coords.bond_lengths_.size();

    vec2d<double> b_matrix_bonds(n_bond_lenghts);
    for (size_t ibond = 0; ibond < n_bond_lenghts; ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        std::array<double, 3> bond_vector = (xyz_coords[n] - xyz_coords[m]) / bond_length;
        const auto &[x, y, z] = bond_vector;
        b_matrix_bonds[ibond] = {x, y, z, -x, -y, -z};
    }

    return b_matrix_bonds;
}

lgopt::vec2d<double> lgopt::buildBMatrixBondAngles(const xyz_coords_t &xyz_coords,
                                                   const RedIntCoords &red_int_coords,
                                                   const double tol)
{
    size_t n_bond_angles = red_int_coords.bond_angles_.size();

    vec2d<double> b_matrix_angles(n_bond_angles);
    for (size_t iangle = 0; iangle < n_bond_angles; iangle++)
    {
        const auto &[m, o, n, bond_angle] = red_int_coords.bond_angles_[iangle];

        // u, v
        std::array<double, 3> u = xyz_coords[m] - xyz_coords[o];
        std::array<double, 3> v = xyz_coords[n] - xyz_coords[o];
        double lambda_u = norm(u);
        double lambda_v = norm(v);
        u = u / lambda_u;
        v = v / lambda_v;

        // w
        bool u_v_parallel = areParallel(u, v, tol);
        std::array<double, 3> w{};
        if (u_v_parallel == false)
            w = cross(u, v);
        else
        {
            std::array<double, 3> tmp{1, -1, 1};
            if (areParallel(u, tmp, tol) == false)
                w = cross(u, tmp);
            else
            {
                tmp = {-1, 1, -1};
                w = cross(u, tmp);
            }
        }

        w = w / norm(w);

        // derivatives
        std::array<double, 3> u_x_w = cross(u, w);
        std::array<double, 3> w_x_v = cross(w, v);

        auto calcContrib = [&](const size_t a) -> std::array<double, 3>
        {
            std::array<double, 3> result{};
            result[0] = zeta(a, m, o) * u_x_w[0] / lambda_u + zeta(a, n, o) * w_x_v[0] / lambda_v;
            result[1] = zeta(a, m, o) * u_x_w[1] / lambda_u + zeta(a, n, o) * w_x_v[1] / lambda_v;
            result[2] = zeta(a, m, o) * u_x_w[2] / lambda_u + zeta(a, n, o) * w_x_v[2] / lambda_v;

            return result;
        };

        std::array<double, 3> d_m = calcContrib(m);
        std::array<double, 3> d_o = calcContrib(o);
        std::array<double, 3> d_n = calcContrib(n);

        b_matrix_angles[iangle] = {
            d_m[0], d_m[1], d_m[2], d_o[0], d_o[1], d_o[2], d_n[0], d_n[1], d_n[2]
        };
    }

    return b_matrix_angles;
}

lgopt::vec2d<double> lgopt::buildBMatrixDihedralAngles(const xyz_coords_t &xyz_coords,
                                                       const RedIntCoords &red_int_coords,
                                                       const double tol)
{
    size_t n_dihedral_angles = red_int_coords.dihedral_angles_.size();

    using std::numbers::pi;

    vec2d<double> b_matrix_dihedrals(n_dihedral_angles);
    for (size_t idihedral = 0; idihedral < n_dihedral_angles; idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        // construct u, v, w and their relative angles
        std::array<double, 3> u = xyz_coords[m] - xyz_coords[o];
        std::array<double, 3> v = xyz_coords[n] - xyz_coords[p];
        std::array<double, 3> w = xyz_coords[p] - xyz_coords[o];
        double lambda_u = norm(u);
        double lambda_v = norm(v);
        double lambda_w = norm(w);
        u = u / lambda_u;
        v = v / lambda_v;
        w = w / lambda_w;
        double phi_u = angle(u, w);
        double phi_v = angle(v, w);
        if (std::fabs(phi_u) < tol || std::fabs(phi_v) < tol ||
            std::fabs(pi - phi_u) < tol || std::fabs(pi - phi_v) < tol)
            throw std::runtime_error(
                std::format("buildBMatrixDihedralAngles(): close to 180 degrees angle between "
                            "atoms: ({}, {}, {}, {})", m, o, p, n));

        double cos_phi_u = std::cos(phi_u);
        double cos_phi_v = std::cos(phi_v);
        double sin_phi_u = std::sin(phi_u);
        double sin_phi_v = std::sin(phi_v);
        double sin2_phi_u = sin_phi_u * sin_phi_u;
        double sin2_phi_v = sin_phi_v * sin_phi_v;
        std::array<double, 3> u_x_w = cross(u, w);
        std::array<double, 3> v_x_w = cross(v, w);

        // Lambda for calculating the eq. (34) with varying atomic index 'a'.
        auto calcContrib = [&](const size_t a) -> std::array<double, 3>
        {
            std::array<double, 3> result{};
            result[0] = zeta(a, m, o) * u_x_w[0] / (lambda_u * sin2_phi_u) +
                        zeta(a, p, n) * v_x_w[0] / (lambda_v * sin2_phi_v) +
                        zeta(a, o, p) * u_x_w[0] * cos_phi_u / (lambda_w * sin2_phi_u) -
                        zeta(a, o, p) * v_x_w[0] * cos_phi_v / (lambda_w * sin2_phi_v);

            result[1] = zeta(a, m, o) * u_x_w[1] / (lambda_u * sin2_phi_u) +
                        zeta(a, p, n) * v_x_w[1] / (lambda_v * sin2_phi_v) +
                        zeta(a, o, p) * u_x_w[1] * cos_phi_u / (lambda_w * sin2_phi_u) -
                        zeta(a, o, p) * v_x_w[1] * cos_phi_v / (lambda_w * sin2_phi_v);

            result[2] = zeta(a, m, o) * u_x_w[2] / (lambda_u * sin2_phi_u) +
                        zeta(a, p, n) * v_x_w[2] / (lambda_v * sin2_phi_v) +
                        zeta(a, o, p) * u_x_w[2] * cos_phi_u / (lambda_w * sin2_phi_u) -
                        zeta(a, o, p) * v_x_w[2] * cos_phi_v / (lambda_w * sin2_phi_v);

            return result;
        };

        std::array<double, 3> d_m = calcContrib(m);
        std::array<double, 3> d_o = calcContrib(o);
        std::array<double, 3> d_p = calcContrib(p);
        std::array<double, 3> d_n = calcContrib(n);

        b_matrix_dihedrals[idihedral] = {
            d_m[0], d_m[1], d_m[2], d_o[0], d_o[1], d_o[2], d_p[0], d_p[1], d_p[2], d_n[0], d_n[1],
            d_n[2]
        };
    }

    return b_matrix_dihedrals;
}
