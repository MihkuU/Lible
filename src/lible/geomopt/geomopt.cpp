#include <lible/geomopt/geomopt.hpp>
#include <lible/geomopt/utils.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <numbers>
#include <stdexcept>

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
    vecvec<size_t> bonding_partners = bondingPartners(atomic_nrs, xyz_coords);

    return {
        bondLengths(bonding_partners, xyz_coords),
        bondAngles(bonding_partners, xyz_coords),
        dihedralAngles(bonding_partners, xyz_coords)
    };
}

std::vector<lgopt::BondLength> lgopt::bondLengths(const vecvec<size_t> &bonding_partners,
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

std::vector<lgopt::BondAngle> lgopt::bondAngles(const vecvec<size_t> &bonding_partners,
                                                const xyz_coords_t &xyz_coords)
{
    std::vector<BondAngle> bond_angles;
    for (size_t o = 0; o < bonding_partners.size(); o++)
        for (size_t m : bonding_partners.at(o))
            for (size_t n : bonding_partners.at(o))
                if (m < n)
                {
                    double bond_angle = bondAngle(xyz_coords[m], xyz_coords[o], xyz_coords[n]);
                    bond_angles.push_back({m, o, n, bond_angle});
                }

    return bond_angles;
}

std::vector<lgopt::DihedralAngle> lgopt::dihedralAngles(const vecvec<size_t> &bonding_partners,
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
                    dihedral_angles.push_back({m, o, p, n, dihedral_angle});
                }
    }

    return dihedral_angles;
}

size_t lgopt::numRedIntCoords(const RedIntCoords &red_int_coords)
{
    const auto &[bonds, angles, dihedrals] = red_int_coords;
    return bonds.size() + angles.size() + dihedrals.size();
}

lible::vecvec<size_t> lgopt::bondingPartners(const std::vector<int> &atomic_nrs,
                                             const xyz_coords_t &xyz_coords)
{
    size_t n_atoms = atomic_nrs.size();

    vecvec<size_t> bonding_partners(n_atoms);
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

    double arg = dot(cross(u, w), cross(v, w)) / (sin_u * sin_v);
    arg = std::clamp(arg, -1.0, 1.0);

    return std::acos(arg);

    // ///////////////////////////////////////////
    // std::array<double, 3> b1 = xyz_o - xyz_m;
    // std::array<double, 3> b2 = xyz_p - xyz_o;
    // std::array<double, 3> b3 = xyz_n - xyz_p;
    //
    // std::array<double, 3> n1 = cross(b1, b2);
    // std::array<double, 3> n2 = cross(b2, b3);
    // n1 = n1 / norm(n1);
    // n2 = n2 / norm(n2);
    //
    // std::array<double, 3> m1 = cross(n1, b2 / norm(b2));
    //
    // double x = dot(n1, n2);
    // double y = dot(m1, n2);
    //
    // // printf("y = %16.12lf, x = %16.12lf, norm(b2) = %16.12lf, std::atan2(y, x) = %16.12lf\n",
    // //     y, x, norm(b2), std::atan2(y, x));
    //
    // return std::atan2(y, x);
}

lible::vecvec<double> lgopt::builBMatrix(const xyz_coords_t &xyz_coords,
                                         const RedIntCoords &red_int_coords,
                                         const double tol)
{
    vecvec<double> b_matrix_bonds = buildBMatrixBondLengths(xyz_coords, red_int_coords);

    vecvec<double> b_matrix_angles = buildBMatrixBondAngles(xyz_coords, red_int_coords, tol);

    vecvec<double> b_matrix_dihedrals = buildBMatrixDihedralAngles(xyz_coords, red_int_coords, tol);

    vecvec<double> b_matrix = b_matrix_bonds;
    b_matrix.insert(b_matrix.end(), b_matrix_angles.begin(), b_matrix_angles.end());
    b_matrix.insert(b_matrix.end(), b_matrix_dihedrals.begin(), b_matrix_dihedrals.end());

    return b_matrix;
}

lible::vecvec<double> lgopt::buildBMatrixBondLengths(const xyz_coords_t &xyz_coords,
                                                     const RedIntCoords &red_int_coords)
{
    size_t n_bond_lenghts = red_int_coords.bond_lengths_.size();

    vecvec<double> b_matrix_bonds(n_bond_lenghts);
    for (size_t ibond = 0; ibond < n_bond_lenghts; ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        std::array<double, 3> bond_vector = (xyz_coords[m] - xyz_coords[n]) / bond_length;
        const auto &[x, y, z] = bond_vector;
        b_matrix_bonds[ibond] = {x, y, z, -x, -y, -z};
    }

    return b_matrix_bonds;
}

lible::vecvec<double> lgopt::buildBMatrixBondAngles(const xyz_coords_t &xyz_coords,
                                                    const RedIntCoords &red_int_coords,
                                                    const double tol)
{
    size_t n_bond_angles = red_int_coords.bond_angles_.size();

    vecvec<double> b_matrix_angles(n_bond_angles);
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
                tmp = {-1, 1, 1};
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

lible::vecvec<double> lgopt::buildBMatrixDihedralAngles(const xyz_coords_t &xyz_coords,
                                                        const RedIntCoords &red_int_coords,
                                                        const double tol)
{
    size_t n_dihedral_angles = red_int_coords.dihedral_angles_.size();

    using std::numbers::pi;

    vecvec<double> b_matrix_dihedrals(n_dihedral_angles);
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

lible::vecvec<double> lgopt::buildBMatrixBondLengthsFD(const double dx,
                                                       const xyz_coords_t &xyz_coords,
                                                       const RedIntCoords &red_int_coords)
{
    size_t n_bond_lenghts = red_int_coords.bond_lengths_.size();

    vecvec<double> b_matrix_bonds(n_bond_lenghts);
    for (size_t ibond = 0; ibond < n_bond_lenghts; ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::vector<double> b_matrix_row(6, 0);
        // m
        for (int im = 0; im < 3; im++)
        {
            std::array<double, 3> coords_m_dx = coords_m;
            coords_m_dx[im] = coords_m[im] + dx;

            double bond_length_dx = bondLength(coords_m_dx, coords_n);

            double deriv = (bond_length_dx - bond_length) / dx;

            b_matrix_row[im] = deriv;
        }

        // n
        for (int in = 0; in < 3; in++)
        {
            std::array<double, 3> coords_n_dx = coords_n;
            coords_n_dx[in] = coords_n[in] + dx;

            double bond_length_dx = bondLength(coords_m, coords_n_dx);

            double deriv = (bond_length_dx - bond_length) / dx;

            b_matrix_row[3 + in] = deriv;
        }
        b_matrix_bonds[ibond] = b_matrix_row;
    }

    return b_matrix_bonds;
}

lible::vecvec<double> lgopt::buildBMatrixBondAnglesFD(const double dx,
                                                      const xyz_coords_t &xyz_coords,
                                                      const RedIntCoords &red_int_coords)
{
    size_t n_bond_angles = red_int_coords.bond_angles_.size();

    vecvec<double> b_matrix_angles(n_bond_angles);
    for (size_t iangle = 0; iangle < n_bond_angles; iangle++)
    {
        const auto &[m, o, n, bond_angle] = red_int_coords.bond_angles_[iangle];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::vector<double> b_matrix_row(9, 0);
        for (int im = 0; im < 3; im++)
        {
            std::array<double, 3> coords_m_dx_up = coords_m;
            std::array<double, 3> coords_m_dx_dn = coords_m;
            coords_m_dx_up[im] = coords_m[im] + dx;
            coords_m_dx_dn[im] = coords_m[im] - dx;

            double angle_dx_up = bondAngle(coords_m_dx_up, coords_o, coords_n);
            double angle_dx_dn = bondAngle(coords_m_dx_dn, coords_o, coords_n);

            double deriv = (angle_dx_up - angle_dx_dn) / (2 * dx);

            b_matrix_row[im] = deriv;
        }

        for (int io = 0; io < 3; io++)
        {
            std::array<double, 3> coords_o_dx_up = coords_o;
            std::array<double, 3> coords_o_dx_dn = coords_o;
            coords_o_dx_up[io] = coords_o[io] + dx;
            coords_o_dx_dn[io] = coords_o[io] - dx;

            double angle_dx_up = bondAngle(coords_m, coords_o_dx_up, coords_n);
            double angle_dx_dn = bondAngle(coords_m, coords_o_dx_dn, coords_n);

            double deriv = (angle_dx_up - angle_dx_dn) / (2 * dx);

            b_matrix_row[3 + io] = deriv;
        }

        for (int in = 0; in < 3; in++)
        {
            std::array<double, 3> coords_n_dx_up = coords_n;
            std::array<double, 3> coords_n_dx_dn = coords_n;
            coords_n_dx_up[in] = coords_n[in] + dx;
            coords_n_dx_dn[in] = coords_n[in] - dx;

            double angle_dx_up = bondAngle(coords_m, coords_o, coords_n_dx_up);
            double angle_dx_dn = bondAngle(coords_m, coords_o, coords_n_dx_dn);

            double deriv = (angle_dx_up - angle_dx_dn) / (2 * dx);

            b_matrix_row[6 + in] = deriv;
        }

        b_matrix_angles[iangle] = b_matrix_row;
    }

    return b_matrix_angles;
}

lible::vecvec<double> lgopt::buildBMatrixDihedralAnglesFD(double dx,
                                                          const xyz_coords_t &xyz_coords,
                                                          const RedIntCoords &red_int_coords)
{
    size_t n_dihedral_angles = red_int_coords.dihedral_angles_.size();

    vecvec<double> b_matrix_dihedrals(n_dihedral_angles);
    for (size_t idihedral = 0; idihedral < n_dihedral_angles; idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_p = xyz_coords[p];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::vector<double> b_matrix_row(12, 0);
        for (int im = 0; im < 3; im++)
        {
            std::array<double, 3> coords_m_dx_up = coords_m;
            std::array<double, 3> coords_m_dx_dn = coords_m;
            coords_m_dx_up[im] += dx;
            coords_m_dx_dn[im] -= dx;

            double dihedral_dx_up = dihedralAngle(coords_m_dx_up, coords_o, coords_p, coords_n);
            double dihedral_dx_dn = dihedralAngle(coords_m_dx_dn, coords_o, coords_p, coords_n);

            double deriv = (dihedral_dx_up - dihedral_dx_dn) / (2 * dx);

            b_matrix_row[im] = deriv;
        }

        for (int io = 0; io < 3; io++)
        {
            std::array<double, 3> coords_o_dx_up = coords_o;
            std::array<double, 3> coords_o_dx_dn = coords_o;
            coords_o_dx_up[io] += dx;
            coords_o_dx_dn[io] -= dx;

            double dihedral_dx_up = dihedralAngle(coords_m, coords_o_dx_up, coords_p, coords_n);
            double dihedral_dx_dn = dihedralAngle(coords_m, coords_o_dx_dn, coords_p, coords_n);

            double deriv = (dihedral_dx_up - dihedral_dx_dn) / (2 * dx);

            b_matrix_row[3 + io] = deriv;
        }

        for (int ip = 0; ip < 3; ip++)
        {
            std::array<double, 3> coords_p_dx_up = coords_p;
            std::array<double, 3> coords_p_dx_dn = coords_p;
            coords_p_dx_up[ip] += dx;
            coords_p_dx_dn[ip] -= dx;

            double dihedral_dx_up = dihedralAngle(coords_m, coords_o, coords_p_dx_up, coords_n);
            double dihedral_dx_dn = dihedralAngle(coords_m, coords_o, coords_p_dx_dn, coords_n);

            double deriv = (dihedral_dx_up - dihedral_dx_dn) / (2 * dx);

            b_matrix_row[6 + ip] = deriv;
        }

        for (int in = 0; in < 3; in++)
        {
            std::array<double, 3> coords_n_dx_up = coords_n;
            std::array<double, 3> coords_n_dx_dn = coords_n;
            coords_n_dx_up[in] += dx;
            coords_n_dx_dn[in] -= dx;

            double dihedral_dx_up = dihedralAngle(coords_m, coords_o, coords_p, coords_n_dx_up);
            double dihedral_dx_dn = dihedralAngle(coords_m, coords_o, coords_p, coords_n_dx_dn);

            double deriv = (dihedral_dx_up - dihedral_dx_dn) / (2 * dx);

            b_matrix_row[9 + in] = deriv;
        }

        b_matrix_dihedrals[idihedral] = b_matrix_row;
    }

    return b_matrix_dihedrals;
}

lible::vec2d lgopt::buildKMatrix(const vecvec<double> &b_matrix,
                                 const std::vector<double> &grad_redint,
                                 const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords)
{
}

lible::vec2d lgopt::buildKMatrixBondLengths(const std::vector<double> &grad_redint,
                                            const xyz_coords_t &xyz_coords,
                                            const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_bonds(Fill(0), n_coords, n_coords);
    for (size_t ibond = 0; ibond < red_int_coords.bond_lengths_.size(); ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        std::array<double, 3> bond_vector = xyz_coords[m] - xyz_coords[n];
        bond_vector = bond_vector / norm(bond_vector);

        size_t ofs_m = 3 * m;
        size_t ofs_n = 3 * n;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double val = grad_redint[ibond] * (bond_vector[i] * bond_vector[j] - delta(i, j)) /
                             bond_length;

                k_matrix_bonds(ofs_m + i, ofs_m + j) -= val;
                k_matrix_bonds(ofs_m + i, ofs_n + j) += val;
                k_matrix_bonds(ofs_n + i, ofs_m + j) += val;
                k_matrix_bonds(ofs_n + i, ofs_n + j) -= val;
            }
    }

    return k_matrix_bonds;
}

lible::vec2d lgopt::buildKMatrixBondAngles(const vecvec<double> &b_matrix,
                                           const std::vector<double> &grad_redint,
                                           const xyz_coords_t &xyz_coords,
                                           const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_angles(Fill(0), n_coords, n_coords);
    for (size_t iangle = 0; iangle < red_int_coords.bond_angles_.size(); iangle++)
    {
        const auto &[m, o, n, bond_angle] = red_int_coords.bond_angles_[iangle];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::array<double, 3> u = coords_m - coords_o;
        std::array<double, 3> v = coords_n - coords_o;

        double lambda_u = norm(u);
        double lambda_v = norm(v);

        u = u / norm(u);
        v = v / norm(v);

        double cos_q = dot(u, v);
        double sin_q = std::sqrt(1 - std::pow(cos_q, 2));

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_n = 3 * n;
        size_t idx = red_int_coords.bond_lengths_.size() + iangle;

        std::array<size_t, 3> mon{m, o, n};
        std::array<size_t, 3> ofs_mon{ofs_m, ofs_o, ofs_n};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double term1 = (u[i] * v[j] + u[j] * v[i] - 3 * u[i] * u[j] * cos_q +
                                delta(i, j) * cos_q) / (std::pow(lambda_u, 2) * sin_q);

                double term2 = (v[i] * u[j] + v[j] * u[i] - 3 * v[i] * v[j] * cos_q +
                                delta(i, j) * cos_q) / (std::pow(lambda_v, 2) * sin_q);

                double term3 = (u[i] * u[j] + v[j] * v[i] - u[i] * v[j] * cos_q - delta(i, j)) /
                               (lambda_u * lambda_v * sin_q);

                double term4 = (v[i] * v[j] + u[j] * u[i] - v[i] * u[j] * cos_q - delta(i, j)) /
                               (lambda_u * lambda_v * sin_q);

                for (int ia = 0; ia < 3; ia++)
                    for (int ib = 0; ib < 3; ib++)
                    {
                        size_t a = mon[ia];
                        size_t b = mon[ib];
                        size_t ofs_a = ofs_mon[ia];
                        size_t ofs_b = ofs_mon[ib];

                        // zeta terms
                        k_matrix_angles(ofs_a + i, ofs_b + j) +=
                            grad_redint[idx] * zeta(a, m, o) * zeta(b, m, o) * term1;
                        k_matrix_angles(ofs_a + i, ofs_b + j) +=
                            grad_redint[idx] * zeta(a, n, o) * zeta(b, n, o) * term2;
                        k_matrix_angles(ofs_a + i, ofs_b + j) +=
                            grad_redint[idx] * zeta(a, m, o) * zeta(b, n, o) * term3;
                        k_matrix_angles(ofs_a + i, ofs_b + j) +=
                            grad_redint[idx] * zeta(a, n, o) * zeta(b, m, o) * term4;

                        // B-matrix term
                        k_matrix_angles(ofs_a + i, ofs_b + j) -=
                                grad_redint[idx] * (cos_q / sin_q) *
                                b_matrix[idx][ia] * b_matrix[idx][ib];
                    }
            }
    }

    return k_matrix_angles;
}

lible::vec2d lgopt::buildKMatrixBondLengthsFD(const double dx, const double dy,
                                              const std::vector<double> &grad_redint,
                                              const xyz_coords_t &xyz_coords,
                                              const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_bonds(Fill(0), n_coords, n_coords);
    for (size_t ibond = 0; ibond < red_int_coords.bond_lengths_.size(); ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        size_t ofs_m = 3 * m;
        size_t ofs_n = 3 * n;

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::array<double, 6> coords_mn{
            coords_m[0], coords_m[1], coords_m[2],
            coords_n[0], coords_n[1], coords_n[2]
        };

        std::array<size_t, 2> ofs_mn{ofs_m, ofs_n};
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        std::array<double, 6> uu = coords_mn;
                        std::array<double, 6> ud = coords_mn;
                        std::array<double, 6> du = coords_mn;
                        std::array<double, 6> dd = coords_mn;

                        size_t i_ = a * 3 + i;
                        size_t j_ = b * 3 + j;
                        uu[i_] += dx;
                        uu[j_] += dy;

                        ud[i_] += dx;
                        ud[j_] -= dy;

                        du[i_] -= dx;
                        du[j_] += dy;

                        dd[i_] -= dx;
                        dd[j_] -= dy;

                        std::array<double, 3> uu_a{uu[0], uu[1], uu[2]};
                        std::array<double, 3> uu_b{uu[3], uu[4], uu[5]};

                        std::array<double, 3> ud_a{ud[0], ud[1], ud[2]};
                        std::array<double, 3> ud_b{ud[3], ud[4], ud[5]};

                        std::array<double, 3> du_a{du[0], du[1], du[2]};
                        std::array<double, 3> du_b{du[3], du[4], du[5]};

                        std::array<double, 3> dd_a{dd[0], dd[1], dd[2]};
                        std::array<double, 3> dd_b{dd[3], dd[4], dd[5]};

                        double deriv = (bondLength(uu_a, uu_b) - bondLength(ud_a, ud_b) -
                                        bondLength(du_a, du_b) + bondLength(dd_a, dd_b)) /
                                       (4 * dx * dy);

                        k_matrix_bonds(ofs_mn[a] + i, ofs_mn[b] + j) += grad_redint[ibond] * deriv;
                    }
    }

    return k_matrix_bonds;
}

lible::vec2d lgopt::buildKMatrixBondAnglesFD(const double dx, const double dy,
                                             const std::vector<double> &grad_redint,
                                             const xyz_coords_t &xyz_coords,
                                             const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_angles(Fill(0), n_coords, n_coords);
    for (size_t iangle = 0; iangle < red_int_coords.bond_angles_.size(); iangle++)
    {
        const auto &[m, o, n, bond_length] = red_int_coords.bond_angles_[iangle];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_n = xyz_coords[n];

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_n = 3 * n;
        size_t idx = red_int_coords.bond_lengths_.size() + iangle;

        std::array<double, 9> coords_mon{
            coords_m[0], coords_m[1], coords_m[2],
            coords_o[0], coords_o[1], coords_o[2],
            coords_n[0], coords_n[1], coords_n[2]
        };

        std::array<size_t, 3> ofs_mon{ofs_m, ofs_o, ofs_n};
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        std::array<double, 9> uu = coords_mon;
                        std::array<double, 9> ud = coords_mon;
                        std::array<double, 9> du = coords_mon;
                        std::array<double, 9> dd = coords_mon;

                        size_t i_ = a * 3 + i;
                        size_t j_ = b * 3 + j;
                        uu[i_] += dx;
                        uu[j_] += dy;

                        ud[i_] += dx;
                        ud[j_] -= dy;

                        du[i_] -= dx;
                        du[j_] += dy;

                        dd[i_] -= dx;
                        dd[j_] -= dy;

                        std::array<double, 3> uu_m{uu[0], uu[1], uu[2]};
                        std::array<double, 3> uu_o{uu[3], uu[4], uu[5]};
                        std::array<double, 3> uu_n{uu[6], uu[7], uu[8]};

                        std::array<double, 3> ud_m{ud[0], ud[1], ud[2]};
                        std::array<double, 3> ud_o{ud[3], ud[4], ud[5]};
                        std::array<double, 3> ud_n{ud[6], ud[7], ud[8]};

                        std::array<double, 3> du_m{du[0], du[1], du[2]};
                        std::array<double, 3> du_o{du[3], du[4], du[5]};
                        std::array<double, 3> du_n{du[6], du[7], du[8]};

                        std::array<double, 3> dd_m{dd[0], dd[1], dd[2]};
                        std::array<double, 3> dd_o{dd[3], dd[4], dd[5]};
                        std::array<double, 3> dd_n{dd[6], dd[7], dd[8]};

                        double deriv = (bondAngle(uu_m, uu_o, uu_n) - bondAngle(ud_m, ud_o, ud_n) -
                                        bondAngle(du_m, du_o, du_n) + bondAngle(dd_m, dd_o, dd_n)) /
                                       (4 * dx * dy);

                        size_t ofs_a = ofs_mon[a];
                        size_t ofs_b = ofs_mon[b];
                        k_matrix_angles(ofs_a + i, ofs_b + j) += grad_redint[idx] * deriv;
                    }
    }

    return k_matrix_angles;
}
