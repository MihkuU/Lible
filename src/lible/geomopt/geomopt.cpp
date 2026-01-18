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
    static constexpr double bonding_factor = 1.3;
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

        double cos_u = std::cos(phi_u);
        double cos_v = std::cos(phi_v);
        double sin_u = std::sin(phi_u);
        double sin_v = std::sin(phi_v);
        double sin2_u = sin_u * sin_u;
        double sin2_v = sin_v * sin_v;
        std::array<double, 3> u_x_w = cross(u, w);
        std::array<double, 3> v_x_w = cross(v, w);

        // Lambda for calculating the eq. (34) with varying atomic index 'a'.
        auto calcContrib = [&](const size_t a) -> std::array<double, 3>
        {
            std::array<double, 3> result{};
            result[0] = zeta(a, m, o) * u_x_w[0] / (lambda_u * sin2_u) +
                        zeta(a, p, n) * v_x_w[0] / (lambda_v * sin2_v) +
                        zeta(a, o, p) * u_x_w[0] * cos_u / (lambda_w * sin2_u) -
                        zeta(a, o, p) * v_x_w[0] * cos_v / (lambda_w * sin2_v);

            result[1] = zeta(a, m, o) * u_x_w[1] / (lambda_u * sin2_u) +
                        zeta(a, p, n) * v_x_w[1] / (lambda_v * sin2_v) +
                        zeta(a, o, p) * u_x_w[1] * cos_u / (lambda_w * sin2_u) -
                        zeta(a, o, p) * v_x_w[1] * cos_v / (lambda_w * sin2_v);

            result[2] = zeta(a, m, o) * u_x_w[2] / (lambda_u * sin2_u) +
                        zeta(a, p, n) * v_x_w[2] / (lambda_v * sin2_v) +
                        zeta(a, o, p) * u_x_w[2] * cos_u / (lambda_w * sin2_u) -
                        zeta(a, o, p) * v_x_w[2] * cos_v / (lambda_w * sin2_v);

            return result;
        };

        std::array<double, 3> d_m = calcContrib(m);
        std::array<double, 3> d_o = calcContrib(o);
        std::array<double, 3> d_p = calcContrib(p);
        std::array<double, 3> d_n = calcContrib(n);

        b_matrix_dihedrals[idihedral] = {
            d_m[0], d_m[1], d_m[2],
            d_o[0], d_o[1], d_o[2],
            d_p[0], d_p[1], d_p[2],
            d_n[0], d_n[1], d_n[2]
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

lible::vecvec<double> lgopt::buildBMatrixDihedralAnglesFD(const double dx,
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

lible::vecvec<double> lgopt::buildBMatrixBondLengthsHD(const xyz_coords_t &xyz_coords,
                                                       const RedIntCoords &red_int_coords)
{
    size_t n_bond_lengths = red_int_coords.bond_lengths_.size();

    vecvec<double> b_matrix_bonds(n_bond_lengths);
    for (size_t ibond = 0; ibond < n_bond_lengths; ibond++)
    {
        const auto &[m, n, bond_length] = red_int_coords.bond_lengths_[ibond];

        std::array<double, 6> derivs = bondLengthGradient(xyz_coords[m], xyz_coords[n]);

        b_matrix_bonds[ibond] = {derivs[0], derivs[1], derivs[2], derivs[3], derivs[4], derivs[5]};
    }

    return b_matrix_bonds;
}

lible::vecvec<double> lgopt::buildBMatrixBondAnglesHD(const xyz_coords_t &xyz_coords,
                                                      const RedIntCoords &red_int_coords)
{
    size_t n_bond_angles = red_int_coords.bond_angles_.size();

    vecvec<double> b_matrix_bonds(n_bond_angles);
    for (size_t iangle = 0; iangle < n_bond_angles; iangle++)
    {
        const auto &[m, o, n, bond_angle] = red_int_coords.bond_angles_[iangle];

        std::array<double, 9> derivs = bondAngleGradient(xyz_coords[m], xyz_coords[o],
                                                         xyz_coords[n]);

        b_matrix_bonds[iangle] = {
            derivs[0], derivs[1], derivs[2], derivs[3], derivs[4], derivs[5],
            derivs[6], derivs[7], derivs[8]
        };
    }

    return b_matrix_bonds;
}

lible::vecvec<double> lgopt::buildBMatrixDihedralAnglesHD(const xyz_coords_t &xyz_coords,
                                                          const RedIntCoords &red_int_coords)
{
    size_t n_dihedral_angles = red_int_coords.dihedral_angles_.size();

    vecvec<double> b_matrix_dihedrals(n_dihedral_angles);
    for (size_t idihedral = 0; idihedral < n_dihedral_angles; idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        std::array<double, 12> derivs = dihedralAngleGradient(xyz_coords[m], xyz_coords[o],
                                                              xyz_coords[p], xyz_coords[n]);

        b_matrix_dihedrals[idihedral] = {
            derivs[0], derivs[1], derivs[2], derivs[3], derivs[4],
            derivs[5], derivs[6], derivs[7], derivs[8], derivs[9],
            derivs[10], derivs[11]
        };
    }

    return b_matrix_dihedrals;
}

lible::vec2d lgopt::buildKMatrix(const vecvec<double> &b_matrix,
                                 const std::vector<double> &grad_redint,
                                 const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords)
{
    return {};
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

        u = u / lambda_u;
        v = v / lambda_v;

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
                                b_matrix[idx][3 * ia + i] * b_matrix[idx][3 * ib + j];
                    }
            }
    }

    return k_matrix_angles;
}

lible::vec2d lgopt::buildKMatrixDihedralAngles(const std::vector<double> &grad_redint,
                                               const xyz_coords_t &xyz_coords,
                                               const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_dihedrals(Fill(0), n_coords, n_coords);
    for (size_t idihedral = 0; idihedral < red_int_coords.dihedral_angles_.size(); idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_p = xyz_coords[p];
        std::array<double, 3> coords_n = xyz_coords[n];

        std::array<double, 3> u = coords_m - coords_o;
        std::array<double, 3> v = coords_n - coords_p;
        std::array<double, 3> w = coords_p - coords_o;

        double lambda_u = norm(u);
        double lambda_v = norm(v);
        double lambda_w = norm(w);
        u = u / lambda_u;
        v = v / lambda_v;
        w = w / lambda_w;

        double cos_u = dot(u, w);
        double sin_u = std::sqrt(1 - std::pow(cos_u, 2));
        double cos_v = -dot(v, w);
        double sin_v = std::sqrt(1 - std::pow(cos_v, 2));

        std::array<double, 3> u_x_w = cross(u, w);
        std::array<double, 3> v_x_w = cross(v, w);

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_p = 3 * p;
        size_t ofs_n = 3 * n;
        size_t idx = red_int_coords.bond_lengths_.size() + red_int_coords.bond_angles_.size() +
                     idihedral;

        printf("m = %zu, o = %zu, p = %zu, n = %zu\n", m, o, p, n);

        std::array<size_t, 4> mopn{m, o, p, n};
        std::array<size_t, 4> ofs_mopn{ofs_m, ofs_o, ofs_p, ofs_n};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                // First 6 terms.
                double term1 = (u_x_w[i] * (w[j] * cos_u - u[j]) +
                                u_x_w[j] * (w[i] * cos_u - u[i])) /
                               (std::pow(lambda_u, 2) * std::pow(sin_u, 4));

                // TODO: comment change
                double term2 = (v_x_w[i] * (-w[j] * cos_v + v[j]) +
                                v_x_w[j] * (-w[i] * cos_v + v[i])) /
                               (std::pow(lambda_v, 2) * std::pow(sin_v, 4));

                double term3 = (u_x_w[i] * (w[j] - 2 * u[j] * cos_u + w[j] * std::pow(cos_u, 2)) +
                                u_x_w[j] * (w[i] - 2 * u[i] * cos_u + w[i] * std::pow(cos_u, 2))) /
                               (2 * lambda_u * lambda_w * std::pow(sin_u, 4));

                // TODO: comment change
                double term4 = (v_x_w[i] * (w[j] + 2 * v[j] * cos_v - w[j] * std::pow(cos_v, 2)) +
                                v_x_w[j] * (w[i] + 2 * v[i] * cos_v - w[i] * std::pow(cos_v, 2))) /
                               (2 * lambda_v * lambda_w * std::pow(sin_v, 4));

                double term5 = (u_x_w[i] * (u[j] + u[j] * std::pow(cos_u, 2) -
                                            3 * w[j] * cos_u + w[j] * std::pow(cos_u, 3)) +
                                u_x_w[j] * (u[i] + u[i] * std::pow(cos_u, 2) -
                                            3 * w[i] * cos_u + w[i] * std::pow(cos_u, 3))) /
                               (2 * std::pow(lambda_w, 2) * std::pow(sin_u, 4));

                // TODO: comment change
                double term6 = (v_x_w[i] * (-v[j] - v[j] * std::pow(cos_v, 2) +
                                            3 * w[j] * cos_v - w[j] * std::pow(cos_v, 3)) +
                                v_x_w[j] * (-v[i] - v[i] * std::pow(cos_v, 2) +
                                            3 * w[i] * cos_v - w[i] * std::pow(cos_v, 3))) /
                               (2 * std::pow(lambda_w, 2) * std::pow(sin_v, 4));

                for (int ia = 0; ia < 4; ia++)
                    for (int ib = 0; ib < 4; ib++)
                    {
                        size_t a = mopn[ia];
                        size_t b = mopn[ib];
                        size_t ofs_a = ofs_mopn[ia];
                        size_t ofs_b = ofs_mopn[ib];

                        // printf("a = %zu, b = %zu\n", a, b);
                        double deriv = zeta(a, m, o) * zeta(b, m, o) * term1;
                        deriv += zeta(a, n, p) * zeta(b, n, p) * term2;
                        deriv += (zeta(a, m, o) * zeta(b, o, p) +
                                  zeta(a, p, o) * zeta(b, o, m)) * term3;
                        deriv += (zeta(a, n, p) * zeta(b, p, o) +
                                  zeta(a, p, o) * zeta(b, n, p)) * term4;
                        deriv += zeta(a, o, p) * zeta(b, p, o) * term5;
                        deriv += zeta(a, p, o) * zeta(b, o, p) * term6; // TODO: comment change

                        // printf("deriv = %16.12lf, grad_redint[idx] = %16.12lf\n",
                        //        deriv, grad_redint[idx]);

                        k_matrix_dihedrals(ofs_a + i, ofs_b + j) += grad_redint[idx] * deriv;
                    }

                // Last two remaining terms.
                if (i != j)
                {
                    // std::set<size_t> ijk{0, 1, 2};
                    // ijk.erase(i);
                    // ijk.erase(j);
                    // size_t k = *ijk.begin();
                    int k = 3 - i - j;
                    // printf("k = %d, i = %d, j = %d\n", k, i, j);

                    // The sin's are squared unlike in the paper. Based on OptKing:
                    // https://github.com/psi-rking/optking.

                    double term7 = (j - i) * std::pow(-0.5, std::abs(j - i)) * (-w[k] * cos_u + u[k]) /
                                   (lambda_u * lambda_w * std::pow(sin_u, 2));

                    double term8 = (j - i) * std::pow(-0.5, std::abs(j - i)) * (w[k] * cos_v - v[k]) /
                                   (lambda_v * lambda_w * std::pow(sin_v, 2));

                    for (int ia = 0; ia < 4; ia++)
                        for (int ib = 0; ib < 4; ib++)
                        {
                            size_t a = mopn[ia];
                            size_t b = mopn[ib];
                            size_t ofs_a = ofs_mopn[ia];
                            size_t ofs_b = ofs_mopn[ib];

                            k_matrix_dihedrals(ofs_a + i, ofs_b + j) +=
                                    grad_redint[idx] * (1 - delta(a, b)) *
                                    (zeta(a, m, o) * zeta(b, p, o) + zeta(a, p, o) * zeta(b, o, m)) *
                                    term7; // TODO: comment change

                            k_matrix_dihedrals(ofs_a + i, ofs_b + j) +=
                                    grad_redint[idx] * (1 - delta(a, b)) *
                                    (zeta(a, n, o) * zeta(b, p, o) + zeta(a, p, o) * zeta(b, o, n)) *
                                    term8; // TODO: comment change
                        }
                }
            }
    }

    return k_matrix_dihedrals;
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
        for (int ia = 0; ia < 2; ia++)
            for (int ib = 0; ib < 2; ib++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        std::array<double, 6> uu = coords_mn;
                        std::array<double, 6> ud = coords_mn;
                        std::array<double, 6> du = coords_mn;
                        std::array<double, 6> dd = coords_mn;

                        size_t i_ = ia * 3 + i;
                        size_t j_ = ib * 3 + j;
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

                        k_matrix_bonds(ofs_mn[ia] + i, ofs_mn[ib] + j) += grad_redint[ibond] * deriv;
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
        for (int ia = 0; ia < 3; ia++)
            for (int ib = 0; ib < 3; ib++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        std::array<double, 9> uu = coords_mon;
                        std::array<double, 9> ud = coords_mon;
                        std::array<double, 9> du = coords_mon;
                        std::array<double, 9> dd = coords_mon;

                        size_t i_ = ia * 3 + i;
                        size_t j_ = ib * 3 + j;
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

                        size_t ofs_a = ofs_mon[ia];
                        size_t ofs_b = ofs_mon[ib];
                        k_matrix_angles(ofs_a + i, ofs_b + j) += grad_redint[idx] * deriv;
                    }
    }

    return k_matrix_angles;
}

lible::vec2d lgopt::buildKMatrixDihedralAnglesFD(const double dx, const double dy,
                                                 const std::vector<double> &grad_redint,
                                                 const xyz_coords_t &xyz_coords,
                                                 const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_dihedrals(Fill(0), n_coords, n_coords);
    for (size_t idihedral = 0; idihedral < red_int_coords.dihedral_angles_.size(); idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        std::array<double, 3> coords_m = xyz_coords[m];
        std::array<double, 3> coords_o = xyz_coords[o];
        std::array<double, 3> coords_p = xyz_coords[p];
        std::array<double, 3> coords_n = xyz_coords[n];

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_p = 3 * p;
        size_t ofs_n = 3 * n;
        size_t idx = red_int_coords.bond_lengths_.size() + red_int_coords.bond_angles_.size() +
                     idihedral;

        std::array<double, 12> coords_mopn{
            coords_m[0], coords_m[1], coords_m[2],
            coords_o[0], coords_o[1], coords_o[2],
            coords_p[0], coords_p[1], coords_p[2],
            coords_n[0], coords_n[1], coords_n[2]
        };

        std::array<size_t, 4> ofs_mopn{ofs_m, ofs_o, ofs_p, ofs_n};

        for (int ia = 0; ia < 4; ia++)
            for (int ib = 0; ib < 4; ib++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        std::array<double, 12> uu = coords_mopn;
                        std::array<double, 12> ud = coords_mopn;
                        std::array<double, 12> du = coords_mopn;
                        std::array<double, 12> dd = coords_mopn;

                        size_t i_ = ia * 3 + i;
                        size_t j_ = ib * 3 + j;
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
                        std::array<double, 3> uu_p{uu[6], uu[7], uu[8]};
                        std::array<double, 3> uu_n{uu[9], uu[10], uu[11]};

                        std::array<double, 3> ud_m{ud[0], ud[1], ud[2]};
                        std::array<double, 3> ud_o{ud[3], ud[4], ud[5]};
                        std::array<double, 3> ud_p{ud[6], ud[7], ud[8]};
                        std::array<double, 3> ud_n{ud[9], ud[10], ud[11]};

                        std::array<double, 3> du_m{du[0], du[1], du[2]};
                        std::array<double, 3> du_o{du[3], du[4], du[5]};
                        std::array<double, 3> du_p{du[6], du[7], du[8]};
                        std::array<double, 3> du_n{du[9], du[10], du[11]};

                        std::array<double, 3> dd_m{dd[0], dd[1], dd[2]};
                        std::array<double, 3> dd_o{dd[3], dd[4], dd[5]};
                        std::array<double, 3> dd_p{dd[6], dd[7], dd[8]};
                        std::array<double, 3> dd_n{dd[9], dd[10], dd[11]};

                        double deriv = (dihedralAngle(uu_m, uu_o, uu_p, uu_n) -
                                        dihedralAngle(ud_m, ud_o, ud_p, ud_n) -
                                        dihedralAngle(du_m, du_o, du_p, du_n) +
                                        dihedralAngle(dd_m, dd_o, dd_p, dd_n)) /
                                       (4 * dx * dy);

                        size_t ofs_a = ofs_mopn[ia];
                        size_t ofs_b = ofs_mopn[ib];
                        k_matrix_dihedrals(ofs_a + i, ofs_b + j) += grad_redint[idx] * deriv;
                    }
    }

    return k_matrix_dihedrals;
}

lible::vec2d lgopt::buildKMatrixBondLengthsHD(const std::vector<double> &grad_redint,
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
        std::array<size_t, 6> idxs{ofs_m, ofs_m + 1, ofs_m + 2, ofs_n, ofs_n + 1, ofs_n + 2};

        size_t idx_grad = ibond;

        arr2d<double, 6, 6> derivs = bondLengthHessian(xyz_coords[m], xyz_coords[n]);

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                k_matrix_bonds(idxs[i], idxs[j]) += derivs[i][j] * grad_redint[idx_grad];
    }

    return k_matrix_bonds;
}

lible::vec2d lgopt::buildKMatrixBondAnglesHD(const std::vector<double> &grad_redint,
                                             const xyz_coords_t &xyz_coords,
                                             const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_bonds(Fill(0), n_coords, n_coords);
    for (size_t iangle = 0; iangle < red_int_coords.bond_angles_.size(); iangle++)
    {
        const auto &[m, o, n, bond_angle] = red_int_coords.bond_angles_[iangle];

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_n = 3 * n;
        std::array<size_t, 9> idxs{
            ofs_m, ofs_m + 1, ofs_m + 2, ofs_o, ofs_o + 1, ofs_o + 2,
            ofs_n, ofs_n + 1, ofs_n + 2
        };

        size_t idx_grad = red_int_coords.bond_lengths_.size() + iangle;

        arr2d<double, 9, 9> derivs = bondAngleHessian(xyz_coords[m], xyz_coords[o], xyz_coords[n]);

        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 9; j++)
                k_matrix_bonds(idxs[i], idxs[j]) += derivs[i][j] * grad_redint[idx_grad];
    }

    return k_matrix_bonds;
}

lible::vec2d lgopt::buildKMatrixDihedralAnglesHD(const std::vector<double> &grad_redint,
                                                 const xyz_coords_t &xyz_coords,
                                                 const RedIntCoords &red_int_coords)
{
    size_t n_coords = 3 * xyz_coords.size();

    vec2d k_matrix_dihedrals(Fill(0), n_coords, n_coords);
    for (size_t idihedral = 0; idihedral < red_int_coords.dihedral_angles_.size(); idihedral++)
    {
        const auto &[m, o, p, n, dihedral_angle] = red_int_coords.dihedral_angles_[idihedral];

        size_t ofs_m = 3 * m;
        size_t ofs_o = 3 * o;
        size_t ofs_p = 3 * p;
        size_t ofs_n = 3 * n;
        std::array<size_t, 12> idxs{
            ofs_m, ofs_m + 1, ofs_m + 2, ofs_o, ofs_o + 1, ofs_o + 2,
            ofs_p, ofs_p + 1, ofs_p + 2, ofs_n, ofs_n + 1, ofs_n + 2
        };

        size_t idx_grad = red_int_coords.bond_lengths_.size() +
                          red_int_coords.bond_angles_.size() + idihedral;

        arr2d<double, 12, 12> derivs = dihedralAngleHessian(xyz_coords[m], xyz_coords[o],
                                                            xyz_coords[p], xyz_coords[n]);

        for (int i = 0; i < 12; i++)
            for (int j = 0; j < 12; j++)
                k_matrix_dihedrals(idxs[i], idxs[j]) += derivs[i][j] * grad_redint[idx_grad];
    }

    return k_matrix_dihedrals;
}
