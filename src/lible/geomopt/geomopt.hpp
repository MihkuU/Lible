#pragma once

#include <lible/geomopt/utils.hpp> // TODO: remove this inclusion?

#include <vector>

namespace lible::geomopt
{
    /// Type alias for a list of xyz-coordinates.
    using xyz_coords_t = std::vector<std::array<double, 3>>;

    /// Structure representing the bond length coordinate between atoms m_ and n_.
    struct BondLength
    {
        size_t m_;
        size_t n_;
        double val_;
    };

    /// Structure representing the bond angle coordinate between atoms m_, o_, n_.
    struct BondAngle
    {
        size_t m_;
        size_t o_;
        size_t n_;
        double val_;
    };

    /// Structure representing the dihedral angle coordinate between atoms m_, o_, p_ and n_.
    struct DihedralAngle
    {
        size_t m_;
        size_t o_;
        size_t p_;
        size_t n_;
        double val_;
    };

    /// Structure representing the redundant internal coordinates: bond lengths, bond angles
    /// and dihedral angles.
    struct RedIntCoords
    {
        std::vector<BondLength> bond_lengths_;
        std::vector<BondAngle> bond_angles_;
        std::vector<DihedralAngle> dihedral_angles_;
    };

    /// Constructs redundant internal coordinates from the given Cartesian coordinates. Expects
    /// coordinates in Bohr (a.u.). The angles are returned in radians.
    RedIntCoords redIntCoords(const std::vector<int> &atomic_nrs, const xyz_coords_t &xyz_coords);

    /// Returns the total number of redundant internal coordinates.
    size_t numRedIntCoords(const RedIntCoords &red_int_coords);

    /// Finds the bonding partner atoms for each atom. Expects coordinates in Bohr (a.u.).
    vecvec<size_t> bondingPartners(const std::vector<int> &atomic_nrs,
                                   const xyz_coords_t &xyz_coords_au);

    /// Calculates the bond length for two atoms. Based on FIG 1 from
    /// https://doi.org/10.1063/1.1515483.
    double bondLength(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_n);

    /// Calculates a bond angle between three atoms. Using eq. (23) from
    /// https://doi.org/10.1063/1.1515483.
    double bondAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                     const std::array<double, 3> &xyz_n);

    /// Calculates a dihedral angle formed by four atoms. Using eq. (31) from
    /// https://doi.org/10.1063/1.1515483.
    double dihedralAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                         const std::array<double, 3> &xyz_p,
                         const std::array<double, 3> &xyz_n); // TODO: contains an angle error?

    /// Constructs and returns a list of bond length coordinates. Expects coordinates in Bohr
    /// (a.u.).
    std::vector<BondLength> bondLengths(const vecvec<size_t> &bonding_partners,
                                        const xyz_coords_t &xyz_coords_au);

    /// Constructs and returns a list of bond angle coordinates (in radians). Expects coordinates
    /// in Bohr (a.u.).
    std::vector<BondAngle> bondAngles(const vecvec<size_t> &bonding_partners,
                                      const xyz_coords_t &xyz_coords_au);

    /// Constructs and returns a list of dihedral angle coordinates (in radians). Expects
    /// coordinates in Bohr (a.u.).
    std::vector<DihedralAngle> dihedralAngles(const vecvec<size_t> &bonding_partners,
                                              const xyz_coords_t &xyz_coords_au);

    /// Calculates the first derivatives of a bond length using eq. (17) from
    /// https://doi.org/10.1063/1.1515483.
    std::array<double, 6> bondLengthGradient(const std::array<double, 3> &xyz_m,
                                             const std::array<double, 3> &xyz_n);

    /// Calculates the first derivatives of a bond angle using eq. (25) from
    /// https://doi.org/10.1063/1.1515483.
    std::array<double, 9> bondAngleGradient(const std::array<double, 3> &xyz_m,
                                            const std::array<double, 3> &xyz_o,
                                            const std::array<double, 3> &xyz_n,
                                            double paralellity_tol = tolerance);

    /// Calculates the first derivatives of a dihedral angle using eq. (34) from
    /// https://doi.org/10.1063/1.1515483.
    std::array<double, 12> dihedralAngleGradient(const std::array<double, 3> &xyz_m,
                                                 const std::array<double, 3> &xyz_o,
                                                 const std::array<double, 3> &xyz_p,
                                                 const std::array<double, 3> &xyz_n,
                                                 double linearity_tol = tolerance);

    /// Calculates the second derivatives of a bond length using eq. (20) from
    /// https://doi.org/10.1063/1.1515483.
    arr2d<double, 6, 6> bondLengthHessian(const std::array<double, 3> &xyz_m,
                                          const std::array<double, 3> &xyz_n);

    /// Calculates the second derivatives of a bond angle using eq. (27) from
    /// https://doi.org/10.1063/1.1515483.
    arr2d<double, 9, 9> bondAngleHessian(const std::array<double, 9> &bond_angle_gradient,
                                         const std::array<double, 3> &xyz_m,
                                         const std::array<double, 3> &xyz_o,
                                         const std::array<double, 3> &xyz_n);

    /// Calculates the second derivatives of a dihedral angle using hyper-dual numbers because
    /// the eq. (35) in https://doi.org/10.1063/1.1515483 contains errors.
    arr2d<double, 12, 12> dihedralAngleHessian(const std::array<double, 3> &xyz_m,
                                               const std::array<double, 3> &xyz_o,
                                               const std::array<double, 3> &xyz_p,
                                               const std::array<double, 3> &xyz_n);

    // B-matrix

    /// Structure representing the Wilson's B-matrix.
    struct BMatrix
    {
        /// Bond length Cartesian derivatives for m, n.
        std::vector<std::array<double, 6>> bmat_bonds_;
        /// Bond angle Cartesian derivatives for m, o, n.
        std::vector<std::array<double, 9>> bmat_angles_;
        /// Dihedral angle Cartesian derivatives for m, o, p, n.
        std::vector<std::array<double, 12>> bmat_dihedrals_;
    };

    /// Constructs the Wilson's B-matrix using formulas from https://doi.org/10.1063/1.1515483.
    BMatrix builBMatrix(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords,
                        double tol = tolerance);

    /// Calculates the B-matrix for bond lenghts using eq. (17) from
    /// https://doi.org/10.1063/1.1515483.
    std::vector<std::array<double, 6>>
    buildBMatrixBondLengths(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond angles using eq. (25) from
    /// https://doi.org/10.1063/1.1515483.
    std::vector<std::array<double, 9>>
    buildBMatrixBondAngles(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords,
                           double tol = tolerance);

    /// Constructs the B-matrix for dihedral angles using eq. (34) from
    /// https://doi.org/10.1063/1.1515483.
    std::vector<std::array<double, 12>>
    buildBMatrixDihedralAngles(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords,
                               double tol = tolerance);

    /// Calculates the B-matrix for bond lengths using finite differences (FD).
    std::vector<std::array<double, 6>>
    buildBMatrixBondLengthsFD(double dx, const xyz_coords_t &xyz_coords,
                              const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond angles using finite differences (FD).
    std::vector<std::array<double, 9>>
    buildBMatrixBondAnglesFD(double dx, const xyz_coords_t &xyz_coords,
                             const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for dihedral angles using finite differences (FD).
    std::vector<std::array<double, 12>>
    buildBMatrixDihedralAnglesFD(double dx, const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond lengths using hyper-dual numbers.
    std::vector<std::array<double, 6>>
    buildBMatrixBondLengthsHD(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond angles using hyper-dual numbers.
    std::vector<std::array<double, 9>>
    buildBMatrixBondAnglesHD(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for dihedral angles using hyper-dual numbers.
    std::vector<std::array<double, 12>>
    buildBMatrixDihedralAnglesHD(const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords);

    // K-matrix

    /// Constructs the K matrix defined in eq. (8) from https://doi.org/10.1063/1.1515483.
    vec2d buildKMatrix(const BMatrix &b_matrix, const std::vector<double> &grad_redint,
                       const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond lengths using eq. (20) from
    /// https://doi.org/10.1063/1.1515483.
    vec2d buildKMatrixBondLengths(const std::vector<double> &grad_redint,
                                  const xyz_coords_t &xyz_coords,
                                  const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond angles using eq. (27) from
    /// https://doi.org/10.1063/1.1515483.
    vec2d buildKMatrixBondAngles(const BMatrix &b_matrix, const std::vector<double> &grad_redint,
                                 const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for dihedral angles using hyper-dual numbers.
    vec2d buildKMatrixDihedralAngles(const std::vector<double> &grad_redint,
                                     const xyz_coords_t &xyz_coords,
                                     const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond lengths using finite differences.
    vec2d buildKMatrixBondLengthsFD(double dx, double dy, const std::vector<double> &grad_redint,
                                    const xyz_coords_t &xyz_coords,
                                    const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond angles using finite differences.
    vec2d buildKMatrixBondAnglesFD(double dx, double dy, const std::vector<double> &grad_redint,
                                   const xyz_coords_t &xyz_coords,
                                   const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for dihedral angles using finite differences.
    vec2d buildKMatrixDihedralAnglesFD(double dx, double dy, const std::vector<double> &grad_redint,
                                       const xyz_coords_t &xyz_coords,
                                       const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond lengths using hyper-dual numbers.
    vec2d buildKMatrixBondLengthsHD(const std::vector<double> &grad_redint,
                                    const xyz_coords_t &xyz_coords,
                                    const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for bond angles using hyper-dual numbers.
    vec2d buildKMatrixBondAnglesHD(const std::vector<double> &grad_redint,
                                   const xyz_coords_t &xyz_coords,
                                   const RedIntCoords &red_int_coords);

    /// Constructs the K matrix for dihedral angles using hyper-dual numbers.
    vec2d buildKMatrixDihedralAnglesHD(const std::vector<double> &grad_redint,
                                       const xyz_coords_t &xyz_coords,
                                       const RedIntCoords &red_int_coords);
}
