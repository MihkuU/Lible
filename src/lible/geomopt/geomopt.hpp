#pragma once

#include <lible/types.hpp>

#include <array>
#include <vector>

namespace lible::geomopt
{
    /// Type alias for a list of xyz-coordinates.
    using xyz_coords_t = std::vector<std::array<double, 3>>;

    /// Tolerance used for judging various things: vector zero length, vector parallelity, etc.
    constexpr double tolerance = 1e-10;

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

    /// Constructs redundant internal coordinates from the given Cartesian coordinates.
    RedIntCoords redIntCoords(const std::vector<int> &atomic_nrs, const xyz_coords_t &xyz_coords);

    /// Returns the total number of redundant internal coordinates.
    size_t numRedIntCoords(const RedIntCoords &red_int_coords);

    /// Finds the bonding partner atoms for each atom.
    vecvec<size_t> bondingPartners(const std::vector<int> &atomic_nrs,
                                   const xyz_coords_t &xyz_coords);

    /// Constructs and returns a list of bond length coordinates.
    std::vector<BondLength> bondLengths(const vecvec<size_t> &bonding_partners,
                                        const xyz_coords_t &xyz_coords);

    /// Constructs and returns a list of bond angle coordinates.
    std::vector<BondAngle> bondAngles(const vecvec<size_t> &bonding_partners,
                                      const xyz_coords_t &xyz_coords);

    /// Constructs and returns a list of dihedral angle coordinates.
    std::vector<DihedralAngle> dihedralAngles(const vecvec<size_t> &bonding_partners,
                                              const xyz_coords_t &xyz_coords);

    // TODO: do this instead?
    struct BMatrix
    {
        std::vector<std::array<double, 6>> bmat_bonds_;
        std::vector<std::array<double, 9>> bmat_angles_;
        std::vector<std::array<double, 12>> bmat_dihedrals_;
    };

    /// Constructs the Wilson's B-matrix using formulas from https://doi.org/10.1063/1.1515483.
    vecvec<double> builBMatrix(const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords,
                               double tol = tolerance);

    /// Calculates the B-matrix for bond lenghts using eq. (17) from
    /// https://doi.org/10.1063/1.1515483.
    vecvec<double> buildBMatrixBondLengths(const xyz_coords_t &xyz_coords,
                                           const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond angles using eq. (25) from
    /// https://doi.org/10.1063/1.1515483.
    vecvec<double> buildBMatrixBondAngles(const xyz_coords_t &xyz_coords,
                                          const RedIntCoords &red_int_coords,
                                          double tol = tolerance);

    /// Constructs the B-matrix for dihedral angles using eq. (34) from
    /// https://doi.org/10.1063/1.1515483.
    vecvec<double> buildBMatrixDihedralAngles(const xyz_coords_t &xyz_coords,
                                              const RedIntCoords &red_int_coords,
                                              double tol = tolerance);

    /// Calculates the B-matrix for bond lengths using finite differences (FD).
    vecvec<double> buildBMatrixBondLengthsFD(double dx, const xyz_coords_t &xyz_coords,
                                             const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for bond angles using finite differences (FD).
    vecvec<double> buildBMatrixBondAnglesFD(double dx, const xyz_coords_t &xyz_coords,
                                            const RedIntCoords &red_int_coords);

    /// Calculates the B-matrix for dihedral angles using finite differences (FD).
    vecvec<double> buildBMatrixDihedralAnglesFD(double dx, const xyz_coords_t &xyz_coords,
                                                const RedIntCoords &red_int_coords);

    /// Constructs the K matrix defined in eq. (8) from https://doi.org/10.1063/1.1515483.
    vec2d buildKMatrix(const vecvec<double> &b_matrix, const std::vector<double> &grad_redint,
                       const xyz_coords_t &xyz_coords, const RedIntCoords &red_int_coords);

    // TODO:
    vec2d buildKMatrixBondLengths(const std::vector<double> &grad_redint,
                                  const xyz_coords_t &xyz_coords,
                                  const RedIntCoords &red_int_coords);

    // TODO:
    vec2d buildKMatrixBondAngles(const vecvec<double> &b_matrix,
                                 const std::vector<double> &grad_redint,
                                 const xyz_coords_t &xyz_coords,
                                 const RedIntCoords &red_int_coords);

    // TODO:
    vec2d buildKMatrixBondLengthsFD(double dx, double dy, const std::vector<double> &grad_redint,
                                    const xyz_coords_t &xyz_coords,
                                    const RedIntCoords &red_int_coords);

    // TODO:
    vec2d buildKMatrixBondAnglesFD(double dx, double dy, const std::vector<double> &grad_redint,
                                   const xyz_coords_t &xyz_coords,
                                   const RedIntCoords &red_int_coords);

    /// Calculates the bond length for two atoms.
    double bondLength(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_n);

    /// Calculates a bond angle between three atoms.
    double bondAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                     const std::array<double, 3> &xyz_n);

    /// Calculates a dihedral angle formed by four atoms. TODO: contains an angle error
    double dihedralAngle(const std::array<double, 3> &xyz_m, const std::array<double, 3> &xyz_o,
                         const std::array<double, 3> &xyz_p, const std::array<double, 3> &xyz_n);
}
