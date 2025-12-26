#pragma once

#include <array>
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

    /// Constructs redundant internal coordinates from the given Cartesian coordinates.
    RedIntCoords redIntCoords(const std::vector<int> &atomic_nrs, const xyz_coords_t &xyz_coords);

    /// Returns the total number of redundant internal coordinates.
    size_t numRedIntCoords(const RedIntCoords &red_int_coords);
}
