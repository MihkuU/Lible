#pragma once

#include <lible/ints/shell.hpp>

namespace lible::ints
{
    /// Structure containing contiguous data from shells for calculating integrals.
    struct ShellData
    {
        /// Constructor for the shell data. Expects that the angular momentum in `shells`
        /// equals `l`. Cannot be called in an OMP parallel region.
        ShellData(int l, const std::vector<Shell> &shells);

        /// Angular momentum.
        int l_{};
        /// Number of shells.
        size_t n_shells_{};
        /// Number of Gaussian primitives.
        size_t n_primitives_{};

        /// Contraction coefficients with the primitive norms multiplied into.
        std::vector<double> coeffs_;
        /// Coordinates of atoms corresponding to shells.
        std::vector<double> coords_;
        /// Exponents of the Gaussian primitives.
        std::vector<double> exps_;
        /// Normalization constants of the atomic orbitals in the spherical basis.
        std::vector<double> norms_;

        /// Indices of atoms involved in the shells.
        std::vector<size_t> atomic_idxs_;
        /// Contraction depth for each shell.
        std::vector<size_t> cdepths_;
        /// Offsets of the contraction data (coeffs and exps) for each shell.
        std::vector<size_t> coffsets_;
        /// Offsets of the spherical Hermite expansion coefficients.
        std::vector<size_t> offsets_ecoeffs_;
        /// Offsets of the shell atomic orbital norms.
        std::vector<size_t> offsets_norms_;
        /// Offsets of the atomic orbital positions in the list of all atomic orbitals.
        std::vector<size_t> offsets_sph_;
        /// Indices of the involved shells in the list of all shells.
        std::vector<size_t> shell_idxs_;
    };

    /// Structure containing contiguous data from shell pairs for calculating integrals. By
    /// default, the Gaussian primitives are screened based on the values of the exponential
    /// pre-factor and contraction coefficients. Screening can be disabled by setting the
    /// threshold to zero.
    struct ShellPairData
    {
        /// Constructor for the shell pair data. Expects that the angular momentum in `shells_a`
        /// equals `la` and the same for `shells_b` and `lb`. Cannot be called in an OMP parallel
        /// region.
        ShellPairData(bool use_symm, int la, int lb, const std::vector<Shell> &shells_a,
                      const std::vector<Shell> &shells_b,
                      double primitives_thrs = 1e-15); // TODO: add this number as a constant somewhere

        /// Flag indicating whether symmetry is used or not.
        bool uses_symm_{};

        /// Threshold for screening the primitive Gaussian pairs.
        double primitives_thrs_{};

        /// Angular momentum in bra-shell.
        int la_{};
        /// Angular momentum in ket-shell.
        int lb_{};
        /// Number of shell pairs after screening.
        size_t n_pairs_{};
        /// Total number of shell pairs.
        size_t n_pairs_total_{};
        /// Number of primitive Gaussian pairs after screening.
        size_t n_ppairs_{};
        /// Total number of primitive Gaussian pairs.
        size_t n_ppairs_total_{};

        /// Exponents of the Gaussian primitives.
        std::vector<double> exps_;
        /// Contraction coefficients including the primitive norms.
        std::vector<double> coeffs_;
        /// Coordinates of shells in the shell pairs.
        std::vector<double> coords_;
        /// Normalization constants of the atomic orbitals in the spherical basis.
        std::vector<double> norms_;

        /// Number of Gaussian primitive pairs for each shell pair.
        std::vector<size_t> nrs_ppairs_;
        /// Offsets of the Gaussians for each shell pair.
        std::vector<size_t> offsets_primitives_;
        /// Offsets of the atomic orbital (spherical) positions in the list of all atomic orbitals.
        std::vector<size_t> offsets_sph_;
        /// Offsets of the atomic orbital (cartesian) positions in the list of all atomic orbitals.
        std::vector<size_t> offsets_cart_;
        /// Offsets of the shell atomic orbital norms.
        std::vector<size_t> offsets_norms_;
        /// Offsets of the spherical Hermite expansion coefficients.
        std::vector<size_t> offsets_ecoeffs_;
        /// Offsets of the 1st derivative spherical Hermite expansion coefficients.
        std::vector<size_t> offsets_ecoeffs_deriv1_;
        /// Offsets of the 2nd derivative spherical Hermite expansion coefficients.
        std::vector<size_t> offsets_ecoeffs_deriv2_;
        /// Indices of atoms involved in the shell pairs.
        std::vector<size_t> atomic_idxs_;
        /// Indices of the involved shells in the list of all shells.
        std::vector<size_t> shell_idxs_;

        /// Returns the angular momentum pair.
        std::pair<int, int> getLPair() const
        {
            return {la_, lb_};
        }

    private:
        /// Counts the numbers of total and screened shell and primitive Gaussian pairs.
        void countPairs(const std::vector<Shell> &shells_a, const std::vector<Shell> &shells_b,
                        size_t &n_pairs, size_t &n_pairs_total, size_t &n_ppairs,
                        size_t &n_ppairs_total) const;
    };


    struct ShellPairDataTmp
    {
    private:
    };
}
