#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/structure.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <array>
#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

#ifdef LIBLE_MAIN_BASIS_DIR
#define path_to_basis_sets LIBLE_MAIN_BASIS_DIR
#endif

#ifdef LIBLE_AUX_BASIS_DIR
#define path_to_aux_basis_sets LIBLE_AUX_BASIS_DIR
#endif

namespace lible::ints
{
    vec2d overlap(const Structure &structure);

    vec2d overlapKernel(size_t ipair, const ShellPairData &sp_data);

    std::array<vec2d, 6> overlapD1Kernel(size_t ipair, const ShellPairData &sp_data);

    vec2d kineticEnergy(const Structure &structure);

    vec2d kineticEnergyKernel(size_t ipair, const ShellPairData &sp_data);

    std::array<vec2d, 6> kineticEnergyD1Kernel(size_t ipair, const ShellPairData &sp_data);

    vec2d nuclearAttraction(const Structure &structure);

    vec2d nuclearAttractionErf(const Structure &structure, const std::vector<double> &omegas);

    vec2d externalCharges(const std::vector<std::array<double, 4>> &point_charges,
                          const Structure &structure);

    vec2d externalChargesErf(const std::vector<std::array<double, 4>> &point_charges,
                             const std::vector<double> &omegas,
                             const Structure &structure);

    vec2d externalChargesKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                const BoysGrid &boys_grid, const ShellPairData &sp_data);

    vec2d externalChargesErfKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                   const std::vector<double> &omegas,
                                   const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::array<vec2d, 6>
    externalChargesD1Kernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                            const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::vector<std::array<vec2d, 3>>
    externalChargesOperatorD1Kernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::vector<vec2d>
    potentialAtExternalChargesKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                     const BoysGrid &boys_grid, const ShellPairData &sp_data);
    std::vector<vec2d>
    potentialAtExternalChargesErfKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                        const std::vector<double> &omegas,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                      const Structure &structure);

    std::array<vec2d, 3> dipoleMomentKernel(size_t ipair, const std::array<double, 3> &origin,
                                            const ShellPairData &sp_data);

    std::array<vec2d, 3> spinOrbitCoupling1El(const Structure &structure);

    std::array<vec2d, 3>
    spinOrbitCoupling1ElKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                               const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::array<vec2d, 3> momentumKernel(size_t ipair, const ShellPairData &sp_data);

    std::array<vec2d, 3> angularMomentumKernel(size_t ipair, const ShellPairData &sp_data);

    arr2d<vec2d, 3, 3> pVpIntegrals(const Structure &structure);

    arr2d<vec2d, 3, 3> pVpKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                 const BoysGrid &boys_grid, const ShellPairData &sp_data);

    std::vector<double> eri2Diagonal(const Structure &structure);

    vec2d eri2(const Structure &structure);

    vec2d eri4Diagonal(const Structure &structure);

    vec3d eri3(const Structure &structure);

    vec4d eri4(const Structure &structure);

    ///
    void eri4Benchmark(const Structure &structure);

    /// Returns the main basis set for an atom.
    basis_atom_t basisForAtom(int atomic_nr, const std::string &basis_set);

    /// Returns the auxiliary basis set for an atom.
    basis_atom_t basisForAtomAux(int atomic_nr, const std::string &aux_basis_set);

    /// Returns the main basis set for a list of atoms.
    basis_atoms_t basisForAtoms(const std::vector<int> &atomic_nrs, const std::string &basis_set);

    /// Returns the auxiliary basis set for a list of atoms.
    basis_atoms_t basisForAtomsAux(const std::vector<int> &atomic_nrs,
                                   const std::string &aux_basis_set);
    ///
    std::vector<std::tuple<int, int, double>> sphericalTrafo(int l);

    /// Constructs the shells from the given basis and coordinates of the atoms (a.u.).
    std::vector<Shell> constructShells(const basis_atoms_t &basis_atoms,
                                       const std::vector<std::array<double, 3>> &coords_atoms);

    /// Calculates the normalization coefficients of the atomic orbitals in a shell.
    std::vector<double> calcShellNorms(int l, const std::vector<double> &coeffs,
                                       const std::vector<double> &exps,
                                       const std::vector<double> &primitive_norms);

    /// Constructs the shell for the auxiliary basis set up to `l_max_aux` (included).
    std::vector<ShellData> shellDataAux(const Structure &structure);

    /// Constructs the shell pair data for the main basis set up to `l_max` (included)
    /// Returns data for (la, lb) pairs. If symmetry is used, returns data for (la >= lb).
    std::vector<ShellPairData> shellPairData(bool use_symm, const Structure &structure);

    ///  Constructs the shell pair data from the given shells. Assumes the angular momentum
    ///  within `shells_a` is the same as in `shells_b`.
    std::vector<ShellPairData> shellPairData(const std::vector<Shell> &shells_a,
                                             const std::vector<Shell> &shells_b);

    ///
    vec2d ecoeffsRecurrence1(double one_o_2a, int l);

    ///
    vec3d ecoeffsRecurrence2(double a, double b, int la, int lb, double PA, double PB, double Kab);

    ///
    vec3d ecoeffsRecurrence2_n1(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0);

    ///
    std::array<vec2d, 3> ecoeffsPrimitive(double a, int l);

    ///
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    ///
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0);

    ///
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    ///
    class BoysGrid;

    /// Calculates the boys function for the given `n`, `x` and `boys_grid`. The returned vector
    /// has the length `n + 1`. At `x = 0`, the function uses eq. (9.8.6) from
    /// https://doi.org/10.1007/s10008-001-0256-1. At large x values of above 30, eq. (9.8.13)
    /// gets used, omitting the exponential terms. At values between 0 and 30, first, eq. (9.8.12)
    /// is used with seven evaluations and then, the downward recursion from 9.8.14 is used.
    std::vector<double> calcBoysF(int n, double x, const BoysGrid &boys_grid);

    // /** TODO: */
    vec3d calcRInts3D(int l, double p, const double *xyz_ab, const double *fnx);

    ///
    double purePrimitiveNorm(int l, double exp);

    ///
    std::vector<std::array<int, 3>> cartExps(int l);

    ///
    std::vector<std::pair<int, int>> getLPairsSymm(int l);

    ///
    std::vector<std::pair<int, int>> getLPairsNoSymm(int l);

    /// Returns the available main basis sets.
    std::set<std::string> availableBasisSets();

    /// Returns the available auxiliary basis sets.
    std::set<std::string> availableBasisSetsAux();

    ///
    class BasisPaths
    {
    public:
        static std::string getMainBasisSetsPath()
        {
            return main_basis_sets_path;
        }

        static std::string getAuxBasisSetsPath()
        {
            return aux_basis_sets_path;
        }

        static void setMainBasisSetsPath(const std::string &path)
        {
            main_basis_sets_path = path;
        }

        static void setAuxBasisSetsPath(const std::string &path)
        {
            aux_basis_sets_path = path;
        }

    private:
        /// Absolute path to the main basis sets.
        static inline std::string main_basis_sets_path{path_to_basis_sets};

        /// Absolute path to the auxiliary basis sets.
        static inline std::string aux_basis_sets_path{path_to_aux_basis_sets};
    };
}
