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
    /// Calculates the overlap integrals. OMP parallelized.
    vec2d overlap(const Structure &structure);

    /// Calculates a batch of overlap integrals.
    vec2d overlapKernel(size_t ipair, const ShellPairData &sp_data);

    /// Calculates a batch of overlap integral derivatives.
    std::array<vec2d, 6> overlapD1Kernel(size_t ipair, const ShellPairData &sp_data);

    /// Calculates the kinetic energy integrals. OMP parallelized.
    vec2d kineticEnergy(const Structure &structure);

    /// Calculates a batch of kinetic energy integrals.
    vec2d kineticEnergyKernel(size_t ipair, const ShellPairData &sp_data);

    /// Calculates a batch of kinetic energy integral derivatives.
    std::array<vec2d, 6> kineticEnergyD1Kernel(size_t ipair, const ShellPairData &sp_data);

    /// Calculates nuclear attraction integrals. OMP parallelized.
    vec2d nuclearAttraction(const Structure &structure);

    /// Calculates attenuated nuclear attraction integrals. OMP parallelized.
    vec2d nuclearAttractionErf(const Structure &structure, const std::vector<double> &omegas);

    /// Calculates one-electron Coulomb integrals with given point charges {x, y, z, q}.
    /// OMP parallelized.
    vec2d externalCharges(const std::vector<std::array<double, 4>> &point_charges,
                          const Structure &structure);

    /// Calculates attenuated one-electron Coulomb integrals with point charges {x, y, z, q}.
    /// OMP parallelized.
    vec2d externalChargesErf(const std::vector<std::array<double, 4>> &point_charges,
                             const std::vector<double> &omegas,
                             const Structure &structure);

    /// Calculates a batch of one-electron Coulomb integrals for the given charges, {x, y, z, q}.
    /// The Boys function must be initialized with l = la + lb.
    vec2d externalChargesKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates a batch of attenuated one-electron Coulomb integrals for the given charges
    /// {x, y, z, q}. The Boys function must be initialized with l = la + lb.
    vec2d externalChargesErfKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                   const std::vector<double> &omegas,
                                   const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates a batch of one-electron Coulomb integral derivatives for the given charges,
    /// {x, y, z, q}. The Boys function must be initialized with l = la + lb + 1. Calculates
    /// only the orbital part.
    std::array<vec2d, 6>
    externalChargesD1Kernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                            const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates a batch of one-electron Coulomb integral derivatives for the given charges,
    /// {x, y, z, q}. The Boys function must be initialized with l = la + lb + 1. Calculates
    /// only the operator part.
    std::vector<std::array<vec2d, 3>>
    externalChargesOperatorD1Kernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates a batch of one-electron Coulomb integrals at the given charges, {x, y, z, q}.
    /// The Boys function must be initialized with l = la + lb.
    std::vector<vec2d>
    potentialAtExternalChargesKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                     const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates a batch of attenuated one-electron Coulomb integrals at the given charges,
    /// {x, y, z, q}. The Boys function must be initialized with l = la + lb.
    std::vector<vec2d>
    potentialAtExternalChargesErfKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                        const std::vector<double> &omegas,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates dipole moment integrals in the three Cartesian directions. OMP parallelized.
    std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                      const Structure &structure);

    /// Calculates a batch of dipole moment integrals in the three Cartesian directions.
    std::array<vec2d, 3> dipoleMomentKernel(size_t ipair, const std::array<double, 3> &origin,
                                            const ShellPairData &sp_data);

    /// Calculates the one-electron spin-orbit coupling integrals in three Cartesian directions.
    /// OMP parallelized.
    std::array<vec2d, 3> spinOrbitCoupling1El(const Structure &structure);

    /// Calculates a batch of one-electron spin-orbit coupling integrals in three Cartesian
    /// directions. The Boys function must be initialized with l = la + lb + 1.
    std::array<vec2d, 3>
    spinOrbitCoupling1ElKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                               const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates the linear momentum integrals in three Cartesian directions. OMP parallelized.
    std::array<vec2d, 3> momentum(const Structure &structure);

    /// Calculates a batch of linear momentum integrals in three Cartesian directions.
    std::array<vec2d, 3> momentumKernel(size_t ipair, const ShellPairData &sp_data);

    /// Calculates the angular momentum integrals in three Cartesian directions. OMP parallelized.
    std::array<vec2d, 3> angularMomentum(const std::array<double, 3> &origin,
                                         const Structure &structure);

    /// Calculates a batch of angular momentum integrals in three Cartesian directions.
    std::array<vec2d, 3> angularMomentumKernel(size_t ipair, const std::array<double, 3> &origin,
                                               const ShellPairData &sp_data);

    /// Calculates the momentum-potential-momentum integrals, used for example in X2C. OMP
    /// parallelized. For reference, see eqs. (6) https://doi.org/10.1063/1.4803693. Returned
    /// for the 2D Cartesian directions as p_i V p_j, with i, j = x, y, z.
    arr2d<vec2d, 3, 3> pVpIntegrals(const Structure &structure);

    /// Calculates a batch of momentum-potential-momentum integrals, used, for example, in X2C.
    /// The charges are expected as {x, y, z, q}. The Boys function must be initialized with
    /// l = la + lb + 2. For reference, see eqs. (6) https://doi.org/10.1063/1.4803693. Returned
    /// for the 2D Cartesian directions as p_i V p_j, with i, j = x, y, z.
    arr2d<vec2d, 3, 3> pVpKernel(size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                 const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /// Calculates the diagonal of the ERI2 over the auxiliary basis set. OMP parallelized.
    std::vector<double> eri2Diagonal(const Structure &structure);

    /// Calculates the ERI2 over the auxiliary basis set. OMP parallelized.
    vec2d eri2(const Structure &structure);

    /// Calculates the ERI4 diagonal, (ab|ab). OMP parallelized.
    vec2d eri4Diagonal(const Structure &structure);

    /// Calculates the ERI3 tensor, (ab|P) where a and b are main basis AOs, P is an auxiliary
    /// basis AO. OMP parallelized. OMP parallelized.
    vec3d eri3(const Structure &structure);

    /// Calculates the ERI4 tensor. OMP parallelized.
    vec4d eri4(const Structure &structure);

    /// Function for doing an ERI4 benchmark.
    void eri4Benchmark(const Structure &structure);

    /// Returns the main basis set for an atom.
    BasisAtom basisForAtom(int atomic_nr, const std::string &basis_set);

    /// Returns the auxiliary basis set for an atom.
    BasisAtom basisForAtomAux(int atomic_nr, const std::string &aux_basis_set);

    /// Returns the main basis set for a list of atoms.
    basis_atoms_t basisForAtoms(const std::vector<int> &atomic_nrs, const std::string &basis_set);

    /// Returns the auxiliary basis set for a list of atoms.
    basis_atoms_t basisForAtomsAux(const std::vector<int> &atomic_nrs,
                                   const std::string &aux_basis_set);

    /// Returns the Cartesian to spherical transformation for given angular momentum. Returned as
    /// a list, {mu, mu_, val}, where mu and mu_ are spherical and Cartesian indices, respectively,
    /// and val is the value of the transformation coefficient.
    std::vector<std::tuple<int, int, double>> sphericalTrafo(int l);

    /// Constructs the shells from the given basis and coordinates of the atoms (a.u.).
    std::vector<Shell> constructShells(const basis_atoms_t &basis_atoms,
                                       const std::vector<std::array<double, 3>> &coords_atoms);

    /// Constructs the shells from the given basis on the given atom.
    std::vector<Shell> constructShells(int atomic_nr, const basis_shells_t &basis_shells,
                                       const std::array<double, 3> &coords_atom);

    /// Calculates the normalization coefficients of the atomic orbitals in a shell.
    std::vector<double> calcShellNorms(int l, const std::vector<double> &coeffs,
                                       const std::vector<double> &exps,
                                       const std::vector<double> &primitive_norms);

    /// Constructs the shell for the auxiliary basis set up to `l_max_aux` (included).
    std::vector<ShellData> shellDataAux(const Structure &structure);

    /// Constructs the shell pair data for the main basis set up to `l_max` (included)
    /// Returns data for (la, lb) pairs. If symmetry is used, returns data for (la >= lb).
    std::vector<ShellPairData> shellPairData(bool use_symm, const Structure &structure);

    ///  Constructs the shell pair data from the given shells.
    std::vector<ShellPairData> shellPairData(const std::vector<Shell> &shells_a,
                                             const std::vector<Shell> &shells_b);

    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// in one Cartesian direction. Based on eq. (8) from https://doi.org/10.1063/5.0217001.
    /// Returned as E^i_t -> (i, t).
    vec2d ecoeffsRecurrence1(double one_o_2a, int l);

    /// Calculates the Hermite expansion coefficients for a primitive Gaussian function product in
    /// one Cartesian direction. Based on eqs. (9.5.6) and (9.5.7) from
    /// https://doi.org/10.1002/9781119019572. Returned as E^{ij}_t -> (i, j, t).
    vec3d ecoeffsRecurrence2(double a, double b, int la, int lb, double PA, double PB,
                             double Kab);

    /// Calculates the first derivative of the Hermite expansion coefficients. Based on eq. (20)
    /// from https://doi.org/10.1007/BF01132826.
    vec3d ecoeffsRecurrence2_n1(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0);

    /// Calculates the second derivative of the Hermite expansion coefficients. Based on eq. (23)
    /// from https://doi.org/10.1007/BF01132826
    vec3d ecoeffsRecurrence2_n2(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0, const vec3d &ecoeffs1);

    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// in the three Cartesian directions.
    std::array<vec2d, 3> ecoeffsPrimitive(double a, int l);

    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    /// Calculates the first derivative of the Hermite expansion coefficients for a single
    /// primitive Gaussian function product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0);

    /// Calculates the second derivative of the Hermite expansion coefficients for a single
    /// primitive Gaussian function product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair_n2(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0,
                                                 const std::array<vec3d, 3> &ecoeffs1);

    /// Calculates the Hermite expansion coefficients for a single shell with given primitive
    /// Gaussian exponents. Returns E^{ii'}_0 E^{jj'}_0 E^{kk'}_0 for each primitive pair.
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    /// Calculates the Hermite expansion coefficients, E_{munu, tuv}, for each primitive Gaussian
    /// pair in each shell pair, where \mu \nu refer to the AO indices. The coefficients have AO
    /// norms and primitive contraction coefficients multiplied into.
    std::vector<double> ecoeffsSHARK(const ShellPairData &sp_data, bool transpose = false);

    /// Calculates the Hermite expansion coefficients, E_{mu, tuv}, for each primitive Gaussian
    /// pair in each shell pair. The coefficients have AO norms and primitive contraction
    /// coefficients multiplied into.
    std::vector<double> ecoeffsSHARK(const ShellData &sh_data, bool transpose = false);

    /// Calculates the first derivate Hermite expansion coefficients,
    /// The coefficients have AO norms and primitive contraction coefficients multiplied into.
    std::vector<double> ecoeffsD1SHARK(const ShellPairData &sp_data, bool transpose = false);

    class BoysGrid; // Forward declaration

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
