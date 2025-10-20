#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/structure.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <array>
#include <cassert>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef LIBLE_MAIN_BASIS_DIR
#define path_to_basis_sets LIBLE_MAIN_BASIS_DIR
#endif

#ifdef LIBLE_AUX_BASIS_DIR
#define path_to_aux_basis_sets LIBLE_AUX_BASIS_DIR
#endif

namespace lible::ints
{
    /**
     * \defgroup IntsMainInterface
     * Groups the contents from the header file `<lible/ints/ints.hpp>` in the `lible::ints`
     * namespace.
     */

    /**
     * \ingroup IntsMainInterface
     * Computes the entire overlap integral matrix for a given molecular structure.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \return Normalized spherical-basis overlap integrals.
     */
    vec2d overlap(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the overlap integrals for a pair of two shells.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis overlap integrals for the given shell pair.
     */
    vec2d overlapKernel(int ipair, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the first derivative of the overlap integrals for a pair of two shells.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return 6D-array of normalized spherical-basis overlap integral derivatives. The integral
     * derivatives are given for \f$(A_x, A_y, A_z, B_x, B_y, B_z)\f$.
     */
    std::array<vec2d, 6> overlapD1Kernel(int ipair, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the entire kinetic energy integral matrix for a given molecular structure.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \return Normalized spherical-basis kinetic energy integrals.
     */
    vec2d kineticEnergy(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the kinetic energy integrals for a pair of two shells.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis kinetic energy integrals for the given shell pair.
     */
    vec2d kineticEnergyKernel(int ipair, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the first derivative of the kinetic energy integrals for a pair of two shells.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return 6D-array of normalized spherical-basis kinetic energy integral derivatives. The
     * integral derivatives are given for \f$(A_x, A_y, A_z, B_x, B_y, B_z)\f$.
     */
    std::array<vec2d, 6> kineticEnergyD1Kernel(int ipair, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the entire nuclear Coulomb attraction integral matrix for a given molecular
     * structure.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \return Normalized spherical-basis nuclear attraction Coulomb integrals.
     */
    vec2d nuclearAttraction(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the entire nuclear Coulomb attraction integral matrix for a given molecular
     * structure with the nuclear charges represented by the error function.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \param omegas list of width parameters for the Gaussian in the error function. The length
     * of this list has to equal the number of atoms in `structure`.
     * \return Normalized spherical-basis erf-attenuated nuclear attraction Coulomb integrals.
     */
    vec2d nuclearAttractionErf(const Structure &structure, const std::vector<double> &omegas);

    /**
     * \ingroup IntsMainInterface
     * Computes the entire integral matrix for the Coulomb interaction between electrons and the
     * given point charges.
     *
     * @param point_charges list of point charges, \f$(x, y, z, q)\f$, with their coordinates in
     * Bohr.
     * @param structure `Structure` object representing the molecular geometry and basis sets.
     * @return Normalized spherical-basis Coulomb interaction integrals.
     */
    vec2d externalCharges(const std::vector<std::array<double, 4>> &point_charges,
                          const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the entire nuclear Coulomb interaction integral matrix for a given molecular
     * structure with the point charges represented with by the error function.
     *
     * \param point_charges list of point charges given by \f$(x, y, z, q)\f$.
     * \param omegas list of width parameters for the Gaussian in the error function. The length
     * of this list has to equal the number of point charges.
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \return Normalized spherical-basis erf-attenuated nuclear attraction Coulomb integrals.
     */
    vec2d externalChargesErf(const std::vector<std::array<double, 4>> &point_charges,
                             const std::vector<double> &omegas,
                             const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the Coulomb interaction integrals for a pair of shells.
     *
     * \param ipair Index of the shell pair.
     * \param charges list of point charges given by \f$(x, y, z, q)\f$.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l_{ab} = l_a + l_b\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis Coulomb interaction integrals for the given shell pair.
     */
    vec2d externalChargesKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the Coulomb interaction integrals for a pair of shells. The external charges
     * are represented by the error function.
     *
     * \param ipair Index of the shell pair.
     * \param charges List of point charges given by \f$(x, y, z, q)\f$.
     * \param omegas List of width parameters for the Gaussian in the error function. The length
     * of this list has to equal the length of `charges`.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l_{ab} = l_a + l_b\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis Coulomb interaction integrals for the given shell pair.
     */
    vec2d externalChargesErfKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                   const std::vector<double> &omegas,
                                   const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the integrals for the first derivative of the external charges for a pair of two
     * shells. Calculates the AO part of the derivatives,
     * \f[
     *   \left\{(\nabla a| \hat{g}(r) | b), (a| \hat{g}(r) |\nabla b)\right\}
     *   \; \text{with the operator} \;
     *   \hat{g}(r) = -\sum_i \frac{q_i}{|\mathbf{r} - \mathbf{r}_i|}
     * \f]
     * The operator part is calculated by `externalChargesOperatorD1Kernel`.
     *
     * \param ipair Index of the shell pair.
     * @param charges list of point charges, \f$(x, y, z, q)\f$, with their coordinates in
     * Bohr.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l = l_a + l_b + 1\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return 6D-array of normalized spherical-basis external integral derivatives. The integral
     * derivatives correspond to \f$(A_x, A_y, A_z, B_x, B_y, B_z)\f$.
     */
    std::array<vec2d, 6>
    externalChargesD1Kernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                            const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the integrals for the first derivative of the Coulomb operator for given external
     * charges,
     * \f[
     *   \{(a| \nabla_i \hat{g}_i(r) |b)\}
     *   \; \text{with the operator} \;
     *   \hat{g}(r) = -\sum_i \frac{q_i}{|\mathbf{r} - \mathbf{r}_i|}
     * \f]
     *
     * \param ipair Index of the shell pair.
     * @param charges list of point charges, \f$(x, y, z, q)\f$, with their coordinates in
     * Bohr.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l = l_a + l_b + 1\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return A vector of 3D-arrays corresponding to the operator derivative integrals. The
     * returned list is of the length of `charges`.
     * The integrals are returned normalized and in the spherical basis.
     */
    std::vector<std::array<vec2d, 3>>
    externalChargesOperatorD1Kernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of integrals for the Coulomb interaction with electrons at each charge.
     * \f[
     *   \{(a| \frac{-q_i}{|\mathbf{r} - \mathbf{r}_i|} |b)\}
     * \f]
     *
     * \param ipair Index of the shell pair.
     * \param charges List of point charges, \f$(x, y, z, q)\f$, with their coordinates in
     * Bohr.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l_{ab} = l_a + l_b\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return A list of normalized spherical-basis integrals with the length of `charges`.
     */
    std::vector<vec2d>
    potentialAtExternalChargesKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                     const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of integrals for the Coulomb interaction with electrons at each charge.
     * \f[
     *   \{(a| \frac{-q_i \cdot \text{erf}(\omega |\mathbf{r} - \mathbf{r}_i|)}
     *   {|\mathbf{r} - \mathbf{r}_i|} |b)\}
     * \f]
     *
     * \param ipair Index of the shell pair.
     * \param charges List of point charges, \f$\{(x, y, z, q)\}\f$, with their coordinates in
     * Bohr.
     * \param omegas List of width parameters for the Gaussian in the error function. The length
     * of this list has to equal the length of `charges`.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l_{ab} = l_a + l_b\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return A list of normalized spherical-basis integrals with the length of `charges`.
     */
    std::vector<vec2d>
    potentialAtExternalChargesErfKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                        const std::vector<double> &omegas,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes all the dipole moment integral matrices for the three Cartesian directions.
     *
     * \param origin Origin of the Cartesian dipole moments, \f$(O_x, O_y, O_z)\f$. Given in
     * atomic units (Bohr).
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \return Normalized spherical-basis dipole moment integrals for each Cartesian direction,
     * \f$(x, y, z)\f$.
     */
    std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                      const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of dipole moment integrals for the three Cartesian directions.
     *
     * \param ipair Index of the shell pair.
     * \param origin Origin of the Cartesian dipole moments, \f$(O_x, O_y, O_z)\f$. Given in
     * atomic units (Bohr).
     * \param sp_data `ShellPairData` object containing all the information required for calculating
     * integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis dipole moment integrals for each Cartesian direction,
     * \f$(x, y, z)\f$.
     */
    std::array<vec2d, 3> dipoleMomentKernel(int ipair, const std::array<double, 3> &origin,
                                            const ShellPairData &sp_data);

    // /**  TODO: dox. */
    std::array<vec2d, 3> spinOrbitCoupling1El(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of spin-orbit coupling integrals
     * \f[
     *   (a| \frac{\mathbf{r}_C}{r^3_C} \nabla_B |b) = \nabla_A \times \nabla_B (a|r^{-1}_C|b)
     * \f]
     *
     * \param ipair Index of the shell pair.
     * \param charges list of point charges given by \f$(x, y, z, q)\f$.
     * \param boys_grid Pre-initialized grid for calculating the Boys function. Must be initialized
     * with \f$l = l_a + l_b + 1\f$ where the \f$(l_a, l_b)\f$ corresponds to the given `sp_data`.
     * \param sp_data `ShellPairData` object containing all the information required for
     * calculating integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis dipole moment integrals for each Cartesian direction.
     */
    std::array<vec2d, 3>
    spinOrbitCoupling1ElKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                               const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of linear momentum integrals, \f$\{-(a| \mathbf{\nabla} |b)\}\f$. The
     * imaginary unit is omitted.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for
     * calculating integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis momentum integrals for each Cartesian direction.
     */
    std::array<vec2d, 3> momentumKernel(int ipair, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes a batch of linear momentum integrals,
     * \f$\{-(a| \mathbf{r} \times \mathbf{\nabla} |b)\}\f$. The imaginary unit is omitted.
     *
     * \param ipair Index of the shell pair.
     * \param sp_data `ShellPairData` object containing all the information required for
     * calculating integrals for the given \f$(l_a, l_b)\f$.
     * \return Normalized spherical-basis angular momentum integrals for each Cartesian direction.
     */
    std::array<vec2d, 3> angularMomentumKernel(int ipair, const ShellPairData &sp_data);

    /** */
    arr2d<vec2d, 3, 3> pVpKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                 const BoysGrid &boys_grid, const ShellPairData &sp_data);

    /**
     * \ingroup IntsMainInterface
     * Computes the diagonal of the two-center electron repulsion (ERI2) integrals over auxiliary
     * basis functions, \f$\{(P|P)\}\f$.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \note Uses OMP parallelization.
     * \return List of normalized spherical-basis integrals.
     */
    std::vector<double> eri2Diagonal(const Structure &structure);


    /**
     * \ingroup IntsMainInterface
     * Computes the two-center electron repulsion integrals (ERI2) over auxiliary basis functions,
     * \f$\{(P|Q)\}\f$.
     *
     * \param structure `Structure` object representing the molecular geometry and basis sets.
     * \note Uses OMP parallelization.
     * \return 2D array of normalized spherical-basis integrals.
     */
    vec2d eri2(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the diagonal of the four-center electron repulsion integrals (ERI4),
     * \f$(\mu\mu|\nu\nu)\f$.
     *
     * @param structure `Structure` object representing the molecular geometry and basis sets.
     * @note Uses OMP parallelization.
     * @return 2D array of normalized spherical-basis integrals.
     */
    vec2d eri4Diagonal(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the three-center electron repulsion integrals (ERI3), \f$(\mu\nu|P)\f$.
     *
     * @param structure `Structure` object representing the molecular geometry and basis sets.
     * @note Uses OMP parallelization.
     * @return 3D array of normalized spherical-basis integrals.
     */
    vec3d eri3(const Structure &structure);

    /**
     * \ingroup IntsMainInterface
     * Computes the four-center electron repulsion integrals (ERI4) over main basis functions,
     * \f$(\mu\nu|\kappa\tau)\f$.
     *
     * @param structure `Structure` object representing the molecular geometry and basis sets.
     * @note Uses OMP parallelization.
     * @return 4D array of normalized spherical-basis integrals
     */
    vec4d eri4(const Structure &structure);

    // /**
    //      *
    //      */
    void eri4Benchmark(const Structure &structure);

    // /**
    //      * \ingroup ints
    //      * Typedef for bundling Gaussian primitive exponents and contraction coefficients
    //      * of a shell.
    //      */
    /** */
    using shell_exps_coeffs_t = std::pair<std::vector<double>, std::vector<double>>;

    /** */
    using basis_atom_t = std::map<int, std::vector<shell_exps_coeffs_t>>;

    /** */
    using basis_atoms_t = std::map<int, basis_atom_t>;

    // /**
    //      * \ingroup ints
    //      * Returns the basis set for the given atom. The exponents and contraction coefficients
    //      * are listed for every angular momentum.
    //      */
    // std::map<int, std::vector<shell_exps_coeffs_t>>
    basis_atom_t basisForAtom(int atomic_nr, const std::string &basis_set);

    // /**
    //      * \ingroup ints
    //      * Returns the auxiliary basis set the given atom. The exponents and contraction
    //      * coefficients are listed for every angular momentum.
    //      */
    // std::map<int, std::vector<shell_exps_coeffs_t>>

    basis_atom_t basisForAtomAux(int atomic_nr, const std::string &aux_basis_set);

    // /**
    //      * \ingroup ints
    //      * Returns the basis set for the given atoms. The exponents and contraction coefficients
    //      * are listed for every angular momentum per atom.
    //      */
    // std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
    basis_atoms_t basisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);

    // /**
    //      * \ingroup ints
    //      * Returns the auxiliary basis set for the given atoms. The exponents and contraction
    //      * coefficients are listed for every angular momentum per atom.
    //      */
    // std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
    basis_atoms_t basisForAtomsAux(const std::set<int> &atomic_nrs, const std::string &aux_basis_set);

    // /**
    //      * \ingroup ints
    //      * Returns the Cartesial to spherical basis transformation,
    //      * \f$\{(\text{i_spherical}, \text{i_cartesian}, \text{val})\}\f$.
    //      */
    std::vector<std::tuple<int, int, double>> sphericalTrafo(int l);

    // /** */
    ERI4Kernel deployERI4Kernel(const ShellPairData &sp_data_ab,
                                const ShellPairData &sp_data_cd);

    // /** */
    ERI3Kernel deployERI3Kernel(const ShellPairData &sp_data_ab,
                                const ShellData &sp_data_cd);

    // /** */
    ERI2Kernel deployERI2Kernel(const ShellData &sp_data_a,
                                const ShellData &sp_data_b);

    // /** */
    ERI4D1Kernel deployERI4D1Kernel(const ShellPairData &sp_data_ab,
                                    const ShellPairData &sp_data_cd);

    // /** */
    ERI3D1Kernel deployERI3D1Kernel(const ShellPairData &sp_data_ab,
                                    const ShellData &sp_data_b);


    // /** */
    ERI2D1Kernel deployERI2D1Kernel(const ShellData &sp_data_a,
                                    const ShellData &sp_data_b);

    // /** */
    ERI2D2Kernel deployERI2D2Kernel(const ShellData &sh_data_a,
                                    const ShellData &sh_data_b);

    // /** */
    ERI4SOCKernel deployERI4SOCKernel(const ShellPairData &sp_data_ab,
                                      const ShellPairData &sp_data_cd);

    // /** */
    ERI3SOCKernel deployERI3SOCKernel(const ShellPairData &sp_data_ab,
                                      const ShellData &sh_data_c);

    // /**
    //      *
    //      */
    std::vector<Shell> constructShells(const basis_atoms_t &basis_atoms,
                                       const std::vector<int> &atomic_nrs,
                                       const std::vector<std::array<double, 3>> &coords_atoms);

    // /**
    //      *
    //      */
    std::vector<double> calcShellNorms(int l, const std::vector<double> &coeffs,
                                       const std::vector<double> &exps,
                                       const std::vector<double> &primitive_norms);

    // /**
    //      * \ingroup ints
    //      * Constructs the shell datas for the auxiliary basis set, up to l_max.
    //      */
    std::vector<ShellData> shellDataAux(const Structure &structure);

    // /**
    //      *
    //      */
    std::vector<ShellPairData> shellPairData(bool use_symm, const Structure &structure);

    // /**
    //      *
    //      */
    std::vector<ShellPairData> shellPairData(const std::vector<Shell> &shells_a,
                                             const std::vector<Shell> &shells_b);

    // /** */
    vec2d ecoeffsRecurrence1(double one_o_2a, int l);

    // /** */
    vec3d ecoeffsRecurrence2(double a, double b, int la, int lb, double PA, double PB, double Kab);

    // /** */
    vec3d ecoeffsRecurrence2_n1(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0);

    // /** */
    std::array<vec2d, 3> ecoeffsPrimitive(double a, int l);

    // /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    // /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0);

    // /** */
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    // /** TODO: */
    class BoysGrid;

    // /** TODO: */
    std::vector<double> calcBoysF(int max_n, double x, const BoysGrid &boys_grid);

    // /** TODO: */
    vec3d calcRInts3D(int l, double p, const double *xyz_ab, const double *fnx);

    /**
     * \ingroup IntsMainInterface
     * Computes the norm of the spherical Gaussian primitive function from
     * \f[
     *   N = \left(\frac{(2a/\pi)^{3/2}(4a)^l}{(2l - 1)!!}\right)^{1/2}
     * \f]
     *
     * \param l angular momentum
     * \param exp Gaussian exponent
     * \return Norm of the spherical Gaussian primitive
     */
    double purePrimitiveNorm(int l, double exp);

    /**
     * \ingroup IntsMainInterface
     * Computes the number of Cartesian Gaussians for the given angular momentum using
     * \f[
     *   N = (l + 1)(l + 2)/2
     * \f]
     *
     * \param l angular momentum
     * \return Number of Cartesian Gaussians.
     */
    constexpr int numCartesians(int l);

    /**
     * \ingroup IntsMainInterface
     * Computes the number of spherical Gaussians for the given angular momentum using
     * \f[
     *   N = 2l + 1
     * \f]
     *
     * \param l angular momentum
     * \return Number of spherical Gaussians.
     */
    constexpr int numSphericals(int l);

    /**
     * \ingroup IntsMainInterface
     * Computes the total number of Hermite Gaussians up the given angular momentum
     * using
     * \f[
     *   N = (l + 1)(l + 2)(l + 3) / 6
     * \f]
     *
     * \param l angular momentum
     * \return Number of Hermite Gaussians up to total angular momentum `l`.
     */
    constexpr int numHermites(int l);

    /**
     * \ingroup IntsMainInterface
     * Computes the list of Cartesian Gaussian exponents for the given angular momentum.
     * The length of this list is given by \f$(l + 1)(l + 2) / 2\f$. The Cartesian gaussians
     * follow the so-called alphabetic ordering.
     *
     * \param l angular momentum
     * \return List of Cartesian Gaussian exponents \f$\{(i, j, k)\}\f$.
     */
    std::vector<std::array<int, 3>> cartExps(int l);

    /**
     * \ingroup IntsMainInterface
     *  Computes a list of angular momentum pairs up to the given angular momentum such that
     *  \f$l_a >= l_b\f$:
     * \f[
     *   \{(0, 0), (1, 0), (1, 1), \ldots, (l, l) \}
     * \f]
     *
     * \param l angular momentum
     * \return List of angular momentum pairs.
     */
    std::vector<std::pair<int, int>> getLPairsSymm(int l);

    /**
     * \ingroup IntsMainInterface
     *  Computes a list of all angular momentum pairs up to the given angular momentum:
     * \f[
     *   \{(0, 0), (0, 1), (1, 0), (1, 1), \ldots, (l, l) \}
     * \f]
     *
     * \param l angular momentum
     * \return List of angular momentum pairs.
     */
    std::vector<std::pair<int, int>> getLPairsNoSymm(int l);

    /**
     * \ingroup IntsMainInterface
     * \return Names (lower case) of all the available main basis sets in the library.
     */
    std::set<std::string> availableBasisSets();

    /**
     * \ingroup IntsMainInterface
     * \return Names (lower case) of all the available auxiliary basis sets in the library.
     */
    std::set<std::string> availableBasisSetsAux();

    /** TODO: */
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
        /** Absolute path to the main basis sets. */
        static inline std::string main_basis_sets_path{path_to_basis_sets};

        /** Absolute path to the auxiliary basis sets. */
        static inline std::string aux_basis_sets_path{path_to_aux_basis_sets};
    };
}
