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
     * Calculates the Coulomb interaction integrals for a pair of shells.
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
     * Calculates the Coulomb interaction integrals for a pair of shells. The external charges
     * are represented by the error function.
     *
     * \param ipair Index of the shell pair.
     * \param charges list of point charges given by \f$(x, y, z, q)\f$.
     * \param omegas list of width parameters for the Gaussian in the error function. The length
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

    std::array<vec2d, 6>
    externalChargesD1Kernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                            const BoysGrid &boys_grid, const ShellPairData &sp_data);

    // /**
    //      * Calculates a batch of normalized Coulombic operator derivative integrals for the shell
    //      * pair 'ipair'. The derivatives are given for each charge as (Ax, Ay, Az). In spherical basis.
    //      * The charges should be given as a list  {(x, y, z, charge)}, with xyz-coordinates in
    //      * atomic units. The Boys grid should be initialized for lab = la + lb in the given shell
    //      * pair data.
    //      */
    std::vector<std::array<vec2d, 3>>
    externalChargesOperatorD1Kernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);

    // /**
    //      * Calculates a batch of normalized Coulombic operator integrals for the shell pair 'ipair'.
    //      * For every charge a batch of integrals is calculated.
    //      * The charges should be given as a list  {(x, y, z, charge)}, with xyz-coordinates in
    //      * atomic units. The Boys grid should be initialized for lab = la + lb in the given shell
    //      * pair data.
    //      */
    std::vector<vec2d>
    potentialAtExternalChargesKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                     const BoysGrid &boys_grid, const ShellPairData &sp_data);

    // /**
    //      * Calculates a batch of normalizederf-attenuated  Coulombic operator integrals for the shell pair 'ipair'.
    //      * For every charge a batch of integrals is calculated.
    //      * The charges should be given as a list  {(x, y, z, charge)}, with xyz-coordinates in
    //      * atomic units. The Boys grid should be initialized for lab = la + lb in the given shell
    //      * pair data. The screening factor omega should be given as a std::vector<double>
    //      */
    std::vector<vec2d>
    potentialAtExternalChargesErfKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                                        const std::vector<double> &omegas,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data);

    // /**
    //      * \ingroup ints
    //      * Calculates the dipole moment integral matrices for the \f$x,y,z\f$-directions.
    //      */
    std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                      const Structure &structure);

    // /**
    //      * Calculates a batch of normalized dipole moment integrals for the shell pair 'ipair'.
    //      * In spherical basis. The integrals are given as (x, y, z). The origin is expected in
    //      * atomic units (bohr).
    //      */
    std::array<vec2d, 3> dipoleMomentKernel(int ipair, const std::array<double, 3> &origin,
                                            const ShellPairData &sp_data);

    // /**  TODO: dox. */
    std::array<vec2d, 3> spinOrbitCoupling1El(const Structure &structure);

    // /**
    //      * TODO: mention that it requires lab + 1 angular momentum boys_grid.
    //      */
    std::array<vec2d, 3>
    spinOrbitCoupling1ElKernel(int ipair, const std::vector<std::array<double, 4>> &charges,
                               const BoysGrid &boys_grid, const ShellPairData &sp_data);

    // /**
    //      * \ingroup ints
    //      * Calculates the diagonal of the two-center ERI matrix over auxiliary basis functions,
    //      * \f$(P|P)\f$.
    //      */
    std::vector<double> eri2Diagonal(const Structure &structure);

    // /**
    //      * \ingroup ints
    //      * Calculates the two-center ERI matrix over the auxiliary basis functions, \f$(P|Q)\f$.
    //      */
    vec2d eri2(const Structure &structure);

    // /**
    //      * \ingroup ints
    //      * Calculates the diagonal of the four-center ERI matrix, \f$(\mu\nu|\mu\nu)\f$.
    //      */
    vec2d eri4Diagonal(const Structure &structure);

    // /**
    //      * \ingroup ints
    //      * Calculates the three-center ERIs, \f$(\mu\nu|P)\f$.
    //      */
    vec3d eri3(const Structure &structure);

    // /**
    //      * \ingroup ints
    //      * Calculates the four-center ERIs, \f$(\mu\nu|\kappa\tau)\f$.
    //      */
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
                                const vec3d &ecoeffs);

    // /** */
    std::array<vec2d, 3> ecoeffsPrimitive(double a, int l);

    // /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    // /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs);

    // /** */
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    // /** TODO: */
    class BoysGrid;

    // /** TODO: */
    std::vector<double> calcBoysF(int max_n, double x, const BoysGrid &boys_grid);

    // /** TODO: */
    vec3d calcRInts3D(int l, double p, const double *xyz_ab, const double *fnx);

    // /**
    //      * \ingroup ints
    //      * TODO: write dox.
    //      */
    double purePrimitiveNorm(int l, double exp);

    // /**
    //      * \ingroup ints
    //      * Returns the number of Cartesian Gaussians.
    //      */
    constexpr int numCartesians(int l);

    // /**
    //      * \ingroup ints
    //      * Returns the number of spherical Gaussians.
    //      */
    constexpr int numSphericals(int l);

    // /**
    //      * \ingroup ints
    //      * Returns the number of Hermite Gaussians.
    //      */
    constexpr int numHermites(int l);

    // /**
    //      * \ingroup ints
    //      * Returns the exponents of a Cartesian Gaussian \f$(x,y,z)\f$-directions for the given
    //      * angular momentum.
    //      */
    std::vector<std::array<int, 3>> cartExps(int l);

    // /**
    //      * \ingroup ints
    //      * Returns a list of angular momentum pairs such that la >= lb:
    //      *   {(0, 0), (1, 0), (1, 1), ..., (l_max, l_max)}.
    //      */
    std::vector<std::pair<int, int>> getLPairsSymm(int l_max);

    // /**
    //      * \ingroup ints
    //      * Returns a list of angular momentum pairs: {(0, 0), (1, 0), (0, 1), ..., (l_max, l_max)}.
    //      */
    std::vector<std::pair<int, int>> getLPairsNoSymm(int l_max);

    // /**
    //      * \ingroup ints
    //      * Returns the names of all available basis sets in lower case.
    //      */
    std::set<std::string> availableBasisSets();

    // /**
    //      * \ingroup ints
    //      * Returns the names of all available auxiliary basis sets in lower case.
    //      */
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
