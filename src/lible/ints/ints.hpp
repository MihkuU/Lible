#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/structure.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <array>
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         * \defgroup ints
         */

        /**
         * \ingroup ints
         * Calculates the diagonal of the two-center ERI matrix over auxiliary basis functions,
         * \f$(P|P)\f$.
         */
        std::vector<double> eri2Diagonal(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the overlap integral matrix.
         */
        vec2d overlap(const Structure &structure);

        /**
         * Calculates a batch of normalized overlap integrals for the shell pair 'ipair'.
         * In spherical basis.
         */
        vec2d overlapKernel(const int ipair, const ShellPairData &sp_data);

        /**
         * Calculates a batch of normalized overlap integral derivatives for the shell pair
         * 'ipair'. The derivatives are given as (Ax, Ay, Az, Bx, By, Bz). In spherical basis.
         */
        std::array<vec2d, 6> overlapD1Kernel(const int ipair, const ShellPairData &sp_data); // TODO: rename to overlapD1Kernel

        /**
         * \ingroup ints
         * Calculates the kinetic energy integral matrix.
         */
        vec2d kineticEnergy(const Structure &structure);

        /**
         * Calculates a batch of normalized kinetic energy integrals for the shell pair 'ipair'.
         * In spherical basis.
         */
        vec2d kineticEnergyKernel(const int ipair, const ShellPairData &sp_data);

        /**
         * Calculates a batch of normalized kinetic energy integral derivatives for the shell pair
         * 'ipair'. The derivatives are given as (Ax, Ay, Az, Bx, By, Bz). In spherical basis.
         */
        std::array<vec2d, 6> kineticEnergyD1Kernel(const int ipair, const ShellPairData &sp_data); // TODO: rename to kineticEnergyD1Kernel

        /**
         * \ingroup ints
         * Calculates the nuclear attraction integral matrix.
         */
        vec2d nuclearAttraction(const Structure &structure);

        /**
         * Calculates the one-electron Coulombic interaction integral matrix for given point
         * charges.
         */
        vec2d externalCharges(const std::vector<std::array<double, 4>> &point_charges,
                              const Structure &structure);

        /**
         * Calculates a batch of normalized Coulombic interaction energy integrals for the shell
         * pair 'ipair'. In spherical basis. The charges should be given as a list
         * {(x, y, z, charge)}, with xyz-coordinates in atomic units. The Boys grid should be
         * initialized for lab = la + lb in the given shell pair data.
         */
        vec2d externalChargesKernel(const int ipair, const std::vector<std::array<double, 4>> &charges,
                                    const BoysGrid &boys_grid, const ShellPairData &sp_data);

        /**
         * Calculates a batch of normalized Coulombic interaction energy integral derivatives for
         * the shell pair 'ipair'. The derivatives are given as (Ax, Ay, Az, Bx, By, Bz).
         * In spherical basis. The charges should be given as a list  {(x, y, z, charge)}, with
         * xyz-coordinates in atomic units. The Boys grid should be initialized for lab = la + lb
         * in the given shell pair data.
         */
        std::array<vec2d, 6>
        externalChargesD1Kernel(const int ipair, const std::vector<std::array<double, 4>> &charges,
                                const BoysGrid &boys_grid, const ShellPairData &sp_data); // TODO: rename to externalChargesD1Kernel

        /**
         * Calculates a batch of normalized Coulombic operator derivative integrals for the shell
         * pair 'ipair'. The derivatives are given for each charge as (Ax, Ay, Az). In spherical basis.
         * The charges should be given as a list  {(x, y, z, charge)}, with xyz-coordinates in
         * atomic units. The Boys grid should be initialized for lab = la + lb in the given shell
         * pair data.
         */
        std::vector<std::array<vec2d, 3>>
        externalChargesOperatorD1Kernel(const int ipair, const std::vector<std::array<double, 4>> &charges,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data); // TODO: rename to externalChargesOperatorD1Kernel

        /**
         * \ingroup ints
         * Calculates the dipole moment integral matrices for the \f$x,y,z\f$-directions.
         */
        std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                          const Structure &structure);

        /**
         * Calculates a batch of normalized dipole moment integrals for the shell pair 'ipair'.
         * In spherical basis. The integrals are given as (x, y, z). The origin is expected in
         * atomic units (bohr).
         */
        std::array<vec2d, 3> dipoleMomentKernel(const int ipair, const std::array<double, 3> &origin,
                                                const ShellPairData &sp_data);

        /**
         * \ingroup ints
         * Calculates the two-center ERI matrix over the auxiliary basis functions, \f$(P|Q)\f$.
         */
        vec2d eri2(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the diagonal of the four-center ERI matrix, \f$(\mu\nu|\mu\nu)\f$.
         */
        vec2d eri4Diagonal(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the three-center ERIs, \f$(\mu\nu|P)\f$.
         */
        vec3d eri3(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the four-center ERIs, \f$(\mu\nu|\kappa\tau)\f$.
         */
        vec4d eri4(const Structure &structure);

        /**
         * \ingroup ints
         * Runs a benchmark of calculating the four-center ERIs for all shell quartets,
         * \f$(l_a l_b|l_c l_d)\f$ with \f$l_a \geq l_b\f$ and \f$(l_a l_b) \geq (l_c l_d)\f$.
         */
        void eri4Benchmark(const Structure &structure);

        /** */
        void eri4BenchmarkTest(const Structure &structure);

        /**
         * \ingroup ints
         * Returns the names of all available basis sets in lower case.
         */
        std::set<std::string> availableBasisSets();

        /**
         * \ingroup ints
         * Returns the names of all available auxiliary basis sets in lower case.
         */
        std::set<std::string> availableBasisSetsAux();

        /**
         * \ingroup ints
         * Typedef for bundling Gaussian primitive exponents and contraction coefficients
         * of a shell.
         */
        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        /**
         * \ingroup ints
         * Returns the basis set for the given atom. The exponents and contraction coefficients
         * are listed for every angular momentum.
         */
        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtom(const int atomic_nr, const std::string &basis_set);

        /**
         * \ingroup ints
         * Returns the auxiliary basis set the given atom. The exponents and contraction
         * coefficients are listed for every angular momentum.
         */
        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtomAux(const int atomic_nr, const std::string &aux_basis_set);

        /**
         * \ingroup ints
         * Returns the basis set for the given atoms. The exponents and contraction coefficients
         * are listed for every angular momentum per atom.
         */
        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);

        /**
         * \ingroup ints
         * Returns the auxiliary basis set for the given atoms. The exponents and contraction
         * coefficients are listed for every angular momentum per atom.
         */
        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtomsAux(const std::set<int> &atomic_nrs, const std::string &aux_basis_set);

        /**
         * \ingroup ints
         * Returns the Cartesial to spherical basis transformation,
         * \f$\{(\text{i_spherical}, \text{i_cartesian}, \text{val})\}\f$.
         */
        std::vector<std::tuple<int, int, double>> sphericalTrafo(const int l);

        /** */
        ERI4Kernel deployERI4Kernel(const ShellPairData &sp_data_ab,
                                    const ShellPairData &sp_data_cd);

        /** */
        ERI3Kernel deployERI3Kernel(const ShellPairData &sp_data_ab,
                                    const ShellData &sp_data_cd);

        /** */
        ERI2Kernel deployERI2Kernel(const ShellData &sp_data_a,
                                    const ShellData &sp_data_b);

        /** */
        ERI4D1Kernel deployERI4D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellPairData &sp_data_cd);

        /** */
        ERI3D1Kernel deployERI3D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellData &sp_data_b);

        /** */
        ERI2D1Kernel deployERI2D1Kernel(const ShellData &sp_data_a,
                                        const ShellData &sp_data_b);

        /**
         * \ingroup ints
         * Constructs the shell data corresponding to the auxilary basis set.
         */
        ShellData shellDataAux(const int l, const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell pair data corresponding to the main basis set.
         */
        ShellPairData shellPairDataSymm(const int la, const int lb, const Structure &structure);

        /**
         * \ingroup ints
         *
         */
        ShellPairData shellPairDataNoSymm(const int la, const int lb, const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell datas for the auxiliary basis set, up to l_max.
         */
        std::vector<ShellData> shellDatasAux(const int l_max, const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell pair datas for the given l-pairs.
         */
        std::vector<ShellPairData>
        shellPairDatasSymm(const std::vector<std::pair<int, int>> &l_pairs,
                           const Structure &structure);

        /** */
        std::vector<ShellPairData>
        shellPairDatasNoSymm(const std::vector<std::pair<int, int>> &l_pairs,
                             const Structure &structure);

        /** */
        vec3d ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                 const double PA, const double PB, const double Kab);

        /** */
        vec3d ecoeffsRecurrence2_n1(const double a, const double b, const int la, const int lb,
                                    const double A, const double B, const vec3d &ecoeffs);

        /** */
        vec2d ecoeffsRecurrence1(const double one_o_2a, const int l);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                  const int lb, const double *xyz_a,
                                                  const double *xyz_b);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair_n1(const double a, const double b, const int la,
                                                     const int lb, const double *xyz_a,
                                                     const double *xyz_b,
                                                     const std::array<lible::vec3d, 3> &ecoeffs);

        /** */
        std::array<lible::vec2d, 3> ecoeffsPrimitive(const double a, const int l);

        /** */
        std::vector<lible::vec3d>
        ecoeffsShellPair_Eij0(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                              const double *exps_a, const double *exps_b, const double *xyz_a,
                              const double *xyz_b);

        /** */
        std::vector<lible::vec4d>
        ecoeffsShellPair_Eijt(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                              const double *exps_a, const double *exps_b, const double *xyz_a,
                              const double *xyz_b);

        /** */
        std::vector<std::vector<double>> ecoeffsShell(const int l, const std::vector<double> &exps);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives in each shell pair.
         */
        std::vector<std::vector<vec3d>>
        ecoeffsSPData_Eij0(const ShellPairData &sp_data);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives in each shell pair.
         */
        std::vector<std::vector<vec4d>>
        ecoeffsSPData_Eijt(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<double>
        ecoeffsSphericalSPData_Bra(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalSPData_BraKet(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<double>
        ecoeffsSphericalShellData_Bra(const ShellData &sh_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalShellData_BraKet(const ShellData &sh_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<std::vector<double>>
        ecoeffsSphericalSPDatas_Bra(const std::vector<ShellPairData> &sp_datas);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsSphericalSPDatas_BraKet(const std::vector<ShellPairData> &sp_datas);

        /** TODO: */
        class BoysGrid;

        /** TODO: */
        std::vector<double> calcBoysF(const int max_n, const double x, const BoysGrid &boys_grid);

        /** TODO: */
        vec3d calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx);

        /**
         * \ingroup ints
         * Returns the number of Cartesian Gaussians.
         */
        constexpr int numCartesians(const int l);

        /**
         * \ingroup ints
         * Returns the number of spherical Gaussians.
         */
        constexpr int numSphericals(const int l);

        /**
         * \ingroup ints
         * Returns the number of Hermite Gaussians.
         */
        constexpr int numHermites(const int l);

        /**
         * \ingroup ints
         * Returns the exponents of a Cartesian Gaussian \f$(x,y,z)\f$-directions for the given
         * angular momentum.
         */
        std::vector<std::array<int, 3>> cartExps(const int l);

        /**
         * \ingroup ints
         * Returns a list of angular momentum pairs such that la >= lb:
         *   {(0, 0), (1, 0), (1, 1), ..., (l_max, l_max)}.
         */
        std::vector<std::pair<int, int>> getLPairsSymm(const int l_max);

        /**
         * \ingroup ints
         * Returns a list of angular momentum pairs: {(0, 0), (1, 0), (0, 1), ..., (l_max, l_max)}.
         */
        std::vector<std::pair<int, int>> getLPairsNoSymm(const int l_max);

#ifdef _LIBLE_USE_HIP_
        namespace gpu
        {
            /**
             *
             */
            vec2d overlap0(const Structure &structure);

            /**
             *
             */
            vec2d overlap(const Structure &structure);
        }
#endif
    }
}
