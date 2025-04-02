#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/structure.hpp>
#include <lible/ints/utils.hpp>

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
         * \ingroup ints
         * Calculates the kinetic energy integral matrix.
         */
        vec2d kineticEnergy(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the nuclear attraction integral matrix.
         */
        vec2d nuclearAttraction(const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the dipole moment integral matrices for the \f$x,y,z\f$-directions.
         */
        std::array<vec2d, 3> dipoleMoment(const std::array<double, 3> &origin,
                                          const Structure &structure);

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

        /** */
        vec4d eri4New(const Structure &structure);

        /**
         * \ingroup ints
         * Runs a benchmark of calculating the four-center ERIs for all shell quartets,
         * \f$(l_a l_b|l_c l_d)\f$ with \f$l_a \geq l_b\f$ and \f$(l_a l_b) \geq (l_c l_d)\f$.
         */
        void eri4Benchmark(const Structure &structure);

        /** */
        void eri4BenchmarkNew(const Structure &structure);

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
         * Returns the exponents of a Cartesian Gaussian \f$(x,y,z)\f$-directions for the given
         * angular momentum.
         */
        std::vector<std::array<int, 3>> cartExps(const int l);

        /**
         * \ingroup ints
         * Returns the Cartesial to spherical basis transformation,
         * \f$\{(\text{i_spherical}, \text{i_cartesian}, \text{val})\}\f$.
         */
        std::vector<std::tuple<int, int, double>> sphericalTrafo(const int l);

        /**
         * \ingroup ints
         * Type alias for the four-center two-electron repulsion integral kernel function.
         */
        using kernel_eri4_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                 const int cdepth_c, const int cdepth_d,
                                                 const double *exps_a, const double *exps_b,
                                                 const double *exps_c, const double *exps_d,
                                                 const double *coords_a, const double *coords_b,
                                                 const double *coords_c, const double *coords_d,
                                                 const double *ecoeffs_ab,
                                                 const double *ecoeffs_cd_tsp,
                                                 double *eri4_batch)>;

        /**
         * \ingroup ints
         * Type alias for the three-center two-electron repulsion integral kernel function.
         */
        using kernel_eri3_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                 const int cdepth_c, const double *exps_a,
                                                 const double *exps_b, const double *exps_c,
                                                 const double *coords_a, const double *coords_b,
                                                 const double *coords_c, const double *ecoeffs_ab,
                                                 const double *ecoeffs_c, double *eri3_batch)>;

        /**
         * \ingroup ints
         * Type alias for the two-center two-electron repulsion integral kernel function.
         */
        using kernel_eri2_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                 const double *exps_a, const double *exps_b,
                                                 const double *coords_a, const double *coords_b,
                                                 const double *ecoeffs_a,
                                                 const double *ecoeffs_b_tsp,
                                                 double *eri2_batch)>;

        /**
         * \ingroup ints
         * Function for deploying a kernel function for calculating the four-center two-electron
         * repulsion integrals.
         */
        kernel_eri4_t deployERI4Kernel(const int la, const int lb, const int lc, const int ld);

        /**
         * \ingroup ints
         * Function for deploying a kernel function for calculating the three-center two-electron
         * repulsion integrals.
         */
        kernel_eri3_t deployERI3Kernel(const int la, const int lb, const int ld);

        /**
         * \ingroup ints
         * Function for deploying a kernel function for calculating the two-center two-electron
         * repulsion integrals.
         */
        kernel_eri2_t deployERI2Kernel(const int la, const int lb);

        /**
         * \ingroup ints
         * Constructs the shell data corresponding to the auxilary basis set.
         */
        ShellData constructShellDataAux(const int l, const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell pair data corresponding to the main basis set.
         */
        ShellPairData constructShellPairData(const int la, const int lb,
                                             const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell datas for the auxiliary basis set, up to l_max.
         */
        std::vector<ShellData>
        constructShellDatasAux(const int l_max, const Structure &structure);

        /**
         * \ingroup ints
         * Constructs the shell pair datas for the given l-pairs.
         */
        std::vector<ShellPairData>
        constructShellPairDatas(const std::vector<std::pair<int, int>> &l_pairs,
                                const Structure &structure);

        /**
         * \ingroup ints
         * Calculates the Hermite expansion coefficients for the given l-pairs and shell-pair
         * datas. The function assumes the given shell-pair datas correspond to the l-pairs.
         * The E-coefficients in the ket are transposed.
         */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsFromSPDatas(const std::vector<std::pair<int, int>> &l_pairs,
                           const std::vector<ShellPairData> &sp_datas);

        /**
         * \ingroup ints
         * Calculates the Hermite expansion coefficients for the given l-value. It is assumed
         * that the shell datas are ordered as 0,...,l_max_aux. The function assumes shell data
         * for the auxiliary basis set.
         */
        std::vector<std::vector<double>>
        ecoeffsFromShellDatasAux(const int l_max_aux, const std::vector<ShellData> &sh_datas);

        /**
         * \ingroup ints
         * Calculates the Hermite expansion coefficients for the given l value. It is assumed
         * that the shell datas are ordered as 0,...,l_max_aux. The function assumes shell data
         * for the auxiliary basis set. The E-coefficients in the ket are transposed.
         */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsFromShellDatasAuxPairs(const int l_max_aux, const std::vector<ShellData> &sh_datas);

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
         * Returns a list of angular momentum pairs such that la >= lb: {(0, 0), (1, 0), (1, 1), ..., (l_max, l_max)}.
         */
        std::vector<std::pair<int, int>> returnLPairs(const int l_max);

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
