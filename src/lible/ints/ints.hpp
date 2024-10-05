#pragma once

#include <lible/types.hpp>
#include <lible/ints/structure.hpp>

#include <cassert>
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
        std::array<vec2d, 3> dipoleMoment(const Structure &structure);

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
        vec3d eri3();

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

        /**
         * \ingroup ints
         * Returns the names of all available basis sets.
         */
        std::set<std::string> availableBasisSets();

        /**
         * \ingroup ints
         * Returns the names of all available auxiliary basis sets.
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
