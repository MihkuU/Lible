#pragma once

#include <array>
#include <utility>
#include <vector>

namespace lible
{
	/*
	 * TODO: clean it up!
	 * Here we deal with atomic orbital shells. We follow conventions outlined in https://iodata.readthedocs.io/en/latest/basis.html.
	 * Definitions:
	 *   Shell - set of basis functions with the same angular momentum, contraction coefficients and
	 *   exponents of gaussian primitives.
	 *
	 */
	namespace ints
	{
		/**		 
		 * Data structure for representing a shell of atomic orbitals. Groups together various
		 * data that defines a shell: angular momentum, contraction coefficients, contraction
		 * exponents etc. 		
		 */
		struct Shell
		{
			/** The main constructor. */
			Shell(const int l, const int z, const size_t dim_cart, const size_t dim_sph,
				  const size_t pos, const size_t pos_cart, const std::array<double, 3> &xyz_coords,
				  const std::vector<double> &coeffs, const std::vector<double> &coeffs_raw,
				  const std::vector<double> &exps, const std::vector<double> &norms)
				: l(l), z(z), dim_cart(dim_cart), dim_sph(dim_sph), pos(pos), pos_cart(pos_cart),
				  xyz_coords(xyz_coords), coeffs(coeffs), coeffs_raw(coeffs_raw), exps(exps),
				  norms(norms)
			{
			}

			/** Angular momentum. */
			const int l;

			/** Atomic number. */
			const int z;

			/** Number of atomic orbitals in Cartesian basis. */
			const size_t dim_cart;

			/** Number of atomic orbitals in spherical basis. */
			const size_t dim_sph;

			/** Starting position in the list of atomic orbitals. */
			const size_t pos;

			/** Starting position in the list of atomic orbitals in Cartesian basis. */
			const size_t pos_cart;

			/** Coordinates of the atom corresponding to the shell. */ // TODO: specify in angstrom
			const std::array<double, 3> xyz_coords;

			/** Contraction coefficients with primitive norms multiplied into. */
			const std::vector<double> coeffs;

			/** Contraction coefficients without primitive norms. */
			const std::vector<double> coeffs_raw;

			/** Exponents of the Gaussian primitives. */
			const std::vector<double> exps;

			/** Normalization constants of the atomic orbitals in spherical basis. */
			const std::vector<double> norms;
		};

		/**
		 * Calculates the norms of the shell atomic orbitals in spherical basis.
		 */
		std::vector<double> calcShellNorms(const int l, const std::vector<double> &coeffs,
										   const std::vector<double> &exps);
	}
}