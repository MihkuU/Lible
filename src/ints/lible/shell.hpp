#pragma once

#define _USE_MATH_DEFINES

#include <array>
#include <utility>
#include <vector>

namespace lible
{
	/*
	 * Here we deal with atomic orbital shells. We follow conventions outlined in https://iodata.readthedocs.io/en/latest/basis.html.
	 * Definitions:
	 *   Shell - set of basis functions with the same angular momentum, contraction coefficients and
	 *   exponents of gaussian primitives.
	 *
	 */
	namespace ints
	{
		struct Shell
		{
			Shell(const int &angular_momentum,
				  const int &atomic_number,
				  const std::size_t &dim_cartesian,
				  const std::size_t &dim_spherical,
				  const std::size_t &pos,
				  const std::array<double, 3> &xyz_coords,
				  const std::vector<double> &coeffs,
				  const std::vector<double> &coeffs_raw,
				  const std::vector<double> &exps,
				  const std::vector<double> &norms)
				: angular_momentum(angular_momentum),
				  atomic_number(atomic_number),
				  dim_cartesian(dim_cartesian),
				  dim_spherical(dim_spherical),
				  pos(pos),
				  xyz_coords(xyz_coords),
				  coeffs(coeffs),
				  coeffs_raw(coeffs_raw),
				  exps(exps),
				  norms(norms)
			{
			}

			Shell(const int &angular_momentum,
				  const std::size_t &dim_cartesian,
				  const std::size_t &dim_spherical)
				: angular_momentum(angular_momentum),
				  atomic_number(0),
				  dim_cartesian(dim_cartesian),
				  dim_spherical(dim_spherical),
				  pos(0),
				  xyz_coords({}),
				  coeffs({}),
				  coeffs_raw({}),
				  exps({}),
				  norms({})
			{
			}

			const int angular_momentum;
			const int atomic_number;
			const std::size_t dim_cartesian;
			const std::size_t dim_spherical;
			const std::size_t pos;
			const std::array<double, 3> xyz_coords;
			const std::vector<double> coeffs;
			const std::vector<double> coeffs_raw;
			const std::vector<double> exps;
			const std::vector<double> norms;
		};

		std::vector<double> calcShellNorms(const int l, const std::vector<double> &coeffs,
										   const std::vector<double> &exps);
	}
}