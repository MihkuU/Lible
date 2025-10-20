#include <lible/ints/ints.hpp>
#include <lible/ints/shell.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <cmath>
#include <numbers>

namespace lints = lible::ints;

using std::array, std::vector;

vector<double> lints::calcShellNorms(const int l, const vector<double> &coeffs, 
                                     const vector<double> &exps, 
                                     const vector<double> &primitive_norms)
{
    const size_t k = exps.size();
    const vector<vector<double>> e_coeffs = ecoeffsShell(l, exps);

    const int dim_cart = numCartesians(l);
    vec2d ints_cart(Fill(0), dim_cart, dim_cart);
    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            const double a = exps[ia];
            const double b = exps[ib];
            const double da = coeffs[ia] * primitive_norms[ia];
            const double db = coeffs[ib] * primitive_norms[ib];

            const double p = a + b;
            const double da_x_db = da * db;

            const double fac = da_x_db * std::pow(std::numbers::pi / p, 1.5);

            for (int mu = 0, munu = 0; mu < dim_cart; mu++)
                for (int nu = 0; nu < dim_cart; nu++, munu++)
                    ints_cart(mu, nu) += fac * e_coeffs[iab][munu];
        }

    vec2d ints_sph = trafo2Spherical(l, l, ints_cart);

    const int dim_sph = numSphericals(l);
    vector<double> norms(dim_sph);
    for (int i = 0; i < dim_sph; i++)
        norms[i] = 1.0 / std::sqrt(ints_sph(i, i));

    return norms;
}

vector<lints::Shell> lints::constructShells(const basis_atoms_t &basis_atoms,
                                            const vector<int> &atomic_nrs,
                                            const vector<array<double, 3>> &coords_atoms)
{
    if (atomic_nrs.size() != coords_atoms.size())
        throw std::runtime_error("lints::constructShells(): number of atomic numbers and coordinates "
                                 "doesn't match");

    int idx_shell = 0;
    int ofs_sph = 0;
    int ofs_cart = 0;

    vector<Shell> shells;
    for (size_t iatom = 0; iatom < atomic_nrs.size(); iatom++)
    {
        const array<double, 3> coords_iatom = coords_atoms[iatom];

        const int atomic_nr = atomic_nrs[iatom];
        basis_atom_t basis_atom = basis_atoms.at(atomic_nr);

        for (const auto &[l, exps_coeffs_list] : basis_atom)
        {
            const int dim_cart = numCartesians(l);
            const int dim_sph = numSphericals(l);

            for (const auto &[exps, coeffs] : exps_coeffs_list)
            {
                if (exps.size() != coeffs.size())
                    throw std::runtime_error("number of exponents/coefficients doesn't match");

                const size_t cdepth = exps.size();

                vector<double> primitive_norms(cdepth);
                for (size_t i = 0; i < cdepth; i++)
                    primitive_norms[i] = purePrimitiveNorm(l, exps[i]);

                vector<double> norms = calcShellNorms(l, coeffs, exps, primitive_norms);

                Shell shell(l, atomic_nr, iatom, dim_cart, dim_sph, ofs_cart, ofs_sph, idx_shell,
                            coords_iatom, exps, coeffs, norms, primitive_norms);

                shells.push_back(shell);

                idx_shell++;
                ofs_sph += dim_sph;
                ofs_cart += dim_cart;
            }
        }
    }

    return shells;
}