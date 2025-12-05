#include <lible/ints/ints.hpp>
#include <lible/ints/shell.hpp>
#include <lible/ints/ecoeffs.hpp.depr>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <cmath>
#include <numbers>

namespace lints = lible::ints;

std::vector<double> lints::calcShellNorms(const int l, const std::vector<double> &coeffs,
                                          const std::vector<double> &exps,
                                          const std::vector<double> &primitive_norms)
{
    size_t cdepth = exps.size();
    std::vector<std::vector<double>> e_coeffs = ecoeffsShell(l, exps);

    int dim_cart = numCartesians(l);
    vec2d ints_cart(Fill(0), dim_cart, dim_cart);
    for (size_t ia = 0, iab = 0; ia < cdepth; ia++)
        for (size_t ib = 0; ib < cdepth; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];
            double da = coeffs[ia] * primitive_norms[ia];
            double db = coeffs[ib] * primitive_norms[ib];

            double p = a + b;
            double da_x_db = da * db;

            double fac = da_x_db * std::pow(std::numbers::pi / p, 1.5);

            for (int mu = 0, munu = 0; mu < dim_cart; mu++)
                for (int nu = 0; nu < dim_cart; nu++, munu++)
                    ints_cart(mu, nu) += fac * e_coeffs[iab][munu];
        }

    vec2d ints_sph = trafo2Spherical(l, l, ints_cart);

    int dim_sph = numSphericals(l);
    std::vector<double> norms(dim_sph);
    for (int i = 0; i < dim_sph; i++)
        norms[i] = 1.0 / std::sqrt(ints_sph(i, i));

    return norms;
}

std::vector<lints::Shell>
lints::constructShells(const basis_atoms_t &basis_atoms,
                       const std::vector<std::array<double, 3>> &coords_atoms)
{
    if (basis_atoms.size() != coords_atoms.size())
        throw std::runtime_error("constructShells(): number of atoms in the basis set and the list "
            "of coordinates doesn't match");

    size_t idx_shell = 0;
    size_t ofs_sph = 0;
    size_t ofs_cart = 0;

    std::vector<Shell> shells;
    for (size_t iatom = 0; iatom < basis_atoms.size(); iatom++)
    {
        std::array<double, 3> coords_iatom = coords_atoms[iatom];
        const auto &[atomic_nr, basis_atom] = basis_atoms[iatom];

        for (const auto &[l, exps, coeffs] : basis_atom)
        {
            const int dim_cart = numCartesians(l);
            const int dim_sph = numSphericals(l);

            if (exps.size() != coeffs.size())
                throw std::runtime_error("constructShells(): number of exponents/coefficients "
                    "doesn't match");

            const size_t cdepth = exps.size();

            std::vector<double> primitive_norms(cdepth);
            for (size_t i = 0; i < cdepth; i++)
                primitive_norms[i] = purePrimitiveNorm(l, exps[i]);

            std::vector<double> norms = calcShellNorms(l, coeffs, exps, primitive_norms);

            Shell shell(l, atomic_nr, iatom, dim_cart, dim_sph, ofs_cart, ofs_sph, idx_shell,
                        coords_iatom, exps, coeffs, norms, primitive_norms);

            shells.push_back(shell);

            idx_shell++;
            ofs_sph += dim_sph;
            ofs_cart += dim_cart;
        }
    }

    return shells;
}
