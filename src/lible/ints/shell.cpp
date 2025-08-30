#include <lible/ints/shell.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <numbers>

namespace lints = lible::ints;

using std::array, std::vector;

vector<double> lints::calcShellNorms(const int l, const vector<double> &coeffs, 
                                     const vector<double> &exps, 
                                     const vector<double> &primitive_norms)
{
    // TODO: error handling
    size_t k = exps.size();
    vector<vector<double>> e_coeffs = ecoeffsShell(l, exps);

    int dim_cart = numCartesians(l);
    vec2d ints_cart(Fill(0), dim_cart, dim_cart);
    for (size_t ia = 0, iab = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++, iab++)
        {
            double a = exps[ia];
            double b = exps[ib];
            double da = coeffs[ia];
            double db = coeffs[ib];

            double na = primitive_norms[ia];
            double nb = primitive_norms[ib];
            da *= na;
            db *= nb;

            double p = a + b;
            double da_x_db = da * db;

            double fac = da_x_db * std::pow(std::numbers::pi / p, 1.5);

            for (int mu = 0, munu = 0; mu < dim_cart; mu++)
                for (int nu = 0; nu < dim_cart; nu++, munu++)
                    ints_cart(mu, nu) += fac * e_coeffs[iab][munu];
        }

    vec2d ints_sph = trafo2Spherical(l, l, ints_cart);

    int dim_sph = numSphericals(l);
    vector<double> norms(dim_sph);
    for (int i = 0; i < dim_sph; i++)
        norms[i] = 1.0 / std::sqrt(ints_sph(i, i));

    return norms;
}

vector<lints::Shell> lints::constructShells(const basis_atoms_t &basis_atoms,
                                            const vector<int> &atomic_nrs,
                                            const vector<array<double, 3>> &coords_atoms)
{
    // TODO: error handling
    int idx_shell = 0;
    int ofs_sph = 0;
    int ofs_cart = 0;

    vector<Shell> shells;
    for (size_t iatom = 0; iatom < atomic_nrs.size(); iatom++)
    {
        array<double, 3> coords_iatom = coords_atoms[iatom];

        int atomic_nr = atomic_nrs[iatom];
        basis_atom_t basis_atom = basis_atoms.at(atomic_nr);

        for (const auto &[l, exps_coeffs_list] : basis_atom)
        {
            int dim_cart = numCartesians(l);
            int dim_sph = numSphericals(l);

            for (const auto &[exps, coeffs] : exps_coeffs_list)
            {
                if (exps.size() != coeffs.size())
                    throw std::runtime_error("number of exponents/coefficients doesn't match");

                size_t cdepth = exps.size();

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