#include <lible/ints/shell.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#define _USE_MATH_DEFINES

namespace LI = lible::ints;

using std::array, std::vector;

namespace lible::ints
{
	std::ostream &operator<<(std::ostream &os, const LI::Shell &obj)
	{
		os << "  l        : " << obj.l << "\n";
		os << "  z        : " << obj.z << "\n";
		os << "  atm index: " << obj.atom_idx << "\n";
		os << "  dim_cart : " << obj.dim_cart << "\n";
		os << "  dim_sph  : " << obj.dim_sph << "\n";
		os << "  position : " << obj.pos << "\n";
		   << obj.xyz_coords[0] << ", " 
		   << obj.xyz_coords[1] << ", " 
		   << obj.xyz_coords[2] << ")\n";
		os << "  coeffs   :";
		for (const auto &coeff : obj.coeffs) 
		{
		    os << coeff << " ";
		}
		os << "\n";
		os << "  exps     : ";
		for (const auto &exp : obj.exps) 
		{
		    os << exp << " ";
		}
		os << "\n";
		os << "  norms    : ";
		for (const auto &norm : obj.norms) 
		{
		    os << norm << " ";
		}
		os << "\n";
		return os;
	}
}


vector<double> LI::calcShellNorms(const int l, const vector<double> &coeffs,
                                  const vector<double> &exps)
{    
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

            double p = a + b;
            double da_x_db = da * db;

            double fac = da_x_db * std::pow(M_PI / p, 1.5);

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
