#include <cmath>
#include <limits>
// #include <numbers>
#include <stdexcept>

#include "defs_geomopt.h"
#include "geometry.h"
#include <fmt/core.h>

using namespace lible::geomopt;
using std::pair;
using std::size_t;
using std::string;
using std::vector;

// Geometry::Geometry(const vector<double> &coords_cart_, const vector<string> &atoms_) : coords_cart(coords_cart_), atoms(atoms_)
// {
//     assert((coords_cart.size() % 3 == 0));
//     assert((coords_cart.size() / 3 == atoms.size()));

//     for (const string &atom : atoms)
//         if (!geomoptdefs::atomic_numbers.contains(atom))
//             throw std::runtime_error(fmt::format("Undefined element {}!", atom));

//     n_atoms = atoms.size();
//     atom_coords_cart.resize(n_atoms);
//     for (size_t matom = 0; matom < n_atoms; matom++)
//         atom_coords_cart[matom] = arma::vec::fixed<3>({coords_cart[3 * matom],
//                                                        coords_cart[3 * matom + 1],
//                                                        coords_cart[3 * matom + 2]});

//     atom_bonding_partners = findAtomBondingPartners(atom_coords_cart);
//     red_int_coords = constructRedIntCoords(atom_coords_cart, atom_bonding_partners);
// }

// vector<double> Geometry::transformCoordsCartToRedint(const vector<double> &coords_cart)
// {
//     vector<double> coords_redint;

//     return coords_redint;
// }

// vector<double> Geometry::transformCoordsRedintToCart(const vector<double> &coords_redint_in)
// {    
//     // size_t max_trafo_iters = 5; //TODO: where to put it?

//     // b_matrix = constructBMatrix(); // TODO: make a separate function for this?
//     // b_matrix_pinv = arma::pinv(b_matrix); // TODO: make a separrate function for this?      
    
//     // arma::dvec coords_redint = arma::conv_to<arma::dvec>::from(coords_redint_in);

//     // arma::dvec coords_cart_target = b_matrix_pinv * coords_redint;

//     // std::cout << "coords_cart:" << std::endl;
//     // std::cout << arma::conv_to<arma::dvec>::from(coords_cart) << std::endl;
        
//     // for (size_t k = 0; k < max_trafo_iters; k++)
//     // {        
//     //     vector<dvec3> atom_coords_cart(n_atoms);
//     //     for (size_t matom = 0; matom < n_atoms; matom++)
//     //         atom_coords_cart[matom] = arma::vec::fixed<3>({coords_cart_target[3 * matom],
//     //                                                        coords_cart_target[3 * matom + 1],
//     //                                                        coords_cart_target[3 * matom + 2]});

//     //     RedIntCoords red_int_coords = constructRedIntCoords(atom_coords_cart, atom_bonding_partners);
//     //     arma::dvec coords_redint_k = arma::conv_to<arma::dvec>::from(red_int_coords.coords_redint);

//     //     arma::dvec delta_coords_redint = coords_redint - coords_redint_k;
//     //     double norm_delta = arma::norm(delta_coords_redint);
//     //     printf("norm_delta = %16.12lf\n", norm_delta);

//     //     coords_cart_target += b_matrix_pinv * delta_coords_redint;
//     //     std::cout << "coords_cart_target:" << std::endl;
//     //     std::cout << coords_cart_target << std::endl;            
//     // }

//     // vector<double> coords_cart_out;
//     // return coords_cart_out;
// }

// vector<double> Geometry::transformGradCartToRedint(const vector<double> &grad_cart_in)
// {    
//     arma::dmat b_matrix = constructBMatrix(red_int_coords);
//     arma::dmat b_matrix_inv = arma::pinv(b_matrix);
//     arma::dmat b_matrix_t = b_matrix.t();
//     arma::dmat b_matrix_t_pinv = arma::pinv(b_matrix_t);

//     arma::dvec grad_cart = arma::conv_to<arma::dvec>::from(grad_cart_in);
//     arma::dvec grad_redint = b_matrix_t_pinv * grad_cart;
//     vector<double> grad_redint_out = arma::conv_to<vector<double>>::from(grad_redint);    

//     return grad_redint_out;
// }

// vector<double> Geometry::transformStepRedIntToCart(const vector<double> &coords_cart,
//                                                    const vector<double> &coords_redint,
//                                                    const vector<double> &step_redint)
// {       
//     const arma::dvec dq = arma::conv_to<arma::dvec>::from(step_redint);
//     const arma::dvec q0 = arma::conv_to<arma::dvec>::from(coords_redint);
//     const arma::dvec x0 = arma::conv_to<arma::dvec>::from(coords_cart);

//     arma::dmat b_matrix = constructBMatrix(red_int_coords); 
//     arma::dmat b_matrix_t = b_matrix.t();
//     arma::dmat b_matrix_pinv = arma::pinv(b_matrix);
    
//     arma::dvec xk = x0;
//     arma::dvec qk = q0;
//     double conv_thrs = 1e-12;
//     size_t max_trafo_iter = 25; //TODO: figure out where to put?
//     for (size_t k = 0; k < max_trafo_iter; k++)
//     {
//         arma::dvec diff_q = (q0 + dq - qk);

//         arma::dvec xk_ = xk;
//         xk = xk + b_matrix_pinv * diff_q;

//         double norm_diff = arma::norm(xk - xk_);
//         printf("norm_diff = %16.12lf\n", norm_diff);
//         double rmse = sqrt(arma::accu(arma::pow(xk - xk_, 2)) / xk.n_elem);
//         if (rmse < conv_thrs)
//         {
//             printf("converged at k = %d\n", k); // TMP
//             break;
//         }

//         vector<dvec3> xk_atom_coords(n_atoms);
//         for (size_t m = 0; m < n_atoms; m++) 
//             xk_atom_coords[m] = arma::vec::fixed<3>({xk[3 * m], xk[3 * m + 1], xk[3 * m + 2]});

//         RedIntCoords red_int_coords = constructRedIntCoords(xk_atom_coords, atom_bonding_partners); //TODO: simplify - 
//         qk = arma::conv_to<arma::dvec>::from(red_int_coords.coords_redint);

//         // b_matrix = constructBMatrix(red_int_coords);
//         // b_matrix_pinv = arma::pinv(b_matrix);
//     }

//     vector<double> step_cart = arma::conv_to<vector<double>>::from(xk - x0);
//     return step_cart;
// }

// double Geometry::calcDistance(const size_t &matom, const size_t &natom,
//                               const vector<dvec3> &atom_coords_cart)
// {
//     return arma::norm(atom_coords_cart[matom] - atom_coords_cart[natom]);
// }

// double Geometry::calcAngle(const size_t &matom, const size_t &oatom, const size_t &natom, 
//                            const vector<dvec3> &atom_coords_cart)
// {
//     dvec3 u = atom_coords_cart[matom] - atom_coords_cart[oatom];
//     dvec3 v = atom_coords_cart[natom] - atom_coords_cart[oatom];
//     u /= arma::norm(u);
//     v /= arma::norm(v);
//     return acos(arma::dot(u, v));
// }

// double Geometry::calcDihedral(const size_t &matom, const size_t &oatom, const size_t &patom, const size_t &natom,
//                               const vector<dvec3> &atom_coords_cart)
// {
//     dvec3 u = atom_coords_cart[matom] - atom_coords_cart[oatom];    
//     dvec3 v = atom_coords_cart[natom] - atom_coords_cart[patom];
//     dvec3 w = atom_coords_cart[patom] - atom_coords_cart[oatom];
//     u /= arma::norm(u);
//     v /= arma::norm(v);
//     w /= arma::norm(w);
//     double sin_uw = sqrt(1 - pow(arma::dot(u, w), 2));
//     double sin_vw = sqrt(1 - pow(arma::dot(v, w), 2));
//     return acos(arma::dot(arma::cross(u, w), arma::cross(v, w)) / (sin_uw * sin_vw));
// }

// // size_t Geometry::findClosestAtom(const size_t &matom)
// // {
// //     double min_distance = std::numeric_limits<double>::max();
// //     size_t closest_atom;
// //     for (size_t natom = matom + 1; natom < n_atoms; natom++)
// //     {
// //         double distance = calcDistance(matom, natom);
// //         if (distance < min_distance)
// //         {
// //             min_distance = distance;
// //             closest_atom = natom;
// //         }
// //     }
// //     return closest_atom;
// // }

// vector<vector<size_t>> Geometry::findAtomBondingPartners(const vector<dvec3> &atom_coords_cart)
// {
//     vector<double> covalent_radii(n_atoms);
//     for (size_t matom = 0; matom < n_atoms; matom++)
//         covalent_radii[matom] = geomoptdefs::covalent_radii.at(atoms[matom]);

//     vector<vector<size_t>> atom_bonding_partners(n_atoms);
//     for (size_t matom = 0; matom < n_atoms; matom++)
//     {
//         vector<size_t> bonding_partners;
//         for (size_t natom = 0; natom < n_atoms; natom++)
//         {
//             if (matom == natom)
//                 continue; 

//             double sum_cov_radii = covalent_radii[matom] + covalent_radii[natom];
//             double bonding_distance = geomoptdefs::bonding_factor * sum_cov_radii;
//             double distance = calcDistance(matom, natom, atom_coords_cart);
//             if (distance < bonding_distance)
//                 bonding_partners.push_back(natom);
//         }

//         atom_bonding_partners[matom] = bonding_partners;
//     }
//     return atom_bonding_partners;
// }

// Geometry::RedIntCoords Geometry::constructRedIntCoords(const vector<dvec3> &atom_coords_cart, 
//                                                        const vector<vector<size_t>> &atom_bonding_partners)
// {
//     //TODO: remove 'atom_bonding_partners'??    
//     RedIntCoords red_int_coords;

//     // Bonds
//     for (size_t matom = 0; matom < n_atoms; matom++)
//     {        
//         vector<size_t> bonding_partners_matom = atom_bonding_partners[matom];
//         for (const size_t &natom : bonding_partners_matom)
//             if (matom < natom)
//                 red_int_coords.bonds.emplace_back(std::make_pair(std::make_tuple(matom, natom),
//                                                                  calcDistance(matom, natom, 
//                                                                  atom_coords_cart)));
//     }

//     // Angles
//     for (size_t oatom = 0; oatom < n_atoms; oatom++)
//     {
//         vector<size_t> bonding_partners_oatom = atom_bonding_partners[oatom];
//         if (bonding_partners_oatom.size() > 1)
//             for (const size_t &matom : bonding_partners_oatom)
//                 for (const size_t &natom : bonding_partners_oatom)
//                     if (matom < natom)
//                         red_int_coords.angles.emplace_back(std::make_pair(std::make_tuple(matom, oatom, natom),
//                                                                           calcAngle(matom, oatom, natom, 
//                                                                           atom_coords_cart)));
//     }

//     // // Dihedrals
//     // for (size_t oatom = 0; oatom < n_atoms; oatom++)
//     // {
//     //     vector<size_t> bonding_partners_oatom = atom_bonding_partners[oatom];
//     //     if (bonding_partners_oatom.size() > 1)
//     //         for (const size_t &patom : bonding_partners_oatom)
//     //             if (oatom < patom)
//     //             {
//     //                 vector<size_t> bonding_partners_patom = atom_bonding_partners[patom];
//     //                 if (bonding_partners_patom.size() > 1)
//     //                     for (const size_t &matom : bonding_partners_oatom)
//     //                         for (const size_t &natom : bonding_partners_patom)
//     //                             if (matom != patom and natom != oatom)
//     //                                 red_int_coords.dihedrals.emplace_back(std::make_pair(std::make_tuple(matom, oatom, patom, natom),
//     //                                                                                      calcDihedral(matom, oatom, patom, natom, 
//     //                                                                                      atom_coords_cart)));
//     //             }
//     // }

//     red_int_coords.coords_redint.resize(red_int_coords.bonds.size() +
//                                         red_int_coords.angles.size() +
//                                         red_int_coords.dihedrals.size());
//     size_t pos = 0;
//     for (auto &bond : red_int_coords.bonds)
//     {
//         red_int_coords.coords_redint[pos] = bond.second;
//         pos++;
//     }

//     // for (auto &[mn, val] : red_int_coords.bonds)
//     // {
//     //     auto [m, n] = mn;
//     //     printf("(%d, %d) %16.12lf\n", m, n, val);
//     // }    

//     for (auto &angle : red_int_coords.angles)
//     {
//         red_int_coords.coords_redint[pos] = angle.second;
//         pos++;
//     }

//     for (auto &[mon_atom, val] : red_int_coords.angles)
//     {
//         auto [m, o, n] = mon_atom;
//         printf("(%d, %d, %d) %16.12lf\n", m, o, n, 57.2958 * val);
//     }

//     for (auto &dihedral : red_int_coords.dihedrals)
//     {
//         red_int_coords.coords_redint[pos] = dihedral.second;
//         pos++;
//     }

//     return red_int_coords;
// }

// arma::dmat Geometry::constructBMatrix(const RedIntCoords &red_int_coords)
// {
//     arma::dmat bmatrix(red_int_coords.coords_redint.size(), coords_cart.size(), arma::fill::zeros);

//     for (size_t ibond = 0; ibond < red_int_coords.bonds.size(); ibond++)
//     {
//         auto [mn_atom, bond_length] = red_int_coords.bonds[ibond];
//         auto [matom, natom] = mn_atom;

//         dvec3 norm_bond_vector = (atom_coords_cart[matom] - atom_coords_cart[natom]) / bond_length;
//         size_t posm = 3 * matom;
//         size_t posn = 3 * natom;
//         bmatrix(ibond, posm) = norm_bond_vector[0];
//         bmatrix(ibond, posm + 1) = norm_bond_vector[1];
//         bmatrix(ibond, posm + 2) = norm_bond_vector[2];
//         bmatrix(ibond, posn) = -norm_bond_vector[0];
//         bmatrix(ibond, posn + 1) = -norm_bond_vector[1];
//         bmatrix(ibond, posn + 2) = -norm_bond_vector[2];
//     }

//     size_t offset = red_int_coords.bonds.size();
//     for (size_t iangle = 0; iangle < red_int_coords.angles.size(); iangle++)
//     {
//         auto [mon_atom, bond_angle] = red_int_coords.angles[iangle];
//         auto [matom, oatom, natom] = mon_atom;

//         dvec3 u = atom_coords_cart[matom] - atom_coords_cart[oatom];
//         dvec3 v = atom_coords_cart[natom] - atom_coords_cart[oatom];
//         double u_norm = arma::norm(u);
//         double v_norm = arma::norm(v);
//         u /= u_norm;
//         v /= v_norm;

//         dvec3 p1 = arma::dvec3({1, -1, 1});
//         dvec3 p2 = arma::dvec3({-1, 1, 1});
//         dvec3 w;
//         if (arma::dot(u, v) != 1)        
//             w = arma::cross(u, v);   
//         else if (arma::dot(u, v) == 1 and (arma::dot(u, p1) != 1 and arma::dot(v, p1) != 1))        
//             w = arma::cross(u, p1);        
//         else if (arma::dot(u, v) == 1 and (arma::dot(u, p1) == 1 and arma::dot(v, p1) == 1))        
//             w = arma::cross(u, p2);        
//         else
//             throw std::runtime_error("Something went wrong with finding the perpendicular vector, \
//                                       contact devs!\n");
//         w /= arma::norm(w);

//         dvec3 u_x_w = arma::cross(u, w);
//         dvec3 w_x_v = arma::cross(w, v);

//         size_t posm = 3 * matom;
//         size_t poso = 3 * oatom;
//         size_t posn = 3 * natom;
//         bmatrix(offset + iangle, posm) = u_x_w(0) / u_norm;
//         bmatrix(offset + iangle, posm + 1) = u_x_w(1) / u_norm;
//         bmatrix(offset + iangle, posm + 2) = u_x_w(2) / u_norm;
//         bmatrix(offset + iangle, poso) = -1 * (w_x_v(0) / v_norm + u_x_w(0) / u_norm);
//         bmatrix(offset + iangle, poso + 1) = -1 * (w_x_v(1) / v_norm + u_x_w(1) / u_norm);
//         bmatrix(offset + iangle, poso + 2) = -1 * (w_x_v(2) / v_norm + u_x_w(2) / u_norm);
//         bmatrix(offset + iangle, posn) = w_x_v(0) / v_norm;
//         bmatrix(offset + iangle, posn + 1) = w_x_v(1) / v_norm;
//         bmatrix(offset + iangle, posn + 2) = w_x_v(2) / v_norm;
//     }

//     offset += red_int_coords.angles.size();
//     for (size_t idihedral = 0; idihedral < red_int_coords.dihedrals.size(); idihedral++)
//     {
//         auto [mopn_atom, dihedral_angle] = red_int_coords.dihedrals[idihedral];
//         auto [matom, oatom, patom, natom] = mopn_atom;

//         dvec3 u = atom_coords_cart[matom] - atom_coords_cart[oatom];
//         dvec3 w = atom_coords_cart[patom] - atom_coords_cart[oatom];
//         dvec3 v = atom_coords_cart[natom] - atom_coords_cart[patom];
//         double norm_u = arma::norm(u);
//         double norm_w = arma::norm(w);
//         double norm_v = arma::norm(v);
//         u /= norm_u;
//         w /= norm_w;
//         v /= norm_v;

//         dvec3 u_x_w = arma::cross(u, w);
//         dvec3 v_x_w = arma::cross(v, w);
//         double cos_phiu = arma::dot(u, w);
//         double cos_phiv = -1 * arma::dot(v, w);
//         double sin_phiu = sqrt(1 - pow(cos_phiu, 2));
//         double sin_phiv = sqrt(1 - pow(cos_phiv, 2));

//         double fac1 = 1.0 / (pow(sin_phiu, 2) * norm_u);
//         double fac2 = 1.0 / (pow(sin_phiv, 2) * norm_v);
//         double fac3 = cos_phiu / (norm_w * pow(sin_phiu, 2));
//         double fac4 = cos_phiv / (norm_w * pow(sin_phiv, 2));

//         size_t posm = 3 * matom;
//         size_t poso = 3 * oatom;
//         size_t posp = 3 * patom;
//         size_t posn = 3 * natom;
//         bmatrix(offset + idihedral, posm) = u_x_w(0) * fac1;
//         bmatrix(offset + idihedral, posm + 1) = u_x_w(1) * fac1;
//         bmatrix(offset + idihedral, posm + 2) = u_x_w(2) * fac1;
//         bmatrix(offset + idihedral, poso) = -u_x_w(0) * fac1 + (u_x_w(0) * fac3 - v_x_w(0) * fac4);
//         bmatrix(offset + idihedral, poso + 1) = -u_x_w(1) * fac1 + (u_x_w(1) * fac3 - v_x_w(1) * fac4);
//         bmatrix(offset + idihedral, poso + 2) = -u_x_w(2) * fac1 + (u_x_w(2) * fac3 - v_x_w(2) * fac4);
//         bmatrix(offset + idihedral, posp) = v_x_w(0) * fac2 - (u_x_w(0) * fac3 - v_x_w(0) * fac4);
//         bmatrix(offset + idihedral, posp + 1) = v_x_w(1) * fac2 - (u_x_w(1) * fac3 - v_x_w(1) * fac4);
//         bmatrix(offset + idihedral, posp + 2) = v_x_w(2) * fac2 - (u_x_w(2) * fac3 - v_x_w(2) * fac4);
//         bmatrix(offset + idihedral, posn) = -v_x_w(0) * fac2;
//         bmatrix(offset + idihedral, posn + 1) = -v_x_w(1) * fac2;
//         bmatrix(offset + idihedral, posn + 2) = -v_x_w(2) * fac2;
//     }

//     return bmatrix;
// }

// arma::sp_dmat Geometry::constructBMatrixSparse()
// {
//     arma::sp_dmat bmatrix;

//     return bmatrix;
// }