#pragma once

#include <armadillo>
#include <cassert>
#include <map>
#include <vector>
#include <tuple>

#include "geomopt.h"

namespace Lible
{
    namespace GeomOpt
    {
        class Geometry
        {
            /*
             * Class for representing the molecular geometry. Provides utilities for transformation
             * between cartesian and redundant internal coordinates.
             *
             * A lot of what is implemented here follows the paper https://doi.org/10.1063/1.1515483.
             *
             * NB! The coordinates are expected to be in Bohr!
             */
        public:
            Geometry(const std::vector<double> &coords_cart_, const std::vector<std::string> &atoms_);

            std::vector<double> transformCoordsCartToRedint(const std::vector<double> &coords_cart);
            std::vector<double> transformCoordsRedintToCart(const std::vector<double> &coords_redint);
            
            std::vector<double> transformGradCartToRedint(const std::vector<double> &grad_cart);
            std::vector<double> transformStepRedIntToCart(const std::vector<double> &coords_cart,
                                                          const std::vector<double> &coords_redint,
                                                          const std::vector<double> &step_redint);

            std::vector<double> getCoordsCart()
            {
                return coords_cart;
            }

            std::vector<double> getCoordsRedint()
            {
                return red_int_coords.coords_redint;
            }

            std::size_t getNumBonds() const
            {
                return red_int_coords.bonds.size();
            }

            std::size_t getNumAngles() const
            {
                return red_int_coords.angles.size();
            }

            std::size_t getNumDihedrals() const
            {
                return red_int_coords.dihedrals.size();
            }

            typedef std::tuple<std::size_t, std::size_t> doublet;
            typedef std::tuple<std::size_t, std::size_t, std::size_t> triplet;
            typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> quartet;
            struct RedIntCoords
            {
                /*
                 * Implementation of redundant internal coordinates as defined in:
                 * https://doi.org/10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0.
                 */
                std::vector<std::pair<doublet, double>> bonds;
                std::vector<std::pair<triplet, double>> angles;
                std::vector<std::pair<quartet, double>> dihedrals;
                std::vector<double> coords_redint;
            };

        private:
            typedef arma::vec::fixed<3> dvec3;

            std::size_t n_atoms;
            std::vector<double> coords_cart;
            std::vector<std::string> atoms;
            std::vector<dvec3> atom_coords_cart;
            std::vector<std::vector<std::size_t>> atom_bonding_partners;

            // arma::dmat b_matrix;
            // arma::dmat b_matrix_pinv;
            // arma::dmat b_matrix_t;

            RedIntCoords red_int_coords;

            double calcAngle(const std::size_t &matom, const std::size_t &oatom, const std::size_t &natom,
                             const std::vector<dvec3> &atom_coords_cart);
            double calcDihedral(const std::size_t &matom, const std::size_t &oatom, const std::size_t &patom, const std::size_t &natom,
                                const std::vector<dvec3> &atom_coords_cart);
            double calcDistance(const std::size_t &matom, const std::size_t &natom,
                                const std::vector<dvec3> &atom_coords_cart);
            // std::size_t findClosestAtom(const std::size_t &iatom);
            std::vector<std::vector<std::size_t>> findAtomBondingPartners(const std::vector<dvec3> &atom_coords_cart);
            arma::dmat constructBMatrix(const RedIntCoords &red_int_coords);
            arma::sp_dmat constructBMatrixSparse();

            RedIntCoords constructRedIntCoords(const std::vector<dvec3> &atom_coords_cart,
                                               const std::vector<std::vector<std::size_t>> &atom_bonding_partners);
        };
    }
}