#pragma once

#include <armadillo>
#include <cassert>
#include <map>
#include <vector>
#include <tuple>
#include <set>

#include "geomopt.h"

namespace Lible
{
    // class GeomOpt::Geometry
    namespace GeomOpt
    {
        class Geometry
        {
        public:
            Geometry(const std::vector<double> &coords_cart_, const std::vector<std::string> &atoms_);

            std::vector<double> getCoordsCart()
            {
                return coords_cart;
            }

        private:
            std::size_t n_atoms;
            std::vector<double> coords_cart;
            std::vector<std::string> atoms;
            std::vector<arma::dvec> atom_coords_cart;

            typedef std::tuple<std::size_t, std::size_t> doublet;
            typedef std::tuple<std::size_t, std::size_t, std::size_t> triplet;
            typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> quartet;
            /*
             * Implementation of redundant internal coordinates as defined in:
             * https://doi.org/10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0.
             */
            struct RedundantInternalCoordinates
            {
                std::vector<std::pair<doublet, double>> bonds;
                std::vector<std::pair<triplet, double>> angles;
                std::vector<std::pair<quartet, double>> dihedrals;
            };

            RedundantInternalCoordinates red_int_coords;

            double calcAngle(const std::size_t &iatom, const std::size_t &jatom, const std::size_t &katom);
            double calcDihedral(const std::size_t &iatom, const std::size_t &jatom, const std::size_t &katom, const std::size_t &latom);
            double calcDistance(const std::size_t &iatom, const std::size_t &jatom);

            size_t findClosestAtom(const std::size_t &iatom);
            std::vector<std::set<std::size_t>> findAtomBondingPartners();
            RedundantInternalCoordinates constructRedIntCoords();
            // arma::dmat constructBMatrix();

            void constructBMatrix(arma::dmat &bmatrix);
            void constructBMatrix(arma::sp_dmat &bmatrix);
        };
    }
}