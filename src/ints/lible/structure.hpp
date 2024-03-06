#pragma once

#include <map>

#include <lible/shell.hpp>

namespace lible
{
    namespace ints
    {
        struct Structure
        {
            /*
             * Struct for representing the atomic structure. Note that the access is public for
             * several data members. What is meant to be used 'internally' is declared private.
             */

            Structure(const std::string &basis_set,                      
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coordinates_angstroem);

            int max_angular_momentum;
            std::size_t n_atomic_orbitals;
            std::map<int, std::vector<Shell>> shells;

        private:
            // TODO: in the future, different atoms must have the capability to have different basis.
            std::size_t n_atoms;
            std::string basis_set;
            std::vector<double> coordinates;
            std::vector<int> atomic_nrs;
            std::vector<std::string> elements;

            // std::string returnBasisPath(const std::string &basis_set);
            std::vector<Shell> parseBasisJSONFile(const std::string &basis_path); // TODO: remove
            void constructShells(std::map<int, std::vector<Shell>> &shells);
            void readBasis(const std::string &basis_set, std::map<int, std::vector<Shell>> &shells);
        };
    }
}