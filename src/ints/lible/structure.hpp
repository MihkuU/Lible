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

            Structure(const std::string &basis_set,
                      const std::vector<std::string> &elements,
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

            void constructShells(int &max_angular_momentum,
                                 size_t &n_atomic_orbitals,
                                 std::map<int, std::vector<Shell>> &shells);

            void readBasis(const std::string &basis_set,
                           int &max_angular_momentum,
                           size_t &n_atomic_orbitals,
                           std::map<int, std::vector<Shell>> &shells);
        };
    }
}