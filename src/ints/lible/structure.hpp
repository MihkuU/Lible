#pragma once

#include <map>
#include <string>

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

            int getMaxL() const;

            int getZ(const size_t iatom) const;

            size_t getDimAO() const;

            size_t getNAtoms() const;

            std::array<double, 3> getCoordsAtom(const size_t iatom) const;

            std::vector<Shell> getShellsL(const int l) const;

            std::map<int, std::vector<Shell>> getShells() const;

        private:
            // TODO: in the future, different atoms must have the capability to have different basis.
            int max_l;
            size_t dim_ao;
            size_t n_atoms;
            std::string basis_set;
            std::vector<double> coordinates;            
            std::vector<int> atomic_nrs;
            std::vector<std::array<double, 3>> coordinates_xyz;
            std::vector<std::string> elements;
            std::map<int, std::vector<Shell>> shells;

            void constructShells(int &max_l, size_t &dim_ao,
                                 std::map<int, std::vector<Shell>> &shells);

            void readBasis(const std::string &basis_set, int &max_l, size_t &dim_ao,
                           std::map<int, std::vector<Shell>> &shells);
        };
    }
}