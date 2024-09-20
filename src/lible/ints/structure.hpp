#pragma once

#include <map>
#include <string>

#include <lible/ints/shell.hpp>

namespace lible
{
    namespace ints
    {
        /** */
        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        /** */
        typedef std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>> basis_atoms_t;

        /**
         * Struct for representing the atomic structure. Note that the access is public for
         * several data members. What is meant to be used 'internally' is declared private.
         *
         * TODO: in the future, different atoms must have the capability to have different basis.
         */
        struct Structure
        {
            Structure(const std::string &basis_set, const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            Structure(const std::string &basis_set, const std::vector<std::string> &elements,
                      const std::vector<double> &coords_angstrom);

            Structure(const std::string &basis_set, const std::string &aux_basis_set,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            Structure(const std::string &basis_set, const std::string &aux_basis_set,
                      const std::vector<std::string> &elements,
                      const std::vector<double> &coords_angstrom);

            bool getUseRI() const;

            int getMaxL() const;

            int getMaxLAux() const;

            int getZ(const size_t iatom) const;

            size_t getDimAO() const;

            size_t getDimAOAux() const;

            size_t getDimAOCart() const;

            size_t getDimAOCartAux() const;

            size_t getNAtoms() const;

            std::array<double, 3> getCoordsAtom(const size_t iatom) const;

            std::vector<Shell> getShellsL(const int l) const;

            std::vector<Shell> getShellsLAux(const int l) const;

            std::map<int, std::vector<Shell>> getShells() const;

            std::map<int, std::vector<Shell>> getShellsAux() const;

        private:
            bool use_ri = false;
            int max_l;            
            int max_l_aux;
            size_t dim_ao;
            size_t dim_ao_aux;
            size_t dim_ao_cart;
            size_t dim_ao_cart_aux;
            size_t n_atoms;
            std::string aux_basis_set;
            std::string basis_set;
            std::vector<double> coords;
            std::vector<int> atomic_nrs;
            std::vector<std::array<double, 3>> coords_xyz;
            std::vector<std::string> elements;
            std::map<int, std::vector<Shell>> shells;
            std::map<int, std::vector<Shell>> shells_aux;

            void constructShells(const basis_atoms_t &basis_atoms, int &max_l, size_t &dim_ao,
                                 size_t &dim_ao_cart, std::map<int, std::vector<Shell>> &shells);
        };
    }
}