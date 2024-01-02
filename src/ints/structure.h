#pragma once

#include <map>

#include "shell.h"

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
                      const std::vector<double> &coordinates_angstroem,
                      const std::vector<std::string> &elements);

            int max_angular_momentum;      // TODO: make const?
            std::size_t n_atomic_orbitals; // TODO: make const?
            std::map<int, std::vector<Shells::Shell>> shells;
            std::map<std::pair<int, int>, std::vector<Shells::ShellPair>> shell_pairs; // TODO: make const?

        private:
            std::size_t n_atoms;
            std::string basis_set;
            std::vector<double> coordinates;
            std::vector<std::string> elements;

            std::string returnBasisPath(const std::string &basis_set);
            std::vector<Shells::Shell> parseBasisJSONFile(const std::string &basis_path);
            void readBasis(const std::string &basis_set);
        };
    }
}