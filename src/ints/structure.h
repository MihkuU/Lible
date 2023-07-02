#pragma once

#include <map>
#include "ints.h"
#include "shell.h"

namespace Lible
{
    /*
     * Struct for representing the atomic structure. Used 'struct' here because its similar to structure...
     * Note that the access is public for several data members. What is meant to be used 'internally' is declared private.
     */
    struct Lible::Ints::Structure
    {
        Structure(const std::string &basis_set, const std::vector<double> &coordinates, const std::vector<std::string> &elements);

        int max_angular_momentum;
        std::map<std::pair<int, int>, std::vector<Shells::ShellPair>> shell_pairs;

    private:
        std::size_t n_atoms;
        std::string basis_set;
        std::vector<double> coordinates;
        std::vector<std::string> elements;
        std::map<int, std::vector<Shells::Shell>> shells_per_angmom;

        std::string returnBasisPath(const std::string &basis_set);
        std::vector<Shells::Shell> parseBasisJSONFile(const std::string &basis_path);
        void readBasis(const std::string &basis_set);
    };
}