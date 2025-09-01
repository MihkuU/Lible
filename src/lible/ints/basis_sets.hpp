#pragma once

#include <set>
#include <string>
#include <map>
#include <vector>
#include <unordered_map>

namespace lible
{
    namespace ints
    {
        enum class AuxBasisFamily
        {
            ahlrichs_fit,
        };

        enum class BasisFamily
        {
            ahlrichs,
            ano,
            dunning,
            pople,
            sto,
        };

        /** */
        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtom(const int atomic_nr, const std::string &basis_set);

        std::map<int, std::vector<shell_exps_coeffs_t>>
        basisForAtomAux(const int atomic_nr, const std::string &aux_basis_set);

        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtoms(const std::set<int> &atomic_nrs, const std::string &basis_set);

        std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>>
        basisForAtomsAux(const std::set<int> &atomic_nrs, const std::string &aux_basis_set);

        std::set<std::string> availableBasisSets();
        
        std::set<std::string> availableBasisSetsAux();

        std::string returnAuxBasisFamilyString(const std::string &aux_basis_set);

        std::string returnBasisFamilyString(const std::string &basis_set);

        std::string returnAuxBasisPath(const std::string &aux_basis_set);

        std::string returnBasisPath(const std::string &basis_set);

        AuxBasisFamily returnAuxBasisFamily(const std::string &aux_basis_set);

        BasisFamily returnBasisFamily(const std::string &basis_set);
    }
}