#pragma once

#include <set>
#include <string>
#include <map>
#include <vector>
#include <unordered_map>

#ifdef LIBLE_MAIN_BASIS_DIR
#define path_to_basis_sets LIBLE_MAIN_BASIS_DIR
#endif 

#ifdef LIBLE_AUX_BASIS_DIR
#define path_to_aux_basis_sets LIBLE_AUX_BASIS_DIR
#endif

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

        class BasisPaths
        {
        public:
            static std::string getMainBasisSetsPath()
            {
                return main_basis_sets_path;
            }

            static std::string getAuxBasisSetsPath()
            {
                return aux_basis_sets_path;
            }

            static void setMainBasisSetsPath(const std::string &path)
            {
                main_basis_sets_path = path;
            }

            static void setAuxBasisSetsPath(const std::string &path)
            {
                aux_basis_sets_path = path;
            }

        private:
            /** Absolute path to the main basis sets. */
            static inline std::string main_basis_sets_path{path_to_basis_sets};

            /** Absolute path to the auxiliary basis sets. */
            static inline std::string aux_basis_sets_path{path_to_aux_basis_sets};
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