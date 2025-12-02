#pragma once

#include <string>

namespace lible::ints
{
    /// Family of an auxiliary basis set. Used for parsing the JSON file from
    /// https://www.basissetexchange.org/.
    enum class AuxBasisFamily
    {
        ahlrichs_fit,
    };

    /// Family of the main basis set. Used for parsing the JSON file from
    /// https://www.basissetexchange.org/.
    enum class BasisFamily
    {
        ahlrichs,
        ano,
        dunning,
        pople,
        sto,
    };

    /// Returns the main basis set family as a string.
    std::string basisFamilyString(const std::string &basis_set);

    /// Returns the auxiliary basis set family as a string.
    std::string auxBasisFamilyString(const std::string &aux_basis_set);

    /// Returns the path to the main basis set.
    std::string basisPath(const std::string &basis_set);

    /// Returns the path to the auxiliary basis set.
    std::string auxBasisPath(const std::string &aux_basis_set);

    /// Returns the family of the main basis set.
    BasisFamily basisFamily(const std::string &basis_set);

    /// Returns the family of the auxiliary basis set.
    AuxBasisFamily auxBasisFamily(const std::string &aux_basis_set);
}
