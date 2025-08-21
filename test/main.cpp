#include <format>
#include <stdexcept>

#include <tests.hpp>

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::string msg =
            std::format("Invalid number of arguments: {}. Only one argument must be given!", (argc - 1));
        throw std::runtime_error(msg);
    }

    std::string test_name = argv[1];

    bool success = false;
    if (test_name == "numCartesians")
        success = lible::tests::numCartesians();
    else if (test_name == "numSphericals")
        success = lible::tests::numSphericals();
    else if (test_name == "overlap")
        success = lible::tests::overlap();
    else if (test_name == "overlapKernel")
        success = lible::tests::overlapKernel();
    else if (test_name == "overlapD1Kernel")
        success = lible::tests::overlapD1Kernel();
    else if (test_name == "kineticEnergy")
        success = lible::tests::kineticEnergy();
    else if (test_name == "kineticEnergyKernel")
        success = lible::tests::kineticEnergyKernel();
    else if (test_name == "kineticEnergyD1Kernel")
        success = lible::tests::kineticEnergyD1Kernel();
    else if (test_name == "nuclearAttraction")
        success = lible::tests::nuclearAttraction();
    else if (test_name == "nuclearAttractionErf")
        success = lible::tests::nuclearAttractionErf();
    else if (test_name == "externalCharges")
        success = lible::tests::externalCharges();
    else if (test_name == "externalChargesErf")
        success = lible::tests::externalChargesErf();
    else if (test_name == "externalChargesKernel")
        success = lible::tests::externalChargesKernel();
    else if (test_name == "externalChargesErfKernel")
        success = lible::tests::externalChargesErfKernel();
    else if (test_name == "externalChargesD1Kernel")
        success = lible::tests::externalChargesD1Kernel();
    else if (test_name == "externalChargesOperatorD1Kernel")
        success = lible::tests::externalChargesOperatorD1Kernel();
    else if (test_name == "potentialAtExternalChargesKernel")
        success = lible::tests::potentialAtExternalChargesKernel();
    else if (test_name == "potentialAtExternalChargesErfKernel")
        success = lible::tests::potentialAtExternalChargesErfKernel();
    else if (test_name == "ERI2Diagonal")
        success = lible::tests::eri2Diagonal();
    else if (test_name == "availableBasisSets")
        success = lible::tests::availableBasisSets();
    else if (test_name == "availableBasisSetsAux")
        success = lible::tests::availableBasisSetsAux();
    else if (test_name == "purePrimitiveNorm")
        success = lible::tests::purePrimitiveNorm();
    else if (test_name == "numHermites")
        success = lible::tests::numHermites();
    else if (test_name == "cartExps")
        success = lible::tests::cartExps();
    else if (test_name == "getLPairsSymm")
        success = lible::tests::getLPairsSymm();
    else if (test_name == "getLPairsNoSymm")
        success = lible::tests::getLPairsNoSymm();
    else
        throw std::runtime_error(std::format("Invalid test name specified: {}", test_name));

    if (success)
        return 0;
    else
        return 1;
}
