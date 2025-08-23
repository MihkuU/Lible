#pragma once

namespace lible
{
    namespace tests
    {
        bool overlap();

        bool overlapKernel();

        bool overlapD1Kernel();

        bool kineticEnergy();

        bool kineticEnergyKernel();

        bool kineticEnergyD1Kernel();

        bool nuclearAttraction();

        bool nuclearAttractionErf();

        bool externalCharges();

        bool externalChargesErf();

        bool externalChargesKernel();

        bool externalChargesErfKernel();

        bool externalChargesD1Kernel();

        bool externalChargesOperatorD1Kernel();

        bool potentialAtExternalChargesKernel();

        bool potentialAtExternalChargesErfKernel();

        bool dipoleMoment();

        bool dipoleMomentKernel();

        bool spinOrbitCoupling1El();

        bool spinOrbitCoupling1ElKernel();

        bool eri2Diagonal();

        bool eri2();

        bool eri4Diagonal();

        bool eri3();
        
        bool eri4();

        bool basisForAtom();

        bool basisForAtomAux();

        bool basisForAtoms();

        bool basisForAtomsAux();
        
        bool availableBasisSets();

        bool availableBasisSetsAux();

        bool sphericalTrafo();

        bool deployERI4Kernel();

        bool deployERI3Kernel();

        bool deployERI2Kernel();

        bool deployERI4D1Kernel();

        bool deployERI3D1Kernel();

        bool deployERI2D1Kernel();

        bool deployERI2D2Kernel();

        bool deployERI4SOCKernel();        

        bool deployERI3SOCKernel();

        bool ecoeffsRecurrence1();

        bool ecoeffsRecurrence2();

        bool ecoeffsRecurrence2_n1();

        bool ecoeffsPrimitive();

        bool ecoeffsPrimitivePair();

        bool ecoeffsPrimitivePair_n1();

        bool ecoeffsShell();

        bool calcBoysF();

        bool calcRInts3D();

        bool purePrimitiveNorm();

        bool numCartesians();

        bool numSphericals();

        bool numHermites();

        bool cartExps();

        bool getLPairsSymm();

        bool getLPairsNoSymm();
    }
}
