#pragma once

namespace lible
{
	namespace tests
	{
		bool numCartesians();	

		bool numSphericals();

		bool eri2Diagonal();

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
	}
}
