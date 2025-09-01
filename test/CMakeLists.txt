set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(testlible PRIVATE 
	${tests_dir}/main.cpp 
	${tests_dir}/tests.cpp)

target_link_libraries(testlible PRIVATE lible)
target_include_directories(testlible PRIVATE ${tests_dir})
target_include_directories(testlible PRIVATE ${src_dir}/src/)

if(LIBLE_BUILD_INTS)
set(IntsTests
	numCartesians 
	numSphericals 
	overlap 
	overlapKernel 
	overlapD1Kernel 
	kineticEnergy
	kineticEnergyKernel
	kineticEnergyD1Kernel
	nuclearAttraction
	nuclearAttractionErf
	externalCharges
	externalChargesErf
	externalChargesKernel
	externalChargesErfKernel
	externalChargesD1Kernel
	externalChargesOperatorD1Kernel
	potentialAtExternalChargesKernel
	potentialAtExternalChargesErfKernel
	dipoleMoment
	dipoleMomentKernel
	spinOrbitCoupling1El
	spinOrbitCoupling1ElKernel
	eri2Diagonal
	eri2
	eri4Diagonal
	eri3
	eri4
	availableBasisSets
	availableBasisSetsAux
	basisForAtom
	basisForAtomAux
	basisForAtoms
	basisForAtomsAux
	deployERI4Kernel
	deployERI3Kernel
	deployERI2Kernel
	deployERI4D1Kernel
	deployERI3D1Kernel
	deployERI2D1Kernel
	deployERI2D2Kernel
	deployERI4SOCKernel
	deployERI3SOCKernel
	ecoeffsRecurrence1
	ecoeffsRecurrence2
	ecoeffsRecurrence2_n1
	ecoeffsPrimitive
	ecoeffsPrimitivePair
	ecoeffsPrimitivePair_n1
	ecoeffsShell
	calcBoysF
	calcRInts3D
	sphericalTrafo
	purePrimitiveNorm
	numHermites
	cartExps
	getLPairsSymm
	getLPairsNoSymm)

foreach(item ${IntsTests})
	add_test(NAME "lible::ints::${item}" COMMAND ${CMAKE_BINARY_DIR}/testlible ${item})
endforeach()
endif()
