set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(testlible PRIVATE 
	${tests_dir}/main.cpp 
	${tests_dir}/tests.cpp)

target_link_libraries(testlible PRIVATE lible)
target_include_directories(testlible PRIVATE ${tests_dir})
target_include_directories(testlible PRIVATE ${src_dir}/src/)

set(Tests 
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
	ERI2Diagonal
	availableBasisSets
	availableBasisSetsAux
	purePrimitiveNorm
	numHermites
	cartExps
	getLPairsSymm
	getLPairsNoSymm)

foreach(item ${Tests})
	add_test(NAME "lible::ints::${item}" COMMAND ${CMAKE_BINARY_DIR}/testlible ${item})
endforeach()
