set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(testlible PRIVATE
        ${tests_dir}/main.cpp
        ${tests_dir}/tests_ints.cpp
        ${tests_dir}/tests_solver.cpp)

target_compile_options(testlible PRIVATE -Wall -Wno-unused-local-typedefs)
target_link_libraries(testlible PRIVATE lible)
target_include_directories(testlible PRIVATE ${tests_dir})
target_include_directories(testlible PRIVATE ${src_dir}/src/)

target_link_libraries(testlible PRIVATE
    "$<TARGET_PROPERTY:lible,LINK_LIBRARIES>")
target_include_directories(testlible PRIVATE
    "$<TARGET_PROPERTY:lible,INCLUDE_DIRECTORIES>")

## Armadillo
#target_compile_definitions(testlible PRIVATE ARMA_DONT_USE_WRAPPER)
#target_compile_definitions(testlible PRIVATE ARMA_DONT_USE_LAPACK)
#target_compile_definitions(testlible PRIVATE ARMA_DONT_USE_BLAS)
#target_include_directories(testlible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/include)
#if(EXISTS "${CMAKE_BINARY_DIR}/extern/arma/build/lib")
#    target_link_libraries(testlible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/lib/libarmadillo.a)
#elseif(EXISTS "${CMAKE_BINARY_DIR}/extern/arma/build/lib64")
#    target_link_libraries(testlible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/lib64/libarmadillo.a)
#endif()
#
#if(CMAKE_BUILD_TYPE STREQUAL "Release")
#    target_compile_definitions(testlible PRIVATE ARMA_NO_DEBUG)
#endif()

# Tests
if (LIBLE_BUILD_INTS)
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
            getLPairsNoSymm
            structureGhost
            structureGhostRI
            structureGhostCustom
            structureGhostCustomRI
    )

    foreach (item ${IntsTests})
        add_test(NAME "lible::ints::${item}" COMMAND ${CMAKE_BINARY_DIR}/testlible ${item})
    endforeach ()
endif ()

if (LIBLE_BUILD_SOLVER)
    set(SolverTests
            preconditionedCG
    )

    foreach (item ${SolverTests})
        add_test(NAME "lible::solver::${item}" COMMAND ${CMAKE_BINARY_DIR}/testlible ${item})
    endforeach ()
endif ()
