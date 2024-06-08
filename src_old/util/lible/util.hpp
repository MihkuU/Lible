#pragma once

/*
 * Header file for general utilities including MPI parallelization.
 * Please be nice and include this file only in the .cpp files.
 */

#include <iostream>
#include <omp.h> //TODO: move to cpp-file
#include <string>

#ifdef _USE_MPI_
#include <mpl/mpl.hpp> //TODO: move to cpp-file
#endif

namespace lible
{
    // TODO: move impl to cpp-file, can move mpl also to cpp then -> good.
    inline void palPrint(const std::string &message, bool shutup = false)
    {
        if (!shutup)
        {
#ifdef _USE_MPI_
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
                std::cout << message;
#else
            std::cout << message;
#endif
        }
    }

    // TODO: overload a version with indent/padding given before message
}