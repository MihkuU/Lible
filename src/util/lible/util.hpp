#pragma once

/*
 * Header file for general utilities including MPI parallelization.
 * Please be nice and include this file only in the .cpp files.
 */

#include <iostream>
#include <omp.h>
#include <string>

#ifdef _USE_MPI_
#include <mpl/mpl.hpp>
#endif

namespace lible
{

// #ifdef _USE_MPI_
//     struct Para
//     {
//         /*
//          * A word is appropriate here.
//          *
//          */
//         const static inline mpl::communicator &comm_world{mpl::environment::comm_world()};
//         // const mpl::communicator &comm_local;
//     };
// #endif

    static void palPrint(const std::string &message, const bool shutup = false)
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

#ifdef _USE_MPI_
    // static int returnTotalRank(const MPI_Comm &comm)
    // {
    //     int rank;
    //     MPI_Comm_rank(comm, &rank);
    //     int num_threads = omp_get_num_threads();
    //     int thread_num = omp_get_thread_num();

    //     return num_threads * rank + thread_num;
    // }

    // static int returnTotalSize(const MPI_Comm &comm)
    // {
    //     int size;
    //     MPI_Comm_size(comm, &size);
    //     int num_threads = omp_get_num_threads();

    //     return num_threads * size;
    // }
#endif

}