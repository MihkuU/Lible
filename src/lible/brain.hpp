#pragma once

#include <vector>

#ifdef _LIBLE_USE_MPI_
#include <mpl/mpl.hpp>
#endif

namespace lible
{
    // TODO: replace the class with a namespace!
    class Brain
    {
    public:
#ifdef _LIBLE_USE_MPI_
        static std::vector<int> returnNodeRanks();

        static mpl::communicator returnNodesComm();

        const static inline mpl::communicator &comm_world{mpl::environment::comm_world()};

        const static inline mpl::communicator &comm_nodes{returnNodesComm()};

        static int returnTotalRank();

        static int returnTotalSize();

        static std::vector<double> bcastVector(const int &root_rank,
                                               const mpl::communicator &comm,
                                               const std::vector<double> &in);

        /*
         * Some code adopted from https://rabauke.github.io/mpl/html/examples/stl_container.html.         
         */

        // send an STL container
        template <typename T>
        static void isend(const int &rank, const mpl::communicator &comm, const T &x)
        {
            mpl::irequest r{comm.isend(x, rank)};
            r.wait();
        }

        // receive an STL container
        template <typename T>
        static T irecv(const int &rank, const mpl::communicator &comm)
        {
            T x;
            mpl::irequest r{comm.irecv(x, 0)};
            mpl::status_t s{r.wait()};
            return x;
        }

#endif
    };
}