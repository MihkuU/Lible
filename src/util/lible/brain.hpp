#pragma once

#ifdef _USE_MPI_
#include <mpl/mpl.hpp>
#endif

namespace lible
{
    class Brain
    {
    public:
#ifdef _USE_MPI_
        static int returnTotalRank();

        static int returnTotalSize();

        static std::vector<int> returnNodeRanks();

        static mpl::communicator returnNodesComm();

        const static inline mpl::communicator &comm_world{mpl::environment::comm_world()};

        const static inline mpl::communicator &comm_nodes{returnNodesComm()};        
#endif
    };
}