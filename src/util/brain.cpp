#include <lible/brain.hpp>

#include <omp.h>
#include <string>
#include <utility>
#include <vector>

using LB = lible::Brain;

using std::pair;
using std::string;
using std::vector;

#ifdef _USE_MPI_

static pair<string, string> splitStrByLast(const string &str, const char &delimiter)
{
    size_t pos_last;
    for (size_t i = 0; i < str.size(); i++)
        if (str[i] == delimiter)
            pos_last = i;

    string part1 = str.substr(0, pos_last);
    string part2 = str.substr(pos_last + 1, str.size() - (pos_last + 1));

    return {part1, part2};
}

int LB::returnTotalRank()
{
    int num_threads = 1, thread_num = 1;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        thread_num = omp_get_thread_num();
    }

    int rank = comm_nodes.rank();

    return num_threads * rank + thread_num;
}

int LB::returnTotalSize()
{
    int num_threads = 1;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    return num_threads * comm_nodes.size();
}

vector<int> LB::returnNodeRanks()
{
    int rank = comm_world.rank();
    string name = mpl::environment::processor_name();
    string name_rank = name + "." + std::to_string(rank);

    int n_nodes;

    std::map<string, vector<int>> nodes_ranks;
    mpl::irequest r_send(comm_world.isend(name_rank, 0));
    if (rank == 0)
    {
        vector<string> names_ranks(comm_world.size());

        mpl::irequest_pool r_pool;
        for (int i = 0; i < comm_world.size(); i++)
            r_pool.push(comm_world.irecv(names_ranks[i], i));
        r_pool.waitall();

        for (int i = 0; i < comm_world.size(); i++)
        {
            auto [node, rank] = splitStrByLast(names_ranks[i], '.');
            nodes_ranks[node].push_back(stoi(rank));
        }

        n_nodes = nodes_ranks.size();
    }
    r_send.wait();

    comm_world.bcast(0, n_nodes);

    vector<int> node_ranks(n_nodes);
    if (rank == 0)
    {
        int idx = 0;
        for (auto &[node, ranks] : nodes_ranks)
        {
            node_ranks[idx] = *ranks.begin();
            idx++;
        }
    }

    comm_world.bcast(0, node_ranks.data(), mpl::contiguous_layout<int>(n_nodes));

    return node_ranks;
}

mpl::communicator LB::returnNodesComm()
{
    int rank = comm_world.rank();

    vector<int> node_ranks = returnNodeRanks();

    bool is_node_rank = (std::find(node_ranks.begin(), node_ranks.end(), rank) != node_ranks.end())
                            ? true
                            : false;

    mpl::ranks node_ranks_mpl(node_ranks.size());
    for (size_t i = 0; i < node_ranks.size(); i++)
        node_ranks_mpl[i] = node_ranks[i];

    mpl::group group_world{comm_world};
    mpl::group group_nodes(mpl::group::include, group_world, node_ranks_mpl);

    return mpl::communicator(mpl::communicator::comm_collective, comm_world, group_nodes);
}

#endif