#include <lible/brain.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

using LB = lible::Brain;

using std::pair;
using std::string;
using std::vector;

#ifdef _LIBLE_USE_MPI_

// Some code adopted from https://rabauke.github.io/mpl/html/examples/stl_container.html.
// Later put them in some utils file?

static pair<string, string> splitStrByLast(const string &str, const char &delimiter)
{
    size_t pos_last{};
    for (size_t i = 0; i < str.size(); i++)
        if (str[i] == delimiter)
            pos_last = i;

    string part1 = str.substr(0, pos_last);
    string part2 = str.substr(pos_last + 1, str.size() - (pos_last + 1));

    return {part1, part2};
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
    vector<int> node_ranks = returnNodeRanks();

    mpl::ranks node_ranks_mpl(node_ranks.size());
    for (size_t i = 0; i < node_ranks.size(); i++)
        node_ranks_mpl[i] = node_ranks[i];

    mpl::group group_world{comm_world};
    mpl::group group_nodes(mpl::group::include, group_world, node_ranks_mpl);

    return mpl::communicator(mpl::communicator::comm_collective, comm_world, group_nodes);
}

int LB::returnTotalRank()
{
    // TODO: documentation - to be called inside omp parallel regions
    int num_threads = omp_get_num_threads();
    int thread_num = omp_get_thread_num();
    int rank = comm_nodes.rank();

    // if (!omp_in_parallel())
    //     throw std::runtime_error("returnTotalRank() not called inside parallel region!");    

    return num_threads * rank + thread_num;
}

int LB::returnTotalSize()
{
    // TODO: documentation - to be called inside omp parallel regions
    int num_threads = omp_get_num_threads();

    // if (!omp_in_parallel())
    //     throw std::runtime_error("returnTotalSize() not called inside parallel region!");

    return num_threads * comm_nodes.size();
}

vector<double> LB::bcastVector(const int &root_rank,
                               const mpl::communicator &comm,
                               const vector<double> &in)
{
    size_t n_elem;
    if (comm.rank() == root_rank)
        n_elem = in.size();

    comm.bcast(root_rank, n_elem);

    vector<double> out;
    if (comm.rank() == root_rank)
        out = in;
    else
        out.resize(n_elem);

    mpl::contiguous_layout<double> layout(n_elem);
    comm.bcast(root_rank, out.data(), layout);

    return out;
}

// vector<string> LB::scatterStrings(const int &root_rank,
//                                   const mpl::communicator &comm,
//                                   const vector<string> &in)
// {
//     /*
//     TODO: move to gci_para.hpp
//     */

//     string chunk;    
//     vector<size_t> sizes;
//     if (root_rank == comm.rank())
//     {
//         vector<string> chunks(comm.size(), "");
//         for (size_t i = 0; i < in.size(); i++)
//             chunks[i % comm.size()] += in[i];

//         chunk = chunks[root_rank];

//         vector<vector<size_t>> sizes_per_proc(comm.size());
//         for (size_t i = 0; i < in.size(); i++)
//             sizes_per_proc[i % comm.size()].push_back(in[i].size());

//         sizes = sizes_per_proc[root_rank];

//         for (int i = 0; i < comm.size(); i++)
//             if (i != root_rank)
//             {
//                 isend(i, comm, chunks[i]);
//                 isend(i, comm, sizes_per_proc[i]);
//             }
//     }

//     if (root_rank != comm.rank())
//     {
//         chunk = irecv<string>(root_rank, comm);
//         sizes = irecv<vector<size_t>>(root_rank, comm);
//     }

//     vector<string> out(sizes.size());
//     for (size_t i = 0, idx = 0; i < sizes.size(); i++)
//     {
//         string item(sizes[i], ' ');
//         for (size_t j = 0; j < sizes[i]; j++, idx++)
//             item[j] = chunk[idx];

//         out[i] = item;
//     }

//     return out;
// }

#endif
