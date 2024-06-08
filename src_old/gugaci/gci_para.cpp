#ifdef _USE_MPI_

#include <lible/gci_para.hpp>
#include <lible/brain.hpp>

#include <tuple>

namespace LG = lible::guga;

using LB = lible::Brain;

using std::map;
using std::set;
using std::string;
using std::tuple;
using std::vector;

using map1_t = map<size_t, set<string>>;
using map2_t = map<string, vector<int>>;

template <typename T>
static vector<T> returnReceived(const int &rank_dest, const mpl::communicator &comm,
                                const T &data)
{
    mpl::irequest r_send(comm.isend(data, rank_dest));

    vector<T> received;
    if (comm.rank() == rank_dest)
    {
        received.resize(comm.size());
        mpl::irequest_pool r_pool;
        for (int i = 0; i < comm.size(); i++)
            r_pool.push(comm.irecv(received[i], i));
        r_pool.waitall();
    }
    r_send.wait();

    return received;
}

auto serializeMap(const map1_t &in)
{
    vector<size_t> keys;
    vector<size_t> value_sizes;
    vector<size_t> item_sizes;
    string items = "";
    for (auto &[key, val] : in)
    {
        keys.push_back(key);
        value_sizes.push_back(val.size());
        for (auto &item : val)
        {
            items += item;
            item_sizes.push_back(item.size());
        }
    }

    keys.shrink_to_fit();
    value_sizes.shrink_to_fit();
    item_sizes.shrink_to_fit();
    items.shrink_to_fit();

    return std::make_tuple(keys, value_sizes, item_sizes, items);
}

static auto serializeMap(const map2_t &in)
{
    string keys = "";
    vector<size_t> key_sizes;
    vector<size_t> value_sizes;
    vector<int> items;

    for (auto &[key, val] : in)
    {
        keys += key;
        key_sizes.push_back(key.size());
        value_sizes.push_back(val.size());
        for (auto &item : val)
            items.push_back(item);
    }

    keys.shrink_to_fit();
    key_sizes.shrink_to_fit();
    value_sizes.shrink_to_fit();
    items.shrink_to_fit();

    return std::make_tuple(keys, key_sizes, value_sizes, items);
}

static auto deSerializeMap(const vector<size_t> &keys, const vector<size_t> &value_sizes,
                           const vector<size_t> &item_sizes, const string &items)
{
    map1_t out;

    size_t idx_item = 0, offset_item = 0;
    for (size_t ikey = 0; ikey < keys.size(); ikey++)
    {
        set<string> val;
        for (size_t i = 0; i < value_sizes[ikey]; i++)
        {
            size_t len = item_sizes[idx_item];

            string item = items.substr(offset_item, len);
            val.insert(item);

            idx_item++;
            offset_item += len;
        }
        size_t key = keys[ikey];
        out[key] = val;
    }

    return out;
}

static auto deSerializeMap(const string &keys, const vector<size_t> &key_sizes,
                           const vector<size_t> &value_sizes, const vector<int> &items)
{
    map2_t out;
    
    size_t idx_item = 0, offset_key = 0;
    for (size_t ikey = 0; ikey < key_sizes.size(); ikey++)
    {
        vector<int> val;
        for (size_t i = 0; i < value_sizes[ikey]; i++)
        {
            val.push_back(items[idx_item]);
            idx_item++;
        }
        string key = keys.substr(offset_key, key_sizes[ikey]);        
        offset_key += key_sizes[ikey];
        out[key] = val;
    }

    return out;
}

template <>
map1_t LG::allReduceMaps(const mpl::communicator &comm, const map1_t &in)
{
    auto [keys, value_sizes, item_sizes, items] = serializeMap(in);

    auto keys_recv = returnReceived(0, comm, keys);
    auto value_sizes_recv = returnReceived(0, comm, value_sizes);
    auto item_sizes_recv = returnReceived(0, comm, item_sizes);
    auto items_recv = returnReceived(0, comm, items);

    map1_t reduced;
    if (comm.rank() == 0)
    {
        for (int iproc = 0; iproc < comm.size(); iproc++)
        {
            map1_t reduced_proc = deSerializeMap(keys_recv[iproc], value_sizes_recv[iproc],
                                                 item_sizes_recv[iproc], items_recv[iproc]);

            for (auto &[key, val] : reduced_proc)
            {
                if (reduced.find(key) == reduced.end())
                    reduced[key] = val;
                else
                    reduced[key].insert(val.begin(), val.end());
            }
        }
    }

    auto [keys_reduced, value_sizes_reduced, item_sizes_reduced, items_reduced] =
        serializeMap(reduced);

    if (comm.rank() == 0)
    {
        for (int rank = 0; rank < comm.size(); rank++)
            if (rank != 0)
            {
                comm.send(keys_reduced, rank, mpl::tag_t(0));
                comm.send(value_sizes_reduced, rank, mpl::tag_t(1));
                comm.send(item_sizes_reduced, rank, mpl::tag_t(2));
                comm.send(items_reduced, rank, mpl::tag_t(3));
            }
    }
    else
    {
        comm.recv(keys_reduced, 0, mpl::tag_t(0));
        comm.recv(value_sizes_reduced, 0, mpl::tag_t(1));
        comm.recv(item_sizes_reduced, 0, mpl::tag_t(2));
        comm.recv(items_reduced, 0, mpl::tag_t(3));
    }

    map1_t out = deSerializeMap(keys_reduced, value_sizes_reduced,
                                item_sizes_reduced, items_reduced);

    return out;
}

#include <iostream> // TMP

template <>
map2_t LG::allReduceMaps(const mpl::communicator &comm, const map2_t &in)
{
    // for (auto &[key, val] : in)    
    // {
    //     std::cout << "key = " << key << "\n";
    //     for (auto &item : val)
    //         std::cout << "item = " << item << "\n";
    // }

    auto [keys, key_sizes, value_sizes, items] = serializeMap(in);

    auto keys_recv = returnReceived(0, comm, keys);
    auto key_sizes_recv = returnReceived(0, comm, key_sizes);
    auto value_sizes_recv = returnReceived(0, comm, value_sizes);
    auto items_recv = returnReceived(0, comm, items);

    // for (size_t i = 0; i < keys.size(); i++)
    // {
    //     if (keys[i] != keys_recv[0][i])
    //         printf("\nputsis\n");
    // }

    // for (size_t i = 0; i < key_sizes.size(); i++)
    // {
    //     if (key_sizes[i] != key_sizes_recv[0][i])
    //         printf("\nputsis2\n");
    // }

    // for (size_t i = 0; i < value_sizes.size(); i++)
    // {
    //     if (value_sizes[i] != value_sizes_recv[0][i])
    //         printf("\nputsis3\n");
    // }

    // for (size_t i = 0; i < value_sizes.size(); i++)
    // {
    //     if (items[i] != items_recv[0][i])
    //         printf("\nputsis4\n");
    // }    

    map2_t reduced;
    if (comm.rank() == 0)
    {
        for (int iproc = 0; iproc < comm.size(); iproc++)
        {
            map2_t reduced_proc = deSerializeMap(keys_recv[iproc], key_sizes_recv[iproc],
                                                 value_sizes_recv[iproc], items_recv[iproc]);

            for (auto &[key, val] : reduced_proc)
            {
                if (reduced.find(key) == reduced.end())
                    reduced[key] = val;
                else
                    reduced[key].insert(reduced[key].end(), val.begin(), val.end());
            }
        }
    }

    auto [keys_reduced, key_sizes_reduced, value_sizes_reduced, items_reduced] =
        serializeMap(reduced);

    if (comm.rank() == 0)
    {
        for (int rank = 0; rank < comm.size(); rank++)
            if (rank != 0)
            {
                comm.send(keys_reduced, rank, mpl::tag_t(0));
                comm.send(key_sizes_reduced, rank, mpl::tag_t(1));
                comm.send(value_sizes_reduced, rank, mpl::tag_t(2));
                comm.send(items_reduced, rank, mpl::tag_t(3));
            }
    }
    else
    {
        comm.recv(keys_reduced, 0, mpl::tag_t(0));
        comm.recv(key_sizes_reduced, 0, mpl::tag_t(1));
        comm.recv(value_sizes_reduced, 0, mpl::tag_t(2));
        comm.recv(items_reduced, 0, mpl::tag_t(3));
    }

    map2_t out = deSerializeMap(keys_reduced, key_sizes_reduced,
                                value_sizes_reduced, items_reduced);

    return out;
}

vector<string> LG::scatterStrings(const int &root_rank,
                                  const mpl::communicator &comm,
                                  const vector<string> &in)
{
    string chunk;    
    vector<size_t> sizes;
    if (root_rank == comm.rank())
    {
        vector<string> chunks(comm.size(), "");
        for (size_t i = 0; i < in.size(); i++)
            chunks[i % comm.size()] += in[i];

        chunk = chunks[root_rank];

        vector<vector<size_t>> sizes_per_proc(comm.size());
        for (size_t i = 0; i < in.size(); i++)
            sizes_per_proc[i % comm.size()].push_back(in[i].size());

        sizes = sizes_per_proc[root_rank];

        for (int i = 0; i < comm.size(); i++)
            if (i != root_rank)
            {
                LB::isend(i, comm, chunks[i]);
                LB::isend(i, comm, sizes_per_proc[i]);
            }
    }

    if (root_rank != comm.rank())
    {
        chunk = LB::irecv<string>(root_rank, comm);
        sizes = LB::irecv<vector<size_t>>(root_rank, comm);
    }

    vector<string> out(sizes.size());
    for (size_t i = 0, idx = 0; i < sizes.size(); i++)
    {
        string item(sizes[i], ' ');
        for (size_t j = 0; j < sizes[i]; j++, idx++)
            item[j] = chunk[idx];

        out[i] = item;
    }

    return out;
}

#endif