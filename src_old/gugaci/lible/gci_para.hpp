#pragma once

#ifdef _USE_MPI_

#include <map>
#include <set>
#include <string>
#include <vector>

#include <mpl/mpl.hpp>

namespace lible
{
    namespace guga
    {
        template <typename T, typename U>
        std::map<T, U> allReduceMaps(const mpl::communicator &comm,
                                     const std::map<T, U> &in);

        std::vector<std::string> scatterStrings(const int &root_rank,
                                                const mpl::communicator &comm,
                                                const std::vector<std::string> &in);
    }
}

#endif