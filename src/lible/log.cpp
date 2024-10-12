#include <lible/log.hpp>

#ifdef _LIBLE_USE_MPI_
#include <lible/brain.hpp>
#endif

#include <filesystem>

#include <fmt/core.h>
#include <omp.h>

namespace LL = lible::log;

using std::string;

LL::Logger::Logger()
{
#ifdef _LIBLE_USE_MPI_
    int rank = Brain::comm_world.rank();
    log_fname = fmt::format("lible.mpi{}.log", rank);
#else
    log_fname = "lible.log";
#endif

    if (std::filesystem::exists(log_fname))
        std::filesystem::remove(log_fname);

    log_file.open(log_fname);
}

LL::Logger::~Logger()
{
    log_file.close();
}

void LL::Logger::operator<<(const std::string &message)
{
    if (omp_in_parallel())
    {
#pragma omp parallel
        {
#pragma omp single
            log_file << message << std::flush;
        }
    }
    else
        log_file << message << std::flush;
}