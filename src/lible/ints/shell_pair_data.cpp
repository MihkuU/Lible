#include <lible/ints/shell_pair_data.hpp>

namespace LI = lible::ints;

LI::ShellPairData::ShellPairData()
{
}

LI::ShellPairData::ShellPairData(const int la, const int lb, const Structure &structure)
    : la(la), lb(lb), structure(&structure)
{
    if (la == lb)
    {
        auto shells_a = structure.getShellsL(la);

        size_t size_a = shells_a.size();

        n_pairs = size_a * (size_a + 1) / 2;

        offsets.reserve(n_pairs);
        coords.reserve(n_pairs);
        ccoeffs.reserve(n_pairs);
        exps.reserve(n_pairs);

        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell <= ishell; jshell++)
            {
                offsets.push_back(std::make_pair(shells_a[ishell].pos,
                                                 shells_a[jshell].pos));

                coords.push_back(std::make_pair(shells_a[ishell].xyz_coords,
                                                shells_a[jshell].xyz_coords));

                ccoeffs.push_back(std::make_pair(shells_a[ishell].coeffs,
                                                 shells_a[jshell].coeffs));

                exps.push_back(std::make_pair(shells_a[ishell].exps,
                                              shells_a[jshell].exps));

                norms.push_back(std::make_pair(shells_a[ishell].norms,
                                               shells_a[jshell].norms));
            }
    }
    else
    {
        auto shells_a = structure.getShellsL(la);
        auto shells_b = structure.getShellsL(lb);

        size_t size_a = shells_a.size();
        size_t size_b = shells_b.size();

        n_pairs = size_a * size_b;

        offsets.reserve(n_pairs);
        coords.reserve(n_pairs);
        ccoeffs.reserve(n_pairs);
        exps.reserve(n_pairs);

        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell < size_b; jshell++)
            {
                offsets.push_back(std::make_pair(shells_a[ishell].pos,
                                                 shells_b[jshell].pos));

                coords.push_back(std::make_pair(shells_a[ishell].xyz_coords,
                                                shells_b[jshell].xyz_coords));

                ccoeffs.push_back(std::make_pair(shells_a[ishell].coeffs,
                                                 shells_b[jshell].coeffs));

                exps.push_back(std::make_pair(shells_a[ishell].exps,
                                              shells_b[jshell].exps));

                norms.push_back(std::make_pair(shells_a[ishell].norms,
                                               shells_b[jshell].norms));
            }
    }
}