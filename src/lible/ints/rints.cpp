#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

using std::array;

// template <int la, int lb>
// constexpr array<LI::RData> LI::generateRRecurrenceTable()
// {
//     array<RData> rrecurrence_table;

//     return rrecurrence_table;
// }

// template <int la, int lb, int dim>
// std::array<LI::RData, dim> LI::generateRRecurrenceTable()
// {
// }

constexpr int LI::calcCartDimSum(const int l)
{
    return (l + 1) * (l + 2) * (l + 3) / 6;
}

constexpr int LI::calcCartDimSum(const int la, const int lb)
{
    int lab = la + lb;
    return (lab + 1) * (lab + 2) * (lab + 3) / 6;
}

constexpr int LI::calcIdx(const int i, const int j, const int k)
{
    int l = i + j + k;
    // int offset = calcRDim()
    int jk = j + k;
    return jk * (jk + 1) / 2 + k;
}

template <int la, int lb>
constexpr array<LI::RData, LI::calcCartDimSum(la, lb) + 1> LI::generateRRecurrenceTable()
{
    array<RData, calcCartDimSum(la, lb) + 1> r_recurrenc_table;

    return r_recurrenc_table;
}