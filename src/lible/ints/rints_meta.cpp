#include <lible/ints/rints_meta.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/utils.hpp>

#include <array>
#include <tuple>

namespace LI = lible::ints;

using std::array, std::tuple, std::vector;

// ///////////////////////////////////////////////////////////////////////////////////////////////////

namespace lible::ints
{
    // constexpr int indexRR(int n, int i, int j, int k)
    // {
    //     int ijk = i + j + k;
    //     if (ijk == 0)
    //         return n;
    //     else
    //     {
    //         int offset_m = numHermitesSum(n - 1) - n;  // n zeros.
    //         int offset_ijk = numHermites(ijk - 1) - 1; // one zero.
    //         int idx_cart = indexCart(i, j, k);

    //         return offset_m + offset_ijk + idx_cart;
    //     }
    // }

    /** */
    // struct RRItem
    // {
    //     int idx0;
    //     int idx1;
    //     int axis;
    //     int fac;
    // };

    // template <int size>
    // using rr_table_t = std::array<RRItem, size>;

    // template <int l>
    // constexpr rr_table_t<LI::numHermites(l) + l - 1> returnRRTable()
    // {
    //     rr_table_t<numHermites(l) - 1> rr_table;

    //     // int start = l + 1;
    //     int idx = 0;
    //     for (int n = 1; n <= l; n++)
    //     {
    //         for (int m = 1; m <= n; m++)
    //             for (int i = m; i >= 0; i--)
    //                 for (int j = m - i; j >= 0; j--)
    //                 {
    //                     int k = m - i - j;

    //                     RRItem rr_item;
    //                     if (i > 0)
    //                     {
    //                         rr_item.axis = 0;
    //                         rr_item.fac = i - 1;

    //                         if (rr_item.fac == 0)
    //                         {
    //                         }
    //                         else
    //                         {
    //                         }
    //                     }
    //                     else
    //                     {
    //                         if (j > 0)
    //                         {
    //                             rr_item.axis = 1;
    //                             rr_item.fac = j - 1;
    //                             if (rr_item.fac == 0)
    //                             {
    //                             }
    //                             else
    //                             {
    //                             }
    //                         }

    //                         else if (k > 0)
    //                         {
    //                             rr_item.axis = 2;
    //                             rr_item.fac = k - 1;
    //                             if (rr_item.fac == 0)
    //                             {
    //                             }
    //                             else
    //                             {
    //                             }
    //                         }
    //                     }

    //                     rr_table[idx] == rr_item;
    //                 }
    //     }

    //     return rr_table;
    // }



    template <int la, int lb, int dim>
    consteval array<tuple<int, int, int, int>, dim> generateRIntsRecursionTable();
}

// Instantiations
namespace lible::ints
{
    // template void calcRInts<0, 0>(const double, const double, double *, double *);
    // template void calcRInts<1, 0>(const double, const double, double *, double *);
    // template void calcRInts<1, 1>(const double, const double, double *, double *);
    // template void calcRInts<2, 0>(const double, const double, double *, double *);
    // template void calcRInts<2, 1>(const double, const double, double *, double *);
    // template void calcRInts<2, 2>(const double, const double, double *, double *);
    // template void calcRInts<3, 0>(const double, const double, double *, double *);
    // template void calcRInts<3, 1>(const double, const double, double *, double *);
    // template void calcRInts<3, 2>(const double, const double, double *, double *);
    // template void calcRInts<3, 3>(const double, const double, double *, double *);
    // template void calcRInts<4, 0>(const double, const double, double *, double *);
    // template void calcRInts<4, 1>(const double, const double, double *, double *);
    // template void calcRInts<4, 2>(const double, const double, double *, double *);
    // template void calcRInts<4, 3>(const double, const double, double *, double *);
    // template void calcRInts<4, 4>(const double, const double, double *, double *);
    // template void calcRInts<5, 0>(const double, const double, double *, double *);
    // template void calcRInts<5, 1>(const double, const double, double *, double *);
    // template void calcRInts<5, 2>(const double, const double, double *, double *);
    // template void calcRInts<5, 3>(const double, const double, double *, double *);
    // template void calcRInts<5, 4>(const double, const double, double *, double *);
    // template void calcRInts<5, 5>(const double, const double, double *, double *);
    // template void calcRInts<6, 0>(const double, const double, double *, double *);
    // template void calcRInts<6, 1>(const double, const double, double *, double *);
    // template void calcRInts<6, 2>(const double, const double, double *, double *);
    // template void calcRInts<6, 3>(const double, const double, double *, double *);
    // template void calcRInts<6, 4>(const double, const double, double *, double *);
    // template void calcRInts<6, 5>(const double, const double, double *, double *);
    // template void calcRInts<6, 6>(const double, const double, double *, double *);
    // template void calcRInts<7, 0>(const double, const double, double *, double *);
    // template void calcRInts<7, 1>(const double, const double, double *, double *);
    // template void calcRInts<7, 2>(const double, const double, double *, double *);
    // template void calcRInts<7, 3>(const double, const double, double *, double *);
    // template void calcRInts<7, 4>(const double, const double, double *, double *);
    // template void calcRInts<7, 5>(const double, const double, double *, double *);
    // template void calcRInts<7, 6>(const double, const double, double *, double *);
    // template void calcRInts<7, 7>(const double, const double, double *, double *);
    // template void calcRInts<8, 0>(const double, const double, double *, double *);
    // template void calcRInts<8, 1>(const double, const double, double *, double *);
    // template void calcRInts<8, 2>(const double, const double, double *, double *);
    // template void calcRInts<8, 3>(const double, const double, double *, double *);
    // template void calcRInts<8, 4>(const double, const double, double *, double *);
    // template void calcRInts<8, 5>(const double, const double, double *, double *);
    // template void calcRInts<8, 6>(const double, const double, double *, double *);
    // template void calcRInts<8, 7>(const double, const double, double *, double *);
    // template void calcRInts<8, 8>(const double, const double, double *, double *);
    // template void calcRInts<9, 0>(const double, const double, double *, double *);
    // template void calcRInts<9, 1>(const double, const double, double *, double *);
    // template void calcRInts<9, 2>(const double, const double, double *, double *);
    // template void calcRInts<9, 3>(const double, const double, double *, double *);
    // template void calcRInts<9, 4>(const double, const double, double *, double *);
    // template void calcRInts<9, 5>(const double, const double, double *, double *);
    // template void calcRInts<9, 6>(const double, const double, double *, double *);
    // template void calcRInts<9, 7>(const double, const double, double *, double *);
    // template void calcRInts<9, 8>(const double, const double, double *, double *);
    // template void calcRInts<9, 9>(const double, const double, double *, double *);
    // template void calcRInts<10, 0>(const double, const double, double *, double *);
    // template void calcRInts<10, 1>(const double, const double, double *, double *);
    // template void calcRInts<10, 2>(const double, const double, double *, double *);
    // template void calcRInts<10, 3>(const double, const double, double *, double *);
    // template void calcRInts<10, 4>(const double, const double, double *, double *);
    // template void calcRInts<10, 5>(const double, const double, double *, double *);
    // template void calcRInts<10, 6>(const double, const double, double *, double *);
    // template void calcRInts<10, 7>(const double, const double, double *, double *);
    // template void calcRInts<10, 8>(const double, const double, double *, double *);
    // template void calcRInts<10, 9>(const double, const double, double *, double *);
    // template void calcRInts<10, 10>(const double, const double, double *, double *);
    // template void calcRInts<11, 0>(const double, const double, double *, double *);
    // template void calcRInts<11, 1>(const double, const double, double *, double *);
    // template void calcRInts<11, 2>(const double, const double, double *, double *);
    // template void calcRInts<11, 3>(const double, const double, double *, double *);
    // template void calcRInts<11, 4>(const double, const double, double *, double *);
    // template void calcRInts<11, 5>(const double, const double, double *, double *);
    // template void calcRInts<11, 6>(const double, const double, double *, double *);
    // template void calcRInts<11, 7>(const double, const double, double *, double *);
    // template void calcRInts<11, 8>(const double, const double, double *, double *);
    // template void calcRInts<11, 9>(const double, const double, double *, double *);
    // template void calcRInts<11, 10>(const double, const double, double *, double *);
    // template void calcRInts<11, 11>(const double, const double, double *, double *);
    // template void calcRInts<12, 0>(const double, const double, double *, double *);
    // template void calcRInts<12, 1>(const double, const double, double *, double *);
    // template void calcRInts<12, 2>(const double, const double, double *, double *);
    // template void calcRInts<12, 3>(const double, const double, double *, double *);
    // template void calcRInts<12, 4>(const double, const double, double *, double *);
    // template void calcRInts<12, 5>(const double, const double, double *, double *);
    // template void calcRInts<12, 6>(const double, const double, double *, double *);
    // template void calcRInts<12, 7>(const double, const double, double *, double *);
    // template void calcRInts<12, 8>(const double, const double, double *, double *);
    // template void calcRInts<12, 9>(const double, const double, double *, double *);
    // template void calcRInts<12, 10>(const double, const double, double *, double *);
    // template void calcRInts<12, 11>(const double, const double, double *, double *);
    // template void calcRInts<12, 12>(const double, const double, double *, double *);
}