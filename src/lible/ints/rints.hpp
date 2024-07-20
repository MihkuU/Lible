#pragma once

#include <array>
#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         *
         */
        struct RData
        {
            double coeff;
            int idx0;
            int idx1;
            int axis;
        };

        // /**
        //  *
        //  */
        // template <int la, int lb>
        // void calcRInts(const std::array<RData> &rdata_table);

        // /**
        //  *
        //  */
        // template <int la, int lb>
        // constexpr std::array<RData> generateRRecurrenceTable();

        // template <int la, int lb, int dim>
        // std::array<RData, dim> generateRRecurrenceTable();

        /** */
        constexpr int calcCartDimSum(const int l);

        /** */
        constexpr int calcCartDimSum(const int la, const int lb);

        /** */
        constexpr int calcIdx(const int i, const int j, const int k);

        /** */
        template <int la, int lb>
        constexpr std::array<RData, calcCartDimSum(la, lb) + 1> generateRRecurrenceTable();

        // /**
        //  *
        //  */
        //  std::array<RData> returnRRecurrenceTable(const int la, const int lb);

    }
}
