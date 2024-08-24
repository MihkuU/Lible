#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

#include <armadillo>

namespace lible
{
    namespace ints
    {
        /**
         *
         */
        void calcRInts(const int la, const int lb, const double p,
                       const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                       vec4d &rints_tmp, vec3d &rints_out);

        /**
         *
         */
        void calcRInts(const int la, const int lb, const double p,
                       const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                       const std::vector<IdxsTUV> &tuv_idxs_a,
                       const std::vector<IdxsTUV> &tuv_idxs_b,
                       vec4d &rints_tmp, arma::dmat &rints_out);

        /**
         *
         */
        void calcRInts(const int la, const int lb, const double p,
                       const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                       const std::vector<IdxsTUV> &tuv_idxs_a,
                       const std::vector<IdxsTUV> &tuv_idxs_b,
                       vec4d &rints_tmp, std::vector<double> &rints_out);

        /**
         *
         */
        void calcRInts(const int la, const int lb, const double p, const double fac,
                       const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                       const std::vector<IdxsTUV> &tuv_idxs_a,
                       const std::vector<IdxsTUV> &tuv_idxs_b,
                       vec4d &rints_tmp, std::vector<double> &rints_out);

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
