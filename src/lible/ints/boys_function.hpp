#pragma once

#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         *
         */
        class BoysF
        {
        public:
            BoysF(const int max_n) : max_n(max_n + 6)
            {
                n_intervals = large_x / interval_size + 1;
                preEvaluate(fnx_grid);
            }

            /**
             *
             */
            void calcFnx(const int max_n, const double x, std::vector<double> &fnx) const;

        private:
            int max_n;
            int n_intervals;

            double boys_f_threshold = 1e-16;
            double interval_size = 0.01;
            double large_x = 30;

            std::vector<double> fnx_grid;

            void preEvaluate(std::vector<double> &fnx_grid);
        };
    }
}