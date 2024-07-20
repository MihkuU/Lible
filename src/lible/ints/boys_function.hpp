#pragma once

#include <armadillo>

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
                preEvaluate(fnx_grid);
            }

            /**
             *
             */
            void calcFnx(const int max_n, const double x, std::vector<double> &fnx) const;

        private:
            int max_n;

            double boys_f_threshold = 1e-16;
            double interval_size = 0.01;
            double large_x = 30;

            arma::dmat fnx_grid;

            void preEvaluate(arma::dmat &fnx_grid);
        };
    }
}