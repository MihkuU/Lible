#pragma once

#include <array>
#include <cmath>
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
                fnx_grid = preEvaluate();
            }

            /**
             *
             */
            void calcFnx(const int max_n, const double x, std::vector<double> &fnx) const; // TODO: make return

        private:
            int max_n;
            int n_intervals;

            double boys_f_threshold = 1e-16;
            double interval_size = 0.01;
            double large_x = 30;

            std::vector<double> fnx_grid;

            std::vector<double> preEvaluate() const;
        };

        /** */
        class BoysGrid
        {
        public:
            BoysGrid();
            
            BoysGrid(const int max_n);

            /** */
            double getLargeX() const;

            /** */
            double getIntervalSize() const;

            /** */
            int getMaxN() const;

            /** */
            // std::vector<double> getFnxGrid() const;

            const std::vector<double>& getFnxGrid() const;

        private:
            double boys_f_threshold = 1e-16; /** */
            double interval_size = 0.01;     /** */
            double large_x = 30;             /** */

            int max_n;       /** */
            int n_intervals; /** */

            std::vector<double> fnx_grid; /** */

            /** */
            std::vector<double> preEvaluate() const;
        };

        std::vector<double> calcBoysF(const int max_n, const double x, const BoysGrid &boys_grid);

        /**                  
         * The Boys function is represented as a 2D grid (matrix) that is rolled out as a vector. 
         * The number of rows corresponds to each value of x from 0 to 30 with a distance of 0.01. 
         * The rows correspond to different angular momentum plus the number of terms in the 
         * Taylor series (eq. (9.8.12) from HJO).        
         *          
         * TODO: Figure out a better name. 
         * TODO: Move this version of the function to a separate file.
         */
        template <const int l>
        class BoysF2
        {                        
        public:
            /** Calculates the Boys function at value x. */
            void calcFnx(const double x, double *fnx) const
            {
                if (x == 0)
                {
                    // (9.8.6) from the bible.
                    fnx[0] = 1;
                    for (int k = 1; k <= l; k++)
                        fnx[k] = 1.0 / (2 * k + 1);
                }
                else if (x > 30.0)
                {
                    // Adapted from HUMMR, should be (9.8.9) in the book.
                    fnx[0] = 0.5 * std::sqrt(M_PI / x);
                    for (int k = 1; k <= l; k++)
                        fnx[k] = fnx[k - 1] * (k - 0.5) / x;
                }
                else
                {
                    // (9.8.12) from HJO.
                    int ival = x / interval_size;

                    double origin_x = ival * interval_size;
                    double delta_x = x - origin_x;

                    double k_factorial = 1.0;
                    double deltax_k = 1.0;

                    double sum = fnx_grid[ival * n_cols + l];
                    for (int k = 1; k < 7; k++)
                    {
                        deltax_k *= -delta_x;
                        k_factorial *= k;
                        sum += fnx_grid[ival * n_cols + (l + k)] * deltax_k / k_factorial;
                    }

                    fnx[l] = sum;

                    double exp_x = std::exp(-x);
                    for (int k = l - 1; k >= 0; k--)
                        fnx[k] = (2.0 * x * fnx[k + 1] + exp_x) / (2 * k + 1);
                }
            }

        private:
            static constexpr double boys_f_threshold = 1e-16;

            static constexpr double interval_size = 0.01;
            
            static constexpr double large_x = 30;            

            static constexpr int n_terms = 7;

            static constexpr int n_rows = large_x / interval_size + 1;

            static constexpr int n_cols = l + n_terms;            

            static consteval int calcFnxGridSize()
            {
                return n_rows * n_cols;
            }

            /** Evaluates the Boys function at the grid points based on (9.8.11) from HJO. */
            static consteval std::array<double, calcFnxGridSize()> preEvaluateGrid()
            {
                std::array<double, calcFnxGridSize()> fnx_grid;

                double x = 0;
                for (int ival = 0; ival < n_rows; ival++)
                {
                    int k = 1;
                    double exp_x = std::exp(-x);
                    double term = 1.0 / (2.0 * (n_cols - 1) + 1.0);
                    double sum = term;
                    while (true)
                    {
                        term *= x / ((n_cols - 1) + 0.5 + k);
                        sum += term;
                        if (std::fabs(term) < boys_f_threshold)
                            break;
                        k++;
                    }

                    fnx_grid[ival * n_cols + (n_cols - 1)] = exp_x * sum;
                    for (int k = n_cols - 2; k >= 0; k--)
                        fnx_grid[ival * n_cols + k] = (2.0 * x * fnx_grid[ival * n_cols + (k + 1)] + exp_x) / (2 * k + 1);

                    x += interval_size;
                }

                return fnx_grid;
            }

            static constexpr std::array<double, calcFnxGridSize()> fnx_grid{preEvaluateGrid()};
        };
    }
}