#pragma once

#include <array>
#include <cmath>
#include <vector>

namespace lible::ints
{
    /// Class for precalculating and using the Boys function grid.
    class BoysGrid
    {
    public:
        /// Default ctor.
        BoysGrid() = default;

        /// Initializes the Boys function grid for the maximal value of `n`.
        explicit BoysGrid(int max_n);

        /// Returns the value of x after which the Boys function is calculated approximately.
        double getLargeX() const;

        /// Returns the interval size.
        double getIntervalSize() const;

        /// Returns the maximum n-value for which the Boys function grid is constructed.
        int getMaxN() const;

        /// Returns a constant reference to the calculated Boys function grid values.
        const std::vector<double> &getFnxGrid() const;

    private:
        /// Convergence threshold for eq. (9.8.12) from https://doi.org/10.1007/s10008-001-0256-1.
        double boys_f_threshold_ = 1e-16;

        /// Length of the Boys function grid interval.
        double interval_size_ = 0.01;

        /// Value of x after which the boys function is calculated approximately. E.g., using
        /// eq. (9.8.13) from https://doi.org/10.1007/s10008-001-0256-1.
        double large_x_ = 30;

        /// Maximum number of `n` for using the recursion to hit the target n.
        int max_n_{};
        /// Number of grid intervals.
        int n_intervals_{};

        /// The precalculated Boys function grid.
        std::vector<double> fnx_grid_;

        /// First, constructs Fn(0) for `max_n_` using eq. (9.8.6) from
        /// https://doi.org/10.1007/s10008-001-0256-1. Then, uses the downward recursion from
        /// eq. (9.8.14) to calculate the grid values for different n.
        std::vector<double> preEvaluate() const;
    };

    /// Templated class for calculating the Boys function. The Boys function grid is evaluated at
    /// compile time.
    template <const int L>
    class BoysF2
    {
    public:
        /// Calculates the Boys function at x.
        void calcFnx(const double x, double *fnx) const
        {
            if (x == 0)
            {
                // (9.8.6) from the bible.
                fnx[0] = 1;
                for (int k = 1; k <= L; k++)
                    fnx[k] = 1.0 / (2 * k + 1);
            }
            else if (x > 30.0)
            {
                // Adapted from HUMMR, should be (9.8.9) in the book.
                fnx[0] = 0.5 * std::sqrt(M_PI / x);
                for (int k = 1; k <= L; k++)
                    fnx[k] = fnx[k - 1] * (k - 0.5) / x;
            }
            else
            {
                // (9.8.12) from HJO.
                int ival = x / interval_size_;

                double origin_x = ival * interval_size_;
                double delta_x = x - origin_x;

                double k_factorial = 1.0;
                double deltax_k = 1.0;

                double sum = fnx_grid_[ival * n_cols_ + L];
                for (int k = 1; k < 7; k++)
                {
                    deltax_k *= -delta_x;
                    k_factorial *= k;
                    sum += fnx_grid_[ival * n_cols_ + (L + k)] * deltax_k / k_factorial;
                }

                fnx[L] = sum;

                double exp_x = std::exp(-x);
                for (int k = L - 1; k >= 0; k--)
                    fnx[k] = (2.0 * x * fnx[k + 1] + exp_x) / (2 * k + 1);
            }
        }

    private:
        static constexpr double boys_f_threshold_ = 1e-16;

        static constexpr double interval_size_ = 0.01;

        static constexpr double large_x_ = 30;

        static constexpr int n_terms_ = 7;

        static constexpr int n_rows_ = large_x_ / interval_size_ + 1;

        static constexpr int n_cols_ = L + n_terms_;

        static consteval int calcFnxGridSize()
        {
            return n_rows_ * n_cols_;
        }

        /// Evaluates the Boys function at the grid points based on (9.8.11) from HJO.
        static consteval std::array<double, calcFnxGridSize()> preEvaluateGrid()
        {
            std::array<double, calcFnxGridSize()> fnx_grid;

            double x = 0;
            for (int ival = 0; ival < n_rows_; ival++)
            {
                int k = 1;
                double exp_x = std::exp(-x);
                double term = 1.0 / (2.0 * (n_cols_ - 1) + 1.0);
                double sum = term;
                while (true)
                {
                    term *= x / ((n_cols_ - 1) + 0.5 + k);
                    sum += term;
                    if (std::fabs(term) < boys_f_threshold_)
                        break;
                    k++;
                }

                fnx_grid[ival * n_cols_ + (n_cols_ - 1)] = exp_x * sum;
                for (int k = n_cols_ - 2; k >= 0; k--)
                    fnx_grid[ival * n_cols_ + k] = (2.0 * x * fnx_grid[ival * n_cols_ + (k + 1)] + exp_x) / (2 * k + 1);

                x += interval_size_;
            }

            return fnx_grid;
        }

        /// Compile-time calculated Boys function grid.
        static constexpr std::array<double, calcFnxGridSize()> fnx_grid_{preEvaluateGrid()};
    };
}
