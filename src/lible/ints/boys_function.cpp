#include <lible/ints/boys_function.hpp>

#include <cmath>

namespace LI = lible::ints;

using std::vector;

void LI::BoysF::preEvaluate(vector<double> &fnx_grid)
{
    fnx_grid.resize(n_intervals * (max_n + 1), 0);

    double x = 0;
    for (int ival = 0; ival < n_intervals; ival++)
    {
        int k = 1;
        double exp_x = std::exp(-x);
        double term = 1.0 / (2.0 * max_n + 1.0);
        double sum = term;
        while (true)
        {
            term *= x / (max_n + 0.5 + k);
            sum += term;
            if (std::fabs(term) < boys_f_threshold)
                break;
            k++;
        }

        fnx_grid[ival * (max_n + 1) + max_n] = exp_x * sum;
        for (int n = max_n - 1; n >= 0; n--)
            fnx_grid[ival * (max_n + 1) + n] = (2.0 * x * fnx_grid[ival * (max_n + 1) + (n + 1)] + exp_x) / (2 * n + 1);
        
        x += interval_size;
    }
}

void LI::BoysF::calcFnx(const int n, const double x, vector<double> &fnx) const
{    
    if (x == 0)
    {
        fnx[0] = 1;
        for (int k = 1; k <= n; k++)
            fnx[k] = 1.0 / (2 * k + 1);
    }
    else if (x > 30.0)
    {
        // Adapted from HUMMR
        fnx[0] = 0.5 * std::sqrt(M_PI / x);
        for (int k = 1; k <= n; k++)
            fnx[k] = fnx[k - 1] * (k - 0.5) / x;
    }
    else
    {
        int ival = x / interval_size;

        double origin_x = ival * interval_size;
        double delta_x = x - origin_x;

        double k_factorial = 1.0;
        double deltax_k = 1.0;

        double sum = fnx_grid[ival * (max_n + 1) + n];
        for (int k = 1; k < 7; k++)
        {
            deltax_k *= -delta_x;
            k_factorial *= k;
            sum += fnx_grid[ival * (max_n + 1) + (n + k)] * deltax_k / k_factorial;
        }

        fnx[n] = sum;

        double exp_x = std::exp(-x);
        for (int k = n - 1; k >= 0; k--)
            fnx[k] = (2.0 * x * fnx[k + 1] + exp_x) / (2 * k + 1);
    }
}