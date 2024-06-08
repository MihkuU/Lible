#include <lible/boys_function.hpp>

namespace LI = lible::ints;

namespace lible::ints
{
    static double boys_f_threshold = 1e-16;
}

void LI::BoysF::preEvaluate(arma::dmat &fnx_grid)
{
    // Adapted from HUMMR
    size_t n_intervals = large_x / interval_size + 1;
    fnx_grid.zeros(n_intervals, max_n + 1);

    double x = 0;
    for (size_t ival = 0; ival < n_intervals; ival++)
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
        fnx_grid(ival, max_n) = exp_x * sum;
        for (int n = max_n - 1; n >= 0; n--)        
            fnx_grid(ival, n) = (2.0 * x * fnx_grid(ival, n + 1) + exp_x) / (2 * n + 1);        

        x += interval_size;
    }    
}

void LI::BoysF::calcFnx(const int max_n, const double x, std::vector<double> &fnx) const
{
    // Adapted from HUMMR
    if (x > 30.0)
    {
        fnx[0] = 0.5 * std::sqrt(M_PI / x);
        for (int n = 1; n <= max_n; n++)
            fnx[n] = fnx[n - 1] * (n - 0.5) / x;
    }
    else 
    {
        size_t ival = x / interval_size;

        double origin_x = ival * interval_size;
        double delta_x = x - origin_x;

        double k_factorial = 1.0;
        double deltax_k = 1.0;

        double sum = fnx_grid(ival, max_n);
        for (int k = 1; k < 7; k++)
        {
            deltax_k *= -delta_x;
            k_factorial *= k;
            sum += fnx_grid(ival, max_n + k) * deltax_k / k_factorial;
        }

        fnx[max_n] = sum;        

        double exp_x = std::exp(-x);
        for (int n = max_n - 1; n >= 0; n--)
            fnx[n] = (2.0 * x * fnx[n + 1] + exp_x) / (2 * n + 1);
    }
}