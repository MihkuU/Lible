#include <lible/ints/boys_function.hpp>

namespace LI = lible::ints;

using std::vector;

vector<double> LI::BoysF::preEvaluate() const
{
    vector<double> fnx_grid(n_intervals * (max_n + 1), 0);

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

    return fnx_grid;
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

LI::BoysGrid::BoysGrid()
{
}

LI::BoysGrid::BoysGrid(const int max_n) : max_n(max_n + 6)
{
    n_intervals = large_x / interval_size + 1;
    fnx_grid = preEvaluate();
}

double LI::BoysGrid::getLargeX() const
{
    return large_x;
}

double LI::BoysGrid::getIntervalSize() const
{
    return interval_size;
}

int LI::BoysGrid::getMaxN() const
{
    return max_n;
}

const vector<double>& LI::BoysGrid::getFnxGrid() const
{
    return fnx_grid;
}

vector<double> LI::BoysGrid::preEvaluate() const
{
    vector<double> fnx_grid(n_intervals * (max_n + 1), 0);

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
        {
            double val = (2.0 * x * fnx_grid[ival * (max_n + 1) + (n + 1)] + exp_x) / (2 * n + 1);
            fnx_grid[ival * (max_n + 1) + n] = val;
        }

        x += interval_size;
    }

    return fnx_grid;
}

vector<double> LI::calcBoysF(const int n, const double x, const BoysGrid &boys_grid)
{
    double large_x = boys_grid.getLargeX();
    double interval_size = boys_grid.getIntervalSize();
    int max_n = boys_grid.getMaxN();
    const vector<double> &fnx_grid = boys_grid.getFnxGrid();

    vector<double> fnx(n + 1, 0);
    if (x == 0)
    {
        fnx[0] = 1;
        for (int k = 1; k <= n; k++)
            fnx[k] = 1.0 / (2 * k + 1);
    }
    else if (x > large_x)
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

    return fnx;
}