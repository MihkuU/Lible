#include <lible/ints/ints_util.hpp>

namespace LI = lible::ints;

double LI::calcPurePrimitiveNorm(const int l, const double exp)
{
    double partial_norm = std::sqrt(std::pow(2 * exp / M_PI, 1.5) *
                                    std::pow(4 * exp, l) /
                                    doubleFactorial(2 * l - 1));

    return partial_norm;
}

double LI::doubleFactorial(const int n)
{
    /*
     * Definition from DOI:10.1002/9781119019572 in eq. (6.5.10) without the odd negative values.
     */

    double double_factorial = 1;
    if (n == 0)
        return double_factorial;
    else if (n == -1)
        return double_factorial;
    else if (n > 0 and (n % 2 == 0))
        for (int i = n; i >= 2; i -= 2)
            double_factorial *= i;
    else if (n > 0 and (n % 2 == 1))
        for (int i = n; i >= 1; i -= 2)
            double_factorial *= i;
    else
        throw std::runtime_error("False input value for the double factorial!");

    return double_factorial;
}

int LI::dimCartesians(const int l)
{
    return (l + 1) * (l + 2) / 2;
}

int LI::dimSphericals(const int l)
{
    return 2 * l + 1;
}

std::vector<std::array<int, 3>> LI::returnCartesianExps(const int l)
{
    std::size_t dim_cart = dimCartesians(l);

    std::vector<std::array<int, 3>> cartesian_exps(dim_cart);
    for (int x = l, pos = 0; x >= 0; x--)
        for (int y = l - x; y >= 0; y--)
        {
            int z = l - x - y;
            cartesian_exps[pos] = {x, y, z};
            pos++;
        }

    return cartesian_exps;
}

std::vector<std::pair<int, int>> LI::returnLPairs(const int l_max)
{
    int n_pairs = (l_max + 1) * (l_max + 2) / 2;
    std::vector<std::pair<int, int>> l_pairs(n_pairs);

    for (int la = l_max, idx = n_pairs - 1; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--, idx--)
            l_pairs[idx] = std::make_pair(la, lb);

    return l_pairs;
}