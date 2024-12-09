#include <lible/ints/utils.hpp>

namespace LI = lible::ints;

using std::array, std::map, std::pair, std::vector;

double LI::calcPurePrimitiveNorm(const int l, const double exp)
{
    double norm = std::sqrt(std::pow(2 * exp / M_PI, 1.5) * std::pow(4 * exp, l) /
                                    doubleFactorial(2 * l - 1));

    return norm;
}

double LI::doubleFactorial(const int n)
{    
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

vector<array<int, 3>> LI::cartExps(const int l)
{
    size_t dim_cart = numCartesians(l);

    vector<array<int, 3>> cartesian_exps(dim_cart);
    for (int x = l, pos = 0; x >= 0; x--)
        for (int y = l - x; y >= 0; y--)
        {
            int z = l - x - y;
            cartesian_exps[pos] = {x, y, z};
            pos++;
        }

    return cartesian_exps;
}

vector<pair<int, int>> LI::returnLPairs(const int l_max)
{
    int n_pairs = (l_max + 1) * (l_max + 2) / 2;
    vector<pair<int, int>> l_pairs(n_pairs);

    for (int la = 0, idx = 0; la <= l_max; la++)
        for (int lb = 0; lb <= la; lb++, idx++)
            l_pairs[idx] = std::make_pair(la, lb);

    return l_pairs;
}

map<int, vector<pair<int, int>>> LI::returnLPairsMap(const int l_max)
{
    int l_max_pair = 2 * l_max;

    map<int, vector<pair<int, int>>> l_pairs_map;
    for (int lab = 0; lab <= l_max_pair; lab++)
    {
        vector<pair<int, int>> lpairs;
        for (int la = 0; la <= l_max; la++)
            for (int lb = 0; lb <= la; lb++)
                if (la + lb == lab)
                    lpairs.push_back(std::make_pair(la, lb));

        l_pairs_map[lab] = lpairs;
    }

    return l_pairs_map;
}

lible::vec3i LI::returnHermiteGaussianPositions(const int l)
{
    vec3i tuv_poss(l + 1, l + 1, l + 1, -1);
    for (int n = 0, tuv = 0; n <= l; n++)
        for (int t = n; t >= 0; t--)
            for (int u = n - t; u >= 0; u--, tuv++)
            {
                int v = n - t - u;
                tuv_poss(t, u, v) = tuv;
            }

    return tuv_poss;
}

vector<array<int, 3>> LI::returnHermiteGaussianIdxs(const int l)
{
    vector<array<int, 3>> idxs_tuv((l + 1) * (l + 2) * (l + 3) / 6);
    for (int n = 0, tuv = 0; n <= l; n++)
        for (int t = n; t >= 0; t--)
            for (int u = n - t; u >= 0; u--, tuv++)
            {
                int v = n - t - u;
                idxs_tuv[tuv] = {t, u, v};                
            }

    return idxs_tuv;
}