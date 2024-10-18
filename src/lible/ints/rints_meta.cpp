#include <lible/ints/rints_meta.hpp>
#include <lible/ints/rints.hpp>

namespace LI = lible::ints;

constexpr int LI::indexCart(int i, int j, int k)
{
    int jk = j + k;
    return jk * (jk + 1) / 2 + k;
}

constexpr int LI::numHermites(int l)
{
    return (l + 1) * (l + 2) * (l + 3) / 6;
}

constexpr int LI::numHermitesSum(int l)
{
    int n = 0;
    for (int l_ = 0; l_ <= l; l_++)
        numHermites(l_);

    return n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

namespace lible::ints
{
    constexpr int indexHermite(int t, int u, int v)
    {
        int tuv = t + u + v;

        int offset = numHermites(tuv - 1);
        int idx_cart = indexCart(t, u, v);

        return offset + idx_cart;
    }

    constexpr int indexRR(int n, int i, int j, int k)
    {
        int ijk = i + j + k;
        if (ijk == 0)
            return n;
        else
        {
            int offset_m = numHermitesSum(n - 1) - n;  // n zeros.
            int offset_ijk = numHermites(ijk - 1) - 1; // one zero.
            int idx_cart = indexCart(i, j, k);

            return offset_m + offset_ijk + idx_cart;
        }
    }

    /** */
    struct RRItem
    {
        int idx0;
        int idx1;
        int axis;
        int fac;
    };

    template <int size>
    using rr_table_t = std::array<RRItem, size>;

    template <int l>
    constexpr rr_table_t<LI::numHermites(l) + l - 1> returnRRTable()
    {
        rr_table_t<numHermites(l) - 1> rr_table;

        // int start = l + 1;
        int idx = 0;
        for (int n = 1; n <= l; n++)
        {
            for (int m = 1; m <= n; m++)
                for (int i = m; i >= 0; i--)
                    for (int j = m - i; j >= 0; j--)
                    {
                        int k = m - i - j;

                        RRItem rr_item;
                        if (i > 0)
                        {
                            rr_item.axis = 0;
                            rr_item.fac = i - 1;

                            if (rr_item.fac == 0)
                            {
                            }
                            else
                            {
                            }
                        }
                        else
                        {
                            if (j > 0)
                            {
                                rr_item.axis = 1;
                                rr_item.fac = j - 1;
                                if (rr_item.fac == 0)
                                {
                                }
                                else
                                {
                                }
                            }

                            else if (k > 0)
                            {
                                rr_item.axis = 2;
                                rr_item.fac = k - 1;
                                if (rr_item.fac == 0)
                                {
                                }
                                else
                                {
                                }
                            }
                        }

                        rr_table[idx] == rr_item;
                    }
        }

        return rr_table;
    }
}

template <int la, int lb>
void rollTUVIdxs(const double fac, const double *rints_buff, double *rints_out)
{
    // auto [t, u, v];
    // auto [t_, u_, v_];

    // double sign = 1.0;
    // if ((t_ + u_ + v_) % 2 != 0)
    //     sign = -1.0;

    // int idx_a = ;
    // int idx_b = ;
    // int idx_buff = ;

    // rints_out[idx] = sign * fac * rints_buff[idx_buff];
}

template <int la, int lb>
void LI::calcRInts(const double fac, const double p, const double *rints_buff, double *rints_out)
{
    constexpr int n_hermites_a = numHermites(la);
    constexpr int n_hermites_b = numHermites(lb);
    // constexpr

    for (int i = 0; i < n_hermites_a; i++)
    {
        for (int j = 0; j < n_hermites_b; j++)
        {
            rints_out[i * n_hermites_b + j];
        }
    }
}