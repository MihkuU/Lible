#include <lible/ints/rints_meta.hpp>

namespace lible::ints
{
template <int la, int lb>
void calcRInts_ERI(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
                   double *rints_out);

template<>
void calcRInts_ERI<0, 0>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 0;
    constexpr int buff_size = 1;
    std::array<double, buff_size> rints_buff{};

    rints_buff[0] = fnx[0];
    double x = -2 * alpha;
    double y = x;
    for (int n = 1; n <= lab; n++)
    {
        rints_buff[n] = fnx[n] * y;
        y *= x;
    }

    // R-ints recursion

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
}
}

template void lible::ints::calcRInts_ERI2_deriv1<0, 0>(const double, const double, const double*, const double*, double*);
