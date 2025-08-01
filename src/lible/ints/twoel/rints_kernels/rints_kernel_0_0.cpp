#include <lible/ints/rints_meta.hpp>

namespace lible::ints
{
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

template void lible::ints::calcRInts_ERI2D1<0, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI2D2<0, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI3D1<0, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

template void lible::ints::calcRInts_ERISOC<0, 0>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

