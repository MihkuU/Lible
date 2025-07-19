#include <lible/ints/rints_meta.hpp>

namespace lible::ints
{
template<>
void calcRInts_ERI<0, 2>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 2;
    constexpr int buff_size = 12;
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
    rints_buff[3] = xyz_ab[0] * rints_buff[2];
    rints_buff[4] = xyz_ab[1] * rints_buff[2];
    rints_buff[5] = xyz_ab[2] * rints_buff[2];
    rints_buff[6] = xyz_ab[0] * rints_buff[3] + 1 * rints_buff[1];
    rints_buff[7] = xyz_ab[0] * rints_buff[4];
    rints_buff[8] = xyz_ab[0] * rints_buff[5];
    rints_buff[9] = xyz_ab[1] * rints_buff[4] + 1 * rints_buff[1];
    rints_buff[10] = xyz_ab[1] * rints_buff[5];
    rints_buff[11] = xyz_ab[2] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[3] = xyz_ab[0] * rints_buff[1];
    rints_buff[4] = xyz_ab[1] * rints_buff[1];
    rints_buff[5] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[1] = -1.0 * fac * rints_buff[3];
    rints_out[2] = -1.0 * fac * rints_buff[4];
    rints_out[3] = -1.0 * fac * rints_buff[5];
    rints_out[4] = 1.0 * fac * rints_buff[6];
    rints_out[5] = 1.0 * fac * rints_buff[7];
    rints_out[6] = 1.0 * fac * rints_buff[8];
    rints_out[7] = 1.0 * fac * rints_buff[9];
    rints_out[8] = 1.0 * fac * rints_buff[10];
    rints_out[9] = 1.0 * fac * rints_buff[11];
}
}

template void lible::ints::calcRInts_ERI2D1<0, 2>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI2D2<0, 2>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_ab, double *rints);

template void lible::ints::calcRInts_ERI3D1<0, 2>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

template void lible::ints::calcRInts_ERISOC<0, 2>(const double alpha, const double fac, const double *fnx,
                                                  const double *xyz_pc, double *rints);

