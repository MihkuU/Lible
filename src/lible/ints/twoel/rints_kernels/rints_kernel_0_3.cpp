#include <lible/ints/rints_meta.hpp>

namespace lible::ints
{
template <int la, int lb>
void calcRInts_ERI(const double alpha, const double fac, const double *fnx, const double *xyz_ab,
                   double *rints_out);

template<>
void calcRInts_ERI<0, 3>(const double alpha, const double fac, const double *fnx,
                         const double *xyz_ab, double *rints_out)
{
    constexpr int lab = 3;
    constexpr int buff_size = 23;
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
    rints_buff[4] = xyz_ab[0] * rints_buff[3];
    rints_buff[5] = xyz_ab[1] * rints_buff[3];
    rints_buff[6] = xyz_ab[2] * rints_buff[3];
    rints_buff[7] = xyz_ab[0] * rints_buff[4] + 1 * rints_buff[2];
    rints_buff[8] = xyz_ab[0] * rints_buff[5];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[1] * rints_buff[5] + 1 * rints_buff[2];
    rints_buff[11] = xyz_ab[1] * rints_buff[6];
    rints_buff[12] = xyz_ab[2] * rints_buff[6] + 1 * rints_buff[2];
    rints_buff[4] = xyz_ab[0] * rints_buff[2];
    rints_buff[5] = xyz_ab[1] * rints_buff[2];
    rints_buff[6] = xyz_ab[2] * rints_buff[2];
    rints_buff[13] = xyz_ab[0] * rints_buff[7] + 2 * rints_buff[4];
    rints_buff[14] = xyz_ab[0] * rints_buff[8] + 1 * rints_buff[5];
    rints_buff[15] = xyz_ab[0] * rints_buff[9] + 1 * rints_buff[6];
    rints_buff[16] = xyz_ab[0] * rints_buff[10];
    rints_buff[17] = xyz_ab[0] * rints_buff[11];
    rints_buff[18] = xyz_ab[0] * rints_buff[12];
    rints_buff[19] = xyz_ab[1] * rints_buff[10] + 2 * rints_buff[5];
    rints_buff[20] = xyz_ab[1] * rints_buff[11] + 1 * rints_buff[6];
    rints_buff[21] = xyz_ab[1] * rints_buff[12];
    rints_buff[22] = xyz_ab[2] * rints_buff[12] + 2 * rints_buff[6];
    rints_buff[7] = xyz_ab[0] * rints_buff[4] + 1 * rints_buff[1];
    rints_buff[8] = xyz_ab[0] * rints_buff[5];
    rints_buff[9] = xyz_ab[0] * rints_buff[6];
    rints_buff[10] = xyz_ab[1] * rints_buff[5] + 1 * rints_buff[1];
    rints_buff[11] = xyz_ab[1] * rints_buff[6];
    rints_buff[12] = xyz_ab[2] * rints_buff[6] + 1 * rints_buff[1];
    rints_buff[4] = xyz_ab[0] * rints_buff[1];
    rints_buff[5] = xyz_ab[1] * rints_buff[1];
    rints_buff[6] = xyz_ab[2] * rints_buff[1];

    // R-ints rollout
    rints_out[0] = 1.0 * fac * rints_buff[0];
    rints_out[1] = -1.0 * fac * rints_buff[4];
    rints_out[2] = -1.0 * fac * rints_buff[5];
    rints_out[3] = -1.0 * fac * rints_buff[6];
    rints_out[4] = 1.0 * fac * rints_buff[7];
    rints_out[5] = 1.0 * fac * rints_buff[8];
    rints_out[6] = 1.0 * fac * rints_buff[9];
    rints_out[7] = 1.0 * fac * rints_buff[10];
    rints_out[8] = 1.0 * fac * rints_buff[11];
    rints_out[9] = 1.0 * fac * rints_buff[12];
    rints_out[10] = -1.0 * fac * rints_buff[13];
    rints_out[11] = -1.0 * fac * rints_buff[14];
    rints_out[12] = -1.0 * fac * rints_buff[15];
    rints_out[13] = -1.0 * fac * rints_buff[16];
    rints_out[14] = -1.0 * fac * rints_buff[17];
    rints_out[15] = -1.0 * fac * rints_buff[18];
    rints_out[16] = -1.0 * fac * rints_buff[19];
    rints_out[17] = -1.0 * fac * rints_buff[20];
    rints_out[18] = -1.0 * fac * rints_buff[21];
    rints_out[19] = -1.0 * fac * rints_buff[22];
}
}

template void lible::ints::calcRInts_ERI2_deriv1<0, 3>(const double, const double, const double*, const double*, double*);
