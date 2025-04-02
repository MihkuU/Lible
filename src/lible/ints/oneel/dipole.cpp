#include <lible/ints/oneel/oneel_detail.hpp>

namespace LIO = lible::ints::one;

using std::array;

template <>
void LIO::kernel<LIO::Option::dipole_moment>(const int la, const int lb,
                                             const ShellPairData &sp_data,
                                             array<lible::vec2d, 3> &ints_out)
{

}