#include "geomopt_impl.h"

using namespace Lible; 
using GeomOpt::BFGS;

using std::vector;

vector<double> GeomOpt::Optimizer<BFGS>::update(const vector<double> &coords_redint_in, const vector<double> &grad_redint)
{
    vector<double> coords_redint_out;

    return coords_redint_out;
}

