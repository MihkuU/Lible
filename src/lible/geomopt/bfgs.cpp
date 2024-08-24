#include "geomopt_impl.h"

using namespace lible; 
using geomopt::BFGS;

using std::vector;

vector<double> geomopt::Optimizer<BFGS>::update(const vector<double> &coords_redint_in, const vector<double> &grad_redint)
{
    vector<double> coords_redint_out;

    return coords_redint_out;
}

