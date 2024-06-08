#include "geomopt_impl.h"

namespace LG = lible::geomopt;
using namespace LG;

using std::vector;

vector<double> LG::Optimizer<GDESCENT>::update(const vector<double> &coords,
                                               const vector<double> &grad,
                                               const size_t &iter)
{
    size_t dim = coords.size();
    vector<double> step(dim, 0);

    if (iter == 0)
    {
        for (size_t i = 0; i < dim; i++)
            step[i] = -grad[i];
    }
    else
    {
        double num = 0, den = 0;
        for (size_t i = 0; i < dim; i++)
        {        
            num += pow(grad[i], 2);
            den += pow(grad_prev[i], 2);
        }

        double beta = num / den;
        for (size_t i = 0; i < dim; i++)        
            step[i] = -grad[i] + beta * step_prev[i];
    }

    grad_prev = grad;
    step_prev = step;

    // printf("\ncoords:");
    // for (auto &coord : coords)
    //     printf("coord = %16.12lf\n", coord);

    // for (auto &s : step)
    //     printf("s = %16.12lf\n", s);

    return step;
}
