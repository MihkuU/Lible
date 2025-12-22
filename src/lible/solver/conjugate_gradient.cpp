#include <lible/solver/solver.hpp>

#include <armadillo>

namespace lsolver = lible::solver;

namespace lible::solver
{
    arma::dvec conv2dvec(const std::vector<double> &vec);

    std::vector<double> conv2vec(const arma::dvec &vec);

    /// Wrapper function for calculating the sigma vector with arma types.
    arma::dvec calcSigmaWrap(const cg_sigma_t &calc_sigma, const arma::dvec &trial);

    /// Wrapper function for calculating the preconditioned residual with arma types.
    arma::dvec preconditionerWrap(const pcg_preconditioner_t &preconditioner,
                                  const arma::dvec &residual);
}

arma::dvec lsolver::conv2dvec(const std::vector<double> &vec)
{
    return arma::conv_to<arma::dvec>::from(vec);
}

std::vector<double> lsolver::conv2vec(const arma::dvec &vec)
{
    return arma::conv_to<std::vector<double> >::from(vec);
}

arma::dvec lsolver::calcSigmaWrap(const cg_sigma_t &calc_sigma, const arma::dvec &trial)
{
    arma::dvec sigma = conv2dvec(calc_sigma(conv2vec(trial)));
    if (sigma.size() != trial.size())
        throw std::runtime_error("calcSigmaWrap(): sigma must have the same size as trial");

    return sigma;
}

arma::dvec lible::solver::preconditionerWrap(const pcg_preconditioner_t &preconditioner,
                                             const arma::dvec &residual)
{
    arma::dvec z = conv2dvec(preconditioner(conv2vec(residual)));
    if (z.size() != residual.size())
        throw std::runtime_error("preconditionerWrap(): preconditioned residual must have the same "
            "size as residual");

    return z;
}

lsolver::CGResults lsolver::conjugateGradient(const std::vector<double> &guess,
                                              const std::vector<double> &rhs,
                                              const cg_sigma_t &calc_sigma,
                                              const CGSettings &settings)
{
    if (guess.size() != rhs.size())
        throw std::runtime_error("conjugateGradient(): guess and rhs must have the same size");

    const auto &[max_iter, conv_tol] = settings;

    arma::dvec x = conv2dvec(guess);
    arma::dvec sigma = calcSigmaWrap(calc_sigma, x);
    arma::dvec r = conv2dvec(rhs) - sigma;
    arma::dvec p = r;

    std::vector<double> max_abs_residuals;
    bool converged = false;
    size_t iter = 0;
    for (; iter < max_iter; iter++)
    {
        double max_abs_res = arma::max(arma::abs(r));
        max_abs_residuals.push_back(max_abs_res);
        if (max_abs_res < conv_tol)
        {
            converged = true;
            break;
        }

        // Solution update
        sigma = calcSigmaWrap(calc_sigma, p);
        double alpha = arma::dot(r, r) / arma::dot(p, sigma);
        x += alpha * p;

        // Trial/residual vector update
        arma::dvec r_prev = r;
        r = r - alpha * sigma;
        double beta = arma::dot(r, r) / arma::dot(r_prev, r_prev);
        p = r + beta * p;
    }

    return {converged, iter, max_abs_residuals, conv2vec(x)};
}

lsolver::CGResults lsolver::preconditionedCG(const std::vector<double> &guess,
                                             const std::vector<double> &rhs,
                                             const cg_sigma_t &calc_sigma,
                                             const pcg_preconditioner_t &preconditioner,
                                             const CGSettings &settings)
{
    if (guess.size() != rhs.size())
        throw std::runtime_error("conjugateGradient(): guess and rhs must have the same size");

    const auto &[max_iter, conv_tol] = settings;

    arma::dvec x = conv2dvec(guess);
    arma::dvec sigma = calcSigmaWrap(calc_sigma, x);
    arma::dvec r = conv2dvec(rhs) - sigma;
    arma::dvec z = preconditionerWrap(preconditioner, r);
    arma::dvec p = z;

    std::vector<double> max_abs_residuals;
    bool converged = false;
    size_t iter = 0;
    for (; iter < max_iter; iter++)
    {
        double max_abs_res = arma::max(arma::abs(r));
        max_abs_residuals.push_back(max_abs_res);
        if (max_abs_res < conv_tol)
        {
            converged = true;
            break;
        }

        // Solution update
        sigma = calcSigmaWrap(calc_sigma, p);
        double alpha = arma::dot(r, z) / arma::dot(p, sigma);
        x += alpha * p;

        // Trial/residual vector update
        arma::dvec r_prev = r;
        r = r - alpha * sigma;
        arma::dvec z_prev = z;
        z = preconditionerWrap(preconditioner, r);
        double beta = arma::dot(r, z) / arma::dot(r_prev, z_prev);
        p = r + beta * p;
    }

    return {converged, iter, max_abs_residuals, conv2vec(x)};
}
