#include <tests.hpp>

#include <lible/solver/solver.hpp>

#include <armadillo>

namespace lsolver = lible::solver;
namespace ltests = lible::tests;

#include <cstdio> // tmp

bool lible::tests::preconditionedCG()
{
    printf("preconditionedCG\n");
    printf("\n");

    size_t dim = 100;

    arma::dmat A(dim, dim, arma::fill::randu);
    A.diag() = A.diag() * 100; // TODO: A less dipshit way of doing diagonally dominat
    arma::dmat AAt = A.t() * A;
    arma::dvec x(dim, arma::fill::randu);
    // arma::sp_mat A = arma::sprandu(dim, dim, 0.01);
    // arma::dvec A_diag(dim, arma::fill::randu);
    // for (size_t i = 0; i < dim; i++)
    //     A(i, i) = 100 * A_diag[i];
    // arma::sp_mat AAt = A.t() * A;

    arma::dvec b = AAt * x;

    // arma::dvec evals_AAt;
    // arma::dmat evecs_AAt;
    // arma::eig_sym(evals_AAt, evecs_AAt, AAt);

    // printf("\nevals_AAt:");
    // std::cout << evals_AAt << std::endl;


    lsolver::pcg_preconditioner_t precond = [&](const std::vector<double> &residual)
    {
        std::vector<double> precond_res(dim, 0.0);
        for (size_t i = 0; i < dim; i++)
            precond_res[i] = (1.0 / AAt(i, i)) * residual[i];

        return precond_res;
    };

    lsolver::pcg_sigma_t calc_sigma = [&](const std::vector<double> &trial)
    {
        return arma::conv_to<std::vector<double>>::from(
            AAt * arma::conv_to<arma::dvec>::from(trial));
    };

    std::vector<double> guess(dim, 0.0);
    std::vector<double> rhs = arma::conv_to<std::vector<double>>::from(b);

    lsolver::PCGSettings settings;
    settings.print_ = true;
    lsolver::PCGResults pcg_results = lsolver::preconditionedCG(
        guess, rhs, calc_sigma, precond, settings);

    return {};
}
