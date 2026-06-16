#include <lible/solver/solver.hpp>

#include <format>
#include <iostream>

#include <armadillo>

namespace lsolver = lible::solver;

namespace lible::solver
{
    /// Calculates the dot product between two vectors.
    double dotp(const std::vector<double> &a, const std::vector<double> &b);
}

double lsolver::dotp(const std::vector<double> &a, const std::vector<double> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

lsolver::PCGResults lsolver::preconditionedCG(
    const std::vector<double> &guess, const std::vector<double> &rhs, const pcg_sigma_t &calc_sigma,
    const pcg_preconditioner_t &preconditioner, const PCGSettings &settings)
{
    if (guess.size() != rhs.size())
        throw std::runtime_error("conjugateGradient(): `guess` and `rhs` must have the same size");

    const auto &[print, max_iter, conv_tol] = settings;

    size_t size = guess.size();

    // Initialization

    std::vector<double> sigma = calc_sigma(guess);
    std::vector<double> r(size);
    for (size_t i = 0; i < size; i++)
        r[i] = rhs[i] - sigma[i];

    std::vector<double> z = preconditioner(r);
    std::vector<double> p = z;
    std::vector<double> x = guess;

    std::vector<double> r_ = r;
    std::vector<double> z_ = z;

    // Iterations
    if (print == true)
        std::cout << std::format("   {:<14}{:<20}\n", "Iteration", "Max. |Residual|");

    std::vector<double> max_abs_residuals{r};
    bool converged = false;
    size_t iter = 0;
    for (; iter < max_iter; iter++)
    {
        // Check convergence and printout (optionally).
        double max_abs_res = 0;
        for (size_t i = 0; i < size; i++)
            if (std::fabs(r[i]) > max_abs_res)
                max_abs_res = std::fabs(r[i]);

        max_abs_residuals.push_back(max_abs_res);

        if (print == true)
            std::cout << std::format("   {:<14}{:<20.4e}\n", iter, max_abs_res);

        if (max_abs_res < conv_tol)
        {
            converged = true;
            break;
        }

        // Update the trial vector.
        if (iter > 0)
        {
            double tau = dotp(r, z) / dotp(r_, z_);
            for (size_t i = 0; i < size; i++)
                p[i] = z[i] + tau * p[i];
        }

        // Update solution and residual.
        sigma = calc_sigma(p);

        double mu = dotp(r, z) / dotp(p, sigma);

        r_ = r;
        for (size_t i = 0; i < size; i++)
        {
            x[i] -= mu * p[i];
            r[i] -= mu * sigma[i];
        }

        z_ = z;
        z = preconditioner(r);
    }

    return {converged, iter, max_abs_residuals, x};
}

lsolver::PCDResults lsolver::pivotedCD(const std::vector<double> &diagonal,
                                       const cd_matrix_element_t &cd_matrix_element,
                                       const PCDSettings &settings)
{
    // size_t dim = diagonal.size();
    //
    // arma::dvec d = conv2dvec(diagonal);
    // arma::uvec pi = arma::linspace<arma::uvec>(1, dim, dim);
    //
    // double error = arma::max(arma::abs(d));
    //
    // bool converged = false;
    // std::vector<size_t> pivot_indices;
    // std::vector<double> errors;
    // std::vector<std::vector<double>> chol_vecs;
    // size_t m = 0;
    // for (; m < dim; m++)
    // {
    //     if (error < settings.conv_tol_)
    //     {
    //         converged = true;
    //         break;
    //     }
    //
    //     size_t imax = arma::dvec(d.elem(pi)).subvec(m, dim - 1).index_max(); // TODO:
    //     std::swap(pi(m), pi(imax));
    //
    //     std::vector<double> chol_vec(dim, 0);
    //     chol_vec[pi(m)] = std::sqrt(d(pi(m)));
    //     for (size_t i = m + 1; i < dim; i++)
    //     {
    //         double l_m_pii = cd_matrix_element(pi(m), pi(i));
    //         for (size_t j = 0; j < m; j++)
    //             l_m_pii -= chol_vecs[j][pi(m)] * chol_vecs[j][pi(i)] / chol_vec[pi(m)];
    //
    //         d[pi(i)] -= std::pow(l_m_pii, 2);
    //     }
    //
    //     error = arma::max(arma::abs(arma::dvec(d.elem(pi)).subvec(m + 1, dim - 1)));
    //
    //     chol_vecs.push_back(chol_vec);
    // }
    // return {converged, m, pivot_indices, errors, chol_vecs};
}

namespace lible::solver
{
    /// Diagonalizes the Hamiltonian in the basis of the trial vectors. Eqs. (72) and (73) from
    /// https://doi.org/10.1016/S0065-3276(08)60532-8.
    std::pair<std::vector<double>, std::vector<std::vector<double>>>
    diagonalizeSubHam(const std::vector<std::vector<double>> &trial_vecs,
                      const std::vector<std::vector<double>> &sigma_vecs);

    /// Makes the eigenvectors using eq. (76) from
    /// https://doi.org/10.1016/S0065-3276(08)60532-8.
    std::vector<std::vector<double>>
    makeEigvecs(size_t n_roots, const std::vector<std::vector<double>> &eigvecs_sub,
                const std::vector<std::vector<double>> &trial_vecs);

    /// Makes the residual vectors using eq. (75) from
    /// https://doi.org/10.1016/S0065-3276(08)60532-8.
    std::vector<std::vector<double>>
    makeResiduals(size_t n_roots, const std::vector<double> &eigvals_sub,
                  const std::vector<std::vector<double>> &eigvecs_sub,
                  const std::vector<std::vector<double>> &trial_vecs,
                  const std::vector<std::vector<double>> &sigma_vecs);

    /// TODO:
    std::vector<std::vector<double>>
    makeNewTrialVecs(double tol_gs, const std::vector<bool> &converged_roots,
                     const std::vector<double> &eigenvalues,
                     const std::vector<std::vector<double>> &residuals,
                     const dvd_precondition_t &precondition,
                     std::vector<std::vector<double>> &trial_vecs);

    /// TODO:
    std::vector<std::vector<double>> collapseSubspace(std::vector<std::vector<double>> &sigma_vecs);


    std::pair<bool, double> checkConvergence(std::vector<bool> &converged_roots, double conv_tol,
                                             const std::vector<std::vector<double>> &residuals);
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> lsolver::diagonalizeSubHam(
    const std::vector<std::vector<double>> &trial_vecs,
    const std::vector<std::vector<double>> &sigma_vecs)
{
    if (trial_vecs.size() != sigma_vecs.size())
        throw std::runtime_error("diagonalizeSubHam(): Unequal number of trial/sigma vectors.");

    arma::dmat sub_ham(trial_vecs.size(), trial_vecs.size(), arma::fill::zeros);
    for (size_t i = 0; i < trial_vecs.size(); i++)
        for (size_t j = 0; j < trial_vecs.size(); j++)
        {
            const auto &trial = trial_vecs[i];
            const auto &sigma = sigma_vecs[j];
            sub_ham(i, j) = std::inner_product(trial.begin(), trial.end(), sigma.begin(), 0.0);
        }

    arma::dvec eigvals;
    arma::dmat eigvecs;
    arma::eig_sym(eigvals, eigvecs, sub_ham);

    std::vector<double> eigvals_out = arma::conv_to<std::vector<double>>::from(eigvals);
    std::vector<std::vector<double>> eigvecs_out(eigvals.size());
    for (size_t i = 0; i < eigvals.size(); i++)
        eigvecs_out[i] = arma::conv_to<std::vector<double>>::from(eigvecs.col(i));

    return {eigvals_out, eigvecs_out};
}

std::vector<std::vector<double>> lsolver::makeResiduals(
    const size_t n_roots, const std::vector<double> &eigvals_sub,
    const std::vector<std::vector<double>> &eigvecs_sub,
    const std::vector<std::vector<double>> &trial_vecs,
    const std::vector<std::vector<double>> &sigma_vecs)
{
    size_t dim = trial_vecs[0].size();

    std::vector<std::vector<double>> residuals(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        std::vector<double> residual(dim, 0);
        for (size_t itrial = 0; itrial < trial_vecs.size(); itrial++)
            for (size_t mu = 0; mu < dim; mu++)
                residual[mu] += eigvecs_sub[iroot][itrial] *
                        (sigma_vecs[itrial][mu] - eigvals_sub[iroot] * trial_vecs[itrial][mu]);

        residuals[iroot] = residual;
    }

    return residuals;
}

std::vector<std::vector<double>> lsolver::makeEigvecs(
    const size_t n_roots, const std::vector<std::vector<double>> &eigvecs_sub,
    const std::vector<std::vector<double>> &trial_vecs)
{
    size_t dim = trial_vecs[0].size();

    std::vector<std::vector<double>> eigvecs(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        std::vector<double> eigvec(dim, 0);
        for (size_t itrial = 0; itrial < trial_vecs.size(); itrial++)
            for (size_t mu = 0; mu < dim; mu++)
                eigvec[mu] += eigvecs_sub[iroot][itrial] * trial_vecs[itrial][mu];

        eigvecs[iroot] = eigvec;
    }

    return eigvecs;
}

std::vector<std::vector<double>> lsolver::makeNewTrialVecs(
    const double tol_gs, const std::vector<bool> &converged_roots,
    const std::vector<double> &eigenvalues, const std::vector<std::vector<double>> &residuals,
    const dvd_precondition_t &precondition, std::vector<std::vector<double>> &trial_vecs)
{
    std::vector<std::vector<double>> correction_vecs = precondition(eigenvalues, residuals);

    size_t dim = trial_vecs[0].size();
    std::vector<std::vector<double>> trial_vecs_new;
    for (size_t iroot = 0; iroot < correction_vecs.size(); iroot++)
    {
        if (converged_roots[iroot] == true)
            continue;

        std::vector<double> delta = correction_vecs[iroot];
        double norm = std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
        for (auto &val : delta)
            val /= norm;

        for (const auto &trial : trial_vecs)
        {
            double num = std::inner_product(trial.begin(), trial.end(), trial.begin(), 0.0);
            double den = std::inner_product(trial.begin(), trial.end(), trial.begin(), 0.0);

            for (size_t i = 0; i < dim; i++)
                delta[i] -= num / den * trial[i];
        }

        norm = std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
        for (auto &val : delta)
            val /= norm;

        if (norm < tol_gs)
            continue;

        trial_vecs_new.push_back(delta);
        trial_vecs.push_back(delta);
    }

    return trial_vecs_new;
}

std::pair<bool, double> lsolver::checkConvergence(
    std::vector<bool> &converged_roots, const double conv_tol,
    const std::vector<std::vector<double>> &residuals)
{
    double max_abs_res = 0;
    for (size_t iroot = 0; iroot < residuals.size(); iroot++)
    {
        double max_abs_res_iroot = *std::ranges::max_element(
            residuals[iroot].begin(), residuals[iroot].end(),
            [](const double a, const double b)
            {
                return std::fabs(a) < std::fabs(b);
            });

        if (max_abs_res_iroot < conv_tol)
            converged_roots[iroot] = true;

        if (max_abs_res_iroot > max_abs_res)
            max_abs_res = max_abs_res_iroot;
    }

    for (bool is_root_converged : converged_roots)
        if (is_root_converged == false)
            return {false, max_abs_res};

    return {true, max_abs_res};
}

lsolver::DVDResults lsolver::diagonalizeDavidson(
    const size_t n_roots, const std::vector<std::vector<double>> &guess_vecs,
    const dvd_precondition_t &precondition, const dvd_sigma_t &calc_sigma,
    const DVDSettings &settings)
{
    if (guess_vecs.size() < n_roots)
        throw std::runtime_error(
            "diagonalize(): `guess_vecs` must have at least `n_roots` vectors");
    // TODO: subspace collapse
    const auto &[print, tol_conv, tol_gs, max_iter, max_dim_subspace] = settings;

    std::vector<std::vector<double>> trial_vecs = guess_vecs;
    std::vector<std::vector<double>> trial_vecs_new = guess_vecs;
    std::vector<std::vector<double>> sigma_vecs;

    double max_abs_res = 0;
    bool converged = false;
    std::vector<bool> converged_roots(n_roots, false);
    std::vector<double> eigenvalues(n_roots);
    std::vector<std::vector<double>> eigenvectors;
    std::vector<std::vector<double>> eigenvalues_history;
    std::vector<double> residuals_history;

    size_t iter = 0;
    for (; iter < max_iter; iter++)
    {
        // Make new sigma vectors, diagonalize subspace Hamiltonian, make eigenvalues/eigenvectors.
        std::vector<std::vector<double>> sigma_vecs_new = calc_sigma(trial_vecs_new);
        sigma_vecs.insert(sigma_vecs.end(), sigma_vecs_new.begin(), sigma_vecs_new.end());

        const auto &[eigvals_sub, eigvecs_sub] = diagonalizeSubHam(trial_vecs, sigma_vecs);

        for (size_t iroot = 0; iroot < n_roots; iroot++)
            eigenvalues[iroot] = eigvals_sub[iroot];
        eigenvalues_history.push_back(eigenvalues);

        eigenvectors = makeEigvecs(n_roots, eigvecs_sub, trial_vecs);

        // Make residuals. Stop if converged.
        std::vector<std::vector<double>> residuals = makeResiduals(
            n_roots, eigvals_sub, eigvecs_sub, trial_vecs, sigma_vecs);

        std::tie(converged, max_abs_res) = checkConvergence(converged_roots, tol_conv, residuals);

        residuals_history.push_back(max_abs_res);

        if (converged)
            break;

        // Make new trial vectors. Collapse if the subspace size is exceeded.
        trial_vecs_new = makeNewTrialVecs(
            tol_gs, converged_roots, eigenvalues, residuals, precondition, trial_vecs);

        // TODO: collapse
    }

    return {
        max_abs_res, converged, iter, eigenvalues, eigenvectors, eigenvalues_history,
        residuals_history
    };
}
