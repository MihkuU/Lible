#include <lible/davidson.hpp>
#include <lible/davidson_settings.hpp>
#include <lible/util.hpp>

#include <armadillo>
#include <fmt/core.h>
#include <omp.h>
#include <stdexcept>

#ifdef _USE_MPI_
#endif

namespace LD = lible::davidson;

using namespace lible;

using std::function;
using std::pair;
using std::string;
using std::vector;

bool checkConvergence(const double &max_res_norm,
                      const arma::dvec &eigenvalues_G,
                      const arma::dmat &eigenvectors_G,
                      const vector<double> &res_norms,
                      vector<double> &eigvals_previous)
{
    size_t n_roots = eigvals_previous.size();
    double max_eigval_diff = 0;
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        double diff_eigval = eigenvalues_G(iroot) - eigvals_previous[iroot];
        eigvals_previous[iroot] = eigenvalues_G(iroot);
        if (std::abs(diff_eigval) > max_eigval_diff)
            max_eigval_diff = diff_eigval;
        if (iroot >= 1)
            palPrint(fmt::format("{:12}\n", " "), LD::Settings::getQuiet());

        string msg = fmt::format("{:>6}  {:>18.10f}  {:>18.10f}  {:>18.10f}",
                                 iroot, eigenvalues_G(iroot), diff_eigval,
                                 res_norms[iroot]);
        palPrint(msg, LD::Settings::getQuiet());
    }

    if (max_res_norm < LD::Settings::getTolResidual())
        return true;
    

    return false;
}

auto calcResiduals(const size_t &n_roots,
                   const arma::dvec &eigenvalues_G,
                   const arma::dmat &eigenvectors_G,
                   const vector<arma::dvec> &sigma_vectors,
                   const vector<arma::dvec> &trial_vectors)
{
    size_t dim = sigma_vectors.at(0).n_elem;

    double max_res_norm = 0;
    vector<double> res_norms(n_roots);
    vector<arma::dvec> res_vectors(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        arma::dvec residual(dim, arma::fill::zeros);
        for (size_t itrial = 0; itrial < trial_vectors.size(); itrial++)
        {
            // std::cout << "\n\n" << sigma_vectors[itrial] << "\n";
            // std::cout << "\n\n" << eigenvalues_G(iroot) * trial_vectors[itrial] << "\n";

            residual += eigenvectors_G(itrial, iroot) *
                        (sigma_vectors[itrial] - eigenvalues_G(iroot) * trial_vectors[itrial]); // TODO preconditioner comes into play here?
        }

        // printf("\nresidual:\n");
        // std::cout << residual << std::endl;

        res_vectors[iroot] = residual;
        double residual_norm = arma::norm(residual);
        res_norms[iroot] = residual_norm;

        if (residual_norm > max_res_norm)
            max_res_norm = residual_norm;
    }

    // printf("\nmax_res_norm = %16.12lf\n", max_res_norm);
    // printf("\nmax_res_norm = %16.12lf\n", max_res_norm);

    return std::make_tuple(max_res_norm, res_norms, res_vectors);
}

auto diagonalizeSubHam(const vector<arma::dvec> &sigma_vectors,
                       const vector<arma::dvec> &trial_vectors)
{
    size_t n_trial = trial_vectors.size();
    arma::dmat G(n_trial, n_trial, arma::fill::zeros);
    for (size_t i = 0; i < n_trial; i++)
        for (size_t j = 0; j <= i; j++)
            G(i, j) = arma::dot(trial_vectors[i], sigma_vectors[j]);
    G += G.t();
    G.diag() *= 0.5;

    arma::dvec eigenvalues_G;
    arma::dmat eigenvectors_G;
    arma::eig_sym(eigenvalues_G, eigenvectors_G, G);

    return std::make_tuple(eigenvalues_G, eigenvectors_G);
}

auto initialize(const size_t &n_roots,
                const function<vector<double>()> &calcDiag,
                const function<vector<vector<double>>(const vector<double> &diag)> &calcGuess,
                const function<arma::dvec(const arma::dvec)> &calcSigma)
{
    palPrint(fmt::format("      Calculating the diagonal...                             "),
             LD::Settings::getQuiet());
    auto start{std::chrono::steady_clock::now()};

    vector<double> diag_raw = calcDiag();
    arma::dvec diag = arma::conv_to<arma::dvec>::from(diag_raw);

    auto end(std::chrono::steady_clock::now());
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format("done {:.2e} s\n", duration.count()),
             LD::Settings::getQuiet());

    palPrint(fmt::format("      Calculating the guess...                                "),
             LD::Settings::getQuiet());
    start = std::chrono::steady_clock::now();

    vector<vector<double>> trial_vectors_raw = calcGuess(diag_raw);
    vector<arma::dvec> trial_vectors(trial_vectors_raw.size());
    if (trial_vectors.size() < n_roots)
        throw std::runtime_error("\nLible::The number of guess vectors is less "
                                 "than the number of roots, aborting!\n");

    for (size_t itrial = 0; itrial < trial_vectors_raw.size(); itrial++)
        trial_vectors[itrial] = arma::conv_to<arma::dvec>::from(trial_vectors_raw[itrial]);

    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration<double>(end - start);
    palPrint(fmt::format("done {:.2e} s\n", duration.count()),
             LD::Settings::getQuiet());

    palPrint(fmt::format("      Calculating initial sigma vectors...                    "),
             LD::Settings::getQuiet());
    start = std::chrono::steady_clock::now();

    vector<arma::dvec> sigma_vectors(trial_vectors.size());
    for (size_t itrial = 0; itrial < trial_vectors.size(); itrial++)
        sigma_vectors[itrial] = calcSigma(trial_vectors[itrial]);

    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration<double>(end - start);
    palPrint(fmt::format("done {:.2e} s\n", duration.count()), LD::Settings::getQuiet());

    return std::make_tuple(diag, trial_vectors, sigma_vectors);
}

auto returnEigenValsVecs(const size_t &n_roots,
                         const arma::dvec &eigenvalues_G,
                         const arma::dmat &eigenvectors_G,
                         const vector<arma::dvec> &trial_vectors)
{
    size_t dim = trial_vectors.at(0).n_elem;

    vector<double> eigenvalues(n_roots);
    vector<vector<double>> eigenvectors(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        arma::dvec eigenvector(dim, arma::fill::zeros);
        for (size_t itrial = 0; itrial < trial_vectors.size(); itrial++)
        {
            double alphaji = eigenvectors_G(itrial, iroot);
            eigenvector = eigenvector + alphaji * trial_vectors[itrial];
        }
        eigenvalues[iroot] = eigenvalues_G(iroot);
        eigenvectors[iroot] = arma::conv_to<vector<double>>::from(eigenvector);
    }

    return std::make_pair(eigenvalues, eigenvectors);
}

auto createNewTrialVecs(const size_t &n_roots,
                        const arma::dvec &diag,
                        const vector<bool> &converged_roots,
                        const vector<double> &eigenvalues,
                        const vector<arma::dvec> &res_vectors,
                        const vector<arma::dvec> &trial_vectors,
                        const vector<vector<double>> &eigenvectors)
{
    size_t dim = trial_vectors.at(0).n_elem;

    int max_n_trial = LD::Settings::getMaxNTrial();
    int max_n_trial_root = LD::Settings::getMaxNTrialRoot();
    if (max_n_trial_root * n_roots > max_n_trial)
        max_n_trial = max_n_trial_root * n_roots;

    // davidson correction vector, apply the preconditioner
    vector<arma::dvec> delta_vectors(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        if (converged_roots[iroot])
            continue;

        delta_vectors[iroot] = -arma::ones(dim) / (diag - eigenvalues[iroot]) % res_vectors[iroot];
    }

    size_t n_trial = trial_vectors.size();
    vector<arma::dvec> trial_vectors_new;
    bool collapsed = false;
    if (trial_vectors.size() >= max_n_trial)
    {
        // Collapse the expansion space according to the idea of Pulay that I found from
        // "Systematic Study of Selected Diagonalization Methods for Conﬁguration Interaction Matrices" (2001, Sherrill).
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            arma::dvec trial_new = arma::conv_to<arma::dvec>::from(eigenvectors[iroot]);
            for (size_t itrial = 0; itrial < n_trial; itrial++)
            {
                arma::dvec trial = trial_vectors[itrial];
                trial_new -= dot(trial_new, trial) / dot(trial, trial) * trial;
            }
            double norm_trial_new = norm(trial_new);
            if (norm_trial_new > LD::Settings::getTolDiscard())
                trial_vectors_new.push_back((1.0 / norm_trial_new) * trial_new);
        }
        collapsed = true;
    }
    else
    {
        /*
         * GS orthonormalize and append the delta vectors to the list of trial vectors
         */
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            if (converged_roots[iroot])
                continue;

            arma::dvec trial_new = delta_vectors[iroot];
            arma::dvec trial_new_(dim, arma::fill::zeros);
            for (size_t itrial = 0; itrial < n_trial; itrial++)
            {
                arma::dvec trial = trial_vectors[itrial];
                trial_new_ -= arma::dot(trial_new, trial) / arma::dot(trial, trial) * trial;
            }
            trial_new += trial_new_;

            double norm_trial_new = arma::norm(trial_new);
            if (norm_trial_new > LD::Settings::getTolDiscard())
            {
                trial_new = (1.0 / norm_trial_new) * trial_new;
                trial_vectors_new.push_back(trial_new);
                n_trial++;
            }
            else
                palPrint(fmt::format("   Trial vector is discarded, norm ({:.2e}) below threshold ({:.2e})!",
                                     norm_trial_new, LD::Settings::getTolDiscard()), LD::Settings::getQuiet());
        }
    }

    return std::pair(trial_vectors_new, collapsed);
}

pair<vector<double>, vector<vector<double>>>
LD::diagonalize(const size_t &n_roots,
                const function<vector<double>()> &calcDiag,
                const function<vector<vector<double>>(const vector<double> &diag)> &calcGuess,
                const function<vector<double>(const vector<double> &trial)> &calcSigma)
{
    /*
     * TODO: test multiroot
     */
    auto start{std::chrono::steady_clock::now()};

    /*
     * Implementation of the Davidson algorithm, based on Section 3.2.1 in
     * https://doi.org/10.1016/S0065-3276(08)60532-8.
     */
    palPrint(fmt::format("\n   Lible::Davidson diagonalization\n\n"), LD::Settings::getQuiet());

    std::function<arma::dvec(const arma::dvec)> calcSigmaAux = [&](const arma::dvec &trial)
    {
        return arma::conv_to<arma::dvec>::from(calcSigma(arma::conv_to<vector<double>>::from(trial)));
    };

    auto [diag, trial_vectors, sigma_vectors] = initialize(n_roots,
                                                           calcDiag,
                                                           calcGuess,
                                                           calcSigmaAux);

    if (Settings::getMaxIter() == 0)
    {
        string msg = fmt::format("Lible::Davidson max_iter is set to {}! "
                                 "Increase the maximum number of iterations:\n"
                                 "      lible::davidson::setMaxIter()\n",
                                 Settings::getMaxIter(), LD::Settings::getQuiet());

        throw std::runtime_error(msg);
    }

    bool converged = false;
    vector<bool> converged_roots(n_roots, false);
    vector<double> eigvals_previous(n_roots, 0);
    std::array<vector<vector<double>>, 2> last2_eigenvectors; // For subspace collapse
    std::pair<vector<double>, vector<vector<double>>> eig_vals_vecs;

    string msg = fmt::format("\n{:6}{:>6}{:>6}{:>20}{:>20}{:>20}{:>20}\n",
                             " ", "Iter", "Root", "E [Ha]",
                             "E-diff", "R-norm", "t [s]");
    palPrint(msg, LD::Settings::getQuiet());

    size_t iter = 0;
    for (; iter < Settings::getMaxIter(); iter++)
    {
        auto start{std::chrono::steady_clock::now()};
        palPrint(fmt::format("{:6}{:>6}", " ", (iter + 1)), LD::Settings::getQuiet());

        /* Diagonalize H in the basis of trial vectors */
        auto [eigenvalues_G, eigenvectors_G] = diagonalizeSubHam(sigma_vectors,
                                                                 trial_vectors);
        /* Calculate residuals */
        auto [max_res_norm, res_norms, res_vectors] = calcResiduals(n_roots,
                                                                    eigenvalues_G,
                                                                    eigenvectors_G,
                                                                    sigma_vectors,
                                                                    trial_vectors);

        /* Check for converged roots */
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            if (res_norms[iroot] < Settings::getTolResidual())
                converged_roots[iroot] = true;

        /* Check for overall convergence */
        converged = checkConvergence(max_res_norm,
                                     eigenvalues_G,
                                     eigenvectors_G,
                                     res_norms,
                                     eigvals_previous);

        /*
         * Create the eigenvectors and return if converged.
         * The CI vectors of the last 2 iterations are kept for the expansion space collapse.
         */
        auto [eigenvalues, eigenvectors] = returnEigenValsVecs(n_roots,
                                                               eigenvalues_G,
                                                               eigenvectors_G,
                                                               trial_vectors);

        if (converged)
        {
            eig_vals_vecs = std::make_pair(eigenvalues, eigenvectors);

            auto end(std::chrono::steady_clock::now());
            std::chrono::duration<double> duration{end - start};
            msg = fmt::format("{:>20.2e}\n", duration.count());
            palPrint(msg, LD::Settings::getQuiet());

            palPrint(fmt::format("{:6}{:^92s}\n", " ",
                                 "*** Convergence of Residuals reached ***"),
                     LD::Settings::getQuiet());

            break;
        }

        last2_eigenvectors[0] = last2_eigenvectors[1];
        last2_eigenvectors[1] = eigenvectors;

        /*
         * Gram-Schmidt orthonormalize and append the correction vectors to the trial vectors.
         * Reset the space of trial vectors if collapsing.
         */
        auto [trial_vectors_new, collapsed] = createNewTrialVecs(n_roots,
                                                                 diag,
                                                                 converged_roots,
                                                                 eigenvalues,
                                                                 res_vectors,
                                                                 trial_vectors,
                                                                 eigenvectors);

        if (collapsed)
            trial_vectors = trial_vectors_new;
        else
            for (auto &trial : trial_vectors_new)
                trial_vectors.push_back(trial);

        /* Calculate new sigma vectors */
        for (size_t itrial = 0; itrial < trial_vectors_new.size(); itrial++)
            sigma_vectors.push_back(calcSigmaAux(trial_vectors_new[itrial]));

        auto end(std::chrono::steady_clock::now());
        std::chrono::duration<double> duration{end - start};
        msg = fmt::format("{:>20.2e}\n", duration.count());
        palPrint(msg, LD::Settings::getQuiet());
    }

    auto end(std::chrono::steady_clock::now());

    if (converged)
    {
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format("\n   Lible::Davidson complete {:.2e} s\n", duration.count()),
                 LD::Settings::getQuiet());
        
        Aux::n_iters = iter + 1;
        return eig_vals_vecs;
    }
    else
    {
        msg = fmt::format("Lible::Davidson didn't converge in {} iterations! "
                          "Increase the maximum number of iterations:\n"
                          "      lible::davidson::setMaxIter()\n",
                          iter);

        Aux::n_iters = iter + 1;

        throw std::runtime_error(msg);
    }
}

// pair<vector<double>, vector<vector<double>>>
// LD::diagonalize(const size_t &n_roots,
//                 const function<vector<double>()> &calcDiag,
//                 const function<vector<vector<double>>(const vector<double> &diag)> &calcGuess,
//                 const function<vector<double>(const vector<double> &trial)> &calcSigma,
//                 const function<vector<double>()> &preConditioner)
// {
// }
