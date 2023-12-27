#include "davidson.h"
#include "davidson_settings.h"
#include "lible_util.h"

#include <armadillo>
#include <fmt/core.h>
#include <omp.h>
#include <stdexcept>

#ifdef _USE_MPI_
#endif

namespace LD = Lible::Davidson;

using namespace Lible;

using std::function;
using std::pair;
using std::vector;

pair<vector<double>, vector<vector<double>>>
LD::diagonalize(const size_t &n_roots,
                const function<vector<double>()> &calcDiag,
                const function<vector<vector<double>>(const vector<double> &diag)> &calcGuess,
                const function<vector<double>(const vector<double> &trial)> &calcSigma)
{
    // TODO: decompose it nicely into functions.
    /*
     * Implementation of the Davidson algorithm, based on Section 3.2.1 in
     * https://doi.org/10.1016/S0065-3276(08)60532-8.
     *
     */
    palPrint(fmt::format("\n   Lible::Davidson diagonalization\n\n"));

    palPrint(fmt::format("      Calculating the diagonal...                             "));
    auto start{std::chrono::steady_clock::now()};
    vector<double> diag_raw = calcDiag();
    arma::dvec diag = arma::conv_to<arma::dvec>::from(diag_raw);
    auto end(std::chrono::steady_clock::now());

    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format("done ({:.4} s)\n", duration.count()));

    palPrint(fmt::format("      Calculating the guess...                                "));
    start = std::chrono::steady_clock::now();
    vector<vector<double>> trial_vectors_raw = calcGuess(diag_raw);
    vector<arma::dvec> trial_vectors(trial_vectors_raw.size());
    for (size_t itrial = 0; itrial < trial_vectors_raw.size(); itrial++)
        trial_vectors[itrial] = arma::conv_to<arma::dvec>::from(trial_vectors_raw[itrial]);    
    end = std::chrono::steady_clock::now();    
    
    diag_raw.clear();
    trial_vectors_raw.clear();

    if (trial_vectors.size() < n_roots)
        throw std::runtime_error("\nLible::The number of guess vectors is less \
                                  than the number of roots, aborting!\n");

    duration = std::chrono::duration<double>(end - start);
    palPrint(fmt::format("done ({:.4} s)\n", duration.count()));

    std::function<arma::dvec(const arma::dvec)> calcSigmaAux = [&](const arma::dvec &trial)
    {
        return arma::conv_to<arma::dvec>::from(calcSigma(arma::conv_to<vector<double>>::from(trial)));
    };

    // TODO: Add some print here!
    vector<arma::dvec> sigma_vectors(trial_vectors.size());
    for (size_t itrial = 0; itrial < trial_vectors.size(); itrial++)
        sigma_vectors[itrial] = calcSigmaAux(trial_vectors[itrial]);

    int max_n_trial = Settings::getMaxNTrial();
    int max_n_trial_root = Settings::getMaxNTrialRoot();
    if (max_n_trial_root * n_roots > max_n_trial)
        max_n_trial = max_n_trial_root * n_roots;

    bool converged = false;
    size_t dim = sigma_vectors[0].n_elem;
    size_t n_trial = trial_vectors.size();
    vector<bool> converged_roots(n_roots, false);
    vector<double> eigvals_previous(n_roots, 0);
    std::array<vector<vector<double>>, 2> last2_eigenvectors; // For subspace collapse
    vector<arma::dvec> trial_vectors_new;

    vector<double> eigenvalues(n_roots);
    vector<vector<double>> eigenvectors(n_roots);
    for (size_t iter = 0; iter < Settings::getMaxIter(); iter++)
    {
        start = std::chrono::steady_clock::now();
        palPrint(fmt::format("      Iter {:3}\n", iter));

        if (iter > 0)
            for (size_t itrial = 0; itrial < trial_vectors_new.size(); itrial++)
                sigma_vectors.push_back(calcSigmaAux(trial_vectors_new[itrial]));

        /* Diagonalize H in the basis of trial vectors */
        int ipal = 0;
        arma::dmat G(n_trial, n_trial, arma::fill::zeros);
        for (size_t i = 0; i < n_trial; i++)
            for (size_t j = 0; j <= i; j++)
                G(i, j) = arma::dot(trial_vectors[i], sigma_vectors[j]);
        G += G.t();
        G.diag() *= 0.5;

        arma::dmat eigenvectors_G;
        arma::dvec eigenvalues_G;
        arma::eig_sym(eigenvalues_G, eigenvectors_G, G);

        /* Residuals */
        double max_residual_norm = 0;
        vector<double> residual_norms(n_roots);
        vector<arma::dvec> residual_vectors(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            arma::dvec residual(dim, arma::fill::zeros);
            for (size_t itrial = 0; itrial < n_trial; itrial++)
                residual += eigenvectors_G(itrial, iroot) *
                            (sigma_vectors[itrial] - eigenvalues_G(iroot) * trial_vectors[itrial]); // TODO preconditioner comes into play here?

            double residual_norm = arma::norm(residual);
            residual_vectors[iroot] = residual;
            residual_norms[iroot] = residual_norm;

            if (residual_norm > max_residual_norm)
                max_residual_norm = residual_norm;

            if (residual_norms[iroot] < Settings::getTolResidual())
                converged_roots[iroot] = true;
        }

        /* Check for convergence */
        double max_eigval_diff = 0;
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            double diff_eigval = eigenvalues_G(iroot) - eigvals_previous[iroot];
            eigvals_previous[iroot] = eigenvalues_G(iroot);
            if (abs(diff_eigval) > max_eigval_diff)
                max_eigval_diff = diff_eigval;
            if (iroot >= 1)
                palPrint(fmt::format("{:15}\n", " "));

            palPrint(fmt::format("   State {:3}: E = {:14.10}  (DE = {:14.10}  N(R) = {:2.10})",
                                 iroot, eigenvalues_G(iroot), diff_eigval, residual_norms[iroot]));
        }

        if (abs(max_residual_norm) < Settings::getTolResidual())
        {
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> duration{end - start};
            palPrint(fmt::format(" ({:.4} s\n)", duration.count()));
            // palPrint(fmt::format(" ({:.6 s})\n", duration.count()));
            palPrint(fmt::format("                  *** Convergence of Residuals reached ***\n"));
            // palPrint(boost::format("  (%6.6lf s)") % (t2 - t1), silent);
            // palPrint(boost::format("%1%") % "\n                  *** Convergence of Residuals reached ***\n", silent);
            converged = true;
        }

        // Make the eigenvectors and return if converged
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            arma::dvec eigenvector(dim, arma::fill::zeros);
            for (size_t itrial = 0; itrial < n_trial; itrial++)
            {
                double alphaji = eigenvectors_G(itrial, iroot);
                eigenvector = eigenvector + alphaji * trial_vectors[itrial];
            }
            eigenvalues[iroot] = eigenvalues_G(iroot);
            eigenvectors[iroot] = arma::conv_to<vector<double>>::from(eigenvector);
        }
        // We keep the CI-vectors of last 2 iterations for the expansion space collapse
        last2_eigenvectors[0] = last2_eigenvectors[1];
        last2_eigenvectors[1] = eigenvectors;

        if (converged)
            return std::make_pair(eigenvalues, eigenvectors);

        // Davidson correction vector, apply the preconditioner
        vector<arma::dvec> delta_vectors(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            if (converged_roots[iroot])
                continue;

            delta_vectors[iroot] = -arma::ones(dim) / (diag - eigenvalues_G(iroot)) % residual_vectors[iroot];
        }

        // GS orthonormalize and append the correction vectors to the trial vectors
        if (trial_vectors.size() >= max_n_trial)
        {
            // Collapse the expansion space according to the idea of Pulay that I found from
            // "Systematic Study of Selected Diagonalization Methods for Conﬁguration Interaction Matrices" (2001, Sherrill).
            trial_vectors_new.clear();
            for (size_t iroot = 0; iroot < n_roots; iroot++)
            {
                arma::dvec trial_new = eigenvectors[iroot];
                for (size_t itrial = 0; itrial < n_trial; itrial++)
                {
                    arma::dvec trial = trial_vectors[itrial];
                    trial_new -= dot(trial_new, trial) / dot(trial, trial) * trial;
                }
                double norm_trial_new = norm(trial_new);
                if (norm_trial_new > Settings::getTolDiscard())
                    trial_vectors_new.push_back((1.0 / norm_trial_new) * trial_new);
            }
            trial_vectors.clear();
            trial_vectors = trial_vectors_new;
            n_trial = trial_vectors.size();
        }
        else
        {
            trial_vectors_new.clear();
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
                if (norm_trial_new > Settings::getTolDiscard())
                {
                    trial_new = (1.0 / norm_trial_new) * trial_new;
                    trial_vectors_new.push_back(trial_new);
                    trial_vectors.push_back(trial_new);
                    n_trial++;
                }
                else
                    palPrint(fmt::format("   Trial vector is discard, norm ({:.2e}) below threshold ({:.2e})!",
                                         norm_trial_new, Settings::getTolDiscard()));
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format("done ({:.4} s)\n", duration.count()));
        // palPrint(boost::format("  (%6.6lf s)\n") % (t2 - t1), silent);
    }

    palPrint(fmt::format("\nLible::Davidson solver didn't converge in {:} iterations, aborting!\n",
                         Settings::getMaxIter()));

    // #ifdef _USE_MPI_
    //     }
    // #endif
}

pair<vector<double>, vector<vector<double>>>
LD::diagonalize(const size_t &n_roots,
                const function<vector<double>()> &calcDiag,
                const function<vector<vector<double>>(const vector<double> &diag)> &calcGuess,
                const function<vector<double>(const vector<double> &trial)> &calcSigma,
                const function<vector<double>()> &preConditioner)
{
}
