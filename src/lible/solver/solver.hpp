#pragma once

#include <functional>
#include <vector>

namespace lible::solver
{
    /* Preconditioned conjugate gradient */

    /// Settings for the conjugate and preconditioned conjugate gradient (CG and PCG) solver.
    struct PCGSettings
    {
        /// Flag for printing the progression.
        bool print_{false};
        /// Maximum number of iterations.
        size_t max_iter_{100};
        /// Maximum absolute value of the residual to signal convergence.
        double conv_tol_{1e-6};
    };

    /// Results returned by the conjugate and preconditioned conjugate gradient (CG and PCG) solver.
    struct PCGResults
    {
        /// Flag for whether convergence was achieved.
        bool converged_{false};
        /// Number of iterations that were run.
        size_t n_iter_{0};
        /// Maximum absolute values of the residuals at each iteration.
        std::vector<double> max_abs_residuals_;
        /// Obtained solution from the conjugate gradient solver.
        std::vector<double> solution_;
    };

    /// Type alias for the CG and PCG sigma vector function object.
    using pcg_sigma_t = std::function<std::vector<double>(
        const std::vector<double> &trial)>;
    /// Type alias for the PCG preconditioner function object.
    using pcg_preconditioner_t = std::function<std::vector<double>(
        const std::vector<double> &residual)>;

    /// Runs the preconditioned conjugate gradient algorithm. Based on Algorithm 11.5.1 from
    /// "Matrix Computations 4th edition".
    PCGResults preconditionedCG(const std::vector<double> &guess, const std::vector<double> &rhs,
                               const pcg_sigma_t &calc_sigma,
                               const pcg_preconditioner_t &preconditioner,
                               const PCGSettings &settings = PCGSettings());

    /* Pivoted Cholesky decomposition */

    /// Settings for the pivoted Cholesky Decomposition.
    struct PCDSettings
    {
        /// Value of the error below which convergence is signaled.
        double conv_tol_{};
    };

    /// Results from the pivoted Cholesky Decomposition.
    struct PCDResults
    {
        /// Flag for whether convergence was achieved.
        bool converged_{false};
        /// Number of iterations that were run
        size_t n_iter_{0};

        /// Pivot indices from the Cholesky decomposition.
        std::vector<size_t> pivot_indices_;
        /// Values of errors at each iteration.
        std::vector<double> errors_;
        /// Cholesky vectors in the lower triangular representation.
        std::vector<std::vector<double>> cholesky_vectors_;
    };

    /// Type alias for a function object that calculates/returns a matrix element A_{i,j}.
    using cd_matrix_element_t = std::function<double(size_t i, size_t j)>;

    /// Runs the pivoted Cholesky decomposition method (pCD). Based on Algorithm 1 from
    /// https://doi.org/10.1016/j.apnum.2011.10.001.
    PCDResults pivotedCD(const std::vector<double> &diagonal,
                         const cd_matrix_element_t &cd_matrix_element,
                         const PCDSettings &settings = PCDSettings());

    /* Davidson diagonalization */

    /// Settings for the Davidson diagonalization.
    struct DVDSettings
    {
        /// Flag for printing out the progression.
        bool print_{false};
        /// Convergence tolerance (max. abs. residual).
        double tol_conv_{1e-6};
        /// Tolerance for omitting new trial vectors after the Gram-Schmidt procedure.
        double tol_gs_{1e-3};
        /// Maximum number of iterations.
        size_t max_iter_{100};
        /// Maximum number of expansion vectors. The subspace is collapsed when exceeded.
        size_t max_dim_subspace_{50};
    };

    /// Results from the Davidson diagonalization.
    struct DVDResults
    {
        /// Final maximum absolute value of the residuals.
        double max_abs_res_{};
        /// Flag for whether the diagonalization converged.
        bool converged_{false};
        /// Number of iterations done.
        size_t n_iter_{};
        /// Eigenvalues.
        std::vector<double> eigenvalues_;
        /// Eigenvectors.
        std::vector<std::vector<double>> eigenvectors_;
        /// Eigenvalues from all the iterations. For each root.
        std::vector<std::vector<double>> eigenvalues_history_;
        /// Maximum absolute values of residuals from all the iterations. Largest value given
        /// among roots per iteration.
        std::vector<double> residuals_history_;
    };

    /// Type alias for the sigma vector function. Returns the sigma vectors.
    using dvd_sigma_t = std::function<std::vector<std::vector<double>>
        (const std::vector<std::vector<double>> &trial_vecs)>;

    /// Type alias for the preconditioner function. Returns the preconditioned residual vectors.
    using dvd_precondition_t = std::function<std::vector<std::vector<double>>
        (const std::vector<double> &eigenvalues,
         const std::vector<std::vector<double>> &residual_vecs)>;

    /// Runs the Davidson-Liu diagonalization algorithm. Based on the section 3.2.1 from
    /// https://doi.org/10.1016/S0065-3276(08)60532-8 and Appendix A from
    /// https://doi.org/10.1002/jcc.1111.
    DVDResults
    diagonalizeDavidson(size_t n_roots, const std::vector<std::vector<double>> &guess_vecs,
                        const dvd_precondition_t &precondition, const dvd_sigma_t &calc_sigma,
                        const DVDSettings &settings = DVDSettings());
}
