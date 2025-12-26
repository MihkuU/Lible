#pragma once

#include <functional>
#include <vector>

namespace lible::solver
{
    // Conjugate gradient

    /// Settings for the conjugate and preconditioned conjugate gradient (CG and PCG) solver.
    struct CGSettings
    {
        /// Maximum number of iterations.
        size_t max_iter_ = 100;
        /// Maximum absolute value of the residual to signal convergence.
        double conv_tol_ = 1e-6;
    };

    /// Results returned by the conjugate and preconditioned conjugate gradient (CG and PCG) solver.
    struct CGResults
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
    using cg_sigma_t = std::function<std::vector<double>(
        const std::vector<double> &trial)>;
    /// Type alias for the PCG preconditioner function object.
    using pcg_preconditioner_t = std::function<std::vector<double>(
        const std::vector<double> &residual)>;

    /// Runs the conjugate gradient method. Based on Section 2.7.6 from "Numerical Recipes 3rd. ed."
    CGResults conjugateGradient(const std::vector<double> &guess, const std::vector<double> &rhs,
                                const cg_sigma_t &calc_sigma,
                                const CGSettings &settings = CGSettings());

    /// Runs the preconditioned conjugate gradient method. Based on Section 2.7.6 from
    /// "Numerical Recipes 3rd. ed."
    CGResults preconditionedCG(const std::vector<double> &guess, const std::vector<double> &rhs,
                               const cg_sigma_t &calc_sigma,
                               const pcg_preconditioner_t &preconditioner,
                               const CGSettings &settings = CGSettings());

    // Cholesky decomposition

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
}
