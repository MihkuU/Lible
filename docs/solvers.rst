.. _Linear-algebra-solvers:

Linear algebra solvers
======================

The solvers module (namespace ``lible::solver``) provides various algorithms to solve linear 
algebra problems in quantum chemistry. These methods are well known and are probably ubiquitous in 
already existing quantum chemistry. The solvers here are implemented for the specific purpose where 
the involved matrices are too large to be kept in memory. Instead, a function calculating 
a specific matrix-vector product is required

Conjugate gradient
------------------

.. cpp:struct:: CGSettings

    Structure for the convergence settings of CG and PCG solvers.

    .. cpp:var:: size_t max_iter_ = 100

       Maximum number of iterations.

    .. cpp:var:: double conv_tol_ = 1e-6

      Maximum absolute value of the residual to signal convergence.

.. cpp:struct:: CGResults 

   Structure for the results returned from CG and PCG solvers.          

   .. cpp:var:: bool converged_
      
      Flag for whether convergence was reached. 

   .. cpp:var:: size_t n_iter_

      Number of iterations that were run.

   .. cpp:var:: std::vector<double> max_abs_residuals_

      Maximum absolute values of the residuals at each iteration.

   .. cpp:var:: std::vector<double> solution_

      Obtained solution from the conjugate gradient solver. 

.. cpp:type:: cg_sigma_t = std::function<std::vector<double>(const std::vector<double> &trial)>;

    Type alias for the CG and PCG sigma vector function object.

.. cpp:type:: pcg_preconditioner_t = std::function<std::vector<double>(const std::vector<double> &residual)>;

    Type alias for the PCG preconditioner function object.

.. cpp:function:: CGResults conjugateGradient(const std::vector<double> &guess,\
     const std::vector<double> &rhs, const cg_sigma_t &calc_sigma, const CGSettings &settings = CGSettings())

   Function for running the conjugate gradient (CG) algorithm. ``guess`` is the initial guess for 
   the solution, :math:`\mathbf{x}_0`; ``rhs`` is the right-hand side of the system 
   :math:`\mathbf{A}\mathbf{x}=\mathbf{b}`; ``calc_sigma`` is a function object for calculating the 
   sigma vector; ``settings`` is a ``CGSettings`` object.

.. cpp:function:: CGResults preconditionedCG(const std::vector<double> &guess,\
    const std::vector<double> &rhs, const cg_sigma_t &calc_sigma, const pcg_preconditioner_t &preconditioner,\
    const CGSettings &settings = CGSettings())

   Function for running the preconditioned conjugate gradient (CG) algorithm. ``guess`` is the 
   initial guess for the solution, :math:`\mathbf{x}_0`; ``rhs`` is the right-hand side of the 
   system, :math:`\mathbf{A}\mathbf{x}=\mathbf{b}`; ``calc_sigma`` is a function object for 
   calculating the sigma vector; ``preconditioner`` is a function object for calculating the 
   preconditioned residual; ``settings`` is a ``CGSettings`` object.


