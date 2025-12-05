#pragma once

namespace lible::ints
{
    /// Specialized SHARK matrix multiplication kernel for R_x_ET with R(tuv, t'u'v') and
    /// ET(t'u'v', \nu).
    template <int lbra, int lket>
    void shark_mm_ket1(const double *R, const double *ET, double *R_x_ET);

    /// Specialized SHARK matrix multiplication kernel for R_x_ET with R(tuv, t'u'v') and
    /// ET(t'u'v', \mu\nu).
    template <int lbra, int lc, int ld>
    void shark_mm_ket2(const double *R, const double *ET, double *R_x_ET);

    /// Specialized SHARK matrix multiplication kernel for E_x_(R_x_ET) with E(\mu, tuv) and
    /// (R_x_ET)_{tuv, \nu}
    template <int lbra, int lket>
    void shark_mm_bra1(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

    /// Specialized SHARK matrix multiplication kernel for E_x_(R_x_ET) with E\mu\nu, tuv) and
    /// (R_x_ET)_{tuv, \ka}.
    template <int la, int lb, int lket>
    void shark_mm_bra2(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

    /// Specialized SHARK matrix multiplication kernel for E_x_(R_x_ET) with E(\mu\nu, tuv) and
    /// (R_x_ET)_{tuv, \ka\ta}.
    template <int la, int lb, int lc, int ld>
    void shark_mm_bra2(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

    /// Generic SHARK matrix multiplication kernel for E_x_(R_x_ET) where E is an m-by-k matrix and
    /// R_x_ET is a k_x_n matrix.
    void shark_mm_bra(int m, int n, int k, const double *E, const double *R_x_ET,
                      double *E_x_R_x_ET);

    /// Generic SHARK matrix multiplication kernel for R_x_ET where R is an m-by-k matrix and ET is
    /// a k_x_n matrix.
    void shark_mm_ket(int m, int n, int k, const double *R, const double *ET, double *R_x_ET);
}
