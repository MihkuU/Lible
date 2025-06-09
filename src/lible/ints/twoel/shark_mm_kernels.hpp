#pragma once

namespace lible
{
    namespace ints
    {
        template <int lbra, int lket>
        void shark_mm_ket1(const double *R, const double *ET, double *R_x_ET);

        template <int lbra, int lc, int ld>
        void shark_mm_ket2(const double *R, const double *ET, double *R_x_ET);

        template <int lbra, int lket>
        void shark_mm_bra1(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

        template <int la, int lb, int lket>
        void shark_mm_bra2(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

        template <int la, int lb, int lc, int ld>
        void shark_mm_bra2(const double *E, const double *R_x_ET, double *E_x_R_x_ET);

        void shark_mm_bra(int m, int n, int k, const double *E, const double *R_x_ET,
                          double *E_x_R_x_ET);

        void shark_mm_ket(int m, int n, int k, const double *R, const double *E, 
                          double *R_X_ET);
    }
}