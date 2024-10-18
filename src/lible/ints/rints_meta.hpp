#pragma once

namespace lible
{
    namespace ints
    {        
        ///// This header file might be temporary!

        /** Calculaces the index of a Cartesian Gaussian with given Cartesian exponents i, j, k.*/
        constexpr int indexCart(int i, int j, int k);

        /** Calculates the number of Hermite Gaussians for given angular momentum l. */
        constexpr int numHermites(int l);

        /**
         * Calculates the numbers of Hermite Gaussians for angular momentum 0,...,l and sums
         * them up.
         */
        constexpr int numHermitesSum(int l);

        /** */
        template <int la, int lb>
        void calcRInts(const double fac, const double p, const double *rints_buff,
                       double *rints_out);
    }
}