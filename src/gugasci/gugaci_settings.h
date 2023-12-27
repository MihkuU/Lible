#pragma once

#include <string>

namespace Lible
{
    namespace GUGA
    {
        class Settings
        {
            /*
             *
             */
        public:

            static double getEnergyTol()
            {
                return energy_tol;
            }

            static double getEpsilonGen()
            {
                return epsilon_gen;
            }

            static double getEpsilonVar()
            {
                return epsilon_var;
            }            

            static int getGuessDim()
            {
                return guess_dim;
            }

            static int getIterContinueEigvec()
            {
                return iter_continue_eigvec;
            }

            static int getMaxIter()
            {
                return max_iter;
            }            

            static void setEnergyTol(const double &e)
            {
                energy_tol = e;
            }

            static void setEpsilonGen(const double &e)
            {
                epsilon_gen = e;
            }

            static void setEpsilonVar(const double &e)
            {
                epsilon_var = e;
            }            

            static void setGuessDim(const int &g)
            {
                guess_dim = g;
            }

            static void setIterContinueEigvec(const int &i)
            {
                iter_continue_eigvec = i;
            }

            static void setMaxIter(const int &m)
            {
                max_iter = m;
            }

        private:
            static inline bool do_pt2 = false;
            static inline double energy_tol = 1e-5;
            static inline double epsilon_gen = 1e-2;
            static inline double epsilon_var = 1e-5;
            static inline double epsilon_pt2 = 1e-6;
            static inline int guess_dim = 512;
            static inline int iter_continue_eigvec = 4;
            static inline int max_iter = 20;
        };
    }
}