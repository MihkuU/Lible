#pragma once

namespace lible
{
    namespace davidson
    {
        class Settings
        {
        public:
            static bool getQuiet()
            {
                return quiet;
            }

            static double getTolDiscard()
            {
                return tol_discard;
            }

            static double getTolResidual()
            {
                return tol_residual;
            }

            static int getMaxIter()
            {
                return max_iter;
            }

            static int getMaxNTrial()
            {
                return max_n_trial;
            }

            static int getMaxNTrialRoot()
            {
                return max_n_trial;
            }

            static void setQuiet(const bool &q)
            {
                quiet = q;
            }

            static void setTolDiscard(const double &t)
            {
                tol_discard = t;
            }

            static void setTolResidual(const double &t)
            {
                tol_residual = t;
            }

            static void setMaxIter(const int &m)
            {
                max_iter = m;
            }

            static void setMaxNTrial(const int &m)
            {
                max_n_trial = m;
            }

            static void setMaxNTrialRoot(const int &m)
            {
                max_n_trial = m;
            }

            Settings() = delete;

        private:
            static inline bool quiet = false;
            static inline double tol_discard = 1e-7;
            static inline double tol_residual = 1e-5;
            static inline int max_iter = 50;
            static inline int max_n_trial = 500;
            static inline int max_n_trial_root = 10;
        };
    }
}