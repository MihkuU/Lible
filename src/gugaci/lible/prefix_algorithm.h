#pragma once

#include <armadillo>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <lible/gci_impl.hpp>

namespace lible
{
    namespace guga
    {
        class GCI::Impl::PrefixAlgorithm
        {
        public:
            /*
             * Wrapper struct for hiding the ugly details of the prefix algorithm.
             */
            PrefixAlgorithm(Impl *impl_) : impl(impl_)
            {
                spin = impl->spin;
                min_nue = impl->min_nue;
                n_orbs = impl->n_orbs;

                int n_doubly_occ = (impl->n_els - min_nue) / 2;
                prefix_size = min_nue + n_doubly_occ;

                constructSortedIntegralLists(impl->one_el_ints, impl->two_el_ints,
                                             max_abs_1el_element, max_abs_2el_element);
            }

            std::vector<std::string> prefixBonanza(const std::set<std::string> &cfgs);

            void generateCFGsAndConnections(const std::vector<std::string> &prefixes,
                                            const std::vector<std::vector<
                                                std::pair<std::string, arma::dvec>>> &generators_by_roots,
                                            const wfn_ptr &wfn_right,
                                            DataFOIS &data_fois, DataVar &data_var);

            std::vector<std::string> findConnectedSFs1El(const std::vector<std::string> &sfs_right,
                                                         const quintet &info_cc);

        private:
            using ints_11_map = std::map<int, std::vector<std::tuple<int, double>>>;
            using ints_12_map = std::map<int, std::vector<std::tuple<int, int, double>>>;
            using ints_21_map = std::map<std::pair<int, int>, std::vector<std::tuple<int, double>>>;
            using ints_13_map = std::map<int, std::vector<std::tuple<int, int, int, double>>>;
            using ints_22_map = std::map<std::pair<int, int>, std::vector<std::tuple<int, int, double>>>;
            using ints_31_map = std::map<std::tuple<int, int, int>, std::vector<std::tuple<int, double>>>;

            Impl *impl;

            double max_abs_1el_element;
            double max_abs_2el_element;
            double spin;
            size_t min_nue;
            size_t n_orbs;
            size_t prefix_size;

            // Epq
            ints_11_map intlist_Epq_pqOut;
            ints_11_map intlist_Epq_pIn_qOut;
            ints_11_map intlist_Epq_qIn_pOut;
            // EpqErr
            ints_12_map intlist_EpqErr_pqrOut;
            ints_12_map intlist_EpqErr_rIn_pqOut;
            ints_12_map intlist_EpqErr_qIn_pOut;
            ints_12_map intlist_EpqErr_pIn_qOut;
            ints_21_map intlist_EpqErr_pqIn;
            // ErrEpq
            ints_12_map intlist_ErrEpq_pqrOut;
            ints_12_map intlist_ErrEpq_rIn_pqOut;
            ints_12_map intlist_ErrEpq_qIn_pOut;
            ints_12_map intlist_ErrEpq_pIn_qOut;
            ints_21_map intlist_ErrEpq_pqIn;
            // EpqEqr
            ints_12_map intlist_EpqEqr_pqrOut;
            ints_12_map intlist_EpqEqr_rIn_pOut;
            ints_12_map intlist_EpqEqr_pIn_rOut;
            ints_21_map intlist_EpqEqr_prIn;
            // EpqErp
            ints_12_map intlist_EpqErp_pqrOut;
            ints_12_map intlist_EpqErp_qIn_rOut;
            ints_12_map intlist_EpqErp_rIn_qOut;
            ints_21_map intlist_EpqErp_qrIn;
            // EpqEpq
            ints_11_map intlist_EpqEpq_pqOut;
            ints_11_map intlist_EpqEpq_qIn_pOut;
            ints_11_map intlist_EpqEpq_pIn_qOut;
            // EpqEpr
            ints_12_map intlist_EpqEpr_pqrOut;
            ints_12_map intlist_EpqEpr_qIn_prOut;
            ints_12_map intlist_EpqEpr_pIn_qrOut;
            ints_21_map intlist_EpqEpr_qrIn_pOut;
            ints_21_map intlist_EpqEpr_pqIn_rOut;
            // EpqErq
            ints_12_map intlist_EpqErq_pqrOut;
            ints_12_map intlist_EpqErq_pIn_qrOut;
            ints_12_map intlist_EpqErq_qIn_prOut;
            ints_21_map intlist_EpqErq_prIn_qOut;
            ints_21_map intlist_EpqErq_pqIn_rOut;
            // EpqErs
            ints_13_map intlist_EpqErs_sIn_pqrOut;
            ints_13_map intlist_EpqErs_qIn_prsOut;
            ints_13_map intlist_EpqErs_pIn_qrsOut;
            ints_22_map intlist_EpqErs_pqrsOut;
            ints_22_map intlist_EpqErs_qsIn_prOut;
            ints_22_map intlist_EpqErs_prIn_qsOut;
            ints_22_map intlist_EpqErs_pqIn_rsOut;
            ints_22_map intlist_EpqErs_psIn_qrOut;
            ints_31_map intlist_EpqErs_pqsIn_rOut;
            ints_31_map intlist_EpqErs_prsIn_qOut;
            ints_31_map intlist_EpqErs_pqrIn_sOut;

            void constructSortedIntegralLists(const vec2d &one_el_ints,
                                              const vec4d &two_el_ints,
                                              double &max_abs_1el_element_out,
                                              double &max_abs_2el_element_out);

            void integralListsHelper_Epq(const vec2d &one_el_ints,
                                         ints_11_map &intlist_Epq_pqOut,
                                         ints_11_map &intlist_Epq_pIn_qOut,
                                         ints_11_map &intlist_Epq_qIn_pOut);

            void integralListsHelper_EpqErr(const vec4d &two_el_ints,
                                            ints_12_map &intlist_EpqErr_pqrOut,
                                            ints_12_map &intlist_EpqErr_rIn_pqOut,
                                            ints_12_map &intlist_EpqErr_qIn_pOut,
                                            ints_12_map &intlist_EpqErr_pIn_qOut,
                                            ints_21_map &intlist_EpqErr_pqIn);

            void integralListsHelper_ErrEpq(const vec4d &two_el_ints,
                                            ints_12_map &intlist_ErrEpq_pqrOut,
                                            ints_12_map &intlist_ErrEpq_rIn_pqOut,
                                            ints_12_map &intlist_ErrEpq_qIn_pOut,
                                            ints_12_map &intlist_ErrEpq_pIn_qOut,
                                            ints_21_map &intlist_ErrEpq_pqIn);

            void integralListsHelper_EpqEqr(const vec4d &two_el_ints,
                                            ints_12_map &intlist_EpqEqr_pqrOut,
                                            ints_12_map &intlist_EpqEqr_rIn_pOut,
                                            ints_12_map &intlist_EpqEqr_pIn_rOut,
                                            ints_21_map &intlist_EpqEqr_prIn);

            void integralListsHelper_EpqErp(const vec4d &two_el_ints,
                                            ints_12_map &intlist_EpqErp_pqrOut,
                                            ints_12_map &intlist_EpqErp_qIn_rOut,
                                            ints_12_map &intlist_EpqErp_rIn_qOut,
                                            ints_21_map &intlist_EpqErp_qrIn);

            void integralListsHelper_EpqEpq(const vec4d &two_el_ints,
                                            ints_11_map &intlist_EpqEpq_pqOut,
                                            ints_11_map &intlist_EpqEpq_qIn_pOut,
                                            ints_11_map &intlist_EpqEpq_pIn_qOut);

            void integralListsHelper_EpqEpr(const vec4d &two_el_ints,
                                            ints_12_map &intlist_EpqEpr_pqrOut,
                                            ints_12_map &intlist_EpqEpr_qIn_prOut,
                                            ints_12_map &intlist_EpqEpr_pIn_qrOut,
                                            ints_21_map &intlist_EpqEpr_qrIn_pOut,
                                            ints_21_map &intlist_EpqEpr_pqIn_rOut);

            void integralListsHelper_EpqErq(const vec4d &two_el_ints,
                                            ints_12_map &intlist_EpqErq_pqrOut,
                                            ints_12_map &intlist_EpqErq_pIn_qrOut,
                                            ints_12_map &intlist_EpqErq_qIn_prOut,
                                            ints_21_map &intlist_EpqErq_prIn_qOut,
                                            ints_21_map &intlist_EpqErq_pqIn_rOut);

            void integralListsHelper_EpqErs(const vec4d &two_el_ints,
                                            ints_13_map &intlist_EpqErs_sIn_pqrOut,
                                            ints_13_map &intlist_EpqErs_qIn_prsOut,
                                            ints_13_map &intlist_EpqErs_pIn_qrsOut,
                                            ints_22_map &intlist_EpqErs_pqrsOut,
                                            ints_22_map &intlist_EpqErs_qsIn_prOut,
                                            ints_22_map &intlist_EpqErs_prIn_qsOut,
                                            ints_22_map &intlist_EpqErs_pqIn_rsOut,
                                            ints_22_map &intlist_EpqErs_psIn_qrOut,
                                            ints_31_map &intlist_EpqErs_pqsIn_rOut,
                                            ints_31_map &intlist_EpqErs_prsIn_qOut,
                                            ints_31_map &intlist_EpqErs_pqrIn_sOut);

            void innerPrefixHelper1El(const int &p, const int &q, const int &icfg_right,
                                      const size_t &nue_right, const std::string &onv_right,
                                      const std::vector<std::string> &sfs_right,
                                      const wfn_ptr &wfn_right,
                                      DataFOIS &data_fois, DataVar &data_var);

            void innerPrefixHelper2El_EpqErr(const int &p, const int &q, const int &r,
                                             const int &icfg_right, const size_t &nue_right,
                                             const std::string &onv_right,
                                             const std::vector<std::string> &sfs_right,
                                             const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                             DataVar &data_var);

            void innerPrefixHelper2El_ErrEpq(const int &p, const int &q, const int &r,
                                             const int &icfg_right, const size_t &nue_right,
                                             const std::string &onv_right,
                                             const std::vector<std::string> &sfs_right,
                                             const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                             DataVar &data_var);

            void innerPrefixHelper2El(const int &p, const int &q, const int &r,
                                      const int &s, const int &icfg_right,
                                      const size_t &nue_right, const std::string &onv_right,
                                      const std::vector<std::string> &sfs_right,
                                      const wfn_ptr &wfn_right,
                                      DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_Epq(const double &max_ci_coeff, const int &na, const int &nc,
                                  const int &occ_diff_sum, const int &p, const int &q,
                                  const int &icfg_right, const std::string &onv_right,
                                  const std::vector<std::string> &sfs_right,
                                  const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                  DataVar &data_var);

            void prefixHelper_EpqErr(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right, const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_ErrEpq(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right, const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqEqr(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right, const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqErp(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right, const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqEpq(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right, const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqEpr(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right,
                                     const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqErq(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right,
                                     const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);

            void prefixHelper_EpqErs(const double &max_ci_coeff, const int &na, const int &nc,
                                     const int &occ_diff_sum, const int &icfg_right,
                                     const std::string &onv_right,
                                     const std::vector<int> &a_idxs,
                                     const std::vector<int> &c_idxs,
                                     const std::vector<std::string> &sfs_right,
                                     const wfn_ptr &wfn_right,
                                     DataFOIS &data_fois, DataVar &data_var);
        };
    }
}