#pragma once

#include <armadillo>

#include "gugaci.h"

namespace Lible
{
    namespace GUGA
    {
        class GCI::CIPSI
        {
        public:
            CIPSI(GCI *gci_) : gci(gci_)
            {
                spin = gci->spin;
                n_roots = gci->n_roots;
            }

            std::set<std::string> selectCFGsAndCSFs(wfn_ptr &wave_function);

            size_t getNCSFsNew()
            {
                return n_csfs_new;
            }

            size_t getNGenerators()
            {
                return n_generators;
            }

            size_t getNPrefixes()
            {
                return n_prefixes;
            }

        private:
            GCI *gci;
            double spin;
            int n_roots;
            size_t n_csfs_new;
            size_t n_generators;
            size_t n_prefixes;

            std::vector<std::vector<std::pair<std::string, arma::dvec>>>
            generateGenerators(size_t &n_generators);

            // WaveFunction setUpCIPSIWaveFunctionDia(const map<string, set<string>> &onvs_sfs_dia);

            // pair<map<string, set<string>>, connection_map_dia> findONVsCSFsDia(const vector<vector<pair<string, dvec>>> &generators_by_roots);

            std::map<std::string, std::vector<int>>
            doCIPSIPruning(const std::vector<wfn_ptr> &wfn_cipsi_para,
                           const std::vector<DataFOIS> &data_fois_para);

            // map<string, vector<size_t>> doCIPSIPruning(const WaveFunction &wfn_cipsi_dia,
            //                                            const connection_map_dia &connections_dia);

            std::vector<wfn_ptr> setUpCIPSIWaveFunction(const std::vector<DataFOIS> &data_fois_para);

            void appendCFGsAndCSFs(const std::map<std::string, std::vector<int>> &selected_cfgs_sfs,
                                   std::set<std::string> &selected_confs,
                                   wfn_ptr &wave_function);

            void appendSpinFunctions(const std::map<std::string, std::set<std::string>> &onvs_sfs_reduced,
                                     std::map<int, std::map<int, std::string>> &sfs_map__idx_to_sf,
                                     std::map<int, std::map<std::string, int>> &sfs_map__sf_to_idx);

            void constructSFPairs(const wfn_ptr &wave_function,
                                  const std::vector<DataFOIS> &data_fois_para,
                                  const std::vector<DataVar> &data_var_para,
                                  const std::vector<wfn_ptr> &wfn_cipsi_para,
                                  sf_pair_map_1el &sf_pairs_1el_cipsi,
                                  sf_pair_map_2el &sf_pairs_2el_cipsi);
        };
    }
}