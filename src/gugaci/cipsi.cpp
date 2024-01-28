#include <lible/cipsi.hpp>
#include <lible/coupling_coeffs.hpp>
#include <lible/prefix_algorithm.hpp>
#include <lible/gci_settings.hpp>

#include <omp.h>

#ifdef _USE_MPI_
#include <lible/brain.hpp>
#include <lible/gci_para.hpp>
#endif

namespace LG = lible::guga;

using namespace lible::guga;
using namespace lible::guga::util;

using arma::dvec;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

set<string> GCI::Impl::CIPSI::selectCFGsAndCSFs(wfn_ptr &wave_function)
{
    auto start_total{std::chrono::steady_clock::now()}; // tmp
    auto start{std::chrono::steady_clock::now()}; // tmp

    vector<vector<pair<string, dvec>>> generators_by_roots = generateGenerators(n_generators);
    set<string> generators;
    for (auto &item1 : generators_by_roots)
    {
        for (auto &item2 : item1)
            generators.insert(item2.first);
    }

    vector<string> prefixes_scattered = impl->prefix_algorithm->prefixBonanza(generators);
    n_prefixes = prefixes_scattered.size();

#ifdef _USE_MPI_
    Brain::comm_nodes.allreduce(mpl::plus<size_t>(), n_prefixes, n_prefixes);
#endif

    int num_threads;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    vector<vector<string>> prefixes_para(num_threads);
    for (size_t iprefix = 0; iprefix < prefixes_scattered.size(); iprefix++)            
        prefixes_para[iprefix % num_threads].push_back(prefixes_scattered[iprefix]);    

    auto end{std::chrono::steady_clock::now()};          // tmp
    std::chrono::duration<double> duration{end - start}; // tmp
    impl->t_selection_setup += duration.count();         // tmp

    // TODO: improve wording and move the comment upwards.
    /*
     * The following parallelization deserves some explanation perhaps...
     * In between HCI-selection and CIPSI-pruning, we need to collect the new spin functions and calculate
     * the coupling coefficients for subsequent steps.
     * Parallel reduction of connections here is not done because each OMP-thread has
     * its own local wave function due to the parallel prefix-algorithm.
     * That being said, we would need to have the connections and wave functions for the parallel CIPSI-pruning.
     * Therefore, we store the data in containers in 'outer' scope with the size of 'nr. of threads' so that they wont be wiped out after
     * the first 'omp parallel' region ends. The OMP threads will pick up their individual data from these containers
     * in the next parallel region for CIPSI-pruning.
     */
    vector<DataFOIS> data_fois_para(num_threads);
    vector<DataVar> data_var_para(num_threads);
    map<string, set<string>> onvs_sfs_reduced;

    start = std::chrono::steady_clock::now(); // tmp
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        DataFOIS data_fois_omp;
        DataVar data_var_omp;
        impl->prefix_algorithm->generateCFGsAndConnections(prefixes_para[thread_num],
                                                           generators_by_roots,
                                                           wave_function,
                                                           data_fois_omp,
                                                           data_var_omp);

#pragma omp critical
        {
            for (const auto &[iconf, sfs] : data_fois_omp.cfgs_sfs)
            {
                string onv = data_fois_omp.wfn->getONV(iconf);
                onvs_sfs_reduced[onv].insert(sfs.begin(), sfs.end());
            }
            data_fois_para[thread_num] = std::move(data_fois_omp);
        }
    }

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_generateCFGsAndConnections += duration.count(); // tmp    

    // auto [onvs_sfs_dia, connections_dia] = findONVsCSFsDia(generators_by_roots);
    // onvs_sfs_reduced.insert(onvs_sfs_dia.begin(), onvs_sfs_dia.begin());

    /*
     * Next is a crucial section where the new spin-functions are appended to the current ones.
     */
    start = std::chrono::steady_clock::now(); // tmp

    appendSpinFunctions(onvs_sfs_reduced,
                        impl->sfs_map__idx_to_sf,
                        impl->sfs_map__sf_to_idx);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_appendSpinFunctions += duration.count(); // tmp                            

    /*
     * Setting up the wave function for CIPSI-pruning. Essentially we transfer the HCI-wave function
     * by adding appropriate spin function positions from the global branching-diagrams map.
     */
    // WaveFunction wfn_cipsi_dia = setUpCIPSIWaveFunctionDia(onvs_sfs_dia);
    start = std::chrono::steady_clock::now(); // tmp

    vector<wfn_ptr> wfn_cipsi_para = setUpCIPSIWaveFunction(data_fois_para);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_setUpCIPSIWaveFunction += duration.count(); // tmp      

    /*
     * Constructing SF-pairs for calculating CCs
     */
    start = std::chrono::steady_clock::now(); // tmp

    sf_pair_map_1el sf_pairs_1el_cipsi;
    sf_pair_map_2el sf_pairs_2el_cipsi;
    constructSFPairs(wave_function, data_fois_para, data_var_para,
                     wfn_cipsi_para, sf_pairs_1el_cipsi, sf_pairs_2el_cipsi);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_constructSFPairs += duration.count(); // tmp                          

    /*
     * Almost there.. now calculating CCs.
     */
    start = std::chrono::steady_clock::now(); // tmp

    impl->coupling_coeffs->constructCCs(sf_pairs_1el_cipsi, sf_pairs_2el_cipsi,
                                       impl->ccs_2el, impl->ccs_1el);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_sel_constructCCs += duration.count(); // tmp                                          

    /*
     * Finally do the CIPSI-pruning....
     */
    start = std::chrono::steady_clock::now(); // tmp

    map<string, vector<int>> selected_cfgs_sfs = doCIPSIPruning(wfn_cipsi_para, data_fois_para);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_doCIPSIPruning += duration.count(); // tmp    

    n_csfs_new = 0;
    for (auto &item : selected_cfgs_sfs)
        n_csfs_new += item.second.size();

    set<string> selected_cfgs;    

    start = std::chrono::steady_clock::now(); // tmp

    appendCFGsAndCSFs(selected_cfgs_sfs, selected_cfgs, wave_function);

    end = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end - start); // tmp
    impl->t_appendCFGsAndCSFs += duration.count(); // tmp   

    auto end_total = std::chrono::steady_clock::now(); // tmp
    duration = std::chrono::duration<double>(end_total - start_total); // tmp
    impl->t_selection_total += duration.count(); // tmp   

    return selected_cfgs;
}

void GCI::Impl::CIPSI::appendCFGsAndCSFs(const map<string, vector<int>> &selected_cfgs_sfs,
                                         set<string> &selected_cfgs,
                                         wfn_ptr &wave_function)
{
    for (auto &item : selected_cfgs_sfs)
    {
        string onv = item.first;
        selected_cfgs.insert(item.first);
        CFG cfg(spin, onv);
        int nue = cfg.getNUE();
        map<int, string> sfs_map = impl->sfs_map__idx_to_sf.at(nue);
        map<string, int> sfs_map_new;
        for (const int &sf_idx : item.second)
            sfs_map_new[sfs_map.at(sf_idx)] = sf_idx;
        cfg.createCSFsFromSFs(sfs_map_new);
        wave_function->insertCFG(cfg);
    }
}

vector<vector<pair<string, dvec>>>
GCI::Impl::CIPSI::generateGenerators(size_t &n_generators)
{
    vector<vector<pair<string, dvec>>> generators_by_roots(n_roots);
    n_generators = 0;
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        vector<pair<string, dvec>> generators_iroot;

        dvec coeffs_iroot = arma::conv_to<arma::dvec>::from(impl->ci_vectors[iroot]);

        for (size_t iconf = 0; iconf < impl->wave_function->getNumCFGs(); iconf++)
        {
            size_t pos = impl->wave_function->getPos(iconf);
            size_t dim = impl->wave_function->getDim(iconf);
            string cfg = impl->wave_function->getONV(iconf);

            dvec coeffs(dim);
            for (size_t mu = 0; mu < dim; mu++)
                coeffs(mu) = coeffs_iroot(pos + mu);

            double max_coeff = arma::max(arma::abs(coeffs));

            if (max_coeff > Settings::getEpsilonGen())
            {
                generators_iroot.push_back(std::make_pair(cfg, coeffs));
                n_generators++;
            }
        }

        std::sort(generators_iroot.begin(), generators_iroot.end(),
                  [](const pair<string, dvec> &lhs, const pair<string, dvec> &rhs)
                  { return arma::max(arma::abs(lhs.second)) > arma::max(arma::abs(rhs.second)); });

        generators_by_roots[iroot] = generators_iroot;
    }

    return generators_by_roots;
}

void GCI::Impl::CIPSI::appendSpinFunctions(const map<string, set<string>> &onvs_sfs_reduced,
                                           map<int, map<int, string>> &sfs_map__idx_to_sf,
                                           map<int, map<string, int>> &sfs_map__sf_to_idx)
{
    map<size_t, set<string>> sfs_reduced;
    for (auto &item : onvs_sfs_reduced)
    {
        string cfg = item.first;
        size_t nue = count(cfg.begin(), cfg.end(), '1');
        sfs_reduced[nue].insert(item.second.begin(), item.second.end());
    }
    
#ifdef _USE_MPI_
    sfs_reduced = LG::allReduceMaps(Brain::comm_nodes, sfs_reduced);
#endif

    for (auto &item : sfs_reduced)
    {
        size_t nue = item.first;
        impl->spin_functions[nue].insert(item.second.begin(), item.second.end());

        if (impl->sfs_map__sf_to_idx.find(nue) == impl->sfs_map__sf_to_idx.end())
        {
            size_t pos = 0;
            for (auto &sf : item.second)
            {
                impl->sfs_map__sf_to_idx[nue][sf] = pos;
                impl->sfs_map__idx_to_sf[nue][pos] = sf;
                pos++;
            }
        }
        else
        {
            for (auto &sf : item.second)
            {
                if (impl->sfs_map__sf_to_idx.at(nue).find(sf) == impl->sfs_map__sf_to_idx.at(nue).end())
                {
                    size_t pos = impl->sfs_map__sf_to_idx.at(nue).size();
                    impl->sfs_map__sf_to_idx[nue][sf] = pos;
                    impl->sfs_map__idx_to_sf[nue][pos] = sf;
                }
                else
                    continue;
            }
        }
    }
}

vector<wfn_ptr> GCI::Impl::CIPSI::setUpCIPSIWaveFunction(const vector<DataFOIS> &data_fois_para)
{
    int num_threads;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    vector<wfn_ptr> wfn_cipsi_para(num_threads);

#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        const WaveFunction *wfn_hci_omp = data_fois_para[thread_num].wfn.get();
        wfn_ptr wfn_cipsi_omp = std::make_unique<WaveFunction>(spin);
        for (const auto &[icfg, sfs] : data_fois_para[thread_num].cfgs_sfs)
        {
            const CFG *cfg = wfn_hci_omp->getCFGPtr(icfg);
            int nue = cfg->getNUE();
            map<string, int> sf_map = returnSFMap(impl->sfs_map__sf_to_idx.at(nue), sfs);           

            CFG cfg_new(spin, cfg->getONV());
            cfg_new.createCSFsFromSFs(sf_map);
            wfn_cipsi_omp->insertCFG(cfg_new);
        }

#pragma omp critical
        {
            wfn_cipsi_para[thread_num] = std::move(wfn_cipsi_omp);
        }
    }

    return wfn_cipsi_para;
}

void GCI::Impl::CIPSI::constructSFPairs(const wfn_ptr &wave_function,
                                        const vector<DataFOIS> &data_fois_para,
                                        const vector<DataVar> &data_var_para,
                                        const vector<wfn_ptr> &wfn_cipsi_para,
                                        sf_pair_map_1el &sf_pairs_1el_cipsi,
                                        sf_pair_map_2el &sf_pairs_2el_cipsi)
{
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        const WaveFunction *wfn_cipsi_omp = wfn_cipsi_para[thread_num].get();
        const DataFOIS *data_fois_omp = &data_fois_para[thread_num];

        sf_pair_map_1el sf_pairs_1el_omp;
        sf_pair_map_2el sf_pairs_2el_omp;
        for (const auto &[key, connections] : data_fois_omp->connections_1el)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection);
                size_t icfg_right = get<1>(connection);
                const CFG *cfg_new = wfn_cipsi_omp->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_omp[key] = sfs_pair;
        }

        for (const auto &[key, connections] : data_fois_omp->connections_EpqErr)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection.first);
                size_t icfg_right = get<1>(connection.first);
                const CFG *cfg_new = wfn_cipsi_omp->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_omp[key].first.insert(sfs_pair.first.begin(), sfs_pair.first.end());
            sf_pairs_1el_omp[key].second.insert(sfs_pair.second.begin(), sfs_pair.second.end());
        }

        for (const auto &[key, connections] : data_fois_omp->connections_ErrEpq)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection.first);
                size_t icfg_right = get<1>(connection.first);
                const CFG *cfg_new = wfn_cipsi_omp->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_omp[key].first.insert(sfs_pair.first.begin(), sfs_pair.first.end());
            sf_pairs_1el_omp[key].second.insert(sfs_pair.second.begin(), sfs_pair.second.end());
        }

        for (const auto &[key, connections] : data_fois_omp->connections_2el)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection);
                size_t icfg_right = get<1>(connection);
                const CFG *cfg_new = wfn_cipsi_omp->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_2el_omp[key] = sfs_pair;
        }

#pragma omp critical
        {
            mergeSFPairs(sf_pairs_1el_omp, sf_pairs_1el_cipsi);
            mergeSFPairs(sf_pairs_2el_omp, sf_pairs_2el_cipsi);
        }
    }
}

map<string, vector<int>>
GCI::Impl::CIPSI::doCIPSIPruning(const vector<wfn_ptr> &wfn_cipsi_para,
                                 const vector<DataFOIS> &data_fois_para)
{
    map<string, vector<int>> selected_cfgs_sfs;
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        const WaveFunction *wfn_cipsi_omp = wfn_cipsi_para[thread_num].get();
        
        dvec diag = arma::conv_to<dvec>::from(impl->calcDiag(wfn_cipsi_omp));

        vector<arma::dvec> sigma_vectors(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            sigma_vectors[iroot] = arma::conv_to<dvec>::from(impl->calcSigma(data_fois_para[thread_num],
                                                                            impl->ci_vectors[iroot],
                                                                            wfn_cipsi_omp));

        map<string, vector<int>> selected_cfgs_sfs_omp;
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            dvec energy_minus_diag(diag.n_elem);
            energy_minus_diag.fill(impl->ci_energies[iroot]);
            energy_minus_diag -= diag;

            dvec cipsi_importance_function = arma::abs(sigma_vectors[iroot] / energy_minus_diag);
            for (size_t iconf = 0; iconf < wfn_cipsi_omp->getNumCFGs(); iconf++)
            {
                const CFG *conf_p = wfn_cipsi_omp->getCFGPtr(iconf);
                vector<int> sf_idxs = conf_p->getSFIdxs();
                size_t dim = sf_idxs.size();
                size_t pos = wfn_cipsi_omp->getPos(iconf);

                dvec cipsi_importance_function_conf = cipsi_importance_function.rows(pos, pos + dim - 1);
                if (max(cipsi_importance_function_conf) > Settings::getEpsilonVar())
                {
                    string selected_cfg = conf_p->getONV();
                    vector<int> selected_sfs;
                    for (size_t mu = 0; mu < dim; mu++)
                    {
                        size_t sf_idx = sf_idxs[mu];
                        if (cipsi_importance_function_conf(mu) > Settings::getEpsilonVar())
                            selected_sfs.push_back(sf_idx);
                    }
                    selected_cfgs_sfs_omp[selected_cfg].insert(selected_cfgs_sfs_omp[selected_cfg].end(),
                                                                 selected_sfs.begin(), selected_sfs.end());
                }
            }
        }

#pragma omp critical
        {
            selected_cfgs_sfs.merge(selected_cfgs_sfs_omp);
        }
    }

#ifdef _USE_MPI_
    selected_cfgs_sfs = LG::allReduceMaps(Brain::comm_nodes, selected_cfgs_sfs);
#endif

    return selected_cfgs_sfs;
}