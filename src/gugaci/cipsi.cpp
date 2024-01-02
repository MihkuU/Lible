#include <lible/cipsi.h>
#include <lible/coupling_coeffs.h>
#include <lible/prefix_algorithm.h>
#include <lible/gci_settings.h>

#include <omp.h>

#ifdef _USE_MPI_
#endif

using namespace lible::guga;
using namespace lible::guga::util;

using arma::dvec;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

set<string> SCI::CIPSI::selectCFGsAndCSFs(wfn_ptr &wave_function)
{
    vector<vector<pair<string, dvec>>> generators_by_roots = generateGenerators(n_generators);
    set<string> generators;
    for (auto &item1 : generators_by_roots)
    {
        for (auto &item2 : item1)
            generators.insert(item2.first);
    }

    vector<string> prefixes_scattered = sci->prefix_algorithm->prefixBonanza(generators);
    n_prefixes = prefixes_scattered.size();

#ifdef _USE_MPI_
    // mpi::all_reduce(mpi::communicator(sci->world, mpi::comm_duplicate),
    //                 n_prefixes, n_prefixes, std::plus<size_t>());
#endif                    

    int num_threads;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    vector<vector<string>> prefixes_para(num_threads);
    for (size_t iprefix = 0; iprefix < prefixes_scattered.size(); iprefix++)
    {
        size_t ipal = iprefix % num_threads;
        prefixes_para[ipal].push_back(prefixes_scattered[iprefix]);
    }

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

#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        DataFOIS data_fois_local;
        DataVar data_var_local;
        sci->prefix_algorithm->generateConfsAndConnections(prefixes_para[thread_num],
                                                           generators_by_roots,
                                                           wave_function,
                                                           data_fois_local,
                                                           data_var_local);

#pragma omp critical
        {
            for (const auto &[iconf, sfs] : data_fois_local.cfgs_sfs)
            {
                string onv = data_fois_local.wfn->getONV(iconf);
                onvs_sfs_reduced[onv].insert(sfs.begin(), sfs.end());
            }
            data_fois_para[thread_num] = std::move(data_fois_local);
        }
    }

    // auto [onvs_sfs_dia, connections_dia] = findONVsCSFsDia(generators_by_roots);
    // onvs_sfs_reduced.insert(onvs_sfs_dia.begin(), onvs_sfs_dia.begin());

    /*
     * Next is a crucial section where the new spin-functions are appended to the current ones.
     */
    appendSpinFunctions(onvs_sfs_reduced, sci->sfs_map__idx_to_sf,
                        sci->sfs_map__sf_to_idx);

    /*
     * Setting up the wave function for CIPSI-pruning. Essentially we transfer the HCI-wave function
     * by adding appropriate spin function positions from the global branching-diagrams map.
     */
    // WaveFunction wfn_cipsi_dia = setUpCIPSIWaveFunctionDia(onvs_sfs_dia);
    vector<wfn_ptr> wfn_cipsi_para = setUpCIPSIWaveFunction(data_fois_para);

    /*
     * Constructing SF-pairs for calculating CCs
     */
    sf_pair_map_1el sf_pairs_1el_cipsi;
    sf_pair_map_2el sf_pairs_2el_cipsi;
    constructSFPairs(wave_function, data_fois_para, data_var_para,
                     wfn_cipsi_para, sf_pairs_1el_cipsi, sf_pairs_2el_cipsi);

    /*
     * Almost there.. now calculating CCs.
     */
    sci->coupling_coeffs->constructCCs(sf_pairs_1el_cipsi, sf_pairs_2el_cipsi,
                                       sci->ccs_2el, sci->ccs_1el);

    /*
     * Finally do the CIPSI-pruning....
     */
    map<string, vector<int>> selected_cfgs_sfs = doCIPSIPruning(wfn_cipsi_para, data_fois_para);
    n_csfs_new = 0;
    for (auto &item : selected_cfgs_sfs)
        n_csfs_new += item.second.size();

    //     // map<string, vector<size_t>> selected_cfgs_sfs_var = doCIPSIPruningOnVarSpace();

    set<string> selected_cfgs;
    appendCFGsAndCSFs(selected_cfgs_sfs, selected_cfgs, wave_function);

    return selected_cfgs;
}

void SCI::CIPSI::appendCFGsAndCSFs(const map<string, vector<int>> &selected_cfgs_sfs,
                                   set<string> &selected_cfgs,
                                   wfn_ptr &wave_function)
{
    for (auto &item : selected_cfgs_sfs)
    {
        string onv = item.first;
        selected_cfgs.insert(item.first);
        CFG cfg(spin, onv);
        int nue = cfg.getNUE();
        map<int, string> sfs_map = sci->sfs_map__idx_to_sf.at(nue);
        map<string, int> sfs_map_new;
        for (const int &sf_idx : item.second)
            sfs_map_new[sfs_map.at(sf_idx)] = sf_idx;
        cfg.createCSFsFromSFs(sfs_map_new);
        wave_function->insertCFG(cfg);
    }
}

vector<vector<pair<string, dvec>>>
SCI::CIPSI::generateGenerators(size_t &n_generators)
{    
    vector<vector<pair<string, dvec>>> generators_by_roots(n_roots);
    n_generators = 0;
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        vector<pair<string, dvec>> generators_iroot;

        dvec coeffs_iroot = arma::conv_to<arma::dvec>::from(sci->ci_vectors[iroot]);

        for (size_t iconf = 0; iconf < sci->wave_function->getNumCFGs(); iconf++)
        {
            size_t pos = sci->wave_function->getPos(iconf);
            size_t dim = sci->wave_function->getDim(iconf);
            string cfg = sci->wave_function->getONV(iconf);

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

void SCI::CIPSI::appendSpinFunctions(const map<string, set<string>> &onvs_sfs_reduced,
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
    // mpi::all_reduce(mpi::communicator(sci->world, mpi::comm_duplicate),
    //                 sfs_reduced, sfs_reduced, addAndReturnMaps());
#endif                    

    for (auto &item : sfs_reduced)
    {
        size_t nue = item.first;
        sci->spin_functions[nue].insert(item.second.begin(), item.second.end());

        if (sci->sfs_map__sf_to_idx.find(nue) == sci->sfs_map__sf_to_idx.end())
        {
            size_t pos = 0;
            for (auto &sf : item.second)
            {
                sci->sfs_map__sf_to_idx[nue][sf] = pos;
                sci->sfs_map__idx_to_sf[nue][pos] = sf;
                pos++;
            }
        }
        else
        {
            for (auto &sf : item.second)
            {
                if (sci->sfs_map__sf_to_idx.at(nue).find(sf) == sci->sfs_map__sf_to_idx.at(nue).end())
                {
                    size_t pos = sci->sfs_map__sf_to_idx.at(nue).size();
                    sci->sfs_map__sf_to_idx[nue][sf] = pos;
                    sci->sfs_map__idx_to_sf[nue][pos] = sf;
                }
                else
                    continue;
            }
        }
    }
}

vector<wfn_ptr> SCI::CIPSI::setUpCIPSIWaveFunction(const vector<DataFOIS> &data_fois_para)
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
        const WaveFunction *wfn_hci_local = data_fois_para[thread_num].wfn.get();
        wfn_ptr wfn_cipsi_local = std::make_unique<WaveFunction>(spin);
        for (const auto &[icfg, sfs] : data_fois_para[thread_num].cfgs_sfs)
        {
            const CFG *cfg = wfn_hci_local->getCFGPtr(icfg);
            int nue = cfg->getNUE();
            map<string, int> sf_map = returnSFMap(sci->sfs_map__sf_to_idx.at(nue), sfs);

            CFG cfg_new(spin, cfg->getONV());
            cfg_new.createCSFsFromSFs(sf_map);
            wfn_cipsi_local->insertCFG(cfg_new);
        }

#pragma omp critical
        {
            wfn_cipsi_para[thread_num] = std::move(wfn_cipsi_local);
        }
    }

    return wfn_cipsi_para;
}

void SCI::CIPSI::constructSFPairs(const wfn_ptr &wave_function,
                                  const vector<DataFOIS> &data_fois_para,
                                  const vector<DataVar> &data_var_para,
                                  const vector<wfn_ptr> &wfn_cipsi_para,
                                  sf_pair_map_1el &sf_pairs_1el_cipsi,
                                  sf_pair_map_2el &sf_pairs_2el_cipsi)
{
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        const WaveFunction *wfn_cipsi_local = wfn_cipsi_para[thread_num].get();
        const DataFOIS *data_fois_local = &data_fois_para[thread_num];

        sf_pair_map_1el sf_pairs_1el_local;
        sf_pair_map_2el sf_pairs_2el_local;
        for (const auto &[key, connections] : data_fois_local->connections_1el)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection);
                size_t icfg_right = get<1>(connection);
                const CFG *cfg_new = wfn_cipsi_local->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_local[key] = sfs_pair;
        }

        for (const auto &[key, connections] : data_fois_local->connections_EpqErr)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection.first);
                size_t icfg_right = get<1>(connection.first);
                const CFG *cfg_new = wfn_cipsi_local->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_local[key].first.insert(sfs_pair.first.begin(), sfs_pair.first.end());
            sf_pairs_1el_local[key].second.insert(sfs_pair.second.begin(), sfs_pair.second.end());
        }

        for (const auto &[key, connections] : data_fois_local->connections_ErrEpq)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection.first);
                size_t icfg_right = get<1>(connection.first);
                const CFG *cfg_new = wfn_cipsi_local->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_1el_local[key].first.insert(sfs_pair.first.begin(), sfs_pair.first.end());
            sf_pairs_1el_local[key].second.insert(sfs_pair.second.begin(), sfs_pair.second.end());
        }

        for (const auto &[key, connections] : data_fois_local->connections_2el)
        {
            sfs_pair_t sfs_pair;
            for (const auto &connection : connections)
            {
                size_t icfg_new = get<0>(connection);
                size_t icfg_right = get<1>(connection);
                const CFG *cfg_new = wfn_cipsi_local->getCFGPtr(icfg_new);
                const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
                vector<int> sfs_new = cfg_new->getSFIdxs();
                vector<int> sfs_right = cfg_right->getSFIdxs();
                sfs_pair.first.insert(sfs_new.begin(), sfs_new.end());
                sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
            }
            sf_pairs_2el_local[key] = sfs_pair;
        }

#pragma omp critical
        {
            mergeSFPairs(sf_pairs_1el_local, sf_pairs_1el_cipsi);
            mergeSFPairs(sf_pairs_2el_local, sf_pairs_2el_cipsi);
        }
    }
}

map<string, vector<int>>
SCI::CIPSI::doCIPSIPruning(const vector<wfn_ptr> &wfn_cipsi_para,
                           const vector<DataFOIS> &data_fois_para)
{
    map<string, vector<int>> selected_cfgs_sfs;
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        const WaveFunction *wfn_cipsi_local = wfn_cipsi_para[thread_num].get();

        dvec diag = arma::conv_to<dvec>::from(sci->calcDiag(wfn_cipsi_local));

        vector<dvec> sigma_vectors(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            sigma_vectors[iroot] = arma::conv_to<dvec>::from(sci->calcSigma(data_fois_para[thread_num],
                                                                            sci->ci_vectors[iroot],
                                                                            wfn_cipsi_local));

        map<string, vector<int>> selected_cfgs_sfs_local;
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            dvec energy_minus_diag(diag.n_elem);
            energy_minus_diag.fill(sci->ci_energies[iroot]);
            energy_minus_diag -= diag;

            dvec cipsi_importance_function = abs(sigma_vectors[iroot] / energy_minus_diag);
            for (size_t iconf = 0; iconf < wfn_cipsi_local->getNumCFGs(); iconf++)
            {
                const CFG *conf_p = wfn_cipsi_local->getCFGPtr(iconf);
                vector<int> sf_idxs = conf_p->getSFIdxs();
                size_t dim = sf_idxs.size();
                size_t pos = wfn_cipsi_local->getPos(iconf);

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
                    selected_cfgs_sfs_local[selected_cfg].insert(selected_cfgs_sfs_local[selected_cfg].end(),
                                                                 selected_sfs.begin(), selected_sfs.end());
                }
            }
        }

#pragma omp critical
        {
            selected_cfgs_sfs.merge(selected_cfgs_sfs_local);
        }
    }

#ifdef _USE_MPI_    
    // mpi::all_reduce(mpi::communicator(sci->world, mpi::comm_duplicate),
    //                 selected_cfgs_sfs, selected_cfgs_sfs, mergeAndReturnMaps());
#endif                    

    return selected_cfgs_sfs;
}