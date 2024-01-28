#include <lible/brain.hpp>
#include <lible/cipsi.hpp>
#include <lible/connections.hpp>
#include <lible/coupling_coeffs.hpp>
#include <lible/prefix_algorithm.hpp>
#include <lible/gci_impl.hpp>
#include <lible/davidson.h>
#include <lible/util.h>

#include <iostream>
#include <chrono>

#include <fmt/core.h>

namespace LD = lible::davidson;

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;
using namespace std::chrono;

using std::map;
using std::set;
using std::string;
using std::unordered_map;
using std::vector;

GCI::Impl::Impl(const int &n_orbs_, const int &n_els_,
                const int &n_roots_, const int &multiplicity_,
                const vec2d &one_el_ints_,
                const vec4d &two_el_ints_,
                const double &core_energy_)
{
#ifdef _USE_MPI_
    if (Brain::comm_nodes.is_valid())
    {
#endif
        n_orbs = n_orbs_;
        n_els = n_els_;
        n_roots = n_roots_;
        multiplicity = multiplicity_;
        one_el_ints = one_el_ints_;
        two_el_ints = two_el_ints_;
        core_energy = core_energy_;

        spin = 0.5 * (multiplicity - 1);
        min_nue = 2 * spin;
        iter_sci = 0;

        one_el_ham = one_el_ints;
        for (int p = 0; p < n_orbs; p++)
            for (int q = 0; q < n_orbs; q++)
            {
                for (int r = 0; r < n_orbs; r++)
                    one_el_ints(p, q) -= 0.5 * two_el_ints(p, r, r, q);
            }

        wave_function = std::make_unique<WaveFunction>(spin);

        cipsi = std::make_unique<CIPSI>(this);
        connections = std::make_unique<Connections>(this);
        coupling_coeffs = std::make_unique<CouplingCoeffs>(this);
        prefix_algorithm = std::make_unique<PrefixAlgorithm>(this);

        LD::Settings::setQuiet(Settings::getQuiet());

#ifdef _USE_MPI_
    }
#endif
}

GCI::Impl::~Impl() = default;

GCI::Impl::Impl(Impl &&) noexcept = default;
GCI::Impl &GCI::Impl::operator=(Impl &&) noexcept = default;

void GCI::Impl::run(vector<double> &ci_energies_out,
                    vector<vector<double>> &ci_vectors_out)
{
#ifdef _USE_MPI_
    if (Brain::comm_nodes.is_valid())
    {
#endif
        if (Settings::getQuiet())
            palPrint("\nLible::GCI calculation (Starting from HF configuration)\n");
        else 
            palPrint("\nLible::GCI calculation\n");

        auto start{std::chrono::steady_clock::now()};

        palPrint("   Starting from the HF configuration. \n", Settings::getQuiet());

        createHFConf(cfgs_new, wave_function, spin_functions,
                     sfs_map__idx_to_sf, sfs_map__sf_to_idx);

        runDriver(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies; // bcast these
        ci_vectors_out = ci_vectors; // bcast these

        auto end{std::chrono::steady_clock::now()};
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format("\nLible::GCI complete {:.2e} s\n", duration.count()));        

#ifdef _USE_MPI_        
    }
    Brain::comm_world.barrier();//tmp
#endif
}

// void GCI::Impl::runFromCFGs(const vector<string> &cfgs,
//                             vector<double> &ci_energies_out,
//                             vector<vector<double>> &ci_vectors_out)
// {
//     palPrint("\nLible::guga-Impl calculation\n\n");

//     palPrint("   Starting from the given CFGs. \n");
//     palPrint("      All CSFs per CFG will be created\n");

// #ifdef _USE_MPI_
//     if (world != MPI_COMM_NULL)
//     {
// #endif

//         runDriver(ci_energies, ci_vectors,
//                   energies_per_iter,
//                   previous_ci_coeffs_map);

//         ci_energies_out = ci_energies;
//         ci_vectors_out = ci_vectors;

// #ifdef _USE_MPI_
//     }
// #endif

//     palPrint("\nLible::guga-Impl finished\n");
// }

// void GCI::Impl::runFromCSFs(const vector<string> &csfs,
//                             vector<double> &ci_energies_out,
//                             vector<vector<double>> &ci_vectors_out)
// {
//     palPrint("\nLible::guga-Impl calculation\n\n");

//     palPrint("   Starting from the given CSFs. \n");

// #ifdef _USE_MPI_
//     if (world != MPI_COMM_NULL)
//     {
// #endif

//         runDriver(ci_energies, ci_vectors,
//                   energies_per_iter,
//                   previous_ci_coeffs_map);

//         ci_energies_out = ci_energies;
//         ci_vectors_out = ci_vectors;

// #ifdef _USE_MPI_
//     }
// #endif

//     palPrint("\nLible::guga-Impl finished\n");
// }

// void GCI::Impl::runFromCFGsFile(const string &cfgs_fname,
//                                 vector<double> &ci_energies_out,
//                                 vector<vector<double>> &ci_vectors_out)
// {
//     palPrint("\nLible::guga-Impl calculation\n\n");

//     palPrint(fmt::format("   Reading starting CFGs from a file: {}\n", cfgs_fname));

// #ifdef _USE_MPI_
//     if (world != MPI_COMM_NULL)
//     {
// #endif

//         runDriver(ci_energies, ci_vectors,
//                   energies_per_iter,
//                   previous_ci_coeffs_map);

//         ci_energies_out = ci_energies;
//         ci_vectors_out = ci_vectors;

// #ifdef _USE_MPI_
//     }
// #endif

//     palPrint("\nLible::guga-Impl finished\n");
// }

void GCI::Impl::runFromCSFsFile(const string &csfs_fname,
                                vector<double> &ci_energies_out,
                                vector<vector<double>> &ci_vectors_out)
{
#ifdef _USE_MPI_
    if (Brain::comm_nodes.is_valid())
    {
#endif
        if (Settings::getQuiet())
            palPrint(fmt::format("\nLible::GCI calculation (Starting CSFs from {})\n", csfs_fname));
        else
            palPrint("\nLible::GCI calculation\n");

        auto start{std::chrono::steady_clock::now()};

        palPrint(fmt::format("   Reading starting CSFs from a file: {}\n", csfs_fname), Settings::getQuiet());

        readCSFs(csfs_fname);

        runDriver(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

        auto end{std::chrono::steady_clock::now()};
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format("\nLible::GCI complete {:.2e} s\n", duration.count()));

        exit(1); // tmp

#ifdef _USE_MPI_
    }
#endif
}

vector<unordered_map<string, double>>
GCI::Impl::mapPreviousCIVector(const vector<vector<double>> &ci_vectors)
{
    vector<unordered_map<string, double>> previous_ci_coeffs_map(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        vector<double> ci_vector = ci_vectors[iroot];
        unordered_map<string, double> previous_ci_coeffs_map_iroot;
        for (size_t iconf = 0, idx = 0; iconf < wave_function->getNumCFGs(); iconf++)
        {
            CFG *cfg = wave_function->getCFGPtr(iconf);
            for (size_t icsf = 0; icsf < cfg->getNumCSFs(); icsf++, idx++)
            {
                string csf = cfg->getCSF(icsf);
                previous_ci_coeffs_map_iroot[csf] = ci_vector[idx];
            }
        }
        previous_ci_coeffs_map[iroot] = previous_ci_coeffs_map_iroot;
    }

    return previous_ci_coeffs_map;
}

void GCI::Impl::createHFConf(set<string> &cfgs_new,
                             wfn_ptr &wave_function,
                             map<int, set<string>> &sfs,
                             map<int, map<int, string>> &sfs_map__idx_to_sf,
                             map<int, map<string, int>> &sfs_map__sf_to_idx)
{
    string onv(n_orbs, '0');

    size_t n = n_els;
    for (size_t p = 0; p < n_orbs; p++)
        if (n > 0)
        {
            if (n > min_nue + 1)
            {
                onv[p] = '2';
                n -= 2;
            }
            else
            {
                onv[p] = '1';
                n -= 1;
            }
        }

    cfgs_new.insert(onv);

    int nue = count(onv.begin(), onv.end(), '1');
    CFG cfg(spin, onv);

    createAllSFs(nue, spin_functions, sfs_map__idx_to_sf, sfs_map__sf_to_idx);

    cfg.createCSFsFromSFs(sfs_map__sf_to_idx[nue]);
    wave_function->insertCFG(cfg);
}

void GCI::Impl::createAllSFs(const int &nue_max,
                             map<int, set<string>> &sfs,
                             map<int, map<int, string>> &sfs_map__idx_to_sf,
                             map<int, map<string, int>> &sfs_map__sf_to_idx)
{
    sfs.clear();
    sfs_map__idx_to_sf.clear();
    sfs_map__sf_to_idx.clear();

    int start;
    if (min_nue == 0)
    {
        sfs[0] = set<string>({""});
        start = 2;
    }
    else
        start = min_nue;

    for (int nue = start; nue <= nue_max;)
    {
        string blank_sf(nue, '0');
        createAllSFsRecursive(nue, '+', 0.5, 0, blank_sf, sfs);
        createAllSFsRecursive(nue, '-', -0.5, 0, blank_sf, sfs);
        nue += 2;
    }

    for (auto &item : sfs)
    {
        int nue = item.first;
        int pos = 0;
        for (const string &diagram : item.second)
        {
            sfs_map__sf_to_idx[nue][diagram] = pos;
            sfs_map__idx_to_sf[nue][pos] = diagram;
            pos++;
        }
    }
}

void GCI::Impl::createAllSFsRecursive(const int &nue, char step,
                                      double s, int i, string sf,
                                      map<int, set<string>> &sfs)
{
    sf[i] = step;

    if (s < 0)
        return;

    if ((i + 1) == nue and s == spin)
    {
        sfs[nue].insert(sf);
        return;
    }

    if (i < nue)
    {
        createAllSFsRecursive(nue, '+', s + 0.5, i + 1, sf, sfs);
        createAllSFsRecursive(nue, '-', s - 0.5, i + 1, sf, sfs);
    }
}

void GCI::Impl::readCSFs(const string &csfs_fname)
{
    std::ifstream file(csfs_fname);
    // TODO: check for redundancy in the file, also catch some runtime errors.
    // TODO: check that the spin-coupling b value is never below 0 for given CSFs

    map<string, set<string>> cfgs_csfs;
    map<int, string> cfgs;
    string csf;
    while (getline(file, csf))
    {
        auto [cfg, sf] = extractCFGandSF(csf);

        size_t nue = count(cfg.begin(), cfg.end(), '1');
        spin_functions[nue].insert(sf);
        if (sfs_map__idx_to_sf.find(nue) == sfs_map__idx_to_sf.end())
        {
            sfs_map__sf_to_idx[nue][sf] = 0;
            sfs_map__idx_to_sf[nue][0] = sf;
        }
        else
        {
            if (sfs_map__sf_to_idx[nue].find(sf) == sfs_map__sf_to_idx[nue].end())
            {
                size_t pos = sfs_map__sf_to_idx[nue].size();
                sfs_map__sf_to_idx[nue][sf] = pos;
                sfs_map__idx_to_sf[nue][pos] = sf;
            }
        }

        if (cfgs_csfs.find(cfg) == cfgs_csfs.end())
        {
            size_t pos = cfgs_csfs.size();
            cfgs_csfs[cfg].insert(csf);
            cfgs[pos] = cfg;
        }
        else
            cfgs_csfs[cfg].insert(csf);
    }

    for (auto &[pos, cfg] : cfgs)
    {
        set<string> csfs = cfgs_csfs.at(cfg);
        size_t nue = count(cfg.begin(), cfg.end(), '1');
        CFG conf(spin, cfg);
        for (const string &csf : csfs)
        {
            string sf = extractSF(csf);
            size_t sf_idx = sfs_map__sf_to_idx.at(nue).at(sf);
            conf.insertCSF(sf_idx, csf);
        }
        wave_function->insertCFG(conf);
        cfgs_new.insert(cfg);
    }
}

void GCI::Impl::runDriver(vector<double> &ci_energies,
                          vector<vector<double>> &ci_vectors,
                          vector<vector<double>> &energies_per_iter,
                          vector<unordered_map<string, double>> &previous_ci_coeffs_map)

{
    // TODO: think about splitting the function to quiet and regular
    
    auto start = std::chrono::steady_clock::now();
    auto start_quiet = start;
    string msg;

    solveCI(ci_energies,
            ci_vectors,
            energies_per_iter,
            previous_ci_coeffs_map);

    auto end_quiet = std::chrono::steady_clock::now();

    if (Settings::getQuiet())
    {
        std::chrono::duration<double> duration_quiet{end_quiet - start_quiet};
        msg = fmt::format("   {:>10}{:>14}{:>14}{:>14}{:>18}{:>10}\n",
                          "Iteration", "N (CFGs)", "N (CSFs)",
                          "N (D-iters)", "E_0 [Ha]", "t [s]"); // TODO: print out over roots in the end at quiet mode
        palPrint(msg);
        msg = fmt::format("{:>13}  {:>12}  {:>12}  {:>12}  {:>16.10f}  {:>.2e}\n",
                          0, wave_function->getNumCFGs(), wave_function->getNumCSFs(),
                          LD::Aux::n_iters, ci_energies[0], duration_quiet.count());
        palPrint(msg);
    }
    else
    {
        palPrint(fmt::format("\n      #CFGs: {:10}\n", wave_function->getNumCFGs()));
        palPrint(fmt::format("      #CSFs: {:10}\n", wave_function->getNumCSFs()));
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            palPrint(fmt::format("      E_var[root {:3}] = {:16.12}\n", iroot, ci_energies[iroot]));
    }

    for (iter_sci = 1; iter_sci <= Settings::getMaxIter(); iter_sci++)
    {
        start_quiet = std::chrono::steady_clock::now();

        palPrint(fmt::format("\n   Lible::GCI iteration {:>2}\n", iter_sci), Settings::getQuiet());

        palPrint(fmt::format("\n   Selecting CFGs and CSFs (HCI+CIPSI)...                "), Settings::getQuiet());

        auto start = std::chrono::steady_clock::now();        
        cfgs_new = cipsi->selectCFGsAndCSFs(wave_function);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format(" done {:.2e} s\n", duration.count()), Settings::getQuiet());

        size_t n_cfgs_new = cfgs_new.size(); // TODO: handle no new CSFs or CFGs
        size_t n_csfs_new = cipsi->getNCSFsNew();
        size_t n_generators = cipsi->getNGenerators();
        size_t n_prefixes = cipsi->getNPrefixes();
        msg = fmt::format("{:3}{:>6}{:>16}{:>16}{:>16}{:>16}\n", " ", "Nr:",
                          "generators", "prefixes",
                          "selected CFGs", "selected CSFs");
        palPrint(msg, Settings::getQuiet());

        msg = fmt::format("{:9}{:>16}{:>16}{:>16}{:>16}\n\n", " ",
                          n_generators, n_prefixes, n_cfgs_new, n_csfs_new);
        palPrint(msg, Settings::getQuiet());

        // palPrint(fmt::format("      Nr. of generator CFGs: {:12}\n", n_generators));
        // palPrint(fmt::format("      Nr. of CFG prefixes: {:14}\n", n_prefixes));
        // palPrint(fmt::format("      Nr. of selected CFGs: {:13}\n", n_cfgs_new));
        // palPrint(fmt::format("      Nr. of selected CSFs: {:13}\n", n_csfs_new));

        if (n_csfs_new == 0)
        {
            palPrint(fmt::format("\n   GCI converged, no new CSFs were found!\n"), Settings::getQuiet());

            if (Settings::getQuiet())
            {
                auto end_quiet = start_quiet = std::chrono::steady_clock::now();
                std::chrono::duration<double> duration_quiet{end_quiet - start_quiet};
                msg = fmt::format("{:>13}  {:>12}  {:>12}  {:>12}  {:>16.10f}  {:>.2e}\n",
                                  iter_sci, wave_function->getNumCFGs(), wave_function->getNumCSFs(),
                                  LD::Aux::n_iters, ci_energies[0], duration_quiet.count());
                palPrint(msg);
            }

            break;
        }

        solveCI(ci_energies,
                ci_vectors,
                energies_per_iter,
                previous_ci_coeffs_map);

        size_t n_cfgs = wave_function->getNumCFGs();
        size_t n_csfs = wave_function->getNumCSFs();
        palPrint(fmt::format("\n      #CFGs: {:9}\n", n_cfgs), Settings::getQuiet());
        palPrint(fmt::format("      #CSFs: {:9}\n", n_csfs), Settings::getQuiet());
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            palPrint(fmt::format("      E_var[root {:3}] = {:16.12}\n", iroot, ci_energies[iroot]),
                     Settings::getQuiet());

        arma::dvec energy_diff(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            energy_diff[iroot] = energies_per_iter[iter_sci][iroot] - energies_per_iter[iter_sci - 1][iroot];

        bool converged = false;
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            if (std::abs(energy_diff[iroot]) < Settings::getEnergyTol())
                converged = true;

        if (Settings::getQuiet())
        {
            end_quiet = std::chrono::steady_clock::now();
            std::chrono::duration<double> duration_quiet{end_quiet - start_quiet};
            msg = fmt::format("{:>13}  {:>12}  {:>12}  {:>12}  {:>16.10f}  {:>.2e}\n",
                              iter_sci, wave_function->getNumCFGs(), wave_function->getNumCSFs(),
                              LD::Aux::n_iters, ci_energies[0], duration_quiet.count());
            palPrint(msg);
        }

        if (converged)
        {
            palPrint(fmt::format("\n   GCI converged - E_diff = {:.2e} less than E_tol = {:.2e}\n",
                                 arma::min(arma::abs(energy_diff)), Settings::getEnergyTol()),
                     Settings::getQuiet());
            break;
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration{end - start};

    palPrint(fmt::format("\n      Final CFG-dimension: {:9d}\n", wave_function->getNumCFGs()), Settings::getQuiet());
    palPrint(fmt::format("      Final CSF-dimension: {:9d}\n", wave_function->getNumCSFs()), Settings::getQuiet());

    /*
     * TODO: final printout in quiet mode.
     */
    // But first, massive timing printout bonanza

    palPrint("   Performance measurements:\n");
    // palPrint(fmt::format("       t_connections_construction: {:.3e}\n", t_connections_construction));
    // palPrint(fmt::format("              t_connections_merge: {:.3e}\n", t_connections_merge));

    // palPrint(fmt::format("\n                    t_ccs_total: {:.3e}\n", t_ccs_total));
    // palPrint(fmt::format("                 t_ccs_critical: {:.3e}\n", t_ccs_critical));
    // palPrint(fmt::format("       t_returnCFGPrototypes1El: {:.3e}\n", t_returnCFGPrototypes1El));
    // palPrint(fmt::format("       t_returnCFGPrototypes2El: {:.3e}\n", t_returnCFGPrototypes2El));
    // palPrint(fmt::format("       t_returnCFGPrototypesDia: {:.3e}\n", t_returnCFGPrototypesDia));
    // palPrint(fmt::format("       t_calcCCs1El: {:.3e}\n", t_calcCCs1El));
    // palPrint(fmt::format("       t_calcCCs2El: {:.3e}\n", t_calcCCs2El));
    // palPrint(fmt::format("       t_calcCCsDia: {:.3e}\n", t_calcCCsDia));



    palPrint(fmt::format("\n                        t_solveCI: {:.3e}\n", t_solveCI));
    palPrint(fmt::format("                   t_constructCCs: {:.3e}\n", t_constructCCs));
    palPrint(fmt::format("       t_connections_construction: {:.3e}\n", t_connections_construction));
    palPrint(fmt::format("              t_connections_merge: {:.3e}\n", t_connections_merge));
    palPrint(fmt::format("                    t_diagonalize: {:.3e}\n", t_diagonalize));

    palPrint(fmt::format("\n       t_constructCCs1: {:.3e}\n", t_constructCCs1));
    palPrint(fmt::format("             t_ccs_1el: {:.3e}\n", t_ccs_1el));
    palPrint(fmt::format("             t_ccs_2el: {:.3e}\n", t_ccs_2el));
    palPrint(fmt::format("             t_ccs_dia: {:.3e}\n", t_ccs_dia));

    palPrint(fmt::format("\n                  t_selection_total: {:.3e}\n", t_selection_total));
    palPrint(fmt::format("                  t_selection_setup: {:.3e}\n", t_selection_setup));
    palPrint(fmt::format("       t_generateCFGsAndConnections: {:.3e}\n", t_generateCFGsAndConnections));
    palPrint(fmt::format("              t_appendSpinFunctions: {:.3e}\n", t_appendSpinFunctions));
    palPrint(fmt::format("           t_setUpCIPSIWaveFunction: {:.3e}\n", t_setUpCIPSIWaveFunction));
    palPrint(fmt::format("                 t_constructSFPairs: {:.3e}\n", t_constructSFPairs));
    palPrint(fmt::format("                 t_sel_constructCCs: {:.3e}\n", t_sel_constructCCs));
    palPrint(fmt::format("                   t_doCIPSIPruning: {:.3e}\n", t_doCIPSIPruning));
    palPrint(fmt::format("                t_appendCFGsAndCSFs: {:.3e}\n", t_appendCFGsAndCSFs));

    size_t n_bytes_1el_conn = 0, n_bytes_2el_conn = 0, n_bytes_dia_conn = 0;

    for (auto &[key, connections] : connections_1el)
    {
        n_bytes_1el_conn += sizeof(key);
        for (auto &item : connections)
            n_bytes_1el_conn += sizeof(item);
    }

    for (auto &[key, connections] : connections_2el)
    {
        n_bytes_2el_conn += sizeof(key);
        for (auto &item : connections)
            n_bytes_2el_conn += sizeof(item);
    }
    
    for (auto &[key, connections] : connections_dia)
    {
        n_bytes_dia_conn += sizeof(key);
        for (auto &item : connections)
            n_bytes_dia_conn += sizeof(item);
    }

    size_t n_bytes_1el_ccs = 0, n_bytes_2el_ccs = 0, n_bytes_dia_ccs = 0;

    for (auto &[key1, map1] : ccs_1el)
    {
        n_bytes_1el_ccs += sizeof(key1);
        n_bytes_1el_ccs += sizeof(map1);
        for (auto &[key2, map2] : map1)
        {
            n_bytes_1el_ccs += sizeof(key2);
            n_bytes_1el_ccs += sizeof(map2);
            for (auto &[key3, val] : map2)
            {
                n_bytes_1el_ccs += sizeof(key3);
                n_bytes_1el_ccs += sizeof(val);
            }
        }
    }

    for (auto &[key1, map1] : ccs_2el)
    {
        n_bytes_2el_ccs += sizeof(key1);
        n_bytes_2el_ccs += sizeof(map1);
        for (auto &[key2, map2] : map1)
        {
            n_bytes_2el_ccs += sizeof(key2);
            n_bytes_2el_ccs += sizeof(map2);
            for (auto &[key3, val] : map2)
            {
                n_bytes_2el_ccs += sizeof(key3);
                n_bytes_2el_ccs += sizeof(val);
            }
        }
    }

    for (auto &[key1, map1] : ccs_dia)
    {
        n_bytes_dia_ccs += sizeof(key1);
        n_bytes_dia_ccs += sizeof(map1);
        for (auto &[key2, map2] : map1)
        {
            n_bytes_dia_ccs += sizeof(key2);
            n_bytes_dia_ccs += sizeof(map2);
            for (auto &[key3, val] : map2)
            {
                n_bytes_dia_ccs += sizeof(key3);
                n_bytes_dia_ccs += sizeof(val);
            }
        }
    }

    double n_gb_1el_conn = double(n_bytes_1el_conn) / std::pow(1024, 3);
    double n_gb_2el_conn = double(n_bytes_2el_conn) / std::pow(1024, 3);
    double n_gb_dia_conn = double(n_bytes_dia_conn) / std::pow(1024, 3);

    double n_gb_1el_ccs = double(n_bytes_1el_ccs) / std::pow(1024, 3);
    double n_gb_2el_ccs = double(n_bytes_2el_ccs) / std::pow(1024, 3);
    double n_gb_dia_ccs = double(n_bytes_dia_ccs) / std::pow(1024, 3);    

    palPrint("\n   Memory footprint\n");
    palPrint(fmt::format("      n_gb_1el_conn: {:.3e} GB\n", n_gb_1el_conn));
    palPrint(fmt::format("      n_gb_2el_conn: {:.3e} GB\n", n_gb_2el_conn));
    palPrint(fmt::format("      n_gb_dia_conn: {:.3e} GB\n", n_gb_dia_conn));
    palPrint(fmt::format("      n_gb_1el_ccs: {:.3e} GB\n", n_gb_1el_ccs));
    palPrint(fmt::format("      n_gb_2el_ccs: {:.3e} GB\n", n_gb_2el_ccs));
    palPrint(fmt::format("      n_gb_dia_ccs: {:.3e} GB\n", n_gb_dia_ccs));    
}

void GCI::Impl::solveCI(vector<double> &ci_energies,
                        vector<vector<double>> &ci_vectors,
                        vector<vector<double>> &energies_per_iter,
                        vector<unordered_map<string, double>> &previous_ci_coeffs_map)
{
    auto start_total{std::chrono::steady_clock::now()}; // tmp
    palPrint("   Constructing configuration connections...              ", Settings::getQuiet());

    auto start{std::chrono::steady_clock::now()};
    connection_map_2el connections_2el_new;
    connection_map_1el connections_1el_new;
    connection_map_dia connections_dia_new;
    connections->constructConnections(cfgs_new, wave_function, wave_function,
                                      connections_1el_new, connections_2el_new,
                                      connections_dia_new);

    auto end(std::chrono::steady_clock::now());
    std::chrono::duration<double> duration{end - start}; 
    t_connections_construction += duration.count(); // TMP

    palPrint(fmt::format("done {:.2e} s\n", duration.count()), Settings::getQuiet());

    palPrint("   Constructing coupling coefficients...                  ", Settings::getQuiet());
    start = std::chrono::steady_clock::now();

    coupling_coeffs->constructCCs(connections_1el_new, connections_2el_new,
                                  connections_dia_new, wave_function, ccs_2el,
                                  ccs_1el, ccs_dia);

    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration<double>(end - start);
    t_constructCCs += duration.count(); // tmp


    palPrint(fmt::format("done {:.2e} s\n", duration.count()), Settings::getQuiet());

    start = steady_clock::now(); // TMP

    mergeConnections(connections_1el_new, connections_1el); // TODO: incorporate timing measurement
    mergeConnections(connections_2el_new, connections_2el); // TODO: incorporate timing measurement
    mergeConnections(connections_dia_new, connections_dia); // TODO: incorporate timing measurement
    
    end = steady_clock::now(); // TMP
    duration = std::chrono::duration<double>(end - start); // tmp
    t_connections_merge += duration.count(); // TMP

    // printf("\nconnections_1el.size() = %d", connections_1el.size());
    // printf("\nconnections_2el.size() = %d", connections_2el.size());
    // printf("\nconnections_dia.size() = %d", connections_dia.size());

    // TODO: printMemoryEstimate();

    /*
     * Diagonalizing the Hamiltonian
     */
    start = std::chrono::steady_clock::now();
    auto [eigenvalues, eigenvectors] = lible::davidson::diagonalize(
        n_roots,
        [this]()
        { return Impl::calcDiag(); },
        [this](const vector<double> &diag)
        { return Impl::calcGuess(diag); },
        [this](const vector<double> &trial)
        { return Impl::calcSigma(trial); });

    ci_energies = std::move(eigenvalues);
    ci_vectors = std::move(eigenvectors);

    energies_per_iter.push_back(ci_energies);

    previous_ci_coeffs_map = mapPreviousCIVector(ci_vectors);
    end = steady_clock::now(); // TMP
    duration = std::chrono::duration<double>(end - start); // tmp
    t_diagonalize += duration.count();

    auto end_total{std::chrono::steady_clock::now()}; // tmp
    duration = std::chrono::duration<double>(end_total - start_total);
    t_solveCI += duration.count();
}
