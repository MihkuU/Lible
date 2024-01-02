#include <lible/cipsi.h>
#include <lible/connections.h>
#include <lible/coupling_coeffs.h>
#include <lible/prefix_algorithm.h>

#include <lible/guga_sci.h>

#include <davidson.h>
#include <lible/util.h>

#include <chrono>

#include <fmt/core.h>

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::map;
using std::set;
using std::string;
using std::unordered_map;
using std::vector;

SCI::SCI(const int &n_orbs_, const int &n_els_,
         const int &n_roots_, const int &multiplicity_,
         const vector<double> &one_el_ints_,
         const vector<double> &two_el_ints_)
    : n_orbs(n_orbs_), n_els(n_els_),
      n_roots(n_roots_), multiplicity(multiplicity_),
      one_el_ints(one_el_ints_), two_el_ints(two_el_ints_)
{
    spin = 0.5 * (multiplicity - 1);
    min_nue = 2 * spin;
}

SCI::~SCI() = default;

SCI::SCI(SCI &&) noexcept = default;
SCI &SCI::operator=(SCI &&) noexcept = default;

void SCI::run(vector<double> &ci_energies_out,
              vector<vector<double>> &ci_vectors_out)
{
    palPrint("\nLible::guga-SCI calculation\n\n");

    palPrint("   Starting from the HF-configuration... \n");

    createHFConf(cfgs_new, wave_function, spin_functions,
                 sfs_map__idx_to_sf, sfs_map__sf_to_idx);

#ifdef _USE_MPI_
    if (world != MPI_COMM_NULL)
    {
#endif

        runKernel(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

#ifdef _USE_MPI_
    }
#endif

    palPrint("\nLible::guga-SCI finished\n");
}

void SCI::runFromCFGs(const vector<string> &cfgs,
                      vector<double> &ci_energies_out,
                      vector<vector<double>> &ci_vectors_out)
{
    palPrint("\nLible::guga-SCI calculation\n\n");

    palPrint("   Starting from the given CFGs... \n");
    palPrint("      All CSFs per CFG will be created\n");

#ifdef _USE_MPI_
    if (world != MPI_COMM_NULL)
    {
#endif

        runKernel(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

#ifdef _USE_MPI_
    }
#endif

    palPrint("\nLible::guga-SCI finished\n");
}

void SCI::runFromCSFs(const vector<string> &csfs,
                      vector<double> &ci_energies_out,
                      vector<vector<double>> &ci_vectors_out)
{
    palPrint("\nLible::guga-SCI calculation\n\n");

    palPrint("   Starting from the given CSFs... \n");

#ifdef _USE_MPI_
    if (world != MPI_COMM_NULL)
    {
#endif

        runKernel(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

#ifdef _USE_MPI_
    }
#endif

    palPrint("\nLible::guga-SCI finished\n");
}

void SCI::runFromCFGsFile(const string &cfgs_fname,
                          vector<double> &ci_energies_out,
                          vector<vector<double>> &ci_vectors_out)
{
    palPrint("\nLible::guga-SCI calculation\n\n");

    palPrint(fmt::format("   Reading starting CFGs from a file: {}\n", cfgs_fname));

#ifdef _USE_MPI_
    if (world != MPI_COMM_NULL)
    {
#endif

        runKernel(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

#ifdef _USE_MPI_
    }
#endif

    palPrint("\nLible::guga-SCI finished\n");
}

void SCI::runFromCSFsFile(const string &csfs_fname,
                          vector<double> &ci_energies_out,
                          vector<vector<double>> &ci_vectors_out)
{
    palPrint("\nLible::guga-SCI calculation\n\n");

    palPrint(fmt::format("   Readinf starting CSFs from a file: {}\n", csfs_fname));

#ifdef _USE_MPI_
    if (world != MPI_COMM_NULL)
    {
#endif

        runKernel(ci_energies, ci_vectors,
                  energies_per_iter,
                  previous_ci_coeffs_map);

        ci_energies_out = ci_energies;
        ci_vectors_out = ci_vectors;

#ifdef _USE_MPI_
    }
#endif

    palPrint("\nLible::guga-SCI finished\n");
}

vector<unordered_map<string, double>>
SCI::mapPreviousCIVector(const vector<vector<double>> &ci_vectors)
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

void SCI::createHFConf(set<string> &cfgs_new,
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

void SCI::createAllSFs(const int &nue_max,
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

void SCI::createAllSFsRecursive(const int &nue, char step,
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

void SCI::runKernel(vector<double> &ci_energies,
                    vector<vector<double>> &ci_vectors,
                    vector<vector<double>> &energies_per_iter,
                    vector<unordered_map<string, double>> &previous_ci_coeffs_map)

{
    auto start = std::chrono::steady_clock::now();

    solveCI(ci_energies,
            ci_vectors,
            energies_per_iter,
            previous_ci_coeffs_map);

    palPrint(fmt::format("      #CFGs: {:10}\n", wave_function->getNumCFGs()));
    palPrint(fmt::format("      #CSFs: {:10}\n", wave_function->getNumCSFs()));
    // TODO: print energies

    for (iter_sci = 1; iter_sci <= Settings::getMaxIter(); iter_sci++)
    {
        palPrint(fmt::format("\n============================|  SCI iteration %3d  |============================\n", iter_sci));

        palPrint(fmt::format("\n   Selecting CFGs and CSFs (HCI+CIPSI)...           "));

        auto start = std::chrono::steady_clock::now();
        cipsi->selectCFGsAndCSFs(wave_function);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration{end - start};
        palPrint(fmt::format("\n   selection done ({:6.4} s)\n", duration.count()));

        size_t n_cfgs_new = cfgs_new.size(); // TODO: handle no new CSFs or CFGs
        size_t n_csfs_new = cipsi->getNCSFsNew();
        size_t n_generators = cipsi->getNGenerators();
        size_t n_prefixes = cipsi->getNPrefixes();
        palPrint(fmt::format("      Nr. of generator CFGs: {:12}\n", n_generators));
        palPrint(fmt::format("      Nr. of CFG prefixes: {:14}\n", n_prefixes));
        palPrint(fmt::format("      Nr. of selected CFGs: {:13}\n", n_cfgs_new));
        palPrint(fmt::format("      Nr. of selected CSFs: {:13}\n", n_csfs_new));

        if (n_csfs_new == 0)
        {
            palPrint(fmt::format("\n   SCI converged, no new CSFs were found!\n"));
            break;
        }

        solveCI(ci_energies, ci_vectors, energies_per_iter,
                previous_ci_coeffs_map);

        size_t n_cfgs = wave_function->getNumCFGs();
        size_t n_csfs = wave_function->getNumCSFs();
        palPrint(fmt::format("         #CFGs: {:9}\n", n_cfgs));
        palPrint(fmt::format("         #CSFs: {:9}\n", n_csfs));
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            palPrint(fmt::format("      E_var[root {:3}] = {:16.12}\n", iroot, ci_energies[iroot]));

        arma::dvec energy_diff(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            energy_diff[iroot] = energies_per_iter[iter_sci][iroot] - energies_per_iter[iter_sci - 1][iroot];

        bool converged = true;
        for (size_t iroot = 0; iroot < n_roots; iroot++)
            if (abs(energy_diff[iroot]) > Settings::getEnergyTol())
                converged = false;

        if (converged)
        {
            palPrint(fmt::format("\n   SCI converged - E_diff = {:.2e} less than E_tol = {:.2e}\n",
                                 arma::min(arma::abs(energy_diff)), Settings::getEnergyTol()));
            break;
        }
    }
    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration{end - start};

    palPrint(fmt::format("\n      Final CFG-dimension: {:9d}\n", wave_function->getNumCFGs()));
    palPrint(fmt::format("      Final CSF-dimension: {:9d}\n", wave_function->getNumCSFs()));

    palPrint(fmt::format("\n\n    Duration of variational part: {:6.4} s.\n", duration.count()));
}

void SCI::solveCI(vector<double> &ci_energies,
                  vector<vector<double>> &ci_vectors,
                  vector<vector<double>> &energies_per_iter,
                  vector<unordered_map<string, double>> &previous_ci_coeffs_map)
{
    palPrint("\n\n   Constructing configuration connections...              ");

    auto start{std::chrono::steady_clock::now()};
    connection_map_2el connections_2el_new;
    connection_map_1el connections_1el_new;
    connection_map_dia connections_dia_new;
    connections->constructConnections(cfgs_new, wave_function, wave_function,
                                      connections_1el_new, connections_2el_new,
                                      connections_dia_new);

    auto end(std::chrono::steady_clock::now());
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format("done ({:.4} s)", duration.count()));

    palPrint("\n   Constructing coupling coefficients...");

    start = std::chrono::steady_clock::now();
    coupling_coeffs->constructCCs(connections_1el_new, connections_2el_new,
                                  connections_dia_new, wave_function, ccs_2el,
                                  ccs_1el, ccs_dia);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration<double>(end - start);
    palPrint(fmt::format("done ({:.4} s)", duration.count()));

    mergeConnections(connections_1el_new, connections_1el);
    mergeConnections(connections_2el_new, connections_2el);
    mergeConnections(connections_dia_new, connections_dia);

    // TODO: printMemoryEstimate();

    /*
     * Diagonalizing the Hamiltonian
     */
    auto [eigenvalues, eigenvectors] = lible::davidson::diagonalize(
        n_roots,
        [this]()
        { return SCI::calcDiag(); },
        [this](const vector<double> &diag)
        { return SCI::calcGuess(diag); },
        [this](const vector<double> &trial)
        { return SCI::calcSigma(trial); });

    ci_energies = std::move(eigenvalues);
    ci_vectors = std::move(eigenvectors);

    energies_per_iter.push_back(ci_energies);

    previous_ci_coeffs_map = mapPreviousCIVector(ci_vectors);
}
