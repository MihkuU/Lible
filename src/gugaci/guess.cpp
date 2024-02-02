#include <lible/gci_impl.hpp>
#include <lible/gci_settings.hpp>
#include <lible/connections.hpp>
#include <lible/coupling_coeffs.hpp>
#include <lible/util.h>

#include <armadillo>
#include <omp.h>

#ifdef _USE_MPI_
#include <lible/brain.hpp>
#endif

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::map;
using std::pair;
using std::set;
using std::string;
using std::tuple;
using std::vector;

vector<vector<double>> GCI::Impl::calcGuess(const vector<double> &diag)
{
    size_t dim_wfn = wave_function->getNumCSFs();

    int guess_dim = Settings::getGuessDim();
    if (guess_dim > dim_wfn)
        guess_dim = dim_wfn;

    /*
     * After some iterations, start using the previous CI-vector as the initial
     * guess.
     */

    if (iter_sci >= Settings::getIterContinueEigvec())
    {
        vector<vector<double>> guess(n_roots);
        for (size_t iroot = 0; iroot < n_roots; iroot++)
        {
            vector<double> guess_vector(wave_function->getNumCSFs(), 0);
            size_t idx = 0;
            for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
            {
                CFG *cfg = wave_function->getCFGPtr(icfg);
                for (size_t icsf = 0; icsf < cfg->getNumCSFs(); icsf++)
                {
                    string csf = cfg->getCSF(icsf);
                    if (previous_ci_coeffs_map[iroot].find(csf) == previous_ci_coeffs_map[iroot].end())
                    {
                        idx++;
                        continue;
                    }
                    else
                    {
                        guess_vector[idx] = previous_ci_coeffs_map[iroot][csf];
                        idx++;
                    }
                }
            }
            guess[iroot] = guess_vector;
        }
        return guess;
    }

    /*
     * Gather CSF-s and configurations corresponding to the lowest H-diagonal
     * values.
     */
    using tuple5 = tuple<string, string, size_t, size_t, double>;
    vector<tuple5> sorted_diagonal_elements;
    for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
    {
        CFG *cfg = wave_function->getCFGPtr(icfg);
        string onv = cfg->getONV();

        int dim = wave_function->getDim(icfg);
        int pos = wave_function->getPos(icfg);

        for (int icsf = 0; icsf < dim; icsf++)
        {
            string csf = cfg->getCSF(icsf);
            double diag_val = diag[pos + icsf];
            sorted_diagonal_elements.push_back(std::make_tuple(onv, csf, pos, icsf, diag_val));
        }
    }

    std::sort(sorted_diagonal_elements.begin(), sorted_diagonal_elements.end(),
              [](const tuple5 &lhs, const tuple5 &rhs)
              { return std::abs(std::get<4>(lhs)) > std::abs(std::get<4>(rhs)); });

    map<string, pair<size_t, map<string, size_t>>> sorted_diagonal_elements_map;
    for (size_t i = 0; i < guess_dim; i++)
    {
        tuple5 bundle = sorted_diagonal_elements[i];

        string onv = std::get<0>(bundle);
        string csf = std::get<1>(bundle);
        size_t pos = std::get<2>(bundle);
        size_t icsf = std::get<3>(bundle);

        sorted_diagonal_elements_map[onv].first = pos;
        sorted_diagonal_elements_map[onv].second[csf] = icsf;
    }

    /*
     * Guess-space connections and CC-s
     */
    wfn_ptr wfn_guess = std::make_unique<WaveFunction>(spin);
    set<string> onvs_guess;
    for (auto &item1 : sorted_diagonal_elements_map)
    {
        string onv = item1.first;

        onvs_guess.insert(onv);

        CFG cfg(spin, onv);
        int nue = cfg.getNUE();
        for (auto &item2 : item1.second.second)
        {
            string csf = item2.first;
            string sf = extractSF(csf);
            size_t pos = sfs_map__sf_to_idx.at(nue).at(sf);
            cfg.insertCSF(pos, csf);
        }
        wfn_guess->insertCFG(cfg);
    }

    connection_map_1el connections_1el;
    connection_map_2el connections_2el;
    connection_map_dia connections_dia;

    connections->constructConnections(onvs_guess, wfn_guess, wfn_guess,
                                      connections_1el, connections_2el,
                                      connections_dia);

    coupling_coeffs->constructCCs(connections_1el, connections_2el,
                                  connections_dia, wfn_guess, ccs_2el,
                                  ccs_1el, ccs_dia);                                 

    /*
     * Constructing and diagonalizing the guess-space Hamiltonian
     */
    size_t dim_guess = wfn_guess->getNumCSFs();
    arma::dmat guess_hamiltonian(dim_guess, dim_guess, arma::fill::zeros);

#pragma omp parallel
    {
        int rank_total, size_total;
#ifdef _USE_MPI_
        rank_total = Brain::returnTotalRank();
        size_total = Brain::returnTotalSize();
#else
        rank_total = omp_get_thread_num();
        size_total = omp_get_num_threads();
#endif
        int thread_num = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        arma::dmat guess_hamiltonian_omp(dim_guess, dim_guess, arma::fill::zeros);

        /* No excitation (diagonal) */
        int ipal = 0;
        for (size_t icfg = 0; icfg < wfn_guess->getNumCFGs(); icfg++)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            CFG *cfg = wfn_guess->getCFGPtr(icfg);
            string onv = cfg->getONV();

            double val = 0;
            for (size_t p = 0; p < n_orbs; p++)
            {
                int p_occ = onv[p] - '0';
                if (p_occ == 0)
                    continue;

                val += p_occ * (one_el_ints(p, p) + 0.5 * p_occ * two_el_ints(p, p, p, p));

                for (size_t q = p + 1; q < n_orbs; q++)
                {
                    int q_occ = onv[q] - '0';
                    if (q_occ == 0)
                        continue;

                    val += p_occ * q_occ * two_el_ints(p, p, q, q);
                }
            }

            int dim = wfn_guess->getDim(icfg);
            int pos = wfn_guess->getPos(icfg);

            for (int mu = 0; mu < dim; mu++)
                guess_hamiltonian_omp(pos + mu, pos + mu) += val;
        }

        ipal = 0;
        for (const auto &[key, connections] : connections_1el)
        {
            if (ipal % num_threads != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_1el.at(key);

            for (const auto &[icfg_left, icfg_right, pq, phase] : connections)
            {
                auto [p, q] = pq1DTo2D(pq, n_orbs);

                CFG *cfg_left = wfn_guess->getCFGPtr(icfg_left);
                CFG *cfg_right = wfn_guess->getCFGPtr(icfg_right);

                string onv_left = cfg_left->getONV();
                string onv_right = cfg_right->getONV();

                double contrib = one_el_ints(p, q);
                for (size_t r = 0; r < n_orbs; r++)
                {
                    int r_occ = onv_right[r] - '0';
                    if (r_occ == 0)
                        continue;

                    contrib += 0.5 * r_occ * two_el_ints(p, q, r, r);
                }

                for (size_t r = 0; r < n_orbs; r++)
                {
                    int r_occ = onv_left[r] - '0';
                    if (r_occ == 0)
                        continue;

                    contrib += 0.5 * r_occ * two_el_ints(r, r, p, q);
                }

                double fac = 1;
                if (phase)
                    fac = -1;
                contrib *= fac;

                int pos_left = wfn_guess->getPos(icfg_left);
                int pos_right = wfn_guess->getPos(icfg_right);

                vector<int> pos_csf_left = cfg_left->getSFIdxs();
                vector<int> pos_csf_right = cfg_right->getSFIdxs();

                for (size_t mu = 0; mu < pos_csf_left.size(); mu++)
                {
                    for (size_t nu = 0; nu < pos_csf_right.size(); nu++)
                    {
                        // double cc = ccs[pos_csf_left[mu]][pos_csf_right[nu]];
                        // double cc = ccs.at(pos_csf_left[mu]).at(pos_csf_right[nu]);
                        double cc = ccs.at(std::make_pair(pos_csf_left[mu], pos_csf_right[nu]));
                        guess_hamiltonian_omp(pos_left + mu, pos_right + nu) += contrib * cc;
                        guess_hamiltonian_omp(pos_right + nu, pos_left + mu) += contrib * cc;
                    }
                }
            }
        }

        /* Excitations to RI-space and back (cfg-diagonal) */
        ipal = 0;
        for (const auto &[key, connections] : connections_dia)
        {
            if (ipal % num_threads != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_dia.at(key);

            for (const auto &[icfg, pqqp] : connections)
            {
                auto [p, q, r, s] = pqrs1DTo4D(pqqp, n_orbs);
                double two_el_int = 0.5 * two_el_ints(p, q, r, s);

                int pos = wfn_guess->getPos(icfg);

                CFG *cfg = wfn_guess->getCFGPtr(icfg);

                vector<int> pos_csf = cfg->getSFIdxs();
                for (size_t mu = 0; mu < pos_csf.size(); mu++)
                    for (size_t nu = 0; nu < pos_csf.size(); nu++)
                    {
                        // double cc = ccs[pos_csf[mu]][pos_csf[nu]];
                        // double cc = ccs.at(pos_csf[mu]).at(pos_csf[nu]);
                        double cc = ccs.at(std::make_pair(pos_csf[mu], pos_csf[nu]));
                        guess_hamiltonian_omp(pos + mu, pos + nu) += two_el_int * cc;
                    }
            }
        }

        /* Two-electron excitations */
        ipal = 0;
        for (const auto &[key, connections] : connections_2el)
        {
            if (ipal % num_threads != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_2el.at(key);

            for (const auto &[icfg_left, icfg_right, pqrs, phase, two_el_ex] : connections)
            {
                size_t pos_left = wfn_guess->getPos(icfg_left);
                size_t pos_right = wfn_guess->getPos(icfg_right);

                double fac = 0.5;
                if (two_el_ex)
                    fac = 1.0;
                if (phase)
                    fac *= -1;

                auto [p, q, r, s] = pqrs1DTo4D(pqrs, n_orbs);
                double contrib = two_el_ints(p, q, r, s);
                contrib *= fac;

                CFG *cfg_left = wfn_guess->getCFGPtr(icfg_left);
                CFG *cfg_right = wfn_guess->getCFGPtr(icfg_right);
                vector<int> pos_csf_left = cfg_left->getSFIdxs();
                vector<int> pos_csf_right = cfg_right->getSFIdxs();

                for (size_t mu = 0; mu < pos_csf_left.size(); mu++)
                    for (size_t nu = 0; nu < pos_csf_right.size(); nu++)
                    {
                        // double cc = ccs[pos_csf_left[mu]][pos_csf_right[nu]];
                        // double cc = ccs.at(pos_csf_left[mu]).at(pos_csf_right[nu]);
                        double cc = ccs.at(std::make_pair(pos_csf_left[mu], pos_csf_right[nu]));
                        guess_hamiltonian_omp(pos_left + mu, pos_right + nu) += contrib * cc;
                        guess_hamiltonian_omp(pos_right + nu, pos_left + mu) += contrib * cc;
                    }
            }
        }

#pragma omp critical
        {
            guess_hamiltonian += guess_hamiltonian_omp;
        }
    }

#ifdef _USE_MPI_
    mpl::contiguous_layout<double> layout(guess_hamiltonian.n_elem);

    Brain::comm_nodes.allreduce([](auto a, auto b)
                                { return a + b; },
                                guess_hamiltonian.memptr(), layout);
#endif

    arma::dmat eigenvectors_guess(dim_guess, dim_guess, arma::fill::zeros);
    arma::dvec eigenvalues_guess(dim_guess, arma::fill::zeros);
    arma::eig_sym(eigenvalues_guess, eigenvectors_guess, guess_hamiltonian);

    /* Mapping the guess eigenvector to CI-vector */
    vector<vector<double>> guess(n_roots);
    for (size_t iroot = 0; iroot < n_roots; iroot++)
    {
        arma::dvec guess_subspace = eigenvectors_guess.col(iroot);
        arma::dvec guess_ci(dim_wfn, arma::fill::zeros);

        for (size_t icfg_guess = 0; icfg_guess < wfn_guess->getNumCFGs(); icfg_guess++)
        {
            CFG *cfg_guess = wfn_guess->getCFGPtr(icfg_guess);
            string onv_guess = cfg_guess->getONV();

            int dim_guess = wfn_guess->getDim(icfg_guess);
            int pos_guess = wfn_guess->getPos(icfg_guess);

            size_t pos_CI = sorted_diagonal_elements_map[onv_guess].first;
            map<string, size_t> csf_idxs_map = sorted_diagonal_elements_map[onv_guess].second;

            for (int icsf = 0; icsf < dim_guess; icsf++)
            {
                string csf = cfg_guess->getCSF(icsf);
                size_t icsf_CI = csf_idxs_map[csf];
                guess_ci(pos_CI + icsf_CI) += guess_subspace(pos_guess + icsf);
            }
        }
        guess[iroot] = arma::conv_to<vector<double>>::from(guess_ci);
    }

    return guess;
}