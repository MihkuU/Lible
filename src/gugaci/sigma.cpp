#include <lible/brain.hpp>
#include <lible/gci_impl.hpp>
#include <lible/util.h>

#ifdef _USE_MPI_
#include <mpl/mpl.hpp>
#endif

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::string;
using std::vector;

vector<double> GCI::Impl::calcSigma(const vector<double> &trial)
{
    arma::dvec sigma(trial.size(), arma::fill::zeros);

#pragma omp parallel
    {
        arma::dvec sigma_omp(wave_function->getNumCSFs(), arma::fill::zeros);
        arma::dvec trial_omp = arma::conv_to<arma::dvec>::from(trial);

        int thread_num = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        int rank_total, size_total;
#ifdef _USE_MPI_
        // int rank_total = returnTotalRank(world);
        // int size_total = returnTotalSize(world);
        rank_total = omp_get_thread_num(); // TMP
        size_total = omp_get_num_threads(); // TMP
#else
        rank_total = omp_get_thread_num();
        size_total = omp_get_num_threads();
#endif

        /* No excitation (diagonal) */
        int ipal = 0;
        for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            CFG *cfg = wave_function->getCFGPtr(icfg);
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

            size_t dim = wave_function->getDim(icfg);
            size_t pos = wave_function->getPos(icfg);

            for (size_t mu = 0; mu < dim; mu++)
                sigma_omp(pos + mu) += val * trial_omp(pos + mu);
        }

        /* One-electron excitations */
        ipal = 0;
        for (const auto &[key, connections] : connections_1el)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_1el.at(key);

            // int nue_right = get<2>(key);
            // int max_sf_dim = weylFormula(spin, nue_right);
            // flat_cc_map ccs_flat = ccs_1el_flat.at(key);

            for (const auto &[icfg_left, icfg_right, pq, phase] : connections)
            {
                auto [p, q] = pq1DTo2D(pq, n_orbs);

                CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
                CFG *cfg_right = wave_function->getCFGPtr(icfg_right);

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

                size_t pos_left = wave_function->getPos(icfg_left);
                size_t pos_right = wave_function->getPos(icfg_right);

                vector<int> sf_idxs_left = cfg_left->getSFIdxs();
                vector<int> sf_idxs_right = cfg_right->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        // double cc = ccs_flat.at(sf_idxs_left[mu] * max_sf_dim + sf_idxs_right[nu]);
                        sigma_omp(pos_left + mu) += contrib * cc * trial_omp(pos_right + nu);
                        sigma_omp(pos_right + nu) += contrib * cc * trial_omp(pos_left + mu);
                    }
            }
        }

        /* Excitations to RI-space and back (conf-diagonal) */
        ipal = 0;
        for (const auto &[key, connections] : connections_dia)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_dia.at(key);

            // int nue_right = get<2>(key);
            // int max_sf_dim = weylFormula(spin, nue_right);
            // flat_cc_map ccs_flat = ccs_dia_flat.at(key);

            for (const auto &[icfg, pqqp] : connections)
            {
                auto [p, q, r, s] = pqrs1DTo4D(pqqp, n_orbs);
                double two_el_int = 0.5 * two_el_ints(p, q, q, p);

                size_t pos = wave_function->getPos(icfg);

                CFG *cfg = wave_function->getCFGPtr(icfg);

                vector<int> sf_idxs = cfg->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs[mu]).at(sf_idxs[nu]);
                        // double cc = ccs_flat.at(sf_idxs[mu] * max_sf_dim + sf_idxs[nu]);
                        sigma_omp(pos + mu) += two_el_int * cc * trial_omp(pos + nu);
                    }
            }
        }

        /* Two-electron excitations */
        ipal = 0;
        for (const auto &[key, connections] : connections_2el)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            cc_map ccs = ccs_2el.at(key);

            // int nue_right = get<4>(key);
            // int max_sf_dim = weylFormula(spin, nue_right);
            // flat_cc_map ccs_flat = ccs_2el_flat.at(key);

            for (const auto &[icfg_left, icfg_right, pqrs, phase, two_el_ex] : connections)
            {
                size_t pos_left = wave_function->getPos(icfg_left);
                size_t pos_right = wave_function->getPos(icfg_right);

                double fac = 0.5;
                if (two_el_ex)
                    fac = 1.0;
                if (phase)
                    fac *= -1;


                auto [p, q, r, s] = pqrs1DTo4D(pqrs, n_orbs);
                double contrib = two_el_ints(p, q, r, s);
                contrib *= fac;

                CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
                CFG *cfg_right = wave_function->getCFGPtr(icfg_right);

                vector<int> sf_idxs_left = cfg_left->getSFIdxs();
                vector<int> sf_idxs_right = cfg_right->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        // double cc = ccs_flat.at(sf_idxs_left[mu] * max_sf_dim + sf_idxs_right[nu]);
                        sigma_omp(pos_left + mu) += contrib * cc * trial_omp(pos_right + nu);
                        sigma_omp(pos_right + nu) += contrib * cc * trial_omp(pos_left + mu);
                    }
            }
        }

#pragma omp critical
        {
            sigma += sigma_omp;
        }
    }

#ifdef _USE_MPI_
    // mpl::contiguous_layout<double> layout(sigma.n_elem);
    // // TODO: make this use the local communicator.
    // Para::comm_world.allreduce([](auto a, auto b)
    //                            { return a + b; },
    //                            sigma.memptr(), layout);
#endif

    return arma::conv_to<vector<double>>::from(sigma);
}

vector<double> GCI::Impl::calcSigma(const vec2d &aux_1el_ints,
                                    const vec4d &aux_2el_ints,
                                    const vector<double> &trial)
{
    /*
     * Calculates the sigma-vector with arbitrary Hamiltonian matrix elements
     *  - symmetry of the matrix elements is not assumed.
     *
     * Can be used in response calculations (e.g. SA-MCSCF gradient) or in CI-coupled
     * second order orbital optimization.
     */
    arma::dvec sigma(trial.size(), arma::fill::zeros);

#pragma omp parallel
    {
        arma::dvec sigma_omp(wave_function->getNumCSFs(), arma::fill::zeros);
        arma::dvec trial_omp = arma::conv_to<arma::dvec>::from(trial);

        int thread_num = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

#ifdef _USE_MPI_
        // int rank_total = returnTotalRank(world);
        // int size_total = returnTotalSize(world);
        int rank_total;
        int size_total;
#else
        int rank_total = thread_num;
        int size_total = num_threads;
#endif

        /* No excitation (diagonal) */
        int ipal = 0;
        for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            CFG *cfg = wave_function->getCFGPtr(icfg);
            string onv = cfg->getONV();

            double val = 0;
            for (size_t p = 0; p < n_orbs; p++)
            {
                int p_occ = onv[p] - '0';
                if (p_occ == 0)
                    continue;

                val += p_occ * (aux_1el_ints(p, p) + 0.5 * p_occ * aux_2el_ints(p, p, p, p));
                for (size_t q = 0; q < n_orbs; q++)
                {
                    if (p == q)
                        continue;

                    int q_occ = onv[q] - '0';
                    if (q_occ == 0)
                        continue;

                    val += 0.5 * p_occ * q_occ * aux_2el_ints(p, p, q, q);
                }
            }

            size_t dim = wave_function->getDim(icfg);
            size_t pos = wave_function->getPos(icfg);

            for (int mu = 0; mu < dim; mu++)
                sigma_omp(pos + mu) += val * trial_omp(pos + mu);
        }

        /* One-electron excitations */
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
                CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
                CFG *cfg_right = wave_function->getCFGPtr(icfg_right);

                string onv_left = cfg_left->getONV();
                string onv_right = cfg_right->getONV();

                // auto [p, q] = pq1DTo2D(pq);
                // size_t qp = pq2DTo1D(q, p);
                auto [p, q] = pq1DTo2D(pq, n_orbs);
                // size_t qp = pq2DTo1D(q, p, n_orbs);
                double contrib1 = aux_1el_ints(p, q);
                double contrib2 = aux_1el_ints(q, p);
                for (size_t r = 0; r < n_orbs; r++)
                {
                    int r_occ = onv_right[r] - '0';
                    if (r_occ == 0)
                        continue;

                    // size_t pqrr = pqrs4DTo1D(p, q, r, r);
                    // size_t rrqp = pqrs4DTo1D(r, r, q, p);
                    // size_t pqrr = pqrs4DTo1D(p, q, r, r, n_orbs);
                    // size_t rrqp = pqrs4DTo1D(r, r, q, p, n_orbs);
                    contrib1 += 0.5 * r_occ * aux_2el_ints(p, q, r, r);
                    contrib2 += 0.5 * r_occ * aux_2el_ints(r, r, q, p);
                }

                for (size_t r = 0; r < n_orbs; r++)
                {
                    int r_occ = onv_left[r] - '0';
                    if (r_occ == 0)
                        continue;

                    // size_t rrpq = pqrs4DTo1D(r, r, p, q);
                    // size_t qprr = pqrs4DTo1D(q, p, r, r);
                    // size_t rrpq = pqrs4DTo1D(r, r, p, q, n_orbs);
                    // size_t qprr = pqrs4DTo1D(q, p, r, r, n_orbs);
                    contrib1 += 0.5 * r_occ * aux_2el_ints(r, r, p, q);
                    contrib2 += 0.5 * r_occ * aux_2el_ints(q, p, r, r);
                }

                double fac = 1;
                if (phase)
                    fac = -1;
                contrib1 *= fac;
                contrib2 *= fac;

                size_t pos_left = wave_function->getPos(icfg_left);
                size_t pos_right = wave_function->getPos(icfg_right);

                vector<int> sf_idxs_left = cfg_left->getSFIdxs();
                vector<int> sf_idxs_right = cfg_right->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        sigma_omp(pos_left + mu) += contrib1 * cc * trial_omp(pos_right + nu);
                        sigma_omp(pos_right + nu) += contrib2 * cc * trial_omp(pos_left + mu);
                    }
            }
        }

        /* Excitations to RI-space and back (conf-diagonal) */
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
                double two_el_int = 0.5 * aux_2el_ints(p, q, q, p);

                size_t pos = wave_function->getPos(icfg);

                CFG *cfg = wave_function->getCFGPtr(icfg);

                vector<int> sf_idxs = cfg->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs[mu]).at(sf_idxs[nu]);
                        sigma_omp(pos + mu) += two_el_int * cc * trial_omp(pos + nu);
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
                // auto [p, q, r, s] = pqrs1DTo4D(pqrs);
                // size_t srqp = pqrs4DTo1D(s, r, q, p);
                auto [p, q, r, s] = pqrs1DTo4D(pqrs, n_orbs);
                size_t srqp = pqrs4DTo1D(s, r, q, p, n_orbs);

                size_t pos_left = wave_function->getPos(icfg_left);
                size_t pos_right = wave_function->getPos(icfg_right);

                double fac = 0.5;
                if (phase)
                    fac *= -1;

                double contrib1 = aux_2el_ints(p, q, r, s);
                double contrib2 = aux_2el_ints(s, r, q, p);

                if (two_el_ex)
                {
                    /*
                     * Due to using symmetry of Hamiltonian
                     */
                    // size_t rspq = pqrs4DTo1D(r, s, p, q);
                    // size_t qpsr = pqrs4DTo1D(q, p, s, r);
                    // contrib1 += eff_2el_part(rspq);
                    // contrib2 += eff_2el_part(qpsr);
                    // TODO: do the above commented method instead
                    double contrib1_ = contrib1;
                    contrib1 += contrib2;
                    contrib2 += contrib1_;
                }
                contrib1 *= fac;
                contrib2 *= fac;

                CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
                CFG *cfg_right = wave_function->getCFGPtr(icfg_right);

                vector<int> sf_idxs_left = cfg_left->getSFIdxs();
                vector<int> sf_idxs_right = cfg_right->getSFIdxs();
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        sigma_omp(pos_left + mu) += contrib1 * cc * trial_omp(pos_right + nu);
                        sigma_omp(pos_right + nu) += contrib2 * cc * trial_omp(pos_left + mu);
                    }
            }
        }

#pragma omp critical
        {
            sigma += sigma_omp;
        }
    }
    // mpi::all_reduce(mpi::communicator(world, mpi::comm_duplicate),
    //                 sigma, sigma, std::plus<arma::dvec>());

    return arma::conv_to<vector<double>>::from(sigma);
}

vector<double> GCI::Impl::calcSigma(const DataFOIS &data_fois,
                                    const std::vector<double> &ci_vector,
                                    const WaveFunction *wfn_cipsi)
{
    vector<double> sigma(wfn_cipsi->getNumCSFs(), 0);

    for (const auto &[key, connections] : data_fois.connections_1el)
    {
        cc_map ccs = ccs_1el.at(key);
        for (const auto &[icfg_left, icfg_right, pq, phase] : connections)
        {
            auto [p, q] = pq1DTo2D(pq, n_orbs);
            double contrib = one_el_ints(p, q);

            double fac = 1;
            if (phase)
                fac = -1;
            contrib *= fac;

            size_t pos_left = wfn_cipsi->getPos(icfg_left);
            size_t pos_right = wave_function->getPos(icfg_right);

            const CFG *cfg_left = wfn_cipsi->getCFGPtr(icfg_left);
            const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
            vector<int> sf_idxs_left = cfg_left->getSFIdxs();
            vector<int> sf_idxs_right = cfg_right->getSFIdxs();
            for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                {
                    double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                    sigma[pos_left + mu] += contrib * cc * ci_vector[pos_right + nu];
                }
        }
    }

    for (const auto &[key, connections] : data_fois.connections_EpqErr)
    {
        cc_map ccs = ccs_1el.at(key);
        for (const auto &connection : connections)
        {
            size_t icfg_left = get<0>(connection.first);
            size_t icfg_right = get<1>(connection.first);
            bool phase = get<2>(connection.first);

            double fac = 1;
            if (phase)
                fac = -1;

            size_t pos_left = wfn_cipsi->getPos(icfg_left);
            size_t pos_right = wave_function->getPos(icfg_right);

            const CFG *cfg_left = wfn_cipsi->getCFGPtr(icfg_left);
            const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
            vector<int> sf_idxs_left = cfg_left->getSFIdxs();
            vector<int> sf_idxs_right = cfg_right->getSFIdxs();
            string onv_right = cfg_right->getONV();

            for (auto &pqrr_r : connection.second)
            {
                size_t pqrr = get<0>(pqrr_r);
                auto [p, q, r, s] = pqrs1DTo4D(pqrr, n_orbs);
                double two_el_int = two_el_ints(p, q, r, s);
                two_el_int *= fac;
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        sigma[pos_left + mu] += two_el_int * cc * ci_vector[pos_right + nu];
                    }
            }
        }
    }

    for (const auto &[key, connections] : data_fois.connections_ErrEpq)
    {
        cc_map ccs = ccs_1el.at(key);
        for (const auto &connection : connections)
        {
            size_t icfg_left = get<0>(connection.first);
            size_t icfg_right = get<1>(connection.first);
            bool phase = get<2>(connection.first);

            double fac = 1;
            if (phase)
                fac = -1;

            size_t pos_left = wfn_cipsi->getPos(icfg_left);
            size_t pos_right = wave_function->getPos(icfg_right);

            const CFG *cfg_left = wfn_cipsi->getCFGPtr(icfg_left);
            const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
            vector<int> sf_idxs_left = cfg_left->getSFIdxs();
            vector<int> sf_idxs_right = cfg_right->getSFIdxs();
            string onv_left = cfg_left->getONV();

            for (auto &rrpq_r : connection.second)
            {
                size_t rrpq = get<0>(rrpq_r);
                auto [p, q, r, s] = pqrs1DTo4D(rrpq, n_orbs);
                double two_el_int = two_el_ints(p, q, r, s);
                two_el_int *= fac;
                for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                    for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                    {
                        double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                        sigma[pos_left + mu] += two_el_int * cc * ci_vector[pos_right + nu];
                    }
            }
        }
    }

    for (const auto &[key, connections] : data_fois.connections_2el)
    {
        cc_map ccs = ccs_2el.at(key);
        for (const auto &[icfg_left, icfg_right, pqrs, phase, two_el_ex] : connections)
        {
            double fac = 0.5;
            if (two_el_ex)
                fac = 1.0;
            if (phase)
                fac *= -1;

            auto [p, q, r, s] = pqrs1DTo4D(pqrs, n_orbs);
            double contrib = two_el_ints(p, q, r, s);
            contrib *= fac;

            size_t pos_left = wfn_cipsi->getPos(icfg_left);
            size_t pos_right = wave_function->getPos(icfg_right);

            const CFG *cfg_left = wfn_cipsi->getCFGPtr(icfg_left);
            CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
            vector<int> sf_idxs_left = cfg_left->getSFIdxs();
            vector<int> sf_idxs_right = cfg_right->getSFIdxs();

            for (size_t mu = 0; mu < sf_idxs_left.size(); mu++)
                for (size_t nu = 0; nu < sf_idxs_right.size(); nu++)
                {
                    double cc = ccs.at(sf_idxs_left[mu]).at(sf_idxs_right[nu]);
                    sigma[pos_left + mu] += contrib * cc * ci_vector[pos_right + nu];
                }
        }
    }

    return sigma;
}
