#include <lible/prefix_algorithm.h>
#include <lible/gci_util.h>

#ifdef _USE_MPI_
#include <mpi.h>
#endif

namespace LG = lible::guga;

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::pair;
using std::set;
using std::string;
using std::tuple;
using std::vector;

bool compByScnd(const tuple<int, double> &a, const tuple<int, double> &b)
{
    return get<1>(a) > get<1>(b);
}

bool compByThrd(const tuple<int, int, double> &a, const tuple<int, int, double> &b)
{
    return get<2>(a) > get<2>(b);
}

bool compByFrth(const tuple<int, int, int, double> &a, const tuple<int, int, int, double> &b)
{
    return get<3>(a) > get<3>(b);
}

vector<string> GCI::Impl::PrefixAlgorithm::prefixBonanza(const set<string> &cfgs)
{
    set<string> prefixes_wfn;
    set<string> prefixes_excited;
    vector<string> prefixes_vec;

    vector<string> prefixes_scattered;

#ifdef _USE_MPI_
    // int rank;
    // MPI_Comm_rank(impl->world, &rank);

    // if (rank == 0)
    // {
#endif
        for (const string &cfg : cfgs)
        {
            string prefix;
            if (prefix_size < n_orbs)
                prefix = cfg.substr(0, prefix_size);
            else if (prefix_size >= n_orbs)
                prefix = cfg.substr(0, n_orbs);
            prefixes_wfn.insert(prefix);
        }

        prefixes_excited = prefixes_wfn;
        if (prefix_size < n_orbs)
        {
            for (const string &prefix : prefixes_wfn)
            {
                // 1x annihilation
                for (size_t q = 0; q < prefix_size; q++)
                {
                    if (prefix[q] == '0')
                        continue;
                    string new_prefix = prefix;
                    new_prefix[q] -= 1;
                    prefixes_excited.insert(new_prefix);
                }

                // 2x annihilation
                for (size_t q = 0; q < prefix_size; q++)
                {
                    if (prefix[q] == '0')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[q] -= 1;
                    for (size_t s = 0; s < prefix_size; s++)
                    {
                        if (new_prefix1[s] == '0')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[s] -= 1;
                        prefixes_excited.insert(new_prefix2);
                    }
                }

                // 1x creation
                for (size_t p = 0; p < prefix_size; p++)
                {
                    if (prefix[p] == '2')
                        continue;
                    string new_prefix = prefix;
                    new_prefix[p] += 1;
                    prefixes_excited.insert(new_prefix);
                }

                // 2x creation
                for (size_t p = 0; p < prefix_size; p++)
                {
                    if (prefix[p] == '2')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[p] += 1;
                    for (size_t r = 0; r < prefix_size; r++)
                    {
                        if (new_prefix1[r] == '2')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[r] += 1;
                        prefixes_excited.insert(new_prefix2);
                    }
                }

                // 2x annihilation and 1x creation
                for (size_t q = 0; q < prefix_size; q++)
                {
                    if (prefix[q] == '0')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[q] -= 1;
                    for (size_t s = 0; s < prefix_size; s++)
                    {
                        if (new_prefix1[s] == '0')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[s] -= 1;
                        for (size_t p = 0; p < prefix_size; p++)
                        {
                            if (new_prefix2[p] == '2')
                                continue;
                            string new_prefix3 = new_prefix2;
                            new_prefix3[p] += 1;
                            prefixes_excited.insert(new_prefix3);
                        }
                    }
                }

                // 1x annihilation and 2x creation
                for (size_t q = 0; q < prefix_size; q++)
                {
                    if (prefix[q] == '2')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[q] += 1;
                    for (size_t s = 0; s < prefix_size; s++)
                    {
                        if (new_prefix1[s] == '2')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[s] += 1;
                        for (size_t p = 0; p < prefix_size; p++)
                        {
                            if (new_prefix2[p] == '0')
                                continue;
                            string new_prefix3 = new_prefix2;
                            new_prefix3[p] -= 1;
                            prefixes_excited.insert(new_prefix3);
                        }
                    }
                }

                // single excitation
                for (size_t q = 0; q < prefix_size; q++)
                {
                    if (prefix[q] == '0')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[q] -= 1;
                    for (size_t p = 0; p < prefix_size; p++)
                    {
                        if (p == q)
                            continue;
                        if (new_prefix1[p] == '2')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[p] += 1;
                        prefixes_excited.insert(new_prefix2);
                    }
                }

                // double excitation
                for (size_t p = 0; p < prefix_size; p++)
                {
                    if (prefix[p] == '2')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[p] += 1;
                    for (size_t q = 0; q < prefix_size; q++)
                    {
                        if (p == q)
                            continue;
                        if (new_prefix1[q] == '0')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[q] -= 1;
                        for (size_t r = 0; r < prefix_size; r++)
                        {
                            if (new_prefix2[r] == '2')
                                continue;
                            string new_prefix3 = new_prefix2;
                            new_prefix3[r] += 1;
                            for (size_t s = 0; s < prefix_size; s++)
                            {
                                if (r == s)
                                    continue;
                                if (new_prefix3[s] == '0')
                                    continue;
                                string new_prefix4 = new_prefix3;
                                new_prefix4[s] -= 1;
                                prefixes_excited.insert(new_prefix4);
                            }
                        }
                    }
                }
            }
        }
        else if (prefix_size >= n_orbs)
        {
            for (const string &prefix : prefixes_wfn)
            {
                // single excitation
                for (size_t p = 0; p < n_orbs; p++)
                {
                    if (prefix[p] == '2')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[p] += 1;
                    for (size_t q = 0; q < n_orbs; q++)
                    {
                        if (p == q)
                            continue;
                        if (new_prefix1[q] == '-')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[q] -= 1;

                        size_t nue = count(new_prefix2.begin(), new_prefix2.end(), '1');
                        if (nue < min_nue)
                            continue;

                        prefixes_excited.insert(new_prefix2);
                    }
                }

                // double excitation
                for (size_t p = 0; p < n_orbs; p++)
                {
                    if (prefix[p] == '2')
                        continue;
                    string new_prefix1 = prefix;
                    new_prefix1[p] += 1;
                    for (size_t q = 0; q < n_orbs; q++)
                    {
                        if (p == q)
                            continue;
                        if (new_prefix1[q] == '0')
                            continue;
                        string new_prefix2 = new_prefix1;
                        new_prefix2[q] -= 1;
                        for (size_t r = 0; r < n_orbs; r++)
                        {
                            if (new_prefix2[r] == '2')
                                continue;
                            string new_prefix3 = new_prefix2;
                            new_prefix3[r] += 1;
                            for (size_t s = 0; s < n_orbs; s++)
                            {
                                if (new_prefix3[s] == '0')
                                    continue;
                                string new_prefix4 = new_prefix3;
                                new_prefix4[s] -= 1;

                                size_t nue = count(new_prefix4.begin(), new_prefix4.end(), '1');
                                if (nue < min_nue)
                                    continue;

                                prefixes_excited.insert(new_prefix4);
                            }
                        }
                    }
                }
            }
        }
        for (const string &prefix : prefixes_excited)
            prefixes_vec.push_back(prefix);
        std::random_shuffle(prefixes_vec.begin(), prefixes_vec.end());
#ifdef _USE_MPI_
    // }
#endif

#ifdef _USE_MPI_
    // int size;
    // MPI_Comm_size(impl->world, &size);
    // vector<vector<string>> prefixes_scatter(size);
    // for (size_t ipref = 0; ipref < prefixes_vec.size(); ipref++)
    // {
    //     int rank = ipref % size;
    //     prefixes_scatter[rank].push_back(prefixes_vec[ipref]);
    // }
#else
    for (size_t ipref = 0; ipref < prefixes_vec.size(); ipref++)
        prefixes_scattered.push_back(prefixes_vec[ipref]);
#endif

#ifdef _USE_MPI__
    // mpi::scatter(impl->world, prefixes_scatter, prefixes_scattered, 0);
#endif

    return prefixes_scattered;
}

void GCI::Impl::PrefixAlgorithm::generateCFGsAndConnections(const vector<string> &prefixes,
                                                             const vector<vector<pair<string, arma::dvec>>> &generators_by_roots,
                                                             const wfn_ptr &wfn_right,
                                                             DataFOIS &data_fois, DataVar &data_var)
{
    for (const string &prefix : prefixes)
    {
        for (size_t iroot = 0; iroot < generators_by_roots.size(); iroot++)
        {
            vector<pair<string, arma::dvec>> generators_iroot = generators_by_roots[iroot];

            /* One-electron excitations */
            for (pair<string, arma::dvec> &onv_coeffs : generators_iroot)
            {
                string onv_right = onv_coeffs.first;
                int iconf_right = wfn_right->findCFGPos(onv_right);
                double max_ci_coeff = arma::max(arma::abs(onv_coeffs.second));

                if (max_ci_coeff * max_abs_1el_element < LG::Settings::getEpsilonVar())                 
                    break;

                int na = 0, nc = 0, occ_diff_sum = 0;
                int p, q;
                for (size_t i = 0; i < prefix_size; i++)
                {
                    int occ_diff = prefix[i] - onv_right[i];
                    if (occ_diff == -1)
                    {
                        q = i;
                        na++;
                    }
                    if (occ_diff == 1)
                    {
                        p = i;
                        nc++;
                    }
                    occ_diff_sum += abs(occ_diff);
                }

                size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');
                const CFG *cfg_right_p = wfn_right->getCFGPtr(iconf_right);
                vector<string> sfs_right = returnSFs(impl->sfs_map__idx_to_sf.at(nue_right), cfg_right_p->getSFIdxs());
                prefixHelper_Epq(max_ci_coeff, na, nc, occ_diff_sum, p, q, iconf_right, onv_right, sfs_right, wfn_right, data_fois, data_var);
            }

            /* Two-electron excitations */
            for (pair<string, arma::dvec> &onv_coeffs : generators_iroot)
            {
                string onv_right = onv_coeffs.first;
                double max_ci_coeff = arma::max(arma::abs(onv_coeffs.second));
                int iconf_right = wfn_right->findCFGPos(onv_right);

                if (max_ci_coeff * max_abs_2el_element < LG::Settings::getEpsilonVar())
                    break;

                int na = 0, nc = 0, occ_diff_sum = 0;
                vector<int> a_idxs, c_idxs;
                for (size_t i = 0; i < prefix_size; i++)
                {
                    int occ_diff = prefix[i] - onv_right[i];
                    occ_diff_sum += std::abs(occ_diff);
                    if (occ_diff == -1)
                    {
                        a_idxs.push_back(i);
                        na++;
                    }
                    else if (occ_diff == -2)
                    {
                        a_idxs.push_back(i);
                        a_idxs.push_back(i);
                        na += 2;
                    }
                    else if (occ_diff == 1)
                    {
                        c_idxs.push_back(i);
                        nc++;
                    }
                    else if (occ_diff == 2)
                    {
                        c_idxs.push_back(i);
                        c_idxs.push_back(i);
                        nc += 2;
                    }
                }

                size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');
                const CFG *cfg_right = wfn_right->getCFGPtr(iconf_right);
                vector<string> sfs_right = returnSFs(impl->sfs_map__idx_to_sf.at(nue_right), cfg_right->getSFIdxs());
                prefixHelper_EpqErr(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_ErrEpq(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqEqr(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqErp(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqEpq(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqEpr(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqErq(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
                prefixHelper_EpqErs(max_ci_coeff, na, nc, occ_diff_sum, iconf_right, onv_right, a_idxs, c_idxs, sfs_right, wfn_right, data_fois, data_var);
            }
        }
    }
}

void GCI::Impl::PrefixAlgorithm::constructSortedIntegralLists(const vec2d &one_el_ints,
                                                              const vec4d &two_el_ints,
                                                              double &max_abs_1el_element_out,
                                                              double &max_abs_2el_element_out)
{
    double max_abs_1el_element = 0;
    double max_abs_2el_element = 0;
    for (size_t p = 0; p < n_orbs; p++)
        for (size_t q = 0; q < n_orbs; q++)
        {
            double abs_1el_element = std::abs(one_el_ints(p, q));
            if (abs_1el_element > max_abs_1el_element)
                max_abs_1el_element = abs_1el_element;
            for (size_t r = 0; r < n_orbs; r++)
                for (size_t s = 0; s < n_orbs; s++)
                {
                    double abs_2el_element = std::abs(two_el_ints(p, q, r, s));
                    if (abs_2el_element > max_abs_2el_element)
                        max_abs_2el_element = abs_2el_element;
                }
        }

    max_abs_1el_element_out = max_abs_1el_element;
    max_abs_2el_element_out = max_abs_2el_element;

    integralListsHelper_Epq(one_el_ints,
                            intlist_Epq_pqOut,
                            intlist_Epq_pIn_qOut,
                            intlist_Epq_qIn_pOut);

    integralListsHelper_EpqErr(two_el_ints,
                               intlist_EpqErr_pqrOut,
                               intlist_EpqErr_rIn_pqOut,
                               intlist_EpqErr_qIn_pOut,
                               intlist_EpqErr_pIn_qOut,
                               intlist_EpqErr_pqIn);

    integralListsHelper_ErrEpq(two_el_ints,
                               intlist_ErrEpq_pqrOut,
                               intlist_ErrEpq_rIn_pqOut,
                               intlist_ErrEpq_qIn_pOut,
                               intlist_ErrEpq_pIn_qOut,
                               intlist_ErrEpq_pqIn);

    integralListsHelper_EpqEqr(two_el_ints,
                               intlist_EpqEqr_pqrOut,
                               intlist_EpqEqr_rIn_pOut,
                               intlist_EpqEqr_pIn_rOut,
                               intlist_EpqEqr_prIn);

    integralListsHelper_EpqErp(two_el_ints,
                               intlist_EpqErp_pqrOut,
                               intlist_EpqErp_qIn_rOut,
                               intlist_EpqErp_rIn_qOut,
                               intlist_EpqErp_qrIn);

    integralListsHelper_EpqEpq(two_el_ints,
                               intlist_EpqEpq_pqOut,
                               intlist_EpqEpq_qIn_pOut,
                               intlist_EpqEpq_pIn_qOut);

    integralListsHelper_EpqEpr(two_el_ints,
                               intlist_EpqEpr_pqrOut,
                               intlist_EpqEpr_qIn_prOut,
                               intlist_EpqEpr_pIn_qrOut,
                               intlist_EpqEpr_qrIn_pOut,
                               intlist_EpqEpr_pqIn_rOut);

    integralListsHelper_EpqErq(two_el_ints,
                               intlist_EpqErq_pqrOut,
                               intlist_EpqErq_pIn_qrOut,
                               intlist_EpqErq_qIn_prOut,
                               intlist_EpqErq_prIn_qOut,
                               intlist_EpqErq_pqIn_rOut);

    integralListsHelper_EpqErs(two_el_ints,
                               intlist_EpqErs_sIn_pqrOut,
                               intlist_EpqErs_qIn_prsOut,
                               intlist_EpqErs_pIn_qrsOut,
                               intlist_EpqErs_pqrsOut,
                               intlist_EpqErs_qsIn_prOut,
                               intlist_EpqErs_prIn_qsOut,
                               intlist_EpqErs_pqIn_rsOut,
                               intlist_EpqErs_psIn_qrOut,
                               intlist_EpqErs_pqsIn_rOut,
                               intlist_EpqErs_prsIn_qOut,
                               intlist_EpqErs_pqrIn_sOut);
}

vector<string> GCI::Impl::PrefixAlgorithm::findConnectedSFs1El(const vector<string> &sfs_right,
                                                               const quintet &info_cc)
{
    auto [exctype, nue_left, nue_right, prel, qrel] = info_cc;

    vector<string> connected_sfs;
    switch (exctype)
    {
        int p, q, norb;
    case (ExcType::DS):
    {
        if (prel >= qrel)
        {
            p = prel + 1;
            q = qrel;
        }
        else
        {
            p = prel;
            q = qrel;
        }
        norb = nue_right + 1;

        string right(norb, '1');
        right[q] = '2';

        CFGProto cfg_right(spin, right, sfs_right);
        connected_sfs = findConnectedSFs_DOMOSOMO(p, q, cfg_right);
        break;
    }
    case (ExcType::DV):
    {
        p = prel;
        q = qrel;

        norb = nue_left;

        string right(norb, '1');
        right[p] = '0';
        right[q] = '2';

        CFGProto cfg_right(spin, right, sfs_right);
        connected_sfs = findConnectedSFs_DOMOVirtual(q, p, cfg_right);
        break;
    }
    case (ExcType::SS):
    {
        p = prel;
        q = qrel;

        int norb = nue_right;

        string right(norb, '1');

        CFGProto cfg_right(spin, right, sfs_right);
        connected_sfs = findConnectedSFs_SOMOSOMO(p, q, cfg_right);        
        break;
    }
    case (ExcType::SV):
    {
        if (prel > qrel)
        {
            p = prel;
            q = qrel;
        }
        else
        {
            p = prel;
            q = qrel + 1;
        }

        norb = nue_right + 1;

        string right(norb, '1');
        right[p] = '0';

        CFGProto cfg_right(spin, right, sfs_right);
        connected_sfs = findConnectedSFs_SOMOVirtual(p, q, cfg_right);
        break;
    }
    }

    return connected_sfs;
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_Epq(const vec2d &one_el_ints,
                                                         ints_11_map &intlist_Epq_pqOut,
                                                         ints_11_map &intlist_Epq_pIn_qOut,
                                                         ints_11_map &intlist_Epq_qIn_pOut)
{
    /* dn = 0 */
    // Epq - pq outside prefix
    for (size_t q = prefix_size; q < n_orbs; q++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            if (p == q)
                continue;

            double abs_eff_h = std::abs(one_el_ints(p, q));
            doublets.push_back(std::make_tuple(p, abs_eff_h));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_Epq_pqOut[q] = doublets;
    }

    /* dn = 1 */
    // Epq - p in prefix, q outside prefix
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            double abs_eff_h = std::abs(one_el_ints(p, q));
            doublets.push_back(std::make_tuple(q, abs_eff_h));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_Epq_pIn_qOut[p] = doublets;
    }

    // Epq - q in prefix, p outside prefix
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            double abs_eff_h = std::abs(one_el_ints(p, q));
            doublets.push_back(std::make_tuple(p, abs_eff_h));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_Epq_qIn_pOut[q] = doublets;
    }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqErr(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_EpqErr_pqrOut,
                                                            ints_12_map &intlist_EpqErr_rIn_pqOut,
                                                            ints_12_map &intlist_EpqErr_qIn_pOut,
                                                            ints_12_map &intlist_EpqErr_pIn_qOut,
                                                            ints_21_map &intlist_EpqErr_pqIn)
{
    /* dn = 0 */
    // EpqErr - pqr outside prefix
    for (size_t r = prefix_size; r < n_orbs; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, r));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErr_pqrOut[r] = triplets;
    }

    // EpqErr - r in prefix, pq outside prefix
    for (size_t r = 0; r < prefix_size; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, r));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErr_rIn_pqOut[r] = triplets;
    }

    /* dn = 1 */
    // EpqErr - q in prefix, p outside prefix, r anywhere
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, r));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErr_qIn_pOut[q] = triplets;
    }

    // EpqErr - p in prefix, q outside prefix, r anywhere
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, r));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErr_pIn_qOut[p] = triplets;
    }

    /* dn = 2 */
    // EpqErr - pq in prefix, r anywhere
    for (size_t p = 0; p < prefix_size; p++)
    {
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (p == q)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, r));
                doublets.push_back(std::make_tuple(r, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqErr_pqIn[std::make_pair(p, q)] = doublets;
        }
    }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_ErrEpq(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_ErrEpq_pqrOut,
                                                            ints_12_map &intlist_ErrEpq_rIn_pqOut,
                                                            ints_12_map &intlist_ErrEpq_qIn_pOut,
                                                            ints_12_map &intlist_ErrEpq_pIn_qOut,
                                                            ints_21_map &intlist_ErrEpq_pqIn)
{
    /* dn = 0 */
    // ErrEpq - pqr outside prefix
    for (size_t r = prefix_size; r < n_orbs; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_ErrEpq_pqrOut[r] = triplets;
    }

    for (size_t r = prefix_size; r < n_orbs; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErr_pqrOut[r] = triplets;
    }

    // ErrEpq - r in prefix, pq outside prefix
    for (size_t r = 0; r < prefix_size; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_ErrEpq_rIn_pqOut[r] = triplets;
    }

    /* dn = 1 */
    // ErrEpq - q in prefix, p outside prefix, r anywhere
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            if (p == q)
                continue;

            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_ErrEpq_qIn_pOut[q] = triplets;
    }

    // ErrEpq - p in prefix, q outside prefix, r anywhere
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_ErrEpq_pIn_qOut[p] = triplets;
    }

    /* dn = 2 */
    // ErrEpq - pq in prefix, r anywhere
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (p == q)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t r = 0; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(r, r, p, q));
                doublets.push_back(std::make_tuple(r, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_ErrEpq_pqIn[std::make_pair(p, q)] = doublets;
        }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqEqr(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_EpqEqr_pqrOut,
                                                            ints_12_map &intlist_EpqEqr_rIn_pOut,
                                                            ints_12_map &intlist_EpqEqr_pIn_rOut,
                                                            ints_21_map &intlist_EpqEqr_prIn)
{
    /* dn = 0 */
    // EpqEqr - pqr out
    for (size_t p = prefix_size; p < n_orbs; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            if (p == q)
                continue;

            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                if (r == p or r == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, q, r));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEqr_pqrOut[p] = triplets;
    }

    /* dn = 1 */
    // EpqEqr - r in, p out, q anywhere
    for (size_t r = 0; r < prefix_size; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t q = 0; q < n_orbs; q++)
            {
                if (q == r or q == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, q, r));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEqr_rIn_pOut[r] = triplets;
    }

    // EpqEqr - p in, r out, q anywhere
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t r = prefix_size; r < n_orbs; r++)
        {
            for (size_t q = 0; q < n_orbs; q++)
            {
                if (q == r or q == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, q, r));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEqr_pIn_rOut[p] = triplets;
    }

    /* dn = 2 */
    // EpqEqr - pr in, q anywhere
    for (size_t p = 0; p < prefix_size; p++)
    {
        for (size_t r = 0; r < prefix_size; r++)
        {
            if (p == r)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t q = 0; q < n_orbs; q++)
            {
                if (q == p or q == r)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, q, r));
                doublets.push_back(std::make_tuple(q, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqEqr_prIn[std::make_pair(p, r)] = doublets;
        }
    }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqErp(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_EpqErp_pqrOut,
                                                            ints_12_map &intlist_EpqErp_qIn_rOut,
                                                            ints_12_map &intlist_EpqErp_rIn_qOut,
                                                            ints_21_map &intlist_EpqErp_qrIn)
{
    /* dn = 0 */
    // EpqErp - pqr out
    for (size_t p = prefix_size; p < n_orbs; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            if (p == q)
                continue;

            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                if (r == p or r == q)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, p));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErp_pqrOut[p] = triplets;
    }

    /* dn = 1 */
    // EpqErp - q in, r out, p anywhere
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t r = prefix_size; r < n_orbs; r++)
            for (size_t p = 0; p < n_orbs; p++)
            {
                if (p == q or r == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, p));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErp_qIn_rOut[q] = triplets;
    }

    // EpqErp - r in, q out, p anywhere
    for (size_t r = 0; r < prefix_size; r++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
            for (size_t p = 0; p < n_orbs; p++)
            {
                if (p == q or r == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, p));
                triplets.push_back(std::make_tuple(p, q, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErp_rIn_qOut[r] = triplets;
    }

    /* dn = 2 */
    // EpqErp - qr in, p anywhere
    for (size_t q = 0; q < prefix_size; q++)
        for (size_t r = 0; r < prefix_size; r++)
        {
            if (q == r)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t p = 0; p < n_orbs; p++)
            {
                if (p == q or r == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, p));
                doublets.push_back(std::make_tuple(p, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqErp_qrIn[std::make_pair(q, r)] = doublets;
        }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqEpq(const vec4d &two_el_ints,
                                                            ints_11_map &intlist_EpqEpq_pqOut,
                                                            ints_11_map &intlist_EpqEpq_qIn_pOut,
                                                            ints_11_map &intlist_EpqEpq_pIn_qOut)
{
    /* dn = 0 */
    // EpqEpq - pq out
    for (size_t q = prefix_size; q < n_orbs; q++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            if (p == q)
                continue;

            double abs_act_eri = std::abs(two_el_ints(p, q, p, q));
            doublets.push_back(std::make_tuple(p, abs_act_eri));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_EpqEpq_pqOut[q] = doublets;
    }

    /* dn = 2 */
    // EpqEpq - q in, p out
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            double abs_act_eri = std::abs(two_el_ints(p, q, p, q));
            doublets.push_back(std::make_tuple(p, abs_act_eri));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_EpqEpq_qIn_pOut[q] = doublets;
    }

    // EpqEpq - p in, q out
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, double>> doublets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            double abs_act_eri = std::abs(two_el_ints(p, q, p, q));
            doublets.push_back(std::make_tuple(q, abs_act_eri));
        }
        std::sort(doublets.begin(), doublets.end(), compByScnd);
        intlist_EpqEpq_pIn_qOut[p] = doublets;
    }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqEpr(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_EpqEpr_pqrOut,
                                                            ints_12_map &intlist_EpqEpr_qIn_prOut,
                                                            ints_12_map &intlist_EpqEpr_pIn_qrOut,
                                                            ints_21_map &intlist_EpqEpr_qrIn_pOut,
                                                            ints_21_map &intlist_EpqEpr_pqIn_rOut)
{
    // Due to symmetry of the excitation, enforcing
    // (pq) < (pr) => q < r

    /* dn = 0 */
    // EpqEpr - pqr out
    for (size_t p = prefix_size; p < n_orbs; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            if (q == p)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r == p)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, p, r));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEpr_pqrOut[p] = triplets;
    }

    /* dn = 1 */
    // EpqEpr - q in, pr out
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                if (p == r)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, p, r));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEpr_qIn_prOut[q] = triplets;
    }

    /* dn = 2 */
    // EpqEpr - qr in, p out
    for (size_t q = 0; q < prefix_size; q++)
    {
        for (size_t r = q + 1; r < prefix_size; r++)
        {
            vector<tuple<int, double>> doublets;
            for (size_t p = prefix_size; p < n_orbs; p++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, p, r));
                doublets.push_back(std::make_tuple(p, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqEpr_qrIn_pOut[std::make_pair(q, r)] = doublets;
        }
    }

    // EpqEpr - p in, qr out
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            for (size_t r = q + 1; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, p, r));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqEpr_pIn_qrOut[p] = triplets;
    }

    /* dn = 3 */
    // EpqEpr - pq in, r out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (p == q)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, p, r));
                doublets.push_back(std::make_tuple(r, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqEpr_pqIn_rOut[std::make_pair(p, q)] = doublets;
        }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqErq(const vec4d &two_el_ints,
                                                            ints_12_map &intlist_EpqErq_pqrOut,
                                                            ints_12_map &intlist_EpqErq_pIn_qrOut,
                                                            ints_12_map &intlist_EpqErq_qIn_prOut,
                                                            ints_21_map &intlist_EpqErq_prIn_qOut,
                                                            ints_21_map &intlist_EpqErq_pqIn_rOut)
{
    // Due to symmetry of the excitation, enforcing
    // (pq) < (rq) => p < r

    /* dn = 0 */
    // EpqErq - pqr out
    for (size_t q = prefix_size; q < n_orbs; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            if (p == q)
                continue;
            for (size_t r = p + 1; r < n_orbs; r++)
            {
                if (q == r)
                    continue;
                double abs_act_eri = std::abs(two_el_ints(p, q, r, q));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErq_pqrOut[q] = triplets;
    }

    /* dn = 1 */
    // EpqErq - p in, qr out
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                if (q == r)
                    continue;

                double abs_act_eri = std::abs(two_el_ints(p, q, r, q));
                triplets.push_back(std::make_tuple(q, r, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErq_pIn_qrOut[p] = triplets;
    }

    /* dn = 2 */
    // EpqErq - q in, pr out
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, double>> triplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t r = p + 1; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, q));
                triplets.push_back(std::make_tuple(p, r, abs_act_eri));
            }
        std::sort(triplets.begin(), triplets.end(), compByThrd);
        intlist_EpqErq_qIn_prOut[q] = triplets;
    }

    // EpqErq - pr in, q out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t r = p + 1; r < prefix_size; r++)
        {
            vector<tuple<int, double>> doublets;
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, q));
                doublets.push_back(std::make_tuple(q, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqErq_prIn_qOut[std::make_pair(p, r)] = doublets;
        }

    /* dn = 3 */
    // EpqErq - pq in, r out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (p == q)
                continue;

            vector<tuple<int, double>> doublets;
            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                double abs_act_eri = std::abs(two_el_ints(p, q, r, q));
                doublets.push_back(std::make_tuple(r, abs_act_eri));
            }
            std::sort(doublets.begin(), doublets.end(), compByScnd);
            intlist_EpqErq_pqIn_rOut[std::make_pair(p, q)] = doublets;
        }
}

void GCI::Impl::PrefixAlgorithm::integralListsHelper_EpqErs(const vec4d &two_el_ints,
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
                                                            ints_31_map &intlist_EpqErs_pqrIn_sOut)
{
    // Due to symmetry of the excitation, enforcing
    // (pq) < (rs) => p < r; q and s undetermined

    /* dn = 0 */
    // EpqErs - pqrs out
    for (size_t p = prefix_size; p < n_orbs; p++)
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            if (p == q)
                continue;

            vector<tuple<int, int, double>> triplets;
            for (size_t r = p + 1; r < n_orbs; r++)
            {
                if (q == r)
                    continue;

                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    if (s == p or s == r or s == q)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    triplets.push_back(std::make_tuple(r, s, abs_act_eri));
                }
            }
            std::sort(triplets.begin(), triplets.end(), compByThrd);
            intlist_EpqErs_pqrsOut[std::make_pair(p, q)] = triplets;
        }

    /* dn = 1 */
    // EpqErs - s in, pqr out
    for (size_t s = 0; s < prefix_size; s++)
    {
        vector<tuple<int, int, int, double>> quadruplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (q == p)
                    continue;

                for (size_t r = p + 1; r < n_orbs; r++)
                {
                    if (r == q)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    quadruplets.push_back(std::make_tuple(p, q, r, abs_act_eri));
                }
            }
        std::sort(quadruplets.begin(), quadruplets.end(), compByFrth);
        intlist_EpqErs_sIn_pqrOut[s] = quadruplets;
    }

    // EpqErs - q in, prs out
    for (size_t q = 0; q < prefix_size; q++)
    {
        vector<tuple<int, int, int, double>> quadruplets;
        for (size_t p = prefix_size; p < n_orbs; p++)
            for (size_t r = p + 1; r < n_orbs; r++)
                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    if (s == p or s == r)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    quadruplets.push_back(std::make_tuple(p, r, s, abs_act_eri));
                }
        std::sort(quadruplets.begin(), quadruplets.end(), compByFrth);
        intlist_EpqErs_qIn_prsOut[q] = quadruplets;
    }

    // EpqErs - p in, qrs out
    for (size_t p = 0; p < prefix_size; p++)
    {
        vector<tuple<int, int, int, double>> quadruplets;
        for (size_t q = prefix_size; q < n_orbs; q++)
            for (size_t r = prefix_size; r < n_orbs; r++)
            {
                if (r == q)
                    continue;

                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    if (s == q or s == r)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    quadruplets.push_back(std::make_tuple(q, r, s, abs_act_eri));
                }
            }
        std::sort(quadruplets.begin(), quadruplets.end(), compByFrth);
        intlist_EpqErs_pIn_qrsOut[p] = quadruplets;
    }

    /* dn = 2 */
    // EpqErs - qs in, pr out
    for (size_t q = 0; q < prefix_size; q++)
        for (size_t s = 0; s < prefix_size; s++)
        {
            if (q == s)
                continue;

            vector<tuple<int, int, double>> triplets;
            for (size_t p = prefix_size; p < n_orbs; p++)
            {
                for (size_t r = p + 1; r < n_orbs; r++)
                {
                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    triplets.push_back(std::make_tuple(p, r, abs_act_eri));
                }
            }
            std::sort(triplets.begin(), triplets.end(), compByThrd);
            intlist_EpqErs_qsIn_prOut[std::make_pair(q, s)] = triplets;
        }

    // EpqErs - pr in, qs out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t r = p + 1; r < prefix_size; r++)
        {
            vector<tuple<int, int, double>> triplets;
            for (size_t q = prefix_size; q < n_orbs; q++)
                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    if (q == s)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    triplets.push_back(std::make_tuple(q, s, abs_act_eri));
                }
            std::sort(triplets.begin(), triplets.end(), compByThrd);
            intlist_EpqErs_prIn_qsOut[std::make_pair(p, r)] = triplets;
        }

    // EpqErs - pq in, rs out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (q == p)
                continue;

            vector<tuple<int, int, double>> triplets;
            for (size_t r = prefix_size; r < n_orbs; r++)
                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    if (s == r)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    triplets.push_back(std::make_tuple(r, s, abs_act_eri));
                }
            std::sort(triplets.begin(), triplets.end(), compByThrd);
            intlist_EpqErs_pqIn_rsOut[std::make_pair(p, q)] = triplets;
        }

    // EpqErs - ps in, qr out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t s = 0; s < prefix_size; s++)
        {
            if (s == p)
                continue;

            vector<tuple<int, int, double>> triplets;
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                for (size_t r = prefix_size; r < n_orbs; r++)
                {
                    if (r == q)
                        continue;

                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    triplets.push_back(std::make_tuple(q, r, abs_act_eri));
                }
            }
            std::sort(triplets.begin(), triplets.end(), compByThrd);
            intlist_EpqErs_psIn_qrOut[std::make_pair(p, s)] = triplets;
        }

    /* dn = 3 */
    // EpqErs - pqs in, r out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (q == p)
                continue;

            for (size_t s = 0; s < prefix_size; s++)
            {
                if (s == p or s == q)
                    continue;

                vector<tuple<int, double>> doublets;
                for (size_t r = prefix_size; r < n_orbs; r++)
                {
                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    doublets.push_back(std::make_tuple(r, abs_act_eri));
                }
                std::sort(doublets.begin(), doublets.end(), compByScnd);
                intlist_EpqErs_pqsIn_rOut[std::make_tuple(p, q, s)] = doublets;
            }
        }

    // EpqErs - prs in, q out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t r = p + 1; r < prefix_size; r++)
            for (size_t s = 0; s < prefix_size; s++)
            {
                if (s == r or s == p)
                    continue;

                vector<tuple<int, double>> doublets;
                for (size_t q = prefix_size; q < n_orbs; q++)
                {
                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    doublets.push_back(std::make_tuple(q, abs_act_eri));
                }
                std::sort(doublets.begin(), doublets.end(), compByScnd);
                intlist_EpqErs_prsIn_qOut[std::make_tuple(p, r, s)] = doublets;
            }

    // EpqErs - pqr in, s out
    for (size_t p = 0; p < prefix_size; p++)
        for (size_t q = 0; q < prefix_size; q++)
        {
            if (q == p)
                continue;

            for (size_t r = p + 1; r < prefix_size; r++)
            {
                if (r == q)
                    continue;

                vector<tuple<int, double>> doublets;
                for (size_t s = prefix_size; s < n_orbs; s++)
                {
                    double abs_act_eri = std::abs(two_el_ints(p, q, r, s));
                    doublets.push_back(std::make_tuple(s, abs_act_eri));
                }
                std::sort(doublets.begin(), doublets.end(), compByScnd);
                intlist_EpqErs_pqrIn_sOut[std::make_tuple(p, q, r)] = doublets;
            }
        }
}

void GCI::Impl::PrefixAlgorithm::innerPrefixHelper1El(const int &p, const int &q, const int &icfg_right,
                                                      const size_t &nue_right, const string &onv_right,
                                                      const vector<string> &sfs_right,
                                                      const wfn_ptr &wfn_right,
                                                      DataFOIS &data_fois, DataVar &data_var)
{
    string onv_new = onv_right;
    if (onv_new[q] == '0')
        return;
    onv_new[q] -= 1;
    if (onv_new[p] == '2')
        return;
    onv_new[p] += 1;

    size_t nue_new = std::count(onv_new.begin(), onv_new.end(), '1');
    if (nue_new < min_nue)
        return;

    int icfg_new_right = wfn_right->findCFGPos(onv_new);
    int icfg_new_left = data_fois.wfn->findCFGPos(onv_new);

    if (icfg_new_right != -1)
        return;

    bool phase = determine1ElPhase(p, q, onv_right);
    int pq = pq2DTo1D(p, q, n_orbs);
    quintet info_cc = returnCCInfo(p, q, nue_new, nue_right, onv_new, onv_right);

    vector<string> sfs_left = findConnectedSFs1El(sfs_right, info_cc);
    if (sfs_left.size() == 0)
        return;

    if (icfg_new_left == -1)
    {
        CFG cfg_new(spin, onv_new);
        icfg_new_left = data_fois.wfn->insertCFGandGetPos(cfg_new);
        data_fois.connections_1el[info_cc].push_back(std::make_tuple(icfg_new_left, icfg_right, pq, phase));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_left.begin(), sfs_left.end());
    }
    else
    {
        data_fois.connections_1el[info_cc].push_back(std::make_tuple(icfg_new_left, icfg_right, pq, phase));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_left.begin(), sfs_left.end());
    }
}

void GCI::Impl::PrefixAlgorithm::innerPrefixHelper2El_EpqErr(const int &p, const int &q, const int &r,
                                                             const int &icfg_right, const size_t &nue_right,
                                                             const string &onv_right,
                                                             const vector<string> &sfs_right,
                                                             const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                                             DataVar &data_var)
{
    string onv_new = onv_right;
    if (onv_new[r] == '0')
        return;
    if (onv_new[q] == '0')
        return;
    onv_new[q] -= 1;
    if (onv_new[p] == '2')
        return;
    onv_new[p] += 1;

    size_t nue_new = std::count(onv_new.begin(), onv_new.end(), '1');
    if (nue_new < min_nue)
        return;

    int icfg_new_right = wfn_right->findCFGPos(onv_new);
    int icfg_new_left = data_fois.wfn->findCFGPos(onv_new);

    if (icfg_new_right != -1)
        return;

    bool phase = determine1ElPhase(p, q, onv_right);
    size_t pqrr = pqrs4DTo1D(p, q, r, r, n_orbs);
    quintet info_cc = returnCCInfo(p, q, nue_new, nue_right, onv_new, onv_right);

    vector<string> sfs_new = findConnectedSFs1El(sfs_right, info_cc);
    if (sfs_new.size() == 0)
        return;

    if (icfg_new_left == -1)
    {
        CFG cfg_new(spin, onv_new);
        icfg_new_left = data_fois.wfn->insertCFGandGetPos(cfg_new);
        data_fois.connections_EpqErr[info_cc][std::make_tuple(icfg_new_left, icfg_right, phase)].push_back(std::make_tuple(pqrr, r));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
    else if (icfg_new_left != -1)
    {
        data_fois.connections_EpqErr[info_cc][std::make_tuple(icfg_new_left, icfg_right, phase)].push_back(std::make_tuple(pqrr, r));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
}

void GCI::Impl::PrefixAlgorithm::innerPrefixHelper2El_ErrEpq(const int &p, const int &q, const int &r,
                                                             const int &icfg_right, const size_t &nue_right,
                                                             const string &onv_right,
                                                             const vector<string> &sfs_right,
                                                             const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                                             DataVar &data_var)
{
    string onv_new = onv_right;
    if (onv_new[q] == '0')
        return;
    onv_new[q] -= 1;
    if (onv_new[p] == '2')
        return;
    onv_new[p] += 1;
    if (onv_new[r] == '0')
        return;

    size_t nue_new = count(onv_new.begin(), onv_new.end(), '1');
    if (nue_new < min_nue)
        return;

    int icfg_new_right = wfn_right->findCFGPos(onv_new);
    int icfg_new_left = data_fois.wfn->findCFGPos(onv_new);

    if (icfg_new_right != -1)
        return;

    bool phase = determine1ElPhase(p, q, onv_right);
    size_t rrpq = pqrs4DTo1D(r, r, p, q, n_orbs);
    quintet info_cc = returnCCInfo(p, q, nue_new, nue_right, onv_new, onv_right);

    vector<string> sfs_new = findConnectedSFs1El(sfs_right, info_cc);
    if (sfs_new.size() == 0)
        return;

    if (icfg_new_left == -1)
    {
        CFG cfg_new(spin, onv_new);
        icfg_new_left = data_fois.wfn->insertCFGandGetPos(cfg_new);
        data_fois.connections_ErrEpq[info_cc][std::make_tuple(icfg_new_left, icfg_right, phase)].push_back(std::make_tuple(rrpq, r));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
    else if (icfg_new_left != -1)
    {
        data_fois.connections_ErrEpq[info_cc][std::make_tuple(icfg_new_left, icfg_right, phase)].push_back(std::make_tuple(rrpq, r));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
}

void GCI::Impl::PrefixAlgorithm::innerPrefixHelper2El(const int &p, const int &q, const int &r,
                                                      const int &s, const int &icfg_right,
                                                      const size_t &nue_right, const string &onv_right,
                                                      const vector<string> &sfs_right,
                                                      const wfn_ptr &wfn_right,
                                                      DataFOIS &data_fois, DataVar &data_var)
{
    string onv_ri = onv_right;
    if (onv_ri[s] == '0')
        return;
    onv_ri[s] -= 1;
    if (onv_ri[r] == '2')
        return;
    onv_ri[r] += 1;

    size_t nue_ri = std::count(onv_ri.begin(), onv_ri.end(), '1');
    if (nue_ri < min_nue)
        return;

    string onv_new = onv_ri;
    if (onv_new[q] == '0')
        return;
    onv_new[q] -= 1;
    if (onv_new[p] == '2')
        return;
    onv_new[p] += 1;

    size_t nue_new = std::count(onv_new.begin(), onv_new.end(), '1');
    if (nue_new < min_nue)
        return;

    int icfg_new_right = wfn_right->findCFGPos(onv_new);
    int icfg_new_left = data_fois.wfn->findCFGPos(onv_new);
    if (icfg_new_right != -1)
        return;

    bool phase = determine2ElPhase(p, q, r, s, onv_ri, onv_right);
    size_t pqrs = pqrs4DTo1D(p, q, r, s, n_orbs);

    quintet info_cc1 = returnCCInfo(p, q, nue_new, nue_ri, onv_new, onv_ri);
    quintet info_cc2 = returnCCInfo(r, s, nue_ri, nue_right, onv_ri, onv_right);

    int exctype1 = get<0>(info_cc1);
    int exctype2 = get<0>(info_cc2);
    int prel1 = get<3>(info_cc1);
    int qrel1 = get<4>(info_cc1);
    int prel2 = get<3>(info_cc2);
    int qrel2 = get<4>(info_cc2);

    vector<string> sfs_ri = findConnectedSFs1El(sfs_right, info_cc2);
    vector<string> sfs_new = findConnectedSFs1El(sfs_ri, info_cc1);
    if (sfs_new.size() == 0)
        return;

    if (icfg_new_left == -1)
    {
        CFG cfg_new(spin, onv_new);
        icfg_new_left = data_fois.wfn->insertCFGandGetPos(cfg_new);
        nonet info_cc1_x_cc2 = std::make_tuple(exctype1, exctype2, nue_new, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
        data_fois.connections_2el[info_cc1_x_cc2].push_back(std::make_tuple(icfg_new_left, icfg_right, pqrs, phase, true));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
    else if (icfg_new_left != -1)
    {
        nonet info_cc1_x_cc2 = std::make_tuple(exctype1, exctype2, nue_new, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
        data_fois.connections_2el[info_cc1_x_cc2].push_back(std::make_tuple(icfg_new_left, icfg_right, pqrs, phase, true));
        data_fois.cfgs_sfs[icfg_new_left].insert(sfs_new.begin(), sfs_new.end());
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_Epq(const double &max_ci_coeff, const int &na, const int &nc,
                                                  const int &occ_diff_sum, const int &p, const int &q,
                                                  const int &icfg_right, const string &onv_right,
                                                  const vector<string> &sfs_right,
                                                  const wfn_ptr &wfn_right, DataFOIS &data_fois,
                                                  DataVar &data_var)
{
    if (occ_diff_sum > 2)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            vector<tuple<int, double>> doublet_list = intlist_Epq_pqOut.at(q);
            for (tuple<int, double> &doublet : doublet_list)
            {                
                int p = get<0>(doublet);
                double one_el_int = get<1>(doublet);                
                if (max_ci_coeff * one_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper1El(p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            vector<tuple<int, double>> doublet_list = intlist_Epq_qIn_pOut.at(q);
            for (tuple<int, double> &doublet : doublet_list)
            {                
                int p = get<0>(doublet);
                double one_el_int = get<1>(doublet);                
                if (max_ci_coeff * one_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper1El(p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            vector<tuple<int, double>> doublet_list = intlist_Epq_pIn_qOut.at(p);
            for (tuple<int, double> &doublet : doublet_list)
            {                
                int q = get<0>(doublet);
                double one_el_int = get<1>(doublet);
                if (max_ci_coeff * one_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper1El(p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 1 and nc == 1)
        {            
            double one_el_int = std::abs(impl->one_el_ints(p, q));
            if (max_ci_coeff * one_el_int > Settings::getEpsilonVar())
                innerPrefixHelper1El(p, q, icfg_right, nue_right, onv_right,
                                     sfs_right, wfn_right, data_fois, data_var);
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqErr(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right, const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 2)
        return;

    size_t nue_right = count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t r = prefix_size; r < n_orbs; r++)
        {
            vector<tuple<int, int, double>> triplet_list = intlist_EpqErr_pqrOut.at(r);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);

                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_EpqErr(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        for (size_t r = 0; r < prefix_size; r++)
        {
            vector<tuple<int, int, double>> triplet_list = intlist_EpqErr_rIn_pqOut.at(r);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_EpqErr(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int q = a_idxs[0];
            vector<tuple<int, int, double>> triplet_list = intlist_EpqErr_qIn_pOut.at(q);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_EpqErr(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplet_list = intlist_EpqErr_pIn_qOut.at(p);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_EpqErr(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 1 and nc == 1)
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublet_list = intlist_EpqErr_pqIn.at(std::make_pair(p, q));
            for (tuple<int, double> &doublet : doublet_list)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_EpqErr(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_ErrEpq(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right, const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 2)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t r = prefix_size; r < n_orbs; r++)
        {
            vector<tuple<int, int, double>> triplet_list = intlist_ErrEpq_pqrOut.at(r);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_ErrEpq(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        for (size_t r = 0; r < prefix_size; r++)
        {
            vector<tuple<int, int, double>> triplet_list = intlist_ErrEpq_rIn_pqOut.at(r);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);

                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_ErrEpq(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int q = a_idxs[0];
            vector<tuple<int, int, double>> triplet_list = intlist_ErrEpq_qIn_pOut.at(q);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);

                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_ErrEpq(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplet_list = intlist_ErrEpq_pIn_qOut.at(p);
            for (tuple<int, int, double> &triplet : triplet_list)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);

                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_ErrEpq(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 1 and nc == 1)
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublet_list = intlist_ErrEpq_pqIn.at(std::make_pair(p, q));
            for (tuple<int, double> &doublet : doublet_list)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El_ErrEpq(p, q, r, icfg_right, nue_right, onv_right,
                                                sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqEqr(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right, const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 2)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            vector<tuple<int, int, double>> triplets = intlist_EpqEqr_pqrOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, q, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int r = a_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqEqr_rIn_pOut.at(r);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, q, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqEqr_pIn_rOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, q, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 1 and nc == 1)
        {
            int r = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqEqr_prIn.at(std::make_pair(p, r));
            for (auto &doublet : doublets)
            {
                int q = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, q, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqErp(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right, const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 2)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            vector<tuple<int, int, double>> triplets = intlist_EpqErp_pqrOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, p, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int q = a_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqErp_qIn_rOut.at(q);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, p, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            int r = c_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqErp_rIn_qOut.at(r);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int q = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, p, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 1 and nc == 1)
        {
            int q = a_idxs[0];
            int r = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqErp_qrIn.at(std::make_pair(q, r));
            for (auto &doublet : doublets)
            {
                int p = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, p, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqEpq(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right, const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 4)
        return;

    size_t nue_right = count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            vector<tuple<int, double>> doublets = intlist_EpqEpq_pqOut.at(q);
            for (auto &doublet : doublets)
            {
                int p = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 2 and (a_idxs[0] == a_idxs[1]))
        {
            int q = a_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqEpq_qIn_pOut.at(q);
            for (auto &doublet : doublets)
            {
                int p = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 2 and (c_idxs[0] == c_idxs[1]))
        {
            int p = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqEpq_pIn_qOut.at(p);
            for (auto &doublet : doublets)
            {
                int q = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (4):
    {
        if (na == 2 and nc == 2)
        {
            if (a_idxs[0] == a_idxs[1] and c_idxs[0] == c_idxs[1])
            {
                int p = c_idxs[0];
                int q = a_idxs[0];
                double two_el_int = std::abs(impl->two_el_ints(p, q, p, q));
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
            }
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqEpr(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right,
                                                     const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 4)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            vector<tuple<int, int, double>> triplets = intlist_EpqEpr_pqrOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int q = a_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqEpr_qIn_prOut.at(q);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 2 and a_idxs[0] != a_idxs[1])
        {
            int q = a_idxs[0];
            int r = a_idxs[1];
            vector<tuple<int, double>> doublets = intlist_EpqEpr_qrIn_pOut.at(std::make_pair(q, r));
            for (auto &doublet : doublets)
            {
                int p = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 2 and c_idxs[0] == c_idxs[1])
        {
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqEpr_pIn_qrOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (3):
    {
        if ((na == 1 and nc == 2) and (c_idxs[0] == c_idxs[1]))
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqEpr_pqIn_rOut.at(std::make_pair(p, q));
            for (auto &doublet : doublets)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (4):
    {
        if ((na == 2 and nc == 2) and (a_idxs[0] != a_idxs[1]) and (c_idxs[0] == c_idxs[1]))
        {
            int q = a_idxs[0];
            int r = a_idxs[1];
            int p = c_idxs[0];

            double two_el_int = std::abs(impl->two_el_ints(p, q, p, r));
            if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                innerPrefixHelper2El(p, q, p, r, icfg_right, nue_right, onv_right,
                                     sfs_right, wfn_right, data_fois, data_var);
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqErq(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right,
                                                     const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 4)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t q = prefix_size; q < n_orbs; q++)
        {
            vector<tuple<int, int, double>> triplets = intlist_EpqErq_pqrOut.at(q);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (1):
    {
        if (nc == 1)
        {
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqErq_pIn_qrOut.at(p);
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 2 and (a_idxs[0] == a_idxs[1]))
        {
            int q = a_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqErq_qIn_prOut.at(q);
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 2 and (c_idxs[0] != c_idxs[1]))
        {
            int p = c_idxs[0];
            int r = c_idxs[1];
            vector<tuple<int, double>> doublets = intlist_EpqErq_prIn_qOut.at(std::make_pair(p, r));
            for (auto &doublet : doublets)
            {
                int q = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (3):
    {
        if ((na == 2 and nc == 1) and (a_idxs[0] == a_idxs[1]))
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqErq_pqIn_rOut.at(std::make_pair(p, q));
            for (auto &doublet : doublets)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (4):
    {
        if ((na == 2 and nc == 2) and (a_idxs[0] == a_idxs[1]) and (c_idxs[0] != c_idxs[1]))
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            int r = c_idxs[1];

            double two_el_int = std::abs(impl->two_el_ints(p, q, r, q));
            if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                innerPrefixHelper2El(p, q, r, q, icfg_right, nue_right, onv_right,
                                     sfs_right, wfn_right, data_fois, data_var);
        }
        break;
    }
    default:
        break;
    }
}

void GCI::Impl::PrefixAlgorithm::prefixHelper_EpqErs(const double &max_ci_coeff, const int &na, const int &nc,
                                                     const int &occ_diff_sum, const int &icfg_right,
                                                     const string &onv_right,
                                                     const vector<int> &a_idxs,
                                                     const vector<int> &c_idxs,
                                                     const vector<string> &sfs_right,
                                                     const wfn_ptr &wfn_right,
                                                     DataFOIS &data_fois, DataVar &data_var)
{
    if (occ_diff_sum > 4)
        return;

    size_t nue_right = std::count(onv_right.begin(), onv_right.end(), '1');

    switch (occ_diff_sum)
    {
    case (0):
    {
        for (size_t p = prefix_size; p < n_orbs; p++)
        {
            for (size_t q = prefix_size; q < n_orbs; q++)
            {
                if (p == q)
                    continue;

                vector<tuple<int, int, double>> triplets = intlist_EpqErs_pqrsOut.at(std::make_pair(p, q));
                for (auto &triplet : triplets)
                {
                    int r = get<0>(triplet);
                    int s = get<1>(triplet);
                    double two_el_int = get<2>(triplet);
                    if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                        innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                             sfs_right, wfn_right, data_fois, data_var);
                    else
                        break;
                }
            }
        }
        break;
    }
    case (1):
    {
        if (na == 1)
        {
            int s = a_idxs[0];
            vector<tuple<int, int, int, double>> quadruplets = intlist_EpqErs_sIn_pqrOut.at(s);
            for (auto &quadruplet : quadruplets)
            {
                int p = get<0>(quadruplet);
                int q = get<1>(quadruplet);
                int r = get<2>(quadruplet);
                double two_el_int = get<3>(quadruplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }

            int q = s;
            quadruplets = intlist_EpqErs_qIn_prsOut.at(q);
            for (auto &quadruplet : quadruplets)
            {
                int p = get<0>(quadruplet);
                int r = get<1>(quadruplet);
                int s = get<2>(quadruplet);
                double two_el_int = get<3>(quadruplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 1)
        {
            int p = c_idxs[0];
            vector<tuple<int, int, int, double>> quadruplets = intlist_EpqErs_pIn_qrsOut.at(p);
            for (auto &quadruplet : quadruplets)
            {
                int q = get<0>(quadruplet);
                int r = get<1>(quadruplet);
                int s = get<2>(quadruplet);
                double two_el_int = get<3>(quadruplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (2):
    {
        if (na == 2 and a_idxs[0] != a_idxs[1])
        {
            int q = a_idxs[0];
            int s = a_idxs[1];
            vector<tuple<int, int, double>> triplets = intlist_EpqErs_qsIn_prOut.at(std::make_pair(q, s));
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }

            triplets = intlist_EpqErs_qsIn_prOut.at(std::make_pair(s, q));
            for (auto &triplet : triplets)
            {
                int p = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, s, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (nc == 2 and c_idxs[0] != c_idxs[1])
        {
            int p = c_idxs[0];
            int r = c_idxs[1];
            vector<tuple<int, int, double>> triplets = intlist_EpqErs_prIn_qsOut.at(std::make_pair(p, r));
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int s = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if (na == 1 and nc == 1)
        {
            int q = a_idxs[0];
            int p = c_idxs[0];
            vector<tuple<int, int, double>> triplets = intlist_EpqErs_pqIn_rsOut.at(std::make_pair(p, q));
            for (auto &triplet : triplets)
            {
                int r = get<0>(triplet);
                int s = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }

            int s = a_idxs[0];
            triplets = intlist_EpqErs_psIn_qrOut.at(std::make_pair(p, s));
            for (auto &triplet : triplets)
            {
                int q = get<0>(triplet);
                int r = get<1>(triplet);
                double two_el_int = get<2>(triplet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (3):
    {
        if ((na == 2 and nc == 1) and (a_idxs[0] != a_idxs[1]))
        {
            int q = a_idxs[0];
            int s = a_idxs[1];
            int p = c_idxs[0];
            vector<tuple<int, double>> doublets = intlist_EpqErs_pqsIn_rOut.at(std::make_tuple(p, q, s));
            for (auto &doublet : doublets)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }

            doublets = intlist_EpqErs_pqsIn_rOut.at(std::make_tuple(p, s, q));
            for (auto &doublet : doublets)
            {
                int r = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, s, r, q, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }

        if ((na == 1 and nc == 2) and (c_idxs[0] != c_idxs[1]))
        {
            int s = a_idxs[0];
            int p = c_idxs[0];
            int r = c_idxs[1];
            vector<tuple<int, double>> doublets = intlist_EpqErs_prsIn_qOut.at(std::make_tuple(p, r, s));
            for (auto &doublet : doublets)
            {
                int q = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }

            int q = s;
            doublets = intlist_EpqErs_pqrIn_sOut.at(std::make_tuple(p, q, r));
            for (auto &doublet : doublets)
            {
                int s = get<0>(doublet);
                double two_el_int = get<1>(doublet);
                if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                    innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                         sfs_right, wfn_right, data_fois, data_var);
                else
                    break;
            }
        }
        break;
    }
    case (4):
    {
        if ((na == 2 and nc == 2) and (a_idxs[0] != a_idxs[1]) and (c_idxs[0] != c_idxs[1]))
        {
            int q = a_idxs[0];
            int s = a_idxs[1];
            int p = c_idxs[0];
            int r = c_idxs[1];

            double two_el_int = std::abs(impl->two_el_ints(p, q, r, s));
            if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                innerPrefixHelper2El(p, q, r, s, icfg_right, nue_right, onv_right,
                                     sfs_right, wfn_right, data_fois, data_var);

            two_el_int = std::abs(impl->two_el_ints(p, s, r, q));
            if (max_ci_coeff * two_el_int > Settings::getEpsilonVar())
                innerPrefixHelper2El(p, s, r, q, icfg_right, nue_right, onv_right,
                                     sfs_right, wfn_right, data_fois, data_var);
        }
        break;
    }
    default:
        break;
    }
}