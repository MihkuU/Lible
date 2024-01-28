#include <lible/connections.hpp>
#include <lible/util.h>

#ifdef _USE_MPI_
#include <lible/brain.hpp>
#endif

#include <algorithm>
#include <omp.h>

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::string;
using std::vector;

void GCI::Impl::Connections::constructConnections(const std::set<std::string> &cfgs_right,
                                                  const wfn_ptr &wfn_left,
                                                  const wfn_ptr &wfn_right,
                                                  connection_map_1el &connections_1el,
                                                  connection_map_2el &connections_2el,
                                                  connection_map_dia &connections_dia)
{
#pragma omp parallel
    {
        connection_map_1el connections_1el_omp;
        connection_map_2el connections_2el_omp;
        connection_map_dia connections_dia_omp;

        int rank_total, size_total;
#ifdef _USE_MPI_
        rank_total = Brain::returnTotalRank();
        size_total = Brain::returnTotalSize();
#else
        rank_total = omp_get_thread_num();
        size_total = omp_get_num_threads();
#endif

        int ipal = 0;
        for (string cfg_right : cfgs_right)
        {
            if (ipal % size_total != rank_total)
            {
                ipal++;
                continue;
            }
            ipal++;

            int icfg_right = wfn_right->findCFGPos(cfg_right);

            /* One-electron excitations */
            // Epq
            constructConnections_Epq(icfg_right, cfg_right, wfn_left, connections_1el_omp);

            // EpqEqp, one-electron excitation to RI-space and back
            constructConnections_EpqEqp(icfg_right, cfg_right, connections_dia_omp);

            // EpqEqr
            constructConnections_EpqEqr(icfg_right, cfg_right, "0ac", wfn_left, connections_2el_omp);
            constructConnections_EpqEqr(icfg_right, cfg_right, "0ca", wfn_left, connections_2el_omp);
            constructConnections_EpqEqr(icfg_right, cfg_right, "a0c", wfn_left, connections_2el_omp);
            constructConnections_EpqEqr(icfg_right, cfg_right, "c0a", wfn_left, connections_2el_omp);
            constructConnections_EpqEqr(icfg_right, cfg_right, "ac0", wfn_left, connections_2el_omp);
            constructConnections_EpqEqr(icfg_right, cfg_right, "ca0", wfn_left, connections_2el_omp);

            // EpqErp
            constructConnections_EpqErp(icfg_right, cfg_right, "0ac", wfn_left, connections_2el_omp);
            constructConnections_EpqErp(icfg_right, cfg_right, "0ca", wfn_left, connections_2el_omp);
            constructConnections_EpqErp(icfg_right, cfg_right, "a0c", wfn_left, connections_2el_omp);
            constructConnections_EpqErp(icfg_right, cfg_right, "c0a", wfn_left, connections_2el_omp);
            constructConnections_EpqErp(icfg_right, cfg_right, "ac0", wfn_left, connections_2el_omp);
            constructConnections_EpqErp(icfg_right, cfg_right, "ca0", wfn_left, connections_2el_omp);

            /* Two-electron excitations */
            // Enforcing here pq < rs

            // EpqEpq
            constructConnections_EpqEpq(icfg_right, cfg_right, wfn_left, connections_2el_omp);

            // EpqEpr
            constructConnections_EpqEpr(icfg_right, cfg_right, "caa", wfn_left, connections_2el_omp);
            constructConnections_EpqEpr(icfg_right, cfg_right, "aca", wfn_left, connections_2el_omp);
            constructConnections_EpqEpr(icfg_right, cfg_right, "aac", wfn_left, connections_2el_omp);

            // EpqErq
            constructConnections_EpqErq(icfg_right, cfg_right, "cca", wfn_left, connections_2el_omp);
            constructConnections_EpqErq(icfg_right, cfg_right, "cac", wfn_left, connections_2el_omp);
            constructConnections_EpqErq(icfg_right, cfg_right, "acc", wfn_left, connections_2el_omp);

            // EpqErs
            constructConnections_EpqErs(icfg_right, cfg_right, "caca", wfn_left, connections_2el_omp);
            constructConnections_EpqErs(icfg_right, cfg_right, "ccaa", wfn_left, connections_2el_omp);
            constructConnections_EpqErs(icfg_right, cfg_right, "caac", wfn_left, connections_2el_omp);
            constructConnections_EpqErs(icfg_right, cfg_right, "acca", wfn_left, connections_2el_omp);
            constructConnections_EpqErs(icfg_right, cfg_right, "aacc", wfn_left, connections_2el_omp);
            constructConnections_EpqErs(icfg_right, cfg_right, "acac", wfn_left, connections_2el_omp);
        }
#pragma omp critical
        {
            mergeConnections(connections_1el_omp, connections_1el);
            mergeConnections(connections_2el_omp, connections_2el);
            mergeConnections(connections_dia_omp, connections_dia);
        }
    }
}

void GCI::Impl::Connections::constructConnections_Epq(const int &icfg_right, const string &cfg_right,
                                                      const wfn_ptr &wfn_left, connection_map_1el &connections_1el)
{
    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t q = 0; q < n_orbs; q++)
    {
        if (testAnnihilation(q, cfg_right))
            continue;

        string cfg = cfg_right;
        annihilation(q, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(q + 1, cfg);
        if (node == nullptr)
            continue;

        size_t q_ = q + 1;
        for (size_t p = q + 1; p < n_orbs; p++)
        {
            if (p > q_)
                node = wfn_left->incrementNode(p - 1, cfg, node);
            if (node == nullptr)
                break;

            if (testCreation(p, cfg))
                continue;

            string cfg_left = cfg;
            creation(p, cfg_left);

            int icfg_left = wfn_left->findCFGPos(p, cfg_left, node);
            if (icfg_left == -1 or icfg_left > icfg_right)
                continue;

            int nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');
            quintet info_cc = returnCCInfo(p, q, nue_left, nue_right, cfg_left, cfg_right);

            size_t pq = pq2DTo1D(p, q, n_orbs);

            bool phase = determine1ElPhase(p, q, cfg_right);

            connections_1el[info_cc].push_back(std::make_tuple(icfg_left, icfg_right, pq, phase));
        }
    }

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (testCreation(p, cfg_right))
            continue;

        string cfg = cfg_right;
        creation(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        size_t p_ = p + 1;
        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p_)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (testAnnihilation(q, cfg))
                continue;

            string cfg_left = cfg;
            annihilation(q, cfg_left);

            int icfg_left = wfn_left->findCFGPos(q, cfg_left, node);
            if (icfg_left == -1 or icfg_left > icfg_right)
                continue;

            int nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');
            quintet info_cc = returnCCInfo(p, q, nue_left, nue_right, cfg_left, cfg_right);

            size_t pq = pq2DTo1D(p, q, n_orbs);
            bool phase = determine1ElPhase(p, q, cfg_right);

            connections_1el[info_cc].push_back(std::make_tuple(icfg_left, icfg_right, pq, phase));
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqEqp(const int &icfg_right, const std::string &cfg_right,
                                                         connection_map_dia &connections_dia)
{
    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (cfg_right[p] == '0')
            continue;

        for (size_t q = 0; q < n_orbs; q++)
        {
            if (p == q)
                continue;
            if (cfg_right[q] == '2')
                continue;

            string cfg_ri = cfg_right;
            cfg_ri[p] -= 1;
            cfg_ri[q] += 1;
            size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');

            if (nue_ri < min_nue)
                continue;

            quintet key = returnCCInfo(q, p, nue_ri, nue_right, cfg_ri, cfg_right);
            size_t pqqp = pqrs4DTo1D(p, q, q, p, n_orbs);

            connections_dia[key].push_back(std::make_tuple(icfg_right, pqqp));
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqEqr(const int &icfg_right, const std::string &cfg_right,
                                                         const std::string &operators, const wfn_ptr &wfn_left,
                                                         connection_map_2el &connections_2el)
{
    char pOp = operators[0];
    char qOp = operators[1];
    char rOp = operators[2];
    std::function<void(const int &p, string &cfg)> pOperate;
    std::function<void(const int &q, string &cfg)> qOperate;
    std::function<void(const int &r, string &cfg)> rOperate;
    if (pOp == '0')
        pOperate = identity;
    else if (pOp == 'a')
        pOperate = annihilation;
    else if (pOp == 'c')
        pOperate = creation;

    if (qOp == '0')
        qOperate = identity;
    else if (qOp == 'a')
        qOperate = annihilation;
    else if (qOp == 'c')
        qOperate = creation;

    if (rOp == '0')
        rOperate = identity;
    else if (rOp == 'a')
        rOperate = annihilation;
    else if (rOp == 'c')
        rOperate = creation;

    std::function<bool(const int &p, const string &cfg)> pTest;
    std::function<bool(const int &q, const string &cfg)> qTest;
    std::function<bool(const int &r, const string &cfg)> rTest;
    if (pOp == '0')
        pTest = testCreation;
    else if (pOp == 'a')
        pTest = testAnnihilation;
    else if (pOp == 'c')
        pTest = testCreation;

    if (qOp == '0')
        qTest = testCreation;
    else if (qOp == 'a')
        qTest = testAnnihilation;
    else if (qOp == 'c')
        qTest = testCreation;

    if (rOp == '0')
        rTest = testCreation;
    else if (rOp == 'a')
        rTest = testAnnihilation;
    else if (rOp == 'c')
        rTest = testCreation;

    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (pTest(p, cfg_right))
            continue;

        string cfg = cfg_right;
        pOperate(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (qTest(q, cfg))
                continue;

            string cfg_ = cfg;
            qOperate(q, cfg_);

            CFGTree::Node *node_ = wfn_left->incrementNode(q, cfg_, node);
            if (node_ == nullptr)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r > q + 1)
                    node_ = wfn_left->incrementNode(r - 1, cfg_, node_);

                if (node_ == nullptr)
                    break;

                if (rTest(r, cfg_))
                    continue;

                string cfg_left = cfg_;
                rOperate(r, cfg_left);

                int icfg_left = wfn_left->findCFGPos(r, cfg_left, node_);
                if (icfg_left == -1 or icfg_left > icfg_right)
                    continue;

                size_t nue_left = count(cfg_left.begin(), cfg_left.end(), '1');

                vector<int> pqrsCanonical = canonicalizeIdxs_EpqEqr(operators, vector<int>({int(p), int(q), int(r)}));
                int p_ = pqrsCanonical[0];
                int q_ = pqrsCanonical[1];
                int r_ = pqrsCanonical[2];
                int s_ = pqrsCanonical[3];

                string cfg_ri = cfg_right;
                cfg_ri[r_] += 1;
                cfg_ri[s_] -= 1;
                size_t nue_ri = count(cfg_ri.begin(), cfg_ri.end(), '1');

                if (nue_ri < min_nue)
                    continue;

                bool phase = determine2ElPhase(p_, q_, r_, s_, cfg_ri, cfg_right);

                quintet info_cc1 = returnCCInfo(p_, q_, nue_left, nue_ri, cfg_left, cfg_ri);
                quintet info_cc2 = returnCCInfo(r_, s_, nue_ri, nue_right, cfg_ri, cfg_right);

                int exctype1 = get<0>(info_cc1);
                int exctype2 = get<0>(info_cc2);
                int prel1 = get<3>(info_cc1);
                int qrel1 = get<4>(info_cc1);
                int prel2 = get<3>(info_cc2);
                int qrel2 = get<4>(info_cc2);

                size_t pqqr = pqrs4DTo1D(p_, q_, r_, s_, n_orbs);

                nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
                connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, pqqr, phase, false));
            }
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqErp(const int &icfg_right, const std::string &cfg_right,
                                                         const std::string &operators, const wfn_ptr &wfn_left,
                                                         connection_map_2el &connections_2el)
{
    char pOp = operators[0];
    char qOp = operators[1];
    char rOp = operators[2];
    std::function<void(const int &p, string &cfg)> pOperate;
    std::function<void(const int &q, string &cfg)> qOperate;
    std::function<void(const int &r, string &cfg)> rOperate;
    if (pOp == '0')
        pOperate = identity;
    else if (pOp == 'a')
        pOperate = annihilation;
    else if (pOp == 'c')
        pOperate = creation;

    if (qOp == '0')
        qOperate = identity;
    else if (qOp == 'a')
        qOperate = annihilation;
    else if (qOp == 'c')
        qOperate = creation;

    if (rOp == '0')
        rOperate = identity;
    else if (rOp == 'a')
        rOperate = annihilation;
    else if (rOp == 'c')
        rOperate = creation;

    std::function<bool(const int &p, const string &cfg)> pTest;
    std::function<bool(const int &q, const string &cfg)> qTest;
    std::function<bool(const int &r, const string &cfg)> rTest;
    if (pOp == '0')
        pTest = testAnnihilation;
    else if (pOp == 'a')
        pTest = testAnnihilation;
    else if (pOp == 'c')
        pTest = testCreation;

    if (qOp == '0')
        qTest = testAnnihilation;
    else if (qOp == 'a')
        qTest = testAnnihilation;
    else if (qOp == 'c')
        qTest = testCreation;

    if (rOp == '0')
        rTest = testAnnihilation;
    else if (rOp == 'a')
        rTest = testAnnihilation;
    else if (rOp == 'c')
        rTest = testCreation;

    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (pTest(p, cfg_right))
            continue;
        string cfg = cfg_right;
        pOperate(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);

            if (node == nullptr)
                break;

            if (qTest(q, cfg))
                continue;

            string cfg_ = cfg;
            qOperate(q, cfg_);

            CFGTree::Node *node_ = wfn_left->incrementNode(q, cfg_, node);
            if (node_ == nullptr)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r > q + 1)
                    node_ = wfn_left->incrementNode(r - 1, cfg_, node_);
                if (node_ == nullptr)
                    break;

                if (rTest(r, cfg_))
                    continue;

                string cfg_left = cfg_;
                rOperate(r, cfg_left);

                int icfg_left = wfn_left->findCFGPos(r, cfg_left, node_);
                if (icfg_left == -1 or icfg_left > icfg_right)
                    continue;

                size_t nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');

                vector<int> pqrsCanonical = canonicalizeIdxs_EpqErp(operators, vector<int>({int(p), int(q), int(r)}));
                int p_ = pqrsCanonical[0];
                int q_ = pqrsCanonical[1];
                int r_ = pqrsCanonical[2];
                int s_ = pqrsCanonical[3];

                string cfg_ri = cfg_right;
                cfg_ri[r_] += 1;
                cfg_ri[s_] -= 1;
                size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');
                if (nue_ri < min_nue)
                    continue;

                bool phase = determine2ElPhase(p_, q_, r_, s_, cfg_ri, cfg_right);

                quintet info_cc1 = returnCCInfo(p_, q_, nue_left, nue_ri, cfg_left, cfg_ri);
                quintet info_cc2 = returnCCInfo(r_, s_, nue_ri, nue_right, cfg_ri, cfg_right);

                int exctype1 = get<0>(info_cc1);
                int exctype2 = get<0>(info_cc2);
                int prel1 = get<3>(info_cc1);
                int qrel1 = get<4>(info_cc1);
                int prel2 = get<3>(info_cc2);
                int qrel2 = get<4>(info_cc2);

                size_t pqrp = pqrs4DTo1D(p_, q_, r_, s_, n_orbs);

                nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
                connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, pqrp, phase, false));
            }
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqEpq(const int &icfg_right, const std::string &cfg_right,
                                                         const wfn_ptr &wfn_left, connection_map_2el &connections_2el)
{
    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t q = 0; q < n_orbs; q++)
    {
        if (testAnnihilation2x(q, cfg_right))
            continue;

        string cfg = cfg_right;
        annihilation2x(q, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(q + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t p = q + 1; p < n_orbs; p++)
        {
            if (p > q + 1)
                node = wfn_left->incrementNode(p - 1, cfg, node);
            if (node == nullptr)
                break;

            if (testCreation2x(p, cfg))
                continue;

            string cfg_left = cfg;
            creation2x(p, cfg_left);

            int icfg_left = wfn_left->findCFGPos(p, cfg_left, node);
            if (icfg_left == -1 or icfg_left > icfg_right)
                continue;

            string cfg_ri = cfg_right;
            cfg_ri[p] += 1;
            cfg_ri[q] -= 1;

            size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');
            if (nue_ri < min_nue)
                continue;

            size_t nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');

            quintet info_ccxcc1 = returnCCInfo(p, q, nue_left, nue_ri, cfg_left, cfg_ri);
            quintet info_ccxcc2 = returnCCInfo(p, q, nue_ri, nue_right, cfg_ri, cfg_right);

            int exctype1 = get<0>(info_ccxcc1);
            int exctype2 = get<0>(info_ccxcc2);
            int prel1 = get<3>(info_ccxcc1);
            int qrel1 = get<4>(info_ccxcc1);
            int prel2 = get<3>(info_ccxcc2);
            int qrel2 = get<4>(info_ccxcc2);

            size_t pqpq = pqrs4DTo1D(p, q, p, q, n_orbs);

            nonet info_ccxcc = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
            connections_2el[info_ccxcc].push_back(std::make_tuple(icfg_left, icfg_right, pqpq, false, false));
        }
    }

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (testCreation2x(p, cfg_right))
            continue;

        string cfg = cfg_right;
        creation2x(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (testAnnihilation2x(q, cfg))
                continue;

            string cfg_left = cfg;
            annihilation2x(q, cfg_left);

            int icfg_left = wfn_left->findCFGPos(q, cfg_left, node);
            if (icfg_left == -1 or icfg_left > icfg_right)
                continue;

            string cfg_ri = cfg_right;
            cfg_ri[p] += 1;
            cfg_ri[q] -= 1;

            size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');
            if (nue_ri < min_nue)
                continue;

            size_t nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');

            quintet info_cc1 = returnCCInfo(p, q, nue_left, nue_ri, cfg_left, cfg_ri);
            quintet info_cc2 = returnCCInfo(p, q, nue_ri, nue_right, cfg_ri, cfg_right);

            int exctype1 = get<0>(info_cc1);
            int exctype2 = get<0>(info_cc2);
            int prel1 = get<3>(info_cc1);
            int qrel1 = get<4>(info_cc1);
            int prel2 = get<3>(info_cc2);
            int qrel2 = get<4>(info_cc2);

            size_t pqpq = pqrs4DTo1D(p, q, p, q, n_orbs);

            nonet info_ccxcc = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
            connections_2el[info_ccxcc].push_back(std::make_tuple(icfg_left, icfg_right, pqpq, false, false));
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqEpr(const int &icfg_right, const std::string &cfg_right,
                                                         const std::string &operators, const wfn_ptr &wfn_left,
                                                         connection_map_2el &connections_2el)
{
    char p_op = operators[0];
    char q_op = operators[1];
    char r_op = operators[2];
    std::function<void(const int &p, string &cfg)> pOperate;
    std::function<void(const int &q, string &cfg)> qOperate;
    std::function<void(const int &r, string &cfg)> rOperate;
    if (p_op == 'c')
        pOperate = creation2x;
    else if (p_op == 'a')
        pOperate = annihilation;

    if (q_op == 'c')
        qOperate = creation2x;
    else if (q_op == 'a')
        qOperate = annihilation;

    if (r_op == 'c')
        rOperate = creation2x;
    else if (r_op == 'a')
        rOperate = annihilation;

    std::function<bool(const int &p, const string &cfg)> pTest;
    std::function<bool(const int &q, const string &cfg)> qTest;
    std::function<bool(const int &r, const string &cfg)> rTest;
    if (p_op == 'c')
        pTest = testCreation2x;
    else if (p_op == 'a')
        pTest = testAnnihilation;

    if (q_op == 'c')
        qTest = testCreation2x;
    else if (q_op == 'a')
        qTest = testAnnihilation;

    if (r_op == 'c')
        rTest = testCreation2x;
    else if (r_op == 'a')
        rTest = testAnnihilation;

    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (pTest(p, cfg_right))
            continue;
        string cfg = cfg_right;
        pOperate(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (qTest(q, cfg))
                continue;
            string cfg_ = cfg;
            qOperate(q, cfg_);

            CFGTree::Node *node_ = wfn_left->incrementNode(q, cfg_, node);
            if (node_ == nullptr)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r > q + 1)
                    node_ = wfn_left->incrementNode(r - 1, cfg_, node_);
                if (node_ == nullptr)
                    break;

                if (rTest(r, cfg_))
                    continue;

                string cfg_left = cfg_;
                rOperate(r, cfg_left);

                int icfg_left = wfn_left->findCFGPos(r, cfg_left, node_);
                if (icfg_left == -1 or icfg_left > icfg_right)
                    continue;

                size_t nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');

                vector<int> pqrsCanonical = canonicalizeIdxs_EpqEpr(operators, vector<int>({int(p), int(q), int(r)}));
                int p_ = pqrsCanonical[0];
                int q_ = pqrsCanonical[1];
                int r_ = pqrsCanonical[2];
                int s_ = pqrsCanonical[3];

                string cfg_ri = cfg_right;
                cfg_ri[r_] += 1;
                cfg_ri[s_] -= 1;
                size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');

                if (nue_ri < min_nue)
                    continue;

                bool phase = determine2ElPhase(p_, q_, r_, s_, cfg_ri, cfg_right);

                quintet info_cc1 = returnCCInfo(p_, q_, nue_left, nue_ri, cfg_left, cfg_ri);
                quintet info_cc2 = returnCCInfo(r_, s_, nue_ri, nue_right, cfg_ri, cfg_right);

                int exctype1 = get<0>(info_cc1);
                int exctype2 = get<0>(info_cc2);
                int prel1 = get<3>(info_cc1);
                int qrel1 = get<4>(info_cc1);
                int prel2 = get<3>(info_cc2);
                int qrel2 = get<4>(info_cc2);

                size_t pqpr = pqrs4DTo1D(p_, q_, r_, s_, n_orbs);

                nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
                connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, pqpr, phase, true));
            }
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqErq(const int &icfg_right, const std::string &cfg_right,
                                                         const std::string &operators, const wfn_ptr &wfn_left,
                                                         connection_map_2el &connections_2el)
{
    char pOp = operators[0];
    char qOp = operators[1];
    char rOp = operators[2];
    std::function<void(const int &p, string &cfg)> pOperate;
    std::function<void(const int &q, string &cfg)> qOperate;
    std::function<void(const int &r, string &cfg)> rOperate;
    if (pOp == 'c')
        pOperate = creation;
    else if (pOp == 'a')
        pOperate = annihilation2x;

    if (qOp == 'c')
        qOperate = creation;
    else if (qOp == 'a')
        qOperate = annihilation2x;

    if (rOp == 'c')
        rOperate = creation;
    else if (rOp == 'a')
        rOperate = annihilation2x;

    std::function<bool(const int &p, const string &cfg)> pTest;
    std::function<bool(const int &q, const string &cfg)> qTest;
    std::function<bool(const int &r, const string &cfg)> rTest;
    if (pOp == 'c')
        pTest = testCreation;
    else if (pOp == 'a')
        pTest = testAnnihilation2x;

    if (qOp == 'c')
        qTest = testCreation;
    else if (qOp == 'a')
        qTest = testAnnihilation2x;

    if (rOp == 'c')
        rTest = testCreation;
    else if (rOp == 'a')
        rTest = testAnnihilation2x;

    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (pTest(p, cfg_right))
            continue;

        string cfg = cfg_right;
        pOperate(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (qTest(q, cfg))
                continue;
            string cfg_ = cfg;
            qOperate(q, cfg_);

            CFGTree::Node *node_ = wfn_left->incrementNode(q, cfg_, node);
            if (node_ == nullptr)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r > q + 1)
                    node_ = wfn_left->incrementNode(r - 1, cfg_, node_);
                if (node_ == nullptr)
                    break;

                if (rTest(r, cfg_))
                    continue;

                string cfg_left = cfg_;
                rOperate(r, cfg_left);

                int icfg_left = wfn_left->findCFGPos(r, cfg_left, node_);
                if (icfg_left == -1 or icfg_left > icfg_right)
                    continue;

                size_t nue_left = std::count(cfg_left.begin(), cfg_left.end(), '1');

                vector<int> pqrsCanonical = canonicalizeIdxs_EpqErq(operators, vector<int>({int(p), int(q), int(r)}));
                int p_ = pqrsCanonical[0];
                int q_ = pqrsCanonical[1];
                int r_ = pqrsCanonical[2];
                int s_ = pqrsCanonical[3];

                string cfg_ri = cfg_right;
                cfg_ri[r_] += 1;
                cfg_ri[s_] -= 1;
                size_t nue_ri = std::count(cfg_ri.begin(), cfg_ri.end(), '1');

                if (nue_ri < min_nue)
                    continue;

                bool phase = determine2ElPhase(p_, q_, r_, s_, cfg_ri, cfg_right);

                quintet info_cc1 = returnCCInfo(p_, q_, nue_left, nue_ri, cfg_left, cfg_ri);
                quintet info_cc2 = returnCCInfo(r_, s_, nue_ri, nue_right, cfg_ri, cfg_right);

                int exctype1 = get<0>(info_cc1);
                int exctype2 = get<0>(info_cc2);
                int prel1 = get<3>(info_cc1);
                int qrel1 = get<4>(info_cc1);
                int prel2 = get<3>(info_cc2);
                int qrel2 = get<4>(info_cc2);

                size_t pqrq = pqrs4DTo1D(p_, q_, r_, s_, n_orbs);

                nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri, nue_right, prel1, qrel1, prel2, qrel2);
                connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, pqrq, phase, true));
            }
        }
    }
}

void GCI::Impl::Connections::constructConnections_EpqErs(const int &icfg_right, const std::string &cfg_right,
                                                         const std::string &operators, const wfn_ptr &wfn_left,
                                                         connection_map_2el &connections_2el)
{
    char pOp = operators[0];
    char qOp = operators[1];
    char rOp = operators[2];
    char sOp = operators[3];
    std::function<void(const int &p, string &cfg)> pOperate;
    std::function<void(const int &q, string &cfg)> qOperate;
    std::function<void(const int &r, string &cfg)> rOperate;
    std::function<void(const int &s, string &cfg)> sOperate;
    if (pOp == 'a')
        pOperate = annihilation;
    else if (pOp == 'c')
        pOperate = creation;

    if (qOp == 'a')
        qOperate = annihilation;
    else if (qOp == 'c')
        qOperate = creation;

    if (rOp == 'a')
        rOperate = annihilation;
    else if (rOp == 'c')
        rOperate = creation;

    if (sOp == 'a')
        sOperate = annihilation;
    else if (sOp == 'c')
        sOperate = creation;

    std::function<bool(const int &p, const string &cfg)> pTest;
    std::function<bool(const int &q, const string &cfg)> qTest;
    std::function<bool(const int &r, const string &cfg)> rTest;
    std::function<bool(const int &s, const string &cfg)> sTest;
    if (pOp == 'a')
        pTest = testAnnihilation;
    else if (pOp == 'c')
        pTest = testCreation;

    if (qOp == 'a')
        qTest = testAnnihilation;
    else if (qOp == 'c')
        qTest = testCreation;

    if (rOp == 'a')
        rTest = testAnnihilation;
    else if (rOp == 'c')
        rTest = testCreation;

    if (sOp == 'a')
        sTest = testAnnihilation;
    else if (sOp == 'c')
        sTest = testCreation;

    size_t nue_right = std::count(cfg_right.begin(), cfg_right.end(), '1');

    for (size_t p = 0; p < n_orbs; p++)
    {
        if (pTest(p, cfg_right))
            continue;

        string cfg = cfg_right;
        pOperate(p, cfg);

        CFGTree::Node *node = wfn_left->searchFromRoot(p + 1, cfg);
        if (node == nullptr)
            continue;

        for (size_t q = p + 1; q < n_orbs; q++)
        {
            if (q > p + 1)
                node = wfn_left->incrementNode(q - 1, cfg, node);
            if (node == nullptr)
                break;

            if (qTest(q, cfg))
                continue;

            string cfg_ = cfg;
            qOperate(q, cfg_);

            CFGTree::Node *node_ = wfn_left->incrementNode(q, cfg_, node);
            if (node_ == nullptr)
                continue;

            for (size_t r = q + 1; r < n_orbs; r++)
            {
                if (r > q + 1)
                    node_ = wfn_left->incrementNode(r - 1, cfg_, node_);
                if (node_ == nullptr)
                    break;

                if (rTest(r, cfg_))
                    continue;

                string cfg__ = cfg_;
                rOperate(r, cfg__);

                CFGTree::Node *node__ = wfn_left->incrementNode(r, cfg__, node_);
                if (node__ == nullptr)
                    continue;

                for (size_t s = r + 1; s < n_orbs; s++)
                {
                    if (s > r + 1)
                        node__ = wfn_left->incrementNode(s - 1, cfg__, node__);
                    if (node__ == nullptr)
                        break;

                    if (sTest(s, cfg__))
                        continue;

                    string cfg_left = cfg__;
                    sOperate(s, cfg_left);

                    int icfg_left = wfn_left->findCFGPos(s, cfg_left, node__);
                    if (icfg_left == -1 or icfg_left > icfg_right)
                        continue;

                    size_t nue_left = count(cfg_left.begin(), cfg_left.end(), '1');

                    vector<int> pqrsCanonical = canonicalizeIdxs_EpqErs(operators, vector<int>({int(p), int(q), int(r), int(s)}));
                    int p_ = pqrsCanonical[0];
                    int q_ = pqrsCanonical[1];
                    int r_ = pqrsCanonical[2];
                    int s_ = pqrsCanonical[3];

                    string cfg_ri1 = cfg_right;
                    cfg_ri1[r_] += 1;
                    cfg_ri1[s_] -= 1;
                    size_t nue_ri1 = std::count(cfg_ri1.begin(), cfg_ri1.end(), '1');

                    if (nue_ri1 >= min_nue)
                    {
                        bool phase = determine2ElPhase(p_, q_, r_, s_, cfg_ri1, cfg_right);

                        quintet info_cc1 = returnCCInfo(p_, q_, nue_left, nue_ri1, cfg_left, cfg_ri1);
                        quintet info_cc2 = returnCCInfo(r_, s_, nue_ri1, nue_right, cfg_ri1, cfg_right);

                        int exctype1 = get<0>(info_cc1);
                        int exctype2 = get<0>(info_cc2);
                        int prel1 = get<3>(info_cc1);
                        int qrel1 = get<4>(info_cc1);
                        int prel2 = get<3>(info_cc2);
                        int qrel2 = get<4>(info_cc2);

                        size_t pqrs = pqrs4DTo1D(p_, q_, r_, s_, n_orbs);

                        nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri1, nue_right, prel1, qrel1, prel2, qrel2);
                        connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, pqrs, phase, true));
                    }

                    string cfg_ri2 = cfg_right;
                    cfg_ri2[r_] += 1;
                    cfg_ri2[q_] -= 1;
                    size_t nue_ri2 = std::count(cfg_ri2.begin(), cfg_ri2.end(), '1');

                    if (nue_ri2 >= min_nue)
                    {
                        bool phase = determine2ElPhase(p_, s_, r_, q_, cfg_ri2, cfg_right);

                        quintet info_cc1 = returnCCInfo(p_, s_, nue_left, nue_ri2, cfg_left, cfg_ri2);
                        quintet info_cc2 = returnCCInfo(r_, q_, nue_ri2, nue_right, cfg_ri2, cfg_right);

                        int exctype1 = get<0>(info_cc1);
                        int exctype2 = get<0>(info_cc2);
                        int prel1 = get<3>(info_cc1);
                        int qrel1 = get<4>(info_cc1);
                        int prel2 = get<3>(info_cc2);
                        int qrel2 = get<4>(info_cc2);

                        size_t psrq = pqrs4DTo1D(p_, s_, r_, q_, n_orbs);

                        nonet key = std::make_tuple(exctype1, exctype2, nue_left, nue_ri2, nue_right, prel1, qrel1, prel2, qrel2);
                        connections_2el[key].push_back(std::make_tuple(icfg_left, icfg_right, psrq, phase, true));
                    }
                }
            }
        }
    }
}