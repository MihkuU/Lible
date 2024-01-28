#include <lible/gci_util.hpp>

#include <stdexcept>

#include <fmt/core.h>

using namespace lible::guga;
using namespace lible::guga::util;

namespace GU = lible::guga::util;

using std::map;
using std::set;
using std::string;
using std::tuple;
using std::vector;

double GU::A(const int &b, const int &x, const int &y)
{
    return sqrt(double(b + x) / double(b + y));
}

double GU::C(const int &b, const int &x)
{
    return sqrt(double((b + x - 1) * (b + x + 1))) / double(b + x);
}

double GU::f(const char &d, const int &b)
{
    /*
     * ACHTUNG: sometimes the following commented code doesn't work. Works in
     * one computer, doesn't in the other. Compiler bug?
     */

    //    if (d == '0' or '3')
    //        return 1;
    //    else if (d == '1')
    //        return (A(b, 2, 0) * A(b, -1, 1));
    //    else if (d == '2')
    //        return (A(b, 0, 2) * A(b, 3, 1));

    if (d == '0' or d == '3')
        return 1;
    else if (d == '1')
        return A(b, 2, 0) * A(b, -1, 1);
    else if (d == '2')
        return A(b, 0, 2) * A(b, 3, 1);
    else
        throw std::runtime_error("Error: wrong step value d!");
}

double GU::fNS(const int &nue, const double &spin)
{
    return (tgamma(double(nue + 1)) / (tgamma(double(0.5 * nue - spin + 1)) * tgamma(double(0.5 * nue + spin + 1))) -
            tgamma(double(nue + 1)) / (tgamma(double(0.5 * nue - spin)) * tgamma(double(0.5 * nue + spin + 2))));
}

int GU::determineStepb(const char &d)
{
    if (d == '0')
        return 0;
    else if (d == '1')
        return 1;
    else if (d == '2')
        return -1;
    else if (d == '3')
        return 0;
    else
        throw std::runtime_error("Error: wrong step value d!");
}

bool GU::determine1ElPhase(const int &p, const int &q, const string &cfg_right)
{
    bool phase = false;

    int counter = 0;
    if (p < q)
        for (int i = p + 1; i < q; i++)
        {
            if (cfg_right[i] == '2')
                counter++;
        }
    else
        for (int i = q + 1; i < p; i++)
        {
            if (cfg_right[i] == '2')
                counter++;
        }

    if (counter % 2 == 1)
        phase = true;

    return phase;
}

bool GU::determine2ElPhase(const int &p, const int &q, const int &r, const int &s,
                           const string &conf_ri, const string &conf_right)
{
    bool total_phase = false;

    int phase1 = 1;
    int phase2 = 1;
    int count_2occ1 = 0;
    if (p < q)
    {
        for (int i = p + 1; i < q; i++)
            if (conf_ri[i] == '2')
                count_2occ1++;
    }
    else
    {
        for (int i = q + 1; i < p; i++)
            if (conf_ri[i] == '2')
                count_2occ1++;
    }

    if (count_2occ1 % 2 == 1)
        phase1 = -1;

    int count_2occ2 = 0;
    if (r < s)
    {
        for (int i = r + 1; i < s; i++)
            if (conf_right[i] == '2')
                count_2occ2++;
    }
    else
    {
        for (int i = s + 1; i < r; i++)
            if (conf_right[i] == '2')
                count_2occ2++;
    }

    if (count_2occ2 % 2 == 1)
        phase2 = -1;

    if (phase1 * phase2 == -1)
        total_phase = not(total_phase);

    return total_phase;
}

size_t GU::pq2DTo1D(const int &p, const int &q, const int &n_orb)
{
    return p * n_orb + q;
}

size_t GU::pqrs4DTo1D(const int &p, const int &q, const int &r, const int &s,
                      const int &n_orb)
{
    return p * pow(n_orb, 3) + q * pow(n_orb, 2) + r * n_orb + s;
}

tuple<int, int> GU::pq1DTo2D(const size_t &pq, const int &n_orbs)
{
    int p = pq / n_orbs;
    int q = pq % n_orbs;

    return {p, q};
}

tuple<int, int, int, int> GU::pqrs1DTo4D(const size_t &pqrs, const int &n_orbs)
{
    size_t x = pqrs;
    int p = x / pow(n_orbs, 3);
    x -= p * pow(n_orbs, 3);

    int q = x / pow(n_orbs, 2);
    x -= q * pow(n_orbs, 2);

    int r = x / n_orbs;
    x -= r * n_orbs;

    int s = x % n_orbs;

    return {p, q, r, s};
}

tuple<string, string> GU::extractCFGandSF(const string &csf)
{
    string cfg, sf;
    for (const char &d : csf)
    {
        if (d == '0')
        {
            cfg.append("0");
        }
        else if (d == '1')
        {
            cfg.append("1");
            sf.append("+");
        }
        else if (d == '2')
        {
            cfg.append("1");
            sf.append("-");
        }
        else if (d == '3')
        {
            cfg.append("2");
        }
        else
        {
            throw std::runtime_error(fmt::format("   Error in reading CSF: %s\
                                                \n      False step-value: %s",
                                                 csf, d));
        }
    }
    return std::make_tuple(cfg, sf);
}

string GU::extractSF(const string &csf)
{
    string sf("");
    for (const char &d : csf)
    {
        if (d == '1')
            sf.append("+");
        else if (d == '2')
            sf.append("-");
    }
    return sf;
}

vector<string> GU::returnSFs(const map<int, string> &sf_map,
                             const vector<int> &sf_idxs)
{
    vector<string> sfs(sf_idxs.size());
    for (size_t i = 0; i < sf_idxs.size(); i++)
        sfs[i] = sf_map.at(sf_idxs.at(i));

    return sfs;
}

map<string, int> GU::returnSFMap(const map<string, int> &sf_map_in, const set<string> &sfs)
{
    map<string, int> sf_map;
    for (auto &sf : sfs)
    {
        size_t sf_idx = sf_map_in.at(sf);
        sf_map[sf] = sf_idx;
    }
    return sf_map;
}

quintet GU::returnCCInfo(const int &p, const int &q, const int &nue_left,
                         const int &nue_right, const string &cfg_left,
                         const string &cfg_right)
{
    ExcType exctype;
    if (cfg_right[p] == '1' and cfg_right[q] == '2')
        exctype = ExcType::DS;
    else if (cfg_right[p] == '0' and cfg_right[q] == '2')
        exctype = ExcType::DV;
    else if (cfg_right[p] == '1' and cfg_right[q] == '1')
        exctype = ExcType::SS;
    else if (cfg_right[p] == '0' and cfg_right[q] == '1')
        exctype = ExcType::SV;
    else
        std::runtime_error("Incompatible orbital occupations, could not deduce exctype!");

    if (exctype == ExcType::DS)
    {
        int prel = 0;
        for (int i = 0; i < p; i++)
            if (cfg_right[i] == '1')
                prel++;
        int qrel = 0;
        for (int i = 0; i < q; i++)
            if (cfg_right[i] == '1')
                qrel++;

        return std::make_tuple(exctype, nue_left, nue_right, prel, qrel);
    }
    else if (exctype == ExcType::DV)
    {
        int prel = 0;
        for (int i = 0; i < p; i++)
            if (cfg_left[i] == '1')
                prel++;
        int qrel = 0;
        for (int i = 0; i < q; i++)
            if (cfg_left[i] == '1')
                qrel++;

        return std::make_tuple(exctype, nue_left, nue_right, prel, qrel);
    }
    else if (exctype == ExcType::SS)
    {
        int prel = 0;
        for (int i = 0; i < p; i++)
            if (cfg_right[i] == '1')
                prel++;
        int qrel = 0;
        for (int i = 0; i < q; i++)
            if (cfg_right[i] == '1')
                qrel++;
        return std::make_tuple(exctype, nue_left, nue_right, prel, qrel);
    }
    else if (exctype == ExcType::SV)
    {
        int prel = 0;
        for (int i = 0; i < p; i++)
            if (cfg_right[i] == '1')
                prel++;
        int qrel = 0;
        for (int i = 0; i < q; i++)
            if (cfg_right[i] == '1')
                qrel++;

        return std::make_tuple(exctype, nue_left, nue_right, prel, qrel);
    }
    else
        throw std::runtime_error("Error: Wrong exctype!");
}

vector<string> GU::findConnectedSFs_DOMOSOMO(const int &p, const int &q, const CFGProto &right)
{
    set<string> connected_sfs;
    if (p < q)
        climbShavittGraphs_DOMOSOMO_R(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);
    else
        climbShavittGraphs_DOMOSOMO_L(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);

    return vector<string>(connected_sfs.begin(), connected_sfs.end());
}

vector<string> GU::findConnectedSFs_DOMOVirtual(const int &p, const int &q,
                                                const CFGProto &right)
{
    set<string> connected_sfs;
    if (p < q)
        climbShavittGraphs_DOMOVirtual_R(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);
    else
        climbShavittGraphs_DOMOVirtual_L(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);

    return vector<string>(connected_sfs.begin(), connected_sfs.end());
}

vector<string> GU::findConnectedSFs_SOMOSOMO(const int &p, const int &q,
                                             const CFGProto &right)
{
    set<string> connected_sfs;
    if (p < q)
        climbShavittGraphs_SOMOSOMO_R(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);
    else
        climbShavittGraphs_SOMOSOMO_L(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);

    return vector<string>(connected_sfs.begin(), connected_sfs.end());
}

vector<string> GU::findConnectedSFs_SOMOVirtual(const int &p, const int &q,
                                                const CFGProto &right)
{
    set<string> connected_sfs;
    if (p < q)
        climbShavittGraphs_SOMOVirtual_R(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);
    else
        climbShavittGraphs_SOMOVirtual_L(p, q, false, 0, 0, "", right.getRoot(), connected_sfs);

    return vector<string>(connected_sfs.begin(), connected_sfs.end());
}

arma::dmat GU::calcCCDOMOSOMO(const int &prel, const int &qrel,
                              const CFGProto &cfg_left,
                              const CFGProto &cfg_right)
{
    int p, q;
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

    arma::dmat ccs(cfg_left.getNumCSFs(), cfg_right.getNumCSFs(),
                   arma::fill::zeros);
    if (p < q)
        walkShavittGraphs_DOMOSOMO_R(p, q, false, 1, 0,
                                     cfg_left.getRoot(),
                                     cfg_right.getRoot(), ccs);
    else
        walkShavittGraphs_DOMOSOMO_L(p, q, false, 1, 0,
                                     cfg_left.getRoot(),
                                     cfg_right.getRoot(), ccs);

    return ccs;
}

arma::dmat GU::calcCCDOMOVirtual(const int &prel, const int &qrel,
                                 const CFGProto &cfg_left,
                                 const CFGProto &cfg_right)
{
    int p, q;
    p = prel;
    q = qrel;

    arma::dmat ccs(cfg_right.getNumCSFs(), cfg_left.getNumCSFs(),
                   arma::fill::zeros);

    if (q < p)
        walkShavittGraphs_SOMOSOMO_R(q, p, false, 1, 0,
                                     cfg_right.getRoot(),
                                     cfg_left.getRoot(), ccs);
    else
        walkShavittGraphs_SOMOSOMO_L(q, p, false, 1, 0,
                                     cfg_right.getRoot(),
                                     cfg_left.getRoot(), ccs);

    return ccs.t();
}

arma::dmat GU::calcCCSOMOSOMO(const int &prel, const int &qrel,
                              const CFGProto &cfg_left,
                              const CFGProto &cfg_right)
{
    int p, q;
    p = prel;
    q = qrel;

    arma::dmat ccs(cfg_left.getNumCSFs(), cfg_right.getNumCSFs(),
                   arma::fill::zeros);
    if (p < q)
        walkShavittGraphs_SOMOSOMO_R(p, q, false, 1, 0,
                                     cfg_left.getRoot(),
                                     cfg_right.getRoot(), ccs);
    else
        walkShavittGraphs_SOMOSOMO_L(p, q, false, 1, 0,
                                     cfg_left.getRoot(),
                                     cfg_right.getRoot(), ccs);

    return ccs;
}

arma::dmat GU::calcCCSOMOVirtual(const int &prel, const int &qrel,
                                 const CFGProto &cfg_left,
                                 const CFGProto &cfg_right)
{
    int p, q;
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

    arma::dmat ccs(cfg_left.getNumCSFs(), cfg_right.getNumCSFs(),
                   arma::fill::zeros);
    if (p < q)
        walkShavittGraphs_SOMOVirtual_R(p, q, false, 1, 0,
                                        cfg_left.getRoot(),
                                        cfg_right.getRoot(), ccs);
    else
        walkShavittGraphs_SOMOVirtual_L(p, q, false, 1, 0,
                                        cfg_left.getRoot(),
                                        cfg_right.getRoot(), ccs);

    return ccs;
}

void GU::climbShavittGraphs_DOMOSOMO_R(const int &p, const int &q, bool flip,
                                       int b_left, int i, string sf,
                                       CSFTree::Node *node_right,
                                       set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left, i + 1, sf,
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_DOMOSOMO_R(p, q, not(flip), b_left, i + 1, sf,
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[3] != nullptr)
        {
            if (b_left + 1 == node_right->b)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[3].get(), connected_sfs);

            if (b_left - 1 == node_right->b)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[3].get(), connected_sfs);
        }
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_R(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_DOMOSOMO_L(const int &p, const int &q, bool flip,
                                       int b_left, int i, string sf,
                                       CSFTree::Node *node_right,
                                       set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[3] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[3].get(), connected_sfs);

        if (node_right->step_values[3] != nullptr)
            climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[3].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[1] != nullptr)
            if (b_left == node_right->step_values[1]->b)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            if (b_left == node_right->step_values[2]->b)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[2].get(), connected_sfs);
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_DOMOSOMO_L(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_DOMOVirtual_R(const int &p, const int &q, bool flip,
                                          int b_right, int i, string sf,
                                          CSFTree::Node *node_left,
                                          set<string> &connected_sfs)
{
    if (b_right < 0)
        return;

    if (node_left->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < p)
    {
        if (node_left->step_values[1] != nullptr)
            climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[1].get(), connected_sfs);

        if (node_left->step_values[2] != nullptr)
            climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[2].get(), connected_sfs);
    }
    else if (i > q)
    {
        if (node_left->step_values[1] != nullptr)
            climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[1].get(), connected_sfs);

        if (node_left->step_values[2] != nullptr)
            climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr)
        {
            climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[3].get(), connected_sfs);

            climbShavittGraphs_DOMOVirtual_R(p, q, not(flip), b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[3].get(), connected_sfs);
        }
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr)
        {
            if (b_right + 1 == node_left->b)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[0].get(), connected_sfs);

            if (b_right - 1 == node_left->b)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[0].get(), connected_sfs);
        }
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[1].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[2].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, not(flip), b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[2].get(), connected_sfs);
        }
        else
        {
            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[1].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[2].get(), connected_sfs);

            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_R(p, q, not(flip), b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[1].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_DOMOVirtual_L(const int &p, const int &q, bool flip,
                                          int b_right, int i, string sf,
                                          CSFTree::Node *node_left,
                                          set<string> &connected_sfs)
{
    if (b_right < 0)
        return;

    if (node_left->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < q)
    {
        if (node_left->step_values[1] != nullptr)
            climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[1].get(), connected_sfs);

        if (node_left->step_values[2] != nullptr)
            climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[2].get(), connected_sfs);
    }
    else if (i > p)
    {
        if (node_left->step_values[1] != nullptr)
            climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[1].get(), connected_sfs);

        if (node_left->step_values[2] != nullptr)
            climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr)
        {
            climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                             node_left->step_values[0].get(), connected_sfs);

            climbShavittGraphs_DOMOVirtual_L(p, q, not(flip), b_right - 1, i + 1, sf + "-",
                                             node_left->step_values[0].get(), connected_sfs);
        }
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr)
        {
            if (b_right + 1 == node_left->b)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[3].get(), connected_sfs);

            if (b_right - 1 == node_left->b)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[3].get(), connected_sfs);
        }
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[1].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[2].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, not(flip), b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[2].get(), connected_sfs);
        }
        else
        {
            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right + 1, i + 1, sf + "+",
                                                 node_left->step_values[1].get(), connected_sfs);

            if (node_left->step_values[2] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, flip, b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[2].get(), connected_sfs);

            if (node_left->step_values[1] != nullptr)
                climbShavittGraphs_DOMOVirtual_L(p, q, not(flip), b_right - 1, i + 1, sf + "-",
                                                 node_left->step_values[1].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_SOMOSOMO_R(const int &p, const int &q, bool flip,
                                       int b_left, int i, string sf,
                                       CSFTree::Node *node_right,
                                       set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left, i + 1, sf,
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_R(p, q, not(flip), b_left, i + 1, sf,
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[1] != nullptr)
            if (b_left == node_right->step_values[1]->b)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            if (b_left == node_right->step_values[2]->b)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[2].get(), connected_sfs);
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_R(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_SOMOSOMO_L(const int &p, const int &q, bool flip,
                                       int b_left, int i, string sf,
                                       CSFTree::Node *node_right,
                                       set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left, i + 1, sf,
                                          node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOSOMO_L(p, q, not(flip), b_left, i + 1, sf,
                                          node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[1] != nullptr)
            if (b_left == node_right->step_values[1]->b)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            if (b_left == node_right->step_values[2]->b)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left, i + 1, sf,
                                              node_right->step_values[2].get(), connected_sfs);
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                              node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOSOMO_L(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                              node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_SOMOVirtual_R(const int &p, const int &q, bool flip,
                                          int b_left, int i, string sf,
                                          CSFTree::Node *node_right,
                                          set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                             node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                             node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                             node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                             node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[0] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                             node_right->step_values[0].get(), connected_sfs);

        if (node_right->step_values[0] != nullptr)
            climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                             node_right->step_values[0].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[1] != nullptr)
            if (b_left == node_right->step_values[1]->b)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left, i + 1, sf,
                                                 node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            if (b_left == node_right->step_values[2]->b)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left, i + 1, sf,
                                                 node_right->step_values[2].get(), connected_sfs);
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, flip, b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_R(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::climbShavittGraphs_SOMOVirtual_L(const int &p, const int &q, bool flip,
                                          int b_left, int i, string sf,
                                          CSFTree::Node *node_right,
                                          set<string> &connected_sfs)
{
    if (b_left < 0)
        return;

    if (node_right->end)
    {
        connected_sfs.insert(sf);
        return;
    }

    if (i < q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                             node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                             node_right->step_values[2].get(), connected_sfs);
    }
    else if (i > p)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                             node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                             node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == q)
    {
        if (node_right->step_values[1] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left, i + 1, sf,
                                             node_right->step_values[1].get(), connected_sfs);

        if (node_right->step_values[2] != nullptr)
            climbShavittGraphs_SOMOVirtual_L(p, q, not(flip), b_left, i + 1, sf,
                                             node_right->step_values[2].get(), connected_sfs);
    }
    else if (i == p)
    {
        if (node_right->step_values[0] != nullptr)
        {
            if (b_left + 1 == node_right->b)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[0].get(), connected_sfs);

            if (b_left - 1 == node_right->b)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[0].get(), connected_sfs);
        }
    }
    else
    {
        if (flip)
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, not(flip), b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[1].get(), connected_sfs);
        }
        else
        {
            if (node_right->step_values[1] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[1].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, flip, b_left - 1, i + 1, sf + "-",
                                                 node_right->step_values[2].get(), connected_sfs);

            if (node_right->step_values[2] != nullptr)
                climbShavittGraphs_SOMOVirtual_L(p, q, not(flip), b_left + 1, i + 1, sf + "+",
                                                 node_right->step_values[2].get(), connected_sfs);
        }
    }
}

void GU::walkShavittGraphs_DOMOSOMO_R(const int &p, const int &q,
                                      bool flip, double cc, int i,
                                      CSFTree::Node *node_left,
                                      CSFTree::Node *node_right,
                                      arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * A(node_right->step_values[1]->b, 1, 0), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[3] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, not(flip), cc * A(node_right->step_values[2]->b, 1, 2), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[3] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * A(node_right->step_values[3]->b, 0, 1), i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[3].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[3] != nullptr)
            walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * A(node_right->step_values[3]->b, 2, 1), i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[3].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * C(node_right->step_values[2]->b, 2), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, not(flip), cc * (1.0 / (node_right->step_values[1]->b)), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * C(node_right->step_values[1]->b, 0), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_R(p, q, not(flip), cc * (-1.0 / (node_right->step_values[2]->b + 2)), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[2].get(), cc_matrix);
        }
    }
}

void GU::walkShavittGraphs_DOMOSOMO_L(const int &p, const int &q,
                                      bool flip, double cc, int i,
                                      CSFTree::Node *node_left,
                                      CSFTree::Node *node_right,
                                      arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[3] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, not(flip), cc * A(node_right->step_values[3]->b, 2, 1), i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[3].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[3] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * A(node_right->step_values[3]->b, 0, 1), i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[3].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * A(node_right->step_values[1]->b, 0, 1), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[3] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * A(node_right->step_values[2]->b, 2, 1), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * C(node_right->step_values[1]->b, 1), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, not(flip), cc * (-1.0 / (node_right->step_values[1]->b + 1)), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, flip, cc * C(node_right->step_values[2]->b, 1), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_DOMOSOMO_L(p, q, not(flip), cc * (1.0 / (node_right->step_values[2]->b + 1)), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[2].get(), cc_matrix);
        }
    }
}

void GU::walkShavittGraphs_SOMOSOMO_R(const int &p, const int &q,
                                      bool flip, double cc, int i,
                                      CSFTree::Node *node_left,
                                      CSFTree::Node *node_right,
                                      arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc * A(node_right->step_values[1]->b, 1, 0), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[3] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, not(flip), cc * A(node_right->step_values[2]->b, 1, 2), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[0].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[0] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc, i + 1,
                                         node_left->step_values[0].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc * C(node_right->step_values[2]->b, 2), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, not(flip), cc * (1.0 / (node_right->step_values[1]->b)), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc * C(node_right->step_values[1]->b, 0), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_R(p, q, not(flip), cc * (-1.0 / (node_right->step_values[2]->b + 2)), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[2].get(), cc_matrix);
        }
    }
}

void GU::walkShavittGraphs_SOMOSOMO_L(const int &p, const int &q,
                                      bool flip, double cc, int i,
                                      CSFTree::Node *node_left,
                                      CSFTree::Node *node_right,
                                      arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[1].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[2].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc, i + 1,
                                         node_left->step_values[0].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[0] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, not(flip), cc, i + 1,
                                         node_left->step_values[0].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[3] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * A(node_right->step_values[1]->b, 0, 1), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[3] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * A(node_right->step_values[2]->b, 2, 1), i + 1,
                                         node_left->step_values[3].get(),
                                         node_right->step_values[2].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * C(node_right->step_values[1]->b, 1), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, not(flip), cc * (-1.0 / (node_right->step_values[1]->b + 1)), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * -1, i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, flip, cc * C(node_right->step_values[2]->b, 1), i + 1,
                                             node_left->step_values[2].get(),
                                             node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOSOMO_L(p, q, not(flip), cc * (1.0 / (node_right->step_values[2]->b + 1)), i + 1,
                                             node_left->step_values[1].get(),
                                             node_right->step_values[2].get(), cc_matrix);
        }
    }
}

void GU::walkShavittGraphs_SOMOVirtual_R(const int &p, const int &q,
                                         bool flip, double cc, int i,
                                         CSFTree::Node *node_left,
                                         CSFTree::Node *node_right,
                                         arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[0] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, not(flip), cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[0].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[0] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[0].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[0].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[0] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc, i + 1,
                                            node_left->step_values[0].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc * -1, i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc * C(node_right->step_values[2]->b, 2), i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, not(flip), cc * (1.0 / (node_right->step_values[1]->b)), i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc * C(node_right->step_values[1]->b, 0), i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, flip, cc * -1, i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_R(p, q, not(flip), cc * (-1.0 / (node_right->step_values[2]->b + 2)), i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[2].get(), cc_matrix);
        }
    }
}

void GU::walkShavittGraphs_SOMOVirtual_L(const int &p, const int &q,
                                         bool flip, double cc, int i,
                                         CSFTree::Node *node_left,
                                         CSFTree::Node *node_right,
                                         arma::dmat &cc_matrix)
{
    if (node_left->end and node_right->end)
    {
        cc_matrix(node_left->pos, node_right->pos) = cc;
        return;
    }

    if (i < q)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else if (i > p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == q)
    {
        if (node_left->step_values[0] != nullptr and node_right->step_values[1] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[0].get(),
                                            node_right->step_values[1].get(), cc_matrix);

        if (node_left->step_values[0] != nullptr and node_right->step_values[2] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, not(flip), cc, i + 1,
                                            node_left->step_values[0].get(),
                                            node_right->step_values[2].get(), cc_matrix);
    }
    else if (i == p)
    {
        if (node_left->step_values[1] != nullptr and node_right->step_values[0] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[1].get(),
                                            node_right->step_values[0].get(), cc_matrix);

        if (node_left->step_values[2] != nullptr and node_right->step_values[0] != nullptr)
            walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc, i + 1,
                                            node_left->step_values[2].get(),
                                            node_right->step_values[0].get(), cc_matrix);
    }
    else
    {
        if (flip)
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc * C(node_right->step_values[1]->b, 1), i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc * -1, i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, not(flip), cc * (-1.0 / (node_right->step_values[1]->b + 1)), i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[1].get(), cc_matrix);
        }
        else
        {
            if (node_left->step_values[1] != nullptr and node_right->step_values[1] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc * -1, i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[1].get(), cc_matrix);

            if (node_left->step_values[2] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, flip, cc * C(node_right->step_values[2]->b, 1), i + 1,
                                                node_left->step_values[2].get(),
                                                node_right->step_values[2].get(), cc_matrix);

            if (node_left->step_values[1] != nullptr and node_right->step_values[2] != nullptr)
                walkShavittGraphs_SOMOVirtual_L(p, q, not(flip), cc * (1.0 / (node_right->step_values[2]->b + 1)), i + 1,
                                                node_left->step_values[1].get(),
                                                node_right->step_values[2].get(), cc_matrix);
        }
    }
}

template <typename T>
void GU::mergeSFPairs(T &sf_pairs_in, T &sf_pairs_out)
{
    for (auto &item : sf_pairs_in)
    {
        auto key = item.first;
        if (sf_pairs_out.find(key) == sf_pairs_out.end())
            sf_pairs_out.emplace(key, move(item.second));
        else
        {
            sf_pairs_out[key].first.merge(sf_pairs_in[key].first);
            sf_pairs_out[key].second.merge(sf_pairs_in[key].second);
        }
    }
}
template void GU::mergeSFPairs(sf_pair_map_1el &, sf_pair_map_1el &);
template void GU::mergeSFPairs(sf_pair_map_2el &, sf_pair_map_2el &);

template <typename T>
void GU::mergeConnections(const T &connections_in, T &connections_out)
{
    for (auto &item : connections_in)
    {
        auto key = item.first;
        // TODO: Take the reference of a map here
        if (connections_out.find(key) == connections_out.end())
            connections_out.emplace(key, move(item.second));
        else
            connections_out[key].insert(connections_out[key].end(),
                                        make_move_iterator(item.second.begin()), make_move_iterator(item.second.end()));
    }
}
template void GU::mergeConnections(const connection_map_1el &, connection_map_1el &);
template void GU::mergeConnections(const connection_map_2el &, connection_map_2el &);
template void GU::mergeConnections(const connection_map_dia &, connection_map_dia &);

template <typename T>
void GU::mergeCCs(const T &ccs_in, T &ccs_out)
{
    for (auto &item : ccs_in)
    {
        auto key = item.first;
        if (ccs_out.find(key) == ccs_out.end())
        {
            ccs_out.emplace(key, std::move(item.second));
        }
        else
        {
            auto map_inner = ccs_out.at(key);
            // for (auto &item1 : item.second)
            // {
            //     size_t idx_left = item1.first;
            //     if (map_inner.find(idx_left) == map_inner.end())
            //     {
            //         map_inner.emplace(idx_left, move(item1.second));
            //     }
            //     else
            //     {
            //         for (auto &item2 : item1.second)
            //         {
            //             size_t idx_right = item2.first;
            //             map_inner[idx_left][idx_right] = item2.second;
            //         }
            //     }
            // }

            for (auto &item1 : item.second)
            {
                size_t idx_left = item1.first;
                for (auto &item2 : item1.second)
                {
                    size_t idx_right = item2.first;
                    map_inner[idx_left][idx_right] = item2.second;
                }
            }

            ccs_out[key] = move(map_inner);
            // ccs_out.emplace(key, move(map_inner));
        }
    }
}
template void GU::mergeCCs(const map<nonet, cc_map> &, map<nonet, cc_map> &);
template void GU::mergeCCs(const map<quintet, cc_map> &, map<quintet, cc_map> &);