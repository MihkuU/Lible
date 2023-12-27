#include "coupling_coeffs.h"

#include <armadillo>
#include <omp.h>

using namespace Lible::GUGA;
using namespace Lible::GUGA::Util;

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

void GCI::CouplingCoeffs::constructCCs(const sf_pair_map_1el &sf_pairs_1el_new,
                                       const sf_pair_map_2el &sf_pairs_2el_new,
                                       std::map<nonet, cc_map> &ccs_2el,
                                       std::map<quintet, cc_map> &ccs_1el)
{
    sf_pair_map_1el sf_pairs_1el_tmp = sf_pairs_1el;
    sf_pair_map_2el sf_pairs_2el_tmp = sf_pairs_2el;

#pragma omp parallel
    {
        map<quintet, cc_map> ccs_1el_omp;
        map<nonet, cc_map> ccs_2el_omp;
        sf_pair_map_1el sf_pairs_1el_omp;
        sf_pair_map_2el sf_pairs_2el_omp;

        int thread_num = omp_get_thread_num();
        int ipal = 0;
        for (auto &[key_new, sfs_pair] : sf_pairs_1el_new)
        {
            if (ipal % omp_get_num_threads() != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            sfs_pair_t sfs_pair_new = returnNewSFPair(sfs_pair, sf_pairs_1el_tmp, key_new);
            if (sf_pairs_1el.find(key_new) == sf_pairs_1el.end())
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes1El(sfs_pair_new, key_new);
                calcCCs1El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_1el_omp);
            }
            else
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes1El(sfs_pair_new, key_new);
                calcCCs1El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_1el_omp);

                sfs_pair_t sfs_pair = sf_pairs_1el.at(key_new);
                sfs_pair_t sfs_pair_leftnew_right = returnSFPair(sfs_pair_new.first, sfs_pair.second);
                proto_1el_tuple cfg_prototypes_leftnew_right = returnCFGPrototypes1El(sfs_pair_leftnew_right, key_new);
                calcCCs1El(sfs_pair_leftnew_right, cfg_prototypes_leftnew_right, key_new, ccs_1el_omp);

                sfs_pair_t sfs_pair_left_rightnew = returnSFPair(sfs_pair.first, sfs_pair_new.second);
                proto_1el_tuple cfg_prototypes_left_rightnew = returnCFGPrototypes1El(sfs_pair_left_rightnew, key_new);
                calcCCs1El(sfs_pair_left_rightnew, cfg_prototypes_left_rightnew, key_new, ccs_1el_omp);
            }
            sf_pairs_1el_omp[key_new] = sfs_pair_new;
        }

        ipal = 0;
        for (auto &[key_new, sfs_pair] : sf_pairs_2el_new)
        {
            if (ipal % omp_get_num_threads() != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            sfs_pair_t sfs_pair_new = returnNewSFPair(sfs_pair, sf_pairs_2el_tmp, key_new);
            if (sf_pairs_2el.find(key_new) == sf_pairs_2el.end())
            {
                proto_2el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes2El(sfs_pair_new, key_new);
                calcCCs2El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_2el_omp);
            }
            else
            {
                proto_2el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes2El(sfs_pair_new, key_new);
                calcCCs2El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_2el_omp);

                sfs_pair_t sfs_pair = sf_pairs_2el.at(key_new);
                sfs_pair_t sfs_pair_leftnew_right = returnSFPair(sfs_pair_new.first, sfs_pair.second);
                proto_2el_tuple cfg_prototypes_leftnew_right = returnCFGPrototypes2El(sfs_pair_leftnew_right, key_new);
                calcCCs2El(sfs_pair_leftnew_right, cfg_prototypes_leftnew_right, key_new, ccs_2el_omp);

                sfs_pair_t sfs_pair_left_rightnew = returnSFPair(sfs_pair.first, sfs_pair_new.second);
                proto_2el_tuple cfg_prototypes_left_rightnew = returnCFGPrototypes2El(sfs_pair_left_rightnew, key_new);
                calcCCs2El(sfs_pair_left_rightnew, cfg_prototypes_left_rightnew, key_new, ccs_2el_omp);
            }
            sf_pairs_2el_omp[key_new] = sfs_pair_new;
        }

#pragma omp critical
        {
            mergeCCs(ccs_1el_omp, ccs_1el);
            mergeCCs(ccs_2el_omp, ccs_2el);
            mergeSFPairs(sf_pairs_1el_omp, sf_pairs_1el);
            mergeSFPairs(sf_pairs_2el_omp, sf_pairs_2el);
        }
    }
}

void GCI::CouplingCoeffs::constructCCs(const connection_map_1el &connections_1el_new,
                                       const connection_map_2el &connections_2el_new,
                                       const connection_map_dia &connections_dia_new,
                                       const wfn_ptr &wave_function,
                                       map<nonet, cc_map> &ccs_2el,
                                       map<quintet, cc_map> &ccs_1el,
                                       map<quintet, cc_map> &ccs_dia)
{
    sf_pair_map_1el sf_pairs_1el_tmp = sf_pairs_1el;
    sf_pair_map_2el sf_pairs_2el_tmp = sf_pairs_2el;
    sf_pair_map_1el sf_pairs_dia_tmp = sf_pairs_dia;

#pragma omp parallel
    {
        map<quintet, cc_map> ccs_1el_omp;
        sf_pair_map_1el sf_pairs_1el_omp;

        int thread_num = omp_get_thread_num();
        int ipal = 0;
        for (auto &[key_new, connections] : connections_1el_new)
        {
            if (ipal % omp_get_num_threads() != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            sfs_pair_t sfs_pair_new = returnNewSFPair(connections, sf_pairs_1el_tmp, key_new, wave_function);
            if (sf_pairs_1el_tmp.find(key_new) == sf_pairs_1el_tmp.end())
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes1El(sfs_pair_new, key_new);
                calcCCs1El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_1el_omp);
            }
            else
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes1El(sfs_pair_new, key_new);
                calcCCs1El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_1el_omp);

                sfs_pair_t sfs_pair = sf_pairs_1el_tmp[key_new];
                sfs_pair_t sfs_pair_leftnew_right = returnSFPair(sfs_pair_new.first, sfs_pair.second);
                proto_1el_tuple cfg_prototypes_leftnew_right = returnCFGPrototypes1El(sfs_pair_leftnew_right, key_new);
                calcCCs1El(sfs_pair_leftnew_right, cfg_prototypes_leftnew_right, key_new, ccs_1el_omp);

                sfs_pair_t sfs_pair_left_rightnew = returnSFPair(sfs_pair.first, sfs_pair_new.second);
                proto_1el_tuple cfg_prototypes_left_rightnew = returnCFGPrototypes1El(sfs_pair_left_rightnew, key_new);
                calcCCs1El(sfs_pair_left_rightnew, cfg_prototypes_left_rightnew, key_new, ccs_1el_omp);
            }
            sf_pairs_1el_omp[key_new] = sfs_pair_new;
        }
#pragma omp critical
        {
            mergeCCs(ccs_1el_omp, ccs_1el);
            mergeSFPairs(sf_pairs_1el_omp, sf_pairs_1el);
        }
    }

#pragma omp parallel
    {
        map<nonet, cc_map> ccs_2el_omp;
        sf_pair_map_2el sf_pairs_2el_omp;

        int thread_num = omp_get_thread_num();
        int ipal = 0;
        for (auto &[key_new, connections] : connections_2el_new)
        {
            if (ipal % omp_get_num_threads() != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            sfs_pair_t sfs_pair_new = returnNewSFPair(connections, sf_pairs_2el_tmp, key_new, wave_function);
            if (sf_pairs_2el_tmp.find(key_new) == sf_pairs_2el_tmp.end())
            {
                proto_2el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes2El(sfs_pair_new, key_new);
                calcCCs2El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_2el_omp);
            }
            else
            {
                proto_2el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypes2El(sfs_pair_new, key_new);
                calcCCs2El(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_2el_omp);

                sfs_pair_t sfs_pair = sf_pairs_2el_tmp[key_new];
                sfs_pair_t sfs_pair_leftnew_right = returnSFPair(sfs_pair_new.first, sfs_pair.second);
                proto_2el_tuple cfg_prototypes_leftnew_right = returnCFGPrototypes2El(sfs_pair_leftnew_right, key_new);
                calcCCs2El(sfs_pair_leftnew_right, cfg_prototypes_leftnew_right, key_new, ccs_2el_omp);

                sfs_pair_t sfs_pair_left_rightnew = returnSFPair(sfs_pair.first, sfs_pair_new.second);
                proto_2el_tuple cfg_prototypes_left_rightnew = returnCFGPrototypes2El(sfs_pair_left_rightnew, key_new);
                calcCCs2El(sfs_pair_left_rightnew, cfg_prototypes_left_rightnew, key_new, ccs_2el_omp);
            }
            sf_pairs_2el_omp[key_new] = sfs_pair_new;
        }
#pragma omp critical
        {
            mergeCCs(ccs_2el_omp, ccs_2el);
            mergeSFPairs(sf_pairs_2el_omp, sf_pairs_2el);
        }
    }

#pragma omp parallel
    {
        map<quintet, cc_map> ccs_dia_omp;
        sf_pair_map_1el sf_pairs_dia_omp;

        int thread_num = omp_get_thread_num();
        int ipal = 0;
        for (auto &[key_new, connections] : connections_dia_new)
        {
            if (ipal % omp_get_num_threads() != thread_num)
            {
                ipal++;
                continue;
            }
            ipal++;

            sfs_pair_t sfs_pair_new = returnNewSFPair(connections, sf_pairs_dia_tmp, key_new, wave_function);
            if (sf_pairs_dia_tmp.find(key_new) == sf_pairs_dia_tmp.end())
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypesDia(sfs_pair_new, key_new);
                calcCCsDia(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_dia_omp);
            }
            else
            {
                proto_1el_tuple cfg_prototypes_leftnew_rightnew = returnCFGPrototypesDia(sfs_pair_new, key_new);
                calcCCsDia(sfs_pair_new, cfg_prototypes_leftnew_rightnew, key_new, ccs_dia_omp);

                sfs_pair_t sfs_pair = sf_pairs_dia_tmp[key_new];
                sfs_pair_t sfs_pair_leftnew_right = returnSFPair(sfs_pair_new.first, sfs_pair.second);
                proto_1el_tuple cfg_prototypes_leftnew_right = returnCFGPrototypesDia(sfs_pair_leftnew_right, key_new);
                calcCCsDia(sfs_pair_leftnew_right, cfg_prototypes_leftnew_right, key_new, ccs_dia_omp);

                sfs_pair_t sfs_pair_left_rightnew = returnSFPair(sfs_pair.first, sfs_pair_new.second);
                proto_1el_tuple cfg_prototypes_left_rightnew = returnCFGPrototypesDia(sfs_pair_left_rightnew, key_new);
                calcCCsDia(sfs_pair_left_rightnew, cfg_prototypes_left_rightnew, key_new, ccs_dia_omp);
            }
        }
#pragma omp critical
        {
            mergeCCs(ccs_dia_omp, ccs_dia);
            mergeSFPairs(sf_pairs_dia_omp, sf_pairs_dia);
        }
    }
}

sfs_pair_t GCI::CouplingCoeffs::returnNewSFPair(const sfs_pair_t &sfs_pair_trial,
                                                const sf_pair_map_1el &sf_pairs_1el,
                                                const quintet &key)
{
    if (sf_pairs_1el.find(key) == sf_pairs_1el.end())
        return sfs_pair_trial;

    sfs_pair_t sfs_pair = sf_pairs_1el.at(key);
    sfs_pair_t sfs_pair_new;
    for (auto &sf_idx : sfs_pair_trial.first)
        if (sfs_pair.first.find(sf_idx) == sfs_pair.first.end())
            sfs_pair_new.first.insert(sf_idx);

    for (auto &sf_idx : sfs_pair_trial.second)
        if (sfs_pair.second.find(sf_idx) == sfs_pair.second.end())
            sfs_pair_new.second.insert(sf_idx);

    return sfs_pair_new;
}

sfs_pair_t GCI::CouplingCoeffs::returnNewSFPair(const sfs_pair_t &sfs_pair_trial,
                                                const sf_pair_map_2el &sf_pairs_2el,
                                                const nonet &key)
{
    if (sf_pairs_2el.find(key) == sf_pairs_2el.end())
        return sfs_pair_trial;

    sfs_pair_t sfs_pair = sf_pairs_2el.at(key);
    sfs_pair_t sfs_pair_new;
    for (auto &sf_idx : sfs_pair_trial.first)
        if (sfs_pair.first.find(sf_idx) == sfs_pair.first.end())
            sfs_pair_new.first.insert(sf_idx);

    for (auto &sf_idx : sfs_pair_trial.second)
        if (sfs_pair.second.find(sf_idx) == sfs_pair.second.end())
            sfs_pair_new.second.insert(sf_idx);

    return sfs_pair_new;
}

sfs_pair_t GCI::CouplingCoeffs::returnNewSFPair(const connection_list_1el &connections,
                                                const sf_pair_map_1el &sf_pairs_1el,
                                                const quintet &key,
                                                const wfn_ptr &wave_function)
{
    sfs_pair_t sf_pair_;
    for (auto &connection : connections)
    {
        int icfg_left = get<0>(connection);
        int icfg_right = get<1>(connection);
        const CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
        const CFG *cfg_p_right = wave_function->getCFGPtr(icfg_right);
        vector<int> sfs_left = cfg_left->getSFIdxs();
        vector<int> sfs_right = cfg_p_right->getSFIdxs();
        sf_pair_.first.insert(sfs_left.begin(), sfs_left.end());
        sf_pair_.second.insert(sfs_right.begin(), sfs_right.end());
    }

    sfs_pair_t sfs_pair_new;
    if (sf_pairs_1el.find(key) == sf_pairs_1el.end())
        return sf_pair_;
    else
    {
        sfs_pair_t sfs_pair = sf_pairs_1el.at(key);
        for (auto &sf_idx : sf_pair_.first)
            if (sfs_pair.first.find(sf_idx) == sfs_pair.first.end())
                sfs_pair_new.first.insert(sf_idx);

        for (auto &sf_idx : sf_pair_.second)
            if (sfs_pair.second.find(sf_idx) == sfs_pair.second.end())
                sfs_pair_new.second.insert(sf_idx);
    }

    return sfs_pair_new;
}

sfs_pair_t GCI::CouplingCoeffs::returnNewSFPair(const connection_list_2el &connections,
                                                const sf_pair_map_2el &sf_pairs_2el,
                                                const nonet &key,
                                                const wfn_ptr &wave_function)
{
    sfs_pair_t sf_pair_;
    for (auto &connection : connections)
    {
        int icfg_left = get<0>(connection);
        int icfg_right = get<1>(connection);
        const CFG *cfg_left = wave_function->getCFGPtr(icfg_left);
        const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
        vector<int> sfs_left = cfg_left->getSFIdxs();
        vector<int> sfs_right = cfg_right->getSFIdxs();
        sf_pair_.first.insert(sfs_left.begin(), sfs_left.end());
        sf_pair_.second.insert(sfs_right.begin(), sfs_right.end());
    }

    sfs_pair_t sfs_pair_new;
    if (sf_pairs_2el.find(key) == sf_pairs_2el.end())
        return sf_pair_;
    else
    {
        sfs_pair_t sfs_pair = sf_pairs_2el.at(key);
        for (auto &sf_idx : sf_pair_.first)
            if (sfs_pair.first.find(sf_idx) == sfs_pair.first.end())
                sfs_pair_new.first.insert(sf_idx);

        for (auto &sf_idx : sf_pair_.second)
            if (sfs_pair.second.find(sf_idx) == sfs_pair.second.end())
                sfs_pair_new.second.insert(sf_idx);
    }

    return sfs_pair_new;
}

sfs_pair_t GCI::CouplingCoeffs::returnNewSFPair(const connection_list_dia &connections,
                                                const sf_pair_map_1el &sf_pairs_dia,
                                                const quintet &key,
                                                const wfn_ptr &wave_function)
{
    sfs_pair_t sf_pair_;
    for (auto &connection : connections)
    {
        int icfg_right = get<0>(connection);
        const CFG *cfg_right = wave_function->getCFGPtr(icfg_right);
        vector<int> sfs_right = cfg_right->getSFIdxs();
        sf_pair_.first.insert(sfs_right.begin(), sfs_right.end());
        sf_pair_.second.insert(sfs_right.begin(), sfs_right.end());
    }

    sfs_pair_t sfs_pair_new;
    if (sf_pairs_dia.find(key) == sf_pairs_dia.end())
    {
        return sf_pair_;
    }
    else
    {
        sfs_pair_t sfs_pair = sf_pairs_dia.at(key);
        for (auto &sf_idx : sf_pair_.first)
            if (sfs_pair.first.find(sf_idx) == sfs_pair.first.end())
                sfs_pair_new.first.insert(sf_idx);

        for (auto &sf_idx : sf_pair_.second)
            if (sfs_pair.second.find(sf_idx) == sfs_pair.second.end())
                sfs_pair_new.second.insert(sf_idx);
    }

    return sfs_pair_new;
}

sfs_pair_t GCI::CouplingCoeffs::returnSFPair(const set<int> &sfs_left,
                                             const set<int> &sfs_right)
{
    sfs_pair_t sfs_pair;
    sfs_pair.first.insert(sfs_left.begin(), sfs_left.end());
    sfs_pair.second.insert(sfs_right.begin(), sfs_right.end());
    return sfs_pair;
}

proto_1el_tuple GCI::CouplingCoeffs::returnCFGPrototypes1El(const sfs_pair_t &sfs_pair,
                                                            const quintet &key)
{
    proto_1el_tuple cfg_prototypes;

    size_t exctype = get<0>(key);
    size_t nue_left = get<1>(key);
    size_t nue_right = get<2>(key);
    size_t prel = get<3>(key);
    size_t qrel = get<4>(key);

    vector<string> sfs_left;
    for (const int &sf_idx : sfs_pair.first)
        sfs_left.push_back(gci->sfs_map__idx_to_sf.at(nue_left).at(sf_idx));

    vector<string> sfs_right;
    for (const int &sf_idx : sfs_pair.second)
        sfs_right.push_back(gci->sfs_map__idx_to_sf.at(nue_right).at(sf_idx));

    int norb, p, q;
    if (exctype == ExcType::DS)
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

        string left(norb, '1');
        string right(norb, '1');
        left[p] = '2';
        right[q] = '2';

        CFGProto cfg_left(spin, left, sfs_left);
        CFGProto cfg_right(spin, right, sfs_right);
        return std::make_tuple(std::move(cfg_left), std::move(cfg_right));
    }
    else if (exctype == ExcType::DV)
    {
        p = prel;
        q = qrel;

        norb = nue_left;

        string left(norb, '1');
        string right(norb, '1');
        right[p] = '0';
        right[q] = '2';

        CFGProto cfg_left(spin, left, sfs_left);
        CFGProto cfg_right(spin, right, sfs_right);
        return std::make_tuple(std::move(cfg_left), std::move(cfg_right));
    }
    else if (exctype == ExcType::SS)
    {
        p = prel;
        q = qrel;

        int norb = nue_right;

        string left(norb, '1');
        string right(norb, '1');
        left[p] = '2';
        left[q] = '0';

        CFGProto cfg_left(spin, left, sfs_left);
        CFGProto cfg_right(spin, right, sfs_right);
        return std::make_tuple(std::move(cfg_left), std::move(cfg_right));
    }
    else if (exctype == ExcType::SV)
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

        string left(norb, '1');
        left[q] = '0';
        string right(norb, '1');
        right[p] = '0';

        CFGProto cfg_left(spin, left, sfs_left);
        CFGProto cfg_right(spin, right, sfs_right);
        return std::make_tuple(std::move(cfg_left), std::move(cfg_right));
    }
    else
        throw std::runtime_error("\nError: Wrong  exctype!");
}

proto_2el_tuple GCI::CouplingCoeffs::returnCFGPrototypes2El(const sfs_pair_t &sfs_pair,
                                                            const nonet &key)
{
    proto_2el_tuple cfg_prototypes;

    vector<string> connected_sfs_right;
    vector<string> connected_sfs_left;

    CFGProto cfg_left;
    CFGProto cfg_ri_left;
    CFGProto cfg_right;
    CFGProto cfg_ri_right;

    size_t exctype1 = get<0>(key);
    size_t exctype2 = get<1>(key);
    size_t nue_left = get<2>(key);
    size_t nue_middle = get<3>(key);
    size_t nue_right = get<4>(key);
    size_t prel1 = get<5>(key);
    size_t qrel1 = get<6>(key);
    size_t prel2 = get<7>(key);
    size_t qrel2 = get<8>(key);

    vector<string> sfs_left;
    vector<string> sfs_right;
    for (const size_t &sf_idx : sfs_pair.first)
        sfs_left.push_back(gci->sfs_map__idx_to_sf.at(nue_left).at(sf_idx));
    for (const size_t &sf_idx : sfs_pair.second)
        sfs_right.push_back(gci->sfs_map__idx_to_sf.at(nue_right).at(sf_idx));

    int norb, p, q;
    // <RI|Ers|right>
    if (exctype2 == ExcType::DS)
    {
        if (prel2 >= qrel2)
        {
            p = prel2 + 1;
            q = qrel2;
        }
        else
        {
            p = prel2;
            q = qrel2;
        }
        norb = nue_right + 1;

        string middle(norb, '1');
        string right(norb, '1');
        middle[p] = '2';
        right[q] = '2';

        cfg_right = CFGProto(spin, right, sfs_right);
        cfg_ri_right = CFGProto(spin, middle);
        connected_sfs_right = findConnectedSFs_DOMOSOMO(p, q, cfg_right);
    }
    else if (exctype2 == ExcType::DV)
    {
        p = prel2;
        q = qrel2;

        norb = nue_middle;

        string middle(norb, '1');
        string right(norb, '1');
        right[p] = '0';
        right[q] = '2';

        cfg_right = CFGProto(spin, right, sfs_right);
        cfg_ri_right = CFGProto(spin, middle);
        connected_sfs_right = findConnectedSFs_DOMOVirtual(q, p, cfg_right);
    }
    else if (exctype2 == ExcType::SS)
    {
        p = prel2;
        q = qrel2;

        norb = nue_right;

        string middle(norb, '1');
        string right(norb, '1');
        middle[p] = '2';
        middle[q] = '0';

        cfg_right = CFGProto(spin, right, sfs_right);
        cfg_ri_right = CFGProto(spin, middle);
        connected_sfs_right = findConnectedSFs_SOMOSOMO(p, q, cfg_right);
    }
    else if (exctype2 == ExcType::SV)
    {
        if (prel2 > qrel2)
        {
            p = prel2;
            q = qrel2;
        }
        else
        {
            p = prel2;
            q = qrel2 + 1;
        }

        norb = nue_right + 1;

        string middle(norb, '1');
        middle[q] = '0';
        string right(norb, '1');
        right[p] = '0';

        cfg_right = CFGProto(spin, right, sfs_right);
        cfg_ri_right = CFGProto(spin, middle);
        connected_sfs_right = findConnectedSFs_SOMOVirtual(p, q, cfg_right);
    }

    // <left|Epq|RI>
    if (exctype1 == ExcType::DS)
    {
        if (prel1 >= qrel1)
        {
            p = prel1 + 1;
            q = qrel1;
        }
        else
        {
            p = prel1;
            q = qrel1;
        }

        norb = nue_middle + 1;

        string left(norb, '1');
        string middle(norb, '1');
        left[p] = '2';
        middle[q] = '2';

        cfg_left = CFGProto(spin, left, sfs_left);
        cfg_ri_left = CFGProto(spin, middle);
        connected_sfs_left = findConnectedSFs_DOMOSOMO(q, p, cfg_left);
    }
    else if (exctype1 == ExcType::DV)
    {
        p = prel1;
        q = qrel1;

        norb = nue_left;

        string left(norb, '1');
        string middle(norb, '1');
        middle[p] = '0';
        middle[q] = '2';

        cfg_left = CFGProto(spin, left, sfs_left);
        cfg_ri_left = CFGProto(spin, middle);
        connected_sfs_left = findConnectedSFs_SOMOSOMO(q, p, cfg_left);
    }
    else if (exctype1 == ExcType::SS)
    {
        p = prel1;
        q = qrel1;

        norb = nue_middle;
        string middle(norb, '1');
        string left(norb, '1');
        left[p] = '2';
        left[q] = '0';

        cfg_left = CFGProto(spin, left, sfs_left);
        cfg_ri_left = CFGProto(spin, middle);
        connected_sfs_left = findConnectedSFs_DOMOVirtual(p, q, cfg_left);
    }
    else if (exctype1 == ExcType::SV)
    {
        if (prel1 > qrel1)
        {
            p = prel1;
            q = qrel1;
        }
        else
        {
            p = prel1;
            q = qrel1 + 1;
        }

        norb = nue_middle + 1;

        string left(norb, '1');
        left[q] = '0';
        string middle(norb, '1');
        middle[p] = '0';

        cfg_left = CFGProto(spin, left, sfs_left);
        cfg_ri_left = CFGProto(spin, middle);
        connected_sfs_left = findConnectedSFs_SOMOVirtual(q, p, cfg_left);
    }

    vector<string> sfs_intersection;
    std::set_intersection(connected_sfs_left.begin(), connected_sfs_left.end(), connected_sfs_right.begin(),
                          connected_sfs_right.end(), inserter(sfs_intersection, sfs_intersection.begin()));
    cfg_ri_right.createCSFsFromSFs(sfs_intersection);
    cfg_ri_left.createCSFsFromSFs(sfs_intersection);

    return std::make_tuple(std::move(cfg_left), std::move(cfg_ri_left),
                           std::move(cfg_ri_right), std::move(cfg_right));
}

proto_1el_tuple GCI::CouplingCoeffs::returnCFGPrototypesDia(const sfs_pair_t &sfs_pair,
                                                            const quintet &key)
{
    proto_1el_tuple cfg_prototypes;

    size_t exctype = get<0>(key);
    size_t nue_ri = get<1>(key);
    size_t nue_right = get<2>(key);
    size_t prel = get<3>(key);
    size_t qrel = get<4>(key);

    vector<string> sfs_right;
    for (const size_t &sf_idx : sfs_pair.second)
        sfs_right.push_back(gci->sfs_map__idx_to_sf[nue_right][sf_idx]);

    int norb, p, q;
    if (exctype == ExcType::DS)
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

        string left(norb, '1');
        string right(norb, '1');
        left[p] = '2';
        right[q] = '2';

        CFGProto cfg_right(spin, right, sfs_right);
        vector<string> sfs_ri = findConnectedSFs_DOMOSOMO(p, q, cfg_right);
        CFGProto cfg_ri(spin, left, sfs_ri);

        return std::make_tuple(std::move(cfg_ri), std::move(cfg_right));
    }
    else if (exctype == ExcType::DV)
    {
        p = prel;
        q = qrel;

        norb = nue_ri;

        string left(norb, '1');
        string right(norb, '1');
        right[p] = '0';
        right[q] = '2';

        CFGProto cfg_right(spin, right, sfs_right);
        vector<string> sfs_ri = findConnectedSFs_DOMOVirtual(q, p, cfg_right);
        CFGProto cfg_ri(spin, left, sfs_ri);

        return std::make_tuple(std::move(cfg_ri), std::move(cfg_right));
    }
    else if (exctype == ExcType::SS)
    {
        p = prel;
        q = qrel;

        int norb = nue_right;

        string left(norb, '1');
        string right(norb, '1');
        left[p] = '2';
        left[q] = '0';

        CFGProto cfg_right(spin, right, sfs_right);
        vector<string> sfs_ri = findConnectedSFs_SOMOSOMO(p, q, cfg_right);
        CFGProto cfg_ri(spin, left, sfs_ri);

        return std::make_tuple(std::move(cfg_ri), std::move(cfg_right));
    }
    else if (exctype == ExcType::SV)
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

        string left(norb, '1');
        left[q] = '0';
        string right(norb, '1');
        right[p] = '0';

        CFGProto cfg_right(spin, right, sfs_right);
        vector<string> sfs_ri = findConnectedSFs_SOMOVirtual(p, q, cfg_right);
        CFGProto cfg_ri(spin, left, sfs_ri);

        return std::make_tuple(std::move(cfg_ri), std::move(cfg_right));
    }
    else
        throw std::runtime_error("Error: Wrong exctype!");
}

void GCI::CouplingCoeffs::calcCCs1El(const sfs_pair_t &sfs_pair,
                                     const proto_1el_tuple &cfg_prototypes,
                                     const quintet &key,
                                     std::map<quintet, cc_map> &ccs_1el)
{
    int exctype = get<0>(key);
    int prel = get<3>(key);
    int qrel = get<4>(key);

    arma::dmat ccs;
    if (exctype == ExcType::DS)
        ccs = calcCCDOMOSOMO(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::DV)
        ccs = calcCCDOMOVirtual(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::SS)
        ccs = calcCCSOMOSOMO(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::SV)
        ccs = calcCCSOMOVirtual(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));

    size_t i = 0;
    for (auto &mu : sfs_pair.first)
    {
        size_t j = 0;
        for (auto &nu : sfs_pair.second)
        {
            ccs_1el[key][mu][nu] = ccs(i, j);
            j++;
        }
        i++;
    }
}

void GCI::CouplingCoeffs::calcCCs2El(const sfs_pair_t &sfs_pair,
                                     const proto_2el_tuple &cfg_prototypes,
                                     const nonet &key,
                                     std::map<nonet, cc_map> &ccs_2el)
{
    int exctype1 = get<0>(key);
    int exctype2 = get<1>(key);
    int prel1 = get<5>(key);
    int qrel1 = get<6>(key);
    int prel2 = get<7>(key);
    int qrel2 = get<8>(key);

    arma::dmat cc1, cc2;
    if (exctype1 == ExcType::DS and exctype2 == ExcType::DS)
    {
        cc1 = calcCCDOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DS and exctype2 == ExcType::DV)
    {
        cc1 = calcCCDOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DS and exctype2 == ExcType::SS)
    {
        cc1 = calcCCDOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DS and exctype2 == ExcType::SV)
    {
        cc1 = calcCCDOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DV and exctype2 == ExcType::DS)
    {
        cc1 = calcCCDOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DV and exctype2 == ExcType::DV)
    {
        cc1 = calcCCDOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DV and exctype2 == ExcType::SS)
    {
        cc1 = calcCCDOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::DV and exctype2 == ExcType::SV)
    {
        cc1 = calcCCDOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SS and exctype2 == ExcType::DS)
    {
        cc1 = calcCCSOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SS and exctype2 == ExcType::DV)
    {
        cc1 = calcCCSOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SS and exctype2 == ExcType::SS)
    {
        cc1 = calcCCSOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SS and exctype2 == ExcType::SV)
    {
        cc1 = calcCCSOMOSOMO(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SV and exctype2 == ExcType::DS)
    {
        cc1 = calcCCSOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SV and exctype2 == ExcType::DV)
    {
        cc1 = calcCCSOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCDOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SV and exctype2 == ExcType::SS)
    {
        cc1 = calcCCSOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOSOMO(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }
    else if (exctype1 == ExcType::SV and exctype2 == ExcType::SV)
    {
        cc1 = calcCCSOMOVirtual(prel1, qrel1, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
        cc2 = calcCCSOMOVirtual(prel2, qrel2, get<2>(cfg_prototypes), get<3>(cfg_prototypes));
    }

    arma::dmat ccs = cc1 * cc2;
    size_t i = 0;
    for (auto &mu : sfs_pair.first)
    {
        size_t j = 0;
        for (auto &nu : sfs_pair.second)
        {
            ccs_2el[key][mu][nu] = ccs(i, j);
            j++;
        }
        i++;
    }
}

void GCI::CouplingCoeffs::calcCCsDia(const sfs_pair_t &sfs_pair,
                                     const proto_1el_tuple &cfg_prototypes,
                                     const quintet &key,
                                     std::map<quintet, cc_map> &ccs_dia)
{
    int exctype = get<0>(key);
    int prel = get<3>(key);
    int qrel = get<4>(key);

    arma::dmat cc;
    if (exctype == ExcType::DS)
        cc = calcCCDOMOSOMO(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::DV)
        cc = calcCCDOMOVirtual(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::SS)
        cc = calcCCSOMOSOMO(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));
    else if (exctype == ExcType::SV)
        cc = calcCCSOMOVirtual(prel, qrel, get<0>(cfg_prototypes), get<1>(cfg_prototypes));

    arma::dmat ccs = cc.t() * cc;
    size_t i = 0;
    for (auto &mu : sfs_pair.first)
    {
        size_t j = 0;
        for (auto &nu : sfs_pair.second)
        {
            ccs_dia[key][mu][nu] = ccs(i, j);
            j++;
        }
        i++;
    }
}