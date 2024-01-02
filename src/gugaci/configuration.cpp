#include <lible/configuration.h>

#include <stdexcept>

using namespace lible::guga;

using std::array;
using std::map;
using std::string;
using std::vector;

CFG::CFG(const double &spin_, const string &occ_nr_vector_)
    : spin(spin_), occ_nr_vector(occ_nr_vector_)
{
    n_orb = occ_nr_vector.size();
    n_el = 0;
    for (int i = 0; i < n_orb; i++)
        n_el += (occ_nr_vector[i] - '0');

    n_unpaired_el = std::count(occ_nr_vector.begin(), occ_nr_vector.end(), '1');

    top[0] = 0.5 * n_el - spin;
    top[1] = 2.0 * spin;
    top[2] = n_orb - top[0] - top[1];
}

vector<string> CFG::extractSFs()
{
    vector<string> sfs;
    if (n_unpaired_el == 0)
        sfs = vector<string>({""});
    else
        for (const string &csf : csfs)
        {
            string sf(n_unpaired_el, '0');
            int idx_sf = 0;
            for (const char &d : csf)
            {
                if (d == '1')
                {
                    sf[idx_sf] = '+';
                    idx_sf++;
                }
                else if (d == '2')
                {
                    sf[idx_sf] = '-';
                    idx_sf++;
                }
            }
            sfs.push_back(sf);
        }
    sfs.shrink_to_fit();

    return sfs;
}

void CFG::createAllCSFs()
{
    resetTree();
    sf_idxs.clear();
    csfs.clear();

    string blank_csf(n_orb, '0');
    array<int, 3> blank_row{0, 0, 0};
    createAllCSFsRecursively('0', 0, blank_row, blank_csf);
}

void CFG::createCSFsFromSFs(const map<string, int> &sfs)
{
    resetTree();
    csfs.clear();
    sf_idxs.clear();

    if (n_unpaired_el == 0)
    {
        string csf(n_orb, '0');
        for (size_t i = 0; i < occ_nr_vector.size(); i++)
        {
            char occ_nr = occ_nr_vector[i];
            if (occ_nr == '0')
                csf[i] = '0';
            else if (occ_nr == '2')
                csf[i] = '3';
        }
        size_t pos = 0;
        csfs.push_back(csf);
        sf_idxs.push_back(pos);
    }
    else
    {
        for (auto &item : sfs)
        {
            string sf = item.first;
            string csf(n_orb, '0');
            int idx_sf = 0;
            for (size_t i = 0; i < occ_nr_vector.size(); i++)
            {
                char occ_nr = occ_nr_vector[i];
                if (occ_nr == '0')
                    csf[i] = '0';
                else if (occ_nr == '2')
                    csf[i] = '3';
                else if (sf[idx_sf] == '+')
                    csf[i] = '1';
                else
                    csf[i] = '2';
                idx_sf++;
            }
            size_t pos = item.second;
            csfs.push_back(csf);
            sf_idxs.push_back(pos);
        }
    }
}

void CFG::insertCSF(const size_t &sf_idx, const string &csf)
{
    csfs.push_back(csf);
    sf_idxs.push_back(sf_idx);
    insertToTree(sf_idx, csf);
}

array<int, 3> CFG::determineStepRow(const char &d)
{
    if (d == '0')
        return {0, 0, 1};
    else if (d == '1')
        return {0, 1, 0};
    else if (d == '2')
        return {1, -1, 1};
    else if (d == '3')
        return {1, 0, 0};
    else
        throw std::runtime_error("Error: wrong step value d!");
}

void CFG::createAllCSFsRecursively(char d, int i, array<int, 3> row, string csf)
{
    if (i != 0)
    {
        array<int, 3> step_row = determineStepRow(d);
        row[0] += step_row[0];
        row[1] += step_row[1];
        row[2] += step_row[2];
        csf[i - 1] = d;
    }

    if (row[1] < 0)
        return;
    if ((2 * row[0] + row[1]) > n_el)
        return;

    if (i == n_orb)
    {
        if (row == top)
        {
            int pos = csfs.size();
            sf_idxs.push_back(pos);
            csfs.push_back(csf);
            insertToTree(pos, csf);
            return;
        }
        else
            return;
    }

    if (occ_nr_vector[i] == '0')
        createAllCSFsRecursively('0', i + 1, row, csf);
    else if (occ_nr_vector[i] == '1')
    {
        createAllCSFsRecursively('1', i + 1, row, csf);
        createAllCSFsRecursively('2', i + 1, row, csf);
    }
    else if (occ_nr_vector[i] == '2')
        createAllCSFsRecursively('3', i + 1, row, csf);
}

CFGProto::CFGProto(const double &spin_, const string &occ_nr_vector_)
    : spin(spin_), occ_nr_vector(occ_nr_vector_)
{
    n_orb = occ_nr_vector.size();
    n_unpaired_el = std::count(occ_nr_vector.begin(), occ_nr_vector.end(), '1');
}

CFGProto::CFGProto(const double &spin_, const string &occ_nr_vector_,
                   const vector<string> &sfs)
    : spin(spin_), occ_nr_vector(occ_nr_vector_)
{
    createCSFsFromSFs(sfs);
}

void CFGProto::createCSFsFromSFs(const vector<string> &sfs)
{
    resetTree();
    csfs.clear();

    if (n_unpaired_el == 0)
    {
        string csf(n_orb, 0);
        for (size_t i = 0; i < occ_nr_vector.size(); i++)
        {
            char occ_nr = occ_nr_vector[i];
            if (occ_nr == '0')
                csf[i] = '0';
            else if (occ_nr == '2')
                csf[i] = '3';

            size_t pos = csfs.size();
            csfs.push_back(csf);
            insertToTree(pos, csf);
        }
    }
    else
    {
        for (const string &sf : sfs)
        {
            string csf(n_orb, '0');
            int sf_idx = 0;
            for (size_t i = 0; i < occ_nr_vector.size(); i++)
            {
                char occ_nr = occ_nr_vector[i];
                if (occ_nr == '0')
                    csf[i] = '0';
                else if (occ_nr == '2')
                    csf[i] = '3';
                else
                {
                    if (sf[sf_idx] == '+')
                        csf[i] = '1';
                    else
                        csf[i] = '2';
                    sf_idx++;
                }
            }

            size_t pos = csfs.size();
            csfs.push_back(csf);
            insertToTree(pos, csf);
        }
    }
    csfs.shrink_to_fit();
}