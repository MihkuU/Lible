#pragma once

#include <map>
#include <set>
#include <vector>

#include <lible/gugaci/csftree.hpp>

namespace lible
{
    namespace guga
    {
        class CFG : public CSFTree
        {
        public:
            CFG(const double &spin_, const std::string &occ_nr_vector_);

            std::vector<std::string> extractSFs();

            void createAllCSFs();

            void createCSFsFromSFs(const std::map<std::string, int> &sfs);

            void insertCSF(const size_t &sf_idx, const std::string &csf);

            int getNUE() const
            {
                return n_unpaired_el;
            }

            size_t getNumCSFs() const
            {
                return csfs.size();
            }

            std::string getONV() const
            {
                return occ_nr_vector;
            }

            std::string getCSF(const int &icsf) const
            {
                return csfs[icsf];
            }

            std::vector<int> getSFIdxs() const
            {
                return sf_idxs;
            }

            std::vector<std::string> getCSFs() const
            {
                return csfs;
            }

        private:
            std::array<int, 3> determineStepRow(const char &d);

            void createAllCSFsRecursively(char d, int i, std::array<int, 3> row, std::string csf);

            double spin;
            int n_el;
            int n_orb;
            int n_unpaired_el;
            int pos;

            std::array<int, 3> top;

            std::string occ_nr_vector;
            std::vector<int> sf_idxs;
            std::vector<std::string> csfs;
        };

        class CFGProto : public CSFTree
        {
        public:
            CFGProto() {}
            CFGProto(const double &spin_, const std::string &occ_nr_vector_);
            CFGProto(const double &spin_, const std::string &occ_nr_vector_,
                     const std::vector<std::string> &sfs);

            void createCSFsFromSFs(const std::vector<std::string> &sfs);

            size_t getNumCSFs() const
            {
                return csfs.size();
            }

            std::string getONV() const
            {
                return occ_nr_vector;
            }

            std::vector<std::string> getCSFs() const
            {
                return csfs;
            }

        private:
            double spin;
            int n_orb;
            int n_unpaired_el;

            std::string occ_nr_vector;
            std::vector<std::string> csfs;
        };
    }
}