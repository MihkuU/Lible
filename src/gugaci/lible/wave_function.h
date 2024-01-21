#pragma once

#include <string>
#include <vector>

#include <lible/configuration.h>
#include <lible/configurationtree.h>

namespace lible
{
    namespace guga
    {
        class WaveFunction : public CFGTree
        {
        public:
            WaveFunction() 
            { 
                spin = 0;
                pos_last_csf = 0;
            }

            WaveFunction(const double &spin) : spin(spin)
            {
                pos_last_csf = 0;
            }

            int insertCFGandGetPos(CFG &cfg)
            {
                int pos_cfg = occ_nr_vectors.size();
                insertToTree(pos_cfg, cfg.getONV());

                size_t dim_csfs = cfg.getNumCSFs();
                occ_nr_vectors.push_back(cfg.getONV());
                dimensions.push_back(cfg.getNumCSFs());
                positions.push_back(pos_last_csf);
                configurations.emplace_back(std::move(cfg));
                pos_last_csf += dim_csfs;

                return pos_cfg;
            }

            void insertCFG(CFG &cfg)
            {
                size_t pos_cfg = occ_nr_vectors.size();
                insertToTree(pos_cfg, cfg.getONV());

                size_t dim_csfs = cfg.getNumCSFs();

                occ_nr_vectors.push_back(cfg.getONV());
                dimensions.push_back(cfg.getNumCSFs());
                positions.push_back(pos_last_csf);
                configurations.emplace_back(std::move(cfg));
                pos_last_csf += dim_csfs;
            }

            size_t getNumCFGs() const
            {
                return configurations.size();
            }

            size_t getNumCSFs() const
            {
                return pos_last_csf;
            }

            size_t getDim(const size_t &icfg) const
            {
                return dimensions[icfg];
            }

            size_t getPos(const size_t &icfg) const
            {
                return positions[icfg];
            }

            std::string getONV(const size_t &icfg) const
            {
                return configurations[icfg].getONV();
            }

            CFG *getCFGPtr(const size_t &icfg)
            {
                return &configurations[icfg];
            }

            const CFG *getCFGPtr(const size_t &icfg) const
            {
                return &configurations[icfg];
            }

        private:
            double spin;
            size_t pos_last_csf;

            std::vector<size_t> dimensions;
            std::vector<size_t> positions;
            std::vector<std::string> occ_nr_vectors;
            std::vector<CFG> configurations;
            std::vector<std::vector<double>> ci_coeffs;
        };
    }
}