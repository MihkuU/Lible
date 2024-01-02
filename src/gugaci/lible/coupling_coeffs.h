#pragma once

#include <map>
#include <set>

#include <lible/guga_sci.h>
#include <lible/gci_util.h>

namespace lible
{
    namespace guga
    {
        class SCI::CouplingCoeffs
        {
        public:
            CouplingCoeffs(SCI *gci_) : sci(gci_)
            {
                spin = sci->spin;
            }

            void constructCCs(const sf_pair_map_1el &sf_pairs_1el_new,
                              const sf_pair_map_2el &sf_pairs_2el_new,
                              std::map<nonet, cc_map> &ccs_2el,
                              std::map<quintet, cc_map> &ccs_1el);

            void constructCCs(const connection_map_1el &connections_1el_new,
                              const connection_map_2el &connections_2el_new,
                              const connection_map_dia &connections_dia_new,
                              const wfn_ptr &wave_function,
                              std::map<nonet, cc_map> &ccs_2el,
                              std::map<quintet, cc_map> &ccs_1el,
                              std::map<quintet, cc_map> &ccs_dia);

        private:
            SCI *sci;
            double spin;

            sf_pair_map_1el sf_pairs_1el;
            sf_pair_map_2el sf_pairs_2el;
            sf_pair_map_1el sf_pairs_dia;

            sfs_pair_t returnNewSFPair(const sfs_pair_t &sfs_pair_trial,
                                       const sf_pair_map_1el &sf_pairs_1el,
                                       const quintet &key);

            sfs_pair_t returnNewSFPair(const sfs_pair_t &sfs_pair,
                                       const sf_pair_map_2el &sf_pairs_2el,
                                       const nonet &key);

            sfs_pair_t returnNewSFPair(const connection_list_1el &connections,
                                       const sf_pair_map_1el &sf_pairs_1el,
                                       const quintet &key,
                                       const wfn_ptr &wave_function);

            sfs_pair_t returnNewSFPair(const connection_list_2el &connections,
                                       const sf_pair_map_2el &sf_pairs_2el,
                                       const nonet &key,
                                       const wfn_ptr &wave_function);

            sfs_pair_t returnNewSFPair(const connection_list_dia &connections,
                                       const sf_pair_map_1el &sf_pairs_dia,
                                       const quintet &key,
                                       const wfn_ptr &wave_function);

            sfs_pair_t returnSFPair(const std::set<int> &sfs_left,
                                    const std::set<int> &sfs_right);

            proto_1el_tuple returnCFGPrototypes1El(const sfs_pair_t &sfs_pair,
                                                   const quintet &key);

            proto_2el_tuple returnCFGPrototypes2El(const sfs_pair_t &sfs_pair,
                                                   const nonet &key);

            proto_1el_tuple returnCFGPrototypesDia(const sfs_pair_t &sfs_pair,
                                                   const quintet &key);

            void calcCCs1El(const sfs_pair_t &sfs_pair,
                            const proto_1el_tuple &cfg_prototypes,
                            const quintet &key,
                            std::map<quintet, cc_map> &ccs_1el_deprecated);

            void calcCCs2El(const sfs_pair_t &sfs_pair,
                            const proto_2el_tuple &cfg_prototypes,
                            const nonet &key,
                            std::map<nonet, cc_map> &ccs_2el);

            void calcCCsDia(const sfs_pair_t &sfs_pair,
                            const proto_1el_tuple &cfg_prototypes,
                            const quintet &key,
                            std::map<quintet, cc_map> &ccs_dia);
        };
    }
}