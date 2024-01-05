#pragma once

#include <set>
#include <string>

#include <lible/guga_sci.h>
#include <lible/gci_util.h>

namespace lible
{
    namespace guga
    {
        class GCI::Impl::Connections
        {
        public:
            Connections(Impl *gci) : impl(gci)
            {
                min_nue = impl->min_nue;
                n_orbs = impl->n_orbs;
            }

            void constructConnections(const std::set<std::string> &cfgs_right,
                                      const wfn_ptr &wfn_left,
                                      const wfn_ptr &wfn_right,
                                      connection_map_1el &connections_1el,
                                      connection_map_2el &connections_2el,
                                      connection_map_dia &connections_dia);

        private:
            Impl *impl;
            
            size_t min_nue;
            size_t n_orbs;            

            void constructConnections_Epq(const int &icfg_right, const std::string &cfg_right,
                                          const wfn_ptr &wfn_left, connection_map_1el &connections_1el);

            void constructConnections_EpqEqp(const int &icfg_right, const std::string &cfg_right,
                                             connection_map_dia &connections_dia);

            void constructConnections_EpqEqr(const int &icfg_right, const std::string &cfg_right,
                                             const std::string &operators, const wfn_ptr &wfn_left,
                                             connection_map_2el &connections_2el);

            void constructConnections_EpqErp(const int &icfg_right, const std::string &cfg_right,
                                             const std::string &operators, const wfn_ptr &wfn_left,
                                             connection_map_2el &connections_2el);

            void constructConnections_EpqEpq(const int &icfg_right, const std::string &cfg_right,
                                             const wfn_ptr &wfn_left, connection_map_2el &connections_2el);

            void constructConnections_EpqEpr(const int &icfg_right, const std::string &cfg_right,
                                             const std::string &operators, const wfn_ptr &wfn_left,
                                             connection_map_2el &connections_2el);

            void constructConnections_EpqErq(const int &icfg_right, const std::string &cfg_right,
                                             const std::string &operators, const wfn_ptr &wfn_left,
                                             connection_map_2el &connections_2el);

            void constructConnections_EpqErs(const int &icfg_right, const std::string &cfg_right,
                                             const std::string &operators, const wfn_ptr &wfn_left,
                                             connection_map_2el &connections_2el);

            static constexpr auto annihilation = [](const int &i, std::string &cfg)
            {
                cfg[i] -= 1;
            };

            static constexpr auto annihilation2x = [](const int &i, std::string &cfg)
            {
                cfg[i] -= 2;
            };

            static constexpr auto creation = [](const int &i, std::string &cfg)
            {
                cfg[i] += 1;
            };

            static constexpr auto creation2x = [](const int &i, std::string &cfg)
            {
                cfg[i] += 2;
            };

            static constexpr auto testAnnihilation = [](const int &i, const std::string &cfg)
            {
                if (cfg[i] == '0')
                    return true;
                else
                    return false;
            };

            static constexpr auto testAnnihilation2x = [](const int &i, const std::string &cfg)
            {
                if (cfg[i] != '2')
                    return true;
                else
                    return false;
            };

            static constexpr auto testCreation = [](const int &i, const std::string &cfg)
            {
                if (cfg[i] == '2')
                    return true;
                else
                    return false;
            };

            static constexpr auto testCreation2x = [](const int &i, const std::string &cfg)
            {
                if (cfg[i] != '0')
                    return true;
                else
                    return false;
            };

            static constexpr auto identity = [](const int &i, std::string &cfg) {
            };

            static constexpr auto canonicalizeIdxs_EpqEqr = [](std::string operators,
                                                               std::vector<int> idxs)
            {
                std::vector<int> idxs_can;

                int idx_c;
                for (int i = 0; i < 3; i++)
                    if (operators[i] == 'c')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 2; i++)
                    if (operators[i] == '0')
                    {
                        idx_c = idxs[i];
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                idxs_can.push_back(idx_c);
                idxs_can.push_back(idxs[0]);
                return idxs_can;
            };

            static constexpr auto canonicalizeIdxs_EpqErp = [](std::string operators, std::vector<int> idxs)
            {
                std::vector<int> idxs_can;

                int idx_a;
                for (int i = 0; i < 3; i++)
                    if (operators[i] == '0')
                    {
                        idx_a = idxs[i];
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 2; i++)

                    if (operators[i] == 'a')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                idxs_can.push_back(idxs[0]);
                idxs_can.push_back(idx_a);
                return idxs_can;
            };

            static constexpr auto canonicalizeIdxs_EpqEpr = [](std::string operators, std::vector<int> idxs)
            {
                std::vector<int> idxs_can;

                int idx_c;
                for (int i = 0; i < 3; i++)
                    if (operators[i] == 'c')
                    {
                        idx_c = idxs[i];
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 2; i++)
                    if (operators[i] == 'a')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                idxs_can.push_back(idx_c);
                idxs_can.push_back(idxs[0]);
                return idxs_can;
            };

            static constexpr auto canonicalizeIdxs_EpqErq = [](std::string operators, std::vector<int> idxs)
            {
                std::vector<int> idxs_can;

                int idx_a;
                for (int i = 0; i < 3; i++)
                    if (operators[i] == 'c')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 2; i++)
                    if (operators[i] == 'a')
                    {
                        idx_a = idxs[i];
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                idxs_can.push_back(idxs[0]);
                idxs_can.push_back(idx_a);
                return idxs_can;
            };

            static constexpr auto canonicalizeIdxs_EpqErs = [](std::string operators, std::vector<int> idxs)
            {
                std::vector<int> idxs_can;

                for (int i = 0; i < 4; i++)                
                    if (operators[i] == 'c')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 3; i++)                
                    if (operators[i] == 'a')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }

                for (int i = 0; i < 2; i++)                
                    if (operators[i] == 'c')
                    {
                        idxs_can.push_back(idxs[i]);
                        idxs.erase(idxs.begin() + i);
                        operators.erase(operators.begin() + i);
                        break;
                    }      

                idxs_can.push_back(idxs[0]);
                return idxs_can;
            };
        };
    }
}