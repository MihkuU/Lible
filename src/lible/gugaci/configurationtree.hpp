#pragma once

#include <memory>
#include <string>

namespace lible
{
    namespace guga
    {
        class CFGTree
        {
        public:
            struct Node
            {   
                bool end = false;
                int pos = -1;
                std::array<std::unique_ptr<Node>, 3> occ_nrs{nullptr, nullptr, nullptr};
            };

            CFGTree()
            {
                root = std::make_unique<Node>();
            }

            int findCFGPos(const std::string &cfg) const
            {
                Node *current = root.get();
                for (size_t i = 0; i < cfg.size(); i++)
                {
                    current = current->occ_nrs[cfg[i] - '0'].get();
                    if (current == nullptr)
                        return -1;
                }
                return current->pos;
            }

            int findCFGPos(const int &start, const std::string &cfg, Node *node) const
            {
                Node *current = node;
                for (size_t i = start; i < cfg.size(); i++)
                {
                    current = current->occ_nrs[cfg[i] - '0'].get();
                    if (current == nullptr)
                        return -1;
                }

                return current->pos;
            }

            Node *incrementNode(const int &start, const std::string &cfg, Node *node) const
            {
                return node->occ_nrs[cfg[start] - '0'].get();
            }

            Node *searchFromRoot(const size_t &stop, const std::string &cfg) const
            {
                Node *current = root.get();
                for (size_t i = 0; i < stop; i++)
                {
                    current = current->occ_nrs[cfg[i] - '0'].get();
                    if (current == nullptr)
                        return nullptr;
                }

                return current;
            }

            void insertToTree(const int &pos, const std::string &cfg)
            {
                Node *current = root.get();
                for (size_t i = 0; i < cfg.size(); i++)
                {
                    int idx = cfg[i] - '0';
                    if (current->occ_nrs[idx] == nullptr)
                        current->occ_nrs[idx] = std::make_unique<Node>();
                    current = current->occ_nrs[idx].get();
                }

                current->end = true;
                current->pos = pos;
            }

        private:
            std::unique_ptr<Node> root;            
        };
    }
}