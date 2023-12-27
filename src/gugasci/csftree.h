#pragma once

#include <array>
#include <memory>
#include <string>

namespace Lible
{
    namespace GUGA
    {
        class CSFTree
        {
        public:
            struct Node
            {
                bool end = false;
                int b = 0;
                int pos = -1;
                std::array<std::unique_ptr<Node>, 4> step_values{nullptr, nullptr, nullptr, nullptr};
            };

            void insertToTree(const int &pos, const std::string &csf)
            {
                Node *current = root.get();
                for (size_t i = 0; i < csf.size(); i++)
                {
                    int idx = csf[i] - '0';
                    if (current->step_values[idx] == nullptr)
                    {
                        current->step_values[idx] = std::make_unique<Node>();
                        if (idx == 1)
                            current->step_values[idx]->b = current->b + 1;
                        else if (idx == 2)
                            current->step_values[idx]->b = current->b - 1;
                        else
                            current->step_values[idx]->b = current->b;
                    }
                    current = current->step_values[idx].get();
                }
            }

            Node *getRoot()
            {
                return root.get();
            }

            Node *getRoot() const
            {
                return root.get();
            }

            void resetTree()
            {
                root.reset();
            }

        private:
            std::unique_ptr<Node> root;
        };
    }
}