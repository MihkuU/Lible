#pragma once

#include <array>
#include <memory>
#include <string>

namespace lible
{
    namespace guga
    {
        class CSFTree
        {
        public:
            struct Node
            {
                Node()
                {
                    end = false;                    
                    b = 0;
                    pos = -1;
                    step_values = {nullptr, nullptr, nullptr, nullptr};
                }

                bool end;
                int b;
                int pos;
                std::array<std::unique_ptr<Node>, 4> step_values;
            };

            CSFTree()
            {                
                root = std::make_unique<Node>();
            }
            ~CSFTree() = default;

            CSFTree(const CSFTree &) = delete;
            CSFTree &operator=(const CSFTree &) = delete;

            CSFTree(CSFTree &&) = default;
            CSFTree &operator=(CSFTree &&) = default;

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

                current->end = true;
                current->pos = pos;
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
                root.reset(new Node());
            }

        private:
            std::unique_ptr<Node> root;
        };
    }
}