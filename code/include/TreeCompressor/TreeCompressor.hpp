#ifndef TREECOMPRESSOR_HPP
#define TREECOMPRESSOR_HPP

#include <vector>
#include <unordered_map>
#include <set>
#include "../DAWG/TreeDAWG.hpp"
#include "../Triplet/Triplet.hpp"

template <typename T>
class TreeCompressor
{
private:
    std::set<Triplet<uint64_t, T, uint64_t>> m_transitions;

public:
    TreeCompressor(const TreeDAWG<T> &TreeDAWG, const std::vector<std::vector<uint64_t>> &node_id_chains, bool verbose = false)
    {
        uint64_t new_node_id = 0;
        std::unordered_map<uint64_t, uint64_t> node_id_map;
        for (const auto &chain : node_id_chains)
        {
            uint64_t prev_class = UINT64_MAX;
            for (uint64_t node_id : chain)
            {
                if (prev_class != TreeDAWG[node_id].get_equivalence_class())
                {
                    node_id_map[node_id] = ++new_node_id;
                    prev_class = TreeDAWG[node_id].get_equivalence_class();
                }
                else
                {
                    node_id_map[node_id] = new_node_id;
                }
            }
        }

        for (const auto node: TreeDAWG)
        {
            uint64_t node_id = node.get_id();
            uint64_t new_node_id = node_id_map[node_id];
            for (const auto &transition : node.get_transitions())
            {
                uint64_t new_dest_id = node_id_map[transition.second];
                m_transitions.emplace(new_node_id, transition.first, new_dest_id);
            }
        }

        if (verbose)
        {
            std::cout << "Compressed transitions: " << m_transitions.size() << std::endl;
        }

    }
};

#endif // TREECOMPRESSOR_HPP
