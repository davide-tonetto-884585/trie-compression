#ifndef TREECOMPRESSOR_HPP
#define TREECOMPRESSOR_HPP

#include <vector>
#include <unordered_map>
#include <set>
#include "../DAWG/TreeDAWG.hpp"
#include "../Triplet/Triplet.hpp"

/**
 * @brief A class for compressing tree structures by merging nodes with equivalent classes.
 *
 * The TreeCompressor takes a TreeDAWG and a set of node chains, then creates a compressed
 * representation by assigning new node IDs based on equivalence classes. Nodes that belong
 * to the same equivalence class within a chain are merged to reduce redundancy.
 *
 * @tparam T The type of labels used in the tree transitions
 */
template <typename T>
class TreeCompressor
{
private:
    uint64_t m_num_new_nodes;
    std::set<Triplet<uint64_t, T, uint64_t>> m_transitions;

public:
    /**
     * @brief Constructs a TreeCompressor and performs the compression.
     *
     * This constructor takes a TreeDAWG and a collection of node chains, then creates
     * a compressed representation by:
     * 1. Assigning new node IDs based on equivalence classes within chains
     * 2. Merging consecutive nodes that belong to the same equivalence class
     * 3. Creating a set of compressed transitions using the new node IDs
     *
     * @param TreeDAWG The input TreeDAWG structure to be compressed
     * @param node_id_chains A vector of chains, where each chain is a vector of node IDs
     *                       representing a path through the tree
     * @param verbose If true, prints compression statistics and transition details
     */
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

        m_num_new_nodes = new_node_id;

        for (const auto node : TreeDAWG)
        {
            uint64_t node_id = node.get_id();
            uint64_t new_node_id = node_id_map[node_id];
            for (const auto &transition : node.get_transitions())
            {
                uint64_t new_dest_id = node_id_map[transition.second];
                m_transitions.emplace(new_node_id, transition.first, new_dest_id);

                // Check if the reversed transition exists in the set
                /* Triplet<uint64_t, T, uint64_t> reversed_transition(new_dest_id, transition.first, new_node_id);
                if (m_transitions.find(reversed_transition) != m_transitions.end())
                {
                    throw std::runtime_error("cycle detected");
                } */
            }
        }

        if (verbose)
        {
            for (const auto &transition : m_transitions)
            {
                std::cout << "(" << transition.first << ", " << transition.second << ", " << transition.third << ")" << std::endl;
            }
        }
    }

    /**
     * @brief Gets the compressed transitions as a set of triplets.
     *
     * Each triplet represents a transition in the compressed tree structure,
     * containing (source_node_id, label, destination_node_id).
     *
     * @return A const reference to the set of compressed transitions
     */
    const std::set<Triplet<uint64_t, T, uint64_t>> &get_transitions() const
    {
        return m_transitions;
    }

    /**
     * @brief Gets the number of new nodes created during compression.
     *
     * @return The total number of unique nodes in the compressed tree
     */
    uint64_t get_num_new_nodes() const
    {
        return m_num_new_nodes;
    }

    /**
     * @brief Checks if the compressed tree has non-deterministic finite automaton (NFA) properties.
     *
     * This method determines if there are two or more transitions with the same label
     * originating from the same node, which would indicate non-deterministic behavior.
     * In a deterministic finite automaton (DFA), each node can have at most one outgoing
     * transition for each label.
     *
     * @return true if the tree has NFA properties (multiple transitions with same label
     *         from the same node), false if it's deterministic
     */
    bool is_nfa() const
    {
        // Group transitions by source node and label
        std::unordered_map<uint64_t, std::unordered_map<T, uint64_t>> node_label_count;
        
        for (const auto& transition : m_transitions)
        {
            uint64_t source_node = transition.first;
            T label = transition.second;
            
            // Count occurrences of each label for each source node
            node_label_count[source_node][label]++;
            
            // If we find more than one transition with the same label from the same node
            if (node_label_count[source_node][label] > 1)
            {
                return true; // NFA detected
            }
        }
        
        return false; // No duplicate labels found, it's deterministic
    }

    
};

#endif // TREECOMPRESSOR_HPP
