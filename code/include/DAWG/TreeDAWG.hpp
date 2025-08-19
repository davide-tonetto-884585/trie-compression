#ifndef TREE_DAWG_HPP
#define TREE_DAWG_HPP

#include <stack>
#include <unordered_map>
#include <algorithm>
#include <string>
#include "DAWG.hpp"
#include "../LabeledTree/LabeledTree.hpp"

/**
 * @brief TreeDAWG is a specialized subclass of DAWG optimized for tree structures.
 *
 * This class overrides the calculate_height_dfs method to provide an iterative,
 * tree-optimized implementation that avoids cycle detection overhead and uses
 * post-order traversal for better performance on tree structures.
 */
template <typename LabelType = char>
class TreeDAWG : public DAWG<TreeDAWG<LabelType>, LabelType>
{
    friend class DAWG<TreeDAWG<LabelType>, LabelType>;

public:
    // Using declarations to access base class members
    using DAWG<TreeDAWG<LabelType>, LabelType>::m_nodes;
    using DAWG<TreeDAWG<LabelType>, LabelType>::m_initial_state_id;
    using DAWG<TreeDAWG<LabelType>, LabelType>::get_num_nodes;

    /**
     * @brief Default constructor for the TreeDAWG class.
     */
    TreeDAWG() : DAWG<TreeDAWG<LabelType>, LabelType>() {}

    /**
     * @brief Constructor from LabeledTree for the TreeDAWG class.
     */
    TreeDAWG(const LabeledTree<LabelType> &tree) : DAWG<TreeDAWG<LabelType>, LabelType>()
    {
        auto nodes = tree.get_nodes();

        // Handle empty tree case
        if (nodes.empty()) {
            return;
        }

        // Create a mapping from tree nodes to DAWG node IDs
        std::unordered_map<std::shared_ptr<Node<LabelType>>, uint64_t> node_to_id;

        // First pass: create all DAWG nodes
        for (const auto &node : nodes)
        {
            uint64_t dawg_node_id = this->add_node();
            node_to_id[node] = dawg_node_id;
        }

        // Second pass: add transitions between nodes and set parent IDs
        for (const auto &tree_node : nodes)
        {
            uint64_t current_id = node_to_id[tree_node];

            // Create transitions for each child
            std::vector<std::pair<LabelType, uint64_t>> transitions;
            for (const auto &child : tree_node->get_children())
            {
                uint64_t child_id = node_to_id[child];
                transitions.emplace_back(child->get_label(), child_id);
                
                // Set the parent ID for the child node
                this->m_nodes[child_id].set_parent_state_id(current_id);
            }

            // Configure the state with its transitions
            if (!transitions.empty())
            {
                typename State<LabelType>::Builder builder;
                for (const auto &[label, target_id] : transitions)
                {
                    builder.add_transition(label, target_id);
                }
                builder.build_into(this->m_nodes[current_id]);
            }
        }

        // Set the root as the initial state
        if (tree.get_root() != nullptr)
        {
            this->set_initial_state(node_to_id[tree.get_root()]);
        }
    }

    /**
     * @brief Stable sorts nodes based on their label path from parent edge transition to root.
     * 
     * This method creates a path of transition labels from each node up to the root,
     * then performs a stable sort based on lexicographic comparison of these paths.
     * The path for each node consists of the labels of transitions from the root down to
     * the parent of that node.
     * 
     * @return A vector of node IDs sorted by their label paths to the root.
     */
    std::vector<uint64_t> stable_sort_nodes_by_label_path(bool verbose = false)
    {
        std::vector<uint64_t> node_ids;
        for (uint64_t i = 0; i < get_num_nodes(); ++i)
        {
            node_ids.push_back(i);
        }

        // Create a map to store label paths for each node
        std::vector<std::string> node_paths(get_num_nodes());
        
        // Build label paths for each node
        for (uint64_t node_id : node_ids)
        {
            std::string path;
            uint64_t current_id = node_id;
            
            // Traverse up to root, collecting transition labels
            while (this->m_nodes[current_id].has_parent())
            {
                uint64_t parent_id = this->m_nodes[current_id].get_parent_state_id();
                
                // Find the transition label from parent to current node
                const auto& parent_transitions = this->m_nodes[parent_id].get_transitions();
                for (const auto& [label, target_id] : parent_transitions)
                {
                    if (target_id == current_id)
                    {
                        path += std::string(1, label);
                        break;
                    }
                }
                current_id = parent_id;
            }
            
            node_paths[node_id] = path;
        }

        // Stable sort based on lexicographic comparison of label paths
        std::stable_sort(node_ids.begin(), node_ids.end(),
            [&node_paths](uint64_t a, uint64_t b) {
                return node_paths[a] < node_paths[b];
            });

        if (verbose)
        {
            std::cout << "Node order: ";
            for (uint64_t i : node_ids)
            {
                std::cout << i << "(" << node_paths[i] << ") ";
            }
            std::cout << std::endl;
        }

        return node_ids;
    }

private:
    /**
     * @brief Validates that a node has no incoming transitions (is a root).
     * @param node_id The ID of the node to validate.
     * @throws std::invalid_argument if the node has incoming transitions.
     */
    void validate_is_root(uint64_t node_id)
    {
        uint64_t num_nodes = get_num_nodes();
        for (uint64_t i = 0; i < num_nodes; ++i)
        {
            if (i == node_id)
                continue; // Skip self

            const auto &transitions = (*this)[i].get_transitions();
            for (const auto &[_, target_id] : transitions)
            {
                if (target_id == node_id)
                {
                    throw std::invalid_argument("Node " + std::to_string(node_id) +
                                                " cannot be set as initial state because it has incoming transitions from node " +
                                                std::to_string(i) + ". In a tree, the root must have no incoming transitions.");
                }
            }
        }
    }
};

#endif // TREE_DAWG_HPP