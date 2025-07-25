#ifndef TREE_DAWG_HPP
#define TREE_DAWG_HPP

#include <stack>
#include "DAWG.hpp"
#include "../LabeledTree/LabeledTree.hpp"

/**
 * @brief TreeDAWG is a specialized subclass of DAWG optimized for tree structures.
 *
 * This class overrides the calculate_height_dfs method to provide an iterative,
 * tree-optimized implementation that avoids cycle detection overhead and uses
 * post-order traversal for better performance on tree structures.
 */
template<typename LabelType = char>
class TreeDAWG : public DAWG<LabelType>
{
public:
    // Using declarations to access base class members
    using DAWG<LabelType>::m_nodes;
    using DAWG<LabelType>::m_initial_state_id;
    using DAWG<LabelType>::get_num_nodes;
    
    /**
     * @brief Default constructor for the TreeDAWG class.
     */
    TreeDAWG() : DAWG<LabelType>() {}

    /**
     * @brief 
     */
    TreeDAWG(const LabeledTree<LabelType> &tree) {
        
    }

    /**
     * @brief Sets the initial state, ensuring it's a valid tree root (no incoming transitions).
     * @param node_id The ID of the node to set as initial state.
     * @throws std::invalid_argument if the node has incoming transitions.
     */
    void set_initial_state(uint64_t node_id) override
    {
        // First validate that the node exists
        if (node_id >= get_num_nodes())
        {
            throw std::out_of_range("Attempt to set non-existent initial state with ID " + std::to_string(node_id));
        }

        // Validate that this node has no incoming transitions (is a root)
        validate_is_root(node_id);

        // Set the initial state
        m_initial_state_id = node_id;
    }

    /**
     * @brief Computes the height of all nodes for a DAWG that is a tree.
     * This version is optimized for trees and does not perform cycle detection.
     * It uses an iterative post-order traversal to avoid stack overflow on deep trees.
     * @return The maximum height of any node in the tree.
     */
    int64_t compute_all_heights() override
    {
        for (auto &node : m_nodes)
            node.m_height = -1;

        int64_t max_h = -1;
        // compute height starting from root
        std::stack<uint64_t> dfs_stack;
        std::stack<uint64_t> processing_stack; // post-order stack

        dfs_stack.push(m_initial_state_id);

        std::vector<bool> visited(m_nodes.size(), false);
        while (!dfs_stack.empty())
        {
            uint64_t current_node_id = dfs_stack.top();
            dfs_stack.pop();

            if (visited[current_node_id])
                continue;

            visited[current_node_id] = true;
            processing_stack.push(current_node_id);

            for (const auto &[_, target_node_id] : m_nodes[current_node_id].m_transitions)
            {
                if (!visited[target_node_id])
                {
                    dfs_stack.push(target_node_id);
                }
            }
        }

        while (!processing_stack.empty())
        {
            uint64_t node_id = processing_stack.top();
            processing_stack.pop();

            State<LabelType> &node = m_nodes[node_id];

            if (node.m_height != -1)
                continue;

            int64_t current_max_h = node.is_final() ? 0 : -2;
            for (const auto &[_, target_node_id] : node.m_transitions)
            {
                int64_t target_height = m_nodes[target_node_id].m_height;

                if (target_height != -2)
                {
                    current_max_h = std::max(current_max_h, 1 + target_height);
                }
            }

            node.m_height = current_max_h;
        }

        // find max height
        for (const auto &node : m_nodes)
        {
            max_h = std::max(max_h, node.m_height);
        }

        return max_h;
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