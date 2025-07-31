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
    TreeDAWG(const LabeledTree<LabelType> &tree)
    {
        auto nodes = tree.get_nodes();

        // Create a mapping from tree nodes to DAWG node IDs
        std::unordered_map<std::shared_ptr<Node<LabelType>>, uint64_t> node_to_id;

        // First pass: create all DAWG nodes
        for (const auto &node : nodes)
        {
            uint64_t dawg_node_id = this->add_node();
            node_to_id[node] = dawg_node_id;
        }

        // Second pass: add transitions between nodes
        for (const auto &tree_node : nodes)
        {
            uint64_t current_id = node_to_id[tree_node];

            // Create transitions for each child
            std::vector<std::pair<LabelType, uint64_t>> transitions;
            for (const auto &child : tree_node->get_children())
            {
                uint64_t child_id = node_to_id[child];
                transitions.emplace_back(child->get_label(), child_id);
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

private:
    /**
     * @brief Sets the initial state, ensuring it's a valid tree root (no incoming transitions).
     * @param node_id The ID of the node to set as initial state.
     * @throws std::invalid_argument if the node has incoming transitions.
     */
    void set_initial_state_impl(uint64_t node_id)
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
    int64_t compute_all_heights_impl()
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