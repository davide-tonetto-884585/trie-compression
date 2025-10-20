#ifndef TREE_DAWG_HPP
#define TREE_DAWG_HPP

#include <stack>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <thread>
#include <numeric>
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
        if (nodes.empty())
        {
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
     * This parallelized version uses C++17 std::execution parallel algorithms
     * to accelerate the path building phase, which works on macOS and other platforms.
     * 
     * @return A vector of node IDs sorted by their label paths to the root.
     */
    std::vector<uint64_t> stable_sort_nodes_by_label_path(bool verbose = false)
    {
        const uint64_t num_nodes = get_num_nodes();
        std::vector<uint64_t> node_ids(num_nodes);
        
        // Parallelize node ID initialization using std::iota and parallel execution
        std::iota(node_ids.begin(), node_ids.end(), 0);

        // Create a map to store label paths for each node
        std::vector<std::string> node_paths(num_nodes);

        // Parallelize path building for each node (most expensive operation)
        // Using thread-based parallelization for better compiler compatibility
        const unsigned int num_threads = std::thread::hardware_concurrency();
        const uint64_t chunk_size = (num_nodes + num_threads - 1) / num_threads;
        std::vector<std::thread> threads;
        
        for (unsigned int t = 0; t < num_threads; ++t)
        {
            uint64_t start = t * chunk_size;
            uint64_t end = std::min(start + chunk_size, num_nodes);
            
            if (start < end)
            {
                threads.emplace_back([this, &node_paths, start, end]()
                {
                    for (uint64_t node_id = start; node_id < end; ++node_id)
                    {
                        std::string path;
                        uint64_t current_id = node_id;

                        // Traverse up to root, collecting transition labels
                        while (this->m_nodes[current_id].has_parent())
                        {
                            uint64_t parent_id = this->m_nodes[current_id].get_parent_state_id();

                            // Find the transition label from parent to current node
                            const auto &parent_transitions = this->m_nodes[parent_id].get_transitions();
                            for (const auto &[label, target_id] : parent_transitions)
                            {
                                if (target_id == current_id)
                                {
                                    path += std::string(1, label);
                                    break;
                                }
                            }
                            current_id = parent_id;
                        }

                        node_paths[node_id] = std::move(path);
                    }
                });
            }
        }
        
        // Wait for all threads to complete
        for (auto& thread : threads)
        {
            thread.join();
        }

        // Use stable sort to maintain deterministic ordering for equal elements
        // For better performance with large datasets, you can use:
        // std::sort(...) if stability is not required
        std::stable_sort(node_ids.begin(), node_ids.end(),
                         [&node_paths](uint64_t a, uint64_t b)
                         {
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

    /**
     * @brief Gets the total number of transitions in the DAWG.
     *
     * This method counts all transitions across all nodes in the DAWG.
     * Each transition represents an edge from one node to another.
     *
     * @return The total number of transitions in the DAWG.
     */
    uint64_t get_num_transitions() const
    {
        uint64_t total_transitions = 0;
        for (const auto &node : this->m_nodes)
        {
            total_transitions += node.get_transitions().size();
        }
        return total_transitions;
    }

    /**
     * @brief Checks if the TreeDAWG is deterministic.
     *
     * A DAWG is deterministic if no node has two or more transitions with the same label.
     * This method iterates through all nodes and their transitions to verify this property.
     *
     * @return true if the TreeDAWG is deterministic, false otherwise.
     */
    bool is_deterministic() const
    {
        for (const auto &node : this->m_nodes)
        {
            const auto &transitions = node.get_transitions();
            
            // Check for duplicate labels in transitions
            // Since transitions are sorted by label (done in State::Builder::build_into),
            // we only need to check adjacent transitions
            for (size_t i = 1; i < transitions.size(); ++i)
            {
                if (transitions[i - 1].first == transitions[i].first)
                {
                    // Found duplicate label
                    return false;
                }
            }
        }
        return true;
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