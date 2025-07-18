#ifndef DAWG_HPP
#define DAWG_HPP

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <list>
#include <functional>

// Special values used in sorting to represent the absence of a transition
// or a transition to a "dead state".
// These values must be outside the range of normal equivalence classes (which are >= 0).
const int NIL_TRANSITION_SORT_VAL = -1;
const int DEAD_TARGET_SORT_VAL = -2;

class DAWG;

/**
 * @brief Represents a single state in the DAWG.
 * 
 * Each state has a unique ID, a flag indicating if it's a final state,
 * and a map of transitions to other states.
 */
struct State
{
public:
    size_t id;
    bool is_final;
    std::map<char, size_t> transitions; // Transitions: character -> destination node id

    /**
     * @brief Get the height of the state.
     * The height is the length of the longest path to a final state.
     * @return The height of the state.
     */
    long long get_height() const { return m_height; }

    /**
     * @brief Get the equivalence class of the state.
     * This is used during the minimization process.
     * @return The equivalence class ID.
     */
    size_t get_equivalence_class() const { return m_equivalence_class; }

    /**
     * @brief Construct a new State object.
     * 
     * @param _id The unique identifier for the state.
     * @param _is_final Whether the state is a final state.
     */
    State(size_t _id, bool _is_final = false)
        : id(_id), is_final(_is_final), m_height(-1), m_equivalence_class(-1) {}

private:
    long long m_height;            // Cached height of the state.
    size_t m_equivalence_class; // Equivalence class ID used in minimization.

    friend class DAWG;
};

/**
 * @brief A class to represent a Deterministic Acyclic Word Graph (DAWG).
 */
class DAWG
{
public:
    std::vector<State> nodes;

private:
    size_t m_initial_state_id;
    bool m_is_acyclic;

    /**
     * @brief Recursively calculates the height of a node using DFS.
     * Also detects cycles in the graph.
     * @param node_id The ID of the node to calculate the height for.
     * @param visited_path_flags A vector to keep track of visited nodes to detect cycles.
     * @return The height of the node, or a special value indicating a cycle or non-final path.
     */
    long long calculate_height_dfs(size_t node_id, std::vector<int> &visited_path_flags)
    {
        State &node = nodes[node_id];
        if (node.m_height != -1)
            return node.m_height;
        visited_path_flags[node_id] = 1; // Mark as visiting

        long long current_max_h = -2; // Represents a non-final state with no path to a final state
        if (node.is_final)
            current_max_h = 0;

        for (const auto &[_, target_node_id] : node.transitions)
        {
            if (target_node_id >= nodes.size())
            {
                std::cerr << "Error: transition from node " << node_id << " to non-existent node " << target_node_id << std::endl;
                continue;
            }
            if (visited_path_flags[target_node_id] == 1)
            {
                m_is_acyclic = false;
                return -3; // Cycle detected
            }
            long long target_height = calculate_height_dfs(target_node_id, visited_path_flags);
            if (!m_is_acyclic)
                return -3;
            if (target_height != -2 && target_height != -3)
            {
                current_max_h = std::max(current_max_h, 1 + target_height);
            }
        }
        visited_path_flags[node_id] = 2; // Mark as visited
        node.m_height = current_max_h;
        return node.m_height;
    }

    /**
     * @brief Performs a counting sort on a block of nodes.
     * 
     * @param block The list of node IDs to sort.
     * @param get_value A function to extract the sorting value from a node ID.
     * @param min_val The minimum possible value for get_value.
     * @param max_val The maximum possible value for get_value.
     * @return A list of new blocks, where each block contains nodes with the same value.
     */
    std::list<std::list<size_t>> counting_sort_block(
        const std::list<size_t> &block,
        std::function<int(size_t)> get_value, // Function to get the sorting value for a nodeId
        int min_val,                          // Minimum possible value for getValue
        int max_val)                          // Maximum possible value for getValue
    {
        std::list<std::list<size_t>> new_blocks;
        if (block.empty())
            return new_blocks;
        if (min_val > max_val)
        { // Can happen if the block is empty or has only unexpected special values
            if (!block.empty())
                new_blocks.push_back(block); // Do not split if the range is not valid
            return new_blocks;
        }

        int range = max_val - min_val + 1;
        std::vector<std::vector<size_t>> buckets(range);

        for (size_t node_id : block)
        {
            int sort_value = get_value(node_id);
            // Ensure the sort value is within the expected range for the buckets
            if (sort_value < min_val || sort_value > max_val)
            {
                // This would indicate a problem with the min_val/max_val logic
                // or with the values returned by getValue. For robustness, one could
                // have an "overflow" bucket or handle the error.
                // For now, we assume that getValue returns values in the range [min_val, max_val].
                // Or, more simply, if it happens, we print an error and do not split.
                std::cerr << "Error: sort_value " << sort_value << " out of range [" << min_val << ", " << max_val << "]" << std::endl;
            }
            buckets[sort_value - min_val].push_back(node_id);
        }

        for (const auto &bucket_content : buckets)
        {
            if (!bucket_content.empty())
            {
                new_blocks.emplace_back(bucket_content.begin(), bucket_content.end());
            }
        }
        return new_blocks;
    }

public:
    /**
     * @brief Default constructor for the DAWG class.
     */
    DAWG() : m_initial_state_id(0), m_is_acyclic(true) {}

    /**
     * @brief Constructor for the DAWG class with a specified initial state.
     * @param set_initial_state_id The ID of the initial state.
     */
    DAWG(size_t set_initial_state_id) : m_initial_state_id(set_initial_state_id), m_is_acyclic(true) {}

    /**
     * @brief Adds a new node to the DAWG.
     * @param is_final Whether the new node is a final state.
     * @return The ID of the newly created node.
     */
    size_t add_node(bool is_final = false)
    {
        size_t id = nodes.size();
        nodes.emplace_back(id, is_final);
        return id;
    }

    /**
     * @brief Adds a transition between two nodes.
     * @param from_node_id The ID of the source node.
     * @param symbol The transition symbol.
     * @param to_node_id The ID of the target node.
     */
    void add_transition(size_t from_node_id, char symbol, size_t to_node_id)
    {
        if (from_node_id < nodes.size() && to_node_id < nodes.size())
        {
            nodes[from_node_id].transitions[symbol] = to_node_id;
        }
        else
        {
            std::cerr << "Error: transition between non-existent nodes." << std::endl;
        }
    }

    /**
     * @brief Sets the initial state of the DAWG.
     * @param node_id The ID of the node to be set as the initial state.
     */
    void set_initial_state(size_t node_id)
    {
        if (node_id < nodes.size())
        {
            m_initial_state_id = node_id;
        }
        else
        {
            std::cerr << "Error: attempt to set a non-existent initial state." << std::endl;
        }
    }

    /**
     * @brief Computes the height of all nodes in the DAWG.
     * The height is the length of the longest path from a node to any final state.
     * This version handles general DAGs and includes cycle detection.
     * @return The maximum height found among all nodes, or -3 if a cycle is detected.
     */
    long long compute_all_heights()
    {
        m_is_acyclic = true;
        for (auto &node : nodes)
            node.m_height = -1;
        long long max_h = -1;
        std::vector<int> visited_path_flags(nodes.size(), 0); // 0: unvisited, 1: visiting, 2: visited
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].m_height == -1)
            {
                long long h_node = calculate_height_dfs(i, visited_path_flags);
                if (!m_is_acyclic)
                    return -3; // Cycle detected
                max_h = std::max(max_h, h_node);
            }
            else
            {
                max_h = std::max(max_h, nodes[i].m_height);
            }
        }
        return max_h;
    }

    /**
     * @brief Computes the height of all nodes for a DAWG that is a tree.
     * This version is optimized for trees and does not perform cycle detection.
     * It uses an iterative post-order traversal to avoid stack overflow on deep trees.
     * @return The maximum height of any node in the tree.
     */
    long long compute_all_heights_for_tree()
    {
        for (auto &node : nodes)
            node.m_height = -1;

        long long max_h = -1;
        // compute height starting from root
        std::stack<size_t> dfs_stack;
        std::stack<size_t> processing_stack; // post-order stack

        dfs_stack.push(m_initial_state_id);

        std::vector<bool> visited(nodes.size(), false);
        while (!dfs_stack.empty())
        {
            size_t current_node_id = dfs_stack.top();
            dfs_stack.pop();

            if (visited[current_node_id])
                continue;

            visited[current_node_id] = true;
            processing_stack.push(current_node_id);

            for (const auto &[_, target_node_id] : nodes[current_node_id].transitions)
            {
                if (!visited[target_node_id])
                {
                    dfs_stack.push(target_node_id);
                }
            }
        }

        while (!processing_stack.empty())
        {
            size_t node_id = processing_stack.top();
            processing_stack.pop();

            State &node = nodes[node_id];

            if (node.m_height != -1)
                continue;

            long long current_max_h = node.is_final ? 0 : -2;
            for (const auto &[_, target_node_id] : node.transitions)
            {
                long long target_height = nodes[target_node_id].m_height;

                if (target_height != -2)
                {
                    current_max_h = std::max(current_max_h, 1 + target_height);
                }
            }

            node.m_height = current_max_h;
        }

        // find max height
        for (const auto &node : nodes)
        {
            max_h = std::max(max_h, node.m_height);
        }

        return max_h;
    }

    /**
     * @brief Minimizes the DAWG using a partition refinement algorithm.
     * This implementation is based on Revuz's algorithm for DAWGs minimization (DOI: 10.1016/0304-3975(92)90142-3).
     * It groups states into equivalence classes.
     * @param is_tree A boolean indicating if the DAWG is a tree, to use the optimized height calculation.
     * @return An unordered_map from original state ID to the representative state ID of its equivalence class.
     */
    std::unordered_map<size_t, size_t> minimize(bool is_tree = true)
    {
        long long max_total_height = is_tree ? compute_all_heights_for_tree() : compute_all_heights();
        if (!m_is_acyclic)
            return {};

        std::vector<std::vector<size_t>> states_by_height(max_total_height + 1);
        for (const auto &node : nodes)
        {
            if (node.m_height >= 0)
            {
                states_by_height[node.m_height].push_back(node.id);
            }
        }

        size_t next_global_eq_class_id = 0;

        if (max_total_height >= 0 && !states_by_height[0].empty())
        {
            for (size_t node_id : states_by_height[0])
            {
                nodes[node_id].m_equivalence_class = next_global_eq_class_id;
            }
            next_global_eq_class_id++;
        }

        // Cache for sorted transition vectors (for efficient k-th access)
        // for nodes in the current height level.
        std::vector<std::vector<std::pair<char, size_t>>> node_sorted_transitions_cache(nodes.size());

        for (size_t h = 1; h <= static_cast<size_t>(max_total_height); ++h)
        {
            if (states_by_height[h].empty())
                continue;

            for (size_t node_id : states_by_height[h])
            {
                node_sorted_transitions_cache[node_id].assign(
                    nodes[node_id].transitions.begin(),
                    nodes[node_id].transitions.end());
            }

            std::list<std::list<size_t>> h_partitions;
            if (!states_by_height[h].empty())
            {
                h_partitions.push_back({});
                for (size_t id : states_by_height[h])
                    h_partitions.back().push_back(id);
            }

            // 1. Partition by is_final
            std::list<std::list<size_t>> p_after_is_final;
            for (auto &block : h_partitions)
            {
                if (block.size() <= 1)
                {
                    p_after_is_final.push_back(block);
                    continue;
                }
                auto get_is_final_val = [&](size_t node_id)
                { return nodes[node_id].is_final ? 1 : 0; };
                auto new_sub_blocks = counting_sort_block(block, get_is_final_val, 0, 1); // min=0, max=1 for bool
                for (auto &sub_block : new_sub_blocks)
                    p_after_is_final.push_back(sub_block);
            }
            h_partitions = p_after_is_final;

            size_t max_k_transitions_for_h = 0;
            for (size_t node_id : states_by_height[h])
            {
                max_k_transitions_for_h = std::max(max_k_transitions_for_h, node_sorted_transitions_cache[node_id].size());
            }

            for (size_t k = 0; k < max_k_transitions_for_h; ++k)
            { // For the k-th transition
                // 2a. Partition by the character of the k-th transition
                std::list<std::list<size_t>> p_after_char;
                int min_val_char = NIL_TRANSITION_SORT_VAL; // -1
                int max_val_char = 255;                     // Max value for unsigned char

                for (auto &block : h_partitions)
                {
                    if (block.size() <= 1)
                    {
                        p_after_char.push_back(block);
                        continue;
                    }
                    auto get_kth_char_val = [&](size_t node_id)
                    {
                        if (k < node_sorted_transitions_cache[node_id].size())
                        {
                            return static_cast<int>(static_cast<unsigned char>(node_sorted_transitions_cache[node_id][k].first));
                        }
                        return NIL_TRANSITION_SORT_VAL;
                    };
                    
                    // Determine the effective range of values in the current block to optimize countingSort
                    int current_block_min_char = max_val_char + 1;
                    int current_block_max_char = min_val_char - 1;
                    bool has_values_char = false;
                    for (size_t node_id : block)
                    {
                        int val = get_kth_char_val(node_id);
                        current_block_min_char = std::min(current_block_min_char, val);
                        current_block_max_char = std::max(current_block_max_char, val);
                        has_values_char = true;
                    }
                    if (!has_values_char)
                    {
                        p_after_char.push_back(block);
                        continue;
                    }

                    auto new_sub_blocks = counting_sort_block(block, get_kth_char_val, current_block_min_char, current_block_max_char);
                    for (auto &sub_block : new_sub_blocks)
                        p_after_char.push_back(sub_block);
                }
                h_partitions = p_after_char;

                // 2b. Partition by the equivalence class of the k-th transition target
                std::list<std::list<size_t>> p_after_target_class;
                // Theoretical range for target classes: DEAD_TARGET_SORT_VAL (-2), NIL_TRANSITION_SORT_VAL (-1), 0 ... next_global_eq_class_id-1
                int min_val_target_overall = DEAD_TARGET_SORT_VAL;
                int max_val_target_overall = (next_global_eq_class_id == 0) ? DEAD_TARGET_SORT_VAL : static_cast<int>(next_global_eq_class_id) - 1;
                if (next_global_eq_class_id == 0)
                { // If Pi_0 was empty and no classes were created
                    max_val_target_overall = std::max(NIL_TRANSITION_SORT_VAL, DEAD_TARGET_SORT_VAL);
                }

                for (auto &block : h_partitions)
                {
                    if (block.size() <= 1)
                    {
                        p_after_target_class.push_back(block);
                        continue;
                    }

                    auto get_kth_target_eq_class_val = [&](size_t node_id)
                    {
                        if (k < node_sorted_transitions_cache[node_id].size())
                        {
                            size_t target_node_id = node_sorted_transitions_cache[node_id][k].second;
                            if (nodes[target_node_id].m_height < 0)
                            { // Target is a "dead state"
                                return DEAD_TARGET_SORT_VAL;
                            }
                            return static_cast<int>(nodes[target_node_id].m_equivalence_class); // Class already finalized
                        }
                        return NIL_TRANSITION_SORT_VAL; // No k-th transition
                    };

                    int current_block_min_target_class = max_val_target_overall + 1; // Initialize to find the effective min
                    int current_block_max_target_class = min_val_target_overall - 1; // Initialize to find the effective max
                    bool has_values_target = false;
                    for (size_t node_id : block)
                    {
                        int val = get_kth_target_eq_class_val(node_id);
                        current_block_min_target_class = std::min(current_block_min_target_class, val);
                        current_block_max_target_class = std::max(current_block_max_target_class, val);
                        has_values_target = true;
                    }
                    if (!has_values_target)
                    {
                        p_after_target_class.push_back(block);
                        continue;
                    }

                    auto new_sub_blocks = counting_sort_block(block, get_kth_target_eq_class_val, current_block_min_target_class, current_block_max_target_class);
                    for (auto &sub_block : new_sub_blocks)
                        p_after_target_class.push_back(sub_block);
                }
                h_partitions = p_after_target_class;
            }

            for (auto &final_block : h_partitions)
            {
                if (!final_block.empty())
                {
                    for (size_t node_id : final_block)
                    {
                        nodes[node_id].m_equivalence_class = next_global_eq_class_id;
                    }
                    next_global_eq_class_id++;
                }
            }
        }

        std::unordered_map<size_t, size_t> result;
        for (const auto &node : nodes)
        {
            if (node.m_height >= 0)
            {
                result[node.id] = node.m_equivalence_class;
            }
            else
            {
                result[node.id] = -1; // Special class for non-useful nodes
            }
        }
        return result;
    }
};

#endif // DAWG_HPP
