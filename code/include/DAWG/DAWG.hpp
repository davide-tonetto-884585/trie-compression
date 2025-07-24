#ifndef DAWG_HPP
#define DAWG_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <stack>
#include <stdexcept>
#include <string>
#include <cstdint>

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
    class Builder
    {
    private:
        std::vector<std::pair<char, size_t>> pending_transitions;

    public:
        Builder &add_transition(char symbol, size_t target_id)
        {
            pending_transitions.emplace_back(symbol, target_id);
            return *this;
        }

        void build_into(State &state)
        {
            state.m_transitions = std::move(pending_transitions);
            // Sort once at the end for efficiency
            std::sort(
                state.m_transitions.begin(),
                state.m_transitions.end(),
                [](const auto &a, const auto &b)
                {
                    return a.first < b.first;
                });
        }
    };

    // Helper method to find transition (replaces map's operator[])
    size_t *find_transition(char symbol)
    {
        auto it = std::lower_bound(
            m_transitions.begin(),
            m_transitions.end(),
            std::make_pair(symbol, size_t(0)),
            [](const auto &a, const auto &b)
            {
                return a.first < b.first;
            }
        );
            
        if (it != m_transitions.end() && it->first == symbol)
        {
            return &(it->second);
        }

        return nullptr;
    }

    // Getter methods for private fields
    size_t get_id() const { return m_id; }
    bool is_final() const { return m_transitions.empty(); }
    const std::vector<std::pair<char, size_t>>& get_transitions() const { return m_transitions; }
    std::vector<std::pair<char, size_t>>& get_transitions() { return m_transitions; }

    /**
     * @brief Get the height of the state.
     * The height is the length of the longest path to a final state.
     * @return The height of the state.
     */
    int64_t get_height() const { return m_height; }

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
    State(size_t _id)
        : m_id(_id), m_height(-1), m_equivalence_class(-1) {}

private:
    size_t m_id;
    std::vector<std::pair<char, size_t>> m_transitions;
    
    int64_t m_height;         // Cached height of the state.
    size_t m_equivalence_class; // Equivalence class ID used in minimization.

    friend class DAWG;
    friend class TreeDAWG;
};

/**
 * @brief A class to represent a Deterministic Acyclic Word Graph (DAWG).
 */
class DAWG
{
protected:
    std::vector<State> m_nodes;
    size_t m_initial_state_id;

private:
    /**
     * @brief Recursively calculates the height of a node using DFS.
     * Also detects cycles in the graph.
     * @param node_id The ID of the node to calculate the height for.
     * @param visited_path_flags A vector to keep track of visited nodes to detect cycles.
     * @return The height of the node, or a special value indicating a cycle or non-final path.
     */
    int64_t calculate_height_dfs(size_t node_id, std::vector<int> &visited_path_flags)
    {
        State &node = m_nodes[node_id];
        if (node.m_height != -1)
            return node.m_height;
        visited_path_flags[node_id] = 1; // Mark as visiting

        int64_t current_max_h = -2; // Represents a non-final state with no path to a final state
        if (node.is_final())
            current_max_h = 0;

        for (const auto &[_, target_node_id] : node.m_transitions)
        {
            if (target_node_id >= m_nodes.size())
            {
                throw std::out_of_range("Transition from node " + std::to_string(node_id) + " to non-existent node " + std::to_string(target_node_id));
            }
            if (visited_path_flags[target_node_id] == 1)
            {
                throw std::runtime_error("Cycle detected in DAWG at node " + std::to_string(target_node_id));
            }
            int64_t target_height = calculate_height_dfs(target_node_id, visited_path_flags);
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
    std::vector<std::vector<size_t>> counting_sort_block(
        const std::vector<size_t> &block,
        std::function<int(size_t)> get_value, // Function to get the sorting value for a nodeId
        int min_val,                          // Minimum possible value for getValue
        int max_val)                          // Maximum possible value for getValue
    {
        std::vector<std::vector<size_t>> new_blocks;
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
                throw std::out_of_range("Sort value " + std::to_string(node_id) + " out of range [" + std::to_string(min_val) + ", " + std::to_string(max_val) + "]");
            }
            buckets[sort_value - min_val].push_back(node_id);
        }

        for (const auto &bucket_content : buckets)
        {
            if (!bucket_content.empty())
            {
                new_blocks.push_back(bucket_content);
            }
        }
        return new_blocks;
    }

public:
    /**
     * @brief Default constructor for the DAWG class.
     */
    DAWG() : m_initial_state_id(0) {}

    /**
     * @brief Constructor for the DAWG class with a specified initial state.
     * @param set_initial_state_id The ID of the initial state.
     */
    DAWG(size_t set_initial_state_id) : m_initial_state_id(set_initial_state_id) {}

    State& operator[](size_t id) {
        if (id >= m_nodes.size()) {
            throw std::out_of_range("State ID " + std::to_string(id) + " is out of range. Maximum ID: " + std::to_string(m_nodes.size() - 1));
        }
        return m_nodes[id];
    }
    
    const State& operator[](size_t id) const {
        if (id >= m_nodes.size()) {
            throw std::out_of_range("State ID " + std::to_string(id) + " is out of range. Maximum ID: " + std::to_string(m_nodes.size() - 1));
        }
        return m_nodes[id];
    }
    
    size_t get_num_nodes() const {
        return m_nodes.size();
    }

    /**
     * @brief Adds a new node to the DAWG.
     * @param is_final Whether the new node is a final state.
     * @return The ID of the newly created node.
     */
    size_t add_node()
    {
        size_t id = m_nodes.size();
        m_nodes.emplace_back(id);
        return id;
    }

    void configure_state(size_t state_id, std::initializer_list<std::pair<char, size_t>> transitions)
    {
        State::Builder builder;
        for (const auto &[symbol, target] : transitions)
        {
            builder.add_transition(symbol, target);
        }

        builder.build_into(m_nodes[state_id]);
    }

    /**
     * @brief Sets the initial state of the DAWG.
     * @param node_id The ID of the node to be set as the initial state.
     */
    virtual void set_initial_state(size_t node_id)
    {
        if (node_id < m_nodes.size())
        {
            m_initial_state_id = node_id;
        }
        else
        {
            throw std::out_of_range("Attempt to set non-existent initial state with ID " + std::to_string(node_id));
        }
    }

    /**
     * @brief Computes the height of all nodes in the DAWG.
     * The height is the length of the longest path from a node to any final state.
     * This version handles general DAGs and includes cycle detection.
     * @return The maximum height found among all nodes, or -3 if a cycle is detected.
     */
    virtual int64_t compute_all_heights()
    {
        for (auto &node : m_nodes)
            node.m_height = -1;
        int64_t max_h = -1;
        std::vector<int> visited_path_flags(m_nodes.size(), 0); // 0: unvisited, 1: visiting, 2: visited
        for (size_t i = 0; i < m_nodes.size(); ++i)
        {
            if (m_nodes[i].m_height == -1)
            {
                int64_t h_node = calculate_height_dfs(i, visited_path_flags);
                max_h = std::max(max_h, h_node);
            }
            else
            {
                max_h = std::max(max_h, m_nodes[i].m_height);
            }
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
    std::unordered_map<size_t, size_t> minimize()
    {
        int64_t max_total_height = compute_all_heights();

        std::vector<std::vector<size_t>> states_by_height(max_total_height + 1);
        for (const auto &node : m_nodes)
        {
            if (node.m_height >= 0)
            {
                states_by_height[node.m_height].push_back(node.m_id);
            }
        }

        size_t next_global_eq_class_id = 0;

        if (max_total_height >= 0 && !states_by_height[0].empty())
        {
            for (size_t node_id : states_by_height[0])
            {
                m_nodes[node_id].m_equivalence_class = next_global_eq_class_id;
            }
            next_global_eq_class_id++;
        }

        // Cache for sorted transition vectors (for efficient k-th access)
        // for nodes in the current height level.
        std::vector<std::vector<std::pair<char, size_t>>> node_sorted_transitions_cache(m_nodes.size());

        for (size_t h = 1; h <= static_cast<size_t>(max_total_height); ++h)
        {
            if (states_by_height[h].empty())
                continue;

            for (size_t node_id : states_by_height[h])
            {
                node_sorted_transitions_cache[node_id].assign(
                    m_nodes[node_id].m_transitions.begin(),
                    m_nodes[node_id].m_transitions.end());
            }

            std::vector<std::vector<size_t>> h_partitions;
            if (!states_by_height[h].empty())
            {
                h_partitions.push_back(states_by_height[h]);
            }

            // 1. Partition by is_final
            std::vector<std::vector<size_t>> p_after_is_final;
            for (auto &block : h_partitions)
            {
                if (block.size() <= 1)
                {
                    p_after_is_final.push_back(block);
                    continue;
                }
                auto get_is_final_val = [&](size_t node_id)
                { return m_nodes[node_id].is_final() ? 1 : 0; };
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
                std::vector<std::vector<size_t>> p_after_char;
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
                std::vector<std::vector<size_t>> p_after_target_class;
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
                            if (m_nodes[target_node_id].m_height < 0)
                            { // Target is a "dead state"
                                return DEAD_TARGET_SORT_VAL;
                            }
                            return static_cast<int>(m_nodes[target_node_id].m_equivalence_class); // Class already finalized
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
                        m_nodes[node_id].m_equivalence_class = next_global_eq_class_id;
                    }
                    next_global_eq_class_id++;
                }
            }
        }

        std::unordered_map<size_t, size_t> result;
        for (const auto &node : m_nodes)
        {
            if (node.m_height >= 0)
            {
                result[node.m_id] = node.m_equivalence_class;
            }
            else
            {
                result[node.m_id] = -1; // Special class for non-useful nodes
            }
        }
        return result;
    }
};

#endif // DAWG_HPP
