#ifndef CHAINS_DIVISION_SOLVER_HPP
#define CHAINS_DIVISION_SOLVER_HPP

#include "../DAWG/TreeDAWG.hpp"
#include "../XBWT/XBWT.hpp"
#include "PerfectMatching.h" // Blossom V header

/**
 * @brief Solver for dividing tree nodes into chains using minimum weight perfect matching
 * 
 * This class implements a chain division algorithm that partitions tree nodes into p chains
 * by solving a minimum weight perfect matching problem on a bipartite graph. The algorithm
 * aims to minimize the number of equivalence class changes across chain boundaries.
 * 
 * The bipartite graph construction follows these rules:
 * - Left side (V1): source nodes s_1...s_p and tree nodes u_0...u_{t-1}
 * - Right side (V2): tree nodes v_0...v_{t-1} and destination nodes d_1...d_p
 * - Edge weights: 0 for same equivalence class, 1 for different equivalence classes
 */
class ChainsDivisionSolver
{
private:
    uint64_t m_num_nodes;                        ///< Total number of tree nodes
    PerfectMatching m_pm;                        ///< Blossom V perfect matching solver instance
    uint16_t m_p;                                ///< Number of chains to create
    std::vector<std::vector<uint64_t>> m_chains; ///< Resulting chains of node indices
    std::vector<uint64_t> m_node_order;          ///< Input node order sequence

    /**
     * @brief Prints the minimum weight perfect matching results
     * 
     * Outputs the matching pairs between left and right nodes in the bipartite graph,
     * showing which source/tree nodes are matched to which tree/destination nodes.
     */
    void print_matching_results()
    {
        std::cout << "Minimum weight perfect matching:" << std::endl;
        uint64_t t = m_num_nodes;
        for (uint64_t i = 0; i < t + m_p; i++)
        {
            uint64_t j = m_pm.GetMatch(i);
            if (i < j)
            {
                std::string left_label, right_label;

                // Determine left node label
                if (i < m_p)
                {
                    left_label = "s_" + std::to_string(i + 1);
                }
                else if (i < m_p + t)
                {
                    left_label = "u_" + std::to_string(i - m_p);
                }

                // Determine right node label
                if (j >= m_p + t && j < m_p + t + t)
                {
                    right_label = "v_" + std::to_string(j - m_p - t);
                }
                else if (j >= m_p + t + t)
                {
                    right_label = "d_" + std::to_string(j - m_p - t - t + 1);
                }

                std::cout << left_label << " -> " << right_label << std::endl;
            }
        }
    }

public:
    /**
     * @brief Constructs the chain division solver and builds the bipartite graph
     * 
     * @tparam T Type of the tree node labels
     * @param treeDAWG The tree DAWG containing nodes with equivalence classes
     * @param node_order Ordered sequence of tree node indices to be divided into chains
     * @param p Number of chains to create (must be > 0)
     * @param verbose If true, prints detailed construction information
     * 
     * The constructor builds a bipartite graph where:
     * - Source nodes s_i connect to the first p distinct equivalence class positions
     * - Tree nodes u_i connect to subsequent nodes with different equivalence classes
     * - Special weight-0 edges preserve equivalence class continuity within chains
     */
    template <typename T>
    ChainsDivisionSolver(const TreeDAWG<T> &treeDAWG, const std::vector<uint64_t> &node_order, uint16_t p, bool verbose = false)
        : m_num_nodes(treeDAWG.get_num_nodes()), m_pm(2 * (treeDAWG.get_num_nodes() + p), p * p + treeDAWG.get_num_nodes() * (p + 1) + treeDAWG.get_num_nodes() * p), m_p(p), m_chains(p), m_node_order(node_order)
    {
        if (verbose)
        {
            std::cout << "Equivalence classes in order:" << std::endl;
            for (uint64_t i = 0; i < node_order.size(); ++i)
            {
                std::cout << "Node " << node_order[i] << ": class "
                          << treeDAWG[node_order[i]].get_equivalence_class() << std::endl;
            }
            std::cout << std::endl;
        }

        uint64_t t = treeDAWG.get_num_nodes();

        // Find positions where equivalence class changes in the node order
        std::vector<std::pair<uint64_t, uint64_t>> class_change_pos;
        uint64_t last_class = UINT64_MAX;
        for (uint64_t i = 0; i < node_order.size(); ++i)
        {
            uint64_t current_class = treeDAWG[node_order[i]].get_equivalence_class();
            if (last_class == UINT64_MAX || last_class != current_class)
            {
                last_class = current_class;
                class_change_pos.push_back({i, current_class});
            }
        }

        // Build bipartite graph edges according to Definition
        // 1. Source edges: Connect source nodes s_1, ..., s_p to first p distinct equivalence classes in V2
        for (uint64_t i = 0; i < p && i < class_change_pos.size(); ++i)
        {
            uint64_t source_node = i;                                 // s_i in V1
            uint64_t target_node = t + p + class_change_pos[i].first; // v_j in V2
            m_pm.AddEdge(source_node, target_node, 1);

            if (verbose)
                std::cout << "Add edge: s_" << i + 1 << " -> v_" << class_change_pos[i].first << " (weight 1)" << std::endl;
        }

        // 2. Tree node edges: For each u_i in V1 (tree nodes)
        uint64_t cp_index = 0;
        for (uint64_t i = 0; i < t; ++i)
        {
            uint64_t u_i = p + i; // u_i position in V1
            uint64_t current_equiv_class = treeDAWG[node_order[i]].get_equivalence_class();

            // Find first p nodes v_j in V2 such that j > i and equiv_class(v_j) != equiv_class(u_i)
            uint64_t edges_added = 0;
            bool found_zero_weight_edge = false;
            // Check if current node i is between two class change positions with his same class
            if (i < class_change_pos[cp_index].first - 1)
            {
                if (cp_index - 1 >= 0 &&
                    i >= class_change_pos[cp_index - 1].first &&
                    current_equiv_class == class_change_pos[cp_index - 1].second)
                {
                    // Add edge to current index + 1 with weight 0
                    uint64_t v_j = t + p + i + 1;
                    m_pm.AddEdge(u_i, v_j, 0);
                    found_zero_weight_edge = true;
                    if (verbose)
                        std::cout << "Add edge: u_" << i << " -> v_" << i + 1 << " (weight 0)" << std::endl;
                }
            }
            else if (i >= class_change_pos[cp_index].first)
            {
                if (cp_index + 1 < class_change_pos.size() &&
                    i < class_change_pos[cp_index + 1].first - 1 &&
                    current_equiv_class == class_change_pos[cp_index].second)
                {
                    // Add edge to current index + 1 with weight 0
                    uint64_t v_j = t + p + i + 1;
                    m_pm.AddEdge(u_i, v_j, 0);
                    found_zero_weight_edge = true;
                    if (verbose)
                        std::cout << "Add edge: u_" << i << " -> v_" << i + 1 << " (weight 0)" << std::endl;
                }
            }

            if (class_change_pos[cp_index].first <= i)
                ++cp_index;

            uint64_t cur_cp_index = cp_index;
            while (cur_cp_index < class_change_pos.size() && edges_added <= p)
            {
                uint64_t j = class_change_pos[cur_cp_index].first;
                uint64_t target_equiv_class = class_change_pos[cur_cp_index].second;
                uint64_t v_j = t + p + j; // v_j position in V2

                if (target_equiv_class != current_equiv_class && edges_added < p)
                {
                    m_pm.AddEdge(u_i, v_j, 1);
                    if (verbose)
                        std::cout << "Add edge: u_" << i << " -> v_" << j << " (weight 1)" << std::endl;

                    edges_added++;
                }
                else if (target_equiv_class == current_equiv_class && !found_zero_weight_edge)
                {
                    m_pm.AddEdge(u_i, v_j, 0);
                    found_zero_weight_edge = true;
                    if (verbose)
                        std::cout << "Add edge: u_" << i << " -> v_" << j << " (weight 0)" << std::endl;
                }

                ++cur_cp_index;
            }

            // If no nodes with different equivalence class found, connect to destinations
            if (!found_zero_weight_edge)
            {
                for (uint64_t j = 0; j < p; ++j)
                {
                    uint64_t d_j = t + p + t + j; // d_j position in V2 (destinations)
                    m_pm.AddEdge(u_i, d_j, 0);
                    if (verbose)
                        std::cout << "Add edge: u_" << i << " -> d_" << j + 1 << " (weight 0, destination)" << std::endl;
                }
            }
        }
    }

    /**
     * @brief Solves the minimum weight perfect matching and constructs the chains
     * 
     * @param verbose If true, prints matching results and final chains
     * @return Vector of p chains, where each chain contains node indices in order
     * 
     * Uses the Blossom V algorithm to find the minimum weight perfect matching,
     * then reconstructs the chains by following the matching edges from source
     * nodes through the tree nodes.
     */
    std::vector<std::vector<uint64_t>> solve(bool verbose = false)
    {
        m_pm.Solve();
        if (verbose)
        {
            print_matching_results();
        }

        std::vector<uint64_t> chain_mem(m_p, 0); // array used to save cur chain node
        for (uint64_t i = 0; i < m_num_nodes; i++)
        {
            uint64_t j = m_pm.GetMatch(i);
            uint64_t real_i = i - m_p;
            uint64_t real_j = j - m_num_nodes - m_p;
            // sources
            if (i < m_p)
            {
                chain_mem[i] = real_j;
                m_chains[i].push_back(m_node_order[real_j]);
            }
            else if (j < m_p + m_num_nodes * 2)
            {
                for (uint64_t k = 0; k < m_p; k++)
                {
                    if (chain_mem[k] == real_i)
                    {
                        m_chains[k].push_back(m_node_order[real_j]);
                        chain_mem[k] = real_j;
                        break;
                    }
                }
            }
        }

        if (verbose)
        {
            for (uint64_t i = 0; i < m_p; i++)
            {
                std::cout << "Chain " << i << ": ";
                for (uint64_t j : m_chains[i])
                {
                    std::cout << j << " ";
                }
                std::cout << std::endl;
            }
        }

        return m_chains;
    }

    /**
     * @brief Returns the constructed chains
     * 
     * @return Const reference to the vector of chains
     * 
     * Note: This method should only be called after solve() has been executed,
     * otherwise it returns empty chains.
     */
    const std::vector<std::vector<uint64_t>>& get_chains() const {
        return m_chains;
    }
};

#endif // CHAINS_DIVISION_SOLVER_HPP