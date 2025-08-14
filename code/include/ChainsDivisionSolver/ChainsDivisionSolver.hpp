#ifndef CHAINS_DIVISION_SOLVER_HPP
#define CHAINS_DIVISION_SOLVER_HPP

#include "../DAWG/TreeDAWG.hpp"
#include "../XBWT/XBWT.hpp"
#include "PerfectMatching.h" // Blossom V header

class ChainsDivisionSolver
{
private:
    uint64_t m_num_nodes;
    PerfectMatching m_pm;
    uint16_t m_p;                                // number of chains
    std::vector<std::vector<uint64_t>> m_chains; // chains of nodes id
    std::vector<uint64_t> m_chain_weights;       // weights of chains
    std::vector<uint64_t> m_node_order;          // store node order for output

public:
    template <typename T>
    ChainsDivisionSolver(const TreeDAWG<T> &treeDAWG, const std::vector<uint64_t> &node_order, uint16_t p, bool verbose = false)
        : m_num_nodes(treeDAWG.get_num_nodes()), m_pm(2 * (treeDAWG.get_num_nodes() + p), p * p + treeDAWG.get_num_nodes() * (p + 1) + treeDAWG.get_num_nodes() * p), m_p(p), m_node_order(node_order)
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

    // Solve minimum weight perfect matching
    std::vector<std::pair<int, int>> solve(bool verbose = false)
    {
        m_pm.Solve();
        if (verbose)
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
                    if (i < m_p) {
                        left_label = "s_" + std::to_string(i + 1);
                    } else if (i < m_p + t) {
                        left_label = "u_" + std::to_string(i - m_p);
                    }
                    
                    // Determine right node label
                    if (j >= m_p + t && j < m_p + t + t) {
                        right_label = "v_" + std::to_string(j - m_p - t);
                    } else if (j >= m_p + t + t) {
                        right_label = "d_" + std::to_string(j - m_p - t - t + 1);
                    }
                    
                    std::cout << left_label << " -> " << right_label << std::endl;
                }
            }
        }

        return {{}};
    }
};

#endif // CHAINS_DIVISION_SOLVER_HPP