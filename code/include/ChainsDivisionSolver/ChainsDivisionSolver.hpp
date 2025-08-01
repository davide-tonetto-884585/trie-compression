#ifndef CHAINS_DIVISION_SOLVER_HPP
#define CHAINS_DIVISION_SOLVER_HPP

#include "../DAWG/TreeDAWG.hpp"
#include "../XBWT/XBWT.hpp"
#include "PerfectMatching.h" // Blossom V header

class ChainsDivisionSolver
{
private:
    PerfectMatching m_pm;
    uint16_t m_p;                                // number of chains
    std::vector<std::vector<uint64_t>> m_chains; // chains of nodes id
    std::vector<uint64_t> m_chain_weights;       // weights of chains

public:
    ChainsDivisionSolver(const TreeDAWG<> &treeDAWG, uint16_t p)
    {
        // compute number of nodes for bipartite graph
        uint64_t num_nodes = treeDAWG.get_num_nodes() + p;
        // compute max number of edges for bipartite graph
        uint64_t max_edges = num_nodes * (p + 1) + p * p + t * p;
        // initialize perfect matching
        m_pm = PerfectMatching(num_nodes, max_edges);

        // order tree nodes using xbwt pathsort and get equivalence classes ordered
        XBWT xbwt = XBWT(treeDAWG);
        std::vector<uint64_t> order = xbwt.get_node_order();
        std::vector<std::pair<uint64_t, uint64_t>> class_change_pos;
        uint64_t last_class = -1;
        for (uint64_t i = 0; i < order.size(); ++i)
        {
            if (last_class == -1)
            {
                last_class = treeDAWG[order[i]].get_equivalence_class();
                class_change_pos.push_back({order[i], last_class});
            }
            else if (last_class != treeDAWG[order[i]].get_equivalence_class())
            {
                last_class = treeDAWG[order[i]].get_equivalence_class();
                class_change_pos.push_back({order[i], last_class});
            }
        }

        for (size_t i = 0; i < num_nodes; ++i)
        {
            // source edges
            if (i < p)
            {
                m_pm.AddEdge(i, treeDAWG.get_num_nodes() + i, 1);
            }
        }
    }

    // Solve minimum weight perfect matching
    std::vector<std::pair<int, int>> solve()
    {
    }
};

#endif // CHAINS_DIVISION_SOLVER_HPP