#include <iostream>
#include <TreeDAWG.hpp>
#include <unordered_map>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>
#include "../../include/XBWT/XBWT.hpp"
#include "../../include/ChainsDivisionSolver/ChainsDivisionSolver.hpp"
#include "../../include/TreeCompressor/TreeCompressor.hpp"

int main()
{
    bool verbose = true, from_file = false;

    std::string filename = "tree_generator/generated_trees/tree_bf3_rp50_sd1-5_letters8_mn10000_s42.txt";
    std::string str = "(1(0(0(0)(1))(1))(1(0(0)(1))(1)))";

    if (from_file)
    { // Read input string from file
        std::ifstream inputFile(filename);
        if (inputFile.is_open())
        {
            getline(inputFile, str);
            inputFile.close();
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
            return 1;
        }
    }

    if (verbose)
        std::cout << "Input string: " << str << std::endl;

    std::cout << "Parsing tree..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    LabeledTree<char> tree(str, [](const std::string &label) -> char
                           { return label[0]; });
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Tree construction time: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " ms" << std::endl;

    std::cout << "Number of nodes: " << tree.get_nodes().size() << std::endl;

    if (verbose)
        std::cout << "Tree: " << tree.to_string() << std::endl;

    // test albero
    std::cout << "Building Tree DAWG..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    TreeDAWG<char> tree_dawg(tree);
    end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Tree DAWG construction completed in " << duration.count() << " ms" << std::endl;

    // Run minimization
    std::cout << "Minimizing Tree DAWG..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    tree_dawg.minimize(); // Assuming it's a tree
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Minimization completed in " << duration.count() << " ms" << std::endl;

    // Print DAWG information using to_string method
    if (verbose)
        std::cout << tree_dawg.to_string() << std::endl;

    // compute node order
    std::cout << "Computing node order..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<uint64_t> node_order = tree_dawg.stable_sort_nodes_by_label_path(verbose);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Node order computed in " << duration.count() << " ms" << std::endl;

    uint16_t p = 2;

    std::cout << "Building bipartite graph..." << std::endl;
    auto start_construction = std::chrono::high_resolution_clock::now();
    ChainsDivisionSolver solver(tree_dawg, node_order, p, verbose);
    auto end_construction = std::chrono::high_resolution_clock::now();
    auto construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_construction - start_construction);
    std::cout << "Construction time: " << construction_time.count() << " ms" << std::endl;

    std::cout << "Solving minimum weight perfect matching..." << std::endl;
    auto start_solving = std::chrono::high_resolution_clock::now();
    auto chains = solver.solve(verbose);
    auto end_solving = std::chrono::high_resolution_clock::now();
    auto solving_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_solving - start_solving);
    std::cout << "Solving time: " << solving_time.count() << " ms" << std::endl;

    if (verbose)
    {
        uint64_t tot_cost = 0;
        for (uint64_t i = 0; i < p; i++)
        {
            std::cout << "Chain " << i << " classes: ";
            uint64_t prev_class = UINT64_MAX;
            for (uint64_t j : chains[i])
            {
                std::cout << tree_dawg[j].get_equivalence_class() << " ";
                if (prev_class != tree_dawg[j].get_equivalence_class())
                    ++tot_cost;

                prev_class = tree_dawg[j].get_equivalence_class();
            }

            std::cout << std::endl;
        }

        std::cout << "Total cost: " << tot_cost << std::endl;
    }

    std::cout << "Compressing tree..." << std::endl;
    TreeCompressor<char> compressor(tree_dawg, chains, verbose);

    std::cout << "--------------------" << std::endl;

    return 0;
}
