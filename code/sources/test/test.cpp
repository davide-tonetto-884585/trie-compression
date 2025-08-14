#include <iostream>
#include <TreeDAWG.hpp>
#include <unordered_map>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>
#include "../../include/XBWT/XBWT.hpp"
#include "../../include/ChainsDivisionSolver/ChainsDivisionSolver.hpp"

// Function to convert tree structure to balanced parentheses
std::string tree_to_balanced_parentheses(const TreeDAWG<char> &tree_dawg)
{
    std::string result;

    // Helper function for DFS traversal
    std::function<void(size_t, char)> dfs = [&](size_t node_id, char label)
    {
        result += '(';
        result += label;

        // Get transitions and sort them for consistent output
        const auto &transitions = tree_dawg[node_id].get_transitions();
        for (const auto &[transition_label, target_id] : transitions)
        {
            dfs(target_id, transition_label);
        }

        result += ')';
    };

    // Start DFS from the initial state with label 'A'
    dfs(tree_dawg.get_initial_state_id(), 'A');
    return result;
}

int main()
{
    bool verbose = false;

    std::string filename = "tree_generator/generated_trees/tree_bf3_rp50_sd1-5_letters8_mn100000_s42.txt";
    std::string str = "(1(0(0(0)(1))(1))(1(0(0)(1))(1)))";

    // Read input string from file
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

    // build an XBWT
    std::vector<unsigned int> intNodesPosSorted{};
    std::vector<unsigned int> testIntNodesPosSorted{};

    // Measure time for XBWT construction
    std::cout << "Building XBWT..." << std::endl;
    /* start = std::chrono::high_resolution_clock::now();
    XBWT<char> xbwt(tree, true, verbose, &intNodesPosSorted);
    end = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "XBWT construction time: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " ms" << std::endl; */

    start = std::chrono::high_resolution_clock::now();
    XBWT<char> xbwt2(tree, false, verbose, &testIntNodesPosSorted);
    end = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "XBWT construction time no path sort: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " ms" << std::endl;

    /* if (intNodesPosSorted != testIntNodesPosSorted)
        throw std::runtime_error("IntNodesPosSorted is not correct"); */

    // test albero
    std::cout << "Building Tree DAWG..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    TreeDAWG<char> tree_dawg(tree);
    end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Tree DAWG construction completed in " << duration.count() << " ms" << std::endl;

    // Print as balanced parentheses
    if (verbose)
    {
        std::string balanced_parentheses = tree_to_balanced_parentheses(tree_dawg);
        std::cout << std::endl
                  << "Tree as balanced parentheses: " << balanced_parentheses << std::endl;
    }

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

    uint16_t p = 2;
    std::vector<uint64_t> node_order(testIntNodesPosSorted.begin(), testIntNodesPosSorted.end());
    
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

    std::cout << "--------------------" << std::endl;

    return 0;
}
