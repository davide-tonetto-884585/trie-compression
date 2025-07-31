#include <iostream>
#include <TreeDAWG.hpp>
#include <unordered_map>
#include <string>
#include <fstream>
#include <chrono>

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
    // test albero
    TreeDAWG<char> tree_dawg;

    // Node mapping from character labels to integer IDs
    std::unordered_map<char, size_t> char_to_id;
    char_to_id['a'] = tree_dawg.add_node();
    char_to_id['b'] = tree_dawg.add_node();
    char_to_id['c'] = tree_dawg.add_node();
    char_to_id['d'] = tree_dawg.add_node();
    char_to_id['e'] = tree_dawg.add_node();
    char_to_id['f'] = tree_dawg.add_node();
    char_to_id['g'] = tree_dawg.add_node();
    char_to_id['h'] = tree_dawg.add_node();
    char_to_id['i'] = tree_dawg.add_node();
    char_to_id['l'] = tree_dawg.add_node();
    char_to_id['m'] = tree_dawg.add_node();

    // Use configure_state for the tree as well
    tree_dawg.configure_state(char_to_id['a'], {{'0', char_to_id['b']}, {'1', char_to_id['c']}});
    tree_dawg.configure_state(char_to_id['b'], {{'0', char_to_id['d']}, {'1', char_to_id['e']}});
    tree_dawg.configure_state(char_to_id['c'], {{'0', char_to_id['f']}, {'1', char_to_id['g']}});
    tree_dawg.configure_state(char_to_id['d'], {{'0', char_to_id['h']}, {'1', char_to_id['i']}});
    tree_dawg.configure_state(char_to_id['f'], {{'0', char_to_id['l']}, {'1', char_to_id['m']}});
    // Final states (e, g, h, i, l, m) already set during add_node

    tree_dawg.set_initial_state(char_to_id['a']);

    // Print as balanced parentheses
    std::string balanced_parentheses = tree_to_balanced_parentheses(tree_dawg);
    std::cout << "Tree as balanced parentheses: " << balanced_parentheses << std::endl;

    // Run minimization
    tree_dawg.minimize(); // Assuming it's a tree

    // Print DAWG information using to_string method
    std::cout << tree_dawg.to_string();

    std::cout << "--------------------" << std::endl;

    // test with generated trees
    std::string str = "(A(B(D(C(C(b))(a)(B(a)))(D(E(c)))(D(B(a))(a)(B(c)))(a))(E(B(D(a))(a)(E(C(C(C(C(b))(a)(B(a)))(D(E(c)))(D(B(a))(a)(B(c)))(b))(a)(B(a)))(D(E(c)))(D(B(a))(a)(B(c)))(b)))(C(D(c))(b)(D(c)))(B(D(b)))(B(D(B(D(a))(a)(E(b)))(C(D(c))(b)(D(c)))(B(D(b)))(a))(a)(E(b)))(C(D(c))(b)(D(c)))(B(D(b)))(b)))(C(D(c))(b)(D(c)))(B(D(b))))";
    // str = "(a(b(c(b(c(x)(y)(z))))))";

    std::string filename = "tree_generator/generated_trees/tree_bf4_t5_rp40_td4-6_letters8_mn1000000_s42.txt";

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

    LabeledTree<char> tree(str, [](const std::string &s)
                           { return s[0]; });

    std::cout << "Building DAWG from tree..." << std::endl;
    TreeDAWG<char> tree_dawg2(tree);

    // Print as balanced parentheses
    // balanced_parentheses = tree_to_balanced_parentheses(tree_dawg2);
    // std::cout << "Tree as balanced parentheses: " << balanced_parentheses << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    tree_dawg2.minimize();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Minimization took: " << duration.count() << " ms" << std::endl;
    // std::cout << tree_dawg2.to_string();

    return 0;
}
