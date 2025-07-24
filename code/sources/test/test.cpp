#include <iostream>
#include <TreeDAWG.hpp> 
#include <unordered_map>

int main()
{
    // test albero
    TreeDAWG tree_dawg;

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

    tree_dawg.set_initial_state(char_to_id['a']);

    // Use configure_state for the tree as well
    tree_dawg.configure_state(char_to_id['a'], {{'0', char_to_id['b']}, {'1', char_to_id['c']}});
    tree_dawg.configure_state(char_to_id['b'], {{'0', char_to_id['d']}, {'1', char_to_id['e']}});
    tree_dawg.configure_state(char_to_id['c'], {{'0', char_to_id['f']}, {'1', char_to_id['g']}});
    tree_dawg.configure_state(char_to_id['d'], {{'0', char_to_id['h']}, {'1', char_to_id['i']}});
    tree_dawg.configure_state(char_to_id['f'], {{'0', char_to_id['l']}, {'1', char_to_id['m']}});
    // Final states (e, g, h, i, l, m) already set during add_node

    // Run minimization
    auto equivalenceClasses_tree = tree_dawg.minimize(); // Assuming it's a tree

    std::cout << "Node equivalence classes:" << std::endl;
    for (const auto &pair : equivalenceClasses_tree)
    {
        // Find the character that corresponds to this node ID
        char node_char = '?';
        for (const auto &mapping : char_to_id)
        {
            if (mapping.second == pair.first)
            {
                node_char = mapping.first;
                break;
            }
        }

        std::cout << "Node " << pair.first
                  << " (char '" << node_char << "', height " << tree_dawg[pair.first].get_height()
                  << "'): Class " << pair.second << std::endl;
    }

    return 0;
}
