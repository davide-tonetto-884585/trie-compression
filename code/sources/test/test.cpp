#include <iostream>
#include <TreeDAWG.hpp> 
#include <unordered_map>

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

    tree_dawg.set_initial_state(char_to_id['a']);

    // Use configure_state for the tree as well
    tree_dawg.configure_state(char_to_id['a'], {{'0', char_to_id['b']}, {'1', char_to_id['c']}});
    tree_dawg.configure_state(char_to_id['b'], {{'0', char_to_id['d']}, {'1', char_to_id['e']}});
    tree_dawg.configure_state(char_to_id['c'], {{'0', char_to_id['f']}, {'1', char_to_id['g']}});
    tree_dawg.configure_state(char_to_id['d'], {{'0', char_to_id['h']}, {'1', char_to_id['i']}});
    tree_dawg.configure_state(char_to_id['f'], {{'0', char_to_id['l']}, {'1', char_to_id['m']}});
    // Final states (e, g, h, i, l, m) already set during add_node

    // Run minimization
    tree_dawg.minimize(); // Assuming it's a tree

    std::cout << "Iterating over all nodes using range-based for loop:" << std::endl;
    for (const auto &node : tree_dawg)
    {
        // Find the character that corresponds to this node ID
        char node_char = '?';
        for (const auto &mapping : char_to_id)
        {
            if (mapping.second == node.get_id())
            {
                node_char = mapping.first;
                break;
            }
        }

        std::cout << "Node " << node.get_id()
                  << " (char '" << node_char 
                  << "', id: " << node.get_id()
                  << ", class: " << node.get_equivalence_class()
                  << ", height: " << node.get_height()
                  << ", is_root: " << (node.get_id() == tree_dawg.get_initial_state_id() ? "yes" : "no")
                  << ", is_final: " << (node.is_final() ? "yes" : "no")
                  << ", transitions: " << node.get_transitions().size() << ")" << std::endl;
    }

    std::cout << "\n=== Testing with different label types ===" << std::endl;
    
    // Test with integer labels
    TreeDAWG<int> int_dawg;
    auto node1 = int_dawg.add_node();
    auto node2 = int_dawg.add_node();
    auto node3 = int_dawg.add_node();
    
    int_dawg.set_initial_state(node1);
    int_dawg.configure_state(node1, {{100, node2}, {200, node3}});
    
    std::cout << "Integer DAWG created with " << int_dawg.get_num_nodes() << " nodes" << std::endl;
    for (const auto &node : int_dawg) {
        std::cout << "Node " << node.get_id() << " has " << node.get_transitions().size() << " transitions" << std::endl;
        for (const auto &[label, target] : node.get_transitions()) {
            std::cout << "  Transition: " << label << " -> " << target << std::endl;
        }
    }
    
    // Test with string labels
    TreeDAWG<std::string> string_dawg;
    auto str_node1 = string_dawg.add_node();
    auto str_node2 = string_dawg.add_node();
    
    string_dawg.set_initial_state(str_node1);
    string_dawg.configure_state(str_node1, {{"hello", str_node2}});
    
    std::cout << "\nString DAWG created with " << string_dawg.get_num_nodes() << " nodes" << std::endl;
    for (const auto &node : string_dawg) {
        std::cout << "Node " << node.get_id() << " has " << node.get_transitions().size() << " transitions" << std::endl;
        for (const auto &[label, target] : node.get_transitions()) {
            std::cout << "  Transition: \"" << label << "\" -> " << target << std::endl;
        }
    }

    return 0;
}
