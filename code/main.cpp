#include <iostream>
#include "DAWG.hpp"

int main()
{
    DAWG dawg;

    // Costruisci l'automa dell'esempio Fig. 1 dell'articolo (adattato a ID da 0 a 14)
    // Gli stati sono numerati da 0 a 14 invece che da 1 a 15 per comodità con i vettori
    for (int i = 0; i < 15; ++i)
    {
        dawg.add_node();
    }

    // Imposta stati finali (originali: 5, 14, 15 -> nostri: 4, 13, 14)
    dawg.nodes[4].is_final = true;  // Stato 5
    dawg.nodes[13].is_final = true; // Stato 14
    dawg.nodes[14].is_final = true; // Stato 15

    dawg.set_initial_state(0); // Stato 1 originale

    // Aggiungi transizioni (adattando gli ID)
    dawg.add_transition(0, 'a', 1); // 1->2
    dawg.add_transition(0, 'b', 2); // 1->3
    dawg.add_transition(0, 'c', 3); // 1->4

    // Livello 2 -> Livello 3
    dawg.add_transition(1, 'a', 4); // 2->5
    dawg.add_transition(1, 'b', 5); // 2->6
    dawg.add_transition(2, 'b', 6); // 3->7
    dawg.add_transition(2, 'a', 9); // 3->10

    dawg.add_transition(3, 'a', 7); // 4->8
    dawg.add_transition(3, 'b', 8); // 4->9

    dawg.add_transition(4, 'b', 9);  // 5->10
    dawg.add_transition(4, 'a', 14); // 5->15
    dawg.add_transition(5, 'a', 9);  // 6->10
    dawg.add_transition(5, 'b', 10); // 6->11

    dawg.add_transition(6, 'a', 9);  // 7->10
    dawg.add_transition(6, 'b', 10); // 7->11

    dawg.add_transition(7, 'a', 11); // 8->12
    dawg.add_transition(8, 'a', 11); // 9->12
    dawg.add_transition(8, 'b', 14); // 9->15

    dawg.add_transition(9, 'b', 12);  // 10->13
    dawg.add_transition(9, 'a', 14);  // 10->15
    dawg.add_transition(10, 'a', 12); // 11->13
    dawg.add_transition(11, 'a', 13); // 12->14
    dawg.add_transition(11, 'c', 14); // 12->15

    dawg.add_transition(12, 'b', 14); // 13->15
    dawg.add_transition(13, 'd', 14); // 14->15

    // Esegui la minimizzazione
    auto equivalenceClasses = dawg.minimize(false);

    std::cout << "Node equivalence classes:" << std::endl;
    for (const auto &pair : equivalenceClasses)
    {
        std::cout << "Nodo " << pair.first << " (originale " << pair.first + 1
                  << ", altezza " << dawg.nodes[pair.first].get_height()
                  << "): Classe " << pair.second << std::endl;
    }

    // L'articolo menziona che gli stati 6 e 7 (nostri 5 e 6) sono equivalenti.
    // Verifichiamo se hanno la stessa classe di equivalenza.
    if (dawg.nodes[5].get_height() >= 0 && dawg.nodes[6].get_height() >= 0 && // assicuriamoci che siano stati validi
        equivalenceClasses.count(5) && equivalenceClasses.count(6) &&
        equivalenceClasses[5] == equivalenceClasses[6])
    {
        std::cout << "\nConferma: Nodo 5 (originale 6) e Nodo 6 (originale 7) sono nella stessa classe di equivalenza: "
                  << equivalenceClasses[5] << std::endl;
    }
    else if (dawg.nodes[5].get_height() >= 0 && dawg.nodes[6].get_height() >= 0)
    {
        std::cout << "\nNota: Nodo 5 (originale 6) e Nodo 6 (originale 7) NON sono nella stessa classe di equivalenza." << std::endl;
        if (equivalenceClasses.count(5))
            std::cout << "Classe nodo 5: " << equivalenceClasses[5] << std::endl;
        if (equivalenceClasses.count(6))
            std::cout << "Classe nodo 6: " << equivalenceClasses[6] << std::endl;
    }

    // test albero
    DAWG tree_dawg;

    // Node mapping from character labels to integer IDs
    std::map<char, size_t> char_to_id;
    char_to_id['a'] = tree_dawg.add_node(false);
    char_to_id['b'] = tree_dawg.add_node(false);
    char_to_id['c'] = tree_dawg.add_node(false);
    char_to_id['d'] = tree_dawg.add_node(false);
    char_to_id['e'] = tree_dawg.add_node(true);
    char_to_id['f'] = tree_dawg.add_node(false);
    char_to_id['g'] = tree_dawg.add_node(true);
    char_to_id['h'] = tree_dawg.add_node(true);
    char_to_id['i'] = tree_dawg.add_node(true);
    char_to_id['l'] = tree_dawg.add_node(true);
    char_to_id['m'] = tree_dawg.add_node(true);

    tree_dawg.set_initial_state(char_to_id['a']);

    // Add transitions based on the diagram
    tree_dawg.add_transition(char_to_id['a'], '0', char_to_id['b']);
    tree_dawg.add_transition(char_to_id['a'], '1', char_to_id['c']);

    tree_dawg.add_transition(char_to_id['b'], '0', char_to_id['d']);
    tree_dawg.add_transition(char_to_id['b'], '1', char_to_id['e']);

    tree_dawg.add_transition(char_to_id['c'], '0', char_to_id['f']);
    tree_dawg.add_transition(char_to_id['c'], '1', char_to_id['g']);

    tree_dawg.add_transition(char_to_id['d'], '0', char_to_id['h']);
    tree_dawg.add_transition(char_to_id['d'], '1', char_to_id['i']);

    tree_dawg.add_transition(char_to_id['f'], '0', char_to_id['l']);
    tree_dawg.add_transition(char_to_id['f'], '1', char_to_id['m']);

    // Run minimization
    auto equivalenceClasses_tree = tree_dawg.minimize(true); // Assuming it's a tree

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
                  << " (char '" << node_char << "', height " << tree_dawg.nodes[pair.first].get_height()
                  << "'): Class " << pair.second << std::endl;
    }

    return 0;
}
