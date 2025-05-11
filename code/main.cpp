#include <iostream>
#include "DAWG.hpp"

int main()
{
    DAWG dawg;

    // Costruisci l'automa dell'esempio Fig. 1 dell'articolo (adattato a ID da 0 a 14)
    // Gli stati sono numerati da 0 a 14 invece che da 1 a 15 per comodità con i vettori
    for (int i = 0; i < 15; ++i)
    {
        dawg.addNode();
    }

    // Imposta stati finali (originali: 5, 14, 15 -> nostri: 4, 13, 14)
    dawg.nodes[4].isFinal = true;  // Stato 5
    dawg.nodes[13].isFinal = true; // Stato 14
    dawg.nodes[14].isFinal = true; // Stato 15

    dawg.setInitialState(0); // Stato 1 originale

    // Aggiungi transizioni (adattando gli ID)
    dawg.addTransition(0, 'a', 1); // 1->2
    dawg.addTransition(0, 'b', 2); // 1->3
    dawg.addTransition(0, 'c', 3); // 1->4

    // Livello 2 -> Livello 3
    dawg.addTransition(1, 'a', 4); // 2->5
    dawg.addTransition(1, 'b', 5); // 2->6
    dawg.addTransition(2, 'b', 6); // 3->7
    dawg.addTransition(2, 'a', 9); // 3->10
    dawg.addTransition(2, 'b', 6); // 3->7

    dawg.addTransition(3, 'a', 7); // 4->8
    dawg.addTransition(3, 'b', 8); // 4->9

    dawg.addTransition(4, 'b', 9);  // 5->10
    dawg.addTransition(4, 'a', 14); // 5->15
    dawg.addTransition(5, 'a', 9);  // 6->10
    dawg.addTransition(5, 'b', 10); // 6->11

    dawg.addTransition(6, 'a', 9);  // 7->10
    dawg.addTransition(6, 'b', 10); // 7->11

    dawg.addTransition(7, 'a', 11); // 8->12
    dawg.addTransition(8, 'a', 11); // 9->12
    dawg.addTransition(8, 'b', 14); // 9->15

    dawg.addTransition(9, 'b', 12);  // 10->13
    dawg.addTransition(9, 'a', 14);  // 10->15
    dawg.addTransition(10, 'a', 12); // 11->13
    dawg.addTransition(11, 'a', 13); // 12->14
    dawg.addTransition(11, 'c', 14); // 12->15

    dawg.addTransition(12, 'b', 14); // 13->15
    dawg.addTransition(13, 'd', 14); // 14->15

    // Esegui la minimizzazione
    std::map<int, int> equivalenceClasses = dawg.minimize();

    if (dawg.isAcyclic)
    {
        std::cout << "Classi di equivalenza dei nodi:" << std::endl;
        for (const auto &pair : equivalenceClasses)
        {
            std::cout << "Nodo " << pair.first << " (originale " << pair.first + 1
                      << ", altezza " << dawg.nodes[pair.first].height
                      << "): Classe " << pair.second << std::endl;
        }

        // L'articolo menziona che gli stati 6 e 7 (nostri 5 e 6) sono equivalenti.
        // Verifichiamo se hanno la stessa classe di equivalenza.
        if (dawg.nodes[5].height >= 0 && dawg.nodes[6].height >= 0 && // assicuriamoci che siano stati validi
            equivalenceClasses.count(5) && equivalenceClasses.count(6) &&
            equivalenceClasses[5] == equivalenceClasses[6])
        {
            std::cout << "\nConferma: Nodo 5 (originale 6) e Nodo 6 (originale 7) sono nella stessa classe di equivalenza: "
                      << equivalenceClasses[5] << std::endl;
        }
        else if (dawg.nodes[5].height >= 0 && dawg.nodes[6].height >= 0)
        {
            std::cout << "\nNota: Nodo 5 (originale 6) e Nodo 6 (originale 7) NON sono nella stessa classe di equivalenza." << std::endl;
            if (equivalenceClasses.count(5))
                std::cout << "Classe nodo 5: " << equivalenceClasses[5] << std::endl;
            if (equivalenceClasses.count(6))
                std::cout << "Classe nodo 6: " << equivalenceClasses[6] << std::endl;
        }
    }
    else
    {
        std::cout << "L'algoritmo non può procedere a causa di cicli rilevati." << std::endl;
    }

    return 0;
}
