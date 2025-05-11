#ifndef DAWG_HPP
#define DAWG_HPP

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <list>       
#include <functional> 

// Valori speciali usati nell'ordinamento per rappresentare
// l'assenza di una transizione o una transizione verso uno "stato morto".
// Questi valori devono essere esterni all'intervallo delle normali classi di equivalenza (che sono >= 0).
const int NIL_TRANSITION_SORT_VAL = -1;
const int DEAD_TARGET_SORT_VAL = -2;

// Struttura per rappresentare uno stato dell'automa
struct State
{
    int id;
    bool isFinal;
    std::map<char, int> transitions; // Transizioni: carattere -> id del nodo destinazione

    int height;
    int equivalenceClass;

    State(int _id, bool _isFinal = false)
        : id(_id), isFinal(_isFinal), height(-1), equivalenceClass(-1) {}
};

class DAWG
{
public:
    std::vector<State> nodes;
    int initialStateId;
    bool isAcyclic;

private:
    int calculateHeightDFS(int nodeId, std::vector<int> &visited_path_flags)
    {
        State &node = nodes[nodeId];
        if (node.height != -1)
            return node.height;
        visited_path_flags[nodeId] = 1;

        int current_max_h = -2;
        if (node.isFinal)
            current_max_h = 0;

        for (const auto &pair : node.transitions)
        {
            int targetNodeId = pair.second;
            if (targetNodeId >= nodes.size())
            {
                std::cerr << "Errore: transizione dal nodo " << nodeId << " al nodo inesistente " << targetNodeId << std::endl;
                continue;
            }
            if (visited_path_flags[targetNodeId] == 1)
            {
                isAcyclic = false;
                return -3;
            }
            int targetHeight = calculateHeightDFS(targetNodeId, visited_path_flags);
            if (!isAcyclic)
                return -3;
            if (targetHeight != -2 && targetHeight != -3)
            {
                current_max_h = std::max(current_max_h, 1 + targetHeight);
            }
        }
        visited_path_flags[nodeId] = 2;
        node.height = current_max_h;
        return node.height;
    }

    // Funzione ausiliaria per l'ordinamento (counting sort) di un blocco di ID di nodi.
    // Si basa su un valore estratto da ciascun nodo.
    // Restituisce una lista di nuovi blocchi (liste di ID di nodi).
    std::list<std::list<int>> countingSortBlock(
        const std::list<int> &block,
        std::function<int(int)> getValue, // Funzione per ottenere il valore di ordinamento per un nodeId
        int min_val,                      // Minimo valore possibile per getValue
        int max_val)                      // Massimo valore possibile per getValue
    {
        std::list<std::list<int>> new_blocks;
        if (block.empty())
            return new_blocks;
        if (min_val > max_val)
        { // Può succedere se il blocco è vuoto o ha solo valori speciali non previsti
            if (!block.empty())
                new_blocks.push_back(block); // Non dividere se il range non è valido
            return new_blocks;
        }

        int range = max_val - min_val + 1;
        std::vector<std::vector<int>> buckets(range);

        for (int nodeId : block)
        {
            int sort_value = getValue(nodeId);
            // Assicura che il valore di sort sia all'interno del range atteso per i bucket
            if (sort_value < min_val || sort_value > max_val)
            {
                // Questo indicherebbe un problema con la logica di min_val/max_val
                // o con i valori restituiti da getValue. Per robustezza, si potrebbe
                // avere un bucket "overflow" o gestire l'errore.
                // Per ora, assumiamo che getValue restituisca valori nel range [min_val, max_val].
                // O, più semplicemente, se capita, stampiamo un errore e non dividiamo.
                std::cerr << "Errore: sort_value " << sort_value << " fuori range [" << min_val << ", " << max_val << "]" << std::endl;
            }
            buckets[sort_value - min_val].push_back(nodeId);
        }

        for (const auto &bucket_content : buckets)
        {
            if (!bucket_content.empty())
            {
                new_blocks.emplace_back(bucket_content.begin(), bucket_content.end());
            }
        }
        return new_blocks;
    }

public:
    DAWG() : initialStateId(0), isAcyclic(true) {}

    int addNode(bool isFinal = false)
    {
        int id = nodes.size();
        nodes.emplace_back(id, isFinal);
        return id;
    }

    void addTransition(int fromNodeId, char symbol, int toNodeId)
    {
        if (fromNodeId < nodes.size() && toNodeId < nodes.size())
        {
            nodes[fromNodeId].transitions[symbol] = toNodeId;
        }
        else
        {
            std::cerr << "Errore: transizione tra nodi inesistenti." << std::endl;
        }
    }

    void setInitialState(int nodeId)
    {
        if (nodeId < nodes.size())
        {
            initialStateId = nodeId;
        }
        else
        {
            std::cerr << "Errore: tentativo di impostare uno stato iniziale inesistente." << std::endl;
        }
    }

    int computeAllHeights()
    {
        isAcyclic = true;
        for (auto &node : nodes)
            node.height = -1;
        int max_h = -1;
        std::vector<int> visited_path_flags(nodes.size(), 0);
        for (int i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].height == -1)
            {
                int h_node = calculateHeightDFS(i, visited_path_flags);
                if (!isAcyclic)
                    return -3;
                max_h = std::max(max_h, h_node);
            }
            else
            {
                max_h = std::max(max_h, nodes[i].height);
            }
        }
        return max_h;
    }

    std::map<int, int> minimize()
    {
        int max_total_height = computeAllHeights();
        if (!isAcyclic)
            return {};

        std::vector<std::vector<int>> statesByHeight(max_total_height + 1);
        for (const auto &node : nodes)
        {
            if (node.height >= 0)
            {
                statesByHeight[node.height].push_back(node.id);
            }
        }

        int nextGlobalEqClassId = 0;

        if (max_total_height >= 0 && !statesByHeight[0].empty())
        {
            for (int nodeId : statesByHeight[0])
            {
                nodes[nodeId].equivalenceClass = nextGlobalEqClassId;
            }
            nextGlobalEqClassId++;
        }

        // Cache per i vettori di transizioni ordinate (per accesso k-esimo efficiente)
        // per i nodi nel livello di altezza corrente.
        std::vector<std::vector<std::pair<char, int>>> node_sorted_transitions_cache(nodes.size());

        for (int h = 1; h <= max_total_height; ++h)
        {
            if (statesByHeight[h].empty())
                continue;

            for (int nodeId : statesByHeight[h])
            {
                node_sorted_transitions_cache[nodeId].assign(
                    nodes[nodeId].transitions.begin(),
                    nodes[nodeId].transitions.end());
            }

            std::list<std::list<int>> H_partitions;
            if (!statesByHeight[h].empty())
            {
                H_partitions.push_back({});
                for (int id : statesByHeight[h])
                    H_partitions.back().push_back(id);
            }

            // 1. Suddivisione per isFinal
            std::list<std::list<int>> P_after_isFinal;
            for (auto &block : H_partitions)
            {
                if (block.size() <= 1)
                {
                    P_after_isFinal.push_back(block);
                    continue;
                }
                auto getIsFinalVal = [&](int nodeId)
                { return nodes[nodeId].isFinal ? 1 : 0; };
                auto new_sub_blocks = countingSortBlock(block, getIsFinalVal, 0, 1); // min=0, max=1 for bool
                for (auto &sub_block : new_sub_blocks)
                    P_after_isFinal.push_back(sub_block);
            }
            H_partitions = P_after_isFinal;

            int max_k_transitions_for_H = 0;
            for (int nodeId : statesByHeight[h])
            {
                max_k_transitions_for_H = std::max(max_k_transitions_for_H, (int)node_sorted_transitions_cache[nodeId].size());
            }

            for (int k = 0; k < max_k_transitions_for_H; ++k)
            { // Per la k-esima transizione
                // 2a. Suddivisione per il carattere della k-esima transizione
                std::list<std::list<int>> P_after_char;
                int min_val_char = NIL_TRANSITION_SORT_VAL; // -1
                int max_val_char = 255;                     // Max valore per unsigned char

                for (auto &block : H_partitions)
                {
                    if (block.size() <= 1)
                    {
                        P_after_char.push_back(block);
                        continue;
                    }
                    auto getKthCharVal = [&](int nodeId)
                    {
                        if (k < node_sorted_transitions_cache[nodeId].size())
                        {
                            return static_cast<int>(static_cast<unsigned char>(node_sorted_transitions_cache[nodeId][k].first));
                        }
                        return NIL_TRANSITION_SORT_VAL;
                    };
                    // Determina il range effettivo dei valori nel blocco corrente per ottimizzare countingSort
                    int current_block_min_char = max_val_char + 1;
                    int current_block_max_char = min_val_char - 1;
                    bool has_values_char = false;
                    for (int nodeId : block)
                    {
                        int val = getKthCharVal(nodeId);
                        current_block_min_char = std::min(current_block_min_char, val);
                        current_block_max_char = std::max(current_block_max_char, val);
                        has_values_char = true;
                    }
                    if (!has_values_char)
                    {
                        P_after_char.push_back(block);
                        continue;
                    }

                    auto new_sub_blocks = countingSortBlock(block, getKthCharVal, current_block_min_char, current_block_max_char);
                    for (auto &sub_block : new_sub_blocks)
                        P_after_char.push_back(sub_block);
                }
                H_partitions = P_after_char;

                // 2b. Suddivisione per la classe di equivalenza del target della k-esima transizione
                std::list<std::list<int>> P_after_target_class;
                // Range teorico per le classi target: DEAD_TARGET_SORT_VAL (-2), NIL_TRANSITION_SORT_VAL (-1), 0 ... nextGlobalEqClassId-1
                int min_val_target_overall = DEAD_TARGET_SORT_VAL;
                int max_val_target_overall = (nextGlobalEqClassId == 0) ? DEAD_TARGET_SORT_VAL : nextGlobalEqClassId - 1;
                if (nextGlobalEqClassId == 0)
                { // Se Pi_0 era vuoto e non sono state create classi
                    max_val_target_overall = std::max(NIL_TRANSITION_SORT_VAL, DEAD_TARGET_SORT_VAL);
                }

                for (auto &block : H_partitions)
                {
                    if (block.size() <= 1)
                    {
                        P_after_target_class.push_back(block);
                        continue;
                    }
                    auto getKthTargetEqClassVal = [&](int nodeId)
                    {
                        if (k < node_sorted_transitions_cache[nodeId].size())
                        {
                            int targetNodeId = node_sorted_transitions_cache[nodeId][k].second;
                            if (nodes[targetNodeId].height < 0)
                            { // Target è uno "stato morto"
                                return DEAD_TARGET_SORT_VAL;
                            }
                            return nodes[targetNodeId].equivalenceClass; // Classe già finalizzata
                        }
                        return NIL_TRANSITION_SORT_VAL; // Nessuna k-esima transizione
                    };

                    int current_block_min_target_class = max_val_target_overall + 1; // Inizializza per trovare il min effettivo
                    int current_block_max_target_class = min_val_target_overall - 1; // Inizializza per trovare il max effettivo
                    bool has_values_target = false;
                    for (int nodeId : block)
                    {
                        int val = getKthTargetEqClassVal(nodeId);
                        current_block_min_target_class = std::min(current_block_min_target_class, val);
                        current_block_max_target_class = std::max(current_block_max_target_class, val);
                        has_values_target = true;
                    }
                    if (!has_values_target)
                    {
                        P_after_target_class.push_back(block);
                        continue;
                    }

                    auto new_sub_blocks = countingSortBlock(block, getKthTargetEqClassVal, current_block_min_target_class, current_block_max_target_class);
                    for (auto &sub_block : new_sub_blocks)
                        P_after_target_class.push_back(sub_block);
                }
                H_partitions = P_after_target_class;
            }

            for (auto &final_block : H_partitions)
            {
                if (!final_block.empty())
                {
                    for (int nodeId : final_block)
                    {
                        nodes[nodeId].equivalenceClass = nextGlobalEqClassId;
                    }
                    nextGlobalEqClassId++;
                }
            }
        }

        std::map<int, int> result;
        for (const auto &node : nodes)
        {
            if (node.height >= 0)
            {
                result[node.id] = node.equivalenceClass;
            }
            else
            {
                result[node.id] = -1; // Classe speciale per nodi non utili
            }
        }
        return result;
    }
};

#endif // DAWG_hpp
