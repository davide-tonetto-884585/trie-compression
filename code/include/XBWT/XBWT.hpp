#ifndef XBWT_HPP
#define XBWT_HPP

#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <iomanip>
#include <unordered_map>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>

#include "../LabeledTree/LabeledTree.hpp"

template <typename T1, typename T2, typename T3>
struct Triplet
{
    T1 first;
    T2 second;
    T3 third;

    Triplet(T1 f, T2 s, T3 t) : first(f), second(s), third(t) {}
    Triplet(const Triplet &other) : first(other.first), second(other.second), third(other.third) {}

    std::string join(const std::string &sep) const
    {
        return first + sep + second + sep + third;
    }
};

std::string padLeft(const std::string &str, unsigned int num, char paddingChar = ' ')
{
    return std::string(num - str.size(), paddingChar) + str;
}

template <typename T>
class XBWT
{
private:
    struct MyImpl
    {
        sdsl::rrr_vector<> SLastCompressed;                // operator[] has complexity O(k), where k is a small constant.
        sdsl::wt_int<sdsl::rrr_vector<>> SAlphaCompressed; // operator[] has complexity O(log |Σ|), where |Σ| denotes the size of the alphabet
        sdsl::rrr_vector<> SAlphaBitCompressed;
        sdsl::rrr_vector<> ACompressed;
        sdsl::rrr_vector<>::rank_1_type SLastCompressedRank;
        sdsl::rrr_vector<>::select_1_type SLastCompressedSelect;
        sdsl::rrr_vector<>::rank_1_type ACompressedRank;
        sdsl::rrr_vector<>::select_1_type ACompressedSelect;

        sdsl::rrr_vector<> SigmaNCompressed;
        sdsl::rrr_vector<>::rank_1_type SigmaNCompressedRank;

        std::unordered_map<T, unsigned int> alphabetMap;
        std::unordered_map<unsigned int, T> alphabetMapInv;

        unsigned int cardSigma;
        unsigned int cardSigmaN;
        unsigned int maxNumDigits;
    };

    std::unique_ptr<MyImpl> pImpl;

    // private methods for building the XBWT
    void createXBWT(const LabeledTree<T> &tree, bool usePathSort = true, bool verbose = false, std::vector<unsigned int> *intNodesPosSorted = nullptr);
    std::vector<unsigned int> pathSort(const std::vector<Triplet<unsigned int, int, int>> &intNodes, bool dummyRoot = false, bool firstIt = true, unsigned int rem = 0, bool j0Encountered = false, unsigned int maxNumDigits = 1) const;
    std::vector<unsigned int> pathSortMerge(std::vector<unsigned int> &intNodesPosNotJSorted, std::vector<unsigned int> &firstIndexIntNodesPosNotJSorted, std::vector<bool> &indexFoundIntNodesPosNotJSorted, std::vector<unsigned int> &intNodesPosJSorted, std::vector<Triplet<unsigned int, int, int>> &tempIntNodes, unsigned int numDummyNodes, short int jv, bool dummyRoot = false, bool firstIt = false, bool j0Encountered = false) const;
    Node<unsigned int> *contractTree(std::vector<Triplet<unsigned int, int, int>> intNodes, short int j, unsigned int *maxNumDigits, std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> *tripletsSorted = nullptr) const;

    template <typename U>
    std::vector<Triplet<unsigned int, int, int>> computeIntNodes(const Node<U> &root, bool encode = false) const;

    // private methods for inverting the XBWT transform
    std::vector<unsigned int> buildF() const;
    std::vector<long int> buildJ() const;

    // utility methods
    void radixSort(std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> &arr) const;

public:
    XBWT(const LabeledTree<T> &tree, bool usePathSort = true, bool verbose = false, std::vector<unsigned int> *intNodesPosSorted = nullptr);
    ~XBWT();

    // rebuild tree
    LabeledTree<T> rebuildTree() const;

    // navigation methods
    std::pair<long int, long int> getChildren(unsigned int i) const;
    long int getRankedChild(unsigned int i, unsigned int k) const;
    long int getCharRankedChild(unsigned int i, T label, unsigned int k) const;
    unsigned int getDegree(unsigned int i) const;
    unsigned int getCharDegree(unsigned int i, T label) const;
    long int getParent(unsigned int i) const;
    std::vector<T> getSubtree(unsigned int i, unsigned int order = 0) const;

    // tree search methods
    std::pair<long int, long int> subPathSearch(const std::vector<T> &path) const;

    // testing methods
    unsigned int getCardSigma() const;
    T getNodeLabel(unsigned int i) const;
    std::vector<T> getUpwardPath(unsigned int i) const;
    std::vector<unsigned int> upwardStableSortConstruction(const LabeledTree<T> &tree, std::vector<Triplet<unsigned int, int, int>> *intNodes = nullptr) const;
};

/**
 * @brief Constructor for the XBWT class.
 *
 * This constructor initializes an XBWT object using a labeled tree,
 * the cardinality of the alphabet, and the cardinality of the internal alphabet.
 *
 * @param tree The labeled tree used to construct the XBWT.
 * @param cardSigma The cardinality of the alphabet.
 * @param cardSigmaN The cardinality of the internal alphabet.
 * @param usePathSort If true, uses path sort to construct the XBWT.
 * @param verbose If true, enables verbose mode for debugging.
 * @param intNodesPosSorted An optional vector of sorted internal node positions.
 */
template <typename T>
XBWT<T>::XBWT(const LabeledTree<T> &tree, bool usePathSort, bool verbose, std::vector<unsigned int> *intNodesPosSorted)
{
    pImpl = std::make_unique<MyImpl>();

    auto nodes = tree.getNodes();
    if (nodes.empty())
    {
        throw std::invalid_argument("The tree is empty");
    }

    unsigned int cardSigmaN = 0;
    std::vector<T> internalLabels;
    std::vector<T> leafLabels;
    
    // Separate internal node labels from leaf labels
    for (const auto &node : nodes)
    {
        if (node->isLeaf())
        {
            auto res = pImpl->alphabetMap.find(node->getLabel());
            if (res == pImpl->alphabetMap.end())
            {
                leafLabels.push_back(node->getLabel());
                pImpl->alphabetMap.emplace(node->getLabel(), 0); // placeholder
            }
        }
        else
        {
            auto res = pImpl->alphabetMap.find(node->getLabel());
            if (res == pImpl->alphabetMap.end())
            {
                internalLabels.push_back(node->getLabel());
                pImpl->alphabetMap.emplace(node->getLabel(), 0); // placeholder
                ++cardSigmaN;
            }
        }
    }

    // Sort both label sets to maintain natural order within each group
    std::sort(internalLabels.begin(), internalLabels.end());
    std::sort(leafLabels.begin(), leafLabels.end());

    // Assign codes: internal nodes get 1, 2, 3, ..., cardSigmaN
    // Leaf nodes get cardSigmaN+1, cardSigmaN+2, ..., cardSigma
    unsigned int maxNumDigits = 1;
    unsigned int code = 1;
    
    // First assign codes to internal nodes
    for (auto label : internalLabels)
    {
        pImpl->alphabetMap[label] = code;
        pImpl->alphabetMapInv.emplace(code, label);
        if (std::to_string(code).size() > maxNumDigits)
            maxNumDigits = std::to_string(code).size();
        ++code;
    }
    
    // Then assign codes to leaf nodes (ensuring they're greater than internal nodes)
    for (auto label : leafLabels)
    {
        pImpl->alphabetMap[label] = code;
        pImpl->alphabetMapInv.emplace(code, label);
        if (std::to_string(code).size() > maxNumDigits)
            maxNumDigits = std::to_string(code).size();
        ++code;
    }
    
    pImpl->cardSigma = code - 1; // Total number of labels
    assert(pImpl->cardSigma == pImpl->alphabetMap.size());
    pImpl->cardSigmaN = cardSigmaN;
    pImpl->maxNumDigits = maxNumDigits;

    if (verbose)
    {
        std::cout << "Cardinality of the alphabet: " << pImpl->cardSigma << std::endl;
        std::cout << "Cardinality of the internal alphabet: " << pImpl->cardSigmaN << std::endl;
    }

    createXBWT(tree, usePathSort, verbose, intNodesPosSorted);
}

template <typename T>
XBWT<T>::~XBWT() = default;

/**
 * @brief Performs radix sort on a vector of pairs.
 *
 * This function sorts a vector of pairs using radix sort. Each pair consists of an unsigned int
 * and a Triplet<unsigned int, int, int>. The sorting is based on the integer values derived from
 * the Triplet.
 *
 * @param arr The vector of pairs to be sorted.
 */
template <typename T>
void XBWT<T>::radixSort(std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> &arr) const
{
    auto countSort = [](std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> &arr, size_t exp)
    {
        std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> output(arr.size(), std::pair<unsigned int, Triplet<std::string, std::string, std::string>>(0, Triplet<std::string, std::string, std::string>("", "", "")));
        std::vector<unsigned int> count(10, 0);

        for (const auto &pair : arr)
        {
            size_t val = static_cast<size_t>(std::stoull(pair.second.join("")));
            ++count[(val / exp) % 10];
        }

        for (unsigned int i = 1; i < 10; ++i)
        {
            count[i] += count[i - 1];
        }

        for (long long int i = arr.size() - 1; i >= 0; --i)
        {
            size_t val = static_cast<size_t>(std::stoull(arr[i].second.join("")));
            output[count[(val / exp) % 10] - 1] = arr[i];
            --count[(val / exp) % 10];
        }

        for (long long int i = 0; i < arr.size(); ++i)
        {
            arr[i] = output[i];
        }
    };

    // TODO: check if this is valid for big integers
    size_t max = static_cast<size_t>(std::stoull(arr[0].second.join("")));
    for (const auto &el : arr)
    {
        size_t val = static_cast<size_t>(std::stoull(el.second.join("")));
        if (val > max)
        {
            max = val;
        }
    }

    for (size_t exp = 1; max / exp > 0; exp *= 10)
    {
        countSort(arr, exp);
    }
}

/**
 * @brief Return an order for the intNodes array of the labeled tree.
 *
 * This function constructs the XBWT from a labeled tree. It uses a stable sort to order the internal nodes
 * upward paths. Not linear time complexity as for pathSort function.
 *
 * @param tree The labeled tree used to construct the XBWT.
 * @param intNodes A vector of Triplets representing the internal nodes of the tree.
 * @return intNodes sorted positions needed for XBWT construction
 */
template <typename T>
std::vector<unsigned int> XBWT<T>::upwardStableSortConstruction(const LabeledTree<T> &tree, std::vector<Triplet<unsigned int, int, int>> *intNodes) const
{
    std::vector<Triplet<unsigned int, int, int>> localIntNodes;
    if (intNodes == nullptr)
    {
        localIntNodes = computeIntNodes(*tree.getRoot(), true);
        intNodes = &localIntNodes;
    }
    else
    {
        *intNodes = computeIntNodes(*tree.getRoot(), true);
    }

    std::vector<std::pair<unsigned int, std::string>> intNodesPaths;
    intNodesPaths.push_back(std::pair<unsigned int, std::string>(0, ""));
    for (unsigned int i = 1; i < (*intNodes).size(); ++i)
    {
        std::ostringstream oss;
        auto parentNode = (*intNodes)[(*intNodes)[i].third - 1];
        oss << std::setfill('0') << std::setw(pImpl->maxNumDigits) << parentNode.first;
        do
        {
            if (parentNode.third != 0)
            {
                parentNode = (*intNodes)[parentNode.third - 1];
                oss << std::setfill('0') << std::setw(pImpl->maxNumDigits) << parentNode.first;
            }
        } while (parentNode.third != 0);

        intNodesPaths.push_back(std::pair<unsigned int, std::string>(i, oss.str()));
    }

    std::stable_sort(intNodesPaths.begin(), intNodesPaths.end(), [](const std::pair<unsigned int, std::string> &a, const std::pair<unsigned int, std::string> &b)
                     { return a.second < b.second; });

    std::vector<unsigned int> intNodesPosSorted(intNodesPaths.size(), 0);
    for (unsigned int i = 0; i < intNodesPaths.size(); ++i)
    {
        intNodesPosSorted[i] = intNodesPaths[i].first;
        // std::cout << intNodesPaths[i].first << " " << intNodesPaths[i].second << std::endl;
    }

    return intNodesPosSorted;
}

/**
 * @brief Computes the internal nodes of a labeled tree.
 *
 * This function traverses a labeled tree in preorder and computes the internal nodes.
 * Each internal node is represented as a Triplet containing the node's label, its level,
 * and the index of its parent node in the array (from 1 to t, root has parent 0).
 *
 * @param root The root node of the labeled tree.
 * @return A vector of Triplets representing the internal nodes of the tree.
 */
template <typename T>
template <typename U>
std::vector<Triplet<unsigned int, int, int>> XBWT<T>::computeIntNodes(const Node<U> &root, bool encode) const
{
    std::vector<Triplet<unsigned int, int, int>> intNodes;

    // stack to visit nodes in reorder traversal
    std::stack<Triplet<const Node<U> *, int, int>> stack;
    stack.push(Triplet<const Node<U> *, int, int>(&root, 0, 0));

    while (!stack.empty())
    {
        auto [node, cur_level, parent_index] = stack.top();
        stack.pop();

        T label = (encode ? pImpl->alphabetMap[node->getLabel()] : node->getLabel());
        intNodes.push_back(Triplet<unsigned int, int, int>(label, cur_level, parent_index));

        // update level and parent index
        parent_index = intNodes.size();
        ++cur_level;

        // visit children in reverse order
        for (auto it = node->getChildren().rbegin(); it != node->getChildren().rend(); ++it)
        {
            stack.push(Triplet<const Node<U> *, int, int>(*it, cur_level, parent_index));
        }
    }

    return intNodes;
}

/**
 * @brief Contracts a labeled tree by removing nodes on levels equal to j (mod 3).
 *
 * This function contracts a labeled tree by removing
 * nodes on levels equal to j (mod 3). If sorted triplets are provided (in the first iteration
 * of pathSort algorithm) then it uses them to rename the nodes.
 *
 * @param intNodes A vector of Triplets representing the internal nodes of the tree.
 * @param j The level used to determine which nodes to remove.
 * @param tripletsSorted An optional vector used to rename the nodes.
 * @return A pointer to the root node of the contracted tree.
 */
template <typename T>
Node<unsigned int> *XBWT<T>::contractTree(std::vector<Triplet<unsigned int, int, int>> intNodes, short int j, unsigned int *maxNumDigits, std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> *tripletsSorted) const
{
    if (tripletsSorted)
    {
        unsigned int newName = 2;
        *maxNumDigits = 1;
        auto &preTriplet = tripletsSorted->front().second;
        for (auto &pair : *tripletsSorted)
        {
            if (pair.second.first != preTriplet.first || pair.second.second != preTriplet.second || pair.second.third != preTriplet.third)
            {
                preTriplet = pair.second;
                ++newName;
            }

            intNodes[pair.first].first = newName;
            if (std::to_string(newName).size() > *maxNumDigits)
                *maxNumDigits = std::to_string(newName).size();
        }
    }

    if (j == 0)
        intNodes[0].first = 1;

    short int jNext[] = {1, 2, 0};
    std::vector<std::pair<unsigned int, unsigned int>> edges;
    if (j == 0)
    {
        for (int i = intNodes.size() - 1; i >= 0; --i) // TODO: complexity can be reduced
        {
            auto &it = intNodes[i];
            if (it.second % 3 != j)
            {
                if (it.second % 3 == jNext[j])
                {
                    // if the parent is the root
                    if (it.third == 1)
                        edges.push_back(std::pair<unsigned int, unsigned int>(it.third - 1, i));
                    else
                        edges.push_back(std::pair<unsigned int, unsigned int>(intNodes[it.third - 1].third - 1, i));
                }
                else
                    edges.push_back(std::pair<unsigned int, unsigned int>(it.third - 1, i));
            }
        }
    }
    else
    {
        for (int i = intNodes.size() - 1; i > 0; --i)
        {
            auto &it = intNodes[i];
            if (it.second % 3 != j)
            {
                if (it.second % 3 == jNext[j])
                    edges.push_back(std::pair<unsigned int, unsigned int>(intNodes[it.third - 1].third - 1, i));
                else
                    edges.push_back(std::pair<unsigned int, unsigned int>(it.third - 1, i));
            }
        }
    }

    // create tree by adding edges in reverse order
    std::vector<Node<unsigned int> *> nodes(intNodes.size(), nullptr);
    for (int i = edges.size() - 1; i >= 0; --i)
    {
        if (nodes[edges[i].first] == nullptr)
            nodes[edges[i].first] = new Node<unsigned int>(intNodes[edges[i].first].first);

        if (nodes[edges[i].second] == nullptr)
            nodes[edges[i].second] = new Node<unsigned int>(intNodes[edges[i].second].first);

        nodes[edges[i].first]->pushBackChild(nodes[edges[i].second]);
    }

    return nodes[0];
}

/**
 * @brief Merges two sorted vectors of node positions based on specific conditions.
 *
 * This function merges two sorted vectors of node positions (`intNodesPosNotJSorted` and `intNodesPosJSorted`)
 * based on specific conditions related to the levels of the nodes and their values. This in order
 * to get a final ordering of the nodes at each recursive step of pathSort algorithm.
 *
 * @param intNodesPosNotJSorted A vector of sorted node positions not at level j.
 * @param firstIndexIntNodesPosNotJSorted A vector containing the first index of elements in `intNodesPosNotJSorted`.
 * @param indexFoundIntNodesPosNotJSorted A vector indicating if a given element is found in `intNodesPosNotJSorted`.
 * @param intNodesPosJSorted A vector of sorted node positions at level j.
 * @param tempIntNodes A vector of internal nodes as described in function `computeIntNodes`.
 * @param numDummyNodes The number of dummy nodes prepended in `tempIntNodes`.
 * @param jv The level j used to determine the merging conditions.
 * @param dummyRoot A boolean indicating if a dummy root is used.
 * @param firstIt A boolean indicating if this is the first iteration of the `pathSort` algorithm.
 * @return A merged vector of node positions sorted.
 */
template <typename T>
std::vector<unsigned int> XBWT<T>::pathSortMerge(std::vector<unsigned int> &intNodesPosNotJSorted, std::vector<unsigned int> &firstIndexIntNodesPosNotJSorted, std::vector<bool> &indexFoundIntNodesPosNotJSorted, std::vector<unsigned int> &intNodesPosJSorted, std::vector<Triplet<unsigned int, int, int>> &tempIntNodes, unsigned int numDummyNodes, short int jv, bool dummyRoot, bool firstIt, bool j0Encountered) const
{
    std::vector<unsigned int> merged(intNodesPosNotJSorted.size() + intNodesPosJSorted.size(), 0);
    std::vector<long int> firstIndexMerged(merged.size(), -1);
    short int jCond[] = {1, 2, 0};
    bool flag = true;
    unsigned int i = 0, j = 0, cont = 0; // i for not j, j for j
    while (flag)
    {
        if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second % 3 == jCond[jv])
        {
            if (!firstIt)
            {
                unsigned int v1 = tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].first;
                unsigned int v2 = tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].first;

                if (v1 == v2)
                {
                    // if same lv
                    // if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second == tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].second)
                    // {
                    if (!indexFoundIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].third - numDummyNodes] ||
                        !indexFoundIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].third - numDummyNodes])
                        throw std::runtime_error("Error: index not found in merge");

                    unsigned int index1 = firstIndexIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].third - numDummyNodes];
                    unsigned int index2 = firstIndexIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].third - numDummyNodes];

                    if (index1 > index2) // prima era >
                        merged[cont++] = intNodesPosNotJSorted[i++];
                    else
                        merged[cont++] = intNodesPosJSorted[j++];
                    /*}
                    else
                    {
                        if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second < tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].second)
                            merged[cont++] = intNodesPosNotJSorted[i++];
                        else
                            merged[cont++] = intNodesPosJSorted[j++];
                    }*/
                }
                else if (v1 < v2)
                    merged[cont++] = intNodesPosNotJSorted[i++];
                else
                    merged[cont++] = intNodesPosJSorted[j++];
            }
            else
            {
                std::pair<unsigned int, unsigned int> pair1{
                    tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].first,
                    tempIntNodes[tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].third].first};

                std::pair<unsigned int, unsigned int> pair2{
                    tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].first,
                    tempIntNodes[tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].third].first};

                if (pair1.first == pair2.first)
                {
                    if (pair1.second == pair2.second)
                    {
                        // must use merged to compare positions of parents
                        if (firstIndexMerged[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes] != -1 ||
                            firstIndexMerged[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes] != -1)
                        {
                            unsigned int index1 = firstIndexMerged[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes];
                            unsigned int index2 = firstIndexMerged[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes];
                            if (index1 > index2)
                                merged[cont++] = intNodesPosNotJSorted[i++];
                            else
                                merged[cont++] = intNodesPosJSorted[j++];
                        }
                        else
                        {
                            unsigned int v1 = tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].first;
                            unsigned int v2 = tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].first;

                            if (v1 == v2)
                            {
                                if (!indexFoundIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].third - numDummyNodes] ||
                                    !indexFoundIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].third - numDummyNodes])
                                    throw std::runtime_error("Error: index not found in merge");

                                unsigned int index1 = firstIndexIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].third - numDummyNodes];
                                unsigned int index2 = firstIndexIntNodesPosNotJSorted[tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].third - numDummyNodes];

                                if (index1 > index2) // prima era >
                                    merged[cont++] = intNodesPosNotJSorted[i++];
                                else
                                    merged[cont++] = intNodesPosJSorted[j++];
                            }
                            else if (v1 < v2)
                                merged[cont++] = intNodesPosNotJSorted[i++];
                            else
                                merged[cont++] = intNodesPosJSorted[j++];
                        }
                    }
                    else if (pair1.second < pair2.second)
                        merged[cont++] = intNodesPosNotJSorted[i++];
                    else
                        merged[cont++] = intNodesPosJSorted[j++];
                }
                else if (pair1.first < pair2.first)
                    merged[cont++] = intNodesPosNotJSorted[i++];
                else
                    merged[cont++] = intNodesPosJSorted[j++];
            }
        }
        else if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second % 3 == jCond[jCond[jv]])
        {
            if (!firstIt)
            {
                unsigned int v1 = tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].first;
                unsigned int v2 = tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].first;

                if (v1 == v2)
                {
                    // if same lv
                    //  if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second == tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].second)
                    //  {
                    if (!indexFoundIntNodesPosNotJSorted[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes] ||
                        !indexFoundIntNodesPosNotJSorted[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes])
                        throw std::runtime_error("Error: index x not found in merge");

                    unsigned int index1 = firstIndexIntNodesPosNotJSorted[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes];
                    unsigned int index2 = firstIndexIntNodesPosNotJSorted[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes];

                    if (index1 > index2) // prima era >
                        merged[cont++] = intNodesPosNotJSorted[i++];
                    else
                        merged[cont++] = intNodesPosJSorted[j++];
                    /*}
                    else
                    {
                        if (tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].second < tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].second)
                            merged[cont++] = intNodesPosNotJSorted[i++];
                        else
                            merged[cont++] = intNodesPosJSorted[j++];
                    } */
                }
                else if (v1 < v2)
                    merged[cont++] = intNodesPosNotJSorted[i++];
                else
                    merged[cont++] = intNodesPosJSorted[j++];
            }
            else
            {
                unsigned int v1 = tempIntNodes[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third].first;
                unsigned int v2 = tempIntNodes[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third].first;

                if (v1 == v2)
                {
                    if (!indexFoundIntNodesPosNotJSorted[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes] ||
                        !indexFoundIntNodesPosNotJSorted[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes])
                        throw std::runtime_error("Error: index not found in merge");

                    unsigned int index1 = firstIndexIntNodesPosNotJSorted[tempIntNodes[intNodesPosJSorted[j] + numDummyNodes].third - numDummyNodes];
                    unsigned int index2 = firstIndexIntNodesPosNotJSorted[tempIntNodes[intNodesPosNotJSorted[i] + numDummyNodes].third - numDummyNodes];

                    if (index1 > index2) // prima era >
                        merged[cont++] = intNodesPosNotJSorted[i++];
                    else
                        merged[cont++] = intNodesPosJSorted[j++];
                }
                else if (v1 < v2)
                    merged[cont++] = intNodesPosNotJSorted[i++];
                else
                    merged[cont++] = intNodesPosJSorted[j++];
            }
        }

        if (j >= intNodesPosJSorted.size())
        {
            for (unsigned int k = i; k < intNodesPosNotJSorted.size(); ++k)
                merged[cont++] = intNodesPosNotJSorted[k];

            flag = false;
        }

        if (i >= intNodesPosNotJSorted.size())
        {
            for (unsigned int k = j; k < intNodesPosJSorted.size(); ++k)
                merged[cont++] = intNodesPosJSorted[k];

            flag = false;
        }

        // update first index merged
        if (firstIndexMerged[merged[cont - 1]] == -1)
            firstIndexMerged[merged[cont - 1]] = cont - 1;
    }

    if (dummyRoot && !firstIt)
    {
        // remove dummy root
        merged.erase(merged.begin());
    }

    /* std::cout << "Merged: " << std::endl;
    for (unsigned int i = 0; i < merged.size(); ++i)
        std::cout << merged[i] << " ";
    std::cout << std::endl; */

    return merged;
}

/**
 * @brief Stable sort the nodes based on the π’s component of its triplets.
 *
 * This function stable sorts the nodes of a labeled tree based on the string obtained by
 * concatenating the labels (π component) on the upward path from the node's parent to the root (the root
 * has empty π component).
 *
 * @param tree The labeled tree to be sorted.
 * @param intNodes A pointer to a vector of Triplets representing the internal nodes of the tree. If nullptr, the internal nodes will be computed using the function `computeIntNodes`.
 * @param dummyRoot A boolean indicating if a dummy root is used.
 * @param firstIt A boolean indicating if this is the first iteration of the pathSort algorithm.
 * @param rem The number of nodes to be removed.
 * @return A vector of sorted node positions.
 */
template <typename T>
std::vector<unsigned int> XBWT<T>::pathSort(const std::vector<Triplet<unsigned int, int, int>> &intNodes, bool dummyRoot, bool firstIt, unsigned int rem, bool j0Encountered, unsigned int maxNumDigits) const
{
    // get number of nodes at level j I {0, 1, 2} mode 3
    long int nodeLevelCounts[] = {0, 0, 0};
    for (const auto &triplet : intNodes)
    {
        ++nodeLevelCounts[triplet.second % 3];
    }

    nodeLevelCounts[0] -= rem;

    // t/3
    long int x = 0;
    if (rem > 0)
    {
        long double realSize = (static_cast<long int>(intNodes.size()) - rem) / 3.;
        x = static_cast<long int>(std::ceil(realSize));
    }
    else
        x = static_cast<long int>(std::ceil(static_cast<long int>(intNodes.size()) / 3.0));

    // compute j value
    short int j = -1;
    for (unsigned int i = 0; i < 3; ++i)
    {
        if (nodeLevelCounts[i] >= x)
        {
            j = i;
            break;
        }
    }

    // std::cout << "j: " << j << std::endl;
    assert(j != -1 && "Error: j value not found");

    unsigned int intNodesPosJSize = j == 0 ? nodeLevelCounts[j] + rem : nodeLevelCounts[j];
    std::vector<unsigned int> intNodesPosJ(intNodesPosJSize, 0);
    std::vector<unsigned int> intNodesPosNotJ(intNodes.size() - intNodesPosJSize, 0);
    for (unsigned int i = 0, k = 0, l = 0; i < intNodes.size(); ++i)
    {
        if (intNodes[i].second % 3 == j)
        {
            intNodesPosJ[k++] = i;
        }
        else
        {
            intNodesPosNotJ[l++] = i;
        }
    }

    // prepemd intNodes with dummy root if needed (2 or 3 based on j level)
    short int numDummyNodes = (j != 0) ? 3 : 2;
    std::vector<Triplet<unsigned int, int, int>> tempIntNodes(intNodes.size() + numDummyNodes, Triplet<unsigned int, int, int>(0, 0, 0));
    if (j != 0)
    {
        tempIntNodes[0] = Triplet<unsigned int, int, int>(0, -3, -1);
        tempIntNodes[1] = Triplet<unsigned int, int, int>(0, -2, 0);
        tempIntNodes[2] = Triplet<unsigned int, int, int>(0, -1, 1);

        for (unsigned int i = 0; i < intNodes.size(); ++i)
        {
            tempIntNodes[i + numDummyNodes] = intNodes[i];
            tempIntNodes[i + numDummyNodes].third += numDummyNodes - 1;
        }
    }
    else
    {
        tempIntNodes[0] = Triplet<unsigned int, int, int>(0, -2, -1);
        tempIntNodes[1] = Triplet<unsigned int, int, int>(0, -1, 0);

        for (unsigned int i = 0; i < intNodes.size(); ++i)
        {
            tempIntNodes[i + numDummyNodes] = intNodes[i];
            tempIntNodes[i + numDummyNodes].third += numDummyNodes - 1;
        }
    }

    std::vector<std::pair<unsigned int, Triplet<std::string, std::string, std::string>>> tripletsToSort(intNodesPosNotJ.size(), std::pair<unsigned int, Triplet<std::string, std::string, std::string>>(0, Triplet<std::string, std::string, std::string>("", "", "")));
    for (unsigned int i = 0; i < intNodesPosNotJ.size(); ++i)
    {
        auto index = intNodesPosNotJ[i] + numDummyNodes;
        tripletsToSort[i].first = intNodesPosNotJ[i];

        tripletsToSort[i].second.first = padLeft(std::to_string(tempIntNodes[tempIntNodes[index].third].first), maxNumDigits, '0');
        index = tempIntNodes[index].third;

        tripletsToSort[i].second.second = padLeft(std::to_string(tempIntNodes[tempIntNodes[index].third].first), maxNumDigits, '0');
        index = tempIntNodes[index].third;

        tripletsToSort[i].second.third = padLeft(std::to_string(tempIntNodes[tempIntNodes[index].third].first), maxNumDigits, '0');
    }

    radixSort(tripletsToSort);

    bool notUnique = false;
    for (unsigned int i = 1; i < tripletsToSort.size(); ++i)
    {
        if (tripletsToSort[i].second.first == tripletsToSort[i - 1].second.first &&
            tripletsToSort[i].second.second == tripletsToSort[i - 1].second.second &&
            tripletsToSort[i].second.third == tripletsToSort[i - 1].second.third)
        {
            notUnique = true;
            break;
        }
    }

    std::vector<unsigned int> intNodesPosNotJSorted(intNodesPosNotJ.size(), 0);
    if (notUnique)
    {
        LabeledTree<unsigned int> cTree;
        if (firstIt)
            cTree.setRoot(contractTree(intNodes, j, &maxNumDigits, &tripletsToSort));
        else
            cTree.setRoot(contractTree(intNodes, j, &maxNumDigits));

        ++rem;
        if (j == 0)
            intNodesPosNotJSorted = pathSort(computeIntNodes<unsigned int>(*cTree.getRoot()), true, false, rem, true, maxNumDigits);
        else
            intNodesPosNotJSorted = pathSort(computeIntNodes<unsigned int>(*cTree.getRoot()), false, false, rem, j0Encountered, maxNumDigits);

        /* std::vector<Triplet<unsigned int, int, int>> localIntNodes = computeIntNodes<unsigned int>(*cTree.getRoot());
        std::vector<std::pair<unsigned int, std::string>> intNodesPaths;
        intNodesPaths.push_back(std::pair<unsigned int, std::string>(0, ""));
        for (unsigned int i = 1; i < localIntNodes.size(); ++i)
        {
            std::ostringstream oss;
            auto parentNode = localIntNodes[localIntNodes[i].third - 1];
            oss << std::setfill('0') << std::setw(maxNumDigits) << localIntNodes[i].first;
            oss << std::setfill('0') << std::setw(maxNumDigits) << parentNode.first;
            do
            {
                if (parentNode.third != 0)
                {
                    parentNode = localIntNodes[parentNode.third - 1];
                    oss << std::setfill('0') << std::setw(maxNumDigits) << parentNode.first;
                }
            } while (parentNode.third != 0);

            intNodesPaths.push_back(std::pair<unsigned int, std::string>(i, oss.str()));
        }

        std::stable_sort(intNodesPaths.begin(), intNodesPaths.end(), [&maxNumDigits](const std::pair<unsigned int, std::string> &a, const std::pair<unsigned int, std::string> &b)
                         {
            if (a.second.substr(0, maxNumDigits) == b.second.substr(0, maxNumDigits))
            {
                std::string aStr = a.second.substr(maxNumDigits);
                std::string bStr = b.second.substr(maxNumDigits);

                return aStr < bStr; 
            } 
            else 
            {
                return a.second.substr(0, maxNumDigits) < b.second.substr(0, maxNumDigits);
            } });

        std::vector<unsigned int> intNodesPosSorted(intNodesPaths.size(), 0);
        std::cout << "Int nodes pos sorted correctly: " << std::endl;
        if (j == 0)
            intNodesPaths.erase(intNodesPaths.begin());
        for (unsigned int i = 0; i < intNodesPaths.size(); ++i)
        {
            // intNodesPosNotJSorted[i] = intNodesPaths[i].first;
            std::cout << intNodesPaths[i].first << " " << intNodesPaths[i].second << std::endl;
        }
        std::cout << std::endl;

        std::cout << "intNodesPosNotJ size: " << intNodesPosNotJ.size() << std::endl;
        std::cout << "cTree size: " << cTree.getNodes().size() << std::endl;
        std::cout << "cTree: " << cTree.toString() << std::endl; */

        if (j == 0)
        {
            for (auto &el : intNodesPosNotJSorted)
                --el;
        }

        /* std::cout << "Int nodes pos not j sorted corr: " << std::endl;
        for (unsigned int i = 0; i < intNodesPosNotJSorted.size(); ++i)
            std::cout << intNodesPosNotJSorted[i] << " ";
        std::cout << std::endl; */
    }
    else if (!notUnique && firstIt)
    {
        std::vector<long int> indexIntNodesPosNotJ(intNodes.size(), -1);
        for (unsigned int i = 0; i < intNodesPosNotJ.size(); ++i)
        {
            if (indexIntNodesPosNotJ[intNodesPosNotJ[i]] == -1)
                indexIntNodesPosNotJ[intNodesPosNotJ[i]] = i;
        }

        for (unsigned int i = 0; i < intNodesPosNotJSorted.size(); ++i)
        {
            assert(indexIntNodesPosNotJ[tripletsToSort[i].first] != -1 && "Error: index not found");
            intNodesPosNotJSorted[i] = indexIntNodesPosNotJ[tripletsToSort[i].first];
        }
    }
    else
    {
        // stable sort intNodesPosNotJ by intNodes.first
        intNodesPosNotJSorted = intNodesPosNotJ;
        std::stable_sort(intNodesPosNotJSorted.begin(), intNodesPosNotJSorted.end(), [&intNodes](unsigned int a, unsigned int b)
                         { return intNodes[a].first < intNodes[b].first; });
    }

    assert(intNodesPosNotJSorted.size() == intNodesPosNotJ.size());

    std::vector<unsigned int> firstIndex(intNodes.size(), 0);
    std::vector<bool> indexFound(intNodes.size(), false);
    for (unsigned int i = 0; i < intNodesPosNotJSorted.size(); ++i)
    {
        if (notUnique || (!notUnique && firstIt))
            intNodesPosNotJSorted[i] = intNodesPosNotJ[intNodesPosNotJSorted[i]];

        if (!indexFound[intNodesPosNotJSorted[i]])
        {
            firstIndex[intNodesPosNotJSorted[i]] = i;
            indexFound[intNodesPosNotJSorted[i]] = true;
        }
    }

    /* std::cout << "Int nodes pos not j sorted: " << std::endl;
    for (unsigned int i = 0; i < intNodesPosNotJSorted.size(); ++i)
        std::cout << intNodesPosNotJSorted[i] << " ";
    std::cout << std::endl; */

    std::vector<Triplet<unsigned int, unsigned int, unsigned int>> pairIntNodesPosJSorted(intNodesPosJ.size(), Triplet<unsigned int, unsigned int, unsigned int>(0, 0, 0));
    short int c = 0;
    if (!firstIt)
    {
        if (j == 0)
        {
            pairIntNodesPosJSorted[0].first = tempIntNodes[intNodesPosJ[0] + numDummyNodes].first;
            pairIntNodesPosJSorted[0].second = 0;
            pairIntNodesPosJSorted[0].third = intNodesPosJ[0];
            c = 1;
        }

        for (unsigned int i = c; i < intNodesPosJ.size(); ++i)
        {
            assert(indexFound[tempIntNodes[intNodesPosJ[i] + numDummyNodes].third - numDummyNodes] && "Error: index not found");

            pairIntNodesPosJSorted[i].first = tempIntNodes[intNodesPosJ[i] + numDummyNodes].first;
            pairIntNodesPosJSorted[i].second = firstIndex[tempIntNodes[intNodesPosJ[i] + numDummyNodes].third - numDummyNodes];
            pairIntNodesPosJSorted[i].third = intNodesPosJ[i];
        }
    }
    else
    {
        if (j == 0)
        {
            pairIntNodesPosJSorted[0].first = tempIntNodes[tempIntNodes[intNodesPosJ[0] + numDummyNodes].third].first;
            pairIntNodesPosJSorted[0].second = 0;
            pairIntNodesPosJSorted[0].third = intNodesPosJ[0];
            c = 1;
        }

        for (unsigned int i = c; i < intNodesPosJ.size(); ++i)
        {
            assert(indexFound[tempIntNodes[intNodesPosJ[i] + numDummyNodes].third - numDummyNodes] && "Error: index not found");

            pairIntNodesPosJSorted[i].first = tempIntNodes[tempIntNodes[intNodesPosJ[i] + numDummyNodes].third].first;
            pairIntNodesPosJSorted[i].second = firstIndex[tempIntNodes[intNodesPosJ[i] + numDummyNodes].third - numDummyNodes];
            pairIntNodesPosJSorted[i].third = intNodesPosJ[i];
        }
    }

    // sort by first and second element
    std::stable_sort(pairIntNodesPosJSorted.begin(), pairIntNodesPosJSorted.end(), [](const Triplet<unsigned int, unsigned int, unsigned int> &a, const Triplet<unsigned int, unsigned int, unsigned int> &b)
                     { return a.first < b.first || (a.first == b.first && a.second < b.second); });

    std::vector<unsigned int> intNodesPosJSorted(intNodesPosJ.size(), 0);
    for (unsigned int i = 0; i < pairIntNodesPosJSorted.size(); ++i)
    {
        intNodesPosJSorted[i] = pairIntNodesPosJSorted[i].third;
    }

    /* std::cout << "Int nodes pos j sorted: " << std::endl;
    for (unsigned int i = 0; i < intNodesPosJSorted.size(); ++i)
        std::cout << intNodesPosJSorted[i] << " ";
    std::cout << std::endl;

    std::cout << "___________________________" << std::endl; */

    return pathSortMerge(intNodesPosNotJSorted, firstIndex, indexFound, intNodesPosJSorted, tempIntNodes, numDummyNodes, j, dummyRoot, firstIt, j0Encountered);
}

/**
 * @brief Creates the XBWT of the given labeled tree.
 *
 * This function creates the XBWT structure from a labeled tree. It computes the internal nodes,
 * sorts them using the pathSort algorithm, and constructs various bit vectors and wavelet trees
 * required for the XBWT representation. Optionally, it can print detailed information about the
 * construction process if verbose mode is enabled.
 *
 * @param tree The labeled tree used to construct the XBWT.
 * @param usePathSort If true, uses the pathSort algorithm to sort the internal nodes.
 * @param verbose If true, enables verbose mode for debugging.
 * @param intNodesPosSorted An optional vector to store the sorted positions of the internal nodes.
 */
template <typename T>
void XBWT<T>::createXBWT(const LabeledTree<T> &tree, bool usePathSort, bool verbose, std::vector<unsigned int> *intNodesPosSorted)
{
    std::vector<unsigned int> posIntNodesSorted;
    std::vector<Triplet<unsigned int, int, int>> intNodes = computeIntNodes(*tree.getRoot(), true);
    if (usePathSort)
        posIntNodesSorted = pathSort(intNodes, false, true, 0, false, pImpl->maxNumDigits);
    else
        posIntNodesSorted = upwardStableSortConstruction(tree, &intNodes);

    if (verbose)
    {
        std::cout << "Sorted pos: " << std::endl;
        for (unsigned int i : posIntNodesSorted)
            std::cout << i << " ";
        std::cout << std::endl;
    }

    if (intNodesPosSorted)
    {
        *intNodesPosSorted = posIntNodesSorted;
    }

    std::vector<unsigned int> posLast(intNodes.size(), 0);
    std::vector<bool> _SLast(intNodes.size(), false);
    std::vector<bool> _SAlphaBit(intNodes.size(), true);
    for (unsigned int i = 1; i < intNodes.size(); ++i)
    {
        posLast[intNodes[i].third - 1] = i;
        _SAlphaBit[intNodes[i].third - 1] = false;
    }

    for (unsigned int i : posLast)
    {
        if (i != 0)
            _SLast[i] = true;
    }

    sdsl::bit_vector SLast(intNodes.size());
    sdsl::int_vector<> tempSAlpha(intNodes.size());
    sdsl::bit_vector SAlphaBit(intNodes.size()); // SAlphaBit[i] = 1 iff the corresponding label is a leaf label
    sdsl::bit_vector A(intNodes.size());         // A[1] = 1, A[j] = 1 iff the first symbol of Sπ[j] differs from the first symbol of Sπ[j − 1]
    sdsl::bit_vector SigmaN(pImpl->cardSigma, 0);
    unsigned int cont = 0;
    unsigned int prev_label = 0;
    for (unsigned int i : posIntNodesSorted)
    {
        SLast[cont] = _SLast[i];
        tempSAlpha[cont] = intNodes[i].first;
        SAlphaBit[cont] = _SAlphaBit[i];
        A[cont] = (cont == 0) ? 0 : (intNodes[intNodes[i].third - 1].first != prev_label); // TODO: in paper A[0] = 1 but it should be 0
        prev_label = (cont == 0) ? 0 : intNodes[intNodes[i].third - 1].first;

        if (SigmaN[intNodes[i].first - 1] == 0 && SAlphaBit[cont] == 0)
            SigmaN[intNodes[i].first - 1] = 1;

        ++cont;
    }

    // Compress the bit vectors
    pImpl->SLastCompressed = sdsl::rrr_vector<>(SLast);
    pImpl->SAlphaBitCompressed = sdsl::rrr_vector<>(SAlphaBit);
    pImpl->ACompressed = sdsl::rrr_vector<>(A);
    pImpl->SigmaNCompressed = sdsl::rrr_vector<>(SigmaN);

    // add rank and select support
    pImpl->SLastCompressedRank = sdsl::rrr_vector<>::rank_1_type(&pImpl->SLastCompressed);
    pImpl->SLastCompressedSelect = sdsl::rrr_vector<>::select_1_type(&pImpl->SLastCompressed);
    pImpl->ACompressedRank = sdsl::rrr_vector<>::rank_1_type(&pImpl->ACompressed);
    pImpl->ACompressedSelect = sdsl::rrr_vector<>::select_1_type(&pImpl->ACompressed);
    pImpl->SigmaNCompressedRank = sdsl::rrr_vector<>::rank_1_type(&pImpl->SigmaNCompressed);

    // Create the wavelet tree
    sdsl::construct_im(pImpl->SAlphaCompressed, tempSAlpha);

    if (verbose)
    {
        // print
        std::cout << "(Index, SLast, SAlpha, SAlphaBit, A) " << std::endl;
        for (unsigned int i = 0; i < intNodes.size(); ++i)
        {
            std::cout << "(" << i << ", " << (pImpl->SLastCompressed[i] ? "1" : "0") << ", " << pImpl->SAlphaCompressed[i] << ", " << (pImpl->SAlphaBitCompressed[i] ? "1" : "0") << ", " << pImpl->ACompressed[i] << ") " << std::endl;
        }

        // compress int_vectors
        sdsl::vlc_vector<> vv(tempSAlpha);
        std::cout << "SAlpha vv size: " << sdsl::size_in_bytes(vv) << " B" << std::endl;

        // build wt_int with non compressed bit vectors
        sdsl::wt_int<> wt;
        sdsl::construct_im(wt, tempSAlpha);
        std::cout << "SAlpha wt size: " << sdsl::size_in_bytes(wt) << " B" << std::endl;        

        std::cout << "SLast size: " << sdsl::size_in_bytes(SLast) << " B" << std::endl;
        std::cout << "SLast compressed size: " << sdsl::size_in_bytes(pImpl->SLastCompressed) << " B" << std::endl;
        std::cout << "SAlpha size: " << sdsl::size_in_bytes(tempSAlpha) << " B" << std::endl;
        std::cout << "SAlpha compressed size: " << sdsl::size_in_bytes(pImpl->SAlphaCompressed) << " B" << std::endl;
        std::cout << "SAlphaBit size: " << sdsl::size_in_bytes(SAlphaBit) << " B" << std::endl;
        std::cout << "SAlphaBit compressed size: " << sdsl::size_in_bytes(pImpl->SAlphaBitCompressed) << " B" << std::endl;
        std::cout << "A size: " << sdsl::size_in_bytes(A) << " B" << std::endl;
        std::cout << "A compressed size: " << sdsl::size_in_bytes(pImpl->ACompressed) << " B" << std::endl;
        std::cout << "SigmaN size: " << sdsl::size_in_bytes(SigmaN) << " B" << std::endl;
        std::cout << "SigmaN compressed size: " << sdsl::size_in_bytes(pImpl->SigmaNCompressed) << " B" << std::endl;
    
        // Compute and print compression ratio
        double originalSize = sdsl::size_in_bytes(wt) + sdsl::size_in_bytes(SAlphaBit);
        double compressedSize = sdsl::size_in_bytes(pImpl->SLastCompressed) + sdsl::size_in_bytes(pImpl->SAlphaCompressed) + sdsl::size_in_bytes(pImpl->SAlphaBitCompressed) + sdsl::size_in_bytes(pImpl->ACompressed);
        double compressionRatio = (1 - (compressedSize / originalSize)) * 100;

        std::cout << "Original size: " << originalSize << " B" << std::endl;
        std::cout << "Compressed size: " << compressedSize << " B" << std::endl;
        std::cout << "Space saving: " << compressionRatio << "%" << std::endl;
    }
}

/**
 * @brief Builds the F array for the XBWT structure.
 *
 * Compute array F such that F[i] = j iff Sπ [ j] is the first entry of S prefixed by i.
 *
 * @return A vector of unsigned integers representing the F array.
 */
template <typename T>
std::vector<unsigned int> XBWT<T>::buildF() const
{
    std::vector<unsigned int> C(pImpl->cardSigma, 0);
    // count occurrences of each symbol only for internal nodes
    for (unsigned int i = 0; i < pImpl->SLastCompressed.size(); ++i)
        if (pImpl->SAlphaBitCompressed[i] == 0)
            ++C[pImpl->SAlphaCompressed[i] - 1];

    std::vector<unsigned int> F(pImpl->cardSigmaN, 0);
    F[0] = 1;
    // iter internal labels cardinality
    for (unsigned int i = 0; i < pImpl->cardSigmaN - 1; ++i)
    {
        unsigned int s = 0, j = F[i];
        while (s != C[i])
        {
            if (pImpl->SLastCompressed[j++] == 1)
                ++s;
        }

        F[i + 1] = j;
    }

    return F;
}

/**
 * @brief Builds the J array for the XBWT structure.
 *
 * Compute array J such that J[i] = j if S[j] is the first child of S[i],
 * and J[i] = −1 if S[i] is a leaf.
 *
 * @return A vector of long integers representing the J array.
 */
template <typename T>
std::vector<long int> XBWT<T>::buildJ() const
{
    std::vector<long int> J(pImpl->SLastCompressed.size(), 0);
    std::vector<long int> FMap(pImpl->cardSigmaN, -1);
    for (unsigned int i = 0; i < J.size(); ++i)
    {
        if (pImpl->SAlphaBitCompressed[i] == 1)
            J[i] = -1;
        else
        {
            unsigned int z = 0, internalNodeIndex = pImpl->SigmaNCompressedRank(pImpl->SAlphaCompressed[i]);
            if (FMap[internalNodeIndex - 1] == -1)
                z = J[i] = FMap[internalNodeIndex - 1] = pImpl->ACompressedSelect(internalNodeIndex);
            else
                z = J[i] = FMap[internalNodeIndex - 1];

            while (pImpl->SLastCompressed[z] != 1)
                ++z;

            FMap[internalNodeIndex - 1] = z + 1;
        }
    }

    return J;
}

/**
 * @brief Rebuilds the labeled tree from the XBWT structure.
 *
 * This function rebuilds the labeled tree from the XBWT structure. It uses the `J` array
 * and reconstructs the tree by traversing the compressed vectors and creating the nodes accordingly.
 *
 * @return A labeled tree reconstructed from the XBWT structure.
 */
template <typename T>
LabeledTree<T> XBWT<T>::rebuildTree() const
{
    /* auto F = buildF(); // TODO: usare direttamente ACompressed
    for (unsigned int i = 0; i < F.size(); ++i)
    {
        std::cout << "F[" << i << "] = " << F[i] << ", ";
    }
    std::cout << std::endl; */

    auto J = buildJ();
    /* for (unsigned int i = 0; i < J.size(); ++i)
    {
        std::cout << "J[" << i << "] = " << J[i] << ", ";
    }
    std::cout << std::endl; */

    Node<T> *root = new Node<T>(pImpl->alphabetMapInv[pImpl->SAlphaCompressed[0]]);
    std::stack<std::pair<unsigned int, Node<T> *>> Q;
    Q.push({0, root});
    while (!Q.empty())
    {
        auto [i, u] = Q.top();
        Q.pop();
        long int j = J[i];
        if (j == -1)
            continue;

        // find first j' >= j such that SLast[j'] = 1
        unsigned int j1 = j;
        while (j1 < J.size() && pImpl->SLastCompressed[j1] != 1)
            ++j1;

        for (unsigned int h = j1; h >= j; --h)
        {
            Node<T> *v = new Node<T>(pImpl->alphabetMapInv[pImpl->SAlphaCompressed[h]]);
            u->prependChild(v);
            Q.push({h, v});
        }
    }

    return LabeledTree<T>(root);
}

/**
 * @brief Retrieves the children of a node from the XBWT structure.
 *
 * This function retrieves the children of a node in the XBWT structure. It uses the compressed
 * vectors to determine the range of positions corresponding to the children of the given node.
 *
 * @param i The index of the node whose children are to be retrieved.
 * @return A pair of long integers representing the range [first, last] of the children.
 *         If the node has no children, returns {-1, -1}.
 * @throws std::runtime_error If the index is out of bounds.
 */
template <typename T>
std::pair<long int, long int> XBWT<T>::getChildren(unsigned int i) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    if (pImpl->SAlphaBitCompressed[i] == 1)
        return {-1, -1};

    unsigned int c = pImpl->SAlphaCompressed[i];
    unsigned int r = pImpl->SAlphaCompressed.rank(i + 1, c);
    unsigned int cIndex = pImpl->SigmaNCompressedRank(c);
    unsigned int y = pImpl->ACompressedSelect(cIndex); // F[c]
    unsigned int z = pImpl->SLastCompressedRank(y);

    long int first = 0;
    if (static_cast<long int>(z + r) - 1 < 1) // TODO: their pseudocode seems to be wrong for i = 0...
        first = 1;
    else
        first = pImpl->SLastCompressedSelect(z + r - 1) + 1;

    long int last = pImpl->SLastCompressedSelect(z + r);

    return {first, last};
}

/**
 * @brief Retrieves the k-th child of a node in the XBWT structure.
 *
 * This function retrieves the k-th child of a node in the XBWT structure. It uses the `getChildren`
 * function to determine the range of positions corresponding to the children of the given node
 * and returns the k-th child within that range.
 *
 * @param i The index of the node whose k-th child is to be retrieved.
 * @param k The rank of the child to be retrieved (1-based index).
 * @return The index of the k-th child, or -1 if the node has fewer than k children or if the node has no children.
 * @throws std::runtime_error If the index is out of bounds or if k is less than 1.
 */
template <typename T>
long int XBWT<T>::getRankedChild(unsigned int i, unsigned int k) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    if (k < 1)
        throw std::runtime_error("Error: k must be greater than 0");

    auto [first, last] = getChildren(i);
    if (first == -1 || last == -1 || (k > last - first + 1))
        return -1;
    else
        return first + k - 1;
}

/**
 * @brief Retrieves the k-th child of a node with a specific label in the XBWT structure.
 *
 * This function retrieves the k-th child of a node with a specific label in the XBWT structure.
 * It uses the `getChildren` function to determine the range of positions corresponding to the children
 * of the given node and returns the k-th child with the specified label within that range.
 *
 * @param i The index of the node whose k-th child with the specified label is to be retrieved.
 * @param c The label of the child to be retrieved.
 * @param k The rank of the child to be retrieved (1-based index).
 * @return The index of the k-th child with the specified label, or -1 if the node has fewer than k children with that label or if the node has no children.
 * @throws std::runtime_error If the index is out of bounds or if k is less than 1.
 */
template <typename T>
long int XBWT<T>::getCharRankedChild(unsigned int i, T label, unsigned int k) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    if (k < 1)
        throw std::runtime_error("Error: k must be greater than 0");

    auto c_it = pImpl->alphabetMap.find(label);
    if (c_it == pImpl->alphabetMap.end())
        return -1;

    auto [first, last] = getChildren(i);
    if (first == -1 or last == -1)
        return -1;

    unsigned int c = c_it->second;
    unsigned int y1 = pImpl->SAlphaCompressed.rank(first, c);
    unsigned int y2 = pImpl->SAlphaCompressed.rank(last + 1, c);
    if (k > y2 - y1)
        return -1;
    else
        return pImpl->SAlphaCompressed.select(y1 + k, c);
}

/**
 * @brief Return the number of children of a node i in the XBWT structure.
 *
 * @param i Index of the node.
 * @return Number of children of the node i.
 */
template <typename T>
unsigned int XBWT<T>::getDegree(unsigned int i) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    if (pImpl->SAlphaBitCompressed[i] == 1)
        return 0;

    auto [first, last] = getChildren(i);
    if (first == -1 || last == -1)
        return 0;

    return last - first + 1;
}

/**
 * @brief Return the number of children of a node i with label c in the XBWT structure.
 *
 * @param i Index of the node.
 * @param c Label of the children.
 * @return Nnumber of children of the node i with label c.
 */
template <typename T>
unsigned int XBWT<T>::getCharDegree(unsigned int i, T label) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    auto c_it = pImpl->alphabetMap.find(label);
    if (c_it == pImpl->alphabetMap.end())
        return 0;

    auto [first, last] = getChildren(i);
    if (first == -1 or last == -1)
        return 0;

    unsigned int c = c_it->second;
    unsigned int y1 = pImpl->SAlphaCompressed.rank(first, c);
    unsigned int y2 = pImpl->SAlphaCompressed.rank(last + 1, c);
    return y2 - y1;
}

/**
 * @brief Retrieves the parent of a node in the XBWT structure.
 *
 * This function retrieves the parent of a node in the XBWT structure. It uses the compressed
 * vectors to determine the position of the parent node.
 *
 * @param i The index of the node whose parent is to be retrieved.
 * @return The index of the parent node, or -1 if the node is the root or if the index is out of bounds.
 * @throws std::runtime_error If the index is out of bounds.
 */
template <typename T>
long int XBWT<T>::getParent(unsigned int i) const
{
    if (i >= pImpl->SLastCompressed.size() || i < 0)
        throw std::runtime_error("Error: index out of bounds");

    if (i == 0)
        return -1;

    unsigned int c = pImpl->ACompressedRank(i + 1);
    unsigned int y = pImpl->ACompressedSelect(pImpl->SigmaNCompressedRank(c));
    unsigned int k = pImpl->SLastCompressedRank(i) - pImpl->SLastCompressedRank(y);
    return pImpl->SAlphaCompressed.select(k + 1, c);
}

/**
 * @brief Return the subtree rooted at node i in the XBWT structure.
 *
 * This function returns the subtree rooted at node i in the XBWT structure. The order parameter
 * specifies the type of traversal to be used to retrieve the subtree. Tse 0 for preorder, 1 for
 * postorder, and 2 for inorder traversal.
 *
 * @param i Index of the root node.
 * @return A vector of unsigned integers representing the labels of nodes in the subtree rooted at node i.
 */
template <typename T>
std::vector<T> XBWT<T>::getSubtree(unsigned int i, unsigned int order) const
{
    std::vector<T> subtree;

    if (order == 0)
    {
        // Pre-order traversal
        std::stack<unsigned int> stack;
        stack.push(i);

        while (!stack.empty())
        {
            unsigned int node = stack.top();
            stack.pop();
            subtree.push_back(pImpl->alphabetMapInv[pImpl->SAlphaCompressed[node]]);

            auto [first, last] = getChildren(node);
            for (unsigned int j = last; j >= first && j != static_cast<unsigned int>(-1); --j)
                stack.push(j);
        }
    }
    else if (order == 1)
    {
        // Post-order traversal
        std::stack<unsigned int> stack;
        std::stack<unsigned int> output;
        stack.push(i);

        while (!stack.empty())
        {
            unsigned int node = stack.top();
            stack.pop();
            output.push(node);

            auto [first, last] = getChildren(node);
            for (unsigned int j = first; j <= last && j != static_cast<unsigned int>(-1); ++j)
                stack.push(j);
        }

        while (!output.empty())
        {
            subtree.push_back(pImpl->alphabetMapInv[pImpl->SAlphaCompressed[output.top()]]);
            output.pop();
        }
    }
    else if (order == 2)
    {
        // In-order traversal
        std::stack<std::pair<unsigned int, bool>> stack;
        stack.push({i, false});

        while (!stack.empty())
        {
            auto [node, visited] = stack.top();
            stack.pop();

            if (visited)
            {
                subtree.push_back(pImpl->alphabetMapInv[pImpl->SAlphaCompressed[node]]);
            }
            else
            {
                auto [first, last] = getChildren(node);
                if (first != static_cast<unsigned int>(-1))
                {
                    for (unsigned int j = last; j >= first && j != static_cast<unsigned int>(-1); --j)
                    {
                        stack.push({j, false});
                    }
                }
                stack.push({node, true});
                if (first != static_cast<unsigned int>(-1))
                {
                    stack.push({first, false});
                }
            }
        }
    }
    else
        throw std::runtime_error("Error: invalid order");

    return subtree;
}

/**
 * @brief Searches for the rannge of nodes whose upward path is prefixed by a given string reversed.
 *
 * This function searches for a subpath in the XBWT structure. It uses the compressed vectors
 * to determine the range of positions corresponding to the nodes whose upward path is prefixed
 * by a given string reversed.
 *
 * @param subPath The subpath to be searched for, represented as a string of labels.
 * @return A pair of long integers representing the range [first, last] of the nodes that match the subpath.
 *         If the subpath is not found, returns {-1, -1}.
 */
template <typename T>
std::pair<long int, long int> XBWT<T>::subPathSearch(const std::vector<T> &path) const
{
    auto label_it = pImpl->alphabetMap.find(path[0]);
    if (label_it == pImpl->alphabetMap.end())
        return {-1, -1};

    long int first = pImpl->ACompressedSelect(pImpl->SigmaNCompressedRank(label_it->second));
    long int last = pImpl->ACompressedSelect(pImpl->SigmaNCompressedRank(label_it->second + 1)) - 1;
    if (first > last)
        return {-1, -1};

    for (unsigned int i = 1; i < path.size(); ++i)
    {
        label_it = pImpl->alphabetMap.find(path[i]);
        if (label_it == pImpl->alphabetMap.end())
            return {-1, -1};

        unsigned int label = label_it->second;
        unsigned int k1 = pImpl->SAlphaCompressed.rank(first, label);
        unsigned int k2 = pImpl->SAlphaCompressed.rank(last + 1, label);

        unsigned int check_for_z1 = pImpl->SAlphaCompressed.rank(pImpl->SAlphaCompressed.size(), label);
        if (k1 + 1 > check_for_z1)
            return {-1, -1};

        unsigned int z1 = pImpl->SAlphaCompressed.select(k1 + 1, label);
        unsigned int z2 = pImpl->SAlphaCompressed.select(k2, label);
        if (z1 > z2)
            return {-1, -1};

        first = getRankedChild(z1, 1);
        last = getChildren(z2).second;
    }

    return {first, last};
}

/**
 * @brief Return the label associated with the node at index i in the XBWT structure.
 *
 * @param i Node index in the XBWT structure
 * @return unsigned int Label associated with the node at index i
 */
template <typename T>
T XBWT<T>::getNodeLabel(unsigned int i) const
{
    return pImpl->alphabetMapInv[pImpl->SAlphaCompressed[i]];
}

/**
 * @brief Return the upward path (from root to parent node) of the node at index i
 * in the XBWT structure.
 *
 * @param i Index of the node in the XBWT structure
 * @return std::string Tpward path from root to the node i parent node
 */
template <typename T>
std::vector<T> XBWT<T>::getUpwardPath(unsigned int i) const
{
    std::vector<T> path;
    auto parent = getParent(i);
    while (parent != -1)
    {
        path.push_back(getNodeLabel(parent));
        parent = getParent(parent);
    }

    // return reversed vector
    std::reverse(path.begin(), path.end());
    return path;
}

template <typename T>
unsigned int XBWT<T>::getCardSigma() const
{
    return pImpl->cardSigma;
}

#endif // XBWT_HPP
