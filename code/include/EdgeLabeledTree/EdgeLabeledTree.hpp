#ifndef EDGELABELEDTREE_HPP
#define EDGELABELEDTREE_HPP

#include <vector>
#include <memory>
#include <string>
#include <stack>
#include <sstream>
#include <iostream>
#include <functional>
#include <algorithm>

/**
 * @brief A template class representing an edge in a labeled tree.
 * 
 * @tparam LabelType The type of the label for each edge.
 */
template <typename LabelType>
class Edge
{
private:
    LabelType m_label;                                      ///< The label of the edge.
    std::shared_ptr<Node<LabelType>> m_target;              ///< The target node of the edge.

public:
    /**
     * @brief Construct a new Edge object.
     *
     * @param label The label of the edge.
     * @param target The target node of the edge.
     */
    Edge(LabelType label, std::shared_ptr<Node<LabelType>> target) : m_label(label), m_target(target) {}

    /**
     * @brief Gets the label of the edge.
     *
     * @return The label of the edge.
     */
    LabelType get_label() const
    {
        return m_label;
    }

    /**
     * @brief Gets the target node of the edge.
     *
     * @return A shared pointer to the target node.
     */
    std::shared_ptr<Node<LabelType>> get_target() const
    {
        return m_target;
    }

    /**
     * @brief Sets the target node of the edge.
     *
     * @param target A shared pointer to the new target node.
     */
    void set_target(std::shared_ptr<Node<LabelType>> target)
    {
        m_target = target;
    }
};

/**
 * @brief A template class representing a node in an edge-labeled tree.
 * 
 * @tparam LabelType The type of the label for each edge.
 */
template <typename LabelType>
class Node : public std::enable_shared_from_this<Node<LabelType>>
{
private:
    std::vector<Edge<LabelType>> m_edges;                   ///< A vector of edges to children nodes.
    std::weak_ptr<Node> m_parent;                           ///< A weak pointer to the parent of the node.

public:
    /**
     * @brief Construct a new Node object.
     *
     * @param prnt A weak pointer to the parent node. Defaults to empty weak_ptr.
     * @param edges A vector of edges to children nodes. Defaults to an empty vector.
     */
    Node(std::weak_ptr<Node> prnt = {}, std::vector<Edge<LabelType>> edges = {}) : m_edges(edges), m_parent(prnt) {}

    /**
     * @brief Checks if the node is a root.
     *
     * @return true if the node is a root, false otherwise.
     */
    bool is_root() const
    {
        return m_parent.expired();
    }

    /**
     * @brief Checks if the node is a leaf.
     *
     * @return true if the node is a leaf, false otherwise.
     */
    bool is_leaf() const
    {
        return m_edges.empty();
    }

    /**
     * @brief Gets the level of the node in the tree.
     *
     * @return The level of the node.
     */
    unsigned int get_level() const
    {
        unsigned int level = 0;
        auto current = m_parent.lock();
        while (current)
        {
            level++;
            current = current->m_parent.lock();
        }
        return level;
    }

    /**
     * @brief Gets the depth of the subtree rooted at this node.
     *
     * @return The depth of the subtree.
     */
    unsigned int get_depth() const
    {
        unsigned int maxDepth = 0;
        for (const auto &edge : m_edges)
        {
            maxDepth = std::max(maxDepth, edge.get_target()->get_depth());
        }
        return maxDepth + 1;
    }

    /**
     * @brief Gets the edges from this node.
     *
     * @return A const reference to the vector of edges.
     */
    const std::vector<Edge<LabelType>> &get_edges() const
    {
        return m_edges;
    }

    /**
     * @brief Gets the children of the node.
     *
     * @return A vector of shared pointers to the children nodes.
     */
    std::vector<std::shared_ptr<Node>> get_children() const
    {
        std::vector<std::shared_ptr<Node>> children;
        for (const auto &edge : m_edges)
        {
            children.push_back(edge.get_target());
        }
        return children;
    }

    /**
     * @brief Gets the parent of the node.
     *
     * @return A shared pointer to the parent node.
     */
    std::shared_ptr<Node> get_parent() const
    {
        return m_parent.lock();
    }

    /**
     * @brief Checks if the node is the rightmost child of its parent.
     *
     * @return true if the node is the rightmost child, false otherwise.
     */
    bool is_rightmost() const
    {
        auto parent = m_parent.lock();
        if (parent)
        {
            const auto &edges = parent->m_edges;
            if (!edges.empty())
            {
                return edges.back().get_target().get() == this;
            }
        }
        return false;
    }

    /**
     * @brief Adds a new child to the node with a given edge label.
     *
     * @param edgeLabel The label of the edge to the new child node.
     */
    void push_back_child(LabelType edgeLabel)
    {
        auto child = std::make_shared<Node>(std::weak_ptr<Node>(this->shared_from_this()));
        m_edges.emplace_back(edgeLabel, child);
    }

    /**
     * @brief Adds an existing node as a child to this node with a given edge label.
     *
     * @param child A shared pointer to the node to be added as a child.
     * @param edgeLabel The label of the edge to the child node.
     */
    void push_back_child(std::shared_ptr<Node> child, LabelType edgeLabel)
    {
        child->m_parent = std::weak_ptr<Node>(this->shared_from_this());
        m_edges.emplace_back(edgeLabel, child);
    }

    /**
     * @brief Prepends an existing node as a child to this node with a given edge label.
     *
     * @param child A shared pointer to the node to be prepended as a child.
     * @param edgeLabel The label of the edge to the child node.
     */
    void prepend_child(std::shared_ptr<Node> child, LabelType edgeLabel)
    {
        child->m_parent = std::weak_ptr<Node>(this->shared_from_this());
        m_edges.insert(m_edges.begin(), Edge<LabelType>(edgeLabel, child));
    }

    /**
     * @brief Clears all children of this node.
     */
    void clear_children()
    {
        m_edges.clear();
    }

    /**
     * @brief Gets the edge label to a specific child node.
     *
     * @param child A shared pointer to the child node.
     * @return The label of the edge to the child, or default LabelType if not found.
     */
    LabelType get_edge_label_to_child(std::shared_ptr<Node> child) const
    {
        for (const auto &edge : m_edges)
        {
            if (edge.get_target() == child)
            {
                return edge.get_label();
            }
        }
        return LabelType{}; // Return default value if not found
    }

    /**
     * @brief Finds a child node by edge label.
     *
     * @param edgeLabel The label of the edge to search for.
     * @return A shared pointer to the child node, or nullptr if not found.
     */
    std::shared_ptr<Node> find_child_by_edge_label(LabelType edgeLabel) const
    {
        for (const auto &edge : m_edges)
        {
            if (edge.get_label() == edgeLabel)
            {
                return edge.get_target();
            }
        }
        return nullptr;
    }

    /**
     * @brief Gets all edge labels from this node.
     *
     * @return A vector of edge labels.
     */
    std::vector<LabelType> get_edge_labels() const
    {
        std::vector<LabelType> labels;
        for (const auto &edge : m_edges)
        {
            labels.push_back(edge.get_label());
        }
        return labels;
    }

    /**
     * @brief Gets the label path from this node to the root.
     * The path is ordered from the node to the root (excluding the root since it has no incoming edge).
     *
     * @return A vector of labels representing the path from this node to root.
     */
    std::vector<LabelType> get_path_to_root() const
    {
        std::vector<LabelType> path;
        auto current = this->shared_from_this();
        auto parent = current->get_parent();
        
        while (parent)
        {
            // Find the edge label from parent to current
            LabelType edgeLabel = parent->get_edge_label_to_child(current);
            path.push_back(edgeLabel);
            current = parent;
            parent = current->get_parent();
        }
        
        return path;
    }
};

/**
 * @brief A template class representing a labeled tree.
 * 
 * @tparam LabelType The type of the label for each node.
 */
template <typename LabelType>
class EdgeLabeledTree
{
private:
    std::shared_ptr<Node<LabelType>> root; ///< A shared pointer to the root node of the tree.

public:
    /**
     * @brief Construct a new Labeled Tree object.
     * Creates an empty tree with no root node.
     */
    EdgeLabeledTree() : root(nullptr) {}

    /**
     * @brief Construct a new Labeled Tree object with a root node.
     * Note: In edge-labeled trees, the root node has no incoming edge and thus no label.
     */
    EdgeLabeledTree(bool createRoot)
    {
        if (createRoot)
        {
            root = std::make_shared<Node<LabelType>>();
        }
        else
        {
            root = nullptr;
        }
    }

    /**
     * @brief Construct a new Labeled Tree object from a string representation.
     *
     * @param str The string representation of the tree.
     * @param strToLabel A function to convert a string to a label.
     */
    EdgeLabeledTree(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        root = fromString(str, strToLabel);
    }

    /**
     * @brief Construct a new Labeled Tree object from a root node.
     *
     * @param root A shared pointer to the root node.
     */
    EdgeLabeledTree(std::shared_ptr<Node<LabelType>> root) : root(root) {}

    /**
     * @brief Copy constructor for the LabeledTree class.
     *
     * @param other The LabeledTree object to copy.
     */
    EdgeLabeledTree(const EdgeLabeledTree &other)
    {
        // recursively copy the tree
        root = copy_tree(other.root);
    }

    /**
     * @brief Destroy the Labeled Tree object.
     */
    ~EdgeLabeledTree()
    {
        // For large trees, we need to manually break the structure to avoid
        // stack overflow during recursive destruction of shared_ptrs
        if (root) {
            std::stack<std::shared_ptr<Node<LabelType>>> stack;
            stack.push(root);
            
            while (!stack.empty()) {
                auto current = stack.top();
                stack.pop();
                
                // Add children to stack before clearing them
                for (auto& child : current->get_children()) {
                    stack.push(child);
                }
                
                // Clear children to break references
                current->clear_children();
            }
            
            root.reset();
        }
    }

    /**
     * @brief Assignment operator for the LabeledTree class.
     *
     * @param other The LabeledTree object to assign from.
     * @return A reference to the assigned LabeledTree object.
     */
    EdgeLabeledTree &operator=(const EdgeLabeledTree &other)
    {
        if (this == &other)
        {
            return *this;
        }

        root = copy_tree(other.root);

        return *this;
    }

    /**
     * @brief Gets the root of the tree.
     *
     * @return A shared pointer to the root node.
     */
    std::shared_ptr<Node<LabelType>> get_root() const
    {
        return root;
    }

    /**
     * @brief Gets the depth of the tree.
     *
     * @return The depth of the tree.
     */
    unsigned int get_depth() const
    {
        return root ? root->get_depth() : 0;
    }

    /**
     * @brief Sets the root of the tree.
     *
     * @param newRoot A shared pointer to the new root node.
     */
    void set_root(std::shared_ptr<Node<LabelType>> newRoot)
    {
        root = newRoot;
    }

    /**
     * @brief Gets all nodes in the tree.
     *
     * @return A vector of shared pointers to all nodes in the tree.
     */
    std::vector<std::shared_ptr<Node<LabelType>>> get_nodes() const
    {
        std::vector<std::shared_ptr<Node<LabelType>>> nodes;
        collect_nodes(root, nodes);
        return nodes;
    }

    /**
     * @brief Converts the tree to a string representation.
     *
     * @return The string representation of the tree.
     */
    std::string to_string() const
    {
        std::ostringstream oss;
        to_string_helper(root, oss);
        return oss.str();
    }

    /**
     * @brief Gets all edge labels in the tree.
     *
     * @return A vector of all edge labels in the tree.
     */
    std::vector<LabelType> get_all_edge_labels() const
    {
        std::vector<LabelType> allLabels;
        if (!root) return allLabels;
        
        std::stack<std::shared_ptr<Node<LabelType>>> stack;
        stack.push(root);
        
        while (!stack.empty())
        {
            auto current = stack.top();
            stack.pop();
            
            for (const auto &edge : current->get_edges())
            {
                allLabels.push_back(edge.get_label());
                stack.push(edge.get_target());
            }
        }
        
        return allLabels;
    }

    /**
     * @brief Creates a simple tree with a root and one child connected by a labeled edge.
     *
     * @param edgeLabel The label for the edge from root to child.
     * @return A LabeledTree with root and one child.
     */
    static EdgeLabeledTree create_simple_tree(LabelType edgeLabel)
    {
        EdgeLabeledTree tree(true); // Create with root
        tree.get_root()->push_back_child(edgeLabel);
        return tree;
    }

    /**
     * @brief Stable sorts nodes based on their label path from node to root.
     * Nodes are sorted lexicographically by their path labels, with shorter paths coming first.
     * For nodes with the same path length, lexicographic comparison is used.
     * The sort is stable, preserving the relative order of nodes with equivalent paths.
     *
     * @return A vector of nodes sorted by their label path to root.
     */
    std::vector<std::shared_ptr<Node<LabelType>>> get_nodes_sorted_by_path() const
    {
        auto nodes = get_nodes();
        
        // Use stable_sort to maintain relative order for equivalent elements
        std::stable_sort(nodes.begin(), nodes.end(), 
            [](const std::shared_ptr<Node<LabelType>>& a, const std::shared_ptr<Node<LabelType>>& b) {
                auto pathA = a->get_path_to_root();
                auto pathB = b->get_path_to_root();
                
                // Compare paths lexicographically
                // Shorter paths come first, then lexicographic comparison
                if (pathA.size() != pathB.size()) {
                    return pathA.size() < pathB.size();
                }
                
                // Same length paths: lexicographic comparison
                return pathA < pathB;
            });
            
        return nodes;
    }

private:
    /**
     * @brief Recursively copies a tree.
     *
     * @param node The root of the tree to copy.
     * @return A shared pointer to the root of the new tree.
     */
    std::shared_ptr<Node<LabelType>> copy_tree(std::shared_ptr<Node<LabelType>> node)
    {
        if (!node)
        {
            return nullptr;
        }

        auto newNode = std::make_shared<Node<LabelType>>();
        for (const auto &edge : node->get_edges())
        {
            auto copiedChild = copy_tree(edge.get_target());
            newNode->push_back_child(copiedChild, edge.get_label());
        }

        return newNode;
    }

    /**
     * @brief Validates the string representation of a tree.
     *
     * @param str The string representation of the tree.
     * @return true if the string is a valid tree, false otherwise.
     */
    bool validate_tree(const std::string &str)
    {
        if (str.empty()) return false;

        enum StateType { START_NODE, CHECK_LABEL, CHECK_CLOSING };

        struct State {
            StateType type;
            State(StateType t) : type(t) {}
        };

        std::stack<State> stack;
        stack.push(State(START_NODE));

        size_t final_pos = 0;
        size_t cur_pos = 0;
        while (!stack.empty()) {
            State current = stack.top();
            stack.pop();

            if (cur_pos >= str.length()) return false;

            switch (current.type) {
                case START_NODE: {
                    if (str[cur_pos] != '(') return false;
                    stack.push(State(CHECK_LABEL));
                    ++cur_pos;
                    break;
                }
                case CHECK_LABEL: {
                    size_t label_start = cur_pos;
                    while (cur_pos < str.length() && std::isalnum(str[cur_pos])) {
                        cur_pos++;
                    }
                    if (label_start == cur_pos) return false;
                    stack.push(State(CHECK_CLOSING));
                    break;
                }
                case CHECK_CLOSING: {
                    if (cur_pos >= str.length()) return false;

                    if (str[cur_pos] == '(') {
                        stack.push(State(CHECK_CLOSING)); 
                        stack.push(State(START_NODE));
                    } else if (str[cur_pos] == ')') {
                        final_pos = ++cur_pos;
                    } else {
                        return false;
                    }

                    break;
                }
            }
        }

        return final_pos == str.length();
    }

    /**
     * @brief Creates a tree from its string representation.
     * In edge-labeled trees, the format is interpreted as: (root_marker(edge_label(...)...))
     * where the first label after '(' represents an edge label to a child.
     *
     * @param str The string representation of the tree.
     * @param strToLabel A function to convert a string to a label.
     * @return A shared pointer to the root of the created tree.
     */
    std::shared_ptr<Node<LabelType>> fromString(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        unsigned int pos = 0;
        if (str.empty() || !validate_tree(str))
        {
            throw std::invalid_argument("Invalid tree string. Error at position: " + std::to_string(pos));
        }

        std::stack<std::shared_ptr<Node<LabelType>>> nodeStack;
        std::istringstream iss(str);
        char ch;
        std::shared_ptr<Node<LabelType>> currentNode, root = nullptr;
        bool expectingEdgeLabel = false;
        
        while (iss >> ch)
        {
            if (ch == '(')
            {
                expectingEdgeLabel = true;
            }
            else if (ch == ')')
            {
                if (!nodeStack.empty())
                {
                    nodeStack.pop();
                }
                expectingEdgeLabel = false;
            }
            else
            {
                std::string labelStr = std::string(1, ch);
                while (iss.peek() != '(' && iss.peek() != ')' && iss >> ch)
                {
                    labelStr += ch;
                }

                if (nodeStack.empty())
                {
                    // First node encountered - this is the root
                    // In edge-labeled trees, root has no incoming edge, so we ignore this label
                    root = std::make_shared<Node<LabelType>>();
                    nodeStack.push(root);
                }
                else if (expectingEdgeLabel)
                {
                    // This label represents an edge to a new child
                    LabelType edgeLabel = strToLabel(labelStr);
                    currentNode = nodeStack.top();
                    currentNode->push_back_child(edgeLabel);
                    // Get the newly created child and push it to stack
                    auto children = currentNode->get_children();
                    nodeStack.push(children.back());
                }
                expectingEdgeLabel = false;
            }
        }

        return root;
    }

    /**
     * @brief Collects all nodes in a subtree.
     *
     * @param node The root of the subtree.
     * @param nodes A reference to a vector to store the nodes.
     */
    void collect_nodes(std::shared_ptr<Node<LabelType>> node, std::vector<std::shared_ptr<Node<LabelType>>> &nodes) const
    {
        if (!node) return;

        std::stack<std::shared_ptr<Node<LabelType>>> stack;
        stack.push(node);

        while (!stack.empty()) {
            auto current = stack.top();
            stack.pop();
            nodes.push_back(current);

            const auto &children = current->get_children();
            for (auto it = children.rbegin(); it != children.rend(); ++it) {
                stack.push(*it);
            }
        }
    }

    /**
     * @brief Helper function to convert a subtree to a string.
     * In edge-labeled trees, the format is: (root_marker(edge_label(...)...))
     *
     * @param node The root of the subtree.
     * @param oss The output string stream.
     */
    void to_string_helper(const std::shared_ptr<Node<LabelType>> node, std::ostringstream &oss) const
    {
        if (!node)
            return;

        oss << '(';
        // For root node, we use a placeholder since it has no incoming edge
        if (node->get_parent() == nullptr)
        {
            oss << "root";
        }
        
        if (!node->is_leaf())
        {
            const auto &edges = node->get_edges();
            for (const auto &edge : edges)
            {
                oss << '(' << edge.get_label();
                to_string_helper(edge.get_target(), oss);
                oss << ')';
            }
        }
        oss << ')';
    }
};

#endif // EDGELABELEDTREE_HPP
