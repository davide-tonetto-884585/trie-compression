#ifndef LABELEDTREE_HPP
#define LABELEDTREE_HPP

#include <vector>
#include <memory>
#include <string>
#include <stack>
#include <sstream>
#include <iostream>
#include <functional>

/**
 * @brief A template class representing a node in a labeled tree.
 * 
 * @tparam LabelType The type of the label for each node.
 */
template <typename LabelType>
class Node
{
private:
    LabelType m_label;              ///< The label of the node.
    std::vector<Node *> m_children; ///< A vector of pointers to the children of the node.
    Node *m_parent;                 ///< A pointer to the parent of the node.

public:
    /**
     * @brief Construct a new Node object.
     *
     * @param lbl The label of the node.
     * @param prnt A pointer to the parent node. Defaults to nullptr.
     * @param children A vector of pointers to the children nodes. Defaults to an empty vector.
     */
    Node(LabelType lbl, Node *prnt = nullptr, std::vector<Node *> children = {}) : m_label(lbl), m_children(children), m_parent(prnt) {}

    /**
     * @brief Checks if the node is a root.
     *
     * @return true if the node is a root, false otherwise.
     */
    bool is_root() const
    {
        return m_parent == nullptr;
    }

    /**
     * @brief Checks if the node is a leaf.
     *
     * @return true if the node is a leaf, false otherwise.
     */
    bool is_leaf() const
    {
        return m_children.empty();
    }

    /**
     * @brief Gets the level of the node in the tree.
     *
     * @return The level of the node.
     */
    unsigned int get_level() const
    {
        unsigned int level = 0;
        Node *current = m_parent;
        while (current)
        {
            level++;
            current = current->m_parent;
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
        for (const auto &child : m_children)
        {
            maxDepth = std::max(maxDepth, child->get_depth());
        }
        return maxDepth + 1;
    }

    /**
     * @brief Gets the label of the node.
     *
     * @return The label of the node.
     */
    LabelType get_label() const
    {
        return m_label;
    }

    /**
     * @brief Gets the children of the node.
     *
     * @return A const reference to the vector of children nodes.
     */
    const std::vector<Node *> &get_children() const
    {
        return m_children;
    }

    /**
     * @brief Gets the parent of the node.
     *
     * @return A pointer to the parent node.
     */
    Node *get_parent() const
    {
        return m_parent;
    }

    /**
     * @brief Checks if the node is the rightmost child of its parent.
     *
     * @return true if the node is the rightmost child, false otherwise.
     */
    bool is_rightmost() const
    {
        if (m_parent)
        {
            const auto &siblings = m_parent->m_children;
            return siblings.back() == this;
        }
        return false;
    }

    /**
     * @brief Adds a new child to the node with a given label.
     *
     * @param lbl The label of the new child node.
     */
    void push_back_child(LabelType lbl)
    {
        auto child = new Node(lbl, this);
        m_children.push_back(child);
    }

    /**
     * @brief Adds an existing node as a child to this node.
     *
     * @param child A pointer to the node to be added as a child.
     */
    void push_back_child(Node *child)
    {
        child->m_parent = this;
        m_children.push_back(child);
    }

    /**
     * @brief Prepends an existing node as a child to this node.
     *
     * @param child A pointer to the node to be prepended as a child.
     */
    void prepend_child(Node *child)
    {
        child->m_parent = this;
        m_children.insert(m_children.begin(), child);
    }
};

/**
 * @brief A template class representing a labeled tree.
 * 
 * @tparam LabelType The type of the label for each node.
 */
template <typename LabelType>
class LabeledTree
{
private:
    Node<LabelType> *root; ///< A pointer to the root node of the tree.

public:
    /**
     * @brief Construct a new Labeled Tree object.
     */
    LabeledTree() : root(nullptr) {}

    /**
     * @brief Construct a new Labeled Tree object with a root node.
     *
     * @param rootLabel The label of the root node.
     */
    LabeledTree(LabelType rootLabel)
    {
        root = new Node<LabelType>(rootLabel);
    }

    /**
     * @brief Construct a new Labeled Tree object from a string representation.
     *
     * @param str The string representation of the tree.
     * @param strToLabel A function to convert a string to a label.
     */
    LabeledTree(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        root = fromString(str, strToLabel);
    }

    /**
     * @brief Construct a new Labeled Tree object from a root node.
     *
     * @param root A pointer to the root node.
     */
    LabeledTree(Node<LabelType> *root) : root(root) {}

    /**
     * @brief Copy constructor for the LabeledTree class.
     *
     * @param other The LabeledTree object to copy.
     */
    LabeledTree(const LabeledTree &other)
    {
        // recursively copy the tree
        root = copy_tree(other.root);
    }

    /**
     * @brief Destroy the Labeled Tree object.
     */
    ~LabeledTree()
    {
        auto nodes = this->get_nodes();
        for (auto node : nodes)
        {
            delete node;
        }
    }

    /**
     * @brief Assignment operator for the LabeledTree class.
     *
     * @param other The LabeledTree object to assign from.
     * @return A reference to the assigned LabeledTree object.
     */
    LabeledTree &operator=(const LabeledTree &other)
    {
        if (this == &other)
        {
            return *this;
        }

        delete root;
        root = copy_tree(other.root);

        return *this;
    }

    /**
     * @brief Gets the root of the tree.
     *
     * @return A pointer to the root node.
     */
    Node<LabelType> *get_root() const
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
     * @param newRoot A pointer to the new root node.
     */
    void set_root(Node<LabelType> *newRoot)
    {
        root = newRoot;
    }

    /**
     * @brief Gets all nodes in the tree.
     *
     * @return A vector of pointers to all nodes in the tree.
     */
    std::vector<Node<LabelType> *> get_nodes() const
    {
        std::vector<Node<LabelType> *> nodes;
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

private:
    /**
     * @brief Recursively copies a tree.
     *
     * @param node The root of the tree to copy.
     * @return A pointer to the root of the new tree.
     */
    Node<LabelType> *copy_tree(Node<LabelType> *node)
    {
        if (!node)
        {
            return nullptr;
        }

        auto newNode = new Node<LabelType>(node->get_label());
        for (const auto &child : node->get_children())
        {
            newNode->push_back_child(copy_tree(child));
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
     *
     * @param str The string representation of the tree.
     * @param strToLabel A function to convert a string to a label.
     * @return A pointer to the root of the created tree.
     */
    Node<LabelType> *fromString(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        unsigned int pos = 0;
        if (str.empty() || !validate_tree(str))
        {
            throw std::invalid_argument("Invalid tree string. Error at position: " + std::to_string(pos));
        }

        std::stack<Node<LabelType> *> nodeStack;
        std::istringstream iss(str);
        char ch;
        Node<LabelType> *currentNode, *root = nullptr;
        while (iss >> ch)
        {
            if (ch == '(')
            {
                // Do nothing, just continue
            }
            else if (ch == ')')
            {
                nodeStack.pop();
            }
            else
            {
                std::string labelStr = std::string(1, ch);
                while (iss.peek() != '(' && iss.peek() != ')' && iss >> ch)
                {
                    labelStr += ch;
                }

                LabelType label = strToLabel(labelStr);
                if (nodeStack.empty())
                {
                    currentNode = root = new Node<LabelType>(label);
                    nodeStack.push(currentNode);
                }
                else
                {
                    currentNode = nodeStack.top();
                    currentNode->push_back_child(label);
                    nodeStack.push(currentNode->get_children().back());
                }
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
    void collect_nodes(Node<LabelType> *node, std::vector<Node<LabelType> *> &nodes) const
    {
        if (!node) return;

        std::stack<Node<LabelType> *> stack;
        stack.push(node);

        while (!stack.empty()) {
            Node<LabelType> *current = stack.top();
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
     *
     * @param node The root of the subtree.
     * @param oss The output string stream.
     */
    void to_string_helper(const Node<LabelType> *node, std::ostringstream &oss) const
    {
        if (!node)
            return;

        oss << '(';
        oss << node->get_label();
        if (!node->is_leaf())
        {
            const auto &children = node->get_children();
            for (size_t i = 0; i < children.size(); ++i)
            {
                to_string_helper(children[i], oss);
            }
        }
        oss << ')';
    }
};

#endif // LABELEDTREE_HPP