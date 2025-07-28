#ifndef LABELEDTREE_HPP
#define LABELEDTREE_HPP

#include <vector>
#include <memory>
#include <string>
#include <stack>
#include <sstream>
#include <iostream>
#include <functional>

template <typename LabelType>
class Node
{
private:
    LabelType m_label;
    std::vector<Node *> m_children;
    Node *m_parent;

public:
    Node(LabelType lbl, Node *prnt = nullptr, std::vector<Node *> children = {}) : m_label(lbl), m_children(children), m_parent(prnt) {}

    ~Node()
    {
        for (auto child : m_children)
            delete child;
    }

    bool is_root() const
    {
        return m_parent == nullptr;
    }

    bool is_leaf() const
    {
        return m_children.empty();
    }

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

    unsigned int get_depth() const
    {
        unsigned int maxDepth = 0;
        for (const auto &child : m_children)
        {
            maxDepth = std::max(maxDepth, child->get_depth());
        }
        return maxDepth + 1;
    }

    LabelType get_label() const
    {
        return m_label;
    }

    const std::vector<Node *> &get_children() const
    {
        return m_children;
    }

    Node * get_parent() const
    {
        return m_parent;
    }

    bool is_rightmost() const
    {
        if (m_parent)
        {
            const auto &siblings = m_parent->m_children;
            return siblings.back() == this;
        }
        return false;
    }

    void push_back_child(LabelType lbl)
    {
        auto child = new Node(lbl, this);
        m_children.push_back(child);
    }

    void push_back_child(Node *child)
    {
        child->m_parent = this;
        m_children.push_back(child);
    }

    void prepend_child(Node *child)
    {
        child->m_parent = this;
        m_children.insert(m_children.begin(), child);
    }
};

template <typename LabelType>
class LabeledTree
{
private:
    Node<LabelType> *root;

public:
    LabeledTree() : root(nullptr) {}

    LabeledTree(LabelType rootLabel)
    {
        root = new Node<LabelType>(rootLabel);
    }

    LabeledTree(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        root = fromString(str, strToLabel);
    }

    LabeledTree(Node<LabelType> *root) : root(root) {}

    LabeledTree(const LabeledTree &other)
    {
        // recursively copy the tree
        root = copy_tree(other.root);
    }

    ~LabeledTree()
    {
        delete root;
    }

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

    Node<LabelType> *get_root() const
    {
        return root;
    }

    unsigned int get_depth() const
    {
        return root ? root->get_depth() : 0;
    }

    void set_root(Node<LabelType> *newRoot)
    {
        root = newRoot;
    }

    std::vector<Node<LabelType> *> get_nodes() const
    {
        std::vector<Node<LabelType> *> nodes;
        collect_nodes(root, nodes);
        return nodes;
    }

    std::string to_string() const
    {
        std::ostringstream oss;
        to_string_helper(root, oss);
        return oss.str();
    }

private:
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

    bool are_parentheses_balanced(const std::string &str)
    {
        int balance = 0;
        for (char ch : str)
        {
            if (ch == '(')
            {
                ++balance;
            }
            else if (ch == ')')
            {
                --balance;
                if (balance < 0)
                {
                    return false; // More closing than opening
                }
            }
        }
        return balance == 0; // Should be zero if balanced
    }

    // Recursive function to validate the tree structure
    bool is_valid_tree(const std::string &str, unsigned int &pos)
    {
        if (pos >= str.length())
        {
            return false;
        }

        // Expecting an opening parenthesis
        if (str[pos] != '(')
        {
            return false;
        }
        ++pos;

        // Expecting a node label consisting of alphanumeric characters
        size_t label_start = pos;
        while (pos < str.length() && std::isalnum(str[pos]))
        {
            ++pos;
        }
        if (label_start == pos)
        {
            return false; // No valid label found
        }

        // Recursively check for child nodes
        while (pos < str.length() && str[pos] == '(')
        {
            if (!is_valid_tree(str, pos))
            {
                return false;
            }
        }

        // Expecting a closing parenthesis
        if (pos >= str.length() || str[pos] != ')')
        {
            return false;
        }
        ++pos;

        return true;
    }

    Node<LabelType> *fromString(const std::string &str, std::function<LabelType(const std::string &)> strToLabel)
    {
        unsigned int pos = 0;
        if (str.empty() || (!are_parentheses_balanced(str) || !is_valid_tree(str, pos)))
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

    void collect_nodes(Node<LabelType> *node, std::vector<Node<LabelType> *> &nodes) const
    {
        nodes.push_back(node);
        for (const auto &child : node->get_children())
        {
            collect_nodes(child, nodes);
        }
    }

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