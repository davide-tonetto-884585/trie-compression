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
    LabelType label;
    std::vector<Node *> children;
    Node *parent;

public:
    Node(LabelType lbl, Node *prnt = nullptr, std::vector<Node *> children = {}) : label(lbl), parent(prnt), children(children) {}

    ~Node()
    {
        for (auto child : children)
            delete child;
    }

    bool isRoot() const
    {
        return parent == nullptr;
    }

    bool isLeaf() const
    {
        return children.empty();
    }

    unsigned int getLevel() const
    {
        unsigned int level = 0;
        Node *current = parent;
        while (current)
        {
            level++;
            current = current->parent;
        }
        return level;
    }

    unsigned int getDepth() const
    {
        unsigned int maxDepth = 0;
        for (const auto &child : children)
        {
            maxDepth = std::max(maxDepth, child->getDepth());
        }
        return maxDepth + 1;
    }

    LabelType getLabel() const
    {
        return label;
    }

    const std::vector<Node *> &getChildren() const
    {
        return children;
    }

    Node * getParent() const
    {
        return parent;
    }

    bool isRightmost() const
    {
        if (parent)
        {
            const auto &siblings = parent->children;
            return siblings.back() == this;
        }
        return false;
    }

    void pushBackChild(LabelType lbl)
    {
        auto child = new Node(lbl, this);
        children.push_back(child);
    }

    void pushBackChild(Node *child)
    {
        child->parent = this;
        children.push_back(child);
    }

    void prependChild(Node *child)
    {
        child->parent = this;
        children.insert(children.begin(), child);
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
        root = copyTree(other.root);
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
        root = copyTree(other.root);

        return *this;
    }

    Node<LabelType> *getRoot() const
    {
        return root;
    }

    unsigned int getDepth() const
    {
        return root ? root->getDepth() : 0;
    }

    void setRoot(Node<LabelType> *newRoot)
    {
        root = newRoot;
    }

    std::vector<Node<LabelType> *> getNodes() const
    {
        std::vector<Node<LabelType> *> nodes;
        collectNodes(root, nodes);
        return nodes;
    }

    std::string toString() const
    {
        std::ostringstream oss;
        toStringHelper(root, oss);
        return oss.str();
    }

private:
    Node<LabelType> *copyTree(Node<LabelType> *node)
    {
        if (!node)
        {
            return nullptr;
        }

        auto newNode = new Node<LabelType>(node->getLabel());
        for (const auto &child : node->getChildren())
        {
            newNode->pushBackChild(copyTree(child));
        }

        return newNode;
    }

    bool areParenthesesBalanced(const std::string &str)
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
    bool isValidTree(const std::string &str, unsigned int &pos)
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
            if (!isValidTree(str, pos))
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
        if (str.empty() || (!areParenthesesBalanced(str) || !isValidTree(str, pos)))
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
                    currentNode->pushBackChild(label);
                    nodeStack.push(currentNode->getChildren().back());
                }
            }
        }

        return root;
    }

    void collectNodes(Node<LabelType> *node, std::vector<Node<LabelType> *> &nodes) const
    {
        nodes.push_back(node);
        for (const auto &child : node->getChildren())
        {
            collectNodes(child, nodes);
        }
    }

    void toStringHelper(const Node<LabelType> *node, std::ostringstream &oss) const
    {
        if (!node)
            return;

        oss << '(';
        oss << node->getLabel();
        if (!node->isLeaf())
        {
            const auto &children = node->getChildren();
            for (size_t i = 0; i < children.size(); ++i)
            {
                toStringHelper(children[i], oss);
            }
        }
        oss << ')';
    }
};

#endif // LABELEDTREE_HPP