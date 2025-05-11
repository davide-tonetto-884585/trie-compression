#ifndef LABELEDTREE_HPP
#define LABELEDTREE_HPP

#include <vector>
#include <optional>
#include <memory> // Required for std::unique_ptr, std::make_unique
#include <string>
#include <stack>
#include <sstream>
#include <iostream>
#include <functional>
#include <algorithm> // Required for std::max
#include <stdexcept> // Required for std::invalid_argument

// Forward declaration of LabeledTree for Node's friend declaration if needed
template <typename LabelType>
class LabeledTree;

template <typename LabelType>
class Node {
private:
    LabelType label;
    std::vector<std::unique_ptr<Node<LabelType>>> children;
    Node<LabelType>* parent; // Non-owning raw pointer

public:
    Node(LabelType lbl, Node<LabelType>* prnt = nullptr) : label(std::move(lbl)), parent(prnt) {}

    bool isRoot() const {
        return parent == nullptr;
    }

    bool isLeaf() const {
        return children.empty();
    }

    unsigned int getLevel() const {
        unsigned int level = 0;
        Node<LabelType>* current = parent;
        while (current) {
            level++;
            current = current->parent;
        }
        return level;
    }

    // Calculates the height of the subtree rooted at this node
    unsigned int getHeight() const {
        unsigned int maxHeight = 0;
        if (isLeaf()) {
            return 0; // Height of a leaf node is 0
        }
        for (const auto& child : children) {
            maxHeight = std::max(maxHeight, child->getHeight());
        }
        return maxHeight + 1;
    }

    const LabelType& getLabel() const {
        return label;
    }

    LabelType& getLabel() { // Non-const version if needed
        return label;
    }

    const std::vector<std::unique_ptr<Node<LabelType>>>& getChildren() const {
        return children;
    }

    // Allows modification of children, e.g., by LabeledTree
    std::vector<std::unique_ptr<Node<LabelType>>>& getChildren_non_const() {
        return children;
    }


    Node<LabelType>* getParent() const {
        return parent;
    }
    
    std::optional<Node<LabelType>*> getParentOptional() const {
        return parent ? std::optional<Node<LabelType>*>(parent) : std::nullopt;
    }


    bool isRightmost() const {
        if (parent) {
            const auto& siblings = parent->children;
            if (siblings.empty()) return false; // Should not happen if parent points to this
            return siblings.back().get() == this;
        }
        return false; // A root node isn't rightmost in the context of siblings
    }

    Node<LabelType>* pushBackChild(LabelType lbl) {
        children.push_back(std::make_unique<Node<LabelType>>(std::move(lbl), this));
        return children.back().get();
    }

    Node<LabelType>* pushBackChild(std::unique_ptr<Node<LabelType>> child) {
        if (child) {
            child->parent = this;
            children.push_back(std::move(child));
            return children.back().get();
        }
        return nullptr;
    }

    Node<LabelType>* prependChild(LabelType lbl) {
        children.insert(children.begin(), std::make_unique<Node<LabelType>>(std::move(lbl), this));
        return children.front().get();
    }
    
    Node<LabelType>* prependChild(std::unique_ptr<Node<LabelType>> child) {
        if (child) {
            child->parent = this;
            children.insert(children.begin(), std::move(child));
            return children.front().get();
        }
        return nullptr;
    }

    // Friend class to allow LabeledTree to manage nodes (e.g. in copy operations or setRoot)
    friend class LabeledTree<LabelType>;
};

template <typename LabelType>
class LabeledTree {
private:
    std::unique_ptr<Node<LabelType>> root;

    // Helper for recursive deep copy
    std::unique_ptr<Node<LabelType>> copyTreeRecursive(const Node<LabelType>* node, Node<LabelType>* parent) {
        if (!node) {
            return nullptr;
        }
        auto newNode = std::make_unique<Node<LabelType>>(node->getLabel(), parent);
        for (const auto& child : node->getChildren()) {
            newNode->pushBackChild(copyTreeRecursive(child.get(), newNode.get()));
        }
        return newNode;
    }

    // Helper for fromString
    static std::unique_ptr<Node<LabelType>> fromStringRecursive(std::istringstream& iss, Node<LabelType>* parent_node, const std::function<LabelType(const std::string&)>& strToLabel) {
        char ch;
        std::string labelStr;

        // Expect '('
        if (!(iss >> ch) || ch != '(') {
            throw std::invalid_argument("Invalid tree string: Expected '('.");
        }

        // Read label until '(' or ')'
        while (iss.peek() != EOF && iss.peek() != '(' && iss.peek() != ')') {
            iss >> ch;
            labelStr += ch;
        }
        if (labelStr.empty()) {
            throw std::invalid_argument("Invalid tree string: Label cannot be empty.");
        }
        
        LabelType label = strToLabel(labelStr);
        auto currentNode = std::make_unique<Node<LabelType>>(std::move(label), parent_node);

        // Read children
        while (iss.peek() == '(') {
            currentNode->pushBackChild(fromStringRecursive(iss, currentNode.get(), strToLabel));
        }

        // Expect ')'
        if (!(iss >> ch) || ch != ')') {
            throw std::invalid_argument("Invalid tree string: Expected ')'.");
        }
        return currentNode;
    }


    void collectNodesRecursive(Node<LabelType>* node, std::vector<Node<LabelType>*>& nodes_vec) const {
        if (!node) return;
        nodes_vec.push_back(node);
        for (const auto& child : node->getChildren()) {
            collectNodesRecursive(child.get(), nodes_vec);
        }
    }

    void toStringHelper(const Node<LabelType>* node, std::ostringstream& oss) const {
        if (!node) return;

        oss << '(';
        // Assuming LabelType is streamable. If not, a specific to_string(LabelType) might be needed.
        oss << node->getLabel(); 
        for (const auto& child : node->getChildren()) {
            toStringHelper(child.get(), oss);
        }
        oss << ')';
    }


public:
    LabeledTree() : root(nullptr) {}

    explicit LabeledTree(LabelType rootLabel) {
        root = std::make_unique<Node<LabelType>>(std::move(rootLabel), nullptr);
    }
    
    // Constructor taking ownership of a raw pointer (use with caution)
    explicit LabeledTree(Node<LabelType>* rawRoot) : root(rawRoot) {
        if (root) {
            root->parent = nullptr; // Ensure root's parent is correctly null
        }
    }

    // Constructor taking ownership of a unique_ptr
    explicit LabeledTree(std::unique_ptr<Node<LabelType>> rootNode) : root(std::move(rootNode)) {
        if (root) {
            root->parent = nullptr;
        }
    }

    LabeledTree(const std::string& str, std::function<LabelType(const std::string&)> strToLabel) : root(nullptr) {
        if (str.empty()) {
            return; // Empty tree
        }
        std::istringstream iss(str);
        root = fromStringRecursive(iss, nullptr, strToLabel);
        
        // Check if the entire string was consumed
        char remainingChar;
        if (iss >> remainingChar) {
            throw std::invalid_argument("Invalid tree string: Extra characters after main tree structure.");
        }
    }

    // Copy constructor
    LabeledTree(const LabeledTree& other) {
        if (other.root) {
            root = copyTreeRecursive(other.root.get(), nullptr);
        } else {
            root = nullptr;
        }
    }

    // Move constructor
    LabeledTree(LabeledTree&& other) noexcept : root(std::move(other.root)) {}

    // Destructor (default is fine due to unique_ptr)
    ~LabeledTree() = default;

    // Copy assignment operator
    LabeledTree& operator=(const LabeledTree& other) {
        if (this == &other) {
            return *this;
        }
        if (other.root) {
            root = copyTreeRecursive(other.root.get(), nullptr);
        } else {
            root.reset(); // Clear current tree
        }
        return *this;
    }

    // Move assignment operator
    LabeledTree& operator=(LabeledTree&& other) noexcept {
        if (this == &other) {
            return *this;
        }
        root = std::move(other.root);
        return *this;
    }

    Node<LabelType>* getRoot() const {
        return root.get();
    }

    // Renamed Node::getDepth() to Node::getHeight() for clarity
    // LabeledTree::getTreeHeight() returns the height of the entire tree
    unsigned int getTreeHeight() const {
        return root ? root->getHeight() : 0; 
        // If root exists, height is root's height. An empty tree or tree with only root has height 0.
        // If definition is number of edges on longest path from root to leaf:
        // return root ? root->getHeight() : -1; // or throw for empty tree
        // If definition is number of nodes on longest path from root to leaf:
        // return root ? root->getHeight() +1 : 0;
        // Current Node::getHeight returns max edges from node to leaf, so this is consistent.
    }


    void setRoot(std::unique_ptr<Node<LabelType>> newRoot) {
        root = std::move(newRoot);
        if (root) {
            root->parent = nullptr;
        }
    }
    
    void setRootLabel(LabelType newLabel) {
        if (root) {
            root->getLabel() = std::move(newLabel);
        } else {
            root = std::make_unique<Node<LabelType>>(std::move(newLabel), nullptr);
        }
    }


    std::vector<Node<LabelType>*> getNodes() const {
        std::vector<Node<LabelType>*> nodes_vec;
        if (root) {
            collectNodesRecursive(root.get(), nodes_vec);
        }
        return nodes_vec;
    }

    std::string toString() const {
        if (!root) return "()"; // Or "" for an empty tree, depending on convention
        std::ostringstream oss;
        toStringHelper(root.get(), oss);
        return oss.str();
    }
};

#endif // LABELEDTREE_HPP
