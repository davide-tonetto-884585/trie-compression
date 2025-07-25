#!/usr/bin/env python3
"""
Repetitive Tree Generator

Generates trees with high levels of repetitiveness by recursively copying subtrees.
Allows control over the number of different subtree patterns and their repetition frequency.
"""

import random
import json
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from copy import deepcopy


@dataclass
class TreeNode:
    """Represents a node in the tree."""
    label: str
    children: List['TreeNode']
    node_id: Optional[int] = None
    
    def __post_init__(self):
        if self.children is None:
            self.children = []
    
    def add_child(self, child: 'TreeNode'):
        """Add a child node."""
        self.children.append(child)
    
    def is_leaf(self) -> bool:
        """Check if this node is a leaf."""
        return len(self.children) == 0
    
    def size(self) -> int:
        """Calculate the total number of nodes in this subtree."""
        return 1 + sum(child.size() for child in self.children)
    
    def depth(self) -> int:
        """Calculate the depth of this subtree."""
        if self.is_leaf():
            return 1
        return 1 + max(child.depth() for child in self.children)
    
    def to_dict(self) -> Dict:
        """Convert tree to dictionary representation."""
        return {
            'label': self.label,
            'node_id': self.node_id,
            'children': [child.to_dict() for child in self.children]
        }
    
    def __str__(self) -> str:
        """String representation of the tree."""
        return self._str_helper(0)
    
    def _str_helper(self, indent: int) -> str:
        """Helper method for string representation."""
        result = '  ' * indent + f"{self.label}\n"
        for child in self.children:
            result += child._str_helper(indent + 1)
        return result


class RepetitiveTreeGenerator:
    """Generates trees with repetitive subtree patterns."""
    
    def __init__(self, 
                 max_depth: int = 6,
                 max_branching_factor: int = 3,
                 num_subtree_templates: int = 3,
                 repetition_probability: float = 0.7,
                 template_depth_range: Tuple[int, int] = (2, 4),
                 alphabet_type: str = "letters",
                 alphabet_size: int = 26,
                 seed: Optional[int] = None):
        """
        Initialize the repetitive tree generator.
        
        Args:
            max_depth: Maximum depth of the generated tree
            max_branching_factor: Maximum number of children per node
            num_subtree_templates: Number of different subtree templates to create
            repetition_probability: Probability of reusing an existing template vs creating new
            template_depth_range: (min, max) depth range for subtree templates
            alphabet_type: "letters" for a-z or "numbers" for 1-n
            alphabet_size: Size of alphabet (max 26 for letters, any positive int for numbers)
            seed: Random seed for reproducible generation
        """
        self.max_depth = max_depth
        self.max_branching_factor = max_branching_factor
        self.num_subtree_templates = num_subtree_templates
        self.repetition_probability = repetition_probability
        self.template_depth_range = template_depth_range
        self.alphabet_type = alphabet_type
        self.alphabet_size = min(alphabet_size, 26) if alphabet_type == "letters" else alphabet_size
        
        if seed is not None:
            random.seed(seed)
        
        # Storage for subtree templates
        self.subtree_templates: List[TreeNode] = []
        self.template_usage_count: List[int] = []
        self.next_node_id = 0
        self.label_counter = 0
        
        # Generate alphabet
        self._generate_alphabet()
        
        # Generate initial subtree templates
        self._generate_subtree_templates()
    
    def _generate_alphabet(self):
        """Generate the alphabet based on the specified type and size."""
        if self.alphabet_type == "letters":
            self.alphabet = [chr(ord('a') + i) for i in range(self.alphabet_size)]
        elif self.alphabet_type == "numbers":
            self.alphabet = [str(i + 1) for i in range(self.alphabet_size)]
        else:
            raise ValueError(f"Invalid alphabet_type: {self.alphabet_type}. Use 'letters' or 'numbers'.")
    
    def _get_next_node_id(self) -> int:
        """Get the next available node ID."""
        node_id = self.next_node_id
        self.next_node_id += 1
        return node_id
    
    def _get_next_label(self) -> str:
        """Generate a label from the alphabet."""
        if self.label_counter < len(self.alphabet):
            label = self.alphabet[self.label_counter]
        else:
            # If we run out of alphabet, use combinations
            base = len(self.alphabet)
            index = self.label_counter
            label = ""
            while True:
                label = self.alphabet[index % base] + label
                index //= base
                if index == 0:
                    break
        
        self.label_counter += 1
        return label
    
    def _generate_label(self, prefix: str = "N") -> str:
        """Generate a label for a node (legacy method for compatibility)."""
        return f"{prefix}{self._get_next_node_id()}"
    
    def _generate_subtree_templates(self):
        """Generate the initial set of subtree templates."""
        for i in range(self.num_subtree_templates):
            template_depth = random.randint(*self.template_depth_range)
            template = self._create_random_subtree(template_depth)
            self.subtree_templates.append(template)
            self.template_usage_count.append(0)
    
    def _create_random_subtree(self, max_depth: int) -> TreeNode:
        """Create a random subtree with the given maximum depth using alphabet labels."""
        root = TreeNode(label=self._get_next_label(), children=[])
        root.node_id = self._get_next_node_id()
        
        if max_depth > 1:
            num_children = random.randint(1, min(self.max_branching_factor, 3))
            for _ in range(num_children):
                child = self._create_random_subtree(max_depth - 1)
                root.add_child(child)
        
        return root
    
    def _copy_subtree_template(self, template: TreeNode) -> TreeNode:
        """Create a copy of a subtree template preserving original labels but with new node IDs."""
        def copy_node(node: TreeNode) -> TreeNode:
            new_node = TreeNode(
                label=node.label,  # Preserve original label
                children=[],
                node_id=self._get_next_node_id()  # New unique node ID
            )
            for child in node.children:
                new_node.add_child(copy_node(child))
            return new_node
        
        return copy_node(template)
    
    def _should_use_template(self, current_depth: int) -> bool:
        """Decide whether to use a template or create a new subtree."""
        # Higher probability of using templates at deeper levels
        depth_factor = min(current_depth / self.max_depth, 1.0)
        adjusted_probability = self.repetition_probability * (0.5 + 0.5 * depth_factor)
        return random.random() < adjusted_probability
    
    def _select_template(self) -> int:
        """Select a template to use, with preference for less-used templates."""
        if not self.template_usage_count:
            return 0
        
        # Calculate weights inversely proportional to usage count
        max_usage = max(self.template_usage_count) + 1
        weights = [max_usage - count for count in self.template_usage_count]
        
        # Weighted random selection
        total_weight = sum(weights)
        if total_weight == 0:
            return random.randint(0, len(self.subtree_templates) - 1)
        
        rand_val = random.uniform(0, total_weight)
        cumulative = 0
        for i, weight in enumerate(weights):
            cumulative += weight
            if rand_val <= cumulative:
                return i
        
        return len(self.subtree_templates) - 1
    
    def generate_tree(self) -> TreeNode:
        """Generate a repetitive tree."""
        # Reset usage counters and label counter for new tree
        self.template_usage_count = [0] * len(self.subtree_templates)
        self.label_counter = 0
        
        root = TreeNode(label=self._get_next_label(), children=[])
        root.node_id = self._get_next_node_id()
        
        self._build_tree_recursive(root, 1)
        return root
    
    def _build_tree_recursive(self, parent: TreeNode, current_depth: int):
        """Recursively build the tree with repetitive patterns."""
        if current_depth >= self.max_depth:
            return
        
        num_children = random.randint(1, self.max_branching_factor)
        
        for _ in range(num_children):
            if (current_depth > 1 and 
                self.subtree_templates and 
                self._should_use_template(current_depth)):
                
                # Use a template
                template_idx = self._select_template()
                template = self.subtree_templates[template_idx]
                child = self._copy_subtree_template(template)
                self.template_usage_count[template_idx] += 1
                parent.add_child(child)
            else:
                # Create a new subtree
                child = TreeNode(label=self._get_next_label(), children=[])
                child.node_id = self._get_next_node_id()
                parent.add_child(child)
                self._build_tree_recursive(child, current_depth + 1)
    
    def get_statistics(self, tree: TreeNode) -> Dict:
        """Get statistics about the generated tree."""
        return {
            'total_nodes': tree.size(),
            'tree_depth': tree.depth(),
            'num_templates': len(self.subtree_templates),
            'template_usage': dict(enumerate(self.template_usage_count)),
            'total_template_usage': sum(self.template_usage_count),
            'repetition_ratio': sum(self.template_usage_count) / tree.size() if tree.size() > 0 else 0
        }
    
    def save_tree_to_json(self, tree: TreeNode, filename: str):
        """Save the tree to a JSON file."""
        with open(filename, 'w') as f:
            json.dump(tree.to_dict(), f, indent=2)
    
    def print_templates(self):
        """Print all subtree templates."""
        print("Subtree Templates:")
        print("=" * 50)
        for i, template in enumerate(self.subtree_templates):
            print(f"Template {i} (used {self.template_usage_count[i]} times):")
            print(template)
            print("-" * 30)


def main():
    """Example usage of the RepetitiveTreeGenerator."""
    print("Repetitive Tree Generator Demo")
    print("=" * 50)
    
    # Demo 1: Tree with letter alphabet
    print("\n=== Demo 1: Letter Alphabet (a-z) ===")
    letter_generator = RepetitiveTreeGenerator(
        max_depth=4,
        max_branching_factor=2,
        num_subtree_templates=2,
        repetition_probability=0.8,
        template_depth_range=(2, 3),
        alphabet_type="letters",
        alphabet_size=5,  # Use only a-e
        seed=42
    )
    
    letter_tree = letter_generator.generate_tree()
    print("Generated Tree with Letter Labels:")
    print(letter_tree)
    
    print("\nTemplates:")
    letter_generator.print_templates()
    
    stats = letter_generator.get_statistics(letter_tree)
    print(f"\nStatistics: {stats['total_nodes']} nodes, {stats['repetition_ratio']:.1%} repetition")
    
    # Demo 2: Tree with number alphabet
    print("\n=== Demo 2: Number Alphabet (1-n) ===")
    number_generator = RepetitiveTreeGenerator(
        max_depth=4,
        max_branching_factor=2,
        num_subtree_templates=2,
        repetition_probability=0.8,
        template_depth_range=(2, 3),
        alphabet_type="numbers",
        alphabet_size=8,  # Use numbers 1-8
        seed=42
    )
    
    number_tree = number_generator.generate_tree()
    print("Generated Tree with Number Labels:")
    print(number_tree)
    
    print("\nTemplates:")
    number_generator.print_templates()
    
    stats = number_generator.get_statistics(number_tree)
    print(f"\nStatistics: {stats['total_nodes']} nodes, {stats['repetition_ratio']:.1%} repetition")
    
    # Demo 3: Demonstrate label preservation in repetitive subtrees
    print("\n=== Demo 3: Label Preservation in Repetitive Subtrees ===")
    preservation_generator = RepetitiveTreeGenerator(
        max_depth=5,
        max_branching_factor=3,
        num_subtree_templates=1,  # Only one template for clear demonstration
        repetition_probability=0.9,  # High repetition
        template_depth_range=(3, 3),  # Fixed template size
        alphabet_type="letters",
        alphabet_size=10,
        seed=123
    )
    
    preservation_tree = preservation_generator.generate_tree()
    print("Tree showing label preservation in repeated subtrees:")
    print(preservation_tree)
    
    # Save examples
    letter_generator.save_tree_to_json(letter_tree, "letter_tree.json")
    number_generator.save_tree_to_json(number_tree, "number_tree.json")
    preservation_generator.save_tree_to_json(preservation_tree, "preservation_tree.json")
    print("\nTrees saved to JSON files.")


if __name__ == "__main__":
    main()