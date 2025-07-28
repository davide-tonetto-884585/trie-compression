import random
import json
import configparser
import os
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class TreeNode:
    """Represents a node in the tree."""

    label: str
    children: List["TreeNode"]
    node_id: Optional[int] = None

    def __post_init__(self):
        if self.children is None:
            self.children = []

    def add_child(self, child: "TreeNode"):
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
            "label": self.label,
            "node_id": self.node_id,
            "children": [child.to_dict() for child in self.children],
        }

    def __str__(self) -> str:
        """String representation of the tree."""
        return self._str_helper(0)

    def _str_helper(self, indent: int) -> str:
        """Helper method for string representation."""
        result = "  " * indent + f"{self.label}\n"
        for child in self.children:
            result += child._str_helper(indent + 1)
        return result


class RepetitiveTreeGenerator:
    """Generates trees with repetitive subtree patterns."""

    def __init__(
        self,
        max_depth: int = 6,
        max_branching_factor: int = 3,
        num_subtree_templates: int = 3,
        repetition_probability: float = 0.7,
        template_depth_range: Tuple[int, int] = (2, 4),
        alphabet_type: str = "letters",
        alphabet_size: int = 26,
        seed: Optional[int] = None,
    ):
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
        self.alphabet_size = (
            min(alphabet_size, 26) if alphabet_type == "letters" else alphabet_size
        )

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
            self.alphabet = [chr(ord("a") + i) for i in range(self.alphabet_size)]
        elif self.alphabet_type == "numbers":
            self.alphabet = [str(i + 1) for i in range(self.alphabet_size)]
        else:
            raise ValueError(
                f"Invalid alphabet_type: {self.alphabet_type}. Use 'letters' or 'numbers'."
            )

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
                node_id=self._get_next_node_id(),  # New unique node ID
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
            if (
                current_depth > 1
                and self.subtree_templates
                and self._should_use_template(current_depth)
            ):

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
            "total_nodes": tree.size(),
            "tree_depth": tree.depth(),
            "num_templates": len(self.subtree_templates),
            "template_usage": dict(enumerate(self.template_usage_count)),
            "total_template_usage": sum(self.template_usage_count),
            "repetition_ratio": (
                sum(self.template_usage_count) / tree.size() if tree.size() > 0 else 0
            ),
        }

    def save_tree_to_json(self, tree: TreeNode, filename: str):
        """Save the tree to a JSON file."""
        with open(filename, "w") as f:
            json.dump(tree.to_dict(), f, indent=2)

    def print_templates(self):
        """Print all subtree templates."""
        print("Subtree Templates:")
        print("=" * 50)
        for i, template in enumerate(self.subtree_templates):
            print(f"Template {i} (used {self.template_usage_count[i]} times):")
            print(template)
            print("-" * 30)


def generate_filename(config_params: Dict) -> str:
    """Generate a descriptive filename based on configuration parameters."""
    params = [
        f"d{config_params['max_depth']}",
        f"bf{config_params['max_branching_factor']}",
        f"t{config_params['num_subtree_templates']}",
        f"rp{int(config_params['repetition_probability'] * 100)}",
        f"td{config_params['template_depth_min']}-{config_params['template_depth_max']}",
        f"{config_params['alphabet_type']}{config_params['alphabet_size']}",
        f"s{config_params['seed']}"
    ]
    return "tree_" + "_".join(params) + ".json"


def load_config(config_file: str = "config.ini") -> Dict:
    """Load configuration from INI file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    
    tree_config = config['TreeGeneration']
    
    return {
        'max_depth': tree_config.getint('max_depth'),
        'max_branching_factor': tree_config.getint('max_branching_factor'),
        'num_subtree_templates': tree_config.getint('num_subtree_templates'),
        'repetition_probability': tree_config.getfloat('repetition_probability'),
        'template_depth_min': tree_config.getint('template_depth_min'),
        'template_depth_max': tree_config.getint('template_depth_max'),
        'alphabet_type': tree_config.get('alphabet_type'),
        'alphabet_size': tree_config.getint('alphabet_size'),
        'seed': tree_config.getint('seed'),
        'output_directory': tree_config.get('output_directory', './generated_trees')
    }


def main():
    """Generate tree using parameters from config.ini file."""
    print("Repetitive Tree Generator - Config-based Generation")
    print("=" * 60)
    
    # Load configuration
    try:
        config_params = load_config()
        print(f"Loaded configuration from config.ini:")
        for key, value in config_params.items():
            print(f"  {key}: {value}")
        print()
    except Exception as e:
        print(f"Error loading config.ini: {e}")
        print("Using default parameters...")
        config_params = {
            'max_depth': 5,
            'max_branching_factor': 3,
            'num_subtree_templates': 2,
            'repetition_probability': 0.8,
            'template_depth_min': 2,
            'template_depth_max': 3,
            'alphabet_type': 'letters',
            'alphabet_size': 10,
            'seed': 42,
            'output_directory': './generated_trees'
        }
    
    # Create output directory if it doesn't exist
    output_dir = config_params['output_directory']
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate tree with loaded parameters
    generator = RepetitiveTreeGenerator(
        max_depth=config_params['max_depth'],
        max_branching_factor=config_params['max_branching_factor'],
        num_subtree_templates=config_params['num_subtree_templates'],
        repetition_probability=config_params['repetition_probability'],
        template_depth_range=(config_params['template_depth_min'], config_params['template_depth_max']),
        alphabet_type=config_params['alphabet_type'],
        alphabet_size=config_params['alphabet_size'],
        seed=config_params['seed']
    )
    
    # Generate the tree
    tree = generator.generate_tree()
    
    # Display tree information
    print("Generated Tree:")
    print(tree)
    print("\nSubtree Templates:")
    generator.print_templates()
    
    # Get and display statistics
    stats = generator.get_statistics(tree)
    print(f"\nTree Statistics:")
    print(f"  Total nodes: {stats['total_nodes']}")
    print(f"  Tree depth: {stats['tree_depth']}")
    print(f"  Number of templates: {stats['num_templates']}")
    print(f"  Template usage: {stats['template_usage']}")
    print(f"  Total template usage: {stats['total_template_usage']}")
    print(f"  Repetition ratio: {stats['repetition_ratio']:.1%}")
    
    # Generate descriptive filename and save
    filename = generate_filename(config_params)
    filepath = os.path.join(output_dir, filename)
    
    generator.save_tree_to_json(tree, filepath)
    print(f"\nTree saved to: {filepath}")


if __name__ == "__main__":
    main()
