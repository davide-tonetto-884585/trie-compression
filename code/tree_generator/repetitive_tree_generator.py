import random
import json
import configparser
import os
from tabnanny import verbose
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
        """Calculate the total number of nodes in this subtree (iterative)."""
        stack = [self]
        total_size = 0
        
        while stack:
            node = stack.pop()
            total_size += 1
            stack.extend(node.children)
        
        return total_size

    def depth(self) -> int:
        """Calculate the depth of this subtree (iterative)."""
        if self.is_leaf():
            return 1
        
        # Stack contains tuples of (node, current_depth)
        stack = [(self, 1)]
        max_depth = 1
        
        while stack:
            node, current_depth = stack.pop()
            max_depth = max(max_depth, current_depth)
            
            for child in node.children:
                stack.append((child, current_depth + 1))
        
        return max_depth

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
        max_branching_factor: int = 3,
        num_subtree_templates: int = 3,
        repetition_probability: float = 0.7,
        template_depth_range: Tuple[int, int] = (2, 4),
        alphabet_type: str = "letters",
        alphabet_size: int = 26,
        max_nodes: Optional[int] = None,
        seed: Optional[int] = None,
    ):
        """
        Initialize the repetitive tree generator.

        Args:
            max_branching_factor: Maximum number of children per node
            num_subtree_templates: Number of different subtree templates to create
            repetition_probability: Probability of reusing an existing template vs creating new
            template_depth_range: (min, max) depth range for subtree templates
            alphabet_type: "letters" for a-z or "numbers" for 1-n
            alphabet_size: Size of alphabet (max 26 for letters, any positive int for numbers)
            max_nodes: Maximum number of nodes in the generated tree (None for no limit)
            seed: Random seed for reproducible generation
        """
        self.max_branching_factor = max_branching_factor
        self.num_subtree_templates = num_subtree_templates
        self.repetition_probability = repetition_probability
        self.template_depth_range = template_depth_range
        self.alphabet_type = alphabet_type
        self.alphabet_size = (
            min(alphabet_size, 26) if alphabet_type == "letters" else alphabet_size
        )
        self.max_nodes = max_nodes

        if seed is not None:
            random.seed(seed)

        # Storage for subtree templates
        self.subtree_templates: List[TreeNode] = []
        self.template_usage_count: List[int] = []
        self.template_sizes: List[int] = []  # Cache template sizes for performance
        self.next_node_id = 0
        self.current_node_count = 0  # Track current number of nodes in the tree

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

    def _get_rand_label(self) -> str:
        """Pick a random label from the alphabet."""
        return random.choice(self.alphabet)

    def _generate_subtree_templates(self):
        """Generate the initial set of subtree templates."""
        for i in range(self.num_subtree_templates):
            template_depth = random.randint(*self.template_depth_range)
            template = self._create_random_subtree(template_depth)
            self.subtree_templates.append(template)
            self.template_usage_count.append(0)
            self.template_sizes.append(template.size())  # Cache the size

    def _create_random_subtree(self, max_depth: int) -> TreeNode:
        """Create a random subtree with the given maximum depth using alphabet labels."""
        root = TreeNode(label=self._get_rand_label(), children=[])
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

    def _add_template_leaves_to_stack(self, node: TreeNode, stack: list):
        """Add leaf nodes of a template subtree to the stack for iterative processing."""
        # Check if we've exceeded max nodes before continuing
        if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
            return
            
        # Use a temporary stack to find all leaf nodes iteratively
        temp_stack = [node]
        
        while temp_stack:
            current = temp_stack.pop()
            
            if current.is_leaf():
                # This is a leaf node, add it to the main stack for processing
                stack.append(current)
            else:
                # Add children to temporary stack to continue searching for leaves
                for child in current.children:
                    temp_stack.append(child)

    def _should_use_template(self) -> bool:
        """Decide whether to use a template or create a new subtree."""
        # Use the base repetition probability without depth adjustment
        return random.random() < self.repetition_probability

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
        # Reset usage counters for new tree
        self.template_usage_count = [0] * len(self.subtree_templates)
        # Reset node counter
        self.current_node_count = 0

        root = TreeNode(label=self._get_rand_label(), children=[])
        root.node_id = self._get_next_node_id()
        self.current_node_count += 1  # Count the root node

        self._build_tree_recursive(root)
        return root

    def _build_tree_recursive(self, parent: TreeNode):
        """Iteratively build the tree with repetitive patterns."""
        # Use a stack to simulate recursion: (node, depth)
        stack = [parent]
        
        while stack:
            current_node = stack.pop()
            
            # Check if we've reached the maximum number of nodes
            if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
                break

            num_children = random.randint(1, self.max_branching_factor)

            for _ in range(num_children):
                # Check node limit before creating each child
                if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
                    break
                    
                if (
                    len(self.subtree_templates) > 0
                    and self._should_use_template()
                ):
                    # Use a template
                    template_idx = self._select_template()
                    template = self.subtree_templates[template_idx]
                    
                    # Check if adding this template would exceed node limit
                    template_size = self.template_sizes[template_idx]  # Use cached size
                    if self.max_nodes is not None and self.current_node_count + template_size > self.max_nodes:
                        # Skip this template and create a single node instead
                        child = TreeNode(label=self._get_rand_label(), children=[])
                        child.node_id = self._get_next_node_id()
                        current_node.add_child(child)
                        self.current_node_count += 1
                        stack.append(child)
                    else:
                        child = self._copy_subtree_template(template)
                        self.template_usage_count[template_idx] += 1
                        current_node.add_child(child)
                        self.current_node_count += template_size
                        
                        # Add template leaves to stack with adjusted depth
                        self._add_template_leaves_to_stack(child, stack)
                else:
                    # Create a new subtree
                    child = TreeNode(label=self._get_rand_label(), children=[])
                    child.node_id = self._get_next_node_id()
                    current_node.add_child(child)
                    self.current_node_count += 1
                    stack.append(child)

    def get_statistics(self, tree: TreeNode) -> Dict:
        """Get statistics about the generated tree."""
        tree_size = tree.size()
        template_nodes_count = sum(usage * self.template_sizes[i] for i, usage in enumerate(self.template_usage_count))
        return {
            "total_nodes": tree_size,
            "tree_depth": tree.depth(),
            "num_templates": len(self.subtree_templates),
            "template_usage": dict(enumerate(self.template_usage_count)),
            "total_template_usage": sum(self.template_usage_count),
            "template_nodes_count": template_nodes_count,
            "repetition_ratio": (
                template_nodes_count / tree_size if tree_size > 0 else 0
            ),
        }

    def save_tree_to_json(self, tree: TreeNode, filename: str):
        """Save the tree to a JSON file."""
        with open(filename, "w") as f:
            json.dump(tree.to_dict(), f, indent=2)
    
    def tree_to_balanced_parentheses(self, tree: TreeNode) -> str:
        """Convert tree to balanced parentheses string with labels (iterative)."""
        result = []
        # Stack contains tuples of (node, action) where action is 'open' or 'close'
        stack = [(tree, 'close'), (tree, 'open')]
        
        while stack:
            node, action = stack.pop()
            
            if action == 'open':
                result.append("(" + node.label)
                # Add children in reverse order so they're processed left-to-right
                for child in reversed(node.children):
                    stack.append((child, 'close'))
                    stack.append((child, 'open'))
            else:  # action == 'close'
                result.append(")")
        
        return "".join(result)
    
    def save_tree_to_parentheses(self, tree: TreeNode, filename: str):
        """Save the tree as a balanced parentheses string with labels."""
        parentheses_str = self.tree_to_balanced_parentheses(tree)
        with open(filename, "w") as f:
            f.write(parentheses_str)

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
        f"bf{config_params['max_branching_factor']}",
        f"t{config_params['num_subtree_templates']}",
        f"rp{int(config_params['repetition_probability'] * 100)}",
        f"td{config_params['template_depth_min']}-{config_params['template_depth_max']}",
        f"{config_params['alphabet_type']}{config_params['alphabet_size']}"
    ]
    
    # Add max_nodes to filename if it's set
    if config_params.get('max_nodes') is not None:
        params.append(f"mn{config_params['max_nodes']}")
    
    params.append(f"s{config_params['seed']}")
    
    return "tree_" + "_".join(params) + ".txt"


def load_config(config_file: str = "config.ini") -> Dict:
    """Load configuration from INI file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    
    tree_config = config['TreeGeneration']
    
    # Handle max_nodes parameter (optional)
    max_nodes = None
    if 'max_nodes' in tree_config:
        max_nodes_value = tree_config.getint('max_nodes')
        if max_nodes_value > 0:
            max_nodes = max_nodes_value
    
    return {
        'max_branching_factor': tree_config.getint('max_branching_factor'),
        'num_subtree_templates': tree_config.getint('num_subtree_templates'),
        'repetition_probability': tree_config.getfloat('repetition_probability'),
        'template_depth_min': tree_config.getint('template_depth_min'),
        'template_depth_max': tree_config.getint('template_depth_max'),
        'alphabet_type': tree_config.get('alphabet_type'),
        'alphabet_size': tree_config.getint('alphabet_size'),
        'seed': tree_config.getint('seed'),
        'output_directory': tree_config.get('output_directory', './generated_trees'),
        'verbose': tree_config.getboolean('verbose'),
        'max_nodes': max_nodes
    }


def main():
    """Generate tree using parameters from config.ini file."""
    print("Repetitive Tree Generator - Config-based Generation")
    print("=" * 60)
    
    # Load configuration
    try:
        config_params = load_config("config_large_repetitive.ini")
        print(f"Loaded configuration from config.ini:")
        for key, value in config_params.items():
            print(f"  {key}: {value}")
        print()
    except Exception as e:
        print(f"Error loading config.ini: {e}")
        print("Using default parameters...")
        config_params = {
            'max_branching_factor': 3,
            'num_subtree_templates': 2,
            'repetition_probability': 0.8,
            'template_depth_min': 2,
            'template_depth_max': 3,
            'alphabet_type': 'letters',
            'alphabet_size': 10,
            'seed': 42,
            'output_directory': './generated_trees',
            'verbose': False,
            'max_nodes': None
        }
    
    # Create output directory if it doesn't exist
    output_dir = config_params['output_directory']
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate tree with loaded parameters
    generator = RepetitiveTreeGenerator(
        max_branching_factor=config_params['max_branching_factor'],
        num_subtree_templates=config_params['num_subtree_templates'],
        repetition_probability=config_params['repetition_probability'],
        template_depth_range=(config_params['template_depth_min'], config_params['template_depth_max']),
        alphabet_type=config_params['alphabet_type'],
        alphabet_size=config_params['alphabet_size'],
        max_nodes=config_params['max_nodes'],
        seed=config_params['seed']
    )
    
    # Generate the tree
    print("Generating tree...")
    tree = generator.generate_tree()
    print("Tree generated.")
    
    if config_params['verbose']:
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
    print(f"  Template nodes count: {stats['template_nodes_count']}")
    print(f"  Repetition ratio: {stats['repetition_ratio']:.1%}")
    
    # Generate balanced parentheses representation
    parentheses_str = generator.tree_to_balanced_parentheses(tree)
    if config_params['verbose']:
        print(f"\nBalanced Parentheses String: {parentheses_str}")
    
    # Generate descriptive filename and save
    filename = generate_filename(config_params)
    filepath = os.path.join(output_dir, filename)
    
    generator.save_tree_to_parentheses(tree, filepath)
    print(f"\nTree saved to: {filepath}")


if __name__ == "__main__":
    main()
