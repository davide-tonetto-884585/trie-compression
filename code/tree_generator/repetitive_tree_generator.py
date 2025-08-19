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
    """Generates trees with repetitive subtree patterns using copy-paste approach."""

    def __init__(
        self,
        max_branching_factor: int = 3,
        repetition_probability: float = 0.7,
        min_subtree_depth: int = 1,
        max_subtree_depth: int = 5,
        alphabet_type: str = "letters",
        alphabet_size: int = 26,
        max_nodes: Optional[int] = None,
        seed: Optional[int] = None,
    ):
        """
        Initialize the repetitive tree generator.

        Args:
            max_branching_factor: Maximum number of children per node
            repetition_probability: Probability of copying an existing subtree vs creating new
            min_subtree_depth: Minimum depth of subtrees to extract for copying
            max_subtree_depth: Maximum depth of subtrees to extract for copying
            alphabet_type: "letters" for a-z or "numbers" for 1-n
            alphabet_size: Size of alphabet (max 26 for letters, any positive int for numbers)
            max_nodes: Maximum number of nodes in the generated tree (None for no limit)
            seed: Random seed for reproducible generation
        """
        self.max_branching_factor = max_branching_factor
        self.repetition_probability = repetition_probability
        self.min_subtree_depth = min_subtree_depth
        self.max_subtree_depth = max_subtree_depth
        self.alphabet_type = alphabet_type
        self.alphabet_size = (
            min(alphabet_size, 26) if alphabet_type == "letters" else alphabet_size
        )
        self.max_nodes = max_nodes

        if seed is not None:
            random.seed(seed)

        # Storage for tracking repetitions
        self.copied_subtrees: List[TreeNode] = []  # Store references to copied subtrees
        self.next_node_id = 0
        self.current_node_count = 0  # Track current number of nodes in the tree

        # Generate alphabet
        self._generate_alphabet()

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
    
    def _get_unique_label(self, used_labels: set) -> Optional[str]:
        """Get a unique label that hasn't been used among siblings."""
        available_labels = [label for label in self.alphabet if label not in used_labels]
        if not available_labels:
            return None  # No more unique labels available
        return random.choice(available_labels)

    def _collect_subtrees(self, root: TreeNode) -> List[TreeNode]:
        """Sample random nodes and collect subtrees rooted at them up to max depth."""
        # Collect all nodes in the tree first
        all_nodes = []
        stack = [root]
        
        while stack:
            node = stack.pop()
            all_nodes.append(node)
            stack.extend(node.children)
        
        # Sample random nodes (limit to avoid too many candidates)
        num_samples = min(50, len(all_nodes))
        if len(all_nodes) <= num_samples:
            sampled_nodes = all_nodes
        else:
            sampled_nodes = random.sample(all_nodes, num_samples)
        
        # For each sampled node, create a subtree with random depth between min and max
        subtrees = []
        for node in sampled_nodes:
            # Sample a random depth between min_subtree_depth and max_subtree_depth
            random_depth = random.randint(self.min_subtree_depth, self.max_subtree_depth)
            subtree = self._extract_subtree_with_max_depth(node, random_depth)
            subtrees.append(subtree)
        
        return subtrees

    def _extract_subtree_with_max_depth(self, root_node: TreeNode, max_depth: int) -> TreeNode:
        """Extract a subtree rooted at the given node up to the specified maximum depth."""
        if max_depth <= 0:
            return None
            
        def copy_with_depth_limit(node: TreeNode, remaining_depth: int) -> TreeNode:
            # Create a copy of the current node
            new_node = TreeNode(
                label=node.label,
                children=[],
                node_id=None  # Will be assigned later when actually used
            )
            
            # If we still have depth remaining, copy children
            if remaining_depth > 1:
                for child in node.children:
                    child_copy = copy_with_depth_limit(child, remaining_depth - 1)
                    if child_copy:
                        new_node.add_child(child_copy)
            
            return new_node
        
        return copy_with_depth_limit(root_node, max_depth)

    def _copy_subtree(self, source: TreeNode) -> TreeNode:
        """Create a deep copy of a subtree with new node IDs but preserving labels."""
        def copy_node(node: TreeNode) -> TreeNode:
            new_node = TreeNode(
                label=node.label,  # Preserve original label
                children=[],
                node_id=self._get_next_node_id(),  # New unique node ID
            )
            for child in node.children:
                new_node.add_child(copy_node(child))
            return new_node

        copied = copy_node(source)
        # Store the size of the copied subtree for accurate statistics
        copied_size = source.size()
        self.copied_subtrees.append((copied, copied_size))
        return copied

    def _should_copy_subtree(self, available_subtrees: List[TreeNode]) -> bool:
        """Decide whether to copy an existing subtree or create a new node."""
        # Only copy if we have available subtrees and meet probability threshold
        return len(available_subtrees) > 0 and random.random() < self.repetition_probability

    def _select_subtree_to_copy(self, available_subtrees: List[TreeNode]) -> TreeNode:
        """Select a random subtree from available options."""
        return random.choice(available_subtrees)

    def generate_tree(self) -> TreeNode:
        """Generate a repetitive tree using copy-paste approach."""
        # Reset counters for new tree
        self.copied_subtrees = []
        self.current_node_count = 0

        root = TreeNode(label=self._get_rand_label(), children=[])
        root.node_id = self._get_next_node_id()
        self.current_node_count += 1

        self._build_tree_iterative(root)
        return root

    def _build_tree_iterative(self, root: TreeNode):
        """Build the tree iteratively with copy-paste repetitive patterns."""
        stack = [root]
        available_subtrees = []  # Cache available subtrees
        subtree_collection_interval = 50  # Collect subtrees every N nodes
        
        while stack:
            current_node = stack.pop()
            
            # Check if we've reached the maximum number of nodes
            if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
                break

            # Periodically collect available subtrees to avoid expensive repeated collection
            if self.current_node_count % subtree_collection_interval == 0:
                available_subtrees = self._collect_subtrees(root)
            
            # Limit the number of children based on remaining node budget
            max_possible_children = self.max_branching_factor
            if self.max_nodes is not None:
                remaining_nodes = self.max_nodes - self.current_node_count
                max_possible_children = min(max_possible_children, remaining_nodes)
            
            if max_possible_children <= 0:
                break
                
            num_children = random.randint(1, max_possible_children)
            
            # Track used labels for this node's children to ensure deterministic behavior
            used_labels = set()

            for _ in range(num_children):
                # Check node limit before creating each child
                if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
                    break
                
                # Decide whether to copy an existing subtree or create a new node
                if self._should_copy_subtree(available_subtrees):
                    # Copy an existing subtree
                    source_subtree = self._select_subtree_to_copy(available_subtrees)
                    subtree_size = source_subtree.size()
                    
                    # Check if adding this subtree would exceed node limit
                    if self.max_nodes is not None and self.current_node_count + subtree_size > self.max_nodes:
                        # Skip copying and create a single node instead
                        unique_label = self._get_unique_label(used_labels)
                        if unique_label is None:
                            break  # No more unique labels available
                        child = TreeNode(label=unique_label, children=[])
                        child.node_id = self._get_next_node_id()
                        current_node.add_child(child)
                        used_labels.add(unique_label)
                        self.current_node_count += 1
                        stack.append(child)
                    else:
                        # Copy the selected subtree and ensure unique root label
                        child = self._copy_subtree(source_subtree)
                        unique_label = self._get_unique_label(used_labels)
                        if unique_label is None:
                            break  # No more unique labels available
                        child.label = unique_label  # Override the copied label to ensure uniqueness
                        current_node.add_child(child)
                        used_labels.add(unique_label)
                        self.current_node_count += subtree_size
                        
                        # Add leaf nodes of copied subtree to stack for further expansion
                        self._add_leaves_to_stack(child, stack)
                else:
                    # Create a new single node
                    unique_label = self._get_unique_label(used_labels)
                    if unique_label is None:
                        break  # No more unique labels available
                    child = TreeNode(label=unique_label, children=[])
                    child.node_id = self._get_next_node_id()
                    current_node.add_child(child)
                    used_labels.add(unique_label)
                    self.current_node_count += 1
                    stack.append(child)

    def _add_leaves_to_stack(self, node: TreeNode, stack: list):
        """Add leaf nodes of a subtree to the stack for further processing."""
        # Check if we've exceeded max nodes before continuing
        if self.max_nodes is not None and self.current_node_count >= self.max_nodes:
            return
            
        temp_stack = [node]
        leaves_added = 0
        max_leaves_to_add = 20  # Limit the number of leaves added at once
        
        while temp_stack and leaves_added < max_leaves_to_add:
            current = temp_stack.pop()
            
            if current.is_leaf():
                # This is a leaf node, add it to the main stack for processing
                stack.append(current)
                leaves_added += 1
            else:
                # Add children to temporary stack to continue searching for leaves
                temp_stack.extend(current.children)

    def get_statistics(self, tree: TreeNode) -> Dict:
        """Get statistics about the generated tree."""
        tree_size = tree.size()
        copied_nodes_count = sum(size for _, size in self.copied_subtrees)
        return {
            "total_nodes": tree_size,
            "tree_depth": tree.depth(),
            "copied_subtrees_count": len(self.copied_subtrees),
            "copied_nodes_count": copied_nodes_count,
            "repetition_ratio": (
                copied_nodes_count / tree_size if tree_size > 0 else 0
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

    def print_copy_statistics(self):
        """Print statistics about copied subtrees."""
        print("Copy-Paste Statistics:")
        print("=" * 50)
        print(f"Number of copied subtrees: {len(self.copied_subtrees)}")
        if self.copied_subtrees:
            print("\nCopied subtrees:")
            for i, (subtree, size) in enumerate(self.copied_subtrees):
                print(f"Copy {i+1} (size: {size}):")
                print(subtree)
                print("-" * 30)


def generate_filename(config_params: Dict) -> str:
    """Generate a descriptive filename based on configuration parameters."""
    params = [
        f"bf{config_params['max_branching_factor']}",
        f"rp{int(config_params['repetition_probability'] * 100)}",
        f"sd{config_params['min_subtree_depth']}-{config_params['max_subtree_depth']}",
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
        'repetition_probability': tree_config.getfloat('repetition_probability'),
        'min_subtree_depth': tree_config.getint('min_subtree_depth', 1),
        'max_subtree_depth': tree_config.getint('max_subtree_depth', 5),
        'alphabet_type': tree_config.get('alphabet_type'),
        'alphabet_size': tree_config.getint('alphabet_size'),
        'seed': tree_config.getint('seed'),
        'output_directory': tree_config.get('output_directory', './generated_trees'),
        'verbose': tree_config.getboolean('verbose'),
        'max_nodes': max_nodes
    }


def main():
    """Generate tree using parameters from config.ini file."""
    print("Repetitive Tree Generator - Copy-Paste Based Generation")
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
            'repetition_probability': 0.8,
            'min_subtree_depth': 1,
            'max_subtree_depth': 5,
            'alphabet_type': 'letters',
            'alphabet_size': 10,
            'seed': 42,
            'output_directory': './generated_trees',
            'verbose': False,
            'max_nodes': 50000
        }
    
    # Create output directory if it doesn't exist
    output_dir = config_params['output_directory']
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate tree with loaded parameters
    generator = RepetitiveTreeGenerator(
        max_branching_factor=config_params['max_branching_factor'],
        repetition_probability=config_params['repetition_probability'],
        min_subtree_depth=config_params['min_subtree_depth'],
        max_subtree_depth=config_params['max_subtree_depth'],
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
        print("\nCopy-Paste Information:")
        generator.print_copy_statistics()
    
    # Get and display statistics
    stats = generator.get_statistics(tree)
    print(f"\nTree Statistics:")
    print(f"  Total nodes: {stats['total_nodes']}")
    print(f"  Tree depth: {stats['tree_depth']}")
    print(f"  Copied subtrees count: {stats['copied_subtrees_count']}")
    print(f"  Copied nodes count: {stats['copied_nodes_count']}")
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
