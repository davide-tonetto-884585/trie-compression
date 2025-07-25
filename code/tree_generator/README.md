# Repetitive Tree Generator

A Python-based tool for generating trees with high levels of repetitiveness by recursively copying subtrees. This generator allows fine-grained control over the repetition patterns and is useful for testing algorithms that work with repetitive tree structures.

## Features

- **Configurable Repetition**: Control the level of repetitiveness through various parameters
- **Template-based Generation**: Uses predefined subtree templates that are copied throughout the tree
- **Flexible Parameters**: Adjust tree depth, branching factor, number of templates, and repetition probability
- **Statistics and Analysis**: Get detailed statistics about the generated trees
- **Export Capabilities**: Save trees to JSON format for further analysis
- **Reproducible Results**: Use seeds for consistent tree generation

## Installation

No external dependencies required. The generator uses only Python standard library modules:

```python
import random
import json
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from copy import deepcopy
```

## Quick Start

```python
from repetitive_tree_generator import RepetitiveTreeGenerator

# Create a generator with default parameters
generator = RepetitiveTreeGenerator()

# Generate a tree
tree = generator.generate_tree()

# Print the tree structure
print(tree)

# Get statistics
stats = generator.get_statistics(tree)
print(f"Tree has {stats['total_nodes']} nodes with {stats['repetition_ratio']:.1%} repetition")
```

## Parameters

### RepetitiveTreeGenerator Parameters

- **`max_depth`** (int, default=6): Maximum depth of the generated tree
- **`max_branching_factor`** (int, default=3): Maximum number of children per node
- **`num_subtree_templates`** (int, default=3): Number of different subtree templates to create
- **`repetition_probability`** (float, default=0.7): Probability of reusing an existing template vs creating new subtree
- **`template_depth_range`** (tuple, default=(2, 4)): (min, max) depth range for subtree templates
- **`alphabet_type`** (str, default="letters"): Type of alphabet to use for node labels ("letters" or "numbers")
- **`alphabet_size`** (int, default=26): Size of the alphabet (max 26 for letters, any positive int for numbers)
- **`seed`** (int, optional): Random seed for reproducible generation

## Usage Examples

### Basic Usage

```python
# Generate a tree with high repetition using letter alphabet
high_rep_generator = RepetitiveTreeGenerator(
    max_depth=6,
    max_branching_factor=2,
    num_subtree_templates=2,
    repetition_probability=0.9,
    alphabet_type="letters",
    alphabet_size=10,
    seed=42
)

tree = high_rep_generator.generate_tree()
print(f"Generated tree with {tree.size()} nodes")
```

### Alphabet Configuration

```python
# Using letter alphabet (a-z)
letter_generator = RepetitiveTreeGenerator(
    alphabet_type="letters",
    alphabet_size=5,  # Use only a, b, c, d, e
    seed=42
)

# Using number alphabet (1-n)
number_generator = RepetitiveTreeGenerator(
    alphabet_type="numbers",
    alphabet_size=10,  # Use numbers 1-10
    seed=42
)

letter_tree = letter_generator.generate_tree()
number_tree = number_generator.generate_tree()
```

### Low Repetition (More Diverse)

```python
low_rep_generator = RepetitiveTreeGenerator(
    max_depth=5,
    num_subtree_templates=5,
    repetition_probability=0.3,
    seed=123
)

tree = low_rep_generator.generate_tree()
```

### Analyzing Templates

```python
generator = RepetitiveTreeGenerator(seed=456)
tree = generator.generate_tree()

# Print all subtree templates
generator.print_templates()

# Get detailed statistics
stats = generator.get_statistics(tree)
print(f"Template usage: {stats['template_usage']}")
print(f"Repetition ratio: {stats['repetition_ratio']:.2%}")
```

### Saving Trees

```python
generator = RepetitiveTreeGenerator()
tree = generator.generate_tree()

# Save to JSON file
generator.save_tree_to_json(tree, "my_tree.json")

# The JSON contains the complete tree structure with labels and node IDs
```

## Tree Structure

### TreeNode Class

Each node in the tree is represented by a `TreeNode` object with the following attributes:

- **`label`**: String identifier for the node
- **`children`**: List of child TreeNode objects
- **`node_id`**: Unique integer identifier (optional)

### Node Methods

- **`add_child(child)`**: Add a child node
- **`is_leaf()`**: Check if the node is a leaf
- **`size()`**: Get total number of nodes in subtree
- **`depth()`**: Get depth of subtree
- **`to_dict()`**: Convert to dictionary representation

## Statistics

The generator provides comprehensive statistics about generated trees:

```python
stats = generator.get_statistics(tree)
# Returns:
# {
#     'total_nodes': 34,
#     'tree_depth': 6,
#     'num_templates': 4,
#     'template_usage': {0: 2, 1: 3, 2: 2, 3: 2},
#     'total_template_usage': 9,
#     'repetition_ratio': 0.26
# }
```

## Configuration Examples

### High Repetition Configuration

```python
# For maximum repetitiveness with letter alphabet
generator = RepetitiveTreeGenerator(
    max_depth=6,
    max_branching_factor=2,
    num_subtree_templates=2,      # Few templates
    repetition_probability=0.95,   # High reuse probability
    template_depth_range=(3, 5),  # Larger templates
    alphabet_type="letters",
    alphabet_size=8,              # Use a-h
    seed=42
)
```

### Balanced Configuration

```python
# For moderate repetitiveness with number alphabet
generator = RepetitiveTreeGenerator(
    max_depth=5,
    max_branching_factor=3,
    num_subtree_templates=3,
    repetition_probability=0.6,
    template_depth_range=(2, 4),
    alphabet_type="numbers",
    alphabet_size=15,             # Use numbers 1-15
    seed=42
)
```

### Low Repetition Configuration

```python
# For minimal repetitiveness (more diverse) with letters
generator = RepetitiveTreeGenerator(
    max_depth=5,
    max_branching_factor=3,
    num_subtree_templates=6,      # Many templates
    repetition_probability=0.3,   # Low reuse probability
    template_depth_range=(2, 3),  # Smaller templates
    alphabet_type="letters",
    alphabet_size=20,             # Use a-t
    seed=42
)
```

## Understanding Repetition

### How Templates Work

1. **Template Generation**: The generator first creates a set of random subtree templates using the configured alphabet
2. **Template Selection**: During tree generation, it probabilistically chooses between:
   - Copying an existing template (creates repetition)
   - Generating a new unique subtree
3. **Template Balancing**: Less-used templates are preferred to ensure balanced usage
4. **Depth-based Probability**: Repetition probability increases with tree depth

### Label Preservation in Repetitive Subtrees

**Key Feature**: When a subtree template is copied to a new position in the tree, **all labels within that subtree remain identical** to the original template. This ensures true structural repetition.

```python
# Example: Template with labels [a, b, c]
# When copied multiple times, each copy will have exactly [a, b, c]
# This creates genuine repetitive patterns in the tree structure

generator = RepetitiveTreeGenerator(
    num_subtree_templates=1,
    repetition_probability=0.9,
    alphabet_type="letters",
    alphabet_size=5
)

tree = generator.generate_tree()
# Result: Multiple subtrees with identical label sequences
```

### Alphabet Systems

- **Letter Alphabet**: Uses lowercase letters a-z (configurable size up to 26)
- **Number Alphabet**: Uses numbers 1-n (configurable size, any positive integer)
- **Label Extension**: If the tree requires more labels than the alphabet size, the generator automatically creates combinations (e.g., aa, ab, ac...)

### Repetition Ratio

The repetition ratio is calculated as:
```
repetition_ratio = total_template_usage / total_nodes
```

A higher ratio indicates more repetitive structure.

## Files

- **`repetitive_tree_generator.py`**: Main generator implementation
- **`example_usage.py`**: Comprehensive usage examples
- **`README.md`**: This documentation

## Running Examples

```bash
# Run the main demo
python repetitive_tree_generator.py

# Run comprehensive examples
python example_usage.py
```

## Use Cases

- **Algorithm Testing**: Test tree algorithms with repetitive structures
- **Compression Research**: Generate trees for testing tree compression algorithms
- **Pattern Analysis**: Study repetitive patterns in tree structures
- **Performance Benchmarking**: Create controlled test cases with varying repetition levels
- **Data Structure Research**: Analyze the impact of repetition on tree operations

## Customization

The generator can be easily extended:

- **Custom Node Types**: Modify the `TreeNode` class for specific data
- **Different Label Schemes**: Customize the `_generate_label` method
- **Advanced Selection**: Modify `_select_template` for different template selection strategies
- **Custom Statistics**: Add new metrics to the `get_statistics` method

## License

This project is open source and available under the MIT License.