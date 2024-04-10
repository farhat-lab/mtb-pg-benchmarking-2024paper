#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Purpose: Update the standard newick format output by IQ-Tree to have internal node names.

# Input Paths: 
### iqtree.treefile

# Output Files:
### iqtree.NodeNamesUpdated.treefile

# Example Usage: IQTree.UpdateNodeNames.V1.py -i iqtree.treefile -o iqtree.NodeNamesUpdated.treefile

import argparse
import ete3
from ete3 import Tree

def update_IQtree_Newick(input_path, output_path):
    # Read the tree output by IQ-tree as Format type #2 (Node name & Support values) 
    iq_tree = Tree(input_path, format=2)

    # Calculate the midpoint node and set as outgroup
    midpoint = iq_tree.get_midpoint_outgroup()
    iq_tree.set_outgroup(midpoint)

    # Rename nodes
    ## A) If an internal node has no name, give it a unique node name
    ## B) If "Bakta" or "PGAP" is found in a leaf name, remove it using split

    node_id = 0
    for node in iq_tree.traverse():
        if not node.is_leaf():
            if node.name == '':
                node.name = f"Node_{node_id}"
                node_id += 1
        else:
            if "Bakta" in node.name:
                node.name = node.name.split(".")[0]
            elif "PGAP" in node.name:
                node.name = node.name.split(".")[0]

    # Write the processed tree to a file
    iq_tree.write(format=1, outfile=output_path)

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description='Process a NEWICK tree file.')
    parser.add_argument('-i', '--input_path', type=str, required=True, help='Path to input NEWICK file output by IQ-tree')
    parser.add_argument('-o', '--output_path', type=str, required=True, help='Path to output updated NEWICK file with internal node names')
    
    # Parse arguments
    args = parser.parse_args()

    # Process the tree
    update_IQtree_Newick(args.input_path, args.output_path)

if __name__ == "__main__":
    main()


