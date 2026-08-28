#!/usr/bin/env python3

from Bio import Phylo
import glob
import os

# Find all IQ-TREE files
files = glob.glob("*.treefile")

trees = []
tip_sets = []

# Read trees and collect taxa
for f in files:
    tree = Phylo.read(f, "newick")
    trees.append((f, tree))

    tips = {t.name for t in tree.get_terminals()}
    tip_sets.append(tips)

# Find shared taxa across all trees
shared_tips = set.intersection(*tip_sets)

print("Number of shared taxa:", len(shared_tips))

# Prune trees
for filename, tree in trees:

    for terminal in list(tree.get_terminals()):
        if terminal.name not in shared_tips:
            tree.prune(terminal)

    outname = filename.replace(
        ".treefile",
        "_pruned.treefile"
    )

    Phylo.write(tree, outname, "newick")

print("Finished pruning all trees.")
