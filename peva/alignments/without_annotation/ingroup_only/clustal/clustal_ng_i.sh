#!/bin/bash

set -e

echo "========================================"
echo "Running PEVA analysis: clustal_ng_i"
echo "========================================"

echo ""
echo "1. Running PEVA model..."
peva -m clustal_ng_i.fasta

echo ""
echo "2. Converting FASTA to DNA format..."
peva convert -m clustal_ng_i.fasta -t dna

echo ""
echo "3. Converting to TNT format..."
peva convert -m clustal_ng_i.fasta -t dna --format tnt -o clustal_ng_i

echo ""
echo "4. Running PEVA search..."
peva search -m clustal_ng_i.nexus --seed 50 -o clustal_ng_i

echo ""
echo "5. Scoring trees..."
peva score -m clustal_ng_i.nexus -g clustal_ng_i_trees.nwk

# echo ""
# echo "6. Running PEVA search with bootstrap..."
# peva search -m clustal_ng_i.nexus --bb -o clustal_ng_i

echo ""
echo "7. Optimizing trees..."
peva optimize -m clustal_ng_i.nexus -g clustal_ng_i_trees.nwk -o clustal_ng_i

echo ""
echo "8. Categorizing results and generating SVG..."
peva categorize -m clustal_ng_i.nexus -g clustal_ng_i_trees.nwk -o clustal_ng_i --svg

echo ""
echo "========================================"
echo "PEVA analysis completed successfully!"
echo "========================================"
