#!/bin/bash

set -e

echo "========================================"
echo "Running PEVA analysis: clustal_g_o"
echo "========================================"

echo ""
echo "1. Running PEVA model..."
peva -m clustal_g_o.fasta

echo ""
echo "2. Converting FASTA to DNA format..."
peva convert -m clustal_g_o.fasta -t dna

echo ""
echo "3. Converting to TNT format..."
peva convert -m clustal_g_o.fasta -t dna --format tnt -o clustal_g_o

echo ""
echo "4. Running PEVA search..."
peva search -m clustal_g_o.nexus --seed 50 -o clustal_g_o

echo ""
echo "5. Scoring trees..."
peva score -m clustal_g_o.nexus -g clustal_g_o_trees.nwk

# echo ""
# echo "6. Running PEVA search with bootstrap..."
# peva search -m clustal_g_o.nexus --bb -o clustal_g_o

echo ""
echo "7. Optimizing trees..."
peva optimize -m clustal_g_o.nexus -g clustal_g_o_trees.nwk -o clustal_g_o

echo ""
echo "8. Categorizing results and generating SVG..."
peva categorize -m clustal_g_o.nexus -g clustal_g_o_trees.nwk -o clustal_g_o --svg

echo ""
echo "========================================"
echo "PEVA analysis completed successfully!"
echo "========================================"
