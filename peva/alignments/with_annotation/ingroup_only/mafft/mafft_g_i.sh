#!/bin/bash

set -e

echo "========================================"
echo "Running PEVA analysis: mafft_g_i"
echo "========================================"

echo ""
echo "1. Running PEVA model..."
peva -m mafft_g_i.fasta

echo ""
echo "2. Converting FASTA to DNA format..."
peva convert -m mafft_g_i.fasta -t dna

echo ""
echo "3. Converting to TNT format..."
peva convert -m mafft_g_i.fasta -t dna --format tnt -o mafft_g_i

echo ""
echo "4. Running PEVA search..."
peva search -m mafft_g_i.nexus --seed 50 -o mafft_g_i

echo ""
echo "5. Scoring trees..."
peva score -m mafft_g_i.nexus -g mafft_g_i_trees.nwk

# echo ""
# echo "6. Running PEVA search with bootstrap..."
# peva search -m mafft_g_i.nexus --bb -o mafft_g_i

echo ""
echo "7. Optimizing trees..."
peva optimize -m mafft_g_i.nexus -g mafft_g_i_trees.nwk -o mafft_g_i

echo ""
echo "8. Categorizing results and generating SVG..."
peva categorize -m mafft_g_i.nexus -g mafft_g_i_trees.nwk -o mafft_g_i --svg

echo ""
echo "========================================"
echo "PEVA analysis completed successfully!"
echo "========================================"
