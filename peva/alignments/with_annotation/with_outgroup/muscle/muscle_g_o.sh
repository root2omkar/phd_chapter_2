#!/bin/bash

set -e

echo "========================================"
echo "Running PEVA analysis: muscle_g_o"
echo "========================================"

echo ""
echo "1. Running PEVA model..."
peva -m muscle_g_o.fasta

echo ""
echo "2. Converting FASTA to DNA format..."
peva convert -m muscle_g_o.fasta -t dna

echo ""
echo "3. Converting to TNT format..."
peva convert -m muscle_g_o.fasta -t dna --format tnt -o muscle_g_o

echo ""
echo "4. Running PEVA search..."
peva search -m muscle_g_o.nexus --seed 50 -o muscle_g_o

echo ""
echo "5. Scoring trees..."
peva score -m muscle_g_o.nexus -g muscle_g_o_trees.nwk

# echo ""
# echo "6. Running PEVA search with bootstrap..."
# peva search -m muscle_g_o.nexus --bb -o muscle_g_o

echo ""
echo "7. Optimizing trees..."
peva optimize -m muscle_g_o.nexus -g muscle_g_o_trees.nwk -o muscle_g_o

echo ""
echo "8. Categorizing results and generating SVG..."
peva categorize -m muscle_g_o.nexus -g muscle_g_o_trees.nwk -o muscle_g_o --svg

echo ""
echo "========================================"
echo "PEVA analysis completed successfully!"
echo "========================================"
