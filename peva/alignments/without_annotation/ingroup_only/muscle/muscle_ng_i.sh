#!/bin/bash

set -e

echo "========================================"
echo "Running PEVA analysis: muscle_ng_i"
echo "========================================"

echo ""
echo "1. Running PEVA model..."
peva -m muscle_ng_i.fasta

echo ""
echo "2. Converting FASTA to DNA format..."
peva convert -m muscle_ng_i.fasta -t dna

echo ""
echo "3. Converting to TNT format..."
peva convert -m muscle_ng_i.fasta -t dna --format tnt -o muscle_ng_i

echo ""
echo "4. Running PEVA search..."
peva search -m muscle_ng_i.nexus --seed 50 -o muscle_ng_i

echo ""
echo "5. Scoring trees..."
peva score -m muscle_ng_i.nexus -g muscle_ng_i_trees.nwk

# echo ""
# echo "6. Running PEVA search with bootstrap..."
# peva search -m muscle_ng_i.nexus --bb -o muscle_ng_i

echo ""
echo "7. Optimizing trees..."
peva optimize -m muscle_ng_i.nexus -g muscle_ng_i_trees.nwk -o muscle_ng_i

echo ""
echo "8. Categorizing results and generating SVG..."
peva categorize -m muscle_ng_i.nexus -g muscle_ng_i_trees.nwk -o muscle_ng_i --svg

echo ""
echo "========================================"
echo "PEVA analysis completed successfully!"
echo "========================================"
