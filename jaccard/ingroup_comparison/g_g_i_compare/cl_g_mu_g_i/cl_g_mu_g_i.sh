#!/bin/bash

function run_j {
        python3 jaccard2.py -v -i \
                clustal_g_i.fasta \
                muscle_g_i.fasta
}

run_j

exit
