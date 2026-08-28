#!/bin/bash

function run_j {
        python3 jaccard2.py -v -i \
                clustal_g_o.fasta \
                muscle_g_o.fasta
}

run_j

exit
