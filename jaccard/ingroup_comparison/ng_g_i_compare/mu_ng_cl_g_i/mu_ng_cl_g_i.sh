#!/bin/bash

function run_j {
        python3 jaccard2.py -v -i \
                muscle_ng_i.fasta \
                clustal_g_i.fasta
}

run_j

exit
