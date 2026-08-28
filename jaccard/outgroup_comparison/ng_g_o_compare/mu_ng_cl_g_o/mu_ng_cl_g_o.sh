#!/bin/bash

function run_j {
        python3 jaccard2.py -v -i \
                muscle_ng_o.fasta \
                clustal_g_o.fasta
}

run_j

exit
