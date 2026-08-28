#!/bin/bash

function run_j {
        python3 jaccard1.py -v -i \
                mafft_ng_i.fasta \
                muscle_ng_i.fasta
}

run_j

exit
