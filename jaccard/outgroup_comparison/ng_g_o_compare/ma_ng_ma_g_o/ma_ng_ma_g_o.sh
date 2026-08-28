#!/bin/bash

function run_j {
        python3 jaccard2.py -v -i \
                mafft_ng_o.fasta \
                mafft_g_o.fasta
}

run_j

exit
