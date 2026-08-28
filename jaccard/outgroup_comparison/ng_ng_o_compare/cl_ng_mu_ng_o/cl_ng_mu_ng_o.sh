#!/bin/bash

function run_j {
        python3 jaccard1.py -v -i \
                clustal_ng_o.fasta \
                muscle_ng_o.fasta
}

run_j

exit
