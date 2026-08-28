#!/bin/bash

function run_j {
        python3 j2.py -v -i \
                muscle.fasta \
                mafft.fasta
}

run_j

exit
