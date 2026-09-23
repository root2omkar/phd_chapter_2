# ReadMe file 

Below is the structure of the alignments directory

```
(base) omarne@CCI3QD19XJALT alignments % tree
.
├── ReadMe.md
├── with_annotation
│   ├── ingroup_only
│   │   ├── clustal
│   │   │   ├── clustal_g_i.fasta
│   │   │   └── partition.nexus
│   │   ├── mafft
│   │   │   ├── mafft_g_i.fasta
│   │   │   └── partition.nexus
│   │   └── muscle
│   │       ├── muscle_g_i.fasta
│   │       └── partition.nexus
│   └── with_outgroup
│       ├── clustal
│       │   ├── clustal_g_o.fasta
│       │   └── partition.nexus
│       ├── mafft
│       │   ├── mafft_g_o.fasta
│       │   └── partition.nexus
│       └── muscle
│           ├── muscle_g_o.fasta
│           └── partition.nexus
└── without_annotation
    ├── ingroup_only
    │   ├── clustal
    │   │   └── clustal_ng_i.fasta
    │   ├── mafft
    │   │   └── mafft_ng_i.fasta
    │   └── muscle
    │       └── muscle_ng_i.fasta
    └── with_outgroup
        ├── clustal
        │   └── clustal_ng_o.fasta
        ├── mafft
        │   └── mafft_ng_o.fasta
        └── muscle
            └── muscle_ng_o.fasta

19 directories, 19 files
(base) omarne@CCI3QD19XJALT alignments % 

```
`with_annotation` folder contains alignments with gene annotation using FLaVi-CoV annotator.

`without_annotation` folder contains alignments without gene annotation (whole genome alignments).

`with_outgroup` folder contains alignments with both ingroup and outgroup (total 69 sequences).

`ingroup_only` folder contains alignments with ingroup sequences only (total 48 sequences).

`clustal`, `mafft` and `muscle` folders contain alignments created using CLUSTAL OMEGA v1.2.4, MUSCLE v5 and MAFFT v7.525, respectively.

`partition.nexus` is a partition file that contains gene partitions.

Note: `without_annotation` alignments do not have a partition file.

`g` - gene annotated

`ng` - no gene annotated

`o` - outgroup and ingroup (total 69 sequences)

`i` - ingroup only (total 48 sequences)

### Commands

Below are the commands which I used to create the alignments.

#### Alignment using CLUSTAL OMEGA v1.2.4, MUSCLE v5 and MAFFT v7.525 with no gene annotation

I used the R code below in the RStudio tool to create the whole-genome alignment of the `small_data. fasta' dataset using CLUSTAL OMEGA v1.2.4, MUSCLE v5, and MAFFT v7.525.

```R

library(Biostrings)
library(msa)
small_data <- readDNAStringSet("small_data.fasta")
clustal_ng_o <- msa(small_data, method = "ClustalOmega", type = "dna")
muscle_ng_o  <- msa(small_data, method = "Muscle", type = "dna")
mafft_ng_o  <- msa(small_data, method = "Mafft", type = "dna")
clustal_ng_o <- DNAStringSet(clustal_ng_o)
muscle_ng_o  <- DNAStringSet(muscle_ng_o)
mafft_ng_o <- DNAStringSet(mafft_ng_o)
writeXStringSet(clustal_ng_o, "clustal_ng_o.fasta", format="fasta")
writeXStringSet(muscle_ng_o, "muscle_ng_o.fasta", format="fasta")
writeXStringSet(mafft_ng_o, "mafft_ng_o.fasta", format="fasta")


```

#### Alignment using Flavi-CoV annotator with CLUSTAL OMEGA v1.2.4, MUSCLE v5 and MAFFT v7.525 to create annotated gene alignment

I used the Flavi-CoV slurm script `run_flavi_cov.slurm` below to create the alignments. This is an example of MAFFT gene annotation alignment. I modified the code for CLUSTAL OMEGA and MUSCLE, respectively.

To download the full Flavi-CoV directory, please go to `https://gitlab.com/MachadoDJ/flavi-cov`


```bash

#!/usr/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=mafft_g_o
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=995gb
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out


# --- Go to working directory ---
cd /users/omarne/small_dataset/ || exit
pwd

# --- Initialize Conda ---
# Source the conda setup script
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh

# Activate your environment
conda activate omkar

# Check Python
which python3
python3 --version

# Test imports
python -c "import parsl; import Bio; print('All OK')"

# --- Set environment variables ---
export WORK_DIR="/users/omarne/small_dataset/"
export FASDATA="small_data.fasta"

# Load modules
# module load clustal-omega/1.2.4 
module load mafft/7.525woe anaconda3 blast/2.15.0+ hmmer
module list

export MAFFT="/apps/pkg/mafft/7.525woe/bin/mafft"
export PY3="/users/omarne/.conda/envs/omkar/bin/python3"
export SIZE=27000

export FLAVI_DIR="${WORK_DIR}/flavi-cov"
export FLAVI="${FLAVI_DIR}/flavi-cov.py"
export DATABASE="${FLAVI_DIR}/database"
export MODULE="${FLAVI_DIR}/sarscov2annotator"
export OUTPUT="output"

export AUX_DIR="${FLAVI_DIR}/auxiliary"
export ALIGNED="${AUX_DIR}/alignments.txt"
export AA2DNA="${AUX_DIR}/aa2dna.py"
export PRE2FA="${AUX_DIR}/predictions2fasta.py"
export MERGER="${AUX_DIR}/mergeFasta.py"

# --- Functions ---
function clean {
    for i in output* extraction* aligned* final_result.fasta hmms unaligned* concat* ; do
        if [ -f "${i}" ]; then rm ${i}; fi
    done
    printf "[ `ls -lh` ]\n"
}

function annotator {
    printf "1. Executing flavi-cov.py pipeline [`date`].\n"
    printf "${PY3} ${FLAVI} -f ${FASDATA} -d ${DATABASE} -m ${MODULE} -o ${OUTPUT} -s ${SIZE}\n"
    ${PY3} ${FLAVI} -f ${FASDATA} -d ${DATABASE} -m ${MODULE} -o ${OUTPUT} -s ${SIZE} > predictions.txt
    wait
    printf "1. Done [`date`].\n"
}

function extract {
    printf "2. Extracting gene fragments from predictions [`date`].\n"
    ${PY3} ${PRE2FA} -f ${FASDATA} -p predictions.txt
    wait
    printf "2. Done [`date`].\n"
}

function align {
    printf "3. Aligning UTRs and ORFs [`date`].\n"
    for i in extraction*_prot.fasta extraction*UTR*_nucl.fasta ; do
        printf "...Aligning: ${i}\n"
        ${MAFFT} --retree 1 --thread 15 ${i} > aligned_${i}
        wait
    done
    printf "3. Done [`date`].\n"
}

function aa2dna {
    printf "4. Transfering gaps from AA alignments to nucleotide sequences [`date`].\n"
    for i in E M N ORF10 ORF1ab ORF1ab-3CLpro ORF1ab-RdRp ORF1ab-endoRNAse ORF1ab-exonuclease ORF1ab-helicase ORF1ab-methyltransferase ORF1ab-nsp10 ORF1ab-nsp1 ORF1ab-nsp2 ORF1ab-nsp3 ORF1ab-nsp4 ORF1ab-nsp6 ORF1ab-nsp7 ORF1ab-nsp8 ORF1ab-nsp9 ORF3a ORF6 ORF7a ORF7b ORF8 S S-E S-FP S-HR1 S-HR2 S-NTD S-RBD S-Region10 S-Region12 S-Region1 S-Region3 S-Region5 S-Region6 S-Region8 S-TM ; do
        printf "...Current partition: ${i}\n"
        ${PY3} ${AA2DNA} -p aligned_extraction_${i}_prot.fasta -n extraction_${i}_nucl.fasta -o aligned
        wait
    done
    printf "4. Done [`date`].\n"
}

function merge {
    printf "5. Merging individual alignments into whole-genome alignment [`date`].\n"
    ${PY3} ${MERGER} -a ${ALIGNED}
    wait
    printf "5. Done [`date`].\n"
}

function organize {
    for i in hmms unaligned_nucl aligned_nucl unaligned_prot aligned_prot ; do
        if [ ! -d ${i} ]; then mkdir ${i}; else for j in ${i}/* ; do if [ -f ${j} ] ; then rm ${j} ; fi ; done ; fi
    done
    for i in aligned*nucl.fasta ; do if [ -f ${i} ] ; then mv ${i} aligned_nucl/ ; fi ; done
    for i in *nucl.fasta ; do if [ -f ${i} ] ; then mv ${i} unaligned_nucl/ ; fi ; done
    for i in aligned*prot.fasta ; do if [ -f ${i} ] ; then mv ${i} aligned_prot/ ; fi ; done
    for i in *prot.fasta ; do if [ -f ${i} ] ; then mv ${i} unaligned_prot/ ; fi ; done
    for i in *.dom.txt *.out.txt ; do if [ -f ${i} ] ; then mv ${i} hmms/ ; fi ; done
}

# --- Execute pipeline ---
clean
annotator
extract
align
aa2dna
merge
organize

# --- Deactivate Conda ---
conda deactivate


```

#### Command to create ingroup only alignments

To create ingroup-only alignments, I directly removed the outgroup sequences from the alignments and kept only the ingroup sequences in the alignments. I did this for all the above types of alignments. I used the below `awk' command to do it.

The command is for the MAFFT alignment file. I modified for CLUSTAL OMEGA and MUSCLE, respectively.

```
awk '
  BEGIN {
    remove["DQ022305_"] = 1
    remove["DQ084199_"] = 1
    remove["GQ153539_"] = 1
    remove["GQ153540_"] = 1
    remove["GQ153541_"] = 1
    remove["GQ153542_"] = 1
    remove["GQ153543_"] = 1
    remove["GQ153544_"] = 1
    remove["GQ153545_"] = 1
    remove["GQ153546_"] = 1
    remove["GQ153547_"] = 1
    remove["GQ153548_"] = 1
    remove["KY352407_"] = 1
    remove["MG772933_"] = 1
    remove["MG772934_"] = 1
    remove["MT040333_"] = 1
    remove["MT040334_"] = 1
    remove["MT040335_"] = 1
    remove["MT040336_"] = 1
    remove["MT072864_"] = 1
    remove["NC025217_"] = 1
    skipping = 0
  }
  /^>/ {
    id = substr($1, 2)
    if (id in remove) {
      skipping = 1
    } else {
      skipping = 0
      print $0
    }
    next
  }
  { if (!skipping) print $0 }
' mafft_g_o.fasta > mafft_g_i.fasta

```