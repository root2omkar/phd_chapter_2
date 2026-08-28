#ReadMe 

Below is the directory structure for the `dataset` directory

```
(base) omarne@CCI3QD19XJALT dataset % tree

├── metadata.csv
├── metadata.xlsx
├── ReadMe.md
├── small_data_fasta
│   ├── DQ022305.fasta
│   ├── DQ084199.fasta
│   ├── GQ153539.fasta
│   ├── GQ153540.fasta
│   ├── GQ153541.fasta
│   ├── GQ153542.fasta
│   ├── GQ153543.fasta
│   ├── GQ153544.fasta
│   ├── GQ153545.fasta
│   ├── GQ153546.fasta
│   ├── GQ153547.fasta
│   ├── GQ153548.fasta
│   ├── HM211101.fasta
│   ├── KY352407.fasta
│   ├── MG772933.fasta
│   ├── MG772934.fasta
│   ├── MN908947.fasta
│   ├── MT040333.fasta
│   ├── MT040334.fasta
│   ├── MT040335.fasta
│   ├── MT040336.fasta
│   ├── MT072864.fasta
│   ├── MW183957.fasta
│   ├── MW494005.fasta
│   ├── MW901580.fasta
│   ├── MZ131148.fasta
│   ├── MZ160298.fasta
│   ├── MZ412499.fasta
│   ├── NC_009021.fasta
│   ├── NC_025217.fasta
│   ├── OD923885.fasta
│   ├── OK177352.fasta
│   ├── OK204161.fasta
│   ├── OK221016.fasta
│   ├── OK433428.fasta
│   ├── OL714154.fasta
│   ├── OL926066.fasta
│   ├── OL965748.fasta
│   ├── ON923618.fasta
│   ├── OP440368.fasta
│   ├── OP578130.fasta
│   ├── OQ558991.fasta
│   ├── OQ896588.fasta
│   ├── OR020535.fasta
│   ├── OR563456.fasta
│   ├── OR689160.fasta
│   ├── OR983979.fasta
│   ├── OU470855.fasta
│   ├── OV948634.fasta
│   ├── OX308112.fasta
│   ├── OX309140.fasta
│   ├── OY784472.fasta
│   ├── OY792307.fasta
│   ├── OZ014103.fasta
│   ├── OZ061529.fasta
│   ├── OZ223682.fasta
│   ├── PP023772.fasta
│   ├── PP030571.fasta
│   ├── PP435502.fasta
│   ├── PP601887.fasta
│   ├── PP615923.fasta
│   ├── PP630509.fasta
│   ├── PQ017244.fasta
│   ├── PQ121865.fasta
│   ├── PQ279056.fasta
│   ├── PQ279092.fasta
│   ├── PQ304301.fasta
│   ├── PQ375771.fasta
│   ├── PQ391889.fasta
│   ├── PQ899058.fasta
│   ├── PQ972523.fasta
│   ├── PV239625.fasta
│   ├── PV497551.fasta
│   ├── PX059619.fasta
│   ├── PX393327.fasta
│   └── PX448754.fasta
├── small_data.fasta
├── small_genbank
│   ├── DQ022305.gb
│   ├── DQ084199.gb
│   ├── GQ153539.gb
│   ├── GQ153540.gb
│   ├── GQ153541.gb
│   ├── GQ153542.gb
│   ├── GQ153543.gb
│   ├── GQ153544.gb
│   ├── GQ153545.gb
│   ├── GQ153546.gb
│   ├── GQ153547.gb
│   ├── GQ153548.gb
│   ├── HM211101.gb
│   ├── KY352407.gb
│   ├── MG772933.gb
│   ├── MG772934.gb
│   ├── MN908947.gb
│   ├── MT040333.gb
│   ├── MT040334.gb
│   ├── MT040335.gb
│   ├── MT040336.gb
│   ├── MT072864.gb
│   ├── MW183957.gb
│   ├── MW494005.gb
│   ├── MW901580.gb
│   ├── MZ131148.gb
│   ├── MZ160298.gb
│   ├── MZ412499.gb
│   ├── NC_009021.gb
│   ├── NC_025217.gb
│   ├── OD923885.gb
│   ├── OK177352.gb
│   ├── OK204161.gb
│   ├── OK221016.gb
│   ├── OK433428.gb
│   ├── OL714154.gb
│   ├── OL926066.gb
│   ├── OL965748.gb
│   ├── ON923618.gb
│   ├── OP440368.gb
│   ├── OP578130.gb
│   ├── OQ558991.gb
│   ├── OQ896588.gb
│   ├── OR020535.gb
│   ├── OR563456.gb
│   ├── OR689160.gb
│   ├── OR983979.gb
│   ├── OU470855.gb
│   ├── OV948634.gb
│   ├── OX308112.gb
│   ├── OX309140.gb
│   ├── OY784472.gb
│   ├── OY792307.gb
│   ├── OZ014103.gb
│   ├── OZ061529.gb
│   ├── OZ223682.gb
│   ├── PP023772.gb
│   ├── PP030571.gb
│   ├── PP435502.gb
│   ├── PP601887.gb
│   ├── PP615923.gb
│   ├── PP630509.gb
│   ├── PQ017244.gb
│   ├── PQ121865.gb
│   ├── PQ279056.gb
│   ├── PQ279092.gb
│   ├── PQ304301.gb
│   ├── PQ375771.gb
│   ├── PQ391889.gb
│   ├── PQ899058.gb
│   ├── PQ972523.gb
│   ├── PV239625.gb
│   ├── PV497551.gb
│   ├── PX059619.gb
│   ├── PX393327.gb
│   ├── PX448754.gb
│   └── small_data.gb
└── terminals.xlsx

```

`metadata.csv` and `metadata.xlsx` are the spreadsheets that include metadata of all the terminals used in the analysis. 

`small_data_fasta` is the folder that includes individual GenBank files of the sequences in the dataset from the NCBI GenBank.

`terminals.xlsx` is the spreadsheet that includes the terminals specifying the outgroup and ingroup. 

`small_genbank` is the dataset that includes outgroup sequences and ingroup SARS-CoV-2 sequences used in the analysis.

`small_data_fasta` is the folder that contains individual fasta files of the sequences in the dataset.


