[toc]

# Comparative Phylogenetics of SARS-CoV-2: Alignment Strategies, Support Metrics, and Variant Delimitation

## Overview
This project investigates the evolution and phylogenetics of coronaviruses, focusing on SARS-CoV-2 (the causative agent of COVID-19). We will retrieve global, open-source SARS-CoV-2 genomic data from **Nextstrain** and **NCBI GenBank**, apply stringent quality filters, and generate curated datasets for downstream analyses.

Two types of datasets will be produced:
- **Small Database**: One representative, high-quality genome (>29,000 bp, low ambiguity, low missing data, available in GenBank) per variant.
- **Large Database**: 10–100 sequences per variant, selected to maximize uniqueness, diversity, and coverage across localities and sequencing conditions.

Using these databases, we will conduct **four sub-projects**, each testing hypotheses about the effects of alignment strategies, phylogenetic support measures, and variant delimitation.

## Sub-Projects

### Sub-Project 1: Alignment Strategy and Phylogenetic Impact
**Lead:** Omkar Marne
**Dataset:** Small Database

This project compares the impact of annotation-based versus whole-genome-based alignments of the SARS-CoV-2 spike glycoprotein on downstream phylogenetic inference.

**Workflow:**
1. Extract spike sequences from GenBank annotations and align with MAFFT, MUSCLE, Geneious, and Clustal.
2. Align complete genomes, extract the spike region based on the Wuhan-Hu-1 reference annotation, and compare results.
3. Normalize alignments to include identical sets of terminals.
4. Conduct phylogenetic analyses using IQ-TREE (version 2 or 3; Omkar to research best choice).
5. Compare alignments (Jaccard distances), trees (YBYRÁ), mutation rates (`--rate`, `--mlrate` in IQ-TREE), and inferred evolutionary models.

**Core Question:**
Does ignoring gene boundaries during whole-genome alignment (and trimming to the spike region afterwards) create biased alignments, mutation rate estimates, and tree topologies compared to annotation-based approaches?

### Sub-Project 2: Bootstrap and Likelihood Ratio Correlation
**Lead:** Anastasiia Duchenko
**Dataset:** Small Database

This project examines how bootstrap proportions (in different flavors) correlate with maximum likelihood support values (likelihood ratios).

**Workflow:**
1. For selected alignments (or all if time allows), calculate branch supports using likelihood ratios and multiple bootstrap methods in IQ-TREE.
2. Write a Python script to extract support values per branch and compile a database of clade-specific support vs. bootstrap measures.
3. Perform statistical analysis to assess correlations.

**Core Question:**
Do any bootstrap methods meaningfully correlate with likelihood ratio support, or do they remain poor proxies, as suggested in prior work (Machado et al. 2022, *Cladistics*)?

### Sub-Project 3: Bootstrap and Alignment Robustness
**Lead:** Brad
**Dataset:** Small Database

This project tests how well bootstrap values predict robustness to perturbations in sequence alignment.

**Workflow:**
1. Use Omkar’s alignments (or replicate them independently).
2. For each alignment, calculate bootstrap values (single chosen flavor) and likelihood ratio supports.
3. Develop a script to compare trees across alignments and calculate clade frequency (how often each clade appears across alignments).
4. Compare clade frequency with bootstrap and likelihood ratio support.

**Core Question:**
To what extent do bootstrap values capture robustness of clades to variation in alignment strategy?

### Sub-Project 4: Variant Delimitation and Synapomorphy Mapping
**Lead:** Omkar Marne
**Dataset:** Large Database

This project applies best-practice phylogenetics to evaluate the reciprocal monophyly and defining mutations of SARS-CoV-2 variants.

**Workflow:**
1. Annotate genomes and generate high-quality alignments using **FLAVi-CoV**.
2. Infer phylogenies using **parsimony**, rather than maximum likelihood.
3. Use **YBYRÁ** to map synapomorphies defining clades.
4. Assess whether variants form reciprocally monophyletic groups and whether inferred synapomorphies match known lineage-defining mutations from Nextstrain, WHO, and Pango.

**Core Question:**
Would rigorous phylogenetic analysis produce a tree and clade-defining mutations that challenge current lineage naming and classification systems?

## Sub-Project Summary Table

| Sub-Project                            | Lead       | Dataset  | Main Goal                                                    |
| -------------------------------------- | ---------- | -------- | ------------------------------------------------------------ |
| 1. Alignment Strategy & Impact         | Omkar      | Small DB | Test annotation-based vs. whole-genome alignment effects on phylogeny. |
| 2. Bootstrap vs. Likelihood            | Anastasiia | Small DB | Assess correlation between bootstrap proportions and likelihood ratios. |
| 3. Bootstrap & Alignment Robustness    | Brad       | Small DB | Test bootstrap’s ability to predict robustness across multiple alignments. |
| 4. Variant Delimitation & Synapomorphy | Omkar      | Large DB | Evaluate reciprocal monophyly of variants and identify lineage-defining mutations. |

## Research Questions

1. Does IQ-TREE v2 (Coronavirus edition) provide advantages over IQ-TREE v3 for coronavirus datasets? This is an important question that will influence all of our work. It should be answered first, via review of the specialized literature and the program's manuals.
2. How do annotation-based alignments differ from whole-genome-based alignments in their impact on alignment quality, mutation rate estimates, tree topologies, and evolutionary models?
3. Do bootstrap proportions (in their various forms) correlate with likelihood ratio support, or do they remain unreliable measures of clade support?
4. How robust are bootstrap and likelihood ratio supports to variation in sequence alignment strategies?
5. Are SARS-CoV-2 variants reciprocally monophyletic when analyzed under rigorous phylogenetic practices?
6. Do inferred synapomorphies match those used by Nextstrain, WHO, and Pango lineage systems to define variants?
7. Would adopting best-practice phylogenetic approaches significantly alter classification and mutation-based lineage definitions?

## Key References

- Jacob Machado, D. (2015). *YBYRÁ facilitates comparison of large phylogenetic trees.* **BMC Bioinformatics**, 16(1), 204.
- Jacob Machado, D., Portella de Luna Marques, F., Jiménez‐Ferbans, L., & Grant, T. (2022). *An empirical test of the relationship between the bootstrap and likelihood ratio support in maximum likelihood phylogenetic analysis.* **Cladistics**, 38(3), 392–401.
- Jacob Machado, D., Castroviejo-Fisher, S., & Grant, T. (2021). *Evidence of absence treated as absence of evidence: effects of missing data on maximum likelihood analysis.* **Mol. Phylogenet. Evol.**, 154, 106966.
- Wenzel, J. (2020). *Origins of SARS‐CoV‐1 and SARS‐CoV‐2 are often poorly explored in leading publications.* **Cladistics**, 36(4), 374–379.
- Jacob Machado, D., de Bernadi Schneider, A., Guirales, S., & Janies, D. A. (2020). *FLAVi: an enhanced annotator for viral genomes of Flaviviridae.* **Viruses**, 12(8), 892.
- Jacob Machado, D., Scott, R., Guirales, S., & Janies, D. A. (2021). *Fundamental evolution of all Orthocoronavirinae including SARS-CoV, MERS-CoV, and SARS-CoV-2.* **Cladistics**, 37(5), 461–488.