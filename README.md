# README

This repository contains code and methods for the publication:

***Nucleotide resolution saturation genome editing of BAP1 to resolve variants of uncertain significance and to functionally classify somatic variants***

## `waters_bap1_sge_analysis.R`

This is code written in [R](https://cran.r-project.org/) to convert next generation sequencing counts into functional scores for BAP1 SGE. 

Broad functions for `waters_bap1_sge_analysis.R` are:

* import SGE counts, which are held as separate files for each sample output from the [QUANTS](https://github.com/cancerit/QUANTS) pipeline, into count matrices for each target region
* annotate variants with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) outputs and [VaLiAnT](https://github.com/cancerit/VaLiAnT) meta-data outputs
* filter counts and create normalisation matrices (synonymous and intronic variants)
* run [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) with adjusted size factors to calculate raw log fold changes (LFCs) and statistics
* produce QC plots for replicate consistency, dispersion, sample correlation and editing efficiency
* median scale and combine library a and library b LFCs to generate a combined LFC
* calculate adjusted z-scores, standard error and p-value 
* FDR correction of p-value
* data frame cleanup and annotation refinement

*Note: PATHS are currently (March 2023) hard coded, with specific steps for BAP1 SGE annotation only.*