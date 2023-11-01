# README

This repository contains code and methods for the publication:

***Comprehensive Nucleotide resolution saturation genome editing of BAP1 to resolve variants of uncertain significance and to functionally classify somatic variants***

<span style="color: red;">**Please note: These data are unpublished and not yet subject to peer review. We provide them as a service to the research community but they are embargoed until publication of our paper. They should not be used as the sole basis for clinical decision making.**</span>

## `waters_bap1_sge_analysis_v2.R`

This is code written in [R](https://cran.r-project.org/) to convert next generation sequencing counts into functional scores for BAP1 SGE. 

Broad functions for waters_bap1_sge_analysis_v2.R` are:

* import SGE counts, which are held as separate files for each sample output from the [QUANTS](https://github.com/cancerit/QUANTS) pipeline, into count matrices for each target region
* annotate variants with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) outputs and [VaLiAnT](https://github.com/cancerit/VaLiAnT) meta-data outputs
* filter counts and create normalisation matrices (synonymous and intronic variants)
* run [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) with adjusted size factors to calculate raw log fold changes (LFCs) and statistics
* produce QC plots for replicate consistency, dispersion, sample correlation and editing efficiency
* median scale and combine Library A and Library B LFCs to generate a functional score
* calculate adjusted p-value and FDR correction of p-value
* compute functional classifications from functional score and FDR values
* data frame clean up and annotation refinement


*Note: PATHS are currently (November 2023) hard coded, with specific steps for BAP1 SGE annotation only.*

## `extended_data_1_(sge_bap1_dataset).xlsx`

The final dataset containing SGE functional scores and functional classifications for 18,108 unique variants. 

## `extended_data_2_(sge_bap1_expanded_dataset).xlsx`

The final dataset containing SGE functional scores and functional classifications for 18,108 unique variants, but with many rows duplicated where the variant has been produced through different VaLiAnT functions or where alternative identifiers for the same editing event are valid, e.g. single base-pair deletions in polynucleotide tracts which have different identifiers but result in a single editing event at the nucleotide sequence level. 23,226 accessions in total. 

## `spliceai.tar.gz`

An updated instance of [spliceAI source code](https://github.com/Illumina/SpliceAI) which allows for multinucleotide score retrieval. 

We adapted the original spliceAI code (https://github.com/Illumina/SpliceAI) to compute the scores of multi-nucleotide variation. We also made necessary changes in the code to extract scores for reference and alternative nucleotides along with acceptor and donor gains and losses. All the spliceAI scores for different variants are generated with a window size of 500, mask value set to 0, and GRCh38 reference genome and its annotations. 

## A5 HAP1 clonal line SNPs from whole genome sequencing

Whole genome sequencing of the HAP1-A5 Cas9+ HAP1- clonal line used in SGE of BAP1 was processed through the GATK pipeline:

[https://gatk.broadinstitute.org/hc/en-us](https://gatk.broadinstitute.org/hc/en-us)

SNPs in this line can be found in the GATK output [here on FigShare](https://figshare.com/articles/dataset/HAP1_variants_called_using_GATK_used_in_SGE_of_BAP1/24480886).

## `Explanation_of_functional_score_and_DESEq2_use_in_SGE.docx`

Word document which explains the functional score and DESEq2 use in SGE.