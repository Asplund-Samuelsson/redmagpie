![alt text](redmagpie.png "Reductive pentose phosphate pathway Machine-Assisted Genomic Pattern Identification and Evaluation")

# Identification of adaptations to the Calvin cycle

## About

This is a method for identification of Calvin cycle-positive microbial genomes and ranking of genetic features (Enzyme Commission numbers and Pfam families) that represent adaptations to the Calvin cycle.

## System requirements

Linux operating system (Tested on Ubuntu 18.04.5 LTS and 20.04.1 LTS)

bash 4.0 (Tested with 4.4.20(1)-release and 5.0.17(1)-release)

Python 3.7 (Tested with 3.7.6)

R ≥ 3.6.3 (Tested with 3.6.3)

GNU parallel 20161222 (Tested with 20161222)

pigz 2.4 (Tested with 2.4)

seqmagick 0.6.2 (Tested with 0.6.2)

ORFfinder 0.4.3 (Tested with 0.4.3)

FastTreeMP 2.1.8 SSE3

[DeepEC](https://bitbucket.org/kaistsystemsbiology/deepec/src/master/) (Tested with commit b7e4546)

Python libraries: sklearn, numpy, pandas

R libraries: doMC, egg, foreach, ggnewscale, ggrepel, ggtree, MidpointRooter, optparse, phytools, RColorBrewer, scales, tidyverse, viridis

## Data

This analysis used 24,706 species representative archaeal and bacterial genomes listed in GTDB release 89 (https://gtdb.ecogenomic.org/). The genomes were downloaded from NCBI in gzip-compressed nucleotide FASTA format followed by ORF annotation with stand-alone ORFfinder v.0.4.3 (https://www.ncbi.nlm.nih.gov/orffinder/) using the following command:

```
ls /path/to/ncbi/genomes/ | parallel --no-notice --jobs 16 '
   Outfile=`echo "data/orf/{}" | sed -e "s/fna.gz$/fasta/"`;
   zcat /path/to/ncbi/genomic/{} > intermediate/{};
   ORFfinder -g 11 -ml 300 -n true -in intermediate/{} -out $Outfile;
   gzip $Outfile;
   rm intermediate/{}
'
```

Note that the input nucleotide files must be named with the GTDB accession IDs, i.e. `<GTDB accession ID>.fna.gz`, which is not the case directly after download from NCBI. If populated correctly, the directory `data/orf/` should contain 24,706 gzip-compressed amino acid FASTA files with ORFs for each of the genomes. Each file must be named `<GTDB accession ID>.fasta.gz` for the analysis to work.

The analysis also relies on the GTDB core protein alignments and Pfam-A HMMs (not included; see steps 02, 03, and 05).

### Example genome ORFs and annotations

Files containing example genome ORFs and their annotations can be [downloaded from Figshare](https://doi.org/10.6084/m9.figshare.13013309) (DOI 10.6084/m9.figshare.13013309).

## Method

These are short descriptions of each step in the analysis. For details, see the corresponding scripts. Some output files have been omitted in the descriptions.

### 00: Extract accession to group association

_Script:_ `source/00_Extract_accession_to_group_association.R`

_Output:_ `intermediate/accession_taxonomy.tab`

Creates a tab-delimited table with taxonomic information for each accession in GTDB.

### 01: Identify Rubisco and Prk

_Script:_ `source/01_Identify_Rubisco_and_Prk.sh`

_Output:_ `data/rubisco.txt`, `data/prk.txt`, `data/positive_genomes.txt`, and more

This script contains a long series of commands used to identify and filter Rubisco and Prk sequences in the GTDB genomes, followed by designating genomes as CBB-positive if they carry both Rubisco and Prk. Manual intervention is required for downloading sequences from [UniProt](https://www.uniprot.org/), submitting sequences to [KEGG BlastKOALA](https://www.kegg.jp/blastkoala/), and selecting the appropriate position of the K174 amino acid in the Rubisco alignment.

### 02: Create phylogenetic trees for Archaea and Bacteria

_Script:_ `source/02_Create_trees.sh`

_Output:_ `intermediate/archaea.tree` and `intermediate/bacteria.tree`

Creates phylogenetic trees for use in the analyses below. Requires the GTDB core protein alignments `ar122_msa.faa` and `bac120_msa.faa` saved as `data/archaea_msa.fasta` and `data/bacteria_msa.fasta`.

### 03: Create distance matrices for Archaea and Bacteria

_Script:_ `source/03_Create_distance_matrices.sh`

_Output:_ `intermediate/archaea.dist` and `intermediate/bacteria.dist`

Creates phylogenetic distance matrices for use in the analyses below. Requires the GTDB core protein alignments `ar122_msa.faa` and `bac120_msa.faa` saved as `data/archaea_msa.fasta` and `data/bacteria_msa.fasta` (not included).

### 04: Select closely related CBB-negative genomes

_Script:_ `source/04_Select_closest_relatives_as_negative.R`

_Output:_ `intermediate/example_genomes.tab`

Selects as many CBB-negative genomes as there are CBB-positive genomes, maintaining a phylogenetic distance that is as short as possible between the two groups, based on the distance matrices generated in the preceding step.

### 05: Identify genetic features in genomes (EC and Pfam)

_Script:_ `source/05_Identify_functions_in_example_genomes.sh`

_Output:_ `intermediate/deepec.tab.gz` and `intermediate/pfam.tab.gz`

Creates a one-column list of all example genomes (CBB-positive and CBB-negative), _i.e._ `intermediate/example_genomes.txt`, and then runs DeepEC and hmmsearch with Pfam on the example genome ORFs to acquire EC and Pfam annotations. Requires the Pfam-A database with HMMs saved as `data/Pfam-A.hmm.gz` (not included).

### 06: Prepare training data for random forest analysis

_Script:_ `source/06_Format_random_forest_input.R`

_Output:_ Files containing EC and Pfam count matrices, CBB-status, and names.

Filters out all Rubisco and Prk EC and Pfam annotations and formats a training dataset for the random forest analysis.

### 07: Random forest analysis

_Script:_ `source/07_Perform_random_forest_analysis.sh`

_Output:_ Log files with accuracies, importance and prediction reports.

Runs the random forest analysis on EC and Pfam features using `source/Rank_features_with_random_forest.py`, yielding importance values and predictions for every genome from 100 random forests. Train:test split is 3:1 and logistic regression is used to select the 600 most promising EC and Pfam features.

### 08: Plot random forest results

_Script:_ `source/08_Create_feature_heatmap.R`

_Output:_ `results/Feature_heatmap.png`

Creates a heatmap with genomes sorted by average random forest prediction accuracy (columns), features sorted by average random forest importance (rows), and feature count represented by fill color.

### 09: Enrichment analysis

_Script:_ `source/09_Feature_enrichment.R`

_Output:_ `intermediate/feature_enrichment.tab`

Compare EC and Pfam counts between CBB-positive and CBB-negative genomes using a Wilcoxon rank sum test.

### 10: Plot taxonomic taxonomic distribution

_Script:_ `source/10_Taxonomic_distribution.R`

_Output:_ `results/taxonomic_distribution.pdf`

Generates a set of colors for taxonomic groups (`intermediate/accession_organism_colours.tab`) and plots the distribution of CBB-positive and CBB-negative genomes.

### 11: ACE analysis in Archaea

_Script:_ `source/11_Correlate_feature_histories_in_Archaea.R`

_Output:_ Feature history correlation results in Archaea and trees.

Performs ancestral character estimation in Archaea, and correlates the emergence of the Calvin cycle to other genetic features. Plots the highest absolute correlations on phylogenetic trees.

### 12: ACE analysis in Bacteria

_Script:_ `source/12_Correlate_feature_histories_in_Bacteria.R`

_Output:_ Subtree table, feature history correlation results in Bacteria, and trees.

Selects bacterial subtrees and plots them on the bacterial phylogenetic tree. Then performs ancestral character estimation in subtrees of Bacteria, and correlates the emergence of the Calvin cycle to other genetic features. Plots the highest absolute correlations on phylogenetic trees.

### 13: Compare methods and create supplementary tables

_Script:_ `source/13_Create_supplementary_and_compare.R`

_Output:_ Tables for enrichment, ACE, and random forest analyses, and comparisons.

Creates supplementary tables for the enrichment, ACE, and random forest analyses, with feature ranks. Correlates ranks between methods. Compares methods with a plot of the ranks.

### 14: Create consensus ranking

_Script:_ `source/14_Create_consensus_ranking.R`

_Output:_ `results/Supplementary_Consensus_rank.tab` and top 20 tables.

Ranks the sum of ranks for all methods to generate a consensus rank for each genetic feature (EC or Pfam). Highlights the top 20 EC and top 20 Pfam features for Table 1 (see step 31).

### 15: Random forest accuracy and genome completeness

_Script:_ `source/15_Compare_accuracy_to_completeness.R`

_Output:_ `results/EC_count_features.acc_vs_comp.png`

Compares average per-genome random forest accuracy and genome completeness in a plot.

### 16: Proximity of genetic features to Rubisco and Prk

_Script:_ `source/16_Gene_proximity.R`

_Output:_ `results/Supplementary_Gene_proximity.tab`

Count the distance in number of ORFs between all genetic features and Rubisco and Prk in CBB-positive genomes.

### 17: Create method comparison plot

_Script:_ `source/17_Comparison_figure.R`

_Output:_ `results/method_feature_rank_comparison_2.pdf`

Plots consensus rank versus method ranks for every genetic feature, with color based on proximity to Rubisco in CBB-positive genomes.

### 18: Create example genomes supplementary dataset

_Script:_ `source/18_Create_Genomes_dataset.R`

_Output:_ `results/Supplementary_Example_genomes.tab.gz`

Creates a supplementary dataset with all feature counts (except for ECs and Pfams representing Rubisco and Prk) as well as information about each example genome; closest selected relative, distance to closest selected relative, subtree affiliation, CBB status, genome completeness, and GTDB taxonomy.

### 19: Random forest importance and abundance difference

_Script:_ `source/19_Correlate_importance_and_difference.R`

_Output:_ A Spearman correlation _r_ value.

Correlate difference in mean count between CBB-positive and CBB-negative genomes for each feature and the random forest importance values of those features.

### 20: Create ORF annotations supplementary dataset

_Script:_ `source/20_Create_ORF_annotations_dataset.sh`

_Output:_ `results/Supplementary_ORF_annotations.tab.gz`

Creates a supplementary dataset with all ORF annotations (EC and Pfam). ORF names can be used to track position of ORFs in official NCBI contigs.

### 21: Investigate completeness of cyanobacterial genomes

_Script:_ `source/21_Cyano_completeness.R`

_Output:_ Test and summary results for cyanobacterial completeness versus CBB status.

Performs a Wilcoxon rank sum test to check whether Cyanobacteria with the CBB-positive and CBB-negative classifications have significantly different genome completeness. Also compares CBB status of highly complete and less complete genomes.

### 22: Generate colors for EC enrichment to use with KEGG maps

_Script:_ `source/22_EC_count_colours_for_KEGG.R`

_Output:_ `results/EC_count_features.KEGG_difference_colours.discrete.txt`

Generate colors based on enrichment/depletion of ECs to use with [KEGG's pathway mapping tool](https://www.genome.jp/kegg/tool/map_pathway2.html).

### 23: Generate colors for EC consensus rank to use with KEGG maps

_Script:_ `source/23_Consensus_rank_colours_for_KEGG.R`

_Output:_ `results/EC_count_features.KEGG_consensus_rank_colours.txt`

Generate colors based on EC consensus ranks to use with [KEGG's pathway mapping tool](https://www.genome.jp/kegg/tool/map_pathway2.html).

### 24: Identify transport features

_Script:_ `source/24_Find_transport_features.R`

_Output:_ Number of transport features that were significant for each method.

Searches for transport-related terms in feature descriptions to count total number of transport features and how many were significant or ranked in the top 1000 (consensus rank).

### 25: Identify photosynthetic genomes

_Script:_ `source/25_Photosynthetic_genomes.R`

_Output:_ `results/Supplementary_Photosynthesis.tab` and `results/genome_photosynthesis_status.tab`

Identify photosynthesis-related features and determine how many CBB-positive and CBB-negative genomes may be photosynthetic.

### 26: Generate a Rubisco phylogenetic tree

_Script:_ `source/26_Make_Rubisco_tree.sh`

_Output:_ `intermediate/rubisco_with_Tabita_examples.tree`

Aligns the identified Rubisco sequences with examples from [Tabita _et al._ (2007)](https://doi.org/10.1128%2FMMBR.00015-07) and constructs a phylogenetic tree to identify Rubisco forms and ensure exclusion of Rubisco-like proteins.

### 27: Plot the Rubisco phylogenetic tree

_Script:_ `source/27_Plot_Rubisco_tree.R`

_Output:_ `results/rubisco_with_Tabita_examples.pdf`

Plots the tree constructed in step 26.

### 28: Count domains of unknown function

_Script:_ `source/28_Count_DUFs.R`

_Output:_ Counts of domains of unknown function among top ranks of the methods.

Counts the domains of unknown function (DUFs/UPFs) in the top 100 of enrichment and ACE analyses, and among all selected 600 Pfams in the random forest analysis.

### 29: Check the ranges of significant enrichment/depletion

_Script:_ `source/29_Range_of_fold_enrichment.R`

_Output:_ Maximum and minimum log<sub>2</sub> fold enrichment.

Calculate the maximum and minimum log<sub>2</sub> fold enrichment of all significant features and also identify the smallest absolute ratio that is still significant.

### 30: Count the number of significant features in ACE in each subtree

_Script:_ `source/30_Count_significant_features_in_subtrees.R`

_Output:_ The number of significant features in ACE in each subtree.

Counts the number of significant features in the ACE analysis for each subtree.

### 31: Create Table 1

_Script:_ `source/31_Create_Table_1.R`

_Output:_ `results/Table_1_Top_consensus_ranks.tab`

Tidy up the top 20 EC and Pfam features output from step 14 to create Table 1.

### 32: List Cyanobacteria that are CBB-positive

_Script:_ `source/32_List_cyanobacteria_with_CBB.R`

_Output:_ `intermediate/cyano_with_CBB.txt`

Make a list of Cyanobacteria that were identified to carry the Calvin cycle.

### 33: Select closest CBB-negative relatives of Cyanobacteria

_Script:_ `source/33_Select_cyano_closest_relatives.R`

_Output:_ `intermediate/cyano_relatives.tab`

Select closest CBB-negative relatives of CBB-positive Cyanobacteria as in step 04.

### 34: Compare distance to relatives in Cyanobacteria and others

_Script:_ `source/34_Cyano_vs_others_relative_dist.R`

_Output:_ `results/cyano_vs_others_relative_dist.tab`

Compare the distance between CBB-positive Cyanobacteria and their CBB-negative relatives to the distance between other CBB-positive and CBB-negative microbes with a Wilcoxon rank sum test.

### Extra: Plotting of ACE results on phylogenetic trees

_Script:_ `source/Plot_ACE.R`

_Output:_ A PDF with a phylogenetic (sub)tree displaying ACE data.

Plots ACE results of any feature on any subtree (as long as the feature is present in the subtree organisms) to allow creating phylogenetic tree figures displaying correlations between features and the CBB status.

## Author
Johannes Asplund-Samuelsson (johannes.aspsam@gmail.com)
