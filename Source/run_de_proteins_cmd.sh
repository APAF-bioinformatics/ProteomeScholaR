#!/bin/bash

# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

BASE_DIR=/home/ignatius/PostDoc/2021/ALPK1_2021
TEMPLATE_DIR=$BASE_DIR/Source/Shared
DATA_DIR=$BASE_DIR/Data/Abundance_Data/P90
RESULTS_DIR="$BASE_DIR/Results/RPE90/Proteins"

Rscript --vanilla $TEMPLATE_DIR/clean_proteins_cmd.R \
 --output-dir=$RESULTS_DIR/DE_Analysis \
 --fasta=$DATA_DIR/ALPK1-set1MusMusculus20201226CanIso.fasta \
 --raw-counts=$DATA_DIR/ALPK1-set1proteinGroups.txt \
 --output-counts=raw_counts_table.tab \
 --ids=cleaned_accession_to_protein_group.tab \
 --column-pattern="Reporter intensity corrected" \
 --group-pattern="RPE" \
 --r-u-count=0 \
 --u-count=1 


# If there are other variables to adjust the linear model (e.g. gender and age), 
# include these as input in the design matrix and formula string.
Rscript --vanilla $TEMPLATE_DIR/de_analysis_cmd.R \
 --max-missing=0 \
 --abundance-thresh=0 \
 --q-value-thresh=0.05 \
 --group-pattern="RPE" \
 --column-pattern="Reporter intensity corrected \\d+ RPE" \
 --ruv-k=4  \
 --num-neg-ctrl=500 \
 --ruv-method="ruv3" \
 --counts=$RESULTS_DIR/DE_Analysis/raw_counts_table.tab \
 --contrasts=$DATA_DIR/contrast_strings_RPE.tab \
 --formula="~ 0 + group " \
 --design-matrix=$RESULTS_DIR/DE_Analysis/design_matrix_cleaned.tab \
 --output-dir=$RESULTS_DIR/DE_Analysis  \
 --sample-id=Sample_ID \
 --group-id=group \
 --row-id=uniprot_acc
 

Rscript --vanilla $TEMPLATE_DIR/annot_proteins_cmd.R \
 --tax-id=10090 \
 --output-dir=$RESULTS_DIR/DE_Analysis \
 --input-wide=$RESULTS_DIR/DE_Analysis/de_proteins_wide.tsv \
 --input-long=$RESULTS_DIR/DE_Analysis/de_proteins_long.tsv \
 --output-wide=$RESULTS_DIR/DE_Analysis/de_proteins_wide_annot.tsv \
 --output-long=$RESULTS_DIR/DE_Analysis/de_proteins_long_annot.tsv \
 --ids=$RESULTS_DIR/DE_Analysis/cleaned_accession_to_protein_group.tab \
 --raw-counts=$DATA_DIR/ALPK1-set1proteinGroups.txt 
