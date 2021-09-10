#!/bin/bash

BASE_DIR=/home/ignatius/PostDoc/2021/MultiPhos2021
DATA_DIR=$BASE_DIR/Data/ALPK1
RESULTS_DIR="$BASE_DIR/Results/ALPK1/RPE/Proteins"

Rscript --vanilla clean_proteins_cmd.R \
 --output-dir=$RESULTS_DIR \
 --fasta=$DATA_DIR/ALPK1-set1MusMusculus20201226CanIso.fasta \
 --raw-counts=$DATA_DIR/ALPK1-set1proteinGroups.txt \
 --output-counts=counts_table_cleaned.tab \
 --column-pattern="Reporter intensity corrected" \
 --r-u-count=0 \
 --u-count=1 

Rscript --vanilla de_proteins_cmd.R \
 --limma-method="contrasts" \
 --min-samples=4 \
 --abundance-thresh=0 \
 --q-value-thresh=0.05 \
 --group-pattern="RPE" \
 --column-pattern="Reporter intensity corrected \\d+ RPE" \
 --ruv-k=4  \
 --num-neg-ctrl=500 \
 --ruv-method="ruv3" \
 --counts=$RESULTS_DIR/counts_table_cleaned.tab \
 --contrasts=$DATA_DIR/contrast_strings.tab \
 --design-matrix=$RESULTS_DIR/design_matrix_cleaned.tab \
 --output-dir=$RESULTS_DIR  \
 --sample-id=Sample_ID \
 --group-id=group \
 --row-id=uniprot_acc
 
 

Rscript --vanilla annot_proteins_cmd.R \
 --output-dir=$RESULTS_DIR \
 --input-wide=$RESULTS_DIR/de_proteins_longer.tsv \
 --input-long=$RESULTS_DIR/de_proteins_wider.tsv \
 --output-wide=$RESULTS_DIR/de_proteins_wider_annot.tsv \
 --output-long=$RESULTS_DIR/de_proteins_longer_annot.tsv
 
 
 
 
 