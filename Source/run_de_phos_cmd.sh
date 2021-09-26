#!/bin/bash

# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

BASE_DIR=/home/ignatius/PostDoc/2021/ALPK1_2021
TEMPLATE_DIR=$BASE_DIR/Source/Shared
DATA_DIR=$BASE_DIR/Data/Abundance_Data/P90
PROTEINS_DIR=$BASE_DIR/Results/RPE90/Proteins
RESULTS_DIR="$BASE_DIR/Results/RPE90/Phosphopeptides"
PHOS_DIR=$RESULTS_DIR/Abundance_Tables
PHOS_RESULTS_DIR=$RESULTS_DIR/DE_Analysis

## Clean phosphorylation data  
Rscript --vanilla $TEMPLATE_DIR/clean_phos_cmd.R \
 --output-dir=$PHOS_DIR \
 --fasta=$DATA_DIR/ALPK1-set1MusMusculus20201226CanIso.fasta \
 --fasta-save=$PHOS_DIR/aa_seq_tbl.RDS  \
 --raw-counts=$DATA_DIR/ALPK1-set1evidence.txt \
 --site-prob=0.75 \
 --recover-prob=0.5 \
 --column-pattern="Reporter intensity corrected" \
 --add-columns="experiment"
 
Rscript --vanilla $TEMPLATE_DIR/de_analysis_cmd.R \
 --max-missing=0 \
 --abundance-thresh=0 \
 --q-value-thresh=0.05 \
 --group-pattern="RPE" \
 --column-pattern="\\d+_RPE" \
 --ruv-k=5  \
 --num-neg-ctrl=500 \
 --ruv-method="ruv3" \
 --counts=$PHOS_DIR/sum_phoshpsites.tsv \
 --contrasts=$DATA_DIR/contrast_strings_RPE.tab \
 --design-matrix=$PROTEINS_DIR/DE_Analysis/design_matrix_cleaned.tab \
 --formula="~ 0 + group " \
 --output-dir=$PHOS_RESULTS_DIR  \
 --sample-id=Sample_ID \
 --group-id=group \
 --row-id=sites_id \
 --prefix=de_phos
 
Rscript --vanilla $TEMPLATE_DIR/annot_phos_cmd.R \
 --tax-id 10090 \
 --output-dir=$PHOS_RESULTS_DIR \
 --raw-counts=$PHOS_DIR/sum_phoshpsites.tsv \
 --input-wide=$PHOS_RESULTS_DIR/de_phos_wide.tsv \
 --input-long=$PHOS_RESULTS_DIR/de_phos_long.tsv \
 --output-wide=$PHOS_RESULTS_DIR/de_phos_wide_annot.tsv \
 --output-long=$PHOS_RESULTS_DIR/de_phos_long_annot.tsv \
 --psp-dir="$BASE_DIR/Data/PhosphositePlus/20210730" \
 --near-ptm=5
 
Rscript --vanilla $TEMPLATE_DIR/norm_phos_by_prot_abundance.R \
 --output-dir=$PHOS_RESULTS_DIR  \
 --protein=$PROTEINS_DIR/DE_Analysis/de_proteins_long_annot.tsv \
 --phospho=$PHOS_RESULTS_DIR/de_phos_long_annot.tsv 
 