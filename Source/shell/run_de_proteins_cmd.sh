#!/bin/bin

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

E_BADARGS=65

#Help string
basename="$(basename $0)"
help_string="Usage: $basename --template_dir TEMPLATE_DIR --data_dir DATA_DIR --results_dir RESULTS_DIR \n
Runs the proteome river pipeline. \n
\n
  -h,  --help\t\t\t Show help options \n
  -t,  --template_dir TEMPLATE_DIR Path to ProteomeRiver's scripts [Default: /opt/ProteomeRiver/Source] \n
  -d,  --data_dir DATA_DIR\t Path to input data (required) \n
  -r,  --results_dir\t\t Path to results [Default: $PWD]\n
 \n"

#Default values:
TEMPLATE_DIR=$BASE_DIR/Source/Shared
RESULTS_DIR=$PWD

#Argument parsing
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-h | --help)
	    HELP=TRUE
	    shift # past argument
	    ;;
	-t | --template_dir)
	    TEMPLATE_DIR="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-d | --data_dir)
	    DATA_DIR="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-r | --results_dir)
	    RESULTS_DIR="$2"
	    shift # past argument
	    shift # past value
	    ;;
	*)    # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
    esac
done

#test help flag
if [  $HELP ];
then
    echo -e "$help_string"
  exit 0
fi

#test mandatory arguments
if [  -z "$DATA_DIR" ];
then
    echo -e "Error: -d,--data_dir is required.\n\n"
    echo -e  $help_string
  exit $E_BADARGS
fi


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
