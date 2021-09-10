## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

options(knitr.duplicate.label = "allow")

#Test if BioManager is installed 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
   BiocManager::install(version = "3.12")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(optparse)

p_load(seqinr)
p_load(ProteomeRiver)
p_load(janitor)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
expt_group <- "RPE" 

## Directories management 
base_dir <- here::here()
data_dir <- file.path( base_dir, "Data", "ALPK1")
results_dir <- file.path(base_dir, "Results", "ALPK1", expt_group,  "Proteins")
source_dir <- file.path(base_dir, "Source")

# Minimum number of samples per experimental group for the protein to be accepted for analysis 

fasta_file <- file.path( data_dir,   "ALPK1-set1MusMusculus20201226CanIso.fasta")
raw_counts_file <- file.path(data_dir, "ALPK1-set1proteinGroups.txt")
column_pattern_input <- "Reporter intensity corrected"
output_counts_file <- "counts_table_cleaned.tab"  


razor_unique_peptides_group_thresh <- 0
unique_peptides_group_thresh <- 1

command_line_options <- commandArgs(trailingOnly = TRUE)
if( length(command_line_options ) > 0 ) {
  parser <- OptionParser(add_help_option =TRUE)
  
  parser <- add_option(parser, c( "--output-dir"), type="character", default="", dest = "results_dir",
                       help="Directory path for all results files.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--fasta"), type="character", default="", dest = "fasta_file",
                       help="Input protein sequence FASTA file with UniProt FASTA header format.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--raw-counts"), type="character", default="", dest = "raw_counts_file",
                       help="Input file with the protein abundance data.",
                       metavar="string")   
  
  parser <- add_option(parser, c( "--output-counts"), type="character", default="", dest = "output_counts_file",
                       help="String representing the name of the output counts data file which will be saved in the directory specified with the --output-dir flag.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--column-pattern"), type="character", default="Reporter intensity corrected", dest = "column_pattern_input",
                       help="String pattern, together with the experimental group pattern, that matches the abundance data columns.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--r-u-count"), type="integer", default=0, dest = "razor_unique_peptides_group_thresh",
                       help="Number of razor and unique peptides for the specified experiemtal group needs to be higher than this threshold for the protein to be included for the analysis.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--u-count"), type="integer", default=1, dest = "unique_peptides_group_thresh",
                       help="Number of unique peptides for the specified experiemtal group needs to be higher than this threshold for the protein to be included for the analysis.",
                       metavar="string")   

  print(command_line_options)
  
  cmd_arguments <- parse_args(parser)
  
  print(cmd_arguments)  

  results_dir <- cmd_arguments$results_dir
  fasta_file <- cmd_arguments$fasta_file
  raw_counts_file <- cmd_arguments$raw_counts_file
  output_counts_file <- cmd_arguments$output_counts_file
  column_pattern_input <- cmd_arguments$column_pattern_input
  razor_unique_peptides_group_thresh <- cmd_arguments$razor_unique_peptides_group_thresh
  unique_peptides_group_thresh <- cmd_arguments$unique_peptides_group_thresh

  
  if( results_dir == "" ) { stop("No value for --output-dir was supplied.") }
  if( fasta_file == "" ) { stop("No value for --fasta was supplied.") }
  if( raw_counts_file == "" ) { stop("No value for -raw-counts was supplied.") }
  if( column_pattern_input == "" ) { stop("No value for --column-pattern was supplied.") }
  if( razor_unique_peptides_group_thresh == "" ) { stop("No value for --r-u-count was supplied.") }
  if( unique_peptides_group_thresh == "" ) { stop("No value for --u-count was supplied.") }

}




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  column_pattern <- paste0(make_clean_names(column_pattern_input),  paste("_\\d+", tolower(expt_group), sep="_")) 

  extract_replicate_group <- paste0(make_clean_names(column_pattern_input),  paste0("_(\\d+)_(", tolower(expt_group), ")")) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
create_dir_if_not_exists( results_dir)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dat_tbl <- vroom::vroom( raw_counts_file )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_cln <-  janitor::clean_names(dat_tbl)

colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids" )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read fasta file 
# colnames(aa_seq_tbl)
print( "Step 1: Reading the fasta file.")

if(!file.exists( file.path( results_dir, "aa_seq_tbl.RDS"))) {
  aa_seq_tbl <- parse_fasta_file( fasta_file)
  saveRDS( aa_seq_tbl,  file.path( results_dir, "aa_seq_tbl.RDS"))
} else {
  aa_seq_tbl <- readRDS( file.path( results_dir, "aa_seq_tbl.RDS"))
}

## Add the row id column and create a column containing the cleaned  peptide
print("Step 2: Get row ID and get cleaned peptide sequence.")
evidence_tbl<- dat_cln %>% 
  mutate( evidence_id = row_number()) 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

razor_unique_peptides_group_col <-  paste0( "razor_unique_peptides_", tolower(expt_group)  ) 
unique_peptides_group_col <-  paste0("unique_peptides_", tolower(expt_group)  )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

peptides_count_helper <- evidence_tbl %>%
  dplyr::select( evidence_id, 
                 protein_ids, 
                 !!rlang::sym(razor_unique_peptides_group_col),
                 !!rlang::sym(unique_peptides_group_col), 
                 reverse,
                 potential_contaminant, 
                 matches(column_pattern)) %>%
  dplyr::filter(  is.na(reverse)	&  
                  is.na(potential_contaminant)) %>%
  dplyr::filter( !str_detect( protein_ids, "CON__") &
                 !str_detect( protein_ids, "REV__")   ) %>%
  dplyr::mutate( protein_ids = str_split( protein_ids, ";" ) ) %>%
  dplyr::mutate( !!rlang::sym(razor_unique_peptides_group_col) := str_split( !!rlang::sym(razor_unique_peptides_group_col), ";" ) )%>%
  dplyr::mutate( !!rlang::sym(unique_peptides_group_col) := str_split( !!rlang::sym(unique_peptides_group_col) , ";" ) ) %>%
  unnest(cols = c( protein_ids, 
                   !!rlang::sym(razor_unique_peptides_group_col), 
                   !!rlang::sym(unique_peptides_group_col)))

evidence_tbl_cleaned <- NA

evidence_tbl_cleaned <- peptides_count_helper %>%
    dplyr::filter( !!rlang::sym(razor_unique_peptides_group_col)  >= razor_unique_peptides_group_thresh &
                   !!rlang::sym(unique_peptides_group_col)  >= unique_peptides_group_thresh  ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  accession_gene_name_tbl <- choose_best_protein_accession( input_tbl = evidence_tbl_cleaned, 
                                                    acc_detail_tab = aa_seq_tbl, 
                                                    accessions_column = protein_ids , 
                                                    row_id_column= uniprot_acc,
                                                    group_id = evidence_id)

evidence_tbl_filt <- evidence_tbl_cleaned %>%
  inner_join( accession_gene_name_tbl %>% 
                dplyr::select(evidence_id, uniprot_acc), by="evidence_id") %>%
  dplyr::select(uniprot_acc, contains(expt_group), -contains(c("razor", "unique"))) %>%
  distinct

colnames( evidence_tbl_filt) <- str_replace_all(colnames(evidence_tbl_filt), extract_replicate_group, "\\1_\\U\\2")  %>% 
  toupper( ) %>%
  str_replace_all( "UNIPROT_ACC", "uniprot_acc") 
  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vroom::vroom_write( evidence_tbl_filt, file.path(results_dir, output_counts_file) ) 

