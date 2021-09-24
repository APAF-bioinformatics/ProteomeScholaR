## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#! /usr/bin/env Rscript

# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
p_load(optparse)
p_load(tictoc)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tic()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(rlang)
p_load( janitor)
p_load(UniProt.ws)
p_load(biomaRt)
p_load(GO.db)
p_load(ProteomeRiver)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list_of_sp_columns <- c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH", "ENSEMBL", "GO-ID", "KEYWORDS")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

group_pattern <- "NR"
## Directories management 
base_dir <- here::here()
data_dir <- file.path( base_dir, "Data", "Abundance_Data", "P90")
results_dir <- file.path(base_dir, "Results",  paste0(group_pattern, "90"),  "Proteins", "DE_Analysis")
source_dir <- file.path(base_dir, "Source")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

input_wide_file <- file.path(results_dir, "de_proteins_wide.tsv" ) 
input_long_file <- file.path(results_dir, "de_proteins_long.tsv" ) 
output_long_file <- file.path(results_dir, "de_proteins_longer_annot.tsv" ) 
output_wide_file <- file.path(results_dir, "de_proteins_wider_annot.tsv" ) 
ids_file <- file.path( results_dir, "cleaned_accession_to_protein_group.tab")
taxonomy_id <- 10090
raw_counts_file <- file.path(data_dir, "ALPK1-set1proteinGroups.txt")


command_line_options <- commandArgs(trailingOnly = TRUE)
if( length(command_line_options ) > 0 ) {
  parser <- OptionParser(add_help_option =TRUE)

  parser <- add_option(parser, c( "--tax-id"), type="integer", default="", dest = "taxonomy_id",
                       help="The NCBI taxonomy ID of the organism being investigated (e.g. M. musculus=10090, H. sapien=9606).",
                       metavar="integer")     
  
  parser <- add_option(parser, c( "--output-dir"), type="character", default="", dest = "results_dir",
                       help="Directory path for all results files.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--input-wide"), type="character", default="", dest = "input_wide_file",
                       help="Results table with values in wider format.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--input-long"), type="character", default="", dest = "input_long_file",
                       help="Results table with values in longer format.",
                       metavar="string") 
  
  parser <- add_option(parser, c("--ids"), type="character", default="", dest = "ids_file",
                       help="File to link the cleaned list of accessions to the original list of protein groups from MaxQuant file. Also, contains the MaxQuant output row ID.",
                       metavar="string") 
  
  parser <- add_option(parser, c("--raw-counts"), type="character", default="", dest = "raw_counts_file",
                       help="Input file with the protein abundance data.",
                       metavar="string")   
  
  parser <- add_option(parser, c( "--output-wide"), type="character", default="", dest = "output_wide_file",
                       help="Results table with values in wider format.",
                       metavar="string")   
  
  parser <- add_option(parser, c( "--output-long"), type="character", default="", dest = "output_long_file",
                       help="Results table with values in longer format.",
                       metavar="string")   


  print(command_line_options)
  
  cmd_arguments <- parse_args(parser)
  
  print(cmd_arguments)  

  input_wide_file <- cmd_arguments$input_wide_file
  input_long_file <- cmd_arguments$input_long_file
  output_wide_file <- cmd_arguments$output_wide_file
  output_long_file <- cmd_arguments$output_long_file
  results_dir <- cmd_arguments$results_dir
  taxonomy_id <- cmd_arguments$taxonomy_id
  ids_file  <- cmd_arguments$ids_file
  raw_counts_file <- cmd_arguments$raw_counts_file

  if( results_dir == "" ) { stop("No value for --output-dir was supplied.") }
  if( group_pattern == "" ) { stop("No value for --group-pattern was supplied.") }
  if( input_wide_file == "" ) { stop("No value for --input-wide was supplied. ") }
  if( input_long_file == "" ) { stop("No value for --input-long was supplied. ") }
  if( output_wide_file == "" ) { stop("No value for --output-wide was supplied. ") }
  if( output_long_file == "" ) { stop("No value for --output-long was supplied. ") }
  if( taxonomy_id == "" ) { stop("No value for --tax_id was supplied. ") }
  if( ids_file == "" ) { stop("No value for --ids was supplied.") }
  if( raw_counts_file == "" ) { stop("No value for -raw-counts was supplied.") }
  
}

create_dir_if_not_exists(results_dir)
de_proteins_wider <- vroom::vroom( input_wide_file ) 
de_proteins_longer <- vroom::vroom( input_long_file ) 
ids_tbl <- vroom::vroom( ids_file ) 
dat_tbl <- vroom::vroom( raw_counts_file )
dat_cln <-  janitor::clean_names(dat_tbl) %>%
  dplyr::rename( gene_names_maxquant = "gene_names")
colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids" )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Download Reactome UniProt to pathways file.")

temp_file <- tempfile()

reactome_file <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=temp_file)

reactome_map <- vroom::vroom( temp_file , 
                              col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )






## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Get the best UniProt accession per row.")

uniprot_acc_tbl <- de_proteins_wider %>%
  mutate( uniprot_acc_copy = uniprot_acc ) %>%
  separate_rows(uniprot_acc_copy, sep=":" ) %>%
  mutate( join_uniprot_acc = clean_isoform_number(uniprot_acc_copy)) %>%
  dplyr::distinct( uniprot_acc, join_uniprot_acc) %>%
  group_by( uniprot_acc) %>%
  mutate( acc_order_id = row_number()) %>% 
  ungroup




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Download information from UniProt.")
uniprot_dat <- NA
up <- NA
if( ! file.exists( file.path(results_dir, "uniprot_data.RDS"))) {

  if( is.na(up)) {
    up <- UniProt.ws(taxId=taxonomy_id )
  } 

 # keytypes(up)
  
  uniprot_dat <- batch_query_evidence(uniprot_acc_tbl, join_uniprot_acc, uniprot_handle=up, 
                                      uniprot_columns = list_of_sp_columns)
  
  saveRDS( uniprot_dat, file.path(results_dir, "uniprot_data.RDS"))
  
} else {
  uniprot_dat <- readRDS( file.path(results_dir, "uniprot_data.RDS"))
}



## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##  # uniprot_dat <- batch_query_evidence(uniprot_acc_tbl %>% head(100), best_uniprot_acc, uniprot_handle=up,
##  #                                      uniprot_columns = list_of_sp_columns)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Merge with Gene Ontology terms.")
goterms <- Term(GOTERM)
gotypes <- Ontology(GOTERM)


tic()
uniprot_dat_cln <- uniprot_go_id_to_term( uniprot_dat, sep="; ", goterms, gotypes  )
toc()



uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
  left_join( uniprot_dat_cln, by=c("join_uniprot_acc" = "UNIPROTKB") )   %>%
  arrange( uniprot_acc, acc_order_id) %>%
  group_by(uniprot_acc ) %>%
  summarise( across( .cols=setdiff( colnames( uniprot_dat_cln), "UNIPROTKB")   , ~paste(., collapse=":"))   ) %>%
  ungroup()   %>%
  dplyr::rename( UNIPROT_GENENAME = "GENENAME")


uniprot_dat_multiple_acc





## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Add reactome pathways annotation.")
reactome_term_tbl <- uniprot_acc_tbl %>%
  left_join( reactome_map, by=c("join_uniprot_acc" = "uniprot_acc") )   %>%
  dplyr::filter(reactome_term != "NA" ) %>%
  group_by(uniprot_acc, join_uniprot_acc) %>%
  summarise( reactome_term = paste(reactome_term, collapse="; ") ) %>%
  ungroup()     %>%
  mutate(reactome_term = str_replace_all( reactome_term , ":", "-")) %>%
  group_by(uniprot_acc ) %>%
  summarise( reactome_term = paste(reactome_term, collapse=":") ) %>%
  ungroup()   




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Output wider format results table with protein annotation.")
de_proteins_wider_annot <- de_proteins_wider %>%
  left_join( ids_tbl, by=c("uniprot_acc" = "uniprot_acc") ) %>% 
  left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
  left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
  left_join( dat_cln, by=c("maxquant_row_id" = "id",
                           "protein_ids" = "protein_ids"))

head( de_proteins_wider_annot ) 

vroom::vroom_write(de_proteins_wider_annot, path=output_wide_file ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Output longer format results table with protein annotation.")
de_proteins_longer_annot <- de_proteins_longer %>%
  left_join( ids_tbl, by=c("uniprot_acc" = "uniprot_acc") ) %>% 
  left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
  left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
  left_join( dat_cln, by=c("maxquant_row_id" = "id",
                           "protein_ids" = "protein_ids"))

head( de_proteins_longer_annot )

vroom::vroom_write(de_proteins_longer_annot, path=output_long_file ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
toc()
sessionInfo()

