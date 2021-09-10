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
p_load(optparse)
p_load(tictoc)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tic()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(rlang)
p_load(UniProt.ws)
p_load(biomaRt)
p_load(GO.db)
p_load(ProteomeRiver)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list_of_sp_columns <- c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH", "ENSEMBL", "GO-ID", "KEYWORDS")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

group_pattern <- "RPE"
## Directories management 
base_dir <- here::here()
data_dir <- file.path( base_dir, "Data", "ALPK1")
results_dir <- file.path(base_dir, "Results", "ALPK1", group_pattern,  "Proteins")
source_dir <- file.path(base_dir, "Source")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

input_wide_file <- file.path(results_dir, "de_proteins_wider.tsv" ) 
input_long_file <- file.path(results_dir, "de_proteins_longer.tsv" ) 
output_long_file <- file.path(results_dir, "de_proteins_longer_annot.tsv" ) 
output_wide_file <- file.path(results_dir, "de_proteins_wider_annot.tsv" ) 

command_line_options <- commandArgs(trailingOnly = TRUE)
if( length(command_line_options ) > 0 ) {
  parser <- OptionParser(add_help_option =TRUE)

  
  parser <- add_option(parser, c( "--output-dir"), type="character", default="", dest = "results_dir",
                       help="Directory path for all results files.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--input-wide"), type="character", default="", dest = "input_wide_file",
                       help="Results table with values in wider format.",
                       metavar="string")   
  
  parser <- add_option(parser, c("--input-long"), type="character", default="", dest = "input_long_file",
                       help="Results table with values in longer format.",
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
  
  if( results_dir == "" ) { stop("No value for --output-dir was supplied.") }
  if( group_pattern == "" ) { stop("No value for --group-pattern was supplied.") }
  if( input_wide_file == "" ) { stop("No value for --input-wide was supplied. ") }
  if( input_long_file == "" ) { stop("No value for --input-long was supplied. ") }
  if( output_wide_file == "" ) { stop("No value for --output-wide was supplied. ") }
  if( output_long_file == "" ) { stop("No value for --output-long was supplied. ") }
  
}

de_proteins_wider <- vroom::vroom( input_wide_file ) 
de_proteins_longer <- vroom::vroom( input_long_file ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Download Reactome UniProt to pathways file.")
temp_file <- tempfile()

reactome_file <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=temp_file)

reactome_map <- vroom::vroom( temp_file , 
                              col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Get the best UniProt accession per row.")
  uniprot_acc_tbl <- de_proteins_wider %>%
    mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(  ., ":" )[[1]][1] )) %>%
    dplyr::select(best_uniprot_acc) %>%
    distinct(best_uniprot_acc)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Download information from UniProt.")
uniprot_dat <- NA
if( ! file.exists( file.path(results_dir, "uniprot_data.RDS"))) {
  
  if(  !"up" %in% ls()) {
    up <- UniProt.ws(taxId=10090 )
  } 

 # keytypes(up)
  
  uniprot_dat <- batch_query_evidence(uniprot_acc_tbl, best_uniprot_acc, uniprot_handle=up, 
                                      uniprot_columns = list_of_sp_columns)
  
  saveRDS( uniprot_dat, file.path(results_dir, "uniprot_data.RDS"))
  
} else {
  uniprot_dat <- readRDS( file.path(results_dir, "uniprot_data.RDS"))
}



## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##  # uniprot_dat <- batch_query_evidence(uniprot_acc_tbl %>% head(100), best_uniprot_acc, uniprot_handle=up,
##  #                                      uniprot_columns = list_of_sp_columns)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("Merge with Gene Ontology terms.")
goterms <- Term(GOTERM)
gotypes <- Ontology(GOTERM)

uniprot_dat_cln <- uniprot_dat %>%
  mutate( go_mapping =  purrr::map( `GO-ID`, ~go_id_to_term(., sep="; ", goterms, gotypes) )) %>%
  unnest(go_mapping) %>%
  dplyr::select(-go_mapping)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Add reactome pathways annotation.")
reactome_term_tbl <- uniprot_acc_tbl %>%
  left_join( reactome_map, by=c("best_uniprot_acc" = "uniprot_acc") )   %>%
  group_by(best_uniprot_acc) %>%
  summarise( reactome_term = paste(reactome_term, collapse="; ") ) %>%
  ungroup()  %>%
  dplyr::filter(reactome_term != "NA" )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Output wider format results table with protein annotation.")
de_proteins_wider_annot <- de_proteins_wider %>% 
    mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(  ., ":" )[[1]][1] )) %>%
  left_join( uniprot_dat_cln, by = c("best_uniprot_acc" = "UNIPROTKB") ) %>%
  left_join( reactome_term_tbl, by = c("best_uniprot_acc" = "best_uniprot_acc"))

head( de_proteins_wider_annot ) 

vroom::vroom_write(de_proteins_wider_annot, path=output_wide_file ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Output longer format results table with protein annotation.")
de_proteins_longer_annot <- de_proteins_longer %>% 
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(  ., ":" )[[1]][1] )) %>%
  left_join( uniprot_dat_cln, by = c("best_uniprot_acc" = "UNIPROTKB") ) %>%
  left_join( reactome_term_tbl, by = c("best_uniprot_acc" = "best_uniprot_acc"))

head( de_proteins_longer_annot )

vroom::vroom_write(de_proteins_longer_annot, path=output_long_file ) 



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
toc()
sessionInfo()

