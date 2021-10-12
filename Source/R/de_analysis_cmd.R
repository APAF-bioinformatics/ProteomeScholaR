#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.12")
}

# load pacman package manager
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(optparse)
p_load(tictoc)
p_load(tidyverse)
p_load(plotly)
p_load(vroom)
p_load(ggplot2)
p_load(ggpubr)
p_load(magrittr)
p_load(knitr)
p_load(rlang)


p_load(limma)
p_load(qvalue)
p_load(ruv)
p_load(mixOmics)

p_load(ProteomeRiver)
p_load(configr)
p_load(logging)
p_load(svglite)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic()


command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option = TRUE)

#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o", "--output_dir"), type = "character", default = "de_analysis", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, "--max_num_samples_miss_per_group", type = "integer", dest = "max_num_samples_miss_per_group",
                     help = "Remove protein if it exceeds this maximum number of samples with missing values per experimental group [default %default]",
                     metavar = "integer")

parser <- add_option(parser, "--abundance_threshold", type = "integer", dest = "abundance_threshold",
                     help = "Abundance threshold above which the protein in the sample is accepted for analysis [default %default]",
                     metavar = "integer")


parser <- add_option(parser, "--group_pattern", type = "character", dest = "group_pattern",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")

parser <- add_option(parser, "--q_val_thresh", type = "double", dest = "q_val_thresh",
                     help = "q-value threshold below which a protein has statistically significant differetial expression [default %default]",
                     metavar = "double")

parser <- add_option(parser, "--ruv_k", type = "integer", dest = "ruv_k",
                     help = "The number of unwanted factors to use.  [default %default]",
                     metavar = "integer")

parser <- add_option(parser, "--num_neg_ctrl", type = "integer", dest = "num_neg_ctrl",
                     help = "The number of negative control proteins to use.  [default %default]",
                     metavar = "number")

parser <- add_option(parser, "--ruv_method", type = "character", dest = "ruv_method",
                     help = "A string representing the ruv3 method to use.  [default %default]",
                     metavar = "character")

parser <- add_option(parser, "--counts_table_file", type = "character", dest = "counts_table_file",
                     help = "Input file with the protein abundance values",
                     metavar = "string")

parser <- add_option(parser, "--test_pairs_file", type = "character", dest = "test_pairs_file",
                     help = "Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.",
                     metavar = "string")

parser <- add_option(parser, "--contrasts_file", type = "character", dest = "contrasts_file",
                     help = "Input file with a table listing all comparisons to be made in string, one comparison per line (e.g. groupB.vs.group_A = groupB - groupA).",
                     metavar = "string")

parser <- add_option(parser, "--formula_string", type = "character", dest = "formula_string",
                     help = "A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group).",
                     metavar = "string")

parser <- add_option(parser, "--design_matrix_file", type = "character", dest = "design_matrix_file",
                     help = "Input file with the design matrix",
                     metavar = "string")

parser <- add_option(parser, "--sample_id", type = "character", dest = "sample_id",
                     help = "A string describing the sample ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--group_id", type = "character", dest = "group_id",
                     help = "A string describing the experimental group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--row_id", type = "character", dest = "row_id",
                     help = "A string describing the row id.",
                     metavar = "string")

parser <- add_option(parser, "--file_prefix", type = "character", dest = "file_prefix",
                     help = "A string to indicate the type of analysis and is used in the file name of the output results table.",
                     metavar = "string")


#parse comand line arguments first.
args <- parse_args(parser)

createOutputDir(args$output_dir, args$no_backup)

## Logger configuration
logReset()
level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
addHandler(writeToConsole, level = ifelse(args$silent, loglevels["ERROR"], level), formatter = cmriFormatter)
addHandler(writeToFile, file = file.path(args$output_dir, args$log_file), level = level, formatter = cmriFormatter)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "de_analysis"), args)
}

cmriWelcome("ProteomeRiver", c("Ignatius Pang", "Pablo Galaviz"))
loginfo("Reading configuration file %s", args$config)
loginfo("Argument: Value")
loginfo("----------------------------------------------------")
for (v in names(args))
{
  loginfo("%s : %s", v, args[v])
}
loginfo("----------------------------------------------------")

testRequiredFiles(c(
  args$counts_table_file
  , args$contrasts_file
  , args$design_matrix_file
))

testRequiredArguments(args, c(
  "group_pattern"
  , "formula_string"
  , "ruv_k"
  , "q_val_thresh"
  , "num_neg_ctrl"
  , "ruv_method"
  , "sample_id"
  , "group_id"
  , "row_id"
  , "file_prefix"
))


args<-parseType(args,
  c("q_val_thresh")
                ,as.double)

args<-parseType(args,
  c("ruv_k"
    ,"num_neg_ctrl"
    ,"max_num_samples_miss_per_group"
    ,"abundance_threshold"
  )
                ,as.integer)

if (args$group_pattern == "") {
  logwarn("Empty group pattern string, using \\d+")
  args$group_pattern <- "\\d+"
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# test_pairs <- NA
# if( limma_method == "pairs" & args$test_pairs_file != "") {
#   print("Read file with pairs of experimental groups to test.")
#   test_pairs <- vroom::vroom( args$test_pairs_file)
# }


loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
contrasts_tbl <- vroom::vroom(args$contrasts_file, delim = "\t")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read file with counts table %s", args$counts_table_file)
evidence_tbl_filt <- vroom::vroom(args$counts_table_file, delim = "\t") %>%
  dplyr::select(one_of(c(args$row_id)), matches(args$group_pattern))


if (!args$row_id %in% colnames(evidence_tbl_filt)) {
  logerror("The row ID specified in --row-id was not found in the input counts file.")
  q()
}
if (length(which(str_detect(setdiff(colnames(evidence_tbl_filt), args$row_id), args$group_pattern))) == 0) {
  logerror("The input counts file did not contain any columns that matches the string provided in --group-pattern=%s", args$group_pattern)
  q()
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Remove empty proteins without abundance data.")
cln_dat_wide_unsorted <- ProteomeRiver::removeEmptyRows(evidence_tbl_filt,
                                                          col_pattern = args$group_pattern,
                                                          row_id = !!rlang::sym(args$row_id))
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

design_mat_cln <- vroom::vroom(args$design_matrix_file) %>%
  as.data.frame() %>%
  dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))

rownames(design_mat_cln) <- design_mat_cln %>% pull(as.name(args$sample_id))

## Check that the sample ID and group ID name is found in the design matrix 
if (length(which(c(args$sample_id, args$group_id) %in% colnames(design_mat_cln))) != 2) {
  logerror("Sample ID and group ID are not matching to the column names used in the design matrix.")
  q()
}

cols_for_analysis <- design_mat_cln %>% pull(as.name(args$sample_id))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cln_dat_wide_cleaned <- NA
if ("max_num_samples_miss_per_group" %in% names(args)) {
  loginfo("Remove row with any missing values")

  cln_dat_wide_cleaned <- removeRowsWithMissingValues(cln_dat_wide_unsorted,
                                                      matches(args$group_pattern),
                                                      design_mat_cln %>%
                                                            dplyr::mutate(Sample_ID = as.character(Sample_ID)),
                                                      !!rlang::sym(args$sample_id),
                                                      !!rlang::sym(args$row_id),
                                                      !!rlang::sym(args$group_id),
                                                      args$max_num_samples_miss_per_group,
                                                      args$abundance_threshold)

} else {
  cln_dat_wide_cleaned <- cln_dat_wide_unsorted
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (length(cln_dat_wide_cleaned) - 1 != length(cols_for_analysis) |
  length(cols_for_analysis) != length(intersect(colnames(cln_dat_wide_cleaned)[-1], cols_for_analysis))) {
  logerror("Column names in design matrix table not consistent with column names in abundance table.")
  q()
}

# Sort the columns in the abundance table according
if (length(which(colnames(cln_dat_wide_cleaned)[-1] == cols_for_analysis)) != length(cols_for_analysis)) {
  logwarn("Sorting the order of the column of the abundance matrix according to the order in the design matrix.")
}

cln_dat_wide <- cln_dat_wide_cleaned[, c(colnames(cln_dat_wide_cleaned)[1], cols_for_analysis)]

loginfo("Assign the indexing row names to the data frame.")
counts_filt <- cln_dat_wide %>%
  column_to_rownames(args$row_id)

vroom::vroom_write(cln_dat_wide, file.path(args$output_dir, "raw_counts_after_removing_empty_rows.tsv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix(design_mat_cln,
                                                     !!rlang::sym(args$sample_id),
                                                     !!rlang::sym(args$group_id))

loginfo(ruvIII_replicates_matrix)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# sample_rows_lists <- get_rows_to_keep_list( cln_dat_wide, 
#                                             matches(args$group_pattern),
#                                             design_mat_cln %>% 
#                                               dplyr::mutate( Sample_ID = as.character(Sample_ID)), 
#                                             !!rlang::sym(args$sample_id),
#                                             !!rlang::sym(args$row_id),
#                                             !!rlang::sym(args$group_id),
#                                             args$min_num_samples_per_group,
#                                             args$abundance_threshold)
# 
# saveRDS(sample_rows_lists, file.path(args$output_dir, "keep_sample_rows_lists.RDS" ))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Count the total number of missing values in total
loginfo("Count the number of missing values for each sample:",table(is.infinite(data.matrix(log2(counts_filt)))))


plot_num_missing_values <- plotNumMissingVales(counts_filt[, cols_for_analysis])


ggsave(filename = file.path(args$output_dir, "num_missing_values.png"), plot = plot_num_missing_values,limitsize = FALSE)
ggsave(filename = file.path(args$output_dir, "num_missing_values.svg"), plot = plot_num_missing_values,limitsize = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Mark entries with zero values as NA, take log of values.")
na_values_marker <- (counts_filt == 0)

counts_na <- counts_filt

counts_na[na_values_marker] <- NA

counts_na.log <- log2(counts_na)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Median centering.")
counts_na.log.quant <- normalizeBetweenArrays(counts_na.log, method = "scale")


vroom::vroom_write(as.data.frame(counts_na.log.quant) %>%
                     rownames_to_column(args$row_id),
                   file.path(args$output_dir, "counts_after_median_scaling_before_imputation.tsv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Missing value imputation")

counts_rnorm.log.quant <- counts_na.log.quant %>%
  as.data.frame %>%
  mutate_all(~{ imputePerCol(., width = 0.3, downshift = 1.8) }) %>%
  as.matrix

counts_rnorm.is_nan <- counts_rnorm.log.quant %>%
  as.data.frame %>%
  mutate_all(~{ !is.finite(.) }) %>%
  as.matrix


vroom::vroom_write(as.data.frame(counts_rnorm.log.quant) %>%
                     rownames_to_column(args$row_id),
                   file.path(args$output_dir, "counts_after_median_scaling_and_imputation.tsv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Run statistical tests without RUV.")

ID <- rownames(counts_rnorm.log.quant)

type_of_grouping <- getTypeOfGrouping(design_matrix = design_mat_cln, group_id = args$group_id, sample_id = args$sample_id)

list_rnorm.log.quant.ruv.r0 <- NA
myRes_rnorm.log.quant <- NA

# if( limma_method == "pairs") {
#   
#   list_rnorm.log.quant.ruv.r0 <- runTests(ID, counts_rnorm.log.quant, test_pairs,  cols_for_analysis, sample_rows_lists, type_of_grouping, 
#                                   design_matrix = design_mat_cln, formula_string= paste0("~ 0 + ", args$group_id),
#                                   contrast_variable=args$group_id,
#                                   weights=NA) 
#   
#   myRes_rnorm.log.quant <- extract_results(list_rnorm.log.quant.ruv.r0)
# 
# } else if ( limma_method == "contrasts") {


list_rnorm.log.quant.ruv.r0 <- runTestsContrasts(counts_rnorm.log.quant,
                                                 contrast_strings = contrasts_tbl[, 1][[1]],
                                                 design_matrix = design_mat_cln,
                                                 formula_string = args$formula_string,
                                                 weights = NA)

myRes_rnorm.log.quant <- list_rnorm.log.quant.ruv.r0$results
# }


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print( "Save the list of differentially expressed proteins, without using RUV.")
# save_de_protein_list(myRes_rnorm.log.quant, 
#                      row_id="protein_id",
#                      sort_by_column =q.mod, 
#                      results_dir = args$output_dir,
#                      file_suffix = "_without_ruv.tsv")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Find the list of negative control genes using ANOVA.")
control_genes_index <- getNegCtrlProtAnova(counts_rnorm.log.quant,
                                           design_matrix = design_mat_cln,
                                           group_column = "group",
                                           num_neg_ctrl = args$num_neg_ctrl,
                                           q_val_thresh = args$q_val_thresh)

vroom::vroom_write(data.frame(temp_col = names(control_genes_index),
                              is_control_genes = control_genes_index) %>%
                     set_colnames(c(args$row_id, "is_control_genes")),
                   file.path(args$output_dir, "ctrl_genes_list_ruv3.tsv"),
                   delim = "\t")

loginfo("Total Num. Genes %d", length(control_genes_index))

loginfo("Num. Control Genes %d", length(which(control_genes_index)))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Draw canonical correlation plot.")
ruv_groups <- data.frame(temp_column = colnames(counts_rnorm.log.quant)) %>%
  dplyr::rename(!!rlang::sym(args$sample_id) := "temp_column") %>%
  left_join(design_mat_cln %>%
              dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id))),
            by = args$sample_id)

cancorplot_r1 <- ruv_cancorplot(t(counts_rnorm.log.quant),
                                X = ruv_groups %>% pull(!!rlang::sym(args$group_id)),
                                ctl = control_genes_index)

#cancorplot_r1

ggsave(plot = cancorplot_r1, filename = file.path(args$output_dir, "cancor_plot_round_1.pdf"),limitsize = FALSE)
ggsave(plot = cancorplot_r1, filename = file.path(args$output_dir, "cancor_plot_round_1.png"),limitsize = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Run RUVIII on the input counts table.")
counts_rnorm.log.ruvIII_v1 <- cmriRUVfit(counts_rnorm.log.quant, X = ruv_groups$group, control_genes_index, Z = 1, k = args$ruv_k,
                                         method = args$ruv_method, M = ruvIII_replicates_matrix) %>% t()

vroom::vroom_write(counts_rnorm.log.ruvIII_v1 %>%
                     as.data.frame %>%
                     rownames_to_column(args$row_id),
                   file.path(args$output_dir, "normalized_counts_after_ruv.tsv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Compare the different experimental groups and obtain lists of differentially expressed proteins.")

list_rnorm.log.quant.ruv.r1 <- NA
myRes_rnorm.log.quant.ruv.r1 <- NA
# if( limma_method == "pairs") {
#   list_rnorm.log.quant.ruv.r1 <- runTests( ID, counts_rnorm.log.ruvIII_v1, test_pairs,  
#                                            cols_for_analysis, sample_rows_lists,
#                                            type_of_grouping, 
#                                            design_matrix = design_mat_cln, 
#                                            formula_string= "~ 0 + group ",
#                                            contrast_variable="group",
#                                            weights=NA)
#   
#   myRes_rnorm.log.quant.ruv.r1 <- extract_results(list_rnorm.log.quant.ruv.r1)
# 
# } else if ( limma_method == "contrasts") {

list_rnorm.log.quant.ruv.r1 <- runTestsContrasts(counts_rnorm.log.ruvIII_v1,
                                                 contrast_strings = contrasts_tbl[, 1][[1]],
                                                 design_matrix = design_mat_cln,
                                                 formula_string = "~ 0 + group ",
                                                 weights = NA)

myRes_rnorm.log.quant.ruv.r1 <- list_rnorm.log.quant.ruv.r1$results

# }


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print("Save the lists of diffierentially expressed proteins (e.g. log fold-change, q-values).")
# save_de_protein_list(myRes_rnorm.log.quant.ruv.r1, 
#                      row_id="protein_id",
#                      sort_by_column =q.mod, 
#                      results_dir = args$output_dir,
#                      file_suffix = "_round_1_ruv.tsv")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Draw the RLE and PCA plots.")
rle_pca_plots_arranged <- rlePcaPlotList(list_of_data_matrix = list(counts_rnorm.log.quant, counts_rnorm.log.ruvIII_v1),
                                         design_matrix = design_mat_cln,
                                         sample_id_column = !!rlang::sym(args$sample_id),
                                         group_column = !!rlang::sym(args$group_id),
                                         list_of_descriptions = list("Before RUVIII", "After RUVIII"))

pdf(file.path(args$output_dir, "rle_pca_plots.pdf"))

rle_pca_plots_arranged

dev.off()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare data for drawing the volcano plots")
selected_data <- getSignificantData(list_of_de_tables = list(myRes_rnorm.log.quant,
                                                             myRes_rnorm.log.quant.ruv.r1),
                                    list_of_descriptions = list("No RUV",
                                                                  "RUV applied"),
                                    row_id = !!sym(args$row_id),
                                    p_value_column = p.mod,
                                    q_value_column = q.mod,
                                    log_q_value_column = lqm,
                                    log_fc_column = logFC,
                                    comparison_column = comparison,
                                    expression_column = expression,
                                    facet_column = analysis_type,
                                    q_val_thresh = 0.05) %>%
  dplyr::rename(log2FC = "logFC")

## Write all the results in one single table 
selected_data %>%
  dplyr:::select(-colour, -lqm) %>%
  vroom::vroom_write(path = file.path(args$output_dir, "lfc_qval_long.tsv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Print the volcano plots")

volplot_gg.all <- plotVolcano(selected_data,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = 0.05,
                              formula_string = "analysis_type ~ comparison")


volplot_gg.all

ggsave(filename = file.path(args$output_dir, "volplot_gg.all.png"), plot = volplot_gg.all, width = 7.29, height = 6)
ggsave(filename = file.path(args$output_dir, "volplot_gg.all.svg"), plot = volplot_gg.all, width = 7.29, height = 6)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Count the number of up or down significnat differentially expressed proteins.")

num_sig_de_genes <- printCountDeGenesTable(list_of_de_tables = list(myRes_rnorm.log.quant,
                                                                    myRes_rnorm.log.quant.ruv.r1),
                                           list_of_descriptions = list("No RUV",
                                                                           "RUV applied"),
                                           formula_string = "analysis_type ~ comparison")

num_sig_de_genes$plot

ggsave(filename = file.path(args$output_dir, "num_de_genes_barplot.png"),
       plot = num_sig_de_genes$plot,
       height = 10,
       width = 7)

ggsave(filename = file.path(args$output_dir, "num_de_genes_barplot.svg"),
       plot = num_sig_de_genes$plot,
       height = 10,
       width = 7)

num_sig_de_genes$table

vroom::vroom_write(num_sig_de_genes$table,
                   path = file.path(args$output_dir,
                                    "num_significant_de_genes_all.tab"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Print p-values distribution figure")
pvalhist <- printPValuesDistribution(selected_data,
                                     p_value_column = p.mod,
                                     formula_string = "analysis_type ~ comparison")

pvalhist

ggsave(filename = file.path(args$output_dir, "p_values_distn.png"),
       plot = pvalhist,
       height = 10,
       width = 7)

ggsave(filename = file.path(args$output_dir, "p_values_distn.svg"),
       plot = pvalhist,
       height = 10,
       width = 7)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

norm_counts <- counts_rnorm.log.ruvIII_v1 %>%
  as.data.frame %>%
  set_colnames(paste0(colnames(counts_rnorm.log.ruvIII_v1), ".log2norm")) %>%
  rownames_to_column(args$row_id)


raw_counts <- counts_filt %>%
  as.data.frame %>%
  set_colnames(paste0(colnames(counts_filt), ".raw")) %>%
  rownames_to_column(args$row_id)

de_proteins_wide <- selected_data %>%
  dplyr::filter(analysis_type == "RUV applied") %>%
  dplyr::select(-lqm, -colour, -analysis_type) %>%
  pivot_wider(id_cols = c(!!sym(args$row_id)),
              names_from = c(comparison),
              names_sep = ":",
              values_from = c(log2FC, q.mod, p.mod)) %>%
  left_join(norm_counts, by = args$row_id) %>%
  left_join(raw_counts, by = args$row_id)

head(de_proteins_wide)

vroom::vroom_write(de_proteins_wide, path = file.path(args$output_dir, paste0(args$file_prefix, "_wide.tsv")))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

norm_counts <- counts_rnorm.log.ruvIII_v1 %>%
  as.data.frame %>%
  rownames_to_column(args$row_id) %>%
  pivot_longer(cols = matches(args$group_pattern),
               names_to = "Sample_ID",
               values_to = "log2norm") %>%
  left_join(design_mat_cln, by = "Sample_ID") %>%
  mutate(temp_id = Sample_ID) %>%
  separate(temp_id, sep = "_", into = c("sample_number", "group_pattern")) %>%
  group_by(!!sym(args$row_id), group) %>%
  arrange(!!sym(args$row_id), group, sample_number) %>%
  mutate(replicate_number = paste0("log2norm.", row_number())) %>%
  ungroup %>%
  pivot_wider(id_cols = c(!!sym(args$row_id), group),
              names_from = replicate_number,
              values_from = log2norm)


raw_counts <- counts_filt %>%
  as.data.frame %>%
  rownames_to_column(args$row_id) %>%
  pivot_longer(cols = matches(args$group_pattern),
               names_to = "Sample_ID",
               values_to = "raw") %>%
  left_join(design_mat_cln, by = "Sample_ID") %>%
  mutate(temp_id = Sample_ID) %>%
  separate(temp_id, sep = "_", into = c("sample_number", "group_pattern")) %>%
  group_by(!!sym(args$row_id), group) %>%
  arrange(!!sym(args$row_id), group, sample_number) %>%
  mutate(replicate_number = paste0("raw.", row_number())) %>%
  ungroup %>%
  pivot_wider(id_cols = c(!!sym(args$row_id), group),
              names_from = replicate_number,
              values_from = raw)

left_join_columns <- rlang::set_names(c(args$row_id, "group"),
                                      c(args$row_id, "left_group"))

right_join_columns <- rlang::set_names(c(args$row_id, "group"),
                                       c(args$row_id, "right_group"))

de_proteins_long <- selected_data %>%
  dplyr::filter(analysis_type == "RUV applied") %>%
  dplyr::select(-lqm, -colour, -analysis_type) %>%
  dplyr::mutate(expression = str_replace_all(expression, "group", "")) %>%
  separate(expression, sep = "-", into = c("left_group", "right_group")) %>%
  left_join(norm_counts, by = left_join_columns) %>%
  left_join(norm_counts, by = right_join_columns,
            suffix = c(".left", ".right")) %>%
  left_join(raw_counts, by = left_join_columns) %>%
  left_join(raw_counts, by = right_join_columns,
            suffix = c(".left", ".right"))


#head(de_proteins_long)


vroom::vroom_write(de_proteins_long, path = file.path(args$output_dir, paste0(args$file_prefix, "_long.tsv")))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

