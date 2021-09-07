## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#'Parse the headers of a Uniprot FASTA file and extract the headers and sequences into a data frame
#'Use seqinr object instead as it seems to be a lot faster to run substring
#' @param path to input faster file with header format described in https://www.uniprot.org/help/fasta-headers
#' @return A table containing the following columns:
#' db  sp for Swiss-Prot, tr for TrEMBL
#' uniprot_acc Uniprot Accession
#' uniprot_id  Uniprot ID
#' species     Species
#' tax_id      Taxonomy ID
#' gene_name   Gene symbol
#' protein_evidence 1 to 5, the lower the value, the more evidence that supports the existence of this protein
#' sequence_version Sequence version
#' is_isoform  Is it a protein isoform (not the canonical form)
#' isoform_num     Isoform number.
#' cleaned_acc Cleaned accession without isoform number.
#' status  Reviewed or unreviewed.
#' seq     Amino acid sequence.
#' seq_length      Sequence length (integer).
#' @export
parse_fasta_file <- function( fasta_file) {

    aa_seqinr <-  read.fasta( file = fasta_file,
                             seqtype="AA",
                             whole.header	=TRUE,
                             as.string=TRUE)

  acc_detail_tab <- parse_fasta_object(aa_seqinr)

  names(aa_seqinr) <- str_match( names(aa_seqinr), "(sp|tr)\\|(.+?)\\|(.*)\\s+" )[,3]

  aa_seq_tbl <- acc_detail_tab %>%
                  mutate(seq = map_chr( aa_seqinr, 1)) %>%
                  mutate(seq_length = purrr::map_int(seq, str_length) )
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'Remove rows in the table where the columns specified by the column regular expression pattern are all zero or NA value.

#'@param input_table Input table with columns recording protein abundances for each sample. The name of these columns matches a regular expression pattern, defined by 'col_pattern'. Remove rows with all samples having no protein abundance.
#'@param col_pattern String representing regular expression pattern that matches the name of columns containing the protein abundance values.
#'@param row_id The column name with the row_id, tidyverse style name.
#'@return A data frame with the rows without abundance values removed.
#' @export
remove_empty_rows <- function(input_table, col_pattern, row_id) {

  sites_to_accept <- input_table %>%
    mutate( across( matches(col_pattern, perl=TRUE), ~.==0 )) %>%
    dplyr::filter( !if_all( matches(col_pattern, perl=TRUE), ~. ==TRUE )) %>%
    dplyr::select( {{row_id}})

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  filtered_table <- input_table %>%
    inner_join(sites_to_accept, by=quo_name(enquo(row_id)))

  return( filtered_table )
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
get_plot_num_missing_vales <- function(input_table) {

  plot_num_missing_values <- apply( data.matrix(log2(input_table)), 2,
       function(x) { length(which(is.infinite(x))) } )  %>%
  t %>% t %>%
  set_colnames( "No. of Missing Values") %>%
  as.data.frame %>%
  rownames_to_column("Samples ID")  %>%
  ggplot( aes( x=`Samples ID`, y=`No. of Missing Values` )) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that does have enough number of samples with abundance values.
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param group_column The name of the column in design_matrix table that has the experimental group.
#'@param min_num_samples_per_group An integer representing the minimum number of samples per group.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@return A list, each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
get_rows_to_keep_list <- function( input_table, cols, design_matrix, sample_id, row_id, group_column, min_num_samples_per_group, abundance_threshold) {

  abundance_long <- input_table %>%
    pivot_longer( cols={{cols}},
                  names_to = quo_name(enquo(sample_id)),
                  values_to = "Abundance") %>%
    left_join( design_matrix, by=quo_name(enquo(sample_id))   )


  count_values_per_group <- abundance_long %>%
    mutate(has_value = ifelse( !is.na(Abundance) & Abundance > abundance_threshold, 1, 0 )) %>%
    group_by({{row_id}}, {{group_column}}) %>%
    summarise( num_values = sum( has_value)) %>%
    ungroup()


  kept_rows_temp <- count_values_per_group %>%
    dplyr::filter( num_values >= min_num_samples_per_group) %>%
    dplyr::select(-num_values) %>%
    group_by( {{group_column}}) %>%
    nest( data = c({{row_id}})) %>%
    ungroup() %>%
    mutate( data  = purrr::map( data, ~{   .[,quo_name(enquo(row_id))][[1]]     }   ))


  sample_rows_lists <- kept_rows_temp$data
  names( sample_rows_lists) <- kept_rows_temp[, quo_name(enquo(group_column))][[1]]

  return( sample_rows_lists)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Data imputation function
#'@param df Data matrix
#'@param width Adjustment factor to the observed standard deviation
#'@param downshift Downshift the mean value by this downshift factor multiplied by the observed standard deviation.
#'@return Data matrix with the missing values from each column replaced with a value randomly sampled from the normal distribution with adjusted mean and standard deviation. The normal distribution parameters are based on the observed distribution of the same column.
#'@export
impute_per_col <- function(temp, width = 0.3, downshift = 1.8) {

    temp[!is.finite(temp)] <- NA

  temp.sd <- width * sd(temp, na.rm = TRUE)   # shrink sd width
  temp.mean <- mean(temp, na.rm = TRUE) -
               downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values

  n.missing <- sum(is.na(temp))
  temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  return(temp)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'Converts a design matrix to a biological replicate matrix for use with ruvIII.
#'@param design_matrix The design matrix with the sample ID in one column and the experimental group in another column
#'@param sample_id_column The name of the column with the sample ID, tidyverse style input.
#'@param group_column The name of the column with the experimental group, tidyverse style input.
#'@param temp_column The name of the temporary column that indicates which samples are biological replicates of the same experimental group.
#'@return A numeric matrix with rows as samples, columns as experimental group, and a value of 1 for samples within the same experimental group represented by the same column, and a value of zero otherwise.
#'@export
get_ruvIII_replicate_matrix <- function( design_matrix, sample_id_column, group_column, temp_column=is_replicate_temp) {

  ruvIII_replicates_matrix <- design_mat_cln %>%
  dplyr::select( {{sample_id_column}}, {{group_column}}) %>%
  mutate( {{temp_column}} := 1) %>%
  pivot_wider( id_cols={{sample_id_column}},
               names_from =  {{group_column}},
               values_from = {{temp_column}},
               values_fill = 0) %>%
  column_to_rownames(quo_name(enquo(sample_id_column))) %>%
  as.matrix

  ruvIII_replicates_matrix
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
print_volcano_plot <- function(data,
                    qm.threshold = 0.05,
                    logFC.threshold = 1,
                    logFC = logFC,
                    qvalue=q.mod) {

    comparison <- unique( data$comparison)

     selected_data <-  data %>%
       rownames_to_column( "ID") %>%
       mutate( lqm=-log2({{qvalue}}), qm={{qvalue}}) %>%
       dplyr::select (ID, lqm, qm,  {{logFC}}) %>%
       dplyr::mutate( colour= case_when ( abs({{logFC}}) >= logFC.threshold & qm >= qm.threshold ~ "orange",
                                  abs({{logFC}}) >= logFC.threshold & qm < qm.threshold ~ "purple",
                                  abs({{logFC}}) < logFC.threshold & qm  < qm.threshold ~ "blue",
                                  TRUE ~ "black" )) %>%
      dplyr::mutate( colour = factor( colour, levels=c("black", "orange", "blue", "purple")))

      p <- selected_data %>%
            ggplot( aes(y=lqm, x=logFC))  +
            geom_point(aes(col=colour)) +
            scale_colour_manual(values = c(levels(selected_data$colour)),
                                labels=c(paste0("Not significant, logFC > ",
                                  logFC.threshold),
                                  paste0("Significant, logFC >= ",
                                         logFC.threshold),
                                  paste0("Significant, logFC <",
                                         logFC.threshold),
                                  "Not Significant")) +
            geom_vline(xintercept=1, colour="black", size=0.2) +
            geom_vline(xintercept=-1, colour="black", size=0.2) +
            geom_hline(yintercept=-log2(qm.threshold)) +
            theme_bw() +
            xlab("Log fold changes") +
            ylab("-log2 moderated q-value") +
           # geom_text(data=selected_data[selected_data$colour == "purple",], aes(y=lqm, x=logFC, label=ID), size=3.0) +
            ggtitle( comparison ) +
            theme(legend.position="none")

   return(p)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
plot_PCA <- function( data,
                      design_matrix,
                      sample_id_column = Sample_ID,
                      group_column = group,
                      title, ... ) {

  pca.res <- pca(t(as.matrix(data)) )

  output <- pca.res$variates$X %>%
    as.data.frame %>%
    rownames_to_column(quo_name(enquo(sample_id_column))) %>%
    left_join( design_matrix, by=quo_name(enquo(sample_id_column))) %>%
    ggplot( aes( PC1, PC2, col=group, label=Sample_ID) ) +
    geom_text() +
    labs(title=title) +
    theme(legend.title = element_blank())

  output
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
my_ruv_rle <- function (Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75,
    0.95), ylim = c(-0.5, 0.5))
{
  #  checks = check.ggplot()
   # if (checks) {
        rle = t(apply(t(Y) - apply(Y, 2, median), 2, quantile,
            probs = probs))
        colnames(rle) = c("min", "lower", "middle", "upper",
            "max")
        df = cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))

        if (!is.null(rowinfo)) {
            rowinfo = data.frame(rowinfo= rowinfo)
            df_temp = cbind(df, rowinfo)

            my.x.factor.levels <- df_temp %>%
              arrange( rowinfo) %>%
              distinct( rle.x.factor)  %>%
              pull( rle.x.factor)

            df <- df_temp %>%
              mutate( rle.x.factor= factor(rle.x.factor,
                                              levels=my.x.factor.levels)) %>%
              arrange( rowinfo)
        }

        rleplot = ggplot(df, aes_string(x = "rle.x.factor")) +
            geom_boxplot(aes_string(lower = "lower", middle = "middle",
                upper = "upper", max = "max", min = "min"),
                stat = "identity") +
          theme_bw() +
           theme(axis.title.x = element_blank() ,
             axis.text.x = element_text(angle = 90) #, axis.ticks.x = element_blank()
             ) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
            geom_hline(yintercept = 0) +
          coord_cartesian(ylim = ylim)
        if (!is.null(rowinfo))
            if (ncol(rowinfo) == 1)
                rleplot = rleplot + aes(fill = rowinfo) + labs(fill = "")
        return(rleplot)
   # }
   # else return(FALSE)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
rle_pca_plot_list <- function( list_of_data_matrix, design_matrix,
                               sample_id_column=Sample_ID, group_column=group, list_of_descriptions) {

  rle_list <- purrr::map2( list_of_data_matrix, list_of_descriptions,
                           ~my_ruv_rle(t(as.matrix(  .x)),
                                       rowinfo=design_matrix[colnames(.x), quo_name(enquo(group_column))]) +
                             labs(title=.y) )

  pca_list <- purrr::map2( list_of_data_matrix, list_of_descriptions,
                           ~plot_PCA(.x,
                                     design_matrix = design_matrix,
                                     sample_id_column = {{sample_id_column}},
                                     group_column  = {{group}},
                                     title=.y, cex=7) )

  list_of_plots <- c( rle_list, pca_list)

rle_pca_plots_arranged <- ggarrange(plotlist=list_of_plots, nrow =2 , ncol =length(list_of_descriptions),
  common.legend = FALSE, legend = "bottom", widths = 14, heights = 14)

rle_pca_plots_arranged
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
count_stat_de_genes <- function(data,
                    qm.threshold = 0.05,
                    logFC.threshold = 0,
                    logFC = logFC,
                    qvalue=q.mod ) {

    comparison <- unique( data$comparison)

     selected_data <-  data %>%
       rownames_to_column( "ID") %>%
       mutate( lqm=-log2({{qvalue}}), qm={{qvalue}}) %>%
       dplyr::select (ID, lqm, qm,  {{logFC}}) %>%
       dplyr::mutate( status= case_when (  qm >= qm.threshold ~"Not significant",
                                  {{logFC}} >= logFC.threshold  & qm < qm.threshold ~ "Significant and Up",
                                  {{logFC}} < logFC.threshold  & qm < qm.threshold ~  "Significant and Down",
                                  TRUE ~ "Not significant" ))

     counts <- selected_data %>%
       group_by( status) %>%
       summarise( counts = n()) %>%
       ungroup()

  return(counts)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
print_count_de_genes_table <- function( list_of_de_tables, list_of_descriptions, formula_string="analysis_type ~ comparison", facet_column=analysis_type )  {

  count_stat_de_genes_helper <-  function(de_table, description) {
    purrr::map( de_table, count_stat_de_genes) %>%
      purrr::map2( names(de_table),
                   ~{.x %>% mutate( comparison = .y)}) %>%
      bind_rows %>%
      mutate(   {{facet_column}} := description)
  }

  num_significant_de_genes_all <- purrr::map2( list_of_de_tables,
                                               list_of_descriptions,
                                               function(a, b) { count_stat_de_genes_helper(de_table=a,
                                                                                           description=b) } ) %>%
    bind_rows()

  num_sig_de_genes_barplot <- num_significant_de_genes_all %>%
    dplyr::filter( status != "Not significant" ) %>%
    ggplot(aes( x=status, y=counts )) +
    geom_bar(stat="identity") +
    geom_text(stat='identity', aes(label= counts), vjust=-0.5) +
    theme(axis.text.x = element_text(angle = 90))

  if( !is.na( formula_string)) {

    num_sig_de_genes_barplot <- num_sig_de_genes_barplot +
      facet_grid(  as.formula(formula_string)  )
  }


  return( list(plot=num_sig_de_genes_barplot, table=num_significant_de_genes_all) )
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
get_significant_data <- function( list_of_de_tables,
                                  list_of_descriptions,
                                  row_id = uniprot_acc,
                                  p_value_column = pmod,
                                  q_value_column = q.mod,
                                  log_fc_column = logFC,
                                  comparison_column = comparison,
                                  facet_column = analysis_type,
                                  q_val_thresh = 0.05) {

  get_row_binded_table <- function( de_table_list, description) {
    output <- purrr::map( de_table_list,
                           function(tbl){tbl %>%
                               rownames_to_column(quo_name(enquo(row_id))) %>%
                               dplyr::select({{row_id}},
                                             {{p_value_column}},
                                             {{q_value_column}},
                                             {{log_fc_column}})} ) %>%
      purrr::map2( names( de_table_list), ~{.x %>% mutate( {{comparison_column}} := .y)}) %>%
      bind_rows() %>%
      mutate(   {{facet_column}} := description)

  }

   logfc_tbl_all <- purrr::map2( list_of_de_tables, list_of_descriptions,
                 function(a, b) { get_row_binded_table( de_table_list=a, description=b) } ) %>%
                   bind_rows()

  selected_data <-  logfc_tbl_all %>%
    rownames_to_column( "ID") %>%
    mutate( lqm=-log10(q.mod), qm=q.mod) %>%
    dplyr::select (ID, lqm, qm,  {{p_value_column}}, {{log_fc_column}},
                   {{comparison_column}},
                   {{facet_column}}) %>%
    dplyr::mutate( colour= case_when ( abs({{log_fc_column}}) >= 1 & qm >= q_val_thresh ~ "orange",
                                       abs({{log_fc_column}}) >= 1 & qm < q_val_thresh ~ "purple",
                                       abs({{log_fc_column}}) < 1 & qm  < q_val_thresh ~ "blue",
                                       TRUE ~ "black" )) %>%
    dplyr::mutate( colour = factor( colour, levels=c("black", "orange", "blue", "purple")))

  selected_data


}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
print_volcano_plot <- function(  selected_data,
                       log_q_value_column = lqm,
                       log_fc_column = logFC,
                       q_val_thresh = 0.05,
                       formula_string="analysis_type ~ comparison" ) {

  volplot_gg.all <- selected_data %>%
    ggplot( aes(y={{log_q_value_column}}, x={{log_fc_column}}))  +
    geom_point(aes(col=colour))  +
    scale_colour_manual(values = c(levels(selected_data$colour)),
                        labels=c(paste0("Not significant, logFC > ",
                                        1),
                                 paste0("Significant, logFC >= ",
                                        1),
                                 paste0("Significant, logFC <",
                                        1),
                                 "Not Significant")) +
    geom_vline(xintercept=1, colour="black", size=0.2) +
    geom_vline(xintercept=-1, colour="black", size=0.2) +
    geom_hline(yintercept=-log10(q_val_thresh)) +
    theme_bw() +
    xlab("Log fold changes") +
    ylab("-log10 q-value") +
    theme(legend.position="none")  +
    facet_grid( analysis_type ~ comparison)

  volplot_gg.all
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
print_p_values_distribution <- function (selected_data, p_value_column = pmod, formula_string="is_ruv_applied ~ comparison" ) {

  breaks <- c( 0, 0.001, 0.01, 0.05,
               seq(0.1, 1, by = 0.1))


  pvalhist <- ggplot(selected_data, aes( {{p_value_column}})) +
    theme(axis.title.y = element_blank()) +
    xlab("P-value") +
    geom_histogram( aes_string(y = "..density.."),
                    breaks = breaks,
                    position = "identity",
                    color = "black") +
    geom_histogram( aes_string(y = "..density.."),
                    breaks = breaks,
                    position = "identity")

  if( !is.na( formula_string)) {

    pvalhist <- pvalhist +
      facet_grid(  as.formula(formula_string)  )
  }


  pvalhist

}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
eb.fit <- function( ID, data, design, contr.matrix)
{
    fit <- lmFit(data, design)
    fit.c <- contrasts.fit(fit, contrasts=contr.matrix)

    fit.eb <- suppressWarnings(eBayes(fit.c))

    logFC <- fit.eb$coefficients[, 1]
    df.r <- fit.eb$df.residual
    df.0 <- rep(fit.eb$df.prior, dim(data)[1])
    s2.0 <- rep(fit.eb$s2.prior, dim(data)[1])
    s2 <- (fit.eb$sigma)^2
    s2.post <- fit.eb$s2.post
    t.ord <- fit.eb$coefficients[, 1]/fit.eb$sigma/fit.eb$stdev.unscaled[, 1]
    t.mod <- fit.eb$t[, 1]
    p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
    p.mod <- fit.eb$p.value[, 1]
    q.ord <- qvalue(p.ord)$q
    q.mod <- qvalue(p.mod)$q

    return( list( table =data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post),
                  fit.eb = fit.eb ))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param A String representing the name of experimental group A for pairwise comparison of B - A.
#'@param B String representing the name of experimental group B for pairwise comparison of B - A.
#'@param group_A Names of all the columns / samples that are in experimental group A.
#'@param group_B Names of all the columns / samples that are in experimental group B.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@return A data frame with the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalized log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalized log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' pmod      moderated t-test p-value
#' qval      t-test q-value
#' q.mod     moderated t-test q-value
#'@export
runTest <- function(ID, A, B, group_A, group_B, design_matrix, formula_string,
                    contrast_variable="group",
                    weights=NA ) {


  ff <- as.formula(formula_string)
  mod_frame <- model.frame( ff, design_matrix)
  design_m <- model.matrix( ff, mod_frame)


  # print("My design matrix")
  # print(design_m)
  # print( paste( "nrow(weights)", nrow(weights), "nrow(design_m)", nrow(design_m)))

  if( !is.na(weights)  ) {
    if( nrow(weights) == nrow(design_m) ) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }

  }

  # print(paste("group_A = ", group_A))
  # print(paste("group_B = ", group_B))

  contr.matrix <- makeContrasts( contrasts = paste0( group_B, "vs", group_A , "=", contrast_variable, group_B, "-", contrast_variable, group_A),
                                 levels = colnames(design_m))

  eb_fit_list <- eb.fit(ID, cbind(A, B), design_m, contr.matrix=contr.matrix)

  r <- eb_fit_list$table
  fit.eb <- eb_fit_list$fit.eb

  return ( list(table= data.frame(row.names = row.names(r),
             comparison = paste( "log(", group_B, ") minus log(", group_A, ")",  sep=""),
             meanA     = rowMeans(A),
             meanB     = rowMeans(B),
             logFC     = r$logFC,
             tstats    = r$t.ord,
             tmod      = r$t.mod,
             pval      = r$p.ord,
             pmod      = r$p.mod,
             qval      = r$q.ord,
             q.mod     = r$q.mod),
             fit.eb = fit.eb))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param data
#'@param tests
#'@param sample_columns
#'@param sample_rows_list
#'@param type_of_grouping
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@results A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalized log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalized log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' pmod      moderated t-test p-value
#' qval      t-test q-value
#' q.mod     moderated t-test q-value
#'@export
runTests <- function(ID, data, tests, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable="group",  weights=NA) {
  r <- list()
  for (i in 1:nrow(tests)) {

    rows_to_keep <- rownames(data )


    if( length( sample_rows_list) > 0) {
      if(    !is.na( sample_rows_list ) &
             #  Check that sample group exists as names inside sample_rows_list
             length(which( c( tests[i, "A"],  tests[i, "B"] ) %in% names( sample_rows_list) )) > 0    ) {


        rows_to_keep <- unique( sample_rows_list[[ tests[[i, "A"]]]],
                                sample_rows_list[[ tests[[i, "B"]]]] )
      }
    }

    tmp  <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( tests[i,]$A, tests[i,]$B) )
    A <- tmp[,type_of_grouping[tests[i,]$A][[1]]]
    B <- tmp[,type_of_grouping[tests[i,]$B][[1]]]

    subset_weights <- NA

    if (!is.na(weights)) {
      subset_weights <- weights[ c( colnames(A), colnames(B)),  ]
    }

    # print(colnames(A))
    # print(colnames(B))
    tmp <- unname(cbind(A,B))
    Aname <- paste(tests[i,]$A, 1:max(1,ncol(A)), sep = "_")
    Bname <- paste(tests[i,]$B, 1:max(1,ncol(B)), sep = "_")
    colnames(tmp) <- c(Aname,Bname)

    selected_sample_ids <- c(type_of_grouping[tests[i,]$A][[1]], type_of_grouping[tests[i,]$B][[1]])
    design_matrix_subset <- design_matrix[ selected_sample_ids, , drop=FALSE]

    # print("My design matrix 1")
    # print( selected_sample_ids)
    # print(design_matrix)
    # print( dim(design_matrix))

    group_A <- tests[i,]$A
    group_B <- tests[i,]$B

    x <- runTest(ID, A, B, group_A, group_B, design_matrix=design_matrix_subset,
                 formula_string=formula_string, contrast_variable=contrast_variable,
                 weights=subset_weights)

    comparison <- paste(group_B, " vs ", group_A, sep="")

    r[[comparison]] <- list(results= x$table, counts=t( cbind(A, B)), fit.eb=x$fit.eb )
  }
  r
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## edited from missMethyl source code https://rdrr.io/bioc/missMethyl/src/R/RUVfunctions.R
#'@export
my_RUVfit <- function (Y, X, ctl, Z = 1, k = NULL, method = c("inv", "rinv",
    "ruv4", "ruv3", "ruv2"), M=NULL, ...)
{
    method <- match.arg(method)
    if ((method %in% c("ruv4", "ruv3", "ruv2")) & is.null(k))
        stop("'k' cannot be NULL if method is 'ruv4', 'ruv3' or 'ruv2'.")
    if (mode(ctl) != "logical")
        stop("'ctl' must be a logical vector.")
    if (is.data.frame(Y))
        Y <- data.matrix(Y)
    if (mode(Y) != "numeric")
        stop("'Y' must be a numeric matrix.")
    if ( method == "ruv3" & is.null(M))
        stop("'M' cannot be NULL if method is 'ruv3'")
    Y <- t(Y)
    fit <- switch(method, inv = ruv::RUVinv(Y = Y, X = X,
                                            ctl = ctl, Z = Z, ...),
        rinv = ruv::RUVrinv(Y = Y, X = X, ctl = ctl,
                            Z = Z, k = k, ...),
        ruv4 = ruv::RUV4(Y = Y, X = X, ctl = ctl,
                         k = k, Z = Z, ...),
        ruv2 = ruv::RUV2(Y = Y, X = X, ctl = ctl,
                         k = k, Z = Z, ...),
        ruv3 = ruv::RUVIII(Y= Y, M=M, ctl= ctl, k=k,
                           ...))
    return(fit)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
extract_ruv_results <- function( results_list) {

    extracted <- purrr::map( results_list,  ~{.$results})

    names(extracted) <- names( results_list)

    return(extracted)
}


#'@export
extract_results <- function( results_list) {

    extracted <- purrr::map( results_list,  ~{.$results})

    names(extracted) <- names( results_list)

    return(extracted)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Save the list of output tables from differential expression analysis of proteins or phosphopeptides into a file and in a specific directory.
#'@param list_of_de_tables A list, each element is a table of log fold-change and q-values from differential expression analysis of proteins / phosphopeptides. Each element in the list has a name, usually the name of the pairwise comparison.
#'@param row_id Add row ID to the output table based on the name (protein or phosphopeptid ID) of each row
#'@param sort_by_column Each table in the list_of_de_tables is sorted in ascending order
#'@param results_dir The results directory to store the output file
#'@param file_suffix The file suffix string to aadd to the name of each comparison from the list_of_de_tables.
#'@export
save_de_protein_list <- function( list_of_de_tables, row_id, sort_by_column =q.mod, results_dir, file_suffix) {

  purrr::walk2( list_of_de_tables, names( list_of_de_tables),

~vroom::vroom_write( .x  %>%
    rownames_to_column(row_id) %>%
    arrange( {{sort_by_column}}),
                     path=file.path( results_dir, paste0(.y, file_suffix ) ) ) )

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
analyse_ranking <- function(data, uniprot_acc_column=uniprot_acc) {

  results_tbl <- data %>%
    as.data.frame %>%
  rownames_to_column( quo_name( enquo(uniprot_acc_column)) ) %>%
  dplyr::select( one_of(c( quo_name( enquo(uniprot_acc_column)), "q.mod", "logFC" ))) %>%
  arrange( desc(q.mod)) %>%
  mutate( ctrl_gene_rank  = row_number())

  return(results_tbl)
}

#'@export
get_control_genes <- function( data,
                               q.value =0.05,
                               logFC_threshold = 1,
                               uniprot_acc_column=uniprot_acc) {

  temp <- purrr::map(data,  ~analyse_ranking(., uniprot_acc_column={{uniprot_acc_column}}) )

  ctrl_genes_list <- temp %>%
  bind_rows( .id="Test") %>%
  dplyr::filter( q.mod >= q.value  &
                   abs(logFC) <= logFC_threshold) %>%
  group_by( {{uniprot_acc_column}} ) %>%
  summarise( total_num_set = n()) %>%
  ungroup() %>%
  dplyr::filter( total_num_set == length(data)) %>%
  dplyr::select( {{uniprot_acc_column}} )

  return( ctrl_genes_list)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@param group_column The name of the column with the experimental group, as a string.
#'@export
get_control_genes_anova <- function(data_matrix, design_matrix, group_column= "group", num_neg_ctrl = 500, q_val_thresh = 0.05) {

  ## inspired from matANOVA function from PhosR package

  grps <- design_matrix[colnames(data_matrix), group_column]

  ps <- apply(data_matrix, 1, function(x) {
    summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
  })

  aov <- qvalue(ps)$qvalues

  filtered_list <- aov[aov > q_val_thresh]

  list_size <- ifelse ( num_neg_ctrl > length(filtered_list),  length(filtered_list), num_neg_ctrl)

  control_genes <- names( sort( filtered_list, decreasing=TRUE )[1:list_size]  )

  nrow(data_matrix) - length(control_genes)
  control_genes_index <- rownames( data_matrix)  %in%  control_genes
  names( control_genes_index ) <- rownames( data_matrix)

  return(control_genes_index)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# up <- UniProt.ws(taxId=10090)
# keytypes(up)
# columns(up)
# test <- batch_query_evidence(subset_tbl, Proteins)

clean_isoform_number <- function( string ) {
  # "Q8K4R4-2"
   str_replace( string, "-\\d+$", "")

}


# Filter for a batch and run analysis on that batch of uniprot accession keys only.
subset_query <- function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH")) {


  # print(subset)
  my_keys <- data %>%
             dplyr::filter( round == subset  ) %>%
             pull({{accessions_col_name}})

  # print(head(my_keys))

  UniProt.ws::select( up,
                      keys=my_keys,
                      columns=uniprot_columns,
                      keytype="UNIPROTKB" )
}


# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batch_query_evidence_helper <- function( uniprot_acc_tbl, uniprot_acc_column) {

  all_uniprot_acc <- uniprot_acc_tbl %>%
    dplyr::select( {{uniprot_acc_column}}) %>%
    mutate( Proteins = str_split({{uniprot_acc_column}}, ";")) %>%
    unnest(Proteins) %>%
    distinct %>%
    arrange(Proteins) %>%
    mutate( Proteins = clean_isoform_number(Proteins) )  %>%
    dplyr::mutate(round = ceiling( row_number() / 100 ) )  ## 100 is the maximum number of queries at one time
}

## Run evidence collection online, giving a table of keys (uniprot_acc_tbl) and the column name (uniprot_acc_column)
#'@export
batch_query_evidence <- function(uniprot_acc_tbl, uniprot_acc_column,  uniprot_handle,
                                 uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH")) {

  # uniprot_evidence_levels <- c("Evidence at protein level",
  #                              "Evidence at transcript level",
  #                              "Inferred from homology",
  #                              "Predicted",
  #                              "Uncertain",
  #                              NA)

  all_uniprot_acc <- batch_query_evidence_helper( uniprot_acc_tbl,
                                                  {{uniprot_acc_column}} )

  partial_subset_query <- partial( subset_query,
                                   data=all_uniprot_acc,
                                   accessions_col_name={{uniprot_acc_column}},
                                   uniprot_handle=uniprot_handle,
                                   uniprot_columns=uniprot_columns)

  rounds_list <- all_uniprot_acc %>%
    distinct(round) %>%
    arrange(round) %>%
    pull(round)

  all_uniprot_evidence <- purrr::map( rounds_list, ~{partial_subset_query(subset=.) } ) %>%
    bind_rows

  return(all_uniprot_evidence )
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Choose the best accession
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parse_fasta_file'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#'@export
choose_best_protein_accession <- function( input_tbl, acc_detail_tab, accessions_column, row_id_column = uniprot_acc, group_id) {


  join_condition <-  rlang::set_names( c(quo_name(enquo(row_id_column)), "cleaned_acc" ) ,
                            c(quo_name(enquo(row_id_column)), "cleaned_acc" ) )

  resolve_acc_helper <- input_tbl %>%
    dplyr::select( {{group_id}}, {{accessions_column}}) %>%
    mutate( {{row_id_column}} := str_split( {{accessions_column}}, ";") ) %>%
    unnest( {{row_id_column}} )   %>%
    mutate( cleaned_acc = clean_isoform_number({{accessions_column}})) %>%
    left_join( acc_detail_tab,
               by=join_condition ) %>%
    dplyr::select({{group_id}}, one_of(c( quo_name(enquo(row_id_column)), "gene_name", "cleaned_acc",
                           "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"  ))) %>%
    distinct %>%
    arrange( {{group_id}}, protein_evidence, status, is_isoform,  desc(seq_length), isoform_num )

  score_isoforms <- resolve_acc_helper %>%
    mutate( gene_name = ifelse( is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    group_by( {{group_id}},  gene_name ) %>%
    arrange( {{group_id}},  protein_evidence,
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc )  %>%
    mutate(ranking = row_number()) %>%
    ungroup


  ## For each gene name find the uniprot_acc with the lowest ranking

  my_group_id <- enquo(group_id)

  join_names <- rlang::set_names( c(quo_name(my_group_id), "ranking", "gene_name" ) ,
                            c(quo_name(my_group_id), "ranking", "gene_name" ) )

  group_gene_names_and_uniprot_accs <- score_isoforms  %>%
        distinct( {{group_id}}, gene_name, ranking )  %>%
        dplyr::filter( ranking == 1) %>%
        left_join( score_isoforms %>%
                     dplyr::select( {{group_id}}, ranking, gene_name, uniprot_acc),
                   by = join_names )   %>%
        dplyr::select(-ranking) %>%
    group_by({{group_id}}) %>%
    summarise( num_gene_names = n(),
               gene_names = paste( gene_name, collapse=":"),
               uniprot_acc = paste( uniprot_acc, collapse=":")) %>%
    ungroup() %>%
    mutate( is_unique = case_when( num_gene_names == 1 ~ "Unique",
                                   TRUE ~ "Multimapped"))


  return( group_gene_names_and_uniprot_accs )

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
my_camera <- function(contrast_name, index_name, abundance_mat, replicates_mat, lists_of_contrasts, list_of_gene_sets, min_set_size=4 ) {

    print(paste("contrast_name =", contrast_name))
    print(paste("index_name =", index_name))


    groupA <- str_replace(  contrast_name , "(.*)\\.vs\\.(.*)", "\\1" )
    groupB <- str_replace(  contrast_name , "(.*)\\.vs\\.(.*)", "\\2" )

    design_choose_column <- replicates_mat[,c(groupA, groupB)]

    design_trimmed <- design_choose_column[  rowSums( design_choose_column) > 0,  ]

    # print(design_trimmed)
    #
    # print(head( abundance_mat[[contrast_name]][ , rownames(design_trimmed)]))


    abundance_mat_trimmed <- abundance_mat[[1]][,rownames(design_trimmed)]

    contrast_mat_trimmed <- lists_of_contrasts[[contrast_name]][ colnames(replicates_mat) %in% colnames( design_trimmed)   ]

    index <- list_of_gene_sets[[index_name]]

    msigdb_ids <- geneIds(index)

    #convert gene sets into a list of gene indices
    camera_indices <- ids2indices(msigdb_ids,
                                  rownames(abundance_mat_trimmed))

    ## At least two genes in the gene set
    camera_indices_filt <- camera_indices [ purrr::map ( camera_indices, length)  >= min_set_size ]



    camera_result <- NA
    if( length( camera_indices_filt) > 0 ) {
      camera_result <- camera( y=abundance_mat_trimmed, design=design_trimmed, index=camera_indices_filt, contrast=contrast_mat_trimmed)
    }

    info_list <- list( camera=camera_result,
                       y= abundance_mat_trimmed,
                       design =design_trimmed,

                       index_name = index_name,
                       index=camera_indices_filt,

                       contrast_name= contrast_name,
                       contrast=contrast_mat_trimmed )

    return( info_list)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
run_gsea   <- function( index_name, contrast_name,  list_of_de_proteins, list_of_gene_sets, min_set_size=4) {

      gene_list <-   list_of_de_proteins[[contrast_name]]

      msigdb_gene_set <-  geneIds(list_of_gene_sets[[index_name]])

      query_gene_list <- data.frame(  gene = names( gene_list) )

      term_to_gene_tab <-  tibble( term = names( msigdb_gene_set), gene=msigdb_gene_set ) %>%
        unnest(gene) %>%
        dplyr::inner_join(  query_gene_list, by=c("gene")  )

      terms_to_keep <- term_to_gene_tab %>%
        group_by(term ) %>%
        summarise(counts= n()) %>%
        ungroup() %>%
        dplyr::filter( counts >= min_set_size ) %>%
        dplyr::select(-counts)

      term_to_gene_tab_filt <- term_to_gene_tab %>%
        inner_join( terms_to_keep, by="term") %>%
        mutate( gene = as.character(gene))

      ## Check that there is overlap
      # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length


      gsea_results <- GSEA(geneList=gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt) )

      return( gsea_results)

    }


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
run_enricher   <- function( index_name, contrast_name,  list_of_de_proteins, list_of_gene_sets, min_set_size=4) {

      gene_list <-   list_of_de_proteins[[contrast_name]]

      msigdb_gene_set <-  geneIds(list_of_gene_sets[[index_name]])

      query_gene_list <- data.frame(  gene =  gene_list)

      term_to_gene_tab <-  tibble( term = names( msigdb_gene_set), gene=msigdb_gene_set ) %>%
        unnest(gene) %>%
        dplyr::inner_join( query_gene_list, by=c("gene")  )

      terms_to_keep <- term_to_gene_tab %>%
        group_by(term ) %>%
        summarise(counts= n()) %>%
        ungroup() %>%
        dplyr::filter( counts >= min_set_size ) %>%
        dplyr::select(-counts)

      term_to_gene_tab_filt <- term_to_gene_tab %>%
        inner_join( terms_to_keep, by="term") %>%
        mutate( gene = as.character(gene))

      ## Check that there is overlap
      # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length

       print( intersect(  gene_list ,  unique( term_to_gene_tab_filt$gene )) %>% length )



      gsea_results <- enricher(gene=gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt) )

      return( gsea_results)

    }


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Convert uniprot_acc to entrez_id
## gene with two proteins, choose the one with best q-value
#'@export
merge_with_entrez_id <- function( input_table, lookup_table  ) {

  de_prot_for_camera_helper <- input_table %>%
    inner_join( lookup_table %>% dplyr::filter(!is.na(ENTREZ_GENE)), by=c("uniprot_acc" = "UNIPROTKB") )

  without_entrez_id <- input_table %>%
    left_join( lookup_table %>% dplyr::filter(!is.na(ENTREZ_GENE)), by=c("uniprot_acc" = "UNIPROTKB") ) %>%
    dplyr::filter( is.na(ENTREZ_GENE))

  ## Find rows with duplicated ids
  duplicated_row_id <-  de_prot_for_camera_helper %>%
    dplyr::group_by(ENTREZ_GENE, comparison) %>%
    summarise( counts =n ()) %>%
    ungroup() %>%
    dplyr::filter( counts > 1) %>%
    dplyr::select(-counts)

  to_be_selected <- de_prot_for_camera_helper %>%
    inner_join( duplicated_row_id, by=c("ENTREZ_GENE", "comparison"))

  ## Duplicates with best q-value
  selected_one <- to_be_selected %>%
    inner_join( to_be_selected %>%
                  group_by(comparison) %>%
                  summarise( q.mod = min (q.mod)) %>%
                  ungroup,  by=c( "comparison" = "comparison",
                                  "q.mod" = "q.mod"))

  ## List of rows that were discarded as they are duplicates
  excluded_duplicates <- to_be_selected %>%
    anti_join( selected_one, by=c( "protein_id" = "protein_id",
                                   "comparison" = "comparison"))

  ## Combine those that are not duplicates
  de_prot_for_camera_no_dup <- de_prot_for_camera_helper %>%
    anti_join( duplicated_row_id, by=c("ENTREZ_GENE", "comparison"))

  de_prot_for_camera_cleaned <- selected_one %>%
    bind_rows(de_prot_for_camera_no_dup ) %>%
    dplyr::arrange( comparison, q.mod)

  de_prot_for_camera_tab <- de_prot_for_camera_cleaned %>%
    dplyr::select(-protein_id, -uniprot_acc) %>%
    relocate( ENTREZ_GENE, .before=comparison) %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup()

  return( list( de_proteins=de_prot_for_camera_tab,
                without_entrez_id=without_entrez_id,
                excluded_duplicates=excluded_duplicates ) )
}

