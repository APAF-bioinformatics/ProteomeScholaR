#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",repos = "https://cran.csiro.au/")
  BiocManager::install()
}

# load pacman package manager
if (!require(pacman)) {
  install.packages("pacman",repos = "https://cran.csiro.au/")
  library(pacman)
}

p_load(devtools)
p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(optparse)

p_load(seqinr)
p_load(janitor)
p_load(tictoc)
p_load(configr)
p_load(logging)

p_load(ggpubr)
p_load(plotly)
p_load(limma)
p_load(qvalue)
p_load(ruv)
p_load(UniProt.ws)

devtools::document()