# Tutorial: DIANN Proteomics Data Analysis with ProteomeScholaR

## Table of Contents
1. [Introduction](#introduction)
   - [What is DIA-NN?](#what-is-diann)
   - [Workflow Overview](#workflow-overview)
   - [Prerequisites](#prerequisites)
2. [Setup and Installation](#setup-and-installation)
   - [System Requirements](#system-requirements)
   - [R Environment Setup](#r-environment)
   - [Package Installation](#package-installation)
   - [Troubleshooting Installation](#installation-troubleshooting)
3. [Project Directory Structure](#project-directory-structure)
4. [Data Preparation](#data-preparation)
5. [Experimental Design](#experimental-design)
6. [Protein Identification and Annotation](#protein-identification)
7. [Data Quality Control and Filtering](#quality-control)
8. [Quality Control and Normalization](#quality-control-and-normalization)
9. [Advanced Analysis](#advanced-analysis)

## Introduction <a name="introduction"></a>

### What is DIA-NN? <a name="what-is-diann"></a>

DIA-NN (Data-Independent Acquisition by Neural Networks) is a software suite for the analysis of data-independent acquisition (DIA) proteomics data. This tutorial focuses on processing DIA-NN output using ProteomeScholaR, a comprehensive R package designed for downstream analysis of proteomics data.

Key features of DIA-NN:
- Deep learning-based peptide identification
- Accurate quantification across large-scale experiments
- Support for library-free and spectral library-based analysis
- Advanced retention time prediction
- Built-in quality control metrics

### Workflow Overview <a name="workflow-overview"></a>

The ProteomeScholaR workflow for DIA-NN data analysis consists of several key steps:

1. Data Import and Validation
   - Loading DIA-NN report files
   - Quality assessment of raw data
   - Validation of experimental design

2. Peptide Processing
   - Filtering based on q-values
   - Proteotypic peptide identification
   - Intensity normalization
   - Peptide-to-protein quantification

3. Statistical Analysis
   - Differential expression analysis
   - Multiple testing correction
   - Pathway enrichment
   - Visualization and reporting

### Prerequisites <a name="prerequisites"></a>

Before starting this workflow, ensure you have:

1. Raw Data Requirements:
   - DIA-NN output file (report.tsv)
   - Protein FASTA file for your organism
   - Experimental design information

2. Technical Requirements:
   - R version 4.1.0 or higher
   - At least 8GB RAM (16GB recommended for large datasets)
   - Sufficient disk space (at least 3x the size of your raw data)

3. Knowledge Requirements:
   - Basic R programming
   - Understanding of proteomics experimental design
   - Familiarity with mass spectrometry data analysis

## Setup and Installation <a name="setup-and-installation"></a>

### System Requirements <a name="system-requirements"></a>

#### Hardware Recommendations:
- CPU: Multi-core processor (4+ cores recommended)
- RAM: Minimum 8GB, 16GB+ recommended for large datasets
- Storage: SSD recommended for better performance
- Operating System: Windows 10/11, macOS 10.15+, or Linux

#### R Environment Requirements:
- R version 4.1.0 or higher
- RStudio version 1.4.0 or higher
- Bioconductor 3.13 or higher

### R Environment Setup <a name="r-environment"></a>

1. Install Base R:
   - Download from [CRAN](https://cran.r-project.org/)
   - Windows users: Choose the appropriate installer
   - macOS users: Install R-patched for better performance
   - Linux users: Use package manager or compile from source

2. Install RStudio:
   - Download from [RStudio website](https://www.rstudio.com/products/rstudio/download/)
   - Choose the appropriate version for your OS
   - Verify installation by launching RStudio

3. Configure R Environment:
   ```R
   # Set CRAN mirror
   options(repos = c(CRAN = "https://cloud.r-project.org"))
   
   # Set BiocManager repository
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   ```

### Package Installation <a name="package-installation"></a>

ProteomeScholaR uses an automated installation function that handles all dependencies:

```R
# Install and load ProteomeScholaR
installProteomeScholaR <- function(verbose = TRUE) {
    # Helper function explanation:
    # This function checks for package existence and installs if missing
    # Parameters:
    #   verbose: Boolean, controls installation messages
    # Returns:
    #   None, installs packages as side effect
    
    # Install helper function
    install_if_missing <- function(pkg, source = "CRAN") {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (verbose) message(sprintf("Installing %s from %s...", pkg, source))
            switch(source,
                "CRAN" = install.packages(pkg),
                "Bioc" = BiocManager::install(pkg, update = FALSE),
                "Github" = devtools::install_github(pkg, force = TRUE)
            )
        } else if (verbose) message(sprintf("%s already installed, skipping...", pkg))
    }

    # Core dependencies installation
    install_if_missing("devtools")
    install_if_missing("BiocManager")
    
    # Bioconductor packages
    bioc_packages <- c("clusterProfiler", "GO.db", "UniProt.ws", "mixOmics")
    invisible(sapply(bioc_packages, function(pkg) install_if_missing(pkg, "Bioc")))
}

# Execute installation
installProteomeScholaR()
loadDependencies()
```

### Troubleshooting Installation <a name="installation-troubleshooting"></a>

#### Common Installation Issues:

1. Package Installation Failures:
   ```
   Error: package 'X' is not available (for R version x.y.z)
   ```
   Solutions:
   - Update R to latest version
   - Install package dependencies manually
   - Check R version compatibility
   
2. BiocManager Issues:
   ```
   Error in BiocManager::install()
   ```
   Solutions:
   - Reset BiocManager:
     ```R
     remove.packages("BiocManager")
     install.packages("BiocManager")
     BiocManager::install(version = "3.14")
     ```
   - Check internet connection
   - Verify proxy settings

3. Dependency Conflicts:
   ```
   Error: package 'X' conflicts with package 'Y'
   ```
   Solutions:
   - Create new R environment
   - Update all packages:
     ```R
     update.packages(ask = FALSE)
     ```
   - Remove conflicting packages and reinstall

4. Memory Issues During Installation:
   ```
   Error: cannot allocate vector of size X Mb
   ```
   Solutions:
   - Increase R memory limit:
     ```R
     memory.limit(size = 8000)  # Windows only
     ```
   - Close other applications
   - Install packages individually

#### Installation Verification:

After installation, verify the setup:
```R
# Check package versions
sessionInfo()

# Verify key functionality
library(ProteomeScholaR)
checkDependencies()  # Custom function to verify all dependencies

# Test data loading
testDataLoad()  # Custom function to test data import
```

#### Getting Help:

If you encounter persistent issues:
1. Check the [ProteomeScholaR GitHub Issues](https://github.com/APAF-bioinformatics/ProteomeScholaR/issues)
2. Search the [Bioconductor Support Site](https://support.bioconductor.org/)
3. Contact package maintainers with:
   - Complete error message
   - R session information
   - Operating system details
   - Steps to reproduce the issue

## Project Directory Structure <a name="project-directory-structure"></a>

### Overview

ProteomeScholaR implements a standardized directory structure to ensure reproducible analysis and organized data management. This structure is crucial for:
- Maintaining data integrity
- Ensuring reproducibility
- Facilitating collaboration
- Managing analysis outputs effectively

### Directory Setup Functions

#### 1. setupDirectories()

```R
setupDirectories()
```

**Function Purpose:**
- Creates a standardized directory structure for your analysis
- Sets up necessary subdirectories for different data types
- Establishes paths for results and visualization outputs

**Created Directory Structure:**
```
project_root/
├── data/
│   ├── proteomics/      # DIA-NN output files
│   ├── UniProt/         # FASTA and protein annotation files
│   └── metadata/        # Experimental design and metadata
├── results/
│   ├── qc/             # Quality control outputs
│   ├── normalization/  # Normalization results
│   └── analysis/       # Final analysis results
├── graphs/
│   ├── qc/             # Quality control visualizations
│   ├── exploratory/    # Exploratory data analysis plots
│   └── publication/    # Publication-ready figures
└── logs/               # Analysis logs and reports
```

**Function Parameters:**
```R
setupDirectories(
    base_dir = NULL,          # Base directory path (default: current working directory)
    create_missing = TRUE,    # Create directories if they don't exist
    verbose = TRUE           # Print progress messages
)
```

**Data Considerations:**
1. Directory Permissions:
   - Ensure write permissions in the base directory
   - Check available disk space (recommend 10GB minimum)
   - Verify path length limitations (Windows max path < 260 characters)

2. Directory Naming:
   - Avoid special characters in directory names
   - Use consistent naming conventions
   - Keep paths as short as practical

#### 2. showDirectories()

```R
showDirectories()
```

**Function Purpose:**
- Displays the current directory structure
- Verifies directory creation success
- Shows file paths for reference

**Output Example:**
```
Project Directories:
  Base Directory: /path/to/project
  Data Directory: /path/to/project/data
    - Proteomics: /path/to/project/data/proteomics
    - UniProt: /path/to/project/data/UniProt
  Results Directory: /path/to/project/results
  Graphs Directory: /path/to/project/graphs
```

### Directory Management Best Practices

1. Data Organization:
   ```R
   # Example of organizing input files
   file.copy("path/to/diann/report.tsv", 
            file.path(data_dir, "proteomics", "report.tsv"))
   file.copy("path/to/fasta/file.fasta", 
            file.path(data_dir, "UniProt", "file.fasta"))
   ```

2. Working Directory Management:
   ```R
   # Set and verify working directory
   setwd(project_dir)
   current_dir <- getwd()
   if (!dir.exists(file.path(current_dir, "data"))) {
       stop("Project directory not properly set up")
   }
   ```

3. Path Management:
   ```R
   # Use file.path for cross-platform compatibility
   results_path <- file.path(results_dir, "analysis")
   graphs_path <- file.path(graphs_dir, "qc")
   ```

### Troubleshooting Directory Setup

#### Common Issues and Solutions

1. Permission Errors:
   ```
   Error: cannot create directory 'X': Permission denied
   ```
   Solutions:
   - Check directory permissions:
     ```R
     # Windows
     shell("icacls path/to/directory")
     # Unix
     system("ls -la path/to/directory")
     ```
   - Request admin rights or change directory location
   - Use tempdir() for temporary analysis:
     ```R
     temp_project_dir <- file.path(tempdir(), "proteomics_analysis")
     setupDirectories(base_dir = temp_project_dir)
     ```

2. Path Length Issues (Windows):
   ```
   Error: cannot create directory: path too long
   ```
   Solutions:
   - Enable long paths in Windows:
     ```R
     # Check if path is too long
     if (nchar(file.path(base_dir, "very/long/path")) > 260) {
         warning("Path may be too long for Windows")
     }
     ```
   - Use shorter directory names
   - Move project closer to root directory

3. Directory Already Exists:
   ```
   Warning: Directory X already exists
   ```
   Solutions:
   - Verify directory contents:
     ```R
     list.files(path = "existing_directory", recursive = TRUE)
     ```
   - Archive old directory:
     ```R
     if (dir.exists("old_directory")) {
         file.rename("old_directory", 
                    paste0("old_directory_", format(Sys.time(), "%Y%m%d")))
     }
     ```
   - Remove and recreate if empty:
     ```R
     if (length(list.files("empty_directory")) == 0) {
         unlink("empty_directory", recursive = TRUE)
         dir.create("empty_directory")
     }
     ```

4. Disk Space Issues:
   ```
   Error: cannot create directory: No space left on device
   ```
   Solutions:
   - Check available space:
     ```R
     # Windows
     shell("dir", intern = TRUE)
     # Unix
     system("df -h", intern = TRUE)
     ```
   - Clean up temporary files:
     ```R
     unlink(file.path(tempdir(), "*"), recursive = TRUE)
     ```
   - Move to drive with more space:
     ```R
     new_base_dir <- file.path("D:", "proteomics_analysis")
     setupDirectories(base_dir = new_base_dir)
     ```

### Directory Validation

After setup, validate the directory structure:

```R
validateDirectoryStructure <- function() {
    required_dirs <- c("data", "data/proteomics", "data/UniProt",
                      "results", "results/qc", "graphs")
    
    missing_dirs <- character(0)
    for (dir in required_dirs) {
        if (!dir.exists(file.path(getwd(), dir))) {
            missing_dirs <- c(missing_dirs, dir)
        }
    }
    
    if (length(missing_dirs) > 0) {
        warning("Missing directories: ", 
                paste(missing_dirs, collapse = ", "))
        return(FALSE)
    }
    
    return(TRUE)
}

# Run validation
if (!validateDirectoryStructure()) {
    message("Directory structure incomplete. Rerunning setup...")
    setupDirectories()
}
```

### Best Practices for Directory Usage

1. File Organization:
   - Keep raw data separate from processed data
   - Use consistent file naming conventions
   - Document directory contents in README files

2. Version Control:
   - Include directory structure in version control
   - Exclude large data files
   - Document directory changes

3. Backup Considerations:
   - Regular backups of critical directories
   - Version control for configuration files
   - Document backup procedures

4. Performance Optimization:
   - Use local drives for active analysis
   - Consider SSD for improved performance
   - Implement cleanup procedures for temporary files

## Data Preparation <a name="data-preparation"></a>

### Overview

The data preparation phase is critical for successful analysis. This section covers:
- Required input file formats and specifications
- Data validation and quality checks
- Configuration setup for analysis parameters
- Common data preparation issues and solutions

### Required Input Files

#### 1. DIA-NN Output File (report.tsv)

**File Specifications:**
- Format: Tab-separated values (.tsv)
- Default filename: "report.tsv"
- Required columns:
  - Protein.Ids: Protein identifiers
  - Stripped.Sequence: Peptide sequences
  - Q.Value: Peptide-level q-values
  - Global.Q.Value: Global q-values
  - Precursor.Quantity: Raw intensity values
  - Run: Sample identifiers
  - Proteotypic: Proteotypic peptide indicators

**Data Validation:**
```R
validateDIANNReport <- function(file_path) {
    # Read the first few lines to check format
    data_sample <- read.delim(file_path, nrows = 5)
    
    # Required columns
    required_cols <- c("Protein.Ids", "Stripped.Sequence", "Q.Value",
                      "Global.Q.Value", "Precursor.Quantity", "Run")
    
    # Check for required columns
    missing_cols <- setdiff(required_cols, colnames(data_sample))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", 
             paste(missing_cols, collapse = ", "))
    }
    
    # Check data types
    if (!is.numeric(data_sample$Q.Value) || 
        !is.numeric(data_sample$Global.Q.Value)) {
        stop("Q.Value and Global.Q.Value must be numeric")
    }
    
    return(TRUE)
}
```

**Common Issues:**
1. Missing Columns:
   ```R
   # Check column names
   data_cols <- colnames(read.delim("report.tsv", nrows = 1))
   print(data_cols)
   ```

2. Invalid Data Types:
   ```R
   # Check data types
   str(read.delim("report.tsv"))
   ```

3. Empty or Corrupted Files:
   ```R
   # Check file size and content
   file.info("report.tsv")
   ```

#### 2. FASTA File

**File Requirements:**
- Format: Standard FASTA format
- Content: Protein sequences for your organism
- Header format: >UniProt_ID|Additional_Info
- Source: UniProt or custom database

**Validation:**
```R
validateFastaFile <- function(fasta_path) {
    # Check file existence
    if (!file.exists(fasta_path)) {
        stop("FASTA file not found")
    }
    
    # Read first few lines
    con <- file(fasta_path, "r")
    first_lines <- readLines(con, n = 10)
    close(con)
    
    # Check header format
    if (!grepl("^>", first_lines[1])) {
        stop("Invalid FASTA format: First line must start with '>'")
    }
    
    # Check sequence lines
    seq_lines <- first_lines[!grepl("^>", first_lines)]
    if (any(grepl("[^ACDEFGHIKLMNPQRSTVWY]", seq_lines))) {
        warning("Non-standard amino acids found in sequences")
    }
    
    return(TRUE)
}
```

### Configuration Setup

The workflow uses a configuration file (config.ini) to manage analysis parameters:

```R
# Read configuration file
config_list <- readConfigFile(file = file.path(source_dir, "config.ini"))

# Essential file specifications
DIANN_filename <- "report.tsv"
fasta_filename <- "fasta.fasta"

# Organism information
taxon_id <- your_organism_id_here
organism_name <- "your_organism_name_here"
```

**Configuration Parameters:**

1. Global Parameters:
   ```ini
   [globalParameters]
   fasta_file = fasta.fasta
   peptides_input_file = report.tsv
   output_directory = results
   threads = 4
   ```

2. Quality Control Parameters:
   ```ini
   [qcParameters]
   q_value_threshold = 0.01
   min_peptides_per_protein = 2
   min_samples_detected = 3
   ```

3. Analysis Parameters:
   ```ini
   [analysisParameters]
   normalization_method = median
   imputation_method = knn
   ```

### Data Loading and Validation

```R
# Load DIA-NN data
data_tbl <- vroom::vroom(
    file.path(data_dir, "proteomics", DIANN_filename),
    .name_repair = "check_unique"
)

# Load FASTA file path
fasta_file_path <- file.path(data_dir, "UniProt", fasta_filename)

# Update configuration
config_list[["globalParameters"]][["fasta_file"]] <- fasta_filename
config_list[["globalParameters"]][["peptides_input_file"]] <- DIANN_filename
```

### Data Quality Checks

1. Missing Value Assessment:
   ```R
   checkMissingValues <- function(data) {
       # Calculate missing value percentage per column
       missing_pct <- colMeans(is.na(data)) * 100
       
       # Check for columns with high missing values
       high_missing <- names(missing_pct[missing_pct > 50])
       if (length(high_missing) > 0) {
           warning("Columns with >50% missing values: ",
                   paste(high_missing, collapse = ", "))
       }
       
       return(missing_pct)
   }
   ```

2. Intensity Distribution Check:
   ```R
   checkIntensityDistribution <- function(data) {
       # Log transform intensities
       log_intensities <- log2(data$Precursor.Quantity)
       
       # Check for potential issues
       if (any(is.infinite(log_intensities))) {
           warning("Zero or negative intensities detected")
       }
       
       # Basic statistics
       stats <- summary(log_intensities)
       return(stats)
   }
   ```

### Troubleshooting Data Preparation

#### 1. File Format Issues

Problem: Invalid file format or encoding
```
Error: 'report.tsv' does not appear to be a tab-delimited file
```

Solutions:
```R
# Check file encoding
readLines("report.tsv", n = 1)

# Convert file encoding if needed
if (Sys.info()["sysname"] == "Windows") {
    raw_data <- readLines("report.tsv")
    writeLines(raw_data, "report_fixed.tsv", useBytes = TRUE)
}
```

#### 2. Memory Issues with Large Files

Problem: Insufficient memory for data loading
```
Error: cannot allocate vector of size X Gb
```

Solutions:
```R
# Use vroom for memory-efficient reading
data_tbl <- vroom::vroom(
    file.path(data_dir, "proteomics", DIANN_filename),
    .name_repair = "check_unique",
    num_threads = parallel::detectCores() - 1
)

# Or read in chunks
readLargeFile <- function(file_path, chunk_size = 1e6) {
    con <- file(file_path, "r")
    chunks <- list()
    
    while (TRUE) {
        chunk <- read.table(con, nrows = chunk_size)
        if (nrow(chunk) == 0) break
        chunks[[length(chunks) + 1]] <- chunk
    }
    close(con)
    
    return(do.call(rbind, chunks))
}
```

#### 3. Missing or Incorrect Identifiers

Problem: Protein IDs don't match FASTA file
```
Warning: X protein IDs not found in FASTA file
```

Solutions:
```R
# Check ID formats
checkProteinIDs <- function(data, fasta_path) {
    # Extract protein IDs from data
    protein_ids <- unique(unlist(strsplit(data$Protein.Ids, ";")))
    
    # Extract IDs from FASTA
    fasta_ids <- readFastaHeaders(fasta_path)
    
    # Compare
    missing_ids <- setdiff(protein_ids, fasta_ids)
    if (length(missing_ids) > 0) {
        warning("Missing protein IDs in FASTA: ",
                length(missing_ids), " IDs")
    }
    
    return(missing_ids)
}
```

#### 4. Configuration File Issues

Problem: Invalid configuration parameters
```
Error: Invalid parameter X in config file
```

Solutions:
```R
# Validate configuration
validateConfig <- function(config) {
    required_params <- c("fasta_file", "peptides_input_file")
    
    # Check required parameters
    missing_params <- setdiff(required_params,
                            names(config[["globalParameters"]]))
    if (length(missing_params) > 0) {
        stop("Missing required parameters: ",
             paste(missing_params, collapse = ", "))
    }
    
    # Validate numeric parameters
    numeric_params <- c("q_value_threshold", "min_peptides_per_protein")
    for (param in numeric_params) {
        value <- config[["qcParameters"]][[param]]
        if (!is.numeric(value) || value < 0) {
            stop("Invalid numeric parameter: ", param)
        }
    }
    
    return(TRUE)
}
```

### Best Practices

1. Data Organization:
   - Keep raw files backed up
   - Use consistent file naming
   - Document any data preprocessing steps

2. Quality Control:
   - Run validation checks before analysis
   - Document quality metrics
   - Keep track of filtering decisions

3. Configuration Management:
   - Version control configuration files
   - Document parameter choices
   - Maintain parameter change log

## Experimental Design <a name="experimental-design"></a>

### Overview

The experimental design is a crucial component that defines:
- Sample grouping and relationships
- Technical and biological replicates
- Batch information
- Contrasts for differential analysis

### Design Matrix Setup

ProteomeScholaR provides two methods for creating the experimental design matrix:

#### 1. Interactive Method (Recommended for New Users)

The interactive method uses a Shiny application for intuitive design matrix creation:

```R
if (exists("design_matrix", envir = .GlobalEnv)) {
    print("Design matrix already set :) No need to run app again!")
} else {
    RunApplet("designMatrix")
}
```

**Interactive App Features:**
- Visual sample organization
- Drag-and-drop interface
- Real-time validation
- Export functionality

**Using the Interactive App:**
1. Launch the app:
   ```R
   RunApplet("designMatrix")
   ```
2. Import sample names:
   - Automatically from DIA-NN output
   - Manual entry
   - File upload

3. Define groups:
   - Assign samples to conditions
   - Specify replicates
   - Add batch information

4. Validate design:
   - Check balance
   - Verify replication
   - Review contrasts

#### 2. Manual Method (For Reproducibility)

Create or load a pre-defined design matrix:

```R
design_matrix <- read.table(
    file = file.path(source_dir, "design_matrix.tab"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)
```

**Required Design Matrix Format:**
```R
# Example design matrix structure
design_matrix_template <- data.frame(
    Run = c("Sample1", "Sample2", "Sample3", "Sample4"),
    group = c("Control", "Control", "Treatment", "Treatment"),
    replicates = c(1, 2, 1, 2),
    batch = c(1, 1, 2, 2),
    stringsAsFactors = FALSE
)
```

### Design Matrix Validation

```R
validateDesignMatrix <- function(design_matrix, data) {
    # Check required columns
    required_cols <- c("Run", "group", "replicates")
    missing_cols <- setdiff(required_cols, colnames(design_matrix))
    if (length(missing_cols) > 0) {
        stop("Missing required columns in design matrix: ",
             paste(missing_cols, collapse = ", "))
    }
    
    # Check sample matching
    data_samples <- unique(data$Run)
    design_samples <- design_matrix$Run
    unmatched_samples <- setdiff(data_samples, design_samples)
    if (length(unmatched_samples) > 0) {
        stop("Samples in data not found in design matrix: ",
             paste(unmatched_samples, collapse = ", "))
    }
    
    # Check group balance
    group_counts <- table(design_matrix$group)
    if (any(group_counts < 2)) {
        warning("Some groups have less than 2 replicates: ",
                paste(names(group_counts[group_counts < 2]),
                      collapse = ", "))
    }
    
    # Check replicate numbering
    replicate_check <- tapply(design_matrix$replicates,
                            design_matrix$group,
                            function(x) {
                                if (any(duplicated(x))) {
                                    return(FALSE)
                                }
                                if (any(x > length(x))) {
                                    return(FALSE)
                                }
                                return(TRUE)
                            })
    if (!all(replicate_check)) {
        stop("Invalid replicate numbering in groups: ",
             paste(names(replicate_check[!replicate_check]),
                   collapse = ", "))
    }
    
    return(TRUE)
}
```

### Design Considerations

1. Replication:
   - Minimum 3 replicates per group recommended
   - Balance technical and biological replication
   - Account for batch effects

2. Grouping:
   - Clear group definitions
   - Avoid confounding factors
   - Consider time points

3. Batch Effects:
   - Record processing batches
   - Balance groups across batches
   - Include batch controls

### Troubleshooting Design Setup

#### 1. Sample Matching Issues

Problem: Samples in data don't match design matrix
```
Error: Samples in data not found in design matrix
```

Solutions:
```R
# Check sample names
checkSampleNames <- function(data, design) {
    # Get unique sample names
    data_samples <- unique(data$Run)
    design_samples <- design$Run
    
    # Check for mismatches
    missing_in_design <- setdiff(data_samples, design_samples)
    missing_in_data <- setdiff(design_samples, data_samples)
    
    # Report issues
    if (length(missing_in_design) > 0) {
        warning("Samples in data missing from design: ",
                paste(missing_in_design, collapse = ", "))
    }
    if (length(missing_in_data) > 0) {
        warning("Samples in design missing from data: ",
                paste(missing_in_data, collapse = ", "))
    }
    
    # Suggest fixes
    if (length(missing_in_design) > 0 || length(missing_in_data) > 0) {
        message("Consider:")
        message("1. Check for typos in sample names")
        message("2. Update design matrix")
        message("3. Filter data to match design")
    }
}
```

#### 2. Replicate Issues

Problem: Invalid replicate numbering
```
Error: Invalid replicate numbering in groups
```

Solutions:
```R
# Fix replicate numbering
fixReplicateNumbers <- function(design) {
    # Reset replicate numbers within each group
    design_fixed <- design %>%
        group_by(group) %>%
        mutate(replicates = 1:n()) %>%
        ungroup()
    
    return(design_fixed)
}

# Check replicate balance
checkReplicateBalance <- function(design) {
    # Count replicates per group
    rep_counts <- table(design$group)
    
    # Check balance
    if (length(unique(rep_counts)) > 1) {
        warning("Unbalanced design detected: ",
                paste(names(rep_counts), "=", rep_counts,
                      collapse = ", "))
    }
    
    return(rep_counts)
}
```

#### 3. Group Definition Issues

Problem: Poorly defined or confounded groups
```
Warning: Potential confounding between group and batch
```

Solutions:
```R
# Check for confounding
checkConfounding <- function(design) {
    if ("batch" %in% colnames(design)) {
        # Create contingency table
        conf_table <- table(design$group, design$batch)
        
        # Check for complete confounding
        if (any(rowSums(conf_table == 0) == ncol(conf_table) - 1)) {
            warning("Complete confounding detected between group and batch")
        }
        
        return(conf_table)
    }
}
```

### Best Practices

1. Documentation:
   ```R
   # Save design matrix with metadata
   saveDesign <- function(design, file_path) {
       # Add metadata
       attr(design, "creation_date") <- Sys.time()
       attr(design, "version") <- "1.0"
       
       # Save as RDS for metadata preservation
       saveRDS(design, file_path)
       
       # Also save as tab-delimited for sharing
       write.table(design,
                  gsub("\\.rds$", ".tab", file_path),
                  sep = "\t", row.names = FALSE)
   }
   ```

2. Quality Control:
   ```R
   # Design quality metrics
   assessDesignQuality <- function(design) {
       metrics <- list(
           n_samples = nrow(design),
           n_groups = length(unique(design$group)),
           min_reps = min(table(design$group)),
           balanced = length(unique(table(design$group))) == 1
       )
       return(metrics)
   }
   ```

3. Version Control:
   ```R
   # Track design changes
   logDesignChange <- function(design, change_description) {
       if (is.null(attr(design, "change_log"))) {
           attr(design, "change_log") <- list()
       }
       
       log_entry <- list(
           timestamp = Sys.time(),
           description = change_description,
           user = Sys.info()["user"]
       )
       
       attr(design, "change_log") <- c(attr(design, "change_log"),
                                     list(log_entry))
       
       return(design)
   }
   ```

### Advanced Design Considerations

1. Complex Designs:
   - Nested factors
   - Time series
   - Crossover designs

2. Power Analysis:
   - Sample size estimation
   - Effect size consideration
   - Variance components

3. Randomization:
   - Processing order
   - Batch assignment
   - Technical replication

## Protein Identification and Annotation <a name="protein-identification"></a>

### Overview

The protein identification and annotation phase is critical for:
- Converting protein identifiers to standard formats
- Mapping peptides to proteins
- Creating structured data objects for analysis
- Ensuring data quality and consistency

### FASTA Processing and Protein ID Conversion

#### 1. processFastaFile Function

```R
fasta_meta_file <- "parsed_fasta_data.rds"
aa_seq_tbl_final <- processFastaFile(
    fasta_file_path,
    uniprot_search_results,
    uniparc_search_results,
    fasta_meta_file,
    organism_name
)
```

**Function Purpose:**
- Parses FASTA file headers and sequences
- Extracts protein identifiers and metadata
- Maps identifiers to UniProt/UniParc accessions
- Creates standardized protein annotation table

**Parameters:**
```R
processFastaFile(
    fasta_path,              # Path to FASTA file
    uniprot_results = NULL,  # UniProt mapping results (optional)
    uniparc_results = NULL,  # UniParc mapping results (optional)
    meta_file = NULL,        # Output file for parsed metadata
    organism = NULL          # Organism name for annotation
)
```

**Data Considerations:**
1. FASTA Header Format:
   ```
   >sp|P12345|GENE_HUMAN Protein description
   >tr|A0A123|GENE2_HUMAN_2 Isoform description
   ```
   - Must contain unique identifiers
   - Should follow UniProt format when possible
   - Can include additional metadata

2. Sequence Requirements:
   - Standard amino acid letters
   - No special characters
   - No ambiguous residues

3. File Size Management:
   - Large FASTA files may require chunked processing
   - Consider memory limitations
   - Use caching for repeated analyses

#### 2. Protein ID Mapping

```R
# Update protein IDs in data
data_cln <- updateProteinIDs(data_cln, aa_seq_tbl_final)
```

**Function Purpose:**
- Maps protein IDs from DIA-NN to standard accessions
- Updates protein annotations
- Validates ID consistency
- Handles multiple protein assignments

**ID Mapping Process:**
1. Extract unique protein IDs:
   ```R
   extractProteinIDs <- function(data) {
       # Split multiple protein assignments
       ids <- unlist(strsplit(data$Protein.Ids, ";"))
       # Remove duplicates
       unique_ids <- unique(ids)
       return(unique_ids)
   }
   ```

2. Map to reference database:
   ```R
   mapProteinIDs <- function(ids, reference_db) {
       # Match IDs to reference
       mapped_ids <- match(ids, reference_db$original_id)
       # Get standard accessions
       standard_ids <- reference_db$uniprot_id[mapped_ids]
       return(standard_ids)
   }
   ```

3. Update data with mapped IDs:
   ```R
   updateIDs <- function(data, mapping) {
       # Replace old IDs with new ones
       data$Protein.Ids <- mapply(function(x) {
           old_ids <- unlist(strsplit(x, ";"))
           new_ids <- mapping[old_ids]
           paste(new_ids, collapse = ";")
       }, data$Protein.Ids)
       return(data)
   }
   ```

### PeptideQuantitativeData Object Creation

#### 1. Object Initialization

```R
peptide_data <- new("PeptideQuantitativeData",
    # Data components
    peptide_data = data_cln,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    proteotypic_peptide_sequence_column = "Proteotypic",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = FALSE,
    
    # Experimental design
    design_matrix = design_matrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    
    # Configuration
    args = config_list
)
```

**Object Structure:**
1. Data Components:
   - Peptide-level quantification
   - Protein assignments
   - Quality metrics
   - Normalization status

2. Experimental Design:
   - Sample information
   - Group assignments
   - Replicate structure

3. Configuration:
   - Analysis parameters
   - Processing options
   - QC thresholds

**Data Validation:**
```R
validatePeptideData <- function(object) {
    # Check required columns
    required_cols <- c(
        object@protein_id_column,
        object@peptide_sequence_column,
        object@q_value_column
    )
    
    missing_cols <- setdiff(required_cols,
                           colnames(object@peptide_data))
    
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }
    
    # Validate data types
    if (!is.numeric(object@peptide_data[[object@q_value_column]])) {
        stop("Q-values must be numeric")
    }
    
    if (!is.numeric(object@peptide_data[[object@raw_quantity_column]])) {
        stop("Quantities must be numeric")
    }
    
    # Check sample matching
    data_samples <- unique(object@peptide_data[[object@sample_id]])
    design_samples <- object@design_matrix[[object@sample_id]]
    
    unmatched <- setdiff(data_samples, design_samples)
    if (length(unmatched) > 0) {
        stop("Samples in data not found in design: ",
             paste(unmatched, collapse = ", "))
    }
    
    return(TRUE)
}
```

### Troubleshooting

#### 1. FASTA Processing Issues

Problem: Invalid FASTA format
```
Error: Invalid FASTA header format
```

Solutions:
```R
# Check FASTA format
checkFastaFormat <- function(fasta_path) {
    # Read first few entries
    con <- file(fasta_path, "r")
    lines <- readLines(con, n = 100)
    close(con)
    
    # Check header format
    headers <- lines[grep("^>", lines)]
    if (length(headers) == 0) {
        stop("No FASTA headers found")
    }
    
    # Validate format
    valid_format <- all(grepl("^>[a-z]{2}\\|[A-Z0-9]+\\|", headers))
    if (!valid_format) {
        warning("Non-standard header format detected")
    }
    
    # Check sequences
    seqs <- lines[!grepl("^>", lines)]
    if (any(grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs))) {
        warning("Non-standard amino acids detected")
    }
    
    return(list(headers = headers, sequences = seqs))
}
```

#### 2. ID Mapping Issues

Problem: Unmapped protein IDs
```
Warning: X% of protein IDs could not be mapped
```

Solutions:
```R
# Analyze unmapped IDs
analyzeUnmappedIDs <- function(original_ids, mapped_ids) {
    # Find unmapped
    unmapped <- setdiff(original_ids, names(mapped_ids))
    
    # Analyze patterns
    patterns <- list(
        uniprot = grep("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", unmapped),
        ensembl = grep("^ENS[A-Z]*[PG]\\d+$", unmapped),
        refseq = grep("^[NX]P_\\d+$", unmapped)
    )
    
    # Report findings
    report <- list(
        total_unmapped = length(unmapped),
        by_pattern = sapply(patterns, length),
        example_ids = head(unmapped)
    )
    
    return(report)
}

# Attempt alternative mappings
alternativeIDMapping <- function(unmapped_ids) {
    # Try UniProt API
    uniprot_mapped <- mapToUniprot(unmapped_ids)
    
    # Try RefSeq
    refseq_mapped <- mapToRefSeq(unmapped_ids)
    
    # Combine results
    all_mapped <- c(uniprot_mapped, refseq_mapped)
    
    return(all_mapped)
}
```

#### 3. Data Object Issues

Problem: Invalid object initialization
```
Error: Object initialization failed
```

Solutions:
```R
# Check data structure
checkDataStructure <- function(data, required_cols) {
    # Verify column presence
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Check data types
    expected_types <- list(
        Protein.Ids = "character",
        Q.Value = "numeric",
        Precursor.Quantity = "numeric"
    )
    
    for (col in names(expected_types)) {
        if (col %in% colnames(data)) {
            actual_type <- class(data[[col]])
            if (actual_type != expected_types[[col]]) {
                warning(sprintf("Column %s has type %s, expected %s",
                              col, actual_type, expected_types[[col]]))
            }
        }
    }
    
    return(TRUE)
}

# Fix common issues
fixDataStructure <- function(data) {
    # Convert character Q-values to numeric
    if (is.character(data$Q.Value)) {
        data$Q.Value <- as.numeric(data$Q.Value)
    }
    
    # Handle missing values
    data[is.na(data$Precursor.Quantity), "Precursor.Quantity"] <- 0
    
    # Ensure protein IDs are character
    data$Protein.Ids <- as.character(data$Protein.Ids)
    
    return(data)
}
```

### Best Practices

1. Data Quality:
   ```R
   # Quality checks for protein identification
   checkProteinIdentification <- function(data) {
       metrics <- list(
           total_proteins = length(unique(unlist(
               strsplit(data$Protein.Ids, ";")))),
           multi_protein_peptides = sum(grepl(";", data$Protein.Ids)),
           missing_proteins = sum(is.na(data$Protein.Ids))
       )
       return(metrics)
   }
   ```

2. ID Mapping:
   ```R
   # Document ID mapping process
   documentIDMapping <- function(original_ids, mapped_ids) {
       report <- list(
           timestamp = Sys.time(),
           total_ids = length(original_ids),
           mapped_ids = length(mapped_ids),
           mapping_rate = length(mapped_ids) / length(original_ids),
           unmapped_examples = head(setdiff(original_ids,
                                          names(mapped_ids)))
       )
       return(report)
   }
   ```

3. Object Management:
   ```R
   # Save intermediate results
   saveIntermediateResults <- function(object, step_name) {
       # Create timestamp
       timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
       
       # Save object
       saveRDS(object,
               file = sprintf("results/%s_%s.rds",
                            step_name, timestamp))
       
       # Log saving
       message(sprintf("Saved %s results at %s",
                      step_name, timestamp))
   }
   ```

### Advanced Considerations

1. Protein Groups:
   - Handle shared peptides
   - Resolve protein inference
   - Quantify protein groups

2. Isoform Analysis:
   - Identify unique peptides
   - Map to specific isoforms
   - Quantify isoform-specific expression

3. Quality Metrics:
   - Protein coverage
   - Unique peptide counts
   - Identification confidence scores

## Quality Control and Normalization

### Pre-normalization Quality Control

Quality control before normalization is crucial for ensuring data reliability and guiding normalization strategy. This section covers essential QC steps and their interpretation.

### RLE Plot Analysis

```R
# Generate initial QC composite figure
QC_composite_figure <- InitialiseGrid()

# Create RLE plot
QC_composite_figure@rle_plots$rle_plot_before_cyclic_loess <- plotRle(
    remove_proteins_with_only_one_rep, 
    "group",
    yaxis_limit=c(-4, 4)
)
```

**What is an RLE Plot?**
- Relative Log Expression plots show the distribution of protein expression ratios
- Each box represents a sample's expression relative to the median across all samples
- Ideal plot shows boxes centered at zero with similar spread

**Interpretation Guidelines:**
1. Box Centers:
   - Should be centered around zero
   - Systematic shifts indicate normalization needs
   - Large deviations suggest batch effects

2. Box Spread:
   - Similar spread indicates consistent variance
   - Different spreads suggest heteroscedasticity
   - Extreme outliers need investigation

3. Common Issues:
   - Shifted medians: Systematic bias present
   - Variable box sizes: Inconsistent variance
   - Extreme outliers: Potential sample problems

### PCA Analysis

```R
# Generate PCA plot
QC_composite_figure@pca_plots$pca_plot_before_cyclic_loess_group <- plotPca(
    remove_proteins_with_only_one_rep,
    grouping_variable = "group",
    label_column = "",
    title = "",
    font_size = 8
)

# Get PCA matrix for additional analysis
pca_mixomics_before_cyclic_loess <- getPcaMatrix(remove_proteins_with_only_one_rep)
```

**PCA Plot Interpretation:**
1. Sample Clustering:
   - Look for expected group separation
   - Technical replicates should cluster tightly
   - Batch effects often visible in PC1/PC2

2. Variance Explained:
   - First PCs should capture biological variation
   - Technical variation in later components
   - Check scree plot for variance distribution

3. Quality Indicators:
   - Clear group separation: Strong biological signal
   - Mixed groups: Possible technical issues
   - Outliers: Samples needing investigation

### Correlation Analysis

```R
# Generate correlation plots
QC_composite_figure@pearson_plots$pearson_correlation_pair_before_cyclic_loess <- plotPearson(
    remove_proteins_with_only_one_rep,
    tech_rep_remove_regex = "pool"
)

# Analyze technical replicate correlation
frozen_protein_matrix_tech_rep <- proteinTechRepCorrelation(
    remove_proteins_with_only_one_rep,
    tech_rep_num_column = "group",
    tech_rep_remove_regex = "pool"
)

# Check correlation thresholds
frozen_protein_matrix_tech_rep |> dplyr::filter(pearson > 0.8) |> nrow()
frozen_protein_matrix_tech_rep |> dplyr::filter(spearman > 0.8) |> nrow()
```

**Expected Correlations:**
1. Technical Replicates:
   - Pearson/Spearman > 0.9 expected
   - Lower values indicate technical issues
   - Consider sample removal if consistently low

2. Biological Replicates:
   - Typically > 0.8 correlation
   - Values vary by sample type
   - Consider biological variation

## Normalization Strategy

### 1. Cyclic Loess Normalization

```R
# Apply cyclic loess normalization
normalised_frozen_protein_matrix_obj <- normaliseBetweenSamples(
    remove_proteins_with_only_one_rep,
    normalisation_method = "cyclicloess"
)

# Generate post-normalization QC
QC_composite_figure@rle_plots$rle_plot_before_ruvIIIc_group <- plotRle(
    normalised_frozen_protein_matrix_obj, 
    "group",
    yaxis_limit=c(-4, 4)
)
```

**Why Cyclic Loess?**
- Non-parametric normalization method
- Handles non-linear intensity dependencies
- Robust to outliers
- Preserves biological variation

**Key Considerations:**
1. Advantages:
   - Effective for intensity-dependent bias
   - Robust to outliers
   - Preserves biological differences

2. Limitations:
   - Computationally intensive
   - May not handle batch effects
   - Requires sufficient sample size

### 2. RUV-III Normalization

```R
# Select negative controls
percentage_as_neg_ctrl <- 10
control_genes_index <- getNegCtrlProtAnova(
    normalised_frozen_protein_matrix_obj,
    percentage_as_neg_ctrl = percentage_as_neg_ctrl
)

# Optimize k parameter
cancorplot_r1 <- ruvCancor(
    normalised_frozen_protein_matrix_obj,
    ctrl = control_genes_index,
    num_components_to_impute = 5,
    ruv_grouping_variable = "group"
)

best_k <- findBestK(cancorplot_r1)
```

**RUV-III Strategy:**
1. Control Selection:
   - Uses least variable proteins
   - Typically 5-15% of total proteins
   - ANOVA-based selection

2. Parameter Optimization:
   - k determines factors to remove
   - Balance correction vs. signal
   - Use canonical correlations

3. Implementation:
   ```R
   ruv_normalised_results_temp_obj <- ruvIII_C_Varying(
       normalised_frozen_protein_matrix_obj,
       ruv_grouping_variable = "group",
       ruv_number_k = best_k,
       ctrl = control_genes_index
   )
   ```

### Post-Normalization QC

```R
# Generate final QC plots
QC_composite_figure@rle_plots$rle_plot_after_ruvIIIc_group <- plotRle(
    ruv_normalised_results_cln_obj,
    group="group",
    yaxis_limit=c(-4, 4)
)

QC_composite_figure@pca_plots$pca_plot_after_ruvIIIc_group <- plotPca(
    ruv_normalised_results_cln_obj,
    grouping_variable = "group",
    label_column = "",
    title = "",
    font_size = 8
)
```

**Success Criteria:**
1. RLE Improvements:
   - Centered boxes around zero
   - Similar box sizes
   - Reduced outliers

2. PCA Changes:
   - Better group separation
   - Reduced batch effects
   - Maintained biological signal

3. Correlation Improvements:
   - Higher replicate correlations
   - More consistent patterns
   - Reduced technical variation

## Best Practices and Recommendations

1. QC Workflow:
   - Always perform pre-normalization QC
   - Document all QC metrics
   - Save intermediate results
   - Compare pre/post normalization

2. Parameter Selection:
   - Start conservative
   - Adjust based on QC
   - Document choices
   - Consider sample size

3. Troubleshooting:
   - Investigate outliers early
   - Consider batch effects
   - Document removals
   - Validate decisions

4. Documentation:
   - Save all QC plots
   - Record parameters
   - Note special cases
   - Track protein numbers

## Advanced Analysis <a name="advanced-analysis"></a>

### Overview

After completing the basic workflow, you can proceed with advanced analyses:
1. Differential Expression Analysis
2. Pathway Analysis
3. Protein-Protein Interaction Networks
4. Gene Ontology Enrichment
5. Custom Visualizations

### Differential Expression Analysis

#### 1. Statistical Testing Setup

```R
# Prepare data for differential expression analysis
de_analysis <- performDifferentialExpression(
    protein_obj = remove_proteins_with_only_one_rep,
    design_matrix = design_matrix,
    contrasts = c("Treatment-Control"),  # Specify your contrasts
    min_samples_per_group = 2
)

# Save results
vroom::vroom_write(
    de_analysis$results,
    file.path(results_dir, "de_analysis", "differential_expression_results.tsv")
)
```

**Key Parameters:**
- `min_samples_per_group`: Minimum samples required per condition
- `contrasts`: Comparison groups for analysis
- `protein_obj`: Filtered protein quantification data

#### 2. Results Processing

```R
# Process and filter significant results
significant_proteins <- processDEResults(
    de_results = de_analysis$results,
    p_value_cutoff = 0.05,
    log2fc_cutoff = 1
)

# Generate summary statistics
de_summary <- summarizeDEResults(significant_proteins)
```

**Important Considerations:**
- P-value adjustment method
- Fold change thresholds
- Multiple testing correction

**Function Details:**
1. Data Preparation:
   - Log2 transformation
   - Missing value handling
   - Sample normalization

2. Statistical Testing:
   - Moderated t-tests
   - Empirical Bayes methods
   - Multiple comparison correction

3. Results Filtering:
   - Significance thresholds
   - Effect size cutoffs
   - Quality control metrics

### Pathway Analysis

#### 1. KEGG Pathway Enrichment

```R
# Perform KEGG pathway analysis
kegg_results <- enrichKEGGPathways(
    gene_list = significant_proteins$gene_ids,
    organism = "hsa",  # Use appropriate organism code
    p_value_cutoff = 0.05
)

# Visualize results
plotKEGGEnrichment(
    kegg_results,
    output_dir = file.path(graphs_dir, "pathway_analysis"),
    top_n = 20
)
```

**Function Features:**
- Automated KEGG database queries
- Multiple visualization options
- Statistical significance testing

**Analysis Steps:**
1. Gene ID Conversion:
   - Map protein IDs to KEGG identifiers
   - Handle multiple mappings
   - Validate gene symbols

2. Enrichment Testing:
   - Hypergeometric testing
   - False discovery rate control
   - Background set definition

3. Results Visualization:
   - Enrichment maps
   - Pathway networks
   - Interactive displays

#### 2. Gene Ontology Analysis

```R
# Perform GO enrichment analysis
go_results <- performGOAnalysis(
    gene_list = significant_proteins$gene_ids,
    ontology = c("BP", "MF", "CC"),
    p_value_cutoff = 0.05
)

# Create enrichment plots
createGOEnrichmentPlots(
    go_results,
    output_dir = file.path(graphs_dir, "go_analysis")
)
```

**Analysis Types:**
- Biological Process (BP)
- Molecular Function (MF)
- Cellular Component (CC)

**Key Considerations:**
1. Ontology Selection:
   - Choose relevant aspects
   - Consider specificity levels
   - Account for redundancy

2. Statistical Methods:
   - Term enrichment tests
   - Multiple testing correction
   - Evidence code filtering

3. Visualization Options:
   - Semantic similarity networks
   - Hierarchical layouts
   - Term clustering

### Advanced Visualization

#### 1. Interactive Volcano Plot

```R
# Create interactive volcano plot
createInteractiveVolcano(
    de_results = de_analysis$results,
    p_value_col = "adj.P.Val",
    fc_col = "logFC",
    labels = TRUE,
    output_file = file.path(graphs_dir, "volcano_plot.html")
)
```

**Plot Features:**
- Interactive tooltips
- Customizable thresholds
- Downloadable format
- Protein annotation display

**Customization Options:**
1. Visual Parameters:
   - Point size and opacity
   - Color schemes
   - Axis scales and labels
   - Legend position

2. Interactivity:
   - Hover information
   - Click events
   - Zoom capabilities
   - Selection tools

3. Export Settings:
   - Resolution control
   - File format options
   - Size specifications

#### 2. Protein Expression Heatmap

```R
# Generate expression heatmap
createExpressionHeatmap(
    protein_data = remove_proteins_with_only_one_rep,
    significant_only = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    filename = file.path(graphs_dir, "expression_heatmap.pdf")
)
```

**Heatmap Features:**
- Hierarchical clustering
- Sample annotations
- Color scaling options
- Dendrograms

**Advanced Options:**
1. Clustering Parameters:
   - Distance metrics
   - Linkage methods
   - Row/column ordering
   - Split points

2. Annotation Tracks:
   - Sample metadata
   - Protein features
   - Statistical results
   - Custom categories

3. Visual Customization:
   - Color gradients
   - Border styles
   - Text formatting
   - Legend placement

#### 3. Network Visualization

```R
# Create protein interaction network plot
createNetworkPlot(
    network_data = ppi_network,
    node_color = "degree",
    node_size = "betweenness",
    layout = "force-directed",
    output_file = file.path(graphs_dir, "network_plot.pdf")
)
```

**Network Features:**
- Force-directed layouts
- Node/edge attributes
- Community detection
- Subnetwork highlighting

**Visualization Options:**
1. Layout Algorithms:
   - Force-directed
   - Circular
   - Hierarchical
   - Grid-based

2. Node/Edge Styling:
   - Size mapping
   - Color schemes
   - Shape options
   - Label placement

3. Interactive Elements:
   - Node selection
   - Edge filtering
   - Zoom controls
   - Path highlighting

### Results Integration

#### 1. Multi-omics Integration

```R
# Integrate multiple data types
integrated_results <- integrateOmicsData(
    proteomics = de_analysis$results,
    pathways = kegg_results,
    networks = ppi_network,
    go_terms = go_results
)
```

**Integration Approaches:**
1. Data Layer Integration:
   - Cross-platform normalization
   - Missing data handling
   - Batch effect correction
   - Scale harmonization

2. Statistical Integration:
   - Meta-analysis methods
   - Joint modeling
   - Network fusion
   - Pathway mapping

3. Biological Integration:
   - Functional synthesis
   - Mechanistic insights
   - Regulatory patterns
   - System-level understanding

#### 2. Biological Interpretation

```R
# Generate biological summary
biological_summary <- interpretResults(
    integrated_results,
    output_dir = file.path(results_dir, "interpretation")
)
```

**Interpretation Framework:**
1. Functional Analysis:
   - Pathway activities
   - Process enrichment
   - Regulatory networks
   - Cellular responses

2. Mechanistic Insights:
   - Key regulators
   - Causal relationships
   - Feedback loops
   - System dynamics

3. Clinical Relevance:
   - Biomarker potential
   - Therapeutic targets
   - Disease mechanisms
   - Patient stratification

### Best Practices for Advanced Analysis

1. Statistical Considerations:
   - Use appropriate multiple testing correction
   - Validate assumptions
   - Consider effect sizes
   - Account for technical variation

2. Visualization Guidelines:
   - Maintain consistent styling
   - Provide clear labels
   - Use appropriate scales
   - Include necessary controls

3. Documentation Requirements:
   - Record all parameters
   - Track filtering decisions
   - Document workflow steps
   - Maintain analysis logs

4. Quality Control:
   - Validate results
   - Check data quality
   - Monitor error rates
   - Assess reproducibility