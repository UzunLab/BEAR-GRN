# GRN Inference Method Benchmark Dataset Collection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17087863.svg)](https://doi.org/10.5281/zenodo.17087863)

A comprehensive collection of Gene Regulatory Network (GRN) inference results from multiple state-of-the-art methods across various single-cell multiomics datasets, designed for benchmarking and comparative analysis of new GRN inference methods. Available as both downloadable datasets and an R package (BEAR-GRN).

## 📊 Dataset Overview

This collection contains GRN inference results from **8 established methods** across **7 single-cell datasets**, along with corresponding ground truth regulatory networks for rigorous benchmarking. The dataset includes multiple analysis scripts for comprehensive evaluation of GRN methods and original multiomics data for generating new GRN inferences.

### Available Data Repositories on Zenodo

- **INFERRED.GRNS**: Pre-computed GRN inference results from all methods
- **GROUND.TRUTHS**: Experimentally validated regulatory networks for evaluation
- **STABILITY_GRNS**: Multiple network replicates for stability analysis
- **INPUT.DATA**: Original multiomics data used to infer GRNs (for new method development)
- **INPUT.DATA.STABILITY**: Multiomics data for generating stability analysis networks

### R Package: BEAR-GRN
The complete analysis pipeline is available as the **BEAR-GRN** R package, which includes:
- All analysis functions
- GRN inference scripts for existing methods
- Data preprocessing pipelines
- Comprehensive documentation and examples

### Included Methods
- **CellOracle**: Dynamic GRN inference from scRNA-seq
- **SCENIC+**: Multi-omics GRN inference (extension of Single-Cell rEgulatory Network Inference and Clustering)
- **Pando**: Lineage-informed GRN reconstruction
- **LINGER**: Lifelong learning for GRN inference
- **FigR**: Functional inference of gene regulation
- **TRIPOD**: Three-step regulatory network inference
- **GRaNIE**: Gene regulatory network inference including enhancers
- **DIRECTNET**: Discover cis-Regulatory Elements and Construct TF regulatory NETwork

### Included Datasets
- **K562**: Human chronic myelogenous leukemia cell line
- **Macrophage_S1**: Mouse macrophage stimulation condition 1
- **Macrophage_S2**: Mouse macrophage stimulation condition 2  
- **mESC_E7.5_rep1**: Mouse embryonic stem cells E7.5 replicate 1
- **mESC_E7.5_rep2**: Mouse embryonic stem cells E7.5 replicate 2
- **mESC_E8.5_rep1**: Mouse embryonic stem cells E8.5 replicate 1
- **mESC_E8.5_rep2**: Mouse embryonic stem cells E8.5 replicate 2

## 📁 Directory Structure

### Zenodo Data Repositories
```
INFERRED.GRNS/                          # Pre-computed GRN results
├── K562/
│   ├── CellOracle/
│   │   └── inference_results.tsv
│   ├── SCENIC+/
│   │   └── inference_results.tsv
│   └── [other methods...]
├── Macrophage_S1/
├── Macrophage_S2/
├── mESC_E7.5_rep1/
├── mESC_E7.5_rep2/
├── mESC_E8.5_rep1/
└── mESC_E8.5_rep2/

GROUND.TRUTHS/                          # Experimental validation data
├── filtered_RN117_K562.tsv
├── filtered_RN204_Buffer1.tsv
├── filtered_RN111_E7.5_rep1.tsv
├── filtered_RN111_E7.5_rep2.tsv
├── filtered_RN111_E8.5_rep1.tsv
└── filtered_RN111_E8.5_rep2.tsv

STABILITY_GRNS/                         # Multiple replicates for stability
├── K562/
├── Macrophage_S1/
└── [other samples with replicates...]

INPUT.DATA/                             # Original multiomics data
├── K562/
├── Macrophage_S1/
└── [preprocessed scRNA-seq + scATAC-seq data...]

INPUT.DATA.STABILITY/                   # Data for stability analysis
├── K562/
├── Macrophage_S1/
└── [multiple data splits/replicates...]
```

### BEAR-GRN R Package Structure
```
BEAR-GRN/
├── R/                                  # Analysis functions
│   ├── benchmark_new_method_early_metrics.R
│   ├── benchmark_new_method_filtered.R
│   ├── benchmark_new_method_stability.R
│   ├── reproduce_ROC_PR.R
│   ├── reproduce_early_metrics.R
│   └── reproduce_stability.R
├── man/                                # Function documentation
│   ├── benchmark_new_method.Rd
│   ├── analyze_network_stability.Rd
│   └── [other .Rd files...]
├── inst/scripts/                       # GRN inference pipelines
│   ├── DATA.PREPROCESSING/
│   │   ├── Step1.QC.R
│   │   ├── Step2.filter_cells_genes.R
│   │   ├── Step3.select_common_cells.R
│   │   └── Step4.subsample_cells.R
│   └── GRN.INFERENCE/
│       ├── CellOracle/
│       │   ├── SCRIPT.CELLORACLE
│       │   └── run_ATAC_RNA.sh
│       ├── DIRECTNET/
│       │   └── New.R
│       ├── FigR/
│       │   └── FigR.R
│       ├── GRaNIE/
│       │   └── GRaNIE_singleCell.R
│       ├── LINGER/
│       ├── Pando/
│       │   └── Pando.R
│       ├── SCENIC+/
│       ├── STREAM/
│       │   └── STREAM2.R
│       ├── TRIPOD/
│       │   ├── TRIPOD_Final.R
│       │   └── make_uniq_grn.R
│       └── scGLUE/
│           ├── run_all.sh
│           ├── step1_scGLUE.py
│           ├── step2_scGLUE.py
│           ├── step3_scGLUE.py
│           └── step4_scGLUE.py
└── tests/
    └── testthat/
        └── test-grn_evaluation.R
```

## 📋 Analysis Functions Overview

This benchmark collection includes four main analysis functions:

### 1. `benchmark_new_method()` - New Method Benchmarking
Evaluates new GRN inference methods against existing methods and ground truth.

### 2. `reproduce_ROC_PR_plots()` - ROC/PR Curve Analysis  
Generates comprehensive ROC and Precision-Recall curves for all methods across datasets.

### 3. `reproduce_early_metrics()` - Top-K Edge Analysis
Computes precision, recall, and F1-scores using filtered networks (top 10K edges).

### 4. `analyze_network_stability()` - Stability Analysis
Analyzes network stability across replicates using Jaccard Index.

## 📋 File Formats and Column Structure

### Method-Specific Column Names

| Method | TF Column | Target Column | Score Column | File Format |
|--------|-----------|---------------|--------------|-------------|
| CellOracle | `source` | `target` | `coef_mean` | CSV |
| SCENIC+ | `Source` | `Target` | `Score` | TSV |
| Pando | `tf` | `target` | `estimate` | CSV |
| LINGER | `Source` | `Target` | `Score` | TSV |
| FigR | `Motif` | `DORC` | `Score` | CSV |
| TRIPOD | `TF` | `gene` | `abs_coef` | CSV |
| GRaNIE | `TF.name` | `gene.name` | `TF_gene.r` | CSV |
| DIRECTNET | `TF` | `Gene` | `NA` | CSV |

### Ground Truth Files
- **Format**: TSV (tab-separated values)
- **Columns**: `Source` (TF), `Target` (target gene)
- **Content**: Experimentally validated regulatory interactions

## 🚀 Installation and Setup

### Option 1: R Package Installation (Recommended)

```r
# Install BEAR-GRN R package (if available from repository)
# install.packages("devtools")
# devtools::install_github("your-repo/BEAR-GRN")
# library(BEAR.GRN)

# Or install from local source
# install.packages("path/to/BEAR-GRN", repos = NULL, type = "source")
```

### Option 2: Manual Download and Setup

#### 1. Download Data from Zenodo

```r
# Download all required datasets
download_grn_data <- function(data_dir = "GRN_BENCHMARK_DATA") {
  
  # Create main directory
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Zenodo repository URLs (update with actual DOIs)
  zenodo_files <- list(
    "INFERRED.GRNS" = "https://zenodo.org/record/XXXXX/files/INFERRED.GRNS.zip",
    "GROUND.TRUTHS" = "https://zenodo.org/record/XXXXX/files/GROUND.TRUTHS.zip", 
    "STABILITY_GRNS" = "https://zenodo.org/record/XXXXX/files/STABILITY_GRNS.zip",
    "INPUT.DATA" = "https://zenodo.org/record/XXXXX/files/INPUT.DATA.zip",
    "INPUT.DATA.STABILITY" = "https://zenodo.org/record/XXXXX/files/INPUT.DATA.STABILITY.zip"
  )
  
  # Download and extract each dataset
  for (dataset_name in names(zenodo_files)) {
    cat("Downloading", dataset_name, "...\n")
    
    zip_file <- file.path(data_dir, paste0(dataset_name, ".zip"))
    
    # Download
    tryCatch({
      download.file(zenodo_files[[dataset_name]], zip_file, mode = "wb")
      cat("Downloaded:", dataset_name, "\n")
      
      # Extract
      unzip(zip_file, exdir = data_dir)
      cat("Extracted:", dataset_name, "\n")
      
      # Remove zip file
      file.remove(zip_file)
      
    }, error = function(e) {
      cat("Error downloading", dataset_name, ":", e$message, "\n")
    })
  }
  
  cat("Data download complete! Files saved to:", data_dir, "\n")
  return(data_dir)
}

# Download all data
data_dir <- download_grn_data("GRN_BENCHMARK_DATA")
```

#### 2. Manual Download Alternative

```bash
# Download directly using wget or curl
wget https://zenodo.org/record/XXXXX/files/INFERRED.GRNS.zip
wget https://zenodo.org/record/XXXXX/files/GROUND.TRUTHS.zip
wget https://zenodo.org/record/XXXXX/files/STABILITY_GRNS.zip
wget https://zenodo.org/record/XXXXX/files/INPUT.DATA.zip
wget https://zenodo.org/record/XXXXX/files/INPUT.DATA.STABILITY.zip

# Extract all files
for file in *.zip; do unzip "$file" && rm "$file"; done
```

### 3. Load Analysis Functions

```r
# If using R package
library(BEAR.GRN)

# Or load individual functions from downloaded scripts
source("path/to/BEAR-GRN/R/benchmark_new_method_early_metrics.R")
source("path/to/BEAR-GRN/R/benchmark_new_method_filtered.R")
source("path/to/BEAR-GRN/R/benchmark_new_method_stability.R")
source("path/to/BEAR-GRN/R/reproduce_ROC_PR.R")
source("path/to/BEAR-GRN/R/reproduce_early_metrics.R")
source("path/to/BEAR-GRN/R/reproduce_stability.R")

# Load required libraries
required_packages <- c("readr", "dplyr", "pROC", "PRROC", "ggplot2", 
                      "stringr", "tibble", "RColorBrewer", "plotrix",
                      "tidyr", "purrr", "viridis", "scales")

lapply(required_packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})
```

## 🧬 Developing New Methods with INPUT.DATA

### Using Original Multiomics Data

The `INPUT.DATA` repository contains the preprocessed multiomics datasets used to generate all GRN inferences. This allows you to:

1. **Train your new GRN inference method**
2. **Generate comparable results** using the same input data
3. **Ensure fair comparison** with existing methods

#### Data Format and Contents

Each dataset directory contains:
- **scRNA-seq data**: Gene expression count matrices
- **scATAC-seq data**: Chromatin accessibility peak matrices  
- **Cell metadata**: Cell type annotations and sample information
- **Gene metadata**: Gene symbols, coordinates, and regulatory annotations
- **Peak metadata**: Genomic coordinates and regulatory region annotations

#### Example: Training a New Method

```r
# Load multiomics data for K562
load_multiomics_data <- function(dataset_name = "K562", data_dir = "INPUT.DATA") {
  
  dataset_path <- file.path(data_dir, dataset_name)
  
  # Load scRNA-seq data
  rna_counts <- readRDS(file.path(dataset_path, "rna_counts.rds"))
  rna_metadata <- read.csv(file.path(dataset_path, "rna_metadata.csv"))
  
  # Load scATAC-seq data  
  atac_counts <- readRDS(file.path(dataset_path, "atac_counts.rds"))
  atac_metadata <- read.csv(file.path(dataset_path, "atac_metadata.csv"))
  
  # Load regulatory annotations
  gene_info <- read.csv(file.path(dataset_path, "gene_info.csv"))
  peak_info <- read.csv(file.path(dataset_path, "peak_info.csv"))
  
  return(list(
    rna = list(counts = rna_counts, metadata = rna_metadata),
    atac = list(counts = atac_counts, metadata = atac_metadata),
    genes = gene_info,
    peaks = peak_info
  ))
}

# Example: Use data with your method
k562_data <- load_multiomics_data("K562")

# Your GRN inference method here
# my_grn <- your_method(
#   rna_data = k562_data$rna$counts,
#   atac_data = k562_data$atac$counts,
#   gene_info = k562_data$genes,
#   peak_info = k562_data$peaks
# )

# Save results in standard format for benchmarking
# write_tsv(my_grn, "my_method_k562_results.tsv")
```

#### Data Preprocessing Pipeline

The `inst/scripts/DATA.PREPROCESSING/` directory contains the preprocessing steps:

```r
# Reproduce the preprocessing pipeline
source("inst/scripts/DATA.PREPROCESSING/Step1.QC.R")          # Quality control
source("inst/scripts/DATA.PREPROCESSING/Step2.filter_cells_genes.R")  # Filtering
source("inst/scripts/DATA.PREPROCESSING/Step3.select_common_cells.R")  # Cell matching
source("inst/scripts/DATA.PREPROCESSING/Step4.subsample_cells.R")      # Subsampling
```

### Generating Stability Data

Use `INPUT.DATA.STABILITY` to create multiple network replicates:

```r
# Example: Generate multiple GRN replicates for stability analysis
generate_stability_networks <- function(dataset_name = "K562", 
                                       n_replicates = 5,
                                       data_dir = "INPUT.DATA.STABILITY") {
  
  results <- list()
  
  for (i in 1:n_replicates) {
    cat("Generating replicate", i, "of", n_replicates, "\n")
    
    # Load replicate data (different cell subsets or bootstrap samples)
    replicate_data <- load_stability_replicate(dataset_name, replicate = i, data_dir)
    
    # Run your method
    # grn_replicate <- your_method(replicate_data)
    
    # Save replicate
    # output_file <- paste0("my_method_", dataset_name, "_rep", i, ".tsv")
    # write_tsv(grn_replicate, output_file)
    
    # results[[i]] <- grn_replicate
  }
  
  return(results)
}

# Generate replicates
# k562_replicates <- generate_stability_networks("K562", n_replicates = 10)
```

### Reference Method Implementations

The `inst/scripts/GRN.INFERENCE/` directory provides reference implementations for all methods:

#### Available Method Scripts
- **CellOracle**: `SCRIPT.CELLORACLE`, `run_ATAC_RNA.sh`
- **DIRECTNET**: `New.R`  
- **FigR**: `FigR.R`
- **GRaNIE**: `GRaNIE_singleCell.R`
- **Pando**: `Pando.R`
- **SCENIC+**: Available in subdirectory
- **STREAM**: `STREAM2.R`
- **TRIPOD**: `TRIPOD_Final.R`, `make_uniq_grn.R`
- **scGLUE**: `run_all.sh`, `step1_scGLUE.py` to `step4_scGLUE.py`

#### Example: Run Existing Method

```r
# Example: Run FigR on your data
source("inst/scripts/GRN.INFERENCE/FigR/FigR.R")

# The script contains the complete pipeline:
# 1. Data loading and preprocessing
# 2. Peak-gene linking
# 3. TF-target prediction  
# 4. Network scoring and filtering

# Modify the script to use your own data or different parameters
```

### Complete Workflow: New Method Development

```r
# Step 1: Develop and test your method
dataset_name <- "K562"
input_data <- load_multiomics_data(dataset_name)

# Your method development here
# new_grn <- develop_new_method(input_data)

# Step 2: Generate network for benchmarking  
# write_tsv(new_grn, paste0("MyMethod_", dataset_name, ".tsv"))

# Step 3: Generate stability replicates (optional)
# stability_networks <- generate_stability_networks(dataset_name, n_replicates = 10)

# Step 4: Benchmark against existing methods
benchmark_results <- benchmark_new_method_early_metrics(
  new_grn_file = paste0("MyMethod_", dataset_name, ".tsv"),
  dataset_name = dataset_name,
  tf_column = "TF",
  target_column = "Gene", 
  score_column = "Score",
  method_name = "MyMethod",
  input_dir = "INFERRED.GRNS",
  ground_truth_dir = "GROUND.TRUTHS"
)

# Step 5: Assess stability (if replicates available)
# stability_results <- benchmark_new_method_stability(
#   base_dir = "STABILITY_GRNS",
#   sample_dirs = dataset_name,
#   methods = c("CellOracle", "FigR"),
#   method_name = "MyMethod", 
#   grn_dir = "path/to/stability/replicates"
# )

print(benchmark_results$ranking)
```

## 📊 Analysis Examples

### New Method Benchmarking Examples

#### Example 1: ROC/PR Curve Benchmarking

```r
# Benchmark your new method using ROC and PR curves
results <- benchmark_new_method(
  dataset_name = "K562",
  new_grn_file = "path/to/your/method_results.tsv",
  new_method_name = "YourNewMethod",
  tf_column = "TF",              # Your TF column name
  target_column = "Gene",        # Your target gene column name  
  score_column = "Confidence",   # Your confidence/score column name
  input_dir = "DATASETS",
  output_dir = "benchmark_results",
  ground_truth_dir = "GROUND_TRUTHS"
)

# View performance summary
print(results$results)
print(results$new_method_performance)

# Generated files:
# - ROC_benchmark.pdf/png: ROC curves with your method highlighted
# - PR_benchmark.pdf: Precision-Recall curves
# - benchmark_results.csv: Detailed performance metrics
```

#### Example 2: Precision/Recall/F1 Benchmarking

```r
# Compare your method using precision, recall, and F1-score on top 10K edges
results <- benchmark_new_method_early_metrics(
  new_grn_file = "my_method_K562.csv",
  dataset_name = "K562",
  tf_column = "TF",
  target_column = "Gene",
  score_column = "Weight",       # Can be NULL if no scores
  method_name = "MyMethod",
  input_dir = "DATASETS",
  output_dir = "my_evaluation",
  max_edges = 10000
)

# View ranking results
print(results$ranking)
print(paste("Your method ranks #", 
            results$ranking[results$ranking$Method == "MyMethod", "F1_Rank"],
            " out of ", nrow(results$ranking), " methods"))

# Generated files:
# - K562_MyMethod_Precision_comparison.pdf/png
# - K562_MyMethod_Recall_comparison.pdf/png  
# - K562_MyMethod_F1_score_comparison.pdf/png
# - K562_MyMethod_detailed_results.tsv
# - K562_MyMethod_ranking.tsv
```

#### Example 3: Stability Analysis for New Method

```r
# Analyze stability of your new method (requires multiple network files)
# Option 1: Single dataset with files in one directory
results <- benchmark_new_method_stability(
  base_dir = "STABILITY_GRNS",
  sample_dirs = c("K562"),
  methods = c("CellOracle", "FigR"),  # Existing methods to compare
  method_name = "MyMethod",
  grn_dir = "path/to/my_method/k562_networks",  # Directory with multiple GRN files
  custom_sort_columns = list(MyMethod = "Score"),
  custom_source_target_columns = list(MyMethod = c("Source", "Target")),
  custom_delim = list(MyMethod = "\t"),
  verbose = TRUE,
  save_plot = TRUE,
  save_results = TRUE
)

# Option 2: Multiple datasets with subdirectories
results <- benchmark_new_method_stability(
  base_dir = "STABILITY_GRNS",
  sample_dirs = c("K562", "Macrophage_S1"),
  methods = c("CellOracle", "FigR", "TRIPOD"),
  method_name = "MyNewMethod",
  grn_dir = "path/to/new_method/grns",  # Should contain K562/ and Macrophage_S1/ subdirs
  custom_sort_columns = list(MyNewMethod = "confidence_score"),
  custom_source_target_columns = list(MyNewMethod = c("TF", "Gene")),
  custom_delim = list(MyNewMethod = ","),
  verbose = TRUE
)

# Generated files:
# - median_ji_heatmap.png: Stability heatmap
# - median_ji_results.csv: Detailed Jaccard Index results
```

### Existing Method Analysis Examples

#### Example 4: Generate ROC and PR Curves for All Methods

```r
# Generate comprehensive ROC and PR curves
roc_pr_results <- reproduce_ROC_PR_plots(
  input_dir = "DATASETS",
  output_dir = "ROC_PR_Results",
  ground_truth_dir = "GROUND_TRUTHS",
  use_filtered_approach = TRUE
)

# Access results for specific dataset
k562_results <- roc_pr_results$results[["K562"]]
print(k562_results$results)  # Performance metrics
```

#### Example 5: Analyze Top-K Edge Performance

```r
# Compute precision/recall/F1 for filtered networks (top 10K edges)
metrics_results <- reproduce_early_metrics(
  input_dir = "DATASETS",
  output_dir = "filtered_metrics_results",
  ground_truth_dir = "GROUND_TRUTHS",
  max_edges = 10000
)

# View summary across all datasets
print(metrics_results$combined_summary)

# Access lollipop plots
precision_plot <- metrics_results$plots[["Precision"]]
recall_plot <- metrics_results$plots[["Recall"]]
f1_plot <- metrics_results$plots[["F1-score"]]
```

#### Example 6: Network Stability Analysis

```r
# Analyze network stability using Jaccard Index
# Note: This requires multiple replicates per method/dataset
stability_results <- analyze_network_stability(
  base_dir = "STABILITY_GRNS",  # Directory with multiple replicates
  sample_dirs = c("K562", "Macrophage_S1", "Macrophage_S2"),
  methods = c("CellOracle", "FigR", "LINGER", "Pando"),
  top_percent = 0.1,    # Top 10% of edges
  random_percent = 0.1, # Random 10% of edges
  create_plot = TRUE,
  save_plot = TRUE,
  save_results = TRUE
)

# View results
print(stability_results$results)
print(stability_results$plot)  # Jaccard Index heatmap
```

#### Example 7: Comprehensive New Method Evaluation

```r
# Complete evaluation workflow for a new method
dataset <- "K562"
new_method <- "DeepGRN"
new_grn_file <- "deepgrn_k562_results.tsv"

# Step 1: ROC/PR Analysis
cat("=== Step 1: ROC/PR Analysis ===\n")
roc_pr_results <- benchmark_new_method(
  dataset_name = dataset,
  new_grn_file = new_grn_file,
  new_method_name = new_method,
  tf_column = "TF",
  target_column = "Target", 
  score_column = "Confidence",
  input_dir = "DATASETS",
  output_dir = paste0(new_method, "_ROC_PR_Benchmark")
)

# Step 2: Precision/Recall/F1 Analysis  
cat("=== Step 2: P/R/F1 Analysis ===\n")
prf_results <- benchmark_new_method_early_metrics(
  new_grn_file = new_grn_file,
  dataset_name = dataset,
  tf_column = "TF",
  target_column = "Target",
  score_column = "Confidence",
  method_name = new_method,
  input_dir = "DATASETS",
  output_dir = paste0(new_method, "_PRF_Evaluation")
)

# Step 3: Stability Analysis (if you have multiple network files)
cat("=== Step 3: Stability Analysis ===\n") 
stability_results <- benchmark_new_method_stability(
  base_dir = "STABILITY_GRNS",
  sample_dirs = c(dataset),
  methods = c("CellOracle", "FigR", "LINGER"),
  method_name = new_method,
  grn_dir = "path/to/deepgrn/multiple_runs",
  custom_sort_columns = list(DeepGRN = "Confidence"),
  custom_source_target_columns = list(DeepGRN = c("TF", "Target")),
  save_plot = TRUE,
  output_dir = paste0(new_method, "_Stability_Analysis")
)

# Summary
cat("=== EVALUATION SUMMARY ===\n")
cat("ROC/PR Performance:\n")
print(roc_pr_results$new_method_performance)

cat("\nP/R/F1 Ranking:\n") 
method_rank <- prf_results$ranking[prf_results$ranking$Method == new_method, ]
print(method_rank)

cat("\nStability Results:\n")
stability_subset <- stability_results$results[stability_results$results$Method == new_method, ]
print(stability_subset)
```

## 🔧 Advanced Usage

### Custom Method Configuration

```r
# For methods with non-standard column names
custom_method_info <- tribble(
  ~method,        ~tf_col,     ~target_col,   ~score_col,
  "MyMethod",     "regulator", "regulated",   "weight",
  "CustomGRN",    "TF_name",   "target_name", "score"
)

# Use in ROC/PR analysis
results <- reproduce_ROC_PR_plots(
  input_dir = "DATASETS",
  output_dir = "custom_results", 
  method_info = custom_method_info
)
```

### Debug Network Files

```r
# Debug file structure issues
debug_network_files(
  base_dir = "DATASETS",
  sample = "K562", 
  method = "DIRECTNET"
)
```

### Custom Visualization

```r
# Create custom stability analysis with different parameters
stability_custom <- analyze_network_stability(
  base_dir = "STABILITY_GRNS",
  sample_dirs = c("K562"),
  methods = c("CellOracle", "FigR"),
  top_percent = 0.05,     # Top 5% edges
  random_percent = 0.05,  # Random 5% edges
  plot_width = 10,
  plot_height = 8,
  verbose = TRUE          # Enable detailed logging
)
```

## 📈 Output Files

### New Method Benchmarking Outputs

#### ROC/PR Benchmarking (`benchmark_new_method`)
- **`benchmark_results.csv`**: Comparative AUROC, AUPRC metrics
- **`ROC_benchmark.pdf/png`**: ROC curves with new method highlighted
- **`PR_benchmark.pdf`**: Precision-Recall curves with gap visualization

#### Early Metrics Benchmarking (`benchmark_new_method_early_metrics`)
- **`[Dataset]_[Method]_detailed_results.tsv`**: Comprehensive performance metrics
- **`[Dataset]_[Method]_ranking.tsv`**: Method rankings by P/R/F1
- **`[Dataset]_[Method]_plot_data.tsv`**: Plot data for visualization
- **`[Dataset]_[Method]_Precision_comparison.pdf/png`**: Precision comparison plots
- **`[Dataset]_[Method]_Recall_comparison.pdf/png`**: Recall comparison plots
- **`[Dataset]_[Method]_F1_score_comparison.pdf/png`**: F1-score comparison plots

#### Stability Benchmarking (`benchmark_new_method_stability`)
- **`median_ji_results.csv`**: Jaccard Index stability results
- **`median_ji_heatmap.png`**: Stability heatmap with new method included

### Existing Method Analysis Outputs

#### ROC/PR Analysis Outputs
- **`evaluation_results.csv`**: AUROC, AUPRC metrics for each method
- **`ROC_plot.pdf/png`**: ROC curve comparison plots
- **`PR_gap_plot.pdf`**: Precision-Recall curves with gap visualization
- **`all_datasets_summary.csv`**: Combined results across datasets

#### Filtered Metrics Outputs  
- **`[Dataset]_filtered_metrics.tsv`**: Detailed metrics per dataset
- **`all_datasets_filtered_metrics.tsv`**: Combined summary
- **`Precision_filtered_lollipop_plot.pdf/png`**: Precision lollipop plots
- **`Recall_filtered_lollipop_plot.pdf/png`**: Recall lollipop plots  
- **`F1_score_filtered_lollipop_plot.pdf/png`**: F1-score lollipop plots

#### Stability Analysis Outputs
- **`median_ji_results.csv`**: Jaccard Index results
- **`median_ji_heatmap.png`**: Stability heatmap visualization

## 📊 Performance Metrics

The analysis functions calculate:

### ROC/PR Analysis
- **AUROC**: Area Under the ROC Curve
- **AUPRC**: Area Under the Precision-Recall Curve  
- **AUPRC_random**: Random baseline AUPRC
- **Edge statistics**: Total edges, evaluable edges, positives, negatives

### Filtered Network Analysis
- **Precision**: TP / (TP + FP) for top-K edges
- **Recall**: TP / (TP + FN) for top-K edges
- **F1-Score**: Harmonic mean of precision and recall
- **Network statistics**: Original/filtered/final network sizes

### Stability Analysis
- **Jaccard Index**: |A ∩ B| / |A ∪ B| between network pairs
- **Median JI**: Median Jaccard Index across all pairwise comparisons
- **Top vs Random**: Comparison of stability in top-ranked vs random edges

## 🔧 Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages
```r
install.packages(c(
  "readr",        # File reading
  "dplyr",        # Data manipulation  
  "pROC",         # ROC analysis
  "PRROC",        # PR analysis
  "ggplot2",      # Plotting
  "stringr",      # String operations
  "tibble",       # Data frames
  "RColorBrewer", # Color palettes
  "plotrix",      # Gap plots
  "tidyr",        # Data reshaping
  "purrr",        # Functional programming
  "viridis",      # Color scales
  "scales"        # Scale functions
))
```

## 📊 Ground Truth Information

Ground truth regulatory networks are derived from:
- ChIP-seq experiments
- Perturbation studies  
- Literature-curated databases
- Experimental validation studies

Each ground truth file contains high-confidence regulatory interactions specific to the cell type and experimental condition.

## 🎯 Analysis Workflow Recommendations

### For New Method Developers

#### Complete Evaluation Workflow
1. **Start with `benchmark_new_method_early_metrics()`** to get practical performance metrics (P/R/F1) and method ranking
2. **Use `benchmark_new_method()`** for comprehensive ROC/PR curve analysis and AUROC/AUPRC comparison
3. **Apply `benchmark_new_method_stability()`** if you have multiple network replicates to assess reproducibility
4. **Use `debug_network_files()`** if you encounter file reading issues

#### Quick Assessment Workflow  
1. **`benchmark_new_method_early_metrics()`** - Fast evaluation with clear ranking
2. **`benchmark_new_method()`** - Detailed performance curves if needed

### For Comparative Studies

#### Comprehensive Analysis
1. **Use `reproduce_ROC_PR_plots()`** for comprehensive method comparison across all datasets
2. **Apply `reproduce_early_metrics()`** for practical performance assessment with lollipop visualizations
3. **Combine with `analyze_network_stability()`** for robustness evaluation across methods

#### Focused Analysis
1. **Start with filtered metrics** (`reproduce_early_metrics()`) for practical insights
2. **Examine ROC/PR curves** (`reproduce_ROC_PR_plots()`) for detailed performance characteristics

### For Method Validation

#### New Method Validation
1. **Benchmark against existing methods** using `benchmark_new_method_early_metrics()` for quick validation
2. **Detailed curve analysis** with `benchmark_new_method()` for publication-quality results  
3. **Stability assessment** with `benchmark_new_method_stability()` for reproducibility claims

#### Existing Method Assessment
1. **Start with filtered metrics** to understand practical performance on top edges
2. **Examine ROC/PR curves** for detailed performance characteristics
3. **Assess stability** across replicates to ensure reproducibility

### Function Selection Guide

| Analysis Goal | Recommended Functions | Key Outputs |
|---------------|----------------------|-------------|
| **New Method: Quick Ranking** | `benchmark_new_method_early_metrics()` | Method ranking, P/R/F1 comparison |
| **New Method: Publication Analysis** | `benchmark_new_method()` + `benchmark_new_method_early_metrics()` | ROC/PR curves + rankings |
| **New Method: Full Validation** | All 3 benchmarking functions | Complete assessment |
| **Literature Comparison** | `reproduce_ROC_PR_plots()` + `reproduce_early_metrics()` | Method comparisons |
| **Stability Assessment** | `analyze_network_stability()` or `benchmark_new_method_stability()` | Jaccard Index analysis |
| **Debugging Issues** | `debug_network_files()` | File format inspection |

## 📝 Citation

If you use this benchmark collection or the BEAR-GRN package in your research, please cite:

### Main Dataset Citation
```bibtex
@dataset{grn_benchmark_2024,
  title={GRN Inference Method Benchmark Dataset Collection},
  author={Your Name and Collaborators},
  year={2024},
  publisher={Zenodo},
  doi={10.5281/zenodo.17087863},
  url={https://doi.org/10.5281/zenodo.17087863}
}
```

### R Package Citation (BEAR-GRN)
```bibtex
@software{bear_grn_2024,
  title={BEAR-GRN: Benchmarking and Evaluation of Automated Regulatory Gene Networks},
  author={Your Name and Collaborators},
  year={2024},
  url={https://github.com/your-repo/BEAR-GRN}
}
```

### Individual Data Repository Citations

If using specific data repositories, please cite:

```bibtex
@dataset{inferred_grns_2024,
  title={Pre-computed Gene Regulatory Network Inference Results},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXX}
}

@dataset{input_data_2024,
  title={Multiomics Input Data for GRN Inference},
  publisher={Zenodo}, 
  doi={10.5281/zenodo.XXXXX}
}

@dataset{stability_data_2024,
  title={Gene Regulatory Network Stability Analysis Data},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXX}
}
```

## 🌟 Complete Ecosystem Overview

The BEAR-GRN ecosystem provides a comprehensive platform for GRN method development and evaluation:

### 📊 **Data Resources**
- **INFERRED.GRNS**: Results from 8+ state-of-the-art methods
- **INPUT.DATA**: Original multiomics data for method training  
- **STABILITY_GRNS**: Multiple replicates for robustness assessment
- **GROUND.TRUTHS**: Experimentally validated regulatory networks

### 🔧 **Analysis Tools**  
- **Existing Method Analysis**: Compare published methods
- **New Method Benchmarking**: Evaluate your novel approaches
- **Stability Assessment**: Measure method reproducibility
- **Visualization Suite**: Publication-ready plots and heatmaps

### 📈 **Method Development**
- **Reference Implementations**: Scripts for all existing methods
- **Preprocessing Pipelines**: Standardized data preparation
- **Training Data**: Consistent input for fair comparison
- **Evaluation Metrics**: ROC/PR curves, precision/recall/F1, stability indices

### 💻 **Software Components**
- **R Package (BEAR-GRN)**: Complete analysis suite
- **Standalone Scripts**: Individual analysis functions  
- **Method Implementations**: Ready-to-use inference pipelines
- **Documentation**: Comprehensive guides and examples

This integrated approach ensures that:
1. **Method developers** can train and evaluate new approaches fairly
2. **Researchers** can compare methods systematically  
3. **Results** are reproducible and standardized
4. **Community** has access to both data and tools

### 🎯 **Key Benefits**
- **Standardized Evaluation**: Same data, metrics, and protocols
- **Fair Comparison**: Identical preprocessing and evaluation criteria
- **Reproducible Science**: All code and data openly available
- **Community Resource**: Shared benchmark for the field
- **Continuous Development**: Easy to add new methods and datasets

## 📄 License

This dataset is released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to:
- **Share**: Copy and redistribute the material
- **Adapt**: Remix, transform, and build upon the material

Under the following terms:
- **Attribution**: You must give appropriate credit

## 🤝 Contributing

We welcome contributions to improve this benchmark collection:

1. **New Methods**: Submit results from additional GRN inference methods
2. **New Datasets**: Contribute results on additional single-cell datasets
3. **Analysis Improvements**: Enhance existing analysis functions
4. **Bug Reports**: Report issues with the analysis pipeline
5. **Documentation**: Improve usage examples and documentation

## 📞 Support

For questions, issues, or contributions:
- **GitHub Issues**: [Link to repository issues]
- **Email**: [yuzun@pennstatehealth.psu.edu ; kkaramveer@pennstatehealth.psu.edu]
- **Documentation**: See function documentation and examples for detailed usage

## 🔄 Version History

- **v1.0.0**: Initial release with 8 methods across 7 datasets
  - Added comprehensive analysis functions
  - ROC/PR curve generation
  - Filtered network metrics analysis  
  - Network stability assessment
  - Enhanced benchmarking capabilities

## 🙏 Acknowledgments

This benchmark collection builds upon the excellent work of the developers of:
CellOracle, SCENIC+, Pando, LINGER, FigR, TRIPOD, DIRECTNET, and GRaNIE methods.

Special thanks to the single-cell genomics community for providing high-quality datasets and ground truth regulatory networks.

## 🚨 Troubleshooting

### Common Issues

#### File Reading Errors
```r
# Debug file structure for existing methods
debug_network_files(base_dir = "DATASETS", sample = "K562", method = "DIRECTNET")

# For new method files, manually inspect:
df <- read_delim("your_method_file.tsv", delim = "\t")
colnames(df)  # Check column names
head(df)      # Check data format
```

#### New Method Column Mapping Issues
- Ensure your TF and target column names are specified correctly
- Check that score column exists (use `NULL` if no scores available)
- Verify file format matches expected delimiter (CSV vs TSV)

#### Empty Results for New Methods
- Verify ground truth files exist for your chosen dataset
- Check that TF and target gene names overlap with ground truth
- Ensure score column contains numeric values (not text)

#### Stability Analysis Issues
- Confirm you have multiple network files (at least 2) per method
- Check that file naming is consistent across replicates
- Verify directory structure matches expected format

#### Memory Issues
For large networks, consider:
- Increasing R memory limits: `options(java.parameters = "-Xmx8g")`
- Processing datasets individually rather than in batch
- Using smaller `max_edges` values in filtered analysis
- Subsampling large networks before analysis

#### Method Performance Issues

##### Poor ROC/PR Performance
- Check if TF-target pairs in your network overlap with ground truth
- Verify scoring is meaningful (higher scores = more confident predictions)

##### Unexpected Ranking Results  
- Examine the `detailed_results.tsv` file for edge statistics
- Check TP, FP, FN counts to understand performance bottlenecks
- Compare network sizes - very small networks may have inflated precision

##### Stability Analysis Problems
- Ensure replicate files contain similar edge sets
- Check that file reading parameters are consistent across files
- Verify that custom column mappings are correct

### Debug Workflows

#### Debug New Method Integration
```r
# Step 1: Check file structure
cat("=== File Check ===\n")
your_file <- "path/to/your/method.tsv"
if (file.exists(your_file)) {
  df <- read_delim(your_file, delim = "\t", n_max = 5)
  print(colnames(df))
  print(df)
} else {
  cat("File not found!\n")
}

# Step 2: Test column mapping
tf_col <- "your_tf_column"
target_col <- "your_target_column"  
score_col <- "your_score_column"

if (all(c(tf_col, target_col, score_col) %in% colnames(df))) {
  cat("All columns found!\n")
} else {
  cat("Missing columns. Available:", paste(colnames(df), collapse = ", "), "\n")
}

# Step 3: Check data quality
print(summary(df[[score_col]]))  # Should be numeric
print(table(is.na(df[[score_col]])))  # Check for NAs
```

#### Debug Stability Analysis Setup
```r
# Check directory structure
base_path <- "your/stability/directory"
sample_name <- "K562"
method_name <- "YourMethod"

# List available files
files <- list.files(file.path(base_path, sample_name, method_name), 
                   pattern = "\\.(csv|tsv)$", full.names = TRUE)
cat("Found files:", length(files), "\n")
cat("File names:", paste(basename(files), collapse = ", "), "\n")

# Test file reading
if (length(files) >= 2) {
  df1 <- read_delim(files[1], delim = "\t", n_max = 3)
  df2 <- read_delim(files[2], delim = "\t", n_max = 3)
  
  cat("File 1 columns:", paste(colnames(df1), collapse = ", "), "\n")
  cat("File 2 columns:", paste(colnames(df2), collapse = ", "), "\n")
  
  cat("Columns match:", identical(colnames(df1), colnames(df2)), "\n")
} else {
  cat("Need at least 2 files for stability analysis\n")
}
```

### File Format Requirements

#### New Method File Formats
Your GRN files should follow these guidelines:

**CSV Format Example:**
```csv
TF,Gene,Score
STAT1,IRF1,0.85
STAT1,IRF7,0.72
NF1,MYC,0.91
```

**TSV Format Example:**
```tsv
Source	Target	Score
STAT1	IRF1	0.85
STAT1	IRF7	0.72
NF1	MYC	0.91
```

#### Directory Structure for Stability Analysis

**Option 1: Single dataset**
```
your_method_networks/
├── network_rep1.tsv
├── network_rep2.tsv
├── network_rep3.tsv
└── ...
```

**Option 2: Multiple datasets**
```
your_method_networks/
├── K562/
│   ├── network_rep1.tsv
│   ├── network_rep2.tsv
│   └── network_rep3.tsv
├── Macrophage_S1/
│   ├── network_rep1.tsv
│   ├── network_rep2.tsv
│   └── network_rep3.tsv
└── ...
```

### Performance Interpretation Guidelines

#### ROC/PR Metrics
- **AUROC > 0.7**: Generally considered good performance
- **AUROC > 0.8**: Strong performance  
- **AUROC > 0.9**: Excellent performance
- **AUPRC**: Compare against random baseline (AUPRC_random)

#### Precision/Recall/F1 Metrics
- **High Precision, Low Recall**: Conservative predictions, fewer false positives
- **Low Precision, High Recall**: Liberal predictions, more false positives
- **Balanced F1**: Good compromise between precision and recall
- **Context matters**: Consider the specific biological system and application

#### Stability Metrics (Jaccard Index)
- **JI > 0.5**: Moderate stability between network replicates
- **JI > 0.7**: Good stability
- **JI > 0.8**: High stability
- **Compare Top vs Random**: Top edges should have higher JI than random edges

### Getting Help

If you encounter issues not covered here:

1. **Check the function documentation**: Each function has detailed parameter descriptions
2. **Use verbose mode**: Set `verbose = TRUE` for detailed logging
3. **Test with small datasets first**: Debug with a subset before full analysis
4. **Verify file formats**: Use the debug functions to inspect your data structure
5. **Contact support**: Email the addresses provided in the support section structure for existing methods
debug_network_files(base_dir = "DATASETS", sample = "K562", method = "DIRECTNET")

# For new method files, manually inspect:
df <- read_delim("your_method_file.tsv", delim = "\t")
colnames(df)  # Check column names
head(df)      # Check data format
```

#### New Method Column Mapping Issues
- Ensure your TF and target column names are specified correctly
- Check that score column exists (use `NULL` if no scores available)
- Verify file format matches expected delimiter (CSV vs TSV)

#### Empty Results for New Methods
- Check the path for new method GRNs 
- Ensure you provide correct TF and target column names
- Ensure score column contains numeric values (not text)

#### Stability Analysis Issues
- Confirm you have multiple network files (at least 2) per method
- Check that file naming is consistent across replicates
- Verify directory structure matches expected format


#### Method Performance Issues

##### Poor ROC/PR Performance
- Check if TF-target pairs in your network overlap with ground truth
- Verify scoring 

##### Unexpected Ranking Results  
- Examine the `detailed_results.tsv` file for edge statistics
- Check TP, FP, FN counts to understand performance bottlenecks
- Compare network sizes - very small networks may have inflated precision

##### Stability Analysis Problems
- Check that file reading parameters are consistent across files
- Verify that custom column mappings are correct

### Debug Workflows

#### Debug New Method Integration
```r
# Step 1: Check file structure
cat("=== File Check ===\n")
your_file <- "path/to/your/method.tsv"
if (file.exists(your_file)) {
  df <- read_delim(your_file, delim = "\t", n_max = 5)
  print(colnames(df))
  print(df)
} else {
  cat("File not found!\n")
}

# Step 2: Test column mapping
tf_col <- "your_tf_column"
target_col <- "your_target_column"  
score_col <- "your_score_column"

if (all(c(tf_col, target_col, score_col) %in% colnames(df))) {
  cat("All columns found!\n")
} else {
  cat("Missing columns. Available:", paste(colnames(df), collapse = ", "), "\n")
}

# Step 3: Check data quality
print(summary(df[[score_col]]))  # Should be numeric
print(table(is.na(df[[score_col]])))  # Check for NAs


#### Missing Columns
Check that your method's column names match the expected format in the method-specific table above.

#### Empty Results  
Ensure ground truth files exist and contain the expected `Source` and `Target` columns.

