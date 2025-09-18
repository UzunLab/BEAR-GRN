# Global variables for this function
utils::globalVariables(c(
  "target", "source_upper", "target_upper", "Metric", "Method", 
  "Dataset", "Score", "tf", "pair", "label", "precision", "recall", "f1_score",
  "F1-score", "Precision", "Recall", "F1_Rank", "Precision_Rank", "Recall_Rank",
  "Sample", "Median_JI", "pivot_wider"
))
#' Reproduce Early GRN Metrics with Lollipop Plots
#'
#' Computes precision, recall, and F1-score metrics for GRN methods using filtered networks
#' (top 10K edges) and generates lollipop plots for visualization.
#'
#' @param input_dir Path to input directory containing dataset subdirectories (e.g., "INFERRED.GRNS")
#' @param output_dir Path to output directory for results and plots (default: "filtered_grn_metrics")
#' @param ground_truth_dir Path to directory containing ground truth files (default: "GROUND.TRUTHS")
#' @param max_edges Maximum number of top edges to consider (default: 10000)
#' @param method_colors Optional named vector of colors for methods. If NULL, uses default colors.
#' @return List containing summary data and plot objects
#' @export
#' @importFrom readr read_tsv write_csv
#' @importFrom dplyr filter select mutate arrange distinct bind_rows left_join group_by summarise ungroup
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs theme_minimal theme element_text coord_flip facet_grid scale_color_manual ggsave
#' @importFrom stringr str_detect str_extract str_replace str_to_upper str_trim
#' @importFrom tibble tibble
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pretty_breaks percent
#' @importFrom grDevices pdf dev.off
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom tidyr pivot_wider pivot_longer
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- reproduce_early_metrics(
#'   input_dir = "INFERRED.GRNS",
#'   output_dir = "filtered_metrics_results"
#' )
#' 
#' # With custom parameters
#' results <- reproduce_early_metrics(
#'   input_dir = "my_data",
#'   ground_truth_dir = "ground_truth_dir",
#'   output_dir = "my_results",
#'   max_edges = 10000
#' )
#' }
reproduce_early_metrics <- function(input_dir, 
                                   output_dir = "filtered_grn_metrics",
                                   ground_truth_dir = "GROUND.TRUTHS",
                                   max_edges = 10000,
                                   method_colors = NULL) {
  
  
  cat("=== Compute Filtered GRN Metrics ===\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Max edges per method:", max_edges, "\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Method-specific column mappings
  column_mappings <- list(
    "CellOracle" = c("source", "target"),
    "FigR" = c("Motif", "DORC"),
    "GRaNIE" = c("TF.name", "gene.name"),
    "LINGER" = c("Source", "Target"),
    "Pando" = c("tf", "target"),
    "SCENIC+" = c("Source", "Target"),
    "Scenic+" = c("Source", "Target"),  # Alternative spelling
    "TRIPOD" = c("TF", "gene"),
    "DIRECTNET" = c("TF", "Gene")
  )
  
  score_mappings <- list(
    "CellOracle" = "coef_abs",
    "FigR" = "Score",
    "GRaNIE" = "TF_gene.r",
    "LINGER" = "Score",
    "Pando" = "estimate",
    "SCENIC+" = "Score",
    "Scenic+" = "Score",  # Alternative spelling
    "TRIPOD" = "abs_coef",
    "DIRECTNET" = "Score"
  )
  
  # Default color palette for methods
  if (is.null(method_colors)) {
    method_colors <- c(
      "CellOracle" = "#1B9E77",  # Teal green
      "FigR"       = "#66A61E",  # Green
      "GRaNIE"     = "#A6761D",  # Brown
      "LINGER"     = "#E7298A",  # Pink/magenta
      "Pando"      = "#7570B3",  # Purple
      "SCENIC+"    = "#D95F02",  # Orange
      "Scenic+"    = "#D95F02",  # Orange (alternative spelling)
      "TRIPOD"     = "#E6AB02",  # Mustard yellow
      "DIRECTNET"  = "#666666"   # Grey
    )
  }
  
  # Dataset to ground truth mapping
  datasets_mapping <- list(
    "K562" = "filtered_RN117_K562.tsv",
    "Macrophage_S1" = "filtered_RN204_Buffer1.tsv",
    "Macrophage_S2" = "filtered_RN204_Buffer1.tsv",
    "mESC_E7.5_rep1" = "filtered_RN111_E7.5_rep1.tsv",
    "mESC_E7.5_rep2" = "filtered_RN111_E7.5_rep2.tsv",
    "mESC_E8.5_rep1" = "filtered_RN111_E8.5_rep1.tsv",
    "mESC_E8.5_rep2" = "filtered_RN111_E8.5_rep2.tsv",
    "E7.5_rep1" = "filtered_RN111_E7.5_rep1.tsv",
    "E7.5_rep2" = "filtered_RN111_E7.5_rep2.tsv",
    "E8.5_rep1" = "filtered_RN111_E8.5_rep1.tsv",
    "E8.5_rep2" = "filtered_RN111_E8.5_rep2.tsv"
  )
  
  # Auto-detect datasets from input directory
  dataset_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  dataset_dirs <- dataset_dirs[dataset_dirs != ""]
  
  if (length(dataset_dirs) == 0) {
    stop("No dataset directories found in input_dir: ", input_dir)
  }
  
  cat("Found datasets:", paste(dataset_dirs, collapse = ", "), "\n")
  
  # Initialize storage
  plot_data <- tibble(Metric = character(), Score = numeric(), Dataset = character(), Method = character())
  all_summaries <- list()
  failed_datasets <- c()
  
  # Main analysis loop
  for (dataset in dataset_dirs) {
    cat("\nProcessing dataset:", dataset, "\n")
    
    # Determine ground truth file
    ground_truth_file <- if (dataset %in% names(datasets_mapping)) {
      file.path(ground_truth_dir, datasets_mapping[[dataset]])
    } else {
      # Try to find a matching ground truth file
      gt_files <- list.files(ground_truth_dir, 
                            pattern = paste0(".*", dataset, ".*\\.tsv$"), 
                            full.names = TRUE)
      if (length(gt_files) == 0) {
        # Fallback: try any .tsv file
        gt_files <- list.files(ground_truth_dir, pattern = "\\.tsv$", full.names = TRUE)
        if (length(gt_files) > 0) gt_files[1] else NULL
      } else {
        gt_files[1]
      }
    }
    
    if (is.null(ground_truth_file) || !file.exists(ground_truth_file)) {
      cat("  Warning: Ground truth file not found for", dataset, "\n")
      failed_datasets <- c(failed_datasets, dataset)
      next
    }
    
    # Load ground truth
    gt <- tryCatch({
      read_tsv(ground_truth_file, col_types = cols(), show_col_types = FALSE)
    }, error = function(e) {
      cat("  Error reading ground truth file:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(gt)) {
      failed_datasets <- c(failed_datasets, dataset)
      next
    }
    
    # Process ground truth
    gt_edges <- unique(tibble(
      source = toupper(gt$Source),
      target = toupper(gt$Target)
    ))
    gt_pairs <- with(gt_edges, paste(source, target, sep = "_"))
    gt_tfs <- unique(gt_edges$source)
    gt_targets <- unique(gt_edges$target)
    
    cat("  Ground truth:", length(gt_pairs), "edges,", length(gt_tfs), "TFs,", length(gt_targets), "targets\n")
    
    # Get available methods in dataset directory
    dataset_input_dir <- file.path(input_dir, dataset)
    method_dirs <- list.dirs(dataset_input_dir, full.names = FALSE, recursive = FALSE)
    methods <- method_dirs[method_dirs != ""]
    
    if (length(methods) == 0) {
      cat("  No method directories found in", dataset_input_dir, "\n")
      failed_datasets <- c(failed_datasets, dataset)
      next
    }
    
    summary_list <- list()
    
    # Process each method
    for (method in methods) {
      if (!(method %in% names(column_mappings))) {
        cat("    ", method, ": Not in column mappings, skipping\n")
        next
      }
      
      method_dir <- file.path(dataset_input_dir, method)
      files <- list.files(method_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)
      
      if (length(files) == 0) {
        cat("    ", method, ": No data files found, skipping\n")
        next
      }
      
      file <- files[1]
      cols <- column_mappings[[method]]
      score_col <- score_mappings[[method]]
      
      # Read and process inferred network
      df <- tryCatch({
        # Determine separator
        sep <- if (grepl("\\.csv$", file)) "," else "\t"
        
        temp_df <- read_delim(file, delim = sep, col_types = cols(), show_col_types = FALSE)
        
        # Rename columns to standard names
        if (length(cols) >= 2 && all(cols %in% names(temp_df))) {
          temp_df <- temp_df %>% rename(source = all_of(cols[1]), target = all_of(cols[2]))
        } else {
          stop("Required columns not found in file")
        }
        
        # Handle scoring column
        if (!is.null(score_col) && score_col %in% names(temp_df)) {
          temp_df[[score_col]] <- abs(temp_df[[score_col]])
          temp_df <- temp_df %>% arrange(desc(.data[[score_col]]))
        } else if (!is.null(score_col)) {
          cat("    Warning: Score column '", score_col, "' not found for ", method, "\n")
        }
        
        temp_df
      }, error = function(e) {
        cat("    ", method, ": Error reading file -", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(df)) next
      
      cat("    ", method, ": Original network size:", nrow(df), "edges\n")
      
      # Filter by ground truth TFs and targets
      df_filtered <- df %>%
        mutate(
          source_upper = toupper(source),
          target_upper = toupper(target)
        ) %>%
        filter(
          source_upper %in% gt_tfs,
          target_upper %in% gt_targets
        )
      
      cat("    ", method, ": Filtered network size:", nrow(df_filtered), "edges\n")
      
      # Take top max_edges edges (or all if fewer than max_edges)
      n_edges <- min(max_edges, nrow(df_filtered))
      if (n_edges == 0) {
        cat("    ", method, ": No edges after filtering, skipping\n")
        next
      }
      
      df_final <- df_filtered %>% head(n_edges)
      cat("    ", method, ": Final network size:", nrow(df_final), "edges\n")
      
      # Create edge pairs for comparison
      df_pairs <- df_final %>%
        transmute(pair = paste(source_upper, target_upper, sep = "_"))
      
      # Calculate metrics
      tp <- intersect(df_pairs$pair, gt_pairs)
      fp <- setdiff(df_pairs$pair, gt_pairs)
      fn <- setdiff(gt_pairs, df_pairs$pair)
      
      precision <- ifelse((length(tp) + length(fp)) > 0, length(tp) / (length(tp) + length(fp)), 0)
      recall <- ifelse((length(tp) + length(fn)) > 0, length(tp) / (length(tp) + length(fn)), 0)
      f1 <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
      
      cat("    ", method, ": Precision =", round(precision, 4), ", Recall =", round(recall, 4), ", F1 =", round(f1, 4), "\n")
      
      # Add to plot data
      plot_data <- bind_rows(plot_data, tibble(
        Metric = c("Precision", "Recall", "F1-score"),
        Score = c(precision, recall, f1),
        Dataset = dataset,
        Method = method
      ))
      
      # Get unique TFs and targets from filtered network
      inferred_tfs <- unique(df_final$source_upper)
      inferred_targets <- unique(df_final$target_upper)
      
      # Create summary statistics
      summary_list[[method]] <- tibble(
        Method = method,
        GT_total_edges = length(gt_pairs),
        GT_TFs = length(gt_tfs),
        GT_targets = length(gt_targets),
        Original_network_size = nrow(df),
        Filtered_network_size = nrow(df_filtered),
        Final_network_size = nrow(df_final),
        Inferred_TFs = length(inferred_tfs),
        Inferred_targets = length(inferred_targets),
        TP = length(tp),
        FP = length(fp),
        FN = length(fn),
        Precision = round(precision, 4),
        Recall = round(recall, 4),
        F1_score = round(f1, 4),
        Common_TFs = length(intersect(gt_tfs, inferred_tfs)),
        Common_targets = length(intersect(gt_targets, inferred_targets))
      )
    }
    
    # Save summary statistics for this dataset
    if (length(summary_list) > 0) {
      df_summary <- bind_rows(summary_list) %>% mutate(Dataset = dataset)
      write_tsv(df_summary, file.path(output_dir, paste0(dataset, "_filtered_metrics.tsv")))
      all_summaries[[dataset]] <- df_summary
      cat("  Saved summary to:", file.path(output_dir, paste0(dataset, "_filtered_metrics.tsv")), "\n")
    }
  }
  
  # Check if we have any data to plot
  if (nrow(plot_data) == 0) {
    cat("No data available for plotting.\n")
    return(list(
      plot_data = plot_data,
      summaries = all_summaries,
      failed_datasets = failed_datasets,
      plots = list()
    ))
  }
  
  # Clean and prepare plot data
  plot_data$Dataset <- str_trim(plot_data$Dataset)
  
  # Set dataset order
  #dataset_order <- unique(plot_data$Dataset)
  #plot_data$Dataset <- factor(plot_data$Dataset, levels = dataset_order)
  # Set custom dataset order
  desired_order <- c("Macrophage_S1", "Macrophage_S2", "K562", "E7.5_rep1", "E7.5_rep2", "E8.5_rep1", "E8.5_rep2")
  available_datasets <- unique(plot_data$Dataset)
  dataset_order <- intersect(desired_order, available_datasets)
  plot_data$Dataset <- factor(plot_data$Dataset, levels = dataset_order)
  
  # Define method order (only include methods that appear in the data)
  all_method_order <- c("CellOracle", "FigR", "GRaNIE", "LINGER", "Pando", "SCENIC+", "Scenic+", "TRIPOD", "DIRECTNET")
  available_methods <- intersect(all_method_order, unique(plot_data$Method))
  plot_data$Method <- factor(plot_data$Method, levels = available_methods)
  
  # Create lollipop plots for each metric
  plots <- list()
  
  for (metric in c("Precision", "Recall", "F1-score")) {
    cat("Creating lollipop plot for:", metric, "\n")
    
    df_metric <- plot_data %>% 
      filter(Metric == metric) %>%
      filter(Method %in% available_methods)  # Only include available methods
    
    if (nrow(df_metric) == 0) {
      cat("  No data available for", metric, "\n")
      next
    }
    
    # Compute dynamic y-axis limit with 10% buffer, capped at 1
    max_score <- max(df_metric$Score, na.rm = TRUE)
    y_limit <- min(1, max(0.1, round(max_score * 1.1, 2)))
    
    p_lollipop <- ggplot(df_metric, aes(x = Dataset, y = Score)) +
      geom_segment(
        aes(x = Dataset, xend = Dataset, y = 0, yend = Score),
        color = "white",
        position = position_dodge(width = 0.8),
        size = 0.8,
        alpha = 0.7
      ) +
      geom_point(
        aes(color = Method),
        shape = 16,
        position = position_dodge(width = 0.8),
        size = 4.5,
        alpha = 0.9
      ) +
      scale_color_manual(values = method_colors) +
      scale_y_continuous(
        limits = c(0, y_limit),
        expand = expansion(mult = c(0, 0.02)),
        labels = scales::number_format(accuracy = 0.01)
      ) +
      scale_x_discrete(expand = expansion(add = c(1.5, 1.5))) +
      labs(
        title = paste0(metric, " Comparison (Filtered Networks, Top ", format(max_edges, scientific = FALSE), " Edges)"),
        subtitle = "Networks filtered by ground truth TFs and targets",
        x = "Dataset", 
        y = metric,
        color = "Method"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        axis.line.x = element_line(color = "black", size = 0.6),
        axis.line.y = element_line(color = "black", size = 0.6),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    plots[[metric]] <- p_lollipop
    
    # Save plot
    filename <- file.path(output_dir, paste0(gsub("-", "_", metric), "_filtered_lollipop_plot.pdf"))
    ggsave(
      filename = filename,
      plot = p_lollipop,
      width = 12,
      height = 7,
      device = "pdf",
      dpi = 300
    )
    
    # Also save PNG version
    filename_png <- file.path(output_dir, paste0(gsub("-", "_", metric), "_filtered_lollipop_plot.png"))
    ggsave(
      filename = filename_png,
      plot = p_lollipop,
      width = 12,
      height = 7,
      device = "png",
      dpi = 300
    )
    
    cat("  Saved plot to:", filename, "\n")
  }
  
  # Create combined summary across all datasets
  if (length(all_summaries) > 0) {
    combined_summary <- bind_rows(all_summaries)
    write_tsv(combined_summary, file.path(output_dir, "all_datasets_filtered_metrics.tsv"))
    cat("Combined summary saved to:", file.path(output_dir, "all_datasets_filtered_metrics.tsv"), "\n")
  }
  
  # Final summary
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("Total datasets found:", length(dataset_dirs), "\n")
  cat("Successfully processed:", length(all_summaries), "\n")
  cat("Failed:", length(failed_datasets), "\n")
  
  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  cat("Analysis completed! Check the '", output_dir, "' folder for results.\n")
  
  return(list(
    plot_data = plot_data,
    summaries = all_summaries,
    combined_summary = if (exists("combined_summary")) combined_summary else NULL,
    failed_datasets = failed_datasets,
    plots = plots
  ))
}
