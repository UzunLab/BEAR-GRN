utils::globalVariables(c(
  "target", "source_upper", "target_upper", "Metric", "Method", 
  "Dataset", "Score", "tf", "pair", "label", "precision", "recall", "f1_score",
  "F1-score", "Precision", "Recall", "F1_Rank", "Precision_Rank", "Recall_Rank",
  "Sample", "Median_JI", "pivot_wider"
))
#' Benchmark New Method on early metrics GRN Performance
#'
#' Evaluates a new method's GRN against ground truth and compares it with existing methods
#' from a specific dataset. Generates comparison plots and detailed metrics.
#'
#' @param new_grn_file Path to the new method's GRN file (CSV or TSV)
#' @param dataset_name Name of the dataset to compare against (e.g., "K562", "Macrophage_S1", "E7.5_rep1")
#' @param tf_column Name of the TF/source column in the new GRN file
#' @param target_column Name of the target column in the new GRN file
#' @param score_column Name of the score column in the new GRN file (optional, can be NULL)
#' @param method_name Name for the new method (for plotting)
#' @param input_dir Path to input directory containing existing methods (default: "INFERRED.GRNS")
#' @param ground_truth_dir Path to directory containing ground truth files (default: "GROUND.TRUTHS")
#' @param output_dir Path to output directory for results (default: "new_method_evaluation")
#' @param max_edges Maximum number of top edges to consider (default: 10000)
#' @param method_color Color for the new method in plots (default: "#FF0000" - red)
#' @return List containing evaluation results, comparison data, and plots
#' @export
#' @importFrom readr read_tsv write_csv
#' @importFrom dplyr filter select mutate arrange distinct bind_rows left_join group_by summarise ungroup
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal theme element_text
#' @importFrom stringr str_detect str_extract str_replace str_to_upper str_trim
#' @importFrom tibble tibble
#' @importFrom stats setNames
#' @importFrom utils head
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- benchmark_new_method_early_metrics(
#'   new_grn_file = "my_method_K562.csv",
#'   dataset_name = "K562",
#'   tf_column = "TF",
#'   target_column = "Gene",
#'   score_column = "Weight",
#'   method_name = "MyMethod"
#' )
#' 
#' # Without score column
#' results <- benchmark_new_method_early_metrics(
#'   new_grn_file = "my_method_results.tsv",
#'   dataset_name = "Macrophage_S1",
#'   tf_column = "source",
#'   target_column = "target",
#'   score_column = NULL,
#'   method_name = "NewApproach"
#' )
#' }
benchmark_new_method_early_metrics <- function(new_grn_file,
                               dataset_name,
                               tf_column,
                               target_column,
                               score_column = NULL,
                               method_name = "NewMethod",
                               input_dir = "INFERRED.GRNS",
                               ground_truth_dir = "GROUND.TRUTHS",
                               output_dir = "new_method_evaluation",
                               max_edges = 10000,
                               method_color = "#FF0000") {
  

  
  cat("=== Evaluating New Method GRN ===\n")
  cat("New GRN file:", new_grn_file, "\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Method name:", method_name, "\n")
  cat("TF column:", tf_column, "\n")
  cat("Target column:", target_column, "\n")
  cat("Score column:", ifelse(is.null(score_column), "None", score_column), "\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Method-specific column mappings for existing methods
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
  
  # Default color palette for existing methods
  existing_method_colors <- c(
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
  
  # Add new method color
  all_method_colors <- c(existing_method_colors, setNames(method_color, method_name))
  
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
  
  # Function to evaluate a single method
  evaluate_single_method <- function(df, method_name, gt_pairs, gt_tfs, gt_targets, max_edges) {
    cat("    ", method_name, ": Original network size:", nrow(df), "edges\n")
    
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
    
    cat("    ", method_name, ": Filtered network size:", nrow(df_filtered), "edges\n")
    
    # Take top max_edges edges (or all if fewer than max_edges)
    n_edges <- min(max_edges, nrow(df_filtered))
    if (n_edges == 0) {
      cat("    ", method_name, ": No edges after filtering, skipping\n")
      return(NULL)
    }
    
    df_final <- df_filtered %>% head(n_edges)
    cat("    ", method_name, ": Final network size:", nrow(df_final), "edges\n")
    
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
    
    cat("    ", method_name, ": Precision =", round(precision, 4), ", Recall =", round(recall, 4), ", F1 =", round(f1, 4), "\n")
    
    # Get unique TFs and targets from filtered network
    inferred_tfs <- unique(df_final$source_upper)
    inferred_targets <- unique(df_final$target_upper)
    
    return(list(
      metrics = c(precision = precision, recall = recall, f1 = f1),
      summary = tibble(
        Method = method_name,
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
    ))
  }
  
  # Check if new GRN file exists
  if (!file.exists(new_grn_file)) {
    stop("New GRN file not found: ", new_grn_file)
  }
  
  # Load and process new method's GRN
  cat("\nLoading new method GRN...\n")
  new_grn <- tryCatch({
    # Determine separator
    sep <- if (grepl("\\.csv$", new_grn_file)) "," else "\t"
    
    temp_df <- read_delim(new_grn_file, delim = sep, col_types = cols(), show_col_types = FALSE)
    
    # Check if required columns exist
    if (!(tf_column %in% names(temp_df))) {
      stop("TF column '", tf_column, "' not found in file. Available columns: ", paste(names(temp_df), collapse = ", "))
    }
    if (!(target_column %in% names(temp_df))) {
      stop("Target column '", target_column, "' not found in file. Available columns: ", paste(names(temp_df), collapse = ", "))
    }
    
    # Rename columns to standard names
    temp_df <- temp_df %>% rename(source = all_of(tf_column), target = all_of(target_column))
    
    # Handle scoring column if provided
    if (!is.null(score_column)) {
      if (score_column %in% names(temp_df)) {
        temp_df[[score_column]] <- abs(temp_df[[score_column]])
        temp_df <- temp_df %>% arrange(desc(.data[[score_column]]))
        cat("  Sorted by score column:", score_column, "\n")
      } else {
        cat("  Warning: Score column '", score_column, "' not found. Proceeding without sorting.\n")
      }
    }
    
    temp_df
  }, error = function(e) {
    stop("Error reading new GRN file: ", e$message)
  })
  
  # Determine ground truth file
  ground_truth_file <- if (dataset_name %in% names(datasets_mapping)) {
    file.path(ground_truth_dir, datasets_mapping[[dataset_name]])
  } else {
    # Try to find a matching ground truth file
    gt_files <- list.files(ground_truth_dir, 
                          pattern = paste0(".*", dataset_name, ".*\\.tsv$"), 
                          full.names = TRUE)
    if (length(gt_files) == 0) {
      stop("Ground truth file not found for dataset: ", dataset_name)
    } else {
      gt_files[1]
    }
  }
  
  if (!file.exists(ground_truth_file)) {
    stop("Ground truth file not found: ", ground_truth_file)
  }
  
  # Load ground truth
  cat("Loading ground truth from:", ground_truth_file, "\n")
  gt <- tryCatch({
    read_tsv(ground_truth_file, col_types = cols(), show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading ground truth file: ", e$message)
  })
  
  # Process ground truth
  gt_edges <- unique(tibble(
    source = toupper(gt$Source),
    target = toupper(gt$Target)
  ))
  gt_pairs <- with(gt_edges, paste(source, target, sep = "_"))
  gt_tfs <- unique(gt_edges$source)
  gt_targets <- unique(gt_edges$target)
  
  cat("Ground truth:", length(gt_pairs), "edges,", length(gt_tfs), "TFs,", length(gt_targets), "targets\n")
  
  # Evaluate new method
  cat("\nEvaluating new method...\n")
  new_method_result <- evaluate_single_method(new_grn, method_name, gt_pairs, gt_tfs, gt_targets, max_edges)
  
  if (is.null(new_method_result)) {
    stop("New method evaluation failed - no valid edges after filtering")
  }
  
  # DEBUG: Check new method result
  cat("\n=== DEBUG INFO ===\n")
  cat("Method name:", method_name, "\n")
  cat("New method result summary:\n")
  print(new_method_result$summary)
  cat("New method metrics:\n")
  print(new_method_result$metrics)
  
  # Load and evaluate existing methods from the specified dataset
  cat("\nLoading existing methods for comparison...\n")
  dataset_input_dir <- file.path(input_dir, dataset_name)
  
  if (!dir.exists(dataset_input_dir)) {
    cat("Warning: No existing methods found for dataset", dataset_name, "in", input_dir, "\n")
    existing_results <- list()
  } else {
    method_dirs <- list.dirs(dataset_input_dir, full.names = FALSE, recursive = FALSE)
    methods <- method_dirs[method_dirs != ""]
    
    existing_results <- list()
    
    for (method in methods) {
      if (!(method %in% names(column_mappings))) {
        cat("  ", method, ": Not in column mappings, skipping\n")
        next
      }
      
      method_dir <- file.path(dataset_input_dir, method)
      files <- list.files(method_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)
      
      if (length(files) == 0) {
        cat("  ", method, ": No data files found, skipping\n")
        next
      }
      
      file <- files[1]
      cols <- column_mappings[[method]]
      score_col <- score_mappings[[method]]
      
      # Read and process existing method's network
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
        }
        
        temp_df
      }, error = function(e) {
        cat("  ", method, ": Error reading file -", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(df)) next
      
      result <- evaluate_single_method(df, method, gt_pairs, gt_tfs, gt_targets, max_edges)
      if (!is.null(result)) {
        existing_results[[method]] <- result
      }
    }
  }
  
  # Combine all results for comparison - FIXED VERSION
  all_results <- existing_results
  
  # Ensure new method is added properly
  if (!is.null(new_method_result) && !is.null(new_method_result$summary)) {
    all_results[[method_name]] <- new_method_result
    cat("Successfully added", method_name, "to results\n")
  } else {
    stop("Failed to add new method to results")
  }
  
  # DEBUG: Check all results
  cat("\n=== ALL RESULTS DEBUG ===\n")
  cat("All results names:", paste(names(all_results), collapse = ", "), "\n")
  cat("Number of results:", length(all_results), "\n")
  
  # Create comparison data
  plot_data <- tibble()
  summary_data <- tibble()
  
  for (method_name_loop in names(all_results)) {
    result <- all_results[[method_name_loop]]
    
    if (!is.null(result) && !is.null(result$metrics) && !is.null(result$summary)) {
      # Add to plot data
      plot_data <- bind_rows(plot_data, tibble(
        Metric = c("Precision", "Recall", "F1-score"),
        Score = c(result$metrics["precision"], result$metrics["recall"], result$metrics["f1"]),
        Dataset = dataset_name,
        Method = method_name_loop
      ))
      
      # Add to summary data
      summary_data <- bind_rows(summary_data, result$summary)
    }
  }
  
  # Remove any rows with NA method names or scores
  plot_data <- plot_data %>% 
    filter(!is.na(Method), !is.na(Score)) %>%
    filter(Method != "")
  
  summary_data <- summary_data %>% 
    filter(!is.na(Method)) %>%
    filter(Method != "")
  
  # DEBUG: Check plot data
  cat("\n=== PLOT DATA DEBUG ===\n")
  cat("Unique methods in plot_data:", paste(unique(plot_data$Method), collapse = ", "), "\n")
  cat("Number of rows in plot_data:", nrow(plot_data), "\n")
  cat("New method data in plot_data:\n")
  print(plot_data %>% filter(Method == method_name))
  
  # Set method order (new method first, then existing methods)
  all_method_order <- c(method_name, names(existing_method_colors))
  available_methods <- intersect(all_method_order, unique(plot_data$Method))
  plot_data$Method <- factor(plot_data$Method, levels = available_methods)
  
  # Create comparison plots
  plots <- list()
  
  for (metric in c("Precision", "Recall", "F1-score")) {
    cat("Creating comparison plot for:", metric, "\n")
    
    df_metric <- plot_data %>% 
      filter(Metric == metric) %>%
      filter(Method %in% available_methods)
    
    if (nrow(df_metric) == 0) {
      cat("  No data available for", metric, "\n")
      next
    }
    
    # Compute dynamic y-axis limit
    max_score <- max(df_metric$Score, na.rm = TRUE)
    y_limit <- min(1, max(0.1, round(max_score * 1.1, 2)))
    
    p_comparison <- ggplot(df_metric, aes(x = Method, y = Score)) +
      geom_col(
        aes(fill = Method),
        width = 0.7,
        alpha = 0.8
      ) +
      scale_fill_manual(values = all_method_colors) +
      scale_y_continuous(
        limits = c(0, y_limit),
        expand = expansion(mult = c(0, 0.02)),
        labels = scales::number_format(accuracy = 0.01)
      ) +
      geom_text(
        aes(label = round(Score, 3)),
        vjust = -0.5,
        size = 3.5,
        fontface = "bold"
      ) +
      labs(
        title = paste0(metric, " Comparison - ", dataset_name, " Dataset"),
        subtitle = paste0("New method (", method_name, ") vs Existing methods (Top ", format(max_edges, scientific = FALSE), " edges)"),
        x = "Method", 
        y = metric,
        fill = "Method"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",  # Remove legend since colors are self-explanatory
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
    
    plots[[metric]] <- p_comparison
    
    # Save plot
    filename <- file.path(output_dir, paste0(dataset_name, "_", gsub("-", "_", metric), "_comparison.pdf"))
    ggsave(
      filename = filename,
      plot = p_comparison,
      width = 10,
      height = 7,
      device = "pdf",
      dpi = 300
    )
    
    # Also save PNG version
    filename_png <- file.path(output_dir, paste0(dataset_name, "_", gsub("-", "_", metric), "_comparison.png"))
    ggsave(
      filename = filename_png,
      plot = p_comparison,
      width = 10,
      height = 7,
      device = "png",
      dpi = 300
    )
    
    cat("  Saved plot to:", filename, "\n")
  }
  
  # Create ranking table
  ranking_data <- plot_data %>%
    filter(!is.na(Score) & !is.na(Method)) %>%  # Remove NA values
    select(Method, Metric, Score) %>%
    pivot_wider(names_from = Metric, values_from = Score) %>%
    filter(!is.na(`F1-score`)) %>%  # Ensure F1-score is not NA
    arrange(desc(`F1-score`)) %>%
    mutate(
      F1_Rank = rank(-`F1-score`, ties.method = "min"),
      Precision_Rank = rank(-Precision, ties.method = "min"),
      Recall_Rank = rank(-Recall, ties.method = "min")
    ) %>%
    select(Method, `F1-score`, F1_Rank, Precision, Precision_Rank, Recall, Recall_Rank)
  
  # Add highlighting for new method
  if (method_name %in% ranking_data$Method) {
    new_method_rank <- ranking_data %>% filter(Method == method_name) %>% pull(F1_Rank)
    cat("\nRanking Results:\n")
    cat("================\n")
    cat(sprintf("Your method '%s' ranks #%d out of %d methods based on F1-score\n", 
                method_name, new_method_rank, nrow(ranking_data)))
  } else {
    cat("\nWarning: New method not found in ranking data\n")
    cat("Available methods:", paste(ranking_data$Method, collapse = ", "), "\n")
  }
  print(ranking_data)
  
  # Save all results
  write_tsv(summary_data, file.path(output_dir, paste0(dataset_name, "_", method_name, "_detailed_results.tsv")))
  write_tsv(ranking_data, file.path(output_dir, paste0(dataset_name, "_", method_name, "_ranking.tsv")))
  write_tsv(plot_data, file.path(output_dir, paste0(dataset_name, "_", method_name, "_plot_data.tsv")))
  
  cat("\nAnalysis completed!\n")
  cat("Results saved to:", output_dir, "\n")
  cat("- Detailed results:", paste0(dataset_name, "_", method_name, "_detailed_results.tsv"), "\n")
  cat("- Method ranking:", paste0(dataset_name, "_", method_name, "_ranking.tsv"), "\n")
  cat("- Plot data:", paste0(dataset_name, "_", method_name, "_plot_data.tsv"), "\n")
  
  return(list(
    new_method_metrics = new_method_result$metrics,
    all_summaries = summary_data,
    ranking = ranking_data,
    plot_data = plot_data,
    plots = plots,
    dataset = dataset_name,
    method_name = method_name,
    ground_truth_stats = list(
      total_edges = length(gt_pairs),
      total_tfs = length(gt_tfs),
      total_targets = length(gt_targets)
    )
  ))
}