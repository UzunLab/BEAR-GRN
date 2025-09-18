# Global variable declarations to satisfy R CMD check
utils::globalVariables(c(
  # Column names used in dplyr operations
  "method_curve", "legend_label", "curve_type", "Source", "Target", 
  "score", "tf", "target", "pair", "label", "method", "AUROC", "AUPRC",
  "fpr", "tpr", "AUPRC_random", "source_upper", "target_upper",
  "Metric", "Method", "Dataset", "Score"
))
#' Benchmark New GRN Method Against Existing Methods
#'
#' Allows users to benchmark their new GRN method against existing methods by providing
#' their inferred GRN file along with method specifications. The function will compare
#' the new method with existing methods using ROC and PR curve analysis.
#'
#' @param dataset_name Name of the dataset for which the GRN was constructed
#' @param new_grn_file Path to the inferred GRN file for the new method
#' @param new_method_name Name of the new method to be benchmarked
#' @param tf_column Name of the transcription factor column in the new GRN file
#' @param target_column Name of the target gene column in the new GRN file  
#' @param score_column Name of the score/confidence column in the new GRN file
#' @param input_dir Path to input directory containing existing method results (optional)
#' @param output_dir Path to output directory for benchmark results
#' @param ground_truth_dir Path to directory containing ground truth files (default: "GROUND.TRUTHS")
#' @param use_filtered_approach Logical, whether to filter to tested TF-target space (default: TRUE)
#' @param existing_method_info Optional data frame with existing method information. If NULL, uses default mapping.
#' @return List containing benchmark results and plots
#' @export
#' @importFrom readr read_tsv write_csv
#' @importFrom dplyr filter select mutate arrange distinct bind_rows left_join tribble
#' @importFrom pROC roc auc
#' @importFrom PRROC pr.curve
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal ggsave
#' @importFrom tibble tibble
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotrix gap.plot
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par lines legend
#' @importFrom stats runif setNames
#' @importFrom utils head
#' @importFrom stringr str_detect str_extract str_ends str_starts str_replace str_to_upper str_trim
#' @examples
#' \dontrun{
#' # Benchmark a new method
#' results <- benchmark_new_method(
#'   dataset_name = "K562",
#'   new_grn_file = "my_method_results.tsv",
#'   new_method_name = "MyNewMethod",
#'   tf_column = "TF",
#'   target_column = "Gene",
#'   score_column = "Confidence",
#'   output_dir = "Benchmark_Results"
#' )
#' 
#' # Benchmark against specific existing methods
#' results <- benchmark_new_method(
#'   dataset_name = "mESC_E7.5_rep1",
#'   new_grn_file = "my_grn_inference.csv",
#'   new_method_name = "NovelGRN",
#'   tf_column = "source_tf",
#'   target_column = "target_gene", 
#'   score_column = "edge_weight",
#'   input_dir = "INFERRED.GRNS",
#'   output_dir = "my_benchmark"
#' )
#' }
benchmark_new_method <- function(dataset_name,
                                new_grn_file,
                                new_method_name,
                                tf_column,
                                target_column,
                                score_column,
                                input_dir = NULL,
                                output_dir,
                                ground_truth_dir = "GROUND.TRUTHS",
                                use_filtered_approach = TRUE,
                                existing_method_info = NULL) {
  
  
  cat("=== Benchmark New Method Against Existing Methods ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("New method:", new_method_name, "\n")
  cat("New GRN file:", new_grn_file, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Validate new GRN file
  if (!file.exists(new_grn_file)) {
    stop("New GRN file not found: ", new_grn_file)
  }
  
  # Default existing method information
  if (is.null(existing_method_info)) {
    existing_method_info <- tribble(
      ~method,      ~tf_col,     ~target_col, ~score_col,
      "CellOracle", "source",    "target",    "coef_mean",
      "SCENIC+",    "Source",    "Target",    "Score",
      "Pando",      "tf",        "target",    "estimate",
      "LINGER",     "Source",    "Target",    "Score",
      "FigR",       "Motif",     "DORC",      "Score",
      "TRIPOD",     "TF",        "gene",      "abs_coef",
      "GRaNIE",     "TF.name",   "gene.name", "TF_gene.r"
    )
  }
  
  # Add new method to method info
  new_method_row <- data.frame(
    method = new_method_name,
    tf_col = tf_column,
    target_col = target_column,
    score_col = score_column,
    stringsAsFactors = FALSE
  )
  
  method_info <- bind_rows(new_method_row, existing_method_info)
  
  # Define dataset to ground truth mapping
  datasets <- list(
    "K562" = "filtered_RN117_K562.tsv",
    "Macrophage_S1" = "filtered_RN204_Buffer1.tsv",
    "Macrophage_S2" = "filtered_RN204_Buffer1.tsv",
    "mESC_E7.5_rep1" = "filtered_RN111_E7.5_rep1.tsv",
    "mESC_E7.5_rep2" = "filtered_RN111_E7.5_rep2.tsv",
    "mESC_E8.5_rep1" = "filtered_RN111_E8.5_rep1.tsv",
    "mESC_E8.5_rep2" = "filtered_RN111_E8.5_rep2.tsv"
  )
  
  # Determine ground truth file
  ground_truth_file <- if (dataset_name %in% names(datasets)) {
    file.path(ground_truth_dir, datasets[[dataset_name]])
  } else {
    # Try to find a matching ground truth file
    gt_files <- list.files(ground_truth_dir, pattern = paste0(".*", dataset_name, ".*\\.tsv$"), full.names = TRUE)
    if (length(gt_files) == 0) {
      gt_files <- list.files(ground_truth_dir, pattern = "\\.tsv$", full.names = TRUE)
      if (length(gt_files) > 0) {
        cat("Warning: Using first available ground truth file:", basename(gt_files[1]), "\n")
        gt_files[1]
      } else {
        stop("No ground truth file found for dataset: ", dataset_name)
      }
    } else {
      gt_files[1]
    }
  }
  
  if (!file.exists(ground_truth_file)) {
    stop("Ground truth file not found: ", ground_truth_file)
  }
  
  cat("Using ground truth file:", ground_truth_file, "\n")
  
  # Helper functions
  subsample_data <- function(df, max_points = 10000) {
    int.length <- nrow(df)
    if (int.length <= max_points) return(df)
    
    int.subsamp_size <- max_points
    vec.subsamp_index <- sample(1:int.length, int.subsamp_size)
    vec.subsamp_index.padded <- c(seq(1, min(10000, int.length), 10), 
                                 vec.subsamp_index, 
                                 seq(max(1, int.length-10000), int.length, 10))
    vec.subsamp_index.sorted <- unique(sort(vec.subsamp_index.padded))
    return(df[vec.subsamp_index.sorted, ])
  }
  
  create_pr_gap_plot <- function(combined_pr, results, method_colors, output_path,
                                gap_region = c(0.3, 0.9), plot_width = 8, plot_height = 8) {
    
    unique_methods <- unique(results$method)
    ordered_method_curves <- unlist(lapply(unique_methods, function(m) {
      c(paste0(m, "_Normal"), paste0(m, "_Randomized"))
    }))
    
    legend_df <- combined_pr %>%
      distinct(method_curve, legend_label, curve_type) %>%
      filter(method_curve %in% ordered_method_curves)
    
    legend_df$col <- sapply(legend_df$method_curve, function(curve) {
      base_method <- gsub("_(Normal|Randomized)$", "", curve)
      return(method_colors[base_method])
    })
    
    legend_df$lty <- ifelse(legend_df$curve_type == "Randomized", 2, 1)
    
    base_methods <- unique(gsub("_(Normal|Randomized)", "", legend_df$method_curve))
    ordered_df <- do.call(rbind, lapply(base_methods, function(method) {
      normal <- legend_df[legend_df$method_curve == paste0(method, "_Normal"), ]
      random <- legend_df[legend_df$method_curve == paste0(method, "_Randomized"), ]
      rbind(normal, random)
    }))
    
    pdf(output_path, width = plot_height, height = plot_height)
    par(mfrow = c(1, 1), tck = -0.02, cex.axis = 0.8, mgp = c(2, 1, 0))
    
    first_label <- ordered_df$legend_label[1]
    first_color <- ordered_df$col[1]
    first_lty <- ifelse(grepl("Random", first_label), 2, 1)
    
    df_first <- subset(combined_pr, legend_label == first_label)
    df_first <- df_first[is.finite(df_first$precision) & df_first$precision <= 1, ]
    df_first_subsampled <- subsample_data(df_first, max_points = 10000)
    
    plot(x = df_first_subsampled$recall, 
             y = df_first_subsampled$precision,
             type = "l", 
             col = first_color, lwd = 2, lty = first_lty,
             ylim = c(0, 1), xlim = c(0, 1),
             xlab = "Recall", ylab = "Precision", 
             main = paste("PR Curve Benchmark -", dataset_name)
             )
    
    for (i in 2:nrow(ordered_df)) {
      df_i <- subset(combined_pr, legend_label == ordered_df$legend_label[i])
      df_i <- df_i[is.finite(df_i$precision) & df_i$precision <= 1, ]
      df_i_subsampled <- subsample_data(df_i, max_points = 10000)
      
      this_lty <- ifelse(grepl("Random", ordered_df$legend_label[i]), 2, 1)
      
          lines(x = df_i_subsampled$recall, y = df_i_subsampled$precision,
               type = "l", 
               col = ordered_df$col[i], lwd = 2, lty = this_lty
               )
    }
    
    legend("topright", legend = ordered_df$legend_label,
           col = ordered_df$col, lwd = 2,
           lty = ifelse(grepl("Random", ordered_df$legend_label), 2, 1), cex = 1)
    
    dev.off()
    cat("PR gap plot saved to:", output_path, "\n")
  }
  
  # Load and preprocess ground truth
  if (str_ends(ground_truth_file, ".csv")) {
    ground_truth <- read_csv(ground_truth_file, col_types = cols())
  } else {
    ground_truth <- read_tsv(ground_truth_file, col_types = cols())
  }
  
  ground_truth <- ground_truth %>% 
    select(Source, Target) %>%
    mutate(Source = str_to_upper(Source), Target = str_to_upper(Target)) %>%
    distinct()
  
  gt_pairs <- paste(ground_truth$Source, ground_truth$Target, sep = "_")
  gt_pairs1 <- unique(gt_pairs)
  tested_tfs <- unique(ground_truth$Source)
  tested_targets <- unique(ground_truth$Target)
  
  cat("Ground truth contains", length(gt_pairs1), "unique interactions\n")
  
  # Initialize results storage
  results <- data.frame(
    method = character(), AUROC = numeric(), AUPRC = numeric(),
    AUPRC_random = numeric(), n_total_edges = integer(),
    n_evaluable_edges = integer(), n_positives = integer(),
    n_negatives = integer(), stringsAsFactors = FALSE
  )
  
  roc_list <- list()
  pr_list <- list()
  random_pr_list <- list()
  
  # Process each method
  for (i in 1:nrow(method_info)) {
    set.seed(42 + i)
    
    method_name <- method_info$method[i]
    tf_col <- method_info$tf_col[i]
    target_col <- method_info$target_col[i]
    score_col <- method_info$score_col[i]
    
    cat("Processing:", method_name, "\n")
    
    # Determine file path
    if (method_name == new_method_name) {
      # Use the provided new GRN file
      file <- new_grn_file
    } else {
      # Look for existing method file
      if (is.null(input_dir)) {
        cat("Skipping existing method", method_name, "- no input directory provided\n")
        next
      }
      
      method_path <- file.path(input_dir, dataset_name, method_name)
      if (!dir.exists(method_path)) {
        cat("Method directory not found:", method_path, ". Skipping.\n")
        next
      }
      
      file_list <- list.files(method_path, full.names = TRUE)
      if (length(file_list) == 0) {
        cat("No files found for method:", method_name, ". Skipping.\n")
        next
      }
      file <- file_list[1]
    }
    
    # Load method data
    df <- if (str_ends(file, ".csv")) {
      read_csv(file, col_types = cols())
    } else {
      read_tsv(file, col_types = cols())
    }
    
    # Check if required columns exist
    if (!tf_col %in% colnames(df)) {
      cat("Warning: TF column '", tf_col, "' not found for", method_name, ". Skipping.\n")
      next
    }
    if (!target_col %in% colnames(df)) {
      cat("Warning: Target column '", target_col, "' not found for", method_name, ". Skipping.\n")
      next
    }
    if (is.na(score_col) || !score_col %in% colnames(df)) {
      cat("Warning: Score column '", score_col, "' not found for", method_name, ". Skipping.\n")
      next
    }
    
    # Preprocess inferred data
    inferred <- df %>%
      select(tf = all_of(tf_col), target = all_of(target_col), score = all_of(score_col)) %>%
      mutate(score = abs(score), tf = str_to_upper(tf), target = str_to_upper(target)) %>%
      group_by(tf, target) %>%
      summarise(score = max(score), .groups = "drop")
    
    n_total_edges <- nrow(inferred)
    
    if (use_filtered_approach) {
      inferred <- inferred %>% filter(tf %in% tested_tfs & target %in% tested_targets)
      cat("  Filtered to tested space:", nrow(inferred), "edges\n")
    }
    
    # Add labels
    inferred <- inferred %>%
      mutate(pair = paste(tf, target, sep = "_"), label = ifelse(pair %in% gt_pairs1, 1, 0))
    
    positives <- inferred %>% filter(label == 1)
    negatives <- inferred %>% filter(label == 0)
    
    if (nrow(positives) == 0 || nrow(negatives) == 0) {
      cat("  Insufficient positive or negative examples for", method_name, ". Skipping.\n")
      next
    }
    
    cat("  Positives:", nrow(positives), "Negatives:", nrow(negatives), "\n")
    
    # AUROC calculation
    if (nrow(negatives) < nrow(positives)) {
      pos_sampled <- positives %>% sample_n(nrow(negatives))
      neg_sampled <- negatives
    } else {
      pos_sampled <- positives
      neg_sampled <- negatives %>% sample_n(nrow(positives))
    }
    
    auroc_data <- bind_rows(pos_sampled, neg_sampled)
    roc_obj <- roc(response = auroc_data$label, predictor = auroc_data$score, 
                   levels = c(0, 1), direction = "<")
    auroc_value <- as.numeric(auc(roc_obj))
    
    roc_df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities, method = method_name)
    roc_df <- roc_df[order(roc_df$fpr), ]
    roc_list[[method_name]] <- roc_df
    
    # AUPRC calculation
    neg_sampled_10x <- negatives %>% sample_n(size = min(nrow(positives) * 10, nrow(negatives)))
    pr_obj <- pr.curve(scores.class0 = positives$score, scores.class1 = neg_sampled_10x$score, curve = TRUE)
    auprc_value <- pr_obj$auc.integral
    
    pr_df <- data.frame(recall = pr_obj$curve[, 1], precision = pr_obj$curve[, 2], method = method_name)
    pr_list[[method_name]] <- pr_df
    
    # Random AUPRC
    inferred_random <- inferred %>% mutate(score = runif(n()))
    positives_r <- inferred_random %>% filter(label == 1)
    negatives_r <- inferred_random %>% filter(label == 0)
    neg_sampled_10x_r <- negatives_r %>% sample_n(size = min(nrow(positives_r) * 10, nrow(negatives_r)))
    
    pr_obj_random <- pr.curve(scores.class0 = positives_r$score, scores.class1 = neg_sampled_10x_r$score, curve = TRUE)
    auprc_random <- pr_obj_random$auc.integral
    
    pr_df_random <- data.frame(recall = pr_obj_random$curve[, 1], precision = pr_obj_random$curve[, 2], method = method_name)
    random_pr_list[[method_name]] <- pr_df_random
    
    # Store results
    results <- rbind(results, data.frame(
      method = method_name, AUROC = auroc_value, AUPRC = auprc_value,
      AUPRC_random = auprc_random, n_total_edges = n_total_edges,
      n_evaluable_edges = nrow(inferred), n_positives = nrow(positives),
      n_negatives = nrow(negatives), stringsAsFactors = FALSE
    ))
    
    cat("  AUROC:", round(auroc_value, 4), "AUPRC:", round(auprc_value, 4), "\n")
  }
  
  if (nrow(results) == 0) {
    stop("No methods could be evaluated. Please check your inputs.")
  }
  
  # Ensure new method appears first in results
  new_method_result <- results[results$method == new_method_name, ]
  other_results <- results[results$method != new_method_name, ]
  results <- bind_rows(new_method_result, other_results)
  
  # Prepare plot data
  results <- results %>%
    mutate(method_auroc = paste0(method, " (", round(AUROC, 3), ")"),
           method_auprc = paste0(method, " (AUPRC=", round(AUPRC, 3), ")"))
  
  roc_combined <- bind_rows(roc_list, .id = "method") %>% left_join(results, by = "method")
  pr_combined <- bind_rows(pr_list, .id = "method") %>% left_join(results, by = "method")
  random_pr_combined <- bind_rows(random_pr_list, .id = "method") %>% left_join(results, by = "method")
  
  # Color scheme - highlight new method
  unique_methods <- unique(results$method)
  num_methods <- length(unique_methods)
  
  # Use red for new method, other colors for existing methods
  method_colors <- setNames(rep("gray60", num_methods), unique_methods)
  method_colors[new_method_name] <- "red"
  
  if (num_methods > 1) {
    other_methods <- unique_methods[unique_methods != new_method_name]
    palette_colors <- brewer.pal(min(max(length(other_methods), 3), 8), "Dark2")
    method_colors[other_methods] <- palette_colors[seq_len(length(other_methods))]
  }
  
  # Custom theme
  custom_theme <- theme_classic(base_size = 14) +
    theme(panel.grid = element_blank(), axis.line = element_line(linewidth = 0.5),
          plot.title = element_text(hjust = 0.5), legend.title = element_blank())
  
  # ROC Plot
  roc_combined$method <- factor(roc_combined$method, levels = unique_methods)
  roc_plot <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = method)) +
    geom_line(linewidth = ifelse(roc_combined$method == new_method_name, 1.2, 0.6)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
    labs(x = "False Positive Rate", y = "True Positive Rate", 
         title = paste("ROC Curve Benchmark -", dataset_name), color = "Method (AUROC)") +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_color_manual(values = method_colors, labels = results$method_auroc) +
    custom_theme
  
  # PR Plot Data Preparation
  random_pr_combined <- random_pr_combined %>%
    mutate(curve_type = "Randomized", legend_label = paste0(method, " Random AUPRC=", round(AUPRC_random, 3)))
  
  pr_combined <- pr_combined %>%
    mutate(curve_type = "Normal", legend_label = paste0(method, " AUPRC=", round(AUPRC, 3)))
  
  combined_pr <- bind_rows(random_pr_combined, pr_combined) %>%
    mutate(method_curve = paste0(method, "_", curve_type))
  
  # Save plots
  ggsave(file.path(output_dir, "ROC_benchmark.png"), plot = roc_plot, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "ROC_benchmark.pdf"), plot = roc_plot, width = 8, height = 6)
  
  # Create PR gap plot
  create_pr_gap_plot(combined_pr = combined_pr, results = results, method_colors = method_colors,
                    output_path = file.path(output_dir, "PR_benchmark.pdf"))
  
  # Save results CSV
  write_csv(results, file.path(output_dir, "benchmark_results.csv"))
  
  # Print summary
  cat("\n=== BENCHMARK SUMMARY ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("New method:", new_method_name, "\n")
  cat("Methods compared:", nrow(results), "\n")
  cat("\nPerformance Summary:\n")
  print(results %>% 
        select(method, AUROC, AUPRC) %>% 
        mutate(AUROC = round(AUROC, 4), AUPRC = round(AUPRC, 4)) %>%
        arrange(desc(AUPRC)))
  
  cat("\nResults saved to:", output_dir, "\n")
  
  return(list(
    results = results, 
    roc_data = roc_combined, 
    pr_data = combined_pr, 
    plots = list(roc = roc_plot),
    new_method_performance = results[results$method == new_method_name, ]
  ))
}
