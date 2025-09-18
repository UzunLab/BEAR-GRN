# Global variable declarations to satisfy R CMD check
utils::globalVariables(c(
  # Column names used in dplyr operations
  "method_curve", "legend_label", "curve_type", "Source", "Target", 
  "score", "tf", "target", "pair", "label", "method", "AUROC", "AUPRC",
  "fpr", "tpr", "AUPRC_random", "source_upper", "target_upper",
  "Metric", "Method", "Dataset", "Score"
))
#' Reproduce ROC and PR Plots for Multiple GRN Methods and Datasets
#'
#' Automatically processes all datasets and methods in the input directory,
#' generates ROC and PR plots, and saves results to CSV files.
#'
#' @param input_dir Path to input directory containing dataset subdirectories
#' @param output_dir Path to output directory for results and plots
#' @param ground_truth_dir Path to directory containing ground truth files (default: "GROUND.TRUTHS")
#' @param use_filtered_approach Logical, whether to filter to tested TF-target space (default: TRUE)
#' @param method_info Optional data frame with method information. If NULL, uses default mapping.
#' @return List containing all results and plots
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
#' # Basic usage - auto-detect everything
#' results <- reproduce_ROC_PR_plots(
#'   input_dir = "INFERRED.GRNS",
#'   output_dir = "Results"
#' )
#' 
#' # With custom ground truth directory
#' results <- reproduce_ROC_PR_plots(
#'   input_dir = "my_data",
#'   output_dir = "my_results",
#'   ground_truth_dir = "my_ground_truths"
#' )
#' }

reproduce_ROC_PR_plots <- function(input_dir, 
                                  output_dir, 
                                  ground_truth_dir = "GROUND.TRUTHS",
                                  use_filtered_approach = TRUE,
                                  method_info = NULL) {
  

  
  cat("=== Reproduce ROC and PR Plots ===\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Create main output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Default method information
  if (is.null(method_info)) {
    method_info <- tribble(
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
  
  # Auto-detect datasets from input directory
  dataset_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  dataset_dirs <- dataset_dirs[dataset_dirs != ""]
  
  if (length(dataset_dirs) == 0) {
    stop("No dataset directories found in input_dir: ", input_dir)
  }
  
  cat("Found datasets:", paste(dataset_dirs, collapse = ", "), "\n")
  
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
             xlab = "Recall", ylab = "Precision", main = "PR Curve"
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
  
  # Main evaluation function for single dataset
  evaluate_dataset <- function(dataset_name, dataset_input_dir, ground_truth_file, dataset_output_dir) {
    
    cat("\n=== Processing Dataset:", dataset_name, "===\n")
    
    # Load ground truth
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
      
      method_path <- file.path(dataset_input_dir, method_name)
      if (!dir.exists(method_path)) {
        cat("Method directory not found:", method_path, "\n")
        next
      }
      
      file_list <- list.files(method_path, full.names = TRUE)
      if (length(file_list) == 0) {
        cat("No files found for method:", method_name, "\n")
        next
      }
      
      file <- file_list[1]
      cat("Processing:", method_name, "\n")
      
      # Load method data
      df <- if (str_ends(file, ".csv")) {
        read_csv(file, col_types = cols())
      } else {
        read_tsv(file, col_types = cols())
      }
      
      if (is.na(score_col) || !score_col %in% colnames(df)) {
        cat("Warning: Score column not found for", method_name, ". Skipping.\n")
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
      cat("No methods could be evaluated for dataset:", dataset_name, "\n")
      return(NULL)
    }
    
    # Prepare plot data
    results <- results %>%
      mutate(method_auroc = paste0(method, " (", round(AUROC, 3), ")"),
             method_auprc = paste0(method, " (AUPRC=", round(AUPRC, 3), ")"))
    
    roc_combined <- bind_rows(roc_list, .id = "method") %>% left_join(results, by = "method")
    pr_combined <- bind_rows(pr_list, .id = "method") %>% left_join(results, by = "method")
    random_pr_combined <- bind_rows(random_pr_list, .id = "method") %>% left_join(results, by = "method")
    
    # Color scheme
    unique_methods <- unique(results$method)
    num_methods <- length(unique_methods)
    palette_colors <- brewer.pal(min(max(num_methods, 3), 8), "Dark2")
    method_colors <- setNames(palette_colors[seq_len(num_methods)], unique_methods)
    
    # Custom theme
    custom_theme <- theme_classic(base_size = 14) +
      theme(panel.grid = element_blank(), axis.line = element_line(linewidth = 0.5),
            plot.title = element_text(hjust = 0.5), legend.title = element_blank())
    
    # ROC Plot
    roc_combined$method <- factor(roc_combined$method, levels = unique_methods)
    roc_plot <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = method)) +
      geom_line(linewidth = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
      labs(x = "False Positive Rate", y = "True Positive Rate", 
           title = paste("ROC Curve -", dataset_name), color = "Method (AUROC)") +
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
    ggsave(file.path(dataset_output_dir, "ROC_plot.png"), plot = roc_plot, width = 8, height = 6, dpi = 300)
    ggsave(file.path(dataset_output_dir, "ROC_plot.pdf"), plot = roc_plot, width = 8, height = 6)
    
    # Create PR gap plot
    create_pr_gap_plot(combined_pr = combined_pr, results = results, method_colors = method_colors,
                      output_path = file.path(dataset_output_dir, "PR_gap_plot.pdf"))
    
    # Save results CSV
    write_csv(results, file.path(dataset_output_dir, "evaluation_results.csv"))
    
    cat("Results saved to:", dataset_output_dir, "\n")
    
    return(list(results = results, roc_data = roc_combined, pr_data = combined_pr, plots = list(roc = roc_plot)))
  }
  
  # Initialize storage for all results
  all_results <- list()
  failed_datasets <- c()
  
  # Process each dataset
  for (dataset_name in dataset_dirs) {
    #cat("\n" , "="*50, "\n")
    cat("\n", strrep("=", 50), "\n")
    cat("Processing dataset:", dataset_name, "\n")
    
    # Determine ground truth file
    ground_truth_file <- if (dataset_name %in% names(datasets)) {
      file.path(ground_truth_dir, datasets[[dataset_name]])
    } else {
      # Try to find a matching ground truth file
      gt_files <- list.files(ground_truth_dir, pattern = paste0(".*", dataset_name, ".*\\.tsv$"), full.names = TRUE)
      if (length(gt_files) == 0) {
        gt_files <- list.files(ground_truth_dir, pattern = "\\.tsv$", full.names = TRUE)
        if (length(gt_files) > 0) gt_files[1] else NULL
      } else {
        gt_files[1]
      }
    }
    
    if (is.null(ground_truth_file) || !file.exists(ground_truth_file)) {
      cat("Warning: Ground truth file not found for", dataset_name, "\n")
      failed_datasets <- c(failed_datasets, dataset_name)
      next
    }
    
    # Set paths
    dataset_input_dir <- file.path(input_dir, dataset_name)
    dataset_output_dir <- file.path(output_dir, paste0(dataset_name, "_Results"))
    
    # Create output directory
    dir.create(dataset_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Check if input directory exists
    if (!dir.exists(dataset_input_dir)) {
      cat("Warning: Input directory does not exist:", dataset_input_dir, "\n")
      failed_datasets <- c(failed_datasets, dataset_name)
      next
    }
    
    # Run evaluation
    tryCatch({
      result <- evaluate_dataset(dataset_name, dataset_input_dir, ground_truth_file, dataset_output_dir)
      if (!is.null(result)) {
        all_results[[dataset_name]] <- result
        cat("Successfully processed:", dataset_name, "\n")
      } else {
        failed_datasets <- c(failed_datasets, dataset_name)
      }
    }, error = function(e) {
      cat("Error processing", dataset_name, ":", e$message, "\n")
      failed_datasets <- c(failed_datasets, dataset_name)
    })
  }
  
  # Create summary results
  if (length(all_results) > 0) {
    summary_results <- bind_rows(lapply(names(all_results), function(dataset) {
      all_results[[dataset]]$results %>% mutate(dataset = dataset)
    }))
    
    write_csv(summary_results, file.path(output_dir, "all_datasets_summary.csv"))
    cat("\nSummary results saved to:", file.path(output_dir, "all_datasets_summary.csv"), "\n")
  }
  
  # Final summary
  cat("\n=== EVALUATION SUMMARY ===\n")
  cat("Total datasets found:", length(dataset_dirs), "\n")
  cat("Successfully processed:", length(all_results), "\n")
  cat("Failed:", length(failed_datasets), "\n")
  
  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  return(list(
    results = all_results,
    summary = if (exists("summary_results")) summary_results else NULL,
    failed = failed_datasets
  ))
}
