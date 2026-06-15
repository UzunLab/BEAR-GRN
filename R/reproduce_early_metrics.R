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
                                    output_dir        = "filtered_grn_metrics",
                                    ground_truth_dir  = "GROUND.TRUTHS",
                                    max_edges         = 10000,
                                    method_colors     = NULL) {

  cat("=== Compute Filtered GRN Metrics ===\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Max edges per method:", max_edges, "\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }

  # ── Method-specific column mappings ────────────────────────────────────────
  column_mappings <- list(
    "CellOracle" = c("source", "target"),
    "FigR"       = c("Motif",  "DORC"),
    "GRaNIE"     = c("TF.name","gene.name"),
    "LINGER"     = c("Source", "Target"),
    "Pando"      = c("tf",     "target"),
    "Pando_xgb"  = c("tf",     "target"),
    "SCENIC+"    = c("Source", "Target"),
    "Scenic+"    = c("Source", "Target"),   # alternative spelling
    "TRIPOD"     = c("TF",     "gene"),
    "DIRECTNET"  = c("TF",     "Gene")
  )

  score_mappings <- list(
    "CellOracle" = "coef_abs",
    "FigR"       = "Score",
    "GRaNIE"     = "TF_gene.r",
    "LINGER"     = "Score",
    "Pando"      = "estimate",
    "Pando_xgb"  = "corr",
    "SCENIC+"    = "Score",
    "Scenic+"    = "Score",
    "TRIPOD"     = "abs_coef",
    "DIRECTNET"  = "Score"
  )

  # ── Default colour palette ──────────────────────────────────────────────────
  if (is.null(method_colors)) {
    method_colors <- c(
      "CellOracle" = "#1B9E77",
      "FigR"       = "#66A61E",
      "GRaNIE"     = "#A6761D",
      "LINGER"     = "#E7298A",
      "Pando"      = "#7570B3",
      "Pando_xgb"  = "#666666",
      "SCENIC+"    = "#D95F02",
      "Scenic+"    = "#D95F02",
      "TRIPOD"     = "#E6AB02",
      "DIRECTNET"  = "#377EB8"
    )
  }

  # ── Dataset → ground-truth file mapping ────────────────────────────────────
  datasets_mapping <- list(
    "K562"           = "filtered_RN117_K562.tsv",
    "Macrophage_S1"  = "filtered_RN204_Buffer1.tsv",
    "Macrophage_S2"  = "filtered_RN204_Buffer2.tsv",
    "mESC_E7.5_rep1" = "filtered_RN111_E7.5_rep1.tsv",
    "mESC_E7.5_rep2" = "filtered_RN111_E7.5_rep2.tsv",
    "mESC_E8.5_rep1" = "filtered_RN111_E8.5_rep1.tsv",
    "mESC_E8.5_rep2" = "filtered_RN111_E8.5_rep2.tsv",
    "E7.5_rep1"      = "filtered_RN111_E7.5_rep1.tsv",
    "E7.5_rep2"      = "filtered_RN111_E7.5_rep2.tsv",
    "E8.5_rep1"      = "filtered_RN111_E8.5_rep1.tsv",
    "E8.5_rep2"      = "filtered_RN111_E8.5_rep2.tsv",
    "iPS"            = "filtered_RN000_iPS.tsv",
    "Naive_mESC"     = "filtered_RN111_Naive_mESC.tsv"
  )

  # ── Auto-detect datasets ────────────────────────────────────────────────────
  dataset_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  dataset_dirs <- dataset_dirs[dataset_dirs != ""]

  if (length(dataset_dirs) == 0) {
    stop("No dataset directories found in input_dir: ", input_dir)
  }

  cat("Found datasets:", paste(dataset_dirs, collapse = ", "), "\n")

  # ── Storage ─────────────────────────────────────────────────────────────────
  plot_data       <- tibble(Metric = character(), Score = numeric(),
                            Dataset = character(), Method = character())
  all_summaries   <- list()
  failed_datasets <- c()

  # ── Main analysis loop ──────────────────────────────────────────────────────
  for (dataset in dataset_dirs) {
    cat("\nProcessing dataset:", dataset, "\n")

    # Determine ground truth file
    ground_truth_file <- if (dataset %in% names(datasets_mapping)) {
      file.path(ground_truth_dir, datasets_mapping[[dataset]])
    } else {
      gt_files <- list.files(ground_truth_dir,
                             pattern    = paste0(".*", dataset, ".*\\.tsv$"),
                             full.names = TRUE)
      if (length(gt_files) == 0) {
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

    if (is.null(gt)) { failed_datasets <- c(failed_datasets, dataset); next }

    gt_edges   <- unique(tibble(source = toupper(gt$Source), target = toupper(gt$Target)))
    gt_pairs   <- with(gt_edges, paste(source, target, sep = "_"))
    gt_tfs     <- unique(gt_edges$source)
    gt_targets <- unique(gt_edges$target)

    cat("  Ground truth:", length(gt_pairs), "edges,",
        length(gt_tfs), "TFs,", length(gt_targets), "targets\n")

    # Get available methods
    dataset_input_dir <- file.path(input_dir, dataset)
    method_dirs <- list.dirs(dataset_input_dir, full.names = FALSE, recursive = FALSE)
    methods     <- method_dirs[method_dirs != ""]

    if (length(methods) == 0) {
      cat("  No method directories found in", dataset_input_dir, "\n")
      failed_datasets <- c(failed_datasets, dataset)
      next
    }

    summary_list <- list()

    # ── Per-method processing ────────────────────────────────────────────────
    for (method in methods) {
      if (!(method %in% names(column_mappings))) {
        cat("    ", method, ": Not in column mappings, skipping\n")
        next
      }

      method_dir <- file.path(dataset_input_dir, method)
      files      <- list.files(method_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)

      if (length(files) == 0) {
        cat("    ", method, ": No data files found, skipping\n")
        next
      }

      file      <- files[1]
      cols      <- column_mappings[[method]]
      score_col <- score_mappings[[method]]

      df <- tryCatch({
        sep     <- if (grepl("\\.csv$", file)) "," else "\t"
        temp_df <- read_delim(file, delim = sep, col_types = cols(), show_col_types = FALSE)

        if (length(cols) >= 2 && all(cols %in% names(temp_df))) {
          temp_df <- temp_df %>% rename(source = all_of(cols[1]), target = all_of(cols[2]))
        } else {
          stop("Required columns not found in file")
        }

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

      # Filter by GT TFs and targets
      df_filtered <- df %>%
        mutate(source_upper = toupper(source), target_upper = toupper(target)) %>%
        filter(source_upper %in% gt_tfs, target_upper %in% gt_targets)

      cat("    ", method, ": Filtered network size:", nrow(df_filtered), "edges\n")

      n_edges <- min(max_edges, nrow(df_filtered))
      if (n_edges == 0) {
        cat("    ", method, ": No edges after filtering, skipping\n")
        next
      }

      df_final <- df_filtered %>% head(n_edges)
      cat("    ", method, ": Final network size:", nrow(df_final), "edges\n")

      df_pairs <- df_final %>%
        transmute(pair = paste(source_upper, target_upper, sep = "_"))

      # Compute metrics
      tp <- intersect(df_pairs$pair, gt_pairs)
      fp <- setdiff(df_pairs$pair,  gt_pairs)
      fn <- setdiff(gt_pairs,        df_pairs$pair)

      precision <- ifelse((length(tp) + length(fp)) > 0,
                          length(tp) / (length(tp) + length(fp)), 0)
      recall    <- ifelse((length(tp) + length(fn)) > 0,
                          length(tp) / (length(tp) + length(fn)), 0)
      f1        <- ifelse((precision + recall) > 0,
                          2 * precision * recall / (precision + recall), 0)

      cat("    ", method,
          ": Precision =", round(precision, 4),
          ", Recall =",    round(recall,    4),
          ", F1 =",        round(f1,        4), "\n")

      plot_data <- bind_rows(plot_data, tibble(
        Metric  = c("Precision", "Recall", "F1-score"),
        Score   = c(precision, recall, f1),
        Dataset = dataset,
        Method  = method
      ))

      inferred_tfs     <- unique(df_final$source_upper)
      inferred_targets <- unique(df_final$target_upper)

      summary_list[[method]] <- tibble(
        Method                = method,
        GT_total_edges        = length(gt_pairs),
        GT_TFs                = length(gt_tfs),
        GT_targets            = length(gt_targets),
        Original_network_size = nrow(df),
        Filtered_network_size = nrow(df_filtered),
        Final_network_size    = nrow(df_final),
        Inferred_TFs          = length(inferred_tfs),
        Inferred_targets      = length(inferred_targets),
        TP                    = length(tp),
        FP                    = length(fp),
        FN                    = length(fn),
        Precision             = round(precision, 4),
        Recall                = round(recall,    4),
        F1_score              = round(f1,        4),
        Common_TFs            = length(intersect(gt_tfs,     inferred_tfs)),
        Common_targets        = length(intersect(gt_targets, inferred_targets))
      )
    }

    # Save per-dataset summary
    if (length(summary_list) > 0) {
      df_summary <- bind_rows(summary_list) %>% mutate(Dataset = dataset)
      write_tsv(df_summary,
                file.path(output_dir, paste0(dataset, "_filtered_metrics.tsv")))
      all_summaries[[dataset]] <- df_summary
      cat("  Saved summary to:",
          file.path(output_dir, paste0(dataset, "_filtered_metrics.tsv")), "\n")
    }
  }

  # ── Early exit if nothing to plot ───────────────────────────────────────────
  if (nrow(plot_data) == 0) {
    cat("No data available for plotting.\n")
    return(list(
      plot_data       = plot_data,
      summaries       = all_summaries,
      failed_datasets = failed_datasets,
      plots           = list()
    ))
  }

  # ── Prepare plot data ───────────────────────────────────────────────────────
  plot_data$Dataset <- str_trim(plot_data$Dataset)

  # Custom dataset order — reversed so first item appears at TOP of Y-axis
  desired_order      <- c("Macrophage_S1", "Macrophage_S2", "K562", "iPS",
                          "E7.5_rep1", "E7.5_rep2", "E8.5_rep1", "E8.5_rep2", "Naive_mESC")
  available_datasets <- unique(plot_data$Dataset)
  dataset_order      <- intersect(desired_order, available_datasets)
  plot_data$Dataset  <- factor(plot_data$Dataset, levels = rev(dataset_order))

  # Method order
  all_method_order  <- c("CellOracle", "FigR", "GRaNIE", "LINGER",
                          "Pando", "Pando_xgb", "SCENIC+", "Scenic+", "TRIPOD", "DIRECTNET")
  available_methods <- intersect(all_method_order, unique(plot_data$Method))
  plot_data$Method  <- factor(plot_data$Method, levels = available_methods)

  # Metric order for facets
  plot_data$Metric  <- factor(plot_data$Metric,
                              levels = c("Precision", "Recall", "F1-score"))

  plots <- list()

  # ── Smart x-axis label formatter ───────────────────────────────────────────
  # Called once per facet panel by ggplot; picks decimal places based on the
  # magnitude of the largest value visible in that panel.
  #
  #   max >= 0.1   →  3 decimal places  e.g. 0.123   (Precision-like range)
  #   max >= 0.01  →  4 decimal places  e.g. 0.0123  (mid-range Recall / F1)
  #   max <  0.01  →  scientific        e.g. 1.23e-03 (very small values)
  #
  smart_label <- function(x) {
    max_val <- max(x, na.rm = TRUE)
    if (is.na(max_val) || max_val == 0) return(as.character(x))

    if (max_val >= 0.1) {
      scales::number_format(accuracy = 0.001)(x)
    } else if (max_val >= 0.01) {
      scales::number_format(accuracy = 0.0001)(x)
    } else {
      scales::scientific_format(digits = 3)(x)
    }
  }

  # ── Combined lollipop plot ──────────────────────────────────────────────────
  cat("\nCreating combined lollipop plot (Datasets on Y, free X scale per metric)\n")

  df_plot_all <- plot_data %>% filter(Method %in% available_methods)

  if (nrow(df_plot_all) > 0) {

    p_combined <- ggplot(df_plot_all, aes(y = Dataset, x = Score)) +
      # Lollipop stick — horizontal, from x = 0 to x = Score
      #geom_segment(
      #  aes(y = Dataset, yend = Dataset, x = 0, xend = Score),
      #  color     = "grey70",
      #  linewidth = 0.6,
      #  alpha     = 0.7,
      #  position  = position_dodge(width = 0.8)
      #) +
      # Dots
      geom_point(
        aes(color = Method),
        shape    = 16,
        size     = 5,           # kept at 6 as requested
        alpha    = 0.9,
        position = position_dodge(width = 0.8)
      ) +
      # ── scales = "free_x": each metric panel gets its own x range ──────────
      facet_grid(. ~ Metric, scales = "free_x") +
      scale_color_manual(values = method_colors) +
      # Lower bound fixed at 0; upper bound auto per panel via free_x.
      # smart_label picks the right number of decimal places per panel range.
      scale_x_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.08)),
        labels = smart_label
      ) +
      labs(
        title    = paste0("GRN Metrics Comparison (Top ",
                          format(max_edges, scientific = FALSE), " Edges)"),
        subtitle = "Networks filtered by ground truth TFs and targets",
        y        = "Dataset",
        x        = "Score",
        color    = "Method"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        strip.text         = element_text(size = 12, face = "bold"),
        strip.background   = element_rect(fill = "grey95", color = NA),
        legend.position    = "right",
        legend.title       = element_text(size = 11, face = "bold"),
        legend.text        = element_text(size = 10),
        axis.line.x        = element_line(color = "black", linewidth = 0.6),
        axis.line.y        = element_line(color = "black", linewidth = 0.6),
        axis.title         = element_text(size = 12, face = "bold"),
        axis.text          = element_text(size = 10),
        axis.text.x        = element_text(angle = 45, hjust = 1),
        plot.title         = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(size = 11, hjust = 0.5, color = "gray50"),
        panel.background   = element_rect(fill = "white", color = NA),
        plot.background    = element_rect(fill = "white", color = NA)
      )

    plots[["combined"]] <- p_combined

    # Height scales with number of datasets
    plot_height <- max(6, length(levels(plot_data$Dataset)) * 0.9)

    pdf_path <- file.path(output_dir, "combined_metrics_lollipop.pdf")
    ggsave(filename = pdf_path, plot = p_combined,
           width = 16, height = plot_height, device = "pdf", dpi = 300)

    png_path <- file.path(output_dir, "combined_metrics_lollipop.png")
    ggsave(filename = png_path, plot = p_combined,
           width = 16, height = plot_height, device = "png", dpi = 300)

    cat("Saved combined plot (PDF):", pdf_path, "\n")
    cat("Saved combined plot (PNG):", png_path, "\n")
  }

  # ── Combined summary across all datasets ────────────────────────────────────
  combined_summary <- NULL
  if (length(all_summaries) > 0) {
    combined_summary <- bind_rows(all_summaries)
    write_tsv(combined_summary,
              file.path(output_dir, "all_datasets_filtered_metrics.tsv"))
    cat("Combined summary saved to:",
        file.path(output_dir, "all_datasets_filtered_metrics.tsv"), "\n")
  }

  # ── Final summary ────────────────────────────────────────────────────────────
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("Total datasets found:",   length(dataset_dirs),   "\n")
  cat("Successfully processed:", length(all_summaries),  "\n")
  cat("Failed:",                 length(failed_datasets), "\n")

  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }

  cat("Analysis completed! Check the '", output_dir, "' folder for results.\n")

  return(list(
    plot_data        = plot_data,
    summaries        = all_summaries,
    combined_summary = combined_summary,
    failed_datasets  = failed_datasets,
    plots            = plots
  ))
}
