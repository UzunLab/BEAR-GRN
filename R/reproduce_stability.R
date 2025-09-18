#' Analyze Network Stability Using Jaccard Index
#'
#' This function analyzes the stability of gene regulatory networks across different
#' methods and samples by calculating Jaccard Index between top and random edge sets.
#' Optionally creates a heatmap visualization and saves results.
#'
#' @param base_dir Character. Base directory containing sample directories.
#' @param sample_dirs Character vector. Names of sample directories to analyze.
#' @param methods Character vector. Names of methods to analyze.
#' @param custom_sort_columns Named list. Custom sorting columns for methods (optional).
#' @param custom_source_target_columns Named list. Custom source/target column mappings (optional).
#' @param custom_delim Named list. Custom delimiter mappings for file reading (optional).
#' @param top_percent Numeric. Percentage of top edges to extract (default: 0.1 for 10%).
#' @param random_percent Numeric. Percentage of random edges to extract (default: 0.1 for 10%).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param file_pattern Character. Pattern to match network files (default: "\\.(csv|tsv|txt)$").
#' @param create_plot Logical. Whether to create a heatmap visualization (default: TRUE).
#' @param save_plot Logical. Whether to save the plot (default: FALSE).
#' @param save_results Logical. Whether to save results as CSV (default: FALSE).
#' @param plot_filename Character. Filename for saved plot (default: "median_ji_heatmap.png").
#' @param results_filename Character. Filename for saved results (default: "median_ji_results.csv").
#' @param plot_width Numeric. Plot width in inches (default: 12).
#' @param plot_height Numeric. Plot height in inches (default: 7).
#' @param plot_dpi Numeric. Plot resolution in DPI (default: 300).
#' @param verbose Logical. Whether to print detailed messages (default: FALSE).
#'
#' @return A list containing:
#'   \item{results}{Data.frame with columns: Sample, Method, Category, Median_JI}
#'   \item{plot}{ggplot2 object (if create_plot = TRUE)}
#'
#' @import dplyr
#' @import readr
#' @import purrr
#' @import ggplot2
#' @importFrom stats median
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off png
#' @importFrom viridis viridis
#'
#' @examples
#' \dontrun{
#' # Basic usage with visualization
#' results <- analyze_network_stability(
#'   base_dir = "STABILITY_GRNS",
#'   sample_dirs = c("K562", "Macrophage_S1"),
#'   methods = c("CellOracle", "FigR"),
#'   create_plot = TRUE,
#'   save_plot = TRUE,
#'   save_results = TRUE
#' )
#'
#' # Debug mode
#' results <- analyze_network_stability(
#'   base_dir = "STABILITY_GRNS",
#'   sample_dirs = c("K562"),
#'   methods = c("DIRECTNET"),
#'   verbose = TRUE
#' )
#' }
#'
#' @export
analyze_network_stability <- function(
    base_dir = ".",
    sample_dirs = c("K562", "Macrophage_S1", "Macrophage_S2", "E7.5_rep1", 
                    "E7.5_rep2", "E8.5_rep1", "E8.5_rep2"),
    methods = c("CellOracle", "FigR", "TRIPOD", "Pando", "LINGER", 
                "Scenic", "GRaNIE", "DIRECTNET"),
    custom_sort_columns = NULL,
    custom_source_target_columns = NULL,
    custom_delim = NULL,
    top_percent = 0.1,
    random_percent = 0.1,
    seed = 123,
    file_pattern = "\\.(csv|tsv|txt)$",
    create_plot = TRUE,
    save_plot = FALSE,
    save_results = FALSE,
    plot_filename = "median_ji_heatmap.png",
    results_filename = "median_ji_results.csv",
    plot_width = 12,
    plot_height = 14,
    plot_dpi = 300,
    verbose = FALSE) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Default configurations
  default_sort_columns <- list(
    CellOracle = "coef_abs",
    FigR = "Score",
    TRIPOD = "abs_coef",
    Pando = "estimate",
    LINGER = "Score",
    Scenic = "Score",
    GRaNIE = "TF_gene.r"
    # DIRECTNET has no score column by default
  )
  
  default_source_target_columns <- list(
    CellOracle = c("source", "target"),
    FigR = c("Motif", "DORC"),
    TRIPOD = c("TF", "gene"),
    Pando = c("tf", "target"),
    LINGER = c("Source", "Target"),
    Scenic = c("Source", "Target"),
    GRaNIE = c("TF.name", "gene.name"),
    DIRECTNET = c("TF", "Gene")
  )
  
  default_delim <- list(
    CellOracle = ",",
    FigR = ",",
    TRIPOD = ",",
    Pando = ",",
    LINGER = ",",
    Scenic = "\t",
    GRaNIE = ",",
    DIRECTNET = ","
  )
  
  # Use custom configurations if provided, otherwise use defaults
  sort_columns <- if (!is.null(custom_sort_columns)) custom_sort_columns else default_sort_columns
  source_target_columns <- if (!is.null(custom_source_target_columns)) custom_source_target_columns else default_source_target_columns
  delim <- if (!is.null(custom_delim)) custom_delim else default_delim
  
  # Validate inputs
  if (!is.numeric(top_percent) || top_percent <= 0 || top_percent > 1) {
    stop("top_percent must be a numeric value between 0 and 1")
  }
  if (!is.numeric(random_percent) || random_percent <= 0 || random_percent > 1) {
    stop("random_percent must be a numeric value between 0 and 1")
  }
  
  # Initialize results dataframe
  median_ji_results <- data.frame(
    Sample = character(),
    Method = character(),
    Category = character(),
    Median_JI = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each sample
  for (sample in sample_dirs) {
    sample_path <- file.path(base_dir, sample)
    
    if (!dir.exists(sample_path)) {
      warning(paste("Sample directory does not exist:", sample_path))
      next
    }
    
    if (verbose) cat("Processing sample:", sample, "\n")
    
    for (method in methods) {
      method_path <- file.path(sample_path, method)
      
      if (!dir.exists(method_path)) {
        warning(paste("Method directory does not exist:", method_path))
        next
      }
      
      if (verbose) cat("  Processing method:", method, "\n")
      
      network_files <- list.files(method_path, pattern = file_pattern, full.names = TRUE)
      
      if (length(network_files) < 2) {
        warning(paste("Not enough network files found for", sample, method, "- need at least 2, found", length(network_files)))
        next
      }
      
      if (verbose) cat("    Found", length(network_files), "files\n")
      
      tryCatch({
        # Process each file
        network_lists <- list()
        for (i in seq_along(network_files)) {
          if (verbose) cat("      Processing file", i, ":", basename(network_files[i]), "\n")
          network_lists[[i]] <- read_and_sort_network_internal(
            network_files[i], method, sort_columns, source_target_columns, 
            delim, top_percent, random_percent, verbose
          )
        }
        
        # Separate top and random networks
        top_networks <- purrr::map(network_lists, "top")
        random_networks <- purrr::map(network_lists, "random")
        
        # Compute Jaccard Index for top and random networks
        ji_top <- compute_pairwise_jaccard_internal(top_networks, verbose)
        ji_random <- compute_pairwise_jaccard_internal(random_networks, verbose)
        
        if (verbose) {
          cat("      Jaccard indices - Top:", length(ji_top), "Random:", length(ji_random), "\n")
          cat("      Median JI - Top:", median(ji_top, na.rm = TRUE), "Random:", median(ji_random, na.rm = TRUE), "\n")
        }
        
        # Append results
        new_results <- data.frame(
          Sample = c(
            paste(sample, paste0("(Top ", top_percent * 100, "%)")),
            paste(sample, paste0("(Random ", random_percent * 100, "%)"))
          ),
          Method = c(method, method),
          Category = c(
            paste0("Top ", top_percent * 100, "%"),
            paste0("Random ", random_percent * 100, "%")
          ),
          Median_JI = c(
            median(ji_top, na.rm = TRUE),
            median(ji_random, na.rm = TRUE)
          ),
          stringsAsFactors = FALSE
        )
        
        median_ji_results <- rbind(median_ji_results, new_results)
        
      }, error = function(e) {
        warning(paste("Error processing", sample, method, ":", e$message))
        if (verbose) cat("      ERROR:", e$message, "\n")
      })
    }
  }
  
  # Convert to factors for proper ordering
  if (nrow(median_ji_results) > 0) {
    median_ji_results$Sample <- factor(median_ji_results$Sample, 
                                       levels = rev(unique(median_ji_results$Sample)))
    median_ji_results$Method <- factor(median_ji_results$Method, 
                                       levels = unique(median_ji_results$Method))
  }
  
  # Create plot if requested
  plot_obj <- NULL
  if (create_plot && nrow(median_ji_results) > 0) {
    plot_obj <- create_jaccard_heatmap_internal(median_ji_results)
    
    # Save plot if requested
    if (save_plot) {
      ggplot2::ggsave(plot_filename, plot = plot_obj, width = plot_width, 
                      height = plot_height, dpi = plot_dpi)
      message(paste("Plot saved as:", plot_filename))
    }
  } else if (create_plot && nrow(median_ji_results) == 0) {
    warning("No results to plot")
  }
  
  # Save results if requested
  if (save_results && nrow(median_ji_results) > 0) {
    utils::write.csv(median_ji_results, results_filename, row.names = FALSE)
    message(paste("Results saved as:", results_filename))
  }
  
  # Return results and plot
  result_list <- list(results = median_ji_results)
  if (create_plot && !is.null(plot_obj)) {
    result_list$plot <- plot_obj
  }
  
  return(result_list)
}


# INTERNAL HELPER FUNCTIONS

#' Read and Sort Network Data - Internal
#' @keywords internal
#' @importFrom readr read_tsv
#' @importFrom dplyr filter select mutate arrange group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal
#' @importFrom stats median
read_and_sort_network_internal <- function(file, method, sort_columns, source_target_columns, 
                                          delim, top_percent, random_percent, verbose = FALSE) {
  
  # Check if file exists
  if (!file.exists(file)) {
    stop(paste("File does not exist:", file))
  }
  
  # Try to read the file
  tryCatch({
    network <- readr::read_delim(file, delim = delim[[method]], show_col_types = FALSE)
  }, error = function(e) {
    stop(paste("Failed to read file", basename(file), "with delimiter '", delim[[method]], "' - Error:", e$message))
  })
  
  # Check if network is empty
  if (nrow(network) == 0) {
    stop(paste("Empty network file:", basename(file)))
  }
  
  # Get column names
  source_col <- source_target_columns[[method]][1]
  target_col <- source_target_columns[[method]][2]
  score_col <- sort_columns[[method]]  # may be NULL for DIRECTNET
  
  # Print column info for debugging
  if (verbose) {
    cat("        File:", basename(file), "\n")
    cat("        Rows:", nrow(network), "\n")
    cat("        Expected source column:", source_col, "\n")
    cat("        Expected target column:", target_col, "\n")
    cat("        Expected score column:", ifelse(is.null(score_col), "None", score_col), "\n")
    cat("        Available columns:", paste(colnames(network), collapse = ", "), "\n")
  }
  
  # Check for required columns
  if (!(source_col %in% colnames(network))) {
    stop(paste("Source column '", source_col, "' not found in", basename(file), 
               ". Available columns: ", paste(colnames(network), collapse = ", ")))
  }
  
  if (!(target_col %in% colnames(network))) {
    stop(paste("Target column '", target_col, "' not found in", basename(file), 
               ". Available columns: ", paste(colnames(network), collapse = ", ")))
  }
  
  # Handle methods without score columns or missing score columns
  if (is.null(score_col) || !(score_col %in% colnames(network))) {
    if (verbose) cat("        No score column - treating all edges equally\n")
    
    network_long <- network %>%
      dplyr::select(dplyr::all_of(source_col), dplyr::all_of(target_col)) %>%
      dplyr::rename(Source = dplyr::all_of(source_col), Target = dplyr::all_of(target_col)) %>%
      dplyr::distinct()  # Remove duplicates
    
    network_long$Score <- 1  # Assign equal scores
    
    # For methods without scores, return all edges as "top" and sample for random
    n_total <- nrow(network_long)
    n_random <- max(1, min(round(random_percent * n_total), n_total))
    
    if (n_total > 0) {
      random_sample <- network_long %>% 
        dplyr::slice_sample(n = n_random, replace = FALSE)
    } else {
      random_sample <- network_long
    }
    
    return(list(top = network_long, random = random_sample))
  }
  
  # Process with score column
  network_long <- network %>%
    dplyr::select(dplyr::all_of(source_col), dplyr::all_of(target_col), dplyr::all_of(score_col)) %>%
    dplyr::rename(Source = dplyr::all_of(source_col), Target = dplyr::all_of(target_col), Score = dplyr::all_of(score_col)) %>%
    dplyr::distinct()  # Remove duplicates
  
  # Remove rows with NA or infinite scores
  initial_rows <- nrow(network_long)
  network_long <- network_long %>% 
    dplyr::filter(!is.na(Score) & is.finite(Score))
  
  if (verbose && nrow(network_long) < initial_rows) {
    cat("        Removed", initial_rows - nrow(network_long), "rows with invalid scores\n")
  }
  
  if (nrow(network_long) == 0) {
    stop(paste("No valid scores found in", basename(file)))
  }
  
  # Extract top and random percentages
  n_total <- nrow(network_long)
  n_top <- max(1, min(round(top_percent * n_total), n_total))
  n_random <- max(1, min(round(random_percent * n_total), n_total))
  
  if (verbose) {
    cat("        Selecting top", n_top, "and random", n_random, "edges from", n_total, "total\n")
  }
  
  top_subset <- network_long %>% 
    dplyr::arrange(dplyr::desc(Score)) %>% 
    dplyr::slice_head(n = n_top)
  
  random_subset <- network_long %>% 
    dplyr::slice_sample(n = n_random, replace = FALSE)
  
  return(list(top = top_subset, random = random_subset))
}

#' Calculate Jaccard Index - Internal
#' @keywords internal
calculate_jaccard_internal <- function(network1, network2) {
  # Create edge identifiers
  edges1 <- paste(network1$Source, network1$Target, sep = "->")
  edges2 <- paste(network2$Source, network2$Target, sep = "->")
  
  # Calculate intersection and union
  common_edges <- intersect(edges1, edges2)
  union_edges <- union(edges1, edges2)
  
  # Handle empty sets
  if (length(union_edges) == 0) return(NA)
  
  return(length(common_edges) / length(union_edges))
}

#' Compute Pairwise Jaccard Indices - Internal
#' @keywords internal
compute_pairwise_jaccard_internal <- function(network_list, verbose = FALSE) {
  ji_values <- c()
  n_networks <- length(network_list)
  
  if (n_networks < 2) {
    if (verbose) cat("        Less than 2 networks, skipping pairwise comparisons\n")
    return(ji_values)
  }
  
  # Compute all pairwise combinations
  for (i in 1:(n_networks - 1)) {
    for (j in (i + 1):n_networks) {
      ji <- calculate_jaccard_internal(network_list[[i]], network_list[[j]])
      ji_values <- c(ji_values, ji)
      
      if (verbose) {
        cat("        Pair", i, "vs", j, "- JI:", round(ji, 4), "\n")
      }
    }
  }
  
  return(ji_values)
}

#' Create Jaccard Index Heatmap - Internal
#' @keywords internal
create_jaccard_heatmap_internal <- function(median_ji_results) {
  
  # Create the heatmap
  ji_heatmap <- ggplot2::ggplot(median_ji_results, ggplot2::aes(x = Method, y = Sample, fill = Median_JI)) +
    ggplot2::geom_tile(color = "grey80", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Median_JI), 
                   color = ifelse(Median_JI > 0.5, "white", "black")),  # threshold at 0.5
      size = 4.2, fontface = "bold", show.legend = FALSE
    ) + 
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_viridis_c(name = "Median JI", option = "D", direction = -1, na.value = "white") +
    ggplot2::coord_fixed(ratio = 0.5) +  # fixes tile to square shape
    ggplot2::labs(
      title = "Median Jaccard Index Heatmap",
      x = "Method",
      y = "Sample"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = ggplot2::element_text(size = 13, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      panel.grid = ggplot2::element_blank()
    )
  
  return(ji_heatmap)
}


# UTILITY FUNCTIONS

#' Debug Network Files
#'
#' Helper function to inspect network files and their structure.
#'
#' @param base_dir Character. Base directory containing sample directories.
#' @param sample Character. Sample name to inspect.
#' @param method Character. Method name to inspect.
#' @param max_files Integer. Maximum number of files to inspect (default: 3).
#'
#' @export
debug_network_files <- function(base_dir, sample, method, max_files = 3) {
  method_path <- file.path(base_dir, sample, method)
  
  cat("=== DEBUGGING NETWORK FILES ===\n")
  cat("Sample:", sample, "\n")
  cat("Method:", method, "\n")
  cat("Directory:", method_path, "\n")
  
  if (!dir.exists(method_path)) {
    cat("ERROR: Directory does not exist!\n")
    return(invisible(NULL))
  }
  
  network_files <- list.files(method_path, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  
  cat("Files found:", length(network_files), "\n")
  
  if (length(network_files) == 0) {
    cat("ERROR: No network files found!\n")
    return(invisible(NULL))
  }
  
  # Default configurations
  default_delim <- list(
    CellOracle = ",", FigR = ",", TRIPOD = ",", Pando = ",",
    LINGER = ",", Scenic = "\t", GRaNIE = ",", DIRECTNET = ","
  )
  
  default_source_target_columns <- list(
    CellOracle = c("source", "target"),
    FigR = c("Motif", "DORC"),
    TRIPOD = c("TF", "gene"),
    Pando = c("tf", "target"),
    LINGER = c("Source", "Target"),
    Scenic = c("Source", "Target"),
    GRaNIE = c("TF.name", "gene.name"),
    DIRECTNET = c("TF", "Gene")
  )
  
  # Show file names
  cat("\nFiles:\n")
  for (i in seq_along(network_files)) {
    cat("  ", i, ":", basename(network_files[i]), "\n")
    if (i >= max_files) {
      cat("  ... (showing first", max_files, "files only)\n")
      break
    }
  }
  
  # Try to read the first file
  cat("\n--- Analyzing first file ---\n")
  cat("File:", basename(network_files[1]), "\n")
  
  expected_delim <- default_delim[[method]]
  expected_cols <- default_source_target_columns[[method]]
  
  cat("Expected delimiter:", ifelse(expected_delim == "\t", "TAB", expected_delim), "\n")
  cat("Expected source column:", expected_cols[1], "\n")
  cat("Expected target column:", expected_cols[2], "\n")
  
  # Try reading with expected delimiter
  tryCatch({
    first_file <- readr::read_delim(network_files[1], delim = expected_delim, 
                                    show_col_types = FALSE, n_max = 5)
    cat("\nSUCCESS: File read successfully!\n")
    cat("Dimensions:", nrow(first_file), "rows,", ncol(first_file), "columns\n")
    cat("Column names:", paste(colnames(first_file), collapse = ", "), "\n")
    
    # Check for expected columns
    source_found <- expected_cols[1] %in% colnames(first_file)
    target_found <- expected_cols[2] %in% colnames(first_file)
    
    cat("Source column found:", source_found, "\n")
    cat("Target column found:", target_found, "\n")
    
    if (!source_found || !target_found) {
      cat("\nERROR: Expected columns not found!\n")
      cat("This is likely the cause of your error.\n")
    }
    
    cat("\nFirst few rows:\n")
    print(first_file)
    
  }, error = function(e) {
    cat("\nERROR reading with expected delimiter:", e$message, "\n")
    
    # Try different delimiters
    cat("\nTrying different delimiters...\n")
    delims_to_try <- c(",", "\t", ";", " ", "|")
    names(delims_to_try) <- c("comma", "tab", "semicolon", "space", "pipe")
    
    for (delim_name in names(delims_to_try)) {
      tryCatch({
        test_read <- readr::read_delim(network_files[1], delim = delims_to_try[[delim_name]], 
                                       show_col_types = FALSE, n_max = 3)
        cat("  \u2713", delim_name, "delimiter works! Columns:", ncol(test_read), "\n")
        if (delim_name == "comma") {
          cat("    Column names:", paste(colnames(test_read), collapse = ", "), "\n")
        }
      }, error = function(e2) {
        cat("  \u2717", delim_name, "delimiter failed\n")
      })
    }
  })
  
  cat("\n=== END DEBUG ===\n")
  
  invisible(network_files)
}

#' Create Custom ggplot2 Theme
#'
#' Creates a custom theme for network stability plots.
#'
#' @param base_size Numeric. Base font size (default: 16).
#'
#' @return ggplot2 theme object.
#'
#' @import ggplot2
#'
#' @export
theme_network_stability <- function(base_size = 16) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.title.y = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = base_size - 2, face = "bold"),
      axis.text.y = ggplot2::element_text(size = base_size - 2, face = "bold"),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(15, 15, 15, 25, "pt")
    )
}
