# Required libraries
library(dplyr)
library(ROCR)
library(PRROC)

# Load your ground truth and inferred networks
ground_truth <- read.csv("filtered_ground_truth_56TFs_3036TGs.csv", header = TRUE)
inferred_network <- read.csv("../output_GRN/5000cells_E7.5_rep1_final_GRN.csv", header = TRUE)

# Prepare the ground truth and inferred network data
ground_truth$interaction <- paste(tolower(ground_truth$Source), tolower(ground_truth$Target), sep = "_")
inferred_network$interaction <- paste(tolower(inferred_network$source), tolower(inferred_network$target), sep = "_")

# Step 1: Subsample ground truth to match inferred network size
set.seed(123)  # For reproducibility
ground_truth_sampled <- ground_truth %>%
  slice_sample(n = nrow(inferred_network))

# Step 2: Merge inferred network with sampled ground truth
merged_network <- inferred_network %>%
  mutate(true_interaction = ifelse(interaction %in% ground_truth_sampled$interaction, 1, 0))

# Step 3: Calculate Accuracy
accuracy <- mean(merged_network$true_interaction == 1)
cat("Accuracy:", accuracy, "\n")

# Step 4: Calculate AUROC
merged_network$predicted_scores <- abs(merged_network$coef_mean)  # Use absolute values for AUROC
pred <- prediction(merged_network$predicted_scores, merged_network$true_interaction)
perf <- performance(pred, "tpr", "fpr")
auroc <- performance(pred, "auc")
cat("AUROC:", auroc@y.values[[1]], "\n")

# Step 5: Calculate AUPRC
pr <- pr.curve(scores.class0 = merged_network$predicted_scores[merged_network$true_interaction == 1],
               scores.class1 = merged_network$predicted_scores[merged_network$true_interaction == 0],
               curve = TRUE)
cat("AUPRC:", pr$auc.integral, "\n")

# Step 6: Save ROC and PR curves
png("figures/1000cells_ROC_curve1.png")
plot(perf, colorize = TRUE, main = "ROC Curve")
dev.off()

png("figures/1000_cells_PR_curve1.png")
plot(pr, main = "Precision-Recall Curve")
dev.off()

# Define metric functions
calculate_f1 <- function(tp, fp, fn) {
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  return(f1)
}

calculate_ji <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

calculate_wjs <- function(ground_truth, inferred_network) {
  common_edges <- intersect(ground_truth$interaction, inferred_network$interaction)
  weighted_common <- inferred_network %>%
    filter(interaction %in% common_edges) %>%
    summarize(weight_sum = sum(abs(coef_mean))) %>%
    pull(weight_sum)
  total_weight <- sum(abs(inferred_network$coef_mean))
  wjs <- weighted_common / total_weight
  return(wjs)
}

calculate_epr <- function(top_network) {
  tp <- sum(top_network$true_interaction == 1)
  total <- nrow(top_network)
  epr <- tp / total
  return(epr)
}

# Modified top_edges function
top_edges <- function(top_k) {
  top_network <- merged_network %>%
    arrange(desc(abs(coef_mean))) %>%
    head(top_k)
  
  tp <- sum(top_network$true_interaction == 1)
  fp <- sum(top_network$true_interaction == 0)
  fn <- sum(ground_truth_sampled$interaction %in% top_network$interaction == FALSE)
  
  f1_score <- calculate_f1(tp, fp, fn)
  ji <- calculate_ji(ground_truth_sampled$interaction, top_network$interaction)
  wjs <- calculate_wjs(ground_truth_sampled, top_network)
  epr <- calculate_epr(top_network)
  
  return(list(f1_score = f1_score, ji = ji, wjs = wjs, epr = epr))
}

# Evaluate for top edges and save results
edge_counts <- c(1000, 2000, 3000, 4000, 5000)
results <- data.frame()

for (k in edge_counts) {
  res <- top_edges(k)
  results <- rbind(results, data.frame(
    Top_K = k,
    F1_Score = res$f1_score,
    Jaccard_Index = res$ji,
    Weighted_Jaccard_Similarity = res$wjs,
    Early_Precision_Rate = res$epr
  ))
}

write.csv(results, "Top_edges_5000cells_evaluation_filtered.csv", row.names = FALSE)
cat("Top edge evaluation completed and saved.\n")
