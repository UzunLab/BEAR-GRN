# Required libraries
library(dplyr)
library(ROCR)
library(PRROC)

# Load your ground truth and inferred networks
ground_truth <- read.delim("ground_truth/RN111.tsv", header = TRUE)
inferred_network <- read.csv("output_GRN/1000cells_E7.5_rep1_final_GRN.csv", header = TRUE)
dim(inferred_network)
dim(ground_truth)
inferred_network
# Step 1: Prepare the ground truth data
ground_truth$interaction <- paste(tolower(ground_truth$Source), tolower(ground_truth$Target), sep = "_")

# Step 2: Prepare the inferred network data
inferred_network$interaction <- paste(tolower(inferred_network$source), tolower(inferred_network$target), sep = "_")
ground_truth
inferred_network
# Step 3: Merge inferred network with ground truth
merged_network <- inferred_network %>%
  mutate(true_interaction = ifelse(interaction %in% ground_truth$interaction, 1, 0))
merged_network
# Step 4: Balance the dataset by equalizing positive and negative samples
positive_samples <- merged_network %>% filter(true_interaction == 1)
negative_samples <- merged_network %>% filter(true_interaction == 0)
dim(positive_samples)
dim(negative_samples)
# Randomly sample the same number of negatives as positives
set.seed(123)  # For reproducibility
############################
#sampled_negatives <- negative_samples %>% sample_n(nrow(positive_samples))

# Combine the balanced positives and negatives
#balanced_network <- bind_rows(positive_samples, sampled_negatives)
datasize=2*(nrow(positive_samples))
datasize
balanced_network <- sample_n(merged_network, size = datasize)
dim(balanced_network)
balanced_network
############################


# Step 5: Calculate Accuracy
accuracy <- mean(balanced_network$true_interaction == 1)
cat("Accuracy:", accuracy, "\n")

# Step 6: Calculate AUROC
balanced_network$predicted_scores <- abs(balanced_network$coef_mean)  # Use absolute values for AUROC
pred <- prediction(balanced_network$predicted_scores, balanced_network$true_interaction)
perf <- performance(pred, "tpr", "fpr")
auroc <- performance(pred, "auc")
cat("AUROC:", auroc@y.values[[1]], "\n")

# Step 7: Calculate AUPRC
pr <- pr.curve(scores.class0 = balanced_network$true_interaction,
               scores.class1 = 1 - balanced_network$true_interaction,
               curve = TRUE)
cat("AUPRC:", pr$auc.integral, "\n")

# Step 8: Save ROC and PR curves
png("output/ROC_curve1.png")
plot(perf, colorize = TRUE, main = "ROC Curve")
dev.off()

png("output/PR_curve1.png")
plot(pr, main = "Precision-Recall Curve")
dev.off()

# Step 9: Evaluation metrics for top edges
# Reusing the previous functions for F1, Jaccard, WJS, and EPR
# Function to calculate F1-Score
calculate_f1 <- function(tp, fp, fn) {
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  return(f1)
}

# Function to calculate Jaccard Index (JI)
calculate_ji <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Function to calculate Weighted Jaccard Similarity (WJS)
calculate_wjs <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  return(intersection / min(length(set1), length(set2)))
}

# Function to calculate Early Precision Rate (EPR)
calculate_epr <- function(top_network) {
  tp <- sum(top_network$true_interaction == 1)
  total <- nrow(top_network)
  epr <- tp / total
  return(epr)
}

# Top edge evaluations
top_edges <- function(top_k) {
  top_network <- balanced_network %>%
    arrange(desc(abs(coef_mean))) %>%  # Rank by absolute coef_mean
    head(top_k)

  tp <- sum(top_network$true_interaction == 1)
  fp <- sum(top_network$true_interaction == 0)
  fn <- sum(ground_truth$interaction %in% top_network$interaction == FALSE)

  # Calculate F1-score
  f1_score <- calculate_f1(tp, fp, fn)

  # Calculate Jaccard Index (JI)
  ji <- calculate_ji(ground_truth$interaction, top_network$interaction)

  # Calculate Weighted Jaccard Similarity (WJS)
  wjs <- calculate_wjs(ground_truth$interaction, top_network$interaction)

  # Calculate Early Precision Rate (EPR)
  epr <- calculate_epr(top_network)

  return(list(f1_score = f1_score, ji = ji, wjs = wjs, epr = epr))
}

# Evaluate for top 1K, 2K, 3K, 5K edges
edge_counts <- c(1000, 2000, 3000, 5000)
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

# Save the results
write.csv(results, "output/Top_edges_evaluation1.csv", row.names = FALSE)

cat("Top edge evaluation completed and saved.\n")

