# Required libraries
library(dplyr)
library(ROCR)
library(PRROC)

# Load your ground truth and inferred networks
ground_truth <- read.csv("ground_truth/filtered_ground_truth_67TFs_8607TGs.csv", header = TRUE)
inferred_network <- read.csv("ground_truth/1000cells_filtered_inferred_grn_67TFs_8607TGs.csv", header = TRUE)
dim(inferred_network)
dim(ground_truth)
#View(inferred_network)
# Step 1: Prepare the ground truth data
ground_truth$interaction <- paste(tolower(ground_truth$Source), tolower(ground_truth$Target), sep = "_")
#View(ground_truth)
# Step 2: Prepare the inferred network data
inferred_network$interaction <- paste(tolower(inferred_network$source), tolower(inferred_network$target), sep = "_")
#View(ground_truth)
inferred_network
# Step 3: Merge inferred network with ground truth
merged_network <- inferred_network %>%
  mutate(true_interaction = ifelse(interaction %in% ground_truth$interaction, 1, 0))
##########################################################
library(ggplot2)
merged_network1 <- merged_network
# Convert true_interaction to a factor for better visualization
merged_network1$true_interaction <- as.factor(merged_network1$true_interaction)
merged_network1$predicted_scores <- merged_network1$coef_abs
# Create the boxplot
ggplot(merged_network1, aes(x = true_interaction, y = predicted_scores, fill = true_interaction)) +
  geom_boxplot() +
  labs(x = "True Interaction", y = "Predicted Scores", title = "Boxplot of Predicted Scores by True Interaction") +
  theme_minimal()

merged_network$coef_abs
###########################################################
#View(merged_network)
# Step 4: Balance the dataset by equalizing positive and negative samples
#positive_samples <- merged_network %>% filter(true_interaction == 1)
#negative_samples <- merged_network %>% filter(true_interaction == 0)
#dim(positive_samples)
#dim(negative_samples)
# Randomly sample the same number of negatives as positives
set.seed(123)  # For reproducibility

# Step 5: Calculate Accuracy
accuracy <- mean(merged_network$true_interaction == 1)
cat("Accuracy:", accuracy, "\n")
# Step 6: Calculate AUROC
merged_network$predicted_scores <- abs(merged_network$coef_mean)  # Use absolute values for AUROC
pred <- prediction(merged_network$predicted_scores, merged_network$true_interaction)
perf <- performance(pred, "tpr", "fpr")
auroc <- performance(pred, "auc")
cat("AUROC:", auroc@y.values[[1]], "\n")


# Step 7: Calculate AUPRC

pr <- pr.curve(scores.class0 = merged_network$true_interaction,
               scores.class1 = 1 - merged_network$true_interaction,
               curve = TRUE)
cat("AUPRC:", pr$auc.integral, "\n")


# Step 8: Save ROC and PR curves
png("../figures/1000cells_ROC_curve_filtered.png")
plot(perf, colorize = TRUE, main = "ROC Curve")
dev.off()

png("../figures/1000_cells_PR_curve_filtered.png")
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
calculate_wjs <- function(ground_truth, inferred_network) {
  # Merge ground truth and inferred network to find common edges
  common_edges <- intersect(ground_truth$interaction, inferred_network$interaction)
  
  # Calculate the sum of weights for common edges
  weighted_common <- inferred_network %>%
    filter(interaction %in% common_edges) %>%
    summarize(weight_sum = sum(abs(coef_mean))) %>%
    pull(weight_sum)
  
  # Calculate the sum of weights for all edges in the inferred network (union of edges)
  total_weight <- sum(abs(inferred_network$coef_mean))
  
  # Calculate WJS
  wjs <- weighted_common / total_weight
  return(wjs)
}



# Function to calculate Early Precision Rate (EPR)
calculate_epr <- function(top_network) {
  tp <- sum(top_network$true_interaction == 1)
  total <- nrow(top_network)
  epr <- tp / total
  return(epr)
}

# Modify top_edges function to use the updated WJS calculation
top_edges <- function(top_k) {
  top_network <- merged_network %>%
    arrange(desc(abs(coef_mean))) %>%  # Rank by absolute coef_mean
    head(top_k)
  
  tp <- sum(top_network$true_interaction == 1)
  fp <- sum(top_network$true_interaction == 0)
  fn <- sum(ground_truth$interaction %in% top_network$interaction == FALSE)
  
  # Calculate F1-score
  f1_score <- calculate_f1(tp, fp, fn)
  
  # Calculate Jaccard Index (JI)
  ji <- calculate_ji(ground_truth$interaction, top_network$interaction)
  
  # Calculate Weighted Jaccard Similarity (WJS) with modified function
  wjs <- calculate_wjs(ground_truth, top_network)
  
  # Calculate Early Precision Rate (EPR)
  epr <- calculate_epr(top_network)
  
  return(list(f1_score = f1_score, ji = ji, wjs = wjs, epr = epr))
}

# Evaluate for top 1K, 2K, 3K, 5K edges
edge_counts <- c(10000, 20000, 30000, 40000, 50000)
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
View(results)
# Save the results
write.csv(results, "Top_edges_1000cells_evaluation_filtered.csv", row.names = FALSE)

cat("Top edge evaluation completed and saved.\n")

