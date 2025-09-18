# Load required libraries
library(dplyr)

# Get the list of input files
file_list <- c(
  "xymats.tX1.neg.df.csv", "xymats.tX1.pos.df.csv", 
  "xymats.tX2.neg.df.csv", "xymats.tX2.pos.df.csv", 
  "xymats.tY1.neg.df.csv", "xymats.tY1.pos.df.csv", 
  "xymats.tY2.neg.df.csv", "xymats.tY2.pos.df.csv"
)

# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop through each file and read the data
for (file in file_list) {
  temp_data <- read.csv(file)
  combined_data <- bind_rows(combined_data, temp_data)
}

# Select only relevant columns (gene, TF, and coef)
filtered_data <- combined_data %>%
  select(gene, TF, coef)

# For each gene-TF pair, select the entry with the highest absolute coef
result <- filtered_data %>%
  mutate(abs_coef = abs(coef)) %>%
  group_by(gene, TF) %>%
  slice_max(order_by = abs_coef, n = 1) %>%
  ungroup() %>%
  select(gene, TF, abs_coef)

# Write the result to a CSV file
write.csv(result, "gene_TF_highest_abs_coef.csv", row.names = FALSE)

