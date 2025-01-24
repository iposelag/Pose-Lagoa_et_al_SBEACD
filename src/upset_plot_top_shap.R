#!/bin/Rscript
#########################################################################################
# Upset plot with the top30 shap genes
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(dplyr)
library(purrr)
source("analysis/processate_metrics.R")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Input parameters
classifiers <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")
input <- "data_driven"
if(input != "disease_related"){
  input_lists <- c("data_driven",
                   "dea",
                   "mrmr",
                   "omnipath_data_driven",
                   "guildify_data_driven")
}else{
  input_lists <- c("disease_related",
                       "disease_related_entire_list",
                       "omnipath_disease_related",
                       "guildify_disease_related")}

threshold <- 0
metric_estimate <- "normMCC"

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. Filter input list-classifier combinations that are greater than a given threshold
# Initialize an empty data frame to store all metrics
all_results_cross_validation <- data.frame()
all_results_test <- data.frame()
all_results_times <- data.frame()

# Loop through each results_models.rds file and process it
for (input_list in input_lists) {
  cat("Adding:", input_list, "\n")
  file_path <- paste0("../COPD/results_", input_list,"/results_models.rds")
  # Read the results_models object
  results_models <- readRDS(file_path)
  
  # Process cross-validation and test metrics
  cross_validation_metrics <- process_cross_validation_metrics(results_models, classifiers)
  test_metrics <- process_test_metrics(results_models, classifiers)
  
  # Process time of classifiers
  classifier_times <- process_times(results_models) 
  
  # Add input_name column to each metrics data frame
  cross_validation_metrics$input_list <- input_list
  test_metrics$input_list <- input_list
  classifier_times$input_list <- input_list
  
  # Append to the overall results
  all_results_cross_validation <- rbind(all_results_cross_validation, cross_validation_metrics)
  all_results_test<- rbind(all_results_test, test_metrics)
  all_results_times <- rbind(all_results_times, classifier_times)
}

# Extract pairs input list-classifier
unique_pairs <- bind_rows(
  all_results_cross_validation %>%
    select(input_list, classifier),
  all_results_test %>%
    select(input_list, classifier)
) %>%
  distinct()

# Actualize input lists present in the pairs
input_lists <- unique(unique_pairs$input_list)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. Load shap data
# Iterate through input lists
results <- map(input_lists, function(input_name) {
  # Read the SHAP results for the given input
  shap_data <- readRDS(paste0("../COPD/results_", input_name, "/results_shap.rds"))
})

# Combine results into a single list by classifier and input
names(results) <- input_lists

# Define a function to filter results based on the unique pairs
filter_results <- function(results, unique_pairs) {
  filtered_results <- list()
  
  for (i in seq_len(nrow(unique_pairs))) {
    # Extract the input_list and classifier for the current pair
    input_list <- unique_pairs$input_list[i]
    classifier <- unique_pairs$classifier[i]
    
    # Check if the pair exists in the results list
    if (!is.null(results[[input_list]][[classifier]])) {
      # Add the SHAP values to the filtered_results list
      filtered_results[[paste(input_list, classifier, sep = "_")]] <- 
        results[[input_list]][[classifier]]$shap_values
    }
  }
  
  return(filtered_results)
}

# Apply the function to filter the results
filtered_results <- filter_results(results, unique_pairs)

# Check the filtered results
print(names(filtered_results))

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. Normalize SHAP values by sample for each combination
normalized_shap <- map(filtered_results, function(res) {
  res %>%
    group_by(sample) %>%
    mutate(
      norm_shap = mean_shap / sum(abs(mean_shap)), # Normalize SHAP values
    ) %>%
    ungroup() %>%
    arrange(desc(norm_shap)) %>% # Optional: sort by normalized_shap
    filter(norm_shap != 0)       # Remove zero normalized SHAP values
})

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Aggregated shap by gene for each combination
aggregated_shap <- map(normalized_shap, function(res){
  res %>%
    group_by(variable_name) %>%
    summarize(mean_norm_shap = mean(abs(norm_shap)), .groups = "drop") %>%
    arrange(desc(mean_norm_shap)) %>% filter(mean_norm_shap !=0)
})

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. Top 30 genes 
top_genes <- map(aggregated_shap, function(res){
  res %>%
    slice_max(mean_norm_shap, n=30) %>% select(variable_name)
})

# Combine all `variable_name` columns into a single vector
all_genes <- unlist(lapply(top_genes, function(df) df$variable_name))
# Count occurrences of each gene
gene_counts <- sort(table(all_genes))
print(gene_counts)
saveRDS(gene_counts, paste0("../COPD/gene_counts_",input,".Rds"))

# Extract the list of genes from the tibbles
gene_list <- lapply(top_genes, function(tibble) tibble$variable_name)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot

# Convert the list of genes into a presence/absence matrix
upset_data <- fromList(gene_list)
comb_mat <- make_comb_mat(gene_list)
comb_degree(comb_mat)

# Open a PDF device
png(paste0("../COPD/plots/upset_plot_",input,".png"), width = 14, height = 8)
upset(
  upset_data, 
  sets = colnames(upset_data), # All columns except the "gene" column
  nintersects = NA,            # Display all intersections
  keep.order = TRUE,            # Keep the order of sets as in the data
  order.by = "freq",           # Order intersections by frequency
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.2, 1.2) # Scale text size (axes, labels, etc.)
)
dev.off()

names <- map(top_genes, function(res){
  res %>%
    filter(variable_name == "TGFBR3")
})