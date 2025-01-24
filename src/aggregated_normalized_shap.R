## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Libraries
library(tidyverse)
library(dplyr)
library(purrr)
source("analysis/processate_metrics.R")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Input parameters
classifiers <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")

input_lists <- c("data_driven",
                 "dea",
                 "mrmr",
                 "disease_related",
                 "disease_related_entire_list",
                 "omnipath_data_driven",
                 "omnipath_disease_related",
                 "omnipath_intersection",
                 "guildify_disease_related",
                 "guildify_data_driven",
                 "omnipath_union")
threshold <- 0.75
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

# Filter cross-validation and test results with estimate > threshold
all_results_cross_validation <- all_results_cross_validation %>% filter(metric == metric_estimate) %>% filter(mean > 0.75)
all_results_test <- all_results_test %>% filter(metric == metric_estimate) %>% filter(estimate>0.75)

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

# saveRDS(normalized_shap, "../COPD/normalized_shap_all_combinations_list.Rds")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Aggregated shap by gene for each combination
aggregated_shap <- map(normalized_shap, function(res){
  res %>%
    group_by(variable_name) %>%
    summarize(mean_norm_shap = mean(abs(norm_shap)), .groups = "drop") %>%
    arrange(desc(mean_norm_shap)) %>% filter(mean_norm_shap !=0)
})

results_shap <- aggregated_shap

results_shap_dfr <- imap_dfr(results_shap, function(data, combination_name) {
  data %>%
    arrange(desc(mean_norm_shap)) %>%                  # Rank genes by mean_shap
    mutate(
      rank = dense_rank(desc(mean_norm_shap)),        # Rank within combination (tie-aware)
      total_genes = n(),                         # Total genes in this list
      norm_rank = 1 - (rank - 1) / (total_genes - 1), # Normalize rank
      input_classifier = combination_name        # Combination name (input and classifier)
    ) 
})

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. Aggregate normalized SHAP results across input list for obtaining the candidate list

# Use the normalized shap
norm_results <- results_shap_dfr %>% 
    group_by(variable_name) %>%
    summarize(
      aggregated_max = max(mean_norm_shap, na.rm =TRUE),
      aggregated_mean = mean(mean_norm_shap, na.rm = TRUE),
      aggregated_median = median(mean_norm_shap, na.rm = TRUE),
      sd = sd(mean_norm_shap, na.rm = TRUE),# Weight by total genes
      occurrence_count = n() # Number of combinations the gene appears in
    ) %>%
    arrange(desc(aggregated_max)) # Sort by aggregated rank

write.csv(norm_results, "../COPD/norm_results.csv", row.names = FALSE)

# Use the rank of the shap values
rank_results <- results_shap_dfr %>%
  group_by(variable_name) %>%
  summarize(
    aggregated_max = max(norm_rank, na.rm = TRUE),
    aggregated_mean = mean(norm_rank, na.rm = TRUE),
    aggregated_median = median(norm_rank, na.rm = TRUE),
    sd = sd(norm_rank, na.rm = TRUE),# Weight by total genes
    occurrence_count = n() # Number of combinations the gene appears in
  ) %>%
  arrange(desc(aggregated_max)) # Sort by aggregated rank

write.csv(rank_results, "../COPD/rank_results.csv", row.names = FALSE)


# ## USE THE RANK
# 
# # Step 1: Flatten and rank genes in each input_list-classifier combination
# ranked_results <- imap_dfr(results_shap, function(data, combination_name) {
#   data %>%
#     arrange(desc(mean_shap)) %>%                  # Rank genes by mean_shap
#     mutate(
#       rank = dense_rank(desc(mean_shap)),        # Rank within combination (tie-aware)
#       total_genes = n(),                         # Total genes in this list
#       normalized_rank = 1 - (rank - 1) / (total_genes - 1), # Normalize rank
#       input_classifier = combination_name        # Combination name (input and classifier)
#     ) 
# })
# 
# # Step 3: Assign final ranking (with frequency penalty)
# aggregated_ranks <- aggregated_ranks %>%
#   mutate(
#     final_rank = rank(-aggregated_rank)
#     # adjusted_rank = aggregated_rank / occurrence_count # Optional: penalize by occurrence
#   )
# 
# 
# # Step 2: Aggregate ranks across all combinations with weighting
# aggregated_ranks <- ranked_results %>%
#   group_by(variable_name) %>%
#   summarize(
#     aggregated_rank = median(normalized_rank, na.rm = TRUE),
#     sd = sd(normalized_rank, na.rm = TRUE),# Weight by total genes
#     occurrence_count = n() # Number of combinations the gene appears in
#   ) %>%
#   arrange(desc(aggregated_rank)) # Sort by aggregated rank
# 
# # Step 3: Assign final ranking (with frequency penalty)
# aggregated_ranks <- aggregated_ranks %>%
#   mutate(
#     final_rank = rank(-aggregated_rank)
#     # adjusted_rank = aggregated_rank / occurrence_count # Optional: penalize by occurrence
#   )
# 
# # View final results
# print(aggregated_ranks)
# 
# write.csv(aggregated_ranks, "../COPD/ranked_genes_from_shap_men.csv", row.names = FALSE)
# 
