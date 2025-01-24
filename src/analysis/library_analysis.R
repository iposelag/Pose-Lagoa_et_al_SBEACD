
## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to compute the intersection among lists
intersection_lists <- function(gene_list, expression){
  # Calculate intersection
  intersection <- intersect(colnames(expression), gene_list)
  return(intersection)
}


## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to categorize genes based on gene list memberships
categorize_genes <- function(ranked_data, gene_lists) {
  # Create a long format data frame of all gene list memberships
  membership_data <- do.call(rbind, lapply(names(gene_lists), function(list_name) {
    data.frame(variable_name = gene_lists[[list_name]], list = list_name, stringsAsFactors = FALSE)
  }))
  
  pairwise_intersections <- combn(names(gene_lists), 2, function(x) {
    list1 <- gene_lists[[x[1]]]
    list2 <- gene_lists[[x[2]]]
    intersect_genes <- intersect(list1, list2)
    list(pair = paste(x[1], x[2], sep = " & "), intersection = intersect_genes, count = length(intersect_genes))
  }, simplify = FALSE)
  
  # Combine membership and intersection data
  all_gene_pairs <- bind_rows(
    membership_data %>% mutate(pair = list),
    pairwise_intersections
  )
  
  # Deduplicate and consolidate into a single string for each gene
  consolidated_pairs <- all_gene_pairs %>%
    group_by(variable_name) %>%
    summarise(pair = paste(unique(pair), collapse = " & "), .groups = "drop")
  
  # Merge consolidated pairs with ranked data
  ranked_data <- ranked_data %>%
    left_join(consolidated_pairs, by = "variable_name")
  
  return(ranked_data)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Process shap value for any classifier (one shap value per gene)
process_classifier <- function(classifier_name, shap_data, top_n = 50) {
  shap_data[[classifier_name]][["shap_values"]] %>%
    group_by(variable_name) %>%
    summarize(mean_shap = mean(abs(mean_shap)), .groups = "drop") %>%
    arrange(desc(mean_shap))
  # slice_head(n = top_n)  # Select top N genes
}

