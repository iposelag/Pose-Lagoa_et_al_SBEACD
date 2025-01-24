#########################################################################################
# Venn diagram and upset seed lists
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
source("analysis/library_analysis.R")
source("feature_selection/library_feature_selection.R")
source("support/library_help.R")
library(UpSetR)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load input parameters
# Input parameters
input_lists <- c("dea", "mrmr", "data_driven", "guildify_data_driven", "omnipath_data_driven",
                 "disease_related", "guildify_disease_related", "omnipath_disease_related",
                 "disease_related_entire_list")
input_lists <- c("dea", "mrmr", "disease_related",
                 "disease_related_entire_list")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

# Load expression data
load("../COPD/results_mrmr/ML_models/expression_train.Rda")

# Load all gene lists into a named list
gene_lists <- lapply(input_lists, load_gene_list)
names(gene_lists) <- input_lists
gene_lists <- lapply(gene_lists, intersection_lists, expression=expression_train)

# Compute pairwise intersections
pairwise_intersections <- combn(names(gene_lists), 2, function(x) {
  list1 <- gene_lists[[x[1]]]
  list2 <- gene_lists[[x[2]]]
  intersect_genes <- intersect(list1, list2)
  list(pair = paste(x[1], x[2], sep = " & "), intersection = intersect_genes, count = length(intersect_genes))
}, simplify = FALSE)

# Convert pairwise intersections to a data frame
pairwise_df <- do.call(rbind, lapply(pairwise_intersections, function(x) {
  data.frame(pair = x$pair, count = x$count, stringsAsFactors = FALSE)
}))

# Print results
print(pairwise_df)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot upset
upset(fromList(gene_lists_intersection), nsets = 9, nintersects = 55 )

# Function to rename the gene lists
rename_gene_lists <- function(gene_lists) {
  # Mapping input_list names to more descriptive names
  new_names <- case_when(
    names(gene_lists) == "dea" ~ "DEA",
    names(gene_lists) == "mrmr" ~ "mRMR",
    names(gene_lists) == "disease_related" ~ "COPD-related curated",
    names(gene_lists) == "disease_related_entire_list" ~ "COPD-related entire list",
    TRUE ~ names(gene_lists)  # Retain original names if no match
  )
  
  # Assign the new names to the gene_lists object
  names(gene_lists) <- new_names
  return(gene_lists)
}

gene_lists <- rename_gene_lists(gene_lists)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot Venn Diagram
# Create a Venn diagram
venn_plot <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = NULL, # Set to NULL to render to the R device
  output = TRUE,
  imagetype = "png",
  fill = genes_list_colors[names(gene_lists)],
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Venn Diagram of Gene Lists"
)

# Plot
grid::grid.draw(venn_plot)
