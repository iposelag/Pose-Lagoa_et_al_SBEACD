#########################################################################################
# Plot input genes lists
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
source("feature_selection/library_feature_selection.R")
source("support/library_help.R")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Input parameters
input_lists<- c("dea", "mrmr", "data_driven", "guildify_data_driven", "omnipath_data_driven",
                "disease_related", "guildify_disease_related", "omnipath_disease_related",
                "disease_related_entire_list", "omnipath_intersection", "omnipath_union")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load expression data
load("../COPD/results_mrmr/ML_models/expression_train.Rda")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Prepare data for plotting
# Initialize a data frame to store results
fsm_genes <- data.frame(
  input_list = character(),
  intersection_count = numeric(),
  stringsAsFactors = FALSE
)

# Loop through input lists and calculate intersections
for (genes_list in input_lists) {
  cat(genes_list, "\n")
  # Load the gene list for the current input
  current_genes <- load_gene_list(genes_list)
  # Calculate intersection
  intersection <- length(intersect(colnames(expression_train), current_genes))
  cat("Intersection:", intersection, "\n")
  # Append results to fsm_genes
  fsm_genes <- rbind(fsm_genes, data.frame(input_list = genes_list, intersection_count = intersection))
  cat("---\n")
}

# View the final data frame
fsm_genes <- fsm_genes %>% convert_input_column()
print(fsm_genes)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot
output_file <- "../COPD/plots/input_lists_barplot.pdf"
barplot_feature_selection <- function(feature_selection_data, output_file){
  
  p <-ggplot(feature_selection_data) +
    aes(
      x = input_list,
      fill = input_list,
      weight = intersection_count
    ) +
    labs(x = "Input selection")+labs(y = "")+
    geom_bar(position = "dodge") +
    geom_text(aes(label = intersection_count, y = intersection_count), 
              position = position_dodge(width = 0.9), vjust = -0.5) +  # Add labels on top of bars
    scale_fill_manual(
      values = genes_list_colors
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.title.x = element_blank()) +  # Remove x-axis title 
    labs(
      title = "Feature Selection Methods Number of Genes"
    )+
    labs(
      fill = "Feature Selection Approach"  # Change this to your desired legend title
    )
  
  # Save the plot
  ggsave(filename = output_file, plot = p, device = "pdf", width = 14, height = 8)
  
}
barplot_feature_selection(fsm_genes, output_file)
