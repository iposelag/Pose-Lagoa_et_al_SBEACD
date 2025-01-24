
## ----------------------------------------------------------------------------------------------------------------------------------------
# Load libraries
source("analysis/processate_metrics.R")
source("analysis/plots.R")
source("support/library_help.R")
source("feature_selection/library_feature_selection.R")
source("analysis/library_analysis.R")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. Processate metrics

classifiers <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")

input_lists <- c("data_driven",
                 "dea",
                 "mrmr",
                 "disease_related",
                 "disease_related_entire_list", 
                 "omnipath_data_driven", 
                 "omnipath_disease_related",
                 "guildify_data_driven",
                 "omnipath_intersection",
                 "guildify_disease_related", 
                 "guildify_data_driven",
                 "omnipath_union")
# input_lists <- c("data_driven", 
#                  "dea", 
#                  "mrmr", 
#                  "omnipath_data_driven", 
#                  "guildify_data_driven")
input_lists <- c("mrmr_30",
                 "disease_related",
                 "mrmr_76", "dea")

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

# Save the combined metrics to a CSV file if needed
write.csv(all_results_cross_validation, "../COPD/all_metrics_cross_validation.csv", row.names = FALSE)
write.csv(all_results_test, "../COPD/all_metrics_test.csv", row.names = FALSE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. Generate Radar Charts
# Read test and cross validation data
# test <- read.csv2("../COPD/all_metrics_test.csv", stringsAsFactors = FALSE, sep = ",")
# cross_validation <- read.csv2("../COPD/all_metrics_cross_validation.csv", stringsAsFactors = FALSE, sep = ",")
test <- all_results_test
cross_validation <- all_results_cross_validation

# Example usage:
test <- convert_input_column(test)
test <- convert_classifiers_column(test)
cross_validation <- convert_input_column(cross_validation)
cross_validation <- convert_classifiers_column(cross_validation)

metrics_to_plot <- c("accuracy")
output_dir <- "../COPD/plots"
minmax <- c(0.5,0.9,0.1)

# Loop through metrics
for (metric in metrics_to_plot) {

  # Test
  data_for_plot <- prepare_data_for_radarchart(test, metric, "estimate", minmax[1], minmax[2])
  plot_radarchart(data_for_plot, metric, "radarchart_test_mrm_comparison", output_dir, minmax[1], minmax[2], minmax[3])
  # summary_barplot(test, metric, ml_models_colors, paste0(output_dir,"/test_summary_barplot"), "estimate")

  # Cross validation
  data_for_plot <- prepare_data_for_radarchart(cross_validation, metric, "mean", minmax[1], minmax[2])
  plot_radarchart(data_for_plot, metric, "radarchart_cross_validation_mrmr_comparison", output_dir, minmax[1], minmax[2], minmax[3])
  # summary_barplot(cross_validation, metric, ml_models_colors, paste0(output_dir,"/cross_validation_summary_barplot"), "mean")

}

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. Anlysis of times
load("../COPD/results_mrmr/ML_models/expression_train.Rda")

input_lists <- c("dea", "mrmr", "data_driven", "guildify_data_driven", "omnipath_data_driven",
                 "disease_related", "guildify_disease_related", "omnipath_disease_related",
                 "disease_related_entire_list", "omnipath_union", "omnipath_intersection")
# Load gene lists
gene_lists <- lapply(input_lists, load_gene_list)
names(gene_lists) <- input_lists
# Intersect gene lists with genes from my dataset
gene_lists_intersection <- lapply(gene_lists, intersection_lists, expression=expression_train)
# Count the number of genes in each gene list
num_genes <- sapply(gene_lists_intersection, length)

# Create a new column in all_results_times for the number of genes
all_results_times$num_genes <- num_genes[match(all_results_times$input_list, names(num_genes))]
all_results_times <- convert_classifiers_column(all_results_times)
# Ensure all classifiers in all_results_times are present in ml_models_colors
all_results_times$classifier <- factor(all_results_times$classifier, levels = names(ml_models_colors))

# Plot comparing time vs number of genes across classifiers
p <- ggplot(all_results_times, aes(x = time, y = num_genes, color = classifier)) +
  geom_point() +
  geom_line(aes(group = classifier)) +
  labs(title = "Comparison of Time vs Number of Genes by Classifier",
       x = "Time (seconds)",
       y = "Number of Genes") +
  scale_color_manual(values = ml_models_colors) +
  theme_minimal()
# Save the plot
ggsave(filename = paste0(output_dir,"/times.pdf"), plot = p, device = "pdf")
