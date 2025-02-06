source("support/library_help.R")
source("analysis/library_analysis.R")
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load expression data
load("../COPD/results_mrmr/ML_models/expression_train.Rda")
expression <- expression_train
# Load phenotypic data
load("../COPD/raw_data/phenotypic.Rda")

# Load normalized shap data
norm <- read.table("../COPD/results/norm_results.csv", header=T, sep=",", stringsAsFactors = F)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Prepare data for plotting 

# Filter normalized shap data
norm_filtered <- norm %>% filter(aggregated_max>0.0109) 
write.csv(norm_filtered, "../COPD/candidate_genes.csv")

top_genes <- norm_filtered$variable_name

# ---------------------
# Prepare expression data
# expression_df_filtered <- expression %>% select(all_of(top_genes))
# expression_spread_matrix <- as.matrix(expression_df_filtered)
# samples_expression <- rownames(expression_spread_matrix)

# ---------------------
# Prepare shap data: aggregated of all normalized scenarios matrix 
normalized_shap <- readRDS("../COPD/results/normalized_shap_all_combinations_list.Rds")

results_shap_dfr <- imap_dfr(normalized_shap, function(data, combination_name) {
  data %>%
    arrange(desc(norm_shap)) %>%                  # Rank genes by mean_shap
    mutate(
      rank = dense_rank(desc(norm_shap)),        # Rank within combination (tie-aware)
      total_genes = n(),                         # Total genes in this list
      norm_rank = 1 - (rank - 1) / (total_genes - 1), # Normalize rank
      input_classifier = combination_name        # Combination name (input and classifier)
    ) 
})

results_shap_dfr <- results_shap_dfr %>% group_by(variable_name) %>% group_by(sample)
# Compute the mean norm_shap for each gene and sample
mean_norm_shap <- results_shap_dfr %>%
  group_by(variable_name, sample) %>%
  summarise(mean_norm_shap = mean(norm_shap, na.rm = TRUE), .groups = "drop")

# Check the resulting data frame
mean_norm_shap <- mean_norm_shap %>% filter(variable_name %in% top_genes)

# Reshape the data into a wide format
heatmap_norm_shap <- mean_norm_shap %>%
  pivot_wider(names_from = sample, values_from = mean_norm_shap, values_fill = 0)
samples_norm_shap <- colnames(heatmap_norm_shap)[-1]

samples_norm_shap <- gsub("^(COPD_|CTRL_)", "", samples_norm_shap)

# Extract the matrix for the heatmap
heatmap_matrix <- as.matrix(heatmap_norm_shap[,-1]) # Remove the first column (gene names) for the matrix
rownames(heatmap_matrix) <- heatmap_norm_shap$variable_name # Set gene names as row names
colnames(heatmap_matrix) <- samples_norm_shap

# ---------------------
# Prepare annotation data for plotting
pdata_interest <- c("dis_condition","GOLD_stage", "sex", "smoker", "age", 
                     "emphysema", "pred_dlco")
phenotypic_ctics_interest <- phenotypic_ctics[, pdata_interest]
phenotypic_ctics_interest <- phenotypic_ctics_interest %>%
  mutate(dis_condition = factor(dis_condition),
         sex = factor(sex),
         # pneumocystis_colonization = factor(pneumocystis_colonization),
         smoker = factor(smoker),
         GOLD_stage = factor(GOLD_stage),
         # platform_id = factor(platform_id),
         age = as.numeric(age),
         emphysema = as.numeric(emphysema),
         pred_dlco = as.numeric(pred_dlco))
phenotypic_ctics_interest$sample_id <- rownames(phenotypic_ctics_interest)
# Filter samples
phenotypic_ctics_interest <- subset(phenotypic_ctics_interest, sample_id %in%
                                      samples_norm_shap)
# Rename sample_id to sample in phenotypic_ctics_interest_train
phenotypic_ctics_interest <- phenotypic_ctics_interest %>%
  rename(sample_id = sample_id)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot
spread_matrix <- heatmap_matrix
samples <- samples_norm_shap
# Annotations plot
# Colors of the annotations
annotation_colors <- list(GOLD_stage = GOLD_stage, dis_condition = dis_condition, sex = sex, smoker = smoker, age = age )
# Match samples
annotation <- phenotypic_ctics_interest[match(samples, rownames(phenotypic_ctics_interest)),]
# Annotations heatmap
ha <- HeatmapAnnotation(
  disease_condition = annotation$dis_condition, 
  GOLD_stage = annotation$GOLD_stage,
  sex = annotation$sex,
  smoker = annotation$smoker,
  age  = annotation$age,
  emphysema = annotation$emphysema,
  pred_DLCO = annotation$pred_dlco,
  col = list(GOLD_stage = GOLD_stage, 
             disease_condition = dis_condition,
             sex = sex,
             # age = age,
             smoker = smoker),
  annotation_name_side = "left",
  simple_anno_size = unit(2, "mm"),
  annotation_name_gp= gpar(fontsize = 6),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 6), 
                                 title_gp = gpar(fontsize = 8))
)

column_dist <- dist(t(spread_matrix))      # Compute distance matrix for columns
column_hclust <- hclust(column_dist)       # Perform hierarchical clustering
# Cut the dendrogram into desired number of clusters (e.g., 5 clusters)
column_clusters <- cutree(column_hclust, k = 5)  # k is the number of clusters
# Split the sample names by their cluster assignment
clusters <- split(names(column_clusters), column_clusters)
# View the results
clusters

# Manually recreate the default palette
# default_colors <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
default_colors <- colorRamp2(seq(min(t(scale(t(heatmap_matrix)))), max(t(scale(t(heatmap_matrix)))), length = 3), c("blue", "#EEEEEE", "red"))
p1 <- Heatmap(
  t(scale(t(heatmap_matrix))),
  na_col = "white",
  col = default_colors,
  cluster_columns = column_hclust,  # Pass the precomputed clustering
  column_split = 5,  # Optional: split columns into 5 clusters
  row_names_gp = gpar(fontsize = 3.8),
  column_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title = "shap value",
    title_gp = gpar(fontsize = 8)
  )
)

p1
# Open a PDF device
pdf("../COPD/plots/heatmap_candidate_genes.pdf", width = 14, height = 14) # Adjust width and height as needed

# Draw the heatmap
draw(p1)

# Close the device
dev.off()

# PV clust significances of clusters
# install.packages("pvclust")
# library(pvclust)
# pv_result <- pvclust(
#   spread_matrix, # Transposed matrix for clustering samples
#   method.hclust = "complete", # Use the same method as `hclust`
#   method.dist = "euclidean", # Use the same distance metric as `dist`
#   nboot = 1000 # Number of bootstrap replicates (increase for better accuracy)
# )
# 
# pdf("../COPD/plots/pv_clust.pdf", width = 16, height = 10) # Adjust width and height as needed
# plot(pv_result, cex = 0.4, cex.pv = 0.6)
# pvrect(pv_result, alpha = 0.9)
# dev.off()

## ----------------------------------------------------------------------------------------------------------------------------------------
# Enrichment of clusters

get_unique_labels <- function(row) {
  unique(colnames(observed)[row])
}

# phenotypic_ctics_interest <- phenotypic_ctics_interest %>% filter(dis_condition != "CTRL")

# Initialize a dataframe to store results
results <- data.frame(
  Variable = character(),
  Test = character(),
  Statistic = numeric(),
  p.value = numeric(),
  Enriched_Class = character(),
  cluster = character(),
  stringsAsFactors = FALSE
)
# Iterate over each cluster
for(cluster in 1:5){
  cluster_samples <- clusters[[cluster]]
  phenotypic_cluster <- phenotypic_ctics_interest
  phenotypic_cluster$cluster <- ifelse(phenotypic_cluster$sample_id %in%
                                         cluster_samples, 1, 0)
  
  pdata <- phenotypic_cluster
  
  # Loop through each column excluding the last two
  for(i in 1:(ncol(pdata)-2)){
    print("-------------------")
    print(names(pdata[i]))
    
    if(class(pdata[[i]]) == "numeric"){
      normality <- by(pdata[[i]], pdata$cluster, shapiro.test)
      if(normality[[1]]$p.value < 0.05 | normality[[2]]$p.value < 0.05){
        test <- wilcox.test(pdata[[i]]~pdata$cluster)
        result <- c(names(pdata[i]), "Wilcox", test$statistic, test$p.value, NA, cluster)
      } else if(var.test(pdata[[i]]~pdata$cluster)$p.value < 0.05){
        test <- wilcox.test(pdata[[i]]~pdata$cluster)
        result <- c(names(pdata[i]), "Wilcox", test$statistic, test$p.value, NA, cluster)
      } else {
        test <- t.test(pdata[[i]]~pdata$cluster)
        result <- c(names(pdata[i]), "T-test", test$statistic, test$p.value, NA, cluster)
      }
    } else {
      tab <- table(pdata$cluster, pdata[[i]])
      if(length(levels(pdata[[i]])) == 2){
        test <- chisq.test(tab, correct = TRUE)
        observed <- test$observed; expected <- test$expected
        comparison <- observed > expected
        enriched_class <- apply(comparison, 1, get_unique_labels)[2]
        enriched_class <- paste(enriched_class, collapse = ", ")
        result <- c(names(pdata[i]), "Chi-squared", test$statistic, test$p.value, enriched_class, cluster)
      } else {
        test <- chisq.test(tab, simulate.p.value = TRUE)
        observed <- test$observed; expected <- test$expected
        comparison <- observed > expected
        enriched_class <- apply(comparison, 1, get_unique_labels)[2]
        enriched_class <- paste(enriched_class, collapse = ", ")
        result <- c(names(pdata[i]), "Chi-squared", test$statistic, test$p.value, enriched_class, cluster)
      }
    }
    
    
    results <- rbind(results, result)
  }
}

# Convert to data frame
results <- as.data.frame(results)
colnames(results) <- c("Variable", "Test", "Statistic", "p.value", "Enriched_Class", "cluster")
print(results)
write.csv(results, "../COPD/enrichment_samples_clusters.csv")
results$p.value <- as.numeric(results$p.value)
# Filter rows with p-value < 0.05
filtered_results <- results[results$p.value < 0.05, ]

# Display the filtered results
print(filtered_results)

## ----------------------------------------------------------------------------------------------------------------------------------------
# # Correlations of shap by phenotypic variables
# shap_df <- as.data.frame(t(heatmap_matrix))
# # Ensure rownames of SHAP data are in a column
# shap_df <- shap_df %>% rownames_to_column(var = "sample_id")

# # Merge dataframes
# merged_df <- merge(shap_df, phenotypic_ctics_interest, by = "sample_id")

# # Correlation between SHAP values and a continuous variable
# cor_results <- merged_df %>%
#   select(-sample_id, -dis_condition, -GOLD_stage, -sex, -smoker, -pneumocystis_colonization, -emphysema) %>%
#   summarise(across(everything(), ~cor(., merged_df$emphysema, method = "spearman", use = "complete.obs")))
