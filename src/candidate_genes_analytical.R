#########################################################################################
# Selection of candidate genes analytical method
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(GGally)
library(reshape2)
library(ggrepel)
library(ggExtra)
library(enrichR)
library(gridExtra)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Input parameters
input_data = "normalized_shaps"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load the data
if (input_data == "rank"){
  data <- read.csv("SHAP/results/rank_results.csv")
}else if(input_data == "normalized_shaps"){
  data <- read.csv("../COPD/results/norm_results.csv")
}

# write.csv(data_filtered, "../COPD/norm_candidate_results.csv", row.names = FALSE)
# data <- data %>% filter(aggregated_max>0.0107)


## ----------------------------------------------------------------------------------------------------------------------------------------
# Load gene lists
# Disease related
disease_related = read.csv("../COPD/raw_data/disease_related.txt", header = FALSE)
disease_related_all = read.csv("../COPD/raw_data/disease_related_entire_list.txt", header = FALSE)
disease_related_guildify = read.csv("../COPD/raw_data/guildify_disease_related.txt", header = FALSE) # Las expansiones están hechas sobre la curada. 
disease_related_omnipath = read.csv("../COPD/raw_data/omnipath_disease_related.txt", header = FALSE)  # Las expansiones están hechas sobre la curada. 

# Data driven
minmax = read.csv("../COPD/raw_data/mrmr.txt", header = FALSE)
dea= read.csv("../COPD/raw_data/dea.txt", header = FALSE)
data_driven = read.csv("../COPD/raw_data/data_driven.txt", header = FALSE)
data_driven_guildify = read.csv("../COPD/raw_data/guildify_data_driven.txt", header = FALSE)
data_driven_omnipath = read.csv("../COPD/raw_data/omnipath_data_driven.txt", header = FALSE)

# Rename the column to a consistent name (e.g., "Gene")
colnames(disease_related) <- colnames(disease_related_guildify) <- colnames(disease_related_omnipath) <- "Gene"
colnames(disease_related_all) <- colnames(minmax) <- colnames(dea) <- colnames(data_driven) <- "Gene"
colnames(data_driven_guildify) <- colnames(data_driven_omnipath) <- "Gene"

# Create a union of all genes
all_genes <- unique(c(
  disease_related$Gene,
  disease_related_all$Gene,
  disease_related_guildify$Gene,
  disease_related_omnipath$Gene,
  minmax$Gene,
  dea$Gene,
  data_driven$Gene,
  data_driven_guildify$Gene,
  data_driven_omnipath$Gene
))
length(all_genes) # 3114

# Create the combined data frame
feature_lists <- data.frame(
  Gene = all_genes,
  Disease_Related = all_genes %in% disease_related$Gene,
  Disease_Related_All = all_genes %in% disease_related_all$Gene,
  Disease_Related_Guildify = all_genes %in% disease_related_guildify$Gene,
  Disease_Related_Omnipath = all_genes %in% disease_related_omnipath$Gene,
  Minmax = all_genes %in% minmax$Gene,
  DEA = all_genes %in% dea$Gene,
  Data_Driven = all_genes %in% data_driven$Gene,
  Data_Driven_Guildify = all_genes %in% data_driven_guildify$Gene,
  Data_Driven_Omnipath = all_genes %in% data_driven_omnipath$Gene
)
 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Analyze relation between aggregated metrics
 
# Pairwise scatter plot matrix to explore relationships between variables
pdf("../COPD/plots/scatter_with_outlier_aggregated_values.pdf", width = 12, height = 8)
ggpairs(data[, c("aggregated_max", "aggregated_mean", "aggregated_median", "sd", "occurrence_count")],
        title = "Pairwise Scatter Plots of SHAP Metrics")
dev.off()

# Outlier removal, gene FANCC with max = 0.18, shown to be implicated in the diseases (part of the 15-genes gene signature) (https://pmc.ncbi.nlm.nih.gov/articles/PMC6422831/)
data_outlier = data[data$aggregated_max < 0.1, ]

pdf("../COPD/plots/scatter_without_outlier_aggregated_values.pdf", width = 12, height = 8)
# Pairwise scatter plot matrix to explore relationships between variables
ggpairs(data_outlier[, c("aggregated_max", "aggregated_mean", "aggregated_median", "sd", "occurrence_count")],
        title = "Pairwise Scatter Plots of SHAP Metrics (withou outlier)")
dev.off()

# Correlation between aggregated_max and occurrence_count
spearman_corr <- cor(data_outlier$aggregated_max, data_outlier$occurrence_count, method = "spearman")
pearson_corr <- cor(data_outlier$aggregated_max, data_outlier$occurrence_count, method = "pearson")

spearman_corr <- cor.test(data_outlier$aggregated_max, data_outlier$occurrence_count, method = "spearman")
pearson_corr <- cor.test(data_outlier$aggregated_max, data_outlier$occurrence_count, method = "pearson")

# Reshape data to long format for easier plotting
data_long <- melt(data_outlier, id.vars = "variable_name")

# Histogram for each variable
p1 <- ggplot(data_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Histograms of SHAP Metrics and Occurrence Count",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal()

# Density plot for each variable
p2 <- ggplot(data_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Density Plots of SHAP Metrics and Occurrence Count",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display plots
print(p1)
print(p2)

p3 <- ggplot(data, aes(x = occurrence_count, y = aggregated_max, color = aggregated_mean)) +
  geom_point(alpha = 0.7) +
  # scale_color_gradient(low = "green", high = "purple", name = "Mean") +
  labs(
    title = "Max SHAP vs. Occurrence Count",
    x = "Occurrence Count",
    y =  "Max SHAP Value",
    color = "Mean"
  ) +
  theme_minimal()
p3

# Filter genes based on criteria for annotation
annotated_genes <- data_outlier[data_outlier$occurrence_count < 10 & data_outlier$aggregated_max > 0.0109, ]
p4 <- ggplot(data_outlier, aes(x = occurrence_count, y = aggregated_max, color = aggregated_mean)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = annotated_genes, 
    aes(label = variable_name),
    size = 3,
    color = "black",
    max.overlaps = 20
  ) +
  geom_hline(yintercept = 0.01, color = "darkgrey")+ #  linetype = "dashed", 
  labs(
    title = "Max SHAP vs. Occurrence Count (Annotated Genes)",
    x = "Occurrence Count",
    y = "Max SHAP Value",
    color = "Mean"
  ) +
  theme_minimal()
p4

# Add marginal density plots
p4_with_density <- ggMarginal(
  p4,
  type = "density",
  margins = "both",    # Add density plots to both x and y axes
  fill = "lightgrey",  # Fill color for density plots
  color = "grey"      # Border color for density plots
)

# Display the plot
print(p4_with_density)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Stablish a threshold for selecting candidate genes

# --------------------Analyse percentils
nrow(data[data$aggregated_max > 0.01, ]) # 193
nrow(data) # 2426
(nrow(data[data$aggregated_max > 0.01, ]) / nrow(data))*100 # 7% of the genes

percentile_10 <- quantile(data$aggregated_max, probs = 0.10); percentile_10 # 0.0004262628
percentile_90 <- quantile(data$aggregated_max, probs = 0.90); percentile_90 # 0.007623987
percentile_95 <- quantile(data$aggregated_max, probs = 0.95); percentile_95 # 0.01618515
nrow(data[data$aggregated_max > percentile_90, ]) # 243
243 - 192 # 51 genes

# Density plot with percentiles
p5 <- ggplot(data, aes(x = aggregated_max)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = percentile_90) +
  geom_vline(xintercept = percentile_95) +
  geom_vline(xintercept = 0.0106) +
  # facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Density Plots of aggregated SHAP max and percentiles",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p5

# --------------------Knee/Elbow Method to find a cut-off
# Compute density
density_data <- density(data$aggregated_max)

# Extract x (values) and y (density) coordinates
x_vals <- density_data$x
y_vals <- density_data$y

# Compute the first derivative (slope)
slope <- diff(y_vals) / diff(x_vals)

# Compute the second derivative (curvature)
curvature <- diff(slope) / diff(x_vals[-length(x_vals)])

# Find the index where the curvature is closest to zero
flatten_index <- which.min(abs(curvature)) # 192

# Get the x-value (threshold) at the start of the tail
tail_start <- x_vals[flatten_index] # 0.06868414

# Create data frames for plotting
density_df <- data.frame(x = x_vals, y = y_vals)
slope_df <- data.frame(x = x_vals[-length(x_vals)], slope = slope)
curvature_df <- data.frame(x = x_vals[-c(length(x_vals), length(x_vals) - 1)], curvature = curvature)

tail_start = 0.0106

# Plot density
density_plot <- ggplot(density_df, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_vline(xintercept = tail_start, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.0109, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Density Plot",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal()

# Plot first derivative (slope)
slope_plot <- ggplot(slope_df, aes(x = x, y = slope)) +
  geom_line(color = "green") +
  geom_vline(xintercept = tail_start, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.0109, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "First Derivative (Slope)",
    x = "Value",
    y = "Slope"
  ) +
  theme_minimal()

# Plot second derivative (curvature)
curvature_plot <- ggplot(curvature_df, aes(x = x, y = curvature)) +
  geom_line(color = "purple") +
  geom_vline(xintercept = tail_start, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.0109, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Second Derivative (Curvature)",
    x = "Value",
    y = "Curvature"
  ) +
  theme_minimal()

# Print all plots one below the other
grid.arrange(density_plot, slope_plot, curvature_plot, ncol = 1)

# --------------------Window sizes
# Compute density
density_data <- density(data$aggregated_max)
x_vals <- density_data$x  # X-axis values
y_vals <- density_data$y  # Density values

# Compute first derivative (slope)
slope <- diff(y_vals) / diff(x_vals)

# Compute second derivative (slope)
curvature <- diff(slope) / diff(x_vals[-length(x_vals)])

# Define the window size (continuous)
window_size <- 0.001  # Adjust this based on your needs

# Initialize an empty data frame for results
results <- data.frame(
  x_start = numeric(),
  x_end = numeric(),
  y_diff = numeric(),
  slope_diff = numeric()
)

# Compute differences in y for each continuous window
start_index <- 1
while (start_index <= length(x_vals)) {
  # Find the indices within the current window
  window_indices <- which(x_vals >= x_vals[start_index] & x_vals < (x_vals[start_index] + window_size))
  
  # Compute the y difference for the window
  if (length(window_indices) > 1) {
    y_diff <- max(y_vals[window_indices]) - min(y_vals[window_indices])
    slope_diff <- max(slope[window_indices]) - min(slope[window_indices])
  } else {
    y_diff <- NA  # Not enough points to compute difference
    slope_diff <- NA  # Not enough points to compute difference
  }
  
  # Store the results
  results <- rbind(
    results,
    data.frame(
      x_start = x_vals[start_index],
      x_end = x_vals[start_index] + window_size,
      y_diff = y_diff,
      slope_diff = slope_diff
    )
  )
  
  # Move to the next window
  start_index <- max(window_indices) + 1
}

# Display the results
head(results)

# Plot the original density curve
density_plot <- ggplot(data = data.frame(x = x_vals, y = y_vals), aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_vline(xintercept = 0.0109, color = "black", linetype = "dashed", size = 0.2) +
  labs(
    title = "Density Plot with Window-Based Y Differences",
    x = "X Value",
    y = "Density"
  ) +
  theme_minimal()

# Plot the y_diff values
diff_plot <- ggplot(results, aes(x = x_start, y = y_diff)) +
  geom_line(color = "red") +
  geom_point(color = "darkred") +
  geom_vline(xintercept = 0.0109, color = "black", linetype = "dashed", size = 0.2) +
  labs(
    title = "Window-Based Y Differences",
    x = "X Start of Window",
    y = "Y Difference"
  ) +
  theme_minimal()

# Plot the slope_diff values
slope_diff_plot <- ggplot(results, aes(x = x_start, y = slope_diff)) +
  geom_line(color = "red") +
  geom_point(color = "darkred") +
  geom_vline(xintercept = 0.0109, color = "black", linetype = "dashed", size = 0.2) +
  labs(
    title = "Window-Based Slope Differences",
    x = "X Start of Window",
    y = "Slope Difference"
  ) +
  theme_minimal()

# Combine the plots
pdf("../COPD/plots/density_aggregated_SHAP.pdf", width = 14 , height =8)
grid.arrange(density_plot, slope_diff_plot, ncol = 1)
dev.off()

# --------------------Selection of the threshold
threshold = 0.0109

sgenes = data[data$aggregated_max > threshold, ]
nrow(sgenes) # 172 genes

(nrow(sgenes) / nrow(data))*100 # 7.08986

# Sort the data frame in descending order of aggregated_max
sgenes_sorted <- sgenes[order(-sgenes$aggregated_max), ]
head(sgenes_sorted)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Enrichment of the selected genes

databases <- enrichR::listEnrichrDbs()
print(head(databases)) # View available databases

# Perform enrichment analysis using selected databases
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "WikiPathways_2023_Human", "GWAS_Catalog_2023",
         "ClinVar_2019", "COVID-19_Related_Gene_Sets_2021", "DisGeNET",
         "Drug_Perturbations_from_GEO_2014", "DrugMatrix", "Jensen_DISEASES", "Reactome_Pathways_2024")
enrichment_results <- enrichR::enrichr(sgenes_sorted$variable_name, dbs)

# Extract results for GO Biological Process
go_results <- enrichment_results[["GO_Biological_Process_2023"]]
kegg_results <- enrichment_results[["KEGG_2021_Human"]]
wiki_results <- enrichment_results[["WikiPathways_2023_Human"]]
gwas_results <- enrichment_results[["GWAS_Catalog_2023"]]

clinvar_results <- enrichment_results[["ClinVar_2019"]]
covid_results <- enrichment_results[["COVID-19_Related_Gene_Sets_2021"]]
disgenet_results <- enrichment_results[["DisGeNET"]]
drugs_perturbations_results <- enrichment_results[["Drug_Perturbations_from_GEO_2014"]]
drugmatrix_results <- enrichment_results[["DrugMatrix"]]
jensen_diseases_results <- enrichment_results[["Jensen_DISEASES"]]
reactome_results <- enrichment_results[["Reactome_Pathways_2024"]]

sgo_results = go_results[go_results$Adjusted.P.value <= 0.05, ]; nrow(sgo_results)
skegg_results = kegg_results[kegg_results$Adjusted.P.value <= 0.05, ]; nrow(skegg_results)
swiki_results = wiki_results[wiki_results$Adjusted.P.value <= 0.05, ]; nrow(swiki_results)
sgwas_results = gwas_results[gwas_results$Adjusted.P.value <= 0.05, ]; nrow(sgwas_results)

sclinvar_results = clinvar_results[clinvar_results$Adjusted.P.value <= 0.05, ]; nrow(sclinvar_results)
scovid_results = covid_results[covid_results$Adjusted.P.value <= 0.05, ]; nrow(scovid_results)
sdisgenet_results = disgenet_results[disgenet_results$Adjusted.P.value <= 0.05, ]; nrow(sdisgenet_results)
sdrugs_perturbations_results = drugs_perturbations_results[drugs_perturbations_results$Adjusted.P.value <= 0.05, ]; nrow(sdrugs_perturbations_results)
sdrugmatrix_results = drugmatrix_results[drugmatrix_results$Adjusted.P.value <= 0.05, ]; nrow(sdrugmatrix_results)
sjensen_diseases_results = jensen_diseases_results[jensen_diseases_results$Adjusted.P.value <= 0.05, ]; nrow(sjensen_diseases_results)
sreactome_results = reactome_results[reactome_results$Adjusted.P.value <= 0.05, ]; nrow(sreactome_results)

# Plot the top 10 enriched pathways
ggplot(go_results[1:20, ], aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched GO Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

# Plot the top 10 enriched pathways
ggplot(kegg_results[1:20, ], aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched GO Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

# Plot the top 10 enriched pathways
ggplot(wiki_results[1:20, ], aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched GO Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

# Plot the top 10 enriched pathways
ggplot(gwas_results[1:20, ], aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched GO Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

sgo_results = go_results[go_results$Adjusted.P.value <= 0.05, ]
ggplot(go_results[1:30, ], aes(x = reorder(Term, Odds.Ratio), y = Odds.Ratio, fill=Adjusted.P.value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched GO Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

sgo_results = go_results[go_results$Adjusted.P.value <= 0.05, ]
ggplot(sreactome_results[1:30, ], aes(x = reorder(Term, Odds.Ratio), y = Odds.Ratio, fill=Adjusted.P.value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched Reactome Biological Processes",
    x = "Biological Process",
    y = "Combined Score"
  ) +
  theme_minimal()

## ----------------------------------------------------------------------------------------------------------------------------------------
# Upset plot origin of selected genes
feature_lists_sgens <- feature_lists %>% filter(Gene %in% sgenes$variable_name)

# Save the original list where each candidate gene came from
candidate_genes <- sgenes[,c("variable_name", "aggregated_max")]
candidate_genes <- candidate_genes %>% rename(Gene = variable_name)
feature_lists_sgens <- merge(feature_lists_sgens,candidate_genes, by = "Gene")
# Order by aggregated max
feature_lists_sgens <- feature_lists_sgens[order(-feature_lists_sgens$aggregated_max), ] 
write.csv(feature_lists_sgens, "../COPD/results/candidate_genes_vs_feature_lists.csv", row.names = FALSE)

# Check the structure of the data
str(feature_lists_sgens)

# Convert logical values (TRUE/FALSE) to binary (1/0)
upset_data <- feature_lists_sgens[, -1]  # Exclude the 'Gene' column
upset_data <- as.data.frame(lapply(upset_data, as.integer))  # Convert to binary (1/0)

# Verify the converted data
print(upset_data)

# Generate the UpSet plot
library(UpSetR)
pdf("../COPD/plots/upset_selected_genes.pdf",  width = 12, height = 8)
upset(
  data = upset_data,
  sets = colnames(upset_data),  # Use all remaining columns as sets
  nintersects = NA,
  order.by = "freq",
  keep.order = TRUE,
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)
)
dev.off()
