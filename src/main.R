#########################################################################################
# SBEACD
#
# Iria Pose Lagoa 2025
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
## Libraries
library(argparse)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read command line arguments 
parser <- ArgumentParser(description = "Run ML pipeline with optional arguments")

# Required arguments
parser$add_argument("--directory_to_load", required = TRUE, help = "Raw data folder for expression data")
parser$add_argument("--file_name", required = TRUE, help = "Name of the expression data .Rda file")
parser$add_argument("--target_var", required = TRUE, help = "Name of the target variable")
parser$add_argument("--procedure", required = TRUE, help = "Feature selection procedure")
parser$add_argument("--ml_models_to_run", nargs = "+", default = c("rf", "knn", "svm_r", "svm_p", "glm", "xgb"),
                   help = "List of ML models to run")

# Optional arguments
parser$add_argument("--directory_to_save", default = NULL, help = "Folder to save results")
parser$add_argument("--threshold_value", type = "numeric", default = NULL, help = "Threshold for top MRMR genes")
parser$add_argument("--disease_code", default = NULL, help = "Disgenet code for the disease")
parser$add_argument("--dea_genes", default = NULL, help = "List of differentially expressed genes")
parser$add_argument("--mrmr_genes", default = NULL, help = "List of mrmr genes")
parser$add_argument("--alternative_genes", default = NULL, help = "List of alternative selected genes")

args <- parser$parse_args()

directory_to_load <- args$directory_to_load
file_name <- args$file_name
target_var <- args$target_var
procedure <- args$procedure
ml_models_to_run_vector <- args$ml_models_to_run
disease_code <- args$disease_code
directory_to_save <- args$directory_to_save
threshold_value <- args$threshold_value
dea_genes <- args$dea_genes
mrmr_genes <- args$mrmr_genes
alternative_genes <- args$alternative_genes

## ----------------------------------------------------------------------------------------------------------------------------------------
# Start
set.seed(1234)
ini <- Sys.time()

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load
# Load libraries
source("feature_selection/library_feature_selection.R")
source("ML_models/library_ml_models.R") # NEEDED for computing the explainer!!!
source("explainability/library_explainability.R")
source("support/library_help.R")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. CREATE TRAIN/TEST SET
## Load expression .Rda file and split the data into training and test set
# 'directory_to_load: raw data folder where you have the expression data
# 'file_name: the naming of you expression data .Rda file
# 'target_var: naming of the target variable
# 'directory_to_save: where you want to save the data
# 'save_objects: TRUE/FALSE if you want to save the splitted train/test data
split_data <- obtain_split_data(directory_to_load = directory_to_load, 
                file_name = file_name, 
                target_var = target_var, 
                directory_to_save = directory_to_save)

expression_train <- split_data$train
expression_test <- split_data$test

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. FEATURE SELECTION
## Call function to perform feature selection
# 'procedure: mrmr, disease_related, disease_related_entire_list, data_driven, omnipath_disease_related, omnipath_data_driven, omnipath_intersection, omnipath_union
# 'expression_data: expression data (samples as rows genes as columns)
# 'target_var: naming of the target variable
# 'threshold_value: numeric thresold for select the top mrmr genes if you already know it
# 'disease_code: disease disgenet code
# 'dea_genes: list of differentially expressed genes
# 'directory_to_load: raw data folder where you have the expression data
# 'directory_to_save: where you want to save the data
# 'save_list: TRUE/FALSE if you want to save the list of genes
genes_list <- run_feature_selection(procedure,
                expression_data = expression_train, 
                target_var = target_var, 
                threshold_value = threshold_value, 
                disease_code = disease_code, 
                dea_genes = dea_genes, 
                mrmr_genes = mrmr_genes,
                alternative_genes = alternative_genes,
                directory_to_load = directory_to_load, 
                directory_to_save = directory_to_save)

## Filter expression data
col_selection <- c(target_var, genes_list)
expression_train <- expression_train[, colnames(expression_train) %in% col_selection]

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. CLASSIFICATION PERFORMANCE

## Call function to run models
# 'data: expression train data
# 'target_var: target variable to perform classification
# 'models_to_run: ml models to run (rf, knn, svm_r, svm_p, glm, xgb)
# 'multi_binary: indicate if the classification is binary or not (binary, FALSE)
# 'directory_to_save_models: where you want to save the results of the models
# 'save_models: TRUE/FALSE if you want to save the models and folds data
results_models <- ML_models(data = expression_train, 
                    target_var = target_var, 
                    models_to_run = ml_models_to_run_vector, 
                    directory_to_save = directory_to_save)

# Extract all models results
# 'models_ro_run: vector with the models you want to run
# 'results_models: list with results for each of the models
# 'expression_train: expression train data
# 'expression_test: expression test data
# 'target_var: name of the target variable
extracted_results <- extract_models_results(
  models_to_run = ml_models_to_run_vector,
  results_models = results_models,
  expression_train = expression_train,
  expression_test = expression_test,
  target_var = target_var
)

# Save list of results
saveRDS(extracted_results, file = paste0(directory_to_save,"results_models.rds"))

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. EXPLAINABILITY
# ' 

shap_results <- calculate_all_shap(ml_models_to_run = ml_models_to_run_vector, 
                                   results_models = results_models, 
                                   expression_data = expression_train, 
                                   target_var = target_var, 
                                   directory_to_save = directory_to_save)


saveRDS(shap_results, file = paste0(directory_to_save,"results_shap.rds"))

fin <- Sys.time()
time <- fin-ini
print("Done!")
print(time)