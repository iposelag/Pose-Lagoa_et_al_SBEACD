#!/bin/Rscript
###############################################################################
######## library w/ functions for classical ml models using tidymodels ########
###############################################################################
## command: 
## output: 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library("DALEXtra")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to generate predictions using the explainer
generate_shap_contributions <- function(explainer, data_to_shap, sample, directory_to_save) {

  data = data_to_shap[sample, -1]
  sample = paste0(data_to_shap[sample, 1],"_",rownames(data_to_shap)[sample])
  shap_values <- predict_parts(explainer = explainer, new_observation = data, 
                                type = "shap", B = 1, random_state =1)
  results <- list(sample = sample,
                  shap_values = shap_values)   
  if(!is.null(directory_to_save)){
    if ("explainability"%in%list.files(paste0(directory_to_save)) == FALSE){
      dir.create(paste0(directory_to_save, "explainability"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"explainability/")
    if (as.character(explainer$label)%in%list.files(paste0(directory)) == FALSE){
      dir.create(paste0(directory, as.character(explainer$label)), showWarnings = FALSE)}
    directory <- paste0(directory,as.character(explainer$label), "/")
    saveRDS(results, file = paste0(directory,sample,".rds"))
  }

  return(results)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate variable importance for a specific batch
calculate_variable_importance_batch <- function(explainer, train_data, samples, directory_to_save) {
  
  ini <- Sys.time()
  pp <- data.frame()
  results <- list()
  for (sample in 1:samples) {
    results_shap <- generate_shap_contributions(explainer, train_data, sample, directory_to_save)
    p <- results_shap$shap_values %>%
      mutate(sample = results_shap$sample)
    pp <- rbind(pp, p)
  }

  df <- as.data.frame(pp)
  df <- df %>% group_by(variable_name, sample) %>%
    summarize(mean_shap = mean((contribution)), .groups = "drop")
  fin <- Sys.time()

  results$shap_values <- df
  results$time <- fin-ini   

  return(results)
}

calculate_all_shap <- function(ml_models_to_run, results_models, expression_data, target_var, directory_to_save = NULL) {
  
  # Initialize list to store results
  all_results_shap <- list()
  
  # Loop through each model in ml_models
  for(classif in ml_models_to_run) {
    # Rename classifier
    classif_name <- rename_classifier(classif)
    
    # Extract the fitted model for the classifier
    fit <- results_models[[classif]]$fit
    
    # Convert target variable to numeric
    target_var_numeric <- as.numeric(expression_data[[target_var]]) - 1
    
    # Create an explainer for SHAP values
    explainer <- explain_tidymodels(fit, 
                                    data = expression_data %>% select(-all_of(target_var)), 
                                    y = target_var_numeric, 
                                    label = classif_name)
    
    # Set batch size to number of rows in training data
    batch <- nrow(expression_data)
    
    # Calculate SHAP values in batch
    results_shap <- calculate_variable_importance_batch(explainer, expression_data, samples = batch, directory_to_save)
    
    # Store the results for the classifier in the list
    all_results_shap[[classif]] <- results_shap
  }
  
  return(all_results_shap)
}
