#!/bin/Rscript
###############################################################################
######## library w/ functions for perform the feature selection step ########
###############################################################################
## command: 
## output: 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(dplyr); library(OmnipathR); 


# Load the feature selection (list of genes)
load_gene_list <- function(genes_list) {
  file_path <- paste0("../COPD/results_",genes_list,"/feature_selection/",genes_list,".txt")
  genes <- scan(file_path, what = "character", sep = ",")
  return(genes)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# ALTERNATIVE: functions to obtain an alternative selection of genes 
obtain_alternative_genes <- function(alternative_genes, directory_to_load,  directory_to_save = NULL) {
  
  cat("Obtaining alternative list of genes:", alternative_genes, "\n")
  # Load gene-disease associations data
  file_path <- paste0(directory_to_load, paste0(alternative_genes,".txt"))
  genes_list <- scan(file_path, what = "character", sep = ",")
  
  # Save the list of genes if save_list is TRUE
  if (!is.null(directory_to_save)) {
    if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
      dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write(genes_list, file = paste0(directory, paste0(alternative_genes,".txt")), ncolumns = 1)
  }
  cat(paste0(alternative_genes, " genes:\n"), genes_list, "\n")
  return(genes_list)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# DISGENET: functions to obtain disease_related genes using disgenet database

# Download data tables from disgenet
########### PENDIENTE #############

# Generate disgenet gene-disease associations curated database
obtain_disease_related_curated_genes <- function(disease_code, directory_to_load,  directory_to_save = NULL) {
  
  cat("Obtaining disease related curated genes\n")
  # Load gene-disease associations data
  disgenet <- read.csv(paste0(directory_to_load, "disgenet_tables/disgenet_curated_gene_disease_associations.tsv"), sep = "\t")
  disgenet_disease_specific_associations <- filter(disgenet, diseaseId == disease_code)
  disgenet_disease_specific_gene_list <- disgenet_disease_specific_associations$geneSymbol
  
  # Save the list of genes if save_list is TRUE
  if (!is.null(directory_to_save)) {
    if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
      dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write(disgenet_disease_specific_gene_list, file = paste0(directory, "disease_related.txt"), ncolumns = 1)
    # Save list of genes: input for guildify
    write(disgenet_disease_specific_gene_list, file = paste0(directory, "disease_related_guildify_input.txt"), 
          sep = ";", ncolumns = length(disgenet_disease_specific_gene_list))
  }

  cat("Disease Related genes:\n", disgenet_disease_specific_gene_list, "\n")
  return(disgenet_disease_specific_gene_list)
}

# Generate disgenet gene-disease associations entire list
obtain_disease_related_entire_genes <- function(disease_code, directory_to_load, directory_to_save = NULL){

    cat("Obtaining disease related entire genes\n")
    # Load gene-disease associations entire data
    disgenet_entire_associations <- read.csv(paste0(directory_to_load, "disgenet_tables/", disease_code, "_disease_gda_summary.tsv"), sep ="\t")
    disgenet_entire_list <- disgenet_entire_associations$Gene
    # Save the list of genes if save_list is TRUE
    if (!is.null(directory_to_save)) {
      if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"feature_selection/")
      # Save list of genes: input for models
      write(disgenet_entire_list,file=paste0(directory,"disease_related_entire_list.txt"), ncolumns = 1)
    }
    cat("Disease Related entire list of genes:\n", disgenet_entire_list, "\n")
    return(disgenet_entire_list)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# DATA-DRIVEN: functions to obtain disease_related genes using disgenet database

## MRMR
obtain_mrmr <- function(expression_data, target_var, threshold_value = NULL, directory_to_save = NULL){

    cat("Obtaining mrmr genes\n")
    if(is.null(threshold_value)){
      cat("Start of bootstrapping process\n")
      all_mrmr_scores_df <- run_mrmr_bootstrapping(expression_data, target_var, n_iterations = 1000, directory_to_save = NULL)
      threshold_value <- extract_score_threshold(all_mrmr_scores_df, "score")
      cat("Bootstrapping done!\n")
    }else{
      threshold_value <- as.numeric(threshold_value)
    }
    cat("Threshold value:\n", threshold_value, "\n")
    mrmr_results <- run_mrmr(expression_data, target_var, variables = ncol(expression_data) - 1, samples = nrow(expression_data))  # nolint: line_length_linter.
    mrmr_results <- extract_mrmr_features(mrmr_results)
    selected_genes <- mrmr_results[mrmr_results[["Score"]] >= threshold_value, ]
    mrmr_list <- selected_genes$Name
    # Remove leading and trailing whitespace
    mrmr_list <- trimws(mrmr_list)

    # Display the cleaned output
    # print(mrmr_list)

    # Optionally save the genes list
    if (!is.null(directory_to_save)) {
      if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"feature_selection/")
      # Save list of genes: input for models
      write(mrmr_list, paste0(directory,"mrmr.txt"), ncolumns = 1)
    }
    cat("MRMR genes:", mrmr_list, "\n")
    return(mrmr_list)
}

shuffling_sample_disease_categories <- function(data, target_var, directory_to_save = NULL){

  shuffled_names <- sample(data[[target_var]])
  data[[target_var]] <- shuffled_names
  # Optionally save the shuffled dataset
  if (!is.null(directory_to_save)) {
    if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
      dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write.csv(data, paste0(directory,"shuffled_dataset.csv"), row.names = FALSE)
  }

  return(data)
}

# Function to run bootstrapping command
run_mrmr <- function(data, target_var, variables, samples) {

  data[[target_var]] <- as.integer(data[[target_var]])

  # Save the shuffled data to a temporary file
  temp_file <- tempfile(fileext = ".csv")
  write.csv(data, temp_file, row.names = FALSE)

  # Define the command template with placeholders for file, variables, and samples
  command_template <- "./mrmr -i INPUT_FILE -n 500 -t 1 -v VARIABLES -s SAMPLES"
  
  # Create the command by replacing placeholders with actual values
  command <- gsub("INPUT_FILE", temp_file, command_template)
  command <- gsub("VARIABLES", as.character(variables), command)
  command <- gsub("SAMPLES", as.character(samples), command)

  # Run the command and capture the output
  output <- system(command, intern = TRUE)

  return(output)
}

# Function to filter mRMR output and extract mRMR features into a dataframe
extract_mrmr_features <- function(results) {

  # Find the index of the line containing "*** mRMR features ***"
  mrmr_start_index <- which(grepl("\\*\\*\\* mRMR features \\*\\*\\*", results)) + 2
  
  # Extract lines until the end of the output
  feature_lines <- results[mrmr_start_index:length(results)]
  
  # Remove empty lines and any other non-data lines
  feature_lines <- feature_lines[feature_lines != ""]
  # Identify the end of the feature lines based on footer detection
  end_index <- which(grepl("^\\s*\\*{3}", feature_lines))
  feature_lines <- feature_lines[1:(end_index[1] - 1)]
  
  # Convert lines to a data frame by splitting on tab character
  features_df <- read.table(text = feature_lines, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Set column names
  colnames(features_df) <- c("Order", "Feature", "Name", "Score")
  
  return(features_df)
}

run_mrmr_bootstrapping <- function(expression_data, target_var, n_iterations = 1000, directory_to_save = NULL) {
  
  # Initialize an empty data frame to store results across iterations
  all_mrmr_results_df <- data.frame()
  
  for (i in 1:n_iterations) {
    
    # Shuffle the sample disease categories
    shuffled_data <- shuffling_sample_disease_categories(expression_data, target_var, directory_to_save)
    
    # Run bootstrapping on the shuffled data with mRMR
    bootstrapping_results <- run_mrmr(data = shuffled_data, target_var = target_var, 
                                      variables = ncol(shuffled_data) - 1, samples = nrow(shuffled_data))
    
    # Extract mRMR features
    mrmr_results_df <- extract_mrmr_features(bootstrapping_results)
    iteration <- data.frame(iteration = rep(i, nrow(mrmr_results_df)))
    mrmr_results_df <- cbind(mrmr_results_df, iteration)
    all_mrmr_results_df <- rbind(all_mrmr_results_df, mrmr_results_df)

    # Save the results if save_results is TRUE
    if (!is.null(directory_to_save)) {
      if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    output_file <- paste0(directory, "all_mrmr_results.csv")
    write.csv(all_mrmr_results_df, output_file, row.names = FALSE)
    }

    # Append iteration number and scores to results
    score <- data.frame(score = mrmr_results_df$Score)
    iteration_score_df <- cbind(iteration, score)
    
    # Combine results from each iteration
    all_scores_df <- cbind(all_mrmr_results_df, iteration_score_df)
  }
  
  return(all_scores_df)
}

# Function to establish a threshold and select top genes
# The function extract the genes that has the top 1% score
extract_score_threshold <- function(data, score_col, threshold_fraction = 0.01) {
  
  # Check if the score column is normally distributed using the Kolmogorov-Smirnov test
  ks_test_result <- ks.test(data[[score_col]], "pnorm", mean = mean(data[[score_col]]), sd = sd(data[[score_col]]))
  
  # If the p-value is low, we assume the distribution is not normal
  if (ks_test_result$p.value < 0.05) {

    # Sort the scores in decreasing order
    sorted_data <- sort(data[[score_col]], decreasing = TRUE)
    
    # Determine the number of genes to select based on the specified fraction
    n <- length(sorted_data)
    index <- ceiling(threshold_fraction * n)  # Get the index for the threshold
    
    # Establish the threshold value
    threshold_value <- sorted_data[index]  # This is the score that acts as the threshold

    # Return the selected genes and the threshold value
    return(threshold_value)
  } else {

    # Calculate the mean and standard deviation of the scores
    mean_score <- mean(data[[score_col]])
    sd_score <- sd(data[[score_col]])

    z_threshold <- qnorm(1 - threshold_fraction) # the value of 2.33 specifically relates to the upper tail (99th percentile) of a normal distribution when selecting the top 1%

    # Calculate the threshold based on Z-score
    threshold_value <- mean_score + z_threshold * sd_score

  return(threshold_value)
  }

}

## DEA
##### PENDIENTE

## DATA-DRIVEN
obtain_data_driven_genes <- function(dea_genes = NULL, mrmr_genes = NULL, expression_data, target_var, threshold_value = NULL, directory_to_load, directory_to_save = NULL){

    cat("Obtaining data driven genes\n")
    # Obtain mrmr_genes
    if(is.null(mrmr_genes)){
      mrmr_genes <- obtain_mrmr(expression_data, target_var, threshold_value, directory_to_save)
    }
    # Obtain dea genes
    if(is.null(dea_genes)){
      file_path <- paste0(directory_to_load, "dea.txt")
      dea_genes <- scan(file_path, what = "character", sep = ",")
      cat("DEA genes:\n", dea_genes, "\n")
    }
    # data-driven
    data_driven <- union(dea_genes, mrmr_genes)
    if(!is.null(directory_to_save)){
      if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
       dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write(data_driven, file=paste0(directory, "data_driven.txt"), ncolumns = 1)
    }
    cat("Data Driven genes:\n", data_driven, "\n")
    return(data_driven)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# EXPANSIONS: functions to expand the seed lists

# Generate Omnipath expansions
obtain_omnipath_expansion <- function(genes, directory_to_save = NULL){

    cat("Obtaining omnipath expansion genes\n")
    # Get the name of the variable as a string
    variable_name <- deparse(substitute(genes))

    # Run omnipath: first interaction genes
    omnipath <- import_all_interactions(
        organism = 9606,
        directed  = 'no'
        ) %>%
        filter(source_genesymbol %in% genes |
           target_genesymbol %in% genes ) %>%
        filter(curation_effort > 1)
    # Target and souce genes
    omnipath_target_genes <- omnipath$target_genesymbol
    omnipath_source_genes <- omnipath$source_genesymbol
    # Union of seed genes and expansion
    expansion_omnipath <- union(omnipath_target_genes, omnipath_source_genes)
    # Add seed genes that are not in the expansion
    expansion_omnipath <- union(expansion_omnipath, genes)
    # Save list of expansion genes: input for models
    if (!is.null(directory_to_save)) {
      if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write(expansion_omnipath,file= paste0(directory, "omnipath_",variable_name ,".txt"), 
          ncolumns = 1)
    }
    cat("Omnipath expansion genes:\n", expansion_omnipath, "\n")
    return(expansion_omnipath)
}

# Generate intersection or union of Omnipath expansions
obtain_omnipath_combined <- function(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save = NULL, operation = "intersection") {
  
  # Obtain the disease-related and data-driven gene lists
  disease_related <- obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save)
  data_driven <- obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save)
  
  # Expand using OmniPath
  expansion_disease_related <- obtain_omnipath_expansion(disease_related, directory_to_save)
  expansion_data_driven <- obtain_omnipath_expansion(data_driven, directory_to_save)
  
  # Perform intersection or union based on the operation parameter
  if (operation == "intersection") {
    cat("Obtaining omnipath intersection genes\n")
    genes_list <- intersect(expansion_data_driven, expansion_disease_related)
    file_name <- "omnipath_intersection.txt"
  } else if (operation == "union") {
    cat("Obtaining omnipath union genes\n")
    genes_list <- union(expansion_data_driven, expansion_disease_related)
    file_name <- "omnipath_union.txt"
  } else {
    stop("Invalid operation. Please choose either 'intersection' or 'union'.")
  }
  
  # Save the results if save_list is TRUE
  if (!is.null(directory_to_save)) {
    if ("feature_selection"%in%list.files(paste0(directory_to_save)) == FALSE){
      dir.create(paste0(directory_to_save, "feature_selection"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"feature_selection/")
    # Save list of genes: input for models
    write(genes_list, file = paste0(directory, file_name), ncolumns = 1)
  }
  cat(genes_list,"\n")
  return(genes_list)
}


# Generate GUILDify expansions
##### PENDIENTE


## ----------------------------------------------------------------------------------------------------------------------------------------
# Define the main function for feature selection
run_feature_selection <- function(procedure, expression_data, target_var, threshold_value =NULL, disease_code = NULL, dea_genes = NULL, mrmr_genes = NULL, alternative_genes = NULL, directory_to_load = NULL, directory_to_save = NULL) {
  
  # Switch case for selecting feature selection procedure
  switch(procedure,
         
         "mrmr" = {
           genes_list <- obtain_mrmr(expression_data, target_var, threshold_value, directory_to_save)
         },

         "disease_related" = {
           genes_list <- obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save)
         },
         
         "disease_related_entire_list" = {
           genes_list <- obtain_disease_related_entire_genes(disease_code, directory_to_load, directory_to_save)
         },
         
         "data_driven" = {
           genes_list <- obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save)

         },
         
         "omnipath_disease_related" = {
           disease_related <- obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save)
           genes_list <- obtain_omnipath_expansion(disease_related, directory_to_save)
         },
         
         "omnipath_data_driven" = {
           data_driven <- obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save)
           genes_list <- obtain_omnipath_expansion(data_driven,  directory_to_save)
         },
         
         "omnipath_intersection" = {
           genes_list <- obtain_omnipath_combined(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save, operation = "intersection")
         },
         
         "omnipath_union" = {
           genes_list <- obtain_omnipath_combined(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save, operation = "union")
         },

         "alternative" = {
           genes_list <- obtain_alternative_genes(alternative_genes, directory_to_load, directory_to_save)
         },
         
        #  "guildify_disease_related" = {
        #    # Example function for guildify expansion
        #    # Replace with the actual function or code for Guildify expansion if available
        #    all_results <- obtain_guildify_expansion(disease_code, directory_to_load)
        #  },
                  
        #  "guildify_data_driven" = {
        #    # Example function for guildify expansion
        #    # Replace with the actual function or code for Guildify expansion if available
        #    all_results <- obtain_guildify_expansion(disease_code, directory_to_load)
        #  },
         
         stop("Invalid procedure specified.")
  )
  
  return(genes_list)
}