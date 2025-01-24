#!/bin/Rscript
###############################################################################
######## library w/ functions for processate ML models results metrics ########
###############################################################################
## command: 
## output: 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to process cross-validation metrics
process_cross_validation_metrics <- function(results_models, classifiers) {
  all_metrics <- data.frame()
  
  for (classif in classifiers) {
    # Extract and structure metrics for cross-validation
    metrics <- as.data.frame(results_models[["cross_validation"]][[classif]][["model_metrics"]])
    colnames(metrics) <- c("metric", "mean", "sd")
    metrics$classifier <- classif
    metrics <- metrics[, c("classifier", "metric", "mean", "sd")]
    
    # Add norMCC as a new metric row
    if ("mcc" %in% metrics$metric) {
      mcc_row <- metrics[metrics$metric == "mcc", ]
      norMCC <- (mcc_row$mean + 1) / 2  # Compute norMCC
      metrics <- rbind(metrics, data.frame(
        classifier = classif,
        metric = "normMCC",
        mean = norMCC,
        sd = NA  # No standard deviation for derived metric
      ))
    }
    
    # Combine into the overall metrics data frame
    all_metrics <- rbind(all_metrics, metrics)
  }
  
  # Add percentage column
  all_metrics$"%" <- all_metrics$mean * 100
  
  return(all_metrics)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to process test metrics
process_test_metrics <- function(results_models, classifiers) {
  all_metrics <- data.frame()
  
  for (classif in classifiers) {
    # Extract and structure metrics for test
    metrics <- as.data.frame(results_models[["test"]][[classif]][["model_metrics"]])
    metrics <- metrics[, c(1, 3, 2)]
    colnames(metrics) <- c("metric", "estimate", "estimator")
    metrics$classifier <- classif
    metrics <- metrics[, c("classifier", "metric", "estimate", "estimator")]
    
    # Add norMCC as a new metric row
    if ("mcc" %in% metrics$metric) {
      mcc_row <- metrics[metrics$metric == "mcc", ]
      norMCC <- (mcc_row$estimate + 1) / 2  # Compute norMCC
      metrics <- rbind(metrics, data.frame(
        classifier = classif,
        metric = "normMCC",
        estimate = norMCC,
        estimator = NA  # No estimator for derived metric
      ))
    }
    
    # Combine into the overall metrics data frame
    all_metrics <- rbind(all_metrics, metrics)
  }
  
  # Add percentage column
  all_metrics$"%" <- all_metrics$estimate * 100
  
  return(all_metrics)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Procesate times
process_times <- function(results_models, unit = "mins") {
  # Ensure the unit is valid for as.numeric()
  valid_units <- c("secs", "mins", "hours", "days")
  if (!(unit %in% valid_units)) {
    stop("Invalid unit. Choose from 'secs', 'mins', 'hours', or 'days'.")
  }
  
  # Initialize an empty data frame
  classifier_time <- data.frame(
    classifier = character(), # Classifier names
    time = numeric()          # Times in specified units
  )
  
  # Loop through the times in results_models and populate the data frame
  for (classif in names(results_models$times)) {
    classifier_time <- rbind(classifier_time, 
                            data.frame(
                                classifier = classif,
                                time = as.numeric(results_models$times[[classif]], units = unit)))
  }
  
  # Return the processed data frame
  return(classifier_time)
}

  