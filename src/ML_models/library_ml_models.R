#!/bin/Rscript
###############################################################################
######## library w/ functions for classical ml models using tidymodels ########
###############################################################################
## command: 
## output: 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(tidyverse)
library(tidymodels)
library(Biobase)
library(vip)
library(sva)
library(themis)
library(glmnet)
library(stacks)
library(bundle)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to split data into training and test sets
obtain_split_data <- function(directory_to_load, file_name, target_var, directory_to_save = NULL) {
  
  # Load the dataset
  set.seed(1234)
  expression <- get(load(paste0(directory_to_load, file_name, ".Rda")))
  
  expression[[target_var]] <- as.factor(expression[[target_var]])

  # Split data into train and test sets with stratification on target_var
  expression_split <- initial_split(expression, strata = target_var)
  
  # Training set
  expression_train <- training(expression_split)
  
  # Optionally save the training set
  if (!is.null(directory_to_save)) {
    if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"ML_models")
    # Save data
    save(expression_train, file = paste0(directory, "/expression_train.Rda"))
  }
  
  # Print class imbalance in the training set
  print("Class imbalance: train set")
  print(round(prop.table(table(expression_train[[target_var]])) * 100, 3))
  
  # Testing set
  expression_test <- testing(expression_split)
  
  # Optionally save the testing set
  if (!is.null(directory_to_save)) {
    if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
    directory <- paste0(directory_to_save,"ML_models/")
    # Save data
    save(expression_test, file = paste0(directory, "/expression_test.Rda"))
  }
  
  # Print class imbalance in the test set
  print("Class imbalance: test set")
  print(round(prop.table(table(expression_test[[target_var]])) * 100, 3))
  
  return(list(train = expression_train, test = expression_test))
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# RECIPES: function to preprocess data: 

# Function to create preprocessing recipe with flexible target variable and processing steps
create_recipe <- function(data, target_var, transformation = NULL, normalize = FALSE) {
  # Define base recipe
  recipe_obj <- recipe(as.formula(paste(target_var, "~ .")), data = data)
  
  # Apply BoxCox transformation if specified
  if (!is.null(transformation) && transformation == "boxcox") {
    recipe_obj <- recipe_obj %>% step_BoxCox(all_numeric_predictors())
  }
  
  # Remove highly correlated features
  recipe_obj <- recipe_obj %>% step_corr(all_numeric_predictors(), threshold = 0.85)
  
  # Apply downsampling for class imbalance
  recipe_obj <- recipe_obj %>% step_downsample(all_outcomes())
  
  # Normalize predictors if specified
  if (normalize) {
    recipe_obj <- recipe_obj %>% step_normalize(all_predictors())
  }
  
  return(recipe_obj)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Functions to define models

## Random Forest
rf_model <- function(data, target_var, folds, predictors, metrics, control, directory_to_save = NULL) {

    # Set seed for reproducibility
    set.seed(1234)

    ini_rf <- Sys.time()
    # Define the model
    rf_spec <- rand_forest(
        mtry = tune(), # num of predictor vbles to sample in each split of the tree
        trees = 1000,  # num of trees
        min_n = tune() # minimum node size
    ) %>%
    set_mode("classification") %>%
    set_engine("ranger", oob.error = TRUE, importance = "impurity", seed=1234)

    # Workflow
    rf_wf <- workflow() %>%
     add_recipe(create_recipe(data, target_var)) %>%
     add_model(rf_spec)

    # Collect model parameters
    parameters_rf <- parameters(rf_spec)
    #' It is necessary to finalize the range of values of mtry() that depend on the
    #' number of columns
    parameters_rf <- rf_spec %>% parameters() %>% finalize(predictors)

    # Tuning
    rf_res <- tune_bayes(
        rf_wf,  # add workflow
        resamples = folds,
        param_info = parameters_rf,
        initial = 7, # This could also be a grid search
        iter = 25,
        metrics = metrics,
        control = control
    )

    # Select best selection of parameters
    best_mcc_rf <- select_best(rf_res, metric = "mcc") # Best based on mcc

    # Finalize models
    final_rf <- finalize_model(
        rf_spec,
        best_mcc_rf
    )

    # Update workflow
    rf_wf_updated <- rf_wf %>%
        update_model(final_rf)

    # Fit the model
    rf_fit <- rf_wf_updated %>%
        fit(data = data)

    # Fit resamples
    rf_rs <- rf_fit %>%
        fit_resamples(
            resamples = folds,
            metrics = metrics,
            control = control 
        )

    # Save the final rf model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(rf_wf, rf_res, final_rf, rf_fit, rf_rs, file = paste0(directory,"/rf_model.Rda"))
    }   

    fin_rf <- Sys.time()

    model <- list(
        time = fin_rf-ini_rf,
        resamples = rf_rs,
        fit = rf_fit
    )

    return(model)
}

## Support Vector Machines: radial
svm_r_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

    # Set seed for reproducibility
    set.seed(1234)

    ini_svm_r <- Sys.time()

    # Define the models 
    svm_r_spec <- svm_rbf(
        cost = tune(),
        rbf_sigma = tune()
    ) %>%
        set_engine("kernlab") %>%
        set_mode("classification")

    # Workflow
    svm_r_wf <- workflow() %>%
        add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
        add_model(svm_r_spec)

    # Collect model parameters
    parameters_r_svm <- parameters(svm_r_spec)

    # Tuning
    svm_r_res <- tune_bayes(
        svm_r_wf,  # add workflow
        resamples = folds,
        param_info = parameters_r_svm,
        initial = 7, # This could also be a grid search
        iter = 25,
        metrics = metrics,
        control = control
    )

    # Select best selection of parameters
    best_mcc_svm_r <- select_best(svm_r_res, metric = "mcc") # Best based on mcc

    # Finalize models
    final_svm_r <- finalize_model(
        svm_r_spec,
        best_mcc_svm_r
    )

    # Update workflow
    svm_r_wf_updated <- svm_r_wf %>%
        update_model(final_svm_r)

    # Fit the model
    svm_r_fit <- svm_r_wf_updated %>%
        fit(data = data)

    # Fit resamples
    svm_r_rs <- svm_r_fit %>%
        fit_resamples(
        resamples = folds,
        metrics = metrics,
        control = control # to obatin the prediction column
        )

    # Save the final svm_r model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(svm_r_fit, svm_r_rs, file = paste(directory,"/svm_r_model.Rda", sep=""))
    }
    
    fin_svm_r <- Sys.time()

    model <- list(
        time = fin_svm_r-ini_svm_r,
        resamples = svm_r_rs,
        fit = svm_r_fit
    )
    
    return(model)
}

## Support Vector Machines: polynomial
svm_p_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL){

    # Set seed for reproducibility
    set.seed(1234)

    ini_svm_p <- Sys.time()

    # Define the model
    svm_p_spec <-
        svm_poly(cost = tune(), degree = tune()) %>%
        set_engine("kernlab") %>%
        set_mode("classification")

    # Workflow    
    svm_p_wf <- workflow() %>%
        add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
        add_model(svm_p_spec)

    # Collect model parameters
    parameters_p_svm <- parameters(svm_p_spec)

    # Tuning
    svm_p_res <- tune_bayes(
        svm_p_wf,  # add workflow
        resamples = folds,
        param_info = parameters_p_svm,
        initial = 7, # This could also be a grid search
        iter = 25,
        metrics = metrics,
        control = control
    )

    # Select best selection of parameters
    best_mcc_svm_p <- select_best(svm_p_res, metric = "mcc") # Best based on mcc

    # Finalize models
    final_svm_p <- finalize_model(
        svm_p_spec,
        best_mcc_svm_p
    )

    # Update workflow
    svm_p_wf_updated <- svm_p_wf %>%
        update_model(final_svm_p)

    # Fit the model
    svm_p_fit <- svm_p_wf_updated %>%
        fit(data = data)

    # Fit resamples
    svm_p_rs <- svm_p_fit %>%
        fit_resamples(
        resamples = folds,
        metrics = metrics,
        control = control # to obatin the prediction column
    )

    fin_svm_p <- Sys.time()

    # Save the final svm_p model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(svm_p_fit, svm_p_rs, file = paste(directory,"/svm_p_model.Rda", sep=""))
    }
    
    model <- list(
        time = fin_svm_p-ini_svm_p,
        resamples = svm_p_rs,
        fit = svm_p_fit
    )

    return(model)
}

## Penalized Regression Model for binary classification
glm_binary_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

    # Set seed for reproducibility
    print("binary classification")
    set.seed(1234)

    ini_glm <- Sys.time()

    # Define the model
    glm_spec <- logistic_reg(penalty = tune(),
                          mixture = tune()) %>%
        set_engine("glmnet") %>%
        set_mode("classification")

    # Workflow
    glm_wf <- workflow() %>%
        add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
        add_model(glm_spec)

    # Collect model parameters
    parameters_glm <- parameters(glm_spec)

    # Tuninggggg
    glm_res <- tune_bayes(
        glm_wf,  # add workflow
        resamples = folds,
        param_info = parameters_glm,
        initial = 7, # This could also be a grid search
        iter = 25,
        metrics = metrics,
        control = control
    )

    # Select best selection of parameters
    best_mcc_glm <- select_best(glm_res, metric = "mcc") # Best based on mcc

    # Finalize models 
    final_glm <- finalize_model(
        glm_spec,
        best_mcc_glm
    )

    # Update workflow
    glm_wf_updated <- glm_wf %>%
        update_model(final_glm)

    # Fit the model
    glm_fit <- glm_wf_updated %>%
        fit(data = data)

    # Fit resamples
    glm_rs <- glm_fit %>%
        fit_resamples(
        resamples = folds,
        metrics = metrics,
        control = control # to obatin the prediction column
    )

    fin_glm <- Sys.time()

    # Save the final glm model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(glm_fit, glm_rs, file = paste(directory,"/glm_model.Rda", sep=""))
    }
    
    model <- list(
        time = fin_glm-ini_glm,
        resamples = glm_rs,
        fit = glm_fit
    )

    return(model)
}

## Penalized Regression Model for multiclassification
glm_multiclass_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {
  
  # Set seed for reproducibility
  set.seed(1234)
  
  ini_glm <- Sys.time()
  
  # Define the model
  glm_spec <- multinom_reg(penalty = tune(),
                               mixture = tune()) %>%
    set_engine("glmnet") %>%
    set_mode("classification")
  
  # Workflow
  glm_wf <- workflow() %>%
    add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    add_model(glm_spec)
  
  # Collect model parameters
  parameters_glm <- parameters(glm_spec)
  
  # Tuninggggg
  glm_res <- tune_bayes(
    glm_wf,  # add workflow
    resamples = folds,
    param_info = parameters_glm,
    initial = 7, # This could also be a grid search
    iter = 25,
    metrics = metrics,
    control = control
  )
  
  # Select best selection of parameters
  best_mcc_glm <- select_best(glm_res, metric = "mcc") # Best based on mcc
  
  # Finalize models 
  final_glm <- finalize_model(
    glm_spec,
    best_mcc_glm
  )
  
  # Update workflow
  glm_wf_updated <- glm_wf %>%
    update_model(final_glm)
  
  # Fit the model
  glm_fit <- glm_wf_updated %>%
    fit(data = data)
  
  # Fit resamples
  glm_rs <- glm_fit %>%
    fit_resamples(
      resamples = folds,
      metrics = metrics,
      control = control # to obatin the prediction column
    )
  
  fin_glm <- Sys.time()
  
  # Save the final glm model
  if(!is.null(directory_to_save)){
    if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
    # Save models
    save(glm_fit, glm_rs, file = paste(directory,"/glm_model.Rda", sep=""))
  }

  model <- list(
      time = fin_glm-ini_glm,
      resamples = glm_rs,
      fit = glm_fit
  )
  
  return(model)
}

## k-Nearest Neighbours
knn_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

    # Set seed for reproducibility
    set.seed(1234)

    ini_knn <- Sys.time()

    # Define the model
    knn_spec <- nearest_neighbor(
        neighbors = tune(),
        dist_power = tune(),
        weight_func = tune()) %>%
        set_engine("kknn") %>%
        set_mode("classification")

    # Workflow
    knn_wf <- workflow() %>%
     add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
     add_model(knn_spec)

    # Collect parameters
    parameters_knn <- parameters(knn_spec)

    # Tuning
    knn_res <- tune_bayes(
        knn_wf,  # add workflow
        resamples = folds,
        param_info = parameters_knn,
        initial = 7, # This could also be a grid search
        iter = 25,
        metrics = metrics,
        control = control
    )

    # Select best selection of parameters
    best_mcc_knn <- select_best(knn_res, metric = "mcc") # Best based on mcc

    # Finalize models
    final_knn <- finalize_model(
        knn_spec,
        best_mcc_knn
    )

    # Update the workflow
    knn_wf_updated <- knn_wf %>%
        update_model(final_knn)

    # Fit the model
    knn_fit <- knn_wf_updated %>%
        fit(data = data)

    # Fit resamples
    knn_rs <- knn_fit %>%
        fit_resamples(
        resamples = folds,
        metrics = metrics,
        control = control # to obatin the prediction column
    )

    fin_knn <- Sys.time()

    # Save the final knn model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(knn_fit, knn_rs, file = paste(directory,"/knn_model.Rda", sep=""))
    }

    model <- list(
        time = fin_knn-ini_knn,
        resamples = knn_rs,
        fit = knn_fit
    )

  return(model)
}

## XGBoost
xgb_model <- function(data, target_var, folds, predictors, metrics, control, directory_to_save = NULL) {

    # Set seed for reproducibility
    set.seed(1234)

    ini_xgb <- Sys.time()

    # Define the model
    xgb_spec <- boost_tree(
      trees = 1000,
      tree_depth = tune(), min_n = tune(),
      loss_reduction = tune(),
      sample_size = tune(), mtry = tune(),
      learn_rate = tune()
    ) %>%
      set_engine("xgboost") %>%
      set_mode("classification")

    # Workflow
    xgb_wf <- workflow() %>%
        add_recipe(create_recipe(data, target_var)) %>%
        add_model(xgb_spec)

    # Collect model parameters
    parameters_xgb <- parameters(xgb_spec)
    #' It is necessary to finalize the range of values of mtry() that depend on the
    #' number of columns
    parameters_xgb <- xgb_spec %>% parameters() %>% finalize(predictors)

    # Tuning
    xgb_res <- tune_bayes(
        xgb_wf,
        resamples = folds,
        param_info = parameters_xgb,
        initial = 7,
        iter = 25,
        control = control,
        metrics = metrics
    )

    # Select best selection of parameters
    best_mcc_xgb <- select_best(xgb_res, metric = "mcc") # Best based on mcc

    # Finalize models
    final_xgb <- finalize_model(
        xgb_spec,
        best_mcc_xgb
    )

    # Update workflow
    xgb_wf_updated <- xgb_wf %>%
        update_model(final_xgb)

    # Fit the model
    xgb_fit <- xgb_wf_updated %>%
        fit(data = data)

    # Fit resamples
    xgb_rs <- xgb_fit %>%
        fit_resamples(
        resamples = folds,
        metrics = metrics,
        control = control # to obatin the prediction column
    )

    # Bundle xgb to save
    xgb_bundle_wf <- bundle(xgb_wf)
    xgb_bundle_res <- bundle(xgb_res)
    final_xgb_bundle <- bundle(final_xgb)
    xgb_bundle_fit <- bundle(xgb_fit)
    xgb_bundle_rs <- bundle(xgb_rs)
    # Save the final xgb model
    if(!is.null(directory_to_save)){
      if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
      # Save models
      save(xgb_bundle_fit, xgb_bundle_rs, file = paste0(directory,"/xgb_bundle_model.Rda"))
    }
    
    fin_xgb <- Sys.time()

    model <- list(
        time = fin_xgb-ini_xgb,
        resamples = xgb_rs,
        fit = xgb_fit
    )

  return(model)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Run ML models: function to run models

ML_models <- function(data, target_var, models_to_run = c("rf", "svm_r", "svm_p", "glm", "knn", "xgb"), directory_to_save = NULL) {
  
  # Set seed for reproducibility
  set.seed(1234)
  
  # Metrics
  metrics_to_save <- metric_set(mcc, bal_accuracy, accuracy, yardstick::sens, yardstick::spec,
                                precision, roc_auc, mn_log_loss)
  # Predictors
  expression_predictors <- select(data, -all_of(target_var))
  # Folds
  # expression_folds <- vfold_cv(data, strata = !!sym(target_var), repeats = 10)
  expression_folds <- vfold_cv(data, repeats = 10)
  
  ## Save the folds object
  if(!is.null(directory_to_save)){
    if ("ML_models"%in%list.files(paste0(directory_to_save)) == FALSE){
        dir.create(paste0(directory_to_save, "ML_models"), showWarnings = FALSE)}
      directory <- paste0(directory_to_save,"ML_models/")
     # Save folds
    saveRDS(expression_folds, file = paste0(directory,"/folds_object.rds"))
  }

  # Control
  ctrl_bayes <- control_bayes(no_improve = 15, verbose = TRUE, verbose_iter = TRUE, 
                              save_pred = TRUE, save_workflow = TRUE, seed = 1234)
  
  # Lists to save results
  results <- list()

  # Run models
  if ("rf" %in% models_to_run){
    results$rf <- rf_model(data, target_var, expression_folds, expression_predictors, metrics_to_save, ctrl_bayes, directory_to_save)
    print("RF done")
  }
  if("svm_r" %in% models_to_run){
    results$svm_r <- svm_r_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    print("SVM_r done")
  }
  if("svm_p" %in% models_to_run){
    results$svm_p <- svm_p_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    print("SVM_p done")
  }
  if("glm" %in% models_to_run){
    if(length(levels(data[[target_var]]))== 2){
      results$glm <- glm_binary_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
      print("GLM done")
    }else{
      results$glm <- glm_multiclass_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
      print("GLM done")
    }
  }
  if("knn" %in% models_to_run){
    results$knn <- knn_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    print("kNN done")
  }
  if("xgb" %in% models_to_run){
    results$xgb <- xgb_model(data, target_var, expression_folds, expression_predictors, metrics_to_save, ctrl_bayes, directory_to_save)
    print("XGB done")
  }

  # Computational time of models:
  return(results)
  
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# POSTPROCESSING: functions to obtain a .rds file with the results
# the results will be given into a R list

## Function to extract misclassified samples
extract_misclassified <- function(predictions, truth_col, pred_col) {
  misclassified <- predictions %>%
    filter(!!sym(truth_col) != !!sym(pred_col)) %>%
    rownames()
  return(misclassified)
}

extract_misclassified_cross_validation <- function(predictions, truth_col, pred_col, data) {
  sample_identifier <- predictions %>%
    filter(!!sym(truth_col) != !!sym(pred_col)) %>%
    select(c(fold,.row))
 
  misclassified <- as.data.frame(table(sample_identifier$.row))
  misclassified <- misclassified %>%
                        rename(row_id = Var1)
  rownames(misclassified) <- rownames(data[
                as.integer(levels(misclassified$row_id)),]) 
  return(misclassified)
}

## Function to extract fitted train and test results
extract_train_test_results <- function(fit, data, target_var) {
  
  # Predict using the fit object
  pred <- fit %>%
    predict(data) %>%
    bind_cols(data)
  
  # Confusion matrix
  matrix <- pred %>%
    conf_mat(truth = target_var, estimate = .pred_class)

  # Metrics
  metrics <- summary(matrix) # we can change the event_level if we want the metrics of the second level look https://yardstick.tidymodels.org/reference/summary.conf_mat.html
  
  # Misclassified samples
  misclassified_samples <- extract_misclassified(pred, target_var, ".pred_class")

  # Return results as a list
  results <- list(
    model_metrics = metrics,
    confusion_matrix = matrix,
    misclassified_samples = misclassified_samples
  )

  return(results)
}

## Function to extract cross validation results
extract_cross_validation_results <- function(classifier, data, target_var){

    cross_validation <- classifier
    pred<- classifier$.predictions %>%
     bind_rows(cross_validation, .id = "column_label")

    pred <-  pred %>% 
        rename(fold = column_label) %>%
        select(-ncol(pred))
    
    # Missclassified samples
    misclassified_samples <- extract_misclassified_cross_validation(pred, target_var, ".pred_class", data)
    
    metrics <- pred %>%
      group_by(fold) %>%
      conf_mat(!!sym(target_var), .pred_class) %>%
      mutate(summary_tbl = lapply(conf_mat, summary)) %>%
      unnest(summary_tbl) %>%
      group_by(.metric) %>%
      summarise(mean = mean(.estimate, na.rm = TRUE), 
                sd = sd(.estimate, na.rm = TRUE))
    
    matrix <- conf_mat_resampled(cross_validation)

      # Return results as a list
    results <- list(
        model_metrics = metrics,
        confusion_matrix = matrix,
        misclassified_samples = misclassified_samples
    )

    return(results)

  }

# Function to extract cross-validation, training, test results, and timing for each model
extract_models_results <- function(models_to_run, results_models, expression_train, expression_test, target_var) {
  
  # Initialize lists for storing results
  times <- list()
  train <- list()
  test <- list()
  cross_validation <- list()
  extracted_results <- list()

  # Loop through each model in the vector and extract results
  for (classif_name in models_to_run) {
    # Extract timing, resampling, and fit objects
    times[[classif_name]] <- results_models[[classif_name]]$time
    resample <- results_models[[classif_name]]$resamples
    fit <- results_models[[classif_name]]$fit

    # Extract cross-validation, train, and test results
    cross_validation[[classif_name]] <- extract_cross_validation_results(resample, expression_train, target_var)
    train[[classif_name]] <- extract_train_test_results(fit, expression_train, target_var)
    test[[classif_name]] <- extract_train_test_results(fit, expression_test, target_var)
  }

  # Store results in a list
  extracted_results$cross_validation <- cross_validation
  extracted_results$train <- train
  extracted_results$test <- test
  extracted_results$times <- times

  return(extracted_results)
}

