# Systems Biology and Explainable AI Ensembles to Uncover Complex Diseases
Iria Pose-Lagoa, Beatriz Urda-García, Jon Sánchez-Valle, Jose Carbonell-Caballero, Alfonso Valencia

### Abstract
Human complex diseases, such as Chronic Obstructive Pulmonary Disease (COPD), pose significant challenges due to their inherent heterogeneity. To enhance the identification of novel candidate genes associated with COPD, this study
presents a novel computational pipeline based on an ensemble of explainable AIs. By integrating different feature selection methods, AI classifiers, and explainability techniques, the study analyzed gene expression data from the Lung Tissue Research Consortium to stratify patients and identify biologically meaningful patterns. The pipeline combines data-driven signals and curated COPD-related gene sets, expanded with known biological interactors that exhibit lower signal-to-noise ratios due to molecular heterogeneity but could play significant roles in the disease. To ensure specificity, explainability scores were used to assess gene contributions across multiple classifier configurations and select the final list of candidate genes. Our results demonstrated superior accuracy compared to previous studies in discriminating between COPD and control samples. Additionally, our final list of candidate genes accurately recapitulated key biological processes in the disease, including immune responses and extracellular matrix degradation, while also highlighting other lesser-known properties, such as the relationships between COPD and cancer, and hormone regulation in lung function. Finally, aggregated explainability scores further identified the presence of potential new molecular subtypes, paving the way for better characterization of the underlying molecular mechanisms of COPD and laying a foundation for improved diagnostics and personalized treatments.

<!-- CODE -->
## Code structure

The code is structure as follows:

```bash
├── COPD/raw_data/
    │
    ├──dea.txt
    ├──expression.Rda
    └──disgenet_tables/
            ├── C0024117_disease_gda_summary.tsv
            └── disgenet_curated_gene_disease_associations.tsv

└── src/   
    │
    ├── main.R/  
    │
    ├── feature_selection/
    │   │  
    │   └── library_feature_selection.R    <- library with functions for the feature selection process
    │ 
    ├── ML_models/   
    │   │  
    |   └──	library_feature_selection.R   <- library with function for the trainining of the models
    │    
    ├── explainability/   
    │   │  
    |   └──	library_explainanility.R   <- library with function for the explainability of the models
    |
    ├── analysis/
    │   ├── library_analysis.R
    │   ├── plots.R
    |   └──	processate_metrics.R  
    |     
    ├── support/
    |   ├── merge_resulst_models.R
    |
    ├── analysis_results_ML_models.R
    ├── aggregated_normalized_shap.R
    ├── candidate_genes_analytical.R
    ├── candidate_genes_mclust.R
    ├── input_lists_barplot.R
    ├── subgroups_analysis_shap.R
    ├── upset_plot_top_shap.R
    └── venn_diagram_seed_lists.R
```

## main.R
Example of command line:  
```
Rscript main.R --directory_to_load "../COPD/raw_data" \
               --file_name "expression" \
               --target_var "dis_condition" \
               --directory_to_save "../COPD/results" \
               --procedure "disease_related" \
               --disease_code "C0024117" \
               --ml_models_to_run "rf" "glm
```
The main script gives as results two .rds files: 

- results_models.rds
- results_shap.rds

> NOTE: Other data results from intermidiate steps (train and test sets, genes lsits, folds object, ml models objects, extended shap samples results …) can be saved providing `directory_to_save` argument. See below.
> 

The script can be run from command line. There are required arguments and optional arguments. The arguments are the following ones:

- Required:
    
    `--directory_to_load` Raw data folder for expression data. *(ex. "../COPD/raw_data/")*
    `--file_name`   Name of the expression data. It should be an .Rda file. It is only required the name without the extension. *(ex. "expression")*
    `--target_var`  Name of the target variable. *(ex. “dis_condition”)*
    `--procedure` Feature selection procedure. *(ex. “disease_related”)*
    `--ml_models_to_run` List of models to run. *(ex. “rf” “glm”)*
    
- Optional arguments:
    
    `--directory_to_save` Folder to save results. *(ex. “../COPD/results/”)*
    
    `--threshold_value` numeric threshold for selecting the top mrmr genes if you already know it, and avoid makin the randomizations (ex. 0.031)
    
    `--disease_code` Disgenet code for the disease *(ex. "C0024117")*
    
    `--dea_genes` Directory for your differential expressed genes list. It should be a .txt files with one gene per line. *(ex. ../COPD/raw_data/dea.txt)*
    

The main script has 4 functions:

### 0. SPLIT DATA

It loads the expression file and split it into training and test. The expression file has to be saved as a data.frame with samples as rows and genes as columns. One column should have the condition of the samples. The extension of the files should be .Rda (ex. “expression.Rda”). This function takes as input:

```r
# REQUIRED:
# 'directory_to_load: raw data folder where you have the expression data (ex. "../COPD/raw_data/")
# 'file_name: the naming of your expression file (ex. "expression")
# 'target_var: name of the target variable (ex. "disease_condition")
# OPTIONAL:
# 'directory_to_save: where you want to save the data (ex. "../COPD/results")
```

### 1. Feature selection step

It perform the feature selection process of interest: “*mrmr*” , “*disease_related*”, “*disease_related_entire_list*”, “*data_driven*”, “*omnipath_disease_related*”, “*omnipath_data_driven”*, “*omnipath_intersection”*, “*omnipath_union”*

The function takes as input:

```r
# REQUIRED:
# 'procedure: ex. "disease_related"
# 'expression_data: expression train data (ex. "expression_train")
# 'target_var: name of the target variable (ex. "disease_condition")
# OPTIONAL: 
# 'threshold_value: numeric thresold for select the top mrmr genes if you already know it
# 'disease_code: disease disgenet code, if needed (ex. "C0024117")
# 'dea_genes: a vector with the differentially expressed genes
# 'directory_to_load: raw data folder where you have the expression data (ex. "../COPD/raw_data/")
# 'directory_to_save: where you want to save the data (ex. "../COPD/results")
```

Notes on procedures:

- MRMR
    
    Note that the randomization *mrmr* performs 1000 iterations of the algorithm to extract the `threshold_value`. If you already know the `threshold_value` for selecting the top genes you can avoid performing the randomization given a value to the variable `threshold_value` .
    
- Disease related
    
    It is needed that the curated table of disease is saved inside the `directory_to_load` folder as: ”disgenet_tables/disgenet_curated_gene_disease_associations.tsv”
    
- Disease related entire list
    
    It is needed to have the entire disgenet table related with the disease saved inside the folder `directory_to_load` folder as: “disgenet_tables/", `disease_code`, "_disease_gda_summary.tsv”
    
- Data-driven
    
    It is required to load the list of differentially expressed genes as a vector
    
- Omnipath expansions
    
    Omnipath requires internet connection. Moreover, each expansion will run the required seed list first, and then, compute the expansion.
    

⚠️It is still pending the integration of the Differential Expression Analysis and the expansion made with GUILDify.

### 2. Classification performance

This function perform the training of the models. It takes as input the following variables:

```r
# REQUIRED:
# 'data: expression train data (ex. "expression_train")
# 'target_var: name of the target variable (ex. "disease_condition")
# 'models_to_run: vector os strings with ml models to run (rf, knn, svm_r, svm_p, glm, xgb) (ex. c("rf","knn"))
# OPTIONAL:
# 'directory_to_save: where you want to save the data (ex. "../COPD/results")
```

There are another extra function `extracted_results` for extract the cross_validation, train, test results of the models. The results is a list of lists as follows:

- Cross_validation
    - Metrics
    - Confusion matrix
    - Misclassified samples
- Train
    - Metrics
    - Confusion matrix
    - Misclassified samples
- Test
    - Metrics
    - Confusion matrix
    - Misclassified samples

### 3. Explainability

This function extract the importance value for each sample and each variable by model. It takes as input the following variables:

```bash
# REQUIRED
# 'ml_models_to_run_vector: models to run vector (ex. "rf" "glm")
# 'results_models: the results from the classification performance function
# 'expression_data: expression train data (ex. "expression_train")
# 'target_var: name of the target variable (ex. "disease_condition")
# OPTIONAL
# 'directory_to_save
```

The results are given as a list of lists as follows:

- Model 1 (RF)
    - shap_results: the shap results are given as a unique importance value for sample and gene.

## 4. ML analysis and generation of candidate list of genes

We create RadarCharts with the prediction performance across all scenarios: `analysis_results_ML_models.R`.
We analyse the shap values, create aggregated metrics and obtain the candidate list of genes: `aggregated_normalized_shap.R`, `candidate_genes_analytical.R`, `candidate_genes_mclust.R`

