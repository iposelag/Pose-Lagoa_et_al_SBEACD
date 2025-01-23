# Systems Biology and Explainable AI Ensembles to Uncover Complex Diseases
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
    │  
    └── support/ 
        │
        └── library_help.R     <- library with other functions
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
