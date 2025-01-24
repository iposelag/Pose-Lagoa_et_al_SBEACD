
# Function to merge lists
merge_results_models <- function(...) {
  # Capture all lists passed as arguments
  list_of_lists <- list(...)
  
  # Get all unique keys from all the lists
  keys <- unique(unlist(lapply(list_of_lists, names)))
  
  # Initialize an empty list to store merged results
  merged <- list()
  
  # Merge the sub-elements under each key
  for (key in keys) {
    merged[[key]] <- unlist(lapply(list_of_lists, `[[`, key), recursive = FALSE)
  }
  
  return(merged)
}

# Read the results_models object
list <- "disease_related_entire_list"
results <- "shap"
list1 <- readRDS(paste0("../COPD/results_",list,"_1/results_",results,".rds"))
list2 <- readRDS(paste0("../COPD/results_",list,"_2/results_",results,".rds"))
list3 <- readRDS(paste0("../COPD/results_",list,"_3/results_",results,".rds"))
list4 <- readRDS(paste0("../COPD/results_",list,"_4/results_",results,".rds"))
list5 <- readRDS(paste0("../COPD/results_",list,"_5/results_",results,".rds"))
list6 <- readRDS(paste0("../COPD/results_",list,"_6/results_",results,".rds"))

# Merge the two lists
results_merged <- merge_results_models(list1, list2,list3, list4,list5, list6)

# Print the merged structure
print(names(results_merged))

# Save the results
directory_to_save <- paste0("../COPD/results_",list,"/")
saveRDS(results_merged, file = paste0(directory_to_save,"results_",results,".rds"))
