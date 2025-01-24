#########################################################################################
# Library with functions for plotting ML results (RadarCharts)
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library("fmsb"); library(dplyr); library(ggplot2)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to prepare data for RadarChart
prepare_data_for_radarchart <- function(results, metric, estimate_column, min = 0, max = 1) {

  # Subset the data for the specific metric
  results_metric <- results[which(results$metric == metric), ]
  classifiers <- unique(results_metric$classifier)
  mlinput <- unique(results_metric$input_list)

  # Initialize matrix to store values for RadarChart
  values_radarchart <- matrix(0, nrow = length(mlinput) + 2, ncol = length(classifiers))
  rownames(values_radarchart) <- c("max", "min", mlinput)
  colnames(values_radarchart) <- classifiers

  #
  best_selection <- rep(7, length(classifiers))
  best_value <- numeric(length(classifiers))

  # Iterate over classifiers
  for (a in seq_along(classifiers)) {

    values <- numeric(length(mlinput))

    # Iterate over ML inputs
    for (b in seq_along(mlinput)) {

      # Save classifier row numbers
      classifier_row_position <- which(results_metric$classifier == classifiers[a])
      # Save input_list row numbers
      mlinput_row_position <- which(results_metric$input_list == mlinput[b])
      # Save the value for the corresponding classifier and ml input
      values[b] <- as.numeric(results_metric[[estimate_column]][
        intersect(classifier_row_position, mlinput_row_position)])

    }

    # Identify the ML input with the maximum value for the corresponding classifier
    values <- replace(values, is.na(values), 0) # NA metric values
    mejorvalor <- which(values == max(values))

    # Chech if the value is the same for several classifiers and update best 
    # classifier and best value selections if needed
    if (length(mejorvalor) == 1) {

          best_selection[a] <- mlinput[mejorvalor]
          best_value[a] <- values[mejorvalor]

        } else if (length(mejorvalor) > 1) {

          best_value[a] <- values[mejorvalor[1]]

    }

    # Fill in the RadarChart matrix
    values_radarchart[, a] <- c(max, min, values)

  }

  return(list(values_radarchart = as.data.frame(values_radarchart), best_selection = best_selection, best_value = best_value))
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to plot RadarChart
plot_radarchart <- function(data, metric, plot_name, output_dir, min, max, seq) {

  values_radarchart <- data$values_radarchart
  best_selection <- data$best_selection
  best_value <- data$best_value
  # Needed becausethe radarchart selects the colors from a vector of colors in order from the first one
  radar_colors <- genes_list_colors[rownames(values_radarchart)[3:nrow(values_radarchart)]]

  # Plot RadarChart
  pdf(file = file.path(output_dir, paste(metric,"_", plot_name, ".pdf", sep = "")), width = 12, height = 10)
  radarchart(values_radarchart, axistype = 5, 
              # Custom the polygons
              pcol = radar_colors, plwd = 3, plty = 1,
              # Custom the grid
              cglcol = "grey", 
              cglwd = 2, 
              axislabcol = radar_colors[best_selection], 
              paxislabels = round(best_value, 3), 
              caxislabels=seq(min,max,seq), 
              calcex = 0.5,
              # Custom labels
              vlcex = 1)
  legend(x = 0.8, y = 1.3, legend = names(radar_colors), bty = "n", pch = 20, col = radar_colors, 
      text.col = "black", cex = 0.6, pt.cex = 1.5)
  dev.off()

}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Function to generate summary barplots
summary_barplot <- function(data, metric, ml_models_colors, output_file, mean_estimate) {

  # Convert necessary columns to numeric
  data <- convert_columns_to_numeric(data, c(mean_estimate))
  
  # Generate results summary
  results_summary <- data %>% 
    filter(metric == metric) %>%
    group_by(classifier) %>% summarise(mean_value = mean(!!sym(mean_estimate)), sd = sd(!!sym(mean_estimate)), median = median(!!sym(mean_estimate)))

  # Reorder the levels of 'classifier' based on 'median'
  results_summary$classifier <- factor(results_summary$classifier, 
                                       levels = results_summary$classifier[order(results_summary$mean_value, decreasing = TRUE)]) 

  # Plotting
  p <- ggplot(results_summary, aes(x = classifier)) +
    geom_bar(aes(y = mean_value, fill = classifier), stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_value - sd, ymax = mean_value + sd), width = 0.2, position = position_dodge(width = 0.8)) +
    geom_point(aes(y = median), color = "#669bbc", size = 3) +
    geom_text(aes(y = mean_value + sd, label = round(sd, 2)), 
              vjust = -0.5, size = 3.5, color = "black") + # Add SD annotations
    labs(title = "Classifier Comparison",
         y = "Mean Value",
         x = "Classifier") +
    theme_minimal() + 
    scale_fill_manual(values = ml_models_colors)
  
  # Save the plot
  ggsave(filename = paste0(output_file,".pdf"), plot = p, device = "pdf")
}
