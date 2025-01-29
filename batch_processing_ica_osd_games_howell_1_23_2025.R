

#install.packages(c("raster","tidyverse","ALS","summarytools","ggplot2","userfriendlyscience", "car", "dplyr"), dependencies = TRUE)
# 
# if (!require("osd")) {
#   install.packages("osd")
# }
# if (!require("raster")) {
#   install.packages("raster")
# }
# 
# if (!require("tidyverse")) {
#   install.packages("tidyverse")
# }
# 
# if (!require("ALS")) {
#   install.packages("ALS")
# }
# 
# if (!require("summarytools")) {
#   install.packages("summarytools")
# }
# if (!require("ggplot2")) {
#   install.packages("ggplot2")
# }
# 
# if (!require("userfriendlyscience")) {
#   install.packages("userfriendlyscience")
# }
# if (!require("car")) {
#   install.packages("car")
# }
# if (!require("dplyr")) {
#   install.packages("dply")
# }

# if (!require("NMF")){
#   install.packages("NMF")
# }

library(raster)
# library(NMF)
library(tidyverse)
library(ALS)
library(summarytools)
library(ggplot2)
library(userfriendlyscience)
library(car)
library(dplyr)
library(osd)

grouping_VNIR_files <- function() {
  spectral_files <-
    list.files(
      path = 'E:/Shared With Me/rdcrlfjl/9K0WWPUH7btUp86f/RS_GIS_CCR_Project/VNIR',
      pattern = 'raw_rd_rf$',
      recursive = TRUE,
      full.names = TRUE
    )
  reference_dir_txt <-
    read.csv(
      'C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/reference_sheet_12_20024.csv')
  group_summary <- list('Dir', 'Group')
  # Creates a summary chart with each files directory and group ID
  for (file in spectral_files) {
    row <-
      reference_dir_txt[sapply(reference_dir_txt$Partial.directory, function(x)
        grepl(x, file, fixed = TRUE)),]
    appending <- list(file, row$GroupAssigned)
    group_summary = rbind(group_summary, appending)
  }
  # Groups by group ID
  dfgroup_summary <- as.data.frame(group_summary[-1, ])
  colnames(dfgroup_summary) <- group_summary[1, ]
  
  print('grouping_VNIR_files-- complete')
  assigning_control_v_treated(dfgroup_summary)
}

generate_treatment_ids <- function(control_id, grouped) {
  # print('check1')
  control_4_dirs<-grouped$Dir[grepl(control_id,grouped$Group)]
  # print('check2')
  treatment5_4_dirs<-grouped$Dir[grepl(sub("_0_", "_5_", control_id), grouped$Group)]
  # print('check3')
  treatment15_4_dirs<-grouped$Dir[grepl(sub("_0_", "_15_", control_id), grouped$Group)]
  # print('check4')
  treatment25_4_dirs<-grouped$Dir[grepl(sub("_0_", "_25_", control_id), grouped$Group)]
  # print('check5')
  return(list(control=control_4_dirs, treatment5=treatment5_4_dirs, treatment15=treatment15_4_dirs, treatment25=treatment25_4_dirs))
}

assigning_control_v_treated<-function(grouped){
  control_dirs<-grouped$Dir[grep("_0_", grouped$Group)]
  # Apply the function to all control IDs
  control_ids<-unique(grouped$Group[grep("_0_", grouped$Group)])
  i=1
  for (control_id in control_ids){
    control_dir<-grouped$Dir[grepl(control_id, grouped$Group)]
    mapped_data<-generate_treatment_ids(control_id, grouped)
    control<-mapped_data$control
    treatment5<-mapped_data$treatment5
    treatment15<-mapped_data$treatment15
    treatment25<-mapped_data$treatment25
    type5<-paste0(control_id,'_5ppm')
    type15<-paste0(control_id,'_15ppm')
    type25<-paste0(control_id,'_25ppm')
    launching_analysis(control_dir, treatment5, type5)
    launching_analysis(control_dir, treatment15, type15)
    launching_analysis(control_dir, treatment25, type25)
    i=i+1
  }
  print('assigning_control_v_treated--debugged')
}

launching_analysis <- function(group1, group2, type) {
  raster_group_summary <- matrix(1, nrow = 1, ncol = 59)
  # group= list of directories, raster_group_summary = blank summary matrix, type="control_vs_treatment_concentration"
  control <-
    open_combine_rasters(group1, raster_group_summary, 'control', type)
  treated <-
    open_combine_rasters(group2, raster_group_summary, 'treatment', type)
  print('Rasters have been combined')
  initial_comp <- complete_pca(control)
  cat('Initial component is: ', initial_comp)
  print('PCA success')
  control_mcr <- complete_mcr_als(control, initial_comp, type)
  cat('MCR component is: ', control_mcr$comp)
  print('Control MCR success')
  control_spectra <- control_mcr$spectra
  #write.csv(control_spectra,'D:/VNIR/control_spectral.csv')
  final_comp <- control_mcr$comp
  cat('final component is: ', final_comp)
  treated_mcr <- complete_mcr_als(treated, final_comp, '')
  print('Treated MCR success')
  treated_spectra <- treated_mcr$spectra

  output <-
    complete_howells_game(control_spectra, treated_spectra, type)

}

open_combine_rasters <-
  function(summary_per_group,
           raster_group_summary,
           type, type_group) {
    # Creates empty matrix that will follow the converted 2D matrix form of (col*row, nbands)
    # Since the sampled matrix will be 100*100*59, we will create 1,10000
    for (dir in summary_per_group) {
      # print(dir)
      raster <- stack(dir)
      mask <- creating_conditions(raster)
      # print('mask created')
      mask_red <- mask$red
      mask_ob_avg <- mask$ob
      mask_r_bg_avg <- mask$rbg
      
      reclassified_raster <-
        classify_raster_based_on_conditions(mask_red, mask_ob_avg, mask_r_bg_avg)
      subset <- stacking_selected_bands(raster)
      masked_subset <-
        applying_mask_on_subset(subset, reclassified_raster)
      # cropped raster should have dimensions of [100,100,59]--> [10000,59]
      cropped_raster <- applying_spatial_subset(masked_subset)
      cat('Dimension of 2D cropped df: ',dim(cropped_raster),'\n')
      #print(raster_group_summary)
      raster_group_summary <- rbind(raster_group_summary, cropped_raster)
      cat('Dimension of raster group summary df:', dim(raster_group_summary),'\n')
    }
    
    raster_group_summary <- raster_group_summary[-1, ]
    cat('Dimension of raster group summary df:',
        dim(raster_group_summary),
        '\n')
    # write.csv(
    #   raster_group_summary,
    # paste('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/20241220-reflectance/one_matrix', "_", type, ".txt"))
    write.csv(raster_group_summary, paste0("C:/Users/RDCRLAAU/Downloads/input_matrix_safeguard",type,"_",type_group, ".csv"))
    # matrix <- scale(raster_group_summary, center = TRUE, scale = TRUE) 
    # write.csv(matrix, paste0("C:/Users/RDCRLAAU/Downloads/scaled_matrix_", type,"_",type_group, ".csv"))
    return(raster_group_summary)
  }

applying_spatial_subset <- function(raster) {
  
  # Initialize the cumulative matrix
  cummulative_for_raster <- matrix(nrow = 0, ncol = 59)
  
  # Get the raster dimensions
  large_matrix <- raster[[1]]$layer.1
  nrow_large <- dim(raster)[1]
  ncol_large <- dim(raster)[2]
  
  # Calculate row and column sums
  rowSum <- rowSums(large_matrix$layer.1)
  colSum <- colSums(large_matrix$layer.1)
  
  # Identify the rows and columns with the largest sums
  largest_row <- order(rowSum, decreasing = TRUE)[1:1000]
  largest_col <- order(colSum, decreasing = TRUE)[1:1000]
  
  for (i in 1:1000) {
    # Current row and column indices
    col <- largest_col[i]
    row <- largest_row[i]
    
    # Define the desired extent
    xmin <- max(row - 10, 1)                 # Ensure extent does not go below 1
    xmax <- min(row + 9, nrow_large)         # Ensure extent does not exceed raster rows
    ymin <- max(col - 10, 1)                 # Ensure extent does not go below 1
    ymax <- min(col + 9, ncol_large)         # Ensure extent does not exceed raster columns
    
    # Adjust extent if it doesn't fit the raster borders
    extent_to_crop <- extent(raster, xmin, xmax, ymin, ymax)
    
    # Crop the raster using the adjusted extent
    cropped_raster <- raster::crop(raster, extent_to_crop)
    
    # Convert the cropped raster to 2D and add it to the cumulative matrix
    cropped_2d <- convert_to_2d(cropped_raster)
    cummulative_for_raster <- rbind(cummulative_for_raster, cropped_2d)
  }
  return(cummulative_for_raster)
}

creating_conditions <- function(Hyp_As_col3_16_3) {
  cat('creating conditions')
  #Averaging red between 750 and 850nm
  Red_avg <- mean(Hyp_As_col3_16_3[[180:204]])
  #Averaging orange between 550 and 600
  Orange_avg <- mean(Hyp_As_col3_16_3[[75:92]])
  #Averaging blue between 410 and 450nm
  Blue_avg <- mean(Hyp_As_col3_16_3[[10:25]])
  #Averaging for a blue-green band, in the paper it claimes between 3 and 450nm which makes no sense.
  #So I'm taking the liberty to make it between 450 and 550
  Blue_Green_avg <- mean(Hyp_As_col3_16_3[[30:69]])
  #Orange Average to Blue Average
  O_B_avg <- Orange_avg / Blue_avg
  
  #Red average to Blue-Green Average
  R_BG_avg <- Red_avg / Blue_Green_avg
  #Red Cutoff Threshold
  # Step 1: Calculate the histogram
  hist_values <-
    hist(getValues(Red_avg), breaks = 256, plot = FALSE)
  
  # Step 2: Implement Otsu's method
  # Compute probability distribution
  p <- hist_values$counts / sum(hist_values$counts)
  omega <- cumsum(p)
  mu <- cumsum(p * hist_values$mids)
  mu_t <- mu[length(mu)]
  
  # Compute between-class variance
  sigma_b_squared <- (mu_t * omega - mu) ^ 2 / (omega * (1 - omega))
  
  # Find the maximum value of between-class variance
  max_sigma_idx <- which.max(sigma_b_squared)
  otsu_thresh <- hist_values$mids[max_sigma_idx]
  
  # Step 3: Apply the threshold
  Red_cutoff <-
    calc(
      Red_avg,
      fun = function(x) {
        ifelse(x > otsu_thresh, 1, 0)
      }
    )
  # Plot to visualize
  # plot(Red_cutoff)
  return(list(red = Red_cutoff, ob = O_B_avg, rbg = R_BG_avg))
}


classify_raster_based_on_conditions <-
  function(red_cutoff_raster,
           orange_blue_avg_raster,
           red_bg_avg_raster) {
    # Define a function to apply to each cell
    reclassify_function <-
      function(red_cutoff,
               orange_blue_avg,
               red_bg_avg) {
        # Check the conditions and assign 1 or 0
        condition1 <- red_cutoff == 1
        condition2 <- (orange_blue_avg >= 1.75) | (red_bg_avg > 2.5)
        
        result <- ifelse(condition1 & condition2, 1, 0)
        return(result)
      }
    
    # Stack the rasters for easier cell-wise operations
    raster_stack <-
      stack(red_cutoff_raster,
            orange_blue_avg_raster,
            red_bg_avg_raster)
    
    # Apply the reclassification function to each cell
    reclassified_raster <- calc(
      raster_stack,
      fun = function(x) {
        reclassify_function(x[1], x[2], x[3])
      }
    )
    # plot(reclassified_raster)
    return(reclassified_raster)
  }

stacking_selected_bands <- function(raster) {
  #range<-list([10:25],[30:69], [75:92],[180:204])
  red_stack <- raster[[180:204]]
  orange_stack <- raster[[75:92]]
  blue_green_stack <- raster[[30:69]]
  blue_stack <- raster[[10:25]]
  subset_raster <- stack(red_stack, orange_stack, blue_stack)
  print('stacking_selected_bands--complete')
  return(subset_raster)
}

applying_mask_on_subset <- function(mask, raster) {
  masked <- mask * raster
  print('applying_mask_on_subset--complete')
  return(masked)
}

convert_to_2d <- function(hypercube) {
  matrix2D <-
    matrix(
      hypercube,
      nrow = dim(hypercube)[1] * dim(hypercube)[2],
      ncol = dim(hypercube)[3]
    )
  return(matrix2D)
}

complete_pca <- function(matrix) { 
  pca_result <- prcomp(matrix)
  
  # Extract eigenvalues
  eigenvalues <- pca_result$sdev ^ 2
  
  # # Calculate first-order derivatives (differences between consecutive eigenvalues)
  # first_order_derivative <- diff(eigenvalues)
  # 
  # number_components<-length(first_order_derivative[log(abs(first_order_derivative))>0])
  n <- length(eigenvalues)
  # Coordinates for the first and last points
  first_point <- c(1, eigenvalues[1])
  last_point <- c(n, eigenvalues[n])
  
  # Compute distances of each point to the line
  distances <- sapply(1:n, function(i) {
    point <- c(i, eigenvalues[i])
    # Distance from point to line
    abs((last_point[2] - first_point[2]) * point[1] -
          (last_point[1] - first_point[1]) * point[2] +
          last_point[1] * first_point[2] -
          last_point[2] * first_point[1]) /
      sqrt((last_point[2] - first_point[2])^2 + (last_point[1] - first_point[1])^2)
  })
  
  # Identify the index of the maximum distance
  number_components <- which.max(distances)
  return(number_components)
}


complete_mcr_als<-function(matrix2D, component, type){
  # write.csv(matrix2D, paste0("C:/Users/RDCRLAAU/Downloads/input_matrix_", type, ".csv"))
  # matrix <- scale(matrix2D, center = TRUE, scale = TRUE) 
  # write.csv(matrix, paste0("C:/Users/RDCRLAAU/Downloads/scaled_matrix_", type, ".csv"))
  cat("NA values?: ", any(is.na(matrix)), '\n')
  matrix[matrix==0]<-1e-6
  osd_result <- osd(D = matrix, k = component, res.method = "ica.osd")
  resolved_spectra<-osd_result$S
  resolved_concentration<-osd_result$C
  cat("Dimensions of resolved spectra: ", dim(resolved_spectra), '\n')
  cat("Dimensions of resolved concentration: ", dim(resolved_concentration), '\n')
  # endtime = Sys.time()
  # # print(endtime)
  # write.csv(resolved_concentration, 'C:/Users/RDCRLAAU/Downloads/resolved_concentration.csv')
  # write.csv(resolved_spectra, 'C:/Users/RDCRLAAU/Downloads/resolved_spectra.csv')
  reconstructed<-resolved_concentration%*%t(resolved_spectra)
  reconstructed[reconstructed==0]<-1e-6
  # write.csv(reconstructed, 'C:/Users/RDCRLAAU/Downloads/reconstructed_og.csv')
  reconstructed <- scale(reconstructed, center = TRUE, scale = TRUE)
  # write.csv(reconstructed, 'C:/Users/RDCRLAAU/Downloads/reconstructed_scaled.csv'
  significant_components<-complete_pca(reconstructed)
  # print(significant_components)
  return(list(comp = significant_components, spectra = resolved_spectra))
}

complete_howells_game<-function(control, treated, name){
  cat('-----Games-Howell Test-----', '\n')
  
  # Ensure that both final_spectra_avg and treated_spectra are vectors of the same length
  # Combine them into a single dataframe with a group label
  
  control_spectra <- data.frame(Spectra = rowMeans(control), Group = "Control")
  treated_spectra_df <- data.frame(Spectra = rowMeans(treated), Group = "Treated")
  
  # Combine the two datasets into one dataframe
  compiled_spectra <- rbind(control_spectra, treated_spectra_df)
  write.csv(compiled_spectra, paste('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/20241220-reflectance/', name,'_compiled_spectra.csv'))
  # Perform the Games-Howell test
  gh_results <- posthocTGH(y = compiled_spectra$Spectra, 
                           x = compiled_spectra$Group, 
                           method = "games-howell")
  
  # Extract mean and confidence intervals for plotting
  means <- tapply(compiled_spectra$Spectra, compiled_spectra$Group, mean)
  ci <- tapply(compiled_spectra$Spectra, compiled_spectra$Group, function(x) {
    c(
      mean(x) - qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x)),
      mean(x) + qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
    )
  })
  # 
  # # Prepare data for plotting
  # plot_data <- data.frame(
  #   Group = names(means),
  #   Mean = as.numeric(means),
  #   CI_Lower = sapply(ci, `[`, 1),
  #   CI_Upper = sapply(ci, `[`, 2)
  # )
  # 
  # # Create the plot
  # ggplot(plot_data, aes(x = Group, y = Mean)) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  #   labs(title = "Mean Spectral Values by Group with 95% Confidence Intervals",
  #        x = "Group", y = "Mean Spectral Value") +
  #   theme_minimal()
  # 
  
  # 4. Prepare a summary dataframe
  summary_df <- data.frame(
    Control = means[1],
    Treated = means[2],
    Diff = gh_results$output$games.howell$diff,
    CI_Lower = gh_results$output$games.howell$ci.lo,
    CI_Upper = gh_results$output$games.howell$ci.hi,
    P_Value = gh_results$output$games.howell$p
  )
  
  # 5. Save the summary to a CSV file
  write.csv(summary_df, paste('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/20241220-reflectance/', name, '_summary.txt'), row.names = FALSE)
}
start = Sys.time()
print(start)
grouping_VNIR_files()
end = Sys.time()
time = end-start
print(time)
