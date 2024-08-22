# Usage of MCR-ALS for creating unique signatures for bioindicating plants
# Main script that will call on individual ones
source('mcr_als_processing_PCA.R')
source('mcr_als_processing_MCR_ALS.R')
source('mcr_als_processing_convert2D.R')

library(raster)
library(ALS)

start<-Sys.time()

launch_main<-function(control, treated){
  # 1.Create mask to separate plant vs. non-plant based on the averaged reflectance
  # values of the red/orange/blue/green bands and their respective ratios
  # *As of 8/14/2024, task has been delegated to Franz Lichtner
  
  
  
  # 2. Apply mask onto the hypercube, thereby classifying which plant vs. background
  # *As of 8/14/2024, task has been delegated to Franz Lichtner
  # - Input: hypercube
  # - Output: masked "unfolded" matrix (ei. converted to a 2D structure)
  
  
  # 3. Perform PCA to find the number of spectral components# 
  # *As of 8/14/2024, task has been delegated to Anna Ulanova
  # - Input: control matrix
  # - OUtput: number of spectral components
  control_pca_output<-complete_pca(control_2D)
  num_components<-control_pca_output[1]
  control_loadings<-control_pca_output[2]
  control_pca<- control_pca_output[3]
  # 4. Performs MCR-ALS on the control until the datasets converge 
  # *As of 8/14/2024, task has been delegated to Anna Ulanova
  # - Input: number of components
  # - Output: calibrated number of components
  control_mcr_als_output<-complete_mcr_als(control_2D, control_loadings, num_components)
  control_sig_components<-control_mcr_als_output[1]
  control_final_pure_spectra<-control_mcr_als_output[2]


  # Performs PCa on the treated spectral dataset to find the pure spectra 
  treated_loadings<-complete_pca(treated_2D)[2]
  
  # Performs MCR-ALS on the treated dataset using the calibrated number of components
  treated_pca_output<-complete_mcr_als(treated_2D, treated_loadings, control_sig_components)
  treated_components<-treated_pca_output[1]
  treated_control_sig_components <- treated_pca_output[2]
  print('Execution time: ', sys.time()-start)
}


#-----------------CONTROL-----------------#
control <- stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/Test Bands/14_1_20230724_2023_07_24_12_33_04/raw_rd_rf')

  
  
#-----------------TREATED-----------------#
treated <- stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/Test Bands/4_1_20230724_2023_07_24_09_00_14/raw_rd_rf')


