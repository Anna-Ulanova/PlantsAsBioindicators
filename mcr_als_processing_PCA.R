# Usage of MCR-ALS for creating unique signatures for bioindicating plants
# This script will complete PCA and report the eigenvalues of the hyperspectral 
# data
install.packages(raster)
library(raster)


complete_pca<-function(convert2D){
  pca_dataset<-prcomp(convert2D)
  print(pca_dataset)
  # Uses Kaiser-Guttman rule which suggest that "the elbow" in scree plots occurs
  # when the eigenvalues dip below 1. 
  # https://datasciencewiki.net/kaisers-rule/
  eigenvalues<-pca_dataset$sdev^2
  num_components<-sum(eigenvalues>1)
  # Use the number of components determined by PCA or manually adjust based on domain knowledge
  num_components <- ifelse(num_components > 1, num_components, 2)  # Ensure at least 2 components
  loadings <-
    pca_dataset$rotation[, 1:num_components]  # Loadings (corresponding to pure spectra)
  return(length(num_components), loadings, pca_dataset)
}


#------------Test Cases------------#

vnir<-stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/PlantsAsBioindicators-Python/VNIR/5_2_20231114_2023_11_14_07_51_18/raw_rd_rf')
length<-dim(vnir)[1]
width<-dim(vnir)[2]
num_bands<-dim(vnir)[3]
convert2D<-matrix(vnir, nrow=length*width,ncol=num_bands)
num_comp<-complete_pca(convert2D)[1]
loadings<-complete_pca(convert2D)[2]
pca_dataset<-complete_pca(convert2D)[3]