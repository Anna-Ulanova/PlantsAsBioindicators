# Usage of MCR-ALS for creating unique signatures for bioindicating plants
# This script will complete MCR-ALS
# https://cran.r-project.org/web/packages/ALS/ALS.pdf
# als(CList, PsiList, S=matrix(), WList=list(), thresh =.001, maxiter=100,
    # forcemaxiter = FALSE,optS1st=TRUE, x=1:nrow(CList[[1]]), x2=1:nrow(S),
    # baseline=FALSE, fixed=vector("list", length(PsiList)),
    # uniC=FALSE, uniS=FALSE, nonnegC = TRUE, nonnegS = TRUE,
    # normS=0, closureC=list())
# CList: A list of concentration profiles. This is where you input your initial 
# estimates for the concentration profiles (usually set to random values initially 
# if you do not have better estimates).
 
# PsiList: A list of spectra matrices for different datasets. If you have multiple
# datasets, each matrix should be one element of the list.
 
# S: The initial estimate for the pure component spectra. This is where you would 
# input your initial spectra estimates, which can be generated randomly or from 
# other methods like PCA.

# WList: A list of weights for different datasets, which can be used to emphasize
# certain datasets more than others during optimization.


library(ALS)

# Requires 2D input data used for initial PCA, PCA data, number of components, dimensions of the 
complete_mcr_als<-function(vnir_2d, loadings, num_components){
  # 1. Estimate number of spectral components using randomly generated spectra and 
  # additional baseline component (constant offset) by using MCR-ALS. Estimate the
  # required number of components
  
  
  # Create initial concentration profile estimates
  initial_concentration_profiles <- matrix(
    runif(nrow(vnir_2d) * num_components), 
    nrow = nrow(vnir_2d), 
    ncol = num_components)
  
  initial_pure_spectra<-loadings
  
  
  mcr_result <- als(
    CList<-list(initial_concentration_profiles), 
    PList<-list(vnir_2d),
    S<-initial_pure_spectra,
    nonnegC = TRUE,
    # Non-negativity constraint on concentrations
    nonnegS = TRUE,
    # Non-negativity constraint on spectra
    maxiter = 100,
    # Maximum number of iterations
    thresh = 1e-6
  )
  
  # mcr_result$C contains the estimated concentration profiles for each component.
  # [n_samples, n_components]
  # Each column in this matrix corresponds to the concentration profile of one
  # component across all samples. These profiles indicate how the concentration of
  # each pure component varies across the different samples.
  # final_concentrations <- mcr_result$C[[1]]
  # mcr_result$S contains the estimated spectral profiles (or pure spectra) for each component.
  # [n_bands, n_components]
  # Purpose: Each column in this matrix corresponds to the spectral profile of one pure component. These profiles represent the characteristic spectra of each pure component across all the spectral bands.
  final_pure_spectra <- mcr_result$S
  
  # reconstructed_spectra <- final_concentrations %*% t(final_pure_spectra)
  # 
  # # Compares observed_spectra, input to the MCR-ALS algorithm [n_samples, n_bands], to reconstructed results
  # # Example for plotting the original vs. reconstructed data
  # n_samples <- nrow(vnir_2d)
  # n_bands <- ncol(vnir_2d)
  # # Plot for the first sample
  # plot(1:n_bands, vnir_2d[1, ], type = "l", col = "red", 
  #      xlab = "Wavelength/Band", ylab = "Intensity",
  #      main = "Comparison of Original and Reconstructed Spectra")
  # lines(1:n_bands, reconstructed_spectra[1, ], col = "blue")
  # legend("topright", legend = c("Original", "Reconstructed"), 
  #        col = c("red", "blue"), lty = 1)
  # 
  # Complete PCA to find if the number of components has changed
  # Eigenvalue calculation (using concentration profiles as an example)
  eigenvalues <- eigen(cov(final_pure_spectra))$values
  
  # Kaiser’s rule: Select components with eigenvalues > 1
  significant_components <- which(eigenvalues > 1)
  return(significant_components, c)
}

sig_components<-complete_pca(pca_dataset)[1]
final_spectra<-complete_pca(pca_dataset)[2]
