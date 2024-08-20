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
complete_mcr_als<-function(convert2D, pca_dataset, num_comp, length, width, num_bands){
  # 1. Estimate number of spectral components using randomly generated spectra and 
  # additional baseline component (constant offset) by using MCR-ALS. Estimate the
  # required number of components
  
  # Initial estimate of pure spectra using the PCA components satisfying Kaiser's rule
  initial_spectra<-pca_dataset$rotation[,1:num_comp]
  
  # Create initial concentration profile estimates
  initial_concentration_profiles <- matrix(
    runif(nrow(convert2D) * num_comp), 
    nrow = nrow(convert2D), 
    ncol = num_comp
  )
  # 2. Rerun MCR-ALS with adjusted number of spectral components until model
  # developed appropriate number of spectral components. Convergence is 
  
  # Define convergence criteria
  tolerance <- 1e-6
  converged <- FALSE
  iteration <- 0
  
  # Initialize previous concentration profiles for convergence checking
  previous_concentration_profiles <- initial_concentration_profiles
  

  while (!converged) {
    # Run MCR-ALS
    mcr_result <- als(
      CList = list(previous_concentration_profiles),  # Use the previous iteration's concentration profiles
      PsiList = list(convert2D),                        # The actual data matrix
      S = initial_concentration_profiles,             # Initial spectra estimate
      nonnegS = TRUE,                                 # Enforce non-negativity for spectra
      nonnegC = TRUE,                                 # Enforce non-negativity for concentration profiles
      maxiter = 10,                                   # Maximum number of iterations per ALS call
      thresh = tolerance                              # Convergence threshold
    )
    # Extract resolved spectra and concentration profiles
    current_spectra <- mcr_result$S
    concentration_profiles <- mcr_result$CList[1]
    
    # Check for convergence: compute the Frobenius norm of the difference
    spectral_diff <- sqrt(sum((current_spectra - previous_spectra)^2))
    
    # Check convergence criteria
    if (spectral_diff < tolerance) {
      converged=TRUE
    }
    else {
      previous_spectra <- current_spectra
    }
  }    
  
  # Plot resolved pure spectra after convergence
  matplot(t(current_spectra), type = "l", lty = 1, col = 1:num_components_kaiser, xlab = "Wavelength Index", ylab = "Intensity")
  legend("topright", legend = paste("Component", 1:num_components_kaiser), col = 1:num_components_kaiser, lty = 1)
  
  # Plot concentration profiles for a subset of samples
  matplot(concentration_profiles[1:100, ], type = "l", lty = 1, col = 1:num_components_kaiser, xlab = "Sample Index", ylab = "Concentration")
  legend("topright", legend = paste("Component", 1:num_components_kaiser), col = 1:num_components_kaiser, lty = 1)
}

complete_pca(pca_dataset)
