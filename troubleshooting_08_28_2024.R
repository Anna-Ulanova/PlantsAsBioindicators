# Load required package for MCR-ALS
if (!requireNamespace("ALS", quietly = TRUE)) {
  install.packages("ALS")
}
library(ALS)
library(raster)
start = Sys.time()
#______________________________________________________________________________#
# CONTROL- find components, find spectra, find reconstructive df
# This script will complete MCR-ALS
# https://cran.r-project.org/web/packages/ALS/ALS.pdf

# TEST CASE control- control and treated
control<-stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/Test Bands/14_1_20230724_2023_07_24_12_33_04/raw_rd_rf')
treated<-stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/MCR-ALS R Processing/Test Bands/4_1_20230724_2023_07_24_09_00_14/raw_rd_rf')

# Reshape the 3D data into a 2D matrix (spatial dimension flattened into rows, spectral bands as columns)
control_2d <- matrix(control, nrow = dim(control)[1] * dim(control)[2], ncol = dim(control)[3])

treated_2d <-matrix(treated, nrow = dim(treated)[1] * dim(treated)[2], ncol = dim(treated)[3])
# Perform PCA on the reshaped data
pca_result <- prcomp(control_2d, scale. = TRUE)

# Extract eigenvalues
eigenvalues <- pca_result$sdev^2

# Create a scree plot to visualize eigenvalues and determine the number of components
plot(eigenvalues, type = "b", xlab = "Component", ylab = "Eigenvalue", main = "Control Scree Plot")
abline(h = 1, col = "red", lty = 2)  # Kaiser's rule threshold

# Uses Kaiser-Guttman rule which suggest that "the elbow" in scree plots occurs
# when the eigenvalues dip below 1.
# https://datasciencewiki.net/kaisers-rule/
num_components <- sum(eigenvalues > 1)

# Load necessary libraries
# library(mcr) # Ensure you have the appropriate package
#-----------------------------------------------------------------------------#
# Step 1: Perform MCR-ALS on control_2d
# Initialize random matrices for MCR-ALS
initial_concentration_profiles_control <- matrix(runif(nrow(control_2d) * num_components), 
                                                 nrow = nrow(control_2d), 
                                                 ncol = num_components)
initial_spectra_control <- matrix(runif(num_components * ncol(control_2d)), 
                                  nrow = ncol(control_2d), 
                                  ncol = num_components)

# Run MCR-ALS on control_2d
mcr_result_control <- als(
  CList = list(initial_concentration_profiles_control), 
  S = initial_spectra_control, 
  PsiList = list(control_2d), 
  nonnegC = TRUE, 
  nonnegS = TRUE, 
  maxiter = 20, 
  thresh = 1e-6
)
# Extract the resolved spectral profiles from control_2d
resolved_spectra_control <- mcr_result_control$S
resolved_concentration_control<-mcr_result_control$C[[1]]
# Step 2: Use the number of components obtained to analyze treated_2d


control_reconstructed_spectra <- resolved_concentration_control %*% t(resolved_spectra_control)
#
# Compares observed_spectra, input to the MCR-ALS algorithm [n_samples, n_bands], to reconstructed results
# Example for plotting the original vs. reconstructed data
n_samples <- nrow(control_2d)
n_bands <- ncol(control_2d)

# Plot for the first sample
plot(1:n_bands, control_2d[1, ], type = "l", col = "red",
     xlab = "Wavelength/Band", ylab = "Intensity",
     main = "CONTROL Comparison of Original and Reconstructed Spectra")
lines(1:n_bands, control_reconstructed_spectra[1, ], col = "blue")
legend("topright", legend = c("Original", "Reconstructed"),
       col = c("red", "blue"))

#Complete PCA to find if the number of components has changed
# # Eigenvalue calculation (using concentration profiles as an example)
eigenvalues <- eigen(cov(control_reconstructed_spectra))$values

# Kaiser’s rule: Select components with eigenvalues > 1
significant_components <- which(eigenvalues > 1)


# Initialize random matrices for MCR-ALS
initial_concentration_profiles_treated <- matrix(runif(nrow(treated_2d) * significant_components), 
                                                 nrow = nrow(treated_2d), 
                                                 ncol = significant_components)
initial_spectra_treated <- matrix(runif(significant_components * ncol(treated_2d)),
                                  nrow = ncol(treated_2d),
                                  ncol = significant_components)

# Run MCR-ALS on treated_2d
mcr_result_treated <- als(
  CList = list(initial_concentration_profiles_treated), 
  S = initial_spectra_treated, 
  PsiList = list(treated_2d), 
  nonnegC = TRUE, 
  nonnegS = TRUE, 
  thresh = 1e-8
)

# Extract the resolved concentration profiles from treated_2d
resolved_concentrations_treated <- mcr_result_treated$CList[[1]]
resolved_spectra_treated<-mcr_result_treated$S
treated_reconstructed_spectra <- resolved_concentrations_treated %*% t(resolved_spectra_treated)

n_samples <- nrow(treated_2d)
n_bands <- ncol(treated_2d)

# Plot for the first sample
plot(1:n_bands, treated_2d[1, ], type = "l", col = "red",
     xlab = "Wavelength/Band", ylab = "Intensity",
     main = "TREATED Comparison of Original and Reconstructed Spectra")
lines(1:n_bands, treated_reconstructed_spectra[1, ], col = "blue")
legend("topright", legend = c("Original", "Reconstructed"),
       col = c("red", "blue"), lty = 1)


