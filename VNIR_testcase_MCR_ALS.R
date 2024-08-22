# Load required package for MCR-ALS
if (!requireNamespace("ALS", quietly = TRUE)) {
  install.packages("ALS")
}
library(ALS)
library(raster)
start = Sys.time()
# This script will complete MCR-ALS
# https://cran.r-project.org/web/packages/ALS/ALS.pdf

# TEST CASE VNIR
vnir<-stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/PlantsAsBioindicators-Python/VNIR/16_3_20230724_2023_07_24_13_10_38/raw_rd_rf')

# Check dimensions of vnir to ensure it's a 3D array
dim(vnir)

# Reshape the 3D data into a 2D matrix (spatial dimension flattened into rows, spectral bands as columns)
vnir_2d <- matrix(vnir, nrow = dim(vnir)[1] * dim(vnir)[2], ncol = dim(vnir)[3])

# Perform PCA on the reshaped data
pca_result <- prcomp(vnir_2d, scale. = TRUE)

# Extract eigenvalues
eigenvalues <- pca_result$sdev^2

# Create a scree plot to visualize eigenvalues and determine the number of components
plot(eigenvalues, type = "b", xlab = "Component", ylab = "Eigenvalue", main = "Scree Plot")
abline(h = 1, col = "red", lty = 2)  # Kaiser's rule threshold

# Uses Kaiser-Guttman rule which suggest that "the elbow" in scree plots occurs
# when the eigenvalues dip below 1. 
# https://datasciencewiki.net/kaisers-rule/
num_components <- sum(eigenvalues > 1)

loadings <-
  pca_result$rotation[, 1:num_components]  # Loadings (corresponding to pure spectra)

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
    # Non-negativity constraint on concentrations

  nonnegC = TRUE,
  nonnegS = TRUE,
  # Non-negativity constraint on spectra
  maxiter = 200,
  # Maximum number of iterations
  thresh = 1e-6
)

# mcr_result$C contains the estimated concentration profiles for each component.
# [n_samples, n_components]
# Each column in this matrix corresponds to the concentration profile of one
# component across all samples. These profiles indicate how the concentration of
# each pure component varies across the different samples.
final_concentrations <- mcr_result$C[[1]]
# mcr_result$S contains the estimated spectral profiles (or pure spectra) for each component.
# [n_bands, n_components]
# Purpose: Each column in this matrix corresponds to the spectral profile of one pure component. These profiles represent the characteristic spectra of each pure component across all the spectral bands.
final_pure_spectra <- mcr_result$S

reconstructed_spectra <- final_concentrations %*% t(final_pure_spectra)
# 
# Compares observed_spectra, input to the MCR-ALS algorithm [n_samples, n_bands], to reconstructed results
# Example for plotting the original vs. reconstructed data
n_samples <- nrow(vnir_2d)
n_bands <- ncol(vnir_2d)

# Plot for the first sample
plot(1:n_bands, vnir_2d[1, ], type = "l", col = "red",
     xlab = "Wavelength/Band", ylab = "Intensity",
     main = "Comparison of Original and Reconstructed Spectra")
lines(1:n_bands, reconstructed_spectra[1, ], col = "blue")
legend("topright", legend = c("Original", "Reconstructed"),
       col = c("red", "blue"), lty = 1)


# Complete PCA to find if the number of components has changed
# Eigenvalue calculation (using concentration profiles as an example)
eigenvalues <- eigen(cov(final_pure_spectra))$values

# Kaiser’s rule: Select components with eigenvalues > 1
significant_components <- which(eigenvalues > 1)

print(Sys.time()-start)