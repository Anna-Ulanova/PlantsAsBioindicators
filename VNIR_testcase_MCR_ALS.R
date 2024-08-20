# Load required package for MCR-ALS
if (!requireNamespace("ALS", quietly = TRUE)) {
  install.packages("ALS")
}
library(ALS)
library(ggplot2)
# This script will complete MCR-ALS
# https://cran.r-project.org/web/packages/ALS/ALS.pdf

# TEST CASE VNIR
vnir<-stack('C:/Users/RDCRLAAU/Desktop/Plant_as_bioindicators/PlantsAsBioindicators-Python/VNIR/5_2_20231114_2023_11_14_07_51_18/raw_rd_rf')

# Check dimensions of vnir to ensure it's a 3D array
dim(vnir)

# Reshape the 3D data into a 2D matrix (spatial dimension flattened into rows, spectral bands as columns)
vnir_2d <- matrix(vnir, nrow = dim(vnir)[1] * dim(vnir)[2], ncol = dim(vnir)[3])

# Perform PCA on the reshaped data
pca_result <- prcomp(vnir_2d, center = TRUE, scale. = TRUE)

# Extract eigenvalues
eigenvalues <- pca_result$sdev^2

# Create a scree plot to visualize eigenvalues and determine the number of components
plot(eigenvalues, type = "b", xlab = "Component", ylab = "Eigenvalue", main = "Scree Plot")
abline(h = 1, col = "red", lty = 2)  # Kaiser's rule threshold

# Uses Kaiser-Guttman rule which suggest that "the elbow" in scree plots occurs
# when the eigenvalues dip below 1. 
# https://datasciencewiki.net/kaisers-rule/
num_components_kaiser <- sum(eigenvalues > 1)
cat("Number of components determined by Kaiser's rule:", num_components_kaiser, "\n")

# Use the number of components determined by PCA or manually adjust based on domain knowledge
num_components <- num_components_kaiser  

# Initial spectra estimate using PCA
initial_spectra <- pca_result$rotation[, 1:num_components]

# Create initial concentration profile estimates
initial_concentration_profiles <- matrix(
  runif(nrow(vnir_2d) * num_components), 
  nrow = nrow(vnir_2d), 
  ncol = num_components
)


# Define convergence criteria
tolerance <- 1 # Initial tolerance was 1e-6
converged <- FALSE
iteration <- 0
max_iteration<-50
# Initialize previous concentration profiles for convergence checking
previous_concentration_profiles <- initial_concentration_profiles

# While loop for convergence
while ((!converged)|(iteration<max_iteration)) {
  # Run MCR-ALS
  # https://cran.r-project.org/web/packages/ALS/ALS.pdf
  mcr_result <- als(
    CList = list(previous_concentration_profiles),  # Use the previous iteration's concentration profiles
    PsiList = list(vnir_2d),                        # The actual data matrix
    S = initial_spectra,                            # Initial spectra estimate
    nonnegS = TRUE,                                 # Enforce non-negativity for spectra
    nonnegC = TRUE,                                 # Enforce non-negativity for concentration profiles
    maxiter = 10,                                   # Maximum number of iterations per ALS call
    thresh = tolerance                              # Convergence threshold
  )
  
  # Extract resolved spectra and updated concentration profiles
  current_spectra <- mcr_result$S
  current_concentration_profiles <- mcr_result$CList[[1]]
  
  # Check for convergence by comparing concentration profiles
  concentration_diff <- sqrt(sum((current_concentration_profiles - previous_concentration_profiles)^2))
  
  # Check convergence criteria
  if (concentration_diff < tolerance) {
    converged <- TRUE
  } else {
    previous_concentration_profiles <- current_concentration_profiles
    iteration<-iteration+1
  }
}

if (converged) {
  cat("Converged in", iteration, "iterations.\n")
} else {
  cat("Did not converge within the maximum number of iterations.\n")
}
final_conc<-mcr_result$C
final_pure_conc<-mcr_result$S

# Assuming concentrations and final_concentrations are matrices
# with dimensions [n_samples, n_components]

n_components <- ncol(initial_concentration_profiles)  # Number of components
n_samples <- nrow(initial_concentration_profiles)     # Number of samples

# Assuming concentrations is the matrix of true concentrations with dimensions [n_samples, n_components]

n_components <- length(initial_concentration_profiles)  # Number of components
n_samples <- nrow(initial_concentration_profiles)   # Number of samples

# Plot for each component
par(mfrow = c(n_components, 1))  # Set up a plotting area with multiple plots

for (component_index in 1:n_components) {
  # Extract the final concentration matrix for this component
  final_component_conc <- final_conc[[component_index]]
  
  # Plot the true concentrations
  plot(1:n_samples, initial_concentration_profiles[, component_index], type = "l", col = "red",
       lty = 2, lwd = 2, xlab = "Sample Index", ylab = "Concentration",
       main = paste("Component", component_index, "Concentration Profile"))
  
  # Add the modeled concentrations to the plot
  lines(1:n_samples, final_component_conc, col = "blue", lty = 1, lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("True Concentration", "Modeled Concentration"),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2)
}

