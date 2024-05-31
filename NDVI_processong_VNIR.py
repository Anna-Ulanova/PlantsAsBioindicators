import numpy as np 
from sklearn.decomposition import PCA
from spectral import imshow, view_cube
import spectral.io.envi as envi
import matplotlib.pyplot as plt
import matplotlib
import os
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from scipy.linalg import eigh
import pywt
import csv

def open_hyperspectral_data(target, metadata, hyperspectral): 
    vnir = envi.open(os.path.join(target, metadata), os.path.join(target, hyperspectral))
    vnir_array = vnir.load()
    return vnir


'''
1. Spectral subset of VNIR files for atmospheric correction
- Adjusting for the effects of the atmosphere on the captured data, ensuring the spectral information accurately represents the ground.
'''

def subset_wavelengths(vnir): 
    # Find indices of subsetted wavelengths
# Extract wavelengths
    wavelengths = np.array(vnir.metadata['wavelength'])

    # Find bands with wavelengths less than or equal to max_wavelength
    subset_indices = np.where(wavelengths <= 920)[0]

    # Subset the hypercube
    subset_hypercube = vnir.read_bands(subset_indices)

    return [subset_hypercube, wavelengths]

def convert_2D(subsetted_hypercube):
    
    # Assume `hyperspectral_data` is a 3D numpy array of shape (height, width, num_bands)
    # Reshape the data to 2D (num_pixels, num_bands)
    num_pixels = subsetted_hypercube.shape[0]*subsetted_hypercube.shape[1]
    reshaped_data = subsetted_hypercube.reshape((num_pixels, subsetted_hypercube.shape[2]))
    return reshaped_data 

'''

2. Feature Extraction

- Given the large number of bands, it's often helpful to reduce the data's dimensionality while preserving important information.

- Principal Component Analysis (PCA): A technique that transforms the data into a set of orthogonal components, reducing redundancy.
PCA also reduces the noise my data. The output increases the contrast between the grass/soil/tray, resulting in the soil core to be 
darker and the tray to be bright. 
'''
def perform_PCA(reshaped_data, subset_hypercube):
    # Apply PCA
    pca = PCA(n_components=10)  # Reduce to 10 principal components
    pca_data = pca.fit_transform(reshaped_data)

    # Reshape back to 3D (height,width, number of components)
    tpca_data_3D = pca_data.reshape((subset_hypercube.shape[0], subset_hypercube.shape[1], 10))
    # Visualize the first principal component
    # Increases brightness of the soil sample
    plt.imshow(pca_data_3D[:, :, 0], cmap='gray')
    plt.title('First Principal Component')
    plt.show()
    return pca_data_3D


'''
4. Finds the indices of bands that have wavelengths corresponding to red/NIR
- Indices for red/NIR are used to determine the NDVI
- The color red occurs between 625 nm and 740 nm 
- NIR occurs between 780 nm and 2500 nm
'''
def find_indices_for_wavelength(wavelengths, wavelength_range):
    """
    Find the indices of "red" bands in hyperspectral data.

    Args:
        wavelengths (ndarray): Array of wavelengths.
        wavelength_range (tuple): Tuple containing the lower and upper bounds of the "red" wavelength range.

    Returns:
        list: List of indices corresponding to "red" bands.
    """
    lower_bound, upper_bound = wavelength_range
    band_indices = np.where((wavelengths >= lower_bound) & (wavelengths <= upper_bound))[0]
    return band_indices
''' 
5. Calculate NDIV using the red and NIR band indices
- Uses the first red/NIR band 
- Define the indices for the red and NIR bands 
'''
def calculate_NDVI(red_reflectance, NIR_reflectance, hyperspectral_data):

    # Compute NDVI
    NDVI = (NIR_reflectance - red_reflectance) / (NIR_reflectance + red_reflectance)
    
    # Display NDVI image
    plt.imshow(NDVI, cmap='RdYlGn')
    plt.colorbar()
    plt.title('NDVI Image')
    plt.show()
    return NDVI

'''
6. Classifies grass/dirt based on threshold input
- Classifies live/healthy grass based on threshold 0.6, reference: https://www.nature.com/articles/s41597-023-02255-3#Sec1 
'''

def classify_dirt_grass(NDVI): 

    return grass_mask 

# Example usage
hyperspectral_data = np.random.rand(100, 100, 200)  # Example 3D hyperspectral array
wavelengths = np.linspace(400, 1000, 200)  # Example array of wavelengths
red_wavelength_range = (600, 700)  # Example "red" wavelength range

red_band_indices = find_red_band_indices(wavelengths, red_wavelength_range)
print("Indices of red bands:", red_band_indices)

target = r'C:\Users\RDCRLAAU\Desktop\Backup\overlapped_SWIR_VNIR\VNIR\GH_20231213\1_1_20231213_2023_12_13_07_44_08'
metadata = 'raw_rd_rf.hdr'
hyperspectral = 'raw_rd_rf'
os.chdir(target)

# Example usage
metadata_file = 'metadata.csv'  # Path to your metadata file
vnir_data = open_hyperspectral_data(target, metadata, hyperspectral)

[subsetted_hypercube, wavelengths] = subset_wavelengths(vnir_data)

reshaped_data = convert_2D(subsetted_hypercube)

pca_data_3D = perform_PCA(reshaped_data, subsetted_hypercube)

red_wavelength_range = [625, 740]

NIR_wavelength_range = [780, 2500]

red_indices = find_indices_for_wavelength(wavelengths, red_wavelength_range)

NIR_indices = find_indices_for_wavelength(wavelengths, NIR_wavelength_range)

# Uses the first red/NIR band 
red_band = red_indices[0]
NIR_band = NIR_indices[0]

# Extract the red and NIR bands
red = hyperspectral_data[:, :, red_band]
nir = hyperspectral_data[:, :, NIR_band]

NDVI = calculate_NDVI(red_indices, NIR_indices, pca_data_3D)
# # Apply ICA
# ica = FastICA(n_components=10)  # Reduce to 10 independent components
# ica_data = ica.fit_transform(reshaped_data)

# # Reshape back to 3D
# ica_data_3d = ica_data.reshape((vnir_array.shape[0], vnir_array.shape[1], 10))

# # Visualize the first independent component
# plt.imshow(ica_data_3d[:, :, 0], cmap='gray')
# plt.title('First Independent Component')
# plt.show()

# # Band selection 
# # How do we choose the selected bands? 
# selected_bands = [10, 20, 30]
# selected_data = vnir_array[:, :, selected_bands]

# # Visualize the first selected band
# plt.imshow(selected_data[:, :, 0], cmap='gray')
# plt.title('Selected Band 10')
# plt.show()

# # Standardize the data
# scaler = StandardScaler()
# scaled_data = scaler.fit_transform(reshaped_data)

# # Compute the covariance matrix
# cov_matrix = np.cov(scaled_data, rowvar=False)


# # Compute the eigenvalues and eigenvectors
# eigenvalues, eigenvectors = eigh(cov_matrix)

# # Select the components with the highest eigenvalues
# num_components = 10
# selected_indices = np.argsort(eigenvalues)[-num_components:]
# mnf_data = np.dot(scaled_data, eigenvectors[:, selected_indices])

# # Reshape back to 3D
# mnf_data_3d = mnf_data.reshape((vnir_array.shape[0], vnir_array.shape[1], num_components))

# # Visualize the first minimum noise fraction transformation component
# # MNF transformation reduces noise and extracts features by decorrelating and rescaling the
# plt.imshow(mnf_data_3d[:, :, 0], cmap='gray')
# plt.title('First MNF Component')
# plt.show()
