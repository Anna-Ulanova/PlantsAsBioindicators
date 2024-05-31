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
#-------------------------------------------------------------------------------------------------------------------------------#
'''
Helper functions used inside other functions 
'''
def closest(list, k): 
     list = np.asarray(list)
     idx = (np.abs(list - k)).argmin()
     print(list[idx])
     print('func closest -- success')
     return idx

def find_wavelength(input):
    wavelengths = np.array(input.metadata['wavelength'])
    wavelengths = wavelengths.astype(float)
    return wavelengths
#-------------------------------------------------------------------------------------------------------------------------------#
'''
Opens the hyperspectral data
'''
def open_hyperspectral_data(target, metadata, hyperspectral): 
    vnir = envi.open(os.path.join(target, metadata), os.path.join(target, hyperspectral))
    vnir_array = vnir.load()
    print('func open_hyperspectral_data -- success')
    return vnir_array
#-------------------------------------------------------------------------------------------------------------------------------#

'''
1. Spectral subset of VNIR files for atmospheric correction
- Adjusting for the effects of the atmosphere on the captured data, ensuring the spectral information accurately represents the ground.
'''
def subset_wavelengths(vnir): 
    # Find indices of subsetted wavelengths
    # Extract wavelengths
    wavelengths = find_wavelength(vnir)
    # Find bands with wavelengths less than or equal to max_wavelength
    subset_indices = np.where(wavelengths <=920.00)[0]
    # Subset the hypercube
    subset_hypercube = vnir.read_bands(subset_indices)
    print('func subset_wavelengths -- success')
    return [subset_hypercube, wavelengths]

def convert_2D(subsetted_hypercube):
    
    # Assume `hyperspectral_data` is a 3D numpy array of shape (height, width, num_bands)
    # Reshape the data to 2D (num_pixels, num_bands)
    num_pixels = subsetted_hypercube.shape[0]*subsetted_hypercube.shape[1]
    reshaped_data = subsetted_hypercube.reshape((num_pixels, subsetted_hypercube.shape[2]))
    print('func convert_2D -- success')
    return reshaped_data 

'''
2. Feature Extraction
- Given the large number of bands, it's often helpful to reduce the data's dimensionality while preserving important information.
- Principal Component Analysis (PCA): A technique that transforms the data into a set of orthogonal components, reducing redundancy.
PCA also reduces the noise my data. The output increases the contrast between the grass/soil/tray, resulting in the soil core to be 
darker and the tray to be bright. 
'''
def perform_PCA(reshaped_data, subset_hypercube, wavelengths):
    # Apply PCA
    pca = PCA(n_components=150)  # Reduce to 10 principal components
    pca_data = pca.fit_transform(reshaped_data)
    # Reshape back to 3D (height,width, number of components)
    pca_data_3D = pca_data.reshape((subset_hypercube.shape[0], subset_hypercube.shape[1], 150))
    # Associate Wavelengths with Bands in PCA-transformed Data
    # Create an array to hold wavelengths associated with each band in PCA-transformed data
    pca_wavelengths = np.empty(pca_data.shape[1])
    # Assign wavelengths to each principal component based on their original bands
    # Example: Assign first 10 wavelengths to the first principal component, next 10 wavelengths to the second principal component, and so on
    for i in range(pca_data.shape[1]):
        pca_wavelengths[i] = wavelengths[i]
    # Visualize the first principal component
    # Increases brightness of the soil sample
    plt.imshow(pca_data_3D[:, :, 0], cmap='gray')
    plt.title('First Principal Component')
    plt.show()
    print('func perform_PCA -- success')
    return [pca_data_3D, pca_wavelengths]

'''
4. Finds the indices of bands that have wavelengths corresponding to red/NIR
- Indices for red/NIR are used to determine the NDVI
- The color red occurs between 625 nm and 740 nm, the average value of the two wavelengths was used. This avoids edge cases and is more representative of how ENVI displays the red in RBG view. 
- NIR occurs between 780 nm and 2500 nm
- Uses helper function closest to find the wavelength value closest to the average of the range
'''
def find_indices_for_wavelength(wavelengths, wavelength):
    """
    Find the indices of "red" bands in hyperspectral data.

    Args:
        wavelengths (ndarray): Array of wavelengths.
        wavelength_range (tuple): Tuple containing the lower and upper bounds of the "red" wavelength range.

    Returns:
        list: List of indices corresponding to "red" bands.
    """
    close_in_list = closest(wavelengths, wavelength)
    # Use list comprehension to return index where the wavelength is found in list
    print('func find_indices_for_wavelength -- success')
    return close_in_list


''' 
5. Calculate NDIV using the red and NIR band indices
- Uses the first red/NIR band 
- Define the indices for the red and NIR bands 
'''
def calculate_NDVI(red_reflectance, NIR_reflectance):
    # Compute NDVI
    NDVI = (NIR_reflectance - red_reflectance) / (NIR_reflectance + red_reflectance)
    print(NDVI)
    # Display NDVI image
    plt.imshow(NDVI, cmap='RdYlGn')
    plt.colorbar()
    plt.title('NDVI Image')
    plt.show()
    print('func calculate_NDVI -- success')
    return NDVI

'''
6. Classifies grass/dirt based on threshold input
- Classifies live/healthy grass based on threshold 0.6, reference: https://www.nature.com/articles/s41597-023-02255-3#Sec1 
'''
def classify_dirt_grass(NDVI): 
    # Threshold NDVI to classify grass and dirt
    ndvi_threshold = 0.6  # Example threshold; values above 0.2 are considered grass
    grass_mask = NDVI > ndvi_threshold
    # Visualize the classification result
    plt.imshow(grass_mask, cmap='gray')
    plt.title('Grass Classification')
    plt.show()
    print('func classify_dirt_grass -- success')
    return grass_mask 
target = r'C:\Users\RDCRLAAU\Desktop\Backup\overlapped_SWIR_VNIR\VNIR\GH_20231213\1_1_20231213_2023_12_13_07_44_08'
metadata = 'raw_rd_rf.hdr'
hyperspectral = 'raw_rd_rf'

os.chdir(target)

vnir_data = open_hyperspectral_data(target, metadata, hyperspectral)

[subsetted_hypercube, wavelength] = subset_wavelengths(vnir_data)

reshaped_data = convert_2D(subsetted_hypercube)

[pca_data_3D, new_wavelength] = perform_PCA(reshaped_data, subsetted_hypercube, wavelength)


red_wavelength_range = (625+740)/2
NIR_wavelength_range = 920

red_band = find_indices_for_wavelength(new_wavelength, red_wavelength_range)
print(red_band)
NIR_band = find_indices_for_wavelength(new_wavelength, NIR_wavelength_range)

# Extract the red and NIR bands
print(pca_data_3D.shape)
red = pca_data_3D[:, :, red_band]
print(red)
nir = pca_data_3D[:, :, NIR_band]
print(nir)
NDVI = calculate_NDVI(red, nir)

grass_classified = classify_dirt_grass(NDVI)