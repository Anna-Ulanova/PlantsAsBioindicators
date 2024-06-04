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
import cv2
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

def crop_raster(raster):
    height, width, num_bands = raster.shape
    top_h = round(height*0.8)
    bottom_h = round(height*0.2)
    print(top_h)
    top_w = round(width*0.8)
    bottom_w = round(width*0.2)
    spatial_subset = raster[bottom_h:top_h, bottom_w:top_w, 0:num_bands]  
    return spatial_subset 
#-------------------------------------------------------------------------------------------------------------------------------#
'''
Opens the hyperspectral data
'''
def open_hyperspectral_data(target, metadata, hyperspectral): 
    vnir = envi.open(os.path.join(target, metadata), os.path.join(target, hyperspectral))
    vnir_array = vnir.load()
    vnir_subset = crop_raster(vnir_array)
    wavelength = np.array(vnir.metadata['wavelength'])
    wavelengths_array = wavelength.astype(float)
    print('func open_hyperspectral_data -- success')
    return [vnir_subset, vnir, wavelengths_array]
#-------------------------------------------------------------------------------------------------------------------------------#

'''
1. Spectral subset of VNIR files for atmospheric correction
- Adjusting for the effects of the atmosphere on the captured data, ensuring the spectral information accurately represents the ground.
'''
def spectral_subset(median_blurred, wavelengths): 
    # Find bands with wavelengths less than or equal to max_wavelength
    subset_indices = np.where(wavelengths <=920.00)[0]
    print(subset_indices)
    # Subset the hypercube
    subset_hypercube = median_blurred[:, :, subset_indices]
    wavelengths = wavelengths[subset_indices]
    print(wavelengths)
    print('func subset_wavelengths -- success')
    return [subset_hypercube, wavelengths]


'''
3. Smoothing- Median filtering to remove the "salt and pepper" noise. 
- Draws a 9 x 9 matrix around each element and determines the median value 
- Output matrix has the same dimensions as the input matrix 
- Moved up earlier because the dimensions remain, unsure how to handle data types
'''
def smoothing_and_blurring(vnir_array):
    height, width, num_bands = vnir_array.shape
    blurred_matrices = [cv2.medianBlur(vnir_array[:,:,i],3) for i in range(num_bands)]
    filtered3D = np.dstack(blurred_matrices)
    print('Filtered Array ', filtered3D.shape)
    return filtered3D
    

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
    ndvi_threshold = 0.2  # Example threshold; values above 0.2 are considered grass
    grass_mask = NDVI > ndvi_threshold
    # Visualize the classification result
    plt.imshow(grass_mask, cmap='gray')
    plt.title('Grass Classification')
    plt.show()
    print('func classify_dirt_grass -- success')
    return grass_mask 

print("Processing")
target = r'C:\Users\RDCRLAAU\Desktop\Backup\overlapped_SWIR_VNIR\VNIR\GH_20231213\1_1_20231213_2023_12_13_07_44_08'
metadata = 'raw_rd_rf.hdr'
hyperspectral = 'raw_rd_rf'

os.chdir(target)

[vnir_data, vnir, wavelengths] = open_hyperspectral_data(target, metadata, hyperspectral)

plt.imshow(vnir_data[:, :, 100], cmap='gray')
plt.title('Original')
plt.colorbar(label='Intensity')
plt.show()

blurred = smoothing_and_blurring(vnir_data)

plt.imshow(blurred[:,:,100], cmap='gray')
plt.title('Median Filtering')
plt.colorbar(label='Intensity')
plt.show()


[subsetted_hypercube, wavelengths] = spectral_subset(blurred, wavelengths)


red_wavelength = (625+740)/2
NIR_wavelength = (780+1000)/2

red_band = find_indices_for_wavelength(wavelengths, red_wavelength)
print('Red band number: ', red_band)
# Extract the red and NIR bands
red = subsetted_hypercube[:, :, red_band]
# Displays red band
plt.imshow(red, cmap='gray')
plt.title('Red Band')
plt.colorbar(label='Intensity')
plt.show()


NIR_band = find_indices_for_wavelength(wavelengths, NIR_wavelength)
print('NIR band number: ', NIR_band)
nir = subsetted_hypercube[:, :, NIR_band]
# Displays red band
plt.imshow(nir, cmap='gray')
plt.title('NIR Band')
plt.colorbar(label='Intensity')
plt.show()


NDVI = calculate_NDVI(red, nir)

grass_classified = classify_dirt_grass(NDVI)