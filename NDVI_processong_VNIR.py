import numpy as np 
from sklearn.decomposition import PCA
from spectral import imshow, view_cube
import spectral.io.envi as envi
import matplotlib.pyplot as plt
import matplotlib
import os
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from scipy.ndimage import gaussian_filter
from scipy.linalg import eigh
import pywt
import csv
import cv2
#-------------------------------------------------------------------------------------------------------------------------------#
'''
Helper functions used inside other functions 
'''
def convert_2D(subsetted_hypercube):

    # Assume `hyperspectral_data` is a 3D numpy array of shape (height, width, num_bands)
    # Reshape the data to 2D (num_pixels, num_bands)
    num_pixels = subsetted_hypercube.shape[0]*subsetted_hypercube.shape[1]
    reshaped_data = subsetted_hypercube.reshape((num_pixels, subsetted_hypercube.shape[2]))
    print('func convert_2D -- success')
    return reshaped_data 

def closest(list, k): 
     list = np.asarray(list)
     idx = (np.abs(list - k)).argmin()
    #  print(list[idx])
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
    # print(top_h)
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
    # print(wavelengths)
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
    #NDVI = (NIR_reflectance - red_reflectance) / (NIR_reflectance + red_reflectance)
    #NDVI = (1- red_reflectance/NIR_reflectance)*(1+red_reflectance/NIR_reflectance)
    print(red_reflectance.shape)
    NDVI = (NIR_reflectance - red_reflectance) / (NIR_reflectance + red_reflectance)
    # np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\PlantsAsBioindicators-Python\ndvi_raw.txt', NDVI)
    print('NDVI shape is: ', NDVI.shape)
    NDVI_minus_background = post_processing(NDVI)
    # np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\PlantsAsBioindicators-Python\ndvi_no_background.txt', NDVI_minus_background)
    print(NDVI_minus_background)
    # Display raw NDVI image
    # plt.imshow(NDVI, cmap='RdYlGn')
    # plt.colorbar()
    # plt.title('NDVI Image Raw')
    # plt.show()
    # Display NDVI image without background
    plt.imshow(NDVI_minus_background, cmap='RdYlGn')
    plt.colorbar()
    plt.title('NDVI Image Remove Background')
    plt.show()
    print('func calculate_NDVI -- success')
    return NDVI_minus_background

def calculate_EVI(red, nir, blue):
    G = 2.5
    C1 = 6
    C2 = 7.5
    L = 1
    evi = G * (nir - red) / (nir + C1 * red - C2 * blue + L)
    # print(evi)
    evi_post_processed = post_processing(evi)
    plt.imshow(evi_post_processed, cmap='RdYlGn')
    plt.colorbar()
    plt.title('EVI Image')
    plt.show()
    return evi_post_processed

'''
Post-processing
-- Masking non-vegetated Areas
-- Gaussian filter: spatial filtering 
'''
def post_processing(array):
    vegetated = array > 0.2 
    masked_vegetated = array * vegetated
    gaussian = gaussian_filter(masked_vegetated, sigma=1)
    return gaussian

'''
Bootstrapping output matrix
-- Resamples a 10 by 10 matrix that does not contain zeros and 
'''
def bootstrapping(indices): 
    num_resamples = 1000
    square_side = 10 
    row, col = indices.shape
    row_lim = row-10 
    col_lim = col -10
    i= 0
    overall_mat = np.zeros((1, 10))
    while i  <num_resamples:
        random_row_point = np.random.randint(10, row_lim)
        random_col_point = np.random.randint(10, col_lim)
        sub_array = indices[random_row_point:(random_row_point + square_side), random_col_point:(random_col_point + square_side)]
        if 0 not in sub_array:
            i=i+1
            overall_mat= np.concatenate([overall_mat,sub_array], axis=0)
    print('sampled bootstrapped array: ', overall_mat.shape)
    NDVI_average = np.concatenate(overall_mat).sum()/(square_side*square_side*num_resamples)
    return NDVI_average
'''
6. Classifies grass/dirt based on threshold input
- Classifies live/healthy grass based on threshold 0.6, reference: https://www.nature.com/articles/s41597-023-02255-3#Sec1 
'''
def classify_dirt_grass(NDVI): 
    # Threshold NDVI to classify grass and dirt
    ndvi_threshold = 0.25 # Example threshold; values above 0.2 are considered grass
    grass_mask = NDVI > ndvi_threshold
    # Visualize the classification result
    # plt.imshow(grass_mask, cmap='gray')
    # plt.title('Grass Classification')
    # plt.show()
    print('func classify_dirt_grass -- success')
    return grass_mask 

print("Processing")

# Dead grass
target = 'C:\\Users\\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\\GH_20231213_1_1_20231213_2023_12_13_07_44_08'
# target = 'C:\\Users\\RDCRLAAU\\Desktop\\Plant as bioindicators\\VNIR\\GH_20230724_12_1_20230724_2023_07_24_12_07_56'
# target = 'C:\\Users\\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\GH_20230726_3_2_20230726_2023_07_26_09_41_37'
# target = 'C:\\Users\\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\GH_20230822_9_2_20230822_2023_08_22_07_09_53'
# target = 'C:\\Users\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\GH_20231116_6_1_20231116_2023_11_16_07_07_28'
# target = 'C:\\Users\\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\GH_20231116_15_2_20231116_2023_11_16_10_50_03'
# target = 'C:\Users\\RDCRLAAU\Desktop\\Plant as bioindicators\\VNIR\\GH_20231121_3_3_20231121_2023_11_21_08_52_18'
# Live grass 
#target = r'D:\VNIR\GH_20231116\15_2_20231116_2023_11_16_10_50_03'
metadata = 'raw_rd_rf.hdr'
hyperspectral = 'raw_rd_rf'

os.chdir(target)

[vnir_data, vnir, wavelengths] = open_hyperspectral_data(target, metadata, hyperspectral)
print('vnir_data shape: ', vnir_data.shape)
# plt.imshow(vnir_data[:, :, 100], cmap='gray')
# plt.title('Original')
# plt.colorbar(label='Intensity')
# plt.show()

blurred = smoothing_and_blurring(vnir_data)
print('median blurring: ', blurred.shape)
# plt.imshow(blurred[:,:,100], cmap='gray')
# plt.title('Median Filtering')
# plt.colorbar(label='Intensity')
# plt.show()


[subsetted_hypercube, wavelengths] = spectral_subset(blurred, wavelengths)


red_wavelength = (625+740)/2
NIR_wavelength = (780+1000)/2
blue_wavelength = (450+495)/2

red_band = find_indices_for_wavelength(wavelengths, red_wavelength)
print('Red band number: ', red_band)
# Extract the red and NIR bands
red = subsetted_hypercube[:, :, red_band]
# Displays red band
# plt.imshow(red, cmap='gray')
# plt.title('Red Band')
# plt.colorbar(label='Intensity')
# plt.show()


NIR_band = find_indices_for_wavelength(wavelengths, NIR_wavelength)
print('NIR band number: ', NIR_band)
nir = subsetted_hypercube[:, :, NIR_band]
# Displays red band
# plt.imshow(nir, cmap='gray')
# plt.title('NIR Band')
# plt.colorbar(label='Intensity')
# plt.show()


blue_band = find_indices_for_wavelength(wavelengths, blue_wavelength)
print('Blue band number: ', blue_band)
blue = subsetted_hypercube[:, :, blue_band]
# Displays red band
# plt.imshow(blue, cmap='gray')
# plt.title('Blue Band')
# plt.colorbar(label='Intensity')
# plt.show()


NDVI = calculate_NDVI(red, nir)
EVI = calculate_EVI(red, nir, blue)
grass_classified = classify_dirt_grass(NDVI)
grass_classified = classify_dirt_grass(EVI)

NDVI_averaged = bootstrapping(NDVI)
EVI_averaged = bootstrapping(EVI)

results_doc_dir = r'C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\PlantsAsBioindicators-Python\summary_NDVI_EVI_averages.txt'
append_new = target +'\t'+str(NDVI_averaged)+'\t'+str(EVI_averaged)
# Open the file in append & read mode ('a+')
with open(results_doc_dir, "a+") as file_object:
    # Move read cursor to the start of file.
    file_object.seek(0)
    # If file is not empty then append '\n'
    data = file_object.read(100)
    if len(data) > 0 :
        file_object.write("\n")
    # Append text at the end of file
    file_object.write(append_new)