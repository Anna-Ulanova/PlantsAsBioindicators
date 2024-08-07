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
#-----------------------------------------Helper functions used inside other functions-----------------------------------------#
'''

'''
def convert_2D(subsetted_hypercube):

    # Assume `hyperspectral_data` is a 3D numpy array of shape (height, width, num_bands)
    # Reshape the data to 2D (num_pixels, num_bands)
    num_pixels = subsetted_hypercube.shape[0]*subsetted_hypercube.shape[1]
    reshaped_data = subsetted_hypercube.reshape((num_pixels, subsetted_hypercube.shape[2]))
    print('func convert_2D -- success')
    return reshaped_data 

def crop_raster(raster):
    height, width, num_bands = raster.shape
    percent = 0
    top_h = round(height*(1-percent))
    bottom_h = round(height*percent)
    # print(top_h)
    top_w = round(width*(1-percent))
    bottom_w = round(width*percent)
    spatial_subset = raster[bottom_h:top_h, bottom_w:top_w, 0:num_bands]  
    return spatial_subset 
#-------------------------------------------------------------------------------------------------------------------------------#
'''
1. Opens the hyperspectral data, retrieves wavelengths associated with each band, crops the raster files by a 20% margin. 
- This function returns a hypercube array, raw raster file, and an array of wavelength per band
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
2. Spectral subset of VNIR files for atmospheric correction
- Adjusting for the effects of the atmosphere on the captured data, ensuring the spectral information accurately represents the ground. 
- Returns the subsetting hypercube and subsetted array of the wavelegnth/band
'''
def spectral_subset(median_blurred, wavelengths): 
    # Find bands with wavelengths less than or equal to max_wavelength
    subset_indices = np.where(wavelengths <=920.00)[0]
    # Subset the hypercube
    subset_hypercube = median_blurred[:, :, subset_indices]
    # Subset the wavelengths
    wavelengths = wavelengths[subset_indices]
    print('func subset_wavelengths -- success')
    return [subset_hypercube, wavelengths]


'''
3. Smoothing- Median filtering to remove the "salt and pepper" noise. 
- Draws a 3 x 3 matrix around each element and determines the median value 
- Output matrix has the same dimensions as the input matrix, returns a hypercube
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
    list = np.asarray(wavelengths)
    idx = (np.abs(list - wavelength)).argmin()
    print('func find_indices_for_wavelength -- success')
    return idx

''' 
5A. Calculate NDVI using the red and NIR band indices
- Uses conventional NDVI formula, sends raw NDVI for background filtering and Gaussian blurring: 
- Unsave regions to view raw and processed NDVI for data
- NDVI highlights:
    - Chlorophyll sensitive 
    - Quantifies plant density and health 

'''
def calculate_NDVI(red_reflectance, NIR_reflectance):
    # Compute NDVI
    NDVI = (NIR_reflectance - red_reflectance) / (NIR_reflectance + red_reflectance)

    # Saves indices associated with each pixel of the raw NDVI matrix 
    np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\ndvi_raw.txt', NDVI)
    
    # Sends for background processing where the NDVI values that are 0.2 or below are converted to zero, and then Gaussian blur is applied
    NDVI_minus_background = post_processing(NDVI)

    # Saves indices associated with each pixel post processing, meaning that background noise will be predominantly zero
    np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\ndvi_no_background.txt', NDVI_minus_background)

    # Display raw NDVI image
    plt.imshow(NDVI, cmap='RdYlGn')
    plt.colorbar()
    plt.title('NDVI Image Raw')
    plt.show()

    # Display NDVI image without background
    plt.imshow(NDVI_minus_background, cmap='RdYlGn')
    plt.colorbar()
    plt.title('NDVI Image Remove Background')
    plt.show()
    print('func calculate_NDVI -- success')
    soil_element_isolation(NDVI_minus_background)
    return NDVI_minus_background
''' 
5B Calculate EVI using the red, blue, and NIR band indices
- Uses constants outlined in: https://doi.org/10.3390/s7112636
- EVI highlights:
    - Responsive to canopy structural variations leaf area index
    - Canopy type and architecture
    - Plant physiognomy
    - Uses blue band to correct for atmosphere
    - Reduces effects of soil background
'''
def calculate_EVI(red, nir, blue):
    G = 2.5
    C1 = 6
    C2 = 7.5
    L = 1
    evi = G * (nir - red) / (nir + C1 * red - C2 * blue + L)
    evi_post_processed = post_processing(evi)
    # Uncomment to see EVI output
    plt.imshow(evi_post_processed, cmap='RdYlGn')
    plt.colorbar()
    plt.title('EVI Image')
    plt.show()
    return evi_post_processed

'''
6. Post-processing
- Masking non-vegetated Areas by assigning any index values above 0.2 as vegetative areas
- Gaussian filter: spatial filtering 
'''
def post_processing(array):
    # Change this value if analyzing overall grass/healthy grass. 
    vegetated = array > 0.10
    masked_vegetated = array * vegetated
    gaussian = gaussian_filter(masked_vegetated, sigma=1)
    return masked_vegetated


'''
6a. NDVI Separation of dead grass vs soil vs live grass
soil = 0.2-0.35
dead grass = 0.35-0.5
live grass = 0.5-1
'''
def soil_element_isolation(array):
    print(array.shape)
    soil1 = 0.2<array
    soil2 = soil1<0.35
    masked_soil = array*soil2
    print(masked_soil.shape)

    plt.imshow(masked_soil, cmap='RdYlGn')
    plt.colorbar()
    plt.title('Soil only?')
    plt.show()

    dead_grass1 = array>0.35  
    dead_grass2 = dead_grass1<0.5
    masked_dead_grass = dead_grass2*array

    plt.imshow(masked_dead_grass, cmap='RdYlGn')
    plt.colorbar()
    plt.title('Dead grass only?')
    plt.show()

    alive_grass1 = array>0.5 
    alive_grass2 = alive_grass1<1
    masked_alive_grass = alive_grass2*array

    plt.imshow(masked_alive_grass, cmap='RdYlGn')
    plt.colorbar()
    plt.title('Alive grass?')
    plt.show()
    print('Isolated soil core elements--work in progress')
'''
7. Bootstrapping output matrix
- Resamples a 10 by 10 matrix that does not contain zeros 1000 times
- Appends all resampled matrices together, finds the total sum of the matrix: 10000 by 10 added up to one number, and then divided by 100000
- Returns the average index value
'''
def bootstrapping(indices): 
    # Number of resamples completed
    num_resamples = 1000
    # Assigns dimensions of the sample matrix, specifies the length of one matrix
    square_side = 10 
    # Finds the number of rows and columns of the input index matrix.
    row, col = indices.shape
    # Cuts off the maximum row number by 10 for picking random points within the index matrix. 
    row_lim = row- 10 
    col_lim = col - 10
    # Loop iteration counter
    i= 0
    # Creates an empty zeros 1 by 10 matrix
    overall_mat = np.zeros((1, 10))
    # Creates a random row number from 10, to row maximux number and a random column number from the same range, creating a random point within the index matrix
    # Uses the random point in matrix to create the resampling sub-matrix. Treats the random point as the bottom left corner, see below:
    #    ----------
    #   |          |
    #   |          |
    #   |          |
    #   |          |
    #  (*)---------
    while i < num_resamples:
        # Finds random row number
        random_row_point = np.random.randint(10, row_lim)
        # Finds random column number
        random_col_point = np.random.randint(10, col_lim)
        # Creates sampling sub-matrix
        sub_array = indices[random_row_point:(random_row_point + square_side), random_col_point:(random_col_point + square_side)]
        # Makes sure that none of the elements inside the matrix are zero
        if 0 not in sub_array:
            i=i+1
            # If matrix does not contain zeroes, the counter increases by one and the sampled matrix is vertically appended to 
            # the overall matrix
            overall_mat= np.concatenate([overall_mat,sub_array], axis=0)
    print('sampled bootstrapped array: ', overall_mat.shape)
    # Calculates the sum of every sampled cell and divides by the number of sampled cells (1000*100)
    NDVI_average = np.concatenate(overall_mat).sum()/(square_side*square_side*num_resamples)
    # Returns a single number
    return NDVI_average
'''
8. Classifies grass/dirt based on threshold input
- Classifies live/healthy grass based on threshold 0.6, reference: https://www.nature.com/articles/s41597-023-02255-3#Sec1 
- Classification shows exactly how much of the soil sample is getting classified as total/healthy grass. 
'''
def classify_dirt_grass(NDVI): 
    # Threshold NDVI to classify grass and dirt
    ndvi_threshold = 0.25 # Example threshold; values above 0.2 are considered grass
    grass_mask = NDVI > ndvi_threshold
    #Visualize the classification result
    plt.imshow(grass_mask, cmap='gray')
    plt.title('Grass Classification')
    plt.show()
    print('func classify_dirt_grass -- success')
    return grass_mask 

def main_launch(target):

    print("---------Processing---------")
    metadata = 'raw_rd_rf.hdr'
    hyperspectral = 'raw_rd_rf'

    # Opens hyperspectral data, crops 20% from all four sides, finds the wavelength list associated with each band
    [vnir_data, vnir, wavelengths] = open_hyperspectral_data(target, metadata, hyperspectral)

    # Uncomment to display the hyperspectral image
    # plt.imshow(vnir_data[:, :, 100], cmap='gray')
    # plt.title('Original')
    # plt.colorbar(label='Intensity')
    # plt.show()

    # Completes median filtering using a 3 x 3 matrix
    blurred = smoothing_and_blurring(vnir_data)

    # Uncomment to display the hyperspectral image
    # plt.imshow(blurred[:,:,100], cmap='gray')
    # plt.title('Median Filtering')
    # plt.colorbar(label='Intensity')
    # plt.show()

    # Completes spectral subset which removes any bands that have a wavelength greater than 920 nm associated with it 
    [subsetted_hypercube, wavelengths] = spectral_subset(blurred, wavelengths)

    # Finds average wavelengths for each band type
    red_wavelength = (625+740)/2
    NIR_wavelength = (780+1000)/2
    blue_wavelength = (450+495)/2


    # Finds the index for red wavelength. Uses helper functions to determine the nearest wavelength for the band. 
    red_band = find_indices_for_wavelength(wavelengths, red_wavelength)
    print('Red band number: ', red_band)
    # Extract the red and NIR bands
    red = subsetted_hypercube[:, :, red_band]
    ## Displays red band
    # plt.imshow(red, cmap='gray')
    # plt.title('Red Band')
    # plt.colorbar(label='Intensity')
    # plt.show()


#######DELETE
    plt.hist(red, bins='auto')
    plt.title('Spread of Red Reflectance Intensities')
    plt.show()
##############
    NIR_band = find_indices_for_wavelength(wavelengths, NIR_wavelength)
    print('NIR band number: ', NIR_band)
    nir = subsetted_hypercube[:, :, NIR_band]
    ## Displays red band
    # plt.imshow(nir, cmap='gray')
    # plt.title('NIR Band')
    # plt.colorbar(label='Intensity')
    # plt.show()

#######DELETE
    plt.hist(nir, bins='auto')
    plt.title('Spread of NIR Reflectance Intensities')
    plt.show()
##############
    blue_band = find_indices_for_wavelength(wavelengths, blue_wavelength)
    print('Blue band number: ', blue_band)
    blue = subsetted_hypercube[:, :, blue_band]
    # Displays red band
    # plt.imshow(blue, cmap='gray')
    # plt.title('Blue Band')
    # plt.colorbar(label='Intensity')
    # plt.show()

    np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\red_band.txt', red)
    np.savetxt(r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\nir_band.txt', nir)
    NDVI = calculate_NDVI(red, nir)
    EVI = calculate_EVI(red, nir, blue)

    NDVI_averaged = bootstrapping(NDVI)
    EVI_averaged = bootstrapping(EVI)
    # CHANGE DIRECTORY TO LOCAL MACHINE
    results_doc_dir = r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\test_20240807.txt'
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

#####################################################################################################################################

# Change to the location of the folder that has the VNIR files
print('Script Initiated...')
# CHANGE DIRECTORY TO LOCAL MACHINE
# Directory where VNIR files would be located.
# - Testing 8/7/2024: copied items:   D:\VNIR\GH_20230724\16_3_20230724_2023_07_24_13_10_38
                                    # D:\VNIR\GH_20230726\9_1_20230726_2023_07_26_10_33_15
                                    # D:\VNIR\GH_20230801\11_3_20230801_2023_08_01_13_07_17
                                    # D:\VNIR\GH_20230822\7_2_20230822_2023_08_22_06_49_43
                                    # D:\VNIR\GH_20231114\5_2_20231114_2023_11_14_07_51_18
directory = r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\VNIR'

vnir_files = []
summary_file = [('File_Name', 'NDVI Summary', 'EVI Summary')]
# CHANGE DIRECTORY TO LOCAL MACHINE
filename = r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\test_20240807.txt'

# Creates a summary document which saves the average NDVI/EVI for hyperspectral data. 
with open(filename, 'w') as file:
    for row in summary_file:
        # Join the elements of the tuple with tabs and write to the file
        file.write('\t'.join(map(str, row)) + '\n')
print('Summary document for NDVI averages created')
# Finds the directories of all VNIR files, returns only the folder where the raw_rd_rf files are saved
for (root, dirs, files) in os.walk(directory, topdown=True):
    for dir in dirs:
        file_path = os.path.join(root, dir)
        if "VNIR" in file_path:
            print(file_path)
            vnir_files.append(file_path)
print('VNIR directories recovered: ', '\n', vnir_files)

main_launch(r'C:\Users\RDCRLAAU\Desktop\Plant_as_bioindicators\PlantsAsBioindicators-Python\VNIR\16_3_20230724_2023_07_24_13_10_38')
# # Launches main function. The main function calls on other subfunctions that complete specific analyses. 
# for vnir in vnir_files: 
#     main_launch(vnir)