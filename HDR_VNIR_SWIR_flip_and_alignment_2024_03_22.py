import os
import cv2
from matplotlib import pyplot as plt
import rasterio
import numpy as np

# AU Note: matplotlib as plt can allow you to give the images as an  output which is helpful 
# helpful during debugging. Uncomment the plt.imshow() lines to see image outputs
 
####################################### List of Functions #######################################
# Function takes in the main directory where the SWIR and VNIR folders are located, an empty and 
# initiated list that will contain the directories of the HDR files, as well as a search key for
# the files you want to isolate (EX:raw_rd_rf.hdr)
# Returns a list of directories where HDR files are found
def directories_of_hdr_files(script_path, directory_list, search_key):
    for root, _ , files in os.walk(script_path, topdown=True):
    # Change search key 
        if search_key in files: 
            hdr_file_directory = os.path.join(root,search_key)
            directory_list.append(hdr_file_directory)
    print('directories_of_hdr_files complete')
    return directory_list
##-------------------------------------------------------------------------------------------###
# Function takes in the directories of hdr files for SWIR and VNIR types, 
# couples hdr directories if their sampling time matches (e.i. identical last 2 folders)
# Returns the two dimensation list that lists the "coupled" HDR files
def matching_vnir_swir_files(swir_directories_with_hdr, vnir_directories_with_hdr):
    # Pairs HDR files if the last 2 folders match, EX: both are in GH_20230724/1_1_20230724_2023_07_24_08_16_25/raw_rd_rf.hdr
    # Creates a list pairs_for_overlap which can be overlapped/processed
    blank =[]
    for directory in swir_directories_with_hdr:
        # Splits lines based on the common delimiter: \ 
        sub_key = directory.split('\\')
        # Joins the last three elements of the split, essentially limiting the search only for files within the last 2 folders
        # This search finds only the hyperspectral SWIR/VNIR files if they were collected for the same plant, during the same time and date
        sub = '\\'.join(sub_key[-2:])
        # Searches for directories that have a partial match for the common folders
        match_found_in_VNIR = [x for x in vnir_directories_with_hdr if sub in x]
        if directory != match_found_in_VNIR[0]:
            row = [directory, match_found_in_VNIR[0]]
            blank.append([row])
    # Turn on if you want to save the coupled list, useful for sanity checks. 
    with open("file.txt", "w") as output:
        output.write(str(blank))
    print('matching_vnir_swir_files complete')
    return blank
##-------------------------------------------------------------------------------------------###
# This function reads the hyperspectral files using rasterio. Function also extracts the first 
# band of the image for downstream processing. 
def opening_hyperspectral_file(file_dir):
    print(file_dir)
    with rasterio.open(file_dir) as src: 
        # Read the first band
        band = src.read(1)  # Assuming bands are indexed starting from 1
        plt.imshow(band)
        plt.title('First Band')
        plt.show()
    print('opening_hyperspectral_file complete')
    return band
##-------------------------------------------------------------------------------------------###
# This function flips the SWIR image horizontally to make a "mirrored" image of itself
# Inputs are the directory of the SWIR image, as well as the new folder that was created 
# for saving the augmented SWIR images
# Returns the pathway of the augmented SWIR image, as well as the "unique" name created for the image
# (unique name is the combination of the two most recent folders which represent the sampling time and date)
def horizontal_flipping(band, flipped_image_dir):
    # Changes directory to the folder where the flipped images are saved
    os.chdir(flipped_image_dir) 
    flipped_band = np.fliplr(band)
    plt.imshow(flipped_band)
    plt.title('Flipped SWIR')
    plt.show()
    # Returns the flipped image directory for overlapping
    print('horizontal_flipping complete')
    return flipped_image_dir, flipped_band
    
##-------------------------------------------------------------------------------------------###
# Function overlaps images based on edges found within the images. First the function detects edges
# within the HDR file, and contours the edges. Then it finds bounding boxes for contours
# Helpful links for Anna: https://datacarpentry.org/image-processing/edge-detection.html
# Explanaition of Canny and edge detection: https://docs.opencv.org/3.4/da/d22/tutorial_py_canny.html
def overlapping_image(flipped_SWIR, VNIR, overlapped_images_dir, unique_image_name):
    os.chdir(overlapped_images_dir)
    # Function to align images based on edge features

    # Detect edges in the images
    edges1 = cv2.Canny(flipped_SWIR, 100, 200)
    edges2 = cv2.Canny(VNIR, 100, 200)

    # Find contours in the edges
    contours1, _ = cv2.findContours(edges1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours2, _ = cv2.findContours(edges2, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Calculate bounding boxes for contours
    x1, y1, w1, h1 = cv2.boundingRect(contours1[0])
    x2, y2, w2, h2 = cv2.boundingRect(contours2[0])

    # Calculate translation offsets
    dx = x2 - x1
    dy = y2 - y1

    # Translate the input image
    rows, cols, _ = im1.shape
    M = np.float32([[1, 0, dx], [0, 1, dy]])
    translated_im = cv2.warpAffine(im1, M, (cols, rows))

    # Resize SWIR image to match VNIR image size
    resized_swir_image = cv2.resize(translated_im, (vnir_image.shape[1], vnir_image.shape[0]))

    # Save the resized SWIR image
    cv2.imwrite(unique_image_name, resized_swir_image)
    print('overlapping_image complete')

##-------------------------------------------------------------------------------------------###

####################################### Input variables #######################################

# Change target directory after testing 
target_directory = "D:\AU"

# Changes to directory to target path and reformats the delimiters to system's operating system
os.chdir(target_directory)
script_path = os.getcwd()

empty_swir_directories_with_hdr = []
empty_vnir_directories_with_hdr = []

swir_directories_with_hdr = directories_of_hdr_files(script_path, empty_swir_directories_with_hdr, search_key='raw_rd_rf')

vnir_directories_with_hdr = directories_of_hdr_files(script_path, empty_vnir_directories_with_hdr, search_key='raw_rd_rf')

# Creates a list of "coupled" SWIR and VNIR files based on their sampling time and date. If the sampling times are the same, 
# directories of the files are extracted and entered into a 2D list (one column for SWIR directories, second column for VNIR directories)
pairs_for_overlap = matching_vnir_swir_files(swir_directories_with_hdr,vnir_directories_with_hdr)
# Removes the first row from the list since it contains empty values
pairs_for_overlap=pairs_for_overlap.pop(0)

# Makes a new directory that saves the SWIR flipped image as backup
flipped_image_dir = os.path.join(script_path, 'flipped_SWIR')
# Creates a folder called "flipped_SWIR" off the main pathway 
#os.mkdir(flipped_image_dir)

# Makes a new directory that saves the overlapped images as backup
overlapped_images_dir = os.path.join(script_path,'overlapped_SWIR_VNIR')
#os.mkdir(overlapped_images_dir)
for pair in pairs_for_overlap:
    print(pair)
    my_list = pair
    first_element_inner_list = my_list[0][0]
    print("First element of the inner list:", first_element_inner_list)
    SWIR = pair[1]
    print('SWIR? \t', SWIR)
    VNIR = pair[0]
    print('VNIR? \t', VNIR)
    swir_band = opening_hyperspectral_file(SWIR)
    vnir_band = opening_hyperspectral_file(VNIR)
    [flipped_SWIR, unique_image_name]= horizontal_flipping(swir_band,flipped_image_dir)
    # Align SWIR image with VNIR image based on object edges
    overlapping_image(swir_band, vnir_band, overlapped_images_dir, unique_image_name)


