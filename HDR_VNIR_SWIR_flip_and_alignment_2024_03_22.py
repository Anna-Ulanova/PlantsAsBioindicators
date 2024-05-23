import os
import numpy as np
import pandas as pd
import re
#import cv2
#from matplotlib import pyplot as plt

# AU Note: matplotlib as plt can allow you to give the images as an  output which is helpful 
# helpful during debugging. Uncomment the plt.imshow() lines to see image outputs
 
####################################### List of Functions #######################################
# Function takes in the main directory where the SWIR and VNIR folders are located, an empty and 
# initiated list that will contain the directories of the HDR files, as well as a search key for
# the files you want to isolate (EX:raw_rd_rf.hdr)
# Returns a list of directories where HDR files are found
def directories_of_hdr_files(script_path, search_key):
    swir_dir = []
    vnir_dir = []
    for root, _ , files in os.walk(script_path, topdown=True):
    # Change search key 
        if search_key in files: 
            hdr_file_directory = os.path.join(root,search_key)
            if 'SWIR' in hdr_file_directory: 
                swir_dir.append(hdr_file_directory)
            if 'VNIR' in hdr_file_directory: 
                vnir_dir.append(hdr_file_directory)        
    return [swir_dir, vnir_dir]
##-------------------------------------------------------------------------------------------###
# Function takes in the directories of hdr files for SWIR and VNIR types, 
# couples hdr directories if their sampling time matches (e.i. identical last 2 folders)
# Returns the two dimensation list that lists the "coupled" HDR files
def matching_vnir_swir_files(swir_directories_with_hdr, vnir_directories_with_hdr):
    # Pairs HDR files if the last 2 folders match, EX: both are in GH_20230724/1_1_20230724_2023_07_24_08_16_25/raw_rd_rf.hdr
    # Creates a list pairs_for_overlap which can be overlapped/processed
    blank = []
    for directory in swir_directories_with_hdr: 
        sub_key = directory.split('\\')
        sub = '\\'.join(sub_key[-3:])
        for vnir in vnir_directories_with_hdr:
            if sub in vnir: 
                blank.append([directory, vnir])
    # Turn on if you want to save the coupled list, useful for sanity checks. 
    df = pd.DataFrame(blank, columns=['SWIR', 'VNIR'])
    df.columns = ['SWIR', 'VNIR']
    df.to_excel('D:/reference.xlsx')
    return df
##-------------------------------------------------------------------------------------------###
def time_offset(perfect_pairs, all_swir, all_vnir):
    oof = 0
    improper = []
    time_offset_pairs = []
    for swir in all_swir: 
        # print(swir)
        if  (swir not in perfect_pairs.loc[:,'SWIR']):
            formatted_swir = r'{}'.format(swir)
            sp = re.split(r'[_]',formatted_swir)
            no_second = '_'.join(sp[3:9])
           # print(sp)
           # print(no_second)
            if any(no_second in x for x in all_vnir):
                    weird_bracket = [x for x in all_vnir if no_second in x]
                    vnir = weird_bracket[0]
                    # print(vnir)
                    time_offset_pairs.append([swir, vnir])
            else:
                improper.append(swir)
        else: 
            oof = oof+1
    print(oof)    
    return [time_offset_pairs, improper]

##-------------------------------------------------------------------------------------------###

def unique_name(directory):
    print(directory)
    parts = directory.split('\\')
    desired_substring = '_'.join(parts[2:4])
    overlapped_images_dir = target_dir+'\\'+desired_substring+'_stacked'
    return overlapped_images_dir
####################################### Input variables #######################################

# Change target directory after testing 
gen_directory = "D:/"
target_directory = 'D:/overlapped_SWIR_VNIR'
# Changes to directory to target path and reformats the delimiters to system's operating system
os.chdir(gen_directory)
script_path = os.getcwd()

[swir_directories_with_hdr, vnir_directories_with_hdr] = directories_of_hdr_files(script_path, search_key='raw_rd_rf')

#Creates a list of "coupled" SWIR and VNIR files based on their sampling time and date. If the sampling times are the same, 
#directories of the files are extracted and entered into a 2D list (one column for SWIR directories, second column for VNIR directories)

pairs_for_overlap = matching_vnir_swir_files(swir_directories_with_hdr,vnir_directories_with_hdr)

[time_offset_pairs, improper] = time_offset(pairs_for_overlap, swir_directories_with_hdr, vnir_directories_with_hdr)

time_offset_pairs_df = pd.DataFrame(time_offset_pairs, columns=['SWIR','VNIR'])    
overall_pairs = pd.concat([pairs_for_overlap,time_offset_pairs_df], ignore_index=True)
os.chdir(target_directory)
target_dir = os.getcwd()
time_offset_pairs_df.to_csv('D:/time_offset_pairs.txt', sep='\t')
overall_pairs.to_csv('D:/overall_pairs.txt', sep='\t')
save_name = []
overall_save_names = overall_pairs['SWIR'].apply(unique_name)
save_names = pd.DataFrame(np.array(overall_save_names), columns=['Save Names'])
print(save_names.shape)
save_names.to_excel('D:/save_name_nomenclature.xlsx')
new = pd.concat([overall_pairs, save_names],axis=1, ignore_index=True)
new = new.drop_duplicates(ignore_index=True)
new.to_csv('D:/wow.csv', sep=',')
