import os
import subprocess
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
    df.to_excel('D:/reference.xlsx')
    return df
##-------------------------------------------------------------------------------------------###
def time_offset(perfect_pairs, all_swir, all_vnir):
    i = 0
    print(perfect_pairs)
    improper = []
    time_offset_pairs = []
    for swir in all_swir: 
        print(swir)
        for i in range(len(perfect_pairs.index)):
            print(perfect_pairs.iloc[i,0])
            if  (swir not in perfect_pairs.iloc[i,0]):
                formatted_swir = r'{}'.format(swir)
                sp = re.split(r'[_]',formatted_swir)
                no_second = '_'.join(sp[1:11])
                print(sp)
                print(no_second)
                if no_second in all_vnir:
                        time_offset_pairs.append([swir, no_second])
                else:
                    improper.append(swir)
            else: 
                print('oof')
            i=i+1
    return [time_offset_pairs, improper]

##-------------------------------------------------------------------------------------------###

def unique_name(directory):
    parts = directory.split('\\')
    desired_substring = '_'.join(parts[2:4])

    return desired_substring
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
# Removes the first row from the list since it contains empty values
#pairs_for_overlap=pairs_for_overlap.pop(0)

[time_offset_pairs, improper] = time_offset(pairs_for_overlap, swir_directories_with_hdr, vnir_directories_with_hdr)
print(time_offset_pairs)
print('###########################')
time_offset_pairs_df = pd.DataFrame(time_offset_pairs, columns=['SWIR','VNIR'])    
time_offset_pairs_df = time_offset_pairs_df.drop_duplicates()
overall_pairs = pd.concat([pairs_for_overlap,time_offset_pairs_df], ignore_index=True)
os.chdir(target_directory)
target_dir = os.getcwd()
time_offset_pairs_df.to_excel('D:/time_offset_pairs.xlsx')
overall_pairs.to_excel('D:/overall_pairs.xlsx')
print(overall_pairs)
print(overall_pairs.shape)
# Wait for ENVI to open before proceeding
# input("Press Enter to continue once ENVI has launched...")
# Once ENVI is open, you can execute IDL commands or scripts within the ENVI/IDL environment
# For example, executing an IDL script:
#os.chdir(r'C:\Program Files\Harris\ENVI56\IDL88\help\online_help\Subsystems\idl')
#idl_script_path = r"C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\IDL Image Processing\main_operation.pro"
save_name = []
# for pair in overall_pairs:
#     for subpair in pair: 
#         SWIR = subpair.iloc[:,1]
#         VNIR = subpair[:,2]
#         #print(SWIR, '\n\n\n\n')
#         #print(VNIR, '\n\n\n\n')
#         collection_date = unique_name(SWIR)
#         # Makes a new directory that saves the overlapped images as backup
#         overlapped_images_dir = target_dir+'\\'+collection_date+'\\stacked'
#         print(overlapped_images_dir)
#         save_name.append(overlapped_images_dir)
#         #subprocess.call(["idl", "-e", r"C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\IDL Image Processing\main_operation.pro , '{}', '{}', '{}'".format(SWIR, VNIR, overlapped_images_dir)])
#         #idl_command = ["idl", "-e", "@" + idl_script_path, SWIR, VNIR,overlapped_images_dir]
#         # Execute IDL script within the ENVI/IDL environment
#         #subprocess.call(idl_command)
print(type(overall_pairs))
print(overall_pairs.head())
overall_pairs= overall_pairs.iloc[1:len(overall_pairs.shape)]
for index, row in overall_pairs.iterrows():
    print(pair)
    SWIR = pair[0]
    VNIR = pair[1]
    #print(SWIR, '\n\n\n\n')
    #print(VNIR, '\n\n\n\n')
    collection_date = unique_name(SWIR, target_dir)
    # Makes a new directory that saves the overlapped images as backup
    overlapped_images_dir = target_dir+'\\'+collection_date+'\\stacked'
    print(overlapped_images_dir)
    save_name.append(overlapped_images_dir)
    #subprocess.call(["idl", "-e", r"C:\Users\RDCRLAAU\Desktop\Plant as bioindicators\IDL Image Processing\main_operation.pro , '{}', '{}', '{}'".format(SWIR, VNIR, overlapped_images_dir)])
    #idl_command = ["idl", "-e", "@" + idl_script_path, SWIR, VNIR,overlapped_images_dir]
    # Execute IDL script within the ENVI/IDL environment
    #subprocess.call(idl_command)
new = pairs_for_overlap.insert(save_name)
new = pd.DataFrame(new)
new.to_excel('D:/wow.xlsx')
