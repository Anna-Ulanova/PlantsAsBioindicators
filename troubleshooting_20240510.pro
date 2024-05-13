;
;swir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\SWIR\1_1_20230724_2023_07_24_08_16_25\EXTRA\raw_rd_rf'
;vnir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\VNIR\1_1_20230724_2023_07_24_08_16_24\EXTRA\raw_rd_rf'
;
;;###################################################################################################################
;e = ENVI()
;; Opening VNIR input file
;swir_raster = e.OpenRaster(swir_dir)
;
;; VNIR dimensions
;vnir_columns = swir_raster.NCOLUMNS
;
;vnir_rows = swir_raster.NROWS
;
;
;; Open Metadata file to get wavelengths
;s_metadata = swir_raster.METADATA
;wavelengths = s_metadata['WAVELENGTH']
;lambda_min = 1350
;lambda_max = 1450
;l1 = 1800 
;l2 = 1900
;indices_subset_1 = WHERE(wavelengths LT lambda_min OR wavelengths GT lambda_max AND wavelengths LT l1 OR wavelengths GT l2, count)
;;i1 = LIST(indices_subset_1)
;;i2 = LIST(indices_subset_2)
;;combined_indices = i1+i2
;;PRINT , combined_indices
;;unique_list = combined_indices[UNIQ(combined_indices, SORT(combined_indices))]
;SUBSET = ENVISubsetRaster(swir_raster, BANDS=indices_subset_1)
;PRINT, 'HI'
;
;View1 = e.GetView()
;
;Layer1 = View1.CreateLayer(SUBSET)
;
;END

swir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\SWIR\1_1_20230724_2023_07_24_08_16_25\EXTRA\raw_rd_rf'
vnir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\VNIR\1_1_20230724_2023_07_24_08_16_24\EXTRA\raw_rd_rf'

;###################################################################################################################
e = ENVI()
; Opening VNIR input file
vnir_raster = e.OpenRaster(vnir_dir)

; VNIR dimensions
vnir_columns = vnir_raster.NCOLUMNS

vnir_rows = vnir_raster.NROWS


; Open Metadata file to get wavelengths
v_metadata = vnir_raster.METADATA
wavelengths = v_metadata['WAVELENGTH']

; Subset VNIR to exclude any bands that have a wavelength more than 920 nm
; file:///C:/Program%20Files/Harris/ENVI56/IDL88/help/online_help/Subsystems/envi/Content/ExtendCustomize/ENVISubsetRaster/ENVISubsetRaster.htm
lambda_max = 920
indices_subset = WHERE(wavelengths LT lambda_max, count)
SUBSET_vnir = ENVISubsetRaster(vnir_raster, BANDS= indices_subset)

PRINT, vnir_columns
; View raster
View1 = e.GetView()

Layer1 = View1.CreateLayer(SUBSET_vnir)

PRINT, SIZE(SUBSET_vnir, /TYPE)
;###################################################################################################################

; Open SWIR input file
swir_raster = e.OpenRaster(swir_dir)
; View raster
;View2 = e.GetView()
;
;Layer2 = View2.CreateLayer(swir_raster)

; SWIR dimensions
swir_columns = swir_raster.NCOLUMNS

swir_rows = swir_raster.NROWS


; Open Metadata file to get wavelengths
s_metadata = swir_raster.METADATA
wavelengths = s_metadata['WAVELENGTH']


; Subset SWIR to exclude bands that have wavelengths 1350-1450, and 1800-1900
lambda_min_1 = 1350
lambda_max_1 = 1450
lambda_min_2 = 1800
lambda_max_2 = 1900
;wavelengths_2 = SUBSET
indices_subset_1 = WHERE(wavelengths LT lambda_min_1 OR wavelengths GT lambda_max_1 AND wavelengths LT lambda_min_2 OR wavelengths GT lambda_max_2, count)
SUBSET_swir = ENVISubsetRaster(swir_raster, BANDS= indices_subset_1)

; Calculate the overlap region
overlap_rows = MIN(swir_rows, vnir_rows)
overlap_cols = MIN(swir_columns, vnir_columns)

; Get the task from the catalog of ENVITasks
Task = ENVITask('MirrorRaster')

; Define inputs
Task.INPUT_RASTER = SUBSET_swir

Task.HORIZONTAL = 'true'

; Run the task
Task.Execute

; Get the collection of objects currently in the Data Manager
DataColl = e.Data

; Add the output to the Data Manager
DataColl.Add, Task.OUTPUT_RASTER

mirror_swir = e.OpenRaster(Task.OUTPUT_RASTER_URI)
View2 = e.GetView()

Layer2 = View2.CreateLayer(mirror_swir)
;###################################################################################################################


sized = ENVISubsetRaster(mirror_swir, SUB_RECT=[0,0,vnir_columns, vnir_rows])

; View raster
View3 = e.GetView()

Layer3 = View3.CreateLayer(sized)


;###################################################################################################################

; Create an ENVIRasterSubset object for Raster1
subset_Raster1 = ENVIRasterSubset(mirror_swir)

; Define the subset parameters based on the dimensions of Raster2
subset_params = {cols: vnir_columns, rows: vnir_rows, x_start: 0, y_start: 0}

; Subset Raster1 to match the dimensions of Raster2
subset_Raster1_subset = mirror_swir.Subset(subset_params)

; Retrieve the subsetted raster data
subsetted_data = subset_Raster1_subset.GetData()


View5 = e.GetView()

Layer5 = View5.CreateLayer(subsetted_data)



Task = ENVITask('BuildingBandStack')
Task.INPUT_RASTERS = [SUBSET_vnir, subsetted_data]
Task.Execute

; Add the output to the Data Manager
e.Data.Add, Task.OUTPUT_RASTER

; Display the result
View5 = e.GetView()
Layer5 = View5.CreateLayer(Task.OUTPUT_RASTER)



;;###################################################################################################################
;
;Task = ENVITask('GenerateTiePointsByCrossCorrelation')
;
;; Define inputs
;Task.INPUT_RASTER1 = mirror_swir
;Task.INPUT_RASTER2 = SUBSET_vnir
;
;; Run the task
;Task.Execute
;
;; Get the output tie points
;TiePoints = Task.OUTPUT_TIEPOINTS
;
;
;;###################################################################################################################
;; Get the image-to-image registration task from the catalog of ENVITasks
;RegistrationTask = ENVITask('ImageToImageRegistration')
;
;; Define inputs
;RegistrationTask.INPUT_TIEPOINTS = TiePoints
;
;RegistrationTask.WARPING = 'RST'
;
;; Run the task
;RegistrationTask.Execute
;
;; Get the output raster
;WarpedRaster = RegistrationTask.OUTPUT_RASTER
;
;View3 = e.GetView()
;
;Layer3 = View3.CreateLayer(WarpedRaster)
;
;bandData = WarpedRaster.GetData()
;newRaster = ENVIRaster(bandData, URI = 'C:\Users\RDCRLAAU\Desktop\Backup\overlapped_SWIR_VNIR\subsetted_reg')
;newRaster.Save

END