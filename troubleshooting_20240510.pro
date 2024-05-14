swir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\SWIR\1_1_20230724_2023_07_24_08_16_25\EXTRA\raw_rd_rf'
vnir_dir = 'C:\Users\RDCRLAAU\Desktop\Backup\VNIR\1_1_20230724_2023_07_24_08_16_24\EXTRA\raw_rd_rf'

;###################################################################################################################
e = ENVI()
; Opening VNIR input file
vnir_raster = e.OpenRaster(vnir_dir)
; Open Metadata file to get wavelengths
v_metadata = vnir_raster.METADATA
wavelengths = v_metadata['WAVELENGTH']

; Subset VNIR to exclude any bands that have a wavelength more than 920 nm
; file:///C:/Program%20Files/Harris/ENVI56/IDL88/help/online_help/Subsystems/envi/Content/ExtendCustomize/ENVISubsetRaster/ENVISubsetRaster.htm
lambda_max = 920
indices_subset = WHERE(wavelengths LT lambda_max, count)
SUBSET_vnir = ENVISubsetRaster(vnir_raster, BANDS=indices_subset)

;View raster
View1 = e.GetView()
Layer1 = View1.CreateLayer(SUBSET_vnir)

;###################################################################################################################

; Open SWIR input file
swir_raster = e.OpenRaster(swir_dir)
; View raster
View2 = e.GetView()
Layer2 = View2.CreateLayer(swir_raster)
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
;https://www.nv5geospatialsoftware.com/docs/envisubsetrastertask.html#SUB_RECT
;;https://www.nv5geospatialsoftware.com/docs/envipixelscaleresamplerastertask.html
; VNIR dimensions
vnir_columns = SUBSET_vnir.NCOLUMNS
PRINT, vnir_columns
vnir_rows = SUBSET_vnir.NROWS
PRINT, vnir_rows
; SWIR dimensions
swir_columns = mirror_swir.NCOLUMNS
PRINT, swir_columns
swir_rows = mirror_swir.NROWS
PRINT, swir_rows
col_factor = float(swir_columns)/float(vnir_columns)
row_factor = float(swir_rows)/float(vnir_rows)
PRINT, col_factor
PRINT, row_factor
Task = ENVITask('PixelScaleResampleRaster')
Task.INPUT_RASTER = mirror_swir
; Pixels cannot be scaled dow
Task.PIXEL_SCALE = [col_factor, row_factor]
Task.Execute 
DataColl = e.Data
DataColl.Add, Task.OUTPUT_RASTER

View3 = e.GetView()
Layer3 = View3.CreateLayer(Task.OUTPUT_RASTER) 
swir_resized = Task.OUTPUT_RASTER

;###################################################################################################################
; file:///C:/Program%20Files/Harris/ENVI56/IDL88/help/online_help/Subsystems/envi/Content/ExtendCustomize/ENVITasks/ENVIBuildBandStackTask.htm
VNIR = ENVISubsetRaster(SUBSET_vnir)
Task = ENVITask('BuildBandStack')
Task.INPUT_RASTERS = [SUBSET_vnir, swir_resized]
; Run the task
Task.Execute

; Add the output to the Data Manager

e.Data.Add, Task.OUTPUT_RASTER
 
; Display the result
View = e.GetView()
Layer1 = View.CreateLayer(Task.OUTPUT_RASTER)
Layer2 = View.CreateLayer(Task.OUTPUT_RASTER, BANDS=[1])

PRINT, "Analysis Complete"
END