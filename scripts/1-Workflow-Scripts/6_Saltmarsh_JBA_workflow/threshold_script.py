"""-------------------------------------------------------------
    Script Name:   threshold_script.py
    Description:   Workflow for LE production of saltmarsh binary raster map developed with JBA Consulting. Developed with accompanying 'input.ini' and 'open_cmd.bat' files.
    Created By:    JBA Consulting, prepared for Living England project, Natural England
    Licence: MIT Licence
    Date Created:  2021-01-15
    Date Modified: 2024-08-13
    Versioning:
    GDAL/OGR version 3.5.1, QGIS LTR v.3.22.10, SAGA v.7.8.2, Python 3.9.5
-------------------------------------------------------------"""

import subprocess
import sys
import os
from osgeo import gdal, osr
import time
import configparser

'''
Useage notes and limitations:
OSGeo4W installer LTR should be installed prior to use on the users machine. This workflow has been developed for Windows use. 
The 'open_cmd.bat' file needs to be placed inside the 'C:\OSGeo4W' folder. 
This allows SAGA tools to be called within the OSGeo4W Shell and should be run to set up OSgeo4W and SAGA environment.
The `input.ini` file contains all the filepaths and variables for running the script. 
Please input a ZONE parameter to match the biogeographic zone to run this on.

'''

# Check file location of GDAL_CALC
GDAL_CALC = r'C:\\OSGeo4W\\apps\\Python39\\Scripts\\gdal_calc.py'

# Inputs - user defined look at input.ini file
config = configparser.ConfigParser()
config.read('C:\OSGeo4W\input.ini')
ZONE = host = config['filepath']['ZONE']
WORKING_DIR = host = config['filepath']['WORKING_DIR']
DTM = host = config['filepath']['DTM']
SPRING = host = config['filepath']['SPRING']
AUTUMN = host = config['filepath']['AUTUMN']
AUTUMN_S1_ASC = host = config['filepath']['AUTUMN_S1_ASC']
SPRING_S1_DESC = host = config['filepath']['SPRING_S1_DESC']
AUTUMN_S1_DESC = host = config['filepath']['AUTUMN_S1_DESC']
SPRING_S1_ASC = host = config['filepath']['SPRING_S1_ASC']
SEA_DEFENCES = host = config['filepath']['SEA_DEFENCES']
LIDAR_2m = host = config['filepath']['LIDAR_2m']
OS_DATA_FORESHORE = host = config['filepath']['OS_DATA_FORESHORE']
OS_DATA_TIDALWATER = host = config['filepath']['OS_DATA_TIDALWATER']
SEA_LEVELS = host = config['filepath']['SEA_LEVELS']

# Outputs - stored in a working dir
SPRING_NDVI = f'{ZONE}_Spring_S2_NDVI.tif'
AUTUMN_NDVI = f'{ZONE}_Autumn_S2_NDVI.tif'
SPRING_NDII = f'{ZONE}_Spring_S2_NDII.tif'
AUTUMN_NDII = f'{ZONE}_Autumn_S2_NDII.tif'
AUTUMN_S1_ABS_ASC = f'{ZONE}_s1_GRD_Autumn2019_Asc_median_ABS_Byte.tif'
SPRING_S1_ABS_DESC = f'{ZONE}_s1_GRD_Spring2020_Desc_median_ABS_Byte.tif'
AUTUMN_S1_ABS_DESC = f'{ZONE}_s1_GRD_Autumn2019_Desc_median_ABS_Byte.tif'
SPRING_S1_ABS_ASC = f'{ZONE}_s1_GRD_Spring2020_Asc_median_ABS_Byte.tif'
SEA_DEFENCES_BURN = f'{ZONE}_EA_SeaDefence_BGZ07_AT_Burn10_Byte.tif'
SEA_DEFENECES_EX_SHR = f'{ZONE}_EA_SeaDefence_BGZ07_AT_Burn10_ExpandAndShrink.tif'
SEA_DEFENECES_EX_SHR_BYTE = f'{ZONE}_EA_SeaDefence_BGZ07_AT_Burn10_ExpandAndShrink_byte.tif'
LIDAR_CLIP = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_Complete_Clipped_Max10m.tif'
LIDAR_COMPLETE_CLIP = f'{ZONE}_EA_SeaDefence_AT_Burn10_ExpandAndShrink_LIDAR_2m_DTM_Composite_2020_Complete_Clipped_Max10m.tif'
CLIPPED_SHORE = f'main_Foreshore_clip.shp' # do not change these
BUFFERED_SHORE = f'main_Foreshore_clip_buffer.shp' # do not change these
DESOLVED_SHORE = f'main_Foreshore_clip_buffer_dissolved.shp' # do not change these
LIDAR_COMPLETE_CLIP_SHORE = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_Complete_Clipped_Max10m_OS_Foreshore_2km_Clipped.tif'
HAT_SEEDS = f'{ZONE}_HAT_seeds_10m.tif'
MHWS_SEEDS = f'{ZONE}_MHWS_seeds_10m.tif'
HAT_DEPTH = f'{ZONE}_HAT_depth.tif'
MHWS_DEPTH = f'{ZONE}_MHWS_depth.tif'
HAT_SURFACE =f'{ZONE}_LIDAR_2m_DTM_Composite_2020_Complete_Clipped_Max10m_OS_Foreshore_2km_Clipped_HAT.tif'
MHWS_SURFACE = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_Complete_Clipped_Max10m_OS_Foreshore_2km_Clipped_MHWS.tif'
HAT_THRES = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_HAT_4p72_mask_byte.tif'
MHWS_THRES = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_MHWS_3p93_mask_byte.tif'
HAT_SIEVE = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_HAT_4p72_mask_byte_sieved.tif'
MHWS_SIEVE = f'{ZONE}_LIDAR_2m_DTM_Composite_2020_MHWS_3p93_mask_byte_sieved.tif'
FORESHORE = f'main_Foreshore_10m.tif'
TIDALWATER = f'main_TidalWater_10m.tif'
HAT_ONLY = f'{ZONE}_HAT_only.tif'
MHWS_ONLY = f'{ZONE}_MHWS_JBA_only.tif'
MLWS_OS_ONLY = f'{ZONE}_MLWS_OS_only.tif'
TIDALSTAGES = f'{ZONE}_TidalStages.tif'
DTM_SLOPE = f'{ZONE}_EA_LiDAR_DTM_Slope.tif'
THRESHOLD_UNION = f'{ZONE}_Thresholding_byte_union.tif'

# For GDAL Calc
NDVI_SPRING_outfile = f'--outfile={SPRING_NDVI}'
NDVI_SPRING_scaled_outfile = f'--outfile={SPRING_NDVI[:-4]}_scaled.tif'
NDVI_AUTUMN_outfile = f'--outfile={AUTUMN_NDVI}'
NDII_SPRING_outfile = f'--outfile={SPRING_NDII}'
NDII_AUTUMN_outfile = f'--outfile={AUTUMN_NDII}'
AUTUMN_S1_ABS_ASC_outfile = f'--outfile={AUTUMN_S1_ABS_ASC}'
SPRING_S1_ABS_DESC_outfile = f'--outfile={SPRING_S1_ABS_DESC}'
AUTUMN_S1_ABS_DESC_outfile = f'--outfile={AUTUMN_S1_ABS_DESC}'
SPRING_S1_ABS_ASC_outfile = f'--outfile={SPRING_S1_ABS_ASC}'
LIDAR_COMPLETE_CLIP_outfile = f'--outfile={LIDAR_COMPLETE_CLIP}'
HAT_THRES_outfile = f'--outfile={HAT_THRES}'
MHWS_THRES_outfile = f'--outfile={MHWS_THRES}'
HAT_ONLY_outfile = f'--outfile={HAT_ONLY}'
MHWS_ONLY_outfile = f'--outfile={MHWS_ONLY}'
MLWS_OS_ONLY_outfile = f'--outfile={MLWS_OS_ONLY}'
TIDALSTAGES_outfile = f'--outfile={TIDALSTAGES}'
THRESHOLD_UNION_outfile = f'--outfile={THRESHOLD_UNION}'

# Set working directory
# check to see that the working directory exists
if not os.path.exists(WORKING_DIR):
   os.makedirs(WORKING_DIR)
os.chdir(WORKING_DIR)

# function to calculate extents of DTM 
def read_raster(raster_in):
    # Function to get raster extents for input into -te flag
    ## Open raster_in as ds (abbreviation for dataset)
    ds=gdal.Open(raster_in)

    ## using GetGeoTransform we can get the upper left X and upper left y coordinates
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    
    ## Calculate lower right x and lower right y we have the coordinates to build a polygon
    ## these values will be returned and be inputs into build bounds
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    coords = [str(ulx), str(lry), str(lrx), str(uly)]

    return coords

#function to run through main workflow
def run():
    
    coords = read_raster(DTM)

    # All the GDAL calls
    # NDVI - calculate NDVI index for both Sentinel-2 composites
    print(f"NDVI being built for {SPRING}")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Float32', '-A',  SPRING, '--A_band', '8', '-B',  SPRING, '--B_band', '4', '--calc=(A.astype(float)-B.astype(float))/(A.astype(float)+B.astype(float))', NDVI_SPRING_outfile])
    
    print(f"NDVI scaled being built for {SPRING}")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--NoDataValue', '-9999', '--type=Int16', '-A', SPRING,  '--A_band', '8', '-B', SPRING, '--B_band', '4', '--calc=(A.astype(float)-B.astype(float))/(A.astype(float)+B.astype(float))*100', NDVI_SPRING_scaled_outfile]) 

    print(f"NDVI being built for {AUTUMN}")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Float32', '-A', AUTUMN, '--A_band', '8', '-B', AUTUMN, '--B_band', '4', '--calc=(A.astype(float)-B.astype(float))/(A.astype(float)+B.astype(float))', NDVI_AUTUMN_outfile])

    # NDII - calculate NDII (Normalised difference Infrared index) index for both Sentinel-2 composites
    print(f"NDII being created - Spring")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Float32', '-A', SPRING, '--A_band', '12', '-B', SPRING, '--B_band', '11', '--calc=(A.astype(float)-B.astype(float))/(A.astype(float)+B.astype(float))', NDII_SPRING_outfile])

    print(f"NDII being created - Autumn")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Float32', '-A', AUTUMN, '--A_band', '12', '-B', AUTUMN, '--B_band', '11', '--calc=(A.astype(float)-B.astype(float))/(A.astype(float)+B.astype(float))', NDII_AUTUMN_outfile])

    # S1 - extract VH backscatter (Band 2) in Absolute value to allow for saving in byte format
    print(f"S1 - extract VH backscatter in Absolute Value")
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', AUTUMN_S1_ASC, '--A_band', '2', '--calc=abs(A)', AUTUMN_S1_ABS_ASC_outfile])

    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', SPRING_S1_DESC, '--A_band', '2', '--calc=abs(A)', SPRING_S1_ABS_DESC_outfile])

    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', AUTUMN_S1_DESC, '--A_band', '2', '--calc=abs(A)', AUTUMN_S1_ABS_DESC_outfile])

    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', SPRING_S1_ASC, '--A_band', '2', '--calc=abs(A)', SPRING_S1_ABS_ASC_outfile])

    # Creation of the DTM with EA Sea defences 
    # Sea Defence (using AIMS Spatial Flood Defences (inc. standardised attributes) ( https://environment.data.gov.uk/dataset/8e5be50f-d465-11e4-ba9a-f0def148f590))
    print("creating DTM with Sea Defence")
    subprocess.call(['gdal_rasterize', '-at', '-l',  'Spatial_Flood_Defences_Including_Standardised_Attributes', '-burn', '10', '-tr', '10', '10', '-te', coords[0], coords[1], coords[2], coords[3], '-ot', 'Byte', '-of', 'GTiff', SEA_DEFENCES, SEA_DEFENCES_BURN])

    # SAGA Expand and Shrink command
    print(f'Using SAGA tools to Expand and Shrink')
    subprocess.call(['saga_cmd', 'grid_tools', 'Shrink and Expand', '-INPUT', SEA_DEFENCES_BURN, '-OPERATION', '3', '-CIRCLE', '1', '-RADIUS', '1', '-EXPAND', '3', '-RESULT', SEA_DEFENECES_EX_SHR])
    
    # SAGA SGRD output has two issues, 1) Float32 with No-data of '-99999', and 2) no CRS. Therefore, need to translate SAGA output to binary (0-1) byte:
    print(f'Convert Expand and Shrink output to Byte and resolve CRS issue')
    subprocess.call(['gdal_translate', '-ot', 'Byte', '-a_srs', 'EPSG:27700', '-a_nodata', '255', SEA_DEFENECES_EX_SHR, SEA_DEFENECES_EX_SHR_BYTE]) 

    # Using 2m composite Lidar - extract 10m maximum for the same extent as original 10m DTM 
    print(f'extract 10m max for same extent as orginal 10m DTM')
    subprocess.call(['gdalwarp', '-r', 'max', '-tr', '10', '10', '-te', coords[0], coords[1], coords[2], coords[3], LIDAR_2m, LIDAR_CLIP])

    # Merge Max 10m Lidar with AIMS coastal defence 10m raster
    print(f'merging max Lidar with AIMS')
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '-A', SEA_DEFENECES_EX_SHR_BYTE, '-B', LIDAR_CLIP, '--calc=(A+B)', LIDAR_COMPLETE_CLIP_outfile])
    
    # clipping
    
    ## Note: the OS Foreshore data was previously merged for this BZ zone to create Foreshore.shp
    print(f'clipping OS OpenMap to DTM extents')
    subprocess.call(['ogr2ogr', '-f', 'ESRI Shapefile', '-clipsrc', coords[0], coords[1], coords[2], coords[3], CLIPPED_SHORE, OS_DATA_FORESHORE])
    
    # buffering
    print(f'buffering OS clip')
    subprocess.call(['ogr2ogr', '-f', 'ESRI Shapefile', BUFFERED_SHORE, CLIPPED_SHORE, '-dialect', 'sqlite', '-sql', 'select ST_buffer(geometry, 2500) as geometry FROM main_Foreshore_clip'])
    
    # dissolving the buffer
    print(f'dissolving the buffer')
    subprocess.call(['ogr2ogr', '-f', 'ESRI Shapefile', DESOLVED_SHORE, BUFFERED_SHORE, '-dialect', 'sqlite', '-sql', 'select ST_union(Geometry) from main_Foreshore_clip_buffer'])
    
    # clip to OS tidal water buffered by 2km (to avoid inland flooding through river networks)
    print(f'clipping the merge of max Lidar with AIMS')
    subprocess.call(['gdalwarp', '-tr', '10', '10', '-of', 'GTiff', '-cutline', DESOLVED_SHORE, '-cl', 'main_Foreshore_clip_buffer_dissolved', LIDAR_COMPLETE_CLIP, LIDAR_COMPLETE_CLIP_SHORE])
    
    # rasterize HAT estimates per point - generate HAT 'seeds' for Lake flooding tool. Important no-data set to 0
    print(f'Rasterize extreme sea levels')
    subprocess.call(['gdal_rasterize', '-l', 'Coastal_Design_Sea_Levels_Coastal_Flood_Boundary_Extreme_Sea_Levels', '-a', 'hat_od', '-tr', '10', '10', '-a_nodata', '0', '-te', coords[0], coords[1], coords[2], coords[3], '-ot', 'Float32', '-of', 'GTiff', SEA_LEVELS, HAT_SEEDS])
  
    
    # rasterize MHWS estimates per point - generate MHWS 'seeds' for Lake flooding tool. Important no-data set to 0
    print(f'rasterize MHWS')
    subprocess.call(['gdal_rasterize', '-l', 'Coastal_Design_Sea_Levels_Coastal_Flood_Boundary_Extreme_Sea_Levels', '-a', 'mhws_od', '-tr', '10', '10', '-a_nodata', '0', '-te', coords[0], coords[1], coords[2], coords[3], '-ot', 'Float32', '-of', 'GTiff', SEA_LEVELS, MHWS_SEEDS])
   
  
    # Run SAGA Lake Flood (Absolute levels), once for HAT and once for MHWS
    
    ###CHECK LIDAR_COMPLETE_CLIP_SHORE is correct input
    
    print(f'Abs levels HAT')
    subprocess.call(['saga_cmd', 'ta_hydrology', 'Lake Flood', '-ELEV', LIDAR_COMPLETE_CLIP_SHORE, '-SEEDS', HAT_SEEDS, '-LEVEL', 'true', '-OUTDEPTH', HAT_DEPTH, '-OUTLEVEL', HAT_SURFACE])
       
    print(f'Abs levels MHWS')
    subprocess.call(['saga_cmd', 'ta_hydrology', 'Lake Flood', '-ELEV', LIDAR_COMPLETE_CLIP_SHORE, '-SEEDS', MHWS_SEEDS, '-LEVEL', 'true', '-OUTDEPTH', MHWS_DEPTH, '-OUTLEVEL', MHWS_SURFACE])
    
    
    ## Having run for both HAT and HMWS, save the output for both as GTiffs, and then run raster calc to extract a binary mask extent for the value for each level (e.g. in BGZ07, the HAT was 4.72m, and the MHWS was 3.93m)
    print(f'extract binary mask for HAT level')
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', HAT_DEPTH, '--calc=(numpy.greater(A,0))', HAT_THRES_outfile])
    
    print(f'extract binary mask for MHWS level')
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', MHWS_DEPTH, '--calc=(numpy.greater(A,0))', MHWS_THRES_outfile])
    
    # Sieve output to remove elevation values outside inter-tidal
    print(f'Sieve output for HAT')
    subprocess.call(['gdal_sieve.bat', '-st', '10', '-4', '-of', 'GTiff', HAT_THRES, HAT_SIEVE])
    
    print(f'Sieve output for MHWS')
    subprocess.call(['gdal_sieve.bat', '-st', '10', '-4', '-of', 'GTiff', MHWS_THRES, MHWS_SIEVE])
    
    ## Tidal stages

    # Using OS Local Data - "M:\OS OpenData\GB\OS OpenMap - Local\Vector\OpenMapLocal_Vector.gdb" - available from https://osdatahub.os.uk/downloads/open/OpenMapLocal
    ## Note: the OS Foreshore data was merged to create Foreshore.shp, as was the TidalWater data. need to ensure this covers the DTM bounds.
        
    print(f'rasterize Foreshore data')
    subprocess.call(['gdal_rasterize',  '-burn', '1', '-tr', '10', '10', '-a_nodata', '-9999', '-te', coords[0], coords[1], coords[2], coords[3], '-ot', 'Byte', '-of', 'GTiff', OS_DATA_FORESHORE, FORESHORE])
      
    print(f'rasterize Tidal Water data')
    subprocess.call(['gdal_rasterize',  '-burn', '1', '-tr', '10', '10', '-a_nodata', '-9999', '-te', coords[0], coords[1], coords[2], coords[3], '-ot', 'Byte', '-of', 'GTiff', OS_DATA_TIDALWATER, TIDALWATER])
    
    #extract single masks for each tidal stage
    print(f'extracting single masks for each tidal stage')
    # remove MWHS extent from HAT so that only HAT extent remains
    print(f'removing MWHS extent from HAT')
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A',  HAT_SIEVE, '-B',  MHWS_SIEVE, '--calc=(A-B)', HAT_ONLY_outfile])
        
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', MHWS_SIEVE, '-B', TIDALWATER, '--calc=numpy.where(numpy.greater(A,0),A-B,0)', MHWS_ONLY_outfile])
    
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', TIDALWATER, '-B', FORESHORE, '--calc=A-B', MLWS_OS_ONLY_outfile])
    
    # Calculate Tidal stages (Combination of all 4 stages)
    print(f'Calculate Tidal stages, combination of all 4 stages)')
    subprocess.call([sys.executable, GDAL_CALC, '--co=COMPRESS=LZW', '--type=Byte', '-A', HAT_ONLY, '-B', MHWS_ONLY, '-C', FORESHORE, '-D', MLWS_OS_ONLY, '--calc=(A*4)+(B*3)+(C*2)+(D*1)', TIDALSTAGES_outfile])
    
    # Create slope from 10m DTM
    print(f'create slope from 10m DTM')
    subprocess.call(['gdaldem', 'slope', DTM, DTM_SLOPE, '-of', 'GTiff', '-b', '1', '-s', '1.0'])
       
    # Using updated Tidal stages (derived from AIMS coastal defences)
    print(f'Using updated Tidal stage')
    print(f'Final step - thresholding')
    subprocess.call([sys.executable, GDAL_CALC, '--extent=union', '--co=COMPRESS=LZW', '--type=Byte', '-A', DTM_SLOPE, '-B', TIDALSTAGES, '-C', AUTUMN_NDVI, '-D', SPRING_NDVI, '-E', AUTUMN_S1_ABS_ASC, '-F', AUTUMN_NDII, '--calc=numpy.less(A,5)*numpy.greater(B,1)*numpy.greater(C,0.2)*numpy.greater(D,0.15)*numpy.less(E,26)*numpy.less(F,-0.18)', THRESHOLD_UNION_outfile]) 
    
    print(f'process completed')
    
    
    print(f'clean up - add projection')
    lsinputs = [MHWS_SURFACE, HAT_SURFACE, HAT_THRES, MHWS_THRES, HAT_SIEVE, MHWS_SIEVE, HAT_DEPTH, MHWS_DEPTH, SEA_DEFENECES_EX_SHR, HAT_ONLY]
    
    for layer in lsinputs:
        outlayer = f'{layer[:-4]}_projected.tif'
        subprocess.call(['gdal_translate', '-of', 'GTiff', '-co', 'COMPRESS=LZW', '-a_srs', 'EPSG:27700', layer, outlayer])
        # delete file
        os.remove(layer)
        # rename file
        os.rename(outlayer, layer)
    
if __name__ == "__main__":
    start_time = time.perf_counter()
    run()
    print(f'time taken {time.perf_counter() - start_time} seconds')