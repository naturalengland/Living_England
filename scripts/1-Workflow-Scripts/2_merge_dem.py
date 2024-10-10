"""-------------------------------------------------------------
    Script Name: merge_dem.py

    Description: Combines DTMs, resamples to 10x10m, and exports in desired projection

    Author: Chris Moore, Natural England

    Licence: MIT Licence
    Date Created: 2023-02-06

    Date Last Modified: 2023-06-16

    Versioning: 1.0
    Python version: 3.7
    Dependencies: 
    arcpy       v.3.1
-------------------------------------------------------------"""
# Import libraries
import os
import time
import shutil
import glob
import datetime
import arcpy

# Select dems that will be merged
lidar = True  # EA LiDAR 2022 - https://www.data.gov.uk/dataset/f0db0249-f17b-4036-9e65-309148c97ce4/national-lidar-programme
ihm = True  # EA Integrated Height Model 2019 - https://www.blueskymapshop.com/products/height-data
srtm = True  # NASA SRTM 30m v3 2000 - https://www.earthdata.nasa.gov/sensors/srtm

# Global Variables
mosaic_srtm = True  # Mosaic SRTM tiles and save as geotiff
keep_gdb = False  # Use existing temporary geodatabase (re-use projected lidar/ihm/srtm dtms)
projection = "BNG"  # Select output projection: "WGS84"/"BNG" (default: wgs84)
wdir = "D:\\DEMs\\" # Set Working Directory

# Output DEM name
out_dem = "mergedDEM_LiDAR_IHM_SRTM" + "_" + projection
print("merge_dem")
print("> loading script params ...")

# Start timer
tic = time.time()

# Reset temp geodatabase
if not keep_gdb:
    if os.path.exists(wdir + "tmp.gdb"):
        shutil.rmtree(wdir + "tmp.gdb")
    arcpy.management.CreateFileGDB(wdir, "tmp.gdb")
tdir = wdir + "tmp.gdb\\"

# Reset output DEM
for f in glob.glob(wdir + out_dem + ".*"):
    os.remove(f)

# Set processing environments
arcpy.env.workspace = wdir + "tmp.gdb"
arcpy.env.overwriteOutput = True
arcpy.env.cellSize = "MINOF"
arcpy.env.resamplingMethod = "BILINEAR"
arcpy.env.nodata = "MINIMUM"

# Set spatial projection
if projection == "BNG":
    proj = arcpy.SpatialReference(27700)  # BNG: 27700
else:
    proj = arcpy.SpatialReference(4326)  # WGS 1984: 4326 | WGS 1984 Web Mercator: 3857
arcpy.env.outputCoordinateSystem = proj

try:
    # Load EA lidar 20222
    if lidar:
        # Get path to raster
        print("> loading EA LiDAR 2022 DTM...")
        lidar_gdb = wdir + "EA_LiDAR_NationalProgramme_2022\\LiDAR_10m_DTM_2022.gdb"
        walk = arcpy.da.Walk(lidar_gdb)
        for dirpath, dirnames, filenames in walk:
            lidar_dem = filenames[0]

        # Project to proj if needed
        d = arcpy.Describe(lidar_gdb + "\\" + lidar_dem)
        if not d.spatialReference.GCS.GCSCode == proj.GCS.GCSCode and not keep_gdb:
            print("> projecting EA LiDAR 2022 to proj...")
            lidar_proj = lidar_dem + "_" + projection
            arcpy.management.ProjectRaster(lidar_gdb + "\\" + lidar_dem, lidar_proj, proj, "CUBIC")
        else:
            lidar_proj = lidar_gdb + "\\" + lidar_dem
    else:
        lidar_proj = "lidar_dem_false"

    # Load ea ihm 2019
    if ihm:
        # Get path to raster
        print("> loading EA IHM 2019 DTM...")
        ihm_gdb = wdir + "EA_Integrated-Height-Model_2019\\EA_IHM_2019_DTM_Resampled10m_National.gdb"
        walk = arcpy.da.Walk(ihm_gdb)
        for dirpath, dirnames, filenames in walk:
            ihm_dem = filenames[0]
        
        # Project to proj
        d = arcpy.Describe(ihm_gdb + "\\" + ihm_dem)
        if not d.spatialReference.GCS.GCSCode == proj.GCS.GCSCode and not keep_gdb:
            print("> projecting EA IHM 2019 to " + projection + "...")
            ihm_proj = ihm_dem + "_" + projection
            arcpy.management.ProjectRaster(ihm_gdb + "\\" + ihm_dem, ihm_proj, proj, "CUBIC")
        else:
            ihm_proj = ihm_gdb + "\\" + ihm_dem
    else:
        ihm_proj = "ihm_dem_false"

    # Load nasa srtm 30m
    if srtm:
        # Set path to raster
        print("> loading NASA SRTM 30m...")
        srtm_dir = wdir + "NASA_SRTM_30m_v3\\"
        srtm_dem = "NASA_SRTM_30m_v3_Eng"
        # Mosaic srtm tiles
        if mosaic_srtm:
            print("> mosaicing SRTM tiles...")
            # Reset output tif
            for f in glob.glob(srtm_dir + srtm_dem + ".*"):
                os.remove(f)
            # Get tile list
            tile_list = []
            for f in glob.glob(srtm_dir + "*.hgt"):
                tile_list.append(f)
            # Create output raster
            arcpy.management.CreateRasterDataset(srtm_dir, srtm_dem + ".tif", number_of_bands = 1,
                                                 raster_spatial_reference=proj, pixel_type = "16_BIT_SIGNED")
            # Mosaic raster tiles
            arcpy.management.Mosaic(tile_list, srtm_dir + srtm_dem + ".tif", "BLEND")
        
        # Project and resample
        srtm_proj10 = srtm_dem + "_" + projection + "_10m"
        if not keep_gdb:
            # Check projection
            d = arcpy.Describe(srtm_dir + srtm_dem + ".tif")
            if not d.spatialReference.GCS.GCSCode == proj.GCS.GCSCode:
                # Reproject
                print("> projecting SRTM to " + projection + "...")
                srtm_proj = srtm_dem + "_" + projection
                arcpy.management.ProjectRaster(srtm_dir + srtm_dem + ".tif", srtm_proj, proj, "CUBIC")
            else:
                srtm_proj = srtm_dir + srtm_dem + ".tif"
            # Resample to 10m
            print("> resampling SRTM to 10m pixels...")
            if lidar:
                cs = arcpy.Describe(lidar_proj).meanCellWidth
            elif ihm:
                cs = arcpy.Describe(ihm_proj).meanCellWidth
            else:
                cs = 10
            arcpy.management.Resample(srtm_proj, srtm_proj10, resampling_type="BILINEAR", cell_size=cs)
    else:
        srtm_proj10 = "srtm_dem_false"

    print("> merging DEMs...")
    # Set reverse order of dems to merge
    dem_rank_tmp = [srtm_proj10, ihm_proj, lidar_proj]

    # Remove elements not being used
    dem_rank = [d for d, tf in zip(dem_rank_tmp, [srtm, ihm, lidar]) if tf is True]

    # Set default reference raster
    arcpy.env.snapRaster = dem_rank[-1]

    # Create output raster
    arcpy.management.CreateRasterDataset(tdir, out_dem, number_of_bands=1, raster_spatial_reference=proj,
                                         pixel_type="16_BIT_SIGNED")
    
    # Mosaic dems
    arcpy.management.Mosaic(dem_rank, out_dem, "LAST")
    
    # Find null pixels
    print("> setting null values to zeros...")
    dem_merge_null = arcpy.ia.IsNull(out_dem)
    
    # change null to zeros
    dem_merge_zeros = arcpy.ia.Con(dem_merge_null, 0, out_dem, "Value = 1")

    print("> exporting merged DEM to geotiff...")
    
    # create output raster
    arcpy.management.CreateRasterDataset(wdir, out_dem + ".tif", number_of_bands=1, raster_spatial_reference=proj,
                                         pixel_type="16_BIT_SIGNED")
    
    # Save dem as geotiff
    dem_merge_zeros.save(wdir + out_dem + ".tif")

# Catch errors
except Exception as e:
    print("\nArcPy Warnings:\n" + arcpy.GetMessages(1) + "\nErrors:\n" + str(e) + "\n")

# End
print("> done")
toc = time.time()
print(datetime.timedelta(seconds=round(toc - tic)))
