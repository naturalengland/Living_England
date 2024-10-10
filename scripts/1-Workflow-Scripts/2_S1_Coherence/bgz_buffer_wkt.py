"""-------------------------------------------------------------
    Script Name:   bgz_buffer_wkt.py
    Description:   Buffers BGZs and exports as WKT polygon and shapefile
    Created By:    Chris Moore
    Date Created:  16.12.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

import time
import datetime
import arcpy
import csv
import os

# buffer size
buffer_km = 20  # km

print("bgz_buffer_wkt")
print("> loading script params ...")
# start timer
tic0 = time.time()

# set processing environments
arcpy.env.overwriteOutput = True
arcpy.env.cellSize = "MAXOF"
arcpy.env.cellSizeProjectionMethod = "CONVERT_UNITS"

# set directories
wdir = "D:\\BioGeographicZones\\"  # working
tdir = wdir + "tmp\\"  # temp
export_shp = wdir + "BGZ_Buffer" + str(buffer_km) + "km_WGS84_Shp\\"  # shp export
export_wkt = wdir + "BGZ_Buffer" + str(buffer_km) + "km_WGS84_Points\\"  # wkt export
# create or clear directories
for cdir in [tdir, export_shp, export_wkt]:
    if not os.path.exists(cdir):
        os.makedirs(cdir)
    else:
        for f in os.listdir(cdir):
            os.remove(cdir + f)
# set scratch
arcpy.env.scratchWorkspace = tdir

# raw BGZs
bgzs_bng = wdir + "BGZ_Phase4\\LivingEngland_BioGeographicZones_Phase4_DRAFT.shp"

print("> projecting from BNG to WGS84...")
# project from BNG to WGS84
bgzs_wgs84 = tdir + "bgzs_wgs84.shp"
wgs84 = arcpy.SpatialReference(4326)  # WGS 1984: 4326 | WGS 1984 Web Mercator: 3857
arcpy.management.Project(bgzs_bng, bgzs_wgs84, wgs84)

print("> applying " + str(buffer_km) + " km buffer...")
# buffer BGZs
bgzs_buffer = tdir + "bgzs_buffer.shp"
arcpy.analysis.Buffer(bgzs_wgs84, bgzs_buffer, str(buffer_km) + " Kilometers", "FULL", "ROUND", "NONE")

print("> exporting BGZ shps...")
# export individual shapefiles
arcpy.analysis.SplitByAttributes(bgzs_buffer, export_shp, "Zone")

print("> calculating coordinates of vertices...")
# convert polygons to vertices
bgzs_points = tdir + "bgzs_points.shp"
arcpy.management.FeatureVerticesToPoints(bgzs_buffer, bgzs_points, "ALL")
# add xy
arcpy.management.AddXY(bgzs_points)

print("> exporting BGZ csvs...")
# keep only necessary fields
arcpy.management.DeleteField(bgzs_points, ["Zone", "POINT_X", "POINT_Y"], "KEEP_FIELDS")
# loop for bgzs
for bgz_short in range(1, 15):
    bgz = "{:02d}".format(bgz_short)
    # split & export csv
    arcpy.conversion.TableToTable(bgzs_points, tdir, "BGZ" + str(bgz) + ".csv", "Zone = " + str(bgz))

print("> format wkt polygons...")
# loop for bgzs
for bgz_short in range(1, 15):
    bgz = "{:02d}".format(bgz_short)
    # initialise WKT file
    wkt_file = export_wkt + "BGZ" + bgz + "_wkt.txt"
    with open(wkt_file, "w") as f:
        f.write("POLYGON((")
    # read csv
    with open(tdir + "BGZ" + bgz + ".csv", newline='') as csvfile:
        csvdata = csv.reader(csvfile)
        # extract data and write to WKT file
        i = 0
        # try:
        for row in csvdata:
            i += 1
            # store first coordinate
            if i == 2:
                start_row = row
            # write coordinates
            if i != 1:
                with open(wkt_file, "a") as f:
                    f.write(row[-2] + " " + row[-1] + ",")
                    # except csv.Error as e:
        #    sys.exit('file {}, line {}: {}'.format(filename,csvdata.line_num,e))
    # finish WKT file
    with open(wkt_file, "a") as f:
        f.write(start_row[-2] + " " + start_row[-1] + "))")

    # end
print("> done")
toc0 = time.time()
print(datetime.timedelta(seconds=round(toc0 - tic0)))
