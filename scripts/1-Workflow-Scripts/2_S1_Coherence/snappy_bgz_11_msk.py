"""-------------------------------------------------------------
    Script Name:   snappy_bgz_11_msk.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 11
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 28.02.23
-------------------------------------------------------------"""

# step 11 (Msk): Mask and clip coherence map to biogeographic zone

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_mrg_ml_tc = sys.argv[1]
shp = sys.argv[2]
bgz = sys.argv[3]
wkt_polygon = sys.argv[4]

# read
s1_mrg_ml_tc = sp.ProductIO.readProduct(io_mrg_ml_tc + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('vectorFile', shp)
parameters.put('separateShapes', False)
# read vector shapefile
s1_mrg_ml_tc_vec = sp.GPF.createProduct('Import-Vector', parameters, s1_mrg_ml_tc)
# set up params
parameters = sp.HashMap()
parameters.put('geometry', str(int(bgz)))
parameters.put('shorelineExtension', 0)
# apply mask
s1_mrg_ml_tc_msk = sp.GPF.createProduct('Land-Sea-Mask', parameters, s1_mrg_ml_tc_vec)
# write
# sp.ProductIO.writeProduct(s1_mrg_ml_tc_msk,io_mrg_ml_tc+'_msk','BEAM-DIMAP')

# read wkt aoi
f = open(wkt_polygon)
aoi = f.read()
f.close()
geometry = sp.WKTReader().read(aoi)
# set up params
parameters = sp.HashMap()
parameters.put('copyMetadata', True)
parameters.put('geoRegion', geometry)
# subset
s1_mrg_ml_tc_sub = sp.GPF.createProduct('Subset', parameters, s1_mrg_ml_tc_msk)
# write
sp.ProductIO.writeProduct(s1_mrg_ml_tc_sub, io_mrg_ml_tc + '_msk', 'BEAM-DIMAP')
# write geotiff
# sp.ProductIO.writeProduct(s1_mrg_ml_tc_sub, io_mrg_ml_tc + '_msk', 'GeoTIFF')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
