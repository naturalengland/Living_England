"""-------------------------------------------------------------
    Script Name:   snappy_bgz_11-2_sub.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 11-2
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 11-2 (Sub): Split coherence map to reduce data size

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_mrg_ml_tc_msk = sys.argv[1]
io_mrg_ml_tc_sub = sys.argv[2]
wkt_polygon = sys.argv[3]

# read wkt aoi
f = open(wkt_polygon)
aoi = f.read()
f.close()
geometry = sp.WKTReader().read(aoi)

# read
s1_mrg_ml_tc_msk = sp.ProductIO.readProduct(io_mrg_ml_tc_msk + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('copyMetadata', True)
parameters.put('geoRegion', geometry)
# subset
s1_mrg_ml_tc_sub = sp.GPF.createProduct('Subset', parameters, s1_mrg_ml_tc_msk)
# write
sp.ProductIO.writeProduct(s1_mrg_ml_tc_sub, io_mrg_ml_tc_sub, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
