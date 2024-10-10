"""-------------------------------------------------------------
    Script Name:   snappy_bgz_02_slc.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 02
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 2 (SLC): Select swaths and bursts

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_slc = sys.argv[1]
io_IWslc = sys.argv[2]
iw = sys.argv[3]
wkt_polygon = sys.argv[4]

# read bgz aoi
f = open(wkt_polygon)
aoi = f.read()
f.close()

# read
s1_slc = sp.ProductIO.readProduct(io_slc + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('subswath', 'IW' + iw)
parameters.put('firstBurstIndex', 1)
parameters.put('lastBurstIndex', 999)
parameters.put('wktAoi', aoi)
# catch error
try:
    # select bursts
    s1_iw = sp.GPF.createProduct('TOPSAR-Split', parameters, s1_slc)
    # write
    sp.ProductIO.writeProduct(s1_iw, io_IWslc, 'BEAM-DIMAP')
except Exception:
    # sys.exc_clear()
    pass

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
