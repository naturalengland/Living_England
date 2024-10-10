"""-------------------------------------------------------------
    Script Name:   snappy_bgz_09_ml.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 09
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 9 (ML): Multilook merged image

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_mrg = sys.argv[1]
io_mrg_ml = sys.argv[2]

# read
s1_mrg = sp.ProductIO.readProduct(io_mrg + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('grSquarePixel', True)
parameters.put('nRgLooks', 4)
parameters.put('nAzLooks', 1)
parameters.put('outputIntensity', False)
# multilook coherence image
s1_mrg_ml = sp.GPF.createProduct('Multilook', parameters, s1_mrg)
# write
sp.ProductIO.writeProduct(s1_mrg_ml, io_mrg_ml, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
