"""-------------------------------------------------------------
    Script Name:   snappy_bgz_03_orb.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 03
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 3 (Orb): Apply orbit correction to SLCs

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_slc = sys.argv[1]

# read
s1_slc = sp.ProductIO.readProduct(io_slc + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('orbitType', 'Sentinel Precise (Auto Download)')
parameters.put('polyDegree', 3)
parameters.put('continueOnFail', True)
# apply orbit files
s1_slc_orb = sp.GPF.createProduct('Apply-Orbit-File', parameters, s1_slc)
# write
sp.ProductIO.writeProduct(s1_slc_orb, io_slc + '_Orb', 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
