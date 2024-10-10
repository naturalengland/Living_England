"""-------------------------------------------------------------
    Script Name:   snappy_bgz_06_cc.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 06
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 6 (Coh): Calculate interferogram coherence

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_rslc_esd = sys.argv[1]
io_cc = sys.argv[2]

# read
s1_rslc_esd = sp.ProductIO.readProduct(io_rslc_esd + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('cohWinRg', 10)  # Coherence Range Window Size
parameters.put('cohWinAz', 2)  # Coherence Azimuth Window Size
parameters.put('subtractFlatEarthPhase', False)  # Subtract flat-earth phase
parameters.put('srpPolynomialDegree', 5)  # Degree of Flat Earth polynomial
parameters.put('srpNumberPoints', 501)  # Number of Flat Earth estimation points
parameters.put('squarePixel', False)  # Independent Window Size
parameters.put('subtractTopographicPhase', False)  # Subtract topographic phase
parameters.put('singleMaster', True)  # Single Master
# calculate ifg coherence
s1_cc = sp.GPF.createProduct('Coherence', parameters, s1_rslc_esd)
# write
sp.ProductIO.writeProduct(s1_cc, io_cc, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
