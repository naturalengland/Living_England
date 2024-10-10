"""-------------------------------------------------------------
    Script Name:   snappy_bgz_07_deb.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 07
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 7 (Deb): Deburst coherence image

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_cc = sys.argv[1]
io_cc_deb = sys.argv[2]

# read
s1_cc = sp.ProductIO.readProduct(io_cc + '.dim')
# set up params
parameters = sp.HashMap()
# deburst coherence maps
s1_cc_deb = sp.GPF.createProduct('TOPSAR-Deburst', parameters, s1_cc)
# write
sp.ProductIO.writeProduct(s1_cc_deb, io_cc_deb, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
