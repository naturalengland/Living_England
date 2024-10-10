"""-------------------------------------------------------------
    Script Name:   snappy_bgz_08_mrg.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 08
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 8 (Mrg): Merge subswaths

import time
import datetime
import sys
import os
import snappy as sp

# start timer
tic = time.time()
# get inputs
wdir = sys.argv[1]
ifg = sys.argv[2]
io_mrg = sys.argv[3]

# initialise
merge_in = []
# loop through subswaths
for iw in [1, 2, 3]:
    # get ifg name
    io_cc_deb = wdir + '\\' + ifg + '_IW' + str(iw) + '_cc_deb'
    # if have ifg
    if os.path.exists(io_cc_deb + '.dim'):
        # read
        s1_cc_deb = sp.ProductIO.readProduct(io_cc_deb + '.dim')
        # store
        merge_in.append(s1_cc_deb)

# merge if have 2/3 subswaths
if len(merge_in) > 1:
    # set up params
    parameters = sp.HashMap()
    # merge subswaths
    s1_mrg = sp.GPF.createProduct('TOPSAR-Merge', parameters, merge_in)
else:
    # get band names
    s1_mrg = s1_cc_deb
    band_names = list(s1_mrg.getBandNames())
    # edit band names
    for name in band_names:
        new_name = name[0:3] + name[7::]  # remove "_IW1"
        s1_mrg.getBand(name).setName(new_name)

# write
sp.ProductIO.writeProduct(s1_mrg, io_mrg, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
