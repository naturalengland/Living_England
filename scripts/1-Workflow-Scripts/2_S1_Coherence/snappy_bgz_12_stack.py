"""-------------------------------------------------------------
    Script Name:   snappy_bgz_12_stack.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 12
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 30.05.23
-------------------------------------------------------------"""

# step 12 (Stack): Stack coherence maps

import time
import datetime
import sys
import os
import glob
import snappy as sp

# start timer
tic = time.time()
# get inputs
wdir = sys.argv[1]
min_btemp = sys.argv[2]
max_btemp = sys.argv[3]
io_stack = sys.argv[4]

# get dates
dates_tmp = []
for s1_zip in glob.glob(wdir + '\\*.zip'):
    timestamp = s1_zip.split("_IW_SLC__1SDV_")[1]
    dates_tmp.append(timestamp[:8])
dates = list(set(dates_tmp))
dates.sort()
# get date pairs
ifgs = []
for d1 in range(len(dates)):
    for d2 in range(d1 + 1, len(dates)):
        ifgs.append(dates[d1] + '_' + dates[d2])

# get dates
dates = [s1_zip.split("_IW_SLC__1SDV_")[1][:8] for s1_zip in glob.glob(wdir + '\\*.zip')]
# sort unique dates
dates = [*set(dates)]
dates.sort()
# get date pairs within temporal baseline limits
ifgs = []
for d1 in range(len(dates)):
    for d2 in range(d1 + 1, len(dates)):
        btemp = abs((datetime.datetime.strptime(dates[d1], "%Y%m%d") -
                     datetime.datetime.strptime(dates[d2], "%Y%m%d")).days)
        if int(min_btemp) <= btemp <= int(max_btemp):
            ifgs.append(dates[d1] + '_' + dates[d2])

# initialise
stack_in = []
# loop through ifgs
for ifg in ifgs:
    # get ifg name
    io_cc_mrg_ml_tc_msk = wdir + '\\' + ifg + '_cc_ML_TC_msk'
    # if have ifg
    if os.path.exists(io_cc_mrg_ml_tc_msk + '.dim'):
        # read
        s1_cc_mrg_ml_tc_msk = sp.ProductIO.readProduct(io_cc_mrg_ml_tc_msk + '.dim')
        # store
        stack_in.append(s1_cc_mrg_ml_tc_msk)

# set up params
parameters = sp.HashMap()
# parameters.put('masterBandNames', stack_in[0].getBandNames()[0])
parameters.put('resamplingType', 'NEAREST_NEIGHBOUR')  # None, 'NEAREST_NEIGHBOUR', 'BILINEAR_INTERPOLATION'
parameters.put('initialOffsetMethod', 'Product Geolocation')  # 'Orbit', 'Product Geolocation'
parameters.put('extent', 'Master')
# stack ifgs
s1_stack = sp.GPF.createProduct('CreateStack', parameters, stack_in)
# write
sp.ProductIO.writeProduct(s1_stack, io_stack, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
