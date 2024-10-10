"""-------------------------------------------------------------
    Script Name:   snappy_bgz_13_bm.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 13
    Created By:    Chris Moore
    Date Created:  11.08.21
    Date Modified: 28.07.22
-------------------------------------------------------------"""

# step 13 (BM): Calculate average temporal coherence

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_stack = sys.argv[1]
pol = sys.argv[2]
min_btemp = sys.argv[3]
max_btemp = sys.argv[4]

# define plus string
plus = ' + '
# get ifg temporal baselines
btemps = range(int(min_btemp), int(max_btemp) + 1, int(min_btemp))
# store new bands
targetBands = sp.jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', len(btemps))
n = 0
# initialise new bands
sp.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
BandDescriptor = sp.jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')
# read
s1_stack = sp.ProductIO.readProduct(io_stack + '.dim')
band_names = list(s1_stack.getBandNames())

# loop through ifg temporal baselines
for btemp in btemps:
    # get bands to use
    band_calc = []
    for name in band_names:
        # get time difference between dates
        date1 = datetime.date(int(name[12:16]), time.strptime(name[9:12], '%b').tm_mon, int(name[7:9]))
        date2 = datetime.date(int(name[22:26]), time.strptime(name[19:22], '%b').tm_mon, int(name[17:19]))
        diff = date2 - date1
        # store ifgs
        if diff.days == btemp and name[4:6] == pol:
            band_calc.append(name)
    # define new band
    cc_mean = BandDescriptor()
    cc_mean.name = 'cc_' + pol + '_' + str(btemp) + 'day_mean'
    cc_mean.type = 'float32'
    cc_mean.expression = '(' + plus.join(band_calc) + ') / ' + str(len(band_calc))
    # store new bands
    targetBands[n] = cc_mean
    n += 1

# set up params
parameters = sp.HashMap()
parameters.put('targetBands', targetBands)
# calculate bands
s1_stack_BM = sp.GPF.createProduct('BandMaths', parameters, s1_stack)
# write
sp.ProductIO.writeProduct(s1_stack_BM, io_stack + '_' + pol + '_BM', 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
