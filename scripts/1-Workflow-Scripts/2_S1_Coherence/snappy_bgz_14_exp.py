"""-------------------------------------------------------------
    Script Name:   snappy_bgz_14_exp.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 14
    Created By:    Chris Moore
    Date Created:  11.08.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 14 (Exp): Export average coherence as geotiff

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_base = sys.argv[1]
pol = sys.argv[2]

# define plus string
plus = ' + '
# read
s1_stack_BM = sp.ProductIO.readProduct(io_base + '_stack_' + pol + '_BM.dim')
band_names = list(s1_stack_BM.getBandNames())

# get bands to use
band_calc = []
for name in band_names:
    if name[3:5] == pol:
        band_calc.append(name)

# store new bands
targetBands = sp.jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
# initialise new bands
sp.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
BandDescriptor = sp.jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

# define new band
cc_mean = BandDescriptor()
cc_mean.name = 'cc_' + pol + '_all_mean'
cc_mean.type = 'float32'
cc_mean.expression = '(' + plus.join(band_calc) + ') / ' + str(len(band_calc))
# store new bands
targetBands[0] = cc_mean

# set up params
parameters = sp.HashMap()
parameters.put('targetBands', targetBands)
# calculate bands
s1_mean = sp.GPF.createProduct('BandMaths', parameters, s1_stack_BM)
# write as geotiff
sp.ProductIO.writeProduct(s1_mean, io_base + '_' + pol + '_cc', 'GeoTIFF-BigTIFF')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
