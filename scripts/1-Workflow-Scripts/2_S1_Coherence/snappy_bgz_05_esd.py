"""-------------------------------------------------------------
    Script Name:   snappy_bgz_05_esd.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 05
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

# step 5 (Esd): Refine cogregistration

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_rslc = sys.argv[1]
io_rslc_esd = sys.argv[2]

# read
s1_rslc = sp.ProductIO.readProduct(io_rslc + '.dim')
# set up params
parameters = sp.HashMap()
# parameters.put('fineWinWidthStr',512) # Registration Window Width
# parameters.put('fineWinHeightStr',512) # Registration Window Height
# parameters.put('fineWinAccAzimuth',16) # Search Window Accuracy in Azimuth Diection
# parameters.put('fineWinAccRange',16) # Search Window Accuracy in Range Diection
# parameters.put('fineWinOversampling',128) # Window oversmapling factor
# parameters.put('xCorrThreshold',0.1) # Cross-Correlation Threshold
# parameters.put('cohThreshold',0.3) # Coherence Threshold for Outlier Removal
# parameters.put('numBlocksPerOverlap',10) # Number of Windows Per Overlap for ESD
# parameters.put('esdEstimator','Periodogram') # ESD Estimator
# parameters.put('weightFunc','Inv Quadratic') # Weight function
# parameters.put('temporalBaselineType','Number of images') # Temporal baseline type
# parameters.put('maxTemporalBaseline',4) # Maximum temporal baseline
# parameters.put('integrationMethod','L1 and L2') # Integration method
# parameters.put('doNotWriteTargetBands',False) # Do not write target bands
# parameters.put('useSuppliedRangeShift',False) # Use user supplied range shift
# parameters.put('useSuppliedAziumuthShift',False) # Use user supplied aziumuth shift
# enhanced spectral diversity
s1_rslc_esd = sp.GPF.createProduct('Enhanced-Spectral-Diversity', parameters, s1_rslc)
# write
sp.ProductIO.writeProduct(s1_rslc_esd, io_rslc_esd, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
