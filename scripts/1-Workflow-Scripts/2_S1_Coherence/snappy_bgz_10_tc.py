"""-------------------------------------------------------------
    Script Name:   snappy_bgz_10_tc.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 10
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 10.02.23
-------------------------------------------------------------"""

# step 10 (TC): Apply terrain correction and geocoding

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_mrg_ml = sys.argv[1]
dem = sys.argv[2]
dem_path = sys.argv[3]

# read
s1_mrg_ml = sp.ProductIO.readProduct(io_mrg_ml + '.dim')
# set up params
parameters = sp.HashMap()
parameters.put('demName', dem)  # Digital Elevation Model
parameters.put('demResamplingMethod', 'BILINEAR_INTERPOLATION')  # DEM Resampling Method
parameters.put('imgResamplingMethod', 'BILINEAR_INTERPOLATION')  # Image Resampling Method
parameters.put('mapProjection', '27700')  # Map Projection
# EPSG:27700 = OSGB 1936 / British National Grid | WGS84(DD)
parameters.put('nodataValueAtSea', False)  # Mask out areas without elevation
parameters.put('outputComplex', False)  # Output complex data
parameters.put('applyRadiometricNormalization', False)  # Apply radiometric normalization
if dem == 'External DEM':
    parameters.put('externalDEMFile', dem_path)  # External DEM to use
#    parameters.put('externalDEMNoDataValue', 32767)  # Null data value in DEM (must be in double format)
    parameters.put('externalDEMApplyEGM', True)  # Apply earth gravitational model to DEM
# terrain correction
s1_mrg_ml_tc = sp.GPF.createProduct('Terrain-Correction', parameters, s1_mrg_ml)
# write
sp.ProductIO.writeProduct(s1_mrg_ml_tc, io_mrg_ml + '_TC', 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
