"""-------------------------------------------------------------
    Script Name:   snappy_bgz_04_rslc.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 04
    Created By:    Chris Moore
    Date Created:  18.02.21
    Date Modified: 10.02.23
-------------------------------------------------------------"""

# step 4 (RSLC): Coregister SLC pairs

import time
import datetime
import sys
import snappy as sp

# start timer
tic = time.time()
# get inputs
io_slc1_orb = sys.argv[1]
io_slc2_orb = sys.argv[2]
io_rslc = sys.argv[3]
dem = sys.argv[4]
dem_path = sys.argv[5]

# read
s1_slc1_orb = sp.ProductIO.readProduct(io_slc1_orb + '.dim')
s1_slc2_orb = sp.ProductIO.readProduct(io_slc2_orb + '.dim')
rslc_in = [s1_slc2_orb, s1_slc1_orb]  # requires reverse order
# set up params
parameters = sp.HashMap()
parameters.put('demName', dem)  # Digital Elevation Model
parameters.put('demResamplingMethod', 'BICUBIC_INTERPOLATION')  # DEM Resampling Method
parameters.put('resamplingType', 'BISINC_5_POINT_INTERPOLATION')  # Resampling Type (BILINEAR_INTERPOLATION,
# BICUBIC_INTERPOLATION, BISINC_5_POINT_INTERPOLATION, BISINC_21_POINT_INTERPOLATION)
parameters.put('maskOutAreaWithoutElevation', False)  # Mask out areas with no elevation
if dem == 'External DEM':
    parameters.put('externalDEMFile', dem_path)  # External DEM to use
#    parameters.put('externalDEMNoDataValue', 32767)  # Null data value in DEM (must be in double format)
    parameters.put('externalDEMApplyEGM', True)  # Apply earth gravitational model to DEM
# back geocode slcs
s1_rslc = sp.GPF.createProduct('Back-Geocoding', parameters, rslc_in)
# write
sp.ProductIO.writeProduct(s1_rslc, io_rslc, 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
