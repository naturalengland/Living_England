"""-------------------------------------------------------------
    Script Name:   snappy_bgz_01_asm.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 01
    Created By:    Chris Moore
    Date Created:  04.08.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""
# step 1 (Asm): Assemble along track scenes

import time
import datetime
import sys
import glob
import sys
# sys.path.insert(1, "C:\\Users\\ne.user\\.snap\\snap-python")
import snappy as sp
from snappy import ProductIO

# start timer
tic = time.time()
# get inputs
wdir = sys.argv[1]
d1 = sys.argv[2]

# read
s1_zips = []
for zip_path in glob.glob(wdir + '\\*' + d1 + '*.zip'):
    s1_zip = sp.ProductIO.readProduct(zip_path)
    s1_zips.append(s1_zip)
# skip if have only 1 zip
if len(s1_zips) == 1:
    # write
    sp.ProductIO.writeProduct(s1_zip, wdir + '\\' + d1 + '_SLC', 'BEAM-DIMAP')
else:
    # set up params
    parameters = sp.HashMap()
    # assemble scenes
    s1_slc = sp.GPF.createProduct('SliceAssembly', parameters, s1_zips)
    # write
    sp.ProductIO.writeProduct(s1_slc, wdir + '\\' + d1 + '_SLC', 'BEAM-DIMAP')

# output
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
