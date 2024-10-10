"""-------------------------------------------------------------
    Script Name:   snappy_mv_outputs.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - move coherence maps to output folder
    Created By:    Chris Moore
    Date Created:  25.07.23
    Date Modified: 25.07.23
-------------------------------------------------------------"""

import os
import shutil
import time
import datetime

# season
season = 'aut2022'  # 'sssYYYY' (sss = spr/sum/aut)

# print script info
print('> snappy_mv_outputs')
# start timer
tic = time.time()
# main directory
mdir = 'D:\\S1_coherence\\phase5\\data'
# create output directory
wdir = mdir + '\\output'
if not os.path.exists(wdir):
    os.mkdir(wdir)

# biogeographic zones
bgzs = [str(i).zfill(2) for i in range(1, 13)] + ['1314']
# loop for bgzs, orbits & polarisations
for bgz in bgzs:
    for ad in ['A', 'D']:
        for pol in ['VH', 'VV']:
            # get data directory
            ddir = mdir + '\\bgz' + bgz + '\\' + season + '_' + ad
            # get file name
            io_cc = 'bgz' + bgz + '_' + season + '_' + ad + '_' + pol + '_cc.tif'
            # copy file
            print('> copying ' + io_cc)
            shutil.copyfile(ddir + '\\' + io_cc, wdir + '\\' + io_cc)

# end
toc = time.time()
print('> done\n> time elapsed: ' + str(datetime.timedelta(seconds=round(toc - tic))))
