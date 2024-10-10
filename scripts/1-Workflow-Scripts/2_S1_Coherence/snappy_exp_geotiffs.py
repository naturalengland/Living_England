"""-------------------------------------------------------------
    Script Name:   snappy_exp_geotiffs.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - export ifg coherence maps
    Created By:    Chris Moore
    Date Created:  30.05.23
    Date Modified: 01.06.23
-------------------------------------------------------------"""

import time
import datetime
import glob
import snappy as sp

# biogeographic zone
bgzs = ['06', '13']  # '01'-'14'
# season
season = 'aut2022'  # 'sssYYYY' (sss = spr/sum/aut)
# asc/desc
ads = ['A', 'D']  # 'A'/'D'

print('> snappy_exp_geotiffs')
# start timer
tic = time.time()
# minimum and maximum temporal baselines that have been processed
min_btemp = 12  # days
max_btemp = 12  # days
# master directory containing processing and data folders
mdir = 'D:\\S1_coherence\\phase5\\data'

# try:
# begin loop for bgzs and asc/desc
for bgz in bgzs:
    for ad in ads:
        print('> exporting geotiffs for bgz' + bgz + ', ' + season + ', ' + ad + '...')

        # set processing directory
        wdir = mdir + '\\bgz' + bgz + '\\' + season + '_' + ad

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
                if min_btemp <= btemp <= max_btemp:
                    ifgs.append(dates[d1] + '_' + dates[d2])

        # loop through ifgs
        for ifg in ifgs:
            # read ifg from BEAM-DIMAP
            ifg_in = wdir + '\\' + ifg + '_cc_ML_TC_msk'
            s1_in = sp.ProductIO.readProduct(ifg_in + '.dim')
            # write to geotiff
            sp.ProductIO.writeProduct(s1_in, ifg_in, 'GeoTIFF')

# # catch errors
# except Exception as e:
#     print(e)
# end
toc = time.time()
print('> done\n> time elapsed: ' + str(datetime.timedelta(seconds=round(toc - tic))))
