"""-------------------------------------------------------------
    Script Name:   snappy_bgz_mosaic.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - mosaic average coherence maps
    Created By:    Chris Moore
    Date Created:  24.09.21
    Date Modified: 26.07.23
-------------------------------------------------------------"""

import time
import datetime
import gc
import snappy as sp

# biogeographic zones and asc/desc to mosaic: ['01'-'14', ['A'/'D'/'A2'/'D2'/'A14'/'D14']] (use A14/D14 to get bgz 1314)
bgzs_ads = [['10', ['D', 'D2']], ['12', ['A', 'A2']], ['12', ['D', 'D2']], ['13', ['A', 'A14']], ['13', ['D', 'D2', 'D14']]]
# season: 'sssYYYY' (sss = spr/sum/aut)
season = 'aut2022'

# print script info
print('> snappy_bgz_mosaic')
# start timer
tic = time.time()
# enable garbage collection
gc.enable()
gc.collect()

# loop for bgzs
for bgz_ads in bgzs_ads:
    # get bgz and asc/desc
    bgz = bgz_ads[0]
    ads = bgz_ads[1]
    print('> averaging bgz ' + str(bgz) + ': ' + ', '.join(ads))

    # set processing directory
    wdir = 'D:\\S1_coherence\\phase5\\data\\bgz' + bgz
    # set output base file name
    out_base = 'bgz' + bgz + '_' + season

    # loop for polarisations
    for pol in ['VH', 'VV']:

        # get coherence maps
        s1_tiffs = []
        for ad in ads:
            if ad == 'D14' or ad == 'A14':
                # read tiff from bgz 14
                s1_tiff = sp.ProductIO.readProduct('D:\\S1_coherence\\phase5\\data\\bgz14\\' + season + '_' + ad[0] +
                                                   '\\bgz14_' + season + '_' + ad[0] + '_' + pol + '_cc.tif')
                # adjust output base file name
                out_base = 'bgz1314_' + season
            else:
                # set data directory
                ddir = wdir + '\\' + season + '_' + ad
                # read tiff
                s1_tiff = sp.ProductIO.readProduct(ddir + '\\bgz' + bgz + '_' + season + '_' + ad + '_' + pol + '_cc.tif')
            # store tiffs
            s1_tiffs.append(s1_tiff)

        # set up params
        parameters = sp.HashMap()
        parameters.put('resamplingMethod', 'BILINEAR_INTERPOLATION')
        parameters.put('average', True)
        parameters.put('normalizeByMean', False)
        # mosaic
        s1_mean = sp.GPF.createProduct('SAR-Mosaic', parameters, s1_tiffs)

        # write as geotiff
        sp.ProductIO.writeProduct(s1_mean, wdir + '\\' + out_base + '_' + ads[0][0] + '_' + pol + '_cc.tif',
                                  'GeoTIFF-BigTIFF')
        print('> processed ' + out_base + '_' + ads[0][0] + '_' + pol + '_cc.tif')

# end
toc = time.time()
print('> done\n> time elapsed: ' + str(datetime.timedelta(seconds=round(toc - tic))))
