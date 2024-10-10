"""-------------------------------------------------------------
    Script Name:   snappy_bgz.py
    Description:   SNAP processing chain for VV/VH S1 coherence data
    Created By:    Chris Moore
    Date Created:  15.01.21
    Date Modified: 10.08.23
-------------------------------------------------------------"""

import sys
import os
import gc
import glob
import datetime
import subprocess

# biogeographic zone
bgzs = ['09', '10', '11', '12', '13', '14']  # '01'-'14'
# bgzs_ads = [['07', ['D']], ['10', ['D2']], ['12', ['A2', 'D2']], ['13', ['D2']]]  # tmp
# season
season = 'spr2023'  # 'sprYYYY'/'sumYYYY'/'autYYYY'
# asc/desc
ads = ['A', 'D']  # 'A'/'D'
# start & end steps
start = 0
end = 14
# clean step outputs? (step 11 excluded for qc)
step_clean = True  # True/False
end_clean = False  # True/False

# minimum and maximum temporal baseline to process
min_btemp = 12  # days
max_btemp = 24  # days
# master directory containing processing and data folders
mdir = 'D:\\S1_coherence\\phase5\\data'
# directory containing running scripts
sdir = 'D:\\S1_coherence\\phase5\\scripts'
# directory containing BGZ WKT polygons
wktdir = 'D:\\S1_coherence\\BGZ_Buffer10km_WGS84_Points'
# directory containing BGZ ESRI shapefiles
shpdir = 'D:\\S1_coherence\\BGZ_Buffer1km_WGS84_Shp'
# dem to use in steps 4 & 10 (External DEM, SRTM 1Sec HGT, SRTM 3Sec, Copernicus 30m Global DEM, Copernicus 90m Global DEM)
dem = 'External DEM'
# path to dem if using external DEM
dem_path = 'D:\\DTM\\mergedDEM_LiDAR_IHM_SRTM.tif'

# step 0 (Zip): Download/copy zip files to directory
# step 1 (Asm): Assemble along-track scenes
# step 2 (SLC): Select swaths and bursts
# step 3 (Orb): Apply orbit correction
# step 4 (RSLC): Coregister SLC pairs
# step 5 (Esd): Refine coregistration
# step 6 (Coh): Calculate interferogram coherence
# step 7 (Deb): De-burst coherence image
# step 8 (Mrg): Merge sub-swaths
# step 9 (ML): Multi-look merged image
# step 10 (TC): Apply terrain correction and geocoding
# step 11 (Msk): Mask and clip to biogeographic zone
# step 12 (Stack): Stack coherence maps
# step 13 (BM): Calculate temporal average coherence
# step 14 (Exp): Export average coherence as geotiff
# END (Clean): Remove temporary processing files

# -------------------- setup -------------------- #

# start timer
tic = datetime.datetime.now().replace(microsecond=0)
# enable garbage collection
gc.enable()
gc.collect()

# script info
print('> snappy_bgz')
run_steps = range(start, end + 1)
# run_steps = []
print('> running steps: ' + ', '.join([str(s) for s in run_steps]))

# begin loop for bgzs and asc/desc
for bgz in bgzs:
    # for bgz_ads in bgzs_ads:  # tmp
    # bgz = bgz_ads[0]  # tmp
    # ads = bgz_ads[1]  # tmp
    for ad in ads:
        # catch errors
        try:
            # start processing chain
            now = datetime.datetime.now()
            print('> [' + now.strftime("%d/%m %H:%M") + '] START: bgz' + bgz + ', ' + season + ', ' + ad)
            # set processing directory
            wdir = mdir + '\\bgz' + bgz + '\\' + season + '_' + ad
            if not os.path.exists(wdir):
                os.mkdir(wdir)
            # get basename for coherence outputs
            io_base = wdir + '\\bgz' + bgz + '_' + season + '_' + ad

            # --------------------- step 0 --------------------- #

            if 0 in run_steps:
                # read urls
                url_file = mdir + '\\bgz' + bgz + '\\bgz' + bgz + '_' + season + '_' + ad + '_urls.txt'
                with open(url_file, 'r') as f:
                    url_list = f.read().splitlines()

                # screen dates for weather effects
                print('> 0 (QC): checking acquisition dates...')
                # create api download directory
                if not os.path.exists(wdir + "\\era5_api"):
                    os.mkdir(wdir + "\\era5_api")
                # get biogeographic zone polygon
                csv_polygon = wktdir + '\\BGZ' + bgz + '.csv'
                # get unique acquisition dates with times (YYYYMMDDhhmmss)
                s1_dates = []
                for u in range(len(url_list)):
                    if url_list[u][56:64] != url_list[u - 1][56:64]:
                        s1_dates.append(url_list[u][56:64] + url_list[u][65:71])
                # get era-5 weather data for 6 days prior to date
                s1_dates_flag = []
                for s1_date in s1_dates:
                    # download era-5 data and apply thresholds
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_00_qc.py', s1_date, csv_polygon, bgz,
                                                        wdir + "\\era5_api"], stderr=subprocess.STDOUT)
                    print('> 0 (QC): ' + s1_date[:8] + ' check ' + pipe_out.decode(sys.stdout.encoding).split()[-3] +
                          ' in ' + pipe_out.decode(sys.stdout.encoding).split()[-2])
                    # store flagged dates
                    if pipe_out.decode(sys.stdout.encoding).split()[-3] == "passed":
                        s1_dates_flag.append(True)
                    else:
                        s1_dates_flag.append(False)
                # filter urls using flags
                url_list_qc = []
                for url_name in url_list:
                    url_i = [d[:8] for d in s1_dates].index(url_name[56:64])
                    if s1_dates_flag[url_i]:
                        url_list_qc.append(url_name)

                # crop url list if have too many dates
                # get dates
                dates_qc = [url_name.split("_IW_SLC__1SDV_")[1][:8] for url_name in url_list_qc]
                # sort unique dates
                dates_qc = [*set(dates_qc)]
                dates_qc.sort()
                # get date pairs within temporal baseline limits
                ifgs_qc = []
                btemp_qc = []
                for d1 in range(len(dates_qc)):
                    for d2 in range(d1 + 1, len(dates_qc)):
                        btemp = abs((datetime.datetime.strptime(dates_qc[d1], "%Y%m%d") -
                                     datetime.datetime.strptime(dates_qc[d2], "%Y%m%d")).days)
                        if min_btemp <= btemp <= max_btemp:
                            ifgs_qc.append(dates_qc[d1] + '_' + dates_qc[d2])
                            btemp_qc.append(btemp)
                # check number of ifgs for each btemp
                crop_flag = True
                for b in [*range(min_btemp, max_btemp+1, min_btemp)]:
                    if btemp_qc.count(b) <= 3:
                        crop_flag = False
                # only crop if all btemps have > 3 ifgs
                if crop_flag:
                    print('> 0 (Crop): reducing number of dates to process...')
                    n = 0
                    while crop_flag:
                        # remove the furthest away date (alternating last and first element of list)
                        n += 1
                        date_pop = dates_qc.pop((n % 2) * -1)
                        print('> 0 (Crop): removed ' + date_pop)
                        # get date pairs within temporal baseline limits
                        ifgs_qc = []
                        btemp_qc = []
                        for d1 in range(len(dates_qc)):
                            for d2 in range(d1 + 1, len(dates_qc)):
                                btemp = abs((datetime.datetime.strptime(dates_qc[d1], "%Y%m%d") -
                                             datetime.datetime.strptime(dates_qc[d2], "%Y%m%d")).days)
                                if min_btemp <= btemp <= max_btemp:
                                    ifgs_qc.append(dates_qc[d1] + '_' + dates_qc[d2])
                                    btemp_qc.append(btemp)
                        # continue crop if all btemps have > 3 ifgs
                        for b in [*range(min_btemp, max_btemp+1, min_btemp)]:
                            if btemp_qc.count(b) <= 3:
                                crop_flag = False
                    # adjust url list
                    url_list_qc_crop = []
                    for url_name in url_list_qc:
                        if url_name.split("_IW_SLC__1SDV_")[1][:8] in dates_qc:
                            url_list_qc_crop.append(url_name)
                else:
                    url_list_qc_crop = url_list_qc

                # download data from asf
                print('> 0 (Zip): downloading zip files...')
                # check if files exist and are not empty (zip size > 1kb)
                for url_name in url_list_qc_crop:
                    if os.path.exists(wdir + '\\' + url_name[39:]) and os.path.getsize(wdir+'\\'+url_name[39:]) > 1024:
                        print('> 0 (Zip): skipping ' + url_name[39:])
                    else:
                        # if exists delete empty zip file
                        if os.path.exists(wdir + '\\' + url_name[39:]):
                            os.remove(wdir + '\\' + url_name[39:])
                        # download data (use python 3 conda env with asf_search module installed)
                        pipe_out = subprocess.check_output(
                            ['C:\\Users\\ne.user\\miniconda3\\envs\\asf_download\\python.exe',
                             sdir + '\\snappy_bgz_00_zip\\snappy_bgz_00_zip.py', url_name, wdir],
                            stderr=subprocess.STDOUT)
                        print('> 0 (Zip): downloaded ' + url_name[39:] + ' in ' +
                              pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- setup --------------------- #

            # get dates
            dates = [s1_zip.split("_IW_SLC__1SDV_")[1][:8] for s1_zip in glob.glob(wdir + '\\*.zip')]
            # sort unique dates
            dates = [*set(dates)]
            dates.sort()
            # convert date format (20210810 -> 10Aug2021)
            dates2 = [d1[6:8] + datetime.date(1900, int(d1[4:6]), 1).strftime('%b') + d1[0:4] for d1 in dates]
            # get date pairs within temporal baseline limits
            ifgs = []
            ifgs2 = []
            for d1 in range(len(dates)):
                for d2 in range(d1 + 1, len(dates)):
                    btemp = abs((datetime.datetime.strptime(dates[d1], "%Y%m%d") -
                                 datetime.datetime.strptime(dates[d2], "%Y%m%d")).days)
                    if min_btemp <= btemp <= max_btemp:
                        ifgs.append(dates[d1] + '_' + dates[d2])
                        ifgs2.append(dates2[d1] + '_' + dates2[d2])

            # --------------------- step 1 --------------------- #

            if 1 in run_steps:
                print('> 1 (Asm): assembling along-track scenes...')
                # loop through dates
                for d1 in dates:
                    # slice assembly
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_01_asm.py', wdir, d1],
                                                       stderr=subprocess.STDOUT)
                    print('> 1 (Asm): processed ' + d1 + '_SLC in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir,
                                                        str(min_btemp), str(max_btemp), io_base, '0'],
                                                       stderr=subprocess.STDOUT)
                    print('> 1 (Asm): deleted step 0 (Zip) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 2 --------------------- #

            if 2 in run_steps:
                print('> 2 (SLC): extracting sub-swath bursts...')
                # get aoi
                wkt_polygon = wktdir + '\\BGZ' + bgz + '_wkt.txt'
                # loop through dates
                for d1 in dates:
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get slc
                        io_slc = wdir + '\\' + d1 + '_SLC'
                        io_IWslc = wdir + '\\' + d1 + '_IW' + str(iw) + '_SLC'
                        # topsar split
                        pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_02_slc.py', io_slc, io_IWslc,
                                                            str(iw), wkt_polygon], stderr=subprocess.STDOUT)
                        if os.path.exists(io_IWslc + '.dim'):
                            print('> 2 (SLC): processed ' + d1 + '_IW' + str(iw) + '_SLC in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir,
                                                        str(min_btemp), str(max_btemp), io_base, '1'],
                                                       stderr=subprocess.STDOUT)
                    print('> 2 (SLC): deleted step 1 (Asm) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 3 --------------------- #

            if 3 in run_steps:
                print('> 3 (Orb): applying precise orbits...')
                # loop through dates
                for d1 in dates:
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get slc
                        io_IWslc = wdir + '\\' + d1 + '_IW' + str(iw) + '_SLC'
                        # if have slc
                        if os.path.exists(io_IWslc + '.dim'):
                            # refine orbits
                            pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_03_orb.py', io_IWslc],
                                                               stderr=subprocess.STDOUT)
                            print('> 3 (Orb): processed ' + d1 + '_IW' + str(iw) + '_SLC_Orb in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '2'], stderr=subprocess.STDOUT)
                    print('> 3 (Orb): deleted step 2 (SLC) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 4 --------------------- #

            if 4 in run_steps:
                print('> 4 (RSLC): co-registering...')
                # loop through ifgs
                for ifg in ifgs:
                    # get ifg dates
                    d1 = ifg[0:8]
                    d2 = ifg[9:19]
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get slcs
                        io_slc1_orb = wdir + '\\' + d1 + '_IW' + str(iw) + '_SLC_Orb'
                        io_slc2_orb = wdir + '\\' + d2 + '_IW' + str(iw) + '_SLC_Orb'
                        # if have slc pair
                        if os.path.exists(io_slc1_orb + '.dim') and os.path.exists(io_slc2_orb + '.dim'):
                            # get output
                            io_rslc = wdir + '\\' + ifg + '_IW' + str(iw) + '_RSLC'
                            # coregister
                            pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_04_rslc.py', io_slc1_orb,
                                                                io_slc2_orb, io_rslc, dem, dem_path],
                                                               stderr=subprocess.STDOUT)
                            print('> 4 (RSLC): processed ' + ifg + '_IW' + str(iw) + '_RSLC in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '3'], stderr=subprocess.STDOUT)
                    print('> 4 (RSLC): deleted step 3 (Orb) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 5 --------------------- #

            if 5 in run_steps:
                print('> 5 (Esd): refining coregistration...')
                # loop through ifgs
                for ifg in ifgs:
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get rslc
                        io_rslc = wdir + '\\' + ifg + '_IW' + str(iw) + '_RSLC'
                        # if have rslc
                        if os.path.exists(io_rslc + '.dim'):
                            # get output
                            io_rslc_esd = io_rslc + '_esd'
                            # esd
                            pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_05_esd.py', io_rslc,
                                                                io_rslc_esd], stderr=subprocess.STDOUT)
                            print('> 5 (Esd): processed ' + ifg + '_IW' + str(iw) + '_RSLC_esd in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '4'], stderr=subprocess.STDOUT)
                    print('> 5 (Esd): deleted step 4 (RSLC) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 6 --------------------- #

            if 6 in run_steps:
                print('> 6 (Coh): calculating coherence...')
                # loop through ifgs
                for ifg in ifgs:
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get rslc
                        io_rslc_esd = wdir + '\\' + ifg + '_IW' + str(iw) + '_RSLC_esd'
                        # if have rslc
                        if os.path.exists(io_rslc_esd + '.dim'):
                            # get output
                            io_cc = wdir + '\\' + ifg + '_IW' + str(iw) + '_cc'
                            # coherence
                            pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_06_cc.py', io_rslc_esd,
                                                                io_cc], stderr=subprocess.STDOUT)
                            print('> 6 (Coh): processed ' + ifg + '_IW' + str(iw) + '_cc in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '5'], stderr=subprocess.STDOUT)
                    print('> 6 (Coh): deleted step 5 (Esd) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 7 --------------------- #

            if 7 in run_steps:
                print('> 7 (Deb): de-bursting...')
                # loop through ifgs
                for ifg in ifgs:
                    # loop through sub-swaths
                    for iw in [1, 2, 3]:
                        # get cc
                        io_cc = wdir + '\\' + ifg + '_IW' + str(iw) + '_cc'
                        # if have rslc
                        if os.path.exists(io_cc + '.dim'):
                            # get output
                            io_cc_deb = io_cc + '_deb'
                            # de-burst
                            pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_07_deb.py', io_cc, io_cc_deb],
                                                               stderr=subprocess.STDOUT)
                            print('> 7 (Deb): processed ' + ifg + '_IW' + str(iw) + '_cc_deb in ' +
                                  pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '6'], stderr=subprocess.STDOUT)
                    print('> 7 (Deb): deleted step 6 (Coh) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 8 --------------------- #

            if 8 in run_steps:
                print('> 8 (Mrg): merging sub-swaths...')
                # loop through ifgs
                for ifg in ifgs:
                    # get output
                    io_mrg = wdir + '\\' + ifg + '_cc'
                    # merge
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_08_mrg.py', wdir, ifg, io_mrg],
                                                       stderr=subprocess.STDOUT)
                    print('> 8 (Mrg): processed ' + ifg + '_cc in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '7'], stderr=subprocess.STDOUT)
                    print('> 8 (Mrg): deleted step 7 (Deb) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 9 --------------------- #

            if 9 in run_steps:
                print('> 9 (ML): multi-looking...')
                # loop through ifgs
                for ifg in ifgs:
                    # get merged cc
                    io_mrg = wdir + '\\' + ifg + '_cc'
                    # if have cc
                    if os.path.exists(io_mrg + '.dim'):
                        # get output
                        io_mrg_ml = io_mrg + '_ML'
                        # multi-look
                        pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_09_ml.py', io_mrg, io_mrg_ml],
                                                           stderr=subprocess.STDOUT)
                        print('> 9 (ML): processed ' + ifg + '_cc_ML in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '8'], stderr=subprocess.STDOUT)
                    print('> 9 (ML): deleted step 8 (Mrg) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 10 --------------------- #

            if 10 in run_steps:
                print('> 10 (TC): applying terrain correction...')
                # loop through ifgs
                for ifg in ifgs:
                    # get merged cc
                    io_mrg_ml = wdir + '\\' + ifg + '_cc_ML'
                    # if have cc
                    if os.path.exists(io_mrg_ml + '.dim'):
                        # geocode
                        pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_10_tc.py', io_mrg_ml, dem,
                                                            dem_path], stderr=subprocess.STDOUT)
                        print('> 10 (TC): processed ' + ifg + '_cc_ML_TC in ' +
                              pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '9'], stderr=subprocess.STDOUT)
                    print('> 10 (TC): deleted step 9 (ML) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 11 --------------------- #

            if 11 in run_steps:
                print('> 11 (Msk): masking and clipping to biogeographic zone...')
                # get biogeographic zone shapefile and wkt polygon
                shp = shpdir + '\\' + str(int(bgz)) + '.shp'
                wkt_polygon = wktdir + '\\BGZ' + bgz + '_wkt.txt'
                # loop through ifgs
                for ifg in ifgs:
                    # get geocoded cc
                    io_mrg_ml_tc = wdir + '\\' + ifg + '_cc_ML_TC'
                    # if have cc
                    if os.path.exists(io_mrg_ml_tc + '.dim'):
                        # mask
                        pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_11_msk.py', io_mrg_ml_tc, shp,
                                                            bgz, wkt_polygon], stderr=subprocess.STDOUT)
                        print('> 11 (Msk): processed ' + ifg + '_cc_ML_TC_msk in ' +
                              pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '10'], stderr=subprocess.STDOUT)
                    print('> 11 (Msk): deleted step 10 (TC) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

                # ### ONLY USE IF NOT ENOUGH RAM FOR STEP 12 ###
                # 1. make new wkt polygons for subset regions
                # 2. make new working directories (e.g. sum2021_A1)
                # 3. copy zip file names into new directories
                # 4. comment out normal step 11 (Msk)
                # 5. run steps 11-14 using e.g. ads=['A1','A2']
                # 6. mosaic geotiffs to re-merge image
                # ###
                # print('> 11-2 (Sub): splitting coherence map...')
                # # get wkt polygon
                # wkt_polygon = wktdir+'\BGZ'+bgz+'_wkt_'+ad[1]+'.txt'
                # # loop through ifgs
                # for ifg in ifgs:
                #    # get input masked cc
                #    io_mrg_ml_tc_msk = wdir[0:-1]+'\\'+ifg+'_cc_ML_TC_msk'
                #    # get output subset cc
                #    io_mrg_ml_tc_sub = wdir+'\\'+ifg+'_cc_ML_TC_msk'
                #    # if have cc
                #    if os.path.exists(io_mrg_ml_tc_msk+'.dim'):
                #        # subset
                #        pipe_out = subprocess.check_output(['python','snappy_bgz_11-2_sub.py',io_mrg_ml_tc_msk,
                #                                            io_mrg_ml_tc_sub,wkt_polygon],stderr=subprocess.STDOUT)
                #        print('> 11-2 (Sub): processed '+ifg+'_cc_ML_TC_msk in '+pipe_out.split()[-1])

            # --------------------- step 12 --------------------- #

            if 12 in run_steps:
                print('> 12 (Stack): stacking coherence images...')
                # stack
                pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_12_stack.py', wdir, str(min_btemp),
                                                    str(max_btemp), io_base + '_stack'], stderr=subprocess.STDOUT)
                print('> 12 (Stack): processed bgz' + bgz + '_' + season + '_' + ad + '_stack in ' +
                      pipe_out.decode(sys.stdout.encoding).split()[-1])
                # # clean previous step
                # if step_clean:
                #     pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir,
                #                                         str(min_btemp), str(max_btemp), io_base, '11'],
                #                                        stderr=subprocess.STDOUT)
                #     print('> 12 (Stack): deleted step 11 (Msk) files in ' +
                #           pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 13 --------------------- #

            if 13 in run_steps:
                print('> 13 (BM): calculating average temporal coherence...')
                # loop through polarisations
                for pol in ['VH', 'VV']:
                    # band maths
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_13_bm.py', io_base + '_stack', pol,
                                                        str(min_btemp), str(max_btemp)], stderr=subprocess.STDOUT)
                    print('> 13 (BM): processed bgz' + bgz + '_' + season + '_' + ad + '_stack_' + pol + '_BM in ' +
                          pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '12'], stderr=subprocess.STDOUT)
                    print('> 13 (BM): deleted step 12 (Stack) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- step 14 --------------------- #

            if 14 in run_steps:
                print('> 14 (Exp): exporting average coherence as geotiff...')
                # loop through polarisations
                for pol in ['VH', 'VV']:
                    # band maths and export
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_14_exp.py', io_base, pol],
                                                       stderr=subprocess.STDOUT)
                    print('> 14 (Exp): processed bgz' + bgz + '_' + season + '_' + ad + '_' + pol + '_cc in ' +
                          pipe_out.decode(sys.stdout.encoding).split()[-1])
                # clean previous step
                if step_clean:
                    pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_step_clean.py', wdir, str(min_btemp),
                                                        str(max_btemp), io_base, '13'], stderr=subprocess.STDOUT)
                    print('> 14 (Exp): deleted step 13 (BM) files in ' + pipe_out.decode(sys.stdout.encoding).split()[-1])

            # --------------------- end --------------------- #

            if end_clean and range(start, end):
                print('> END (Clean): removing processing files...')
                # delete data
                pipe_out = subprocess.check_output(['python', sdir + '\\snappy_bgz_END_clean.py', wdir, str(min_btemp),
                                                    str(max_btemp), str(start), str(end), io_base], stderr=subprocess.STDOUT)
                print('> END (Clean): deleted steps ' + ', '.join([str(s) for s in range(start, end)]) + ' in ' +
                      pipe_out.decode(sys.stdout.encoding).split()[-1])

            # finish processing timer
            toc = datetime.datetime.now().replace(microsecond=0)
            print('> DONE: ' + str(toc - tic))

        # catch errors
        except Exception as e:
            # print error message
            print('> ERROR: skipping remaining processing steps...')
            print(e)
            # finish processing timer
            toc = datetime.datetime.now().replace(microsecond=0)
            print('> DONE: ' + str(toc - tic))
            # continue script
            pass
# end
print('> finished')
