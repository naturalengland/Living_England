"""-------------------------------------------------------------
    Script Name:   snappy_bgz_00_download.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 00b
    Created By:    Chris Moore
    Date Created:  26.07.22
    Date Modified: 24.09.24

    Notes: Separated out step 0 of coherence code to be run separately. 
    Needs to use asf_download conda environment.
    
-------------------------------------------------------------"""

# step 0 (Zip): Download zip files from ASF urls

import sys
import os
import gc
import datetime
import subprocess

# biogeographic zone
# bgzs = ['01']  # '01'-'14'
bgzs = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14']  # '01'-'14' # doing 02, A & D - if finishes, remove 02
# bgzs_ads = [['07', ['D']], ['10', ['D2']], ['12', ['A2', 'D2']], ['13', ['D2']]]  # tmp
# season
season = 'sum2021'  # 'sprYYYY'/'sumYYYY'/'autYYYY'
# asc/desc
ads = ['A', 'D'],  # 'A'/'D'
# start & end steps
start = 0
end = 0
# clean step outputs? (step 11 excluded for qc)
step_clean = True  # True/False
end_clean = True  # True/False

# minimum and maximum temporal baseline to process
min_btemp = 12  # days
max_btemp = 24  # days
# master directory containing processing and data folders
mdir = 'D:\\S1_Coherence\\LE2021\\data'
# directory containing running scripts
sdir = 'D:\\le_code\\LivingEngland-Private\\1-Workflow-Scripts\\2_S1_Coherence'
# directory containing BGZ WKT polygons
wktdir = 'D:\\S1_Coherence\\LE2021\\BGZ_Buffer10km_WGS84_Points'
# directory containing BGZ ESRI shapefiles
shpdir = 'D:\\S1_Coherence\\LE2021\\BGZ_Buffer1km_WGS84_Shp'


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
                print("url_list = ", url_list)

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
                    if pipe_out.decode(sys.stdout.encoding).split()[-2] == "passed":
                        s1_dates_flag.append(True)
                        print(f"number of dates passed: {len(s1_dates_flag)}")
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
print('> downloaded qc urls')
