"""-------------------------------------------------------------
    Script Name:   s1_bgz_weather_qc.py
    Description:   Check S1 acquisition dates with weather quality check (step 0 QC)
    Created By:    Chris Moore
    Date Created:  11.08.23
    Date Modified: 15.08.23
-------------------------------------------------------------"""

import sys
import os
import datetime
import subprocess
import pandas as pd

# start timer
now = datetime.datetime.now()
tic = now.replace(microsecond=0)
print(now.strftime("%d/%m %H:%M") + ' s1_bgz_weather_qc')

# seasons
seasons = ["sum2023"]  # "aut2022"/"spr2023"/"sum2023"
# biogeographic zones
bgzs = [str(b).zfill(2) for b in range(1, 15)]
# bgzs = ["01"]
# ascending/descending
ads = ["A", "D", "A2", "D2"]

# output directory
out_dir = "C:\\Users\\Chris Moore\\Documents\\Living_England\\SAR_rainfall"
# master directory containing processing and data folders
mdir = "D:\\S1"
# directory containing qc script
sdir = "C:\\Users\\Chris Moore\\Documents\\Living_England\\S1_coherence"
# directory containing BGZ WKT polygons
wktdir = "D:\\BioGeographicZones\\BGZ_Buffer10km_WGS84_Points"

# start processing
try:
    # read file containing S1 tracks for each bgz
    tracks = pd.read_csv(mdir + "\\bgz_tracks.txt", sep="\t")
    # initialise store
    cols = ["bgz", "track", "date", "qc"]
    out_df = pd.DataFrame(columns=cols)

    # loop through bgzs, seasons & orbits
    for season in seasons:
        for bgz in bgzs:
            for ad in ads:
                # check if have url list
                url_file = mdir + '\\bgz' + bgz + '\\bgz' + bgz + '_' + season + '_' + ad + '_urls.txt'
                if os.path.exists(url_file):
                    print("> checking acquisition dates for " + season + ", bgz" + bgz + ", " + ad)

                    # set processing directory
                    wdir = mdir + '\\bgz' + bgz + '\\' + season + '_' + ad
                    if not os.path.exists(wdir):
                        os.mkdir(wdir)
                    # create api download directory
                    if not os.path.exists(wdir + "\\era5_api"):
                        os.mkdir(wdir + "\\era5_api")

                    # get biogeographic zone polygon
                    csv_polygon = wktdir + '\\BGZ' + bgz + '.csv'

                    # read urls
                    with open(url_file, 'r') as f:
                        url_list = f.read().splitlines()
                    # get unique acquisition dates with times (YYYYMMDDhhmmss)
                    s1_dates = []
                    for u in range(len(url_list)):
                        if url_list[u][56:64] != url_list[u - 1][56:64]:
                            s1_dates.append(url_list[u][56:64] + url_list[u][65:71])

                    # get era-5 weather data for 6 days prior to date
                    for s1_date in s1_dates:
                        # download era-5 data and apply thresholds
                        pipe_out = subprocess.check_output(
                            ['python', sdir + '\\snappy_bgz_00_qc.py', s1_date, csv_polygon, bgz,
                             wdir + "\\era5_api"], stderr=subprocess.STDOUT)
                        print('> ' + s1_date[:8] + ' check ' + pipe_out.decode(sys.stdout.encoding).split()[-3] + ' in '
                              + pipe_out.decode(sys.stdout.encoding).split()[-2])
                        # store flagged dates
                        if pipe_out.decode(sys.stdout.encoding).split()[-3] == "passed":
                            s1_date_flag = True
                        else:
                            s1_date_flag = False
                        # store data
                        track = str(int(tracks[ad].loc[tracks["BGZ"] == int(bgz)].values.tolist()[0])).zfill(3)
                        out_df.loc[len(out_df)] = [bgz, track, s1_date[:8], s1_date_flag]

        # export dataframe
        out_df.to_csv(out_dir + "\\s1_bgz_weather_qc_" + season + ".csv", columns=cols, index=False)

# catch errors
except Exception as e:
    print(e)
# end
toc = datetime.datetime.now().replace(microsecond=0)
print('> done in ' + str(toc - tic))
