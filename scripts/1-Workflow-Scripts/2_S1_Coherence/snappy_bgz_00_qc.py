"""-------------------------------------------------------------
    Script Name:   snappy_bgz_00_qc.py
    Description:   SNAP processing chain for VV/VH S1 coherence data - step 00a
    Created By:    Chris Moore
    Date Created:  26.07.22
    Date Modified: 15.08.23
-------------------------------------------------------------"""

# step 0 (QC): Download zip files from ASF urls

import os
import time
import datetime
import sys
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.path as mplp
import cdsapi

# start timer
tic = time.time()
# get inputs
s1_date = sys.argv[1]
bgz_csv = sys.argv[2]
bgz = sys.argv[3]
era5_dir = sys.argv[4]

# thresholds
era5_flag = np.array([
    10,  # wind speed at acquisition date (m/s)
    15,  # rainfall 0-1 days before date (mm/day)
    15,  # rainfall 1-3 days before date (mm/day)
    15,  # rainfall 3-6 days before date (mm/day)
    1,   # snowfall 0-1 days before date (mm/day)
    1,   # snowfall 1-3 days before date (mm/day)
    1    # snowfall 3-6 days before date (mm/day)
])

# read bgz polygon points
bgz_df = pd.read_csv(bgz_csv)
bgz_points = [tuple(row) for row in bgz_df[['POINT_X', 'POINT_Y']].to_numpy()]

# convert YYYYMMDDhhmmss to datetime
date = pd.to_datetime(s1_date, format="%Y%m%d%H%M%S")
# get dates to download
years = []
months = []
days = []
# record all dates up to 6 days before acquisition
for d in range(7):
    # time difference
    t_delta = date - pd.Timedelta(days=d)
    # record date
    years.append(t_delta.year)
    months.append(t_delta.month)
    days.append(t_delta.day)
# remove duplicates, sort and format lists
years = [str(s).zfill(2) for s in set(years)]
months = [str(s).zfill(2) for s in set(months)]
days = [str(s).zfill(2) for s in set(days)]
# netcdf file name
fname = "bgz" + bgz + "_y" + "-".join(years) + "_m" + "-".join(months) + "_d" + "-".join(days) + ".nc"

# catch errors
try:
    # download data from cds api
    if not os.path.exists(era5_dir + "\\" + fname):
        # ERA5 data source
        data_source = "reanalysis-era5-single-levels"
        product = "reanalysis"
        # output format
        out_type = "netcdf"
        # data types to download
        variables = ["10m_u_component_of_wind", "10m_v_component_of_wind", "snowfall", "total_precipitation"]
        # hourly time steps
        times = ["00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                 "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"]
        # area of interest for bgzs [y max, x min, y min, x max]
        bbox_dict = {
            "01": [55.09998782, -3.702568339, 53.90933757, -1.840694843],
            "02": [55.84761397, -2.902215556, 54.3904869, -0.949713418],
            "03": [54.62663541, -3.342593206, 53.18742381, -1.420024592],
            "04": [54.62663192, -1.739141622, 52.90127161, -0.012829053],
            "05": [53.3780367, -3.183219618, 52.17337551, -1.395471542],
            "06": [53.90212085, -2.163717527, 52.38285334, -0.504578137],
            "07": [54.1134966, -0.83716391, 52.21966776, 0.574107403],
            "08": [53.02760598, -0.300906889, 51.68300339, 1.822551742],
            "09": [52.65929243, -1.742799086, 51.44637266, 0.216462399],
            "10": [51.97215252, -1.557572202, 50.73697559, 1.509343873],
            "11": [52.53294537, -3.250251265, 50.92797416, -1.088777339],
            "12": [51.97766502, -2.987131564, 50.47687975, 0.337686488],
            "13": [51.28280131, -5.803202263, 49.9160602, -2.803928282],
            "14": [50.01870611, -6.511325134, 49.82625956, -6.186329462]
        }
        # get api client
        c = cdsapi.Client()
        print("> 0 (QC): getting api client")
        # download data
        data_dict = {"product_type": product, "format": out_type, "variable": variables, "year": years, "month": months,
                     "day": days, "time": times, "area": bbox_dict[bgz]}
        print(f"> 0 (QC): retrieving {data_dict}")
        c.retrieve(data_source, data_dict, era5_dir + "\\" + fname)
        print(f"> 0 (QC): successful  {era5_dir}\\ {fname}")
        

    # import era5 netcdf
    era5_nc = netCDF4.Dataset(era5_dir + "\\" + fname)
    print(f"> 0 (QC): importing netcdf {era5_dir}\\ {fname}")
    # create mask
    mpath = mplp.Path(bgz_points)
    era5_lon_grd, era5_lat_grd = np.meshgrid(era5_nc['longitude'][:].data, era5_nc['latitude'][:].data)
    mpoints = np.array((era5_lon_grd.flatten(), era5_lat_grd.flatten())).T
    mask = mpath.contains_points(mpoints).reshape(era5_lon_grd.shape)
    print(f"> 0 (QC): creating mask")  
    # convert time into datetime
    era5_time = pd.to_datetime(era5_nc['valid_time'][:].data, unit='s')
    # get variable names
    # era5_vars = []
    # for var in era5_nc.variables.values():
    #     if var.name not in list(era5_nc.dimensions.keys()):
    #         era5_vars.append(var.name)
    era5_vars = ['u10', 'v10', 'sf', 'tp']
    # read data
    era5_mean = np.zeros((len(era5_time), len(era5_vars)))
    era5_std = np.zeros((len(era5_time), len(era5_vars)))
    for v in range(len(era5_vars)):
        # get variable data
        var_arr = era5_nc[era5_vars[v]][:].data
        # convert m to mm for tp and sf
        if era5_nc[era5_vars[v]].units == "m" or era5_nc[era5_vars[v]].units == "m of water equivalent":
            var_arr = np.multiply(var_arr, 1000)
        # get timestep data
        for t in range(len(era5_time)):
            var_arr_msk = var_arr[t]
            # apply mask
            var_arr_msk[~mask] = np.nan
            # calculate mean
            era5_mean[t, v] = np.nanmean(var_arr_msk)
    # calculate total wind speed from u and v
    print(f"> 0 (QC): calculating wind speed")
    if "u10" in era5_vars and "v10" in era5_vars:
        # new variable name
        era5_vars.append("w10")
        # calculate mean: w = (u^2 + v^2)^0.5
        era5_mean = np.hstack((era5_mean, np.atleast_2d(np.sqrt(np.add(
            np.square(era5_mean[:, era5_vars.index("u10")]), np.square(era5_mean[:, era5_vars.index("v10")])))).T))
    # convert to df
    print(f"> 0 (QC): converting to dataframe")
    era5_df = pd.DataFrame(era5_mean, columns=[var + "_mean" for var in era5_vars]).set_index(era5_time)

    # get stats for acquisition date
    print(f"> 0 (QC): getting stats for acquisition date")
    s1_stats = np.array([
        # wind at time of pass (m/s)
        era5_df.loc[date.round("h")]["w10_mean"],
        # rainfall in time steps (mm/day)
        era5_df.loc[date.round("h") - pd.Timedelta(days=1):date.round("h")]["tp_mean"].sum(),
        era5_df.loc[date.round("h") - pd.Timedelta(days=3):date.round("h") - pd.Timedelta(days=1)]["tp_mean"].sum() / 2,
        era5_df.loc[date.round("h") - pd.Timedelta(days=6):date.round("h") - pd.Timedelta(days=3)]["tp_mean"].sum() / 3,
        # snowfall in time steps (mm/day)
        era5_df.loc[date.round("h") - pd.Timedelta(days=1):date.round("h")]["sf_mean"].sum(),
        era5_df.loc[date.round("h") - pd.Timedelta(days=3):date.round("h") - pd.Timedelta(days=1)]["sf_mean"].sum() / 2,
        era5_df.loc[date.round("h") - pd.Timedelta(days=6):date.round("h") - pd.Timedelta(days=3)]["sf_mean"].sum() / 3
    ])
    print(f" threshold: {min(era5_flag > s1_stats)}")

    # output threshold check
    if not min(era5_flag > s1_stats):  # if any stat is more than threshold then is false (fails)
        print("failed")
    else:
        print("passed")

# catch errors
except Exception as e:
    print("Error: " + str(e))
    print("failed")

# output time
toc = time.time()
print(str(datetime.timedelta(seconds=round(toc - tic))))
