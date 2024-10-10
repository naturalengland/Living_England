"""-------------------------------------------------------------
    Script Name:   test_s1cc_weather.py
    Description:   Tests sensitivity of S1 coherence to recent rainfall, snowfall & wind
    Created By:    Chris Moore
    Date Created:  14.02.23
    Date Modified: 01.06.23
-------------------------------------------------------------"""

import os
import string
import time
import datetime
import math
import netCDF4
import warnings
import numpy as np
import pandas as pd
import rasterio as rs
# import geopandas as gpd
# import xarray as xr
# import rioxarray as rxr
import matplotlib.path as mplp
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# import matplotlib.units as munits
# from ea_rainfall_api import get_ea_rainfall

# biogeographic zone
bgz = "01"  # '01'-'14'
# season
season = "aut2022"  # sssYYYY

print("test_s1cc_weather\n> loading script params...")
tic = time.time()
# suppress warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
# set working directory
wdir = "C:\\Users\\Chris Moore\\Documents\\Living_England\\SAR_rainfall\\"

# get s1 acquisition dates/times
s1_dir = "D:\\S1\\bgz" + bgz + "\\"
s1_dates_times = []
s1_ifgs_all = []
for ad in ["A"]:
    with open(s1_dir + "bgz" + bgz + "_" + season + "_" + ad + "_urls.txt", "r") as f:
        urls = f.readlines()
    for u in range(len(urls)):
        # get unique dates
        if urls[u][56:64] != urls[u - 1][56:64]:
            s1_dates_times.append([ad, urls[u][56:64], urls[u][65:71]])
            # get ifg pairs
            if u > 0:
                s1_ifgs_all.append([ad, urls[u][56:64] + "_" + urls[u - 1][56:64]])

# # read bgz shapefile
# bgz_shp = "D:\\BioGeographicZones\\BGZ_Buffer1km_WGS84_Shp\\" + str(int(bgz)) + ".shp"
# bgz_gdf = gpd.read_file(bgz_shp)
# read bgz polygon points
bgz_csv = "D:\\BioGeographicZones\\BGZ_Buffer1km_WGS84_Points\\BGZ" + str(int(bgz)).zfill(2) + ".csv"
bgz_df = pd.read_csv(bgz_csv)
bgz_points = [tuple(row) for row in bgz_df[['POINT_X', 'POINT_Y']].to_numpy()]

try:
    # load s1 coherence geotiff data
    print("> loading S1 coherence GeoTIFFs...")
    # initialise arrays
    s1_ifgs_exist = []
    s1_stats = []
    # loop through ifgs
    for ifg in s1_ifgs_all:
        fname = s1_dir + season + "_" + ifg[0] + "\\" + ifg[1] + "_cc_ML_TC_msk.tif"
        if os.path.exists(fname):
            # read geotiff
            with rs.open(fname, "r") as cc:
                arr = cc.read()
            # set 0 to nan
            arr[arr == 0] = np.nan
            # store ifg name and dates/times
            d1i = [s1_dates_times.index(row) for row in s1_dates_times if ifg[1][:8] in row][0]
            d2i = [s1_dates_times.index(row) for row in s1_dates_times if ifg[1][9:] in row][0]
            s1_ifgs_exist.append(ifg + [s1_dates_times[d1i][1] + s1_dates_times[d1i][2],
                                        s1_dates_times[d2i][1] + s1_dates_times[d2i][2]])
            # calculate stats for vh and vv
            s1_stats.append([np.nanmedian(arr[0]), np.nanmedian(arr[1]), np.nanmean(arr[0]), np.nanmean(arr[1]),
                             np.nanstd(arr[0]), np.nanstd(arr[1])])
    # convert ifg info to df and format datetime to the nearest hour
    s1_df = pd.DataFrame(s1_ifgs_exist, columns=["ad", "ifg", "dt1", "dt2"])
    s1_df["dt1"] = pd.to_datetime(s1_df["dt1"], format="%Y%m%d%H%M%S").dt.round("H")
    s1_df["dt2"] = pd.to_datetime(s1_df["dt2"], format="%Y%m%d%H%M%S").dt.round("H")
    # convert stats to df and merge with ifg info df
    s1_stats = pd.DataFrame(s1_stats, columns=["vh_med", "vv_med", "vh_mean", "vv_mean", "vh_std", "vv_std"])
    s1_df = pd.concat([s1_df, s1_stats], axis=1)

    # load gpm hdf5 data
    # print("> loading GPM HDF5...")
    # gpm_dir = "D:\\Rainfall\\NASA_GPM_IMERG_Late_L3_HalfHourly_v6\\" + season + "\\"

    # # load era5 grib data
    # print("> loading ERA5 GRIB...")
    # era5_dir = "D:\\Rainfall\\ESA_ERA5\\" + season + "\\"
    # fname = "adaptor.mars.internal-1674813233.52303-4408-13-21baff49-5052-4044-ae0a-1dd3039aac6c.grib"
    # key_dicts = [{"filter_by_keys": {"shortName": "tp"}},
    #              {"filter_by_keys": {"shortName": "sf"}}]
    # era5_bgz = []
    # for key_dict in key_dicts:
    #     grib = xr.open_dataset(era5_dir + fname, engine="cfgrib", backend_kwargs=key_dict)
    #     era5_bgz.append(grib)

    # load era5 netcdf data
    print("> loading ERA5 NetCDF...")
    era5_dir = "D:\\Rainfall\\ESA_ERA5\\" + season + "\\"
    if season == "spr2021":
        # fname = "adaptor.mars.internal-1678786278.4258015-5031-8-de93e1d7-712f-4d9d-838c-09d0820b8b3d.nc"
        fname = "adaptor.mars.internal-1683714566.2010736-30308-18-a8f8c272-22e7-4eec-918c-d5654bce8831.nc"
    else:
        fname = "bgz" + bgz + "_" + season + ".nc"
    era5_nc = netCDF4.Dataset(era5_dir + fname)
    # create mask
    mpath = mplp.Path(bgz_points)
    era5_lon_grd, era5_lat_grd = np.meshgrid(era5_nc['longitude'][:].data, era5_nc['latitude'][:].data)
    mpoints = np.array((era5_lon_grd.flatten(), era5_lat_grd.flatten())).T
    mask = mpath.contains_points(mpoints).reshape(era5_lon_grd.shape)
    # convert time into datetime
    era5_time = pd.to_datetime(era5_nc['time'][:].data, unit='h', origin=pd.Timestamp('1900-01-01'))
    # get variable names
    era5_vars = []
    for var in era5_nc.variables.values():
        if var.name not in list(era5_nc.dimensions.keys()):
            era5_vars.append(var.name)
    # read data
    era5_med = np.zeros((len(era5_time), len(era5_vars)))
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
            # calculate stats
            era5_med[t, v] = np.nanmedian(var_arr_msk)
            era5_mean[t, v] = np.nanmean(var_arr_msk)
            era5_std[t, v] = np.nanstd(var_arr_msk)
    # calculate total wind speed from u and v
    if "u10" in era5_vars and "v10" in era5_vars:
        # new variable name
        era5_vars.append("w10")
        # calculate for mean and median: w = (u^2 + v^2)^0.5
        era5_med = np.hstack((era5_med, np.atleast_2d(np.sqrt(np.add(
            np.square(era5_med[:, era5_vars.index("u10")]), np.square(era5_med[:, era5_vars.index("v10")])))).T))
        era5_mean = np.hstack((era5_mean, np.atleast_2d(np.sqrt(np.add(
            np.square(era5_mean[:, era5_vars.index("u10")]), np.square(era5_mean[:, era5_vars.index("v10")])))).T))
        # propagate error for std: ew = (((eu^2 * u^2) + (ev^2 * v^2)) / (u^2 + v^2))^0.5
        era5_std = np.hstack((era5_std, np.atleast_2d(np.sqrt(np.divide(
            np.add(np.multiply(np.square(era5_std[:, era5_vars.index("u10")]),
                               np.square(era5_mean[:, era5_vars.index("u10")])),
                   np.multiply(np.square(era5_std[:, era5_vars.index("v10")]),
                               np.square(era5_mean[:, era5_vars.index("v10")]))),
            np.add(np.square(era5_mean[:, era5_vars.index("u10")]), np.square(era5_mean[:, era5_vars.index("v10")]))
        ))).T))
    # convert to df
    era5_df = pd.DataFrame(np.hstack((era5_med, era5_mean, era5_std)),
                           columns=[var + "_med" for var in era5_vars] + [var + "_mean" for var in era5_vars] +
                                   [var + "_std" for var in era5_vars]).set_index(era5_time)

    # get weather stats for start and end ifg dates/times
    print("> extracting weather stats for ifg dates...")
    tp_steps = [0, 1, 3, 6, 12]  # days
    s1_era5_w10 = np.full((len(s1_df), 6), np.nan)
    s1_era5_tp = np.full((len(tp_steps), len(s1_df), 6), np.nan)
    s1_era5_sf = np.full((len(tp_steps), len(s1_df), 6), np.nan)
    # loop through ifgs
    for row in s1_df.iterrows():
        if row[1][2] in era5_df.index and row[1][3] in era5_df.index:
            # wind at time of pass
            s1_era5_w10[row[0], :] = [era5_df.loc[row[1][2]]["w10_med"], era5_df.loc[row[1][2]]["w10_mean"],
                                      era5_df.loc[row[1][2]]["w10_std"], era5_df.loc[row[1][3]]["w10_med"],
                                      era5_df.loc[row[1][3]]["w10_mean"], era5_df.loc[row[1][3]]["w10_std"]]
            # rainfall at time of pass
            s1_era5_tp[0, row[0], :] = [era5_df.loc[row[1][2]]["tp_med"], era5_df.loc[row[1][2]]["tp_mean"],
                                        era5_df.loc[row[1][2]]["tp_std"], era5_df.loc[row[1][3]]["tp_med"],
                                        era5_df.loc[row[1][3]]["tp_mean"], era5_df.loc[row[1][3]]["tp_std"]]
            # snowfall at time of pass
            s1_era5_sf[0, row[0], :] = [era5_df.loc[row[1][2]]["sf_med"], era5_df.loc[row[1][2]]["sf_mean"],
                                        era5_df.loc[row[1][2]]["sf_std"], era5_df.loc[row[1][3]]["sf_med"],
                                        era5_df.loc[row[1][3]]["sf_mean"], era5_df.loc[row[1][3]]["sf_std"]]
            # total rainfall in time steps before each pass
            for i in range(1, len(tp_steps)):
                # get sum precipitation in step: tp = t1 + t2 + ... + tN
                d1_med = era5_df.loc[row[1][2] - pd.Timedelta(tp_steps[i], unit="days"):
                                     row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_med"].sum()
                d1_mean = era5_df.loc[row[1][2] - pd.Timedelta(tp_steps[i], unit="days"):
                                      row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_mean"].sum()
                d2_med = era5_df.loc[row[1][3] - pd.Timedelta(tp_steps[i], unit="days"):
                                     row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_med"].sum()
                d2_mean = era5_df.loc[row[1][3] - pd.Timedelta(tp_steps[i], unit="days"):
                                      row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_mean"].sum()
                # propagate errors in step: etp = (et1^2 + et2^2 + ... + etN^2)^0.5
                d1_std = np.sqrt(np.sum(np.square(
                    era5_df.loc[row[1][2] + pd.Timedelta(tp_steps[i], unit="days"):
                                row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_std"].to_numpy())))
                d2_std = np.sqrt(np.sum(np.square(
                    era5_df.loc[row[1][3] + pd.Timedelta(tp_steps[i], unit="days"):
                                row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["tp_std"].to_numpy())))
                # store
                s1_era5_tp[i, row[0], :] = [d1_med, d1_mean, d1_std, d2_med, d2_mean, d2_std]
            # total snowfall in time steps before each pass
            for i in range(1, len(tp_steps)):
                # get sum precipitation in step: tp = t1 + t2 + ... + tN
                d1_med = era5_df.loc[row[1][2] - pd.Timedelta(tp_steps[i], unit="days"):
                                     row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_med"].sum()
                d1_mean = era5_df.loc[row[1][2] - pd.Timedelta(tp_steps[i], unit="days"):
                                      row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_mean"].sum()
                d2_med = era5_df.loc[row[1][3] - pd.Timedelta(tp_steps[i], unit="days"):
                                     row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_med"].sum()
                d2_mean = era5_df.loc[row[1][3] - pd.Timedelta(tp_steps[i], unit="days"):
                                      row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_mean"].sum()
                # propagate errors in step: esf = (es1^2 + es2^2 + ... + esN^2)^0.5
                d1_std = np.sqrt(np.sum(np.square(
                    era5_df.loc[row[1][2] + pd.Timedelta(tp_steps[i], unit="days"):
                                row[1][2] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_std"].to_numpy())))
                d2_std = np.sqrt(np.sum(np.square(
                    era5_df.loc[row[1][3] + pd.Timedelta(tp_steps[i], unit="days"):
                                row[1][3] - pd.Timedelta(tp_steps[i-1], unit="days")]["sf_std"].to_numpy())))
                # store
                s1_era5_sf[i, row[0], :] = [d1_med, d1_mean, d1_std, d2_med, d2_mean, d2_std]
    # convert stats to df
    w10_cols = ["dt1_w10_med", "dt1_w10_mean", "dt1_w10_std", "dt2_w10_med", "dt2_w10_mean", "dt2_w10_std"]
    tp_cols = ["dt1_tp_med", "dt1_tp_mean", "dt1_tp_std", "dt2_tp_med", "dt2_tp_mean", "dt2_tp_std"]
    sf_cols = ["dt1_sf_med", "dt1_sf_mean", "dt1_sf_std", "dt2_sf_med", "dt2_sf_mean", "dt2_sf_std"]
    s1_df = pd.concat([s1_df, pd.DataFrame(s1_era5_w10, columns=w10_cols)], axis=1)
    for i in range(len(tp_steps)):
        tp_cols_step = [c[:6] + str(tp_steps[i]).zfill(2) + c[6:] for c in tp_cols]
        sf_cols_step = [c[:6] + str(tp_steps[i]).zfill(2) + c[6:] for c in sf_cols]
        s1_df = pd.concat([s1_df, pd.DataFrame(s1_era5_tp[i], columns=tp_cols_step),
                           pd.DataFrame(s1_era5_sf[i], columns=sf_cols_step)], axis=1)

    # remove nan rows from df
    s1_df.dropna(inplace=True)
    s1_df.reset_index(drop=True, inplace=True)
    # add index names
    letters = list(string.ascii_uppercase)
    letters = letters[:len(s1_df)]
    letters.reverse()
    s1_df = s1_df.set_index(pd.Index(letters))

    # find thresholds for wind, rainfall & snowfall
    print("> testing for weather thresholds...")
    # wind (m/s)
    # get sorted unique list of values
    thresholds_w10 = s1_df["dt1_w10_mean"].values.tolist()
    thresholds_w10.append(s1_df["dt2_w10_mean"].values[0])
    thresholds_w10 = sorted(list(set(thresholds_w10)))
    # initialise
    thresholds_w10_cc = np.full((len(thresholds_w10), 3), np.nan)
    n = 0
    # loop through thresholds
    for t in thresholds_w10:
        s1_df_flag = s1_df[(s1_df["dt1_w10_mean"] <= t) & (s1_df["dt2_w10_mean"] <= t)]
        thresholds_w10_cc[n, :] = [t, s1_df_flag["vh_mean"].mean(), s1_df_flag["vv_mean"].mean()]
        n += 1
    # rainfall (mm/day)
    # get sorted unique list of values
    thresholds_tp = []
    for i in range(1, len(tp_steps)-1):
        step = str(tp_steps[i]).zfill(2)
        t_diff = tp_steps[i] - tp_steps[i-1]
        thresholds_tp.extend(np.divide(s1_df["dt1_tp" + step + "_mean"].values.tolist(), t_diff).tolist())
        thresholds_tp.append(s1_df["dt2_tp" + step + "_mean"].values[0] / t_diff)
    thresholds_tp = sorted(list(set(thresholds_tp)))
    # initialise
    thresholds_tp_cc = np.full((len(thresholds_tp), 3), np.nan)
    n = 0
    # loop through thresholds
    for t in thresholds_tp:
        s1_df_flag = s1_df
        for i in range(1, len(tp_steps)-1):
            step = str(tp_steps[i]).zfill(2)
            t_diff = tp_steps[i] - tp_steps[i-1]
            s1_df_flag = s1_df_flag[(s1_df_flag["dt1_tp" + step + "_mean"] <= (t*t_diff)) &
                                    (s1_df_flag["dt2_tp" + step + "_mean"] <= (t*t_diff))]
        thresholds_tp_cc[n, :] = [t, s1_df_flag["vh_mean"].mean(), s1_df_flag["vv_mean"].mean()]
        n += 1
    # snowfall (mm/day)
    # get sorted unique list of values
    thresholds_sf = []
    for i in range(1, len(tp_steps)-1):
        step = str(tp_steps[i]).zfill(2)
        t_diff = tp_steps[i] - tp_steps[i-1]
        thresholds_sf.extend(np.divide(s1_df["dt1_sf" + step + "_mean"].values.tolist(), t_diff).tolist())
        thresholds_sf.append(s1_df["dt2_sf" + step + "_mean"].values[0] / t_diff)
    thresholds_sf = sorted(list(set(thresholds_sf)))
    # initialise
    thresholds_sf_cc = np.full((len(thresholds_sf), 3), np.nan)
    n = 0
    # loop through thresholds
    for t in thresholds_sf:
        s1_df_flag = s1_df
        for i in range(1, len(tp_steps)-1):
            step = str(tp_steps[i]).zfill(2)
            t_diff = tp_steps[i] - tp_steps[i-1]
            s1_df_flag = s1_df_flag[(s1_df_flag["dt1_sf" + step + "_mean"] <= (t*t_diff)) &
                                    (s1_df_flag["dt2_sf" + step + "_mean"] <= (t*t_diff))]
        thresholds_sf_cc[n, :] = [t, s1_df_flag["vh_mean"].mean(), s1_df_flag["vv_mean"].mean()]
        n += 1

    # create plots
    print("> plotting data...")
    # adjust season name
    if season[:3] == "spr":
        season_name = "Spring"
    elif season[:3] == "sum":
        season_name = "Summer"
    elif season[:3] == "aut":
        season_name = "Autumn"
    else:
        season_name = ""

    # overview figure set up
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # plot atmospheric data by time
    h1, = ax2.plot(era5_df.index, era5_df["w10_mean"], "-", color="c", linewidth=0.5, alpha=0.5, label="Wind Speed")
    for i in range(len(era5_mean)):
        h2, = ax2.plot([era5_df.index[i], era5_df.index[i]], [0, era5_df["tp_mean"][i]], "-", color="b", linewidth=0.5,
                       label="Rainfall")
        h21, = ax2.plot([era5_df.index[i], era5_df.index[i]], [0, era5_df["sf_mean"][i]], "-", color="m", linewidth=0.5,
                        label="Snowfall")
    # plot coherence data by time
    s1_fill_dates = []
    s1_fill_vh_minus = []
    s1_fill_vh_plus = []
    s1_fill_vv_minus = []
    s1_fill_vv_plus = []
    for i in range(len(s1_df)):
        s1_fill_dates.extend([s1_df["dt2"][i], s1_df["dt1"][i]])
        s1_fill_vh_minus.extend([s1_df["vh_mean"][i] - s1_df["vh_std"][i], s1_df["vh_mean"][i] - s1_df["vh_std"][i]])
        s1_fill_vh_plus.extend([s1_df["vh_mean"][i] + s1_df["vh_std"][i], s1_df["vh_mean"][i] + s1_df["vh_std"][i]])
        s1_fill_vv_minus.extend([s1_df["vv_mean"][i] - s1_df["vv_std"][i], s1_df["vv_mean"][i] - s1_df["vv_std"][i]])
        s1_fill_vv_plus.extend([s1_df["vv_mean"][i] + s1_df["vv_std"][i], s1_df["vv_mean"][i] + s1_df["vv_std"][i]])
    ax1.fill_between(s1_fill_dates + list(reversed(s1_fill_dates)), s1_fill_vh_minus + list(reversed(s1_fill_vh_plus)),
                     color="r", linewidth=0.5, alpha=0.15)
    ax1.fill_between(s1_fill_dates + list(reversed(s1_fill_dates)), s1_fill_vv_minus + list(reversed(s1_fill_vv_plus)),
                     color="tab:orange", linewidth=0.5, alpha=0.15)
    for i in range(len(s1_df)):
        h3, = ax1.plot([s1_df["dt1"][i], s1_df["dt2"][i]], [s1_df["vh_mean"][i], s1_df["vh_mean"][i]], "-",
                       color="r", linewidth=2, label="VH Coherence")
        h4, = ax1.plot([s1_df["dt1"][i], s1_df["dt2"][i]], [s1_df["vv_mean"][i], s1_df["vv_mean"][i]], "-",
                       color="tab:orange", linewidth=2, label="VV Coherence")
        h5 = ax1.text(s1_df["dt1"][i] + (s1_df["dt2"][i] - s1_df["dt1"][i]) / 2,
                      (s1_df["vv_mean"][i] + s1_df["vh_mean"][i]) / 2, list(s1_df.index)[i],
                      horizontalalignment="center", verticalalignment="center")
    # axis params
    ax1.set_ylabel("Coherence")
    ax2.set_ylabel("Precipitation (mm) / Wind Speed (m/s)")
    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 10)
    # ax1.set_xlim((era5_df.index[0], era5_df.index[-1] + np.timedelta64(1, "h")))
    ax1.set_xlim((s1_df["dt1"][-1] - np.timedelta64(6, "D"), s1_df["dt2"][0] + np.timedelta64(6, "D")))
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    # ax1.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1, 11, 21]))
    ax1.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1]))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m"))
    # legend
    leg = plt.legend(handles=[h3, h4, h21, h2, h1], loc=9, ncols=5)
    for line in leg.get_lines():
        line.set_linewidth(2)
    # figure params
    fig.set_figwidth(9.5)
    fig.set_figheight(4.5)
    plt.grid(visible=True, which="major", alpha=0.5)
    plt.title("BGZ " + bgz + " - " + season_name + " " + season[3:])
    fig.tight_layout()
    # save
    plt.show()
    plt.savefig(wdir + "bgz" + bgz + "_" + season + "_vh_vv_tp_w10_sf.png")
    plt.close()

    # wind analysis figure set up
    fig, (ax1, ax2) = plt.subplots(1, 2)
    # plot wind for each ifg with coherence
    ax1.scatter(s1_df["dt1_w10_mean"], s1_df["dt2_w10_mean"], s=[(cc * 5) ** 5 for cc in s1_df["vv_mean"]], c="b")
    ax1.scatter(s1_df["dt1_w10_mean"], s1_df["dt2_w10_mean"], s=[(cc * 5) ** 5 for cc in s1_df["vh_mean"]], c="r")
    # plot coherence by sum of wind
    # sum_w10 = s1_df[["dt1_w10_mean", "dt2_w10_mean"]].sum(axis=1)
    max_w10 = s1_df[["dt1_w10_mean", "dt2_w10_mean"]].max(axis=1)
    z_vh = np.polyfit(max_w10, s1_df["vh_mean"], 1)
    z_vv = np.polyfit(max_w10, s1_df["vv_mean"], 1)
    p_vh = np.poly1d(z_vh)
    p_vv = np.poly1d(z_vv)
    ax2.plot(max_w10, p_vh(max_w10), "r-", linewidth=2, alpha=0.4)
    ax2.plot(max_w10, p_vv(max_w10), "b-", linewidth=2, alpha=0.4)
    h1 = ax2.errorbar(max_w10, s1_df["vh_mean"], yerr=s1_df["vh_std"], elinewidth=0.5, capsize=2, fmt="ro",
                      markersize=5, label="VH")
    h2 = ax2.errorbar(max_w10, s1_df["vv_mean"], yerr=s1_df["vv_std"], elinewidth=0.5, capsize=2, fmt="bo",
                      markersize=5, label="VV")
    # axis params
    ax1.grid(visible=True, which="major", alpha=0.5)
    ax2.grid(visible=True, which="major", alpha=0.5)
    ax1.set_ylim(0, math.ceil(max_w10.max()))
    ax2.set_ylim(0.2, 0.8)
    ax1.set_xlim(0, math.ceil(max_w10.max()))
    ax2.set_xlim(max_w10.min(), max_w10.max())
    ax1.set_ylabel("Wind Speed at Pass 2 (m/s)")
    ax2.set_ylabel("Coherence")
    ax1.set_xlabel("Wind Speed at Pass 1 (m/s)")
    ax2.set_xlabel("Maximum Wind Speed (m/s)")
    # add labels
    for i in range(len(s1_df)):
        h9 = ax1.text(s1_df["dt1_w10_mean"][i], s1_df["dt2_w10_mean"][i], list(s1_df.index)[i],
                      horizontalalignment="left", verticalalignment="bottom")
        h10 = ax2.text(max_w10[i], (s1_df["vh_mean"][i] + s1_df["vv_mean"][i]) / 2, list(s1_df.index)[i],
                       horizontalalignment="center", verticalalignment="center")
    # legend
    ax2.legend(handles=[h1, h2], loc=1)
    h3, = ax1.plot(-10, -10, "ro", markersize=5, label="VH")
    h4, = ax1.plot(-10, -10, "bo", markersize=5, label="VV")
    h5 = ax1.scatter(-10, -10, s=(0.6 * 5) ** 5, c="k", label="0.6")
    h6 = ax1.scatter(-10, -10, s=(0.5 * 5) ** 5, c="k", label="0.5")
    h7 = ax1.scatter(-10, -10, s=(0.4 * 5) ** 5, c="k", label="0.4")
    h8 = ax1.scatter(-10, -10, s=(0.3 * 5) ** 5, c="k", label="0.3")
    leg = ax1.legend(handles=[h5, h6, h7, h8], loc=1, title="Coherence")
    ax1.legend(handles=[h3, h4], loc=1, bbox_to_anchor=(0.73, 1))
    ax1.add_artist(leg)
    # figure params
    fig.set_figwidth(8.5)
    fig.set_figheight(4.5)
    fig.tight_layout()
    # save
    plt.show()
    plt.savefig(wdir + "bgz" + bgz + "_" + season + "_vh_vv_wind.png")
    plt.close()

    # loop through rainfall steps
    for i in range(len(tp_steps)):
        step = str(tp_steps[i]).zfill(2)
        if i == 0:
            axis_lab = "at"
        else:
            axis_lab = str(tp_steps[i - 1]) + "-" + str(tp_steps[i]) + " days before"
        # rainfall analysis figure set up
        fig, (ax1, ax2) = plt.subplots(1, 2)
        # plot wind for each ifg with coherence
        ax1.scatter(s1_df["dt1_tp" + step + "_mean"], s1_df["dt2_tp" + step + "_mean"],
                    s=[(cc * 5) ** 5 for cc in s1_df["vv_mean"]], c="b")
        ax1.scatter(s1_df["dt1_tp" + step + "_mean"], s1_df["dt2_tp" + step + "_mean"],
                    s=[(cc * 5) ** 5 for cc in s1_df["vh_mean"]], c="r")
        # plot coherence by sum of rainfall
        sum_tp = s1_df[["dt1_tp" + step + "_mean", "dt2_tp" + step + "_mean"]].sum(axis=1)
        max_tp = s1_df[["dt1_tp" + step + "_mean", "dt2_tp" + step + "_mean"]].max(axis=1)
        z_vh = np.polyfit(sum_tp, s1_df["vh_mean"], 1)
        z_vv = np.polyfit(sum_tp, s1_df["vv_mean"], 1)
        p_vh = np.poly1d(z_vh)
        p_vv = np.poly1d(z_vv)
        ax2.plot(sum_tp, p_vh(sum_tp), "r-", linewidth=2, alpha=0.4)
        ax2.plot(sum_tp, p_vv(sum_tp), "b-", linewidth=2, alpha=0.4)
        h1 = ax2.errorbar(sum_tp, s1_df["vh_mean"], yerr=s1_df["vh_std"], elinewidth=0.5, capsize=2, fmt="ro",
                          markersize=5, label="VH")
        h2 = ax2.errorbar(sum_tp, s1_df["vv_mean"], yerr=s1_df["vv_std"], elinewidth=0.5, capsize=2, fmt="bo",
                          markersize=5, label="VV")
        # axis params
        ax1.grid(visible=True, which="major", alpha=0.5)
        ax2.grid(visible=True, which="major", alpha=0.5)
        ax1.set_ylim(0, math.ceil(max_tp.max()))
        ax2.set_ylim(0.2, 0.8)
        ax1.set_xlim(0, math.ceil(max_tp.max()))
        ax2.set_xlim(sum_tp.min(), sum_tp.max())
        ax1.set_ylabel("Precipitation " + axis_lab + " Pass 2 (mm)")
        ax2.set_ylabel("Coherence")
        ax1.set_xlabel("Precipitation " + axis_lab + " Pass 1 (mm)")
        ax2.set_xlabel("Combined Precipitation (mm)")
        # add labels
        for j in range(len(s1_df)):
            h9 = ax1.text(s1_df["dt1_tp" + step + "_mean"][j], s1_df["dt2_tp" + step + "_mean"][j],
                          list(s1_df.index)[j], horizontalalignment="left", verticalalignment="bottom")
            h10 = ax2.text(sum_tp[j], (s1_df["vh_mean"][j] + s1_df["vv_mean"][j]) / 2, list(s1_df.index)[j],
                           horizontalalignment="center", verticalalignment="center")
        # legend
        ax2.legend(handles=[h1, h2], loc=1)
        h3, = ax1.plot(-10, -10, "ro", markersize=5, label="VH")
        h4, = ax1.plot(-10, -10, "bo", markersize=5, label="VV")
        h5 = ax1.scatter(-10, -10, s=(0.6 * 5) ** 5, c="k", label="0.6")
        h6 = ax1.scatter(-10, -10, s=(0.5 * 5) ** 5, c="k", label="0.5")
        h7 = ax1.scatter(-10, -10, s=(0.4 * 5) ** 5, c="k", label="0.4")
        h8 = ax1.scatter(-10, -10, s=(0.3 * 5) ** 5, c="k", label="0.3")
        leg = ax1.legend(handles=[h5, h6, h7, h8], loc=1, title="Coherence")
        ax1.legend(handles=[h3, h4], loc=1, bbox_to_anchor=(0.73, 1))
        ax1.add_artist(leg)
        # figure params
        fig.set_figwidth(8.5)
        fig.set_figheight(4.5)
        fig.tight_layout()
        # save
        plt.show()
        plt.savefig(wdir + "bgz" + bgz + "_" + season + "_vh_vv_tp" + step + ".png")
        plt.close()

    # threshold figure set up
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    # plot wind and rainfall thresholds
    h1, = ax1.plot(thresholds_tp_cc[:, 0], thresholds_tp_cc[:, 1], "o", color="r", ms=5, mec="k", mew=0.5, label="VH")
    h2, = ax1.plot(thresholds_tp_cc[:, 0], thresholds_tp_cc[:, 2], "o", color="b", ms=5, mec="k", mew=0.5, label="VV")
    h3, = ax2.plot(thresholds_w10_cc[:, 0], thresholds_w10_cc[:, 1], "o", color="r", ms=5, mec="k", mew=0.5, label="VH")
    h4, = ax2.plot(thresholds_w10_cc[:, 0], thresholds_w10_cc[:, 2], "o", color="b", ms=5, mec="k", mew=0.5, label="VV")
    h5, = ax3.plot(thresholds_sf_cc[:, 0], thresholds_sf_cc[:, 1], "o", color="r", ms=5, mec="k", mew=0.5, label="VH")
    h6, = ax3.plot(thresholds_sf_cc[:, 0], thresholds_sf_cc[:, 2], "o", color="b", ms=5, mec="k", mew=0.5, label="VV")
    # axis params
    ax1.grid(visible=True, which="major", alpha=0.5)
    ax2.grid(visible=True, which="major", alpha=0.5)
    ax3.grid(visible=True, which="major", alpha=0.5)
    ax1.set_xlim(0, math.ceil(thresholds_tp_cc[:, 0][-1]))
    ax2.set_xlim(0, math.ceil(thresholds_w10_cc[:, 0][-1]))
    ax3.set_xlim(0, math.ceil(thresholds_sf_cc[:, 0][-1]))
    ax1.set_ylim(min([np.nanmin(thresholds_tp_cc[:, 1:3]), np.nanmin(thresholds_w10_cc[:, 1:3]),
                      np.nanmin(thresholds_sf_cc[:, 1:3])]),
                 max([np.nanmax(thresholds_tp_cc[:, 1:3]), np.nanmax(thresholds_w10_cc[:, 1:3]),
                      np.nanmax(thresholds_sf_cc[:, 1:3])]))
    ax2.set_ylim(min([np.nanmin(thresholds_tp_cc[:, 1:3]), np.nanmin(thresholds_w10_cc[:, 1:3]),
                      np.nanmin(thresholds_sf_cc[:, 1:3])]),
                 max([np.nanmax(thresholds_tp_cc[:, 1:3]), np.nanmax(thresholds_w10_cc[:, 1:3]),
                      np.nanmax(thresholds_sf_cc[:, 1:3])]))
    ax3.set_ylim(min([np.nanmin(thresholds_tp_cc[:, 1:3]), np.nanmin(thresholds_w10_cc[:, 1:3]),
                      np.nanmin(thresholds_sf_cc[:, 1:3])]),
                 max([np.nanmax(thresholds_tp_cc[:, 1:3]), np.nanmax(thresholds_w10_cc[:, 1:3]),
                      np.nanmax(thresholds_sf_cc[:, 1:3])]))
    ax1.set_xlabel("Rainfall (mm/day)")
    ax2.set_xlabel("Wind Speed (m/s)")
    ax3.set_xlabel("Snowfall (mm/day)")
    ax1.set_ylabel("Average Coherence")
    ax2.set_ylabel("Average Coherence")
    ax3.set_ylabel("Average Coherence")
    # legend
    ax1.legend(handles=[h1, h2], loc=1)
    ax2.legend(handles=[h3, h4], loc=1)
    ax3.legend(handles=[h3, h4], loc=1)
    # figure params
    fig.set_figwidth(12.5)
    fig.set_figheight(4.5)
    fig.tight_layout()
    # save
    plt.show()
    plt.savefig(wdir + "bgz" + bgz + "_" + season + "_vh_vv_tp_w10_sf_thresholds.png")
    plt.close()

except Exception as e:
    print(e)
# end
print("> done")
toc = time.time()
print(datetime.timedelta(seconds=round(toc - tic)))
